#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# nav_rinex2_to_df.py
"""
Parse a GPS RINEX 2.xx NAV file (.n / .16n) into a DataFrame-like list of dicts,
and provide utilities to compute satellite ECEF positions, azimuth/elevation from
a receiver, and the thin-shell mapping function M(E) used for VTEC.

Outputs (when run as a script):
  - <out_prefix>.gps_nav.csv : one row per ephemeris (PRN, toe, parameters)
  - (Optionally) a demo CSV (<out_prefix>.demo_positions.csv) if you pass
    --demo-time "YYYY-MM-DD HH:MM:SS" and a receiver ECEF via --rx "X,Y,Z"

Usage:
  python nav_rinex2_to_df.py OUCA0010.16n nav_2016-10-01 \      --demo-time "2016-10-01 00:30:00" --rx "5411106.3084,-747628.5974,3286931.3223"

"""
import sys, math, csv, re
from datetime import datetime, timedelta, timezone

MU  = 3.986005e14               # m^3/s^2 (WGS-84 Earth gravitational)
OMEGA_E_DOT = 7.2921151467e-5   # rad/s (Earth rotation rate)
Frel = -4.442807633e-10         # s/m^1/2 (relativistic correction factor)

def _to_float(s):
    s = s.replace('D','E').replace('d','e')
    try:
        return float(s)
    except:
        return math.nan

def _fields(line, widths):
    out = []
    i = 0
    for w in widths:
        seg = line[i:i+w]
        out.append(seg)
        i += w
    return out

def _parse_epoch(line):
    prn = line[0:2].strip()
    if prn and prn[0].isdigit():
        prn = 'G' + prn.zfill(2)
    elif line[0] in ('G','R','E','C'):
        prn = line[0] + line[1:3].strip().zfill(2)
    else:
        tok = re.findall(r'\S+', line[:3])
        if tok:
            t = tok[0]
            if t[0] in ('G','R','E','C'):
                prn = t[0] + t[1:].zfill(2)
            else:
                prn = 'G' + t.zfill(2)
        else:
            prn = 'G00'

    yy = int(line[3:5]); mo = int(line[6:8]); dd = int(line[9:11])
    hh = int(line[12:14]); mi = int(line[15:17]); ss = int(float(line[18:22]))
    year = 1900+yy if yy>=80 else 2000+yy
    t0 = datetime(year, mo, dd, hh, mi, ss, tzinfo=timezone.utc)

    af0 = _to_float(line[22:41])
    af1 = _to_float(line[41:60])
    af2 = _to_float(line[60:79])
    return prn, t0, af0, af1, af2

def parse_rinex2_nav(path):
    with open(path, 'r', errors='ignore') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines) and "END OF HEADER" not in lines[i]:
        i += 1
    i += 1

    ephs = []
    while i < len(lines):
        if i+7 >= len(lines):
            break
        l1 = lines[i].rstrip("\n"); l2 = lines[i+1].rstrip("\n")
        l3 = lines[i+2].rstrip("\n"); l4 = lines[i+3].rstrip("\n")
        l5 = lines[i+4].rstrip("\n"); l6 = lines[i+5].rstrip("\n")
        l7 = lines[i+6].rstrip("\n"); l8 = lines[i+7].rstrip("\n")
        try:
            prn, toc, af0, af1, af2 = _parse_epoch(l1)
        except Exception:
            i += 1
            continue

        def get4(l):
            return [_to_float(l[3+19*k:3+19*(k+1)]) for k in range(4)]

        IODE, Crs, d_n, M0     = get4(l2)
        Cuc,  e,   Cus, sqrtA  = get4(l3)
        Toe,  Cic, Omega0, Cis = get4(l4)
        i0,   Crc, omega, OMEGA_DOT = get4(l5)
        IDOT, L2_codes, week, L2P_flag = get4(l6)
        sv_accuracy, sv_health, Tgd, IODC = get4(l7)
        TxTime, fit_int, spare1, spare2 = get4(l8)

        ephs.append({
            "prn": prn, "toc": toc,
            "af0": af0, "af1": af1, "af2": af2,
            "IODE": IODE, "Crs": Crs, "d_n": d_n, "M0": M0,
            "Cuc": Cuc, "e": e, "Cus": Cus, "sqrtA": sqrtA,
            "Toe": Toe, "Cic": Cic, "Omega0": Omega0, "Cis": Cis,
            "i0": i0, "Crc": Crc, "omega": omega, "OMEGA_DOT": OMEGA_DOT,
            "IDOT": IDOT, "week": week, "sv_accuracy": sv_accuracy,
            "sv_health": sv_health, "Tgd": Tgd, "IODC": IODC,
            "TxTime": TxTime, "fit_int": fit_int
        })
        i += 8

    return ephs

def gpsweek(dt):
    epoch = datetime(1980,1,6, tzinfo=timezone.utc)
    delta = dt - epoch
    sow = delta.total_seconds()
    week = int(sow // 604800)
    tow  = sow - week*604800
    return week, tow

def _wrap_tk(tk):
    half = 302400.0
    if tk >  half: tk -= 2*half
    if tk < -half: tk += 2*half
    return tk

def sat_ecef_from_ephem(eph, t_utc):
    week_t, tow_t = gpsweek(t_utc)
    week_toe = int(eph["week"]) if not math.isnan(eph["week"]) else week_t
    toe = float(eph["Toe"])
    tk = (week_t - week_toe)*604800.0 + (tow_t - toe)
    tk = _wrap_tk(tk)

    A = (eph["sqrtA"])**2
    n0 = math.sqrt(MU / A**3)
    n  = n0 + eph["d_n"]
    Mk = eph["M0"] + n*tk

    Ek = Mk
    for _ in range(10):
        Ek_next = Mk + eph["e"]*math.sin(Ek)
        if abs(Ek_next - Ek) < 1e-12:
            Ek = Ek_next
            break
        Ek = Ek_next

    sin_vk = math.sqrt(1 - eph["e"]**2) * math.sin(Ek) / (1 - eph["e"]*math.cos(Ek))
    cos_vk = (math.cos(Ek) - eph["e"]) / (1 - eph["e"]*math.cos(Ek))
    vk = math.atan2(sin_vk, cos_vk)

    phik = vk + eph["omega"]
    duk = eph["Cus"]*math.sin(2*phik) + eph["Cuc"]*math.cos(2*phik)
    drk = eph["Crs"]*math.sin(2*phik) + eph["Crc"]*math.cos(2*phik)
    dik = eph["Cis"]*math.sin(2*phik) + eph["Cic"]*math.cos(2*phik)

    uk = phik + duk
    rk = A*(1 - eph["e"]*math.cos(Ek)) + drk
    ik = eph["i0"] + eph["IDOT"]*tk + dik

    x_orb = rk * math.cos(uk)
    y_orb = rk * math.sin(uk)

    Omega = eph["Omega0"] + (eph["OMEGA_DOT"] - OMEGA_E_DOT)*tk - OMEGA_E_DOT*eph["Toe"]

    cosO = math.cos(Omega); sinO = math.sin(Omega)
    cosi = math.cos(ik);    sini = math.sin(ik)

    x = x_orb*cosO - y_orb*cosi*sinO
    y = x_orb*sinO + y_orb*cosi*cosO
    z = y_orb*sini
    return x, y, z

def ecef_to_geodetic(x, y, z):
    a = 6378137.0
    f = 1/298.257223563
    e2 = f*(2-f)
    b = a*(1-f)
    ep2 = (a*a - b*b)/(b*b)

    r = math.hypot(x, y)
    E2 = a*a - b*b
    F = 54*b*b*z*z
    G = r*r + (1 - e2)*z*z - e2*E2
    c = (e2*e2*F*r*r)/(G*G*G)
    s = (1 + c + math.sqrt(c*c + 2*c))**(1/3)
    P = F/(3*(s + 1/s + 1)**2 * G*G)
    Q = math.sqrt(1 + 2*e2*e2*P)
    r0 = -(P*e2*r)/(1+Q) + math.sqrt(0.5*a*a*(1+1/Q) - (P*(1-e2)*z*z)/(Q*(1+Q)) - 0.5*P*r*r)
    U = math.sqrt((r - e2*r0)**2 + z*z)
    V = math.sqrt((r - e2*r0)**2 + (1 - e2)*z*z)
    Z0 = b*b*z/(a*V)
    h = U*(1 - b*b/(a*V))
    lat = math.atan2(z + ep2*Z0, r)
    lon = math.atan2(y, x)
    return lat, lon, h

def enu_from_ecef(rx_xyz, sat_xyz):
    x, y, z = rx_xyz
    lat, lon, _ = ecef_to_geodetic(x, y, z)
    slat, clat = math.sin(lat), math.cos(lat)
    slon, clon = math.sin(lon), math.cos(lon)
    R = [
        [-slon,          clon,         0],
        [-slat*clon, -slat*slon,  clat],
        [ clat*clon,  clat*slon,  slat],
    ]
    dx = sat_xyz[0]-x; dy = sat_xyz[1]-y; dz = sat_xyz[2]-z
    e = R[0][0]*dx + R[0][1]*dy + R[0][2]*dz
    n = R[1][0]*dx + R[1][1]*dy + R[1][2]*dz
    u = R[2][0]*dx + R[2][1]*dy + R[2][2]*dz
    return e, n, u

def az_el_from_ecef(rx_xyz, sat_xyz):
    e, n, u = enu_from_ecef(rx_xyz, sat_xyz)
    az = math.atan2(e, n) % (2*math.pi)
    r = math.sqrt(e*e + n*n + u*u)
    el = math.asin(u / r)
    return az, el

def vtec_mapping(el_rad, shell_height=350e3, Re=6378137.0):
    cosE = math.cos(el_rad)
    sin_zp = (Re/(Re + shell_height)) * cosE
    sin_zp = max(min(sin_zp, 1.0), -1.0)
    zp = math.asin(sin_zp)
    M = 1.0 / math.cos(zp)
    return M

def write_csv(rows, path):
    if not rows:
        return
    hdr = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def main():
    if len(sys.argv) < 3:
        print("Usage:\n  python nav_rinex2_to_df.py <input.16n> <out_prefix> [--demo-time 'YYYY-MM-DD HH:MM:SS' --rx 'X,Y,Z']")
        sys.exit(2)
    in_nav = sys.argv[1]
    out_prefix = sys.argv[2]

    ephs = parse_rinex2_nav(in_nav)
    if not ephs:
        print("No ephemerides parsed. Check RINEX version (2.xx) and file.")
        sys.exit(1)

    out_csv = f"{out_prefix}.gps_nav.csv"
    rows = []
    for e in ephs:
        d = dict(e)
        d["toc"] = e["toc"].isoformat()
        rows.append(d)
    write_csv(rows, out_csv)
    print(f"OK -> wrote {out_csv}  (ephemerides: {len(rows)})")

    if "--demo-time" in sys.argv and "--rx" in sys.argv:
        tstr = sys.argv[sys.argv.index("--demo-time")+1]
        rx_str = sys.argv[sys.argv.index("--rx")+1]
        t_demo = datetime.fromisoformat(tstr).replace(tzinfo=timezone.utc)
        rx_xyz = tuple(float(s) for s in rx_str.split(","))

        demo_rows = []
        week_t, tow_t = gpsweek(t_demo)
        prns = sorted({e["prn"] for e in ephs})
        for prn in prns:
            recs = [e for e in ephs if e["prn"]==prn]
            best = None; best_abs = 1e18
            for e in recs:
                wk = int(e["week"]) if not math.isnan(e["week"]) else week_t
                tk = (week_t - wk)*604800.0 + (tow_t - e["Toe"])
                tk = _wrap_tk(tk)
                if abs(tk) < best_abs:
                    best_abs = abs(tk); best = e
            if best is None:
                continue
            sat_xyz = sat_ecef_from_ephem(best, t_demo)
            az, el = az_el_from_ecef(rx_xyz, sat_xyz)
            M = vtec_mapping(el)
            demo_rows.append({
                "time": t_demo.isoformat(), "prn": prn,
                "az_deg": math.degrees(az), "el_deg": math.degrees(el),
                "mapping_M": M
            })
        demo_csv = f"{out_prefix}.demo_positions.csv"
        write_csv(demo_rows, demo_csv)
        print(f"OK -> wrote {demo_csv}  (records: {len(demo_rows)})")

if __name__ == "__main__":
    main()
