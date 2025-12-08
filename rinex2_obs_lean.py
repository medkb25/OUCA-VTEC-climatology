#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rinex2_obs_lean.py — Minimal RINEX 2 observation parser (GPS) for VTEC:
    returns a DataFrame with columns:
        - time (UTC)
        - prn
        - P2_minus_C1_m  (P2 − C1 in meters)

Assumptions:
- RINEX 2.x (2.10 / 2.11) observation files
- Uses C1 (L1 code pseudorange) and P2 (L2 code pseudorange), in meters
- GPS only (PRN starting with "G"); other GNSS constellations are ignored

Usage:
    from rinex2_obs_lean import parse_rinex2_obs_lean
    df = parse_rinex2_obs_lean("RS45_20230517.23o")
"""

from pathlib import Path
import math, re
import pandas as pd

_OBS_LINE_WIDTH = 16   # largeur par observable dans RINEX 2.x
_PER_LINE = 5          # 5 observables par ligne
_epoch_re = re.compile(
    r"^\s*(?P<yy>\d{2})\s+(?P<mo>\d{1,2})\s+(?P<dy>\d{1,2})\s+(?P<hh>\d{1,2})\s+(?P<mm>\d{1,2})\s+(?P<ss>[0-9.]+)\s+(?P<flag>\d)\s+(?P<nsat>\d+)"
)

def _yy_to_year(two):
    return 1900 + two if two >= 80 else 2000 + two

def _parse_header(fh):
    """Lit l'en-tête et retourne (obs_types:list[str], meta:dict)."""
    obs_types, meta = [], {}
    ntypes_expected = None
    pushback = None

    while True:
        line = pushback if pushback is not None else fh.readline()
        pushback = None
        if not line:
            break

        if "APPROX POSITION XYZ" in line:
            try:
                meta["rx_X"] = float(line[0:14])
                meta["rx_Y"] = float(line[14:28])
                meta["rx_Z"] = float(line[28:42])
            except Exception:
                pass

        if "# / TYPES OF OBSERV" in line:
            # nombre déclaré
            try:
                ntypes_expected = int(line[0:6])
            except Exception:
                ntypes_expected = int(line.split()[0])
            # types présents sur cette ligne (dans les 60 premières colonnes)
            part = line[:60]
            toks = [t for t in part.split() if re.fullmatch(r"[A-Z0-9]{2,3}", t)]
            toks = [t for t in toks if not t.isdigit()]
            obs_types.extend(toks)

            # lignes de continuation (# / TYPES OF OBSERV)
            while ntypes_expected is not None and len(obs_types) < ntypes_expected:
                cont = fh.readline()
                if not cont:
                    break
                if "# / TYPES OF OBSERV" not in cont:
                    pushback = cont
                    break
                part = cont[:60]
                toks = [t for t in part.split() if re.fullmatch(r"[A-Z0-9]{2,3}", t)]
                toks = [t for t in toks if not t.isdigit()]
                obs_types.extend(toks)

        if "END OF HEADER" in line:
            break

    return obs_types, meta

def _parse_epoch_header(line):
    m = _epoch_re.match(line)
    if not m:
        return None
    yy = int(m.group("yy")); mo = int(m.group("mo")); dy = int(m.group("dy"))
    hh = int(m.group("hh")); mm = int(m.group("mm")); ss = float(m.group("ss"))
    nsat = int(m.group("nsat"))
    year = _yy_to_year(yy)
    t = pd.Timestamp(year=year, month=mo, day=dy, hour=hh, minute=mm, second=0, tz="UTC") + pd.to_timedelta(ss, unit="s")
    remainder = line[m.end():].strip()
    return t, nsat, remainder

def _collect_epoch_prns(first_remainder, nsat, fh):
    # Extrait PRN concaténés type "G05G07..."; peut déborder sur lignes suivantes
    prns = re.findall(r"[A-Z]\d{2}", first_remainder)
    while len(prns) < nsat:
        cont = fh.readline()
        if not cont:
            break
        prns.extend(re.findall(r"[A-Z]\d{2}", cont.strip()))
    return prns[:nsat]

def _read_obs_values(fh, ntypes):
    """Lit ntypes observables (RINEX2: blocs de 16 colonnes; 5 par ligne)."""
    vals = []
    remaining = ntypes
    nlines = math.ceil(ntypes / _PER_LINE)
    for _ in range(nlines):
        raw = fh.readline()
        if not raw:
            break
        line = raw.rstrip("\n")
        for i in range(min(_PER_LINE, remaining)):
            chunk = line[i*_OBS_LINE_WIDTH:(i+1)*_OBS_LINE_WIDTH]
            # Valeur sur colonnes 1..14 (LLI/SNR en 15..16)
            val_str = chunk[:14].strip().replace("D", "E")
            if val_str:
                try:
                    vals.append(float(val_str))
                except Exception:
                    vals.append(None)
            else:
                vals.append(None)
        remaining -= _PER_LINE
        if remaining <= 0:
            break
    while len(vals) < ntypes:
        vals.append(None)
    return vals

def parse_rinex2_obs_lean(path):
    """
    Retourne DataFrame: time (UTC), prn, P2_minus_C1_m
    - Ignore les systèmes ≠ GPS (PRN non 'Gxx')
    - Nécessite que C1 et P2 soient déclarés dans l'en-tête et présents en data
    """
    rows = []
    with open(path, "r", errors="ignore") as fh:
        obs_types, _meta = _parse_header(fh)
        if not obs_types:
            return pd.DataFrame(columns=["time","prn","P2_minus_C1_m"])

        # indices C1 / P2 (ordre exact déclaré dans l'en-tête)
        try:
            idx_C1 = obs_types.index("C1")
            idx_P2 = obs_types.index("P2")
        except ValueError:
            # pas exploitable pour P2-C1
            return pd.DataFrame(columns=["time","prn","P2_minus_C1_m"])

        ntypes = len(obs_types)

        # lecture du corps
        while True:
            raw = fh.readline()
            if not raw:
                break
            ep = _parse_epoch_header(raw.rstrip("\n"))
            if not ep:
                continue
            t_utc, nsat, remainder = ep
            prns = _collect_epoch_prns(remainder, nsat, fh)

            for prn in prns:
                # ne garder que GPS
                if not prn.startswith("G"):
                    # consommer quand même les lignes d'observables
                    _ = _read_obs_values(fh, ntypes)
                    continue

                vals = _read_obs_values(fh, ntypes)
                v_c1 = vals[idx_C1] if idx_C1 is not None else None
                v_p2 = vals[idx_P2] if idx_P2 is not None else None
                if v_c1 is None or v_p2 is None:
                    continue
                rows.append((t_utc, prn, float(v_p2 - v_c1)))

    df = pd.DataFrame(rows, columns=["time","prn","P2_minus_C1_m"])
    if not df.empty:
        df["time"] = pd.to_datetime(df["time"], utc=True)
        df["prn"] = df["prn"].astype("category")
        df["P2_minus_C1_m"] = df["P2_minus_C1_m"].astype("float32")
    return df

# Petit test manuel:
if __name__ == "__main__":
    import sys
    if len(sys.argv) >= 2:
        df = parse_rinex2_obs_lean(sys.argv[1])
        print(df.head())
        print(df.describe(include="all"))
