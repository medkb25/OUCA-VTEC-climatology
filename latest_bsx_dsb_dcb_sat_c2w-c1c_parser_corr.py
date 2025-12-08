#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re, math, gzip, os, argparse, csv

C = 299_792_458.0  # m/s
_OBS_ALIAS = {"C1P":"C1W","C2P":"C2W"}  # GPS: P-code→W

# ---------- utils ----------
def _unit_to_m(v: float, unit: str) -> float:
    u = (unit or "ns").strip().lower()
    if u in ("m","meter","metre"): return v
    return v * 1e-9 * C  # défaut ns

def _norm_obs(code: str) -> str:
    return _OBS_ALIAS.get(code.strip().upper(), code.strip().upper())

def _is_prn(tok: str) -> bool: return bool(re.match(r"^G\d{2}$", tok))
def _is_obs(tok: str) -> bool: return bool(re.match(r"^[A-Z0-9]{3}$", tok))
def _open_text(path: str):
    return gzip.open(path, "rt", errors="ignore") if path.lower().endswith(".gz") else open(path, "r", errors="ignore")

# ---------- SINEX bias ----------
def _extract_line_fields(tokens, default_unit: str):
    if len(tokens) < 8: return None
    tp = tokens[0]
    if tp not in ("DCB","DSB"): return None
    prn_idx = next((i for i,t in enumerate(tokens[:8]) if _is_prn(t)), None)
    if prn_idx is None: return None
    obs = []
    for t in tokens[prn_idx+1:prn_idx+7]:
        if _is_obs(t):
            obs.append(_norm_obs(t))
            if len(obs)==2: break
    if len(obs)!=2: return None
    unit = tokens[-3] if re.match(r"^[A-Za-z]+$", tokens[-3]) else (default_unit or "ns")
    try:
        val = float(tokens[-2])
    except Exception:
        return None
    return tp, tokens[prn_idx], obs[0], obs[1], unit, _unit_to_m(val, unit)

def parse_sinex_bias(path: str, debug=False):
    rows, inside, default_unit = [], False, "ns"
    with _open_text(path) as f:
        for raw in f:
            s = raw.strip()
            if not s or s.startswith(("*","%")): continue
            if "BIAS_UNIT" in s and "+BIAS/SOLUTION" not in s and "+FILE/REFERENCE" not in s:
                parts = s.split()
                for i,t in enumerate(parts[:-1]):
                    if t.upper()=="BIAS_UNIT": default_unit = parts[i+1].lower(); break
            if "+BIAS/SOLUTION" in s: inside=True; continue
            if "-BIAS/SOLUTION" in s: inside=False; continue
            if not inside: continue
            toks = s.split()
            got = _extract_line_fields(toks, default_unit)
            if not got: continue
            tp, prn, o1, o2, unit, val = got
            if prn[0]!="G": continue
            rows.append((prn, f"{o1}-{o2}", val, tp))
    if debug: print(f"[SINEX] {os.path.basename(path)}: {len(rows)} lignes")
    return rows

# ---------- CODE mensuels ----------
def parse_code_monthly_table(path: str, debug=False):
    """
    Lit un fichier texte CODE mensuel décompressé (ex: P1C12201.DCB ou P1P22201.DCB).
    Retour: dict { 'Gxx': value_m } ; la nature (P1-C1 ou P1-P2) dépend du fichier.
    Hypothèse format 'G01  <val_ns>  <rms_ns>'.
    """
    out = {}
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if not line.startswith("G"): continue
            parts = line.split()
            if len(parts) < 2: continue
            prn = parts[0]
            try:
                val_ns = float(parts[1])
            except Exception:
                continue
            out[prn] = _unit_to_m(val_ns, "ns")
    if debug: print(f"[CODE] {os.path.basename(path)}: {len(out)} PRN")
    return out

# ---------- fusion logique ----------
def build_table(rows, prefer=("DSB","DCB"), monthly_p1c1=None, monthly_p1p2=None):
    """
    rows: [(prn, pair, value_m, bias_type)] depuis SINEX journalier
    monthly_p1c1: dict PRN -> (P1-C1) = (C1W-C1C) [m]
    monthly_p1p2: dict PRN -> (P1-P2) = (C1W-C2W) [m]
    Retour: [(prn, C2W-C1C_m, source)]
    """
    data = {}
    for prn, pair, val, tp in rows:
        data.setdefault(prn, {}).setdefault(tp, {})[pair] = val

    out = []
    for prn, by_type in sorted(data.items()):
        p2c1, src = None, "missing"

        # 1) direct / inversé
        for tp in prefer:
            d = by_type.get(tp, {})
            if "C2W-C1C" in d:      p2c1, src = d["C2W-C1C"], f"direct:{tp}"; break
            if "C1C-C2W" in d:      p2c1, src = -d["C1C-C2W"], f"direct:inverted:{tp}"; break

        # 2) dérivé même type:  -( (C1C-C1W) + (C1W-C2W) )
        if p2c1 is None:
            for tp in prefer:
                d = by_type.get(tp, {})
                a = d.get("C1C-C1W"); b = d.get("C1W-C2W")
                if a is not None and b is not None:
                    p2c1, src = -(a+b), f"derived:{tp}"; break

        # 3) dérivé mix SINEX
        if p2c1 is None:
            d_dsb, d_dcb = by_type.get("DSB", {}), by_type.get("DCB", {})
            a = d_dsb.get("C1C-C1W", d_dcb.get("C1C-C1W"))
            b = d_dsb.get("C1W-C2W", d_dcb.get("C1W-C2W"))
            if a is not None and b is not None:
                p2c1 = -(a+b)
                tag_a = "DSB" if "C1C-C1W" in d_dsb else ("DCB" if "C1C-C1W" in d_dcb else "NA")
                tag_b = "DSB" if "C1W-C2W" in d_dsb else ("DCB" if "C1W-C2W" in d_dcb else "NA")
                src = "derived:mixed" if tag_a!=tag_b else f"derived:{tag_a}"

        # 4) SINEX(C1W-C2W) + mensuel (P1-C1)
        if p2c1 is None and monthly_p1c1:
            # prendre C1W-C2W ou son inverse si dispo
            c1w_c2w = None
            for d in by_type.values():
                if "C1W-C2W" in d: c1w_c2w = d["C1W-C2W"]; break
                if "C2W-C1W" in d: c1w_c2w = -d["C2W-C1W"]; break
            c1w_c1c = monthly_p1c1.get(prn)
            if c1w_c2w is not None and c1w_c1c is not None:
                # C2W-C1C = (C1W-C1C) - (C1W-C2W)
                p2c1 = c1w_c1c - c1w_c2w
                src = "derived:sinex(C1W-C2W)+monthly(P1-C1)"

        # 5) Mensuels seuls: (P1-C1) - (P1-P2)
        if p2c1 is None and monthly_p1c1 and monthly_p1p2:
            a = monthly_p1c1.get(prn)   # C1W-C1C
            b = monthly_p1p2.get(prn)   # C1W-C2W
            if a is not None and b is not None:
                p2c1 = a - b            # = C2W-C1C
                src  = "derived:monthly(P1-C1)-(P1-P2)"

        out.append((prn, float('nan') if p2c1 is None else float(p2c1), src))
    return out

def save_table_csv(rows, out_csv: str):
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["prn","C2W-C1C_sat_m","source"]); w.writerows(rows)

# ---------- CLI ----------
def main():
    p = argparse.ArgumentParser(description="C2W−C1C from SINEX-Bias with CODE monthly fallbacks")
    p.add_argument("in_bsx", help="Input *.BIA|*.BSX[.gz]")
    p.add_argument("out_csv", help="Output CSV")
    p.add_argument("--prefer", default="DSB,DCB", help="Order e.g. 'DSB,DCB' or 'DCB,DSB'")
    p.add_argument("--monthly_p1c1", help="CODE monthly P1-C1 (decompressed .DCB)")
    p.add_argument("--monthly_p1p2", help="CODE monthly P1-P2 (decompressed .DCB)")
    p.add_argument("--debug", action="store_true")
    args = p.parse_args()

    prefer = tuple(s.strip().upper() for s in args.prefer.split(",") if s.strip())
    rows   = parse_sinex_bias(args.in_bsx, debug=args.debug)

    m_p1c1 = parse_code_monthly_table(args.monthly_p1c1, debug=args.debug) if args.monthly_p1c1 else None
    m_p1p2 = parse_code_monthly_table(args.monthly_p1p2, debug=args.debug) if args.monthly_p1p2 else None

    tbl = build_table(rows, prefer=prefer, monthly_p1c1=m_p1c1, monthly_p1p2=m_p1p2)
    save_table_csv(tbl, args.out_csv)
    if args.debug:
        n_good = sum(1 for _,v,_ in tbl if not math.isnan(v))
        print(f"OK PRN: {n_good}/{len(tbl)}")

if __name__ == "__main__":
    main()
