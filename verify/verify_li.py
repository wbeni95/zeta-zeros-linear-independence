"""verify_li.py: rigorous verification of the m-independence certificates for the first n zeros of the Riemann zeta function
(THEOREM.md).  Uses FLINT/Arb through python-flint.

Claim (C_M) for a basis of size n: there is no (m_0, m_1, ..., m_N) in Z^n, N = n - 1, with (m_1..m_N) != 0, max_{j>=1} |m_j| <= M and
    m_0 gamma_1 + sum_{j=1}^{N} m_j gamma_{j+1} = 0.
Certificate: rows b_i in L (exact check), full rank, UB^2 >= max_i ||b*_i||^2 and E1 >= sum_j |Y_j - 2^nu x_j| give (C_M) whenever
    M^2 (N + E1^2) UB^2 < 2^(2 nu),      L = { (v, z) in Z^(N+1) : v_j = z Y_j mod 2^nu },  x_j = gamma_{j+1}/gamma_1,  nu = 32000.
Arithmetic:
  * D_i = ||b*_i||^2 enclosed by Arb ball arithmetic: every row entry x -> ball [x>>sh, (x>>sh)+1] * 2^sh, ball Gram product,
    right-looking (Schur complement) ball LDL; precision doubled until every relative radius < 2^-40.  Full rank over R <=> all
    D_i certainly > 0; second rank witness: rank mod 2^59-55.
  * zeros: decimal midpoints from data/zeros_1_701_midpoints.txt as exact rationals +- DELTA, with DELTA proved by Hardy-Z sign
    changes and Turing-indexed neighbour balls (mode 'zeros'); E1 = exact rational sum of dyadic upper endpoints.
  * M: exact integer / Fraction arithmetic.
Modes:
  python verify_li.py zeros [part nparts]        prove |gamma_i - mid_i| < DELTA for i = 1..700 (writes output/zero_proofs_*.jsonl)
  python verify_li.py cert <n> <basis.txt.gz> [claimed_M]
  python verify_li.py table                      all bases listed in ../results/LI_TABLE.json
By default 'cert' and 'table' use the zero proofs in output/ if present (all 700), otherwise the shipped ../results/zero_proofs_shipped.jsonl.
Requirements: Python 3.11+, python-flint >= 0.6."""
import glob
import gzip
import hashlib
import json
import os
import sys
import time
from fractions import Fraction
from math import isqrt, log2

from flint import acb, arb, arb_mat, ctx, nmod_mat

sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
OUT = os.path.join(HERE, "output")
ZEROS = os.path.join(PKG, "data", "zeros_1_701_midpoints.txt")
TABLE = os.path.join(PKG, "results", "LI_TABLE.json")
SHIPPED_PROOFS = os.path.join(PKG, "results", "zero_proofs_shipped.jsonl")
NU = 32000
NMAX = 700
DELTA = "2e-9665"           # proved radius around the midpoints
DELTA_EVAL = "1.999e-9665"  # Hardy-Z evaluation offset
ZPREC = 32600
P_RANK = 2 ** 59 - 55


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def dyadic(x):
    """exact Fraction value of an exact arb (zero radius), e.g. x.upper()"""
    m, e = x.man_exp()
    m, e = int(m), int(e)
    return Fraction(m * 2 ** e) if e >= 0 else Fraction(m, 2 ** (-e))


def up(x):
    return dyadic(x.upper())


def lo(x):
    return dyadic(x.lower())


# ------------------------------------------------------------------------------------------------ zero data
def read_mids():
    return [l.strip() for l in open(ZEROS, encoding="ascii") if l.strip() and not l.startswith("#")]


def dec_ball(s):
    """rigorous ball around the exact decimal rational s (current ctx.prec)"""
    ip, _, fp = s.partition(".")
    return arb(int(ip + fp)) / arb(10) ** len(fp)


def sha1s(s):
    return hashlib.sha1(s.encode()).hexdigest()


def hardy_Z(t):
    """Z(t) = exp(i theta(t)) zeta(1/2 + i t), theta(t) = Im lgamma(1/4 + i t/2) - (t/2) log pi; returns an acb ball"""
    th = acb(arb(1) / 4, t / 2).lgamma().imag - t / 2 * arb.pi().log()
    return acb(0, th).exp() * acb(arb(1) / 2, t).zeta()


def prove_zero(i, mid, nb_im):
    """sign change of Hardy Z at mid -+ DELTA_EVAL, and the window [mid -+ DELTA] strictly between the Turing-indexed isolating
    balls of zeros i-1 and i+1 (and meeting the ball of zero i)  =>  the unique zero gamma_i satisfies |gamma_i - mid| < DELTA"""
    old = ctx.prec
    try:
        ctx.prec = ZPREC
        m = dec_ball(mid)
        de = arb(DELTA_EVAL)
        tm, tp = m - de, m + de
        assert tm.rad() < arb("1e-9668") and tp.rad() < arb("1e-9668")
        Zm, Zp = hardy_Z(tm), hardy_Z(tp)
        real_ok = Zm.imag.contains(0) and Zp.imag.contains(0)
        sm = 1 if Zm.real > 0 else (-1 if Zm.real < 0 else 0)
        sp = 1 if Zp.real > 0 else (-1 if Zp.real < 0 else 0)
        sign_change = sm * sp == -1
        ctx.prec = 256
        dd = arb(DELTA)
        wl, wr = dec_ball(mid) - dd, dec_ball(mid) + dd
        left_ok = True if i == 1 else bool(nb_im[i - 2] < wl)
        right_ok = bool(wr < nb_im[i])
        own_overlap = nb_im[i - 1].overlaps(arb(dec_ball(mid).mid(), dd))
        ok = bool(real_ok and sign_change and left_ok and right_ok and own_overlap)
        return {"i": i, "ok": ok, "delta": DELTA, "sign_minus": sm, "sign_plus": sp, "Z_real": bool(real_ok), "left_ok": left_ok,
                "right_ok": right_ok, "own_ball_overlap": bool(own_overlap), "mid_sha1": sha1s(mid), "prec": ZPREC}
    finally:
        ctx.prec = old


def mode_zeros(part=0, nparts=1):
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, f"zero_proofs_part{part}of{nparts}.jsonl")
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)["i"] for l in open(out) if l.strip()}
    mids = read_mids()
    ctx.prec = 128
    nb_im = [b.imag for b in acb.zeta_zeros(1, NMAX + 1)]        # Turing-indexed isolating balls of zeros 1..701
    todo = [i for i in range(1, NMAX + 1) if (i - 1) % nparts == part and i not in done]
    t0 = time.time()
    for k, i in enumerate(todo, 1):
        rec = prove_zero(i, mids[i - 1], nb_im)
        with open(out, "a") as f:
            f.write(json.dumps(rec) + "\n")
        if not rec["ok"]:
            log("ZERO PROOF FAILED", rec)
        if k % 25 == 0:
            log(f"zeros: {k}/{len(todo)} ({time.time() - t0:.0f} s)")
    log(f"zeros part {part}/{nparts}: {len(todo)} new proofs in {time.time() - t0:.0f} s")


def load_proofs():
    own = {}
    for f in glob.glob(os.path.join(OUT, "zero_proofs_part*.jsonl")):
        for l in open(f):
            if l.strip():
                r = json.loads(l); own[r["i"]] = r
    if all(i in own for i in range(1, NMAX + 1)):
        return own, "own (output/)"
    shipped = {}
    for l in open(SHIPPED_PROOFS):
        if l.strip():
            r = json.loads(l); shipped[r["i"]] = r
    return shipped, "shipped (results/zero_proofs_shipped.jsonl)"


def proofs_status(mids, n, recs):
    bad = [i for i in range(1, n + 1) if not (i in recs and recs[i]["ok"] and recs[i]["delta"] == DELTA
                                                and recs[i]["mid_sha1"] == sha1s(mids[i - 1]))]
    return len(bad) == 0, bad[:10], len(bad)


def y_e1(n, nu, mids, delta=DELTA):
    """Y_j (j = 1..N) as Python ints, E1 = exact Fraction upper bound of sum_j |Y_j - 2^nu x_j|, E (max), ambiguity count;
    gamma_i in [mid_i - delta, mid_i + delta]"""
    old = ctx.prec
    ctx.prec = nu + 800
    try:
        dl = arb("[+/- " + delta + "]")
        g = [dec_ball(mids[i]) + dl for i in range(n)]
        two = arb(2) ** nu
        T = [(g[j] * two) / g[0] for j in range(1, n)]
        Y, E1, E, amb = [], Fraction(0), Fraction(0), 0
        for t in T:
            c = dyadic(t.mid())
            y = (2 * c.numerator + c.denominator) // (2 * c.denominator)     # round half up of the exact midpoint
            if (t - (arb(y) + arb(1) / 2)).contains(0) or (t - (arb(y) - arb(1) / 2)).contains(0):
                amb += 1
            e = up(abs(arb(y) - t))
            E1 += e
            E = max(E, e)
            Y.append(y)
        return Y, E1, E, amb
    finally:
        ctx.prec = old


# ------------------------------------------------------------------------------------------------ lattice checks
def read_basis(path):
    rows = []
    with gzip.open(path, "rt", encoding="ascii") as f:
        for l in f:
            if l.startswith("#") or not l.strip():
                continue
            rows.append([int(t, 16) for t in l.split()])
    return rows


def membership(rows, Y, nu):
    """every row (v_1..v_N, z) satisfies v_j = z Y_j mod 2^nu (exact)"""
    N = len(Y)
    mask = (1 << nu) - 1          # x mod 2^nu == x & mask for every Python int (two's complement semantics)
    Ym = [y & mask for y in Y]
    bad = 0
    for r in rows:
        if len(r) != N + 1:
            return False, -1
        z = r[N] & mask
        for j in range(N):
            if (r[j] - z * Ym[j]) & mask:
                bad += 1
                break
    return bad == 0, bad


def rank_mod_p(rows, p=P_RANK):
    return nmod_mat([[x % p for x in r] for r in rows], p).rank()


def ball_gram(rows, prec):
    """ball enclosure of the Gram matrix of the rows, scaled by 2^(-2 sh)"""
    mb = max(abs(x).bit_length() for r in rows for x in r)
    sh = max(0, mb - (prec - 16))
    ctx.prec = prec
    half = arb("[0.5 +/- 0.5]")                         # [0, 1]
    if sh:
        A = arb_mat([[arb(x >> sh) + half for x in r] for r in rows])   # x / 2^sh in [x>>sh, (x>>sh) + 1]
    else:
        A = arb_mat([[arb(x) for x in r] for r in rows])
    G = A * A.transpose()
    d = len(rows)
    S = [[G[i, j] for j in range(i + 1)] for i in range(d)]
    del A, G
    return S, sh


def ball_ldl(S):
    """right-looking LDL^T of the symmetric ball matrix (lower triangle lists, modified in place); returns (D, fail_k)"""
    d = len(S)
    D = []
    for k in range(d):
        Dk = S[k][k]
        if not (Dk > 0):
            return D, k
        D.append(Dk)
        col = [S[i][k] for i in range(k + 1, d)]
        for a in range(len(col)):
            i = k + 1 + a
            la = col[a] / Dk
            Si = S[i]
            Si[k + 1:i + 1] = [s - la * c for s, c in zip(Si[k + 1:i + 1], col[:a + 1])]
            Si[k] = None
        S[k] = None
    return D, None


def gs_enclosure(rows, prec0=1024, maxprec=32768):
    """D_i balls (scaled back exactly as Fractions), precision used, max relative radius"""
    prec = prec0
    while True:
        S, sh = ball_gram(rows, prec)
        D, fail = ball_ldl(S)
        rel = max((float(x.rad() / x.mid()) for x in D), default=1.0) if fail is None else None
        if fail is None and rel < 2.0 ** -40:
            sc = Fraction(2 ** (2 * sh))
            return {"D_lo": [lo(x) * sc for x in D], "D_up": [up(x) * sc for x in D], "prec": prec, "max_rel_rad": rel}
        log(f"  ball LDL at prec {prec}: fail_k={fail} max_rel={rel}; doubling")
        if prec >= maxprec:
            return {"error": f"no certified LDL up to prec {maxprec} (fail_k={fail}, rel={rel})", "prec": prec}
        prec *= 2


def M_max(N, E1, Dup, nu):
    """largest M >= 0 with M^2 (N + E1^2) Dup < 2^(2 nu), exact rationals"""
    A = (N + E1 * E1) * Dup
    R = Fraction(2 ** (2 * nu)) / A
    M = isqrt(R.numerator // R.denominator)
    while M > 0 and M * M >= R:
        M -= 1
    while (M + 1) * (M + 1) < R:
        M += 1
    return M


def flog2(fr):
    return log2(fr.numerator) - log2(fr.denominator)


def certify(path, n, claimed=None):
    t0 = time.time()
    mids = read_mids()
    recs, src = load_proofs()
    rows = read_basis(path)
    N = n - 1
    Y, E1, E, amb = y_e1(n, NU, mids)
    memb, nbad = membership(rows, Y, NU)
    rk = rank_mod_p(rows)
    gs = gs_enclosure(rows)
    pr_ok, pr_bad, pr_nbad = proofs_status(mids, n, recs)
    res = {"basis": os.path.basename(path), "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(), "n": n, "nu": NU,
           "rows": len(rows), "membership": memb, "bad_rows": nbad, "rank_mod_p": rk, "p": P_RANK, "Y_ambiguous": amb,
           "E1": float(E1), "E": float(E), "zero_proofs": src, "zero_proofs_ok": pr_ok, "zero_proofs_missing": pr_nbad}
    if "error" in gs:
        res.update({"gs_error": gs["error"], "full_rank": False, "M": "0", "claim_verified": False})
        return res
    Dup = max(gs["D_up"])
    posdef = all(x > 0 for x in gs["D_lo"])
    full = posdef and rk == len(rows) == N + 1
    ok = memb and full and pr_ok
    M = M_max(N, E1, Dup, NU) if ok else 0
    res.update({"full_rank": full, "ldl_prec": gs["prec"], "ldl_max_rel_rad": gs["max_rel_rad"], "log2_maxGS_upper": 0.5 * flog2(Dup),
                "M": str(M), "log10_M": (log2(M) / log2(10)) if M > 0 else None,
                "claimed_M": str(claimed) if claimed is not None else None,
                "claim_verified": bool(ok and claimed is not None and M >= claimed), "s": round(time.time() - t0, 1)})
    return res


def save(res):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "certificates.jsonl"), "a") as f:
        f.write(json.dumps(res) + "\n")


def mode_table():
    tab = json.load(open(TABLE))
    allok = True
    for e in tab["rows"]:
        path = os.path.join(PKG, "bases", e["basis"])
        if not os.path.exists(path):
            log(f"n={e['n']}: {e['basis']} not found in bases/ (download it from the Zenodo record), skipped")
            allok = False
            continue
        r = certify(path, e["n"], int(e["M"]))
        save(r)
        log(f"n={e['n']}: M={r.get('M')} claimed={e['M']} verified={r['claim_verified']} memb={r['membership']} full_rank={r['full_rank']} "
            f"zero_proofs_ok={r['zero_proofs_ok']} ({r.get('s')} s)")
        allok &= r["claim_verified"]
    print("all claims verified" if allok else "NOT ALL CLAIMS VERIFIED")
    return allok


if __name__ == "__main__":
    md = sys.argv[1] if len(sys.argv) > 1 else ""
    if md == "zeros":
        mode_zeros(int(sys.argv[2]) if len(sys.argv) > 2 else 0, int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    elif md == "cert":
        r = certify(sys.argv[3], int(sys.argv[2]), int(sys.argv[4]) if len(sys.argv) > 4 else None)
        save(r)
        print(json.dumps(r, indent=1))
        sys.exit(0 if r["claim_verified"] or len(sys.argv) <= 4 else 1)
    elif md == "table":
        sys.exit(0 if mode_table() else 1)
    else:
        print(__doc__)
