"""Cross-checks for the two refinements proved in PAPER.en.md.

Theorem G. On Gamma_1(12), the coordinates (A, B, C) of (chi_3, chi_4, 4)
satisfy A ≡ 5 p^2, B ≡ p d, C ≡ d - 1 (mod 9), where p = c/3. Hence
U = (3/2) S-tilde is a square modulo 9. The written argument is the proof.
This file checks the nine generator matrices, their inverses, the four
explicit residues of U, and, as a cross-check, closure of the law on a
complete set of residues modulo 36.

Theorem H. On Gamma_1(15), writing b_n = 3^{n+1} a_n / 4 and r = -c/3,
one has b_0 ≡ 2 r^2, b_1 ≡ 2 r, and b_n ≡ 0 for n >= 2 (mod 3). Hence
(3/4) S-tilde is 0 or 1 modulo 3, never 2. This file checks that law on
every saved Schreier polynomial and on each of those polynomials slashed
by (1, 0, -15, 1).

Remark after Theorem C. Modulo 16 every even residue of S-tilde occurs, and
modulo 36 every residue compatible with S-tilde ≡ 0 or 4 (mod 6) occurs.
The quotients below are the proof of that sharpness statement. A path of
exact slash operations is compared with the reduced slash.

Standard library only. From the repository root:

    python src/verify_refinements.py
"""

from collections import deque
from fractions import Fraction
from math import comb
from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_k7 import stilde

ROOT = Path(__file__).resolve().parents[1]
LOG = []


def say(msg):
    print(msg, flush=True)
    LOG.append(msg)


def matmul(m, g):
    a, b, c, d = m
    e, f, gg, h = g
    return (a * e + b * gg, a * f + b * h, c * e + d * gg, c * f + d * h)


def matinv(m):
    a, b, c, d = m
    return (d, -b, -c, a)


def parse_alphabet(name):
    gens = []
    for line in (ROOT / "certificates" / name).read_text(encoding="utf-8").splitlines():
        if line.startswith("A "):
            parts = [int(x) for x in line.split()[1:]]
            gens.append(tuple(parts[1:]))
    return gens


# ---------- Theorem G ----------

LINES = {
    (1, 1, 0, 1): (0, 0, 0),
    (-11, 1, -12, 1): (8, -4, 0),
    (49, -9, 60, -11): (560, -616, 168),
    (13, -3, 48, -11): (-1087, 1426, -471),
    (13, -4, 36, -11): (-72, 120, -48),
    (61, -25, 144, -59): (-1800, 4422, -2715),
    (13, -6, 24, -11): (-40, 110, -75),
    (25, -16, 36, -23): (-216, 840, -816),
    (37, -27, 48, -35): (1127, -4862, 5247),
}


def slash12(tri, g):
    A, B, C = tri
    a, b, c, d = g
    p = c // 3
    return (
        A * a * a + B * a * p + C * p * p,
        6 * A * a * b + B * (a * d + b * c) + 2 * C * p * d,
        9 * A * b * b + 3 * B * b * d + C * d * d,
    )


def U_of(tri, mat):
    A, B, C = tri
    return A * mat[3] * mat[3] - B * mat[3] * (mat[2] // 3) + C * (mat[2] // 3) ** 2


def law12(tri, mat):
    A, B, C = tri
    p = mat[2] // 3
    d = mat[3]
    return (A - 5 * p * p) % 9 == 0 and (B - p * d) % 9 == 0 and (C - (d - 1)) % 9 == 0


def check_theorem_g():
    polys = dict(LINES)
    polys[(1, 0, 0, 1)] = (0, 0, 0)
    for mat, tri in list(polys.items()):
        inv = matinv(mat)
        hinv = tuple(-x for x in slash12(tri, inv))
        if inv in polys and polys[inv] != hinv:
            raise SystemExit(f"inverse clash {inv}")
        polys[inv] = hinv
        if not law12(hinv, inv):
            raise SystemExit(f"inverse leaves the mod-9 law {inv}")
    for mat, tri in LINES.items():
        if not law12(tri, mat):
            raise SystemExit(f"generator leaves the mod-9 law {mat}")
    witnesses = {
        0: ((1, 0, 0, 1), (0, 0, 0), 0),
        1: ((-11, 1, -12, 1), LINES[(-11, 1, -12, 1)], -8),
        4: ((13, -6, 24, -11), LINES[(13, -6, 24, -11)], 40),
        7: ((49, -9, 60, -11), LINES[(49, -9, 60, -11)], -560),
    }
    for residue, (mat, tri, expected) in witnesses.items():
        got = U_of(tri, mat)
        if got != expected or got % 9 != residue:
            raise SystemExit((residue, got))
    say("theorem G: nine generators and their inverses obey the mod-9 law")
    say("theorem G: U takes 0, -8, 40, -560, hence all four squares modulo 9")

    fail = 0
    checked = 0
    for p in range(9):
        for dM in (1, 13, 25):
            A = 5 * p * p
            B = p * dM
            C = dM - 1
            for al in range(3):
                for be in range(3):
                    a = 1 + 12 * al
                    dg = 1 + 12 * be
                    for q in range(0, 36, 4):
                        c = 3 * q
                        target = (a * dg - 1) % 36
                        for b in range(36):
                            if (b * c - target) % 36:
                                continue
                            sA, sB, sC = slash12((A, B, C), (a, b, c, dg))
                            nA = sA + 5 * q * q
                            nB = sB + q * dg
                            nC = sC + (dg - 1)
                            np = p * a + dM * q
                            nd = 3 * p * b + dM * dg
                            checked += 1
                            if (nA - 5 * np * np) % 9 or (nB - np * nd) % 9 or (nC - (nd - 1)) % 9:
                                fail += 1
    if fail or checked != 26244:
        raise SystemExit(f"closure cross-check checked {checked} fail {fail}")
    say(f"theorem G cross-check: {checked} residue pairs with det ≡ 1 (mod 36) stay in the law")


# ---------- Theorem H ----------

def slash_b(bs, g):
    coef = [Fraction(4 * b, 3 ** (n + 1)) for n, b in enumerate(bs)]
    new_a = [Fraction(0)] * 6
    aa, bb, cc, dd = g
    for n in range(6):
        for s in range(5 - n + 1):
            for t in range(n + 1):
                j = 5 - (s + t)
                new_a[j] += (
                    coef[n]
                    * comb(5 - n, s)
                    * aa ** s
                    * bb ** (5 - n - s)
                    * comb(n, t)
                    * cc ** t
                    * dd ** (n - t)
                )
    out = []
    for j in range(6):
        val = new_a[j] * (3 ** (j + 1)) / 4
        if val.denominator != 1:
            raise SystemExit(f"nonintegral slash {g} {val}")
        out.append(int(val))
    return tuple(out)


def law7(bs, mat):
    r = (-mat[2]) // 3
    if any(bs[n] % 3 for n in range(2, 6)):
        return False
    return bs[0] % 3 == (2 * r * r) % 3 and bs[1] % 3 == (2 * r) % 3


def check_theorem_h():
    text = (ROOT / "certificates" / "verify-k7-log.txt").read_text(encoding="utf-8")
    if not text.strip().endswith("K7 CERTIFICATE OK"):
        raise SystemExit("verify-k7-log.txt is not a finished certificate")
    polys = {(1, 0, 0, 1): (0, 0, 0, 0, 0, 0)}
    for line in text.splitlines():
        if not line.startswith("P "):
            continue
        left, right = line[2:].split("|")
        mat = tuple(int(x) for x in left.split())
        bs = tuple(int(x) for x in right.split())
        polys[mat] = bs
    if len(polys) != 47:
        raise SystemExit(f"expected 47 polynomials, got {len(polys)}")
    for mat, bs in polys.items():
        if mat[0] * mat[3] - mat[1] * mat[2] != 1:
            raise SystemExit(f"not in SL2 {mat}")
        if not law7(bs, mat):
            raise SystemExit(f"mod-3 law fails {mat}")
    anchor = (1, 0, -15, 1)
    if polys[anchor] != (4250, -3200, 900, -120, 0, 0):
        raise SystemExit("specimen polynomial changed")
    if stilde(1, 15) != Fraction(-17000, 3):
        raise SystemExit("finite sum moved")
    # (3/4) S-tilde(anchor inverse) = sum b_n d^{5-n} (-c/3)^n
    b = polys[anchor]
    total = 0
    rr = 1
    step = (-anchor[2]) // 3
    for n, bn in enumerate(b):
        total += bn * (anchor[3] ** (5 - n)) * rr
        rr *= step
    if total != -4250 or total % 3 != 1:
        raise SystemExit(total)
    bad = 0
    for mat, bs in polys.items():
        prod_b = tuple(x + y for x, y in zip(slash_b(bs, anchor), polys[anchor]))
        prod = matmul(mat, anchor)
        if not law7(prod_b, prod):
            bad += 1
    if bad:
        raise SystemExit(f"{bad} products leave the mod-3 law")
    say("theorem H: 46 Schreier polynomials and the identity obey b0 ≡ 2r^2, b1 ≡ 2r, b_n ≡ 0 (n>=2) mod 3")
    say("theorem H: each of those 47 polynomials, slashed by (1,0;-15,1) and added, stays in the law")
    say("theorem H: (3/4) S-tilde(1,0;15,1) = -4250 ≡ 1 (mod 3); the identity gives 0")


# ---------- sharpness modulo 16 and 36 ----------

RAW = {
    (1, 1, 0, 1): (0, 0, 0),
    (-8, 1, -9, 1): (2, -4, 0),
    (19, -4, 81, -17): (522, -1980, 208),
    (10, -3, 27, -8): (-174, 876, -124),
    (-26, 11, -45, 19): (-220, 1672, -352),
    (37, -16, 81, -35): (-612, 4752, -1024),
    (19, -12, 27, -17): (156, -1716, 526),
}


def slash9(tri, g):
    A, B1, C1 = tri
    a, b, c, d = g
    mu = c // 9
    return (
        A * a * a + B1 * a * mu + 9 * C1 * mu * mu,
        18 * A * a * b + B1 * (a * d + 9 * mu * b) + 18 * C1 * mu * d,
        9 * A * b * b + B1 * b * d + C1 * d * d,
    )


def V_of(tri, mat):
    A, B1, C1 = tri
    m = mat[2] // 9
    d = mat[3]
    return A * d * d - B1 * d * m + 9 * C1 * m * m


def pack(vals, mod):
    key = 0
    for v in vals:
        key = key * mod + (v % mod)
    return key


def check_level9_sharpness():
    alph = parse_alphabet("gamma1_9_words.txt")
    polys = dict(RAW)
    polys[(1, 0, 0, 1)] = (0, 0, 0)
    for mat, tri in list(polys.items()):
        inv = matinv(mat)
        polys[inv] = tuple(-x for x in slash9(tri, inv))
    for g in alph:
        if g not in polys:
            raise SystemExit(f"missing polynomial {g}")

    def quotient(mod):
        seen = {pack((0, 0, 0, 0, 1), mod)}
        q = deque([(0, 0, 0, 0, 1 % mod)])
        values = set()
        while q:
            A, B1, C1, m, d = q.popleft()
            values.add((A * d * d - B1 * d * m + 9 * C1 * m * m) % mod)
            for g in alph:
                tri = polys[g]
                s = slash9((A, B1, C1), g)
                a, b, c, dg = g
                nxt = (
                    (s[0] + tri[0]) % mod,
                    (s[1] + tri[1]) % mod,
                    (s[2] + tri[2]) % mod,
                    (m * a + d * (c // 9)) % mod,
                    (9 * m * b + d * dg) % mod,
                )
                key = pack(nxt, mod)
                if key not in seen:
                    seen.add(key)
                    q.append(nxt)
        return seen, values

    for mod, allowed in (
        (16, [n for n in range(16) if n % 2 == 0]),
        (36, [n for n in range(36) if n % 6 in (0, 4)]),
    ):
        seen, values = quotient(mod)
        missing = [n for n in allowed if n not in values]
        extra = [n for n in values if n not in allowed]
        if missing or extra:
            raise SystemExit(f"mod {mod} missing {missing} extra {extra}")
        say(f"level 9 mod {mod}: {len(seen)} states, every allowed residue occurs, no others")

    rng = random.Random(0)
    mismatches = 0
    for mod in (16, 36):
        for _ in range(20):
            mat = (1, 0, 0, 1)
            tri = (0, 0, 0)
            state = (0, 0, 0, 0, 1 % mod)
            for _step in range(6):
                g = alph[rng.randrange(len(alph))]
                tri = tuple(x + y for x, y in zip(slash9(tri, g), polys[g]))
                mat = matmul(mat, g)
                A, B1, C1, m, d = state
                s = slash9((A, B1, C1), g)
                a, b, c, dg = g
                state = (
                    (s[0] + polys[g][0]) % mod,
                    (s[1] + polys[g][1]) % mod,
                    (s[2] + polys[g][2]) % mod,
                    (m * a + d * (c // 9)) % mod,
                    (9 * m * b + d * dg) % mod,
                )
            exact = (
                tri[0] % mod,
                tri[1] % mod,
                tri[2] % mod,
                (mat[2] // 9) % mod,
                mat[3] % mod,
            )
            if exact != state or V_of(tri, mat) % mod != (
                state[0] * state[4] ** 2
                - state[1] * state[4] * state[3]
                + 9 * state[2] * state[3] ** 2
            ) % mod:
                mismatches += 1
    if mismatches:
        raise SystemExit(f"path mismatches {mismatches}")
    say("level 9: 40 exact words of length 6 match the mod-16 and mod-36 reductions")


def main():
    say("python " + sys.version.split()[0])
    check_theorem_g()
    check_theorem_h()
    check_level9_sharpness()
    say("REFINEMENTS CERTIFICATE OK")
    (ROOT / "certificates" / "verify-refinements-log.txt").write_text(
        "\n".join(LOG) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
