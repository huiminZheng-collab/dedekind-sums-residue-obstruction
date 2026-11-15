"""Certificate for the second set-equality failure, (chi_3, chi_4, k=4) on Gamma_1(12).

Finite pieces checked here:
  * author period polynomials agree with our finite sum at >= 3 cusps;
  * formula (12) agrees with our finite sum;
  * the integer coordinates (A, B, C) satisfy the mod-3 invariant;
  * stored Schreier words multiply back to every nontrivial Schreier generator;
  * the two witness values are -16/3 and -2254/3, with rational gcd 2/3.

The cocycle preservation of the invariant is the written proof, not this file.
Standard library only.
"""

from fractions import Fraction
from math import comb, gcd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_words import check_file

ROOT = Path(__file__).resolve().parents[1]
LOG = []


def say(msg):
    print(msg)
    LOG.append(msg)


def matmul(A, B):
    return (
        A[0] * B[0] + A[1] * B[2],
        A[0] * B[1] + A[1] * B[3],
        A[2] * B[0] + A[3] * B[2],
        A[2] * B[1] + A[3] * B[3],
    )


def matinv(A):
    return (A[3], -A[1], -A[2], A[0])


def det(A):
    return A[0] * A[3] - A[1] * A[2]


def bernoulli_coeffs(k):
    coeffs = [Fraction(0)] * (k + 1)
    for m in range(k + 1):
        for n in range(m + 1):
            sign = -1 if n % 2 else 1
            factor = Fraction(sign * comb(m, n), m + 1)
            for r in range(k + 1):
                coeffs[r] += factor * comb(k, r) * (n ** (k - r))
    return coeffs


B1 = bernoulli_coeffs(1)
B3 = bernoulli_coeffs(3)


def bern(coeffs, x):
    x = Fraction(x)
    if x.denominator == 1:
        return Fraction(0)
    t = x - (x.numerator // x.denominator)
    val = Fraction(0)
    for c in reversed(coeffs):
        val = val * t + c
    return val


def chi3(n):
    r = n % 3
    return 0 if r == 0 else (1 if r == 1 else -1)


def chi4(n):
    r = n % 4
    if r % 2 == 0:
        return 0
    return 1 if r == 1 else -1


def stilde(a, c):
    """S-tilde for chi_3, chi_4, k=4. Even weight, so negative c uses -gamma."""
    if c < 0:
        return stilde(-a, -c)
    if c == 0:
        return Fraction(0)
    total = Fraction(0)
    for j in range(c):
        cj = chi4(j)
        if cj == 0:
            continue
        b1 = bern(B1, Fraction(j, c))
        for n in (1, 2):
            total += chi3(n) * cj * b1 * bern(B3, Fraction(a * j, c) + Fraction(n, 3))
    return c * c * total


# Author file hF-k4,3c3;2-1;4c4;3-1, copied and checked below.
LINES = [
    ((1, 1, 0, 1), {2: Fraction(0), 1: Fraction(0), 0: Fraction(0)}),
    ((-11, 1, -12, 1), {2: Fraction(16, 3), 1: Fraction(-8, 9), 0: Fraction(0)}),
    ((49, -9, 60, -11), {2: Fraction(1120, 3), 1: Fraction(-1232, 9), 0: Fraction(112, 9)}),
    ((13, -3, 48, -11), {2: Fraction(-2174, 3), 1: Fraction(2852, 9), 0: Fraction(-314, 9)}),
    ((13, -4, 36, -11), {2: Fraction(-48), 1: Fraction(80, 3), 0: Fraction(-32, 9)}),
    ((61, -25, 144, -59), {2: Fraction(-1200), 1: Fraction(2948, 3), 0: Fraction(-1810, 9)}),
    ((13, -6, 24, -11), {2: Fraction(-80, 3), 1: Fraction(220, 9), 0: Fraction(-50, 9)}),
    ((25, -16, 36, -23), {2: Fraction(-144), 1: Fraction(560, 3), 0: Fraction(-544, 9)}),
    ((37, -27, 48, -35), {2: Fraction(2254, 3), 1: Fraction(-9724, 9), 0: Fraction(1166, 3)}),
]


def abc_of(coef):
    A = coef[2] * Fraction(3, 2)
    B = coef[1] * Fraction(9, 2)
    C = coef[0] * Fraction(27, 2)
    if A.denominator != 1 or B.denominator != 1 or C.denominator != 1:
        raise SystemExit(f"nonintegral coordinates {coef}")
    return int(A), int(B), int(C)


def reduce_pair(num, den):
    if den < 0:
        num, den = -num, -den
    g = gcd(num, den)
    return num // g, den // g


def S_cusp(num, den):
    num, den = reduce_pair(num, den)
    if den % 12 or gcd(num, den) != 1:
        raise ValueError((num, den))
    return stilde(num, den) / (den * den)


def check_polynomials_and_invariant():
    cusps = [(1, 12), (5, 12), (7, 12), (11, 12), (5, 24), (7, 24)]
    for mat, coef in LINES:
        A, B, C = abc_of(coef)
        cpp = mat[2] // 3
        if C % 3 != 0 or B % 3 != cpp % 3 or A % 3 != (2 * cpp * cpp) % 3:
            raise SystemExit(f"invariant fails {mat}")
        if mat == (1, 1, 0, 1):
            continue
        a, b, c, d = mat
        ncmp = 0
        for p, q in cusps:
            alpha = Fraction(p, q)
            factor = c * alpha + d
            if factor == 0:
                continue
            num = a * p + b * q
            den = c * p + d * q
            got = S_cusp(p, q) - factor ** 2 * S_cusp(num, den)
            claimed = coef[2] * alpha * alpha + coef[1] * alpha + coef[0]
            if got != claimed:
                raise SystemExit(f"mismatch {mat} at {alpha}: {got} vs {claimed}")
            ncmp += 1
        if ncmp < 3:
            raise SystemExit(f"too few cusps {mat}")
        # formula (12): S-tilde(inverse) = c2 d^2 - c1 d c + c0 c^2
        pred = coef[2] * d * d - coef[1] * d * c + coef[0] * c * c
        inv = matinv(mat)
        got = stilde(inv[0], inv[2])
        if pred != got:
            raise SystemExit(f"formula (12) {mat}: {pred} vs {got}")
    say("polynomials: 8 generators match our finite sum at >= 3 cusps; (12) matches; mod-3 invariant holds")


def rgcd(x, y):
    x, y = Fraction(x), Fraction(y)
    return abs(Fraction(gcd(x.numerator * y.denominator, y.numerator * x.denominator), x.denominator * y.denominator))


def check_witnesses():
    w1 = (1, -1, 12, -11)
    w2 = (-11, 3, -48, 13)
    for w in (w1, w2):
        if det(w) != 1 or w[0] % 12 != 1 or w[3] % 12 != 1 or w[2] % 12 != 0:
            raise SystemExit(f"illegal witness {w}")
    v1 = stilde(w1[0], w1[2])
    v2 = stilde(w2[0], w2[2])
    if v1 != Fraction(-16, 3) or v2 != Fraction(-2254, 3):
        raise SystemExit((v1, v2))
    if rgcd(v1, v2) != Fraction(2, 3):
        raise SystemExit("gcd")
    # Both lie in the allowed classes: (3/2) S-tilde ≢ 2 (mod 3).
    for v in (v1, v2):
        n = v * Fraction(3, 2)
        if n.denominator != 1 or int(n) % 3 == 2:
            raise SystemExit(f"witness in a forbidden class {v}")
    say("witnesses: S-tilde(1,-1;12,-11)=-16/3, S-tilde(-11,3;-48,13)=-2254/3, gcd 2/3")


def check_cocycle_sample():
    """One product, cocycle evaluation against the finite sum. A sample, not the proof."""
    g1 = (-11, 1, -12, 1)
    g2 = (13, -6, 24, -11)
    coef = {M: C for M, C in LINES}
    # h_{g1 g2} = h_g1 | g2 + h_g2, then formula (12) on the product.
    c1 = coef[g1]
    c2 = coef[g2]
    a, b, c, d = g2

    def slash_poly(coef, g):
        A, B, C = coef[2], coef[1], coef[0]
        ga, gb, gc, gd = g
        nA = A * ga * ga + B * ga * gc + C * gc * gc
        nB = 2 * A * ga * gb + B * (ga * gd + gb * gc) + 2 * C * gc * gd
        nC = A * gb * gb + B * gb * gd + C * gd * gd
        return {2: nA, 1: nB, 0: nC}

    sl = slash_poly(c1, g2)
    h = {k: sl[k] + c2[k] for k in (0, 1, 2)}
    prod = matmul(g1, g2)
    # S-tilde(prod^{-1}) from h_prod
    pa, pb, pc, pd = prod
    pred = h[2] * pd * pd - h[1] * pd * pc + h[0] * pc * pc
    inv = matinv(prod)
    got = stilde(inv[0], inv[2])
    if pred != got:
        raise SystemExit(f"cocycle sample {pred} vs {got}")
    say(f"sample: cocycle on one product agrees with the finite sum, value {got}")


def check_words():
    info = check_file(
        ROOT / "certificates" / "gamma1_12_words.txt",
        expected_cosets=96,
        expected_schreier=27,
    )
    say(
        f"words: {info['cosets']} cosets, {info['schreier']} Schreier generators, "
        f"{info['nontrivial_words']} explicit words multiplied back, max length {info['max_word_length']}"
    )


def main():
    say("python " + sys.version.split()[0])
    if B1 != [Fraction(-1, 2), Fraction(1)]:
        raise SystemExit("B1")
    if B3 != [Fraction(0), Fraction(1, 2), Fraction(-3, 2), Fraction(1)]:
        raise SystemExit("B3")
    check_polynomials_and_invariant()
    check_witnesses()
    check_cocycle_sample()
    check_words()
    say("CHI34 CERTIFICATE OK")
    (ROOT / "certificates" / "verify-chi34-log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
