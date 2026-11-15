"""Certificate for S-tilde_{chi_3, chi_3, 4} on Gamma_1(9).

The group-wide congruence is the coefficient invariant in the proof:
even integral coordinates (A, B1, C1), and
    A ≡ 2 m^2 (mod 3),  B1 ≡ -2 m (mod 3),
preserved by the weight -2 cocycle. This program checks the finite pieces
that the proof quotes:

  * the six quadratics agree with the finite sum at three or more cusps;
  * formula (12) agrees with the finite sum on those matrices;
  * the seven triples satisfy the invariant;
  * the -2 m^2 family, by expanding the eight blocks and summing with
    Faulhaber's formulae (no external computer-algebra system);
  * each nontrivial Schreier generator equals the product of its stored word.

A 1458-state residue automaton is run afterwards as a cross-check of the
same congruence. It is not the argument in the proof.

Standard library only. From the repository root:

    python src/verify.py
    python src/check_words.py certificates/gamma1_9_words.txt
"""

from collections import deque
from fractions import Fraction
from math import comb, gcd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_words import check_file

CHI = (0, 1, -1)
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
    a, b, c, d = A
    return (d, -b, -c, a)


def det(A):
    return A[0] * A[3] - A[1] * A[2]


def bernoulli_coeffs(k):
    """Definition 1.1, expanded as a polynomial in the fractional part."""
    coeffs = [Fraction(0)] * (k + 1)
    for m in range(k + 1):
        for n in range(m + 1):
            sign = -1 if n % 2 else 1
            factor = Fraction(sign * comb(m, n), m + 1)
            for r in range(k + 1):
                coeffs[r] += factor * comb(k, r) * (n ** (k - r))
    return coeffs


def stilde_bernoulli(a, c):
    """Theorem 1.7 times c^2, for chi_3, chi_3, k=4, c>0."""
    b1c = bernoulli_coeffs(1)
    b3c = bernoulli_coeffs(3)
    total = Fraction(0)
    for j in range(c):
        cj = CHI[j % 3]
        if cj == 0:
            continue
        b1 = _eval_bernoulli(b1c, Fraction(j, c))
        for n in (1, 2):
            total += CHI[n] * cj * b1 * _eval_bernoulli(b3c, Fraction(a * j, c) + Fraction(n, 3))
    return c * c * total


def _eval_bernoulli(coeffs, x):
    x = Fraction(x)
    if x.denominator == 1:
        return Fraction(0)
    t = x - (x.numerator // x.denominator)
    val = Fraction(0)
    for c in reversed(coeffs):
        val = val * t + c
    return val


def stilde(a, c):
    """Same value as stilde_bernoulli, cleared to an integer formula. c>0."""
    if c <= 0:
        raise ValueError(c)
    a %= c
    num = 0
    for j in range(1, c):
        ch2 = CHI[j % 3]
        if ch2 == 0:
            continue
        twoj = 2 * j - c
        for n in (1, 2):
            r = (3 * (a * j) + n * c) % (3 * c)
            if r == 0:
                continue
            num += CHI[n] * ch2 * twoj * r * (2 * r - 3 * c) * (r - 3 * c)
    den = 108 * c * c
    if num % den:
        raise ValueError("nonintegral value")
    return num // den


def slash(triple, g):
    """Weight -2 slash on (A, 9B, 9C)."""
    A, B1, C1 = (Fraction(x) for x in triple)
    a, b, c, d = g
    B = B1 / 9
    C = C1 / 9
    nA = A * a * a + B * a * c + C * c * c
    nB = 2 * A * a * b + B * (a * d + b * c) + 2 * C * c * d
    nC = A * b * b + B * b * d + C * d * d
    return (nA, 9 * nB, 9 * nC)


def stilde_of_inverse(tri, mat):
    """Formula (12): tri stores h_mat, return S-tilde(mat^{-1})."""
    A, B1, C1 = (Fraction(x) for x in tri)
    a, b, c, d = mat
    if c == 0:
        return Fraction(0)
    return A * d * d - (B1 * d * c) / 9 + (C1 * c * c) / 9


# h = A x^2 + (B1/9) x + C1/9, copied from the author's hF-k4 file and
# checked against the finite sum below.
RAW = {
    (1, 1, 0, 1): (0, 0, 0),
    (-8, 1, -9, 1): (2, -4, 0),
    (19, -4, 81, -17): (522, -1980, 208),
    (10, -3, 27, -8): (-174, 876, -124),
    (-26, 11, -45, 19): (-220, 1672, -352),
    (37, -16, 81, -35): (-612, 4752, -1024),
    (19, -12, 27, -17): (156, -1716, 526),
}


def poly_of(tri):
    A, B1, C1 = tri
    return lambda x: A * x * x + Fraction(B1, 9) * x + Fraction(C1, 9)


def S_cusp(num, den):
    g = gcd(num, den)
    num //= g
    den //= g
    if den < 0:
        num, den = -num, -den
    if den % 9 or gcd(num, den) != 1:
        raise ValueError((num, den))
    return Fraction(stilde(num, den), den * den)


def h_at(gamma, alpha):
    a, b, c, d = gamma
    factor = c * alpha + d
    if factor == 0:
        raise ZeroDivisionError
    beta_n = a * alpha.numerator + b * alpha.denominator
    beta_d = c * alpha.numerator + d * alpha.denominator
    return S_cusp(alpha.numerator, alpha.denominator) - factor ** 2 * S_cusp(beta_n, beta_d)


def check_bernoulli_agreement():
    b1 = bernoulli_coeffs(1)
    b3 = bernoulli_coeffs(3)
    if b1 != [Fraction(-1, 2), Fraction(1)]:
        raise SystemExit(f"B1 coeffs {b1}")
    if b3 != [Fraction(0), Fraction(1, 2), Fraction(-3, 2), Fraction(1)]:
        raise SystemExit(f"B3 coeffs {b3}")
    samples = [(1, 9), (1, 18), (10, 27), (19, 81)]
    for a, c in samples:
        if stilde(a, c) != stilde_bernoulli(a, c):
            raise SystemExit(f"integer formula disagrees with Definition 1.1 at {(a, c)}")
    say("bernoulli: B1 and B3 match Definition 1.1, and the cleared formula matches on 4 samples")


def check_polynomials():
    cusps = [Fraction(2, 9), Fraction(4, 9), Fraction(5, 9), Fraction(7, 9), Fraction(8, 9)]
    for g, tri in RAW.items():
        if g == (1, 1, 0, 1):
            continue
        claimed = poly_of(tri)
        ncmp = 0
        for alpha in cusps:
            try:
                got = h_at(g, alpha)
            except (ValueError, ZeroDivisionError):
                continue
            if got != claimed(alpha):
                raise SystemExit(f"polynomial mismatch {g} at {alpha}: {got} vs {claimed(alpha)}")
            ncmp += 1
        if ncmp < 3:
            raise SystemExit(f"too few cusps for {g}")
        pred = stilde_of_inverse(tri, g)
        inv = matinv(g)
        if inv[2] > 0:
            got = stilde(inv[0], inv[2])
        elif inv[2] < 0:
            neg = tuple(-x for x in inv)
            # k-2 = 2 is even, so S-tilde(gamma) = S-tilde(-gamma).
            got = stilde(neg[0], neg[2])
        else:
            got = 0
        if pred != got:
            raise SystemExit(f"formula (12) mismatch {g}: {pred} vs {got}")
    say("polynomials: 6 generators match the finite sum at >= 3 cusps, and (12) matches")


def check_invariant_on_generators():
    for g, tri in RAW.items():
        A, B1, C1 = tri
        m = g[2] // 9
        if A % 2 or B1 % 2 or C1 % 2:
            raise SystemExit(f"not even {g} {tri}")
        if (A - 2 * m * m) % 3 or (B1 + 2 * m) % 3:
            raise SystemExit(f"mod 3 fails {g}")
        inv = matinv(g)
        sl = slash(tri, inv)
        inv_tri = tuple(-int(x) for x in sl)
        mi = inv[2] // 9
        Ai, Bi, Ci = inv_tri
        if Ai % 2 or Bi % 2 or Ci % 2 or (Ai - 2 * mi * mi) % 3 or (Bi + 2 * mi) % 3:
            raise SystemExit(f"inverse invariant fails {inv} {inv_tri}")
    say("invariant: 7 generators and their inverses have even (A, B1, C1), A≡2m^2 and B1≡-2m mod 3")


# --- family: bipolynomials in (m, q), then Faulhaber in q ---

def _trim(p):
    return {k: v for k, v in p.items() if v}


def mul_bi(p, q):
    r = {}
    for (i, j), a in p.items():
        for (k, l), b in q.items():
            key = (i + k, j + l)
            r[key] = r.get(key, Fraction(0)) + a * b
    return _trim(r)


def add_bi(p, q, scale=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, Fraction(0)) + scale * v
    return _trim(r)


def mul_uni(p, q):
    r = {}
    for i, a in p.items():
        for j, b in q.items():
            r[i + j] = r.get(i + j, Fraction(0)) + a * b
    return {k: v for k, v in r.items() if v}


def add_uni(p, q, scale=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, Fraction(0)) + scale * v
    return {k: v for k, v in r.items() if v}


def pow_sums_0(N):
    """sum_{q=0}^{N-1} q^k for k=0..4, as polynomials in m. N is a polynomial in m."""
    one = {0: Fraction(-1)}
    Lm1 = add_uni(N, one)
    twoLm1 = add_uni({k: 2 * v for k, v in N.items()}, one)
    L2 = mul_uni(N, N)
    cubic = add_uni(add_uni({k: 3 * v for k, v in L2.items()}, {k: -3 * v for k, v in N.items()}), one)
    half = mul_uni(N, Lm1)
    S0 = N
    S1 = {k: v / 2 for k, v in half.items()}
    S2 = {k: v / 6 for k, v in mul_uni(half, twoLm1).items()}
    S3 = {k: v / 4 for k, v in mul_uni(half, half).items()}
    S4 = {k: v / 30 for k, v in mul_uni(mul_uni(half, twoLm1), cubic).items()}
    return [S0, S1, S2, S3, S4]


def sums_inclusive(q0, q1):
    hi = pow_sums_0(add_uni(q1, {0: Fraction(1)}))
    lo = pow_sums_0(q0)
    return [add_uni(hi[k], lo[k], -1) for k in range(5)]


def block_numerator(r, n, floor):
    """Numerator N(m,q) of the contribution 81 m^2 chi B1 B3 = N / (324 m^2)."""
    # B1 = (6q + 2r - 9m) / (18m)
    # f = (3q + 3 n m + r - 9 m floor) / (9m)
    # f-1/2 = (6q + 6 n m + 2r - 18 m floor - 9m) / (18m)
    # f-1 = (3q + 3 n m + r - 9 m floor - 9m) / (9m)
    def lin(m_coeff, q_coeff, const):
        p = {}
        if m_coeff:
            p[(1, 0)] = Fraction(m_coeff)
        if q_coeff:
            p[(0, 1)] = Fraction(q_coeff)
        if const:
            p[(0, 0)] = Fraction(const)
        return p

    n1 = lin(-9, 6, 2 * r)
    nf = lin(3 * n - 9 * floor, 3, r)
    nh = lin(6 * n - 18 * floor - 9, 6, 2 * r)
    nm = lin(3 * n - 9 * floor - 9, 3, r)
    chi = (1 if r == 1 else -1) * (1 if n == 1 else -1)
    prod = mul_bi(mul_bi(n1, nf), mul_bi(nh, nm))
    if chi == -1:
        prod = {k: -v for k, v in prod.items()}
    return prod


def sum_block(r, n, floor, q0, q1):
    numer = block_numerator(r, n, floor)
    by_q = {}
    for (i, j), coeff in numer.items():
        by_q.setdefault(j, {})
        by_q[j][i] = by_q[j].get(i, Fraction(0)) + coeff
    Sk = sums_inclusive(q0, q1)
    total = {}
    for j, poly in by_q.items():
        total = add_uni(total, mul_uni(poly, Sk[j]))
    # contribution sum = (sum of numerators) / (324 m^2). Leave the m^2 in place
    # by returning the numerator polynomial; the caller divides.
    return total


def check_family():
    # Four blocks with r=1. The proof pairs each with a block of equal sum.
    blocks = [
        (1, 1, 0, {0: Fraction(0)}, {1: Fraction(2), 0: Fraction(-1)}),
        (1, 1, 1, {1: Fraction(2)}, {1: Fraction(3), 0: Fraction(-1)}),
        (1, 2, 0, {0: Fraction(0)}, {1: Fraction(1), 0: Fraction(-1)}),
        (1, 2, 1, {1: Fraction(1)}, {1: Fraction(3), 0: Fraction(-1)}),
    ]
    numer_sum = {}
    for spec in blocks:
        numer_sum = add_uni(numer_sum, sum_block(*spec))
    # symmetry doubles the sum; contribution = numer / (324 m^2)
    # 2 * numer / (324 m^2) = -2 m^2, so numer = -2 m^2 * 324 m^2 / 2 = -324 m^4
    expect = {4: Fraction(-324)}
    if numer_sum != expect:
        raise SystemExit(f"family numerator {numer_sum}")
    # The unsplit double sum, for small m, against the same closed form.
    for m in range(1, 13):
        if stilde(1, 9 * m) != -2 * m * m:
            raise SystemExit(f"direct sum fails m={m}")
    # Pointwise symmetry of summands, checked on a finite grid and proved in the text.
    for m in range(1, 5):
        for q in range(0, 3 * m):
            for r, n, floor in ((1, 1, 0), (1, 1, 1), (1, 2, 0), (1, 2, 1)):
                # only where this (q, floor) is the true floor
                boundary = m * (3 - n)
                true_floor = 1 if q >= boundary else 0
                if floor != true_floor:
                    continue
                left = summand_value(m, q, r, n)
                right = summand_value(m, 3 * m - 1 - q, 3 - r, 3 - n)
                if left != right:
                    raise SystemExit(f"symmetry fail m={m} q={q}")
    say("family: four Faulhaber blocks, doubled by the (j,n)->(c-j,3-n) symmetry, sum to -2*m**2")
    say("family direct: unsplit double sum equals -2*m**2 for m=1..12")


def summand_value(m, q, r, n):
    """One term of S-tilde, Bernoulli form, for the family a=1, c=9m."""
    c = 9 * m
    j = 3 * q + r
    x = Fraction(j, c) + Fraction(n, 3)
    # B3 via the polynomial, B1 via the polynomial; chi included.
    b1 = _eval_bernoulli(bernoulli_coeffs(1), Fraction(j, c))
    b3 = _eval_bernoulli(bernoulli_coeffs(3), x)
    chi = (1 if r == 1 else -1) * (1 if n == 1 else -1)
    return 81 * m * m * chi * b1 * b3


def check_words():
    info = check_file(
        ROOT / "certificates" / "gamma1_9_words.txt",
        expected_cosets=72,
        expected_schreier=19,
    )
    say(
        f"words: {info['cosets']} cosets, {info['schreier']} Schreier generators, "
        f"{info['nontrivial_words']} explicit words multiplied back, max length {info['max_word_length']}"
    )


def h_table():
    table = {}
    for g, tri in RAW.items():
        table[g] = tri
        sl = slash(tri, matinv(g))
        inv_tri = tuple(-x for x in sl)
        if any(x != int(x) for x in inv_tri):
            raise SystemExit(f"inverse polynomial not integral for {g}")
        table[matinv(g)] = tuple(int(x) for x in inv_tri)
    return table


def check_automaton_crosscheck():
    """Cross-check only. The proof does not use this graph."""
    table = h_table()
    gens = list(table.items())
    mod_m, mod_c = 54, 6

    def red_mat(M):
        return tuple(x % mod_m for x in M)

    def red_tri(tri):
        return tuple(int(Fraction(x)) % mod_c for x in tri)

    seen = {(red_mat((1, 0, 0, 1)), (0, 0, 0))}
    qq = deque(seen)
    checked = 0
    while qq:
        rm, rtri = qq.popleft()
        a, b, c, d = rm
        if c % 9 or a % 9 != 1 or d % 9 != 1 or (a * d - b * c - 1) % mod_m:
            raise SystemExit(f"illegal class {rm}")
        A, B1, C1 = rtri
        m6 = (c // 9) % 6
        d6 = d % 6
        val = (A * d6 * d6 - B1 * d6 * m6 + 3 * C1 * m6 * m6) % 6
        if val != (-2 * m6 * m6) % 6:
            raise SystemExit(f"cross-check fail {rm} {rtri}")
        # The hand invariant, reduced modulo 6, on every reached state.
        if A % 2 or B1 % 2 or C1 % 2:
            raise SystemExit(f"parity fail in automaton {rm} {rtri}")
        if (A - 2 * m6 * m6) % 3 or (B1 + 2 * m6) % 3:
            raise SystemExit(f"mod 3 fail in automaton {rm} {rtri}")
        checked += 1
        for g, gtri in gens:
            nm = red_mat(matmul(rm, g))
            sl = slash(rtri, g)
            ntri = red_tri((sl[0] + gtri[0], sl[1] + gtri[1], sl[2] + gtri[2]))
            st = (nm, ntri)
            if st not in seen:
                seen.add(st)
                qq.append(st)
    say(f"cross-check: {checked} residue states, invariant and the mod-6 congruence both hold")


def main():
    say("python " + sys.version.split()[0])
    check_bernoulli_agreement()
    check_polynomials()
    check_invariant_on_generators()
    check_family()
    check_words()
    check_automaton_crosscheck()
    say("CERTIFICATE OK")
    out = ROOT / "certificates" / "verify-log.txt"
    out.write_text("\n".join(LOG) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
