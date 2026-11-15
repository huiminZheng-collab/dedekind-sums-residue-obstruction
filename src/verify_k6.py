"""Finite-sum input for (chi_3, chi_3, 6) on Gamma_1(9).

The periodic Bernoulli functions are the classical polynomials
B_1(x) = x - 1/2 and B_5(x) = x^5 - (5/2)x^4 + (5/3)x^3 - (1/6)x,
evaluated at the fractional part and set to 0 at integers.
The weight-4 value -2 on (1, 0; 9, 1) is a calibration against the
hand block sum, not an input taken from another log.
"""

from collections import deque
from fractions import Fraction
from math import gcd
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
LOG = []


def say(msg):
    print(msg, flush=True)
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


def chi3(n):
    r = n % 3
    return 0 if r == 0 else (1 if r == 1 else -1)


def frac_poly(coeffs, x):
    """coeffs[k] is the coefficient of x^k. Zero at integers."""
    x = Fraction(x)
    if x.denominator == 1:
        return Fraction(0)
    t = x - (x.numerator // x.denominator)
    val = Fraction(0)
    for c in reversed(coeffs):
        val = val * t + c
    return val


# Low-degree coefficient first.
B1 = [Fraction(-1, 2), Fraction(1)]
B3 = [Fraction(0), Fraction(1, 2), Fraction(-3, 2), Fraction(1)]
B5 = [
    Fraction(0),
    Fraction(-1, 6),
    Fraction(0),
    Fraction(5, 3),
    Fraction(-5, 2),
    Fraction(1),
]


def stilde(a, c, weight):
    """c^{weight-2} times the Theorem 1.7 sum for chi_3, chi_3."""
    if c < 0:
        sign = 1 if (weight - 2) % 2 == 0 else -1
        return sign * stilde(-a, -c, weight)
    if c == 0:
        return Fraction(0)
    a %= c
    key = (weight, a, c)
    if key in _SUM:
        return _SUM[key]
    poly = B3 if weight == 4 else B5
    total = Fraction(0)
    for j in range(c):
        cj = chi3(j)
        if cj == 0:
            continue
        b1 = frac_poly(B1, Fraction(j, c))
        base = Fraction(a * j, c)
        for n in (1, 2):
            total += chi3(n) * cj * b1 * frac_poly(poly, base + Fraction(n, 3))
    total *= c ** (weight - 2)
    _SUM[key] = total
    return total


_SUM = {}


def reduce_pair(num, den):
    if den < 0:
        num, den = -num, -den
    g = gcd(num, den)
    return num // g, den // g


def S_cusp(num, den):
    num, den = reduce_pair(num, den)
    if den % 9 or gcd(num, den) != 1:
        raise ValueError((num, den))
    return stilde(num, den, 6) / (den ** 4)


def h_at(gamma, alpha):
    a, b, c, d = gamma
    factor = c * alpha + d
    if factor == 0:
        raise ZeroDivisionError(gamma, alpha)
    num = a * alpha.numerator + b * alpha.denominator
    den = c * alpha.numerator + d * alpha.denominator
    return S_cusp(alpha.numerator, alpha.denominator) - factor ** 4 * S_cusp(num, den)


def interpolate(xs, ys):
    """Return coefficients of x^0, ..., x^{n-1}."""
    n = len(xs)
    A = []
    for x, y in zip(xs, ys):
        row = [Fraction(1)]
        for _ in range(n - 1):
            row.append(row[-1] * x)
        row.append(y)
        A.append(row)
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        inv = 1 / A[col][col]
        for j in range(col, n + 1):
            A[col][j] *= inv
        for r in range(n):
            if r == col or A[r][col] == 0:
                continue
            factor = A[r][col]
            for j in range(col, n + 1):
                A[r][j] -= factor * A[col][j]
    return [A[i][n] for i in range(n)]


def eval_low(coeffs, x):
    val = Fraction(0)
    for c in reversed(coeffs):
        val = val * x + c
    return val


def schreier_targets(path):
    N = None
    alph = []
    words = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if parts[0] == "N":
            N = int(parts[1])
        elif parts[0] == "A":
            alph.append(tuple(int(x) for x in parts[2:6]))
        elif parts[0] == "W":
            bar = parts.index("|")
            mat = tuple(int(x) for x in parts[1:bar])
            word = [int(x) for x in parts[bar + 1 :]]
            words.append((mat, word))
    S = (0, -1, 1, 0)
    T = (1, 1, 0, 1)
    X = [S, T, matinv(S), matinv(T)]
    reps = {(0, 1): (1, 0, 0, 1)}
    qq = deque([(0, 1)])
    while qq:
        r = reps[qq.popleft()]
        for g in X:
            rg = matmul(r, g)
            k = (rg[2] % N, rg[3] % N)
            if k not in reps:
                reps[k] = rg
                qq.append(k)
    sch = set()
    for r in reps.values():
        for g in X:
            rg = matmul(r, g)
            k = (rg[2] % N, rg[3] % N)
            h = matmul(rg, matinv(reps[k]))
            sch.add(h)
    produced = {(1, 0, 0, 1)}
    for mat, word in words:
        acc = (1, 0, 0, 1)
        for i in word:
            acc = matmul(acc, alph[i])
        if acc != mat or det(mat) != 1:
            raise SystemExit(f"word mismatch {mat} {acc}")
        produced.add(acc)
    if produced != sch or len(reps) != 72 or len(sch) != 19:
        raise SystemExit(f"generators {len(reps)} {len(sch)} {len(produced)}")
    return [mat for mat, _word in words]


def usable_cusps(gamma, pool):
    a, b, c, d = gamma
    good = []
    for alpha in pool:
        if c * alpha + d != 0:
            good.append(alpha)
    if len(good) < 7:
        raise SystemExit(f"not enough cusps for {gamma}")
    return good


def poly_of(gamma, pool):
    xs = usable_cusps(gamma, pool)
    ys = [h_at(gamma, x) for x in xs[:5]]
    coeffs = interpolate(xs[:5], ys)
    for x in xs[5:7]:
        if eval_low(coeffs, x) != h_at(gamma, x):
            raise SystemExit(f"interpolation failed at {x} for {gamma}")
    return coeffs


def main():
    say("python " + sys.version.split()[0])
    if stilde(1, 9, 4) != -2:
        raise SystemExit(f"weight-4 calibration failed: {stilde(1, 9, 4)}")
    say("calibration: Stilde_{chi3,chi3,4}(1,0;9,1) = -2")
    words = ROOT / "certificates" / "gamma1_9_words.txt"
    gens = schreier_targets(words)
    say(f"generators: 72 cosets, 19 Schreier elements, {len(gens)} stored words")
    pool = []
    seen = set()
    for t in range(1, 10):
        for num in range(1, 9 * t):
            alpha = Fraction(num, 9 * t)
            if alpha in seen or alpha.denominator % 9 or gcd(alpha.numerator, alpha.denominator) != 1:
                continue
            seen.add(alpha)
            pool.append(alpha)
            if len(pool) == 12:
                break
        if len(pool) == 12:
            break
    rows = []
    content = 0
    for gamma in gens:
        low = poly_of(gamma, pool)
        # a_i is the coefficient of x^{4-i}.
        ai = [low[4 - i] for i in range(5)]
        raw = []
        b = []
        for i, a in enumerate(ai):
            cleared = (9 ** (i + 1)) * a
            if cleared.denominator != 1:
                raise SystemExit(f"not integral after 9^(i+1): {gamma} {i} {cleared}")
            raw.append(int(cleared))
            if int(cleared) % 30:
                raise SystemExit(f"lambda=30 fails: {gamma} {i} {cleared}")
            b.append(int(cleared) // 30)
            content = gcd(content, abs(int(cleared)))
        a, _bb, c, d = gamma
        r = -c // 9
        if any(x % 3 for x in b[2:]):
            raise SystemExit(f"higher coefficient not 0 mod 3: {gamma} {b}")
        if b[0] % 3 != (2 * r * r) % 3 or b[1] % 3 != (-4 * r) % 3:
            raise SystemExit(f"quadratic law fails: {gamma} r={r} b={b}")
        # Formula (12) on the inverse, compared with the finite sum.
        if c == 0:
            predicted = Fraction(0)
        else:
            predicted = sum(ai[i] * (d ** (4 - i)) * ((-c) ** i) for i in range(5))
        direct = stilde(d, -c, 6)
        if predicted != direct:
            raise SystemExit(f"evaluation mismatch {gamma}: {predicted} != {direct}")
        rows.append((gamma, b, r, predicted))
        say(f"g={list(gamma)} b={b} r={r} U={3 * predicted / 10}")
    if content != 30:
        raise SystemExit(f"global content {content} != 30")
    say("content of all 9^(i+1) a_i is 30; lambda=30, alpha=2 on all 18 words")
    p = stilde(1, 9, 6)
    p2 = stilde(1, 18, 6)
    if p != Fraction(310, 3) or p2 != Fraction(4480, 3):
        raise SystemExit(f"witnesses {p} {p2}")
    if gcd(31, 448) != 1:
        raise SystemExit("gcd")
    say("direct sums: Stilde(P)=310/3, Stilde(P^2)=4480/3, gcd(31,448)=1")
    # The polynomial of P^{-1}=(1,0;-9,1) is read back from its b-row.
    target = (1, 0, -9, 1)
    b = next(row[1] for row in rows if row[0] == target)
    if b != [-31, 116, -162, 108, 0]:
        raise SystemExit(f"P inverse coefficients {b}")
    for m in range(1, 31):
        acc = sum(-31 + 116 * t - 162 * t * t + 108 * t ** 3 for t in range(1, m + 1))
        if acc != m * m * (27 * m * m + 4):
            raise SystemExit(f"power sum failed at m={m}")
    for m in range(1, 7):
        got = stilde(1, 9 * m, 6)
        expect = Fraction(10, 3) * m * m * (27 * m * m + 4)
        if got != expect:
            raise SystemExit(f"family m={m}: {got} != {expect}")
    say("power sums m=1..30 and direct sums m=1..6 match (10/3)m^2(27m^2+4)")
    say("K6 CERTIFICATE OK")
    out = ROOT / "certificates" / "verify-k6-log.txt"
    out.write_text("\n".join(LOG) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
