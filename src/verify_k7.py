"""Containment certificate for S-tilde_{chi_3, chi_5, 7}(Gamma_1(15)).

The Schreier generators of Gamma_1(15) are computed from S and T.
For each non-identity generator the period polynomial is interpolated
from the finite sum at six cusps and checked at a seventh. Membership in
P(7; 4, 3) is the input to the slash-closure argument in PAPER.en.md.
The two explicit values with gcd 4/3 are recomputed here.

Standard library only. From the repository root:

    python src/verify_k7.py
"""

from fractions import Fraction
from math import comb, gcd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_words import schreier_generators

ROOT = Path(__file__).resolve().parents[1]
LOG = []
_BERN = {}
_SUM = {}


def say(msg):
    print(msg, flush=True)
    LOG.append(msg)


def bernoulli_coeffs(k):
    coeffs = [Fraction(0)] * (k + 1)
    for m in range(k + 1):
        for n in range(m + 1):
            sign = -1 if n % 2 else 1
            factor = Fraction(sign * comb(m, n), m + 1)
            for r in range(k + 1):
                coeffs[r] += factor * comb(k, r) * (n ** (k - r))
    return coeffs


def B(k, x):
    x = Fraction(x)
    if x.denominator == 1:
        return Fraction(0)
    coeffs = _BERN.get(k)
    if coeffs is None:
        coeffs = bernoulli_coeffs(k)
        _BERN[k] = coeffs
    t = x - (x.numerator // x.denominator)
    val = Fraction(0)
    for c in reversed(coeffs):
        val = val * t + c
    return val


def chi3(n):
    r = n % 3
    return 0 if r == 0 else (1 if r == 1 else -1)


def chi5(n):
    a = n % 5
    if a == 0:
        return 0
    return 1 if pow(a, 2, 5) == 1 else -1


def stilde(a, c):
    """S-tilde for chi_3, chi_5, k=7. k-2 is odd, so S-tilde(-gamma) = -S-tilde(gamma)."""
    if c < 0:
        return -stilde(-a, -c)
    if c == 0:
        return Fraction(0)
    a %= c
    key = (a, c)
    if key in _SUM:
        return _SUM[key]
    total = Fraction(0)
    for j in range(c):
        cj = chi5(j)
        if cj == 0:
            continue
        b1 = B(1, Fraction(j, c))
        aj = Fraction(a * j, c)
        for n in (1, 2):
            total += chi3(n) * cj * b1 * B(6, aj + Fraction(n, 3))
    total *= c ** 5
    _SUM[key] = total
    return total


def rgcd(x, y):
    x, y = Fraction(x), Fraction(y)
    return abs(
        Fraction(
            gcd(x.numerator * y.denominator, y.numerator * x.denominator),
            x.denominator * y.denominator,
        )
    )


def reduce_pair(num, den):
    if den < 0:
        num, den = -num, -den
    g = gcd(num, den)
    return num // g, den // g


def S_cusp(num, den):
    num, den = reduce_pair(num, den)
    if den % 15 or gcd(num, den) != 1:
        raise ValueError((num, den))
    return stilde(num, den) / (den ** 5)


def h_at(gamma, alpha):
    a, b, c, d = gamma
    factor = c * alpha + d
    if factor == 0:
        raise ZeroDivisionError
    num = a * alpha.numerator + b * alpha.denominator
    den = c * alpha.numerator + d * alpha.denominator
    return S_cusp(alpha.numerator, alpha.denominator) - factor ** 5 * S_cusp(num, den)


def cusp_list():
    out = []
    for den in (15, 30, 45, 60, 75, 90):
        for num in range(1, den):
            if gcd(num, den) == 1:
                out.append(Fraction(num, den))
    return out


def interpolate(xs, ys):
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
        inv = Fraction(1, 1) / A[col][col]
        for j in range(col, n + 1):
            A[col][j] *= inv
        for r in range(n):
            if r == col or A[r][col] == 0:
                continue
            factor = A[r][col]
            for j in range(col, n + 1):
                A[r][j] -= factor * A[col][j]
    return [A[i][n] for i in range(n)]


def eval_poly(coeffs, x):
    val = Fraction(0)
    for c in reversed(coeffs):
        val = val * x + c
    return val


def cleared_integers(coeffs):
    """3^{n+1} a_n / 4, where a_n is the coefficient of x^{5-n}."""
    out = []
    for n in range(6):
        an = coeffs[5 - n]
        cleared = (3 ** (n + 1)) * an
        if cleared.denominator != 1 or int(cleared) % 4:
            raise SystemExit(f"not in P(7;4,3): n={n}, 3^(n+1)a_n={cleared}")
        out.append(int(cleared) // 4)
    return out


def admissible_bottom_rows(N):
    return sum(1 for c in range(N) for d in range(N) if gcd(gcd(c, d), N) == 1)


def main():
    say("python " + sys.version.split()[0])
    b6 = bernoulli_coeffs(6)
    expected = [
        Fraction(1, 42),
        Fraction(0),
        Fraction(-1, 2),
        Fraction(0),
        Fraction(5, 2),
        Fraction(-3),
        Fraction(1),
    ]
    if b6 != expected:
        raise SystemExit(f"B6 coeffs {b6}")
    say("bernoulli: B6 matches x^6 - 3x^5 + (5/2)x^4 - (1/2)x^2 + 1/42")

    m1 = (1, -1, 15, -14)
    if m1[0] * m1[3] - m1[1] * m1[2] != 1 or m1[0] % 15 != 1 or m1[3] % 15 != 1 or m1[2] % 15:
        raise SystemExit("witness not in Gamma_1(15)")
    if not (1 <= m1[0] < 750 and 15 <= m1[2] < 750):
        raise SystemExit("witness not in G_50(15)")
    v1 = stilde(1, 15)
    if v1 != Fraction(-17000, 3):
        raise SystemExit(v1)
    say(f"G_50 witness S-tilde(1,-1;15,-14) = {v1}")

    m3 = (16, -9, 105, -59)
    if m3[0] * m3[3] - m3[1] * m3[2] != 1 or m3[0] % 15 != 1 or m3[3] % 15 != 1 or m3[2] % 15:
        raise SystemExit("second witness illegal")
    v3 = stilde(16, 105)
    g = rgcd(v1, v3)
    if g != Fraction(4, 3):
        raise SystemExit(g)
    say(f"second value S-tilde(16,-9;105,-59) = {v3}; gcd = {g}")

    author = (-14, 1, -15, 1)
    author_h = [
        Fraction(0),
        Fraction(0),
        Fraction(-160, 27),
        Fraction(400, 3),
        Fraction(-12800, 9),
        Fraction(17000, 3),
    ]
    got = h_at(author, Fraction(2, 15))
    claimed = eval_poly(author_h, Fraction(2, 15))
    if got != claimed:
        raise SystemExit(f"author polynomial disagrees at 2/15: {got} vs {claimed}")
    say("author polynomial of (-14,1;-15,1) matches the finite sum at 2/15")

    rows = admissible_bottom_rows(15)
    cosets, sch = schreier_generators(15)
    if cosets != rows or cosets != 192 or len(sch) != 47:
        raise SystemExit(f"cosets {cosets}, admissible {rows}, schreier {len(sch)}")
    say(f"generation: {cosets} cosets, {len(sch)} Schreier generators, admissible bottom rows {rows}")

    cusps = cusp_list()
    failures = 0
    max_den = 1
    for mat in sorted(sch):
        if mat == (1, 0, 0, 1):
            continue
        xs, ys = [], []
        for alpha in cusps:
            try:
                val = h_at(mat, alpha)
            except (ValueError, ZeroDivisionError):
                continue
            den = alpha.denominator
            img_den = abs(mat[2] * alpha.numerator + mat[3] * alpha.denominator)
            max_den = max(max_den, den, img_den)
            xs.append(alpha)
            ys.append(val)
            if len(xs) == 7:
                break
        if len(xs) < 7:
            raise SystemExit(f"only {len(xs)} cusps for {mat}")
        coeffs = interpolate(xs[:6], ys[:6])
        if eval_poly(coeffs, xs[6]) != ys[6]:
            failures += 1
            say(f"seventh cusp failed {mat}")
            continue
        cleared = cleared_integers(coeffs)
        say(f"P {mat[0]} {mat[1]} {mat[2]} {mat[3]} | " + " ".join(str(v) for v in cleared))
    if failures:
        raise SystemExit(f"{failures} interpolation failures")
    say(
        f"polynomials: 46 Schreier generators interpolated at 6 cusps, "
        f"checked at a 7th, all lie in P(7;4,3); largest cusp denominator touched {max_den}"
    )
    say("containment input: every Schreier generator has h in P(7;4,3); the identity has h = 0")
    say("K7 CERTIFICATE OK")
    out = ROOT / "certificates" / "verify-k7-log.txt"
    out.write_text("\n".join(LOG) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
