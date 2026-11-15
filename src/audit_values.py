"""Independent finite-sum checks used in the Hydrogen audit.

Recomputes:
  * S-tilde_{chi_3, chi_5, 7}(1, -1; 15, -14) and two further matrices;
  * the same value from the author's stored degree-5 polynomial via formula (12);
  * the Hurwitz (k-1) denominator statement on the finite range Hydrogen stated,
    plus a few conductors and weights outside that range.

Does not call any module under dedekind_research.
"""

from fractions import Fraction
from math import comb, gcd
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
LOG = []


def say(msg):
    print(msg)
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


_CACHE = {}


def B(k, x):
    x = Fraction(x)
    if x.denominator == 1:
        return Fraction(0)
    if k not in _CACHE:
        _CACHE[k] = bernoulli_coeffs(k)
    t = x - (x.numerator // x.denominator)
    val = Fraction(0)
    for c in reversed(_CACHE[k]):
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


def legendre(n, p):
    a = n % p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return -1 if r == p - 1 else 1


def chi4(n):
    r = n % 4
    if r % 2 == 0:
        return 0
    return 1 if r == 1 else -1


def chi8a(n):
    r = n % 8
    if r % 2 == 0:
        return 0
    return 1 if r in (1, 7) else -1


def stilde(k, chi2, q2, a, c):
    """chi_1 = chi_3, conductor 3. Sign for c<0: S-tilde(gamma) = (-1)^{k-2} S-tilde(-gamma)."""
    if c < 0:
        sign = 1 if (k - 2) % 2 == 0 else -1
        return sign * stilde(k, chi2, q2, -a, -c)
    if c == 0:
        return Fraction(0)
    total = Fraction(0)
    for j in range(c):
        cj = chi2(j)
        if cj == 0:
            continue
        b1 = B(1, Fraction(j, c))
        for n in (1, 2):
            total += chi3(n) * cj * b1 * B(k - 1, Fraction(a * j, c) + Fraction(n, 3))
    return (c ** (k - 2)) * total


def rgcd(x, y):
    x, y = Fraction(x), Fraction(y)
    return abs(
        Fraction(
            gcd(x.numerator * y.denominator, y.numerator * x.denominator),
            x.denominator * y.denominator,
        )
    )


def check_k7():
    m1 = (1, -1, 15, -14)
    if m1[0] * m1[3] - m1[1] * m1[2] != 1 or m1[0] % 15 != 1 or m1[3] % 15 != 1 or m1[2] % 15:
        raise SystemExit("m1 not in Gamma_1(15)")
    v1 = stilde(7, chi5, 5, 1, 15)
    if v1 != Fraction(-17000, 3):
        raise SystemExit(v1)
    say(f"k7 direct sum: S-tilde(1,-1;15,-14) = {v1}")
    say(f"k7 in 4Z? {v1.denominator == 1 and int(v1) % 4 == 0}")

    # Author polynomial of gamma = inverse of m1, evaluated by formula (12).
    # h = 17000/3 x^5 - 12800/9 x^4 + 400/3 x^3 - 160/27 x^2
    # gamma = (-14, 1, -15, 1), gamma^{-1} = m1, cusp of the inverse is 1/15.
    x = Fraction(1, 15)
    h = (
        Fraction(17000, 3) * x ** 5
        - Fraction(12800, 9) * x ** 4
        + Fraction(400, 3) * x ** 3
        - Fraction(160, 27) * x ** 2
    )
    pred = (15 ** 5) * h
    if pred != v1:
        raise SystemExit(f"author polynomial gives {pred}, sum gives {v1}")
    say(f"k7 formula (12) on the author's first polynomial also gives {pred}")

    # Positive-c second matrix, so the sign convention is not used.
    m3 = (16, -9, 105, -59)
    if m3[0] * m3[3] - m3[1] * m3[2] != 1 or m3[0] % 15 != 1 or m3[3] % 15 != 1 or m3[2] % 15:
        raise SystemExit("m3 illegal")
    v3 = stilde(7, chi5, 5, m3[0], m3[2])
    g13 = rgcd(v1, v3)
    say(f"k7 second positive-c value S-tilde(16,-9;105,-59) = {v3}")
    say(f"k7 gcd of the two positive-c values = {g13}")

    m2 = (-29, 4, -225, 31)
    if m2[0] * m2[3] - m2[1] * m2[2] != 1 or m2[0] % 15 != 1 or m2[3] % 15 != 1 or m2[2] % 15:
        raise SystemExit("m2 illegal")
    v2 = stilde(7, chi5, 5, m2[0], m2[2])
    g = rgcd(g13, v2)
    say(f"k7 third value, negative c, sign (-1)^{{k-2}}: {v2}")
    say(f"k7 gcd of the three values = {g}")
    say("k7 span statement: these values generate a group containing (4/3)Z; they do not prove every value lies in (4/3)Z")


def hurwitz_bad(q1, chi, k):
    b = bernoulli_coeffs(k - 1)
    r = [Fraction(0)] * k
    for n in range(1, q1):
        cn = chi(n)
        if cn == 0:
            continue
        shift = Fraction(n, q1)
        for m in range(k):
            for j in range(m + 1):
                r[j] += cn * b[m] * comb(m, j) * (shift ** (m - j))
    r = [c / (k - 1) for c in r]
    rad = []
    t = q1
    p = 2
    while p * p <= t:
        if t % p == 0:
            rad.append(p)
            while t % p == 0:
                t //= p
        p += 1
    if t > 1:
        rad.append(t)
    bad = []
    for c in r:
        d = c.denominator
        for rp in rad:
            while d % rp == 0:
                d //= rp
        if d != 1:
            bad.append(c)
    return bad


def check_hurwitz():
    chars = [
        (3, chi3),
        (4, chi4),
        (5, chi5),
        (7, lambda n: legendre(n, 7)),
        (8, chi8a),
    ]
    fails = []
    for q1, chi in chars:
        for k in range(3, 10):
            if hurwitz_bad(q1, chi, k):
                fails.append((q1, k))
    say(f"hurwitz stated range q in {{3,4,5,7,8}}, k=3..9: failures {fails or 'none'}")
    probes = []
    for q1, ks in [(3, (10, 12, 14)), (11, (3, 4, 5, 6, 8, 10)), (13, (3, 4, 6))]:
        chi = chi3 if q1 == 3 else (lambda n, q1=q1: legendre(n, q1))
        for k in ks:
            bad = hurwitz_bad(q1, chi, k)
            probes.append((q1, k, "fail" if bad else "ok"))
    say("hurwitz probes outside the stated range: " + ", ".join(f"({q},{k}) {s}" for q, k, s in probes))
    say("hurwitz status: finite experiment, not a theorem for all conductors and weights")


def main():
    say("python " + sys.version.split()[0])
    check_k7()
    check_hurwitz()
    say("AUDIT VALUES OK")
    (ROOT / "certificates" / "audit-values-log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
