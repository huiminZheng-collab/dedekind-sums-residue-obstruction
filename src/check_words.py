"""Multiply stored Schreier words back to their matrices.

The search that found the words is not repeated. Completeness is the
equality of two sets: the Schreier generators computed from S and T, and
the products of the stored words.

Usage:
    python src/check_words.py certificates/gamma1_9_words.txt
"""

from collections import deque
import sys


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


def schreier_generators(N):
    """Reidemeister-Schreier generators of Gamma_1(N) from {S, T, S^{-1}, T^{-1}}."""
    S = (0, -1, 1, 0)
    T = (1, 1, 0, 1)
    X = [S, T, matinv(S), matinv(T)]
    reps = {(0, 1): (1, 0, 0, 1)}
    qq = deque([(0, 1)])
    while qq:
        key = qq.popleft()
        r = reps[key]
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
            if det(h) != 1 or h[0] % N != 1 or h[3] % N != 1 or h[2] % N != 0:
                raise SystemExit(f"Schreier generator outside Gamma_1({N}): {h}")
            sch.add(h)
    return len(reps), sch


def load_words(path):
    N = None
    alph = []
    words = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0] == "N":
                N = int(parts[1])
            elif parts[0] == "A":
                idx = int(parts[1])
                mat = tuple(int(x) for x in parts[2:6])
                if idx != len(alph):
                    raise SystemExit(f"alphabet index {idx} out of order")
                alph.append(mat)
            elif parts[0] == "W":
                bar = parts.index("|")
                mat = tuple(int(x) for x in parts[1:bar])
                word = [int(x) for x in parts[bar + 1 :]]
                words.append((mat, word))
    if N is None or not alph or not words:
        raise SystemExit(f"incomplete word file {path}")
    return N, alph, words


def check_file(path, expected_cosets=None, expected_schreier=None):
    N, alph, words = load_words(path)
    cosets, sch = schreier_generators(N)
    if expected_cosets is not None and cosets != expected_cosets:
        raise SystemExit(f"cosets {cosets} != {expected_cosets}")
    if expected_schreier is not None and len(sch) != expected_schreier:
        raise SystemExit(f"Schreier count {len(sch)} != {expected_schreier}")
    produced = {(1, 0, 0, 1)}
    max_len = 0
    for mat, word in words:
        acc = (1, 0, 0, 1)
        for i in word:
            acc = matmul(acc, alph[i])
        if acc != mat:
            raise SystemExit(f"word does not multiply to {mat}, got {acc}")
        if acc in produced and acc != (1, 0, 0, 1):
            raise SystemExit(f"duplicate word target {acc}")
        produced.add(acc)
        max_len = max(max_len, len(word))
    if produced != sch:
        missing = sch - produced
        extra = produced - sch
        raise SystemExit(f"word set mismatch missing {len(missing)} extra {len(extra)}")
    nontrivial = len(sch) - (1 if (1, 0, 0, 1) in sch else 0)
    return {
        "N": N,
        "cosets": cosets,
        "schreier": len(sch),
        "nontrivial_words": nontrivial,
        "alphabet": len(alph),
        "max_word_length": max_len,
    }


def main(argv):
    if len(argv) != 2:
        raise SystemExit("usage: check_words.py WORDFILE")
    info = check_file(argv[1])
    print(
        f"words ok: Gamma_1({info['N']}), {info['cosets']} cosets, "
        f"{info['schreier']} Schreier generators, {info['nontrivial_words']} words, "
        f"alphabet {info['alphabet']}, max length {info['max_word_length']}"
    )


if __name__ == "__main__":
    main(sys.argv)
