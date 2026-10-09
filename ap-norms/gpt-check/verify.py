"""Independent, deliberately literal checker for the AP/norm certificates.

Only the exchange request and the binary format described there were used.
The --pairs limit is for a timing sample; omit it for an exhaustive run.
"""
import argparse
import json
import math
import struct
import time
from pathlib import Path


def prime_divisors(m):
    ans = []
    d = 2
    while d * d <= m:
        if m % d == 0:
            ans.append(d)
            while m % d == 0:
                m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        ans.append(m)
    return ans


class Field:
    def __init__(self, digits):
        self.n = len(digits)
        self.f = tuple(map(int, reversed(digits)))
        self.zero = (0,) * self.n
        self.one = (1,) + (0,) * (self.n - 1)
        self.x = (0, 1) + (0,) * (self.n - 2)

    def add(self, a, b):
        return tuple((u + v) % 3 for u, v in zip(a, b))

    def mul(self, a, b):
        n = self.n
        c = [0] * (2 * n - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    c[i + j] += ai * bj
        f = self.f
        for t in range(2 * n - 2, n - 1, -1):
            z = c[t] % 3
            if z:
                for j in range(n):
                    c[t - n + j] -= z * f[j]
        return tuple(z % 3 for z in c[:n])

    def power(self, a, e):
        out = self.one
        while e:
            if e & 1:
                out = self.mul(out, a)
            e >>= 1
            if e:
                a = self.mul(a, a)
        return out

    def normal(self, a, q, k, n, basis):
        # Rank over F_3 of {b_l * a^(q^i): 0<=l<k, 0<=i<n}.
        # F_q-basis b_l is {1,g_q,...,g_q^(k-1)}.
        columns = []
        z = a
        for _ in range(n):
            for b in basis:
                columns.append(self.mul(b, z))
            z = self.power(z, q)
        rows = [[columns[j][i] for j in range(self.n)] for i in range(self.n)]
        rank = 0
        for col in range(self.n):
            pivot = next((r for r in range(rank, self.n) if rows[r][col]), None)
            if pivot is None:
                continue
            rows[rank], rows[pivot] = rows[pivot], rows[rank]
            inv = rows[rank][col]  # 1 and 2 are self-inverse in F_3
            if inv == 2:
                rows[rank] = [(2 * v) % 3 for v in rows[rank]]
            for r in range(rank + 1, self.n):
                z = rows[r][col]
                if z:
                    rows[r] = [(u - z * v) % 3 for u, v in zip(rows[r], rows[rank])]
            rank += 1
            if rank == self.n:
                return True
        return False


def read_header(folder):
    lines = (folder / 'header.txt').read_text(encoding='utf-8').splitlines()
    first = dict(piece.split('=') for piece in lines[0].split())
    digits = next(s.split('=', 1)[1] for s in lines if s.startswith('poly_high_to_low='))
    return *(int(first[key]) for key in ('q', 'n', 'N', 'M')), digits


def orbit_coverage(files, row_bytes, M, N):
    """Mark every exponent reached by Frobenius and sign from each stored j."""
    bits = bytearray((M + 7) // 8)
    marked = 0
    repeated_representatives = 0
    half = M // 2
    for path in files:
        with path.open('rb') as f:
            while True:
                row = f.read(row_bytes)
                if not row:
                    break
                if len(row) != row_bytes:
                    raise ValueError(f'incomplete record in {path}')
                j = struct.unpack_from('<I', row)[0]
                if j >= M:
                    raise ValueError(f'j out of range in {path}: {j}')
                if bits[j >> 3] & (1 << (j & 7)):
                    repeated_representatives += 1
                y = j
                for _ in range(N):
                    for z in (y, (y + half) % M):
                        index = z >> 3
                        mask = 1 << (z & 7)
                        if not bits[index] & mask:
                            bits[index] |= mask
                            marked += 1
                    y = 3 * y % M
    return {'uncovered_orbit_elements': M - marked,
            'repeated_representatives': repeated_representatives}


def check(folder, pair_limit):
    q, n, N, M, digits = read_header(folder)
    k = N // n
    assert q == 3 ** k and M == 3 ** N - 1 and len(digits) == N
    need_indices = [c for c in range(1, q) if math.gcd(c, q - 1) == 1]
    need = len(need_indices) ** 2
    row_bytes = 4 * (need + 1)
    field = Field(digits)
    factors = prime_divisors(M)
    assert field.power(field.x, M) == field.one
    assert all(field.power(field.x, M // p) != field.one for p in factors)
    gq = field.power(field.x, M // (q - 1))
    target = [field.power(gq, c) for c in need_indices]
    basis = [field.power(gq, l) for l in range(k)]
    files = sorted(folder.glob('cert_*.bin'))
    total_bytes = sum(p.stat().st_size for p in files)
    assert all(p.stat().st_size % row_bytes == 0 for p in files)
    expected_records = total_bytes // row_bytes
    start = time.perf_counter()
    failures = []
    failure_count = 0
    checked = 0
    visited_records = 0
    for path in files:
        with path.open('rb') as f:
            while True:
                chunk = f.read(row_bytes)
                if not chunk:
                    break
                assert len(chunk) == row_bytes
                data = struct.unpack('<' + 'I' * (need + 1), chunk)
                j = data[0]
                assert j < M
                beta = field.power(field.x, j)
                visited_records += 1
                for slot, exponent in enumerate(data[1:]):
                    if pair_limit is not None and checked >= pair_limit:
                        break
                    expected_norm1 = target[slot // len(need_indices)]
                    expected_norm2 = target[slot % len(need_indices)]
                    alpha = field.power(field.x, exponent) if exponent < M else field.zero
                    other = field.add(alpha, beta)
                    ok = exponent < M and math.gcd(exponent, M) == 1
                    if ok:
                        ok = all(field.power(alpha, M // p) != field.one for p in factors)
                    if ok:
                        ok = other != field.zero and all(field.power(other, M // p) != field.one for p in factors)
                    if ok:
                        ok = field.power(alpha, M // (q - 1)) == expected_norm1
                    if ok:
                        ok = field.power(other, M // (q - 1)) == expected_norm2
                    if ok:
                        ok = field.normal(alpha, q, k, n, basis) or field.normal(other, q, k, n, basis)
                    if not ok:
                        failure_count += 1
                        if len(failures) < 10:
                            failures.append({'file': path.name, 'record': visited_records, 'j': j, 'slot': slot, 'i': exponent})
                    checked += 1
                if pair_limit is not None and checked >= pair_limit:
                    break
        if pair_limit is not None and checked >= pair_limit:
            break
    elapsed = time.perf_counter() - start
    coverage = orbit_coverage(files, row_bytes, M, N) if pair_limit is None else None
    return {'pair': [q, n], 'M': M, 'prime_divisors': factors,
            'record_bytes': row_bytes, 'records_total': expected_records,
            'pairs_total': expected_records * need, 'records_visited': visited_records,
            'pairs_checked': checked, 'failure_count': failure_count,
            'failures_first_10': failures, 'orbit_coverage': coverage,
            'seconds': elapsed, 'estimated_full_hours_at_sample_rate':
            elapsed / checked * expected_records * need / 3600 if checked else None}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=Path)
    parser.add_argument('--pairs', type=int, default=None)
    args = parser.parse_args()
    print(json.dumps(check(args.folder, args.pairs), ensure_ascii=False, indent=2))
