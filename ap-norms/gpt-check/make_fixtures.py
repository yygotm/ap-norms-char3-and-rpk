"""Generate cross-language arithmetic/rank fixtures from our own Python code."""
import random
from pathlib import Path
from verify import Field, read_header

import os
root = Path(os.environ.get('DATA', '.'))  # directory holding the <q>_<n> certificate folders
out = Path(__file__).with_name('fixtures.txt')
rng = random.Random(202610092056)

def unpack(a, n):
    values = []
    for _ in range(n):
        values.append(a % 3)
        a //= 3
    return tuple(values)

def pack(a):
    return sum(v * 3 ** i for i, v in enumerate(a))

with out.open('w', encoding='ascii', newline='\n') as stream:
    for pair in ('3_6', '3_7', '3_8', '3_10', '3_12', '9_6', '9_8', '27_6'):
        q, n, N, M, poly = read_header(root / pair)
        field = Field(poly)
        k = N // n
        gq = field.power(field.x, M // (q - 1))
        basis = [field.power(gq, ell) for ell in range(k)]
        for a in [0, 1, 2, 3, M] + [rng.randrange(M + 1) for _ in range(24)]:
            b, e = rng.randrange(M + 1), rng.randrange(M + 1)
            aa, bb = unpack(a, N), unpack(b, N)
            values = (N, poly, q, n, a, b, e,
                      pack(field.add(aa, bb)), pack(field.mul(aa, bb)),
                      pack(field.power(aa, e)), int(field.normal(aa, q, k, n, basis)))
            stream.write(' '.join(map(str, values)) + '\n')
print(f'Wrote {8 * 29} fixtures to {out}')
