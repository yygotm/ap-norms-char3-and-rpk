# Two finite-field existence problems: AP pairs with prescribed norms in characteristic 3, and last two coefficients of r-primitive k-normal polynomials

This repository holds computer-search results and independent checks for two papers by
Chatterjee, Sharma, Shukla and Tiwari (first paper) and Chatterjee, Sharma and Tiwari
(second paper).

1. **Candidate 31.** Chatterjee–Sharma–Shukla–Tiwari, *Arithmetic progression in a finite
   field with prescribed norms*, arXiv:2401.01819, Theorem 5.7 (m = 2, q = 3^k, n >= 6)
   leaves ten pairs (q, n) as possible exceptions:
   (3,6) (9,6) (27,6) (3,7) (3,8) (9,8) (3,9) (3,10) (3,12) (3,16).
2. **Candidate 32.** Chatterjee–Sharma–Tiwari, *r-primitive k-normal polynomials over finite
   fields with last two coefficients prescribed*, arXiv:2501.04999 v3, Lemmas 4.3–4.9 leave
   78 pairs (q, m) (r = 3, k = 1, m >= 7) as possible exceptions.

## Results

### 1. arXiv:2401.01819, Theorem 5.7: none of the ten pairs is an exception

For q = 3 the norm condition is automatic (the only primitive element of F_3 is -1), so the
q = 3 cases are the question whether for every beta != 0 there is an alpha with alpha and
alpha + beta both primitive and at least one of them normal.

* (3,9) and (3,16) are not in Table 6 of Lemos–Neumann–Ribas (arXiv:2208.02876), so
  Theorem 1.1(i) there covers them; combined with the q = 3 remark above they lie in the set
  of the theorem. (This step is a combination of a published theorem and a remark, not
  something we computed.)
* The other eight pairs (3,6) (3,7) (3,8) (3,10) (3,12) (9,6) (9,8) (27,6) were settled by
  computation. For every beta (one per orbit of beta -> beta^3 and beta -> -beta) and every
  ordered pair of primitive norms (c1, c2) a certificate alpha is stored (`ap-norms/`), and a
  separate checker verifies all of them. (27,6) has 10,762,819 beta-representatives and
  1,549,845,936 certificates.

An independent re-implementation (OpenAI GPT, written without reading our code) rechecked
all eight sets of certificates: 1,571,757,379 witnesses in total, 0 failures, 0 uncovered
orbits (`reports/20261009-candidate31-independent-check-part2.md`).

### 2. arXiv:2501.04999: a problem in Lemma 3.4 / Theorem 3.7, and one genuine exception

* **Norm range.** If xi is r-primitive, xi = eps^t with gcd(t, q^m - 1) = r, so N(xi) =
  N(eps)^t has order exactly (q-1)/w, w = gcd(r, q-1). Lemma 3.4 says the norm is "w-free",
  which does not match this (for q = 4, r = 3 the only possible norm is 1, which is not
  3-free; for q = 5, r = 3 the possible norms are the elements of order 4, but every
  element is 1-free). Theorem 3.7 asserts the conclusion for every b in F_q^*, which is false
  as stated.
* **A counterexample to the printed statement of Theorem 3.7** (found by GPT, we checked it):
  q = 4, m = 64, r = 3, k = 1, g = X - 1. The sieve inequality
  q^(m/2-k-2) = 4^29 > r W(Q) W((x^m-1)/g) = 3 * 256 * 2 holds, but for b of order 3 in F_4^*
  no 3-primitive element has norm b.
* **With the corrected statement (b of order (q-1)/w), among the 78 listed pairs exactly one
  is an exception:** (q, m) = (2, 8), (a, b) = (0, 1). In F_{2^8} there are 64 3-primitive
  elements, 8 of them are 1-normal, and all eight have Tr(xi^-1) = 1; so no 3-primitive
  1-normal xi has trace of the inverse equal to 0. (A cell (a, b) here means a = Tr(xi^-1) N(xi) and b = N(xi), the unsigned last two elementary symmetric values; the minimal polynomial is X^m - a_1 X^(m-1) + ... + (-1)^(m-1) a X + (-1)^m b, so the signs depend on the parity of m.)
* All other cells of all 78 pairs are covered. For 21 pairs (q^m up to 2.7e8) by exhaustive
  enumeration, for the other 57 by explicit witnesses (462,321 of them, in
  `rpk-normal/witnesses/`).
* Independent check (GPT, own C++ code, not shown ours): the 21 pairs were fully enumerated,
  the 57 sets of witnesses were verified (irreducibility of the modulus, order, 1-normality,
  coefficients, coverage), 171 witnesses were independently rechecked using a second method. 0 mismatches
  (`reports/20261010-candidate32-independent-check.md`).

## Layout

```
ap-norms/
  claude/            apnorm.c (search), checkcert.c (checker), logs, scripts
  gpt-check/         independent C++ checker (check.cpp), Python reference, logs
  certificates-small/    certificates for (3,6) (3,7) (3,8) (3,10) (3,12) (9,6)
  certificates-large.sha256   SHA-256 of the certificates of (9,8) (91.5 MB) and (27,6) (6.24 GB)
rpk-normal/
  claude/            brute.c (exhaustive), rsearch.c (random search), logs, params
  witnesses/         the 57 witness files (.bin/.hdr), 462,321 witnesses
  gpt-check/         independent C++ checker, Python reference, logs
reports/             reports of the independent checks
```

The certificates for (9,8) (91.5 MB) and (27,6) (6.24 GB) are not stored here because of their
size. `bash ap-norms/claude/run_main.sh` (output directory `OUT`, default `~/ap31`) regenerates them (our run: (9,8) in about 30 s, (27,6) in
about 13.5 min on 16 threads with 4.5 GB peak memory). `certificates-large.sha256` lists the
SHA-256 of the files we checked; we have not verified that a rerun gives byte-identical files.

Build examples (gcc 15.2, Ubuntu on WSL):

```
gcc -O2 -fopenmp -o brute rpk-normal/claude/brute.c       # usage: brute p e m   (q = p^e)
gcc -O2 -fopenmp -o rsearch rpk-normal/claude/rsearch.c   # usage: rsearch params budget threads outprefix
gcc -O2 -fopenmp -o apnorm ap-norms/claude/apnorm.c       # see run_main.sh
```

Independent checkers (g++ with OpenMP; DATA is the directory holding the certificate folders):

```
DATA=ap-norms/certificates-small bash ap-norms/gpt-check/run.sh build
DATA=ap-norms/certificates-small bash ap-norms/gpt-check/run.sh full 3_6
bash rpk-normal/gpt-check/run.sh full all       # uses rpk-normal/witnesses
```

The `.log` files and the reports are records of the original runs. For publication, local paths were normalized to repository-relative ones and private workflow descriptions were removed; numerical results and verdicts were preserved. SHA-256 values printed inside the logs and reports are those of the original runs' files, not of the files in this repository (use `input-manifest.json` and `certificates-large.sha256` for the shipped data). The two `finish_report.py` scripts are archived records of how the original reports were written; they are not meant to be rerun from this repository alone: the AP one needs original logs that are not shipped, the RPK one checks hashes of the original (CRLF) input files and prints the date of the original run. They write to `*.regenerated` and never overwrite the archived reports. `apnorm` refuses to run in a directory that already contains `header.txt` or `cert_*.bin`; always give it an empty directory. After those runs the shared-variable accesses in `brute.c`, `apnorm.c` and `rsearch.c` were changed to atomic loads/stores (memory-order cleanup only; small cases were rerun with identical results).

## Scope and honesty about method

These are computer searches with explicit certificates and independent re-checks; they are
not formally verified proofs. The computations were carried out with the help of an AI
assistant (Anthropic's Claude) and checked independently with OpenAI's GPT. Statements
marked as a combination of a published theorem and a remark are not computed here. The
Lemma 3.4 / Theorem 3.7 discussion refers to the arXiv v3 PDF.

## Author

Hirotaka Shimizu
