// Candidate 31 (arXiv:2401.01819, Thm 5.7): for q = 3^k and F_{q^n} = F_{3^N}, N = k*n,
// decide whether for every beta != 0 and every pair (c1,c2) of primitive elements of F_q
// there is alpha with alpha, alpha+beta primitive, N(alpha)=c1, N(alpha+beta)=c2,
// and alpha or alpha+beta normal over F_q.
//
// Elements are handled as logs to a primitive element g (Zech logarithms).
//   primitive  <=> gcd(log, M) = 1,            M = 3^N - 1
//   norm index  = log mod (q-1)                 (N(g) generates F_q^*)
//   alpha+beta  = g^{j + Z[i-j]}  for alpha = g^i, beta = g^j
// beta is reduced modulo x -> 3x (absolute Frobenius) and beta -> -beta (swap the pair).
//
// usage: apnorm k n [threads] [selftest] [certdir]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <dirent.h>

static int K, n, N;
static uint64_t M, Q;        // M = 3^N - 1, Q = q
static uint32_t *expd;       // expd[t] = dense base-3 index of g^t
static uint32_t *zech;       // zech[t] = log(1 + g^t), UINT32_MAX if 1 + g^t = 0
static uint8_t *coprime;     // bitset: gcd(t, M) = 1
static uint8_t *normc;       // 2 bits per log: 0 unknown, 1 normal, 2 not normal
static uint64_t pw3[40];
static const char *certdir;

static uint64_t gcd64(uint64_t a, uint64_t b) { while (b) { uint64_t t = a % b; a = b; b = t; } return a; }
static inline int isprim(uint64_t t) { return (coprime[t >> 3] >> (t & 7)) & 1; }

// F_3 rank of N vectors (rows) of length N.
static int rank3(int8_t m[][40]) {
  int r = 0;
  for (int c = 0; c < N && r < N; c++) {
    int p = -1;
    for (int i = r; i < N; i++) if (m[i][c]) { p = i; break; }
    if (p < 0) continue;
    if (p != r) for (int j = 0; j < N; j++) { int8_t t = m[p][j]; m[p][j] = m[r][j]; m[r][j] = t; }
    int inv = m[r][c];  // 1 -> 1, 2 -> 2 (2*2 = 4 = 1 mod 3)
    for (int j = c; j < N; j++) m[r][j] = (int8_t)((m[r][j] * inv) % 3);
    for (int i = 0; i < N; i++) if (i != r && m[i][c]) {
      int f = m[i][c];
      for (int j = c; j < N; j++) m[i][j] = (int8_t)(((m[i][j] - f * m[r][j]) % 3 + 3) % 3);
    }
    r++;
  }
  return r;
}

// alpha = g^t normal over F_q <=> {g_q^u * alpha^{q^s} : u < K, s < n} is an F_3-basis.
static int normal_raw(uint64_t t) {
  int8_t m[40][40];
  uint64_t sub = M / (Q - 1);  // g_q = g^sub
  int row = 0;
  uint64_t ls = t;
  for (int s = 0; s < n; s++) {
    for (int u = 0; u < K; u++) {
      uint64_t l = (ls + (uint64_t)u * sub) % M;
      uint64_t d = expd[l];
      for (int j = 0; j < N; j++) { m[row][j] = (int8_t)(d % 3); d /= 3; }
      row++;
    }
    ls = (unsigned __int128)ls * Q % M;
  }
  return rank3(m) == N;
}

static inline int isnormal(uint64_t t) {
  uint64_t w = t >> 2; int sh = (int)(t & 3) * 2;
  int v = (__atomic_load_n(&normc[w], __ATOMIC_RELAXED) >> sh) & 3;
  if (v) return v == 1;
  int r = normal_raw(t);
  __atomic_fetch_or(&normc[w], (uint8_t)((r ? 1 : 2) << sh), __ATOMIC_RELAXED);
  return r;
}

// dense index of x * (element with dense index d) modulo f; poly[j] = coefficient of x^j in f (monic)
static uint64_t mulx(uint64_t d, const int *poly) {
  uint64_t top = d / pw3[N - 1];
  uint64_t r = (d % pw3[N - 1]) * 3;
  if (!top) return r;
  // x^N = -sum poly[j] x^j ; add top * (-poly[j]) to trit j
  uint64_t out = 0;
  for (int j = 0; j < N; j++) {
    int tr = (int)(r / pw3[j] % 3);
    int add = (int)((3 - poly[j]) % 3) * (int)top % 3;
    out += (uint64_t)((tr + add) % 3) * pw3[j];
  }
  return out;
}

static int try_poly(const int *poly) {
  uint64_t d = 1;
  for (uint64_t t = 0; t < M; t++) {
    if (t && d == 1) return 0;
    expd[t] = (uint32_t)d;
    d = mulx(d, poly);
  }
  return d == 1;
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: apnorm k n [threads] [selftest]\n"); return 2; }
  K = atoi(argv[1]); n = atoi(argv[2]); N = K * n;
  int threads = argc > 3 ? atoi(argv[3]) : 16;
  int selftest = argc > 4 ? atoi(argv[4]) : 0;
  certdir = argc > 5 ? argv[5] : NULL;
  if (certdir) {  // refuse to touch a directory that already holds certificate output (before any write)
    DIR *dp = opendir(certdir);
    if (!dp) { fprintf(stderr, "cannot open certificate directory %s\n", certdir); exit(2); }
    struct dirent *de;
    while ((de = readdir(dp))) {
      size_t l = strlen(de->d_name);
      if (!strcmp(de->d_name, "header.txt") || (l > 4 && !strncmp(de->d_name, "cert_", 5) && !strcmp(de->d_name + l - 4, ".bin"))) {
        fprintf(stderr, "%s/%s already exists: use an empty certificate directory\n", certdir, de->d_name); exit(2);
      }
    }
    closedir(dp);
  }
  omp_set_num_threads(threads);
  pw3[0] = 1; for (int i = 1; i < 40; i++) pw3[i] = pw3[i - 1] * 3;
  M = pw3[N] - 1; Q = pw3[K];
  double t0 = omp_get_wtime();
  expd = malloc(sizeof(uint32_t) * M);
  // search a primitive polynomial of degree N over F_3
  int poly[40]; int found = 0;
  for (uint64_t c = 1; c < pw3[N] && !found; c++) {
    if (c % 3 == 0) continue;  // constant term nonzero
    uint64_t x = c; for (int j = 0; j < N; j++) { poly[j] = (int)(x % 3); x /= 3; }
    found = try_poly(poly);
  }
  if (!found) { fprintf(stderr, "no primitive polynomial\n"); return 1; }
  fprintf(stderr, "q=%llu n=%d N=%d M=%llu poly:", (unsigned long long)Q, n, N, (unsigned long long)M);
  for (int j = N - 1; j >= 0; j--) fprintf(stderr, "%d", poly[j]);
  fprintf(stderr, "  [%.1fs]\n", omp_get_wtime() - t0);
  uint32_t *logd = malloc(sizeof(uint32_t) * (M + 1));
  logd[0] = UINT32_MAX;
  #pragma omp parallel for schedule(static)
  for (uint64_t t = 0; t < M; t++) logd[expd[t]] = (uint32_t)t;
  zech = malloc(sizeof(uint32_t) * M);
  #pragma omp parallel for schedule(static)
  for (uint64_t t = 0; t < M; t++) {
    uint64_t d = expd[t];
    uint64_t e = (d % 3 == 2) ? d - 2 : d + 1;
    zech[t] = logd[e];
  }
  free(logd);
  coprime = calloc(M / 8 + 1, 1);
  {
    // primes of M
    uint64_t ps[64]; int np = 0; uint64_t x = M;
    for (uint64_t p = 2; p * p <= x; p++) if (x % p == 0) { ps[np++] = p; while (x % p == 0) x /= p; }
    if (x > 1) ps[np++] = x;
    fprintf(stderr, "M primes:"); for (int i = 0; i < np; i++) fprintf(stderr, " %llu", (unsigned long long)ps[i]); fprintf(stderr, "\n");
    #pragma omp parallel for schedule(static)
    for (uint64_t b = 0; b < M / 8 + 1; b++) {
      uint8_t v = 0;
      for (int i = 0; i < 8; i++) {
        uint64_t t = b * 8 + i; if (t >= M) break;
        int ok = 1; for (int j = 0; j < np; j++) if (t % ps[j] == 0) { ok = 0; break; }
        if (ok) v |= (uint8_t)(1 << i);
      }
      coprime[b] = v;
    }
  }
  normc = calloc(M / 4 + 1, 1);
  fprintf(stderr, "tables [%.1fs]\n", omp_get_wtime() - t0);

  if (selftest) {
    uint64_t np = 0, nn = 0;
    #pragma omp parallel for reduction(+:np,nn) schedule(dynamic, 4096)
    for (uint64_t t = 0; t < M; t++) { np += isprim(t); nn += isnormal(t); }
    printf("selftest: primitive=%llu normal(nonzero)=%llu\n", (unsigned long long)np, (unsigned long long)nn);
    return 0;
  }

  // primitive norm indices c in Z/(q-1) with gcd(c, q-1) = 1
  int qm1 = (int)(Q - 1); int cidx[64]; int nc = 0; int pos[64];
  for (int c = 0; c < qm1; c++) { pos[c] = -1; if (gcd64(c, qm1) == 1) { pos[c] = nc; cidx[nc++] = c; } }
  int need = nc * nc;
  // beta orbit representatives under j -> 3j, j -> j + M/2
  uint8_t *seen = calloc(M / 8 + 1, 1);
  uint32_t *reps = malloc(sizeof(uint32_t) * (M / N + 16)); uint64_t nr = 0, cap = M / N + 16;
  for (uint64_t j = 0; j < M; j++) {
    if ((seen[j >> 3] >> (j & 7)) & 1) continue;
    if (nr == cap) { cap *= 2; reps = realloc(reps, sizeof(uint32_t) * cap); }
    reps[nr++] = (uint32_t)j;
    uint64_t a = j;
    for (int s = 0; s < N; s++) {
      uint64_t b = (a + M / 2) % M;
      seen[a >> 3] |= (uint8_t)(1 << (a & 7)); seen[b >> 3] |= (uint8_t)(1 << (b & 7));
      a = a * 3 % M;
    }
  }
  free(seen);
  fprintf(stderr, "beta reps=%llu pairs per beta=%d [%.1fs]\n", (unsigned long long)nr, need, omp_get_wtime() - t0);

  // stride coprime to M for the alpha scan order
  uint64_t stride = 1000003; while (gcd64(stride, M) != 1) stride += 2;
  uint64_t fails = 0, maxscan = 0, totscan = 0;
  FILE *ff = fopen("fails.txt", "w");
  // certificate: header + one file per thread of records {uint32 j, uint32 alpha_log[need]},
  // alpha_log[pos[c1]*nc + pos[c2]] = log of alpha with N(alpha) = g_q^c1, N(alpha+beta) = g_q^c2
  FILE *cf[256] = {0};
  if (certdir) {
    char path[1024];
    snprintf(path, sizeof path, "%s/header.txt", certdir);
    FILE *h = fopen(path, "w");
    fprintf(h, "q=%llu n=%d N=%d M=%llu\npoly_high_to_low=", (unsigned long long)Q, n, N, (unsigned long long)M);
    for (int j = N - 1; j >= 0; j--) fprintf(h, "%d", poly[j]);
    fprintf(h, "\nmonic x^N + sum poly[j] x^j over F_3; g = x; g_q = g^(M/(q-1)); dense element = sum trit_j 3^j\n");
    fprintf(h, "beta = g^j, j over representatives of j -> 3j, j -> j+M/2 (count %llu)\n", (unsigned long long)nr);
    fprintf(h, "primitive norm indices c (gcd(c,q-1)=1):");
    for (int c = 0; c < nc; c++) fprintf(h, " %d", cidx[c]);
    fprintf(h, "\nrecord = uint32 j, uint32 alpha_log[%d] (index pos(c1)*%d+pos(c2)), little endian\n", need, nc);
    fclose(h);
    for (int t = 0; t < threads; t++) {
      snprintf(path, sizeof path, "%s/cert_%03d.bin", certdir, t);
      cf[t] = fopen(path, "wb");
    }
  }
  #pragma omp parallel for schedule(dynamic, 256) reduction(+:fails,totscan) reduction(max:maxscan)
  for (uint64_t r = 0; r < nr; r++) {
    uint64_t j = reps[r];
    uint32_t rec[1 + 144];
    rec[0] = (uint32_t)j;
    uint64_t got[3] = {0, 0, 0}; int have = 0;  // need <= 144 bits
    uint64_t scan = 0, i = (j * 7919 + 1) % M;
    for (; scan < M && have < need; scan++, i = (i + stride) % M) {
      if (!isprim(i)) continue;
      uint64_t dlt = (i + M - j) % M;
      uint32_t z = zech[dlt];
      if (z == UINT32_MAX) continue;
      uint64_t s = (j + z) % M;
      if (!isprim(s)) continue;
      int b = pos[i % qm1] * nc + pos[s % qm1];
      if ((got[b >> 6] >> (b & 63)) & 1) continue;
      if (!isnormal(i) && !isnormal(s)) continue;
      got[b >> 6] |= 1ULL << (b & 63); have++;
      rec[1 + b] = (uint32_t)i;
    }
    totscan += scan; if (scan > maxscan) maxscan = scan;
    if (have == need && certdir) fwrite(rec, sizeof(uint32_t), 1 + need, cf[omp_get_thread_num()]);
    if (have < need) {
      fails++;
      #pragma omp critical
      { fprintf(ff, "beta=g^%llu covered %d of %d\n", (unsigned long long)j, have, need); fflush(ff); }
    }
  }
  fclose(ff);
  for (int t = 0; t < threads; t++) if (cf[t]) fclose(cf[t]);
  printf("q=%llu n=%d: beta reps=%llu fails=%llu mean_scan=%.1f max_scan=%llu time=%.1fs\n",
         (unsigned long long)Q, n, (unsigned long long)nr, (unsigned long long)fails,
         (double)totscan / nr, (unsigned long long)maxscan, omp_get_wtime() - t0);
  return 0;
}
