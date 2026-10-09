// Checker for apnorm certificates (candidate 31). Written separately from apnorm.c:
// no Zech table; alpha+beta is computed by digit-wise addition of base-3 vectors.
//
// For every record {j, alpha_log[need]} it checks, with beta = g^j, alpha = g^i:
//   gcd(i, M) = 1, alpha+beta != 0 and gcd(log(alpha+beta), M) = 1,
//   i mod (q-1) = c1 and log(alpha+beta) mod (q-1) = c2 for the pair the slot stands for,
//   alpha or alpha+beta is normal over F_q (F_3-rank of {g_q^u alpha^{q^s}} is N).
// It also checks that x is primitive for the polynomial (period exactly M) and that
// the set of j covers every class of j -> 3j, j -> j+M/2.
//
// usage: checkcert certdir threads
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

static int K, n, N; static uint64_t M, Q, P3[40];
static uint32_t *ex, *lg;

static int rank_f3(uint8_t a[40][40], int rows, int cols) {
  int r = 0;
  for (int c = 0; c < cols && r < rows; c++) {
    int p = r; while (p < rows && a[p][c] == 0) p++;
    if (p == rows) continue;
    for (int k = 0; k < cols; k++) { uint8_t t = a[r][k]; a[r][k] = a[p][k]; a[p][k] = t; }
    if (a[r][c] == 2) for (int k = 0; k < cols; k++) a[r][k] = (uint8_t)(a[r][k] * 2 % 3);
    for (int i = 0; i < rows; i++) {
      if (i == r || a[i][c] == 0) continue;
      uint8_t f = a[i][c];
      for (int k = 0; k < cols; k++) a[i][k] = (uint8_t)((a[i][k] + 3 * 3 - f * a[r][k]) % 3);
    }
    r++;
  }
  return r;
}

static int normal(uint64_t i) {
  uint8_t a[40][40]; int row = 0;
  uint64_t s0 = i, sub = M / (Q - 1);
  for (int s = 0; s < n; s++) {
    for (int u = 0; u < K; u++) {
      uint32_t d = ex[(s0 + (uint64_t)u * sub) % M];
      for (int t = 0; t < N; t++) { a[row][t] = (uint8_t)(d % 3); d /= 3; }
      row++;
    }
    s0 = (uint64_t)((unsigned __int128)s0 * Q % M);
  }
  return rank_f3(a, N, N) == N;
}

static uint32_t add3(uint32_t x, uint32_t y) {
  uint32_t r = 0;
  for (int t = 0; t < N; t++) { r += (uint32_t)(((x % 3) + (y % 3)) % 3) * (uint32_t)P3[t]; x /= 3; y /= 3; }
  return r;
}

static uint64_t gcdu(uint64_t a, uint64_t b) { while (b) { uint64_t t = a % b; a = b; b = t; } return a; }

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: checkcert certdir threads\n"); return 2; }
  const char *dir = argv[1]; int threads = atoi(argv[2]); omp_set_num_threads(threads);
  char path[1024], polys[64]; unsigned long long qq, mm;
  snprintf(path, sizeof path, "%s/header.txt", dir);
  FILE *h = fopen(path, "r"); if (!h) { perror(path); return 1; }
  if (fscanf(h, "q=%llu n=%d N=%d M=%llu\npoly_high_to_low=%63s", &qq, &n, &N, &mm, polys) != 5) { fprintf(stderr, "bad header\n"); return 1; }
  fclose(h);
  Q = qq; M = mm; K = 0; for (uint64_t t = Q; t > 1; t /= 3) K++;
  P3[0] = 1; for (int t = 1; t < 40; t++) P3[t] = P3[t - 1] * 3;
  if (K * n != N || P3[N] - 1 != M || (int)strlen(polys) != N) { fprintf(stderr, "header inconsistent\n"); return 1; }
  int f[40]; for (int t = 0; t < N; t++) f[t] = polys[N - 1 - t] - '0';
  // powers of x as trit arrays, multiplying by x and reducing with x^N = -sum f[t] x^t
  ex = malloc(sizeof(uint32_t) * M); lg = malloc(sizeof(uint32_t) * (M + 1));
  for (uint64_t t = 0; t <= M; t++) lg[t] = UINT32_MAX;
  int v[40] = {0}; v[0] = 1;
  for (uint64_t e = 0; e < M; e++) {
    uint32_t d = 0; for (int t = 0; t < N; t++) d += (uint32_t)v[t] * (uint32_t)P3[t];
    if (lg[d] != UINT32_MAX) { fprintf(stderr, "x is not primitive (repeat at %llu)\n", (unsigned long long)e); return 1; }
    ex[e] = d; lg[d] = (uint32_t)e;
    int top = v[N - 1];
    for (int t = N - 1; t > 0; t--) v[t] = v[t - 1];
    v[0] = 0;
    for (int t = 0; t < N; t++) v[t] = (v[t] + 3 * 3 - top * f[t]) % 3;
  }
  if (!(v[0] == 1 && ({ int z = 1; for (int t = 1; t < N; t++) if (v[t]) z = 0; z; }))) { fprintf(stderr, "x^M != 1\n"); return 1; }
  int qm1 = (int)(Q - 1), cs[64], nc = 0;
  for (int c = 0; c < qm1; c++) if (gcdu(c, qm1) == 1) cs[nc++] = c;
  int need = nc * nc;
  // class bitset for coverage
  uint8_t *cov = calloc(M / 8 + 1, 1);
  uint64_t recs = 0, bad = 0;
  for (int t = 0; t < 256; t++) {
    snprintf(path, sizeof path, "%s/cert_%03d.bin", dir, t);
    FILE *fp = fopen(path, "rb"); if (!fp) continue;
    size_t B = 1 << 16; uint32_t *buf = malloc(sizeof(uint32_t) * (1 + need) * B);
    size_t got;
    while ((got = fread(buf, sizeof(uint32_t) * (1 + need), B, fp)) > 0) {
      #pragma omp parallel for schedule(dynamic, 64) reduction(+:bad)
      for (size_t r = 0; r < got; r++) {
        uint32_t *rec = buf + r * (1 + need);
        uint64_t j = rec[0]; int ok = j < M;
        for (int a = 0; ok && a < nc; a++) for (int b = 0; ok && b < nc; b++) {
          uint64_t i = rec[1 + a * nc + b];
          if (i >= M || gcdu(i, M) != 1 || (int)(i % qm1) != cs[a]) { ok = 0; break; }
          uint32_t sd = add3(ex[i], ex[j]);
          if (sd == 0) { ok = 0; break; }
          uint64_t s = lg[sd];
          if (gcdu(s, M) != 1 || (int)(s % qm1) != cs[b]) { ok = 0; break; }
          if (!normal(i) && !normal(s)) { ok = 0; break; }
        }
        if (!ok) { bad++;
          #pragma omp critical
          fprintf(stderr, "bad record j=%u\n", rec[0]); }
      }
      for (size_t r = 0; r < got; r++) {
        uint64_t a = buf[r * (1 + need)];
        for (int s = 0; s < N; s++) {
          uint64_t b = (a + M / 2) % M;
          cov[a >> 3] |= (uint8_t)(1 << (a & 7)); cov[b >> 3] |= (uint8_t)(1 << (b & 7));
          a = (uint64_t)((unsigned __int128)a * 3 % M);
        }
      }
      recs += got;
    }
    free(buf); fclose(fp);
  }
  uint64_t miss = 0; for (uint64_t j = 0; j < M; j++) if (!((cov[j >> 3] >> (j & 7)) & 1)) miss++;
  printf("q=%llu n=%d: records=%llu bad=%llu uncovered_beta=%llu -> %s\n", qq, n,
         (unsigned long long)recs, (unsigned long long)bad, (unsigned long long)miss,
         (bad == 0 && miss == 0) ? "CERTIFIED" : "FAILED");
  return (bad == 0 && miss == 0) ? 0 : 1;
}
