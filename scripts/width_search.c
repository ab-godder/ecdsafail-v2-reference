/* Exact-integer port of scripts/width_check.py's `search`; see that file and DESIGN.md,
 * "A width cut at secp256k1". Same nodes, same order, same tests, in fixed-width integers.
 *
 * usage: width_search Q_HEX LOG_UMIN SHARD SHARDS
 * Prints "T a1 a0 b1 b0 k" for each terminal solution with odd u (the caller derives m and x
 * and applies 2m < q), then "N nodes" for this shard's node count.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
#define L 8 /* 512-bit words */
typedef struct { u64 w[L]; } W;

enum { P = 32, SHARD_DEPTH = 12 };

static W Q;            /* q */
static W QP;           /* q << P */
static int QP_BITS;
static int LOG_UMIN, SHARD, SHARDS;
static long long nodes, seen;

static void die(const char *m) { fprintf(stderr, "width_search: %s\n", m); exit(2); }

static W w128(u128 x) { W r; memset(&r, 0, sizeof r); r.w[0] = (u64)x; r.w[1] = (u64)(x >> 64); return r; }

static int bitlen(const W *a) {
  for (int i = L - 1; i >= 0; i--)
    if (a->w[i]) return 64 * i + 64 - __builtin_clzll(a->w[i]);
  return 0;
}

static int cmp(const W *a, const W *b) {
  for (int i = L - 1; i >= 0; i--)
    if (a->w[i] != b->w[i]) return a->w[i] < b->w[i] ? -1 : 1;
  return 0;
}

static W add(const W *a, const W *b) {
  W r; u128 c = 0;
  for (int i = 0; i < L; i++) { c += (u128)a->w[i] + b->w[i]; r.w[i] = (u64)c; c >>= 64; }
  if (c) die("add overflow");
  return r;
}

static W shl(const W *a, int s) {
  if (s < 0) die("negative shift");
  if (bitlen(a) + s > 64 * L - 1) die("shift overflow");
  W r; memset(&r, 0, sizeof r);
  int q = s / 64, b = s % 64;
  for (int i = L - 1; i >= q; i--) {
    u64 v = a->w[i - q] << b;
    if (b && i - q - 1 >= 0) v |= a->w[i - q - 1] >> (64 - b);
    r.w[i] = v;
  }
  return r;
}

static W mul(const W *a, const W *b) {
  if (bitlen(a) + bitlen(b) > 64 * L - 1) die("mul overflow");
  W r; memset(&r, 0, sizeof r);
  for (int i = 0; i < L; i++) {
    if (!a->w[i]) continue;
    u128 c = 0;
    for (int j = 0; i + j < L; j++) {
      c += (u128)a->w[i] * b->w[j] + r.w[i + j];
      r.w[i + j] = (u64)c; c >>= 64;
    }
  }
  return r;
}

static W mul_small(const W *a, u64 k) { W b = w128(k); return mul(a, &b); }

/* a << s compared with b, for s that may push a past 512 bits: then a << s is larger,
 * as every b here is far below 2^511. */
static int cmp_shl(const W *a, int s, const W *b) {
  if (bitlen(b) > 500) die("comparand too large");
  if (bitlen(a) && bitlen(a) + s > 511) return 1;
  W t = shl(a, s);
  return cmp(&t, b);
}

static u64 mod_small(const W *a, u64 m) {
  u128 r = 0;
  for (int i = L - 1; i >= 0; i--) r = ((r << 64) | a->w[i]) % m;
  return (u64)r;
}

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }

static void terminal(u64 a1, u64 a0, u64 b1, u64 b0, const W *IN0, const W *IN1) {
  /* a1 u + a0 4^k = q with IN0 4^k <= q 2^P <= IN1 4^k; odd u. */
  int k_lo = (QP_BITS - bitlen(IN1) - 1) / 2 - 1, k_hi = (QP_BITS - bitlen(IN0) + 1) / 2 + 1;
  if (QP_BITS - bitlen(IN1) - 1 < 0) k_lo = -1;   /* Python floor division */
  if (k_lo < 1) k_lo = 1;
  for (int k = k_lo; k <= k_hi; k++) {
    if (cmp_shl(IN0, 2 * k, &QP) > 0 || cmp_shl(IN1, 2 * k, &QP) < 0) continue;
    /* r = q - a0 4^k > 0, r % a1 == 0, r / a1 odd  <=>  r mod 2 a1 == a1 (r > 0 follows). */
    if (a1 >= (1ULL << 62)) die("a1 too large");
    u64 m2 = 2 * a1;
    u64 p4 = 1 % m2, b4 = 4 % m2;
    for (int e = k; e; e >>= 1) { if (e & 1) p4 = mulmod(p4, b4, m2); b4 = mulmod(b4, b4, m2); }
    u64 qm = mod_small(&Q, m2), t = mulmod(a0 % m2, p4, m2);
    u64 r = (qm + m2 - t) % m2;
    if (r == a1) printf("T %llu %llu %llu %llu %d\n", (unsigned long long)a1, (unsigned long long)a0,
                        (unsigned long long)b1, (unsigned long long)b0, k);
  }
}

static u64 small(u128 v) { if (v >> 62) die("coefficient too large"); return (u64)v; }

static void node(u64 a1, u64 a0, u64 b1, u64 b0, int T, u64 Y0, u64 Y1, int depth) {
  if (depth == SHARD_DEPTH) {
    seen++;
    if (seen % SHARDS != SHARD) return;
  }
  if (depth >= SHARD_DEPTH || SHARD == 0) nodes++;
  u128 IA_ = (u128)a1 * Y0 + ((u128)a0 << P), IB_ = (u128)b1 * Y0 + ((u128)b0 << P);
  if ((IA_ >> 120) || (IB_ >> 120)) die("form too large");
  W IA = w128(IA_), IB = w128(IB_);
  W hi = IA_ >= IB_ ? IA : IB, lo = IA_ >= IB_ ? IB : IA;
  W Y1w = w128(Y1), Y1sq = mul(&Y1w, &Y1w);
  W hh = mul(&hi, &hi), hhl = mul(&hh, &lo);
  if (cmp_shl(&Y1sq, 2 * T + 1 + P, &hhl) <= 0) return;           /* hi^2 lo >= Y1^2 << 2T+1+P */
  W cap = mul(&Q, &Y1w);
  if (cmp_shl(&hi, LOG_UMIN + 1, &cap) > 0) return;
  if (Y1 - Y0 > (Y0 >> 3)) {
    u64 Ym = (Y0 + Y1) >> 1;
    u128 Ma_ = (u128)a1 * Ym + ((u128)a0 << P), Mb_ = (u128)b1 * Ym + ((u128)b0 << P);
    W Ma = w128(Ma_ >= Mb_ ? Ma_ : Mb_), Mb = w128(Ma_ >= Mb_ ? Mb_ : Ma_);
    W lhs = mul(&Ma, &Ma); lhs = mul(&lhs, &Mb); lhs = mul_small(&lhs, 10);
    W Ymw = w128(Ym), ysq = mul(&Ymw, &Ymw);
    ysq = mul_small(&ysq, 13);
    if (Y1 > 3 * Y0 || cmp_shl(&ysq, 2 * T + P, &lhs) <= 0) {
      node(a1, a0, b1, b0, T, Y0, Ym, depth + 1);
      node(a1, a0, b1, b0, T, Ym, Y1, depth + 1);
      return;
    }
  }
  for (int j = 0;; j++) {
    W fIB = shl(&IB, j + 1), IN0 = add(&IA, &fIB);
    if (cmp_shl(&IN0, LOG_UMIN, &cap) > 0) break;
    int T2 = T + 1 + (j + 1) / 2;
    W s = mul(&IN0, &IN0); s = mul(&s, &IA);
    if (cmp_shl(&Y1sq, 2 * T2 + 1 + P, &s) > 0) {                   /* IN0^2 IA < Y1^2 << ... */
      if (j + 1 >= 62) die("j too large");
      u64 n1 = small((u128)a1 + ((u128)b1 << (j + 1))), n0 = small((u128)a0 + ((u128)b0 << (j + 1)));
      W IN1 = w128((u128)n1 * Y1 + ((u128)n0 << P));
      terminal(n1, n0, a1, a0, &IN0, &IN1);
      if (cmp_shl(&IN0, LOG_UMIN + 1, &cap) <= 0) {
        node(n1, n0, a1, a0, T2, Y0, Y1, depth + 1);
        node(a1, a0, n1, n0, T2, Y0, Y1, depth + 1);
      }
    }
  }
}

int main(int argc, char **argv) {
  if (argc != 5) die("usage: width_search Q_HEX LOG_UMIN SHARD SHARDS");
  memset(&Q, 0, sizeof Q);
  const char *h = argv[1];
  if (h[0] == '0' && (h[1] == 'x' || h[1] == 'X')) h += 2;
  for (; *h; h++) {
    int d = (*h >= '0' && *h <= '9') ? *h - '0' : (*h | 32) - 'a' + 10;
    if (d < 0 || d > 15) die("bad hex");
    Q = shl(&Q, 4); Q.w[0] |= (u64)d;
  }
  LOG_UMIN = atoi(argv[2]); SHARD = atoi(argv[3]); SHARDS = atoi(argv[4]);
  QP = shl(&Q, P); QP_BITS = bitlen(&QP);
  u64 Y_LO = ((u64)17 << P) / 100, Y_HI = (((u64)583 << P) + 99) / 100;
  if (SHARD == 0) {
    W I0 = w128((u128)Y_LO + ((u128)1 << P)), I1 = w128((u128)Y_HI + ((u128)1 << P));
    terminal(1, 1, 1, 0, &I0, &I1);
  }
  node(1, 1, 1, 0, 1, Y_LO, Y_HI, 0);
  node(1, 0, 1, 1, 1, Y_LO, Y_HI, 0);
  printf("N %lld\n", nodes);
  return 0;
}
