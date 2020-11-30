#include <cmath>
#include <immintrin.h>

// Jakub's NearbyINT based function
void CalcBondsNINT(const float *coords1, const float *coords2, const float *box,
                 size_t nvals, float *output) {
  for (size_t i = 0; i < nvals; ++i) {
    float dist = 0.0;
    for (size_t j = 0; j < 3; ++j) {
      float r = coords1[3 * i + j] - coords2[3 * i + j];
      float b = box[j];
      float adj = nearbyintf(r / b);
      r -= adj * b;
      dist += r * r;
    }
    output[i] = sqrtf(dist);
  }
}

// Jakub's NINT + FMA based function
void CalcBondsFMA(const float *coords1, const float *coords2, const float *box,
                  size_t nvals, float *output) {
  for (size_t i = 0; i < nvals; ++i) {
    float dist = 0.0;
    for (size_t j = 0; j < 3; ++j) {
      float r = coords1[3 * i + j] - coords2[3 * i + j];
      float b = box[j];
      float adj = nearbyintf(r / b);
      r = fmaf(-adj, b, r);
      dist = j == 0 ? r * r : fmaf(r, r, dist);
    }
    output[i] = sqrtf(dist);
  }
}

#if defined(__AVX2__) && defined(__FMA__)

typedef struct {
    __m256 a;
    __m256 b;
    __m256 c;
} m256_3;

static inline m256_3
__attribute__((__always_inline__))
mm256_transpose_8x3_ps(__m256 a, __m256 b, __m256 c) {
    /* a = x0y0z0x1y1z1x2y2 */
    /* b = z2x3y3z3x4y4z4x5 */
    /* c = y6z6x7y7z7x8y8z8 */
    
    __m256 m1 = _mm256_blend_ps(a, b, 0xf0);
    __m256 m2 = _mm256_permute2f128_ps(a, c, 0x21);
    __m256 m3 = _mm256_blend_ps(b, c, 0xf0);
    /* m1 = x0y0z0x1x4y4z4x5 */
    /* m2 = y1z1x2y2y5z5x6y6 */
    /* m3 = z2x3y3z3z6x7y7z7 */
    
    __m256 t1 = _mm256_shuffle_ps(m2, m3, _MM_SHUFFLE(2,1,3,2));
    __m256 t2 = _mm256_shuffle_ps(m1, m2, _MM_SHUFFLE(1,0,2,1));
    /* t1 = x2y2x3y3x6y6x7y7 */
    /* t2 = y0z0y1z1y4z4y5z5 */
    
    __m256 x = _mm256_shuffle_ps(m1, t1, _MM_SHUFFLE(2,0,3,0));
    __m256 y = _mm256_shuffle_ps(t2, t1, _MM_SHUFFLE(3,1,2,0));
    __m256 z = _mm256_shuffle_ps(t2, m3, _MM_SHUFFLE(3,0,3,1));
    /* x = x0x1x2x3x4x5x6x7 */
    /* y = y0y1y2y3y4y5y6y7 */
    /* z = z0z1z2z3z4z5z6z7 */
    
    m256_3 res = {x, y, z};
    return res;
}

static inline __m256
__attribute__((__always_inline__))
mm256_periodic_boundary_distance_round(
    __m256 m1, __m256 m2,
    __m256 box, __m256 rbox
) {
    __m256 diff = m1 - m2;
    __m256 shift = diff * rbox;
    __m256 int_shift = _mm256_round_ps(shift, _MM_ROUND_NEAREST | _MM_FROUND_NO_EXC);
    return _mm256_fnmadd_ps(int_shift, box, diff);
}

void CalcBonds256(
    const float *arr1,
    const float *arr2,
    const float *box,
    size_t n,
    float *out
) {
    if (n & 0x7) {
        throw "Number of particles must be a multiple of 8";
    }
    n >>= 3;

    const __m256 *arr1_256 = (const __m256 *) arr1;
    const __m256 *arr2_256 = (const __m256 *) arr2;
    __m256 *out_256 = (__m256 *) out;

    __m256 boxv = {box[0], box[1], box[2], NAN, box[1], box[2], box[0], NAN};
    __m256 rboxv = _mm256_set1_ps(1.0f) / boxv;
    __m256 box1 = _mm256_permute_ps(boxv, _MM_SHUFFLE(0,2,1,0));
    __m256 box2 = _mm256_permute_ps(boxv, _MM_SHUFFLE(2,1,0,2));
    __m256 box3 = _mm256_permute_ps(boxv, _MM_SHUFFLE(1,0,2,1));
    __m256 rbox1 = _mm256_permute_ps(rboxv, _MM_SHUFFLE(0,2,1,0));
    __m256 rbox2 = _mm256_permute_ps(rboxv, _MM_SHUFFLE(2,1,0,2));
    __m256 rbox3 = _mm256_permute_ps(rboxv, _MM_SHUFFLE(1,0,2,1));
    
#pragma unroll(2)
#pragma GCC unroll(2)
    for (size_t i = 0; i < n; ++i) {
        size_t j = i * 3;
        
        __m256 m11 = arr1_256[j];
        __m256 m12 = arr1_256[j+1];
        __m256 m13 = arr1_256[j+2];
        __m256 m21 = arr2_256[j];
        __m256 m22 = arr2_256[j+1];
        __m256 m23 = arr2_256[j+2];
        
        __m256 diffm1 = mm256_periodic_boundary_distance_round(m11, m21, box1, rbox1);
        __m256 diffm2 = mm256_periodic_boundary_distance_round(m12, m22, box2, rbox2);
        __m256 diffm3 = mm256_periodic_boundary_distance_round(m13, m23, box3, rbox3);

        m256_3 transpose_res = mm256_transpose_8x3_ps(diffm1, diffm2, diffm3);
        __m256 x_diff = transpose_res.a;
        __m256 y_diff = transpose_res.b;
        __m256 z_diff = transpose_res.c;

        __m256 dist_sq = x_diff * x_diff;
        dist_sq = _mm256_fmadd_ps(y_diff, y_diff, dist_sq);
        dist_sq = _mm256_fmadd_ps(z_diff, z_diff, dist_sq);

        __m256 dist = _mm256_sqrt_ps(dist_sq);
        out_256[i] = dist;
    }
}

#endif // __AVX2__ && __FMA__
