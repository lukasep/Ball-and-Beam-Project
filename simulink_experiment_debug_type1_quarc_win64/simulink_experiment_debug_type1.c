/*
 * simulink_experiment_debug_type1.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "simulink_experiment_debug_type1".
 *
 * Model version              : 1.1
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C source code generated on : Wed Apr 30 14:04:21 2025
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "simulink_experiment_debug_type1.h"
#include "rtwtypes.h"
#include "simulink_experiment_debug_type1_types.h"
#include <math.h>
#include <emmintrin.h>
#include <string.h>
#include "simulink_experiment_debug_type1_private.h"
#include "rt_nonfinite.h"
#include "simulink_experiment_debug_type1_dt.h"

/* Block signals (default storage) */
B_simulink_experiment_debug_t_T simulink_experiment_debug_typ_B;

/* Block states (default storage) */
DW_simulink_experiment_debug__T simulink_experiment_debug_ty_DW;

/* Real-time model */
static RT_MODEL_simulink_experiment__T simulink_experiment_debug_ty_M_;
RT_MODEL_simulink_experiment__T *const simulink_experiment_debug_ty_M =
  &simulink_experiment_debug_ty_M_;

/* Forward declaration for local functions */
static real_T simulink_experiment_debug__norm(const real_T x[16]);
static void simulink_experiment_debu_mpower(const real_T a[16], real_T b, real_T
  c[16]);
static real_T simulink_experiment_debug__log2(real_T x);
static void simulink_expe_padeApproximation(const real_T A[16], const real_T A2
  [16], const real_T A4[16], const real_T A6[16], int32_T m, real_T F[16]);
static void simulink_exp_recomputeBlockDiag(const real_T A[16], const real_T F
  [16], const int32_T blockFormat[3], real_T b_F[16]);
static real_T simulink_experiment_debug_xnrm2(int32_T n, const real_T x[16],
  int32_T ix0);
static void simulink_experiment_debug_xgerc(int32_T m, int32_T n, real_T alpha1,
  int32_T ix0, const real_T y[4], const real_T A[16], int32_T ia0, real_T b_A[16]);
static void simulink_experiment_debu_xgehrd(const real_T a[16], real_T b_a[16],
  real_T tau[3]);
static real_T simulink_experiment_deb_xnrm2_j(int32_T n, const real_T x[3]);
static void simulink_experiment_d_xzlarfg_j(int32_T n, real_T alpha1, real_T x[3],
  real_T *b_alpha1, real_T *tau);
static void simulink_experiment_deb_xdlanv2(real_T a, real_T b, real_T c, real_T
  d, real_T *rt1r, real_T *rt1i, real_T *rt2r, real_T *rt2i, real_T *b_a, real_T
  *b_b, real_T *b_c, real_T *b_d, real_T *cs, real_T *sn);
static void simulink_experiment_debu_xhseqr(const real_T h[16], const real_T z
  [16], real_T b_h[16], int32_T *info, real_T b_z[16]);
static void simulink_experiment_debug_schur(const real_T A[16], real_T V[16],
  real_T T[16]);
static void simulink_experiment_debug__expm(real_T A[16], real_T F[16]);
static void rate_monotonic_scheduler(void);
time_T rt_SimUpdateDiscreteEvents(
  int_T rtmNumSampTimes, void *rtmTimingData, int_T *rtmSampleHitPtr, int_T
  *rtmPerTaskSampleHits )
{
  rtmSampleHitPtr[1] = rtmStepTask(simulink_experiment_debug_ty_M, 1);
  rtmSampleHitPtr[2] = rtmStepTask(simulink_experiment_debug_ty_M, 2);
  UNUSED_PARAMETER(rtmNumSampTimes);
  UNUSED_PARAMETER(rtmTimingData);
  UNUSED_PARAMETER(rtmPerTaskSampleHits);
  return(-1);
}

/*
 *         This function updates active task flag for each subrate
 *         and rate transition flags for tasks that exchange data.
 *         The function assumes rate-monotonic multitasking scheduler.
 *         The function must be called at model base rate so that
 *         the generated code self-manages all its subrates and rate
 *         transition flags.
 */
static void rate_monotonic_scheduler(void)
{
  /* To ensure a deterministic data transfer between two rates,
   * data is transferred at the priority of a fast task and the frequency
   * of the slow task.  The following flags indicate when the data transfer
   * happens.  That is, a rate interaction flag is set true when both rates
   * will run, and false otherwise.
   */

  /* tid 1 shares data with slower tid rate: 2 */
  if (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[1] == 0) {
    simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2 =
      (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2] == 0);

    /* update PerTaskSampleHits matrix for non-inline sfcn */
    simulink_experiment_debug_ty_M->Timing.perTaskSampleHits[5] =
      simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2;
  }

  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2])++;
  if ((simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2]) > 4) {/* Sample time: [0.01s, 0.0s] */
    simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2] = 0;
  }
}

static real_T simulink_experiment_debug__norm(const real_T x[16])
{
  real_T i;
  real_T s;
  real_T y;
  int32_T b_j;
  int32_T j;
  boolean_T b;
  boolean_T exitg1;
  y = 0.0;
  b_j = 0;
  exitg1 = false;
  while ((!exitg1) && (b_j < 4)) {
    j = b_j + 1;
    i = x[(j - 1) << 2];
    i = fabs(i);
    s = i;
    i = x[((j - 1) << 2) + 1];
    i = fabs(i);
    s += i;
    i = x[((j - 1) << 2) + 2];
    i = fabs(i);
    s += i;
    i = x[((j - 1) << 2) + 3];
    i = fabs(i);
    s += i;
    b = rtIsNaN(s);
    if (b) {
      y = (rtNaN);
      exitg1 = true;
    } else {
      if (s > y) {
        y = s;
      }

      b_j++;
    }
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T tmp;
  real_T tmp_0;
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

static void simulink_experiment_debu_mpower(const real_T a[16], real_T b, real_T
  c[16])
{
  __m128d tmp;
  __m128d tmp_0;
  real_T aBuffer[16];
  real_T a_0[16];
  real_T cBuffer[16];
  real_T c_0[16];
  real_T e;
  real_T ed2;
  int32_T exitg1;
  int32_T i;
  int32_T n;
  int32_T n_0;
  int32_T nb;
  int32_T nbitson;
  boolean_T aBufferInUse;
  boolean_T firstmult;
  boolean_T lsb;
  memcpy(&a_0[0], &a[0], sizeof(real_T) << 4U);
  e = fabs(b);
  if (e <= 2.147483647E+9) {
    e = fabs(b);
    n = (int32_T)e;
    n_0 = n;
    nbitson = 0;
    nb = -1;
    while (n_0 > 0) {
      nb++;
      if (((uint32_T)n_0 & 1U) != 0U) {
        nbitson++;
      }

      n_0 >>= 1;
    }

    if (n <= 2) {
      if (b == 2.0) {
        memcpy(&c_0[0], &a_0[0], sizeof(real_T) << 4U);
        for (nb = 0; nb < 4; nb++) {
          for (i = 0; i <= 2; i += 2) {
            _mm_storeu_pd(&c[i + (nb << 2)], _mm_set1_pd(0.0));
            tmp = _mm_loadu_pd(&a_0[i]);
            tmp = _mm_mul_pd(_mm_set1_pd(c_0[nb << 2]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
            tmp = _mm_loadu_pd(&a_0[i + 4]);
            tmp = _mm_mul_pd(_mm_set1_pd(c_0[(nb << 2) + 1]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
            tmp = _mm_loadu_pd(&a_0[i + 8]);
            tmp = _mm_mul_pd(_mm_set1_pd(c_0[(nb << 2) + 2]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
            tmp = _mm_loadu_pd(&a_0[i + 12]);
            tmp = _mm_mul_pd(_mm_set1_pd(c_0[(nb << 2) + 3]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
          }
        }
      } else {
        firstmult = false;
        for (n_0 = 0; n_0 < 16; n_0++) {
          c_0[n_0] = a_0[n_0];
          if (firstmult) {
            firstmult = true;
          } else {
            e = c_0[n_0];
            aBufferInUse = rtIsNaN(e);
            if (aBufferInUse) {
              firstmult = true;
            }
          }
        }

        if (firstmult) {
          for (nb = 0; nb < 16; nb++) {
            c[nb] = (rtNaN);
          }
        } else {
          for (nb = 0; nb < 16; nb++) {
            c[nb] = 0.0;
          }

          c[0] = 1.0;
          c[5] = 1.0;
          c[10] = 1.0;
          c[15] = 1.0;
        }
      }
    } else {
      firstmult = true;
      aBufferInUse = false;
      lsb = (((uint32_T)nbitson & 1U) != 0U);
      lsb = !lsb;
      nbitson = nb - 1;
      for (n_0 = 0; n_0 <= nbitson; n_0++) {
        if (((uint32_T)n & 1U) != 0U) {
          if (firstmult) {
            firstmult = false;
            if (lsb) {
              if (aBufferInUse) {
                memcpy(&cBuffer[0], &aBuffer[0], sizeof(real_T) << 4U);
              } else {
                memcpy(&cBuffer[0], &a_0[0], sizeof(real_T) << 4U);
              }
            } else if (aBufferInUse) {
              memcpy(&c[0], &aBuffer[0], sizeof(real_T) << 4U);
            } else {
              memcpy(&c[0], &a_0[0], sizeof(real_T) << 4U);
            }
          } else {
            if (aBufferInUse) {
              if (lsb) {
                memcpy(&c_0[0], &aBuffer[0], sizeof(real_T) << 4U);
                for (nb = 0; nb < 4; nb++) {
                  for (i = 0; i < 4; i++) {
                    c[nb + (i << 2)] = 0.0;
                    e = c[(i << 2) + nb];
                    e += c_0[i << 2] * cBuffer[nb];
                    c[nb + (i << 2)] = e;
                    e = c[(i << 2) + nb];
                    e += c_0[(i << 2) + 1] * cBuffer[nb + 4];
                    c[nb + (i << 2)] = e;
                    e = c[(i << 2) + nb];
                    e += c_0[(i << 2) + 2] * cBuffer[nb + 8];
                    c[nb + (i << 2)] = e;
                    e = c[(i << 2) + nb];
                    e += c_0[(i << 2) + 3] * cBuffer[nb + 12];
                    c[nb + (i << 2)] = e;
                  }
                }
              } else {
                memcpy(&c_0[0], &aBuffer[0], sizeof(real_T) << 4U);
                for (nb = 0; nb < 4; nb++) {
                  for (i = 0; i < 4; i++) {
                    cBuffer[nb + (i << 2)] = 0.0;
                    e = cBuffer[(i << 2) + nb];
                    e += c_0[i << 2] * c[nb];
                    cBuffer[nb + (i << 2)] = e;
                    e = cBuffer[(i << 2) + nb];
                    e += c_0[(i << 2) + 1] * c[nb + 4];
                    cBuffer[nb + (i << 2)] = e;
                    e = cBuffer[(i << 2) + nb];
                    e += c_0[(i << 2) + 2] * c[nb + 8];
                    cBuffer[nb + (i << 2)] = e;
                    e = cBuffer[(i << 2) + nb];
                    e += c_0[(i << 2) + 3] * c[nb + 12];
                    cBuffer[nb + (i << 2)] = e;
                  }
                }
              }
            } else if (lsb) {
              memcpy(&c_0[0], &a_0[0], sizeof(real_T) << 4U);
              for (nb = 0; nb < 4; nb++) {
                for (i = 0; i < 4; i++) {
                  c[nb + (i << 2)] = 0.0;
                  e = c[(i << 2) + nb];
                  e += c_0[i << 2] * cBuffer[nb];
                  c[nb + (i << 2)] = e;
                  e = c[(i << 2) + nb];
                  e += c_0[(i << 2) + 1] * cBuffer[nb + 4];
                  c[nb + (i << 2)] = e;
                  e = c[(i << 2) + nb];
                  e += c_0[(i << 2) + 2] * cBuffer[nb + 8];
                  c[nb + (i << 2)] = e;
                  e = c[(i << 2) + nb];
                  e += c_0[(i << 2) + 3] * cBuffer[nb + 12];
                  c[nb + (i << 2)] = e;
                }
              }
            } else {
              memcpy(&c_0[0], &a_0[0], sizeof(real_T) << 4U);
              for (nb = 0; nb < 4; nb++) {
                for (i = 0; i < 4; i++) {
                  cBuffer[nb + (i << 2)] = 0.0;
                  e = cBuffer[(i << 2) + nb];
                  e += c_0[i << 2] * c[nb];
                  cBuffer[nb + (i << 2)] = e;
                  e = cBuffer[(i << 2) + nb];
                  e += c_0[(i << 2) + 1] * c[nb + 4];
                  cBuffer[nb + (i << 2)] = e;
                  e = cBuffer[(i << 2) + nb];
                  e += c_0[(i << 2) + 2] * c[nb + 8];
                  cBuffer[nb + (i << 2)] = e;
                  e = cBuffer[(i << 2) + nb];
                  e += c_0[(i << 2) + 3] * c[nb + 12];
                  cBuffer[nb + (i << 2)] = e;
                }
              }
            }

            lsb = !lsb;
          }
        }

        n >>= 1;
        if (aBufferInUse) {
          memcpy(&c_0[0], &aBuffer[0], sizeof(real_T) << 4U);
          for (nb = 0; nb < 4; nb++) {
            for (i = 0; i < 4; i++) {
              a_0[nb + (i << 2)] = 0.0;
              ed2 = a_0[(i << 2) + nb];
              ed2 += c_0[i << 2] * aBuffer[nb];
              a_0[nb + (i << 2)] = ed2;
              ed2 = a_0[(i << 2) + nb];
              ed2 += c_0[(i << 2) + 1] * aBuffer[nb + 4];
              a_0[nb + (i << 2)] = ed2;
              ed2 = a_0[(i << 2) + nb];
              ed2 += c_0[(i << 2) + 2] * aBuffer[nb + 8];
              a_0[nb + (i << 2)] = ed2;
              ed2 = a_0[(i << 2) + nb];
              ed2 += c_0[(i << 2) + 3] * aBuffer[nb + 12];
              a_0[nb + (i << 2)] = ed2;
            }
          }
        } else {
          memcpy(&c_0[0], &a_0[0], sizeof(real_T) << 4U);
          for (nb = 0; nb < 4; nb++) {
            for (i = 0; i < 4; i++) {
              aBuffer[nb + (i << 2)] = 0.0;
              e = aBuffer[(i << 2) + nb];
              e += c_0[i << 2] * a_0[nb];
              aBuffer[nb + (i << 2)] = e;
              e = aBuffer[(i << 2) + nb];
              e += c_0[(i << 2) + 1] * a_0[nb + 4];
              aBuffer[nb + (i << 2)] = e;
              e = aBuffer[(i << 2) + nb];
              e += c_0[(i << 2) + 2] * a_0[nb + 8];
              aBuffer[nb + (i << 2)] = e;
              e = aBuffer[(i << 2) + nb];
              e += c_0[(i << 2) + 3] * a_0[nb + 12];
              aBuffer[nb + (i << 2)] = e;
            }
          }
        }

        aBufferInUse = !aBufferInUse;
      }

      if (firstmult) {
        if (aBufferInUse) {
          memcpy(&c[0], &aBuffer[0], sizeof(real_T) << 4U);
        } else {
          memcpy(&c[0], &a_0[0], sizeof(real_T) << 4U);
        }
      } else if (aBufferInUse) {
        for (nb = 0; nb < 4; nb++) {
          for (i = 0; i <= 2; i += 2) {
            _mm_storeu_pd(&c[i + (nb << 2)], _mm_set1_pd(0.0));
            tmp = _mm_loadu_pd(&cBuffer[i]);
            tmp = _mm_mul_pd(_mm_set1_pd(aBuffer[nb << 2]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
            tmp = _mm_loadu_pd(&cBuffer[i + 4]);
            tmp = _mm_mul_pd(_mm_set1_pd(aBuffer[(nb << 2) + 1]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
            tmp = _mm_loadu_pd(&cBuffer[i + 8]);
            tmp = _mm_mul_pd(_mm_set1_pd(aBuffer[(nb << 2) + 2]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
            tmp = _mm_loadu_pd(&cBuffer[i + 12]);
            tmp = _mm_mul_pd(_mm_set1_pd(aBuffer[(nb << 2) + 3]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
          }
        }
      } else {
        for (nb = 0; nb < 4; nb++) {
          for (i = 0; i <= 2; i += 2) {
            _mm_storeu_pd(&c[i + (nb << 2)], _mm_set1_pd(0.0));
            tmp = _mm_loadu_pd(&cBuffer[i]);
            tmp = _mm_mul_pd(_mm_set1_pd(a_0[nb << 2]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
            tmp = _mm_loadu_pd(&cBuffer[i + 4]);
            tmp = _mm_mul_pd(_mm_set1_pd(a_0[(nb << 2) + 1]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
            tmp = _mm_loadu_pd(&cBuffer[i + 8]);
            tmp = _mm_mul_pd(_mm_set1_pd(a_0[(nb << 2) + 2]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
            tmp = _mm_loadu_pd(&cBuffer[i + 12]);
            tmp = _mm_mul_pd(_mm_set1_pd(a_0[(nb << 2) + 3]), tmp);
            tmp_0 = _mm_loadu_pd(&c[(nb << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&c[i + (nb << 2)], tmp);
          }
        }
      }
    }
  } else {
    aBufferInUse = rtIsInf(b);
    firstmult = !aBufferInUse;
    aBufferInUse = rtIsNaN(b);
    aBufferInUse = !aBufferInUse;
    aBufferInUse = (firstmult && aBufferInUse);
    if (aBufferInUse) {
      e = fabs(b);
      firstmult = true;
      do {
        exitg1 = 0;
        ed2 = floor(e / 2.0);
        if (2.0 * ed2 != e) {
          if (firstmult) {
            memcpy(&c[0], &a_0[0], sizeof(real_T) << 4U);
            firstmult = false;
          } else {
            for (nb = 0; nb < 4; nb++) {
              for (i = 0; i < 4; i++) {
                c_0[nb + (i << 2)] = 0.0;
                e = c_0[(i << 2) + nb];
                e += a_0[i << 2] * c[nb];
                c_0[nb + (i << 2)] = e;
                e = c_0[(i << 2) + nb];
                e += a_0[(i << 2) + 1] * c[nb + 4];
                c_0[nb + (i << 2)] = e;
                e = c_0[(i << 2) + nb];
                e += a_0[(i << 2) + 2] * c[nb + 8];
                c_0[nb + (i << 2)] = e;
                e = c_0[(i << 2) + nb];
                e += a_0[(i << 2) + 3] * c[nb + 12];
                c_0[nb + (i << 2)] = e;
              }
            }

            memcpy(&c[0], &c_0[0], sizeof(real_T) << 4U);
          }
        }

        if (ed2 == 0.0) {
          exitg1 = 1;
        } else {
          e = ed2;
          memcpy(&c_0[0], &a_0[0], sizeof(real_T) << 4U);
          for (nb = 0; nb < 4; nb++) {
            for (i = 0; i < 4; i++) {
              cBuffer[nb + (i << 2)] = 0.0;
              ed2 = cBuffer[(i << 2) + nb];
              ed2 += c_0[i << 2] * a_0[nb];
              cBuffer[nb + (i << 2)] = ed2;
              ed2 = cBuffer[(i << 2) + nb];
              ed2 += c_0[(i << 2) + 1] * a_0[nb + 4];
              cBuffer[nb + (i << 2)] = ed2;
              ed2 = cBuffer[(i << 2) + nb];
              ed2 += c_0[(i << 2) + 2] * a_0[nb + 8];
              cBuffer[nb + (i << 2)] = ed2;
              ed2 = cBuffer[(i << 2) + nb];
              ed2 += c_0[(i << 2) + 3] * a_0[nb + 12];
              cBuffer[nb + (i << 2)] = ed2;
            }
          }

          memcpy(&a_0[0], &cBuffer[0], sizeof(real_T) << 4U);
        }
      } while (exitg1 == 0);
    } else {
      for (nb = 0; nb < 16; nb++) {
        c[nb] = (rtNaN);
      }
    }
  }
}

static real_T simulink_experiment_debug__log2(real_T x)
{
  real_T f;
  real_T fdbl;
  int32_T eint;
  boolean_T b;
  boolean_T c;
  if (x == 0.0) {
    f = (rtMinusInf);
  } else {
    b = rtIsInf(x);
    c = !b;
    b = rtIsNaN(x);
    b = !b;
    b = (c && b);
    if (b) {
      b = rtIsInf(x);
      c = !b;
      b = rtIsNaN(x);
      b = !b;
      b = (c && b);
      if (b) {
        fdbl = frexp(x, &eint);
      } else {
        fdbl = x;
        eint = 0;
      }

      if (fdbl == 0.5) {
        f = (real_T)eint - 1.0;
      } else if ((eint == 1) && (fdbl < 0.75)) {
        f = log(2.0 * fdbl) / 0.69314718055994529;
      } else {
        f = log(fdbl) / 0.69314718055994529 + (real_T)eint;
      }
    } else {
      f = x;
    }
  }

  return f;
}

static void simulink_expe_padeApproximation(const real_T A[16], const real_T A2
  [16], const real_T A4[16], const real_T A6[16], int32_T m, real_T F[16])
{
  __m128d tmp;
  __m128d tmp_0;
  __m128d tmp_1;
  __m128d tmp_2;
  __m128d tmp_3;
  __m128d tmp_4;
  real_T V[16];
  real_T y[16];
  real_T d;
  real_T s;
  real_T x;
  int32_T a;
  int32_T g_k;
  int32_T ix;
  int32_T jBcol;
  int32_T jj;
  int32_T jm1;
  int32_T jp1j;
  int32_T jpiv;
  int32_T jpiv_offset;
  int32_T jrow;
  int32_T mmj;
  int8_T ipiv[4];
  if (m == 3) {
    memcpy(&F[0], &A2[0], sizeof(real_T) << 4U);
    F[0] += 60.0;
    F[5] += 60.0;
    F[10] += 60.0;
    F[15] += 60.0;
    for (ix = 0; ix < 4; ix++) {
      for (g_k = 0; g_k <= 2; g_k += 2) {
        _mm_storeu_pd(&y[g_k + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A[g_k]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[ix << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&y[ix]);
      _mm_storeu_pd(&F[ix], tmp);
      tmp = _mm_loadu_pd(&A2[ix]);
      tmp = _mm_mul_pd(_mm_set1_pd(12.0), tmp);
      _mm_storeu_pd(&V[ix], tmp);
    }

    d = 120.0;
  } else if (m == 5) {
    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&A2[ix]);
      tmp = _mm_mul_pd(_mm_set1_pd(420.0), tmp);
      tmp_0 = _mm_loadu_pd(&A4[ix]);
      tmp = _mm_add_pd(tmp_0, tmp);
      _mm_storeu_pd(&F[ix], tmp);
    }

    F[0] += 15120.0;
    F[5] += 15120.0;
    F[10] += 15120.0;
    F[15] += 15120.0;
    for (ix = 0; ix < 4; ix++) {
      for (g_k = 0; g_k <= 2; g_k += 2) {
        _mm_storeu_pd(&y[g_k + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A[g_k]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[ix << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&y[ix]);
      _mm_storeu_pd(&F[ix], tmp);
      tmp = _mm_loadu_pd(&A4[ix]);
      tmp = _mm_mul_pd(_mm_set1_pd(30.0), tmp);
      tmp_0 = _mm_loadu_pd(&A2[ix]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(3360.0), tmp_0);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&V[ix], tmp);
    }

    d = 30240.0;
  } else if (m == 7) {
    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&A4[ix]);
      tmp = _mm_mul_pd(_mm_set1_pd(1512.0), tmp);
      tmp_0 = _mm_loadu_pd(&A2[ix]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(277200.0), tmp_0);
      tmp_1 = _mm_loadu_pd(&A6[ix]);
      tmp = _mm_add_pd(tmp_1, tmp);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&F[ix], tmp);
    }

    F[0] += 8.64864E+6;
    F[5] += 8.64864E+6;
    F[10] += 8.64864E+6;
    F[15] += 8.64864E+6;
    for (ix = 0; ix < 4; ix++) {
      for (g_k = 0; g_k <= 2; g_k += 2) {
        _mm_storeu_pd(&y[g_k + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A[g_k]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[ix << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&y[ix]);
      _mm_storeu_pd(&F[ix], tmp);
      tmp = _mm_loadu_pd(&A6[ix]);
      tmp = _mm_mul_pd(_mm_set1_pd(56.0), tmp);
      tmp_0 = _mm_loadu_pd(&A4[ix]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(25200.0), tmp_0);
      tmp_1 = _mm_loadu_pd(&A2[ix]);
      tmp_1 = _mm_mul_pd(_mm_set1_pd(1.99584E+6), tmp_1);
      tmp = _mm_add_pd(tmp, tmp_0);
      tmp = _mm_add_pd(tmp, tmp_1);
      _mm_storeu_pd(&V[ix], tmp);
    }

    d = 1.729728E+7;
  } else if (m == 9) {
    for (ix = 0; ix < 4; ix++) {
      for (g_k = 0; g_k <= 2; g_k += 2) {
        _mm_storeu_pd(&V[g_k + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A6[g_k]);
        tmp = _mm_mul_pd(_mm_set1_pd(A2[ix << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A6[g_k + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(A2[(ix << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A6[g_k + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(A2[(ix << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A6[g_k + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(A2[(ix << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&A6[ix]);
      tmp = _mm_mul_pd(_mm_set1_pd(3960.0), tmp);
      tmp_0 = _mm_loadu_pd(&A4[ix]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(2.16216E+6), tmp_0);
      tmp_1 = _mm_loadu_pd(&A2[ix]);
      tmp_1 = _mm_mul_pd(_mm_set1_pd(3.027024E+8), tmp_1);
      tmp_2 = _mm_loadu_pd(&V[ix]);
      tmp = _mm_add_pd(tmp_2, tmp);
      tmp = _mm_add_pd(tmp, tmp_0);
      tmp = _mm_add_pd(tmp, tmp_1);
      _mm_storeu_pd(&F[ix], tmp);
    }

    F[0] += 8.8216128E+9;
    F[5] += 8.8216128E+9;
    F[10] += 8.8216128E+9;
    F[15] += 8.8216128E+9;
    for (ix = 0; ix < 4; ix++) {
      for (g_k = 0; g_k <= 2; g_k += 2) {
        _mm_storeu_pd(&y[g_k + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A[g_k]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[ix << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(F[(ix << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&y[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&y[g_k + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&y[ix]);
      _mm_storeu_pd(&F[ix], tmp);
      tmp = _mm_loadu_pd(&V[ix]);
      tmp = _mm_mul_pd(_mm_set1_pd(90.0), tmp);
      tmp_0 = _mm_loadu_pd(&A6[ix]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(110880.0), tmp_0);
      tmp_1 = _mm_loadu_pd(&A4[ix]);
      tmp_1 = _mm_mul_pd(_mm_set1_pd(3.027024E+7), tmp_1);
      tmp_2 = _mm_loadu_pd(&A2[ix]);
      tmp_2 = _mm_mul_pd(_mm_set1_pd(2.0756736E+9), tmp_2);
      tmp = _mm_add_pd(tmp, tmp_0);
      tmp = _mm_add_pd(tmp, tmp_1);
      tmp = _mm_add_pd(tmp, tmp_2);
      _mm_storeu_pd(&V[ix], tmp);
    }

    d = 1.76432256E+10;
  } else {
    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&A6[ix]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(3.352212864E+10), tmp);
      tmp_1 = _mm_loadu_pd(&A4[ix]);
      tmp_2 = _mm_mul_pd(_mm_set1_pd(1.05594705216E+13), tmp_1);
      tmp_3 = _mm_loadu_pd(&A2[ix]);
      tmp_4 = _mm_mul_pd(_mm_set1_pd(1.1873537964288E+15), tmp_3);
      tmp_0 = _mm_add_pd(tmp_0, tmp_2);
      tmp_0 = _mm_add_pd(tmp_0, tmp_4);
      tmp_1 = _mm_mul_pd(_mm_set1_pd(16380.0), tmp_1);
      tmp_2 = _mm_mul_pd(_mm_set1_pd(4.08408E+7), tmp_3);
      tmp = _mm_add_pd(tmp, tmp_1);
      tmp = _mm_add_pd(tmp, tmp_2);
      _mm_storeu_pd(&F[ix], tmp_0);
      _mm_storeu_pd(&y[ix], tmp);
    }

    F[0] += 3.238237626624E+16;
    F[5] += 3.238237626624E+16;
    F[10] += 3.238237626624E+16;
    F[15] += 3.238237626624E+16;
    for (ix = 0; ix < 4; ix++) {
      for (g_k = 0; g_k <= 2; g_k += 2) {
        _mm_storeu_pd(&V[g_k + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A6[g_k]);
        tmp = _mm_mul_pd(_mm_set1_pd(y[ix << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A6[g_k + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(y[(ix << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A6[g_k + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(y[(ix << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A6[g_k + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(y[(ix << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&V[ix]);
      tmp_0 = _mm_loadu_pd(&F[ix]);
      tmp = _mm_add_pd(tmp, tmp_0);
      _mm_storeu_pd(&V[ix], tmp);
    }

    for (ix = 0; ix < 4; ix++) {
      for (g_k = 0; g_k <= 2; g_k += 2) {
        _mm_storeu_pd(&F[g_k + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A[g_k]);
        tmp = _mm_mul_pd(_mm_set1_pd(V[ix << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&F[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&F[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(V[(ix << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&F[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&F[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(V[(ix << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&F[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&F[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A[g_k + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(V[(ix << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&F[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&F[g_k + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&A6[ix]);
      tmp = _mm_mul_pd(_mm_set1_pd(182.0), tmp);
      tmp_0 = _mm_loadu_pd(&A4[ix]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(960960.0), tmp_0);
      tmp_1 = _mm_loadu_pd(&A2[ix]);
      tmp_1 = _mm_mul_pd(_mm_set1_pd(1.32324192E+9), tmp_1);
      tmp = _mm_add_pd(tmp, tmp_0);
      tmp = _mm_add_pd(tmp, tmp_1);
      _mm_storeu_pd(&y[ix], tmp);
    }

    for (ix = 0; ix < 4; ix++) {
      for (g_k = 0; g_k <= 2; g_k += 2) {
        _mm_storeu_pd(&V[g_k + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A6[g_k]);
        tmp = _mm_mul_pd(_mm_set1_pd(y[ix << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A6[g_k + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(y[(ix << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A6[g_k + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(y[(ix << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A6[g_k + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(y[(ix << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&V[(ix << 2) + g_k]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&V[g_k + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&A6[ix]);
      tmp = _mm_mul_pd(_mm_set1_pd(6.704425728E+11), tmp);
      tmp_0 = _mm_loadu_pd(&A4[ix]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(1.29060195264E+14), tmp_0);
      tmp_1 = _mm_loadu_pd(&A2[ix]);
      tmp_1 = _mm_mul_pd(_mm_set1_pd(7.7717703038976E+15), tmp_1);
      tmp_2 = _mm_loadu_pd(&V[ix]);
      tmp = _mm_add_pd(tmp_2, tmp);
      tmp = _mm_add_pd(tmp, tmp_0);
      tmp = _mm_add_pd(tmp, tmp_1);
      _mm_storeu_pd(&V[ix], tmp);
    }

    d = 6.476475253248E+16;
  }

  V[0] += d;
  V[5] += d;
  V[10] += d;
  V[15] += d;
  for (g_k = 0; g_k < 16; g_k++) {
    d = (real_T)g_k + 1.0;
    V[(int32_T)d - 1] -= F[(int32_T)d - 1];
    F[(int32_T)d - 1] *= 2.0;
  }

  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  ipiv[3] = 4;
  for (g_k = 0; g_k < 3; g_k++) {
    jBcol = g_k;
    jm1 = jBcol;
    mmj = 3 - jBcol;
    jrow = jm1 * 5;
    a = 1;
    jj = jrow + 1;
    jp1j = jj + 1;
    jrow = mmj + 1;
    memcpy(&y[0], &V[0], sizeof(real_T) << 4U);
    ix = jj - 1;
    x = y[jj - 1];
    s = fabs(x);
    d = s;
    for (jpiv_offset = 2; jpiv_offset <= jrow; jpiv_offset++) {
      ix++;
      x = y[ix];
      s = fabs(x);
      if (s > d) {
        a = jpiv_offset;
        d = s;
      }
    }

    jpiv_offset = a - 1;
    jpiv = (jj + jpiv_offset) - 1;
    if (V[jpiv] != 0.0) {
      if (jpiv_offset != 0) {
        jrow = (jBcol + jpiv_offset) + 1;
        ipiv[jBcol] = (int8_T)jrow;
        jrow = jm1;
        ix = jrow + jpiv_offset;
        d = V[jrow];
        V[jrow] = V[ix];
        V[ix] = d;
        jrow += 4;
        ix += 4;
        d = V[jrow];
        V[jrow] = V[ix];
        V[ix] = d;
        jrow += 4;
        ix += 4;
        d = V[jrow];
        V[jrow] = V[ix];
        V[ix] = d;
        jrow += 4;
        ix += 4;
        d = V[jrow];
        V[jrow] = V[ix];
        V[ix] = d;
      }

      jrow = mmj;
      jpiv = (jp1j + jrow) - 1;
      for (ix = jp1j; ix <= jpiv; ix++) {
        x = V[ix - 1];
        s = V[jj - 1];
        d = x / s;
        V[ix - 1] = d;
      }
    }

    jrow = 3 - jBcol;
    jpiv_offset = jj + 3;
    jBcol = jj + 4;
    jj = jBcol;
    jpiv = jrow - 1;
    for (jBcol = 0; jBcol <= jpiv; jBcol++) {
      d = V[jpiv_offset];
      if (d != 0.0) {
        d = -d;
        ix = jp1j - 1;
        jrow = jj;
        jm1 = mmj + jj;
        for (a = jrow + 1; a <= jm1; a++) {
          V[a - 1] += V[ix] * d;
          ix++;
        }
      }

      jpiv_offset += 4;
      jj += 4;
    }
  }

  for (mmj = 0; mmj < 3; mmj++) {
    ix = mmj;
    if (ix + 1 != ipiv[ix]) {
      jpiv_offset = ipiv[ix] - 1;
      d = F[ix];
      F[ix] = F[jpiv_offset];
      F[jpiv_offset] = d;
      d = F[ix + 4];
      F[ix + 4] = F[jpiv_offset + 4];
      F[jpiv_offset + 4] = d;
      d = F[ix + 8];
      F[ix + 8] = F[jpiv_offset + 8];
      F[jpiv_offset + 8] = d;
      d = F[ix + 12];
      F[ix + 12] = F[jpiv_offset + 12];
      F[jpiv_offset + 12] = d;
    }
  }

  memcpy(&y[0], &V[0], sizeof(real_T) << 4U);
  for (g_k = 0; g_k < 4; g_k++) {
    jBcol = g_k;
    jBcol = (jBcol << 2) - 1;
    for (mmj = 0; mmj < 4; mmj++) {
      jpiv_offset = mmj;
      jp1j = (jpiv_offset << 2) - 1;
      if (F[(jpiv_offset + jBcol) + 1] != 0.0) {
        jpiv = jpiv_offset;
        for (ix = jpiv + 2; ix < 5; ix++) {
          F[ix + jBcol] -= F[(jpiv_offset + jBcol) + 1] * y[ix + jp1j];
        }
      }
    }
  }

  for (g_k = 0; g_k < 4; g_k++) {
    jBcol = g_k;
    jBcol <<= 2;
    for (jpiv_offset = 3; jpiv_offset >= 0; jpiv_offset--) {
      jp1j = jpiv_offset << 2;
      if (F[jpiv_offset + jBcol] != 0.0) {
        F[jpiv_offset + jBcol] /= V[jpiv_offset + jp1j];
        jpiv = jpiv_offset;
        jrow = jpiv - 1;
        for (mmj = 0; mmj <= jrow; mmj++) {
          ix = mmj;
          F[ix + jBcol] -= F[jpiv_offset + jBcol] * V[ix + jp1j];
        }
      }
    }
  }

  F[0]++;
  F[5]++;
  F[10]++;
  F[15]++;
}

static void simulink_exp_recomputeBlockDiag(const real_T A[16], const real_T F
  [16], const int32_T blockFormat[3], real_T b_F[16])
{
  real_T a;
  real_T b;
  real_T c;
  real_T coshdelta;
  real_T delta;
  real_T sinchdelta;
  int32_T blockType;
  memcpy(&b_F[0], &F[0], sizeof(real_T) << 4U);
  blockType = blockFormat[0];
  if (blockType != 0) {
    if (blockType == 1) {
      sinchdelta = A[0];
      delta = A[5];
      a = sinchdelta;
      a = exp(a);
      b = delta;
      b = exp(b);
      c = (sinchdelta + delta) / 2.0;
      coshdelta = sinchdelta - delta;
      coshdelta = fabs(coshdelta);
      coshdelta /= 2.0;
      if ((c >= coshdelta) || rtIsNaN(coshdelta)) {
        coshdelta = c;
      }

      if (coshdelta < 709.782712893384) {
        c = exp(c);
        coshdelta = (delta - sinchdelta) / 2.0;
        if (coshdelta == 0.0) {
          coshdelta = 1.0;
        } else {
          sinchdelta = coshdelta;
          sinchdelta = sinh(sinchdelta);
          coshdelta = sinchdelta / coshdelta;
        }

        coshdelta *= A[4] * c;
      } else {
        coshdelta = (b - a) * A[4] / (delta - sinchdelta);
      }

      b_F[0] = a;
      b_F[4] = coshdelta;
      b_F[5] = b;
    } else {
      a = A[0];
      b = A[4];
      c = A[1];
      coshdelta = b * c;
      delta = fabs(coshdelta);
      delta = sqrt(delta);
      a = exp(a);
      coshdelta = delta;
      coshdelta = cos(coshdelta);
      if (delta == 0.0) {
        sinchdelta = 1.0;
      } else {
        sinchdelta = delta;
        sinchdelta = sin(sinchdelta);
        sinchdelta /= delta;
      }

      b_F[0] = a * coshdelta;
      b_F[1] = a * c * sinchdelta;
      b_F[4] = a * b * sinchdelta;
      b_F[5] = b_F[0];
    }
  }

  blockType = blockFormat[1];
  if (blockType != 0) {
    if (blockType == 1) {
      sinchdelta = A[5];
      delta = A[10];
      a = sinchdelta;
      a = exp(a);
      b = delta;
      b = exp(b);
      c = (sinchdelta + delta) / 2.0;
      coshdelta = sinchdelta - delta;
      coshdelta = fabs(coshdelta);
      coshdelta /= 2.0;
      if ((c >= coshdelta) || rtIsNaN(coshdelta)) {
        coshdelta = c;
      }

      if (coshdelta < 709.782712893384) {
        c = exp(c);
        coshdelta = (delta - sinchdelta) / 2.0;
        if (coshdelta == 0.0) {
          coshdelta = 1.0;
        } else {
          sinchdelta = coshdelta;
          sinchdelta = sinh(sinchdelta);
          coshdelta = sinchdelta / coshdelta;
        }

        coshdelta *= A[9] * c;
      } else {
        coshdelta = (b - a) * A[9] / (delta - sinchdelta);
      }

      b_F[5] = a;
      b_F[9] = coshdelta;
      b_F[10] = b;
    } else {
      a = A[5];
      b = A[9];
      c = A[6];
      coshdelta = b * c;
      delta = fabs(coshdelta);
      delta = sqrt(delta);
      a = exp(a);
      coshdelta = delta;
      coshdelta = cos(coshdelta);
      if (delta == 0.0) {
        sinchdelta = 1.0;
      } else {
        sinchdelta = delta;
        sinchdelta = sin(sinchdelta);
        sinchdelta /= delta;
      }

      b_F[5] = a * coshdelta;
      b_F[6] = a * c * sinchdelta;
      b_F[9] = a * b * sinchdelta;
      b_F[10] = b_F[5];
    }
  }

  blockType = blockFormat[2];
  if (blockType != 0) {
    if (blockType == 1) {
      sinchdelta = A[10];
      delta = A[15];
      a = sinchdelta;
      a = exp(a);
      b = delta;
      b = exp(b);
      c = (sinchdelta + delta) / 2.0;
      coshdelta = sinchdelta - delta;
      coshdelta = fabs(coshdelta);
      coshdelta /= 2.0;
      if ((c >= coshdelta) || rtIsNaN(coshdelta)) {
        coshdelta = c;
      }

      if (coshdelta < 709.782712893384) {
        c = exp(c);
        coshdelta = (delta - sinchdelta) / 2.0;
        if (coshdelta == 0.0) {
          coshdelta = 1.0;
        } else {
          sinchdelta = coshdelta;
          sinchdelta = sinh(sinchdelta);
          coshdelta = sinchdelta / coshdelta;
        }

        coshdelta *= A[14] * c;
      } else {
        coshdelta = (b - a) * A[14] / (delta - sinchdelta);
      }

      b_F[10] = a;
      b_F[14] = coshdelta;
      b_F[15] = b;
    } else {
      a = A[10];
      b = A[14];
      c = A[11];
      coshdelta = b * c;
      delta = fabs(coshdelta);
      delta = sqrt(delta);
      a = exp(a);
      coshdelta = delta;
      coshdelta = cos(coshdelta);
      if (delta == 0.0) {
        sinchdelta = 1.0;
      } else {
        sinchdelta = delta;
        sinchdelta = sin(sinchdelta);
        sinchdelta /= delta;
      }

      b_F[10] = a * coshdelta;
      b_F[11] = a * c * sinchdelta;
      b_F[14] = a * b * sinchdelta;
      b_F[15] = b_F[10];
    }
  }

  if (blockFormat[2] == 0) {
    coshdelta = A[15];
    coshdelta = exp(coshdelta);
    b_F[15] = coshdelta;
  }
}

static real_T simulink_experiment_debug_xnrm2(int32_T n, const real_T x[16],
  int32_T ix0)
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  int32_T k;
  int32_T kend;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    absxk = x[ix0 - 1];
    y = fabs(absxk);
  } else {
    scale = 3.3121686421112381E-170;
    kend = ix0 + 1;
    for (k = ix0; k <= kend; k++) {
      absxk = x[k - 1];
      absxk = fabs(absxk);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T a;
  real_T b;
  real_T y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = sqrt(a * a + 1.0) * b;
  } else if (a > b) {
    b /= a;
    y = sqrt(b * b + 1.0) * a;
  } else if (rtIsNaN(b)) {
    y = (rtNaN);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

static void simulink_experiment_debug_xgerc(int32_T m, int32_T n, real_T alpha1,
  int32_T ix0, const real_T y[4], const real_T A[16], int32_T ia0, real_T b_A[16])
{
  real_T yjy;
  int32_T b;
  int32_T d;
  int32_T ijA;
  int32_T ix;
  int32_T j;
  int32_T jA;
  int32_T jy;
  memcpy(&b_A[0], &A[0], sizeof(real_T) << 4U);
  if (!(alpha1 == 0.0)) {
    jA = ia0;
    jy = 0;
    b = (uint8_T)n - 1;
    for (j = 0; j <= b; j++) {
      yjy = y[jy];
      if (yjy != 0.0) {
        yjy *= alpha1;
        ix = ix0 - 1;
        d = (m + jA) - 1;
        for (ijA = jA; ijA <= d; ijA++) {
          b_A[ijA - 1] += b_A[ix] * yjy;
          ix++;
        }
      }

      jy++;
      jA += 4;
    }
  }
}

static void simulink_experiment_debu_xgehrd(const real_T a[16], real_T b_a[16],
  real_T tau[3])
{
  __m128d tmp;
  real_T A[16];
  real_T x[16];
  real_T work[4];
  real_T alpha1;
  real_T b;
  real_T c_0;
  real_T xnorm;
  real_T y;
  int32_T b_i;
  int32_T c;
  int32_T e;
  int32_T exitg1;
  int32_T i;
  int32_T ia;
  int32_T im1;
  int32_T im1n;
  int32_T in;
  int32_T ip1;
  int32_T iy;
  int32_T jA;
  int32_T jy;
  int32_T knt;
  int32_T mm1;
  int32_T nm1;
  int32_T rowright;
  boolean_T exitg2;
  memcpy(&b_a[0], &a[0], sizeof(real_T) << 4U);
  work[0] = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  work[3] = 0.0;
  for (b_i = 0; b_i < 3; b_i++) {
    im1 = b_i;
    ip1 = b_i + 2;
    rowright = 4;
    im1n = im1 << 2;
    in = (b_i + 1) << 2;
    alpha1 = b_a[((b_i << 2) + ip1) - 1];
    c = b_i + 3;
    jA = im1 << 2;
    if (c <= 4) {
      rowright = c;
    }

    im1 = rowright + jA;
    c = 2 - b_i;
    b = 0.0;
    nm1 = c - 1;
    xnorm = simulink_experiment_debug_xnrm2(nm1 + 1, b_a, im1);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(alpha1, xnorm);
      if (alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      y = fabs(xnorm);
      if (y < 1.0020841800044864E-292) {
        knt = -1;
        do {
          knt++;
          c = nm1;
          rowright = im1 + c;
          jA = ((((rowright - im1) + 1) / 2) << 1) + im1;
          jy = jA - 2;
          for (c = im1; c <= jy; c += 2) {
            tmp = _mm_loadu_pd(&b_a[c - 1]);
            tmp = _mm_mul_pd(tmp, _mm_set1_pd(9.9792015476736E+291));
            _mm_storeu_pd(&b_a[c - 1], tmp);
          }

          for (c = jA; c <= rowright; c++) {
            b_a[c - 1] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          alpha1 *= 9.9792015476736E+291;
          y = fabs(xnorm);
        } while ((y < 1.0020841800044864E-292) && (knt + 1 < 20));

        xnorm = simulink_experiment_debug_xnrm2(nm1 + 1, b_a, im1);
        xnorm = rt_hypotd_snf(alpha1, xnorm);
        if (alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        b = xnorm - alpha1;
        b /= xnorm;
        y = alpha1 - xnorm;
        alpha1 = 1.0 / y;
        c = nm1;
        rowright = im1 + c;
        jA = ((((rowright - im1) + 1) / 2) << 1) + im1;
        jy = jA - 2;
        for (c = im1; c <= jy; c += 2) {
          tmp = _mm_loadu_pd(&b_a[c - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(alpha1));
          _mm_storeu_pd(&b_a[c - 1], tmp);
        }

        for (c = jA; c <= rowright; c++) {
          b_a[c - 1] *= alpha1;
        }

        for (c = 0; c <= knt; c++) {
          xnorm *= 1.0020841800044864E-292;
        }

        alpha1 = xnorm;
      } else {
        b = xnorm - alpha1;
        b /= xnorm;
        y = alpha1 - xnorm;
        alpha1 = 1.0 / y;
        c = nm1;
        rowright = im1 + c;
        jA = ((((rowright - im1) + 1) / 2) << 1) + im1;
        jy = jA - 2;
        for (c = im1; c <= jy; c += 2) {
          tmp = _mm_loadu_pd(&b_a[c - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(alpha1));
          _mm_storeu_pd(&b_a[c - 1], tmp);
        }

        for (c = jA; c <= rowright; c++) {
          b_a[c - 1] *= alpha1;
        }

        alpha1 = xnorm;
      }
    }

    tau[b_i] = b;
    b_a[(ip1 + (b_i << 2)) - 1] = 1.0;
    knt = 3 - b_i;
    jy = ip1 + im1n;
    c = in + 1;
    xnorm = tau[b_i];
    if (xnorm != 0.0) {
      jA = knt;
      jA--;
      i = (jy + jA) - 1;
      while ((knt > 0) && (b_a[i] == 0.0)) {
        knt--;
        i--;
      }

      im1 = 4;
      exitg2 = false;
      while ((!exitg2) && (im1 > 0)) {
        jA = im1;
        i = (c + jA) - 1;
        jA = knt;
        jA = (jA - 1) << 2;
        rowright = i + jA;
        do {
          exitg1 = 0;
          if (i <= rowright) {
            if (b_a[i - 1] != 0.0) {
              exitg1 = 1;
            } else {
              i += 4;
            }
          } else {
            im1--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      knt = 0;
      im1 = 0;
    }

    if (knt > 0) {
      for (nm1 = 0; nm1 < 16; nm1++) {
        b = b_a[nm1];
        x[nm1] = b;
        A[nm1] = b;
      }

      if (im1 != 0) {
        mm1 = im1 - 1;
        nm1 = knt;
        jA = mm1;
        rowright = jA + 1;
        rowright--;
        for (i = 0; i <= rowright; i++) {
          iy = i;
          work[iy] = 0.0;
        }

        i = jy - 1;
        jA = (nm1 - 1) << 2;
        nm1 = c + jA;
        for (rowright = c; rowright <= nm1; rowright += 4) {
          c_0 = x[i];
          iy = 0;
          e = rowright + mm1;
          for (ia = rowright; ia <= e; ia++) {
            work[iy] += A[ia - 1] * c_0;
            iy++;
          }

          i++;
        }
      }

      b = -xnorm;
      if (!(b == 0.0)) {
        jA = c;
        rowright = knt - 1;
        for (knt = 0; knt <= rowright; knt++) {
          xnorm = b_a[jy - 1];
          if (xnorm != 0.0) {
            xnorm *= b;
            i = 0;
            c = jA;
            nm1 = (im1 + jA) - 1;
            for (mm1 = c; mm1 <= nm1; mm1++) {
              b_a[mm1 - 1] += work[i] * xnorm;
              i++;
            }
          }

          jy++;
          jA += 4;
        }
      }
    }

    knt = 3 - b_i;
    im1 = 3 - b_i;
    c = ip1 + im1n;
    jA = ip1 + in;
    xnorm = tau[b_i];
    if (xnorm != 0.0) {
      im1n = knt;
      im1n--;
      i = (c + im1n) - 1;
      while ((knt > 0) && (b_a[i] == 0.0)) {
        knt--;
        i--;
      }

      exitg2 = false;
      while ((!exitg2) && (im1 > 0)) {
        im1n = im1;
        im1n = (im1n - 1) << 2;
        in = jA + im1n;
        im1n = knt;
        im1n = (in + im1n) - 1;
        do {
          exitg1 = 0;
          if (in <= im1n) {
            if (b_a[in - 1] != 0.0) {
              exitg1 = 1;
            } else {
              in++;
            }
          } else {
            im1--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      knt = 0;
      im1 = 0;
    }

    if (knt > 0) {
      for (nm1 = 0; nm1 < 16; nm1++) {
        b = b_a[nm1];
        x[nm1] = b;
        A[nm1] = b;
      }

      if (im1 != 0) {
        mm1 = knt;
        nm1 = im1 - 1;
        im1n = nm1;
        rowright = im1n + 1;
        rowright--;
        for (i = 0; i <= rowright; i++) {
          iy = i;
          work[iy] = 0.0;
        }

        iy = 0;
        im1n = nm1 << 2;
        nm1 = jA + im1n;
        for (rowright = jA; rowright <= nm1; rowright += 4) {
          i = c - 1;
          c_0 = 0.0;
          e = (rowright + mm1) - 1;
          for (ia = rowright; ia <= e; ia++) {
            y = A[ia - 1];
            b = x[i];
            b *= y;
            c_0 += b;
            i++;
          }

          work[iy] += c_0;
          iy++;
        }
      }

      memcpy(&x[0], &b_a[0], sizeof(real_T) << 4U);
      simulink_experiment_debug_xgerc(knt, im1, -xnorm, c, work, x, jA, b_a);
    }

    b_a[(ip1 + (b_i << 2)) - 1] = alpha1;
  }
}

static real_T simulink_experiment_deb_xnrm2_j(int32_T n, const real_T x[3])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    absxk = x[1];
    y = fabs(absxk);
  } else {
    scale = 3.3121686421112381E-170;
    absxk = x[1];
    absxk = fabs(absxk);
    if (absxk > 3.3121686421112381E-170) {
      y = 1.0;
      scale = absxk;
    } else {
      t = absxk / 3.3121686421112381E-170;
      y = t * t;
    }

    absxk = x[2];
    absxk = fabs(absxk);
    if (absxk > scale) {
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }

    y = scale * sqrt(y);
  }

  return y;
}

static void simulink_experiment_d_xzlarfg_j(int32_T n, real_T alpha1, real_T x[3],
  real_T *b_alpha1, real_T *tau)
{
  __m128d tmp;
  real_T xnorm;
  real_T y;
  int32_T a;
  int32_T c;
  int32_T knt;
  int32_T nm1;
  int32_T scalarLB;
  int32_T vectorUB;
  *b_alpha1 = alpha1;
  *tau = 0.0;
  if (n > 0) {
    nm1 = n - 2;
    xnorm = simulink_experiment_deb_xnrm2_j(nm1 + 1, x);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(*b_alpha1, xnorm);
      if (*b_alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      y = fabs(xnorm);
      if (y < 1.0020841800044864E-292) {
        knt = -1;
        do {
          knt++;
          c = nm1;
          c += 2;
          scalarLB = (((c - 1) / 2) << 1) + 2;
          vectorUB = scalarLB - 2;
          for (a = 2; a <= vectorUB; a += 2) {
            tmp = _mm_loadu_pd(&x[a - 1]);
            tmp = _mm_mul_pd(tmp, _mm_set1_pd(9.9792015476736E+291));
            _mm_storeu_pd(&x[a - 1], tmp);
          }

          for (a = scalarLB; a <= c; a++) {
            x[a - 1] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          *b_alpha1 *= 9.9792015476736E+291;
          y = fabs(xnorm);
        } while ((y < 1.0020841800044864E-292) && (knt + 1 < 20));

        xnorm = simulink_experiment_deb_xnrm2_j(nm1 + 1, x);
        xnorm = rt_hypotd_snf(*b_alpha1, xnorm);
        if (*b_alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        y = xnorm - *b_alpha1;
        *tau = y / xnorm;
        y = *b_alpha1 - xnorm;
        *b_alpha1 = 1.0 / y;
        c = nm1;
        c += 2;
        scalarLB = (((c - 1) / 2) << 1) + 2;
        vectorUB = scalarLB - 2;
        for (a = 2; a <= vectorUB; a += 2) {
          tmp = _mm_loadu_pd(&x[a - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(*b_alpha1));
          _mm_storeu_pd(&x[a - 1], tmp);
        }

        for (a = scalarLB; a <= c; a++) {
          x[a - 1] *= *b_alpha1;
        }

        for (a = 0; a <= knt; a++) {
          xnorm *= 1.0020841800044864E-292;
        }

        *b_alpha1 = xnorm;
      } else {
        y = xnorm - *b_alpha1;
        *tau = y / xnorm;
        y = *b_alpha1 - xnorm;
        *b_alpha1 = 1.0 / y;
        c = nm1;
        c += 2;
        scalarLB = (((c - 1) / 2) << 1) + 2;
        vectorUB = scalarLB - 2;
        for (a = 2; a <= vectorUB; a += 2) {
          tmp = _mm_loadu_pd(&x[a - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(*b_alpha1));
          _mm_storeu_pd(&x[a - 1], tmp);
        }

        for (a = scalarLB; a <= c; a++) {
          x[a - 1] *= *b_alpha1;
        }

        *b_alpha1 = xnorm;
      }
    }
  }
}

static void simulink_experiment_deb_xdlanv2(real_T a, real_T b, real_T c, real_T
  d, real_T *rt1r, real_T *rt1i, real_T *rt2r, real_T *rt2i, real_T *b_a, real_T
  *b_b, real_T *b_c, real_T *b_d, real_T *cs, real_T *sn)
{
  real_T bcmax;
  real_T bcmis;
  real_T p;
  real_T sab;
  real_T scale;
  real_T temp;
  real_T z;
  int32_T count;
  *b_d = d;
  *b_c = c;
  *b_b = b;
  *b_a = a;
  if (*b_c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b_b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    temp = *b_d;
    *b_d = *b_a;
    *b_a = temp;
    *b_b = -*b_c;
    *b_c = 0.0;
  } else if ((*b_a - *b_d == 0.0) && ((*b_b < 0.0) != (*b_c < 0.0))) {
    *cs = 1.0;
    *sn = 0.0;
  } else {
    temp = *b_a - *b_d;
    p = 0.5 * temp;
    scale = fabs(*b_b);
    z = fabs(*b_c);
    if ((scale >= z) || rtIsNaN(z)) {
      bcmax = scale;
    } else {
      bcmax = z;
    }

    z = 1.0;
    if (!(*b_b < 0.0)) {
      scale = 1.0;
    } else {
      scale = -1.0;
    }

    if (*b_c < 0.0) {
      z = -1.0;
    }

    bcmis = fabs(*b_b);
    sab = fabs(*b_c);
    if ((bcmis <= sab) || rtIsNaN(sab)) {
      sab = bcmis;
    }

    bcmis = sab * scale * z;
    scale = fabs(p);
    if ((!(scale >= bcmax)) && (!rtIsNaN(bcmax))) {
      scale = bcmax;
    }

    z = p / scale * p + bcmax / scale * bcmis;
    if (z >= 8.8817841970012523E-16) {
      scale = sqrt(scale);
      z = sqrt(z);
      z *= scale;
      if (!(p < 0.0)) {
        scale = z;
      } else {
        scale = -z;
      }

      z = p + scale;
      *b_a = *b_d + z;
      *b_d -= bcmax / z * bcmis;
      bcmax = rt_hypotd_snf(*b_c, z);
      *cs = z / bcmax;
      *sn = *b_c / bcmax;
      *b_b -= *b_c;
      *b_c = 0.0;
    } else {
      bcmis = *b_b + *b_c;
      scale = fabs(temp);
      z = fabs(bcmis);
      if ((!(scale >= z)) && (!rtIsNaN(z))) {
        scale = z;
      }

      count = 0;
      while ((scale >= 7.4428285367870146E+137) && (count <= 20)) {
        bcmis *= 1.3435752215134178E-138;
        temp *= 1.3435752215134178E-138;
        scale = fabs(temp);
        z = fabs(bcmis);
        if ((!(scale >= z)) && (!rtIsNaN(z))) {
          scale = z;
        }

        count++;
      }

      while ((scale <= 1.3435752215134178E-138) && (count <= 20)) {
        bcmis *= 7.4428285367870146E+137;
        temp *= 7.4428285367870146E+137;
        scale = fabs(temp);
        z = fabs(bcmis);
        if ((!(scale >= z)) && (!rtIsNaN(z))) {
          scale = z;
        }

        count++;
      }

      p = 0.5 * temp;
      bcmax = rt_hypotd_snf(bcmis, temp);
      scale = fabs(bcmis);
      *cs = (scale / bcmax + 1.0) * 0.5;
      *cs = sqrt(*cs);
      if (!(bcmis < 0.0)) {
        scale = 1.0;
      } else {
        scale = -1.0;
      }

      *sn = -(p / (bcmax * *cs)) * scale;
      temp = *b_a * *cs + *b_b * *sn;
      p = -*b_a * *sn + *b_b * *cs;
      scale = *b_c * *cs + *b_d * *sn;
      z = -*b_c * *sn + *b_d * *cs;
      *b_a = temp * *cs + scale * *sn;
      *b_b = p * *cs + z * *sn;
      *b_c = -temp * *sn + scale * *cs;
      *b_d = -p * *sn + z * *cs;
      temp = (*b_a + *b_d) * 0.5;
      *b_a = temp;
      *b_d = temp;
      if (*b_c != 0.0) {
        if (*b_b != 0.0) {
          if ((*b_b < 0.0) == (*b_c < 0.0)) {
            sab = fabs(*b_b);
            sab = sqrt(sab);
            bcmis = fabs(*b_c);
            bcmis = sqrt(bcmis);
            z = sab * bcmis;
            if (!(*b_c < 0.0)) {
              p = z;
            } else {
              p = -z;
            }

            scale = *b_b + *b_c;
            scale = fabs(scale);
            scale = sqrt(scale);
            bcmax = 1.0 / scale;
            *b_a = temp + p;
            *b_d = temp - p;
            *b_b -= *b_c;
            *b_c = 0.0;
            p = sab * bcmax;
            scale = bcmis * bcmax;
            temp = *cs * p - *sn * scale;
            *sn = *cs * scale + *sn * p;
            *cs = temp;
          }
        } else {
          *b_b = -*b_c;
          *b_c = 0.0;
          temp = *cs;
          *cs = -*sn;
          *sn = temp;
        }
      }
    }
  }

  *rt1r = *b_a;
  *rt2r = *b_d;
  if (*b_c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    scale = fabs(*b_b);
    scale = sqrt(scale);
    z = fabs(*b_c);
    z = sqrt(z);
    *rt1i = scale * z;
    *rt2i = -*rt1i;
  }
}

static void simulink_experiment_debu_xhseqr(const real_T h[16], const real_T z
  [16], real_T b_h[16], int32_T *info, real_T b_z[16])
{
  real_T b_v[3];
  real_T aa;
  real_T ab;
  real_T ba;
  real_T h12;
  real_T h21s;
  real_T h22;
  real_T htmp1;
  real_T rt2i;
  real_T tst;
  real_T y;
  real_T y_0;
  int32_T L;
  int32_T b_k;
  int32_T hoffset;
  int32_T i;
  int32_T ix;
  int32_T j;
  int32_T k;
  int32_T kdefl;
  int32_T nr;
  int32_T t;
  boolean_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T goto150;
  memcpy(&b_z[0], &z[0], sizeof(real_T) << 4U);
  memcpy(&b_h[0], &h[0], sizeof(real_T) << 4U);
  *info = 0;
  b_v[0] = 0.0;
  b_v[1] = 0.0;
  b_v[2] = 0.0;
  b_h[2] = 0.0;
  b_h[3] = 0.0;
  b_h[7] = 0.0;
  kdefl = 0;
  i = 3;
  exitg1 = false;
  while ((!exitg1) && (i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    ix = 0;
    exitg2 = false;
    while ((!exitg2) && (ix < 301)) {
      k = i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > L)) {
        aa = b_h[((k - 1) << 2) + k];
        htmp1 = fabs(aa);
        if (htmp1 <= 4.0083367200179456E-292) {
          exitg3 = true;
        } else {
          aa = b_h[(((k - 1) << 2) + k) - 1];
          htmp1 = fabs(aa);
          aa = b_h[(k << 2) + k];
          ab = fabs(aa);
          tst = htmp1 + ab;
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              aa = b_h[(((k - 2) << 2) + k) - 1];
              tst = fabs(aa);
            }

            if (k + 2 <= 4) {
              aa = b_h[((k << 2) + k) + 1];
              htmp1 = fabs(aa);
              tst += htmp1;
            }
          }

          aa = b_h[((k - 1) << 2) + k];
          htmp1 = fabs(aa);
          if (htmp1 <= 2.2204460492503131E-16 * tst) {
            htmp1 = fabs(aa);
            aa = b_h[((k << 2) + k) - 1];
            tst = fabs(aa);
            if (htmp1 > tst) {
              ab = htmp1;
              ba = tst;
            } else {
              ab = tst;
              ba = htmp1;
            }

            aa = b_h[(k << 2) + k];
            htmp1 = fabs(aa);
            aa = b_h[(((k - 1) << 2) + k) - 1] - b_h[(k << 2) + k];
            tst = fabs(aa);
            if (htmp1 > tst) {
              aa = htmp1;
              htmp1 = tst;
            } else {
              aa = tst;
            }

            tst = aa + ab;
            aa = aa / tst * htmp1 * 2.2204460492503131E-16;
            if ((aa <= 4.0083367200179456E-292) || rtIsNaN(aa)) {
              aa = 4.0083367200179456E-292;
            }

            if (ab / tst * ba <= aa) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }

      L = k + 1;
      if (k + 1 > 1) {
        b_h[k + ((k - 1) << 2)] = 0.0;
      }

      if (k + 1 >= i) {
        goto150 = true;
        exitg2 = true;
      } else {
        kdefl++;
        t = kdefl / 20;
        t *= 20;
        t = kdefl - t;
        if (t == 0) {
          aa = b_h[((i - 1) << 2) + i];
          htmp1 = fabs(aa);
          aa = b_h[(((i - 2) << 2) + i) - 1];
          ab = fabs(aa);
          tst = htmp1 + ab;
          aa = b_h[(i << 2) + i] + 0.75 * tst;
          h12 = -0.4375 * tst;
          ba = tst;
          h22 = aa;
        } else {
          t = kdefl / 10;
          t *= 10;
          t = kdefl - t;
          if (t == 0) {
            aa = b_h[((k << 2) + k) + 1];
            htmp1 = fabs(aa);
            aa = b_h[(((k + 1) << 2) + k) + 2];
            ab = fabs(aa);
            tst = htmp1 + ab;
            aa = b_h[(k << 2) + k] + 0.75 * tst;
            h12 = -0.4375 * tst;
            ba = tst;
            h22 = aa;
          } else {
            aa = b_h[(((i - 1) << 2) + i) - 1];
            ba = b_h[((i - 1) << 2) + i];
            h12 = b_h[((i << 2) + i) - 1];
            h22 = b_h[(i << 2) + i];
          }
        }

        htmp1 = fabs(aa);
        ab = fabs(h12);
        tst = fabs(ba);
        h21s = fabs(h22);
        tst = ((htmp1 + ab) + tst) + h21s;
        if (tst == 0.0) {
          ba = 0.0;
          h22 = 0.0;
          h12 = 0.0;
          rt2i = 0.0;
        } else {
          aa /= tst;
          ba /= tst;
          h12 /= tst;
          h22 /= tst;
          htmp1 = (aa + h22) / 2.0;
          ab = (aa - htmp1) * (h22 - htmp1) - h12 * ba;
          aa = fabs(ab);
          aa = sqrt(aa);
          if (ab >= 0.0) {
            ba = htmp1 * tst;
            h12 = ba;
            h22 = aa * tst;
            rt2i = -h22;
          } else {
            ba = htmp1 + aa;
            h12 = htmp1 - aa;
            aa = ba - h22;
            htmp1 = fabs(aa);
            aa = h12 - h22;
            ab = fabs(aa);
            if (htmp1 <= ab) {
              ba *= tst;
              h12 = ba;
            } else {
              h12 *= tst;
              ba = h12;
            }

            h22 = 0.0;
            rt2i = 0.0;
          }
        }

        t = i - 2;
        exitg3 = false;
        while ((!exitg3) && (t + 1 >= k + 1)) {
          h21s = b_h[((t << 2) + t) + 1];
          aa = b_h[(t << 2) + t] - h12;
          htmp1 = fabs(aa);
          ab = fabs(rt2i);
          tst = fabs(h21s);
          tst += htmp1 + ab;
          h21s = b_h[((t << 2) + t) + 1] / tst;
          b_v[0] = ((b_h[(t << 2) + t] - h12) / tst * (b_h[(t << 2) + t] - ba) +
                    b_h[((t + 1) << 2) + t] * h21s) - rt2i / tst * h22;
          b_v[1] = (((b_h[(((t + 1) << 2) + t) + 1] + b_h[(t << 2) + t]) - ba) -
                    h12) * h21s;
          b_v[2] = b_h[(((t + 1) << 2) + t) + 2] * h21s;
          aa = b_v[0];
          htmp1 = fabs(aa);
          aa = b_v[1];
          ab = fabs(aa);
          aa = b_v[2];
          tst = fabs(aa);
          tst += htmp1 + ab;
          b_v[0] /= tst;
          b_v[1] /= tst;
          b_v[2] /= tst;
          if (t + 1 == k + 1) {
            exitg3 = true;
          } else {
            aa = b_h[t];
            htmp1 = fabs(aa);
            aa = b_v[1];
            ab = fabs(aa);
            aa = b_v[2];
            tst = fabs(aa);
            aa = b_v[0];
            h21s = fabs(aa);
            aa = b_h[0];
            y = fabs(aa);
            aa = b_h[(t << 2) + t];
            y_0 = fabs(aa);
            aa = b_h[(((t + 1) << 2) + t) + 1];
            aa = fabs(aa);
            if ((ab + tst) * htmp1 <= ((y + y_0) + aa) * (2.2204460492503131E-16
                 * h21s)) {
              exitg3 = true;
            } else {
              t--;
            }
          }
        }

        k = i;
        for (b_k = t + 1; b_k <= k; b_k++) {
          nr = (i - b_k) + 2;
          if (nr >= 3) {
            nr = 3;
          }

          if (b_k > t + 1) {
            hoffset = (((b_k - 2) << 2) + b_k) - 2;
            for (j = 0; j < nr; j++) {
              b_v[j] = b_h[(j + hoffset) + 1];
            }
          }

          htmp1 = b_v[0];
          simulink_experiment_d_xzlarfg_j(nr, htmp1, b_v, &aa, &tst);
          b_v[0] = aa;
          if (b_k > t + 1) {
            b_h[(b_k + ((b_k - 2) << 2)) - 1] = b_v[0];
            b_h[b_k + ((b_k - 2) << 2)] = 0.0;
            if (b_k < i) {
              b_h[b_k + 1] = 0.0;
            }
          } else if (t + 1 > L) {
            b_h[b_k - 1] *= 1.0 - tst;
          }

          htmp1 = b_v[1];
          aa = tst * htmp1;
          if (nr == 3) {
            ba = b_v[2];
            h21s = tst * ba;
            for (j = b_k; j < 5; j++) {
              ab = (b_h[(((j - 1) << 2) + b_k) - 1] + b_h[((j - 1) << 2) + b_k] *
                    htmp1) + b_h[(((j - 1) << 2) + b_k) + 1] * ba;
              b_h[(b_k + ((j - 1) << 2)) - 1] -= ab * tst;
              b_h[b_k + ((j - 1) << 2)] -= ab * aa;
              b_h[(b_k + ((j - 1) << 2)) + 1] -= ab * h21s;
            }

            nr = b_k + 3;
            j = i + 1;
            if (nr <= j) {
              j = nr;
            }

            nr = (uint8_T)j - 1;
            for (hoffset = 0; hoffset <= nr; hoffset++) {
              j = hoffset;
              ab = (b_h[((b_k - 1) << 2) + j] + b_h[(b_k << 2) + j] * htmp1) +
                b_h[((b_k + 1) << 2) + j] * ba;
              b_h[j + ((b_k - 1) << 2)] -= ab * tst;
              b_h[j + (b_k << 2)] -= ab * aa;
              b_h[j + ((b_k + 1) << 2)] -= ab * h21s;
            }

            for (nr = 0; nr < 4; nr++) {
              j = nr;
              ab = (b_z[((b_k - 1) << 2) + j] + b_z[(b_k << 2) + j] * htmp1) +
                b_z[((b_k + 1) << 2) + j] * ba;
              b_z[j + ((b_k - 1) << 2)] -= ab * tst;
              b_z[j + (b_k << 2)] -= ab * aa;
              b_z[j + ((b_k + 1) << 2)] -= ab * h21s;
            }
          } else if (nr == 2) {
            for (j = b_k; j < 5; j++) {
              ab = b_h[(((j - 1) << 2) + b_k) - 1] + b_h[((j - 1) << 2) + b_k] *
                htmp1;
              b_h[(b_k + ((j - 1) << 2)) - 1] -= ab * tst;
              b_h[b_k + ((j - 1) << 2)] -= ab * aa;
            }

            nr = (uint8_T)(i + 1) - 1;
            for (hoffset = 0; hoffset <= nr; hoffset++) {
              j = hoffset;
              ab = b_h[((b_k - 1) << 2) + j] + b_h[(b_k << 2) + j] * htmp1;
              b_h[j + ((b_k - 1) << 2)] -= ab * tst;
              b_h[j + (b_k << 2)] -= ab * aa;
            }

            for (nr = 0; nr < 4; nr++) {
              j = nr;
              ab = b_z[((b_k - 1) << 2) + j] + b_z[(b_k << 2) + j] * htmp1;
              b_z[j + ((b_k - 1) << 2)] -= ab * tst;
              b_z[j + (b_k << 2)] -= ab * aa;
            }
          }
        }

        ix++;
      }
    }

    if (!goto150) {
      *info = i + 1;
      exitg1 = true;
    } else {
      if ((i + 1 != L) && (L == i)) {
        kdefl = i - 1;
        simulink_experiment_deb_xdlanv2(b_h[kdefl + (kdefl << 2)], b_h[kdefl +
          (i << 2)], b_h[i + (kdefl << 2)], b_h[i + (i << 2)], &aa, &ab, &ba,
          &h21s, &h12, &h22, &rt2i, &y, &tst, &htmp1);
        b_h[kdefl + (kdefl << 2)] = h12;
        b_h[kdefl + (i << 2)] = h22;
        b_h[i + (kdefl << 2)] = rt2i;
        b_h[i + (i << 2)] = y;
        if (i + 1 < 4) {
          nr = 3 - i;
          ix = ((i + 1) << 2) + kdefl;
          t = ((i + 1) << 2) + i;
          b_k = (uint8_T)nr - 1;
          for (k = 0; k <= b_k; k++) {
            aa = tst * b_h[ix] + htmp1 * b_h[t];
            b_h[t] = tst * b_h[t] - htmp1 * b_h[ix];
            b_h[ix] = aa;
            t += 4;
            ix += 4;
          }
        }

        nr = kdefl;
        ix = (i - 1) << 2;
        t = (kdefl + 1) << 2;
        if (nr >= 1) {
          b_k = (uint8_T)nr - 1;
          for (k = 0; k <= b_k; k++) {
            aa = tst * b_h[ix] + htmp1 * b_h[t];
            b_h[t] = tst * b_h[t] - htmp1 * b_h[ix];
            b_h[ix] = aa;
            t++;
            ix++;
          }
        }

        ix = (i - 1) << 2;
        t = (kdefl + 1) << 2;
        aa = tst * b_z[ix] + htmp1 * b_z[t];
        b_z[t] = tst * b_z[t] - htmp1 * b_z[ix];
        b_z[ix] = aa;
        t++;
        ix++;
        aa = tst * b_z[ix] + htmp1 * b_z[t];
        b_z[t] = tst * b_z[t] - htmp1 * b_z[ix];
        b_z[ix] = aa;
        t++;
        ix++;
        aa = tst * b_z[ix] + htmp1 * b_z[t];
        b_z[t] = tst * b_z[t] - htmp1 * b_z[ix];
        b_z[ix] = aa;
        t++;
        ix++;
        aa = tst * b_z[ix] + htmp1 * b_z[t];
        b_z[t] = tst * b_z[t] - htmp1 * b_z[ix];
        b_z[ix] = aa;
      }

      kdefl = 0;
      i = L - 2;
    }
  }

  b_h[3] = 0.0;
}

static void simulink_experiment_debug_schur(const real_T A[16], real_T V[16],
  real_T T[16])
{
  __m128d tmp;
  real_T A_0[16];
  real_T x[16];
  real_T x_data[16];
  real_T work[4];
  real_T tau[3];
  real_T a;
  real_T b;
  real_T c_0;
  real_T x_0;
  int32_T coltop;
  int32_T e;
  int32_T exitg1;
  int32_T i;
  int32_T ia;
  int32_T iaii;
  int32_T iajm1;
  int32_T istart;
  int32_T itau;
  int32_T ix;
  int32_T iy;
  int32_T lastc;
  int32_T lastv;
  boolean_T b_0;
  boolean_T c;
  boolean_T exitg2;
  memcpy(&x[0], &A[0], sizeof(real_T) << 4U);
  c = true;
  for (istart = 0; istart < 16; istart++) {
    x_data[istart] = x[istart];
    if (c) {
      x_0 = x_data[istart];
      b_0 = rtIsInf(x_0);
      c = !b_0;
      b_0 = rtIsNaN(x_0);
      b_0 = !b_0;
      b_0 = (c && b_0);
      if (b_0) {
        c = true;
      } else {
        c = false;
      }
    } else {
      c = false;
    }
  }

  c = !c;
  if (c) {
    for (i = 0; i < 16; i++) {
      V[i] = (rtNaN);
    }

    istart = 2;
    for (lastc = 0; lastc < 3; lastc++) {
      lastv = lastc;
      for (iajm1 = istart; iajm1 < 5; iajm1++) {
        V[(iajm1 + (lastv << 2)) - 1] = 0.0;
      }

      istart++;
    }

    for (i = 0; i < 16; i++) {
      T[i] = (rtNaN);
    }
  } else {
    simulink_experiment_debu_xgehrd(A, x_data, tau);
    memcpy(&V[0], &x_data[0], sizeof(real_T) << 4U);
    for (lastv = 2; lastv >= 0; lastv--) {
      istart = lastv;
      istart = (istart + 1) << 2;
      ia = istart;
      coltop = lastv;
      istart = coltop;
      for (itau = 0; itau <= istart; itau++) {
        iajm1 = itau;
        i = ia + iajm1;
        V[i] = 0.0;
      }

      iajm1 = ia - 4;
      coltop = lastv + 3;
      for (itau = coltop; itau < 5; itau++) {
        istart = (ia + itau) - 1;
        i = (iajm1 + itau) - 1;
        V[istart] = V[i];
      }
    }

    V[1] = 0.0;
    V[2] = 0.0;
    V[3] = 0.0;
    V[0] = 1.0;
    itau = 2;
    work[0] = 0.0;
    work[1] = 0.0;
    work[2] = 0.0;
    work[3] = 0.0;
    for (iajm1 = 2; iajm1 >= 0; iajm1--) {
      istart = iajm1 + 6;
      i = iajm1;
      i <<= 2;
      iaii = istart + i;
      if (iajm1 + 1 < 3) {
        V[iaii - 1] = 1.0;
        istart = 2 - iajm1;
        lastv = istart + 1;
        lastc = 2 - iajm1;
        istart = iaii + 4;
        x_0 = tau[itau];
        if (x_0 != 0.0) {
          i = lastv;
          i--;
          i = (iaii + i) - 1;
          while ((lastv > 0) && (V[i] == 0.0)) {
            lastv--;
            i--;
          }

          exitg2 = false;
          while ((!exitg2) && (lastc > 0)) {
            i = lastc;
            i = (i - 1) << 2;
            coltop = istart + i;
            i = lastv;
            i = (coltop + i) - 1;
            do {
              exitg1 = 0;
              if (coltop <= i) {
                if (V[coltop - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  coltop++;
                }
              } else {
                lastc--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = 0;
          lastc = 0;
        }

        if (lastv > 0) {
          for (i = 0; i < 16; i++) {
            c_0 = V[i];
            x[i] = c_0;
            A_0[i] = c_0;
          }

          if (lastc != 0) {
            ia = lastc - 1;
            i = ia;
            coltop = i + 1;
            coltop--;
            for (i = 0; i <= coltop; i++) {
              iy = i;
              work[iy] = 0.0;
            }

            iy = 0;
            i = ia << 2;
            coltop = istart + i;
            for (i = istart; i <= coltop; i += 4) {
              ix = iaii - 1;
              c_0 = 0.0;
              e = (i + lastv) - 1;
              for (ia = i; ia <= e; ia++) {
                a = A_0[ia - 1];
                b = x[ix];
                b *= a;
                c_0 += b;
                ix++;
              }

              work[iy] += c_0;
              iy++;
            }
          }

          memcpy(&x[0], &V[0], sizeof(real_T) << 4U);
          simulink_experiment_debug_xgerc(lastv, lastc, -x_0, iaii, work, x,
            istart, V);
        }

        istart = 1 - iajm1;
        i = iaii + 1;
        a = -tau[itau];
        coltop = i + istart;
        lastv = ((((coltop - i) + 1) / 2) << 1) + i;
        lastc = lastv - 2;
        for (istart = i; istart <= lastc; istart += 2) {
          tmp = _mm_loadu_pd(&V[istart - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(a));
          _mm_storeu_pd(&V[istart - 1], tmp);
        }

        for (istart = lastv; istart <= coltop; istart++) {
          V[istart - 1] *= a;
        }
      }

      V[iaii - 1] = 1.0 - tau[itau];
      coltop = iajm1;
      istart = coltop - 1;
      for (lastc = 0; lastc <= istart; lastc++) {
        lastv = lastc;
        i = (iaii - lastv) - 2;
        V[i] = 0.0;
      }

      itau--;
    }

    memcpy(&x[0], &V[0], sizeof(real_T) << 4U);
    simulink_experiment_debu_xhseqr(x_data, x, T, &istart, V);
  }
}

static void simulink_experiment_debug__expm(real_T A[16], real_T F[16])
{
  __m128d tmp;
  __m128d tmp_0;
  real_T a[16];
  real_T b[16];
  real_T scaledT[16];
  real_T x[16];
  real_T x_0[16];
  real_T x_data[16];
  real_T alpha;
  real_T d10;
  real_T d4;
  real_T d6;
  int32_T blockFormat[3];
  int32_T b_i;
  int32_T b_k;
  int32_T b_k_0;
  int32_T eint;
  int32_T exitg1;
  int32_T i;
  boolean_T b_0;
  boolean_T c;
  boolean_T exitg2;
  boolean_T guard1;
  boolean_T guard2;
  boolean_T guard3;
  boolean_T guard4;
  boolean_T recomputeDiags;
  memcpy(&x[0], &A[0], sizeof(real_T) << 4U);
  c = true;
  for (b_k_0 = 0; b_k_0 < 16; b_k_0++) {
    x_data[b_k_0] = x[b_k_0];
    if (c) {
      alpha = x_data[b_k_0];
      b_0 = rtIsInf(alpha);
      c = !b_0;
      b_0 = rtIsNaN(alpha);
      b_0 = !b_0;
      b_0 = (c && b_0);
      if (b_0) {
        c = true;
      } else {
        c = false;
      }
    } else {
      c = false;
    }
  }

  c = !c;
  if (c) {
    for (b_i = 0; b_i < 16; b_i++) {
      F[b_i] = (rtNaN);
    }
  } else {
    c = true;
    b_k_0 = 1;
    exitg2 = false;
    while ((!exitg2) && (b_k_0 - 1 < 4)) {
      b_k = b_k_0;
      b_i = 1;
      do {
        exitg1 = 0;
        if (b_i - 1 < 4) {
          if ((b_i != b_k) && (!(A[(((b_k - 1) << 2) + b_i) - 1] == 0.0))) {
            c = false;
            exitg1 = 1;
          } else {
            b_i++;
          }
        } else {
          b_k_0++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }

    if (c) {
      for (b_i = 0; b_i < 16; b_i++) {
        F[b_i] = 0.0;
      }

      alpha = A[0];
      alpha = exp(alpha);
      F[0] = alpha;
      alpha = A[5];
      alpha = exp(alpha);
      F[5] = alpha;
      alpha = A[10];
      alpha = exp(alpha);
      F[10] = alpha;
      alpha = A[15];
      alpha = exp(alpha);
      F[15] = alpha;
    } else {
      memcpy(&x_data[0], &A[0], sizeof(real_T) << 4U);
      c = true;
      b_k_0 = 0;
      exitg2 = false;
      while ((!exitg2) && (b_k_0 < 4)) {
        d10 = (real_T)b_k_0 + 1.0;
        b_k = (int32_T)d10 - 1;
        b_i = 0;
        do {
          exitg1 = 0;
          if (b_i <= b_k) {
            d6 = (real_T)b_i + 1.0;
            if (!(x_data[((((int32_T)d10 - 1) << 2) + (int32_T)d6) - 1] ==
                  x_data[((((int32_T)d6 - 1) << 2) + (int32_T)d10) - 1])) {
              c = false;
              exitg1 = 1;
            } else {
              b_i++;
            }
          } else {
            b_k_0++;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }

      if (c) {
        simulink_experiment_debug_schur(A, x_data, a);
        for (b_k_0 = 0; b_k_0 < 4; b_k_0++) {
          b_k = b_k_0;
          d10 = a[(b_k << 2) + b_k];
          d10 = exp(d10);
          F[b_k << 2] = x_data[b_k << 2] * d10;
          F[(b_k << 2) + 1] = x_data[(b_k << 2) + 1] * d10;
          F[(b_k << 2) + 2] = x_data[(b_k << 2) + 2] * d10;
          F[(b_k << 2) + 3] = x_data[(b_k << 2) + 3] * d10;
        }

        for (b_i = 0; b_i < 4; b_i++) {
          for (i = 0; i < 4; i++) {
            a[b_i + (i << 2)] = 0.0;
            d10 = a[(i << 2) + b_i];
            d10 += F[b_i] * x_data[i];
            a[b_i + (i << 2)] = d10;
            d10 = a[(i << 2) + b_i];
            d10 += F[b_i + 4] * x_data[i + 4];
            a[b_i + (i << 2)] = d10;
            d10 = a[(i << 2) + b_i];
            d10 += F[b_i + 8] * x_data[i + 8];
            a[b_i + (i << 2)] = d10;
            d10 = a[(i << 2) + b_i];
            d10 += F[b_i + 12] * x_data[i + 12];
            a[b_i + (i << 2)] = d10;
          }
        }

        memcpy(&F[0], &a[0], sizeof(real_T) << 4U);
        for (b_i = 0; b_i < 4; b_i++) {
          a[b_i << 2] = (F[b_i << 2] + F[b_i]) / 2.0;
          a[(b_i << 2) + 1] = (F[(b_i << 2) + 1] + F[b_i + 4]) / 2.0;
          a[(b_i << 2) + 2] = (F[(b_i << 2) + 2] + F[b_i + 8]) / 2.0;
          a[(b_i << 2) + 3] = (F[(b_i << 2) + 3] + F[b_i + 12]) / 2.0;
        }

        memcpy(&F[0], &a[0], sizeof(real_T) << 4U);
      } else {
        recomputeDiags = true;
        b_k = 3;
        while (recomputeDiags && (b_k <= 4)) {
          b_i = b_k;
          while (recomputeDiags && (b_i <= 4)) {
            recomputeDiags = (x_data[(((b_k - 3) << 2) + b_i) - 1] == 0.0);
            b_i++;
          }

          b_k++;
        }

        if (recomputeDiags) {
          b_k = 1;
          exitg2 = false;
          while ((!exitg2) && (b_k - 1 < 3)) {
            b_k_0 = b_k;
            if (A[((b_k_0 - 1) << 2) + b_k_0] != 0.0) {
              if ((b_k_0 != 3) && (A[((b_k_0 << 2) + b_k_0) + 1] != 0.0)) {
                recomputeDiags = false;
                exitg2 = true;
              } else if (A[(((b_k_0 - 1) << 2) + b_k_0) - 1] != A[(b_k_0 << 2) +
                         b_k_0]) {
                recomputeDiags = false;
                exitg2 = true;
              } else {
                alpha = A[((b_k_0 - 1) << 2) + b_k_0];
                if (rtIsNaN(alpha)) {
                  alpha = (rtNaN);
                } else if (alpha < 0.0) {
                  alpha = -1.0;
                } else {
                  alpha = (alpha > 0.0);
                }

                d10 = A[((b_k_0 << 2) + b_k_0) - 1];
                if (rtIsNaN(d10)) {
                  d10 = (rtNaN);
                } else if (d10 < 0.0) {
                  d10 = -1.0;
                } else {
                  d10 = (d10 > 0.0);
                }

                if (alpha * d10 != -1.0) {
                  recomputeDiags = false;
                  exitg2 = true;
                } else {
                  b_k++;
                }
              }
            } else {
              b_k++;
            }
          }
        }

        d10 = 0.0;
        memcpy(&a[0], &A[0], sizeof(real_T) << 4U);
        for (b_i = 0; b_i < 4; b_i++) {
          for (i = 0; i <= 2; i += 2) {
            _mm_storeu_pd(&b[i + (b_i << 2)], _mm_set1_pd(0.0));
            tmp = _mm_loadu_pd(&a[i]);
            tmp = _mm_mul_pd(_mm_set1_pd(x_data[b_i << 2]), tmp);
            tmp_0 = _mm_loadu_pd(&b[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&b[i + (b_i << 2)], tmp);
            tmp = _mm_loadu_pd(&a[i + 4]);
            tmp = _mm_mul_pd(_mm_set1_pd(x_data[(b_i << 2) + 1]), tmp);
            tmp_0 = _mm_loadu_pd(&b[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&b[i + (b_i << 2)], tmp);
            tmp = _mm_loadu_pd(&a[i + 8]);
            tmp = _mm_mul_pd(_mm_set1_pd(x_data[(b_i << 2) + 2]), tmp);
            tmp_0 = _mm_loadu_pd(&b[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&b[i + (b_i << 2)], tmp);
            tmp = _mm_loadu_pd(&a[i + 12]);
            tmp = _mm_mul_pd(_mm_set1_pd(x_data[(b_i << 2) + 3]), tmp);
            tmp_0 = _mm_loadu_pd(&b[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&b[i + (b_i << 2)], tmp);
          }
        }

        memcpy(&x_data[0], &b[0], sizeof(real_T) << 4U);
        for (b_i = 0; b_i < 4; b_i++) {
          for (i = 0; i <= 2; i += 2) {
            _mm_storeu_pd(&a[i + (b_i << 2)], _mm_set1_pd(0.0));
            tmp = _mm_loadu_pd(&b[i]);
            tmp = _mm_mul_pd(_mm_set1_pd(x_data[b_i << 2]), tmp);
            tmp_0 = _mm_loadu_pd(&a[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&a[i + (b_i << 2)], tmp);
            tmp = _mm_loadu_pd(&b[i + 4]);
            tmp = _mm_mul_pd(_mm_set1_pd(x_data[(b_i << 2) + 1]), tmp);
            tmp_0 = _mm_loadu_pd(&a[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&a[i + (b_i << 2)], tmp);
            tmp = _mm_loadu_pd(&b[i + 8]);
            tmp = _mm_mul_pd(_mm_set1_pd(x_data[(b_i << 2) + 2]), tmp);
            tmp_0 = _mm_loadu_pd(&a[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&a[i + (b_i << 2)], tmp);
            tmp = _mm_loadu_pd(&b[i + 12]);
            tmp = _mm_mul_pd(_mm_set1_pd(x_data[(b_i << 2) + 3]), tmp);
            tmp_0 = _mm_loadu_pd(&a[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&a[i + (b_i << 2)], tmp);
          }
        }

        for (b_i = 0; b_i < 4; b_i++) {
          for (i = 0; i <= 2; i += 2) {
            _mm_storeu_pd(&x_data[i + (b_i << 2)], _mm_set1_pd(0.0));
            tmp = _mm_loadu_pd(&a[i]);
            tmp = _mm_mul_pd(_mm_set1_pd(b[b_i << 2]), tmp);
            tmp_0 = _mm_loadu_pd(&x_data[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&x_data[i + (b_i << 2)], tmp);
            tmp = _mm_loadu_pd(&a[i + 4]);
            tmp = _mm_mul_pd(_mm_set1_pd(b[(b_i << 2) + 1]), tmp);
            tmp_0 = _mm_loadu_pd(&x_data[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&x_data[i + (b_i << 2)], tmp);
            tmp = _mm_loadu_pd(&a[i + 8]);
            tmp = _mm_mul_pd(_mm_set1_pd(b[(b_i << 2) + 2]), tmp);
            tmp_0 = _mm_loadu_pd(&x_data[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&x_data[i + (b_i << 2)], tmp);
            tmp = _mm_loadu_pd(&a[i + 12]);
            tmp = _mm_mul_pd(_mm_set1_pd(b[(b_i << 2) + 3]), tmp);
            tmp_0 = _mm_loadu_pd(&x_data[(b_i << 2) + i]);
            tmp = _mm_add_pd(tmp, tmp_0);
            _mm_storeu_pd(&x_data[i + (b_i << 2)], tmp);
          }
        }

        alpha = simulink_experiment_debug__norm(a);
        d4 = rt_powd_snf(alpha, 0.25);
        alpha = simulink_experiment_debug__norm(x_data);
        d6 = rt_powd_snf(alpha, 0.16666666666666666);
        if ((!(d4 >= d6)) && (!rtIsNaN(d6))) {
          d4 = d6;
        }

        guard1 = false;
        guard2 = false;
        guard3 = false;
        guard4 = false;
        if (d4 <= 0.01495585217958292) {
          for (b_k_0 = 0; b_k_0 < 16; b_k_0++) {
            x_0[b_k_0] = x[b_k_0];
            alpha = x_0[b_k_0];
            alpha = fabs(alpha);
            scaledT[b_k_0] = alpha;
            scaledT[b_k_0] *= 0.19285012468241128;
          }

          simulink_experiment_debu_mpower(scaledT, 7.0, x_0);
          alpha = simulink_experiment_debug__norm(x_0);
          alpha /= simulink_experiment_debug__norm(x);
          alpha = simulink_experiment_debug__log2(2.0 * alpha /
            2.2204460492503131E-16) / 6.0;
          alpha = ceil(alpha);
          if (!(alpha >= 0.0)) {
            alpha = 0.0;
          }

          if (alpha == 0.0) {
            b_k_0 = 3;
          } else {
            guard4 = true;
          }
        } else {
          guard4 = true;
        }

        if (guard4) {
          if (d4 <= 0.253939833006323) {
            for (b_k_0 = 0; b_k_0 < 16; b_k_0++) {
              alpha = A[b_k_0];
              x_0[b_k_0] = alpha;
              x[b_k_0] = alpha;
              alpha = x_0[b_k_0];
              alpha = fabs(alpha);
              scaledT[b_k_0] = alpha;
              scaledT[b_k_0] *= 0.12321872304378752;
            }

            simulink_experiment_debu_mpower(scaledT, 11.0, x_0);
            alpha = simulink_experiment_debug__norm(x_0);
            alpha /= simulink_experiment_debug__norm(x);
            alpha = simulink_experiment_debug__log2(2.0 * alpha /
              2.2204460492503131E-16) / 10.0;
            alpha = ceil(alpha);
            if (!(alpha >= 0.0)) {
              alpha = 0.0;
            }

            if (alpha == 0.0) {
              b_k_0 = 5;
            } else {
              guard3 = true;
            }
          } else {
            guard3 = true;
          }
        }

        if (guard3) {
          simulink_experiment_debu_mpower(a, 2.0, x_0);
          alpha = simulink_experiment_debug__norm(x_0);
          d4 = rt_powd_snf(alpha, 0.125);
          if ((!(d6 >= d4)) && (!rtIsNaN(d4))) {
            d6 = d4;
          }

          if (d6 <= 0.95041789961629319) {
            for (b_k_0 = 0; b_k_0 < 16; b_k_0++) {
              alpha = A[b_k_0];
              x_0[b_k_0] = alpha;
              x[b_k_0] = alpha;
              alpha = x_0[b_k_0];
              alpha = fabs(alpha);
              scaledT[b_k_0] = alpha;
              scaledT[b_k_0] *= 0.090475336558796943;
            }

            simulink_experiment_debu_mpower(scaledT, 15.0, x_0);
            alpha = simulink_experiment_debug__norm(x_0);
            alpha /= simulink_experiment_debug__norm(x);
            alpha = simulink_experiment_debug__log2(2.0 * alpha /
              2.2204460492503131E-16) / 14.0;
            alpha = ceil(alpha);
            if (!(alpha >= 0.0)) {
              alpha = 0.0;
            }

            if (alpha == 0.0) {
              b_k_0 = 7;
            } else {
              guard2 = true;
            }
          } else {
            guard2 = true;
          }
        }

        if (guard2) {
          if (d6 <= 2.097847961257068) {
            for (b_k_0 = 0; b_k_0 < 16; b_k_0++) {
              alpha = A[b_k_0];
              x_0[b_k_0] = alpha;
              x[b_k_0] = alpha;
              alpha = x_0[b_k_0];
              alpha = fabs(alpha);
              scaledT[b_k_0] = alpha;
              scaledT[b_k_0] *= 0.071467735648795785;
            }

            simulink_experiment_debu_mpower(scaledT, 19.0, x_0);
            alpha = simulink_experiment_debug__norm(x_0);
            alpha /= simulink_experiment_debug__norm(x);
            alpha = simulink_experiment_debug__log2(2.0 * alpha /
              2.2204460492503131E-16) / 18.0;
            alpha = ceil(alpha);
            if (!(alpha >= 0.0)) {
              alpha = 0.0;
            }

            if (alpha == 0.0) {
              b_k_0 = 9;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }
        }

        if (guard1) {
          for (b_i = 0; b_i < 4; b_i++) {
            for (i = 0; i <= 2; i += 2) {
              _mm_storeu_pd(&x[i + (b_i << 2)], _mm_set1_pd(0.0));
              tmp = _mm_loadu_pd(&a[i]);
              tmp = _mm_mul_pd(_mm_set1_pd(x_data[b_i << 2]), tmp);
              tmp_0 = _mm_loadu_pd(&x[(b_i << 2) + i]);
              tmp = _mm_add_pd(tmp, tmp_0);
              _mm_storeu_pd(&x[i + (b_i << 2)], tmp);
              tmp = _mm_loadu_pd(&a[i + 4]);
              tmp = _mm_mul_pd(_mm_set1_pd(x_data[(b_i << 2) + 1]), tmp);
              tmp_0 = _mm_loadu_pd(&x[(b_i << 2) + i]);
              tmp = _mm_add_pd(tmp, tmp_0);
              _mm_storeu_pd(&x[i + (b_i << 2)], tmp);
              tmp = _mm_loadu_pd(&a[i + 8]);
              tmp = _mm_mul_pd(_mm_set1_pd(x_data[(b_i << 2) + 2]), tmp);
              tmp_0 = _mm_loadu_pd(&x[(b_i << 2) + i]);
              tmp = _mm_add_pd(tmp, tmp_0);
              _mm_storeu_pd(&x[i + (b_i << 2)], tmp);
              tmp = _mm_loadu_pd(&a[i + 12]);
              tmp = _mm_mul_pd(_mm_set1_pd(x_data[(b_i << 2) + 3]), tmp);
              tmp_0 = _mm_loadu_pd(&x[(b_i << 2) + i]);
              tmp = _mm_add_pd(tmp, tmp_0);
              _mm_storeu_pd(&x[i + (b_i << 2)], tmp);
            }
          }

          alpha = simulink_experiment_debug__norm(x);
          d10 = rt_powd_snf(alpha, 0.1);
          if ((d4 >= d10) || rtIsNaN(d10)) {
            d10 = d4;
          }

          if ((d6 <= d10) || rtIsNaN(d10)) {
            d10 = d6;
          }

          alpha = simulink_experiment_debug__log2(d10 / 5.3719203511481517);
          alpha = ceil(alpha);
          if (alpha >= 0.0) {
            d10 = alpha;
          } else {
            d10 = 0.0;
          }

          d6 = rt_powd_snf(2.0, d10);
          for (b_k_0 = 0; b_k_0 < 16; b_k_0++) {
            alpha = A[b_k_0] / d6;
            x_0[b_k_0] = alpha;
            x[b_k_0] = alpha;
            alpha = x_0[b_k_0];
            alpha = fabs(alpha);
            scaledT[b_k_0] = alpha;
            scaledT[b_k_0] *= 0.05031554467093536;
          }

          simulink_experiment_debu_mpower(scaledT, 27.0, x_0);
          alpha = simulink_experiment_debug__norm(x_0);
          alpha /= simulink_experiment_debug__norm(x);
          alpha = simulink_experiment_debug__log2(2.0 * alpha /
            2.2204460492503131E-16) / 26.0;
          alpha = ceil(alpha);
          if (!(alpha >= 0.0)) {
            alpha = 0.0;
          }

          d10 += alpha;
          b_0 = rtIsInf(d10);
          if (b_0) {
            d6 = simulink_experiment_debug__norm(A) / 5.3719203511481517;
            b_0 = rtIsInf(d6);
            c = !b_0;
            b_0 = rtIsNaN(d6);
            b_0 = !b_0;
            b_0 = (c && b_0);
            if (b_0) {
              d6 = frexp(d6, &eint);
              d10 = eint;
            } else {
              d10 = 0.0;
            }

            if (d6 == 0.5) {
              d10--;
            }
          }

          b_k_0 = 13;
        }

        if (d10 != 0.0) {
          alpha = rt_powd_snf(2.0, d10);
          for (b_i = 0; b_i <= 14; b_i += 2) {
            tmp = _mm_loadu_pd(&A[b_i]);
            tmp = _mm_div_pd(tmp, _mm_set1_pd(alpha));
            _mm_storeu_pd(&A[b_i], tmp);
          }

          d6 = 2.0 * d10;
          alpha = rt_powd_snf(2.0, d6);
          for (b_i = 0; b_i <= 14; b_i += 2) {
            tmp = _mm_loadu_pd(&b[b_i]);
            tmp = _mm_div_pd(tmp, _mm_set1_pd(alpha));
            _mm_storeu_pd(&b[b_i], tmp);
          }

          d6 = 4.0 * d10;
          alpha = rt_powd_snf(2.0, d6);
          for (b_i = 0; b_i <= 14; b_i += 2) {
            tmp = _mm_loadu_pd(&a[b_i]);
            tmp = _mm_div_pd(tmp, _mm_set1_pd(alpha));
            _mm_storeu_pd(&a[b_i], tmp);
          }

          d6 = 6.0 * d10;
          alpha = rt_powd_snf(2.0, d6);
          for (b_i = 0; b_i <= 14; b_i += 2) {
            tmp = _mm_loadu_pd(&x_data[b_i]);
            tmp = _mm_div_pd(tmp, _mm_set1_pd(alpha));
            _mm_storeu_pd(&x_data[b_i], tmp);
          }
        }

        if (recomputeDiags) {
          blockFormat[0] = 0;
          blockFormat[1] = 0;
          blockFormat[2] = 0;
          b_k = 0;
          while (b_k + 1 < 3) {
            if (A[((b_k << 2) + b_k) + 1] != 0.0) {
              blockFormat[b_k] = 2;
              blockFormat[b_k + 1] = 0;
              b_k += 2;
            } else if ((A[((b_k << 2) + b_k) + 1] == 0.0) && (A[(((b_k + 1) << 2)
              + b_k) + 2] == 0.0)) {
              blockFormat[b_k] = 1;
              b_k++;
            } else {
              blockFormat[b_k] = 0;
              b_k++;
            }
          }

          if (A[11] != 0.0) {
            blockFormat[2] = 2;
          } else if ((blockFormat[1] == 0) || (blockFormat[1] == 1)) {
            blockFormat[2] = 1;
          }
        }

        simulink_expe_padeApproximation(A, b, a, x_data, b_k_0, F);
        if (recomputeDiags) {
          memcpy(&x[0], &F[0], sizeof(real_T) << 4U);
          simulink_exp_recomputeBlockDiag(A, x, blockFormat, F);
        }

        b_k = (int32_T)d10 - 1;
        for (b_k_0 = 0; b_k_0 <= b_k; b_k_0++) {
          memcpy(&x_data[0], &F[0], sizeof(real_T) << 4U);
          for (b_i = 0; b_i < 4; b_i++) {
            for (i = 0; i < 4; i++) {
              a[b_i + (i << 2)] = 0.0;
              d10 = a[(i << 2) + b_i];
              d10 += x_data[i << 2] * F[b_i];
              a[b_i + (i << 2)] = d10;
              d10 = a[(i << 2) + b_i];
              d10 += x_data[(i << 2) + 1] * F[b_i + 4];
              a[b_i + (i << 2)] = d10;
              d10 = a[(i << 2) + b_i];
              d10 += x_data[(i << 2) + 2] * F[b_i + 8];
              a[b_i + (i << 2)] = d10;
              d10 = a[(i << 2) + b_i];
              d10 += x_data[(i << 2) + 3] * F[b_i + 12];
              a[b_i + (i << 2)] = d10;
            }
          }

          memcpy(&F[0], &a[0], sizeof(real_T) << 4U);
          if (recomputeDiags) {
            for (b_i = 0; b_i <= 14; b_i += 2) {
              tmp = _mm_loadu_pd(&A[b_i]);
              tmp = _mm_mul_pd(_mm_set1_pd(2.0), tmp);
              _mm_storeu_pd(&A[b_i], tmp);
            }

            memcpy(&x[0], &F[0], sizeof(real_T) << 4U);
            simulink_exp_recomputeBlockDiag(A, x, blockFormat, F);
          }
        }
      }
    }
  }
}

/* Model output function for TID0 */
void simulink_experiment_debug_type1_output0(void) /* Sample time: [0.0s, 0.0s] */
{
  studentControllerInterface_si_T *obj;
  studentControllerInterface_si_T *obj_0;
  studentControllerInterface_si_T *obj_1;
  real_T tauVec[100];
  real_T A_pred[16];
  real_T P_mat[16];
  real_T P_new[16];
  real_T P_plus[16];
  real_T b_y[16];
  real_T K[8];
  real_T a[8];
  real_T b_y_0[8];
  real_T Bd[4];
  real_T K_new[4];
  real_T S[4];
  real_T a_0[4];
  real_T prevState[4];
  real_T stateDerivative[4];
  real_T yPred[2];
  real_T amp;
  real_T b_x;
  real_T b_x_0;
  real_T b_x_1;
  real_T b_x_2;
  real_T b_x_3;
  real_T b_x_4;
  real_T b_x_5;
  real_T b_x_6;
  real_T b_x_7;
  real_T c;
  real_T c3;
  real_T c_0;
  real_T c_1;
  real_T c_2;
  real_T c_3;
  real_T c_4;
  real_T c_5;
  real_T c_6;
  real_T deltaTime;
  real_T error_x1;
  real_T error_x3;
  real_T refAccel;
  real_T refPos;
  real_T refVel;
  real_T s3;
  real_T t_ratio;
  real_T t_square;
  real_T u0;
  real_T x;
  real_T x_0;
  int32_T b_k;
  int32_T r1;
  int32_T r2;
  boolean_T p;
  boolean_T p_0;
  boolean_T p_1;
  static const real_T tmp[16] = { 1000.0, 0.0, 0.0, 0.0, 0.0, 400.0, 0.0, 0.0,
    0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 2.0E-5 };

  __m128d tmp_0;
  __m128d tmp_1;
  real_T z_idx_0;
  real_T z_idx_1;
  static const real_T tmp_2[16] = { 0.011, 0.0, 0.0, 0.0, 0.0, 0.055, 0.0, 0.0,
    0.0, 0.0, 0.022, 0.0, 0.0, 0.0, 0.0, 0.22 };

  static const int8_T tmp_3[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    1 };

  static const int8_T tmp_4[8] = { 1, 0, 0, 0, 0, 1, 0, 0 };

  static const int8_T tmp_5[16] = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
    0 };

  boolean_T exitg1;
  boolean_T exitg2;

  {                                    /* Sample time: [0.0s, 0.0s] */
    rate_monotonic_scheduler();
  }

  /* S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_read_encoder
      (simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task, 1,
       &simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Buffer);
    if (result < 0) {
      simulink_experiment_debug_typ_B.HILReadEncoderTimebase = 0;
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    } else {
      simulink_experiment_debug_typ_B.HILReadEncoderTimebase =
        simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Buffer;
    }
  }

  /* S-Function (hil_read_analog_block): '<S1>/HIL Read Analog' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Analog (hil_read_analog_block) */
  {
    t_error result = hil_read_analog
      (simulink_experiment_debug_ty_DW.HILInitialize_Card,
       &simulink_experiment_debug_typ_P.HILReadAnalog_channels, 1,
       &simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    }

    simulink_experiment_debug_typ_B.HILReadAnalog =
      simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer;
  }

  /* Gain: '<S1>/BB01 Sensor  Gain (m//V)' */
  simulink_experiment_debug_typ_B.BB01SensorGainmV =
    simulink_experiment_debug_typ_P.BB01SensorGainmV_Gain *
    simulink_experiment_debug_typ_B.HILReadAnalog;

  /* Gain: '<S1>/Encoder Calibration  (rad//count)' */
  simulink_experiment_debug_typ_B.EncoderCalibrationradcount =
    simulink_experiment_debug_typ_P.EncoderCalibrationradcount_Gain *
    simulink_experiment_debug_typ_B.HILReadEncoderTimebase;

  /* Bias: '<S1>/Bias' */
  simulink_experiment_debug_typ_B.Bias =
    simulink_experiment_debug_typ_B.EncoderCalibrationradcount +
    simulink_experiment_debug_typ_P.Bias_Bias;

  /* Clock: '<Root>/Clock' */
  simulink_experiment_debug_typ_B.Clock =
    simulink_experiment_debug_ty_M->Timing.t[0];

  /* MATLABSystem: '<Root>/MATLAB System' */
  u0 = simulink_experiment_debug_typ_B.Clock;
  z_idx_0 = simulink_experiment_debug_typ_B.BB01SensorGainmV;
  z_idx_1 = simulink_experiment_debug_typ_B.Bias;
  obj = &simulink_experiment_debug_ty_DW.obj;

  /*  stepImpl Run one control step: estimate state, compute control */
  /*   Inputs: */
  /*     currentTime          - current simulation time (s) */
  /*     measuredBallPosition - measured ball position along pendulum (m) */
  /*     measuredAngle        - measured pendulum angle (rad) */
  /*   Outputs: */
  /*     V_servo              - computed servo voltage command */
  /*     estimatedPosition    - Kalman-estimated position (m) */
  /*     estimatedVelocity    - Kalman-estimated velocity (m/s) */
  /*     estimatedAngle       - Kalman-estimated pendulum angle (rad) */
  /*     estimatedAngularVelocity - estimated angular velocity (rad/s) */
  /*     cumulativeErrorOut   - current integral of position error */
  /*  Physical constants for pendulum-ball system */
  /*  ball radius (m) */
  /*  pendulum length (m) */
  /*  gravitational acceleration (m/s^2) */
  /*  gain mapping input voltage to torque */
  /*  motor time constant (s) */
  /*  Initialize constants for feedback linearization if not set */
  if (obj->const_1 == 0.0) {
    obj->const_1 = 0.41828772872251135;
    obj->const_2 = 0.0025453075675320605;
  }

  /*  Determine time step dt, handling first-call or zero increment */
  if (u0 == obj->previousTime) {
    deltaTime = obj->previousDeltaTime;
  } else {
    deltaTime = u0 - obj->previousTime;
  }

  /*  Get desired reference trajectory at current time */
  if (u0 < 5.0) {
    refPos = 0.0;
    refVel = 0.0;
    refAccel = 0.0;
  } else if (u0 < 61.85) {
    t_ratio = u0 - 5.0;
    t_square = t_ratio / 56.85;
    if (t_square < 0.5) {
      amp = t_square / 0.5 * 0.090000000000000011 + 0.05;
      x = 0.11423973285781065 * t_ratio;
      x = sin(x);
      x = 0.83775804095727813 * t_ratio - 0.2094395102393195 * x /
        0.11423973285781065;
      x = sin(x);
      error_x1 = 0.11423973285781065 * t_ratio;
      error_x1 = sin(error_x1);
      error_x1 = 0.83775804095727813 * t_ratio - 0.2094395102393195 * error_x1 /
        0.11423973285781065;
      error_x1 = cos(error_x1);
      s3 = 0.11423973285781065 * t_ratio;
      s3 = cos(s3);
      refVel = (0.83775804095727813 - 0.2094395102393195 * s3) * (amp * error_x1)
        + 0.00316622691292876 * x;
      x = 6.2831853071795862 * t_ratio / 55.0;
      x = cos(x);
      x = 0.83775804095727813 - 3.1415926535897931 * x / 15.0;
      t_square = x * x;
      x = 6.2831853071795862 * t_ratio / 55.0;
      x = sin(x);
      x = 11.0 * x / 6.0 - 12.566370614359172 * t_ratio / 15.0;
      x = cos(x);
      error_x1 = 6.2831853071795862 * t_ratio / 55.0;
      error_x1 = cos(error_x1);
      s3 = 6.2831853071795862 * t_ratio / 55.0;
      s3 = sin(s3);
      s3 = 11.0 * s3 / 6.0 - 12.566370614359172 * t_ratio / 15.0;
      s3 = sin(s3);
      x_0 = 6.2831853071795862 * t_ratio / 55.0;
      x_0 = sin(x_0);
      refAccel = 6.2831853071795862 * t_ratio / 55.0;
      refAccel = sin(refAccel);
      refAccel = 11.0 * refAccel / 6.0 - 12.566370614359172 * t_ratio / 15.0;
      refAccel = cos(refAccel);
      refAccel = ((0.83775804095727813 - 3.1415926535897931 * error_x1 / 15.0) *
                  (12.0 * x) / 1895.0 + (6.0 * t_ratio / 1895.0 + 0.05) * s3 *
                  t_square) + (6.0 * t_ratio / 1895.0 + 0.05) *
        (19.739208802178716 * x_0 * refAccel) / 825.0;
    } else {
      amp = 0.14;
      x = 0.11423973285781065 * t_ratio;
      x = sin(x);
      x = 0.83775804095727813 * t_ratio - 0.2094395102393195 * x /
        0.11423973285781065;
      x = cos(x);
      error_x1 = 0.11423973285781065 * t_ratio;
      error_x1 = cos(error_x1);
      refVel = (0.83775804095727813 - 0.2094395102393195 * error_x1) * (0.14 * x);
      x = 6.2831853071795862 * t_ratio / 55.0;
      x = cos(x);
      x = 0.83775804095727813 - 3.1415926535897931 * x / 15.0;
      t_square = x * x;
      x = 6.2831853071795862 * t_ratio / 55.0;
      x = sin(x);
      x = 11.0 * x / 6.0 - 12.566370614359172 * t_ratio / 15.0;
      x = sin(x);
      error_x1 = 6.2831853071795862 * t_ratio / 55.0;
      error_x1 = sin(error_x1);
      s3 = 6.2831853071795862 * t_ratio / 55.0;
      s3 = sin(s3);
      s3 = 11.0 * s3 / 6.0 - 12.566370614359172 * t_ratio / 15.0;
      s3 = cos(s3);
      refAccel = 7.0 * x * t_square / 50.0 + 69.0872308076255 * error_x1 * s3 /
        20625.0;
    }

    x = 0.11423973285781065 * t_ratio;
    x = sin(x);
    x = 0.83775804095727813 * t_ratio - 0.2094395102393195 * x /
      0.11423973285781065;
    x = sin(x);
    refPos = amp * x;
  } else if (u0 < 65.0) {
    refPos = 0.0;
    refVel = 0.0;
    refAccel = 0.0;
  } else if (u0 < 85.0) {
    t_square = u0 - 65.0;
    t_ratio = t_square / 20.0;
    if (t_ratio < 0.5) {
      t_ratio = 0.05;
    } else {
      t_ratio = 0.1;
    }

    x = 0.62831853071795862 * t_square;
    x = sin(x);
    if (x < 0.0) {
      x = -1.0;
    } else {
      x = (x > 0.0);
    }

    refPos = t_ratio * x;
    refVel = 0.0;
    refAccel = 0.0;
  } else {
    refPos = 0.0;
    refVel = 0.0;
    refAccel = 0.0;
  }

  /*  Define process and measurement noise covariance matrices */
  /*  Retrieve previous state estimate and control input */
  prevState[0] = obj->stateEstimate[0];
  prevState[1] = obj->stateEstimate[1];
  prevState[2] = obj->stateEstimate[2];
  prevState[3] = obj->stateEstimate[3];
  amp = obj->controlInput;

  /*  Precompute sin and cos of estimated angle for speed */
  c3 = prevState[2];
  c3 = cos(c3);
  s3 = prevState[2];
  s3 = sin(s3);

  /*  Continuous-time state derivative under nonlinear dynamics */
  x = prevState[3];
  t_square = x * x;
  t_ratio = c3 * c3;
  stateDerivative[0] = prevState[1];
  stateDerivative[1] = 0.41828772872251135 * s3 - (0.21275 - prevState[0]) *
    0.7142857142857143 * 0.0035634305945448845 * t_square * t_ratio;
  stateDerivative[2] = prevState[3];
  stateDerivative[3] = -prevState[3] / 0.025 + 60.0 * amp;

  /*  Euler integration to predict next state */
  t_ratio = stateDerivative[0];
  t_ratio *= deltaTime;
  t_ratio += prevState[0];
  stateDerivative[0] = t_ratio;
  t_ratio = stateDerivative[1];
  t_ratio *= deltaTime;
  t_ratio += prevState[1];
  stateDerivative[1] = t_ratio;
  t_ratio = stateDerivative[2];
  t_ratio *= deltaTime;
  t_ratio += prevState[2];
  stateDerivative[2] = t_ratio;
  t_ratio = stateDerivative[3];
  t_ratio *= deltaTime;
  t_ratio += prevState[3];
  stateDerivative[3] = t_ratio;

  /*  Linearized discretized A matrix for covariance prediction */
  t_square = x * x;
  t_ratio = c3 * c3;
  amp = x * x;
  s3 = x * x;
  x_0 = c3 * c3;
  x = 2.0 * prevState[2];
  x = sin(x);
  error_x1 = 2.0 * prevState[2];
  error_x1 = sin(error_x1);
  A_pred[0] = 1.0;
  A_pred[4] = deltaTime;
  A_pred[8] = 0.0;
  A_pred[12] = 0.0;
  A_pred[1] = 5.0 * deltaTime * 0.00064516 * t_square * t_ratio / 1.26735175;
  A_pred[5] = 1.0;
  A_pred[9] = ((0.0108077 * amp * x + 8.34831 * c3) - 0.0508 * prevState[0] * s3
               * error_x1) * (5.0 * deltaTime * 0.0254) / 2.5347035;
  A_pred[13] = -(5.0 * deltaTime * 0.00064516 * prevState[3] * x_0 * (0.4255 -
    2.0 * prevState[0])) / 1.26735175;
  A_pred[2] = 0.0;
  A_pred[6] = 0.0;
  A_pred[10] = 1.0;
  A_pred[14] = deltaTime;
  A_pred[3] = 0.0;
  A_pred[7] = 0.0;
  A_pred[11] = 0.0;
  A_pred[15] = 1.0 - deltaTime / 0.025;

  /*  A-priori covariance update */
  for (r2 = 0; r2 < 16; r2++) {
    P_mat[r2] = obj->estimateCovariance[r2];
  }

  for (r2 = 0; r2 < 4; r2++) {
    for (b_k = 0; b_k < 4; b_k++) {
      b_y[r2 + (b_k << 2)] = 0.0;
      t_ratio = b_y[(b_k << 2) + r2];
      t_ratio += P_mat[b_k << 2] * A_pred[r2];
      b_y[r2 + (b_k << 2)] = t_ratio;
      t_ratio = b_y[(b_k << 2) + r2];
      t_ratio += P_mat[(b_k << 2) + 1] * A_pred[r2 + 4];
      b_y[r2 + (b_k << 2)] = t_ratio;
      t_ratio = b_y[(b_k << 2) + r2];
      t_ratio += P_mat[(b_k << 2) + 2] * A_pred[r2 + 8];
      b_y[r2 + (b_k << 2)] = t_ratio;
      t_ratio = b_y[(b_k << 2) + r2];
      t_ratio += P_mat[(b_k << 2) + 3] * A_pred[r2 + 12];
      b_y[r2 + (b_k << 2)] = t_ratio;
    }

    for (b_k = 0; b_k < 4; b_k++) {
      P_plus[r2 + (b_k << 2)] = 0.0;
      t_ratio = P_plus[(b_k << 2) + r2];
      t_ratio += b_y[r2] * A_pred[b_k];
      P_plus[r2 + (b_k << 2)] = t_ratio;
      t_ratio = P_plus[(b_k << 2) + r2];
      t_ratio += b_y[r2 + 4] * A_pred[b_k + 4];
      P_plus[r2 + (b_k << 2)] = t_ratio;
      t_ratio = P_plus[(b_k << 2) + r2];
      t_ratio += b_y[r2 + 8] * A_pred[b_k + 8];
      P_plus[r2 + (b_k << 2)] = t_ratio;
      t_ratio = P_plus[(b_k << 2) + r2];
      t_ratio += b_y[r2 + 12] * A_pred[b_k + 12];
      P_plus[r2 + (b_k << 2)] = t_ratio;
    }
  }

  for (r2 = 0; r2 < 16; r2++) {
    b_y[r2] = tmp_2[r2];
    P_mat[r2] = tmp_3[r2];
  }

  for (r2 = 0; r2 < 4; r2++) {
    for (b_k = 0; b_k < 4; b_k++) {
      A_pred[r2 + (b_k << 2)] = 0.0;
      x = A_pred[(b_k << 2) + r2];
      x += b_y[r2] * P_mat[b_k];
      A_pred[r2 + (b_k << 2)] = x;
      x = A_pred[(b_k << 2) + r2];
      x += b_y[r2 + 4] * P_mat[b_k + 4];
      A_pred[r2 + (b_k << 2)] = x;
      x = A_pred[(b_k << 2) + r2];
      x += b_y[r2 + 8] * P_mat[b_k + 8];
      A_pred[r2 + (b_k << 2)] = x;
      x = A_pred[(b_k << 2) + r2];
      x += b_y[r2 + 12] * P_mat[b_k + 12];
      A_pred[r2 + (b_k << 2)] = x;
    }
  }

  for (r2 = 0; r2 <= 14; r2 += 2) {
    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_loadu_pd(&P_plus[r2]);
    tmp_1 = _mm_loadu_pd(&A_pred[r2]);
    tmp_0 = _mm_add_pd(tmp_0, tmp_1);

    /* MATLABSystem: '<Root>/MATLAB System' */
    _mm_storeu_pd(&P_plus[r2], tmp_0);
  }

  /* MATLABSystem: '<Root>/MATLAB System' */
  /*  Kalman filter update step */
  /*  measurement vector */
  for (r2 = 0; r2 < 8; r2++) {
    a[r2] = tmp_4[r2];
  }

  /*  predicted measurement */
  Bd[0] = 1.0;
  Bd[1] = 0.0;
  Bd[2] = 0.0;
  Bd[3] = 0.09;
  for (r2 = 0; r2 < 2; r2++) {
    yPred[r2] = 0.0;
    for (b_k = 0; b_k < 4; b_k++) {
      amp = yPred[r2];
      amp += a[(b_k << 1) + r2] * stateDerivative[b_k];
      b_y_0[r2 + (b_k << 1)] = 0.0;
      t_ratio = b_y_0[(b_k << 1) + r2];
      t_ratio += P_plus[b_k << 2] * a[r2];
      b_y_0[r2 + (b_k << 1)] = t_ratio;
      t_ratio = b_y_0[(b_k << 1) + r2];
      t_ratio += P_plus[(b_k << 2) + 1] * a[r2 + 2];
      b_y_0[r2 + (b_k << 1)] = t_ratio;
      t_ratio = b_y_0[(b_k << 1) + r2];
      t_ratio += P_plus[(b_k << 2) + 2] * a[r2 + 4];
      b_y_0[r2 + (b_k << 1)] = t_ratio;
      t_ratio = b_y_0[(b_k << 1) + r2];
      t_ratio += P_plus[(b_k << 2) + 3] * a[r2 + 6];
      b_y_0[r2 + (b_k << 1)] = t_ratio;
      yPred[r2] = amp;
    }

    for (b_k = 0; b_k < 2; b_k++) {
      S[r2 + (b_k << 1)] = 0.0;
      s3 = S[(b_k << 1) + r2];
      s3 += b_y_0[r2] * a[b_k];
      S[r2 + (b_k << 1)] = s3;
      s3 = S[(b_k << 1) + r2];
      s3 += b_y_0[r2 + 2] * a[b_k + 2];
      S[r2 + (b_k << 1)] = s3;
      s3 = S[(b_k << 1) + r2];
      s3 += b_y_0[r2 + 4] * a[b_k + 4];
      S[r2 + (b_k << 1)] = s3;
      s3 = S[(b_k << 1) + r2];
      s3 += b_y_0[r2 + 6] * a[b_k + 6];
      S[r2 + (b_k << 1)] = s3;
    }

    prevState[r2] = 0.0;
    prevState[r2] += Bd[r2];
    prevState[r2 + 2] = 0.0;
    prevState[r2 + 2] += Bd[r2 + 2];
  }

  S[0] += prevState[0];
  S[1] += prevState[1];
  S[2] += prevState[2];
  S[3] += prevState[3];
  memcpy(&b_y[0], &P_plus[0], sizeof(real_T) << 4U);
  for (r2 = 0; r2 < 4; r2++) {
    for (b_k = 0; b_k < 2; b_k++) {
      b_y_0[r2 + (b_k << 2)] = 0.0;
      b_y_0[r2 + (b_k << 2)] += b_y[r2] * a[b_k];
      b_y_0[r2 + (b_k << 2)] += b_y[r2 + 4] * a[b_k + 2];
      b_y_0[r2 + (b_k << 2)] += b_y[r2 + 8] * a[b_k + 4];
      b_y_0[r2 + (b_k << 2)] += b_y[r2 + 12] * a[b_k + 6];
    }
  }

  x = S[1];
  t_square = fabs(x);
  t_ratio = t_square;
  x = S[0];
  t_square = fabs(x);
  if (t_ratio > t_square) {
    r1 = 1;
    r2 = 0;
  } else {
    r1 = 0;
    r2 = 1;
  }

  t_square = S[r2] / S[r1];
  t_ratio = S[r2 + 2] - S[r1 + 2] * t_square;
  K[r1 << 2] = b_y_0[0] / S[r1];
  K[r2 << 2] = (b_y_0[4] - K[r1 << 2] * S[r1 + 2]) / t_ratio;
  K[r1 << 2] -= K[r2 << 2] * t_square;
  K[(r1 << 2) + 1] = b_y_0[1] / S[r1];
  K[(r2 << 2) + 1] = (b_y_0[5] - K[(r1 << 2) + 1] * S[r1 + 2]) / t_ratio;
  K[(r1 << 2) + 1] -= K[(r2 << 2) + 1] * t_square;
  K[(r1 << 2) + 2] = b_y_0[2] / S[r1];
  K[(r2 << 2) + 2] = (b_y_0[6] - K[(r1 << 2) + 2] * S[r1 + 2]) / t_ratio;
  K[(r1 << 2) + 2] -= K[(r2 << 2) + 2] * t_square;
  K[(r1 << 2) + 3] = b_y_0[3] / S[r1];
  K[(r2 << 2) + 3] = (b_y_0[7] - K[(r1 << 2) + 3] * S[r1 + 2]) / t_ratio;
  K[(r1 << 2) + 3] -= K[(r2 << 2) + 3] * t_square;

  /*  Kalman gain */
  t_ratio = z_idx_0;
  t_ratio -= yPred[0];
  z_idx_0 = t_ratio;
  t_ratio = z_idx_1;
  t_ratio -= yPred[1];
  z_idx_1 = t_ratio;

  /*  state correction */
  /*  Covariance correction (Joseph form for numerical stability) */
  for (r2 = 0; r2 < 4; r2++) {
    t_ratio = stateDerivative[r2];
    amp = K[r2] * z_idx_0;
    amp += K[r2 + 4] * z_idx_1;
    t_ratio += amp;
    for (b_k = 0; b_k < 4; b_k++) {
      A_pred[r2 + (b_k << 2)] = 0.0;
      x = A_pred[(b_k << 2) + r2];
      x += a[b_k << 1] * K[r2];
      A_pred[r2 + (b_k << 2)] = x;
      x = A_pred[(b_k << 2) + r2];
      x += a[(b_k << 1) + 1] * K[r2 + 4];
      A_pred[r2 + (b_k << 2)] = x;
    }

    stateDerivative[r2] = t_ratio;
  }

  for (r2 = 0; r2 < 16; r2++) {
    P_mat[r2] = 0.0;
  }

  P_mat[0] = 1.0;
  P_mat[5] = 1.0;
  P_mat[10] = 1.0;
  P_mat[15] = 1.0;
  for (r2 = 0; r2 <= 14; r2 += 2) {
    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_loadu_pd(&P_mat[r2]);
    tmp_1 = _mm_loadu_pd(&A_pred[r2]);
    tmp_0 = _mm_sub_pd(tmp_0, tmp_1);

    /* MATLABSystem: '<Root>/MATLAB System' */
    _mm_storeu_pd(&P_mat[r2], tmp_0);
  }

  /* MATLABSystem: '<Root>/MATLAB System' */
  for (r2 = 0; r2 < 4; r2++) {
    for (b_k = 0; b_k <= 2; b_k += 2) {
      _mm_storeu_pd(&A_pred[b_k + (r2 << 2)], _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&K[b_k]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(a[r2 << 1]), tmp_0);
      tmp_1 = _mm_loadu_pd(&A_pred[(r2 << 2) + b_k]);
      tmp_0 = _mm_add_pd(tmp_0, tmp_1);
      _mm_storeu_pd(&A_pred[b_k + (r2 << 2)], tmp_0);
      tmp_0 = _mm_loadu_pd(&K[b_k + 4]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(a[(r2 << 1) + 1]), tmp_0);
      tmp_1 = _mm_loadu_pd(&A_pred[(r2 << 2) + b_k]);
      tmp_0 = _mm_add_pd(tmp_0, tmp_1);
      _mm_storeu_pd(&A_pred[b_k + (r2 << 2)], tmp_0);
    }
  }

  for (r2 = 0; r2 < 16; r2++) {
    P_new[r2] = 0.0;
  }

  P_new[0] = 1.0;
  P_new[5] = 1.0;
  P_new[10] = 1.0;
  P_new[15] = 1.0;
  for (r2 = 0; r2 <= 14; r2 += 2) {
    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_loadu_pd(&P_new[r2]);
    tmp_1 = _mm_loadu_pd(&A_pred[r2]);
    tmp_0 = _mm_sub_pd(tmp_0, tmp_1);

    /* MATLABSystem: '<Root>/MATLAB System' */
    _mm_storeu_pd(&P_new[r2], tmp_0);
  }

  /* MATLABSystem: '<Root>/MATLAB System' */
  for (r2 = 0; r2 < 4; r2++) {
    for (b_k = 0; b_k <= 2; b_k += 2) {
      _mm_storeu_pd(&b_y[b_k + (r2 << 2)], _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&P_new[b_k]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(P_plus[r2 << 2]), tmp_0);
      tmp_1 = _mm_loadu_pd(&b_y[(r2 << 2) + b_k]);
      tmp_0 = _mm_add_pd(tmp_0, tmp_1);
      _mm_storeu_pd(&b_y[b_k + (r2 << 2)], tmp_0);
      tmp_0 = _mm_loadu_pd(&P_new[b_k + 4]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(P_plus[(r2 << 2) + 1]), tmp_0);
      tmp_1 = _mm_loadu_pd(&b_y[(r2 << 2) + b_k]);
      tmp_0 = _mm_add_pd(tmp_0, tmp_1);
      _mm_storeu_pd(&b_y[b_k + (r2 << 2)], tmp_0);
      tmp_0 = _mm_loadu_pd(&P_new[b_k + 8]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(P_plus[(r2 << 2) + 2]), tmp_0);
      tmp_1 = _mm_loadu_pd(&b_y[(r2 << 2) + b_k]);
      tmp_0 = _mm_add_pd(tmp_0, tmp_1);
      _mm_storeu_pd(&b_y[b_k + (r2 << 2)], tmp_0);
      tmp_0 = _mm_loadu_pd(&P_new[b_k + 12]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(P_plus[(r2 << 2) + 3]), tmp_0);
      tmp_1 = _mm_loadu_pd(&b_y[(r2 << 2) + b_k]);
      tmp_0 = _mm_add_pd(tmp_0, tmp_1);
      _mm_storeu_pd(&b_y[b_k + (r2 << 2)], tmp_0);
    }
  }

  for (r2 = 0; r2 < 4; r2++) {
    for (b_k = 0; b_k < 4; b_k++) {
      P_plus[r2 + (b_k << 2)] = 0.0;
      t_ratio = P_plus[(b_k << 2) + r2];
      t_ratio += b_y[r2] * P_mat[b_k];
      P_plus[r2 + (b_k << 2)] = t_ratio;
      t_ratio = P_plus[(b_k << 2) + r2];
      t_ratio += b_y[r2 + 4] * P_mat[b_k + 4];
      P_plus[r2 + (b_k << 2)] = t_ratio;
      t_ratio = P_plus[(b_k << 2) + r2];
      t_ratio += b_y[r2 + 8] * P_mat[b_k + 8];
      P_plus[r2 + (b_k << 2)] = t_ratio;
      t_ratio = P_plus[(b_k << 2) + r2];
      t_ratio += b_y[r2 + 12] * P_mat[b_k + 12];
      P_plus[r2 + (b_k << 2)] = t_ratio;
    }

    a[r2] = 0.0;
    a[r2] += K[r2];
    a[r2] += K[r2 + 4] * 0.0;
    a[r2 + 4] = 0.0;
    a[r2 + 4] += K[r2] * 0.0;
    a[r2 + 4] += K[r2 + 4];
  }

  K_new[0] = 1.0;
  K_new[1] = 0.0;
  K_new[2] = 0.0;
  K_new[3] = 0.09;
  for (r2 = 0; r2 < 2; r2++) {
    for (b_k = 0; b_k <= 2; b_k += 2) {
      _mm_storeu_pd(&b_y_0[b_k + (r2 << 2)], _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&a[b_k]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(K_new[r2 << 1]), tmp_0);
      tmp_1 = _mm_loadu_pd(&b_y_0[(r2 << 2) + b_k]);
      tmp_0 = _mm_add_pd(tmp_1, tmp_0);
      _mm_storeu_pd(&b_y_0[b_k + (r2 << 2)], tmp_0);
      tmp_0 = _mm_loadu_pd(&a[b_k + 4]);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(K_new[(r2 << 1) + 1]), tmp_0);
      tmp_1 = _mm_loadu_pd(&b_y_0[(r2 << 2) + b_k]);
      tmp_0 = _mm_add_pd(tmp_0, tmp_1);
      _mm_storeu_pd(&b_y_0[b_k + (r2 << 2)], tmp_0);
    }
  }

  for (r2 = 0; r2 < 4; r2++) {
    a[r2] = 0.0;
    a[r2] += b_y_0[r2];
    a[r2] += b_y_0[r2 + 4] * 0.0;
    a[r2 + 4] = 0.0;
    a[r2 + 4] += b_y_0[r2] * 0.0;
    a[r2 + 4] += b_y_0[r2 + 4];
    for (b_k = 0; b_k < 4; b_k++) {
      A_pred[r2 + (b_k << 2)] = 0.0;
      x = A_pred[(b_k << 2) + r2];
      x += a[r2] * K[b_k];
      A_pred[r2 + (b_k << 2)] = x;
      x = A_pred[(b_k << 2) + r2];
      x += a[r2 + 4] * K[b_k + 4];
      A_pred[r2 + (b_k << 2)] = x;
    }
  }

  for (r2 = 0; r2 <= 14; r2 += 2) {
    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_loadu_pd(&P_plus[r2]);
    tmp_1 = _mm_loadu_pd(&A_pred[r2]);
    tmp_0 = _mm_add_pd(tmp_0, tmp_1);

    /* MATLABSystem: '<Root>/MATLAB System' */
    _mm_storeu_pd(&P_plus[r2], tmp_0);
  }

  /* MATLABSystem: '<Root>/MATLAB System' */
  /*  Save updated estimates */
  obj->stateEstimate[0] = stateDerivative[0];
  obj->stateEstimate[1] = stateDerivative[1];
  obj->stateEstimate[2] = stateDerivative[2];
  obj->stateEstimate[3] = stateDerivative[3];
  for (r2 = 0; r2 < 16; r2++) {
    obj->estimateCovariance[r2] = P_plus[r2];
  }

  /*  Control calculation - choose between feedback linearization and LQR */
  if (obj->useFeedbackLinearization) {
    /*  Feedback linearization control */
    /*  Compute the error state using Lie derivatives */
    error_x1 = stateDerivative[0] - refPos;
    c3 = stateDerivative[1] - refVel;
    obj_0 = obj;
    z_idx_0 = stateDerivative[0];
    refVel = stateDerivative[2];
    refPos = stateDerivative[3];

    /*  Second Lie derivative */
    t_square = refPos * refPos;
    x = refVel;
    x = cos(x);
    t_ratio = x * x;
    refVel = sin(refVel);
    s3 = obj_0->const_1 * refVel - (0.21275 - z_idx_0) * obj_0->const_2 *
      t_square * t_ratio;
    error_x3 = s3 - refAccel;
    z_idx_1 = stateDerivative[1];
    refVel = stateDerivative[2];

    /*  Third Lie derivative */
    t_square = refPos * refPos;
    x = refVel;
    x = cos(x);
    t_ratio = x * x;
    amp = rt_powd_snf(refPos, 3.0);
    s3 = refPos * refPos;
    x = refVel;
    x = cos(x);
    x_0 = x * x;
    refAccel = refVel;
    refAccel = cos(refAccel);
    x = refVel;
    x = cos(x);
    refVel = sin(refVel);
    t_square = (2.0 * obj_0->const_2 * (0.21275 - z_idx_0) * amp * x * refVel +
                (obj_0->const_2 * z_idx_1 * t_square * t_ratio + obj_0->const_1 *
                 refPos * refAccel)) + 2.0 * obj_0->const_2 / 0.025 * (0.21275 -
      z_idx_0) * s3 * x_0;
    prevState[0] = error_x1;
    prevState[1] = c3;
    prevState[2] = error_x3;
    prevState[3] = t_square;

    /*  Apply feedback linearization control law */
    error_x1 = -23.4521 * prevState[0];
    error_x1 += -29.3025 * prevState[1];
    error_x1 += -17.2402 * prevState[2];
    error_x1 += -5.872 * prevState[3];
    refVel = stateDerivative[2];

    /*  Feedback linearization control law: u = (v - phi(x))/psi(x) */
    obj_1 = obj_0;

    /*  phi(x) part of the control law */
    t_square = rt_powd_snf(refPos, 3.0);
    t_ratio = refPos * refPos;
    x = refVel;
    x = cos(x);
    amp = x * x;
    s3 = refPos * refPos;
    x = refVel;
    x = cos(x);
    x_0 = x * x;
    c3 = refPos * refPos;
    x = refVel;
    x = cos(x);
    error_x3 = x * x;
    c = refPos * refPos;
    c_0 = rt_powd_snf(refPos, 3.0);
    x = refVel;
    x = cos(x);
    c_1 = x * x;
    x = refVel;
    x = sin(x);
    c_2 = x * x;
    c_3 = refPos * refPos;
    x = refVel;
    x = cos(x);
    c_4 = x * x;
    c_5 = refPos * refPos;
    x = refVel;
    x = cos(x);
    c_6 = x * x;
    refAccel = refVel;
    refAccel = cos(refAccel);
    x = refVel;
    x = sin(x);
    b_x = refVel;
    b_x = sin(b_x);
    b_x_0 = refVel;
    b_x_0 = cos(b_x_0);
    b_x_1 = refVel;
    b_x_1 = sin(b_x_1);
    b_x_2 = refVel;
    b_x_2 = sin(b_x_2);
    b_x_3 = refVel;
    b_x_3 = cos(b_x_3);
    b_x_4 = refVel;
    b_x_4 = sin(b_x_4);
    b_x_5 = refVel;
    b_x_5 = cos(b_x_5);
    b_x_6 = refVel;
    b_x_6 = cos(b_x_6);
    b_x_7 = refVel;
    b_x_7 = sin(b_x_7);
    s3 = ((((-2.0 * obj_1->const_2 * z_idx_1 * c * b_x_0 * b_x_1 -
             obj_1->const_1 * refPos * b_x_2) + 2.0 * obj_1->const_2 * (0.21275
             - z_idx_0) * c_0 * (c_1 - c_2)) - 4.0 * obj_1->const_2 / 0.025 *
           (0.21275 - z_idx_0) * c_3 * b_x_3 * b_x_4) * refPos + ((-2.0 *
            obj_1->const_2 * t_square * refAccel * x - 2.0 * obj_1->const_2 /
            0.025 * t_ratio * amp) * z_idx_1 + (obj_1->const_1 * b_x - (0.21275
             - z_idx_0) * obj_1->const_2 * c3 * error_x3) * (obj_1->const_2 * s3
            * x_0))) + (((2.0 * obj_1->const_2 * z_idx_1 * refPos * c_4 +
                          obj_1->const_1 * b_x_5) + 6.0 * obj_1->const_2 *
                         (0.21275 - z_idx_0) * c_5 * b_x_6 * b_x_7) + 4.0 *
                        obj_1->const_2 / 0.025 * (0.21275 - z_idx_0) * refPos *
                        c_6) * (-refPos / 0.025);

    /*  psi(x) part of the control law */
    x = refVel;
    x = cos(x);
    t_square = x * x;
    t_ratio = refPos * refPos;
    x = refVel;
    x = cos(x);
    amp = x * x;
    refAccel = refVel;
    refAccel = cos(refAccel);
    x = refVel;
    x = cos(x);
    refVel = sin(refVel);
    t_square = (((2.0 * obj_0->const_2 * z_idx_1 * refPos * t_square +
                  obj_0->const_1 * refAccel) + 6.0 * obj_0->const_2 * (0.21275 -
      z_idx_0) * t_ratio * x * refVel) + 4.0 * obj_0->const_2 / 0.025 * (0.21275
      - z_idx_0) * refPos * amp) * 60.0;
    t_square = (error_x1 - s3) / t_square;
  } else {
    /*  Original LQR control (from the first script) */
    /*  Nonlinear dynamics feedforward terms */
    x = stateDerivative[3];
    t_square = x * x;
    x = stateDerivative[2];
    x = cos(x);
    error_x1 = stateDerivative[2];
    error_x1 = sin(error_x1);
    t_ratio = -16.731509148900454 * stateDerivative[3] * x - 0.41828772872251141
      * t_square * error_x1;
    x = stateDerivative[2];
    x = cos(x);
    amp = 25.097263723350682 * x;

    /*  Compute discrete-time system matrices for LQR */
    for (r2 = 0; r2 < 16; r2++) {
      x = tmp_5[r2];
      x *= deltaTime;
      A_pred[r2] = x;
    }

    simulink_experiment_debug__expm(A_pred, P_plus);

    /*  Numerical integration fallback for rank-deficient A */
    tauVec[99] = deltaTime;
    tauVec[0] = 0.0;
    if (-deltaTime == 0.0) {
      t_square = deltaTime / 99.0;
      for (b_k = 0; b_k < 98; b_k++) {
        r1 = ((b_k + 2) << 1) - 101;
        tauVec[b_k + 1] = (real_T)r1 * t_square;
      }
    } else if ((deltaTime < 0.0) && (fabs(deltaTime) > 8.9884656743115785E+307))
    {
      s3 = deltaTime / 99.0;
      for (b_k = 0; b_k < 98; b_k++) {
        t_square = (real_T)b_k + 1.0;
        tauVec[(int32_T)t_square] = s3 * t_square;
      }
    } else {
      s3 = deltaTime / 99.0;
      for (b_k = 0; b_k < 98; b_k++) {
        t_square = (real_T)b_k + 1.0;
        tauVec[(int32_T)t_square] = t_square * s3;
      }
    }

    t_square = deltaTime / 99.0;
    Bd[0] = 0.0;
    Bd[1] = 0.0;
    Bd[2] = 0.0;
    Bd[3] = 0.0;
    for (b_k = 0; b_k < 100; b_k++) {
      s3 = (real_T)b_k + 1.0;
      s3 = tauVec[(int32_T)s3 - 1];
      for (r2 = 0; r2 < 16; r2++) {
        x = tmp_5[r2];
        x *= s3;
        A_pred[r2] = x;
      }

      simulink_experiment_debug__expm(A_pred, b_y);
      for (r2 = 0; r2 <= 2; r2 += 2) {
        tmp_0 = _mm_loadu_pd(&b_y[r2]);
        tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(0.0));
        tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
        tmp_1 = _mm_loadu_pd(&b_y[r2 + 4]);
        tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(0.0));
        tmp_0 = _mm_add_pd(tmp_1, tmp_0);
        tmp_1 = _mm_loadu_pd(&b_y[r2 + 8]);
        tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(0.0));
        tmp_0 = _mm_add_pd(tmp_1, tmp_0);
        tmp_1 = _mm_loadu_pd(&b_y[r2 + 12]);
        tmp_0 = _mm_add_pd(tmp_1, tmp_0);
        tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(t_square));
        tmp_1 = _mm_loadu_pd(&Bd[r2]);
        tmp_0 = _mm_add_pd(tmp_1, tmp_0);
        _mm_storeu_pd(&Bd[r2], tmp_0);
      }
    }

    /*  Discrete-time LQR via iterative Riccati loop */
    for (r2 = 0; r2 < 16; r2++) {
      P_mat[r2] = 0.0;
    }

    P_mat[0] = 1.0;
    prevState[0] = 0.0;
    P_mat[5] = 1.0;
    prevState[1] = 0.0;
    P_mat[10] = 1.0;
    prevState[2] = 0.0;
    P_mat[15] = 1.0;
    prevState[3] = 0.0;
    r1 = 0;
    exitg1 = false;
    while ((!exitg1) && (r1 < 100000)) {
      memcpy(&b_y[0], &P_mat[0], sizeof(real_T) << 4U);
      for (r2 = 0; r2 < 4; r2++) {
        for (b_k = 0; b_k < 4; b_k++) {
          A_pred[r2 + (b_k << 2)] = 0.0;
          x = A_pred[(b_k << 2) + r2];
          x += P_plus[r2 << 2] * b_y[b_k << 2];
          A_pred[r2 + (b_k << 2)] = x;
          x = A_pred[(b_k << 2) + r2];
          x += P_plus[(r2 << 2) + 1] * b_y[(b_k << 2) + 1];
          A_pred[r2 + (b_k << 2)] = x;
          x = A_pred[(b_k << 2) + r2];
          x += P_plus[(r2 << 2) + 2] * b_y[(b_k << 2) + 2];
          A_pred[r2 + (b_k << 2)] = x;
          x = A_pred[(b_k << 2) + r2];
          x += P_plus[(r2 << 2) + 3] * b_y[(b_k << 2) + 3];
          A_pred[r2 + (b_k << 2)] = x;
        }

        for (b_k = 0; b_k < 4; b_k++) {
          P_new[r2 + (b_k << 2)] = 0.0;
          z_idx_0 = P_new[(b_k << 2) + r2];
          z_idx_0 += P_plus[b_k << 2] * A_pred[r2];
          P_new[r2 + (b_k << 2)] = z_idx_0;
          z_idx_0 = P_new[(b_k << 2) + r2];
          z_idx_0 += P_plus[(b_k << 2) + 1] * A_pred[r2 + 4];
          P_new[r2 + (b_k << 2)] = z_idx_0;
          z_idx_0 = P_new[(b_k << 2) + r2];
          z_idx_0 += P_plus[(b_k << 2) + 2] * A_pred[r2 + 8];
          P_new[r2 + (b_k << 2)] = z_idx_0;
          z_idx_0 = P_new[(b_k << 2) + r2];
          z_idx_0 += P_plus[(b_k << 2) + 3] * A_pred[r2 + 12];
          P_new[r2 + (b_k << 2)] = z_idx_0;
        }
      }

      memcpy(&b_y[0], &P_mat[0], sizeof(real_T) << 4U);
      for (r2 = 0; r2 < 4; r2++) {
        S[r2] = 0.0;
        for (b_k = 0; b_k < 4; b_k++) {
          s3 = S[r2];
          A_pred[r2 + (b_k << 2)] = 0.0;
          x = A_pred[(b_k << 2) + r2];
          x += P_plus[r2 << 2] * b_y[b_k << 2];
          A_pred[r2 + (b_k << 2)] = x;
          x = A_pred[(b_k << 2) + r2];
          x += P_plus[(r2 << 2) + 1] * b_y[(b_k << 2) + 1];
          A_pred[r2 + (b_k << 2)] = x;
          x = A_pred[(b_k << 2) + r2];
          x += P_plus[(r2 << 2) + 2] * b_y[(b_k << 2) + 2];
          A_pred[r2 + (b_k << 2)] = x;
          x = A_pred[(b_k << 2) + r2];
          x += P_plus[(r2 << 2) + 3] * b_y[(b_k << 2) + 3];
          A_pred[r2 + (b_k << 2)] = x;
          s3 += A_pred[(b_k << 2) + r2] * Bd[b_k];
          S[r2] = s3;
        }
      }

      memcpy(&b_y[0], &P_mat[0], sizeof(real_T) << 4U);
      t_square = 0.0;
      for (r2 = 0; r2 < 4; r2++) {
        z_idx_0 = b_y[r2 << 2] * Bd[0];
        z_idx_0 += b_y[(r2 << 2) + 1] * Bd[1];
        z_idx_0 += b_y[(r2 << 2) + 2] * Bd[2];
        z_idx_0 += b_y[(r2 << 2) + 3] * Bd[3];
        t_square += z_idx_0 * Bd[r2];
      }

      x = t_square + 0.055;
      t_square = 1.0 / x;
      S[0] *= t_square;
      S[1] *= t_square;
      S[2] *= t_square;
      S[3] *= t_square;
      memcpy(&b_y[0], &P_mat[0], sizeof(real_T) << 4U);
      for (r2 = 0; r2 < 4; r2++) {
        z_idx_0 = b_y[r2 << 2] * Bd[0];
        z_idx_0 += b_y[(r2 << 2) + 1] * Bd[1];
        z_idx_0 += b_y[(r2 << 2) + 2] * Bd[2];
        z_idx_0 += b_y[(r2 << 2) + 3] * Bd[3];
        a_0[r2] = z_idx_0;
      }

      for (r2 = 0; r2 < 4; r2++) {
        z_idx_0 = P_plus[r2 << 2] * a_0[0];
        z_idx_0 += P_plus[(r2 << 2) + 1] * a_0[1];
        z_idx_0 += P_plus[(r2 << 2) + 2] * a_0[2];
        z_idx_0 += P_plus[(r2 << 2) + 3] * a_0[3];
        A_pred[r2 << 2] = S[0] * z_idx_0;
        A_pred[(r2 << 2) + 1] = S[1] * z_idx_0;
        A_pred[(r2 << 2) + 2] = S[2] * z_idx_0;
        A_pred[(r2 << 2) + 3] = S[3] * z_idx_0;
      }

      for (r2 = 0; r2 <= 14; r2 += 2) {
        tmp_0 = _mm_loadu_pd(&P_new[r2]);
        tmp_1 = _mm_loadu_pd(&A_pred[r2]);
        tmp_0 = _mm_sub_pd(tmp_0, tmp_1);
        tmp_1 = _mm_loadu_pd(&tmp[r2]);
        tmp_0 = _mm_add_pd(tmp_0, tmp_1);
        tmp_1 = _mm_loadu_pd(&P_mat[r2]);
        _mm_storeu_pd(&b_y[r2], tmp_1);
        _mm_storeu_pd(&P_new[r2], tmp_0);
      }

      t_square = 0.0;
      for (r2 = 0; r2 < 4; r2++) {
        s3 = b_y[r2 << 2] * Bd[0];
        s3 += b_y[(r2 << 2) + 1] * Bd[1];
        s3 += b_y[(r2 << 2) + 2] * Bd[2];
        s3 += b_y[(r2 << 2) + 3] * Bd[3];
        t_square += s3 * Bd[r2];
        s3 = P_mat[r2 << 2] * Bd[0];
        s3 += P_mat[(r2 << 2) + 1] * Bd[1];
        s3 += P_mat[(r2 << 2) + 2] * Bd[2];
        s3 += P_mat[(r2 << 2) + 3] * Bd[3];
        S[r2] = s3;
      }

      t_square += 0.055;
      for (r2 = 0; r2 < 4; r2++) {
        z_idx_0 = P_plus[r2 << 2] * S[0];
        z_idx_0 += P_plus[(r2 << 2) + 1] * S[1];
        z_idx_0 += P_plus[(r2 << 2) + 2] * S[2];
        z_idx_0 += P_plus[(r2 << 2) + 3] * S[3];
        z_idx_0 /= t_square;
        K_new[r2] = z_idx_0;
      }

      S[0] = prevState[0];
      a_0[0] = K_new[0];
      S[1] = prevState[1];
      a_0[1] = K_new[1];
      S[2] = prevState[2];
      a_0[2] = K_new[2];
      S[3] = prevState[3];
      a_0[3] = K_new[3];
      p = false;
      p_0 = true;
      b_k = 0;
      exitg2 = false;
      while ((!exitg2) && (b_k < 4)) {
        t_square = (real_T)b_k + 1.0;
        z_idx_0 = S[(int32_T)t_square - 1];
        z_idx_1 = a_0[(int32_T)t_square - 1];
        p_1 = (z_idx_0 == z_idx_1);
        if (!p_1) {
          p_0 = false;
          exitg2 = true;
        } else {
          b_k++;
        }
      }

      if (p_0) {
        p = true;
      }

      if (p) {
        exitg1 = true;
      } else {
        memcpy(&P_mat[0], &P_new[0], sizeof(real_T) << 4U);
        prevState[0] = K_new[0];
        prevState[1] = K_new[1];
        prevState[2] = K_new[2];
        prevState[3] = K_new[3];
        r1++;
      }
    }

    /*  Form LQR error vector and nominal control */
    x = stateDerivative[2];
    x = sin(x);
    error_x1 = stateDerivative[2];
    error_x1 = cos(error_x1);
    K_new[0] = stateDerivative[0] - refPos;
    K_new[1] = stateDerivative[1] - refVel;
    K_new[2] = 0.41828772872251141 * x - refAccel;
    K_new[3] = 0.41828772872251141 * stateDerivative[3] * error_x1;
    t_square = prevState[0] * K_new[0];
    t_square += prevState[1] * K_new[1];
    t_square += prevState[2] * K_new[2];
    t_square += prevState[3] * K_new[3];
    t_square = 1.0 / amp * (-t_ratio - t_square);

    /*  Integrator anti-windup: accumulate and clamp error */
    obj->cumulativeError += (stateDerivative[0] - refPos) * deltaTime;
    amp = obj->cumulativeError;
    if (!(amp <= 1.0)) {
      amp = 1.0;
    }

    t_ratio = -1.0;
    if (amp >= -1.0) {
      t_ratio = amp;
    }

    obj->cumulativeError = t_ratio;
    t_ratio = 2.0 * obj->cumulativeError;

    /*  integral correction */
    t_square -= t_ratio;
  }

  /*  Store final control input */
  obj->controlInput = t_square;

  /*  Outputs to Simulink */
  /*  servo voltage command */
  t_ratio = stateDerivative[0];

  /*  estimated ball position */
  amp = stateDerivative[1];

  /*  estimated ball velocity */
  s3 = stateDerivative[2];

  /*  estimated pendulum angle */
  x_0 = stateDerivative[3];

  /*  estimated pendulum angular velocity */
  z_idx_0 = obj->cumulativeError;

  /*  integrator state */
  /*  Update time history for next step */
  obj->previousTime = u0;
  obj->previousDeltaTime = deltaTime;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o1 = t_square;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o2 = t_ratio;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o3 = amp;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o4 = s3;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o5 = x_0;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o6 = z_idx_0;

  /* Saturate: '<Root>/+//-10V' */
  u0 = simulink_experiment_debug_typ_B.MATLABSystem_o1;
  deltaTime = simulink_experiment_debug_typ_P.u0V_LowerSat;
  refPos = simulink_experiment_debug_typ_P.u0V_UpperSat;
  if (u0 > refPos) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = refPos;
  } else if (u0 < deltaTime) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = deltaTime;
  } else {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u0;
  }

  /* End of Saturate: '<Root>/+//-10V' */

  /* Gain: '<S1>/Motor  Gain (V//V)' */
  simulink_experiment_debug_typ_B.MotorGainVV =
    simulink_experiment_debug_typ_P.MotorGainVV_Gain *
    simulink_experiment_debug_typ_B.u0V;

  /* S-Function (hil_write_analog_block): '<S1>/HIL Write Analog' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Write Analog (hil_write_analog_block) */
  {
    t_error result;
    result = hil_write_analog(simulink_experiment_debug_ty_DW.HILInitialize_Card,
      &simulink_experiment_debug_typ_P.HILWriteAnalog_channels, 1,
      &simulink_experiment_debug_typ_B.MotorGainVV);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    }
  }

  /* MATLAB Function: '<Root>/MATLAB Function' */
  /* MATLAB Function 'MATLAB Function': '<S2>:1' */
  /* '<S2>:1:3' */
  if (simulink_experiment_debug_typ_B.Clock < 5.0) {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else if (simulink_experiment_debug_typ_B.Clock < 61.85) {
    t_square = (simulink_experiment_debug_typ_B.Clock - 5.0) / 56.85;
    if (t_square < 0.5) {
      amp = t_square / 0.5 * 0.090000000000000011 + 0.05;
      simulink_experiment_debug_typ_B.v_ref = cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * amp *
        (0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
          * 0.11423973285781065) * 0.2094395102393195) + sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * 0.00316622691292876;
      u0 = 0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock -
        5.0) * 6.2831853071795862 / 55.0) * 3.1415926535897931 / 15.0;
      simulink_experiment_debug_typ_B.a_ref = (cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * 12.0 * (0.83775804095727813 - cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 3.1415926535897931 / 15.0) / 1895.0 + sin(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * ((simulink_experiment_debug_typ_B.Clock -
        5.0) * 6.0 / 1895.0 + 0.05) * (u0 * u0)) + cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * (sin((simulink_experiment_debug_typ_B.Clock
        - 5.0) * 6.2831853071795862 / 55.0) * 19.739208802178716) *
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.0 / 1895.0 + 0.05) /
        825.0;
    } else {
      amp = 0.14;
      simulink_experiment_debug_typ_B.v_ref = cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * 0.14 *
        (0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
          * 0.11423973285781065) * 0.2094395102393195);
      x = 0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
        * 6.2831853071795862 / 55.0) * 3.1415926535897931 / 15.0;
      simulink_experiment_debug_typ_B.a_ref = sin(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * 7.0 * (x * x) / 50.0 + cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * (sin((simulink_experiment_debug_typ_B.Clock
        - 5.0) * 6.2831853071795862 / 55.0) * 69.0872308076255) / 20625.0;
    }

    simulink_experiment_debug_typ_B.p_ref = sin
      ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 - sin
       ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065) *
       0.2094395102393195 / 0.11423973285781065) * amp;
  } else if (simulink_experiment_debug_typ_B.Clock < 65.0) {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else if (simulink_experiment_debug_typ_B.Clock < 85.0) {
    if ((simulink_experiment_debug_typ_B.Clock - 65.0) / 20.0 < 0.5) {
      t_ratio = 0.05;
    } else {
      t_ratio = 0.1;
    }

    u0 = sin((simulink_experiment_debug_typ_B.Clock - 65.0) *
             0.62831853071795862);
    if (rtIsNaN(u0)) {
      u0 = (rtNaN);
    } else if (u0 < 0.0) {
      u0 = -1.0;
    } else {
      u0 = (u0 > 0.0);
    }

    simulink_experiment_debug_typ_B.p_ref = t_ratio * u0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  }

  /* End of MATLAB Function: '<Root>/MATLAB Function' */

  /* Gain: '<Root>/m to cm' */
  /* '<S2>:1:3' */
  simulink_experiment_debug_typ_B.mtocm[0] =
    simulink_experiment_debug_typ_P.mtocm_Gain *
    simulink_experiment_debug_typ_B.p_ref;
  simulink_experiment_debug_typ_B.mtocm[1] =
    simulink_experiment_debug_typ_P.mtocm_Gain *
    simulink_experiment_debug_typ_B.BB01SensorGainmV;

  /* Gain: '<S3>/Gain' */
  simulink_experiment_debug_typ_B.Gain =
    simulink_experiment_debug_typ_P.Gain_Gain *
    simulink_experiment_debug_typ_B.Bias;

  /* RateTransition: '<Root>/Rate Transition' */
  if (simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2) {
    simulink_experiment_debug_ty_DW.RateTransition_Buffer =
      simulink_experiment_debug_typ_B.Clock;

    /* RateTransition: '<Root>/Rate Transition1' */
    simulink_experiment_debug_ty_DW.RateTransition1_Buffer =
      simulink_experiment_debug_typ_B.p_ref;

    /* RateTransition: '<Root>/Rate Transition2' */
    simulink_experiment_debug_ty_DW.RateTransition2_Buffer =
      simulink_experiment_debug_typ_B.MATLABSystem_o1;

    /* RateTransition: '<Root>/Rate Transition3' */
    simulink_experiment_debug_ty_DW.RateTransition3_Buffer =
      simulink_experiment_debug_typ_B.BB01SensorGainmV;

    /* RateTransition: '<Root>/Rate Transition4' */
    simulink_experiment_debug_ty_DW.RateTransition4_Buffer =
      simulink_experiment_debug_typ_B.Bias;
  }

  /* End of RateTransition: '<Root>/Rate Transition' */
}

/* Model update function for TID0 */
void simulink_experiment_debug_type1_update0(void) /* Sample time: [0.0s, 0.0s] */
{
  /* Update absolute time */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick0)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH0;
  }

  simulink_experiment_debug_ty_M->Timing.t[0] =
    simulink_experiment_debug_ty_M->Timing.clockTick0 *
    simulink_experiment_debug_ty_M->Timing.stepSize0 +
    simulink_experiment_debug_ty_M->Timing.clockTickH0 *
    simulink_experiment_debug_ty_M->Timing.stepSize0 * 4294967296.0;

  /* Update absolute time */
  /* The "clockTick1" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick1"
   * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick1 and the high bits
   * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick1)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH1;
  }

  simulink_experiment_debug_ty_M->Timing.t[1] =
    simulink_experiment_debug_ty_M->Timing.clockTick1 *
    simulink_experiment_debug_ty_M->Timing.stepSize1 +
    simulink_experiment_debug_ty_M->Timing.clockTickH1 *
    simulink_experiment_debug_ty_M->Timing.stepSize1 * 4294967296.0;
}

/* Model output function for TID2 */
void simulink_experiment_debug_type1_output2(void) /* Sample time: [0.01s, 0.0s] */
{
  /* RateTransition: '<Root>/Rate Transition2' */
  simulink_experiment_debug_typ_B.RateTransition2 =
    simulink_experiment_debug_ty_DW.RateTransition2_Buffer;

  /* RateTransition: '<Root>/Rate Transition1' */
  simulink_experiment_debug_typ_B.RateTransition1 =
    simulink_experiment_debug_ty_DW.RateTransition1_Buffer;

  /* RateTransition: '<Root>/Rate Transition3' */
  simulink_experiment_debug_typ_B.RateTransition3 =
    simulink_experiment_debug_ty_DW.RateTransition3_Buffer;

  /* RateTransition: '<Root>/Rate Transition4' */
  simulink_experiment_debug_typ_B.RateTransition4 =
    simulink_experiment_debug_ty_DW.RateTransition4_Buffer;

  /* RateTransition: '<Root>/Rate Transition' */
  simulink_experiment_debug_typ_B.RateTransition =
    simulink_experiment_debug_ty_DW.RateTransition_Buffer;
}

/* Model update function for TID2 */
void simulink_experiment_debug_type1_update2(void) /* Sample time: [0.01s, 0.0s] */
{
  /* Update absolute time */
  /* The "clockTick2" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick2"
   * and "Timing.stepSize2". Size of "clockTick2" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick2 and the high bits
   * Timing.clockTickH2. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick2)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH2;
  }

  simulink_experiment_debug_ty_M->Timing.t[2] =
    simulink_experiment_debug_ty_M->Timing.clockTick2 *
    simulink_experiment_debug_ty_M->Timing.stepSize2 +
    simulink_experiment_debug_ty_M->Timing.clockTickH2 *
    simulink_experiment_debug_ty_M->Timing.stepSize2 * 4294967296.0;
}

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void simulink_experiment_debug_type1_output(int_T tid)
{
  switch (tid) {
   case 0 :
    simulink_experiment_debug_type1_output0();
    break;

   case 2 :
    simulink_experiment_debug_type1_output2();
    break;

   default :
    /* do nothing */
    break;
  }
}

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void simulink_experiment_debug_type1_update(int_T tid)
{
  switch (tid) {
   case 0 :
    simulink_experiment_debug_type1_update0();
    break;

   case 2 :
    simulink_experiment_debug_type1_update2();
    break;

   default :
    /* do nothing */
    break;
  }
}

/* Model initialize function */
void simulink_experiment_debug_type1_initialize(void)
{
  {
    studentControllerInterface_si_T *b_obj;
    int32_T i;
    static const real_T tmp[16] = { 0.01, 0.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0,
      0.0, 0.0, 0.01, 0.0, 0.0, 0.0, 0.0, 0.01 };

    /* Start for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

    /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
    {
      t_int result;
      t_boolean is_switching;
      result = hil_open("q2_usb", "0",
                        &simulink_experiment_debug_ty_DW.HILInitialize_Card);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      is_switching = false;
      result = hil_set_card_specific_options
        (simulink_experiment_debug_ty_DW.HILInitialize_Card,
         "d0=digital;d1=digital;led=auto;update_rate=normal", 50);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      result = hil_watchdog_clear
        (simulink_experiment_debug_ty_DW.HILInitialize_Card);
      if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AIPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AIPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0] =
          (simulink_experiment_debug_typ_P.HILInitialize_AILow);
        simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[1] =
          (simulink_experiment_debug_typ_P.HILInitialize_AILow);
        simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AIHigh;
        simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AIHigh;
        result = hil_set_analog_input_ranges
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AIChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0],
           &simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AOPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AOPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0] =
          (simulink_experiment_debug_typ_P.HILInitialize_AOLow);
        simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[1] =
          (simulink_experiment_debug_typ_P.HILInitialize_AOLow);
        simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOHigh;
        simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOHigh;
        result = hil_set_analog_output_ranges
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0],
           &simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AOStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AOEnter && is_switching))
      {
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOInitial;
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOInitial;
        result = hil_write_analog
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if (simulink_experiment_debug_typ_P.HILInitialize_AOReset) {
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOWatchdog;
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOWatchdog;
        result = hil_watchdog_set_analog_expiration_state
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      result = hil_set_digital_directions
        (simulink_experiment_debug_ty_DW.HILInitialize_Card, NULL, 0U,
         simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_DOStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_DOEnter && is_switching))
      {
        {
          int_T i1;
          boolean_T *dw_DOBits =
            &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0];
          for (i1=0; i1 < 8; i1++) {
            dw_DOBits[i1] =
              simulink_experiment_debug_typ_P.HILInitialize_DOInitial;
          }
        }

        result = hil_write_digital
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U,
           (t_boolean *) &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if (simulink_experiment_debug_typ_P.HILInitialize_DOReset) {
        {
          int_T i1;
          int32_T *dw_DOStates =
            &simulink_experiment_debug_ty_DW.HILInitialize_DOStates[0];
          for (i1=0; i1 < 8; i1++) {
            dw_DOStates[i1] =
              simulink_experiment_debug_typ_P.HILInitialize_DOWatchdog;
          }
        }

        result = hil_watchdog_set_digital_expiration_state
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U, (const
            t_digital_state *)
           &simulink_experiment_debug_ty_DW.HILInitialize_DOStates[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_EIPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_EIPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[0] =
          simulink_experiment_debug_typ_P.HILInitialize_EIQuadrature;
        simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[1] =
          simulink_experiment_debug_typ_P.HILInitialize_EIQuadrature;
        result = hil_set_encoder_quadrature_mode
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_EIChannels, 2U,
           (t_encoder_quadrature_mode *)
           &simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_EIStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_EIEnter && is_switching))
      {
        simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[0] =
          simulink_experiment_debug_typ_P.HILInitialize_EIInitial;
        simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[1] =
          simulink_experiment_debug_typ_P.HILInitialize_EIInitial;
        result = hil_set_encoder_counts
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_EIChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }
    }

    /* Start for S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

    /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_create_encoder_reader
        (simulink_experiment_debug_ty_DW.HILInitialize_Card,
         simulink_experiment_debug_typ_P.HILReadEncoderTimebase_SamplesI,
         &simulink_experiment_debug_typ_P.HILReadEncoderTimebase_Channels, 1,
         &simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task);
      if (result >= 0) {
        result = hil_task_set_buffer_overflow_mode
          (simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task,
           (t_buffer_overflow_mode)
           (simulink_experiment_debug_typ_P.HILReadEncoderTimebase_Overflow - 1));
      }

      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
      }
    }

    /* Start for MATLABSystem: '<Root>/MATLAB System' */
    b_obj = &simulink_experiment_debug_ty_DW.obj;
    b_obj->previousTime = 0.0;
    b_obj->previousDeltaTime = 0.001;
    b_obj->controlInput = 0.0;
    b_obj->stateEstimate[0] = 0.0;
    b_obj->stateEstimate[1] = 0.0;
    b_obj->stateEstimate[2] = 0.0;
    b_obj->stateEstimate[3] = 0.0;
    for (i = 0; i < 16; i++) {
      b_obj->estimateCovariance[i] = tmp[i];
    }

    b_obj->cumulativeError = 0.0;
    b_obj->const_1 = 0.0;
    b_obj->const_2 = 0.0;
    b_obj->useFeedbackLinearization = true;
    simulink_experiment_debug_ty_DW.objisempty = true;

    /* End of Start for MATLABSystem: '<Root>/MATLAB System' */
  }
}

/* Model terminate function */
void simulink_experiment_debug_type1_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_digital_outputs = 0;
    hil_task_stop_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_monitor_stop_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    is_switching = false;
    if ((simulink_experiment_debug_typ_P.HILInitialize_AOTerminate &&
         !is_switching) || (simulink_experiment_debug_typ_P.HILInitialize_AOExit
         && is_switching)) {
      simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
        simulink_experiment_debug_typ_P.HILInitialize_AOFinal;
      simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
        simulink_experiment_debug_typ_P.HILInitialize_AOFinal;
      num_final_analog_outputs = 2U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((simulink_experiment_debug_typ_P.HILInitialize_DOTerminate &&
         !is_switching) || (simulink_experiment_debug_typ_P.HILInitialize_DOExit
         && is_switching)) {
      {
        int_T i1;
        boolean_T *dw_DOBits =
          &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOBits[i1] = simulink_experiment_debug_typ_P.HILInitialize_DOFinal;
        }
      }

      num_final_digital_outputs = 8U;
    } else {
      num_final_digital_outputs = 0;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_digital_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(simulink_experiment_debug_ty_DW.HILInitialize_Card
                         ,
                         simulink_experiment_debug_typ_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , NULL, 0
                         ,
                         simulink_experiment_debug_typ_P.HILInitialize_DOChannels,
                         num_final_digital_outputs
                         , NULL, 0
                         ,
                         &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages
                         [0]
                         , NULL
                         , (t_boolean *)
                         &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (simulink_experiment_debug_ty_DW.HILInitialize_Card,
             simulink_experiment_debug_typ_P.HILInitialize_AOChannels,
             num_final_analog_outputs,
             &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_digital_outputs > 0) {
          local_result = hil_write_digital
            (simulink_experiment_debug_ty_DW.HILInitialize_Card,
             simulink_experiment_debug_typ_P.HILInitialize_DOChannels,
             num_final_digital_outputs, (t_boolean *)
             &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_monitor_delete_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_close(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    simulink_experiment_debug_ty_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  if (tid == 1)
    tid = 0;
  simulink_experiment_debug_type1_output(tid);
}

void MdlUpdate(int_T tid)
{
  if (tid == 1)
    tid = 0;
  simulink_experiment_debug_type1_update(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  simulink_experiment_debug_type1_initialize();
}

void MdlTerminate(void)
{
  simulink_experiment_debug_type1_terminate();
}

/* Registration function */
RT_MODEL_simulink_experiment__T *simulink_experiment_debug_type1(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)simulink_experiment_debug_ty_M, 0,
                sizeof(RT_MODEL_simulink_experiment__T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&simulink_experiment_debug_ty_M->solverInfo,
                          &simulink_experiment_debug_ty_M->Timing.simTimeStep);
    rtsiSetTPtr(&simulink_experiment_debug_ty_M->solverInfo, &rtmGetTPtr
                (simulink_experiment_debug_ty_M));
    rtsiSetStepSizePtr(&simulink_experiment_debug_ty_M->solverInfo,
                       &simulink_experiment_debug_ty_M->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&simulink_experiment_debug_ty_M->solverInfo,
                          (&rtmGetErrorStatus(simulink_experiment_debug_ty_M)));
    rtsiSetRTModelPtr(&simulink_experiment_debug_ty_M->solverInfo,
                      simulink_experiment_debug_ty_M);
  }

  rtsiSetSimTimeStep(&simulink_experiment_debug_ty_M->solverInfo,
                     MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange
    (&simulink_experiment_debug_ty_M->solverInfo, false);
  rtsiSetSolverName(&simulink_experiment_debug_ty_M->solverInfo,
                    "FixedStepDiscrete");

  /* Initialize timing info */
  {
    int_T *mdlTsMap =
      simulink_experiment_debug_ty_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    mdlTsMap[2] = 2;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "simulink_experiment_debug_ty_M points to
       static memory which is guaranteed to be non-NULL" */
    simulink_experiment_debug_ty_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    simulink_experiment_debug_ty_M->Timing.sampleTimes =
      (&simulink_experiment_debug_ty_M->Timing.sampleTimesArray[0]);
    simulink_experiment_debug_ty_M->Timing.offsetTimes =
      (&simulink_experiment_debug_ty_M->Timing.offsetTimesArray[0]);

    /* task periods */
    simulink_experiment_debug_ty_M->Timing.sampleTimes[0] = (0.0);
    simulink_experiment_debug_ty_M->Timing.sampleTimes[1] = (0.002);
    simulink_experiment_debug_ty_M->Timing.sampleTimes[2] = (0.01);

    /* task offsets */
    simulink_experiment_debug_ty_M->Timing.offsetTimes[0] = (0.0);
    simulink_experiment_debug_ty_M->Timing.offsetTimes[1] = (0.0);
    simulink_experiment_debug_ty_M->Timing.offsetTimes[2] = (0.0);
  }

  rtmSetTPtr(simulink_experiment_debug_ty_M,
             &simulink_experiment_debug_ty_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = simulink_experiment_debug_ty_M->Timing.sampleHitArray;
    int_T *mdlPerTaskSampleHits =
      simulink_experiment_debug_ty_M->Timing.perTaskSampleHitsArray;
    simulink_experiment_debug_ty_M->Timing.perTaskSampleHits =
      (&mdlPerTaskSampleHits[0]);
    mdlSampleHits[0] = 1;
    simulink_experiment_debug_ty_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(simulink_experiment_debug_ty_M, 20.0);
  simulink_experiment_debug_ty_M->Timing.stepSize0 = 0.002;
  simulink_experiment_debug_ty_M->Timing.stepSize1 = 0.002;
  simulink_experiment_debug_ty_M->Timing.stepSize2 = 0.01;

  /* External mode info */
  simulink_experiment_debug_ty_M->Sizes.checksums[0] = (578646490U);
  simulink_experiment_debug_ty_M->Sizes.checksums[1] = (3948803095U);
  simulink_experiment_debug_ty_M->Sizes.checksums[2] = (10994780U);
  simulink_experiment_debug_ty_M->Sizes.checksums[3] = (3569915800U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[3];
    simulink_experiment_debug_ty_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    systemRan[2] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(simulink_experiment_debug_ty_M->extModeInfo,
      &simulink_experiment_debug_ty_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(simulink_experiment_debug_ty_M->extModeInfo,
                        simulink_experiment_debug_ty_M->Sizes.checksums);
    rteiSetTPtr(simulink_experiment_debug_ty_M->extModeInfo, rtmGetTPtr
                (simulink_experiment_debug_ty_M));
  }

  simulink_experiment_debug_ty_M->solverInfoPtr =
    (&simulink_experiment_debug_ty_M->solverInfo);
  simulink_experiment_debug_ty_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&simulink_experiment_debug_ty_M->solverInfo, 0.002);
  rtsiSetSolverMode(&simulink_experiment_debug_ty_M->solverInfo,
                    SOLVER_MODE_MULTITASKING);

  /* block I/O */
  simulink_experiment_debug_ty_M->blockIO = ((void *)
    &simulink_experiment_debug_typ_B);

  {
    simulink_experiment_debug_typ_B.HILReadEncoderTimebase = 0.0;
    simulink_experiment_debug_typ_B.HILReadAnalog = 0.0;
    simulink_experiment_debug_typ_B.BB01SensorGainmV = 0.0;
    simulink_experiment_debug_typ_B.EncoderCalibrationradcount = 0.0;
    simulink_experiment_debug_typ_B.Bias = 0.0;
    simulink_experiment_debug_typ_B.Clock = 0.0;
    simulink_experiment_debug_typ_B.u0V = 0.0;
    simulink_experiment_debug_typ_B.MotorGainVV = 0.0;
    simulink_experiment_debug_typ_B.mtocm[0] = 0.0;
    simulink_experiment_debug_typ_B.mtocm[1] = 0.0;
    simulink_experiment_debug_typ_B.Gain = 0.0;
    simulink_experiment_debug_typ_B.RateTransition2 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition1 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition3 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition4 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o1 = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o2 = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o3 = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o4 = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o5 = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o6 = 0.0;
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  }

  /* parameters */
  simulink_experiment_debug_ty_M->defaultParam = ((real_T *)
    &simulink_experiment_debug_typ_P);

  /* states (dwork) */
  simulink_experiment_debug_ty_M->dwork = ((void *)
    &simulink_experiment_debug_ty_DW);
  (void) memset((void *)&simulink_experiment_debug_ty_DW, 0,
                sizeof(DW_simulink_experiment_debug__T));
  simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_FilterFrequency[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_FilterFrequency[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition1_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition2_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition3_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition4_Buffer = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    simulink_experiment_debug_ty_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 22;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  simulink_experiment_debug_ty_M->Sizes.numContStates = (0);/* Number of continuous states */
  simulink_experiment_debug_ty_M->Sizes.numY = (0);/* Number of model outputs */
  simulink_experiment_debug_ty_M->Sizes.numU = (0);/* Number of model inputs */
  simulink_experiment_debug_ty_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  simulink_experiment_debug_ty_M->Sizes.numSampTimes = (3);/* Number of sample times */
  simulink_experiment_debug_ty_M->Sizes.numBlocks = (32);/* Number of blocks */
  simulink_experiment_debug_ty_M->Sizes.numBlockIO = (24);/* Number of block outputs */
  simulink_experiment_debug_ty_M->Sizes.numBlockPrms = (88);/* Sum of parameter "widths" */
  return simulink_experiment_debug_ty_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
