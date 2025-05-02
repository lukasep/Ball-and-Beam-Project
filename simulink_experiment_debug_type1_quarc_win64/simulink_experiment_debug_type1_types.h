/*
 * simulink_experiment_debug_type1_types.h
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

#ifndef RTW_HEADER_simulink_experiment_debug_type1_types_h_
#define RTW_HEADER_simulink_experiment_debug_type1_types_h_
#include "rtwtypes.h"
#ifndef struct_tag_5b1iKqDGDha33LAulW7mb
#define struct_tag_5b1iKqDGDha33LAulW7mb

struct tag_5b1iKqDGDha33LAulW7mb
{
  real_T previousTime;
  real_T previousDeltaTime;
  real_T controlInput;
  real_T stateEstimate[4];
  real_T estimateCovariance[16];
  real_T cumulativeError;
  real_T const_1;
  real_T const_2;
  boolean_T useFeedbackLinearization;
};

#endif                                 /* struct_tag_5b1iKqDGDha33LAulW7mb */

#ifndef typedef_studentControllerInterface_si_T
#define typedef_studentControllerInterface_si_T

typedef struct tag_5b1iKqDGDha33LAulW7mb studentControllerInterface_si_T;

#endif                             /* typedef_studentControllerInterface_si_T */

/* Parameters (default storage) */
typedef struct P_simulink_experiment_debug_t_T_ P_simulink_experiment_debug_t_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_simulink_experiment_d_T RT_MODEL_simulink_experiment__T;

#endif                 /* RTW_HEADER_simulink_experiment_debug_type1_types_h_ */
