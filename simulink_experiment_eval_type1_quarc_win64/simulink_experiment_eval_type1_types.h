/*
 * simulink_experiment_eval_type1_types.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "simulink_experiment_eval_type1".
 *
 * Model version              : 1.9
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C source code generated on : Wed Apr  9 13:35:22 2025
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_simulink_experiment_eval_type1_types_h_
#define RTW_HEADER_simulink_experiment_eval_type1_types_h_
#include "rtwtypes.h"
#ifndef struct_tag_XNpIO79DGqXtd5S6MMzAiG
#define struct_tag_XNpIO79DGqXtd5S6MMzAiG

struct tag_XNpIO79DGqXtd5S6MMzAiG
{
  int32_T isInitialized;
  real_T t_prev;
  real_T V_servo;
  boolean_T fl_lin;
  real_T ball_rad;
  real_T beam_len;
  real_T g_val;
  real_T servo_gain;
  real_T tau_val;
  real_T const_1;
  real_T const_2;
  real_T state_estimate[4];
  real_T K[4];
  real_T C[8];
  real_T observer_gain[8];
};

#endif                                 /* struct_tag_XNpIO79DGqXtd5S6MMzAiG */

#ifndef typedef_studentControllerInterface_si_T
#define typedef_studentControllerInterface_si_T

typedef struct tag_XNpIO79DGqXtd5S6MMzAiG studentControllerInterface_si_T;

#endif                             /* typedef_studentControllerInterface_si_T */

/* Parameters (default storage) */
typedef struct P_simulink_experiment_eval_ty_T_ P_simulink_experiment_eval_ty_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_simulink_experiment_e_T RT_MODEL_simulink_experiment__T;

#endif                  /* RTW_HEADER_simulink_experiment_eval_type1_types_h_ */
