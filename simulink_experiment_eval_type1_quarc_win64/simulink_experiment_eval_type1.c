/*
 * simulink_experiment_eval_type1.c
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

#include "simulink_experiment_eval_type1.h"
#include "simulink_experiment_eval_type1_types.h"
#include "rtwtypes.h"
#include <math.h>
#include <emmintrin.h>
#include "simulink_experiment_eval_type1_private.h"
#include "rt_nonfinite.h"
#include <string.h>
#include "simulink_experiment_eval_type1_dt.h"

/* Block signals (default storage) */
B_simulink_experiment_eval_ty_T simulink_experiment_eval_type_B;

/* Block states (default storage) */
DW_simulink_experiment_eval_t_T simulink_experiment_eval_typ_DW;

/* Real-time model */
static RT_MODEL_simulink_experiment__T simulink_experiment_eval_typ_M_;
RT_MODEL_simulink_experiment__T *const simulink_experiment_eval_typ_M =
  &simulink_experiment_eval_typ_M_;

/* Forward declaration for local functions */
static studentControllerInterface_si_T *studentControllerInterface_stud
  (studentControllerInterface_si_T *obj);
static void rate_monotonic_scheduler(void);
time_T rt_SimUpdateDiscreteEvents(
  int_T rtmNumSampTimes, void *rtmTimingData, int_T *rtmSampleHitPtr, int_T
  *rtmPerTaskSampleHits )
{
  rtmSampleHitPtr[1] = rtmStepTask(simulink_experiment_eval_typ_M, 1);
  rtmSampleHitPtr[2] = rtmStepTask(simulink_experiment_eval_typ_M, 2);
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
  if (simulink_experiment_eval_typ_M->Timing.TaskCounters.TID[1] == 0) {
    simulink_experiment_eval_typ_M->Timing.RateInteraction.TID1_2 =
      (simulink_experiment_eval_typ_M->Timing.TaskCounters.TID[2] == 0);

    /* update PerTaskSampleHits matrix for non-inline sfcn */
    simulink_experiment_eval_typ_M->Timing.perTaskSampleHits[5] =
      simulink_experiment_eval_typ_M->Timing.RateInteraction.TID1_2;
  }

  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (simulink_experiment_eval_typ_M->Timing.TaskCounters.TID[2])++;
  if ((simulink_experiment_eval_typ_M->Timing.TaskCounters.TID[2]) > 4) {/* Sample time: [0.01s, 0.0s] */
    simulink_experiment_eval_typ_M->Timing.TaskCounters.TID[2] = 0;
  }
}

static studentControllerInterface_si_T *studentControllerInterface_stud
  (studentControllerInterface_si_T *obj)
{
  studentControllerInterface_si_T *b_obj;
  real_T a;
  int32_T i;
  static const int8_T tmp[8] = { 1, 0, 0, 0, 0, 1, 0, 0 };

  static const real_T tmp_0[8] = { 22.9338, 130.7798, 0.957, 10.983, 1.0388,
    12.977, 23.0662, 132.2187 };

  b_obj = obj;
  b_obj->t_prev = -1.0;
  b_obj->V_servo = 0.0;
  b_obj->fl_lin = true;
  b_obj->ball_rad = 0.0254;
  b_obj->beam_len = 0.4255;
  b_obj->g_val = 9.81;
  b_obj->servo_gain = 1.5;
  b_obj->tau_val = 0.025;
  b_obj->state_estimate[0] = 0.0;
  b_obj->state_estimate[1] = 0.0;
  b_obj->state_estimate[2] = 0.0;
  b_obj->state_estimate[3] = 0.0;
  for (i = 0; i < 8; i++) {
    b_obj->C[i] = tmp[i];
  }

  for (i = 0; i < 8; i++) {
    b_obj->observer_gain[i] = tmp_0[i];
  }

  b_obj->isInitialized = 0;

  /*  Constructor to initialize the controller */
  b_obj->const_1 = 5.0 * b_obj->g_val * b_obj->ball_rad / (7.0 * b_obj->beam_len);
  a = b_obj->ball_rad / b_obj->beam_len;
  a *= a;
  b_obj->const_2 = 0.7142857142857143 * a;

  /*  % Set LQR weights and compute gains */
  /*  obj.Q(1,1) = 10;  % Position error weight */
  /*  obj.Q(2,2) = 70;   % Velocity error weight */
  /*  obj.R = 0.1;         % Control effort weight */
  /*  Calculate LQR gain for the Feedback Linearization method */
  /* obj.K = lqr(obj.A, obj.B, obj.Q, obj.R) */
  /* obj.K = [23.4521   29.3025   17.2402    5.8720]; % MATLAB SINE WAVE-- Score: 0.84 */
  b_obj->K[0] = 100.0;
  b_obj->K[1] = 89.3115;
  b_obj->K[2] = 36.3827;
  b_obj->K[3] = 8.5303;

  /*  MATLAB SQUARE WAVE-- Score: 3.33 */
  /*              obj. K = [3.0002   6.0101    8.0129   4.0423]; % SIMULINK SINE WAVE-- Score: 1.5 */
  /* obj.K = [3.0102   6.1143    7.9920   4.0120]; % SIMULINK SQUARE WAVE-- 4.37 */
  /*  Display which controller is active */
  return b_obj;
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

/* Model output function for TID0 */
void simulink_experiment_eval_type1_output0(void) /* Sample time: [0.0s, 0.0s] */
{
  __m128d tmp;
  __m128d tmp_0;
  studentControllerInterface_si_T *obj;
  studentControllerInterface_si_T *obj_0;
  studentControllerInterface_si_T *obj_1;
  real_T a[8];
  real_T error[4];
  real_T state_derivative[4];
  real_T measurements[2];
  real_T a_ball_ref;
  real_T amplitude;
  real_T b_x;
  real_T b_x_0;
  real_T b_x_1;
  real_T b_x_2;
  real_T b_x_3;
  real_T b_x_4;
  real_T b_x_5;
  real_T b_x_6;
  real_T b_x_7;
  real_T ball_velocity;
  real_T beam_angle;
  real_T beam_angular_velocity;
  real_T br;
  real_T c;
  real_T c_0;
  real_T c_1;
  real_T c_2;
  real_T c_3;
  real_T c_4;
  real_T c_5;
  real_T c_6;
  real_T dt;
  real_T error_x1;
  real_T error_x2;
  real_T error_x3;
  real_T omega;
  real_T p_ball_ref;
  real_T u0;
  real_T u1;
  real_T u2;
  real_T v_ball_ref;
  int32_T i;

  {                                    /* Sample time: [0.0s, 0.0s] */
    rate_monotonic_scheduler();
  }

  /* S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

  /* S-Function Block: simulink_experiment_eval_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_read_encoder
      (simulink_experiment_eval_typ_DW.HILReadEncoderTimebase_Task, 1,
       &simulink_experiment_eval_typ_DW.HILReadEncoderTimebase_Buffer);
    if (result < 0) {
      simulink_experiment_eval_type_B.HILReadEncoderTimebase = 0;
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
    } else {
      simulink_experiment_eval_type_B.HILReadEncoderTimebase =
        simulink_experiment_eval_typ_DW.HILReadEncoderTimebase_Buffer;
    }
  }

  /* S-Function (hil_read_analog_block): '<S1>/HIL Read Analog' */

  /* S-Function Block: simulink_experiment_eval_type1/Ball and Beam Hardware Interface/HIL Read Analog (hil_read_analog_block) */
  {
    t_error result = hil_read_analog
      (simulink_experiment_eval_typ_DW.HILInitialize_Card,
       &simulink_experiment_eval_type_P.HILReadAnalog_channels, 1,
       &simulink_experiment_eval_typ_DW.HILReadAnalog_Buffer);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
    }

    simulink_experiment_eval_type_B.HILReadAnalog =
      simulink_experiment_eval_typ_DW.HILReadAnalog_Buffer;
  }

  /* Gain: '<S1>/BB01 Sensor  Gain (m//V)' */
  simulink_experiment_eval_type_B.BB01SensorGainmV =
    simulink_experiment_eval_type_P.BB01SensorGainmV_Gain *
    simulink_experiment_eval_type_B.HILReadAnalog;

  /* Gain: '<S1>/Encoder Calibration  (rad//count)' */
  simulink_experiment_eval_type_B.EncoderCalibrationradcount =
    simulink_experiment_eval_type_P.EncoderCalibrationradcount_Gain *
    simulink_experiment_eval_type_B.HILReadEncoderTimebase;

  /* Bias: '<S1>/Bias' */
  simulink_experiment_eval_type_B.Bias =
    simulink_experiment_eval_type_B.EncoderCalibrationradcount +
    simulink_experiment_eval_type_P.Bias_Bias;

  /* Clock: '<Root>/Clock' */
  simulink_experiment_eval_type_B.Clock =
    simulink_experiment_eval_typ_M->Timing.t[0];

  /* MATLABSystem: '<Root>/MATLAB System' */
  u0 = simulink_experiment_eval_type_B.Clock;
  u1 = simulink_experiment_eval_type_B.BB01SensorGainmV;
  u2 = simulink_experiment_eval_type_B.Bias;
  obj = &simulink_experiment_eval_typ_DW.obj;

  /*  This is the main function called every iteration. You have to implement */
  /*  the controller in this function, bu you are not allowed to */
  /*  change the signature of this function. */
  /*  Input arguments: */
  /*    t: current time */
  /*    ball_pos: position of the ball provided by the ball position sensor (m) */
  /*    ball_vel: velocity of the ball (m/s) */
  /*    beam_angle: servo motor angle provided by the encoder of the motor (rad) */
  /*    beam_angular_vel: angular velocity of the beam (rad/s) */
  /*  Output: */
  /*    V_servo: voltage to the servo input */
  /*  Safety Params */
  /*  Calculate time step */
  dt = u0 - obj->t_prev;
  if (obj->t_prev < 0.0) {
    dt = 0.001;
  }

  /*  Get reference trajectory */
  /*  [p_ref, v_ref, a_ref] = get_ref_traj(t) */
  /*  Reference trajectory at time t: */
  /*    Inputs */
  /*        t: query time */
  /*    Outputs */
  /*        p_ref: reference position of the ball */
  /*        v_ref: reference velocity of the ball */
  /*        a_ref: reference acceleration of the ball */
  /*  m */
  /*  sec */
  /*     %% Sine wave. */
  beam_angle = 0.62831853071795862 * u0;
  beam_angle = sin(beam_angle);
  p_ball_ref = 0.04 * beam_angle;
  beam_angle = 0.62831853071795862 * u0;
  beam_angle = cos(beam_angle);
  v_ball_ref = 0.025132741228718346 * beam_angle;
  beam_angle = 0.62831853071795862 * u0;
  beam_angle = sin(beam_angle);
  a_ball_ref = -0.015791367041742971 * beam_angle;

  /*     %% Square wave. */
  /*      p_ref = amplitude * sign(sin(omega * t)); */
  /*      v_ref = 0; */
  /*      a_ref = 0; */
  /*  Measurements available from sensors */
  measurements[0] = u1;
  measurements[1] = u2;

  /*  Get the previous control input */
  amplitude = obj->V_servo;

  /*  Update state estimate using the Luenberger observer */
  obj_0 = obj;
  error[0] = obj->state_estimate[0];
  error[1] = obj->state_estimate[1];
  error[2] = obj->state_estimate[2];
  error[3] = obj->state_estimate[3];

  /*  Luenberger Observer implementation */
  /*  Extract state variables for clarity */
  omega = error[0];
  ball_velocity = error[1];
  beam_angle = error[2];
  beam_angular_velocity = error[3];

  /*  Calculate the nonlinear system dynamics */
  state_derivative[0] = ball_velocity;
  ball_velocity = beam_angular_velocity * beam_angular_velocity;
  br = beam_angle;
  br = cos(br);
  u1 = br * br;
  beam_angle = sin(beam_angle);
  state_derivative[1] = obj_0->const_1 * beam_angle - (obj_0->beam_len / 2.0 -
    omega) * obj_0->const_2 * ball_velocity * u1;
  state_derivative[2] = beam_angular_velocity;
  state_derivative[3] = (obj_0->servo_gain * amplitude - beam_angular_velocity) /
    obj_0->tau_val;

  /*  Add correction term based on measurement error */
  for (i = 0; i < 8; i++) {
    a[i] = obj_0->C[i];
  }

  for (i = 0; i <= 0; i += 2) {
    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp = _mm_loadu_pd(&a[i]);
    tmp = _mm_mul_pd(tmp, _mm_set1_pd(error[0]));
    tmp = _mm_add_pd(tmp, _mm_set1_pd(0.0));
    tmp_0 = _mm_loadu_pd(&a[i + 2]);
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(error[1]));
    tmp = _mm_add_pd(tmp_0, tmp);

    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_loadu_pd(&a[i + 4]);
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(error[2]));
    tmp = _mm_add_pd(tmp_0, tmp);

    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_loadu_pd(&a[i + 6]);
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(error[3]));
    tmp = _mm_add_pd(tmp_0, tmp);

    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_loadu_pd(&measurements[i]);
    tmp = _mm_sub_pd(tmp_0, tmp);

    /* MATLABSystem: '<Root>/MATLAB System' */
    _mm_storeu_pd(&measurements[i], tmp);
  }

  /* MATLABSystem: '<Root>/MATLAB System' */
  for (i = 0; i < 8; i++) {
    a[i] = obj_0->observer_gain[i];
  }

  /*  Update the state estimate using forward Euler integration */
  for (i = 0; i <= 2; i += 2) {
    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp = _mm_loadu_pd(&a[i]);
    tmp = _mm_mul_pd(tmp, _mm_set1_pd(measurements[0]));
    tmp = _mm_add_pd(tmp, _mm_set1_pd(0.0));
    tmp_0 = _mm_loadu_pd(&a[i + 4]);
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(measurements[1]));
    tmp = _mm_add_pd(tmp_0, tmp);

    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_loadu_pd(&state_derivative[i]);
    tmp = _mm_add_pd(tmp_0, tmp);

    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp = _mm_mul_pd(tmp, _mm_set1_pd(dt));
    tmp_0 = _mm_loadu_pd(&error[i]);
    tmp_0 = _mm_add_pd(tmp_0, tmp);

    /* MATLABSystem: '<Root>/MATLAB System' */
    _mm_storeu_pd(&error[i], tmp_0);
    _mm_storeu_pd(&state_derivative[i], tmp);
  }

  /* MATLABSystem: '<Root>/MATLAB System' */
  obj->state_estimate[0] = error[0];
  obj->state_estimate[1] = error[1];
  obj->state_estimate[2] = error[2];
  obj->state_estimate[3] = error[3];

  /*  Extract state estimates */
  dt = obj->state_estimate[0];
  amplitude = obj->state_estimate[1];
  omega = obj->state_estimate[2];
  beam_angle = obj->state_estimate[3];

  /*  Use estimated states for control */
  if (obj->fl_lin) {
    /*  Feedback linearization controller */
    /*  Compute the error state */
    error_x1 = dt - p_ball_ref;
    error_x2 = amplitude - v_ball_ref;
    obj_0 = obj;

    /*  Second Lie derivative */
    ball_velocity = beam_angle * beam_angle;
    br = omega;
    br = cos(br);
    u1 = br * br;
    br = omega;
    br = sin(br);
    p_ball_ref = obj_0->const_1 * br - (obj_0->beam_len / 2.0 - dt) *
      obj_0->const_2 * ball_velocity * u1;
    error_x3 = p_ball_ref - a_ball_ref;
    obj_0 = obj;

    /*  Third Lie derivative */
    ball_velocity = beam_angle * beam_angle;
    br = omega;
    br = cos(br);
    u1 = br * br;
    v_ball_ref = rt_powd_snf(beam_angle, 3.0);
    a_ball_ref = beam_angle * beam_angle;
    br = omega;
    br = cos(br);
    p_ball_ref = br * br;
    br = omega;
    br = cos(br);
    beam_angular_velocity = omega;
    beam_angular_velocity = cos(beam_angular_velocity);
    b_x = omega;
    b_x = sin(b_x);
    ball_velocity = ((obj_0->beam_len / 2.0 - dt) * (2.0 * obj_0->const_2) *
                     v_ball_ref * beam_angular_velocity * b_x + (obj_0->const_2 *
      amplitude * ball_velocity * u1 + obj_0->const_1 * beam_angle * br)) + 2.0 *
      obj_0->const_2 / obj_0->tau_val * (obj_0->beam_len / 2.0 - dt) *
      a_ball_ref * p_ball_ref;
    error[0] = error_x1;
    error[1] = error_x2;
    error[2] = error_x3;
    error[3] = ball_velocity;

    /*  This K is using the hard coded matrices from the function */
    /*  set below */
    state_derivative[0] = -obj->K[0];
    state_derivative[1] = -obj->K[1];
    state_derivative[2] = -obj->K[2];
    state_derivative[3] = -obj->K[3];
    error_x1 = state_derivative[0] * error[0];
    error_x1 += state_derivative[1] * error[1];
    error_x1 += state_derivative[2] * error[2];
    error_x1 += state_derivative[3] * error[3];

    /*  Apply feedback linearization to get the actual control */
    /*                  V_servo = obj.computeControl(estimated_ball_pos, estimated_ball_vel, estimated_beam_ang, estimated_beam_ang_vel, v); */
    obj_0 = obj;

    /*  Feedback linearization control law: u = (v - phi(x))/psi(x) */
    obj_1 = obj_0;

    /*  phi(x) part of the control law */
    ball_velocity = rt_powd_snf(beam_angle, 3.0);
    u1 = beam_angle * beam_angle;
    br = omega;
    br = cos(br);
    v_ball_ref = br * br;
    a_ball_ref = beam_angle * beam_angle;
    br = omega;
    br = cos(br);
    p_ball_ref = br * br;
    error_x2 = beam_angle * beam_angle;
    br = omega;
    br = cos(br);
    error_x3 = br * br;
    c = beam_angle * beam_angle;
    c_0 = rt_powd_snf(beam_angle, 3.0);
    br = omega;
    br = cos(br);
    c_1 = br * br;
    br = omega;
    br = sin(br);
    c_2 = br * br;
    c_3 = beam_angle * beam_angle;
    br = omega;
    br = cos(br);
    c_4 = br * br;
    c_5 = beam_angle * beam_angle;
    br = omega;
    br = cos(br);
    c_6 = br * br;
    br = omega;
    br = cos(br);
    beam_angular_velocity = omega;
    beam_angular_velocity = sin(beam_angular_velocity);
    b_x = omega;
    b_x = sin(b_x);
    b_x_0 = omega;
    b_x_0 = cos(b_x_0);
    b_x_1 = omega;
    b_x_1 = sin(b_x_1);
    b_x_2 = omega;
    b_x_2 = sin(b_x_2);
    b_x_3 = omega;
    b_x_3 = cos(b_x_3);
    b_x_4 = omega;
    b_x_4 = sin(b_x_4);
    b_x_5 = omega;
    b_x_5 = cos(b_x_5);
    b_x_6 = omega;
    b_x_6 = cos(b_x_6);
    b_x_7 = omega;
    b_x_7 = sin(b_x_7);
    p_ball_ref = ((((-2.0 * obj_1->const_2 * amplitude * c * b_x_0 * b_x_1 -
                     obj_1->const_1 * beam_angle * b_x_2) + (obj_1->beam_len /
      2.0 - dt) * (2.0 * obj_1->const_2) * c_0 * (c_1 - c_2)) - 4.0 *
                   obj_1->const_2 / obj_1->tau_val * (obj_1->beam_len / 2.0 - dt)
                   * c_3 * b_x_3 * b_x_4) * beam_angle + ((obj_1->const_1 * b_x
      - (obj_1->beam_len / 2.0 - dt) * obj_1->const_2 * error_x2 * error_x3) *
      (obj_1->const_2 * a_ball_ref * p_ball_ref) + (-2.0 * obj_1->const_2 *
      ball_velocity * br * beam_angular_velocity - 2.0 * obj_1->const_2 /
      obj_1->tau_val * u1 * v_ball_ref) * amplitude)) + (((obj_1->beam_len / 2.0
      - dt) * (6.0 * obj_1->const_2) * c_5 * b_x_6 * b_x_7 + (2.0 *
      obj_1->const_2 * amplitude * beam_angle * c_4 + obj_1->const_1 * b_x_5)) +
      4.0 * obj_1->const_2 / obj_1->tau_val * (obj_1->beam_len / 2.0 - dt) *
      beam_angle * c_6) * (-beam_angle / obj_1->tau_val);

    /*  psi(x) part of the control law */
    br = omega;
    br = cos(br);
    ball_velocity = br * br;
    u1 = beam_angle * beam_angle;
    br = omega;
    br = cos(br);
    v_ball_ref = br * br;
    br = omega;
    br = cos(br);
    beam_angular_velocity = omega;
    beam_angular_velocity = cos(beam_angular_velocity);
    b_x = omega;
    b_x = sin(b_x);
    ball_velocity = (((obj_0->beam_len / 2.0 - dt) * (6.0 * obj_0->const_2) * u1
                      * beam_angular_velocity * b_x + (2.0 * obj_0->const_2 *
      amplitude * beam_angle * ball_velocity + obj_0->const_1 * br)) + 4.0 *
                     obj_0->const_2 / obj_0->tau_val * (obj_0->beam_len / 2.0 -
      dt) * beam_angle * v_ball_ref) * (obj_0->servo_gain / obj_0->tau_val);
    ball_velocity = (error_x1 - p_ball_ref) / ball_velocity;
  } else {
    /*  this doesn't work rn */
    /*  LQR for system linearized around trajectory */
    /*  Calculate error state */
    error_x1 = dt - p_ball_ref;
    error_x2 = amplitude - v_ball_ref;
    error[0] = error_x1;
    error[1] = error_x2;
    error[2] = omega;
    error[3] = beam_angle;

    /*  Calculate Jacobian at current state */
    /*  Set LQR weights */
    /*  Position error weight */
    /*  Velocity error weight */
    /*  Angle weight */
    /*  Angular velocity weight */
    /*  Control effort weight */
    /*  Calculate LQR gain */
    /*  K = lqr(jacb, B_control, obj.Q, obj.R) */
    /*  Apply LQR to get control */
    ball_velocity = -20.0 * error[0];
    ball_velocity += -25.1525 * error[1];
    ball_velocity += -13.0233 * error[2];
    ball_velocity += -2.6315 * error[3];
  }

  /*  Apply safety limits to servo voltage */
  if (u2 > 0.78539816339744828) {
    u2 = (0.78539816339744828 - u2) * 10.0;
    if (!(ball_velocity <= u2)) {
      ball_velocity = u2;
    }
  } else if (u2 < -0.78539816339744828) {
    u2 = (-0.78539816339744828 - u2) * 10.0;
    if (!(ball_velocity >= u2)) {
      ball_velocity = u2;
    }
  }

  obj->V_servo = ball_velocity;
  obj->t_prev = u0;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_eval_type_B.MATLABSystem_o1 = ball_velocity;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_eval_type_B.MATLABSystem_o2 = dt;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_eval_type_B.MATLABSystem_o3 = amplitude;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_eval_type_B.MATLABSystem_o4 = omega;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_eval_type_B.MATLABSystem_o5 = beam_angle;

  /* Saturate: '<Root>/+//-10V' */
  u0 = simulink_experiment_eval_type_B.MATLABSystem_o1;
  u1 = simulink_experiment_eval_type_P.u0V_LowerSat;
  u2 = simulink_experiment_eval_type_P.u0V_UpperSat;
  if (u0 > u2) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_eval_type_B.u0V = u2;
  } else if (u0 < u1) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_eval_type_B.u0V = u1;
  } else {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_eval_type_B.u0V = u0;
  }

  /* End of Saturate: '<Root>/+//-10V' */

  /* Gain: '<S1>/Motor  Gain (V//V)' */
  simulink_experiment_eval_type_B.MotorGainVV =
    simulink_experiment_eval_type_P.MotorGainVV_Gain *
    simulink_experiment_eval_type_B.u0V;

  /* S-Function (hil_write_analog_block): '<S1>/HIL Write Analog' */

  /* S-Function Block: simulink_experiment_eval_type1/Ball and Beam Hardware Interface/HIL Write Analog (hil_write_analog_block) */
  {
    t_error result;
    result = hil_write_analog(simulink_experiment_eval_typ_DW.HILInitialize_Card,
      &simulink_experiment_eval_type_P.HILWriteAnalog_channels, 1,
      &simulink_experiment_eval_type_B.MotorGainVV);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
    }
  }

  /* RateTransition: '<Root>/Rate Transition' */
  if (simulink_experiment_eval_typ_M->Timing.RateInteraction.TID1_2) {
    simulink_experiment_eval_typ_DW.RateTransition_Buffer =
      simulink_experiment_eval_type_B.Clock;

    /* RateTransition: '<Root>/Rate Transition2' */
    simulink_experiment_eval_typ_DW.RateTransition2_Buffer =
      simulink_experiment_eval_type_B.MATLABSystem_o1;

    /* RateTransition: '<Root>/Rate Transition3' */
    simulink_experiment_eval_typ_DW.RateTransition3_Buffer =
      simulink_experiment_eval_type_B.BB01SensorGainmV;

    /* RateTransition: '<Root>/Rate Transition4' */
    simulink_experiment_eval_typ_DW.RateTransition4_Buffer =
      simulink_experiment_eval_type_B.Bias;
  }

  /* End of RateTransition: '<Root>/Rate Transition' */
}

/* Model update function for TID0 */
void simulink_experiment_eval_type1_update0(void) /* Sample time: [0.0s, 0.0s] */
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
  if (!(++simulink_experiment_eval_typ_M->Timing.clockTick0)) {
    ++simulink_experiment_eval_typ_M->Timing.clockTickH0;
  }

  simulink_experiment_eval_typ_M->Timing.t[0] =
    simulink_experiment_eval_typ_M->Timing.clockTick0 *
    simulink_experiment_eval_typ_M->Timing.stepSize0 +
    simulink_experiment_eval_typ_M->Timing.clockTickH0 *
    simulink_experiment_eval_typ_M->Timing.stepSize0 * 4294967296.0;

  /* Update absolute time */
  /* The "clockTick1" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick1"
   * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick1 and the high bits
   * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_eval_typ_M->Timing.clockTick1)) {
    ++simulink_experiment_eval_typ_M->Timing.clockTickH1;
  }

  simulink_experiment_eval_typ_M->Timing.t[1] =
    simulink_experiment_eval_typ_M->Timing.clockTick1 *
    simulink_experiment_eval_typ_M->Timing.stepSize1 +
    simulink_experiment_eval_typ_M->Timing.clockTickH1 *
    simulink_experiment_eval_typ_M->Timing.stepSize1 * 4294967296.0;
}

/* Model output function for TID2 */
void simulink_experiment_eval_type1_output2(void) /* Sample time: [0.01s, 0.0s] */
{
  /* RateTransition: '<Root>/Rate Transition2' */
  simulink_experiment_eval_type_B.RateTransition2 =
    simulink_experiment_eval_typ_DW.RateTransition2_Buffer;

  /* RateTransition: '<Root>/Rate Transition3' */
  simulink_experiment_eval_type_B.RateTransition3 =
    simulink_experiment_eval_typ_DW.RateTransition3_Buffer;

  /* RateTransition: '<Root>/Rate Transition4' */
  simulink_experiment_eval_type_B.RateTransition4 =
    simulink_experiment_eval_typ_DW.RateTransition4_Buffer;

  /* RateTransition: '<Root>/Rate Transition' */
  simulink_experiment_eval_type_B.RateTransition =
    simulink_experiment_eval_typ_DW.RateTransition_Buffer;
}

/* Model update function for TID2 */
void simulink_experiment_eval_type1_update2(void) /* Sample time: [0.01s, 0.0s] */
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
  if (!(++simulink_experiment_eval_typ_M->Timing.clockTick2)) {
    ++simulink_experiment_eval_typ_M->Timing.clockTickH2;
  }

  simulink_experiment_eval_typ_M->Timing.t[2] =
    simulink_experiment_eval_typ_M->Timing.clockTick2 *
    simulink_experiment_eval_typ_M->Timing.stepSize2 +
    simulink_experiment_eval_typ_M->Timing.clockTickH2 *
    simulink_experiment_eval_typ_M->Timing.stepSize2 * 4294967296.0;
}

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void simulink_experiment_eval_type1_output(int_T tid)
{
  switch (tid) {
   case 0 :
    simulink_experiment_eval_type1_output0();
    break;

   case 2 :
    simulink_experiment_eval_type1_output2();
    break;

   default :
    /* do nothing */
    break;
  }
}

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void simulink_experiment_eval_type1_update(int_T tid)
{
  switch (tid) {
   case 0 :
    simulink_experiment_eval_type1_update0();
    break;

   case 2 :
    simulink_experiment_eval_type1_update2();
    break;

   default :
    /* do nothing */
    break;
  }
}

/* Model initialize function */
void simulink_experiment_eval_type1_initialize(void)
{
  {
    studentControllerInterface_si_T *obj;

    /* Start for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

    /* S-Function Block: simulink_experiment_eval_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
    {
      t_int result;
      t_boolean is_switching;
      result = hil_open("q2_usb", "0",
                        &simulink_experiment_eval_typ_DW.HILInitialize_Card);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
        return;
      }

      is_switching = false;
      result = hil_set_card_specific_options
        (simulink_experiment_eval_typ_DW.HILInitialize_Card,
         "d0=digital;d1=digital;led=auto;update_rate=normal", 50);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
        return;
      }

      result = hil_watchdog_clear
        (simulink_experiment_eval_typ_DW.HILInitialize_Card);
      if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
        return;
      }

      if ((simulink_experiment_eval_type_P.HILInitialize_AIPStart &&
           !is_switching) ||
          (simulink_experiment_eval_type_P.HILInitialize_AIPEnter &&
           is_switching)) {
        simulink_experiment_eval_typ_DW.HILInitialize_AIMinimums[0] =
          (simulink_experiment_eval_type_P.HILInitialize_AILow);
        simulink_experiment_eval_typ_DW.HILInitialize_AIMinimums[1] =
          (simulink_experiment_eval_type_P.HILInitialize_AILow);
        simulink_experiment_eval_typ_DW.HILInitialize_AIMaximums[0] =
          simulink_experiment_eval_type_P.HILInitialize_AIHigh;
        simulink_experiment_eval_typ_DW.HILInitialize_AIMaximums[1] =
          simulink_experiment_eval_type_P.HILInitialize_AIHigh;
        result = hil_set_analog_input_ranges
          (simulink_experiment_eval_typ_DW.HILInitialize_Card,
           simulink_experiment_eval_type_P.HILInitialize_AIChannels, 2U,
           &simulink_experiment_eval_typ_DW.HILInitialize_AIMinimums[0],
           &simulink_experiment_eval_typ_DW.HILInitialize_AIMaximums[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_eval_type_P.HILInitialize_AOPStart &&
           !is_switching) ||
          (simulink_experiment_eval_type_P.HILInitialize_AOPEnter &&
           is_switching)) {
        simulink_experiment_eval_typ_DW.HILInitialize_AOMinimums[0] =
          (simulink_experiment_eval_type_P.HILInitialize_AOLow);
        simulink_experiment_eval_typ_DW.HILInitialize_AOMinimums[1] =
          (simulink_experiment_eval_type_P.HILInitialize_AOLow);
        simulink_experiment_eval_typ_DW.HILInitialize_AOMaximums[0] =
          simulink_experiment_eval_type_P.HILInitialize_AOHigh;
        simulink_experiment_eval_typ_DW.HILInitialize_AOMaximums[1] =
          simulink_experiment_eval_type_P.HILInitialize_AOHigh;
        result = hil_set_analog_output_ranges
          (simulink_experiment_eval_typ_DW.HILInitialize_Card,
           simulink_experiment_eval_type_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_eval_typ_DW.HILInitialize_AOMinimums[0],
           &simulink_experiment_eval_typ_DW.HILInitialize_AOMaximums[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_eval_type_P.HILInitialize_AOStart &&
           !is_switching) ||
          (simulink_experiment_eval_type_P.HILInitialize_AOEnter && is_switching))
      {
        simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[0] =
          simulink_experiment_eval_type_P.HILInitialize_AOInitial;
        simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[1] =
          simulink_experiment_eval_type_P.HILInitialize_AOInitial;
        result = hil_write_analog
          (simulink_experiment_eval_typ_DW.HILInitialize_Card,
           simulink_experiment_eval_type_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
          return;
        }
      }

      if (simulink_experiment_eval_type_P.HILInitialize_AOReset) {
        simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[0] =
          simulink_experiment_eval_type_P.HILInitialize_AOWatchdog;
        simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[1] =
          simulink_experiment_eval_type_P.HILInitialize_AOWatchdog;
        result = hil_watchdog_set_analog_expiration_state
          (simulink_experiment_eval_typ_DW.HILInitialize_Card,
           simulink_experiment_eval_type_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
          return;
        }
      }

      result = hil_set_digital_directions
        (simulink_experiment_eval_typ_DW.HILInitialize_Card, NULL, 0U,
         simulink_experiment_eval_type_P.HILInitialize_DOChannels, 8U);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
        return;
      }

      if ((simulink_experiment_eval_type_P.HILInitialize_DOStart &&
           !is_switching) ||
          (simulink_experiment_eval_type_P.HILInitialize_DOEnter && is_switching))
      {
        {
          int_T i1;
          boolean_T *dw_DOBits =
            &simulink_experiment_eval_typ_DW.HILInitialize_DOBits[0];
          for (i1=0; i1 < 8; i1++) {
            dw_DOBits[i1] =
              simulink_experiment_eval_type_P.HILInitialize_DOInitial;
          }
        }

        result = hil_write_digital
          (simulink_experiment_eval_typ_DW.HILInitialize_Card,
           simulink_experiment_eval_type_P.HILInitialize_DOChannels, 8U,
           (t_boolean *) &simulink_experiment_eval_typ_DW.HILInitialize_DOBits[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
          return;
        }
      }

      if (simulink_experiment_eval_type_P.HILInitialize_DOReset) {
        {
          int_T i1;
          int32_T *dw_DOStates =
            &simulink_experiment_eval_typ_DW.HILInitialize_DOStates[0];
          for (i1=0; i1 < 8; i1++) {
            dw_DOStates[i1] =
              simulink_experiment_eval_type_P.HILInitialize_DOWatchdog;
          }
        }

        result = hil_watchdog_set_digital_expiration_state
          (simulink_experiment_eval_typ_DW.HILInitialize_Card,
           simulink_experiment_eval_type_P.HILInitialize_DOChannels, 8U, (const
            t_digital_state *)
           &simulink_experiment_eval_typ_DW.HILInitialize_DOStates[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_eval_type_P.HILInitialize_EIPStart &&
           !is_switching) ||
          (simulink_experiment_eval_type_P.HILInitialize_EIPEnter &&
           is_switching)) {
        simulink_experiment_eval_typ_DW.HILInitialize_QuadratureModes[0] =
          simulink_experiment_eval_type_P.HILInitialize_EIQuadrature;
        simulink_experiment_eval_typ_DW.HILInitialize_QuadratureModes[1] =
          simulink_experiment_eval_type_P.HILInitialize_EIQuadrature;
        result = hil_set_encoder_quadrature_mode
          (simulink_experiment_eval_typ_DW.HILInitialize_Card,
           simulink_experiment_eval_type_P.HILInitialize_EIChannels, 2U,
           (t_encoder_quadrature_mode *)
           &simulink_experiment_eval_typ_DW.HILInitialize_QuadratureModes[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_eval_type_P.HILInitialize_EIStart &&
           !is_switching) ||
          (simulink_experiment_eval_type_P.HILInitialize_EIEnter && is_switching))
      {
        simulink_experiment_eval_typ_DW.HILInitialize_InitialEICounts[0] =
          simulink_experiment_eval_type_P.HILInitialize_EIInitial;
        simulink_experiment_eval_typ_DW.HILInitialize_InitialEICounts[1] =
          simulink_experiment_eval_type_P.HILInitialize_EIInitial;
        result = hil_set_encoder_counts
          (simulink_experiment_eval_typ_DW.HILInitialize_Card,
           simulink_experiment_eval_type_P.HILInitialize_EIChannels, 2U,
           &simulink_experiment_eval_typ_DW.HILInitialize_InitialEICounts[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
          return;
        }
      }
    }

    /* Start for S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

    /* S-Function Block: simulink_experiment_eval_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_create_encoder_reader
        (simulink_experiment_eval_typ_DW.HILInitialize_Card,
         simulink_experiment_eval_type_P.HILReadEncoderTimebase_SamplesI,
         &simulink_experiment_eval_type_P.HILReadEncoderTimebase_Channels, 1,
         &simulink_experiment_eval_typ_DW.HILReadEncoderTimebase_Task);
      if (result >= 0) {
        result = hil_task_set_buffer_overflow_mode
          (simulink_experiment_eval_typ_DW.HILReadEncoderTimebase_Task,
           (t_buffer_overflow_mode)
           (simulink_experiment_eval_type_P.HILReadEncoderTimebase_Overflow - 1));
      }

      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
      }
    }

    /* Start for MATLABSystem: '<Root>/MATLAB System' */
    studentControllerInterface_stud(&simulink_experiment_eval_typ_DW.obj);
    simulink_experiment_eval_typ_DW.objisempty = true;
    obj = &simulink_experiment_eval_typ_DW.obj;
    obj->isInitialized = 1;
  }
}

/* Model terminate function */
void simulink_experiment_eval_type1_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

  /* S-Function Block: simulink_experiment_eval_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_digital_outputs = 0;
    hil_task_stop_all(simulink_experiment_eval_typ_DW.HILInitialize_Card);
    hil_monitor_stop_all(simulink_experiment_eval_typ_DW.HILInitialize_Card);
    is_switching = false;
    if ((simulink_experiment_eval_type_P.HILInitialize_AOTerminate &&
         !is_switching) || (simulink_experiment_eval_type_P.HILInitialize_AOExit
         && is_switching)) {
      simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[0] =
        simulink_experiment_eval_type_P.HILInitialize_AOFinal;
      simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[1] =
        simulink_experiment_eval_type_P.HILInitialize_AOFinal;
      num_final_analog_outputs = 2U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((simulink_experiment_eval_type_P.HILInitialize_DOTerminate &&
         !is_switching) || (simulink_experiment_eval_type_P.HILInitialize_DOExit
         && is_switching)) {
      {
        int_T i1;
        boolean_T *dw_DOBits =
          &simulink_experiment_eval_typ_DW.HILInitialize_DOBits[0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOBits[i1] = simulink_experiment_eval_type_P.HILInitialize_DOFinal;
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
      result = hil_write(simulink_experiment_eval_typ_DW.HILInitialize_Card
                         ,
                         simulink_experiment_eval_type_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , NULL, 0
                         ,
                         simulink_experiment_eval_type_P.HILInitialize_DOChannels,
                         num_final_digital_outputs
                         , NULL, 0
                         ,
                         &simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages
                         [0]
                         , NULL
                         , (t_boolean *)
                         &simulink_experiment_eval_typ_DW.HILInitialize_DOBits[0]
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (simulink_experiment_eval_typ_DW.HILInitialize_Card,
             simulink_experiment_eval_type_P.HILInitialize_AOChannels,
             num_final_analog_outputs,
             &simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_digital_outputs > 0) {
          local_result = hil_write_digital
            (simulink_experiment_eval_typ_DW.HILInitialize_Card,
             simulink_experiment_eval_type_P.HILInitialize_DOChannels,
             num_final_digital_outputs, (t_boolean *)
             &simulink_experiment_eval_typ_DW.HILInitialize_DOBits[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_eval_typ_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(simulink_experiment_eval_typ_DW.HILInitialize_Card);
    hil_monitor_delete_all(simulink_experiment_eval_typ_DW.HILInitialize_Card);
    hil_close(simulink_experiment_eval_typ_DW.HILInitialize_Card);
    simulink_experiment_eval_typ_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  if (tid == 1)
    tid = 0;
  simulink_experiment_eval_type1_output(tid);
}

void MdlUpdate(int_T tid)
{
  if (tid == 1)
    tid = 0;
  simulink_experiment_eval_type1_update(tid);
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
  simulink_experiment_eval_type1_initialize();
}

void MdlTerminate(void)
{
  simulink_experiment_eval_type1_terminate();
}

/* Registration function */
RT_MODEL_simulink_experiment__T *simulink_experiment_eval_type1(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)simulink_experiment_eval_typ_M, 0,
                sizeof(RT_MODEL_simulink_experiment__T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&simulink_experiment_eval_typ_M->solverInfo,
                          &simulink_experiment_eval_typ_M->Timing.simTimeStep);
    rtsiSetTPtr(&simulink_experiment_eval_typ_M->solverInfo, &rtmGetTPtr
                (simulink_experiment_eval_typ_M));
    rtsiSetStepSizePtr(&simulink_experiment_eval_typ_M->solverInfo,
                       &simulink_experiment_eval_typ_M->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&simulink_experiment_eval_typ_M->solverInfo,
                          (&rtmGetErrorStatus(simulink_experiment_eval_typ_M)));
    rtsiSetRTModelPtr(&simulink_experiment_eval_typ_M->solverInfo,
                      simulink_experiment_eval_typ_M);
  }

  rtsiSetSimTimeStep(&simulink_experiment_eval_typ_M->solverInfo,
                     MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange
    (&simulink_experiment_eval_typ_M->solverInfo, false);
  rtsiSetSolverName(&simulink_experiment_eval_typ_M->solverInfo,
                    "FixedStepDiscrete");

  /* Initialize timing info */
  {
    int_T *mdlTsMap =
      simulink_experiment_eval_typ_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    mdlTsMap[2] = 2;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "simulink_experiment_eval_typ_M points to
       static memory which is guaranteed to be non-NULL" */
    simulink_experiment_eval_typ_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    simulink_experiment_eval_typ_M->Timing.sampleTimes =
      (&simulink_experiment_eval_typ_M->Timing.sampleTimesArray[0]);
    simulink_experiment_eval_typ_M->Timing.offsetTimes =
      (&simulink_experiment_eval_typ_M->Timing.offsetTimesArray[0]);

    /* task periods */
    simulink_experiment_eval_typ_M->Timing.sampleTimes[0] = (0.0);
    simulink_experiment_eval_typ_M->Timing.sampleTimes[1] = (0.002);
    simulink_experiment_eval_typ_M->Timing.sampleTimes[2] = (0.01);

    /* task offsets */
    simulink_experiment_eval_typ_M->Timing.offsetTimes[0] = (0.0);
    simulink_experiment_eval_typ_M->Timing.offsetTimes[1] = (0.0);
    simulink_experiment_eval_typ_M->Timing.offsetTimes[2] = (0.0);
  }

  rtmSetTPtr(simulink_experiment_eval_typ_M,
             &simulink_experiment_eval_typ_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = simulink_experiment_eval_typ_M->Timing.sampleHitArray;
    int_T *mdlPerTaskSampleHits =
      simulink_experiment_eval_typ_M->Timing.perTaskSampleHitsArray;
    simulink_experiment_eval_typ_M->Timing.perTaskSampleHits =
      (&mdlPerTaskSampleHits[0]);
    mdlSampleHits[0] = 1;
    simulink_experiment_eval_typ_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(simulink_experiment_eval_typ_M, 90.0);
  simulink_experiment_eval_typ_M->Timing.stepSize0 = 0.002;
  simulink_experiment_eval_typ_M->Timing.stepSize1 = 0.002;
  simulink_experiment_eval_typ_M->Timing.stepSize2 = 0.01;

  /* External mode info */
  simulink_experiment_eval_typ_M->Sizes.checksums[0] = (1603816804U);
  simulink_experiment_eval_typ_M->Sizes.checksums[1] = (34713326U);
  simulink_experiment_eval_typ_M->Sizes.checksums[2] = (2992390849U);
  simulink_experiment_eval_typ_M->Sizes.checksums[3] = (505437970U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[2];
    simulink_experiment_eval_typ_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(simulink_experiment_eval_typ_M->extModeInfo,
      &simulink_experiment_eval_typ_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(simulink_experiment_eval_typ_M->extModeInfo,
                        simulink_experiment_eval_typ_M->Sizes.checksums);
    rteiSetTPtr(simulink_experiment_eval_typ_M->extModeInfo, rtmGetTPtr
                (simulink_experiment_eval_typ_M));
  }

  simulink_experiment_eval_typ_M->solverInfoPtr =
    (&simulink_experiment_eval_typ_M->solverInfo);
  simulink_experiment_eval_typ_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&simulink_experiment_eval_typ_M->solverInfo, 0.002);
  rtsiSetSolverMode(&simulink_experiment_eval_typ_M->solverInfo,
                    SOLVER_MODE_MULTITASKING);

  /* block I/O */
  simulink_experiment_eval_typ_M->blockIO = ((void *)
    &simulink_experiment_eval_type_B);

  {
    simulink_experiment_eval_type_B.HILReadEncoderTimebase = 0.0;
    simulink_experiment_eval_type_B.HILReadAnalog = 0.0;
    simulink_experiment_eval_type_B.BB01SensorGainmV = 0.0;
    simulink_experiment_eval_type_B.EncoderCalibrationradcount = 0.0;
    simulink_experiment_eval_type_B.Bias = 0.0;
    simulink_experiment_eval_type_B.Clock = 0.0;
    simulink_experiment_eval_type_B.u0V = 0.0;
    simulink_experiment_eval_type_B.MotorGainVV = 0.0;
    simulink_experiment_eval_type_B.RateTransition2 = 0.0;
    simulink_experiment_eval_type_B.RateTransition3 = 0.0;
    simulink_experiment_eval_type_B.RateTransition4 = 0.0;
    simulink_experiment_eval_type_B.RateTransition = 0.0;
    simulink_experiment_eval_type_B.MATLABSystem_o1 = 0.0;
    simulink_experiment_eval_type_B.MATLABSystem_o2 = 0.0;
    simulink_experiment_eval_type_B.MATLABSystem_o3 = 0.0;
    simulink_experiment_eval_type_B.MATLABSystem_o4 = 0.0;
    simulink_experiment_eval_type_B.MATLABSystem_o5 = 0.0;
  }

  /* parameters */
  simulink_experiment_eval_typ_M->defaultParam = ((real_T *)
    &simulink_experiment_eval_type_P);

  /* states (dwork) */
  simulink_experiment_eval_typ_M->dwork = ((void *)
    &simulink_experiment_eval_typ_DW);
  (void) memset((void *)&simulink_experiment_eval_typ_DW, 0,
                sizeof(DW_simulink_experiment_eval_t_T));
  simulink_experiment_eval_typ_DW.HILInitialize_AIMinimums[0] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_AIMinimums[1] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_AIMaximums[0] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_AIMaximums[1] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_AOMinimums[0] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_AOMinimums[1] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_AOMaximums[0] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_AOMaximums[1] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[0] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_AOVoltages[1] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_FilterFrequency[0] = 0.0;
  simulink_experiment_eval_typ_DW.HILInitialize_FilterFrequency[1] = 0.0;
  simulink_experiment_eval_typ_DW.HILReadAnalog_Buffer = 0.0;
  simulink_experiment_eval_typ_DW.RateTransition_Buffer = 0.0;
  simulink_experiment_eval_typ_DW.RateTransition2_Buffer = 0.0;
  simulink_experiment_eval_typ_DW.RateTransition3_Buffer = 0.0;
  simulink_experiment_eval_typ_DW.RateTransition4_Buffer = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    simulink_experiment_eval_typ_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 22;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  simulink_experiment_eval_typ_M->Sizes.numContStates = (0);/* Number of continuous states */
  simulink_experiment_eval_typ_M->Sizes.numY = (0);/* Number of model outputs */
  simulink_experiment_eval_typ_M->Sizes.numU = (0);/* Number of model inputs */
  simulink_experiment_eval_typ_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  simulink_experiment_eval_typ_M->Sizes.numSampTimes = (3);/* Number of sample times */
  simulink_experiment_eval_typ_M->Sizes.numBlocks = (20);/* Number of blocks */
  simulink_experiment_eval_typ_M->Sizes.numBlockIO = (17);/* Number of block outputs */
  simulink_experiment_eval_typ_M->Sizes.numBlockPrms = (86);/* Sum of parameter "widths" */
  return simulink_experiment_eval_typ_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
