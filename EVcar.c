/*
 * TWOWheeled_EV.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "TWOWheeled_EV".
 *
 * Model version              : 4.4
 * Simulink Coder version : 23.2 (R2023b) 01-Aug-2023
 * C source code generated on : Sat Jan 13 13:21:15 2024
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Atmel->AVR
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "TWOWheeled_EV.h"
#include "rtwtypes.h"
#include "TWOWheeled_EV_private.h"
#include <string.h>
#include <stddef.h>
#include "rt_nonfinite.h"
#include <math.h>
#include <float.h>

/* Block signals (default storage) */
B_TWOWheeled_EV_T TWOWheeled_EV_B;

/* Continuous states */
X_TWOWheeled_EV_T TWOWheeled_EV_X;

/* Disabled State Vector */
XDis_TWOWheeled_EV_T TWOWheeled_EV_XDis;

/* Block states (default storage) */
DW_TWOWheeled_EV_T TWOWheeled_EV_DW;

/* Mass Matrices */
MassMatrix_TWOWheeled_EV_T TWOWheeled_EV_MassMatrix;

/* Real-time model */
static RT_MODEL_TWOWheeled_EV_T TWOWheeled_EV_M_;
RT_MODEL_TWOWheeled_EV_T *const TWOWheeled_EV_M = &TWOWheeled_EV_M_;
real_T look1_pbinlcapw(real_T u0, const real_T bp0[], const real_T table[],
  uint32_T prevIndex[], uint32_T maxIndex)
{
  real_T frac;
  real_T y;
  real_T yL_0d0;
  uint32_T bpIdx;
  uint32_T found;
  uint32_T iLeft;
  uint32_T iRght;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'on'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Clip'
     Use last breakpoint for index at or above upper limit: 'on'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'on'
     Use last breakpoint for index at or above upper limit: 'on'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0UL]) {
    bpIdx = 0UL;
    frac = 0.0;
  } else if (u0 < bp0[maxIndex]) {
    /* Binary Search using Previous Index */
    bpIdx = prevIndex[0UL];
    iLeft = 0UL;
    iRght = maxIndex;
    found = 0UL;
    while (found == 0UL) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx - 1UL;
        bpIdx = ((bpIdx + iLeft) - 1UL) >> 1UL;
      } else if (u0 < bp0[bpIdx + 1UL]) {
        found = 1UL;
      } else {
        iLeft = bpIdx + 1UL;
        bpIdx = ((bpIdx + iRght) + 1UL) >> 1UL;
      }
    }

    frac = (u0 - bp0[bpIdx]) / (bp0[bpIdx + 1UL] - bp0[bpIdx]);
  } else {
    bpIdx = maxIndex;
    frac = 0.0;
  }

  prevIndex[0UL] = bpIdx;

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'on'
     Overflow mode: 'portable wrapping'
   */
  if (bpIdx == maxIndex) {
    y = table[bpIdx];
  } else {
    yL_0d0 = table[bpIdx];
    y = (table[bpIdx + 1UL] - yL_0d0) * frac + yL_0d0;
  }

  return y;
}

/* Projection for root system: '<Root>' */
void TWOWheeled_EV_projection(void)
{
  NeslSimulationData *simulationData;
  NeuDiagnosticManager *diagnosticManager;
  NeuDiagnosticTree *diagnosticTree;
  char *msg;
  real_T tmp_0[16];
  real_T time;
  int32_T tmp_2;
  int_T tmp_1[5];
  boolean_T tmp;

  /* Projection for SimscapeExecutionBlock: '<S42>/STATE_1' */
  simulationData = (NeslSimulationData *)TWOWheeled_EV_DW.STATE_1_SimData;
  time = TWOWheeled_EV_M->Timing.t[0];
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time;
  simulationData->mData->mContStates.mN = 27;
  simulationData->mData->mContStates.mX =
    &TWOWheeled_EV_X.TWOWheeled_EVBatterycharge[0];
  simulationData->mData->mDiscStates.mN = 0;
  simulationData->mData->mDiscStates.mX = &TWOWheeled_EV_DW.STATE_1_Discrete;
  simulationData->mData->mModeVector.mN = 6;
  simulationData->mData->mModeVector.mX = (int32_T *)
    &TWOWheeled_EV_DW.STATE_1_Modes[0];
  tmp = false;
  simulationData->mData->mFoundZcEvents = tmp;
  simulationData->mData->mIsMajorTimeStep = rtmIsMajorTimeStep(TWOWheeled_EV_M);
  tmp = false;
  simulationData->mData->mIsSolverAssertCheck = tmp;
  simulationData->mData->mIsSolverCheckingCIC = false;
  tmp = rtsiIsSolverComputingJacobian(&TWOWheeled_EV_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = tmp;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  simulationData->mData->mIsModeUpdateTimeStep = rtsiIsModeUpdateTimeStep
    (&TWOWheeled_EV_M->solverInfo);
  tmp_1[0] = 0;
  tmp_0[0] = TWOWheeled_EV_B.INPUT_1_1_1[0];
  tmp_0[1] = TWOWheeled_EV_B.INPUT_1_1_1[1];
  tmp_0[2] = TWOWheeled_EV_B.INPUT_1_1_1[2];
  tmp_0[3] = TWOWheeled_EV_B.INPUT_1_1_1[3];
  tmp_1[1] = 4;
  tmp_0[4] = TWOWheeled_EV_B.INPUT_2_1_1[0];
  tmp_0[5] = TWOWheeled_EV_B.INPUT_2_1_1[1];
  tmp_0[6] = TWOWheeled_EV_B.INPUT_2_1_1[2];
  tmp_0[7] = TWOWheeled_EV_B.INPUT_2_1_1[3];
  tmp_1[2] = 8;
  tmp_0[8] = TWOWheeled_EV_B.INPUT_4_1_1[0];
  tmp_0[9] = TWOWheeled_EV_B.INPUT_4_1_1[1];
  tmp_0[10] = TWOWheeled_EV_B.INPUT_4_1_1[2];
  tmp_0[11] = TWOWheeled_EV_B.INPUT_4_1_1[3];
  tmp_1[3] = 12;
  tmp_0[12] = TWOWheeled_EV_B.INPUT_3_1_1[0];
  tmp_0[13] = TWOWheeled_EV_B.INPUT_3_1_1[1];
  tmp_0[14] = TWOWheeled_EV_B.INPUT_3_1_1[2];
  tmp_0[15] = TWOWheeled_EV_B.INPUT_3_1_1[3];
  tmp_1[4] = 16;
  simulationData->mData->mInputValues.mN = 16;
  simulationData->mData->mInputValues.mX = &tmp_0[0];
  simulationData->mData->mInputOffsets.mN = 5;
  simulationData->mData->mInputOffsets.mX = (int32_T *)&tmp_1[0];
  diagnosticManager = (NeuDiagnosticManager *)TWOWheeled_EV_DW.STATE_1_DiagMgr;
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  tmp_2 = ne_simulator_method((NeslSimulator *)
    TWOWheeled_EV_DW.STATE_1_Simulator, NESL_SIM_PROJECTION, simulationData,
    diagnosticManager);
  if (tmp_2 != 0L) {
    tmp = error_buffer_is_empty(rtmGetErrorStatus(TWOWheeled_EV_M));
    if (tmp) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(TWOWheeled_EV_M, msg);
    }
  }

  /* End of Projection for SimscapeExecutionBlock: '<S42>/STATE_1' */
}

/* ForcingFunction for root system: '<Root>' */
void TWOWheeled_EV_forcingfunction(void)
{
  NeslSimulationData *simulationData;
  NeuDiagnosticManager *diagnosticManager;
  NeuDiagnosticTree *diagnosticTree;
  XDot_TWOWheeled_EV_T *_rtXdot;
  char *msg;
  real_T tmp[16];
  real_T time;
  int32_T tmp_1;
  int_T tmp_0[5];
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot_TWOWheeled_EV_T *) TWOWheeled_EV_M->derivs);

  /* ForcingFunction for Integrator: '<S34>/Integrator' */
  _rtXdot->Integrator_CSTATE = TWOWheeled_EV_B.Product_d;

  /* ForcingFunction for Integrator: '<S26>/Integrator1' */
  lsat = (TWOWheeled_EV_X.Integrator1_CSTATE <=
          TWOWheeled_EV_P.Integrator1_LowerSat);
  usat = (TWOWheeled_EV_X.Integrator1_CSTATE >=
          TWOWheeled_EV_P.Integrator1_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (TWOWheeled_EV_B.Sum8 > 0.0)) || (usat &&
       (TWOWheeled_EV_B.Sum8 < 0.0))) {
    _rtXdot->Integrator1_CSTATE = TWOWheeled_EV_B.Sum8;
  } else {
    /* in saturation */
    _rtXdot->Integrator1_CSTATE = 0.0;
  }

  /* End of ForcingFunction for Integrator: '<S26>/Integrator1' */

  /* ForcingFunction for SimscapeExecutionBlock: '<S42>/STATE_1' */
  simulationData = (NeslSimulationData *)TWOWheeled_EV_DW.STATE_1_SimData;
  time = TWOWheeled_EV_M->Timing.t[0];
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time;
  simulationData->mData->mContStates.mN = 27;
  simulationData->mData->mContStates.mX =
    &TWOWheeled_EV_X.TWOWheeled_EVBatterycharge[0];
  simulationData->mData->mDiscStates.mN = 0;
  simulationData->mData->mDiscStates.mX = &TWOWheeled_EV_DW.STATE_1_Discrete;
  simulationData->mData->mModeVector.mN = 6;
  simulationData->mData->mModeVector.mX = (int32_T *)
    &TWOWheeled_EV_DW.STATE_1_Modes[0];
  lsat = false;
  simulationData->mData->mFoundZcEvents = lsat;
  simulationData->mData->mIsMajorTimeStep = rtmIsMajorTimeStep(TWOWheeled_EV_M);
  lsat = false;
  simulationData->mData->mIsSolverAssertCheck = lsat;
  simulationData->mData->mIsSolverCheckingCIC = false;
  lsat = rtsiIsSolverComputingJacobian(&TWOWheeled_EV_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = lsat;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  simulationData->mData->mIsModeUpdateTimeStep = rtsiIsModeUpdateTimeStep
    (&TWOWheeled_EV_M->solverInfo);
  tmp_0[0] = 0;
  tmp[0] = TWOWheeled_EV_B.INPUT_1_1_1[0];
  tmp[1] = TWOWheeled_EV_B.INPUT_1_1_1[1];
  tmp[2] = TWOWheeled_EV_B.INPUT_1_1_1[2];
  tmp[3] = TWOWheeled_EV_B.INPUT_1_1_1[3];
  tmp_0[1] = 4;
  tmp[4] = TWOWheeled_EV_B.INPUT_2_1_1[0];
  tmp[5] = TWOWheeled_EV_B.INPUT_2_1_1[1];
  tmp[6] = TWOWheeled_EV_B.INPUT_2_1_1[2];
  tmp[7] = TWOWheeled_EV_B.INPUT_2_1_1[3];
  tmp_0[2] = 8;
  tmp[8] = TWOWheeled_EV_B.INPUT_4_1_1[0];
  tmp[9] = TWOWheeled_EV_B.INPUT_4_1_1[1];
  tmp[10] = TWOWheeled_EV_B.INPUT_4_1_1[2];
  tmp[11] = TWOWheeled_EV_B.INPUT_4_1_1[3];
  tmp_0[3] = 12;
  tmp[12] = TWOWheeled_EV_B.INPUT_3_1_1[0];
  tmp[13] = TWOWheeled_EV_B.INPUT_3_1_1[1];
  tmp[14] = TWOWheeled_EV_B.INPUT_3_1_1[2];
  tmp[15] = TWOWheeled_EV_B.INPUT_3_1_1[3];
  tmp_0[4] = 16;
  simulationData->mData->mInputValues.mN = 16;
  simulationData->mData->mInputValues.mX = &tmp[0];
  simulationData->mData->mInputOffsets.mN = 5;
  simulationData->mData->mInputOffsets.mX = (int32_T *)&tmp_0[0];
  simulationData->mData->mDx.mN = 27;
  simulationData->mData->mDx.mX = &_rtXdot->TWOWheeled_EVBatterycharge[0];
  diagnosticManager = (NeuDiagnosticManager *)TWOWheeled_EV_DW.STATE_1_DiagMgr;
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  tmp_1 = ne_simulator_method((NeslSimulator *)
    TWOWheeled_EV_DW.STATE_1_Simulator, NESL_SIM_FORCINGFUNCTION, simulationData,
    diagnosticManager);
  if (tmp_1 != 0L) {
    lsat = error_buffer_is_empty(rtmGetErrorStatus(TWOWheeled_EV_M));
    if (lsat) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(TWOWheeled_EV_M, msg);
    }
  }

  /* End of ForcingFunction for SimscapeExecutionBlock: '<S42>/STATE_1' */

  /* ForcingFunction for Integrator: '<Root>/Integrator' */
  _rtXdot->Integrator_CSTATE_f = TWOWheeled_EV_B.OUTPUT_1_0[2];

  /* ForcingFunction for Integrator: '<S35>/Integrator2' */
  _rtXdot->Integrator2_CSTATE = TWOWheeled_EV_B.Product;
}

/* MassMatrix for root system: '<Root>' */
void TWOWheeled_EV_massmatrix(void)
{
  NeslSimulationData *simulationData;
  NeuDiagnosticManager *diagnosticManager;
  NeuDiagnosticTree *diagnosticTree;
  char *msg;
  real_T tmp_0[16];
  real_T time;
  real_T *tmp_2;
  real_T *tmp_3;
  int32_T tmp_4;
  int_T tmp_1[5];
  boolean_T tmp;

  /* MassMatrix for SimscapeExecutionBlock: '<S42>/STATE_1' */
  simulationData = (NeslSimulationData *)TWOWheeled_EV_DW.STATE_1_SimData;
  time = TWOWheeled_EV_M->Timing.t[0];
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time;
  simulationData->mData->mContStates.mN = 27;
  simulationData->mData->mContStates.mX =
    &TWOWheeled_EV_X.TWOWheeled_EVBatterycharge[0];
  simulationData->mData->mDiscStates.mN = 0;
  simulationData->mData->mDiscStates.mX = &TWOWheeled_EV_DW.STATE_1_Discrete;
  simulationData->mData->mModeVector.mN = 6;
  simulationData->mData->mModeVector.mX = (int32_T *)
    &TWOWheeled_EV_DW.STATE_1_Modes[0];
  tmp = false;
  simulationData->mData->mFoundZcEvents = tmp;
  simulationData->mData->mIsMajorTimeStep = rtmIsMajorTimeStep(TWOWheeled_EV_M);
  tmp = false;
  simulationData->mData->mIsSolverAssertCheck = tmp;
  simulationData->mData->mIsSolverCheckingCIC = false;
  tmp = rtsiIsSolverComputingJacobian(&TWOWheeled_EV_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = tmp;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  simulationData->mData->mIsModeUpdateTimeStep = rtsiIsModeUpdateTimeStep
    (&TWOWheeled_EV_M->solverInfo);
  tmp_1[0] = 0;
  tmp_0[0] = TWOWheeled_EV_B.INPUT_1_1_1[0];
  tmp_0[1] = TWOWheeled_EV_B.INPUT_1_1_1[1];
  tmp_0[2] = TWOWheeled_EV_B.INPUT_1_1_1[2];
  tmp_0[3] = TWOWheeled_EV_B.INPUT_1_1_1[3];
  tmp_1[1] = 4;
  tmp_0[4] = TWOWheeled_EV_B.INPUT_2_1_1[0];
  tmp_0[5] = TWOWheeled_EV_B.INPUT_2_1_1[1];
  tmp_0[6] = TWOWheeled_EV_B.INPUT_2_1_1[2];
  tmp_0[7] = TWOWheeled_EV_B.INPUT_2_1_1[3];
  tmp_1[2] = 8;
  tmp_0[8] = TWOWheeled_EV_B.INPUT_4_1_1[0];
  tmp_0[9] = TWOWheeled_EV_B.INPUT_4_1_1[1];
  tmp_0[10] = TWOWheeled_EV_B.INPUT_4_1_1[2];
  tmp_0[11] = TWOWheeled_EV_B.INPUT_4_1_1[3];
  tmp_1[3] = 12;
  tmp_0[12] = TWOWheeled_EV_B.INPUT_3_1_1[0];
  tmp_0[13] = TWOWheeled_EV_B.INPUT_3_1_1[1];
  tmp_0[14] = TWOWheeled_EV_B.INPUT_3_1_1[2];
  tmp_0[15] = TWOWheeled_EV_B.INPUT_3_1_1[3];
  tmp_1[4] = 16;
  simulationData->mData->mInputValues.mN = 16;
  simulationData->mData->mInputValues.mX = &tmp_0[0];
  simulationData->mData->mInputOffsets.mN = 5;
  simulationData->mData->mInputOffsets.mX = (int32_T *)&tmp_1[0];
  tmp_2 = TWOWheeled_EV_M->massMatrixPr;
  tmp_3 = double_pointer_shift(tmp_2, TWOWheeled_EV_DW.STATE_1_MASS_MATRIX_PR);
  simulationData->mData->mMassMatrixPr.mN = 8;
  simulationData->mData->mMassMatrixPr.mX = tmp_3;
  diagnosticManager = (NeuDiagnosticManager *)TWOWheeled_EV_DW.STATE_1_DiagMgr;
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  tmp_4 = ne_simulator_method((NeslSimulator *)
    TWOWheeled_EV_DW.STATE_1_Simulator, NESL_SIM_MASSMATRIX, simulationData,
    diagnosticManager);
  if (tmp_4 != 0L) {
    tmp = error_buffer_is_empty(rtmGetErrorStatus(TWOWheeled_EV_M));
    if (tmp) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(TWOWheeled_EV_M, msg);
    }
  }

  /* End of MassMatrix for SimscapeExecutionBlock: '<S42>/STATE_1' */
}

void local_evaluateMassMatrix(RTWSolverInfo *si, real_T *Mdest )
{
  /* Refresh global mass matrix */
  TWOWheeled_EV_massmatrix();

  /* Copy the mass matrix from system to the destination, if needed. */
  if (Mdest != rtsiGetSolverMassMatrixPr(si)) {
    real_T *Msrc = rtsiGetSolverMassMatrixPr(si);
    int_T nzmax = rtsiGetSolverMassMatrixNzMax(si);
    (void) memcpy(Mdest, Msrc,
                  (uint_T)nzmax*sizeof(real_T));
  }
}

/* Simplified version of numjac.cpp, for use with RTW. */
void local_numjac( RTWSolverInfo *si, real_T *y, const real_T *Fty, real_T *fac,
                  real_T *dFdy )
{
  /* constants */
  real_T THRESH = 1e-6;
  real_T EPS = 2.2e-16;                /* utGetEps(); */
  real_T BL = pow(EPS, 0.75);
  real_T BU = pow(EPS, 0.25);
  real_T FACMIN = pow(EPS, 0.78);
  real_T FACMAX = 0.1;
  int_T nx = 31;
  real_T *x = rtsiGetContStates(si);
  boolean_T *xdis = rtsiGetContStateDisabledPtr(si);
  real_T del;
  real_T difmax;
  real_T FdelRowmax;
  real_T temp;
  real_T Fdiff;
  real_T maybe;
  real_T xscale;
  real_T fscale;
  real_T *p;
  int_T rowmax;
  int_T i,j;
  if (x != y)
    (void) memcpy(x, y,
                  (uint_T)nx*sizeof(real_T));
  rtsiSetSolverComputingJacobian(si,true);
  for (p = dFdy, j = 0; j < nx; j++, p += nx) {
    /* Zero column j of dFdy if state j is currently disabled. */
    if (xdis[j]) {
      (void) memset(p, 0,
                    (uint_T)nx*sizeof(p[0]));
      continue;
    }

    /* Select an increment del for a difference approximation to
       column j of dFdy.  The vector fac accounts for experience
       gained in previous calls to numjac. */
    xscale = fabs(x[j]);
    if (xscale < THRESH)
      xscale = THRESH;
    temp = (x[j] + fac[j]*xscale);
    del = temp - y[j];
    while (del == 0.0) {
      if (fac[j] < FACMAX) {
        fac[j] *= 100.0;
        if (fac[j] > FACMAX)
          fac[j] = FACMAX;
        temp = (x[j] + fac[j]*xscale);
        del = temp - x[j];
      } else {
        del = THRESH;                  /* thresh is nonzero */
        break;
      }
    }

    /* Keep del pointing into region. */
    if (Fty[j] >= 0.0)
      del = fabs(del);
    else
      del = -fabs(del);

    /* Form a difference approximation to column j of dFdy. */
    temp = x[j];
    x[j] += del;
    TWOWheeled_EV_step();
    rtsiSetdX(si,p);
    TWOWheeled_EV_forcingfunction();
    x[j] = temp;
    difmax = 0.0;
    rowmax = 0;
    FdelRowmax = p[0];
    temp = 1.0 / del;
    for (i = 0; i < nx; i++) {
      Fdiff = p[i] - Fty[i];
      maybe = fabs(Fdiff);
      if (maybe > difmax) {
        difmax = maybe;
        rowmax = i;
        FdelRowmax = p[i];
      }

      p[i] = temp * Fdiff;
    }

    /* Adjust fac for next call to numjac. */
    if (((FdelRowmax != 0.0) && (Fty[rowmax] != 0.0)) || (difmax == 0.0)) {
      fscale = fabs(FdelRowmax);
      if (fscale < fabs(Fty[rowmax]))
        fscale = fabs(Fty[rowmax]);
      if (difmax <= BL*fscale) {
        /* The difference is small, so increase the increment. */
        fac[j] *= 10.0;
        if (fac[j] > FACMAX)
          fac[j] = FACMAX;
      } else if (difmax > BU*fscale) {
        /* The difference is large, so reduce the increment. */
        fac[j] *= 0.1;
        if (fac[j] < FACMIN)
          fac[j] = FACMIN;
      }
    }
  }

  rtsiSetSolverComputingJacobian(si,false);
}                                      /* end local_numjac */

/*
 * This function updates continuous states using the ODE14X fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static int_T rt_ODE14x_N[4] = { 12, 8, 6, 4 };

  time_T t0 = rtsiGetT(si);
  time_T t1 = t0;
  time_T h = rtsiGetStepSize(si);
  real_T *x1 = rtsiGetContStates(si);
  int_T order = rtsiGetSolverExtrapolationOrder(si);
  int_T numIter = rtsiGetSolverNumberNewtonIterations(si);
  ODE14X_IntgData *id = (ODE14X_IntgData *)rtsiGetSolverData(si);
  real_T *x0 = id->x0;
  real_T *f0 = id->f0;
  real_T *x1start = id->x1start;
  real_T *f1 = id->f1;
  real_T *Delta = id->Delta;
  real_T *E = id->E;
  real_T *fac = id->fac;
  real_T *dfdx = id->DFDX;
  real_T *W = id->W;
  int_T *pivots = id->pivots;
  real_T *xtmp = id->xtmp;
  real_T *ztmp = id->ztmp;
  boolean_T *xdis = rtsiGetContStateDisabledPtr(si);
  int_T *Mpattern_ir = rtsiGetSolverMassMatrixIr(si);
  int_T *Mpattern_jc = rtsiGetSolverMassMatrixJc(si);
  real_T *M = id->M;
  int_T col,row,rowidx;
  int_T *N = &(rt_ODE14x_N[0]);
  int_T i,j,k,iter;
  int_T nx = 31;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(x0, x1,
                (uint_T)nx*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  if (id->isFirstStep) {
    local_evaluateMassMatrix(si,M );
    id->isFirstStep = false;
  }

  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  TWOWheeled_EV_forcingfunction();
  local_numjac(si,x0,f0,fac,dfdx );
  for (j = 0; j < order; j++) {
    real_T *p;
    real_T hN = h/N[j];

    /* Get the iteration matrix and solution at t0 */

    /* [L,U] = lu(M - hN*J) */
    (void) memcpy(W, dfdx,
                  (uint_T)nx*nx*sizeof(real_T));
    for (p = W, i = 0; i < nx*nx; i++, p++) {
      *p *= (-hN);
    }

    for (col = 0, p = W; col < nx; col++, p += nx) {
      if (xdis[col]) {
        (void) memset(p, 0,
                      (uint_T)nx*sizeof(p[0]));
        p[col] = 1.0;
      } else {
        for (rowidx = Mpattern_jc[col]; rowidx < Mpattern_jc[col+1]; rowidx++) {
          real_T m_row_col = M[rowidx];
          row = Mpattern_ir[rowidx];
          p[row] += m_row_col;
        }
      }
    }

    rt_lu_real(W, nx,
               pivots);

    /* First Newton's iteration at t0. */
    /* rhs = hN*f0 */
    for (i = 0; i < nx; i++) {
      Delta[i] = hN*f0[i];
    }

    /* Delta = (U \ (L \ rhs)) */
    rt_ForwardSubstitutionRR_Dbl(W, Delta,
      f1, nx,
      1, pivots,
      1);
    rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
      Delta, nx,
      1, 0);

    /* ytmp = y0 + Delta
       ztmp = (ytmp-y0)/h
     */
    (void) memcpy(x1, x0,
                  (uint_T)nx*sizeof(real_T));
    for (i = 0; i < nx; i++) {
      x1[i] += Delta[i];
      ztmp[i] = Delta[i]/hN;
    }

    /* Additional Newton's iterations, if desired.
       for iter = 2:NewtIter
       rhs = hN*feval(odefun,tn,ytmp,extraArgs{:}) - M*(ytmp - yn);
       if statedepM   % only for state dep. Mdel ~= 0
       Mdel = M - feval(massfun,tn,ytmp);
       rhs = rhs + Mdel*ztmp*h;
       end
       Delta = ( U \ ( L \ rhs ) );
       ytmp = ytmp + Delta;
       ztmp = (ytmp - yn)/h
       end
     */
    rtsiSetT(si, t0);
    rtsiSetdX(si, f1);
    for (iter = 1; iter < numIter; iter++) {
      TWOWheeled_EV_step();
      TWOWheeled_EV_forcingfunction();
      for (i = 0; i < nx; i++) {
        Delta[i] = hN*f1[i];
        xtmp[i] = x1[i] - x0[i];
      }

      /* rhs = hN*f(tn,ytmp) - M*(ytmp-yn) */
      for (col = 0; col < nx; col++) {
        for (rowidx = Mpattern_jc[col]; rowidx < Mpattern_jc[col+1]; rowidx++) {
          real_T m_row_col = M[rowidx];
          row = Mpattern_ir[rowidx];
          Delta[row] -= m_row_col*xtmp[col];
        }
      }

      rt_ForwardSubstitutionRR_Dbl(W, Delta,
        f1, nx,
        1, pivots,
        1);
      rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
        Delta, nx,
        1, 0);

      /* ytmp = ytmp + delta
         ztmp = (ytmp - yn)/h
       */
      for (i = 0; i < nx; i++) {
        x1[i] += Delta[i];
        ztmp[i] = (x1[i] - x0[i])/hN;
      }
    }

    /* Steps from t0+hN to t1 -- subintegration of N(j) steps for extrapolation
       ttmp = t0;
       for i = 2:N(j)
       ttmp = ttmp + hN
       ytmp0 = ytmp;
       for iter = 1:NewtIter
       rhs = (ytmp0 - ytmp) + hN*feval(odefun,ttmp,ytmp,extraArgs{:});
       Delta = ( U \ ( L \ rhs ) );
       ytmp = ytmp + Delta;
       end
       end
     */
    for (k = 1; k < N[j]; k++) {
      t1 = t0 + k*hN;
      (void) memcpy(x1start, x1,
                    (uint_T)nx*sizeof(real_T));
      rtsiSetT(si, t1);
      rtsiSetdX(si, f1);
      for (iter = 0; iter < numIter; iter++) {
        TWOWheeled_EV_step();
        TWOWheeled_EV_forcingfunction();
        if (iter == 0) {
          for (i = 0; i < nx; i++) {
            Delta[i] = hN*f1[i];
          }
        } else {
          for (i = 0; i < nx; i++) {
            Delta[i] = hN*f1[i];
            xtmp[i] = (x1[i]-x1start[i]);
          }

          /* rhs = hN*f(tn,ytmp) - M*(ytmp-yn) */
          for (col = 0; col < nx; col++) {
            for (rowidx = Mpattern_jc[col]; rowidx < Mpattern_jc[col+1]; rowidx
                 ++) {
              real_T m_row_col = M[rowidx];
              row = Mpattern_ir[rowidx];
              Delta[row] -= m_row_col*xtmp[col];
            }
          }
        }

        rt_ForwardSubstitutionRR_Dbl(W, Delta,
          f1, nx,
          1, pivots,
          1);
        rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
          Delta, nx,
          1, 0);

        /* ytmp = ytmp + Delta
           ztmp = (ytmp - ytmp0)/h
         */
        for (i = 0; i < nx; i++) {
          x1[i] += Delta[i];
          ztmp[i] = (x1[i] - x1start[i])/hN;
        }
      }
    }

    /* Extrapolate to order j
       E(:,j) = ytmp
       for k = j:-1:2
       coef = N(k-1)/(N(j) - N(k-1))
       E(:,k-1) = E(:,k) + coef*( E(:,k) - E(:,k-1) )
       end
     */
    (void) memcpy(&(E[nx*j]), x1,
                  (uint_T)nx*sizeof(real_T));
    for (k = j; k > 0; k--) {
      real_T coef = (real_T)(N[k-1]) / (N[j]-N[k-1]);
      for (i = 0; i < nx; i++) {
        x1[i] = E[nx*k+i] + coef*(E[nx*k+i] - E[nx*(k-1)+i]);
      }

      (void) memcpy(&(E[nx*(k-1)]), x1,
                    (uint_T)nx*sizeof(real_T));
    }
  }

  /* x1 = E(:,1); */
  (void) memcpy(x1, E,
                (uint_T)nx*sizeof(real_T));

  /* t1 = t0 + h; */
  rtsiSetT(si,rtsiGetSolverStopTime(si));
  TWOWheeled_EV_step();
  TWOWheeled_EV_projection();
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * System initialize for enable system:
 *    '<S29>/Pass Through'
 *    '<S31>/Pass Through'
 */
void TWOWheeled_EV_PassThrough_Init(real_T rtp_IC, B_PassThrough_TWOWheeled_EV_T
  *localB)
{
  /* SystemInitialize for SignalConversion generated from: '<S30>/u' incorporates:
   *  Outport: '<S30>/y'
   */
  localB->u = rtp_IC;
}

/*
 * Disable for enable system:
 *    '<S29>/Pass Through'
 *    '<S31>/Pass Through'
 */
void TWOWheeled__PassThrough_Disable(DW_PassThrough_TWOWheeled_EV_T *localDW)
{
  localDW->PassThrough_MODE = false;
}

/*
 * Output and update for enable system:
 *    '<S29>/Pass Through'
 *    '<S31>/Pass Through'
 */
void TWOWheeled_EV_PassThrough(RT_MODEL_TWOWheeled_EV_T * const TWOWheeled_EV_M,
  boolean_T rtu_Enable, real_T rtu_u, B_PassThrough_TWOWheeled_EV_T *localB,
  DW_PassThrough_TWOWheeled_EV_T *localDW)
{
  /* Outputs for Enabled SubSystem: '<S29>/Pass Through' incorporates:
   *  EnablePort: '<S30>/Enable'
   */
  if (rtmIsMajorTimeStep(TWOWheeled_EV_M) && rtsiIsModeUpdateTimeStep
      (&TWOWheeled_EV_M->solverInfo)) {
    if (rtu_Enable) {
      localDW->PassThrough_MODE = true;
    } else if (localDW->PassThrough_MODE) {
      TWOWheeled__PassThrough_Disable(localDW);
    }
  }

  if (localDW->PassThrough_MODE) {
    /* SignalConversion generated from: '<S30>/u' */
    localB->u = rtu_u;
  }

  /* End of Outputs for SubSystem: '<S29>/Pass Through' */
}

real_T rt_modd_snf(real_T u0, real_T u1)
{
  real_T q;
  real_T y;
  boolean_T yEq;
  y = u0;
  if (u1 == 0.0) {
    if (u0 == 0.0) {
      y = u1;
    }
  } else if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (u0 == 0.0) {
    y = 0.0 / u1;
  } else if (rtIsInf(u1)) {
    if ((u1 < 0.0) != (u0 < 0.0)) {
      y = u1;
    }
  } else {
    y = fmod(u0, u1);
    yEq = (y == 0.0);
    if ((!yEq) && (u1 > floor(u1))) {
      q = fabs(u0 / u1);
      yEq = !(fabs(q - floor(q + 0.5)) > DBL_EPSILON * q);
    }

    if (yEq) {
      y = u1 * 0.0;
    } else if ((u0 < 0.0) != (u1 < 0.0)) {
      y += u1;
    }
  }

  return y;
}

/* Model step function */
void TWOWheeled_EV_step(void)
{
  /* local block i/o variables */
  real_T rtb_FromWs;
  real_T rtb_Gain1;
  if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
    /* set solver stop time */
    if (!(TWOWheeled_EV_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&TWOWheeled_EV_M->solverInfo,
                            ((TWOWheeled_EV_M->Timing.clockTickH0 + 1) *
        TWOWheeled_EV_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&TWOWheeled_EV_M->solverInfo,
                            ((TWOWheeled_EV_M->Timing.clockTick0 + 1) *
        TWOWheeled_EV_M->Timing.stepSize0 + TWOWheeled_EV_M->Timing.clockTickH0 *
        TWOWheeled_EV_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(TWOWheeled_EV_M)) {
    TWOWheeled_EV_M->Timing.t[0] = rtsiGetT(&TWOWheeled_EV_M->solverInfo);
  }

  {
    NeslSimulationData *simulationData;
    NeuDiagnosticManager *diagnosticManager;
    NeuDiagnosticTree *diagnosticTree;
    NeuDiagnosticTree *diagnosticTree_0;
    char *msg;
    char *msg_0;
    real_T tmp_3[49];
    real_T tmp[16];
    real_T rtb_Add;
    real_T time;
    real_T time_0;
    real_T time_1;
    real_T time_2;
    real_T time_tmp;
    int32_T tmp_1;
    int16_T i;
    int_T tmp_4[6];
    int_T tmp_0[5];
    boolean_T tmp_5[6];
    boolean_T LogicalOperator2_tmp;
    boolean_T LogicalOperator2_tmp_tmp;
    boolean_T tmp_2;
    LogicalOperator2_tmp_tmp = rtmIsMajorTimeStep(TWOWheeled_EV_M);
    LogicalOperator2_tmp = !LogicalOperator2_tmp_tmp;

    /* Logic: '<S27>/Logical Operator2' */
    TWOWheeled_EV_B.LogicalOperator2 = (LogicalOperator2_tmp &&
      TWOWheeled_EV_B.LogicalOperator2);

    /* Clock: '<S14>/Clock' incorporates:
     *  SimscapeExecutionBlock: '<S42>/OUTPUT_1_0'
     *  SimscapeExecutionBlock: '<S42>/STATE_1'
     */
    rtb_Add = TWOWheeled_EV_M->Timing.t[0];

    /* FromWorkspace: '<S7>/FromWs' */
    {
      real_T *pDataValues = (real_T *) TWOWheeled_EV_DW.FromWs_PWORK.DataPtr;
      real_T *pTimeValues = (real_T *) TWOWheeled_EV_DW.FromWs_PWORK.TimePtr;
      int_T currTimeIndex = TWOWheeled_EV_DW.FromWs_IWORK.PrevIndex;
      real_T t = TWOWheeled_EV_M->Timing.t[0];

      /* Get index */
      if (t <= pTimeValues[0]) {
        currTimeIndex = 0;
      } else if (t >= pTimeValues[7]) {
        currTimeIndex = 6;
      } else {
        if (t < pTimeValues[currTimeIndex]) {
          while (t < pTimeValues[currTimeIndex]) {
            currTimeIndex--;
          }
        } else {
          while (t >= pTimeValues[currTimeIndex + 1]) {
            currTimeIndex++;
          }
        }
      }

      TWOWheeled_EV_DW.FromWs_IWORK.PrevIndex = currTimeIndex;

      /* Post output */
      {
        real_T t1 = pTimeValues[currTimeIndex];
        real_T t2 = pTimeValues[currTimeIndex + 1];
        if (t1 == t2) {
          if (t < t1) {
            rtb_FromWs = pDataValues[currTimeIndex];
          } else {
            rtb_FromWs = pDataValues[currTimeIndex + 1];
          }
        } else {
          real_T f1 = (t2 - t) / (t2 - t1);
          real_T f2 = 1.0 - f1;
          real_T d1;
          real_T d2;
          int_T TimeIndex = currTimeIndex;
          d1 = pDataValues[TimeIndex];
          d2 = pDataValues[TimeIndex + 1];
          rtb_FromWs = (real_T) rtInterpolate(d1, d2, f1, f2);
          pDataValues += 8;
        }
      }
    }

    /* MultiPortSwitch: '<Root>/Multiport Switch' incorporates:
     *  Constant: '<Root>/Constant'
     */
    if ((int16_T)TWOWheeled_EV_P.Constant_Value_c == 1) {
      /* Switch: '<S14>/Switch' incorporates:
       *  Constant: '<S14>/repeat'
       */
      if (TWOWheeled_EV_P.repeat_Value > TWOWheeled_EV_P.Switch_Threshold) {
        /* Switch: '<S14>/Switch' incorporates:
         *  Clock: '<S14>/Clock'
         *  Constant: '<S14>/tFinal'
         *  Math: '<S14>/Math Function'
         */
        TWOWheeled_EV_B.Switch_o = rt_modd_snf(rtb_Add,
          TWOWheeled_EV_P.tFinal_Value);
      } else {
        /* Switch: '<S14>/Switch' incorporates:
         *  Clock: '<S14>/Clock'
         */
        TWOWheeled_EV_B.Switch_o = rtb_Add;
      }

      /* End of Switch: '<S14>/Switch' */

      /* MultiPortSwitch: '<Root>/Multiport Switch' incorporates:
       *  Lookup_n-D: '<S14>/1-D Lookup Table'
       *  Switch: '<S14>/Switch'
       */
      TWOWheeled_EV_B.MultiportSwitch = look1_pbinlcapw(TWOWheeled_EV_B.Switch_o,
        TWOWheeled_EV_P.uDLookupTable_bp01Data,
        TWOWheeled_EV_P.uDLookupTable_tableData, &TWOWheeled_EV_DW.m_bpIndex,
        2473UL);
    } else {
      /* MultiPortSwitch: '<Root>/Multiport Switch' */
      TWOWheeled_EV_B.MultiportSwitch = rtb_FromWs;
    }

    /* End of MultiPortSwitch: '<Root>/Multiport Switch' */

    /* Gain: '<S26>/Kff//vnom' */
    TWOWheeled_EV_B.Kffvnom = TWOWheeled_EV_P.LongitudinalDriver_Kff /
      TWOWheeled_EV_P.LongitudinalDriver_vnom * TWOWheeled_EV_B.MultiportSwitch;

    /* Integrator: '<S34>/Integrator' */
    TWOWheeled_EV_B.Integrator = TWOWheeled_EV_X.Integrator_CSTATE;

    /* Gain: '<S26>/Kp//vnom' */
    TWOWheeled_EV_B.Kpvnom = TWOWheeled_EV_P.LongitudinalDriver_Kp /
      TWOWheeled_EV_P.LongitudinalDriver_vnom * TWOWheeled_EV_B.Integrator;

    /* Integrator: '<S26>/Integrator1' */
    /* Limited  Integrator  */
    if (TWOWheeled_EV_X.Integrator1_CSTATE >=
        TWOWheeled_EV_P.Integrator1_UpperSat) {
      TWOWheeled_EV_X.Integrator1_CSTATE = TWOWheeled_EV_P.Integrator1_UpperSat;
    } else if (TWOWheeled_EV_X.Integrator1_CSTATE <=
               TWOWheeled_EV_P.Integrator1_LowerSat) {
      TWOWheeled_EV_X.Integrator1_CSTATE = TWOWheeled_EV_P.Integrator1_LowerSat;
    }

    /* Integrator: '<S26>/Integrator1' */
    TWOWheeled_EV_B.Integrator1 = TWOWheeled_EV_X.Integrator1_CSTATE;
    if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
      /* Gain: '<S26>/Kg' incorporates:
       *  Constant: '<Root>/Grade'
       */
      TWOWheeled_EV_B.Kg = TWOWheeled_EV_P.LongitudinalDriver_Kg *
        TWOWheeled_EV_P.Grade_Value;
    }

    /* Sum: '<S26>/Sum1' */
    TWOWheeled_EV_B.Sum1 = ((TWOWheeled_EV_B.Kffvnom + TWOWheeled_EV_B.Kpvnom) +
      TWOWheeled_EV_B.Integrator1) + TWOWheeled_EV_B.Kg;

    /* Saturate: '<S26>/-1 to 1 ' */
    if (TWOWheeled_EV_B.Sum1 > TWOWheeled_EV_P.uto1_UpperSat) {
      /* Saturate: '<S26>/-1 to 1 ' */
      TWOWheeled_EV_B.uto1 = TWOWheeled_EV_P.uto1_UpperSat;
    } else if (TWOWheeled_EV_B.Sum1 < TWOWheeled_EV_P.uto1_LowerSat) {
      /* Saturate: '<S26>/-1 to 1 ' */
      TWOWheeled_EV_B.uto1 = TWOWheeled_EV_P.uto1_LowerSat;
    } else {
      /* Saturate: '<S26>/-1 to 1 ' */
      TWOWheeled_EV_B.uto1 = TWOWheeled_EV_B.Sum1;
    }

    /* End of Saturate: '<S26>/-1 to 1 ' */

    /* Switch: '<S27>/Switch1' incorporates:
     *  Saturate: '<S27>/0~1'
     *  Switch: '<S27>/Switch'
     */
    if (TWOWheeled_EV_B.LogicalOperator2) {
      /* Switch: '<S27>/Switch1' */
      TWOWheeled_EV_B.Switch1 = 0.0;
    } else {
      if (TWOWheeled_EV_B.uto1 > TWOWheeled_EV_P.u1_UpperSat) {
        /* Saturate: '<S27>/0~1' incorporates:
         *  Switch: '<S27>/Switch'
         */
        TWOWheeled_EV_B.Accel = TWOWheeled_EV_P.u1_UpperSat;
      } else if (TWOWheeled_EV_B.uto1 < TWOWheeled_EV_P.u1_LowerSat) {
        /* Saturate: '<S27>/0~1' incorporates:
         *  Switch: '<S27>/Switch'
         */
        TWOWheeled_EV_B.Accel = TWOWheeled_EV_P.u1_LowerSat;
      } else {
        /* Saturate: '<S27>/0~1' incorporates:
         *  Switch: '<S27>/Switch'
         */
        TWOWheeled_EV_B.Accel = TWOWheeled_EV_B.uto1;
      }

      /* Switch: '<S27>/Switch' incorporates:
       *  Saturate: '<S27>/0~1'
       */
      TWOWheeled_EV_B.ut_c = TWOWheeled_EV_B.Accel;

      /* Switch: '<S27>/Switch1' */
      TWOWheeled_EV_B.Switch1 = TWOWheeled_EV_B.ut_c;
    }

    /* End of Switch: '<S27>/Switch1' */

    /* Saturate: '<S27>/Saturation' */
    if (TWOWheeled_EV_B.Switch1 > TWOWheeled_EV_P.Saturation_UpperSat) {
      /* Saturate: '<S27>/Saturation' */
      TWOWheeled_EV_B.Saturation = TWOWheeled_EV_P.Saturation_UpperSat;
    } else if (TWOWheeled_EV_B.Switch1 < TWOWheeled_EV_P.Saturation_LowerSat) {
      /* Saturate: '<S27>/Saturation' */
      TWOWheeled_EV_B.Saturation = TWOWheeled_EV_P.Saturation_LowerSat;
    } else {
      /* Saturate: '<S27>/Saturation' */
      TWOWheeled_EV_B.Saturation = TWOWheeled_EV_B.Switch1;
    }

    /* End of Saturate: '<S27>/Saturation' */

    /* Logic: '<S29>/NOT' */
    TWOWheeled_EV_B.NOT = (LogicalOperator2_tmp_tmp || TWOWheeled_EV_B.NOT);

    /* Outputs for Enabled SubSystem: '<S29>/Pass Through' */
    TWOWheeled_EV_PassThrough(TWOWheeled_EV_M, TWOWheeled_EV_B.NOT,
      TWOWheeled_EV_B.Saturation, &TWOWheeled_EV_B.PassThrough,
      &TWOWheeled_EV_DW.PassThrough);

    /* End of Outputs for SubSystem: '<S29>/Pass Through' */

    /* SimscapeInputBlock: '<S42>/INPUT_1_1_1' */
    TWOWheeled_EV_B.INPUT_1_1_1[0] = TWOWheeled_EV_B.PassThrough.u;
    TWOWheeled_EV_B.INPUT_1_1_1[1] = 0.0;
    TWOWheeled_EV_B.INPUT_1_1_1[2] = 0.0;
    TWOWheeled_EV_B.INPUT_1_1_1[3] = 0.0;

    /* Logic: '<S28>/Logical Operator2' */
    TWOWheeled_EV_B.LogicalOperator2_f = (LogicalOperator2_tmp &&
      TWOWheeled_EV_B.LogicalOperator2_f);

    /* Switch: '<S28>/Switch1' incorporates:
     *  Saturate: '<S28>/-1~0'
     *  Switch: '<S28>/Switch'
     */
    if (TWOWheeled_EV_B.LogicalOperator2_f) {
      /* Switch: '<S28>/Switch1' */
      TWOWheeled_EV_B.Switch1_e = 0.0;
    } else {
      if (TWOWheeled_EV_B.uto1 > TWOWheeled_EV_P.u0_UpperSat) {
        /* Saturate: '<S28>/-1~0' incorporates:
         *  Switch: '<S28>/Switch'
         */
        TWOWheeled_EV_B.u0 = TWOWheeled_EV_P.u0_UpperSat;
      } else if (TWOWheeled_EV_B.uto1 < TWOWheeled_EV_P.u0_LowerSat) {
        /* Saturate: '<S28>/-1~0' incorporates:
         *  Switch: '<S28>/Switch'
         */
        TWOWheeled_EV_B.u0 = TWOWheeled_EV_P.u0_LowerSat;
      } else {
        /* Saturate: '<S28>/-1~0' incorporates:
         *  Switch: '<S28>/Switch'
         */
        TWOWheeled_EV_B.u0 = TWOWheeled_EV_B.uto1;
      }

      /* Switch: '<S28>/Switch' incorporates:
       *  Saturate: '<S28>/-1~0'
       *  UnaryMinus: '<S28>/Unary Minus'
       */
      TWOWheeled_EV_B.ut = -TWOWheeled_EV_B.u0;

      /* Switch: '<S28>/Switch1' */
      TWOWheeled_EV_B.Switch1_e = TWOWheeled_EV_B.ut;
    }

    /* End of Switch: '<S28>/Switch1' */

    /* Saturate: '<S28>/Saturation' */
    if (TWOWheeled_EV_B.Switch1_e > TWOWheeled_EV_P.Saturation_UpperSat_g) {
      /* Saturate: '<S28>/Saturation' */
      TWOWheeled_EV_B.Saturation_l = TWOWheeled_EV_P.Saturation_UpperSat_g;
    } else if (TWOWheeled_EV_B.Switch1_e < TWOWheeled_EV_P.Saturation_LowerSat_m)
    {
      /* Saturate: '<S28>/Saturation' */
      TWOWheeled_EV_B.Saturation_l = TWOWheeled_EV_P.Saturation_LowerSat_m;
    } else {
      /* Saturate: '<S28>/Saturation' */
      TWOWheeled_EV_B.Saturation_l = TWOWheeled_EV_B.Switch1_e;
    }

    /* End of Saturate: '<S28>/Saturation' */

    /* Logic: '<S31>/NOT' */
    TWOWheeled_EV_B.NOT_c = (LogicalOperator2_tmp_tmp || TWOWheeled_EV_B.NOT_c);

    /* Outputs for Enabled SubSystem: '<S31>/Pass Through' */
    TWOWheeled_EV_PassThrough(TWOWheeled_EV_M, TWOWheeled_EV_B.NOT_c,
      TWOWheeled_EV_B.Saturation_l, &TWOWheeled_EV_B.PassThrough_h,
      &TWOWheeled_EV_DW.PassThrough_h);

    /* End of Outputs for SubSystem: '<S31>/Pass Through' */

    /* SimscapeInputBlock: '<S42>/INPUT_2_1_1' */
    TWOWheeled_EV_B.INPUT_2_1_1[0] = TWOWheeled_EV_B.PassThrough_h.u;
    TWOWheeled_EV_B.INPUT_2_1_1[1] = 0.0;
    TWOWheeled_EV_B.INPUT_2_1_1[2] = 0.0;
    TWOWheeled_EV_B.INPUT_2_1_1[3] = 0.0;

    /* SimscapeInputBlock: '<S42>/INPUT_4_1_1' incorporates:
     *  Constant: '<S11>/Headwind speed'
     */
    TWOWheeled_EV_B.INPUT_4_1_1[0] = TWOWheeled_EV_P.Headwindspeed_Value;
    TWOWheeled_EV_B.INPUT_4_1_1[1] = 0.0;
    TWOWheeled_EV_B.INPUT_4_1_1[2] = 0.0;
    if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
      TWOWheeled_EV_DW.INPUT_4_1_1_Discrete[0] = !(TWOWheeled_EV_B.INPUT_4_1_1[0]
        == TWOWheeled_EV_DW.INPUT_4_1_1_Discrete[1]);
      TWOWheeled_EV_DW.INPUT_4_1_1_Discrete[1] = TWOWheeled_EV_B.INPUT_4_1_1[0];
    }

    TWOWheeled_EV_B.INPUT_4_1_1[0] = TWOWheeled_EV_DW.INPUT_4_1_1_Discrete[1];
    TWOWheeled_EV_B.INPUT_4_1_1[3] = TWOWheeled_EV_DW.INPUT_4_1_1_Discrete[0];

    /* End of SimscapeInputBlock: '<S42>/INPUT_4_1_1' */

    /* SimscapeInputBlock: '<S42>/INPUT_3_1_1' incorporates:
     *  Constant: '<S11>/Grade'
     */
    TWOWheeled_EV_B.INPUT_3_1_1[0] = TWOWheeled_EV_P.Grade_Value_m;
    TWOWheeled_EV_B.INPUT_3_1_1[1] = 0.0;
    TWOWheeled_EV_B.INPUT_3_1_1[2] = 0.0;
    if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
      TWOWheeled_EV_DW.INPUT_3_1_1_Discrete[0] = !(TWOWheeled_EV_B.INPUT_3_1_1[0]
        == TWOWheeled_EV_DW.INPUT_3_1_1_Discrete[1]);
      TWOWheeled_EV_DW.INPUT_3_1_1_Discrete[1] = TWOWheeled_EV_B.INPUT_3_1_1[0];
    }

    TWOWheeled_EV_B.INPUT_3_1_1[0] = TWOWheeled_EV_DW.INPUT_3_1_1_Discrete[1];
    TWOWheeled_EV_B.INPUT_3_1_1[3] = TWOWheeled_EV_DW.INPUT_3_1_1_Discrete[0];

    /* End of SimscapeInputBlock: '<S42>/INPUT_3_1_1' */

    /* SimscapeExecutionBlock: '<S42>/STATE_1' incorporates:
     *  SimscapeExecutionBlock: '<S42>/OUTPUT_1_0'
     */
    simulationData = (NeslSimulationData *)TWOWheeled_EV_DW.STATE_1_SimData;
    time = rtb_Add;
    simulationData->mData->mTime.mN = 1;
    simulationData->mData->mTime.mX = &time;
    simulationData->mData->mContStates.mN = 27;
    simulationData->mData->mContStates.mX =
      &TWOWheeled_EV_X.TWOWheeled_EVBatterycharge[0];
    simulationData->mData->mDiscStates.mN = 0;
    simulationData->mData->mDiscStates.mX = &TWOWheeled_EV_DW.STATE_1_Discrete;
    simulationData->mData->mModeVector.mN = 6;
    simulationData->mData->mModeVector.mX = (int32_T *)
      &TWOWheeled_EV_DW.STATE_1_Modes[0];
    LogicalOperator2_tmp_tmp = false;
    simulationData->mData->mFoundZcEvents = LogicalOperator2_tmp_tmp;
    LogicalOperator2_tmp_tmp = rtmIsMajorTimeStep(TWOWheeled_EV_M);
    simulationData->mData->mIsMajorTimeStep = LogicalOperator2_tmp_tmp;
    LogicalOperator2_tmp = false;
    simulationData->mData->mIsSolverAssertCheck = LogicalOperator2_tmp;
    simulationData->mData->mIsSolverCheckingCIC = false;
    LogicalOperator2_tmp = rtsiIsSolverComputingJacobian
      (&TWOWheeled_EV_M->solverInfo);
    simulationData->mData->mIsComputingJacobian = LogicalOperator2_tmp;
    simulationData->mData->mIsEvaluatingF0 = false;
    simulationData->mData->mIsSolverRequestingReset = false;
    LogicalOperator2_tmp = rtsiIsModeUpdateTimeStep(&TWOWheeled_EV_M->solverInfo);
    simulationData->mData->mIsModeUpdateTimeStep = LogicalOperator2_tmp;
    tmp_0[0] = 0;
    tmp[0] = TWOWheeled_EV_B.INPUT_1_1_1[0];
    tmp[1] = TWOWheeled_EV_B.INPUT_1_1_1[1];
    tmp[2] = TWOWheeled_EV_B.INPUT_1_1_1[2];
    tmp[3] = TWOWheeled_EV_B.INPUT_1_1_1[3];
    tmp_0[1] = 4;
    tmp[4] = TWOWheeled_EV_B.INPUT_2_1_1[0];
    tmp[5] = TWOWheeled_EV_B.INPUT_2_1_1[1];
    tmp[6] = TWOWheeled_EV_B.INPUT_2_1_1[2];
    tmp[7] = TWOWheeled_EV_B.INPUT_2_1_1[3];
    tmp_0[2] = 8;
    tmp[8] = TWOWheeled_EV_B.INPUT_4_1_1[0];
    tmp[9] = TWOWheeled_EV_B.INPUT_4_1_1[1];
    tmp[10] = TWOWheeled_EV_B.INPUT_4_1_1[2];
    tmp[11] = TWOWheeled_EV_B.INPUT_4_1_1[3];
    tmp_0[3] = 12;
    tmp[12] = TWOWheeled_EV_B.INPUT_3_1_1[0];
    tmp[13] = TWOWheeled_EV_B.INPUT_3_1_1[1];
    tmp[14] = TWOWheeled_EV_B.INPUT_3_1_1[2];
    tmp[15] = TWOWheeled_EV_B.INPUT_3_1_1[3];
    tmp_0[4] = 16;
    simulationData->mData->mInputValues.mN = 16;
    simulationData->mData->mInputValues.mX = &tmp[0];
    simulationData->mData->mInputOffsets.mN = 5;
    simulationData->mData->mInputOffsets.mX = (int32_T *)&tmp_0[0];
    simulationData->mData->mOutputs.mN = 33;
    simulationData->mData->mOutputs.mX = &TWOWheeled_EV_B.STATE_1[0];
    simulationData->mData->mTolerances.mN = 0;
    simulationData->mData->mTolerances.mX = NULL;
    simulationData->mData->mCstateHasChanged = false;
    time_tmp = TWOWheeled_EV_M->Timing.t[0];
    time_0 = time_tmp;
    simulationData->mData->mTime.mN = 1;
    simulationData->mData->mTime.mX = &time_0;
    simulationData->mData->mSampleHits.mN = 0;
    simulationData->mData->mSampleHits.mX = NULL;
    simulationData->mData->mIsFundamentalSampleHit = false;
    diagnosticManager = (NeuDiagnosticManager *)TWOWheeled_EV_DW.STATE_1_DiagMgr;
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    tmp_1 = ne_simulator_method((NeslSimulator *)
      TWOWheeled_EV_DW.STATE_1_Simulator, NESL_SIM_OUTPUTS, simulationData,
      diagnosticManager);
    if (tmp_1 != 0L) {
      tmp_2 = error_buffer_is_empty(rtmGetErrorStatus(TWOWheeled_EV_M));
      if (tmp_2) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(TWOWheeled_EV_M, msg);
      }
    }

    /* SimscapeExecutionBlock: '<S42>/OUTPUT_1_0' */
    simulationData = (NeslSimulationData *)TWOWheeled_EV_DW.OUTPUT_1_0_SimData;
    time_1 = rtb_Add;
    simulationData->mData->mTime.mN = 1;
    simulationData->mData->mTime.mX = &time_1;
    simulationData->mData->mContStates.mN = 0;
    simulationData->mData->mContStates.mX = NULL;
    simulationData->mData->mDiscStates.mN = 0;
    simulationData->mData->mDiscStates.mX =
      &TWOWheeled_EV_DW.OUTPUT_1_0_Discrete;
    simulationData->mData->mModeVector.mN = 0;
    simulationData->mData->mModeVector.mX = (int32_T *)
      &TWOWheeled_EV_DW.OUTPUT_1_0_Modes;
    tmp_2 = false;
    simulationData->mData->mFoundZcEvents = tmp_2;
    simulationData->mData->mIsMajorTimeStep = LogicalOperator2_tmp_tmp;
    LogicalOperator2_tmp_tmp = false;
    simulationData->mData->mIsSolverAssertCheck = LogicalOperator2_tmp_tmp;
    simulationData->mData->mIsSolverCheckingCIC = false;
    simulationData->mData->mIsComputingJacobian = false;
    simulationData->mData->mIsEvaluatingF0 = false;
    simulationData->mData->mIsSolverRequestingReset = false;
    simulationData->mData->mIsModeUpdateTimeStep = LogicalOperator2_tmp;
    tmp_4[0] = 0;
    tmp_3[0] = TWOWheeled_EV_B.INPUT_1_1_1[0];
    tmp_3[1] = TWOWheeled_EV_B.INPUT_1_1_1[1];
    tmp_3[2] = TWOWheeled_EV_B.INPUT_1_1_1[2];
    tmp_3[3] = TWOWheeled_EV_B.INPUT_1_1_1[3];
    tmp_4[1] = 4;
    tmp_3[4] = TWOWheeled_EV_B.INPUT_2_1_1[0];
    tmp_3[5] = TWOWheeled_EV_B.INPUT_2_1_1[1];
    tmp_3[6] = TWOWheeled_EV_B.INPUT_2_1_1[2];
    tmp_3[7] = TWOWheeled_EV_B.INPUT_2_1_1[3];
    tmp_4[2] = 8;
    tmp_3[8] = TWOWheeled_EV_B.INPUT_4_1_1[0];
    tmp_3[9] = TWOWheeled_EV_B.INPUT_4_1_1[1];
    tmp_3[10] = TWOWheeled_EV_B.INPUT_4_1_1[2];
    tmp_3[11] = TWOWheeled_EV_B.INPUT_4_1_1[3];
    tmp_4[3] = 12;
    tmp_3[12] = TWOWheeled_EV_B.INPUT_3_1_1[0];
    tmp_3[13] = TWOWheeled_EV_B.INPUT_3_1_1[1];
    tmp_3[14] = TWOWheeled_EV_B.INPUT_3_1_1[2];
    tmp_3[15] = TWOWheeled_EV_B.INPUT_3_1_1[3];
    tmp_4[4] = 16;
    memcpy(&tmp_3[16], &TWOWheeled_EV_B.STATE_1[0], 33U * sizeof(real_T));
    tmp_4[5] = 49;
    simulationData->mData->mInputValues.mN = 49;
    simulationData->mData->mInputValues.mX = &tmp_3[0];
    simulationData->mData->mInputOffsets.mN = 6;
    simulationData->mData->mInputOffsets.mX = (int32_T *)&tmp_4[0];
    simulationData->mData->mOutputs.mN = 3;
    simulationData->mData->mOutputs.mX = &TWOWheeled_EV_B.OUTPUT_1_0[0];
    simulationData->mData->mTolerances.mN = 0;
    simulationData->mData->mTolerances.mX = NULL;
    simulationData->mData->mCstateHasChanged = false;
    time_2 = time_tmp;
    simulationData->mData->mTime.mN = 1;
    simulationData->mData->mTime.mX = &time_2;
    simulationData->mData->mSampleHits.mN = 0;
    simulationData->mData->mSampleHits.mX = NULL;
    simulationData->mData->mIsFundamentalSampleHit = false;
    diagnosticManager = (NeuDiagnosticManager *)
      TWOWheeled_EV_DW.OUTPUT_1_0_DiagMgr;
    diagnosticTree_0 = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    tmp_1 = ne_simulator_method((NeslSimulator *)
      TWOWheeled_EV_DW.OUTPUT_1_0_Simulator, NESL_SIM_OUTPUTS, simulationData,
      diagnosticManager);
    if (tmp_1 != 0L) {
      LogicalOperator2_tmp_tmp = error_buffer_is_empty(rtmGetErrorStatus
        (TWOWheeled_EV_M));
      if (LogicalOperator2_tmp_tmp) {
        msg_0 = rtw_diagnostics_msg(diagnosticTree_0);
        rtmSetErrorStatus(TWOWheeled_EV_M, msg_0);
      }
    }

    if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
      /* UnitDelay: '<S35>/Unit Delay' */
      TWOWheeled_EV_B.UnitDelay[0] = TWOWheeled_EV_DW.UnitDelay_DSTATE[0];
      TWOWheeled_EV_B.UnitDelay[1] = TWOWheeled_EV_DW.UnitDelay_DSTATE[1];
    }

    /* Sum: '<S21>/Sum7' */
    TWOWheeled_EV_B.Sum7 = TWOWheeled_EV_B.MultiportSwitch -
      TWOWheeled_EV_B.OUTPUT_1_0[1];

    /* SignalConversion generated from: '<S35>/Vector Concatenate1' incorporates:
     *  Concatenate: '<S35>/Vector Concatenate1'
     *  UnitConversion: '<S24>/Unit Conversion'
     */
    /* Unit Conversion - from: m/s to: m/s
       Expression: output = (1*input) + (0) */
    TWOWheeled_EV_B.Switch[0] = TWOWheeled_EV_B.Sum7;

    /* SignalConversion generated from: '<S35>/Vector Concatenate1' incorporates:
     *  Concatenate: '<S35>/Vector Concatenate1'
     *  UnitConversion: '<S24>/Unit Conversion'
     */
    TWOWheeled_EV_B.Switch[1] = TWOWheeled_EV_B.Sum7;

    /* Switch: '<S35>/Switch' incorporates:
     *  RelationalOperator: '<S35>/Relational Operator'
     *  RelationalOperator: '<S35>/Relational Operator1'
     *  UnitConversion: '<S24>/Unit Conversion'
     */
    if (!(TWOWheeled_EV_B.Sum7 > TWOWheeled_EV_B.UnitDelay[0])) {
      /* Switch: '<S35>/Switch' incorporates:
       *  Concatenate: '<S35>/Vector Concatenate1'
       */
      TWOWheeled_EV_B.Switch[0] = TWOWheeled_EV_B.UnitDelay[0];
    }

    if (!(TWOWheeled_EV_B.Sum7 < TWOWheeled_EV_B.UnitDelay[1])) {
      /* Switch: '<S35>/Switch' incorporates:
       *  Concatenate: '<S35>/Vector Concatenate1'
       */
      TWOWheeled_EV_B.Switch[1] = TWOWheeled_EV_B.UnitDelay[1];
    }

    /* End of Switch: '<S35>/Switch' */

    /* Product: '<S35>/Product' incorporates:
     *  UnitConversion: '<S24>/Unit Conversion'
     */
    TWOWheeled_EV_B.Product = TWOWheeled_EV_B.Sum7 * TWOWheeled_EV_B.Sum7;
    if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
      /* Logic: '<S21>/Logical Operator' */
      TWOWheeled_EV_B.LogicalOperator = false;
      tmp_5[0] = false;
      tmp_5[1] = false;
      tmp_5[2] = false;
      tmp_5[3] = false;
      tmp_5[4] = false;
      tmp_5[5] = false;
      for (i = 0; i < 5; i++) {
        TWOWheeled_EV_B.LogicalOperator = (TWOWheeled_EV_B.LogicalOperator ||
          tmp_5[i + 1]);
      }

      /* End of Logic: '<S21>/Logical Operator' */
    }

    /* Switch: '<S21>/Switch' */
    if (TWOWheeled_EV_B.LogicalOperator) {
      /* Switch: '<S21>/Switch' incorporates:
       *  Constant: '<S21>/Zero'
       */
      TWOWheeled_EV_B.Switch_c = TWOWheeled_EV_P.Zero_Value;
    } else {
      /* Switch: '<S21>/Switch' */
      TWOWheeled_EV_B.Switch_c = TWOWheeled_EV_B.Sum7;
    }

    /* End of Switch: '<S21>/Switch' */

    /* Sum: '<S34>/Sum' */
    TWOWheeled_EV_B.Sum = TWOWheeled_EV_B.Switch_c - TWOWheeled_EV_B.Integrator;
    if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
      /* Product: '<S33>/Divide' incorporates:
       *  Constant: '<S33>/Constant'
       */
      TWOWheeled_EV_B.Divide = 1.0 / TWOWheeled_EV_P.LongitudinalDriver_tauerr;
    }

    /* Product: '<S34>/Product' */
    TWOWheeled_EV_B.Product_d = TWOWheeled_EV_B.Divide * TWOWheeled_EV_B.Sum;
    if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
      /* Gain: '<S6>/Gain1' */
      rtb_Gain1 = TWOWheeled_EV_P.Gain1_Gain * TWOWheeled_EV_B.OUTPUT_1_0[0];

      /* DigitalClock: '<S14>/Digital Clock' */
      TWOWheeled_EV_B.DigitalClock = (((TWOWheeled_EV_M->Timing.clockTick1+
        TWOWheeled_EV_M->Timing.clockTickH1* 4294967296.0)) * 0.5);
    }

    /* Sum: '<S26>/Sum5' incorporates:
     *  UnaryMinus: '<S28>/Unary Minus1'
     */
    TWOWheeled_EV_B.Sum5 = (TWOWheeled_EV_B.PassThrough.u - TWOWheeled_EV_B.Sum1)
      - (-TWOWheeled_EV_B.PassThrough_h.u);

    /* Gain: '<S26>/Kaw' */
    TWOWheeled_EV_B.Kaw = TWOWheeled_EV_P.LongitudinalDriver_Kaw *
      TWOWheeled_EV_B.Sum5;

    /* Gain: '<S26>/Ki//vnom' */
    TWOWheeled_EV_B.Kivnom = TWOWheeled_EV_P.LongitudinalDriver_Ki /
      TWOWheeled_EV_P.LongitudinalDriver_vnom * TWOWheeled_EV_B.Integrator;

    /* Sum: '<S26>/Sum8' */
    TWOWheeled_EV_B.Sum8 = TWOWheeled_EV_B.Kivnom + TWOWheeled_EV_B.Kaw;

    /* Sum: '<S14>/Add1' incorporates:
     *  Clock: '<S14>/Clock'
     */
    TWOWheeled_EV_B.Add1 = rtb_Add - TWOWheeled_EV_B.DigitalClock;
    if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
      /* Sum: '<S6>/Add' incorporates:
       *  Constant: '<S6>/Constant4'
       *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
       */
      rtb_Add = TWOWheeled_EV_P.Constant4_Value -
        TWOWheeled_EV_DW.DiscreteTimeIntegrator_DSTATE;
    }
  }

  if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
    /* Matfile logging */
    rt_UpdateTXYLogVars(TWOWheeled_EV_M->rtwLogInfo, (TWOWheeled_EV_M->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
    NeslSimulationData *simulationData;
    NeuDiagnosticManager *diagnosticManager;
    NeuDiagnosticTree *diagnosticTree;
    char *msg;
    real_T tmp_0[16];
    real_T time;
    int32_T tmp_2;
    int_T tmp_1[5];
    boolean_T tmp;

    /* Update for SimscapeExecutionBlock: '<S42>/STATE_1' */
    simulationData = (NeslSimulationData *)TWOWheeled_EV_DW.STATE_1_SimData;
    time = TWOWheeled_EV_M->Timing.t[0];
    simulationData->mData->mTime.mN = 1;
    simulationData->mData->mTime.mX = &time;
    simulationData->mData->mContStates.mN = 27;
    simulationData->mData->mContStates.mX =
      &TWOWheeled_EV_X.TWOWheeled_EVBatterycharge[0];
    simulationData->mData->mDiscStates.mN = 0;
    simulationData->mData->mDiscStates.mX = &TWOWheeled_EV_DW.STATE_1_Discrete;
    simulationData->mData->mModeVector.mN = 6;
    simulationData->mData->mModeVector.mX = (int32_T *)
      &TWOWheeled_EV_DW.STATE_1_Modes[0];
    tmp = false;
    simulationData->mData->mFoundZcEvents = tmp;
    simulationData->mData->mIsMajorTimeStep = rtmIsMajorTimeStep(TWOWheeled_EV_M);
    tmp = false;
    simulationData->mData->mIsSolverAssertCheck = tmp;
    simulationData->mData->mIsSolverCheckingCIC = false;
    tmp = rtsiIsSolverComputingJacobian(&TWOWheeled_EV_M->solverInfo);
    simulationData->mData->mIsComputingJacobian = tmp;
    simulationData->mData->mIsEvaluatingF0 = false;
    simulationData->mData->mIsSolverRequestingReset = false;
    simulationData->mData->mIsModeUpdateTimeStep = rtsiIsModeUpdateTimeStep
      (&TWOWheeled_EV_M->solverInfo);
    tmp_1[0] = 0;
    tmp_0[0] = TWOWheeled_EV_B.INPUT_1_1_1[0];
    tmp_0[1] = TWOWheeled_EV_B.INPUT_1_1_1[1];
    tmp_0[2] = TWOWheeled_EV_B.INPUT_1_1_1[2];
    tmp_0[3] = TWOWheeled_EV_B.INPUT_1_1_1[3];
    tmp_1[1] = 4;
    tmp_0[4] = TWOWheeled_EV_B.INPUT_2_1_1[0];
    tmp_0[5] = TWOWheeled_EV_B.INPUT_2_1_1[1];
    tmp_0[6] = TWOWheeled_EV_B.INPUT_2_1_1[2];
    tmp_0[7] = TWOWheeled_EV_B.INPUT_2_1_1[3];
    tmp_1[2] = 8;
    tmp_0[8] = TWOWheeled_EV_B.INPUT_4_1_1[0];
    tmp_0[9] = TWOWheeled_EV_B.INPUT_4_1_1[1];
    tmp_0[10] = TWOWheeled_EV_B.INPUT_4_1_1[2];
    tmp_0[11] = TWOWheeled_EV_B.INPUT_4_1_1[3];
    tmp_1[3] = 12;
    tmp_0[12] = TWOWheeled_EV_B.INPUT_3_1_1[0];
    tmp_0[13] = TWOWheeled_EV_B.INPUT_3_1_1[1];
    tmp_0[14] = TWOWheeled_EV_B.INPUT_3_1_1[2];
    tmp_0[15] = TWOWheeled_EV_B.INPUT_3_1_1[3];
    tmp_1[4] = 16;
    simulationData->mData->mInputValues.mN = 16;
    simulationData->mData->mInputValues.mX = &tmp_0[0];
    simulationData->mData->mInputOffsets.mN = 5;
    simulationData->mData->mInputOffsets.mX = (int32_T *)&tmp_1[0];
    diagnosticManager = (NeuDiagnosticManager *)TWOWheeled_EV_DW.STATE_1_DiagMgr;
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    tmp_2 = ne_simulator_method((NeslSimulator *)
      TWOWheeled_EV_DW.STATE_1_Simulator, NESL_SIM_UPDATE, simulationData,
      diagnosticManager);
    if (tmp_2 != 0L) {
      tmp = error_buffer_is_empty(rtmGetErrorStatus(TWOWheeled_EV_M));
      if (tmp) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(TWOWheeled_EV_M, msg);
      }
    }

    /* End of Update for SimscapeExecutionBlock: '<S42>/STATE_1' */
    if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
      /* Update for UnitDelay: '<S35>/Unit Delay' */
      TWOWheeled_EV_DW.UnitDelay_DSTATE[0] = TWOWheeled_EV_B.Switch[0];
      TWOWheeled_EV_DW.UnitDelay_DSTATE[1] = TWOWheeled_EV_B.Switch[1];

      /* Update for DiscreteIntegrator: '<S6>/Discrete-Time Integrator' */
      TWOWheeled_EV_DW.DiscreteTimeIntegrator_DSTATE +=
        TWOWheeled_EV_P.DiscreteTimeIntegrator_gainval * rtb_Gain1;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(TWOWheeled_EV_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(TWOWheeled_EV_M)!=-1) &&
          !((rtmGetTFinal(TWOWheeled_EV_M)-(((TWOWheeled_EV_M->Timing.clockTick1
               +TWOWheeled_EV_M->Timing.clockTickH1* 4294967296.0)) * 0.5)) >
            (((TWOWheeled_EV_M->Timing.clockTick1+
               TWOWheeled_EV_M->Timing.clockTickH1* 4294967296.0)) * 0.5) *
            (DBL_EPSILON))) {
        rtmSetErrorStatus(TWOWheeled_EV_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&TWOWheeled_EV_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++TWOWheeled_EV_M->Timing.clockTick0)) {
      ++TWOWheeled_EV_M->Timing.clockTickH0;
    }

    TWOWheeled_EV_M->Timing.t[0] = rtsiGetSolverStopTime
      (&TWOWheeled_EV_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.5s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.5, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      TWOWheeled_EV_M->Timing.clockTick1++;
      if (!TWOWheeled_EV_M->Timing.clockTick1) {
        TWOWheeled_EV_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void TWOWheeled_EV_derivatives(void)
{
  NeslSimulationData *simulationData;
  NeuDiagnosticManager *diagnosticManager;
  NeuDiagnosticTree *diagnosticTree;
  XDot_TWOWheeled_EV_T *_rtXdot;
  char *msg;
  real_T tmp[16];
  real_T time;
  int32_T tmp_1;
  int_T tmp_0[5];
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot_TWOWheeled_EV_T *) TWOWheeled_EV_M->derivs);

  /* Derivatives for Integrator: '<S34>/Integrator' */
  _rtXdot->Integrator_CSTATE = TWOWheeled_EV_B.Product_d;

  /* Derivatives for Integrator: '<S26>/Integrator1' */
  lsat = (TWOWheeled_EV_X.Integrator1_CSTATE <=
          TWOWheeled_EV_P.Integrator1_LowerSat);
  usat = (TWOWheeled_EV_X.Integrator1_CSTATE >=
          TWOWheeled_EV_P.Integrator1_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (TWOWheeled_EV_B.Sum8 > 0.0)) || (usat &&
       (TWOWheeled_EV_B.Sum8 < 0.0))) {
    _rtXdot->Integrator1_CSTATE = TWOWheeled_EV_B.Sum8;
  } else {
    /* in saturation */
    _rtXdot->Integrator1_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S26>/Integrator1' */

  /* Derivatives for SimscapeExecutionBlock: '<S42>/STATE_1' */
  simulationData = (NeslSimulationData *)TWOWheeled_EV_DW.STATE_1_SimData;
  time = TWOWheeled_EV_M->Timing.t[0];
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time;
  simulationData->mData->mContStates.mN = 27;
  simulationData->mData->mContStates.mX =
    &TWOWheeled_EV_X.TWOWheeled_EVBatterycharge[0];
  simulationData->mData->mDiscStates.mN = 0;
  simulationData->mData->mDiscStates.mX = &TWOWheeled_EV_DW.STATE_1_Discrete;
  simulationData->mData->mModeVector.mN = 6;
  simulationData->mData->mModeVector.mX = (int32_T *)
    &TWOWheeled_EV_DW.STATE_1_Modes[0];
  lsat = false;
  simulationData->mData->mFoundZcEvents = lsat;
  simulationData->mData->mIsMajorTimeStep = rtmIsMajorTimeStep(TWOWheeled_EV_M);
  lsat = false;
  simulationData->mData->mIsSolverAssertCheck = lsat;
  simulationData->mData->mIsSolverCheckingCIC = false;
  lsat = rtsiIsSolverComputingJacobian(&TWOWheeled_EV_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = lsat;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  simulationData->mData->mIsModeUpdateTimeStep = rtsiIsModeUpdateTimeStep
    (&TWOWheeled_EV_M->solverInfo);
  tmp_0[0] = 0;
  tmp[0] = TWOWheeled_EV_B.INPUT_1_1_1[0];
  tmp[1] = TWOWheeled_EV_B.INPUT_1_1_1[1];
  tmp[2] = TWOWheeled_EV_B.INPUT_1_1_1[2];
  tmp[3] = TWOWheeled_EV_B.INPUT_1_1_1[3];
  tmp_0[1] = 4;
  tmp[4] = TWOWheeled_EV_B.INPUT_2_1_1[0];
  tmp[5] = TWOWheeled_EV_B.INPUT_2_1_1[1];
  tmp[6] = TWOWheeled_EV_B.INPUT_2_1_1[2];
  tmp[7] = TWOWheeled_EV_B.INPUT_2_1_1[3];
  tmp_0[2] = 8;
  tmp[8] = TWOWheeled_EV_B.INPUT_4_1_1[0];
  tmp[9] = TWOWheeled_EV_B.INPUT_4_1_1[1];
  tmp[10] = TWOWheeled_EV_B.INPUT_4_1_1[2];
  tmp[11] = TWOWheeled_EV_B.INPUT_4_1_1[3];
  tmp_0[3] = 12;
  tmp[12] = TWOWheeled_EV_B.INPUT_3_1_1[0];
  tmp[13] = TWOWheeled_EV_B.INPUT_3_1_1[1];
  tmp[14] = TWOWheeled_EV_B.INPUT_3_1_1[2];
  tmp[15] = TWOWheeled_EV_B.INPUT_3_1_1[3];
  tmp_0[4] = 16;
  simulationData->mData->mInputValues.mN = 16;
  simulationData->mData->mInputValues.mX = &tmp[0];
  simulationData->mData->mInputOffsets.mN = 5;
  simulationData->mData->mInputOffsets.mX = (int32_T *)&tmp_0[0];
  simulationData->mData->mDx.mN = 27;
  simulationData->mData->mDx.mX = &_rtXdot->TWOWheeled_EVBatterycharge[0];
  diagnosticManager = (NeuDiagnosticManager *)TWOWheeled_EV_DW.STATE_1_DiagMgr;
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  tmp_1 = ne_simulator_method((NeslSimulator *)
    TWOWheeled_EV_DW.STATE_1_Simulator, NESL_SIM_DERIVATIVES, simulationData,
    diagnosticManager);
  if (tmp_1 != 0L) {
    lsat = error_buffer_is_empty(rtmGetErrorStatus(TWOWheeled_EV_M));
    if (lsat) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(TWOWheeled_EV_M, msg);
    }
  }

  /* End of Derivatives for SimscapeExecutionBlock: '<S42>/STATE_1' */

  /* Derivatives for Integrator: '<Root>/Integrator' */
  _rtXdot->Integrator_CSTATE_f = TWOWheeled_EV_B.OUTPUT_1_0[2];

  /* Derivatives for Integrator: '<S35>/Integrator2' */
  _rtXdot->Integrator2_CSTATE = TWOWheeled_EV_B.Product;
}

/* Model initialize function */
void TWOWheeled_EV_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)TWOWheeled_EV_M, 0,
                sizeof(RT_MODEL_TWOWheeled_EV_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&TWOWheeled_EV_M->solverInfo,
                          &TWOWheeled_EV_M->Timing.simTimeStep);
    rtsiSetTPtr(&TWOWheeled_EV_M->solverInfo, &rtmGetTPtr(TWOWheeled_EV_M));
    rtsiSetStepSizePtr(&TWOWheeled_EV_M->solverInfo,
                       &TWOWheeled_EV_M->Timing.stepSize0);
    rtsiSetdXPtr(&TWOWheeled_EV_M->solverInfo, &TWOWheeled_EV_M->derivs);
    rtsiSetContStatesPtr(&TWOWheeled_EV_M->solverInfo, (real_T **)
                         &TWOWheeled_EV_M->contStates);
    rtsiSetNumContStatesPtr(&TWOWheeled_EV_M->solverInfo,
      &TWOWheeled_EV_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&TWOWheeled_EV_M->solverInfo,
      &TWOWheeled_EV_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&TWOWheeled_EV_M->solverInfo,
      &TWOWheeled_EV_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&TWOWheeled_EV_M->solverInfo,
      &TWOWheeled_EV_M->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&TWOWheeled_EV_M->solverInfo, (boolean_T**)
      &TWOWheeled_EV_M->contStateDisabled);
    rtsiSetErrorStatusPtr(&TWOWheeled_EV_M->solverInfo, (&rtmGetErrorStatus
      (TWOWheeled_EV_M)));
    rtsiSetSolverMassMatrixIr(&TWOWheeled_EV_M->solverInfo,
      TWOWheeled_EV_MassMatrix.ir);
    rtsiSetSolverMassMatrixJc(&TWOWheeled_EV_M->solverInfo,
      TWOWheeled_EV_MassMatrix.jc);
    rtsiSetSolverMassMatrixPr(&TWOWheeled_EV_M->solverInfo,
      TWOWheeled_EV_MassMatrix.pr);
    rtsiSetRTModelPtr(&TWOWheeled_EV_M->solverInfo, TWOWheeled_EV_M);
  }

  rtsiSetSimTimeStep(&TWOWheeled_EV_M->solverInfo, MAJOR_TIME_STEP);
  TWOWheeled_EV_M->intgData.x0 = TWOWheeled_EV_M->odeX0;
  TWOWheeled_EV_M->intgData.f0 = TWOWheeled_EV_M->odeF0;
  TWOWheeled_EV_M->intgData.x1start = TWOWheeled_EV_M->odeX1START;
  TWOWheeled_EV_M->intgData.f1 = TWOWheeled_EV_M->odeF1;
  TWOWheeled_EV_M->intgData.Delta = TWOWheeled_EV_M->odeDELTA;
  TWOWheeled_EV_M->intgData.E = TWOWheeled_EV_M->odeE;
  TWOWheeled_EV_M->intgData.fac = TWOWheeled_EV_M->odeFAC;

  /* initialize */
  {
    int_T i;
    real_T *f = TWOWheeled_EV_M->intgData.fac;
    for (i = 0; i < (int_T)(sizeof(TWOWheeled_EV_M->odeFAC)/sizeof(real_T)); i++)
    {
      f[i] = 1.5e-8;
    }
  }

  TWOWheeled_EV_M->intgData.DFDX = TWOWheeled_EV_M->odeDFDX;
  TWOWheeled_EV_M->intgData.W = TWOWheeled_EV_M->odeW;
  TWOWheeled_EV_M->intgData.pivots = TWOWheeled_EV_M->odePIVOTS;
  TWOWheeled_EV_M->intgData.xtmp = TWOWheeled_EV_M->odeXTMP;
  TWOWheeled_EV_M->intgData.ztmp = TWOWheeled_EV_M->odeZTMP;
  TWOWheeled_EV_M->intgData.M = TWOWheeled_EV_M->odeMASSMATRIX_M;
  TWOWheeled_EV_M->intgData.isFirstStep = true;
  rtsiSetSolverExtrapolationOrder(&TWOWheeled_EV_M->solverInfo, 4);
  rtsiSetSolverNumberNewtonIterations(&TWOWheeled_EV_M->solverInfo, 1);
  TWOWheeled_EV_M->contStates = ((X_TWOWheeled_EV_T *) &TWOWheeled_EV_X);
  TWOWheeled_EV_M->contStateDisabled = ((XDis_TWOWheeled_EV_T *)
    &TWOWheeled_EV_XDis);
  TWOWheeled_EV_M->Timing.tStart = (0.0);
  TWOWheeled_EV_M->massMatrixType = ((ssMatrixType)1);
  TWOWheeled_EV_M->massMatrixNzMax = (12);
  TWOWheeled_EV_M->massMatrixIr = (TWOWheeled_EV_MassMatrix.ir);
  TWOWheeled_EV_M->massMatrixJc = (TWOWheeled_EV_MassMatrix.jc);
  TWOWheeled_EV_M->massMatrixPr = (TWOWheeled_EV_MassMatrix.pr);
  rtsiSetSolverMassMatrixType(&TWOWheeled_EV_M->solverInfo, (ssMatrixType)1);
  rtsiSetSolverMassMatrixNzMax(&TWOWheeled_EV_M->solverInfo, 12);
  rtsiSetSolverData(&TWOWheeled_EV_M->solverInfo, (void *)
                    &TWOWheeled_EV_M->intgData);
  rtsiSetIsMinorTimeStepWithModeChange(&TWOWheeled_EV_M->solverInfo, false);
  rtsiSetSolverName(&TWOWheeled_EV_M->solverInfo,"ode14x");
  rtmSetTPtr(TWOWheeled_EV_M, &TWOWheeled_EV_M->Timing.tArray[0]);
  rtmSetTFinal(TWOWheeled_EV_M, 200.0);
  TWOWheeled_EV_M->Timing.stepSize0 = 0.5;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = (NULL);
    TWOWheeled_EV_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(TWOWheeled_EV_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(TWOWheeled_EV_M->rtwLogInfo, (NULL));
    rtliSetLogT(TWOWheeled_EV_M->rtwLogInfo, "tout");
    rtliSetLogX(TWOWheeled_EV_M->rtwLogInfo, "");
    rtliSetLogXFinal(TWOWheeled_EV_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(TWOWheeled_EV_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(TWOWheeled_EV_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(TWOWheeled_EV_M->rtwLogInfo, 0);
    rtliSetLogDecimation(TWOWheeled_EV_M->rtwLogInfo, 1);
    rtliSetLogY(TWOWheeled_EV_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(TWOWheeled_EV_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(TWOWheeled_EV_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &TWOWheeled_EV_B), 0,
                sizeof(B_TWOWheeled_EV_T));

  /* states (continuous) */
  {
    (void) memset((void *)&TWOWheeled_EV_X, 0,
                  sizeof(X_TWOWheeled_EV_T));
  }

  /* disabled states */
  {
    (void) memset((void *)&TWOWheeled_EV_XDis, 0,
                  sizeof(XDis_TWOWheeled_EV_T));
  }

  /* global mass matrix */
  {
    int_T *ir = TWOWheeled_EV_MassMatrix.ir;
    int_T *jc = TWOWheeled_EV_MassMatrix.jc;
    real_T *pr = TWOWheeled_EV_MassMatrix.pr;
    (void) memset((void *)ir, 0,
                  12*sizeof(int_T));
    (void) memset((void *)jc, 0,
                  (31+1)*sizeof(int_T));
    (void) memset((void *)pr, 0,
                  12*sizeof(real_T));
  }

  /* states (dwork) */
  (void) memset((void *)&TWOWheeled_EV_DW, 0,
                sizeof(DW_TWOWheeled_EV_T));

  /* Root-level init GlobalMassMatrixPr offset */
  {
    TWOWheeled_EV_DW.STATE_1_MASS_MATRIX_PR = 2;/* '<S42>/STATE_1' */
  }

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(TWOWheeled_EV_M->rtwLogInfo, 0.0,
    rtmGetTFinal(TWOWheeled_EV_M), TWOWheeled_EV_M->Timing.stepSize0,
    (&rtmGetErrorStatus(TWOWheeled_EV_M)));

  {
    NeModelParameters modelParameters;
    NeModelParameters modelParameters_0;
    NeslSimulationData *tmp_1;
    NeslSimulator *tmp;
    NeuDiagnosticManager *diagnosticManager;
    NeuDiagnosticTree *diagnosticTree;
    NeuDiagnosticTree *diagnosticTree_0;
    char *msg;
    char *msg_0;
    real_T tmp_2;
    int32_T tmp_3;
    boolean_T tmp_0;

    /* Start for FromWorkspace: '<S7>/FromWs' */
    {
      static real_T pTimeValues0[] = { 0.0, 24.0, 24.0, 135.0, 135.0, 135.5,
        135.5, 200.0 } ;

      static real_T pDataValues0[] = { 0.0, 18.8, 18.8, 29.6, 29.6, 19.5, 19.5,
        29.0 } ;

      TWOWheeled_EV_DW.FromWs_PWORK.TimePtr = (void *) pTimeValues0;
      TWOWheeled_EV_DW.FromWs_PWORK.DataPtr = (void *) pDataValues0;
      TWOWheeled_EV_DW.FromWs_IWORK.PrevIndex = 0;
    }

    /* Start for SimscapeExecutionBlock: '<S42>/STATE_1' */
    tmp = nesl_lease_simulator("TWOWheeled_EV/Solver Configuration_1", 0L, 0L);
    TWOWheeled_EV_DW.STATE_1_Simulator = (void *)tmp;
    tmp_0 = pointer_is_null(TWOWheeled_EV_DW.STATE_1_Simulator);
    if (tmp_0) {
      TWOWheeled_EV_3cd5c335_1_gateway();
      tmp = nesl_lease_simulator("TWOWheeled_EV/Solver Configuration_1", 0L, 0L);
      TWOWheeled_EV_DW.STATE_1_Simulator = (void *)tmp;
    }

    tmp_1 = nesl_create_simulation_data();
    TWOWheeled_EV_DW.STATE_1_SimData = (void *)tmp_1;
    diagnosticManager = rtw_create_diagnostics();
    TWOWheeled_EV_DW.STATE_1_DiagMgr = (void *)diagnosticManager;
    modelParameters.mSolverType = NE_SOLVER_TYPE_DAE;
    modelParameters.mSolverAbsTol = 0.001;
    modelParameters.mSolverRelTol = 0.001;
    modelParameters.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters.mStartTime = 0.0;
    modelParameters.mLoadInitialState = false;
    modelParameters.mUseSimState = false;
    modelParameters.mLinTrimCompile = false;
    modelParameters.mLoggingMode = SSC_LOGGING_OFF;
    modelParameters.mRTWModifiedTimeStamp = 6.27033059E+8;
    modelParameters.mUseModelRefSolver = false;
    modelParameters.mTargetFPGAHIL = false;
    tmp_2 = 0.001;
    modelParameters.mSolverTolerance = tmp_2;
    tmp_2 = 0.5;
    modelParameters.mFixedStepSize = tmp_2;
    tmp_0 = false;
    modelParameters.mVariableStepSolver = tmp_0;
    tmp_0 = false;
    modelParameters.mIsUsingODEN = tmp_0;
    modelParameters.mZcDisabled = true;
    diagnosticManager = (NeuDiagnosticManager *)TWOWheeled_EV_DW.STATE_1_DiagMgr;
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    tmp_3 = nesl_initialize_simulator((NeslSimulator *)
      TWOWheeled_EV_DW.STATE_1_Simulator, &modelParameters, diagnosticManager);
    if (tmp_3 != 0L) {
      tmp_0 = error_buffer_is_empty(rtmGetErrorStatus(TWOWheeled_EV_M));
      if (tmp_0) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(TWOWheeled_EV_M, msg);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S42>/STATE_1' */

    /* Start for SimscapeExecutionBlock: '<S42>/OUTPUT_1_0' */
    tmp = nesl_lease_simulator("TWOWheeled_EV/Solver Configuration_1", 1L, 0L);
    TWOWheeled_EV_DW.OUTPUT_1_0_Simulator = (void *)tmp;
    tmp_0 = pointer_is_null(TWOWheeled_EV_DW.OUTPUT_1_0_Simulator);
    if (tmp_0) {
      TWOWheeled_EV_3cd5c335_1_gateway();
      tmp = nesl_lease_simulator("TWOWheeled_EV/Solver Configuration_1", 1L, 0L);
      TWOWheeled_EV_DW.OUTPUT_1_0_Simulator = (void *)tmp;
    }

    tmp_1 = nesl_create_simulation_data();
    TWOWheeled_EV_DW.OUTPUT_1_0_SimData = (void *)tmp_1;
    diagnosticManager = rtw_create_diagnostics();
    TWOWheeled_EV_DW.OUTPUT_1_0_DiagMgr = (void *)diagnosticManager;
    modelParameters_0.mSolverType = NE_SOLVER_TYPE_DAE;
    modelParameters_0.mSolverAbsTol = 0.001;
    modelParameters_0.mSolverRelTol = 0.001;
    modelParameters_0.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters_0.mStartTime = 0.0;
    modelParameters_0.mLoadInitialState = false;
    modelParameters_0.mUseSimState = false;
    modelParameters_0.mLinTrimCompile = false;
    modelParameters_0.mLoggingMode = SSC_LOGGING_OFF;
    modelParameters_0.mRTWModifiedTimeStamp = 6.27033059E+8;
    modelParameters_0.mUseModelRefSolver = false;
    modelParameters_0.mTargetFPGAHIL = false;
    tmp_2 = 0.001;
    modelParameters_0.mSolverTolerance = tmp_2;
    tmp_2 = 0.5;
    modelParameters_0.mFixedStepSize = tmp_2;
    tmp_0 = false;
    modelParameters_0.mVariableStepSolver = tmp_0;
    tmp_0 = false;
    modelParameters_0.mIsUsingODEN = tmp_0;
    modelParameters_0.mZcDisabled = true;
    diagnosticManager = (NeuDiagnosticManager *)
      TWOWheeled_EV_DW.OUTPUT_1_0_DiagMgr;
    diagnosticTree_0 = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    tmp_3 = nesl_initialize_simulator((NeslSimulator *)
      TWOWheeled_EV_DW.OUTPUT_1_0_Simulator, &modelParameters_0,
      diagnosticManager);
    if (tmp_3 != 0L) {
      tmp_0 = error_buffer_is_empty(rtmGetErrorStatus(TWOWheeled_EV_M));
      if (tmp_0) {
        msg_0 = rtw_diagnostics_msg(diagnosticTree_0);
        rtmSetErrorStatus(TWOWheeled_EV_M, msg_0);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S42>/OUTPUT_1_0' */
  }

  {
    int_T tmp_1;
    int_T tmp_2;
    int_T tmp_3;
    int_T tmp_4;
    int_T tmp_5;
    int_T tmp_6;
    boolean_T tmp;
    boolean_T tmp_0;

    /* InitializeConditions for Integrator: '<S34>/Integrator' */
    TWOWheeled_EV_X.Integrator_CSTATE = TWOWheeled_EV_P.Integrator_IC;

    /* InitializeConditions for Integrator: '<S26>/Integrator1' */
    TWOWheeled_EV_X.Integrator1_CSTATE = TWOWheeled_EV_P.Integrator1_IC;

    /* InitializeConditions for SimscapeExecutionBlock: '<S42>/STATE_1' */
    tmp = false;
    tmp_0 = false;
    if (tmp_0 || tmp) {
      tmp_1 = strcmp(rtsiGetSolverName(&TWOWheeled_EV_M->solverInfo), "daessc");
      tmp_2 = strcmp(rtsiGetSolverName(&TWOWheeled_EV_M->solverInfo), "ode14x");
      tmp_3 = strcmp(rtsiGetSolverName(&TWOWheeled_EV_M->solverInfo), "ode15s");
      tmp_4 = strcmp(rtsiGetSolverName(&TWOWheeled_EV_M->solverInfo), "ode1be");
      tmp_5 = strcmp(rtsiGetSolverName(&TWOWheeled_EV_M->solverInfo), "ode23t");
      tmp_6 = strcmp(rtsiGetSolverName(&TWOWheeled_EV_M->solverInfo), "odeN");
      if ((boolean_T)((tmp_1 != 0) & (tmp_2 != 0) & (tmp_3 != 0) & (tmp_4 != 0)
                      & (tmp_5 != 0) & (tmp_6 != 0))) {
        rtmSetErrorStatus(TWOWheeled_EV_M,
                          "Detected inconsistent solvers in the model reference hierarchy. Model built with ode14x requires one of {daessc, ode14x, ode15s, ode1be, ode23t, odeN} solvers to run. Use one of the required solvers in the top model.");
      }
    }

    /* End of InitializeConditions for SimscapeExecutionBlock: '<S42>/STATE_1' */

    /* InitializeConditions for UnitDelay: '<S35>/Unit Delay' */
    TWOWheeled_EV_DW.UnitDelay_DSTATE[0] =
      TWOWheeled_EV_P.UnitDelay_InitialCondition[0];
    TWOWheeled_EV_DW.UnitDelay_DSTATE[1] =
      TWOWheeled_EV_P.UnitDelay_InitialCondition[1];

    /* InitializeConditions for Integrator: '<Root>/Integrator' */
    TWOWheeled_EV_X.Integrator_CSTATE_f = TWOWheeled_EV_P.Integrator_IC_m;

    /* InitializeConditions for Integrator: '<S35>/Integrator2' */
    TWOWheeled_EV_X.Integrator2_CSTATE = TWOWheeled_EV_P.Integrator2_IC;

    /* InitializeConditions for DiscreteIntegrator: '<S6>/Discrete-Time Integrator' */
    TWOWheeled_EV_DW.DiscreteTimeIntegrator_DSTATE =
      TWOWheeled_EV_P.DiscreteTimeIntegrator_IC;

    /* SystemInitialize for Enabled SubSystem: '<S29>/Pass Through' */
    TWOWheeled_EV_PassThrough_Init(TWOWheeled_EV_P.SignalHold_IC,
      &TWOWheeled_EV_B.PassThrough);

    /* End of SystemInitialize for SubSystem: '<S29>/Pass Through' */

    /* SystemInitialize for Enabled SubSystem: '<S31>/Pass Through' */
    TWOWheeled_EV_PassThrough_Init(TWOWheeled_EV_P.SignalHold_IC_d,
      &TWOWheeled_EV_B.PassThrough_h);

    /* End of SystemInitialize for SubSystem: '<S31>/Pass Through' */

    /* Root-level InitSystemMatrices */
    {
      static int_T modelMassMatrixIr[12] = { 0, 1, 2, 4, 3, 9, 5, 6, 9, 7, 29,
        30 };

      static int_T modelMassMatrixJc[32] = { 0, 1, 2, 3, 4, 6, 7, 9, 10, 10, 10,
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
        10, 11, 12 };

      static real_T modelMassMatrixPr[12] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0 };

      (void) memcpy(TWOWheeled_EV_MassMatrix.ir, modelMassMatrixIr,
                    12*sizeof(int_T));
      (void) memcpy(TWOWheeled_EV_MassMatrix.jc, modelMassMatrixJc,
                    32*sizeof(int_T));
      (void) memcpy(TWOWheeled_EV_MassMatrix.pr, modelMassMatrixPr,
                    12*sizeof(real_T));
    }
  }
}

/* Model terminate function */
void TWOWheeled_EV_terminate(void)
{
  /* Terminate for SimscapeExecutionBlock: '<S42>/STATE_1' */
  neu_destroy_diagnostic_manager((NeuDiagnosticManager *)
    TWOWheeled_EV_DW.STATE_1_DiagMgr);
  nesl_destroy_simulation_data((NeslSimulationData *)
    TWOWheeled_EV_DW.STATE_1_SimData);
  nesl_erase_simulator("TWOWheeled_EV/Solver Configuration_1");
  nesl_destroy_registry();

  /* Terminate for SimscapeExecutionBlock: '<S42>/OUTPUT_1_0' */
  neu_destroy_diagnostic_manager((NeuDiagnosticManager *)
    TWOWheeled_EV_DW.OUTPUT_1_0_DiagMgr);
  nesl_destroy_simulation_data((NeslSimulationData *)
    TWOWheeled_EV_DW.OUTPUT_1_0_SimData);
  nesl_erase_simulator("TWOWheeled_EV/Solver Configuration_1");
  nesl_destroy_registry();
}
