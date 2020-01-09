/*********************************************************************************************
 *  Fungsi Discrete Extended Kalman Filter (yang dihitung secara diskrit)
 *
 *         Formulasi plant yang diestimasi:
 *              x(k+1) = f[x(k), u(k)] + w(t)          ; x=Nx1,    u=Mx1            ...{EKF_1}
 *              y(k)   = h[x(k), u(k)] + v(t)          ; z=Zx1                      ...{EKF_2}
 *
 *
 *
 * bla bla bla lorem ipsum bla bla bla
 * 
 * 
 * See https://github.com/pronenewbits for more!
 ********************************************************************************************/

#include "ekf.h"



EKF::EKF(Matrix &XInit, const float_prec PInit, const float_prec QInit, const float_prec RInit,
               bool (*bNonlinearUpdateX)(Matrix &, Matrix &, Matrix &), bool (*bNonlinearUpdateY)(Matrix &, Matrix &, Matrix &), 
               bool (*bCalcJacobianA)(Matrix &, Matrix &, Matrix &), bool (*bCalcJacobianC)(Matrix &, Matrix &, Matrix &))
{
    /* Initialization:
     *  x(k=0|k=0)  = Expected value of x at time-0 (i.e. x(k=0)), typically set to zero.
     *  P (k=0|k=0) = Identity * covariant(P(k=0)), typically initialized with some big number.
     *  Q, R        = Covariance matrices of process & measurement. As this implementation 
     *                the noise as AWGN (and same value for every variable), this is set
     *                to Q=diag(QInit,...,QInit) and R=diag(RInit,...,RInit).
     */
    X_Est = XInit;
    P.vIsiDiagonal(PInit);
    Q.vIsiDiagonal(QInit);
    R.vIsiDiagonal(RInit);
    this->bNonlinearUpdateX = bNonlinearUpdateX;
    this->bNonlinearUpdateY = bNonlinearUpdateY;
    this->bCalcJacobianA = bCalcJacobianA;
    this->bCalcJacobianC = bCalcJacobianC;
}

void EKF::vReset(Matrix &XInit, const float_prec PInit, const float_prec QInit, const float_prec RInit)
{
    X_Est = XInit;
    P.vIsiDiagonal(PInit);
    Q.vIsiDiagonal(QInit);
    R.vIsiDiagonal(RInit);
}

bool EKF::bUpdate(Matrix &Y, Matrix &U)
{
    /* Run once every sampling time */
    
    
    /* =============== Calculate the Jacobian matrix of f (i.e. A) =============== */
    /* A = d(f(..))/dx |x(k-1|k-1),u(k)                                 ...{EKF_3} */
    if (!bCalcJacobianA(A, X_Est, U)) {
        return false;
    }
    
    
    /* =========================== Prediction of x & P =========================== */
    /* x(k|k-1) = f[x(k-1|k-1), u(k)]                                   ...{EKF_4} */
    if (!bNonlinearUpdateX(X_Est, X_Est, U)) {
        return false;
    }

    /* P(k|k-1)  = A*P(k-1|k-1)*A' + Q                                  ...{EKF_5} */
    P = A*P*(A.Transpose()) + Q;
    
    
    
    
    
    /* =============== Calculate the Jacobian matrix of h (i.e. C) =============== */
    /* C = d(h(..))/dx |x(k|k-1),u(k)                                   ...{EKF_8} */
    if (!bCalcJacobianC(C, X_Est, U)) {
        return false;
    }
    
    
    /* =========================== Correction of x & P =========================== */
    /* e(k)    = y(k) - h[x(k|k-1), u(k)]                               ...{EKF_9} */
    if (!bNonlinearUpdateY(Y_Est, X_Est, U)) {
        return false;
    }
    Err = Y - Y_Est;

    /* S       = C*P(k|k-1)*C' + R                                      ...{EKF_10} */
    S = (C*P*(C.Transpose())) + R;

    /* K       = P(k|k-1)*C'*(S^-1)                                     ...{EKF_11} */
    Gain = P*(C.Transpose())*(S.Invers());

    /* x(k|k) = x(k|k-1) + K*e(k)                                       ...{EKF_12} */
    X_Est = X_Est + (Gain*Err);

    /* P(k|k)  = (I - K*C)*P(k|k-1)                                     ...{EKF_13} */
    Matrix Identitas = Matrix(SS_X_LEN, SS_X_LEN);
    Identitas.vSetIdentitas();
    P = (Identitas - (Gain*C))*P;
    
    
    return true;
}
