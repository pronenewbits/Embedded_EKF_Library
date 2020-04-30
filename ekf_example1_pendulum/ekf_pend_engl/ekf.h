/*************************************************************************************************************
 * Class for Discrete Extended Kalman Filter
 * 
 * 
 * See https://github.com/pronenewbits for more!
 ************************************************************************************************************/
#ifndef EKF_H
#define EKF_H

#include "konfig.h"
#include "matrix.h"

class EKF
{
public:
    EKF(const Matrix& XInit, const Matrix& P, const Matrix& Q, const Matrix& R,
        bool (*bNonlinearUpdateX)(Matrix&, const Matrix&, const Matrix&),
        bool (*bNonlinearUpdateY)(Matrix&, const Matrix&, const Matrix&), 
        bool (*bCalcJacobianF)(Matrix&, const Matrix&, const Matrix&),
        bool (*bCalcJacobianH)(Matrix&, const Matrix&, const Matrix&));
    void vReset(const Matrix& XInit, const Matrix& P, const Matrix& Q, const Matrix& R);
    bool bUpdate(const Matrix& Y, const Matrix& U);
    const Matrix GetX()   const { return X_Est; }
    const Matrix GetY()   const { return Y_Est; }
    const Matrix GetP()   const { return P; }
    const Matrix GetErr() const { return Err; }

protected:
    bool (*bNonlinearUpdateX) (Matrix& X_dot, const Matrix& X, const Matrix& U);
    bool (*bNonlinearUpdateY) (Matrix& Y_Est, const Matrix& X, const Matrix& U);
    bool (*bCalcJacobianF) (Matrix& F, const Matrix& X, const Matrix& U);
    bool (*bCalcJacobianH) (Matrix& H, const Matrix& X, const Matrix& U);

private:
    Matrix X_Est{SS_X_LEN, 1};
    Matrix P{SS_X_LEN, SS_X_LEN};
    Matrix F{SS_X_LEN, SS_X_LEN};
    Matrix H{SS_Z_LEN, SS_X_LEN};
    Matrix Y_Est{SS_Z_LEN, 1};
    Matrix Err{SS_Z_LEN, 1};
    Matrix Q{SS_X_LEN, SS_X_LEN};
    Matrix R{SS_Z_LEN, SS_Z_LEN};
    Matrix S{SS_Z_LEN, SS_Z_LEN};
    Matrix Gain{SS_X_LEN, SS_Z_LEN};
};

#endif // EKF_H
