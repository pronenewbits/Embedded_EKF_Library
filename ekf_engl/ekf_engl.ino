#include <Wire.h>
#include <elapsedMillis.h>
#include "konfig.h"
#include "matrix.h"
#include "ekf.h"


/* Just example */
#define P_INIT      (1000.)
#define Q_INIT      (0.001)
#define R_INIT      (0.001)


bool Main_bUpdateNonlinearX(Matrix &X_Next, Matrix &X, Matrix &U);
bool Main_bUpdateNonlinearY(Matrix &Y, Matrix &X, Matrix &U);
bool Main_bCalcJacobianF(Matrix &F, Matrix &X, Matrix &U);
bool Main_bCalcJacobianH(Matrix &H, Matrix &X, Matrix &U);

Matrix X(SS_X_LEN, 1);
Matrix Y(SS_Z_LEN, 1);
Matrix U(SS_U_LEN, 1);
EKF EKF_IMU(X, P_INIT, Q_INIT, R_INIT, Main_bUpdateNonlinearX, Main_bUpdateNonlinearY, Main_bCalcJacobianF, Main_bCalcJacobianH);

elapsedMillis timerLed, timerEKF;
uint64_t u64compuTime;
char bufferTxSer[100];



void setup() {
    /* serial to display data */
    Serial.begin(115200);
    while(!Serial) {}
    
    X.vIsiNol();
    EKF_IMU.vReset(X, P_INIT, Q_INIT, R_INIT);
}


void loop() {
    if (timerEKF > SS_DT_MILIS) {
        /* ================== Read the sensor data / simulate the system here ================== */
        /* ------------------ Read the sensor data / simulate the system here ------------------ */
        
        
        
        /* ============================= Update the Kalman Filter ============================== */
        u64compuTime = micros();
        if (!EKF_IMU.bUpdate(Y, U)) {
            X.vIsiNol();
            EKF_IMU.vReset(X, P_INIT, Q_INIT, R_INIT);
            Serial.println("Whoop ");
        }
        u64compuTime = (micros() - u64compuTime);
        /* ----------------------------- Update the Kalman Filter ------------------------------ */
        
        
        
        /* =========================== Print to serial (for plotting) ========================== */
        #if (1)
            /* Print: Computation time, X_Est[0] */
            snprintf(bufferTxSer, sizeof(bufferTxSer)-1, "%.3f %.3f ", ((float)u64compuTime)/1000., EKF_IMU.GetX()[0][0]);
            Serial.print(bufferTxSer);
        #endif
        Serial.print('\n');
        /* --------------------------- Print to serial (for plotting) -------------------------- */
        
        
        timerEKF = 0;
    }
}


bool Main_bUpdateNonlinearX(Matrix &X_Next, Matrix &X, Matrix &U)
{
    /* Insert the nonlinear update transformation here
     *          x(k+1) = f[x(k), u(k)]
     */
    
    return true;
}

bool Main_bUpdateNonlinearY(Matrix &Y, Matrix &X, Matrix &U)
{
    /* Insert the nonlinear measurement transformation here
     *          y(k)   = h[x(k), u(k)]
     */
    
    return true;
}

bool Main_bCalcJacobianF(Matrix &F, Matrix &X, Matrix &U)
{
    /* Insert the linearized update transformation here (i.e. Jacobian matrix of f[x(k), u(k)]) */
    
    return true;
}

bool Main_bCalcJacobianH(Matrix &H, Matrix &X, Matrix &U)
{
    /* Insert the linearized measurement transformation here (i.e. Jacobian matrix of h[x(k), u(k)]) */
    
    return true;
}





void SPEW_THE_ERROR(char const * str)
{
    #if (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_PC)
        cout << (str) << endl;
    #elif (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_EMBEDDED_ARDUINO)
        Serial.println(str);
    #else
        /* Silent function */
    #endif
    while(1);
}
