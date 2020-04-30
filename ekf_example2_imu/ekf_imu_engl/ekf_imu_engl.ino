/*************************************************************************************************************
 *  
 * In this example, we'll process IMU data from MPU9250 sensor.
 * 
 * See https://github.com/pronenewbits for more!
 ************************************************************************************************************/
#include <Wire.h>
#include <elapsedMillis.h>
#include "konfig.h"
#include "matrix.h"
#include "ekf.h"
#include "simple_mpu9250.h"


/* ================================================== The AHRS/IMU variables ================================================== */
/* Gravity vector constant (align with global Z-axis) */
#define IMU_ACC_Z0          (1)
/* Magnetic vector constant (align with local magnetic vector) */
float_prec IMU_MAG_B0_data[3] = {cos(0), sin(0), 0.000000};
Matrix IMU_MAG_B0(3, 1, IMU_MAG_B0_data);
/* The hard-magnet bias */
float_prec HARD_IRON_BIAS_data[3] = {8.832973, 7.243323, 23.95714};
Matrix HARD_IRON_BIAS(3, 1, HARD_IRON_BIAS_data);



/* ============================================ EKF variables/function declaration ============================================ */
/* Just example; in konfig.h: 
 *  SS_X_LEN = 4
 *  SS_Z_LEN = 6
 *  SS_U_LEN = 3
 */
/* EKF initialization constant -------------------------------------------------------------------------------------- */
#define P_INIT      (10.)
#define Q_INIT      (1e-6)
#define R_INIT_ACC  (0.0015/10.)
#define R_INIT_MAG  (0.0015/10.)
/* P(k=0) variable -------------------------------------------------------------------------------------------------- */
float_prec EKF_PINIT_data[SS_X_LEN*SS_X_LEN] = {P_INIT, 0,      0,      0,
                                                0,      P_INIT, 0,      0,
                                                0,      0,      P_INIT, 0,
                                                0,      0,      0,      P_INIT};
Matrix EKF_PINIT(SS_X_LEN, SS_X_LEN, EKF_PINIT_data);
/* Q constant ------------------------------------------------------------------------------------------------------- */
float_prec EKF_QINIT_data[SS_X_LEN*SS_X_LEN] = {Q_INIT, 0,      0,      0,
                                                0,      Q_INIT, 0,      0,
                                                0,      0,      Q_INIT, 0,
                                                0,      0,      0,      Q_INIT};
Matrix EKF_QINIT(SS_X_LEN, SS_X_LEN, EKF_QINIT_data);
/* R constant ------------------------------------------------------------------------------------------------------- */
float_prec EKF_RINIT_data[SS_Z_LEN*SS_Z_LEN] = {R_INIT_ACC, 0,          0,          0,          0,          0,
                                                0,          R_INIT_ACC, 0,          0,          0,          0,
                                                0,          0,          R_INIT_ACC, 0,          0,          0,
                                                0,          0,          0,          R_INIT_MAG, 0,          0,
                                                0,          0,          0,          0,          R_INIT_MAG, 0,
                                                0,          0,          0,          0,          0,          R_INIT_MAG};
Matrix EKF_RINIT(SS_Z_LEN, SS_Z_LEN, EKF_RINIT_data);
/* Nonlinear & linearization function ------------------------------------------------------------------------------- */
bool Main_bUpdateNonlinearX(Matrix& X_Next, const Matrix& X, const Matrix& U);
bool Main_bUpdateNonlinearY(Matrix& Y, const Matrix& X, const Matrix& U);
bool Main_bCalcJacobianF(Matrix& F, const Matrix& X, const Matrix& U);
bool Main_bCalcJacobianH(Matrix& H, const Matrix& X, const Matrix& U);
/* EKF variables ---------------------------------------------------------------------------------------------------- */
Matrix quaternionData(SS_X_LEN, 1);
Matrix Y(SS_Z_LEN, 1);
Matrix U(SS_U_LEN, 1);
/* EKF system declaration ------------------------------------------------------------------------------------------- */
EKF EKF_IMU(quaternionData, EKF_PINIT, EKF_QINIT, EKF_RINIT,
            Main_bUpdateNonlinearX, Main_bUpdateNonlinearY, Main_bCalcJacobianF, Main_bCalcJacobianH);



/* ========================================= Auxiliary variables/function declaration ========================================= */
elapsedMillis timerLed, timerEKF;
uint64_t u64compuTime;
char bufferTxSer[100];
/* The command from the PC */
char cmd;
/* A SimpleMPU9250 object with the MPU-9250 sensor on I2C bus 0 with address 0x68 */
SimpleMPU9250 IMU(Wire, 0x68);



void setup() {
    /* Serial initialization -------------------------------------- */
    Serial.begin(115200);
    while(!Serial) {}
    Serial.println("Calibrating IMU bias...");

    /* IMU initialization ----------------------------------------- */
    int status = IMU.begin();   /* start communication with IMU */
    if (status < 0) {
        Serial.println("IMU initialization unsuccessful");
        Serial.println("Check IMU wiring or try cycling power");
        Serial.print("Status: ");
        Serial.println(status);
        while(1) {}
    }
    
    /* EKF initialization ----------------------------------------- */
    /* x(k=0) = [1 0 0 0]' */
    quaternionData.vSetToZero();
    quaternionData[0][0] = 1.0;
    EKF_IMU.vReset(quaternionData, EKF_PINIT, EKF_QINIT, EKF_RINIT);
    
    snprintf(bufferTxSer, sizeof(bufferTxSer)-1, "EKF in Teensy 4.0 (%s)\r\n",
                                                 (FPU_PRECISION == PRECISION_SINGLE)?"Float32":"Double64");
    Serial.print(bufferTxSer);
}


void loop() {
    if (timerEKF >= SS_DT_MILIS) {
        timerEKF = 0;
        
        
        /* ================== Read the sensor data / simulate the system here ================== */
        /* Read the raw data */
        IMU.readSensor();
        float Ax = IMU.getAccelX_mss();
        float Ay = IMU.getAccelY_mss();
        float Az = IMU.getAccelZ_mss();
        float Bx = IMU.getMagX_uT();
        float By = IMU.getMagY_uT();
        float Bz = IMU.getMagZ_uT();
        float p = IMU.getGyroX_rads();
        float q = IMU.getGyroY_rads();
        float r = IMU.getGyroZ_rads();
        /* Input 1:3 = gyroscope */
        U[0][0] = p;  U[1][0] = q;  U[2][0] = r;
        /* Output 1:3 = accelerometer */
        Y[0][0] = Ax; Y[1][0] = Ay; Y[2][0] = Az;
        /* Output 4:6 = magnetometer */
        Y[3][0] = Bx; Y[4][0] = By; Y[5][0] = Bz;
        
        /* Compensating Hard-Iron Bias for magnetometer */
        Y[3][0] = Y[3][0]-HARD_IRON_BIAS[0][0];
        Y[4][0] = Y[4][0]-HARD_IRON_BIAS[1][0];
        Y[5][0] = Y[5][0]-HARD_IRON_BIAS[2][0];
        
        /* Normalizing the output vector */
        float_prec _normG = sqrt(Y[0][0] * Y[0][0]) + (Y[1][0] * Y[1][0]) + (Y[2][0] * Y[2][0]);
        Y[0][0] = Y[0][0] / _normG;
        Y[1][0] = Y[1][0] / _normG;
        Y[2][0] = Y[2][0] / _normG;
        float_prec _normM = sqrt(Y[3][0] * Y[3][0]) + (Y[4][0] * Y[4][0]) + (Y[5][0] * Y[5][0]);
        Y[3][0] = Y[3][0] / _normM;
        Y[4][0] = Y[4][0] / _normM;
        Y[5][0] = Y[5][0] / _normM;
        /* ------------------ Read the sensor data / simulate the system here ------------------ */
        
        
        /* ============================= Update the Kalman Filter ============================== */
        u64compuTime = micros();
        if (!EKF_IMU.bUpdate(Y, U)) {
            quaternionData.vSetToZero();
            quaternionData[0][0] = 1.0;
            EKF_IMU.vReset(quaternionData, EKF_PINIT, EKF_QINIT, EKF_RINIT);
            Serial.println("Whoop ");
        }
        u64compuTime = (micros() - u64compuTime);
        /* ----------------------------- Update the Kalman Filter ------------------------------ */
    }
    
    
    /* The serial data is sent by responding to command from the PC running Processing scipt */
    if (Serial.available()) {
        cmd = Serial.read();
        if (cmd == 'v') {
            snprintf(bufferTxSer, sizeof(bufferTxSer)-1, "EKF in Teensy 4.0 (%s)\r\n",
                                                         (FPU_PRECISION == PRECISION_SINGLE)?"Float32":"Double64");
            Serial.print(bufferTxSer);
            Serial.print('\n');
        } else if (cmd == 'q') {
            /* =========================== Print to serial (for plotting) ========================== */
            quaternionData = EKF_IMU.GetX();

            while (!Serial.available());
            uint8_t count = Serial.read();
            for (uint8_t _i = 0; _i < count; _i++) {
                serialFloatPrint(quaternionData[0][0]);
                Serial.print(",");
                serialFloatPrint(quaternionData[1][0]);
                Serial.print(",");
                serialFloatPrint(quaternionData[2][0]);
                Serial.print(",");
                serialFloatPrint(quaternionData[3][0]);
                Serial.print(",");
                serialFloatPrint((float)u64compuTime);
                Serial.print(",");
                Serial.println("");
            }
            /* --------------------------- Print to serial (for plotting) -------------------------- */
        }
    }
}


/* Function to interface with the Processing script in the PC */
void serialFloatPrint(float f) {
    byte * b = (byte *) &f;
    for (int i = 0; i < 4; i++) {
        byte b1 = (b[i] >> 4) & 0x0f;
        byte b2 = (b[i] & 0x0f);

        char c1 = (b1 < 10) ? ('0' + b1) : 'A' + b1 - 10;
        char c2 = (b2 < 10) ? ('0' + b2) : 'A' + b2 - 10;

        Serial.print(c1);
        Serial.print(c2);
    }
}


bool Main_bUpdateNonlinearX(Matrix& X_Next, const Matrix& X, const Matrix& U)
{
    /* Insert the nonlinear update transformation here
     *          x(k+1) = f[x(k), u(k)]
     *
     * The quaternion update function:
     *  q0_dot = 1/2. * (  0   - p*q1 - q*q2 - r*q3)
     *  q1_dot = 1/2. * ( p*q0 +   0  + r*q2 - q*q3)
     *  q2_dot = 1/2. * ( q*q0 - r*q1 +  0   + p*q3)
     *  q3_dot = 1/2. * ( r*q0 + q*q1 - p*q2 +  0  )
     * 
     * Euler method for integration:
     *  q0 = q0 + q0_dot * dT;
     *  q1 = q1 + q1_dot * dT;
     *  q2 = q2 + q2_dot * dT;
     *  q3 = q3 + q3_dot * dT;
     */
    float_prec q0, q1, q2, q3;
    float_prec p, q, r;
    
    q0 = X[0][0];
    q1 = X[1][0];
    q2 = X[2][0];
    q3 = X[3][0];
    
    p = U[0][0];
    q = U[1][0];
    r = U[2][0];
    
    X_Next[0][0] = (0.5 * (+0.00 -p*q1 -q*q2 -r*q3))*SS_DT + q0;
    X_Next[1][0] = (0.5 * (+p*q0 +0.00 +r*q2 -q*q3))*SS_DT + q1;
    X_Next[2][0] = (0.5 * (+q*q0 -r*q1 +0.00 +p*q3))*SS_DT + q2;
    X_Next[3][0] = (0.5 * (+r*q0 +q*q1 -p*q2 +0.00))*SS_DT + q3;
    
    
    /* ======= Additional ad-hoc quaternion normalization to make sure the quaternion is a unit vector (i.e. ||q|| = 1) ======= */
    if (!X_Next.bNormVector()) {
        /* System error, return false so we can reset the UKF */
        return false;
    }
    
    return true;
}

bool Main_bUpdateNonlinearY(Matrix& Y, const Matrix& X, const Matrix& U)
{
    /* Insert the nonlinear measurement transformation here
     *          y(k)   = h[x(k), u(k)]
     *
     * The measurement output is the gravitational and magnetic projection to the body:
     *     DCM     = [(+(q0**2)+(q1**2)-(q2**2)-(q3**2)),                    2*(q1*q2+q0*q3),                    2*(q1*q3-q0*q2)]
     *               [                   2*(q1*q2-q0*q3), (+(q0**2)-(q1**2)+(q2**2)-(q3**2)),                    2*(q2*q3+q0*q1)]
     *               [                   2*(q1*q3+q0*q2),                    2*(q2*q3-q0*q1), (+(q0**2)-(q1**2)-(q2**2)+(q3**2))]
     * 
     *  G_proj_sens = DCM * [0 0 1]             --> Gravitational projection to the accelerometer sensor
     *  M_proj_sens = DCM * [Mx My Mz]          --> (Earth) magnetic projection to the magnetometer sensor
     */
    float_prec q0, q1, q2, q3;
    float_prec q0_2, q1_2, q2_2, q3_2;

    q0 = X[0][0];
    q1 = X[1][0];
    q2 = X[2][0];
    q3 = X[3][0];

    q0_2 = q0 * q0;
    q1_2 = q1 * q1;
    q2_2 = q2 * q2;
    q3_2 = q3 * q3;
    
    Y[0][0] = (2*q1*q3 -2*q0*q2) * IMU_ACC_Z0;

    Y[1][0] = (2*q2*q3 +2*q0*q1) * IMU_ACC_Z0;

    Y[2][0] = (+(q0_2) -(q1_2) -(q2_2) +(q3_2)) * IMU_ACC_Z0;
    
    Y[3][0] = (+(q0_2)+(q1_2)-(q2_2)-(q3_2)) * IMU_MAG_B0[0][0]
             +(2*(q1*q2+q0*q3)) * IMU_MAG_B0[1][0]
             +(2*(q1*q3-q0*q2)) * IMU_MAG_B0[2][0];

    Y[4][0] = (2*(q1*q2-q0*q3)) * IMU_MAG_B0[0][0]
             +(+(q0_2)-(q1_2)+(q2_2)-(q3_2)) * IMU_MAG_B0[1][0]
             +(2*(q2*q3+q0*q1)) * IMU_MAG_B0[2][0];

    Y[5][0] = (2*(q1*q3+q0*q2)) * IMU_MAG_B0[0][0]
             +(2*(q2*q3-q0*q1)) * IMU_MAG_B0[1][0]
             +(+(q0_2)-(q1_2)-(q2_2)+(q3_2)) * IMU_MAG_B0[2][0];
    
    return true;
}

bool Main_bCalcJacobianF(Matrix& F, const Matrix& X, const Matrix& U)
{
    /* In Main_bUpdateNonlinearX():
     *  q0 = q0 + q0_dot * dT;
     *  q1 = q1 + q1_dot * dT;
     *  q2 = q2 + q2_dot * dT;
     *  q3 = q3 + q3_dot * dT;
     */
    float_prec p, q, r;

    p = U[0][0];
    q = U[1][0];
    r = U[2][0];

    F[0][0] =  1.000;
    F[1][0] =  0.5*p * SS_DT;
    F[2][0] =  0.5*q * SS_DT;
    F[3][0] =  0.5*r * SS_DT;

    F[0][1] = -0.5*p * SS_DT;
    F[1][1] =  1.000;
    F[2][1] = -0.5*r * SS_DT;
    F[3][1] =  0.5*q * SS_DT;

    F[0][2] = -0.5*q * SS_DT;
    F[1][2] =  0.5*r * SS_DT;
    F[2][2] =  1.000;
    F[3][2] = -0.5*p * SS_DT;

    F[0][3] = -0.5*r * SS_DT;
    F[1][3] = -0.5*q * SS_DT;
    F[2][3] =  0.5*p * SS_DT;
    F[3][3] =  1.000;
    
    return true;
}

bool Main_bCalcJacobianH(Matrix& H, const Matrix& X, const Matrix& U)
{
    /* In Main_bUpdateNonlinearY():
     * 
     * The measurement output is the gravitational and magnetic projection to the body:
     *     DCM     = [(+(q0**2)+(q1**2)-(q2**2)-(q3**2)),                    2*(q1*q2+q0*q3),                    2*(q1*q3-q0*q2)]
     *               [                   2*(q1*q2-q0*q3), (+(q0**2)-(q1**2)+(q2**2)-(q3**2)),                    2*(q2*q3+q0*q1)]
     *               [                   2*(q1*q3+q0*q2),                    2*(q2*q3-q0*q1), (+(q0**2)-(q1**2)-(q2**2)+(q3**2))]
     * 
     *  G_proj_sens = DCM * [0 0 -g]            --> Gravitational projection to the accelerometer sensor
     *  M_proj_sens = DCM * [Mx My Mz]          --> (Earth) magnetic projection to the magnetometer sensor
     */
    float_prec q0, q1, q2, q3;

    q0 = X[0][0];
    q1 = X[1][0];
    q2 = X[2][0];
    q3 = X[3][0];
    
    H[0][0] = -2*q2 * IMU_ACC_Z0;
    H[1][0] = +2*q1 * IMU_ACC_Z0;
    H[2][0] = +2*q0 * IMU_ACC_Z0;
    H[3][0] =  2*q0*IMU_MAG_B0[0][0] + 2*q3*IMU_MAG_B0[1][0] - 2*q2*IMU_MAG_B0[2][0];
    H[4][0] = -2*q3*IMU_MAG_B0[0][0] + 2*q0*IMU_MAG_B0[1][0] + 2*q1*IMU_MAG_B0[2][0];
    H[5][0] =  2*q2*IMU_MAG_B0[0][0] - 2*q1*IMU_MAG_B0[1][0] + 2*q0*IMU_MAG_B0[2][0];
    
    H[0][1] = +2*q3 * IMU_ACC_Z0;
    H[1][1] = +2*q0 * IMU_ACC_Z0;
    H[2][1] = -2*q1 * IMU_ACC_Z0;
    H[3][1] =  2*q1*IMU_MAG_B0[0][0]+2*q2*IMU_MAG_B0[1][0] + 2*q3*IMU_MAG_B0[2][0];
    H[4][1] =  2*q2*IMU_MAG_B0[0][0]-2*q1*IMU_MAG_B0[1][0] + 2*q0*IMU_MAG_B0[2][0];
    H[5][1] =  2*q3*IMU_MAG_B0[0][0]-2*q0*IMU_MAG_B0[1][0] - 2*q1*IMU_MAG_B0[2][0];
    
    H[0][2] = -2*q0 * IMU_ACC_Z0;
    H[1][2] = +2*q3 * IMU_ACC_Z0;
    H[2][2] = -2*q2 * IMU_ACC_Z0;
    H[3][2] = -2*q2*IMU_MAG_B0[0][0]+2*q1*IMU_MAG_B0[1][0] - 2*q0*IMU_MAG_B0[2][0];
    H[4][2] =  2*q1*IMU_MAG_B0[0][0]+2*q2*IMU_MAG_B0[1][0] + 2*q3*IMU_MAG_B0[2][0];
    H[5][2] =  2*q0*IMU_MAG_B0[0][0]+2*q3*IMU_MAG_B0[1][0] - 2*q2*IMU_MAG_B0[2][0];
    
    H[0][3] = +2*q1 * IMU_ACC_Z0;
    H[1][3] = +2*q2 * IMU_ACC_Z0;
    H[2][3] = +2*q3 * IMU_ACC_Z0;
    H[3][3] = -2*q3*IMU_MAG_B0[0][0]+2*q0*IMU_MAG_B0[1][0] + 2*q1*IMU_MAG_B0[2][0];
    H[4][3] = -2*q0*IMU_MAG_B0[0][0]-2*q3*IMU_MAG_B0[1][0] + 2*q2*IMU_MAG_B0[2][0];
    H[5][3] =  2*q1*IMU_MAG_B0[0][0]+2*q2*IMU_MAG_B0[1][0] + 2*q3*IMU_MAG_B0[2][0];
    
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
