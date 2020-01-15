/**************************************************************************************************
 * 
 * 
 * See https://github.com/pronenewbits for more!
 *************************************************************************************************/
#include <Wire.h>
#include <elapsedMillis.h>
#include "konfig.h"
#include "matrix.h"
#include "ekf.h"
#include "MPU9250.h"


// 0.01,-0.03,-1.00,    -0.41,-0.73,0.55,   -0.82,0.60,-0.12,
// 0.00,-0.02,-1.00,    -0.38,0.74,-0.55,   -0.05,-0.19,0.14,
// #define IMU_ACC_Z0          (-1)
// #define IMU_MAG_B0x         (-0.38)
// #define IMU_MAG_B0y         (-0.74)
// #define IMU_MAG_B0z         (-0.55)
#define IMU_ACC_Z0          (-1)
#define IMU_MAG_B0x         (1)
#define IMU_MAG_B0y         (1)
#define IMU_MAG_B0z         (0.0)


/* The hard-magnet bias */
float_prec BIAS_MAG[3] = {8.9254665,  8.040476, -25.126487};
// float_prec BIAS_MAG[3] = {12.11, 13.25, -25.83};

/* An MPU9250 object with the MPU-9250 sensor on I2C bus 0 with address 0x68 */
MPU9250 IMU(Wire, 0x68);

/* Just example */
#define P_INIT      (10.)
#define Q_INIT      (1e-6)
#define R_INIT      (0.0015)



bool Main_bUpdateNonlinearX(Matrix &X_Next, Matrix &X, Matrix &U);
bool Main_bUpdateNonlinearY(Matrix &Y, Matrix &X, Matrix &U);
bool Main_bCalcJacobianF(Matrix &F, Matrix &X, Matrix &U);
bool Main_bCalcJacobianH(Matrix &H, Matrix &X, Matrix &U);

Matrix quaternionData(SS_X_LEN, 1);
Matrix Y(SS_Z_LEN, 1);
Matrix U(SS_U_LEN, 1);
EKF EKF_IMU(quaternionData, P_INIT, Q_INIT, R_INIT, Main_bUpdateNonlinearX, Main_bUpdateNonlinearY, Main_bCalcJacobianF, Main_bCalcJacobianH);

elapsedMillis timerLed, timerEKF;
uint64_t u64compuTime;
char bufferTxSer[100];
/* The command from the PC */
char cmd;



void setup() {
    /* serial to display data */
    Serial.begin(115200);
    while(!Serial) {}
    
    /* start communication with IMU */
    int status = IMU.begin();
    if (status < 0) {
        Serial.println("IMU initialization unsuccessful");
        Serial.println("Check IMU wiring or try cycling power");
        Serial.print("Status: ");
        Serial.println(status);
        while(1) {}
    }
    IMU.setAccelRange(MPU9250::ACCEL_RANGE_2G);
    IMU.setGyroRange(MPU9250::GYRO_RANGE_2000DPS);
    IMU.setDlpfBandwidth(MPU9250::DLPF_BANDWIDTH_184HZ);
    IMU.setSrd(19);
    
    
    /* x(k=0) = [1 0 0 0]' */
    quaternionData.vIsiNol();
    quaternionData[0][0] = 1.0;
    EKF_IMU.vReset(quaternionData, P_INIT, Q_INIT, R_INIT);
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


void loop() {
    if (timerEKF > SS_DT_MILIS) {
        /* ================== Read the sensor data / simulate the system here ================== */
        IMU.readSensor();
        Y[0][0] = IMU.getAccelX_mss();
        Y[1][0] = IMU.getAccelY_mss();
        Y[2][0] = IMU.getAccelZ_mss();
        Y[3][0] = IMU.getMagX_uT()-BIAS_MAG[0];
        Y[4][0] = IMU.getMagY_uT()-BIAS_MAG[1];
        Y[5][0] = IMU.getMagZ_uT()-BIAS_MAG[2];
        U[0][0] = IMU.getGyroX_rads();
        U[1][0] = IMU.getGyroY_rads();
        U[2][0] = IMU.getGyroZ_rads();
        /* Normalisasi */
        float_prec _normG = (Y[0][0] * Y[0][0]) + (Y[1][0] * Y[1][0]) + (Y[2][0] * Y[2][0]);
        _normG = sqrt(_normG);
        Y[0][0] = Y[0][0] / _normG;
        Y[1][0] = Y[1][0] / _normG;
        Y[2][0] = Y[2][0] / _normG;
        float_prec _normM = (Y[3][0] * Y[3][0]) + (Y[4][0] * Y[4][0]) + (Y[5][0] * Y[5][0]);
        _normM = sqrt(_normM);
        Y[3][0] = Y[3][0] / _normM;
        Y[4][0] = Y[4][0] / _normM;
        Y[5][0] = Y[5][0] / _normM;
        /* ------------------ Read the sensor data / simulate the system here ------------------ */
        
        
        
        /* ============================= Update the Kalman Filter ============================== */
        u64compuTime = micros();
        if (!EKF_IMU.bUpdate(Y, U)) {
            quaternionData.vIsiNol();
            quaternionData[0][0] = 1.0;
            EKF_IMU.vReset(quaternionData, P_INIT, Q_INIT, R_INIT);
            Serial.println("Whoop ");
        }
        u64compuTime = (micros() - u64compuTime);
        /* ----------------------------- Update the Kalman Filter ------------------------------ */
        
        
        
        /* =========================== Print to serial (for plotting) ========================== */
        #if (0)     /* Enable this to calibrate BIAS_MAG (hard-magnet values) */
            static float_prec MAG_BIAS_X_MAX = 0, MAG_BIAS_X_MIN = 0;
            static float_prec MAG_BIAS_Y_MAX = 0, MAG_BIAS_Y_MIN = 0;
            static float_prec MAG_BIAS_Z_MAX = 0, MAG_BIAS_Z_MIN = 0;
            static float_prec MAGNETO_X = 0, MAGNETO_Y = 0, MAGNETO_Z = 0;
            
            MAGNETO_X = IMU.getMagX_uT();
            MAGNETO_Y = IMU.getMagY_uT();
            MAGNETO_Z = IMU.getMagZ_uT();
            
            (MAGNETO_X > MAG_BIAS_X_MAX)? MAG_BIAS_X_MAX = MAGNETO_X:1;
            (MAGNETO_Y > MAG_BIAS_Y_MAX)? MAG_BIAS_Y_MAX = MAGNETO_Y:1;
            (MAGNETO_Z > MAG_BIAS_Z_MAX)? MAG_BIAS_Z_MAX = MAGNETO_Z:1;
            (MAGNETO_X < MAG_BIAS_X_MIN)? MAG_BIAS_X_MIN = MAGNETO_X:1;
            (MAGNETO_Y < MAG_BIAS_Y_MIN)? MAG_BIAS_Y_MIN = MAGNETO_Y:1;
            (MAGNETO_Z < MAG_BIAS_Z_MIN)? MAG_BIAS_Z_MIN = MAGNETO_Z:1;
            
            /* Print raw magnetometer */
            Serial.print(MAGNETO_X);
            Serial.print(", ");
            Serial.print(MAGNETO_Y);
            Serial.print(", ");
            Serial.print(MAGNETO_Z);
            Serial.print(", ");
            
            /* Print calibrated magnetometer */
            Serial.print(MAGNETO_X - (MAG_BIAS_X_MAX+MAG_BIAS_X_MIN)/2.);
            Serial.print(", ");
            Serial.print(MAGNETO_Y - (MAG_BIAS_Y_MAX+MAG_BIAS_Y_MIN)/2.);
            Serial.print(", ");
            Serial.print(MAGNETO_Z - (MAG_BIAS_Z_MAX+MAG_BIAS_Z_MIN)/2.);
            Serial.print(", ");
            
            /* Print the calculated bias */
            Serial.print((MAG_BIAS_X_MAX+MAG_BIAS_X_MIN)/2.);
            Serial.print(", ");
            Serial.print((MAG_BIAS_Y_MAX+MAG_BIAS_Y_MIN)/2.);
            Serial.print(", ");
            Serial.print((MAG_BIAS_Z_MAX+MAG_BIAS_Z_MIN)/2.);
            Serial.print(", ");
            
            Serial.println("");
        #endif
        
        
        #if (0)     /* Enable this to calibrate IMU_MAG_B0x, IMU_MAG_B0y, IMU_MAG_B0z value (use this __after__ you calibrated the BIAS_MAG!) */
            for (uint8_t _i = 0; _i < SS_Z_LEN; _i++) {
                Serial.print(Y[_i][0]);
                Serial.print(",");
            }
            for (uint8_t _i = 0; _i < SS_U_LEN; _i++) {
                Serial.print(U[_i][0]);
                Serial.print(",");
            }
            Serial.println("");
        #endif
        /* --------------------------- Print to serial (for plotting) -------------------------- */
        
        
        timerEKF = 0;
    }
    
    
    /* The serial data is sent by responding to command from the PC running Processing scipt */
    if (Serial.available()) {
        cmd = Serial.read();
        if (cmd == 'v') {
            snprintf(bufferTxSer, sizeof(bufferTxSer)-1, "EKF in Teensy 4.0 (%s)", (FPU_PRECISION == PRECISION_SINGLE)?"Float32":"Double64");
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


bool Main_bUpdateNonlinearX(Matrix &X_Next, Matrix &X, Matrix &U)
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

bool Main_bUpdateNonlinearY(Matrix &Y, Matrix &X, Matrix &U)
{
    /* Insert the nonlinear measurement transformation here
     *          y(k)   = h[x(k), u(k)]
     *
     * The measurement output is the gravitational and magnetic projection to the body
     *     DCM     = [(+(q0**2)+(q1**2)-(q2**2)-(q3**2)),                        2*(q1*q2+q0*q3),                        2*(q1*q3-q0*q2)]
     *               [                   2*(q1*q2-q0*q3),     (+(q0**2)-(q1**2)+(q2**2)-(q3**2)),                        2*(q2*q3+q0*q1)]
     *               [                   2*(q1*q3+q0*q2),                        2*(q2*q3-q0*q1),     (+(q0**2)-(q1**2)-(q2**2)+(q3**2))]
     * 
     *  G_proj_sens = DCM * [0 0 -g]            --> Gravitational projection to the accelerometer sensor
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

    Y[3][0] = (+(q0_2)+(q1_2)-(q2_2)-(q3_2)) * IMU_MAG_B0x
             +(2*(q1*q2+q0*q3)) * IMU_MAG_B0y
             +(2*(q1*q3-q0*q2)) * IMU_MAG_B0z;

    Y[4][0] = (2*(q1*q2-q0*q3)) * IMU_MAG_B0x
             +(+(q0_2)-(q1_2)+(q2_2)-(q3_2)) * IMU_MAG_B0y
             +(2*(q2*q3+q0*q1)) * IMU_MAG_B0z;

    Y[5][0] = (2*(q1*q3+q0*q2)) * IMU_MAG_B0x
             +(2*(q2*q3-q0*q1)) * IMU_MAG_B0y
             +(+(q0_2)-(q1_2)-(q2_2)+(q3_2)) * IMU_MAG_B0z;
    
    return true;
}

bool Main_bCalcJacobianF(Matrix &F, Matrix &X, Matrix &U)
{
    /* Kode python box_quaternion_tanpa_magnetometer_EKF_vIMU6DOF+HMC.py:
     *  q0_dot = 1/2. * (  0   - p*q1 - q*q2 - r*q3)
     *  q1_dot = 1/2. * ( p*q0 +   0  + r*q2 - q*q3)
     *  q2_dot = 1/2. * ( q*q0 - r*q1 +  0   + p*q3)
     *  q3_dot = 1/2. * ( r*q0 + q*q1 - p*q2 +  0  )
     */
    /* Kode python box_quaternion_tanpa_magnetometer_EKF_vIMU6DOF+HMC.py:
        F     = numpy.array([[    0, -1/2.*p, -1/2.*q, -1/2.*r],
                             [1/2.*p,      0,  1/2.*r, -1/2.*q],
                             [1/2.*q, -1/2.*r,      0,  1/2.*p],
                             [1/2.*r,  1/2.*q, -1/2.*p,      0]],
                             "float64")
                             
                             
    X_Next[0][0] = (0.5 * (+0.00 -p*q1 -q*q2 -r*q3))*SS_DT + q0;
    X_Next[1][0] = (0.5 * (+p*q0 +0.00 +r*q2 -q*q3))*SS_DT + q1;
    X_Next[2][0] = (0.5 * (+q*q0 -r*q1 +0.00 +p*q3))*SS_DT + q2;
    X_Next[3][0] = (0.5 * (+r*q0 +q*q1 -p*q2 +0.00))*SS_DT + q3;
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

bool Main_bCalcJacobianH(Matrix &H, Matrix &X, Matrix &U)
{
    /* Kode python box_quaternion_EKF_vIMU6DOF+HMC.py:
     *     DCM     = numpy.array([[(+(q0**2)+(q1**2)-(q2**2)-(q3**2)),                        2*(q1*q2+q0*q3),                        2*(q1*q3-q0*q2)],
     *                           [                   2*(q1*q2-q0*q3),     (+(q0**2)-(q1**2)+(q2**2)-(q3**2)),                        2*(q2*q3+q0*q1)],
     *                           [                   2*(q1*q3+q0*q2),                        2*(q2*q3-q0*q1),     (+(q0**2)-(q1**2)-(q2**2)+(q3**2))]],
     *                           "float32")
     * 
     *  G_proj_sens = DCM * [0 0 -g]              --> Proyeksi gravitasi ke sensors
     *  M_proj_sens = DCM * [Mx My 0]             --> Proyeksi medan magnet ke sensors
     */

    float_prec q0, q1, q2, q3;

    q0 = X[0][0];
    q1 = X[1][0];
    q2 = X[2][0];
    q3 = X[3][0];

    H[0][0] =  2*q2 * IMU_ACC_Z0;
    H[1][0] = -2*q1 * IMU_ACC_Z0;
    H[2][0] = -2*q0 * IMU_ACC_Z0;
    H[3][0] =  2*q0*IMU_MAG_B0x + 2*q3*IMU_MAG_B0y - 2*q2*IMU_MAG_B0z;
    H[4][0] = -2*q3*IMU_MAG_B0x + 2*q0*IMU_MAG_B0y + 2*q1*IMU_MAG_B0z;
    H[5][0] =  2*q2*IMU_MAG_B0x - 2*q1*IMU_MAG_B0y + 2*q0*IMU_MAG_B0z;

    H[0][1] = -2*q3 * IMU_ACC_Z0;
    H[1][1] = -2*q0 * IMU_ACC_Z0;
    H[2][1] =  2*q1 * IMU_ACC_Z0;
    H[3][1] =  2*q1*IMU_MAG_B0x+2*q2*IMU_MAG_B0y + 2*q3*IMU_MAG_B0z;
    H[4][1] =  2*q2*IMU_MAG_B0x-2*q1*IMU_MAG_B0y + 2*q0*IMU_MAG_B0z;
    H[5][1] =  2*q3*IMU_MAG_B0x-2*q0*IMU_MAG_B0y - 2*q1*IMU_MAG_B0z;

    H[0][2] =  2*q0 * IMU_ACC_Z0;
    H[1][2] = -2*q3 * IMU_ACC_Z0;
    H[2][2] =  2*q2 * IMU_ACC_Z0;
    H[3][2] = -2*q2*IMU_MAG_B0x+2*q1*IMU_MAG_B0y - 2*q0*IMU_MAG_B0z;
    H[4][2] =  2*q1*IMU_MAG_B0x+2*q2*IMU_MAG_B0y + 2*q3*IMU_MAG_B0z;
    H[5][2] =  2*q0*IMU_MAG_B0x+2*q3*IMU_MAG_B0y - 2*q2*IMU_MAG_B0z;

    H[0][3] = -2*q1 * IMU_ACC_Z0;
    H[1][3] = -2*q2 * IMU_ACC_Z0;
    H[2][3] = -2*q3 * IMU_ACC_Z0;
    H[3][3] = -2*q3*IMU_MAG_B0x+2*q0*IMU_MAG_B0y + 2*q1*IMU_MAG_B0z;
    H[4][3] = -2*q0*IMU_MAG_B0x-2*q3*IMU_MAG_B0y + 2*q2*IMU_MAG_B0z;
    H[5][3] =  2*q1*IMU_MAG_B0x+2*q2*IMU_MAG_B0y + 2*q3*IMU_MAG_B0z;
    
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
