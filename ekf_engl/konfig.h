#ifndef KONFIG_H
#define KONFIG_H

#include <stdlib.h>
#include <stdint.h>
#include <math.h>



/* State Space dimension */
#define SS_X_LEN    (4)
#define SS_Z_LEN    (6)
#define SS_U_LEN    (3)
#define SS_DT_MILIS (10)                            /* 10 ms */
#define SS_DT       float_prec(SS_DT_MILIS/1000.)   /* Sampling time */



/* Change this size based on the biggest matrix you will use */
#define MATRIX_MAXIMUM_SIZE     (10)

/* Define this to enable matrix bound checking */
#define MATRIX_PAKAI_BOUND_CHECKING

/* Set this define to choose math precision of the system */
#define PRECISION_SINGLE    1
#define PRECISION_DOUBLE    2
#define FPU_PRECISION       (PRECISION_SINGLE)

#if (FPU_PRECISION == PRECISION_SINGLE)
    #define float_prec      float
    #define float_prec_ZERO (1e-8)
#elif (FPU_PRECISION == PRECISION_DOUBLE)
    #define float_prec      double
    #define float_prec_ZERO (1e-15)
#else
    #error("FPU_PRECISION has not been defined!");
#endif



/* Set this define to choose system implementation (mainly used to define how you print the matrix via the Matrix::vCetak() function) */
#define SYSTEM_IMPLEMENTATION_PC                    1
#define SYSTEM_IMPLEMENTATION_EMBEDDED_NO_PRINT     2
#define SYSTEM_IMPLEMENTATION_EMBEDDED_ARDUINO      3

#define SYSTEM_IMPLEMENTATION                       (SYSTEM_IMPLEMENTATION_EMBEDDED_ARDUINO)



#endif // KONFIG_H
