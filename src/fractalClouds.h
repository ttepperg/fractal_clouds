// --------------------------------------------------------------------------------------
// YOUR COMPILE OPTIONS:
//

//! Set to 1 for debbuging
#define DEBUG 0

//! choose ONE data type (set to 1): float (4 bytes) or double (8 bytes)
#define FLOAT 0
#define DOUBLE 1

//! Number of iterations (min: 1; no more than ~20 should be necessary)
#define MAX_ITER 15

//! Set to 1 (0) to (not) write out binary file
#define BINARYOUT 1

//! Set to 1 to use a random (otherwise fixed to -1) random-generator seed
#define RANDOM 0

//! Set to 1 to include Nyquist frequency in power spectrum
#define NYQUIST 0

/*!
	\f$\beta\f$ < 0 is the power-law index of the target power spectrum  \f$P(k) \propto k^{\beta}\f$. Note that the spectral density \f$|F|^2\f$ will be a power-law \f$k^s\f$ with \f$s = \beta - (n - 1)\f$, where \f$n\f$ is the number of dimensions. 
*/
#define BETA (-1.5)
// #define BETA (-5./3.) // Kolmogorov turbulence spectrum

//! Minimum wavenumber (i.e., largest scale) to be present (not yet used)
//! Minimum value: 1  Maximum: DIM-1
#define KMIN 1

// Input Gaussian parameters
#define SIG2_GAUSS 4.   // if = 0, use target log-normal parameters directly;
								// if != 0, make sure to set MU_GAUSS to the desired value!
#define MU_GAUSS 0.

// Target log-normal parameters (will be set automatically if SIG2_GAUSS != 0)
#define SIG2_LOGN 5.
#define MU_LOGN 1.

// Choose one:
// #define SIM1D 1
#define SIM2D 2
//  #define SIM3D 3

/*!	IMPORTANT: the cube linear size DIM must be a power of 2!
		NOTE: strictly, the size of the useful simulation will be half the size
		define by DIM (due to aliasing effects when FFT-ing), although the full
		cube is still useful.
		MAXIMUM RECOMMENDED sizes for each data type are:
			DOUBLE: 8192 (1D); 1024 (2D); 256 (3D)
			FLOAT: >8192 (1D); 2048 (2D); 512 (3D)
*/
#define DIM (512)

// Type of fit to power spectrum

#define SGFILTER 0	// Savitzky-Golay filter of order:
#define SGORDER 1
#define HALFWINDOW 2 // window size is 2*HALFWINDOW + 1

//OR (default: SGFILTER of order 1)

#define POLYFIT 1	// polynomial fit of degree:
#define POLYDEG 3



// --------------------------------------------------------------------------------------
