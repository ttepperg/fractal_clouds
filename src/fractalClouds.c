/*! 
  \file  fractalClouds.c
  \brief Fractal cloud generator

	This code generates a \e single-point log-normal random density field <b>\e H</b> - via exponentiation of a \e seed Gaussian random field <b>\e G</b> - that has a power spectrum \f$P(k) \propto k^{\beta}\f$, with \f$\beta < 0\f$, with a specified index \f$1 < |\beta| < 3\f$. \f$F\f$ is the Fourier transform of \f$H\f$. It is an implementation of the algorithm by [Lewis & Austin (2002)](http://www.eos.ubc.ca/research/clouds/radjp416), used by, e.g., [Sutherland & Bicknell (2003)](http://stacks.iop.org/0004-637X/591/i=1/a=238), [Sutherland & Bicknell (2007)](http://adsabs.harvard.edu/abs/2007ApJS..173...37S) to generate fractal gas clouds such as the one used by [Bland-Hawthorn et al. (2007)](http://adsabs.harvard.edu/abs/2007ApJ...670L.109B) and [Tepper-Garcia et al. (2015a)](http://adsabs.harvard.edu/abs/2015ApJ...813...94T). For now only available in serial. Note that 1D simulation power spectra are
	very noisy, since their value at each wavenumber includes at most two contributions. In contrast, 8*(k-1) and 2*(12*(k-1)*(k-1) + 1) values enter into the calculation of the power spectrum at k = 1, 2, ... in 2D and 3D simulations, respectively, and appear therefore smoother than 1D power spectra, even more with increasing k.

	\attention 	We make use of a number of routines provided by the Numerical Recipes (NR). These have all (and their appropriate arguments) been modified to be of type Real (float or double depending on the user's choice) precision. For this reason, all the NR routine names have been pre-appended with a 'd', i.e. \c gasdev has been renamed to \c dgasdev. Also, NR routines adopt the "unit-offset" convention, meaning that the lowest subscript of an array is 1 (rather than 0 as in ANSI C). We use the functions provided by \c dnrutil to allocate matrices and tensors for consistency with other NR routines.\n Both the Gaussian ad log-normal random field are considered \e complex (although the imaginary part is set to 0 in the initial Gaussian random field), since their respective Fourier transform will be in general complex. For this reason, we use of the routines \c dfour1 and \c dfourn to compute FFT in 1 and 2,3 dimensions, respectively. Note that C does not provide a complex data type, and thus a complex number consists of a pair of real numbers, each member representing the real and imaginary part, respectively.


	\b Compilation (or replace gcc-mp-4.9 by your gcc compiler)

	$> gcc-mp-4.9 dnrutil.c dlubksb.c dludcmp.c dsavgol.c dfour1.c dfourn.c dran1.c dgasdev.c dgammln.c dgcf.c dgser.c dgammq.c dfit.c dlfit.c dcovsrt.c dgaussj.c alloc_funcs.c fractalClouds.c -o fractalClouds

	<b> Program call </b>

	$> ./fractalClouds

	\b Output
	
	All files are simple binary files of type Real (4 or 8 byte depending on the user's choice) per value. For a simulation with dimensionality N=1,2,3 and size DIM, these files are:
		- \c gaussNd_ini.bdat: File with the initial Gaussian random field, and it contains \f$2({\rm DIM})^{\rm N}\f$ values. Each pair of consecutive values correspond to the real and imaginary part of a complex number representing the value of the field at \f$x\f$, \f$(x,y)\f$, or \f$(x,y,z)\f$ for N=1, 2, or 3, respectively. The values are stored contiguously along the \f$x\f$-coordinate, followed by \f$y\f$, and \f$z\f$.
		- <c> powSpecNd_<field>#.bdat </c>: Power spectrum \f$P(k)\f$ of the <field> = gauss | lognorm, at iteration number #. The file contains DIM numbers, each representing the value of \f$P(k)\f$ at \f$k_n\f$, where \f$k_n = n\f$/DIM, \f$n = 0, 1, ...,\f$DIM-1. 
		- <c> searchDirNd_#.bdat </c>: Search direction, i.e., difference between target and current (at iteration number #) power spectrum of log-normal random field. The data structure is the same as <c> powSpecNd_<field>#.bdat</c>.
		- <c> <field>Nd_filt </c>: Final random field, where <field> = gauss | lognorm. The data structure is the same as \c gaussNd_ini.bdat.

	<b> Accompanying software </b>

		- \c plotField.gp: An editable gnuplot script to visualise the result. Edits include setting the dimensionality of the problem, and the field to visualise. Requires gnuplot 4.6 or higher. 

		Script call:

		$> gnuplot plotField.gp

	\authors T.Tepper Garcia (tepper@physics.usyd.edu.au)

	\date   Nov 9, 2015

	\todo (!) = URGENT
		- Encode convergence condition to terminate iteration
		- Write initial guess for spectral density to disk (binary)
		- Encapsulate computation of target spectral density into a callable function
		- Encapsulate computation of power spectrum into a callable function
		- Encapsulate computation of new spectral density into a callable function
		- Encapsulate computation of search direction into a callable function
		- Encapsulate writing of binary files into a callable function
		- Expand plotField.gp documentation
		- Parallelise code

	<b> Modification history </b>

		- 21 NOV 2016 (TTG): Added binary dump of REAL part of filtered Gaussian / lognormal field
		- 2 DEC 2016 (TTG): Introduced a new parameter (KMIN; minimum wavenumbe) to control the maximum cloudlet size,
		given by DIM / KMIN

	\b References
		- Press, W.H. and Teukolsky, S.A. and Vetterling, W.T. and Flannery, B.P., "Numerical Recipes in C", Cambridge University Press, 1992, Second Edition
		- Lewis, G. M. and Austin, P. H., "An iterative method for generating scaling log-normal simulations", Proceedings of the 11th Conference on Atmospheric Radiation, 2002, American Meteorological Society Conference Series 123-126
		- Sutherland, Ralph S. and Bicknell, Geoffrey V. and Dopita, Michael A., "The Numerical Simulation of Radiative Shocks. II. Thermal Instabilities in Two-dimensional Models", 2003, The Astrophysical Journal, Volume 591, Number 1, Pages, 238-257
		- Sutherland, R.S. and Bicknell, G.V., "Interactions of a Light Hypersonic Jet with a Nonuniform Interstellar Medium", The Astrophysical Journal Supplement, 2007, Volume 173, Pages 37-69

*/
// --------------------------------------------------------------------------------------

#include "fractalCloudsDefs.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "alloc_funcs.h"
#include "dnrutil.h"

#define PI 3.14159265358979
#define SQRTPI sqrt(PI)
#define SQRT2 sqrt(2.)
#define SQRT3 sqrt(3.)

// 2D is default
#ifdef SIM2D
	#undef SIM1D
	#undef SIM3D
// 	#define DIM (1024)
	#define DIM1 1
	#define DIM2 DIM
	#define DIM3 DIM
#endif

#ifdef SIM1D
	#undef SIM2D
	#undef SIM3D
// 	#define DIM (8192)
	#define DIM3 DIM
#endif

#ifdef SIM3D
	#undef SIM1D
	#undef SIM2D
// 	#define DIM (256)
	#define DIM1 DIM
	#define DIM2 DIM
	#define DIM3 DIM
#endif

// ---------------------------------------------------------------------------------------
// FUNCTION PROTOTYPES

Real dgasdev(long *idum);
Real dran1(long *idum);
void dfourn(Real data[], unsigned long nn[], int ndim, int isign);
void dlfit(Real x[], Real y[], Real sig[], int ndat, Real a[], int ia[],
	int ma, Real **covar, Real *chisq, void (*funcs)(Real, Real [], int));
void dfit(Real x[], Real y[], int ndata, Real sig[], int mwt, Real *a,
	Real *b, Real *siga, Real *sigb, Real *chi2, Real *q);
void dsavgol(Real c[], int np, int nl, int nr, int ld, int m);

short wrap(short i, short n, short m);
Real sed(Real k, Real s);
void poly(Real x, Real pn[], int n); // polynomial of degree n
Real min(Real x, Real y);
Real max(Real x, Real y);

Real cmp_mu(Real m, Real s);
Real cmp_sigma(Real m, Real s);
Real cmp_m(Real mu, Real sig);
Real cmp_s(Real mu, Real sig);

#ifdef SIM1D
	Real average(Real *arr, int l);
	Real stddev(Real *arr, int l);
	Real average_log(Real *arr, int l);
	Real stddev_log(Real *arr, int l);
#else
	Real average(Real ***arr, int n, int m, int l);
	Real stddev(Real ***arr, int n, int m, int l);
	Real average_log(Real ***arr, int n, int m, int l);
	Real stddev_log(Real ***arr, int n, int m, int l);
#endif

// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{ /*!
	Main Program
*/

// ---------------------------------------------------------------------------------------
// VARIABLE DECLARATION

	int iterNum = 0; // iteration counter

	long seed = 0; // type long required by dgasdev.c

	Real dimCheck;
	unsigned long* dims = NULL; // type long required by dfourn.c
	short numDims = 1; // will be set to correct value later

	short m; // common dummy indices; used to compute annular average of pow.spec.
	int sumTerms; // number count of points at each annulus

	// Gaussian parameters
	Real muGauss = MU_GAUSS;	// if sigmaGauss2 != 0, use Gaussian parameters to set log-normal
										// parameters
	Real sigmaGauss2 = SIG2_GAUSS;	// if = 0, use target log-normal parameters directly;
										// if != 0, make sure to set muGauss to the desired value!
	Real sigmaGauss;

	Real muAux, sigmaAux;

	// Target log-normal parameters
	Real muLogn = MU_LOGN;
	Real sigmaLogn2 = SIG2_LOGN;
	Real muLogn2;		// don't initialise!
	Real sigmaLogn;	// don't initialise!

	double normFactSpecDensIni; // normalisation of guess spectral density
	double normFactSpecDensFin; // normalisation of final spectral density
	
	Real normFact; // dummy variable used for normalisation in different situations
	
	Real* powSpecGauss = NULL;

	Real* powSpecLogNorm = NULL;
	Real* logPowSpecLogNorm = NULL;	// log10 of power spectrum P(k)
	Real* logWaveNumber = NULL; // log10(k); to fit with a linear model
	Real* sigmaWeight = NULL; // error assigned to P(k) to fit linear model;
										// the larger this value, the lesser the weight

	// linear fitting parameters
	Real normConstLogNorm, slopeLogNorm;
	Real inter = 0., slopeErr = 0., interErr = 0., chi2 = 0., qParam = 0.;

	// polynomial fitting parameters
#if (POLYFIT)
	Real* polyFitLogPowSpecLogNorm = NULL; // polynomal fit to log-normal power spectrum
	int polyDeg = POLYDEG;	// degree of polynomial;; Lewis a& Austin (2002) recommend 5
						  			// both 3 and 5 work well
	int numParams = polyDeg+1; // number of parameters to fit
	Real* polyVal = NULL;	// stores each of the values x^0, x^1, ... , x^polyDeg
	Real* polyCoeff = NULL; 	// list of parameters for polynomial fitting to log-normal P(k)
	int* polyCoeffFlag = NULL;	// list of parameters to be held constant (0) during fit;
										// default for each parameter is 1 (i.e. optimize)
	Real **covar;			// covariance matrix
	Real chisq;				// chi^2 value of fit
#endif // POLYFIT

	// Savitzky-Golary filter (smoothing) parameters
	// alternative to polynomial fit to log-normal P(k)
#if (SGFILTER)
	int sgOrder = SGORDER;
	int derOrder = 0.; // don't change!
	int halfWindow = HALFWINDOW;
	int windowSize = 2*halfWindow + 1;
	Real *sgCoeff;
	Real* sgLogPowSpecLogNorm = NULL; // SG fit to log-normal power spectrum
	Real signal;
#endif // SGFILTER

	Real *searchDir;// difference between target and actual log-normal spectrum
	Real searchDirMean;// maximum relative difference
	Real stepSize, searchStep;// adaptive step to take difference

	FILE* filePtr = NULL;	// generic file pointer
	char fileName[100];		// generic output file name
	size_t fileErr;			// generic file read/write error handle
	char fileExt[6] = ".bdat";	// binary data file extension
	char fileInfix[3];		// file name infix (differs with the dimension of the sim)

#ifdef SIM1D
	Real* noise = NULL;
	Real* noiseFiltered = NULL;
	Real* specDens = NULL;
	Real* logNormal = NULL;
	Real x;
	short i,j,k;
	short i3;
	short ii3;
	Real pointOutX;
#else
	Real*** noise = NULL;
	Real*** noiseFiltered = NULL;
	Real*** specDens = NULL;
	Real*** logNormal = NULL;
	Real x, y, z, w;
	short h, i, j, k;
	short i1, i2, i3;
	short ii1, ii2, ii3;
	Real pointOutX, pointOutY, pointOutZ;
#endif // SIM2D / SIM3D

// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// CHECK CORRECTNESS OF DIMENSIONS (should be an integral power of 2)
{
	dimCheck= log2(DIM3);
	if ( (dimCheck - (floor(dimCheck))) > FLT_EPSILON ) {
		printf("Dimension 1 not a power of 2! log2(DIM3) = %12.8e\n", dimCheck);
		return(0);
	}

#ifdef SIM2D
	dimCheck= log2(DIM1);
	if ( (dimCheck - (floor(dimCheck))) > FLT_EPSILON ) {
		printf("Dimension 1 not a power of 2! log2(DIM1) = %12.8e\n", dimCheck);
		return(0);
	}
	dimCheck= log2(DIM2);
	if ( (dimCheck - (floor(dimCheck))) > FLT_EPSILON ) {
		printf("Dimension 2 not a power of 2! log2(DIM2) = %12.8e\n", dimCheck);
		return(0);
	}
	dimCheck= log2(DIM3);
	if ( (dimCheck - (floor(dimCheck))) > FLT_EPSILON ) {
		printf("Dimension 3 not a power of 2! log2(DIM3) = %12.8e\n", dimCheck);
		return(0);
	}
#endif // SIM2D

#ifdef SIM3D
	dimCheck= log2(DIM1);
	if ( (dimCheck - (floor(dimCheck))) > FLT_EPSILON ) {
		printf("Dimension 1 not a power of 2! log2(DIM1) = %12.8e\n", dimCheck);
		return(0);
	}
	dimCheck= log2(DIM2);
	if ( (dimCheck - (floor(dimCheck))) > FLT_EPSILON ) {
		printf("Dimension 2 not a power of 2! log2(DIM2) = %12.8e\n", dimCheck);
		return(0);
	}
	dimCheck= log2(DIM3);
	if ( (dimCheck - (floor(dimCheck))) > FLT_EPSILON ) {
		printf("Dimension 3 not a power of 2! log2(DIM3) = %12.8e\n", dimCheck);
		return(0);
	}
#endif // SIM3D
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// SANITY CHECKS / DEFAULTS / ETC.

#if (POLYFIT)
#if (SGFILTER) // Savitzky-Golay filter is default
#undef POLYFIT
#define POLYFIT 0
#endif
#else
#if !(SGFILTER)
	printf("\nNo type of fit (POLYFIT | SGFILTER) specified.\n\n");
	return(0);
#endif
#endif

// ---------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------
// ALLOCATE INITAL ARRAYS
{
#ifdef SIM1D
	printf("\n\n# This simulation is in 1D (%d)\n\n", DIM3);
	numDims = SIM1D;
	noise = vector(1,2*DIM3); 				// Gaussian noise vector of size 2N to accomodate N
													// complex numbers
	noiseFiltered = vector(1,2*DIM3);	// filtered FFT of noise
	specDens = vector(1,DIM3);				// power spectrum; of size N to accommodate N reals
	logNormal = vector(1,2*DIM3);			// log-normal random field
	sprintf(fileInfix, "%s", "1d");
#endif // SIM1D

#ifdef SIM2D
	printf("\n\n# This simulation is in 2D (%dx%d)\n\n", DIM2, DIM3);
	numDims = SIM2D;
	sprintf(fileInfix, "%s", "2d");
#endif // SIM2D

#ifdef SIM3D
	printf("\n\n# This simulation is in 3D (%dx%dx%d)\n\n", DIM1, DIM2, DIM3);
	numDims = SIM3D;
	sprintf(fileInfix, "%s", "3d");
#endif // SIM3D

#ifdef SIM1D
	dims = lvector(1,1);
	dims[1] = DIM3;
#else
	dims = lvector(1,3);
	dims[1] = DIM1;
	dims[2] = DIM2;
	dims[3] = DIM3;
	noise = f3tensor(1,DIM1, 1,DIM2, 1,2*DIM3);
	noiseFiltered = f3tensor(1,DIM1, 1,DIM2, 1,2*DIM3);
	specDens = f3tensor(1,DIM1, 1,DIM2, 1,DIM3);
	logNormal = f3tensor(1,DIM1, 1,DIM2, 1,2*DIM3);
#endif // SIM2D / SIM3D
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// Generate Savitzky-Golay filter coefficients (used for smoothing; alternative to fitting
// a polynomial to log-normal P(k)
#if (SGFILTER)
{
	sgCoeff = vector(1,windowSize);
	dsavgol(&sgCoeff[1]-1,windowSize,halfWindow,halfWindow,derOrder,sgOrder);
	// cyclically shift array to transform from wrap-around order to 'right' order
	for(i=1;i<=halfWindow;i++){
		sgCoeff[0] = sgCoeff[windowSize];
		for(j=windowSize;j>0;j--){
			sgCoeff[j] = sgCoeff[j-1];
		}
	}
}
#endif
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// GENERATE 1-POINT UNCORRELATED REAL GAUSSIAN RANDOM FIELD (white noise)
// Note that muGauss = 0 is always assumed, and the field is then shifted by muGauss to generate the corresponding log-normal field
{
	// Compute Gaussian variance from target log-normal parameters
	// See Sutherland & Bicknell (2007), their equations B5, B6
	// Note that muGauss = 0 is used; and muGauss is used later to generate a log-normal field

	// sanity check
	if (muLogn <= 0.) {
		printf("ERROR: mu' must be larger than 0!\n");
		return(0);
	}
	if (sigmaLogn2 <= 0.) {
		printf("ERROR: sigma^2 must be larger than 0!\n");
		return(0);
	}

	if (sigmaGauss2 == 0.) { // compute Gaussian parameters from target log-normal parameters

		muLogn2 = muLogn * muLogn;
		sigmaLogn = sqrt(sigmaLogn2);

		muGauss = cmp_m(muLogn,sigmaLogn);
		sigmaGauss = cmp_s(muLogn,sigmaLogn);
		sigmaGauss2 = sigmaGauss*sigmaGauss;

	} else { // compute target log-normal parameters from input Gaussian parameters

		sigmaGauss = sqrt(sigmaGauss2);

		muLogn = cmp_mu(muGauss, sigmaGauss);
		sigmaLogn = cmp_sigma(muGauss, sigmaGauss);
		sigmaLogn2 = sigmaLogn*sigmaLogn;


	}

#if RANDOM
   srand(time(NULL));
	seed = -1 * (1 + rand() % 100); // seed between -101 and -1
#else
   srand(seed);
	seed = -1;
#endif
printf("# Random seed is %ld\n", seed );

#ifdef SIM1D

	for(k = 1; k <= DIM3; k++) {

		noise[2*k-1] = dgasdev(&seed);// + muGauss; // real part
		noise[2*k] = dgasdev(&seed);// + muGauss; // imaginary part

	}


	// compute actual mean / variance
	muAux = average(noise,DIM3);
	sigmaAux = stddev(noise,DIM3);

	printf("# current mu (GAUSS) is %lf\n", muAux);
	printf("# current sigma (GAUSS) is %lf\n", sigmaAux);

	sigmaAux = 1./sigmaAux;

	for(k = 1; k <= DIM3; k++) {

		noise[2*k-1] = sigmaAux*(noise[2*k-1] - muAux); // real part
		noise[2*k] = sigmaAux*(noise[2*k] - muAux); // imaginary part

	}

	// re-compute actual mean / variance
	muAux = average(noise,DIM3);
	sigmaAux = stddev(noise,DIM3);

#else // 2D / 3D 

	for(i = 1; i <= DIM1; i++)
		for(j = 1; j <= DIM2; j++)
			for(k = 1; k <= DIM3; k++) {

				noise[i][j][2*k-1] = dgasdev(&seed);// + muGauss; // real part
				noise[i][j][2*k] = dgasdev(&seed);// + muGauss; // imaginary part

			}

	// compute actual mean / variance
	muAux = average(noise,DIM1,DIM2,DIM3);
	sigmaAux = stddev(noise,DIM1,DIM2,DIM3);

	printf("# current mu (GAUSS) is %lf\n", muAux);
	printf("# current sigma (GAUSS) is %lf\n", sigmaAux);

	sigmaAux = 1./sigmaAux;

	// renormalise data
	for(i = 1; i <= DIM1; i++)
		for(j = 1; j <= DIM2; j++)
			for(k = 1; k <= DIM3; k++) {

				noise[i][j][2*k-1] = sigmaAux*(noise[i][j][2*k-1] - muAux);// + muGauss; // real part
				noise[i][j][2*k] = sigmaAux*(noise[i][j][2*k] - muAux);// + muGauss; // imaginary part

			}

	// re-compute actual mean / variance
	muAux = average(noise,DIM1,DIM2,DIM3);
	sigmaAux = stddev(noise,DIM1,DIM2,DIM3);

#endif // SIM2D / SIM3D

printf("# adjusted mu (GAUSS) is %lf\n", muAux);
printf("# adjusted sigma (GAUSS) is %lf\n", sigmaAux);



}
// ---------------------------------------------------------------------------------------

// INFO:
	printf("# Data type is %s\n", (FLOAT ? "FLOAT" : "DOUBLE") );
	printf("# KMIN is %d\n", KMIN);
	printf("# target (input) mu (GAUSS) is %lf\n", muGauss);
	printf("# target (input) sigma (GAUSS) is %lf\n", sigmaGauss);
	printf("# target (input) mu (LOGN) is %lf\n", muLogn);
	printf("# target (input) sigma (LOGN) is %lf\n\n", sigmaLogn);


// ---------------------------------------------------------------------------------------
// BINARY OUTPUT: INITIAL GAUSSIAN RANDOM FIELD (WHITE NOISE)
#if BINARYOUT
{
	// set output file name
	sprintf(fileName, "%s%s_ini%s", "gauss", fileInfix, fileExt);

#ifdef SIM1D

	if ( (filePtr = fopen(fileName, "w")) == NULL) {
		printf("ERROR: could not open file '%s'.\n", fileName);
		return(0);
	}

	fileErr = fwrite(&noise[1], sizeof(Real), 2*DIM3, filePtr);

	if ( fileErr != (2*DIM3)) {
		printf("ERROR: could not complete write to '%s'.\n", fileName);
		return(0);
	}

#else // 2D / 3D

	if ( (filePtr = fopen(fileName, "w")) == NULL) {
		printf("ERROR: could not open file '%s'.\n", fileName);
		return(0);
	}

	fileErr = 0;

	for (i = 1; i <= DIM1; i++)
		for (j = 1; j <= DIM2; j++)
			for (k = 1; k <= DIM3; k++) {
				fileErr += fwrite(&noise[i][j][2*k-1], sizeof(Real), 1, filePtr);
				fileErr += fwrite(&noise[i][j][2*k], sizeof(Real), 1, filePtr);
			}

	if ( fileErr != (DIM1*DIM2*2*DIM3)) {
		printf("ERROR: could not complete write to '%s'.\n", fileName);
		return(0);
	}

#endif

	// DON'T FORGET!!!
	fclose(filePtr);

}
#endif // BINARYOUT
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// COMPUTE *forward* Fast Fourier Transform of Gaussian random field
// IMPORTANT: normalisation factor equal to the product of the sampling rate Delta in each
// dimension, and which is necessary to compare the result from the Num.Rec. FFT to the
// corresponding result in the continuous case, is ignored here because is implicitly set
// to 1; see Num.Recipes, their eq. 12.1.18
{
#ifdef SIM1D

	dfour1(noise,DIM3,1);

#else // 2D/3D

	dfourn(&noise[1][1][1] - 1,dims,3,1); // test with complex input data array

#endif
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// Compute target spectral density D(k) ~ k^s of Gaussian random field, with s=beta-(n-1)
// n = 1, 2, 3 for 1D, 2D or 3D
// IMPORTANT: spectrum is folded at DIM/2 to avoid aliasing effects
{
#ifdef SIM1D
{

	ii3 = (short) (0.5 * (Real) DIM3);

	// initialise
	for (k = 1; k <= DIM3; k++){
		specDens[k] = 0.;
	}

	for(i3 = 2; i3 <= ii3+1; i3++) { // up to Nyquist frequency; skip wavenumber k = 0

		k = i3;
		x = (Real) (k-1);

		// take into account cut-off at minimum wavelength KMIN
		specDens[k] = ( (x >= KMIN) ? sed(x, (0.5 * (BETA - numDims + 1))) : 0. );

		if ((k > 1) && (k <= ii3)) { // fold wrt k
			specDens[DIM3+2-k] = specDens[k];
		} // k > 1

	} // i3-loop

	// normalise target spectral density
	normFactSpecDensIni = 0.;

	// compute normalisation constant to preserve variance
	// taking into account cut-off at minimum wavelength KMIN
	for (k = 1; k <= DIM3; k++)
		normFactSpecDensIni += (double)((double)specDens[k]*(double)specDens[k]);

  	normFactSpecDensIni = sqrt((double) (DIM3) / normFactSpecDensIni);

	// normalisation to preserve variance
	for (k = 1; k <= DIM3; k++)
		specDens[k] *= (Real) normFactSpecDensIni;

	// reset value at k = 0
	specDens[1] = 1.0; // condition to preserve mean value of field

}
#else // 2D / 3D

#ifdef SIM2D
{

	ii2 = (short) (0.5 * (Real) DIM2);
	ii3 = (short) (0.5 * (Real) DIM3);

	// initialise
	for (j = 1; j <= DIM2; j++)
		for (k = 1; k <= DIM3; k++){
			specDens[1][j][k] = 0.;
		}

		for (i2 = 1; i2 <= ii2+1; i2++) { // up to Nyquist frequency

			j = i2;
			y = (Real) (j-1);

			for(i3 = 1; i3 <= ii3+1; i3++) { // up to Nyquist frequency

				k = i3;
				z = (Real) (k-1);

				if ( (i2 > 1) || (i3 > 1) ) { // avoid w = 0

					w = (Real) sqrt(y*y + z*z);

					// take into account cut-off at minimum wavelength KMIN
					specDens[1][i2][i3] = ( (w >= KMIN) ? sed(w, (0.5 * (BETA - numDims + 1))) : 0. );

					if ((i2 > 1) && (i2 <= ii2)) { // fold 1 wrt j at frequency below Nyquist

						specDens[1][DIM2+2-i2][i3] = specDens[1][i2][i3];

						if ((i3 > 1) && (i3 <= ii3)) // fold 1 wrt j and k
							specDens[1][DIM2+2-i2][DIM3+2-i3] = specDens[1][i2][i3];

					} // 1 < j < ii2

					if ((i3 > 1) && (i3 <= ii3)){ // fold 1 wrt k
							specDens[1][i2][DIM3+2-i3] = specDens[1][i2][i3];
					}

				} // avoid w = 0

			} // i3-loop

		} // i2-loop

	// normalise target spectral density
	normFactSpecDensIni = 0.;

	// compute normalisation constant to preserve variance
	// taking into account cut-off at minimum wavelength KMIN
 	for (j = 1; j <= DIM2; j++)
		for (k = 1; k <= DIM3; k++)
			normFactSpecDensIni += (double)((double)specDens[1][j][k]*(double)specDens[1][j][k]);

  	normFactSpecDensIni = sqrt((double) (DIM2*DIM3) / normFactSpecDensIni);

	// normalisation to preserve variance
	for (j = 1; j <= DIM2; j++)
		for (k = 1; k <= DIM3; k++)
			specDens[1][j][k] *= (Real) normFactSpecDensIni;

	// reset value at k = 0
	specDens[1][1][1] = 1.; // condition to preserve mean value of field


}
#else // 3D
{

	ii1 = (short) (0.5 * (Real) DIM1);
	ii2 = (short) (0.5 * (Real) DIM2);
	ii3 = (short) (0.5 * (Real) DIM3);

	// initialise
	for (i = 1; i <= DIM1; i++)
		for (j = 1; j <= DIM2; j++)
			for (k = 1; k <= DIM3; k++){
				specDens[i][j][k] = 0.;
			}

		for (i1 = 1; i1 <= ii1+1; i1++) { // up to Nyquist frequency

				i = i1;
				x = (Real) (i-1);

			for (i2 = 1; i2 <= ii2+1; i2++) { // up to Nyquist frequency

				j = i2;
				y = (Real) (j-1);

				for (i3 = 1; i3 <= ii3+1; i3++) { // up to Nyquist frequency

					k = i3;
					z = (Real) (k-1);

					if ( (i1 > 1) || (i2 > 1) || (i3 > 1) ) { // avoid w = 0


						w = (Real) sqrt(x*x + y*y + z*z);

						// take into account cut-off at minimum wavelength KMIN
						specDens[i1][i2][i3] = ( (w >= KMIN) ? sed(w, (0.5 * (BETA - numDims + 1))) : 0. );

						if ((i1 > 1) && (i1 <= ii1)) { // fold 1 wrt i

							specDens[DIM1+2-i1][i2][i3] = specDens[i1][i2][i3];

							if ((i2 > 1) && (i2 <= ii2)) {  // fold 1 wrt i and j

								specDens[DIM1+2-i1][DIM2+2-i2][i3] = specDens[i1][i2][i3];

								if ((i3 > 1) && (i3 <= ii3)){ // fold 1 wrt to i, j, k
									specDens[DIM1+2-i1][DIM2+2-i2][DIM3+2-i3] =
									specDens[i1][i2][i3];
								}

							}

							if ((i3 > 1) && (i3 <= ii3)){ // fold 1 wrt i and k
								specDens[DIM1+2-i1][i2][DIM3+2-i3] = specDens[i1][i2][i3];
							}

						}
						
						if ((i2 > 1) && (i2 <= ii2)) { // fold 1 wrt j

							specDens[i1][DIM2+2-i2][i3] = specDens[i1][i2][i3];

							if ((i3 > 1) && (i3 <= ii3)){ // fold 1 wrt j and k
								specDens[i1][DIM2+2-i2][DIM3+2-i3] = specDens[i1][i2][i3];
							}

						}

						if ((i3 > 1) && (i3 <= ii3)){ // fold 1 wrt k
							specDens[i1][i2][DIM3+2-i3] = specDens[i1][i2][i3];
						}

					} // avoid w = 0

				} // k-loop

			} // j-loop

		} // i-loop

	// normalise target spectral density
	normFactSpecDensIni = 0.;

	// compute normalisation constant to preserve variance
	// taking into account cut-off at minimum wavelength KMIN
	for (i = 1; i <= DIM1; i++)
		for (j = 1; j <= DIM2; j++)
			for (k = 1; k <= DIM3; k++)
				normFactSpecDensIni += (double)((double)specDens[i][j][k]*(double)specDens[i][j][k]);

  	normFactSpecDensIni = sqrt((double) (DIM1*DIM2*DIM3) / normFactSpecDensIni);

	// normalisation to preserve variance
	for (i = 1; i <= DIM1; i++)
		for (j = 1; j <= DIM2; j++)
			for (k = 1; k <= DIM3; k++)
				specDens[i][j][k] *= (Real) normFactSpecDensIni;

	// reset value at k = 0
	specDens[1][1][1] = 1.0; // condition to preserve mean value of field

}
#endif

#endif // 1D / 2D / 3D
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// ASCII OUTPUT: TARGET SPECTRAL DENSITY D(k)
{
// #ifdef SIM1D
// 	
// 	for(k = 1; k <= DIM3; k++)
// 		printf("%d %e %e\n", k, specDens[k], normFactSpecDensIni);
// 
// 	return(0);
// 
// #else // 2D / 3D
// 
// #ifdef SIM2D
// 
// 		for (j = 1; j <= DIM2; j++)
// 			for (k = 1; k <= DIM3; k++)
// 				printf("%d %d %e %e\n",
// 				j, k, specDens[1][j][k], normFactSpecDensIni);
// 
// #else // 3D
// 
// 	for (i = 1; i <= DIM1; i++)
// 		for (j = 1; j <= DIM2; j++)
// 			for (k = 1; k <= DIM3; k++)
// 				printf("%d %d %d %e %e\n",
// 				i, j, k, specDens[i][j][k], normFactSpecDensIni);
// 
// #endif
// 
// 	return(0);
// 
// #endif
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// START ITERATION

	for(iterNum = 0; iterNum <= MAX_ITER ; iterNum++) {

// ---------------------------------------------------------------------------------------
// Filter (apodize) white noise with spectral density D(k)
// taking into account cut-off at KMIN
{
#ifdef SIM1D

	for (k = 1; k <= DIM3; k++) {

		noiseFiltered[2*k-1] = noise[2*k-1] * specDens[k];	// real part
		noiseFiltered[2*k] = noise[2*k] * specDens[k];		// imaginary part

	} // i1 loop

#else // 2D / 3D

	for (i = 1; i <= DIM1; i++)
		for (j = 1; j <= DIM2; j++)
			for (k = 1; k <= DIM3; k++){
				noiseFiltered[i][j][2*k-1] = noise[i][j][2*k-1] * specDens[i][j][k];	// real part
				noiseFiltered[i][j][2*k] = noise[i][j][2*k] * specDens[i][j][k];		// imaginary part
			}

#endif
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// COMPUTE / OUTPUT POWER SPECTRUM P(k) OF FILTERED GAUSSIAN RANDOM FIELD
// ACTUALLY: PERIODOGRAM; BEWARE OF INNACURACIES
// COMPARE THIS TO THE TARGET *ISOTROPIC* POWER SPECTRUM
// HERE, THE POWER SPECTRUM OVER SOLID ANGLE IS COMPUTED AS AN AVERAGE OVER ANNULAR SEGMENTS
// OF INCREASING RADIUS (WAVENUMBER) AND CONSTANT WIDTH
{
#ifdef SIM1D // 1D
{
// IMPORTANT: Here, in contrast to 2D and 3D, may consider to average over a given binsize
// to reduce noise effects when displaying the power spectrum; not necessary in 2D / 3D because
the number of points used to compute P(k) at each k increases like N and N^2

	ii3 = (short) (0.5 * (Real) DIM3);

#if NYQUIST
	powSpecGauss = vector(1,ii3+1); // + 1 to include Nyquist frequency
#else
	powSpecGauss = vector(1,ii3);
#endif

	// normalisation to recover *input* power spectrum; a factor 1 / N * sigmaGauss^2
	// is due to the average spectrum of the initially white noise

	normFact = 1. / ((Real) (DIM3 * sigmaGauss2)); 

#if NYQUIST
	for(k = 1; k <= ii3+1; k++) {//only to half the grid size; INCLUDE Nyquist frequency
#else
	for(k = 1; k <= ii3; k++) {// only to half the grid size; OMIT Nyquist frequency
#endif
										// note that the other half of the computing grid
										// corresponds to negative wavenumbers and these
										// contribute to their positive counterpart
										// we are hence computing the one-sided power
										// spectrum (see Num.Rec. eq. 12.0.14)
		powSpecGauss[k] = 0;

		sumTerms = 0;

		// sum over "annular" element
		powSpecGauss[k] +=
		noiseFiltered[2*k-1] * noiseFiltered[2*k-1];	// real squared; at +wavenumb.
		powSpecGauss[k] +=
		noiseFiltered[2*k] * noiseFiltered[2*k];		// imaginary squared; at +wavenumb.
		sumTerms += 1;

		if ((k > 1) && (k <= ii3)) { // fold at frequencies below Nyquist frequency
			powSpecGauss[k] +=
			noiseFiltered[2*(DIM3+2-k)-1] * noiseFiltered[2*(DIM3+2-k)-1];	// -wavenumb.
			powSpecGauss[k] +=
			noiseFiltered[2*(DIM3+2-k)] * noiseFiltered[2*(DIM3+2-k)];		// -wavenumb.
			sumTerms += 1;
		}

		powSpecGauss[k] *= normFact / sumTerms;	// average over added terms

	} // k-loop

#if BINARYOUT
{
	// set output file name
#if (DEBUG)
	sprintf(fileName, "%s%s_gauss%d%s", "powSpec", fileInfix, iterNum, fileExt);
#else
	sprintf(fileName, "%s%s_gauss%s", "powSpec", fileInfix, fileExt);
#endif // DEBUG

	if ( (filePtr = fopen(fileName, "w")) == NULL) {
		printf("ERROR: could not open file '%s'.\n", fileName);
		return(0);
	}

#if NYQUIST
	fileErr = fwrite(&powSpecGauss[1], sizeof(Real), ii3+1, filePtr); //Nyquist
#else
	fileErr = fwrite(&powSpecGauss[1], sizeof(Real), ii3, filePtr);
#endif

#if NYQUIST
	if ( fileErr != ii3+1) { //Nyquist
#else
	if ( fileErr != ii3) {
#endif
		printf("ERROR: could not complete write to '%s'.\n", fileName);
		return(0);
	}

	// DON'T FORGET!!!
	fclose(filePtr);

}
#endif // BINARYOUT

} // SIM1D
#else // 2D / 3D

#ifdef SIM2D
{
	ii2 = (short) (0.5 * (Real) DIM2);

#if NYQUIST
	powSpecGauss = vector(1,ii2+1); // + 1 to include Nyquist frequency
#else
	powSpecGauss = vector(1,ii2);
#endif

	normFact = 1. / ((Real) (DIM2 * DIM3) * sigmaGauss2); 

#if NYQUIST
	for (m = 1; m <= ii2+1; m++) {// only to half the grid; INCLUDE Nyquist frequency;
#else
	for (m = 1; m <= ii2; m++) {	// only to half the grid; OMIT Nyquist frequency;
#endif
											// recall that the other half corresponds to negative
											// frequencies (or wavenumbers), which will contribute
											// to the same value of P(k) since k = |k|
		powSpecGauss[m] = 0;

		sumTerms = 0;

		for (j = 1; j <= m; j++) {
			for (k = 1; k <= m; k++) {

				x = sqrt( (Real) ((k-1)*(k-1) + (j-1)*(j-1)));

				if ( (j == m) || (k == m) ) { // avoid interior points

					// point 1
					powSpecGauss[m] += 
						noiseFiltered[1][j][2*k-1] *
						noiseFiltered[1][j][2*k-1];
					powSpecGauss[m] += 
						noiseFiltered[1][j][2*k] *
						noiseFiltered[1][j][2*k];
					sumTerms += 1;

					if ((j > 1) && (j <= ii2)) { // fold 1 wrt j at frequency below Nyquist

						powSpecGauss[m] += 
							noiseFiltered[1][DIM2+2-j][2*k-1] *
							noiseFiltered[1][DIM2+2-j][2*k-1];
						powSpecGauss[m] += 
							noiseFiltered[1][DIM2+2-j][2*k] *
							noiseFiltered[1][DIM2+2-j][2*k];
						sumTerms += 1;

						if ((k > 1) && (k <= ii2)) { // fold 1 wrt j and k

							powSpecGauss[m] += 
								noiseFiltered[1][DIM2+2-j][2*(DIM3+2-k)-1] *
								noiseFiltered[1][DIM2+2-j][2*(DIM3+2-k)-1];
							powSpecGauss[m] += 
								noiseFiltered[1][DIM2+2-j][2*(DIM3+2-k)] *
								noiseFiltered[1][DIM2+2-j][2*(DIM3+2-k)];
							sumTerms += 1;

						} // 1< j, k < ii2

					} // 1 < j < ii2

					if ((k > 1) && (k <= ii2)) { // fold 1 wrt k

						powSpecGauss[m] += 
							noiseFiltered[1][j][2*(DIM3+2-k)-1] *
							noiseFiltered[1][j][2*(DIM3+2-k)-1];
						powSpecGauss[m] += 
							noiseFiltered[1][j][2*(DIM3+2-k)] *
							noiseFiltered[1][j][2*(DIM3+2-k)];
						sumTerms += 1;

					} // 1 < k < ii2

				} // avoid interior points

			} // k-loop

		} // j-loop

		// the number 8*(m-1) is the *expected* number of terms to
		// contribute to each m > 1; if equal to sumTerms, everything is 0K!
// 		printf("%d %d %d\n", m, sumTerms, 8*(m-1));

		powSpecGauss[m] *= normFact; // average over ensemble
		
		powSpecGauss[m] *= x / ((Real) sumTerms); // average over all points in annulus

	} // m-loop

#if BINARYOUT
{
	// set output file name
#if (DEBUG)
	sprintf(fileName, "%s%s_gauss%d%s", "powSpec", fileInfix, iterNum, fileExt);
#else
	sprintf(fileName, "%s%s_gauss%s", "powSpec", fileInfix, fileExt);
#endif // DEBUG

	if ( (filePtr = fopen(fileName, "w")) == NULL) {
		printf("ERROR: could not open file '%s'.\n", fileName);
		return(0);
	}

#if NYQUIST
	fileErr = fwrite(&powSpecGauss[1], sizeof(Real), ii2+1, filePtr); // Nyquist
#else
	fileErr = fwrite(&powSpecGauss[1], sizeof(Real), ii2, filePtr);
#endif

#if NYQUIST
	if ( fileErr != ii2+1) { // Nyquist
#else
	if ( fileErr != ii2) {
#endif
		printf("ERROR: could not complete write to '%s'.\n", fileName);
		return(0);
	}

	// DON'T FORGET!!!
	fclose(filePtr);

}
#endif // BINARYOUT

} // SIM2D
#else // 3D
{
	ii1 = (short) (0.5 * (Real) DIM1);

#if NYQUIST
	powSpecGauss = vector(1,ii1+1); // + 1 to include Nyquist frequency
#else
	powSpecGauss = vector(1,ii1);
#endif

	normFact = 1. / ((Real) (DIM1 * DIM2 * DIM3) * sigmaGauss2); 

#if NYQUIST
	for (m = 1; m <= ii1+1; m++) {// only to half the grid; INCLUDE Nyquist frequency;
#else
	for (m = 1; m <= ii1; m++) {	// only to half the grid; OMIT Nyquist frequency;
#endif
											// recall that the other half corresponds to negative
											// frequencies (or wavenumbers), which will contribute
											// to the same value of P(k) since k = |k|
		powSpecGauss[m] = 0;

		sumTerms = 0;

		for (i = 1; i <= m; i++) {
			for (j = 1; j <= m; j++) {
				for (k = 1; k <= m; k++) {

					x = sqrt( (Real) ((k-1)*(k-1) + (j-1)*(j-1) + (i-1)*(i-1)));

					if ( (i == m) || (j == m) || (k == m) ) { // avoid interior points

						// point 1
						powSpecGauss[m] += 
							noiseFiltered[i][j][2*k-1] *
							noiseFiltered[i][j][2*k-1];
						powSpecGauss[m] += 
							noiseFiltered[i][j][2*k] *
							noiseFiltered[i][j][2*k];
						sumTerms += 1;

						if ((i > 1) && (i <= ii1)){ // fold 1 wrt i

							powSpecGauss[m] += 
								noiseFiltered[DIM1+2-i][j][2*k-1] *
								noiseFiltered[DIM1+2-i][j][2*k-1];
							powSpecGauss[m] += 
								noiseFiltered[DIM1+2-i][j][2*k] *
								noiseFiltered[DIM1+2-i][j][2*k];
							sumTerms += 1;

							if ((j > 1) && (j <= ii1)) {  // fold 1 wrt i and j

								powSpecGauss[m] += 
									noiseFiltered[DIM1+2-i][DIM2+2-j][2*k-1] *
									noiseFiltered[DIM1+2-i][DIM2+2-j][2*k-1];
								powSpecGauss[m] += 
									noiseFiltered[DIM1+2-i][DIM2+2-j][2*k] *
									noiseFiltered[DIM1+2-i][DIM2+2-j][2*k];
								sumTerms += 1;

								if ((k > 1) && (k <= ii1)) { // fold 1 wrt to i, j, k

									powSpecGauss[m] += 
										noiseFiltered[DIM1+2-i][DIM2+2-j][2*(DIM3+2-k)-1] *
										noiseFiltered[DIM1+2-i][DIM2+2-j][2*(DIM3+2-k)-1];
									powSpecGauss[m] += 
										noiseFiltered[DIM1+2-i][DIM2+2-j][2*(DIM3+2-k)] *
										noiseFiltered[DIM1+2-i][DIM2+2-j][2*(DIM3+2-k)];
									sumTerms += 1;

								} // 1 < i, j, k < ii1

							} // 1 < i, j < ii1

							if ((k > 1) && (k <= ii1)) {  // fold 1 wrt i and k

								powSpecGauss[m] += 
									noiseFiltered[DIM1+2-i][j][2*(DIM3+2-k)-1] *
									noiseFiltered[DIM1+2-i][j][2*(DIM3+2-k)-1];
								powSpecGauss[m] += 
									noiseFiltered[DIM1+2-i][j][2*(DIM3+2-k)] *
									noiseFiltered[DIM1+2-i][j][2*(DIM3+2-k)];
								sumTerms += 1;

							} // 1 < i, k < ii1

						} // 1 < i < ii1
						
						if ((j > 1) && (j <= ii1)) { // fold 1 wrt j

							powSpecGauss[m] += 
								noiseFiltered[i][DIM2+2-j][2*k-1] *
								noiseFiltered[i][DIM2+2-j][2*k-1];
							powSpecGauss[m] += 
								noiseFiltered[i][DIM2+2-j][2*k] *
								noiseFiltered[i][DIM2+2-j][2*k];
							sumTerms += 1;

							if ((k > 1) && (k <= ii1)){ // fold 1 wrt j and k

								powSpecGauss[m] += 
									noiseFiltered[i][DIM2+2-j][2*(DIM3+2-k)-1] *
									noiseFiltered[i][DIM2+2-j][2*(DIM3+2-k)-1];
								powSpecGauss[m] += 
									noiseFiltered[i][DIM2+2-j][2*(DIM3+2-k)] *
									noiseFiltered[i][DIM2+2-j][2*(DIM3+2-k)];
								sumTerms += 1;

							} // 1 < j, k < ii1

						} // 1 < j < ii1

						if ((k > 1) && (k <= ii1)){ // fold 1 wrt k

							powSpecGauss[m] += 
								noiseFiltered[i][j][2*(DIM3+2-k)-1] *
								noiseFiltered[i][j][2*(DIM3+2-k)-1];
							powSpecGauss[m] += 
								noiseFiltered[i][j][2*(DIM3+2-k)] *
								noiseFiltered[i][j][2*(DIM3+2-k)];
							sumTerms += 1;

						} // 1 < k < ii1

					} // avoid interior points

				} // k-loop

			} // j-loop

		} // i-loop

		// the number 2*(12*(m-1)*(m-1) + 1) is the *expected* number of terms to
		// contribute to each m > 1; if equal to sumTerms, everything is 0K!
// 		printf("%d %d %d\n", m, sumTerms, 2*(12*(m-1)*(m-1) + 1));

		powSpecGauss[m] *= normFact; // average over ensemble
		
		powSpecGauss[m] *= x*x / ((Real) sumTerms); // average over all points in annulus


	} // m-loop

#if BINARYOUT
{
	// set output file name
#if (DEBUG)
	sprintf(fileName, "%s%s_gauss%d%s", "powSpec", fileInfix, iterNum, fileExt);
#else
	sprintf(fileName, "%s%s_gauss%s", "powSpec", fileInfix, fileExt);
#endif // DEBUG

	if ( (filePtr = fopen(fileName, "w")) == NULL) {
		printf("ERROR: could not open file '%s'.\n", fileName);
		return(0);
	}

#if NYQUIST
	fileErr = fwrite(&powSpecGauss[1], sizeof(Real), ii1+1, filePtr); // Nyquist
#else
	fileErr = fwrite(&powSpecGauss[1], sizeof(Real), ii1, filePtr);
#endif

#if NYQUIST
	if ( fileErr != ii1+1) { // Nyquist
#else
	if ( fileErr != ii1) {
#endif
		printf("ERROR: could not complete write to '%s'.\n", fileName);
		return(0);
	}

	// DON'T FORGET!!!
	fclose(filePtr);

}
#endif // BINARYOUT

} // SIM3D

#endif // 2D / 3D

#endif // 1 or 2D / 3D?

}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// COMPUTE *inverse* Fast Fourier Transform of filtered Gaussian noise
{
// IMPORTANT: normalisation factor usually includes the sampling rate Delta in each
// dimension, which is necessary to compare the result from the Num.Rec.FFT routine to the
// corresponding result in the continuous case, but it is ignored here because is implicitly
// set to 1; see Num.Recipes, their eq. 12.1.18

#ifdef SIM1D // 1D

	dfour1(noiseFiltered,DIM3,-1);

	// factor 1/(N) not taken into account by NR routine
	normFact = (1. / (Real) DIM3);

	for (k = 1; k <= DIM3; k++){
		noiseFiltered[2*k-1] *= normFact;	// real part
		noiseFiltered[2*k] *= normFact; 		// imag part
	}

#else // 2D / 3D

	dfourn(&noiseFiltered[1][1][1] - 1,dims,3,-1);

	// factor 1/(N1*N2*N3) not taken into account by NR routine
	normFact = (1. / (Real) (DIM1 * DIM2 * DIM3));

	for(i = 1; i <= DIM1; i++)
		for(j = 1; j <= DIM2; j++)
			for(k = 1; k <= DIM3; k++) {
				noiseFiltered[i][j][2*k-1] *= normFact;	// real part
				noiseFiltered[i][j][2*k] *= normFact; 		// imag part
			}

#endif // end of 1D, 2D, or 3D?
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// RENORMALISE GAUSSIAN RANDOM FIEL (i.e. shift and scale) to enforce target statistics
{
#ifdef SIM1D

	// compute actual mean / variance
	muAux = average(noiseFiltered,DIM3);
	sigmaAux = stddev(noiseFiltered,DIM3);

#if (DEBUG)
	printf("# current mu (GAUSS) is %lf\n", muAux);
	printf("# current sigma (GAUSS) is %lf\n", sigmaAux);
#endif // DEBUG

	sigmaAux = 1./sigmaAux;

	for(k = 1; k <= DIM3; k++) {

		noiseFiltered[2*k-1] = sigmaAux*(noiseFiltered[2*k-1] - muAux); // real part
		noiseFiltered[2*k] = sigmaAux*(noiseFiltered[2*k] - muAux); // imaginary part

	}

	// re-compute actual mean / variance
	muAux = average(noiseFiltered,DIM3);
	sigmaAux = stddev(noiseFiltered,DIM3);

#else // 2D / 3D 

	// compute actual mean / variance
	muAux = average(noiseFiltered,DIM1,DIM2,DIM3);
	sigmaAux = stddev(noiseFiltered,DIM1,DIM2,DIM3);

#if (DEBUG)
	printf("# current mu (GAUSS) is %lf\n", muAux);
	printf("# current sigma (GAUSS) is %lf\n", sigmaAux);
#endif // DEBUG

	sigmaAux = 1./sigmaAux;

	// renormalise data
	for(i = 1; i <= DIM1; i++)
		for(j = 1; j <= DIM2; j++)
			for(k = 1; k <= DIM3; k++) {

				noiseFiltered[i][j][2*k-1] = sigmaAux*(noiseFiltered[i][j][2*k-1] - muAux);// + muGauss; // real part
				noiseFiltered[i][j][2*k] = sigmaAux*(noiseFiltered[i][j][2*k] - muAux);// + muGauss; // imaginary part

			}

	// re-compute actual mean / variance
	muAux = average(noiseFiltered,DIM1,DIM2,DIM3);
	sigmaAux = stddev(noiseFiltered,DIM1,DIM2,DIM3);

#endif // SIM2D / SIM3D

#if (DEBUG)
	printf("# adjusted mu (GAUSS) is %lf\n", muAux);
	printf("# adjusted sigma (GAUSS) is %lf\n", sigmaAux);
#endif // DEBUG



}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// BINARY OUTPUT: FILTERED GAUSSIAN RANDOM FIELD
// NOTE: TWO FILES ARE WRITTEN: ONE CONTAINING BOTH THE REAL AND IMAGINARY PARTS OF THE
// FIELD, AND ONE CONTAINING THE REAL PART ONLY. THE LATTER SHALL PROVE MORE USEFUL FOR
// PRACTICAL APPLICATIONS. NOTE ALSO THAT, BECAUSE THE INITIAL GAUSSIAN FIELD IS REAL,
// THE IMAGINARY PART IS TRIVIALLY ZERO.
// NOTE: will be overwritten at each iteration!
#if BINARYOUT
{

#ifdef SIM1D

	// 1: FULL (REAL AND IMAGINARY PARTS) FIELD

		// set output file name
		sprintf(fileName, "%s%s_filt%s", "gauss", fileInfix, fileExt);

		filePtr = NULL;
		if ( (filePtr = fopen(fileName, "w")) == NULL) {
			printf("ERROR: could not open file '%s'.\n", fileName);
			return(0);
		}

		fileErr = fwrite(&noiseFiltered[1], sizeof(Real), 2*DIM3, filePtr);

		if ( fileErr != (2*DIM3) ) {
			printf("ERROR: could not complete write to '%s'.\n", fileName);
			return(0);
		}

		// DON'T FORGET!!!
		fclose(filePtr);

	// 2: REAL FIELD

		// set output file name
		sprintf(fileName, "%s%s_filt_Real%s", "gauss", fileInfix, fileExt);

		filePtr = NULL;
		if ( (filePtr = fopen(fileName, "w")) == NULL) {
			printf("ERROR: could not open file '%s'.\n", fileName);
			return(0);
		}

		fileErr = 0;
		for (k = 1; k <= DIM3; k++)
			fileErr += fwrite(&noiseFiltered[2*k-1], sizeof(Real), 1, filePtr);

		if ( fileErr != DIM3 ) {
			printf("ERROR: could not complete write to '%s'.\n", fileName);
			return(0);
		}

		// DON'T FORGET!!!
		fclose(filePtr);

#else // 2D / 3D

	// 1: FULL (REAL AND IMAGINARY PARTS) FIELD

		// set output file name
		sprintf(fileName, "%s%s_filt%s", "gauss", fileInfix, fileExt);

		filePtr = NULL;
		if ( (filePtr = fopen(fileName, "w")) == NULL) {
			printf("ERROR: could not open file '%s'.\n", fileName);
			return(0);
		}

		fileErr = 0;

		for (i = 1; i <= DIM1; i++)
			for (j = 1; j <= DIM2; j++)
				for (k = 1; k <= DIM3; k++) {
					fileErr += fwrite(&noiseFiltered[i][j][2*k-1], sizeof(Real), 1, filePtr);
					fileErr += fwrite(&noiseFiltered[i][j][2*k], sizeof(Real), 1, filePtr);
				}

		if ( fileErr != (DIM1*DIM2*2*DIM3) ) {
			printf("ERROR: could not complete write to '%s'.\n", fileName);
			return(0);
		}

		// DON'T FORGET!!!
		fclose(filePtr);

	// 2: REAL FIELD

		// set output file name
		sprintf(fileName, "%s%s_filt_Real%s", "gauss", fileInfix, fileExt);

		filePtr = NULL;
		if ( (filePtr = fopen(fileName, "w") ) == NULL) {
			printf("ERROR: could not open file '%s'.\n", fileName);
			return(0);
		}

		fileErr = 0;

		for (i = 1; i <= DIM1; i++)
			for (j = 1; j <= DIM2; j++)
				for (k = 1; k <= DIM3; k++)
					fileErr += fwrite(&noiseFiltered[i][j][2*k-1], sizeof(Real), 1, filePtr);

		if ( fileErr != (DIM1*DIM2*DIM3) ) {
			printf("ERROR: could not complete write to '%s'.\n", fileName);
			return(0);
		}

		// DON'T FORGET!!!
		fclose(filePtr);

#endif

}
#endif // BINARYOUT
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// NOW WE HAVE A GUESS FOR THE SEED GAUSSIAN RANDOM FIELD
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// COMPUTE 'target' log-normal random field, scaling and shifting by sigma and mu, resp.,
// taking into account cut-off at KMIN
{
#ifdef SIM1D
	
	for(k = 1; k <= DIM3; k++) {
		logNormal[2*k-1] =
			exp((sigmaGauss * noiseFiltered[2*k-1] + muGauss));// real part
		logNormal[2*k] =
			exp((sigmaGauss * noiseFiltered[2*k] + muGauss));	// imaginary part
	}

#else // 2D / 3D

	for (i = 1; i <= DIM1; i++)
		for (j = 1; j <= DIM2; j++)
			for (k = 1; k <= DIM3; k++) {
				logNormal[i][j][2*k-1] =
					exp((sigmaGauss * noiseFiltered[i][j][2*k-1] + muGauss));// real part
				logNormal[i][j][2*k] =
					exp((sigmaGauss * noiseFiltered[i][j][2*k] + muGauss));	// imaginary part
			}
#endif
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// COMPUTE *forward* Fast Fourier Transform of log-normal random field
{
#ifdef SIM1D

	dfour1(logNormal,DIM3,1);

#else // 2D/3D

	dfourn(&logNormal[1][1][1] - 1,dims,3,1); // test with complex input data array

#endif
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// COMPUTE / OUTPUT POWER SPECTRUM P(k) OF FILTERED LOG-NORMAL RANDOM FIELD
// ACTUALLY: PERIODOGRAM; BEWARE OF INNACURACIES
// COMPARE THIS TO THE TARGET *ISOTROPIC* POWER SPECTRUM
// HERE, THE POWER SPECTRUM OVER SOLID ANGLE IS COMPUTED AS AN AVERAGE OVER ANNULAR SEGMENTS OF
// INCREASING RADIUS (WAVENUMBER) AND CONSTANT WIDTH
{
#ifdef SIM1D
{
// IMPORTANT: Here, in contrast to 2D and 3D, may consider to average over a given binsize to reduce noise effects when displaying the power spectrum; not necessary in 2D / 3D because the number of points used to compute P(k) at each k increases like N and N^2

	ii3 = (short) (0.5 * (Real) DIM3);

#if NYQUIST
	powSpecLogNorm = vector(1,ii3+1); // + 1 to include Nyquist frequency
#else
	powSpecLogNorm = vector(1,ii3);
#endif

	// normalisation to recover *input* power spectrum;
	// a factor 1 / N * sigmaLogn^2 is due to the average spectrum of the initially white noise
	normFact = 1. / ((Real) (DIM3 * sigmaLogn2)); 

#if NYQUIST
	for(k = 1; k <= ii3+1; k++) {//only to half the grid size; INCLUDE Nyquist frequency
#else
	for(k = 1; k <= ii3; k++) {	// only to half the grid size; OMIT Nyquist frequency
#endif
										// note that the other half of the computing grid
										// corresponds to negative wavenumbers and these
										// contribute to their positive counterpart
										// we are hence computing the one-sided power
										// spectrum (see Num.Rec. eq. 12.0.14)
		powSpecLogNorm[k] = 0;

		sumTerms = 0;

		// sum over "annular" element
		powSpecLogNorm[k] +=
		logNormal[2*k-1] * logNormal[2*k-1];	// real squared; at +wavenumb.
		powSpecLogNorm[k] +=
		logNormal[2*k] * logNormal[2*k];		// imaginary squared; at +wavenumb.
		sumTerms += 1;

		if ((k > 1) && (k <= ii3)) { // fold at frequencies below Nyquist frequency
			powSpecLogNorm[k] +=
			logNormal[2*(DIM3+2-k)-1] * logNormal[2*(DIM3+2-k)-1];	// -wavenumb.
			powSpecLogNorm[k] +=
			logNormal[2*(DIM3+2-k)] * logNormal[2*(DIM3+2-k)];		// -wavenumb.
			sumTerms += 1;
		}

		powSpecLogNorm[k] *= normFact / sumTerms;	// average over added terms

	} // i-loop

#if BINARYOUT
{
	// set output file name
#if (DEBUG)
	sprintf(fileName, "%s%s_lognorm%d%s", "powSpec", fileInfix, iterNum, fileExt);
#else
// NOTE: will be overwritten at each iteration!
	sprintf(fileName, "%s%s_lognorm%s", "powSpec", fileInfix, fileExt);
#endif // DEBUG

	if ( (filePtr = fopen(fileName, "w")) == NULL) {
		printf("ERROR: could not open file '%s'.\n", fileName);
		return(0);
	}

#if NYQUIST
	fileErr = fwrite(&powSpecLogNorm[1], sizeof(Real), ii3+1, filePtr); // Nyquist
#else
	fileErr = fwrite(&powSpecLogNorm[1], sizeof(Real), ii3, filePtr);
#endif

#if NYQUIST
	if ( fileErr != ii3+1) { // Nyquist
#else
	if ( fileErr != ii3) {
#endif
		printf("ERROR: could not complete write to '%s'.\n", fileName);
		return(0);
	}

	// DON'T FORGET!!!
	fclose(filePtr);

}
#endif // BINARYOUT

} // SIM1D

#else // 2D / 3D

#ifdef SIM2D
{
	ii2 = (short) (0.5 * (Real) DIM2);

#if NYQUIST
	powSpecLogNorm = vector(1,ii2+1); // + 1 to include Nyquist frequency
#else
	powSpecLogNorm = vector(1,ii2);
#endif

	// normalisation to recover *input* power spectrum;
	// a factor 1 / N^2 * sigmaLogn^2 is due to the average spectrum of the initially white noise
	normFact = 1. / ((Real) (DIM2 * DIM3) * sigmaLogn2); 

#if NYQUIST
	for (m = 1; m <= ii2+1; m++) {// only to half the grid; INCLUDE Nyquist frequency;
#else
	for (m = 1; m <= ii2; m++) {// only to half the grid; OMIT Nyquist frequency;
#endif
											// recall that the other half corresponds to negative
											// frequencies (or wavenumbers), which will contribute
											// to the same value of P(k) since k = |k|
		powSpecLogNorm[m] = 0;

		sumTerms = 0;

		for (j = 1; j <= m; j++) {
			for (k = 1; k <= m; k++) {

				x = sqrt( (Real) ((k-1)*(k-1) + (j-1)*(j-1)));

				if ( (j == m) || (k == m) ) { // avoid interior points

					// point 1
					powSpecLogNorm[m] += 
						logNormal[1][j][2*k-1] *
						logNormal[1][j][2*k-1];
					powSpecLogNorm[m] += 
						logNormal[1][j][2*k] *
						logNormal[1][j][2*k];
					sumTerms += 1;

					if ((j > 1) && (j <= ii2)) { // fold 1 wrt j at frequency below Nyquist

						powSpecLogNorm[m] += 
							logNormal[1][DIM2+2-j][2*k-1] *
							logNormal[1][DIM2+2-j][2*k-1];
						powSpecLogNorm[m] += 
							logNormal[1][DIM2+2-j][2*k] *
							logNormal[1][DIM2+2-j][2*k];
						sumTerms += 1;

						if ((k > 1) && (k <= ii2)) { // fold 1 wrt j and k

							powSpecLogNorm[m] += 
								logNormal[1][DIM2+2-j][2*(DIM3+2-k)-1] *
								logNormal[1][DIM2+2-j][2*(DIM3+2-k)-1];
							powSpecLogNorm[m] += 
								logNormal[1][DIM2+2-j][2*(DIM3+2-k)] *
								logNormal[1][DIM2+2-j][2*(DIM3+2-k)];
							sumTerms += 1;

						} // 1< j, k < ii2

					} // 1 < j < ii2

					if ((k > 1) && (k <= ii2)) { // fold 1 wrt k

						powSpecLogNorm[m] += 
							logNormal[1][j][2*(DIM3+2-k)-1] *
							logNormal[1][j][2*(DIM3+2-k)-1];
						powSpecLogNorm[m] += 
							logNormal[1][j][2*(DIM3+2-k)] *
							logNormal[1][j][2*(DIM3+2-k)];
						sumTerms += 1;

					} // 1 < k < ii2

				} // avoid interior points

			} // k-loop

		} // j-loop

		// the number 8*(m-1) is the *expected* number of terms to
		// contribute to each m > 1; if equal to sumTerms, everything is 0K!
// 		printf("%d %d %d\n", m, sumTerms, 8*(m-1));

		powSpecLogNorm[m] *= normFact; // average over ensemble
		
		powSpecLogNorm[m] *= x / ((Real) sumTerms); // average over all points in annulus

	} // m-loop

#if BINARYOUT
{
	// set output file name
#if (DEBUG)
	sprintf(fileName, "%s%s_lognorm%d%s", "powSpec", fileInfix, iterNum, fileExt);
#else
	sprintf(fileName, "%s%s_lognorm%s", "powSpec", fileInfix, fileExt);
#endif // DEBUG

	if ( (filePtr = fopen(fileName, "w")) == NULL) {
		printf("ERROR: could not open file '%s'.\n", fileName);
		return(0);
	}

#if NYQUIST
	fileErr = fwrite(&powSpecLogNorm[1], sizeof(Real), ii2+1, filePtr); // Nyquist
#else
	fileErr = fwrite(&powSpecLogNorm[1], sizeof(Real), ii2, filePtr);
#endif

#if NYQUIST
	if ( fileErr != ii2+1) { // Nyquist
#else
	if ( fileErr != ii2) {
#endif
		printf("ERROR: could not complete write to '%s'.\n", fileName);
		return(0);
	}

	// DON'T FORGET!!!
	fclose(filePtr);

}
#endif // BINARYOUT

} // SIM2D
#else // 3D
{
	ii1 = (short) (0.5 * (Real) DIM1);

#if NYQUIST
	powSpecLogNorm = vector(1,ii1+1); // + 1 to include Nyquist frequency
#else
	powSpecLogNorm = vector(1,ii1);
#endif

	// normalisation to recover *input* power spectrum;
	// a factor 1 / N^3 * sigmaLogn^2 is due to the average spectrum of the initially white noise
	normFact = 1. / ((Real) (DIM1 * DIM2 * DIM3) * sigmaLogn2); 

#if NYQUIST
	for (m = 1; m <= ii1+1; m++) {// only to half the grid; INCLUDE Nyquist frequency;
#else
	for (m = 1; m <= ii1; m++) {	// only to half the grid; OMIT Nyquist frequency;
#endif
											// recall that the other half corresponds to negative
											// frequencies (or wavenumbers), which will contribute
											// to the same value of P(k) since k = |k|
		powSpecLogNorm[m] = 0;

		sumTerms = 0;

		for (i = 1; i <= m; i++) {
			for (j = 1; j <= m; j++) {
				for (k = 1; k <= m; k++) {

					x = sqrt( (Real) ((k-1)*(k-1) + (j-1)*(j-1) + (i-1)*(i-1)));

					if ( (i == m) || (j == m) || (k == m) ) { // avoid interior points

						// point 1
						powSpecLogNorm[m] += 
							logNormal[i][j][2*k-1] *
							logNormal[i][j][2*k-1];
						powSpecLogNorm[m] += 
							logNormal[i][j][2*k] *
							logNormal[i][j][2*k];
						sumTerms += 1;

						if ((i > 1) && (i <= ii1)){ // fold 1 wrt i

							powSpecLogNorm[m] += 
								logNormal[DIM1+2-i][j][2*k-1] *
								logNormal[DIM1+2-i][j][2*k-1];
							powSpecLogNorm[m] += 
								logNormal[DIM1+2-i][j][2*k] *
								logNormal[DIM1+2-i][j][2*k];
							sumTerms += 1;

							if ((j > 1) && (j <= ii1)) {  // fold 1 wrt i and j

								powSpecLogNorm[m] += 
									logNormal[DIM1+2-i][DIM2+2-j][2*k-1] *
									logNormal[DIM1+2-i][DIM2+2-j][2*k-1];
								powSpecLogNorm[m] += 
									logNormal[DIM1+2-i][DIM2+2-j][2*k] *
									logNormal[DIM1+2-i][DIM2+2-j][2*k];
								sumTerms += 1;

								if ((k > 1) && (k <= ii1)) { // fold 1 wrt to i, j, k

									powSpecLogNorm[m] += 
										logNormal[DIM1+2-i][DIM2+2-j][2*(DIM3+2-k)-1] *
										logNormal[DIM1+2-i][DIM2+2-j][2*(DIM3+2-k)-1];
									powSpecLogNorm[m] += 
										logNormal[DIM1+2-i][DIM2+2-j][2*(DIM3+2-k)] *
										logNormal[DIM1+2-i][DIM2+2-j][2*(DIM3+2-k)];
									sumTerms += 1;

								} // 1 < i, j, k < ii1

							} // 1 < i, j < ii1

							if ((k > 1) && (k <= ii1)) {  // fold 1 wrt i and k

								powSpecLogNorm[m] += 
									logNormal[DIM1+2-i][j][2*(DIM3+2-k)-1] *
									logNormal[DIM1+2-i][j][2*(DIM3+2-k)-1];
								powSpecLogNorm[m] += 
									logNormal[DIM1+2-i][j][2*(DIM3+2-k)] *
									logNormal[DIM1+2-i][j][2*(DIM3+2-k)];
								sumTerms += 1;

							} // 1 < i, k < ii1

						} // 1 < i < ii1
						
						if ((j > 1) && (j <= ii1)) { // fold 1 wrt j

							powSpecLogNorm[m] += 
								logNormal[i][DIM2+2-j][2*k-1] *
								logNormal[i][DIM2+2-j][2*k-1];
							powSpecLogNorm[m] += 
								logNormal[i][DIM2+2-j][2*k] *
								logNormal[i][DIM2+2-j][2*k];
							sumTerms += 1;

							if ((k > 1) && (k <= ii1)){ // fold 1 wrt j and k

								powSpecLogNorm[m] += 
									logNormal[i][DIM2+2-j][2*(DIM3+2-k)-1] *
									logNormal[i][DIM2+2-j][2*(DIM3+2-k)-1];
								powSpecLogNorm[m] += 
									logNormal[i][DIM2+2-j][2*(DIM3+2-k)] *
									logNormal[i][DIM2+2-j][2*(DIM3+2-k)];
								sumTerms += 1;

							} // 1 < j, k < ii1

						} // 1 < j < ii1

						if ((k > 1) && (k <= ii1)){ // fold 1 wrt k

							powSpecLogNorm[m] += 
								logNormal[i][j][2*(DIM3+2-k)-1] *
								logNormal[i][j][2*(DIM3+2-k)-1];
							powSpecLogNorm[m] += 
								logNormal[i][j][2*(DIM3+2-k)] *
								logNormal[i][j][2*(DIM3+2-k)];
							sumTerms += 1;

						} // 1 < k < ii1

					} // avoid interior points

				} // k-loop

			} // j-loop

		} // i-loop

		// the number 2*(12*(m-1)*(m-1) + 1) is the *expected* number of terms to
		// contribute to each m > 1; if equal to sumTerms, everything is 0K!
// 		printf("%d %d %d\n", m, sumTerms, 2*(12*(m-1)*(m-1) + 1));

		powSpecLogNorm[m] *= normFact; // average over ensemble
		
		powSpecLogNorm[m] *= x*x / ((Real) sumTerms); // average over all points in annulus

	} // m-loop

#if BINARYOUT
{
	// set output file name
#if (DEBUG)
	sprintf(fileName, "%s%s_lognorm%d%s", "powSpec", fileInfix, iterNum, fileExt);
#else
	sprintf(fileName, "%s%s_lognorm%s", "powSpec", fileInfix, fileExt);
#endif // DEBUG

	if ( (filePtr = fopen(fileName, "w")) == NULL) {
		printf("ERROR: could not open file '%s'.\n", fileName);
		return(0);
	}

#if NYQUIST
	fileErr = fwrite(&powSpecLogNorm[1], sizeof(Real), ii1+1, filePtr); // Nyquist
#else
	fileErr = fwrite(&powSpecLogNorm[1], sizeof(Real), ii1, filePtr);
#endif

#if NYQUIST
	if ( fileErr != ii1+1) { // Nyquist
#else
	if ( fileErr != ii1) {
#endif
		printf("ERROR: could not complete write to '%s'.\n", fileName);
		return(0);
	}

	// DON'T FORGET!!!
	fclose(filePtr);

}
#endif // BINARYOUT

} // SIM3D

#endif // 2D / 3D

#endif // 1 or 2D / 3D?

}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// COMPUTE LOGNORMAL LOG10(P(k)) vs. LOG10(k)
{
#ifdef SIM1D
{
#if NYQUIST
	i3 = (short) (0.5 * (Real) DIM3)+1; // + 1 to include Nyquist frequency
#else
	i3 = (short) (0.5 * (Real) DIM3);
#endif

	// allocate arrays
	logPowSpecLogNorm = vector(1,i3);
	logWaveNumber = vector(1,i3);
	sigmaWeight = vector(1,i3);

	for(k = 2; k <= i3; k++){ // skip wavenumber k = 0

		logWaveNumber[k] = log10((Real) (k-1));

		// take into account cut-off at minimum wavelength KMIN
		logPowSpecLogNorm[k] = ( (k >= KMIN) ? log10(powSpecLogNorm[k]) : 0. );

		sigmaWeight[k] = 1.; // set weight to whatever is appropriate

	}
}
#else
{
	#ifdef SIM2D

#if NYQUIST
		i1 = (short) (0.5 * (Real) DIM2)+1; // + 1 to include Nyquist frequency
#else
		i1 = (short) (0.5 * (Real) DIM2);
#endif

	#else

#if NYQUIST
		i1 = (short) (0.5 * (Real) DIM1)+1; // + 1 to include Nyquist frequency
#else
		i1 = (short) (0.5 * (Real) DIM1);
#endif

	#endif

	// allocate arrays
	logPowSpecLogNorm = vector(1,i1);
	logWaveNumber = vector(1,i1);
	sigmaWeight = vector(1,i1);

	for(i = 2; i <= i1; i++){ // skip k = 0

		logWaveNumber[i] = log10((Real) (i-1));

		// take into account cut-off at minimum wavelength KMIN
		logPowSpecLogNorm[i] = ( (i >= KMIN) ? log10(powSpecLogNorm[i]) : 0. );

		sigmaWeight[i] = 1.; // set weight to whatever is appropriate

	}
}
#endif
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// FIT POLYNOMIAL TO LOGNORMAL LOG10(P(k)) vs. LOG10(k)
#if (POLYFIT)
{
#ifdef SIM1D

#if NYQUIST
	i3 = (short) (0.5 * (Real) DIM3)+1; // + 1 to include Nyquist frequency
#else
	i3 = (short) (0.5 * (Real) DIM3);
#endif

#else // 2D / 3D

#ifdef SIM2D

#if NYQUIST
		i3 = (short) (0.5 * (Real) DIM2)+1; // + 1 to include Nyquist frequency
#else
		i3 = (short) (0.5 * (Real) DIM2);
#endif

#else // 2D

#if NYQUIST
		i3 = (short) (0.5 * (Real) DIM1)+1; // + 1 to include Nyquist frequency
#else
		i3 = (short) (0.5 * (Real) DIM1);
#endif

#endif // 3D

#endif // 1D / 2D / 3D

	// allocate arrays
	polyFitLogPowSpecLogNorm = vector(1,i3);
	polyVal=vector(1,numParams);
	polyCoeff=vector(1,numParams);
	polyCoeffFlag=ivector(1,numParams);
	covar=matrix(1,numParams,1,numParams);

	// initialise arrays
	for (k = 1; k <= i3; k++)
		polyFitLogPowSpecLogNorm[k] = 0.;

	for (k = 1; k <= numParams; k++) // set all parameters to be fitted
		polyCoeffFlag[k] = 1;

	// ignore term at wavenumber k = 0
	// take into account cut-off at minimum wavelength KMIN, i.e. fit only over non-zero values of P(k)
	dlfit(&logWaveNumber[1+KMIN]-1, &logPowSpecLogNorm[1+KMIN]-1, &sigmaWeight[1+KMIN]-1, (int) (i3-(KMIN)), &polyCoeff[1]-1, &polyCoeffFlag[1]-1, numParams, covar, &chisq, poly);


	// compute value of polynomial fit at x = log(k)
	for(k = 2; k <= i3; k++) {
		x = logWaveNumber[k];
		poly(x, polyVal, numParams);		
		for (ii3 = 1; ii3 <= numParams; ii3++)
			polyFitLogPowSpecLogNorm[k] +=  polyCoeff[ii3]*polyVal[ii3];
 	} // i-loop

}
#endif // POLYFIT
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// SMOOTH LOGNORMAL LOG10(P(k)) vs. LOG10(k) USING A SAVITZKY-GOLAY FILTER
#if (SGFILTER)
{
#ifdef SIM1D

#if NYQUIST
	i3 = (short) (0.5 * (Real) DIM3)+1; // + 1 to include Nyquist frequency
#else
	i3 = (short) (0.5 * (Real) DIM3);
#endif

#else // 2D / 3D

	#ifdef SIM2D

#if NYQUIST
		i3 = (short) (0.5 * (Real) DIM2)+1; // + 1 to include Nyquist frequency
#else
		i3 = (short) (0.5 * (Real) DIM2);
#endif

	#else

#if NYQUIST
		i3 = (short) (0.5 * (Real) DIM1)+1; // + 1 to include Nyquist frequency
#else
		i3 = (short) (0.5 * (Real) DIM1);
#endif

	#endif

#endif

	// allocate arrays
	sgLogPowSpecLogNorm = vector(1,i3);

	for(k = KMIN+halfWindow+1; k <= i3-halfWindow; k++) {

		sgLogPowSpecLogNorm[k] = 0.;

		for(i=-halfWindow;i<=halfWindow;i++){
			signal = ( (k+i >= KMIN+1)&&(k+i <= i3) ? logPowSpecLogNorm[k+i] : 0.);
			sgLogPowSpecLogNorm[k] +=  signal * sgCoeff[i+halfWindow+1];
		} // i-loop

// 		printf("%d %e %e\n", k, sgLogPowSpecLogNorm[k], logPowSpecLogNorm[k]);

	} // k-loop

	// reset values at edge to avoid padding issues
	for(i=2;i<=KMIN+halfWindow+1;i++)
		sgLogPowSpecLogNorm[i] = logPowSpecLogNorm[i];

	// reset values at edge to avoid padding issues
	for(i=1;i<=halfWindow;i++)
		sgLogPowSpecLogNorm[i3-halfWindow+i] = logPowSpecLogNorm[i3-halfWindow+i];

}
#endif // SGFILTER
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// FIT LINEAR MODEL TO LOGNORMAL LOG10(P(k)) vs. LOG10(k)
{
#ifdef SIM1D

#if NYQUIST
	i3 = (short) (0.5 * (Real) DIM3)+1; // + 1 to include Nyquist frequency
#else
	i3 = (short) (0.5 * (Real) DIM3);
#endif

#else

	#ifdef SIM2D

#if NYQUIST
		i3 = (short) (0.5 * (Real) DIM2)+1; // + 1 to include Nyquist frequency
#else
		i3 = (short) (0.5 * (Real) DIM2);
#endif

	#else

#if NYQUIST
		i3 = (short) (0.5 * (Real) DIM1)+1; // + 1 to include Nyquist frequency
#else
		i3 = (short) (0.5 * (Real) DIM1);
#endif

	#endif

#endif

	// Fit: 5th parameter set to 1 to use provided weights; otherwise set to 0
	// ignore term at wavenumber k = 0
	// take into account cut-off at minimum wavelength KMIN, i.e. fit only over non-zero values of P(k)
	dfit(&logWaveNumber[1+KMIN]-1, &logPowSpecLogNorm[1+KMIN]-1, (int) (i3-(KMIN)), &sigmaWeight[1+KMIN]-1, 1, &inter, &slopeLogNorm, &interErr, &slopeErr, &chi2, &qParam);

	// set normalisation constant
	normConstLogNorm = inter;

}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// COMPUTE *inverse* Fast Fourier Transform of filtered log-normal field
{
// IMPORTANT: normalisation factor usually includes the sampling rate Delta in each dimension, which is necessary to compare the result from the Num.Rec.FFT routine to the corresponding result in the continuous case, but it is ignored here because is implicitly set to 1; see Num.Recipes, their eq. 12.1.18

#ifdef SIM1D // 1D

	dfour1(logNormal,DIM3,-1);

	// factor 1/(N) not taken into account by NR routine
	normFact = (1. / (Real) DIM3);

	for (k = 1; k <= DIM3; k++){
		logNormal[2*k-1] *= normFact;	// real part
		logNormal[2*k] *= normFact; 		// imag part
	}

#else // 2D / 3D

	dfourn(&logNormal[1][1][1] - 1,dims,3,-1);

	// factor 1/(N1*N2*N3) not taken into account by NR routine
	normFact = (1. / (Real) (DIM1 * DIM2 * DIM3));

	for(i = 1; i <= DIM1; i++)
		for(j = 1; j <= DIM2; j++)
			for(k = 1; k <= DIM3; k++) {
				logNormal[i][j][2*k-1] *= normFact;	// real part
				logNormal[i][j][2*k] *= normFact; 		// imag part
			}

#endif // end of 1D, 2D, or 3D?
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// BINARY OUTPUT: FILTERED LOGNORMAL RANDOM FIELD (THIS IS THE ACTUAL INTENT OF THE CODE)
// NOTE: TWO FILES ARE WRITTEN: ONE CONTAINING BOTH THE REAL AND IMAGINARY PARTS OF THE
// FIELD, AND ONE CONTAINING THE REAL PART ONLY. THE LATTER SHALL PROVE MORE USEFUL FOR
// PRACTICAL APPLICATIONS. NOTE ALSO THAT, BECAUSE THE INITIAL GAUSSIAN FIELD IS REAL,
// THE IMAGINARY PART IS TRIVIALLY CONSTANT.
// NOTE: will be overwritten at each iteration!
#if BINARYOUT
{
#ifdef SIM1D

	// 1: FULL (REAL AND IMAGINARY PARTS) FIELD

		// set output file name
		sprintf(fileName, "%s%s_filt%s", "lognorm", fileInfix, fileExt);

		filePtr = NULL;
		if ( (filePtr = fopen(fileName, "w") ) == NULL) {
			printf("ERROR: could not open file '%s'.\n", fileName);
			return(0);
		}

		fileErr = fwrite(&logNormal[1], sizeof(Real), 2*DIM3, filePtr);

		if ( fileErr != (2*DIM3) ) {
			printf("ERROR: could not complete write to '%s'.\n", fileName);
			return(0);
		}

		// DON'T FORGET!!!
		fclose(filePtr);

	// 2: REAL FIELD

		// set output file name
		sprintf(fileName, "%s%s_filt_Real%s", "lognorm", fileInfix, fileExt);

		filePtr = NULL;
		if ( (filePtr = fopen(fileName, "w") ) == NULL) {
			printf("ERROR: could not open file '%s'.\n", fileName);
			return(0);
		}

		fileErr = 0;
		for(k = 1; k <= DIM3; k++)
			fileErr += fwrite(&logNormal[2*k-1], sizeof(Real), 1, filePtr);

		if ( fileErr != DIM3 ) {
			printf("ERROR: could not complete write to '%s'.\n", fileName);
			return(0);
		}

		// DON'T FORGET!!!
		fclose(filePtr);

#else // 2D / 3D

	// 1: FULL (REAL AND IMAGINARY PARTS) FIELD

		// set output file name
		sprintf(fileName, "%s%s_filt%s", "lognorm", fileInfix, fileExt);

		filePtr = NULL;
		if ( (filePtr = fopen(fileName, "w") ) == NULL) {
			printf("ERROR: could not open file '%s'.\n", fileName);
			return(0);
		}

		fileErr = 0;

		for (i = 1; i <= DIM1; i++)
			for (j = 1; j <= DIM2; j++)
				for (k = 1; k <= DIM3; k++) {
					fileErr += fwrite(&logNormal[i][j][2*k-1], sizeof(Real), 1, filePtr);
					fileErr += fwrite(&logNormal[i][j][2*k], sizeof(Real), 1, filePtr);
				}

		if ( fileErr != (DIM1*DIM2*2*DIM3) ) {
			printf("ERROR: could not complete write to '%s'.\n", fileName);
			return(0);
		}

		// DON'T FORGET!!!
		fclose(filePtr);

	// 1: REAL FIELD

		// set output file name
		sprintf(fileName, "%s%s_filt_Real%s", "lognorm", fileInfix, fileExt);

		filePtr = NULL;
		if ( (filePtr = fopen(fileName, "w") ) == NULL) {
			printf("ERROR: could not open file '%s'.\n", fileName);
			return(0);
		}

		fileErr = 0;

		for (i = 1; i <= DIM1; i++)
			for (j = 1; j <= DIM2; j++)
				for (k = 1; k <= DIM3; k++)
					fileErr += fwrite(&logNormal[i][j][2*k-1], sizeof(Real), 1, filePtr);

		if ( fileErr != (DIM1*DIM2*DIM3) ) {
			printf("ERROR: could not complete write to '%s'.\n", fileName);
			return(0);
		}

		// DON'T FORGET!!!
		fclose(filePtr);

#endif

}
#endif // BINARYOUT
// ---------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------
// ASCII OUTPUT: POLYNOMIAL FIT TO LOGNORMAL LOG10(P(k)) vs. LOG10(k)
// #if (DEBUG)
{
// #ifdef SIM1D
// 
// #if NYQUIST
// 	i3 = (short) (0.5 * (Real) DIM3)+1; // + 1 to include Nyquist frequency
// #else
// 	i3 = (short) (0.5 * (Real) DIM3);
// #endif
// 	
// #else // 2D / 3D
// 
// 	#ifdef SIM2D
// 
// #if NYQUIST
// 		i3 = (short) (0.5 * (Real) DIM2)+1; // + 1 to include Nyquist frequency
// #else
// 	i3 = (short) (0.5 * (Real) DIM2);
// #endif
// 	
// 	#else
// 
// #if NYQUIST
// 	i3 = (short) (0.5 * (Real) DIM1)+1; // + 1 to include Nyquist frequency
// #else
// 	i3 = (short) (0.5 * (Real) DIM1);
// #endif
// 	
// 	#endif
// 
// #endif
// 
// 	for(k = 2; k <= i3; k++)
// 		printf("%d %e %e %e\n", k, 
// 		logWaveNumber[k], logPowSpecLogNorm[k], polyFitLogPowSpecLogNorm[k]);
// 
// 	printf("\n\n");
// 
// 	return(0);
}
// #endif // DEBUG
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
#if !(DEBUG)
{
	printf("# iter: %3d   target/current slope: %8e/%8e (+/- %8e)  norm: %8e\n",
	iterNum, BETA, slopeLogNorm,  fabs((BETA - slopeLogNorm)/BETA), normConstLogNorm);
#ifdef SIM1D
	printf("# target/actual mu (GAUSS) is %lf/%lf\n", muGauss, muGauss+average(noiseFiltered,DIM3));
	printf("# target/actual sigma (GAUSS) is %lf/%lf\n", sigmaGauss, sigmaGauss*stddev(noiseFiltered,DIM3));
	printf("# target/actual mu (LOGN) is %lf/%lf\n", muLogn,
	cmp_mu( average_log(logNormal,DIM3), stddev_log(logNormal,DIM3) ));
	printf("# target/actual sigma (LOGN) is %lf/%lf\n", sigmaLogn,
	cmp_sigma( average_log(logNormal,DIM3), stddev_log(logNormal,DIM3) ));
#else
	printf("# target/actual mu (GAUSS) is %lf/%lf\n", muGauss, muGauss+average(noiseFiltered,DIM1,DIM2,DIM3));
	printf("# target/actual sigma (GAUSS) is %lf/%lf\n", sigmaGauss, sigmaGauss*stddev(noiseFiltered,DIM1,DIM2,DIM3));
	printf("# target/actual mu (LOGN) is %lf/%lf\n", muLogn,
	cmp_mu( average_log(logNormal,DIM1,DIM2,DIM3), stddev_log(logNormal,DIM1,DIM2,DIM3) ));
	printf("# target/actual sigma (LOGN) is %lf/%lf\n", sigmaLogn,
	cmp_sigma( average_log(logNormal,DIM1,DIM2,DIM3), stddev_log(logNormal,DIM1,DIM2,DIM3) ));
#endif
}
#endif
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// COMPUTE 'SEARCH DIRECTION'
{
#ifdef SIM1D

#if NYQUIST
	ii3 = (short) (0.5 * (Real) DIM3)+1; // + 1 to include Nyquist frequency
#else
	ii3 = (short) (0.5 * (Real) DIM3);
#endif

#else // 2D / 3D

	#ifdef SIM2D

#if NYQUIST
		ii3 = (short) (0.5 * (Real) DIM2)+1; // + 1 to include Nyquist frequency
#else
		ii3 = (short) (0.5 * (Real) DIM2);
#endif

	#else

#if NYQUIST
		ii3 = (short) (0.5 * (Real) DIM1)+1; // + 1 to include Nyquist frequency
#else
		ii3 = (short) (0.5 * (Real) DIM1);
#endif

	#endif

#endif

	searchDir = vector(1,ii3);

	searchDirMean = 0.;

	for (i3 = 2; i3 <= ii3; i3++) {

		searchDir[i3] = 0.;

		x = (Real) (i3-1);

#if (POLYFIT)
		// compute value of polynomial fit to log-normal P(k) at log(x)
		polyFitLogPowSpecLogNorm[1] =  0.;
		poly(log10(x), polyVal, numParams);		
		for (m = 1; m <= numParams; m++)
			polyFitLogPowSpecLogNorm[1] +=  polyCoeff[m]*polyVal[m];

		// difference between log[target P(k)] and log[actual log-normal P(k)]
		searchDir[i3] = BETA * log10(x) + normConstLogNorm - polyFitLogPowSpecLogNorm[1];
#endif // POLYFIT

#if (SGFILTER)
		searchDir[i3] = BETA * log10(x) + normConstLogNorm - sgLogPowSpecLogNorm[i3];
#endif // SGFILTER

		searchDirMean += fabs(searchDir[i3] / (BETA * log10(x) + normConstLogNorm));

// #if DEBUG
// 		printf("%e %e %e %e %e\n",
// 		log10(x), polyFitLogPowSpecLogNorm[1],
// 		BETA * log10(x) + normConstLogNorm,
// 		BETA * log10(x) + normConstLogNorm + searchDir[i3],
// 		searchDir[i3]);
// #endif

	} // i3-loop

	// average difference
	searchDirMean /= (Real) ii3;

#if DEBUG
	printf("\n\n");
#endif

}
// ---------------------------------------------------------------------------------------

// #if (DEBUG)
	printf("# Search direction (mean): %e\n\n", searchDirMean);
// #endif

// ---------------------------------------------------------------------------------------
// BINARY OUTPUT: SEARCH DIRECTION VECTOR
#if (DEBUG)
#if BINARYOUT
{	// set output file name
	sprintf(fileName, "%s%s_%d%s", "searchDir", fileInfix, iterNum, fileExt);

	if ( (filePtr = fopen(fileName, "w")) == NULL) {
		printf("ERROR: could not open file '%s'.\n", fileName);
		return(0);
	}

	fileErr = fwrite(&searchDir[1], sizeof(Real), ii3, filePtr);

	if ( fileErr != ii3) {
		printf("ERROR: could not complete write to '%s'.\n", fileName);
		return(0);
	}

}
#endif // BINARYOUT
#endif // DEBUG
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// UPDATE GUESS FOR P(k) USING "SEARCH DIRECTION" [EFFECTIVELY, DIFFERENCE BETWEEN ACTUAL AND TARGET P(k)]
// IMPORTANT: spectrum is folded at 1 + DIM/2 to avoid aliasing effects
{
#ifdef SIM1D
{

	ii3 = (short) (0.5 * (Real) DIM3);

	stepSize = 0.5; // trial-and-error value; smaller / larger values affect convergence

	for (i3 = 2; i3 <= ii3+1; i3++) { // up to Nyquist; skip wavenumber k = 0

		// set step direction
		searchStep = searchDir[i3];
		// scale step
		searchStep *= stepSize;
		// "add" logarithmic (hence 10^x), square root (0.5) step to previous guess
		specDens[i3] *= pow(1.0e1, 0.5*searchStep);

		if ((i3 > 1) && (i3 <= ii3)) // fold wrt k
			specDens[DIM3+2-i3] = specDens[i3];

	} // i3-loop

	normFactSpecDensFin = 0.;

	for (k = 1; k <= DIM3; k++)
		normFactSpecDensFin += (double)((double)specDens[k]*(double)specDens[k]);

  	normFactSpecDensFin = sqrt((double) (DIM3) / normFactSpecDensFin);

	// normalisation to preserve variance
	for (k = 1; k <= DIM3; k++)
		specDens[k] *= (Real) normFactSpecDensFin;

	// reset value at k = 0
	specDens[1] = 1.0; // condition to preserve mean value of field



}
#else // 2D / 3D

#ifdef SIM2D // 2D
{

#if NYQUIST
	m = (short) (0.5 * (Real) DIM2)+1; // Nyquist
#else
	m = (short) (0.5 * (Real) DIM2);
#endif

	ii2 = (short) (0.5 * (Real) DIM2);
	ii3 = (short) (0.5 * (Real) DIM3);

	stepSize = 0.5; // trial-and-error value; smaller / larger values affect convergence

	for (i2 = 1; i2 <= ii2+1; i2++) { // up to Nyquist

		j = i2;

		y = (Real) (j-1);

		for(i3 = 1; i3 <= ii3+1; i3++) { // up to Nyquist

			k = i3; // wrap-around

			z = (Real) (k-1);

			if ( (i2 > 1) || (i3 > 1) ) { // avoid w = 0

				w = (Real) sqrt(y*y + z*z);

				// update guess for spectral density

				// set step direction
				searchStep =
				( (1+(int) floor(w) <= m) ? searchDir[1+(int) floor(w)] : searchDir[m]);

				// scale step
				searchStep *= stepSize;

				// "add" logarithmic (hence 10^x), square root (0.5) step to previous guess
				specDens[1][i2][i3] *= pow(1.0e1, 0.5*searchStep);

				if ((i2 > 1) && (i2 <= ii2)) { // fold 1 wrt j at frequency below Nyquist

					specDens[1][DIM2+2-i2][i3] = specDens[1][i2][i3];

					if ((i3 > 1) && (i3 <= ii3)) // fold 1 wrt j and k
						specDens[1][DIM2+2-i2][DIM3+2-i3] = specDens[1][i2][i3];

				} // 1 < j < ii2

				if ((i3 > 1) && (i3 <= ii3)) // fold 1 wrt k
						specDens[1][i2][DIM3+2-i3] = specDens[1][i2][i3];

			} // avoid w = 0

		} // i3-loop

	} // i2-loop

	normFactSpecDensFin = 0.;

	specDens[1][1][1] = 0.0;

	// compute normalisation constant to preserve variance
	// taking into account cut-off at KMIN
		for (j = 1; j <= DIM2; j++)
			for (k = 1; k <= DIM3; k++)
				normFactSpecDensFin += (double)((double)specDens[1][j][k] * (double)specDens[1][j][k]);

	normFactSpecDensFin = sqrt((double) (DIM2*DIM3) / normFactSpecDensFin);

		for (j = 1; j <= DIM2; j++)
			for (k = 1; k <= DIM3; k++)
				specDens[1][j][k] *= (Real) normFactSpecDensFin;


	// reset value at k = 0
	specDens[1][1][1] = 1.0;	// condition to preserve mean value of field
}
#else // 3D
{

#if NYQUIST
	m = (short) (0.5 * (Real) DIM1)+1; // Nyquist
#else
	m = (short) (0.5 * (Real) DIM1);
#endif

	ii1 = (short) (0.5 * (Real) DIM1);
	ii2 = (short) (0.5 * (Real) DIM2);
	ii3 = (short) (0.5 * (Real) DIM3);

	stepSize = 0.5; // trial-and-error value; smaller / larger values affect convergence

	for (i1 = 1; i1 <= ii1+1; i1++) { // up to Nyquist frequency

		i = i1;
		x = (Real) (i-1);

		for (i2 = 1; i2 <= ii2+1; i2++) { // up to Nyquist frequency

			j = i2;
			y = (Real) (j-1);

			for (i3 = 1; i3 <= ii3+1; i3++) { // up to Nyquist frequency

				k = i3;
				z = (Real) (k-1);

				if ( (i1 > 1) || (i2 > 1) || (i3 > 1) ) { // avoid w = 0


					w = (Real) sqrt(x*x + y*y + z*z);


					// update guess for spectral density

					// set step direction
					searchStep =
					( (1+(int) floor(w) <= m) ? searchDir[1+(int) floor(w)] : searchDir[m]);

					// scale step
					searchStep *= stepSize;

					// "add" logarithmic (hence 10^x), square root (0.5) step to previous guess
					specDens[i1][i2][i3] *= pow(1.0e1, 0.5 *searchStep);

					if ((i1 > 1) && (i1 <= ii1)) { // fold 1 wrt i

						specDens[DIM1+2-i1][i2][i3] = specDens[i1][i2][i3];

						if ((i2 > 1) && (i2 <= ii2)) {  // fold 1 wrt i and j

							specDens[DIM1+2-i1][DIM2+2-i2][i3] = specDens[i1][i2][i3];

							if ((i3 > 1) && (i3 <= ii3)) // fold 1 wrt to i, j, k
								specDens[DIM1+2-i1][DIM2+2-i2][DIM3+2-i3] =
								specDens[i1][i2][i3];

						}

						if ((i3 > 1) && (i3 <= ii3)) // fold 1 wrt i and k
							specDens[DIM1+2-i1][i2][DIM3+2-i3] = specDens[i1][i2][i3];

					}
					
					if ((i2 > 1) && (i2 <= ii2)) { // fold 1 wrt j

						specDens[i1][DIM2+2-i2][i3] = specDens[i1][i2][i3];

						if ((i3 > 1) && (i3 <= ii3)) // fold 1 wrt j and k
							specDens[i1][DIM2+2-i2][DIM3+2-i3] = specDens[i1][i2][i3];

					}

					if ((i3 > 1) && (i3 <= ii3)) // fold 1 wrt k
						specDens[i1][i2][DIM3+2-i3] = specDens[i1][i2][i3];

				} // avoid w = 0

			} // k-loop

		} // j-loop

	} // i-loop

	normFactSpecDensFin = 0.;

	for (i = 1; i <= DIM1; i++)
		for (j = 1; j <= DIM2; j++)
			for (k = 1; k <= DIM3; k++)
				normFactSpecDensFin += (double) ((double)specDens[i][j][k]*(double)specDens[i][j][k]);

  	normFactSpecDensFin = sqrt((double) (DIM1*DIM2*DIM3) / normFactSpecDensFin);

	for (i = 1; i <= DIM1; i++)
		for (j = 1; j <= DIM2; j++)
			for (k = 1; k <= DIM3; k++)
				specDens[i][j][k] *= (Real) normFactSpecDensFin;

	// reset value at k = 0
	specDens[1][1][1] = 1.0; // condition to preserve mean value of field

}
#endif // 2D or 3D?

#endif // 1D?

// 	printf("%e %e %e\n\n",
// 	normFactSpecDensFin, normFactSpecDensIni, normConstLogNorm);

}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// ASCII OUTPUT: UPDATED SPECTRAL DENSITY
// #if (DEBUG)
{
// #ifdef SIM1D
// 	
// 	for(k = 1; k <= DIM3; k++)
// 		printf("%d %e %e\n", k, specDens[k], normFactSpecDensFin);
// 
// 	return(0);
// 
// #else // 2D / 3D
// 
// 	for (i = 1; i <= DIM1; i++)
// 		for (j = 1; j <= DIM2; j++)
// 			for (k = 1; k <= DIM3; k++)
// 				printf("%d %d %d %e %e\n",
// 				i, j, k, specDens[i][j][k], normFactSpecDensFin);
// 
// 	printf("\n\n");
// 
// // 	return(0);
// 
// #endif
// #endif // DEBUG
}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// FREE MEMORY
{
#ifdef SIM1D
	#if NYQUIST
		i3 = (short) (0.5 * (Real) DIM3)+1; // + 1 to include Nyquist frequency
	#else
		i3 = (short) (0.5 * (Real) DIM3);
	#endif
#endif
#ifdef SIM2D
	#if NYQUIST
		i3 = (short) (0.5 * (Real) DIM2)+1; // + 1 to include Nyquist frequency
	#else
		i3 = (short) (0.5 * (Real) DIM2);
	#endif
#endif
#ifdef SIM3D
	#if NYQUIST
		i3 = (short) (0.5 * (Real) DIM1)+1; // + 1 to include Nyquist frequency
	#else
		i3 = (short) (0.5 * (Real) DIM1);
	#endif
#endif

#if (POLYFIT)
	free_vector(polyFitLogPowSpecLogNorm,1,i3);
	free_vector(polyVal,1,numParams);
	free_vector(polyCoeff,1,numParams);
	free_ivector(polyCoeffFlag,1,numParams);
	free_matrix(covar,1,numParams,1,numParams);
#endif
	free_vector(logPowSpecLogNorm,1,i3);
	free_vector(powSpecLogNorm,1,i3);

	free_vector(logWaveNumber,1,i3);
	free_vector(sigmaWeight,1,i3);


}
// ---------------------------------------------------------------------------------------

	printf("\n");

	} // END ITERATION

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
	printf("\n\n");
	return(0);

} // end of main
// ---------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------
// USER function declaration

Real cmp_m(Real mu, Real sig)
{//returns Gaussian expected value for a log-normal mean mu and sigma sig

	return 2.*log(mu) - 0.5 * log(sig*sig + mu*mu);

} // end of cmp_mu

Real cmp_s(Real mu, Real sig)
{//returns Gaussian std. deviation for a log-normal mean mu and sigma sig

	return sqrt(log(sig*sig + mu*mu) - 2.*log(mu));

} // end of cmp_mu

Real cmp_mu(Real m, Real s)
{//returns expected value for a Gaussian mean m and sigma s

	return exp(m + 0.5 * s*s);

} // end of cmp_mu

Real cmp_sigma(Real m, Real s)
{//returns log-normal std. deviation for a Gaussian mean m and sigma s
	Real mu;

	mu = cmp_mu(m,s);
	return mu * sqrt(exp(s*s) - 1.);

} // end of cmp_mu

#ifdef SIM1D

	Real average(Real *arr, int l)
	{
		int k;
		Real avg = 0.;
		for(k=1;k<=l;k++)
				avg += arr[2*k-1];
		avg /= (l);
		return avg;
	} // end of average

	Real stddev(Real *arr, int l)
	{
		int k;
		Real avg = 0., std = 0., val;
		avg = average(arr, l);
		for(k=1;k<=l;k++){
			val = arr[2*k-1];
			std += (avg - val) * (avg - val);
		}
		std = sqrt(std/(l - 1));
		return std;
	} // end of stddev

	Real average_log(Real *arr, int l)
	{
		int k;
		Real avg = 0.;
		for(k=1;k<=l;k++)
				avg += log(arr[2*k-1]);
		avg /= (l);
		return avg;
	} // end of average

	Real stddev_log(Real *arr, int l)
	{
		int k;
		Real avg = 0., std = 0., val;
		avg = average(arr, l);
		for(k=1;k<=l;k++){
			val = log(arr[2*k-1]);
			std += (avg - val) * (avg - val);
		}
		std = sqrt(std/(l - 1));
		return std;
	} // end of stddev

#else

	Real average_log(Real ***arr, int n, int m, int l)
	{
		int i,j,k;
		Real avg = 0.;
		for(i=1;i<=n;i++)
			for(j=1;j<=m;j++)
				for(k=1;k<=l;k++)
					avg += log(arr[i][j][2*k-1]);
		avg /= (n * m * l);
		return avg;
	} // end of average

	Real stddev_log(Real ***arr, int n, int m, int l)
	{
		int i,j,k;
		Real avg = 0., std = 0., val;
		avg = average_log(arr, n, m, l);
		for(i=1;i<=n;i++)
			for(j=1;j<=m;j++)
				for(k=1;k<=l;k++){
					val = log(arr[i][j][2*k-1]);
					std += (avg - val) * (avg - val);
				}
		std = sqrt(std/(n * m * l - 1));
		return std;
	} // end of stddev


	Real average(Real ***arr, int n, int m, int l)
	{
		int i,j,k;
		Real avg = 0.;
		for(i=1;i<=n;i++)
			for(j=1;j<=m;j++)
				for(k=1;k<=l;k++)
					avg += arr[i][j][2*k-1];
		avg /= (n * m * l);
		return avg;
	} // end of average

	Real stddev(Real ***arr, int n, int m, int l)
	{
		int i,j,k;
		Real avg = 0., std = 0., val;
		avg = average(arr, n, m, l);
		for(i=1;i<=n;i++)
			for(j=1;j<=m;j++)
				for(k=1;k<=l;k++){
					val = arr[i][j][2*k-1];
					std += (avg - val) * (avg - val);
				}
		std = sqrt(std/(n * m * l - 1));
		return std;
	} // end of stddev

#endif

Real min(Real x, Real y)
{
	return (x < y ? x : y);
} // end of min

Real max(Real x, Real y)
{
	return (x > y ? x : y);
} // end of mas

Real sed(Real k, Real s)
{ /*!
	\fn
	 Returns the value of a power-law k^s
	*/

	return (Real) pow(k, s);

} // end of sed

short wrap(short i, short n, short m)
{ /*!
	\fn
	Calculates the value of index i wrapped-around (to the right) n places on an array of size m
	*/

	short wi;
	wi = (i + n) % (m);
	return (wi != 0 ? wi : m);

// Alternative (more convoluted but interesting):
// 	wi = ( (wi = (i + n) % (m)) != 0 ? wi : m );
// 	return wi;

} // end of wrap

void poly(Real x, Real pn[], int n)
{ /*!
	\fn
	Computes and returns in order in array pn each of the n+1 terms of a polynomial of degree n evaluated at x
	*/

	short i;

	pn[1] = 1.;
	for (i=2; i<=n+1; i++)
		pn[i] = pn[i-1]*x;

} // end of poly

// ---------------------------------------------------------------------------------------
