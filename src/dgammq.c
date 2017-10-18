#include "fractalCloudsDefs.h"

Real dgammq(Real a, Real x)
{
	void dgcf(Real *gammcf, Real a, Real x, Real *gln);
	void dgser(Real *gamser, Real a, Real x, Real *gln);
	void nrerror(char error_text[]);
	Real gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine dgammq");
	if (x < (a+1.0)) {
		dgser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		dgcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software "!D. */
