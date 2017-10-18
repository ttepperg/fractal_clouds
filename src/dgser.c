#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7
#include "fractalCloudsDefs.h"

void dgser(Real *gamser, Real a, Real x, Real *gln)
{
	Real dgammln(Real xx);
	void nrerror(char error_text[]);
	int n;
	Real sum,del,ap;

	*gln=dgammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine dgser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine dgser");
		return;
	}
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software "!D. */
