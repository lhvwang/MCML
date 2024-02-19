/***********************************************************
 *	Some routines modified from Numerical Recipes in C,
 *	including error report, array or matrix declaration
 *	and releasing, integrations.
 *
 *	Some frequently used routines are also included here.
 ****/

/*
 * #include <stdlib.h> #include <stdio.h> #include <math.h>
 */
#include "conv.h"

/***********************************************************
 *	Report error message to stderr, then exit the program
 *	with signal 1.
 ****/
void
nrerror(char error_text[])
{
  fprintf(stderr, "%s.\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

/***********************************************************
 *	Allocate an array with index from nl to nh inclusive.
 *
 *	Original matrix and vector from Numerical Recipes in C
 *	don't initialize the elements to zero. This will
 *	be accomplished by the following functions.
 ****/
double     *
AllocVector(short nl, short nh)
{
  double     *v;
  short       i;

  v = (double *) malloc((unsigned) (nh - nl + 1) * sizeof(double));
  if (!v)
    nrerror("allocation failure in vector()");

  for (i = nl; i <= nh; i++)
    v[i] = 0.0;			/* init. */
  return v - nl;
}

/***********************************************************
 *	Allocate a matrix with row index from nrl to nrh
 *	inclusive, and column index from ncl to nch
 *	inclusive.
 ****/
double    **
AllocMatrix(short nrl, short nrh,
	    short ncl, short nch)
{
  short       i, j;
  double    **m;

  m = (double **) malloc((unsigned) (nrh - nrl + 1)
			 * sizeof(double *));
  if (!m)
    nrerror("allocation failure 1 in matrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (double *) malloc((unsigned) (nch - ncl + 1)
			     * sizeof(double));
    if (!m[i])
      nrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }

  for (i = nrl; i <= nrh; i++)
    for (j = ncl; j <= nch; j++)
      m[i][j] = 0.0;
  return m;
}

/***********************************************************
 *	Release the memory.
 ****/
void
FreeVector(double *v, short nl, short nh)
{
  free((char *) (v + nl));
}

/***********************************************************
 *	Release the memory.
 ****/
void
FreeMatrix(double **m, short nrl, short nrh,
	   short ncl, short nch)
{
  short       i;

  for (i = nrh; i >= nrl; i--)
    free((char *) (m[i] + ncl));
  free((char *) (m + nrl));
}

/***********************************************************
 *	Trapzoidal integration.
 ****/

#define FUNC(x) ((*func)(x))

float
trapzd(
       float (*func) (float),
       float a, float b, int n)
{
  float       x, tnm, sum, del;
  static float s;
  static int  it;
  int         j;

  if (n == 1) {
    it = 1;
    s = 0.5 * (b - a) * (FUNC(a) + FUNC(b));
  } else {
    tnm = it;
    del = (b - a) / tnm;
    x = a + 0.5 * del;
    for (sum = 0.0, j = 1; j <= it; j++, x += del)
      sum += FUNC(x);
    it *= 2;
    s = 0.5 * (s + (b - a) * sum / tnm);
  }
  return (s);
}
#undef FUNC

/***********************************************************
 *  Allow user to change EPS.
 *  Modifed to do at least three trapzd() in case data is
 *  noisy.  9/26/1994 Lihong Wang.
 ****/
#define JMAX 30
float
qtrap(
      float (*func) (float),
      float a,
      float b,
      float EPS)
{
  int         j;
  float       s, s_old = 0;

  for (j = 1; j <= JMAX; j++) {
    s = trapzd(func, a, b, j);
    if (j <= 3 || fabs(s - s_old) > EPS * fabs(s_old))
      s_old = s;
    else
      break;
  }
  return (s);
}
#undef JMAX

/***********************************************************
 *	Modified Bessel function exp(-x) I0(x), for x >=0.
 *	We modified from the original bessi0(). Instead of
 *	I0(x) itself, it returns I0(x) exp(-x).
 ****/
double
BessI0(double x)
{
  double      ax, ans;
  double      y;

  if ((ax = fabs(x)) < 3.75) {
    y = x / 3.75;
    y *= y;
    ans = exp(-ax) * (1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
	      + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2))))));
  } else {
    y = 3.75 / ax;
    ans = (1 / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1
		+ y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
	    + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
					       + y * 0.392377e-2))))))));
  }
  return ans;
}

/****************************************************************
 ****/
short
GetShort(short Lo, short Hi)
{
  char        in_str[STRLEN];
  short       x;

  gets(in_str);
  sscanf(in_str, "%hd", &x);
  while (x < Lo || x > Hi) {
    printf("...Wrong paramter.  Input again: ");
    gets(in_str);
    sscanf(in_str, "%hd", &x);
  }
  return (x);
}

/****************************************************************
 ****/
float
GetFloat(float Lo, float Hi)
{
  char        in_str[STRLEN];
  float       x;

  gets(in_str);
  sscanf(in_str, "%f", &x);
  while (x < Lo || x > Hi) {
    printf("...Wrong paramter.  Input again: ");
    gets(in_str);
    sscanf(in_str, "%f", &x);
  }
  return (x);
}
