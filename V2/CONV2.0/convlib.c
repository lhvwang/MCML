/**************************************************************************
 *  	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992.
 *
 *	Some routines modified from Numerical Recipes in C,
 *	including error report, array or matrix declaration
 *	and releasing, integrations.
 *
 *	Some frequently used routines are also included here.
 ****/

#include "conv.h"

/**************************************************************************
 *	Report error message to stderr, then exit the program
 *	with signal 1.
 ****/
void
ErrorExit(char error_text[])
{
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

/**************************************************************************
 *	Allocate an array with index from nl to nh inclusive.
 *
 *	Original matrix and vector from Numerical Recipes in C
 *	don't initialize the elements to zero. This will
 *	be accomplished by the following functions.
 *
 *      exit the program when mode = 1. 
 ****/
double     *
AllocArray1D(short nl, short nh, char mode)
{
  double     *v;
  short       i;

  v = (double *) malloc((unsigned) (nh - nl + 1) * sizeof(double));

  if (v != NULL) { 
    v = v - nl;
    for (i = nl; i <= nh; i++)
      v[i] = 0.0;			/* init. */
  } else if (mode == 1) 
    ErrorExit("allocation failure in AllocArray1D()");

  return v;
}

/**************************************************************************
 *	Allocate a matrix with row index from nrl to nrh
 *	inclusive, and column index from ncl to nch
 *	inclusive.
 *
 *      exit the program when mode = 1. 
 ****/
double    **
AllocArray2D(short nrl, short nrh,
	     short ncl, short nch, char mode)
{
  short       i, j;
  double    **m;

  m = (double **) malloc((unsigned) (nrh - nrl + 1)
			 * sizeof(double *));

  if (m != NULL) { 
    m = m - nrl;
    for (i = nrl; i <= nrh; i++) 
      if ((m[i] = AllocArray1D(ncl, nch, mode)) == NULL) { 
	if (mode == 1)
	  ErrorExit("allocation failure in AllocArray2D()");
	else {
	  for (j = nrl; j < i; j ++)
	    FreeArray1D(m[j], ncl, nch);
	  free(m);
	  m = NULL;
	  break;
	}
      }

  } else if (mode == 1) 
    ErrorExit("allocation failure in AllocArray2D()");

  return m;
}

/**************************************************************************
 *	Allocate a 3D array with row index from nrl to nrh
 *	inclusive, column index from ncl to nch
 *	inclusive, and depth index from ndl to ndh
 *	inclusive.
 *
 *      exit the program when mode = 1. 
 ****/
double   ***
AllocArray3D(short nrl, short nrh,
	     short ncl, short nch,
	     short ndl, short ndh, char mode)
{
  short       i, j;
  double   ***m;

  m = (double ***) malloc((unsigned) (nrh - nrl + 1) * sizeof(double **));

  if (m != NULL) { 
    m = m - nrl;
    for (i = nrl; i <= nrh; i++)
      if ((m[i] = AllocArray2D(ncl, nch, ndl, ndh, mode)) == NULL) {
	if (mode == 1)
	  ErrorExit("allocation failure in AllocArray3D()");
	else {
	  for (j = nrl; j < i; j ++)
	    FreeArray2D(m[j], ncl, nch, ndl, ndh);
          free(m);
          m = NULL;
          break;
        }
      }
  } else if (mode == 1) 
    ErrorExit("allocation failure in AllocArray3D()");

  return m;
}

/**************************************************************************
 *	Release the memory.
 ****/
void
FreeArray1D(double *v, short nl, short nh)
{
  if (v != NULL) 
    free((char *) (v + nl));
}

/**************************************************************************
 *	Release the memory.
 ****/
void
FreeArray2D(double **m, short nrl, short nrh,
	    short ncl, short nch)
{
  short       i;

  if (m != NULL) {
    for (i = nrh; i >= nrl; i--)
      free((char *) (m[i] + ncl));
    free((char *) (m + nrl));
  }
}

/**************************************************************************
 *  Release the memory.
 ****/
void
FreeArray3D(double ***m, short nrl, short nrh,
	    short ncl, short nch, short ndl, short ndh)
{
  short       i;

  if (m != NULL) {
    for (i = nrh; i >= nrl; i--)
      FreeArray2D(m[i], ncl, nch, ndl, ndh);
    free((char *) (m + nrl));
  }
}

/**************************************************************************
 *	Trapezoidal integration.
 ****/

#define FUNC(x, y) ((*func)(x, y))
float
trapzd(
       float (*func) (float, ConvStru *),
       float a,
       float b,
       int n,
       ConvStru *Conv_Ptr)
{
  float       x, tnm, sum, del;
  static float s;
  static int  it;
  int         j;

  if (n == 1) {
    it = 1;
    return (s = 0.5 * (b - a) * (FUNC(a, Conv_Ptr) + FUNC(b, Conv_Ptr)));
  } else {
    tnm = it;
    del = (b - a) / tnm;
    x = a + 0.5 * del;
    for (sum = 0.0, j = 1; j <= it; j++, x += del)
      sum += FUNC(x, Conv_Ptr);
    it *= 2;
    s = 0.5 * (s + (b - a) * sum / tnm);
    return s;
  }
}
#undef FUNC

/**************************************************************************
 *  	Returns the integral of the function func from a to b.  EPS is the 
 *	relative error of the integration.  This function is based on
 *	W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.  Flannery, 
 *	"Numerical Recipes in C," Cambridge University Press, 2nd edition, 
 *	(1992).
 ****/
#define JMAX 20
float
qtrap(
      float (*func) (float, ConvStru *), 
      float a, 
      float b, 
      ConvStru *Conv_Ptr)
{
  int         j;
  float       s, s_old;

  s_old = -1.0e30;
  for (j = 1; j <= JMAX; j++) {
    s = trapzd(func, a, b, j, Conv_Ptr);
    if (fabs(s - s_old) <= Conv_Ptr->eps * fabs(s_old))
      break;
    s_old = s;
  }
  return (s);
}
#undef JMAX

/**************************************************************************
 *	Modified Bessel function exp(-x) I0(x), for x >=0.
 *	This function was modified from the original bessi0(). Instead of
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

/**************************************************************************
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

/**************************************************************************
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
