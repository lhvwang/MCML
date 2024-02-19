/**************************************************************************
 *      Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992.
 *      Program name: CONV.
 *
 *	A program used to process the data from MCML -
 *	A Monte Carlo simulation of light transport in
 *	multi-layered turbid media in ANSI Standard C.
 ****
 *      Starting Date:  10/1991.
 *      Current Date:   07/1993.
 *
 *      Lihong Wang, Ph.D.
 *      Steven L. Jacques, Ph.D.
 *      Laser Biology Research Laboratory - 17
 *      M.D. Anderson Cancer Center
 *      University of Texas
 *      1515 Holcombe Blvd.
 *      Houston, Texas 77030
 *      USA
 *
 *      Liqiong Zheng, B.S.
 *      Department of Computer Science
 *      University of Houston
 *      Houston, Texas
 ****
 *	See comments in convmcml.h.
 ****/

#include "convmcml.h"

#define MAX(x, y) ((x)>(y) ? (x):(y))
#define MIN(x, y) ((x)<(y) ? (x):(y))

/*************************************************************************
 *      If IBMPC is 1, file extensions are limited to 3 letters at most.
 ****/
#define IBMPC 0 

/********************************* Stuctures *****************************/

/****
 *	Parameters to describe a photon beam.
 *  	Original:  infinitely narrow beam. This is default for the
 *		   beam from the MCML output.
 *	Flat:	   Flat beam with radius R.
 *	Gaussian:  Gaussian with 1/e2 radius R.
 *	Arbitrary: general beam described by points with interpolation.
 ****/
typedef struct {
  enum {
    original, flat, gaussian, arbitrary
  }           type;		/* beam type. */
  double      P;		/* total power. [J] */
  double      R;		/* radius. [cm]     */

  char        pro_fname[STRLEN];/* profile name of the arbitrary beam. */
  short	      np;		/* number of points describing profile.*/
  double     *r;		/* r for beam profile. [cm] */
  double     *pd;		/* power density for profile. [J/cm2] */
}           BeamStru;

/****
 *	Specify which quantities are available to be extracted.
 *	E.g., Td_r_a_t means Td_r at a specific a and t.
 ****/
typedef struct {		/* use bit field to save space. */
  int         Td_r_a_t:1;	/* Td_r@a@t. */
  int         Td_r_a:1;		/* Td_r@a*t. */
  int         Td_r_t:1;		/* Td_r*a@t. */
  int         Td_r:1;		/* Td_r*a*t. */

  int         Td_a_r_t:1;	/* Td_a@r@t. */
  int         Td_a_r:1;		/* Td_a@r*t. */
  int         Td_a_t:1;		/* Td_a*r@t. */
  int         Td_a:1;		/* Td_a*r*t. */

  int         Td_t_r_a:1;	/* Td_t@r@a. */
  int         Td_t_r:1;		/* Td_t@a*a. */
  int         Td_t_a:1;		/* Td_t*r@a. */
  int         Td_t:1;		/* Td_t*r*a. */

  int         Rd_r_a_t:1;	/* Rd_r@a@t. */
  int         Rd_r_a:1;		/* Rd_r@a*t. */
  int         Rd_r_t:1;		/* Rd_r*a@t. */
  int         Rd_r:1;		/* Rd_r*a*t. */

  int         Rd_a_r_t:1;	/* Rd_a@r@t. */
  int         Rd_a_r:1;		/* Rd_a@r*t. */
  int         Rd_a_t:1;		/* Rd_a*r@t. */
  int         Rd_a:1;		/* Rd_a*r*t. */

  int         Rd_t_r_a:1;	/* Rd_t@r@a. */
  int         Rd_t_r:1;		/* Rd_t@a*a. */
  int         Rd_t_a:1;		/* Rd_t*r@a. */
  int         Rd_t:1;		/* Rd_t*r*a. */

  int         A_rz_t:1;		/* A_rz@t. */
  int         A_rz:1;		/* A_rz*t. */
  int         A_z_t:1;		/* A_z*r@t. */
  int         A_z:1;		/* A_z*r*t. */
  int         A_t_r_z:1;	/* A_t@r@z. */
  int         A_t_z:1;		/* A_t*r@z. */
  int         A_t:1;		/* A_t*r*z. */
  int         A_l:1;		/* A_layer. */
}           ExtractStru;

/*************************************************************************
 *	Data structures for the binary tree used to store evaluations of
 *	part of the integrand.
 ****/
struct Node {
  float       x, y;
  struct Node *left, *right;
};

typedef struct Node *LINK;/* link is a pointer to a node. */
typedef LINK TREE;	/* tree is a link. */

/****
 *	A structure used by CONV.
 *	x may represent z, a, t.
 ****/
typedef struct {
  Boolean     datain;		/* mco data is availalbe. */
  char        fversion;		/* mco file version. */
  ExtractStru extract;		/* extractable quantites. */

  /* input these parameters to the convolution function. */
  BeamStru    beam;		/* incident beam of finite size. */
  float       eps;		/* relative error in convolution. */
  short       nr;		/* nr of InStru. */
  double      dr;		/* dr of InStru. */
  short       nrc;		/* number of r gridlines in convolution. */
  double      drc;		/* r grid size. */
  short       nxc;		/* number of x gridlines in convolution. */
  double      dxc;		/* x grid size. */
  double    **M_rx;		/* impulse response to be convolved. */
  double     *Mb_x;		/* ballistic. */

  /* the following parameters are used by the convolution function. */
  double      r;		/* current r at which doing integration. */
  double      rc;		/* current r for which doing convolution. */
  short       ix;		/* current ix. */
  TREE        tree;		/* Tree to store ITheta() & ExpBessI0(). */
}           ConvStru;

/**************************************************************************
 *	Routine prototypes for dynamic memory allocation and
 *	release of arrays and matrices.
 *	Modified from Numerical Recipes in C.
 ****/
double     *AllocArray1D(short, short, char);
void        FreeArray1D(double *, short, short);
double    **AllocArray2D(short, short, short, short, char);
void        FreeArray2D(double **, short, short, short, short);
double   ***AllocArray3D(short, short, short, short, short, short, char);
void        FreeArray3D(double ***, short, short, short, short, short, short);

