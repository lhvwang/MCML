/*************************************************************************
 *  	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992.
 *
 *	Specify photon beam profile, and convolute over the beam.
 ****/

#include "conv.h"

/*************************************************************************
 *	Use the center of each r grid element if 0, use the optimized
 *	point if 1.
 ****/
#define BESTR 1

/*************************************************************************
 *	The optimized point in the ith grid.
 ****/
#define RBEST(i) ((i+0.5) + 1/(12*(i+0.5)))*dr

/*************************************************************************
 *	Used by convolution over Gaussian beam.  Ignore the value
 *	beyond GAUSSLIMIT*radius.
 ****/
#define GAUSSLIMIT 4

/********************************Declarations****************************/
float       qtrap(float (*) (float, ConvStru *), float, float, ConvStru *);
double      BessI0(double);
short       GetShort(short, short);
float       GetFloat(float, float);
Boolean     CommentLineQ(char *);

/********************************Binary Tree*****************************/
LINK
FillNode(float x, float y)
{
  LINK        l;

  l = (LINK) malloc(sizeof(struct Node));
  if (l) {
    l->x = x;
    l->y = y;
    l->left = l->right = NULL;
  }
  return (l);
}

/*************************************************************************
 *	Assume the (x, y) is not in the tree.
 ****/
void
InsertNode(TREE * TreePtr, float x, float y)
{
  LINK        l1, l2;		/* l1 buffers l2. */

  l1 = NULL;
  l2 = *TreePtr;
  while (l2 != NULL) {
    l1 = l2;
    if (x < l2->x)
      l2 = l2->left;
    else
      l2 = l2->right;
  }

  if (l1 == NULL)		/* Empty tree. */
    *TreePtr = FillNode(x, y);
  else if (x < l1->x)
    l1->left = FillNode(x, y);
  else
    l1->right = FillNode(x, y);
}

LINK
SearchNode(TREE Tree, float x)
{
  LINK        l = Tree;
  Boolean     found = 0;

  while (l != NULL && !found)
    if (x < l->x)
      l = l->left;
    else if (x > l->x)
      l = l->right;
    else
      found = 1;

  return (l);
}

void
FreeTree(TREE Tree)
{
  if (Tree) {
    FreeTree(Tree->left);
    FreeTree(Tree->right);
    free(Tree);
  }
}

/********************************Convolution*****************************/
/*************************************************************************
 *	Specify the parameters for a flat beam.
 ****/
void
GetFlatBeam(BeamStru * Beam_Ptr)
{
  Beam_Ptr->type = flat;

  printf("Total energy of the flat beam [J]: ");
  Beam_Ptr->P = GetFloat(FLT_MIN, FLT_MAX);
  printf("Radius of the flat beam [cm]: ");
  Beam_Ptr->R = GetFloat(FLT_MIN, FLT_MAX);
  printf("Total power: %10.2g J, and radius: %10.2g cm.\n",
	 Beam_Ptr->P, Beam_Ptr->R);
}

/*************************************************************************
 *	Specify the parameters for a Gaussian beam.
 ****/
void
GetGaussianBeam(BeamStru * Beam_Ptr)
{
  Beam_Ptr->type = gaussian;

  printf("Total energy of the Gaussian beam [J]: ");
  Beam_Ptr->P = GetFloat(FLT_MIN, FLT_MAX);
  printf("1/e2 Radius of the Gaussian beam [cm]: ");
  Beam_Ptr->R = GetFloat(FLT_MIN, FLT_MAX);
  printf("Total power: %10.2g J, and radius: %10.2g cm.\n",
	 Beam_Ptr->P, Beam_Ptr->R);
}

/*************************************************************************
 *      Specify the parameters for a arbitrary beam.
 ****/
void
GetArbitraryBeam(BeamStru * Beam_Ptr)
{
  FILE       *file;
  char        Fname[STRLEN], buf[STRLEN];
  long        file_pos;
  short       i;
  Boolean     problem = 0;

  printf("Input filename of two-column beam profile (or . to quit): ");
  scanf("%s", Fname);
  if (strlen(Fname) == 1 && Fname[0] == '.')
    return;
  if (!(file = fopen(Fname, "r")))
    return;

  strcpy(Beam_Ptr->pro_fname, Fname);
  /* skip comment line. */
  do
    if (fgets(buf, 255, file) == NULL) {
      puts("No data found in the profile file.");
      return;
    }
  while (CommentLineQ(buf));

  /* find number of points. */
  file_pos = ftell(file);
  Beam_Ptr->np = 0;
  while (fgets(buf, 255, file) != NULL)
    if (CommentLineQ(buf)) {
      puts("No comment or space line allowed in the middle of data points");
      return;
    } else
      Beam_Ptr->np++;
  fseek(file, file_pos, SEEK_SET);

  /* allocate arrays. */
  if ((Beam_Ptr->r = AllocArray1D(0, Beam_Ptr->np - 1, 0)) == NULL
      || (Beam_Ptr->pd = AllocArray1D(0, Beam_Ptr->np - 1, 0)) == NULL) {
    puts("Not enough memory to get the arbitrary beam.");
    return;
  }

  /* read the points. */
  for (i = 0; i < Beam_Ptr->np; i++) {
    fscanf(file, "%lf%lf", &(Beam_Ptr->r[i]), &(Beam_Ptr->pd[i]));
    if (i >= 1 && Beam_Ptr->r[i] < Beam_Ptr->r[i - 1]) {
      puts("r coordinates not sorted in profile file.");
      problem = 1;
    } else if (Beam_Ptr->pd[i] < 0) {
      puts("Negative power density in profile file.");
      problem = 1;
    }
  }

  fclose(file);
  if (problem)
    Beam_Ptr->type = original;
  else
    Beam_Ptr->type = arbitrary;
}

/*************************************************************************
 *	Specify the beam profile.
 ****/
void
LaserBeam(BeamStru * Beam_Ptr)
{
  char        cmd_str[STRLEN];

  printf("Beam profile:f=flat, g=Gaussian, a=arbitrary, q=quit: ");
  do
    gets(cmd_str);
  while (!strlen(cmd_str));

  switch (toupper(cmd_str[0])) {
  case 'F':
    GetFlatBeam(Beam_Ptr);
    break;
  case 'G':
    GetGaussianBeam(Beam_Ptr);
    break;
  case 'A':
    GetArbitraryBeam(Beam_Ptr);
    break;
  case 'Q':
    break;
  default:;
    puts("Unsupported beam type.");
    break;
  }
}

/*************************************************************************
 *	Specify the relative convolution error. 0.001-0.1 recomended.
 ****/
void
ConvError(float *eps)
{
  printf("Relative convolution error\n");
  printf("Current value is %8.2g (0.001-0.1 recommended): ",
	 *eps);
  *eps = GetFloat(FLT_MIN, 1);
}

/*************************************************************************
 *      Specify the resolution and the number of points in r direction.
 ****/
void
ConvResolution(ConvStru * Conv_Ptr)
{
  char        string[STRLEN];

  printf("Current resolution: %10.2g cm and number of points: %4d\n",
	 Conv_Ptr->drc, Conv_Ptr->nrc);
  do {
    printf("Do you want to change them? (Y/N): ");
    do {
      gets(string);
    } while (!strlen(string));
  } while (toupper(string[0]) != 'Y' && toupper(string[0]) != 'N');

  if (toupper(string[0]) == 'Y') {
    printf("Input resolution in r direction [cm]: ");
    Conv_Ptr->drc = GetFloat(FLT_MIN, FLT_MAX);
    printf("Input number of points in r direction: ");
    Conv_Ptr->nrc = GetShort(1, SHRT_MAX);
    printf("New resolution: %10.2g cm and number of points: %4d\n",
	   Conv_Ptr->drc, Conv_Ptr->nrc);
  }
}

/********************************Integration*****************************/
/*************************************************************************
 *	Compute Itheta shown in the manual.
 ****/
double
ITheta(double r, double r2, double R)
{
  double      temp;

  if (R >= r + r2)
    return (1);
  else if (fabs(r - r2) <= R) {
    temp = (r * r + r2 * r2 - R * R) / (2 * r * r2);
    if (fabs(temp) > 1)
      temp = SIGN(temp);
    return (acos(temp) / PI);
  } else			/* R < fabs(r-r2) */
    return (0);
}

/*************************************************************************
 ****/
double
ExpBessI0(double r, double r2, double R)
{
  double      expbess;
  double      _RR = 1 / (R * R);
  double      x = 4 * r * r2 * _RR;
  double      y = 2 * (r2 * r2 + r * r) * _RR;

  expbess = exp(-y + x) * BessI0(x);
  return (expbess);
}

/*************************************************************************
 *	Interpolate for the arrays A_rz[][iz] with the optimized
 *	points in each grid element.
 *
 *	rb(i) = [(i+0.5) + 1/(12(i+0.5))] dr.
 ****/
double
M_rxInterp(double r2, ConvStru * Conv_Ptr)
{
  double      ir2_dbl, r2lo, M_lo, M_hi, M_at_r2, dr = Conv_Ptr->dr;
  short       ix = Conv_Ptr->ix, ir2lo, nr = Conv_Ptr->nr;

  if (nr <= 2)
    M_at_r2 = Conv_Ptr->M_rx[0][ix];
  else {
    /* find the low index according to r2. */
    ir2_dbl = r2 / dr;
    if (ir2_dbl <= 1)
      ir2lo = 0;
    else
      ir2lo = 0.5 * (ir2_dbl + sqrt(ir2_dbl * ir2_dbl - 1 / 3) - 0.5);
    ir2lo = MIN(ir2lo, nr - 3);

    r2lo = RBEST(ir2lo);
    M_lo = Conv_Ptr->M_rx[ir2lo][ix];
    M_hi = Conv_Ptr->M_rx[ir2lo + 1][ix];
    M_at_r2 = M_lo + (M_hi - M_lo) * (r2 - r2lo)
      / (RBEST(ir2lo + 1) - r2lo);
  }

  return (MAX(0, M_at_r2));
}

/*************************************************************************
 *
 ****/
double
ArbSource(double rho,		/* radius from center of beam. */
	  ConvStru * Conv_Ptr)
{
  BeamStru    beam;
  short       ip = 0;
  double      temp;

  beam = Conv_Ptr->beam;
  while (ip < beam.np)
    if (beam.r[ip++] > rho)
      break;

  if (beam.np == 1)		/* single point. */
    return ((ip < beam.np) ? beam.pd[0] : 0);
  else				/* interpolate. */
    temp = beam.pd[ip] + (rho - beam.r[ip]) *
      (beam.pd[ip - 1] - beam.pd[ip]) / (beam.r[ip - 1] - beam.r[ip]);
  return (MAX(0, temp));
}

/*************************************************************************
 *
 ****/
float
ArbIntegrand(float theta2, ConvStru * Conv_Ptr)
{
  double      r = Conv_Ptr->r, rc = Conv_Ptr->rc;

  return (ArbSource(sqrt(r * r + rc * rc - 2 * r * rc * cos(theta2)),
		    Conv_Ptr));
}

/*************************************************************************
 *	If the integrand at r2 has not been stored in the binary tree,
 *	evaluate it for either flat or gaussian beams, then store it in
 *	the tree.
 *	If the integarnd has bee stored in the binary tree, just fetch
 *	the result.
 *
 *	The integral variable r2 is the r" in the formula shown in the
 *	manual.  When r2 is in the range of recorded array, interpolation
 *	is used to evaluate the A_rz at r2.
 *
 *	Note that since the last grid elements collect all the
 *	photon weight that falls beyond the grid system, we should
 *	avoid using them in the convolution.
 ****/
float
M_rxIntegrand(float r2, ConvStru * Conv_Ptr)
{
  float       integrand;	/* part of the integrand. */
  LINK        link;
  BeamStru    beam;

  beam = Conv_Ptr->beam;
  if ((link = SearchNode(Conv_Ptr->tree, r2)))	/* integrand in tree. */
    integrand = link->y;
  else {
    if (beam.type == flat)
      integrand = ITheta(Conv_Ptr->rc, r2, beam.R);
    else if (beam.type == gaussian)
      integrand = ExpBessI0(Conv_Ptr->rc, r2, beam.R);
    else
      integrand = qtrap(ArbIntegrand, 0, 2 * PI, Conv_Ptr);
    InsertNode(&Conv_Ptr->tree, r2, integrand);
  }

  return (integrand * M_rxInterp(r2, Conv_Ptr) * r2);
}

/*************************************************************************
 *	Integrate Func for a flat beam, where Func is the integrand.
 ****/
double
FlatIntegration(
		float (*Func) (float, ConvStru *),
		ConvStru * Conv_Ptr)
{
  double      R = Conv_Ptr->beam.R;
  double      b_max = (Conv_Ptr->nr - 0.5) * Conv_Ptr->dr;
  double      a = MAX(0, Conv_Ptr->rc - R);
  double      b = MIN(b_max, Conv_Ptr->rc + R);

  if (a >= b)
    return (0);
  else
    return (qtrap(Func, a, b, Conv_Ptr));
}

/*************************************************************************
 *	Integrate Func for a Gaussian beam, where Func is the integrand.
 ****/
double
GaussIntegration(
		 float (*Func) (float, ConvStru *),
		 ConvStru * Conv_Ptr)
{
  double      R = Conv_Ptr->beam.R;
  double      b_max = (Conv_Ptr->nr - 0.5) * Conv_Ptr->dr;
  double      a = MAX(0, Conv_Ptr->rc - GAUSSLIMIT * R);
  double      b = MIN(b_max, Conv_Ptr->rc + GAUSSLIMIT * R);

  if (a >= b)
    return (0);
  else
    return (qtrap(Func, a, b, Conv_Ptr));
}

/*************************************************************************
 *	Integrate Func for an arbitrary beam, where Func is the integrand.
 ****/
double
ArbIntegration(
	       float (*Func) (float, ConvStru *),
	       ConvStru * Conv_Ptr)
{
  double      a = 0;
  double      b = (Conv_Ptr->nr - 0.5) * Conv_Ptr->dr;

  if (a >= b)
    return (0);
  else
    return (qtrap(Func, a, b, Conv_Ptr));
}

/*************************************************************************
 *	Convolve the matrix M_rx inside Conv_Ptr and assign the convolved
 *	matrix to M_rxc.
 *
 *	Parameters needed for the convolution are:
 *	  beam: the beam parameters.
 *	  eps:  convolution error.
 *	  Conv_Ptr->nr
 *	  Conv_Ptr->dr
 *	  Conv_Ptr->nrc
 *	  Conv_Ptr->drc
 *	  Conv_Ptr->nxc
 *	  Conv_Ptr->dxc
 *	  Conv_Ptr->M_rx: the impulse response vs r & x, where x can be
 *			  z, a, or t.
 *	  Conv_Ptr->Mb_x: ballistic component.
 *
 *	The estimate of "ballistic" is limited by 1 grid line separation.
 ****/
void
ConvM_rx(ConvStru * Conv_Ptr,
	 double **M_rxc /* convolved array. */ )
{
  short       irc, ix;
  double      P = Conv_Ptr->beam.P, R = Conv_Ptr->beam.R;
  double      rc, factor, ballistic;

  puts("The convolution may take a little while. Wait...");
  for (irc = 0; irc < Conv_Ptr->nrc; irc++) {
    Conv_Ptr->rc = rc = (irc + 0.5) * Conv_Ptr->drc;
    Conv_Ptr->tree = NULL;	/* init the tree. */

    if (Conv_Ptr->beam.type == flat)
      factor = (rc <= R) ? 0.5 / PI : 0;
    else if (Conv_Ptr->beam.type = gaussian)
      factor = 0.5 / PI * exp(-2 * rc * rc / (R * R));
    else
      factor = ArbSource(0, Conv_Ptr);

    /* convolution for different ix. */
    for (ix = 0; ix < Conv_Ptr->nxc; ix++) {
      Conv_Ptr->ix = ix;

      /* ballistic contribution. */
      if (Conv_Ptr->Mb_x != NULL)
	ballistic = Conv_Ptr->Mb_x[ix] * factor;
      else
	ballistic = 0;

      if (Conv_Ptr->beam.type == flat)
	M_rxc[irc][ix] = 2 * P / (R * R)
	  * (FlatIntegration(M_rxIntegrand, Conv_Ptr) + ballistic);
      else if (Conv_Ptr->beam.type == gaussian)
	M_rxc[irc][ix] = 4 * P / (R * R)
	  * (GaussIntegration(M_rxIntegrand, Conv_Ptr) + ballistic);
      else			/* Arbitrary. */
	M_rxc[irc][ix] = ArbIntegration(M_rxIntegrand, Conv_Ptr) + ballistic;
    }
    FreeTree(Conv_Ptr->tree);
  }
}
