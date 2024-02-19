/****************************************************************
 *	Specify photon beam profile, and convolute over the beam.
 *
 ****/

#include "conv.h"

/****************************************************************
 *	If IBMPC is 1, file extensions are limited to 3 letters
 *	at most.
 ****/
#define IBMPC 0

/****************************************************************
 *	Used by convolution over Gaussian beam.  Ignore the value
 *	beyond GAUSSLIMIT radius.
 ****/
#define GAUSSLIMIT 4


/***********************Declarations****************************/
FILE       *GetWriteFile(char *);
void        AllocConvData(InputStruct *, OutStruct *);
void        FreeConvData(InputStruct *, OutStruct *);
float       qtrap(float (*) (float), float, float, float);
double      BessI0(double);

/****************************************************************
 *	Data structures for the binary tree used to store part of
 *	the integrand evaluation.
 ****/
struct Node {
  float       x, y;
  struct Node *left, *right;
};

typedef struct Node *LINK;
typedef LINK TREE;

/****************************************************************
 *	A global structure to pass the current coordinate of the
 * 	physical quantities being evaluated and the pointers of the
 *	input and output parameters to the integration function.
 ****/
struct {
  double      r;
  short       iz, ia;
  InputStruct *in_ptr;
  OutStruct  *out_ptr;
  TREE        tree;		/* A tree to store ITheta() &*
				 * ExpBessI0(). */
}           ConvVar;


/***********************Convolution*****************************/
/****************************************************************
 *	Specify the parameters for a flat beam.
 ****/
void
GetFlatBeam(BeamStruct * Beam_Ptr)
{
  Beam_Ptr->type = flat;

  printf("Total energy of the flat beam [J]: ");
  Beam_Ptr->P = GetFloat(FLT_MIN, FLT_MAX);
  printf("Radius of the flat beam [cm]: ");
  Beam_Ptr->R = GetFloat(FLT_MIN, FLT_MAX);
  printf("Total power: %8.2lg J, and radius: %8.2lg cm.\n",
	 Beam_Ptr->P, Beam_Ptr->R);
}

/****************************************************************
 *	Specify the parameters for a Gaussian beam.
 ****/
void
GetGaussianBeam(BeamStruct * Beam_Ptr)
{
  Beam_Ptr->type = gaussian;

  printf("Total energy of the Gaussian beam [J]: ");
  Beam_Ptr->P = GetFloat(FLT_MIN, FLT_MAX);
  printf("1/e2 Radius of the Gaussian beam [cm]: ");
  Beam_Ptr->R = GetFloat(FLT_MIN, FLT_MAX);
  printf("Total power: %8.2lg J, and radius: %8.2lg cm.\n",
	 Beam_Ptr->P, Beam_Ptr->R);
}

/****************************************************************
 *	Specify the beam profile.
 ****/
void
LaserBeam(BeamStruct * Beam_Ptr, OutStruct * Out_Ptr)
{
  ConvStruct  null_conved = NULLCONVSTRUCT;
  char        cmd_str[STRLEN];

  printf("Beam profile:f=flat, g=Gaussian. q=quit: ");
  do
    gets(cmd_str);
  while (!strlen(cmd_str));

  switch (toupper(cmd_str[0])) {
  case 'F':
    Out_Ptr->conved = null_conved;
    GetFlatBeam(Beam_Ptr);
    break;
  case 'G':
    Out_Ptr->conved = null_conved;
    GetGaussianBeam(Beam_Ptr);
    break;
  default:;
  }
}


/****************************************************************
 *	Convolution.
 ***************************************************************/

/****************************************************************
 *	Specify the resolution and the number of points in r
 *	direction.
 *	Set the Out_Ptr->conved to null.
 *	Reallocate the arrays for the convolution arrays.
 ****/
void
ConvResolution(InputStruct * In_Ptr, OutStruct * Out_Ptr)
{
  char        in_str[STRLEN];
  ConvStruct  null_conved = NULLCONVSTRUCT;

  if (!Out_Ptr->allocated) {
    puts("...No mcml output data to work with");
    return;
  }
  printf("Current resolution: %8.2lg cm and number of points: %hd\n",
	 In_Ptr->drc, In_Ptr->nrc);
  printf("Input resolution in r direction [cm]: ");
  In_Ptr->drc = GetFloat(FLT_MIN, FLT_MAX);
  printf("Input number of points in r direction: ");
  In_Ptr->nrc = GetShort(1, SHRT_MAX);

  printf("Resolution: %8.2lg cm and number of points: %hd\n",
	 In_Ptr->drc, In_Ptr->nrc);

  Out_Ptr->conved = null_conved;
  FreeConvData(In_Ptr, Out_Ptr);
  AllocConvData(In_Ptr, Out_Ptr);
}

/****************************************************************
 *	Specify the relative convolution error. 0.001-0.1 recomended.
 *	Set the Out_Ptr->conved to null.
 ****/
void
ConvError(InputStruct * In_Ptr, OutStruct * Out_Ptr)
{
  char        in_str[STRLEN];
  ConvStruct  null_conved = NULLCONVSTRUCT;
  float       eps;

  printf("Relative convolution error\n");
  printf("Current value is %8.2g (0.001-0.1 recommended): ",
	 In_Ptr->eps);
  In_Ptr->eps = GetFloat(FLT_MIN, 1);
  Out_Ptr->conved = null_conved;
}

/***********************Binary Tree*****************************/
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

/****************************************************************
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

/***********************Integration*****************************/
/****************************************************************
 *	Compute Itheta shown in the manual.
 ****/
double
ITheta(double r, double r2, double R)
{
  double      temp;

  if (R >= r + r2)
    temp = 1;
  else if (fabs(r - r2) <= R) {
    temp = (r * r + r2 * r2 - R * R) / (2 * r * r2);
    if (fabs(temp) > 1)
      temp = SIGN(temp);
    temp = acos(temp) / PI;
  } else			/* R < fabs(r-r2) */
    temp = 0;

  return (temp);
}

/****************************************************************
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

/****************************************************************
 *	Interpolate for the arrays A_rz[].
 ****/
double
A_rzInterp(double **A_rz, double r2)
{
  double      ir2, A_lo, A_hi, A_at_r2;
  short       iz, ir2lo, nr = ConvVar.in_ptr->nr;

  ir2 = r2 / ConvVar.in_ptr->dr;
  iz = ConvVar.iz;
  if (nr < 3)
    A_at_r2 = A_rz[0][iz];
  else if (ir2 < nr - 1.5) {	/* interpolation. */
    ir2lo = MAX(0, (short) (ir2 - 0.5));	/* truncation. */
    A_lo = A_rz[ir2lo][iz];
    A_hi = A_rz[ir2lo + 1][iz];
    A_at_r2 = A_lo + (A_hi - A_lo) * (ir2 - ir2lo - 0.5);
  } else {			/* extrapolation. */
    ir2lo = nr - 3;
    A_lo = A_rz[ir2lo][iz];
    A_hi = A_rz[ir2lo + 1][iz];
    if (A_lo >= A_hi)		/* Noise test. */
      A_at_r2 = A_lo + (A_hi - A_lo) * (ir2 - ir2lo - 0.5);
    else
      A_at_r2 = 0.0;
  }

  return (MAX(0, A_at_r2));
}

/****************************************************************
 *	Interpolate for the arrays Rd_ra[] or Tt_ra[].
 ****/
double
RT_raInterp(double **RT_ra, double r2)
{
  double      ir2, RT_lo, RT_hi, RT_at_r2;
  short       ia, ir2lo, nr = ConvVar.in_ptr->nr;

  ir2 = r2 / ConvVar.in_ptr->dr;
  ia = ConvVar.ia;
  if (nr < 3)
    RT_at_r2 = RT_ra[0][ia];
  else if (ir2 < nr - 1.5) {	/* interpolation. */
    ir2lo = MAX(0, (short) (ir2 - 0.5));	/* truncation. */
    RT_lo = RT_ra[ir2lo][ia];
    RT_hi = RT_ra[ir2lo + 1][ia];
    RT_at_r2 = RT_lo + (RT_hi - RT_lo) * (ir2 - ir2lo - 0.5);
  } else {			/* extrapolation. */
    ir2lo = nr - 3;
    RT_lo = RT_ra[ir2lo][ia];
    RT_hi = RT_ra[ir2lo + 1][ia];
    if (RT_lo >= RT_hi)		/* Noise test. */
      RT_at_r2 = RT_lo + (RT_hi - RT_lo) * (ir2 - ir2lo - 0.5);
    else
      RT_at_r2 = 0.0;
  }

  return (MAX(0, RT_at_r2));
}

/****************************************************************
 *	Interpolate for the arrays Rd_r[] or Tt_r[].
 ****/
double
RT_rInterp(double *RT_r, double r2)
{
  double      ir2, RT_lo, RT_hi, RT_at_r2;
  short       ir2lo, nr = ConvVar.in_ptr->nr;

  ir2 = r2 / ConvVar.in_ptr->dr;
  if (nr < 3)
    RT_at_r2 = RT_r[0];
  else if (ir2 < nr - 1.5) {	/* interpolation. */
    ir2lo = MAX(0, (short) (ir2 - 0.5));	/* truncation. */
    RT_lo = RT_r[ir2lo];
    RT_hi = RT_r[ir2lo + 1];
    RT_at_r2 = RT_lo + (RT_hi - RT_lo) * (ir2 - ir2lo - 0.5);
  } else {			/* extrapolation. */
    ir2lo = nr - 3;
    RT_lo = RT_r[ir2lo];
    RT_hi = RT_r[ir2lo + 1];
    if (RT_lo >= RT_hi)		/* Noise test. */
      RT_at_r2 = RT_lo + (RT_hi - RT_lo) * (ir2 - ir2lo - 0.5);
    else
      RT_at_r2 = 0.0;
  }

  return (MAX(0, RT_at_r2));
}

/****************************************************************
 *	Convolution integrand for either flat or gaussian beams.
 *	Return the integrand for the convolution integral.
 *	r2 is the r" in the formula shown in the manual.
 *	When r2 is in the range of recorded array, interpolation
 *	is used to evaluate the diffuse reflectance at r2.
 *	Note that since the last grid elements collect all the
 *	photon weight that falls beyond the grid system, we should
 *	avoid using them in the convolution.
 ****/
float
A_rzFGIntegrand(float r2)
{				/* r" in the integration. */
  float       f;
  short       nr = ConvVar.in_ptr->nr;
  double      R, r, A_at_r2;
  LINK        link;

  A_at_r2 = A_rzInterp(ConvVar.out_ptr->A_rz, r2);

  R = ConvVar.in_ptr->beam.R;
  r = ConvVar.r;
  if ((link = SearchNode(ConvVar.tree, r2)))	/* f in tree. */
    f = link->y;
  else {
    if (ConvVar.in_ptr->beam.type == flat)
      f = ITheta(r, r2, R);
    else			/* Gaussian. */
      f = ExpBessI0(r, r2, R);
    InsertNode(&ConvVar.tree, r2, f);
  }

  f *= A_at_r2 * r2;
  return (f);
}

/****************************************************************
 *	Convolution integrand for either flat or gaussian beams.
 *	See comments for A_rzFGIntegrand().
 ****/
float
Rd_raFGIntegrand(float r2)
{				/* r" in the integration. */
  float       f;
  short       nr = ConvVar.in_ptr->nr;
  double      R, r, Rd_at_r2;
  LINK        link;

  Rd_at_r2 = RT_raInterp(ConvVar.out_ptr->Rd_ra, r2);

  R = ConvVar.in_ptr->beam.R;
  r = ConvVar.r;
  if ((link = SearchNode(ConvVar.tree, r2)))	/* f in tree. */
    f = link->y;
  else {
    if (ConvVar.in_ptr->beam.type == flat)
      f = ITheta(r, r2, R);
    else			/* Gaussian. */
      f = ExpBessI0(r, r2, R);
    InsertNode(&ConvVar.tree, r2, f);
  }

  f *= Rd_at_r2 * r2;
  return (f);
}

/****************************************************************
 *	Convolution integrand for either flat or gaussian beams.
 *	See comments for A_rzFGIntegrand().
 ****/
float
Rd_rFGIntegrand(float r2)
{				/* r" in the integration. */
  float       f;
  short       nr = ConvVar.in_ptr->nr;
  double      R, r, Rd_at_r2;

  Rd_at_r2 = RT_rInterp(ConvVar.out_ptr->Rd_r, r2);

  R = ConvVar.in_ptr->beam.R;
  r = ConvVar.r;
  if (ConvVar.in_ptr->beam.type == flat)
    f = Rd_at_r2 * ITheta(r, r2, R) * r2;
  else				/* Gaussian. */
    f = Rd_at_r2 * ExpBessI0(r, r2, R) * r2;

  return (f);
}

/****************************************************************
 *	Convolution integrand for either flat or gaussian beams.
 *	See comments for A_rzFGIntegrand().
 ****/
float
Tt_raFGIntegrand(float r2)
{				/* r" in the integration. */
  float       f;
  short       nr = ConvVar.in_ptr->nr;
  double      R, r, Tt_at_r2;
  LINK        link;

  Tt_at_r2 = RT_raInterp(ConvVar.out_ptr->Tt_ra, r2);

  R = ConvVar.in_ptr->beam.R;
  r = ConvVar.r;
  if ((link = SearchNode(ConvVar.tree, r2)))	/* f in tree. */
    f = link->y;
  else {
    if (ConvVar.in_ptr->beam.type == flat)
      f = ITheta(r, r2, R);
    else			/* Gaussian. */
      f = ExpBessI0(r, r2, R);
    InsertNode(&ConvVar.tree, r2, f);
  }

  f *= Tt_at_r2 * r2;
  return (f);
}

/****************************************************************
 *	Convolution integrand for either flat or gaussian beams.
 *	See comments for A_rzFGIntegrand().
 ****/
float
Tt_rFGIntegrand(float r2)
{				/* r" in the integration. */
  float       f;
  short       nr = ConvVar.in_ptr->nr;
  double      R, r, Tt_at_r2;

  Tt_at_r2 = RT_rInterp(ConvVar.out_ptr->Tt_r, r2);

  R = ConvVar.in_ptr->beam.R;
  r = ConvVar.r;
  if (ConvVar.in_ptr->beam.type == flat)
    f = Tt_at_r2 * ITheta(r, r2, R) * r2;
  else				/* Gaussian. */
    f = Tt_at_r2 * ExpBessI0(r, r2, R) * r2;

  return (f);
}

/****************************************************************
 ****/
double
FlatIntegration(float (*Func) (float))
{
  double      rc = ConvVar.r;
  double      R = ConvVar.in_ptr->beam.R;
  double      b_max = (ConvVar.in_ptr->nr - 0.5) * ConvVar.in_ptr->dr;
  double      a = MAX(0, rc - R);
  double      b = MIN(b_max, rc + R);

  if (a >= b)
    return (0);
  else
    return (qtrap(Func, a, b, ConvVar.in_ptr->eps));
}

/****************************************************************
 ****/
double
GaussIntegration(float (*Func) (float))
{
  double      rc = ConvVar.r;
  double      R = ConvVar.in_ptr->beam.R;
  double      b_max = (ConvVar.in_ptr->nr - 0.5) * ConvVar.in_ptr->dr;
  double      a = MAX(0, rc - GAUSSLIMIT * R);
  double      b = MIN(b_max, rc + GAUSSLIMIT * R);

  if (a >= b)
    return (0);
  else
    return (qtrap(Func, a, b, ConvVar.in_ptr->eps));
}

/****************************************************************
 ****/
void
ConvA_rz(InputStruct * In_Ptr,
	 OutStruct * Out_Ptr)
{
  short       irc, iz;
  double      rc, P = In_Ptr->beam.P, R = In_Ptr->beam.R;

  puts("The convolution may take a little while. Wait...");
  for (irc = 0; irc < In_Ptr->nrc; irc++) {
    rc = (irc + 0.5) * In_Ptr->drc;
    ConvVar.r = rc;
    ConvVar.tree = NULL;	/* init the tree. */
    for (iz = 0; iz < In_Ptr->nz; iz++) {
      ConvVar.iz = iz;
      if (In_Ptr->beam.type == flat)
	Out_Ptr->A_rzc[irc][iz] = 2 * P / (R * R)
	  * FlatIntegration(A_rzFGIntegrand);
      else			/* Gaussian. */
	Out_Ptr->A_rzc[irc][iz] = 4 * P / (R * R)
	  * GaussIntegration(A_rzFGIntegrand);
    }
    FreeTree(ConvVar.tree);
  }

  Out_Ptr->conved.A_rz = 1;
}

/****************************************************************
 ****/
void
ConvRd_ra(InputStruct * In_Ptr,
	  OutStruct * Out_Ptr)
{
  short       irc, ia;
  double      rc, P = In_Ptr->beam.P, R = In_Ptr->beam.R;
  double      b_max = (In_Ptr->nr - 1) * In_Ptr->dr;

  puts("The convolution may take a little while. Wait...");
  for (irc = 0; irc < In_Ptr->nrc; irc++) {
    rc = (irc + 0.5) * In_Ptr->drc;
    ConvVar.r = rc;
    ConvVar.tree = NULL;	/* init the tree. */
    for (ia = 0; ia < In_Ptr->na; ia++) {
      ConvVar.ia = ia;
      if (In_Ptr->beam.type == flat)
	Out_Ptr->Rd_rac[irc][ia] = 2 * P / (R * R)
	  * FlatIntegration(Rd_raFGIntegrand);
      else			/* Gaussian. */
	Out_Ptr->Rd_rac[irc][ia] = 4 * P / (R * R)
	  * GaussIntegration(Rd_raFGIntegrand);
    }
    FreeTree(ConvVar.tree);
  }

  Out_Ptr->conved.Rd_ra = 1;
}

/****************************************************************
 ****/
void
ConvRd_r(InputStruct * In_Ptr, OutStruct * Out_Ptr)
{
  short       irc;
  double      rc, P = In_Ptr->beam.P, R = In_Ptr->beam.R;
  double      b_max = (In_Ptr->nr - 1) * In_Ptr->dr;

  for (irc = 0; irc < In_Ptr->nrc; irc++) {
    rc = (irc + 0.5) * In_Ptr->drc;
    ConvVar.r = rc;
    if (In_Ptr->beam.type == flat)
      Out_Ptr->Rd_rc[irc] = 2 * P / (R * R)
	* FlatIntegration(Rd_rFGIntegrand);
    else			/* Gaussian. */
      Out_Ptr->Rd_rc[irc] = 4 * P / (R * R)
	* GaussIntegration(Rd_rFGIntegrand);
  }

  Out_Ptr->conved.Rd_r = 1;
}

/****************************************************************
 ****/
void
ConvTt_ra(InputStruct * In_Ptr, OutStruct * Out_Ptr)
{
  short       irc, ia;
  double      rc, P = In_Ptr->beam.P, R = In_Ptr->beam.R;
  double      b_max = (In_Ptr->nr - 1) * In_Ptr->dr;

  puts("The convolution may take a little while. Wait...");
  for (irc = 0; irc < In_Ptr->nrc; irc++) {
    rc = (irc + 0.5) * In_Ptr->drc;
    ConvVar.r = rc;
    ConvVar.tree = NULL;	/* init the tree. */
    for (ia = 0; ia < In_Ptr->na; ia++) {
      ConvVar.ia = ia;
      if (In_Ptr->beam.type == flat)
	Out_Ptr->Tt_rac[irc][ia] = 2 * P / (R * R)
	  * FlatIntegration(Tt_raFGIntegrand);
      else			/* Gaussian. */
	Out_Ptr->Tt_rac[irc][ia] = 4 * P / (R * R)
	  * GaussIntegration(Tt_raFGIntegrand);
    }
    FreeTree(ConvVar.tree);
  }

  Out_Ptr->conved.Tt_ra = 1;
}

/****************************************************************
 ****/
void
ConvTt_r(InputStruct * In_Ptr, OutStruct * Out_Ptr)
{
  short       irc;
  double      rc, P = In_Ptr->beam.P, R = In_Ptr->beam.R;
  double      b_max = (In_Ptr->nr - 1) * In_Ptr->dr;

  for (irc = 0; irc < In_Ptr->nrc; irc++) {
    rc = (irc + 0.5) * In_Ptr->drc;
    ConvVar.r = rc;
    if (In_Ptr->beam.type == flat)
      Out_Ptr->Tt_rc[irc] = 2 * P / (R * R)
	* FlatIntegration(Tt_rFGIntegrand);
    else			/* Gaussian. */
      Out_Ptr->Tt_rc[irc] = 4 * P / (R * R)
	* GaussIntegration(Tt_rFGIntegrand);
  }

  Out_Ptr->conved.Tt_r = 1;
}

/****************************************************************
 ****/
void
ShowOutConvMenu(char *in_fname)
{
  printf("Arz = absorption vs r & z [J/cm3]\n");
  printf("Frz = fluence vs r & z [J/cm2]\n");

  printf("Rr  = diffuse reflectance vs radius r [J/cm2]\n");
  printf("Rra = diffuse reflectance vs radius and angle [J/(cm2 sr)]\n");

  printf("Tr  = transmittance vs radius r [J/cm2]\n");
  printf("Tra = transmittance vs radius and angle [J/(cm2 sr)]\n");

  printf("Q   = Quit to main menu\n");
  printf("* input filename: %s \n", in_fname);
}


/****************************************************************
 *	3 numbers each line: r, z, A[r][z].
 ****/
void
WriteA_rzc(InputStruct * In_Ptr,
	   double **A_rzc)
{
  FILE       *file;
  short       ir, iz, nr = In_Ptr->nr, nz = In_Ptr->nz;
  double      r, z, dr = In_Ptr->dr, dz = In_Ptr->dz;
  char        fname[STRLEN];

#if IBMPC
  strcpy(fname, "Arz");
#else
  strcpy(fname, "Arzc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "%-12s\t%-12s\t%-s[J/cm3]\n",
	  "r[cm]", "z[cm]", fname);
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    for (iz = 0; iz < nz; iz++) {
      z = (iz + 0.5) * dz;
      fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n",
	      r, z, A_rzc[ir][iz]);
    }
  }

  fclose(file);
}

/****************************************************************
 *	3 numbers each line: r, z, F[r][z].
 ****/
void
WriteF_rzc(InputStruct * In_Ptr,
	   double **A_rzc)
{
  FILE       *file;
  short       ir, iz, nr = In_Ptr->nr, nz = In_Ptr->nz;
  double      mua, r, z, dr = In_Ptr->dr, dz = In_Ptr->dz;
  char        fname[STRLEN];

#if IBMPC
  strcpy(fname, "Frz");
#else
  strcpy(fname, "Frzc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "%-12s\t%-12s\t%-s[J/cm3]\n",
	  "r[cm]", "z[cm]", fname);
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    for (iz = 0; iz < nz; iz++) {
      z = (iz + 0.5) * dz;
      mua = In_Ptr->layerspecs[IzToLayer(iz, In_Ptr)].mua;
      if (mua > 0.0)
	fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n",
		r, z, A_rzc[ir][iz] / mua);
      else
	fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n", r, z, 0.0);
    }
  }

  fclose(file);
}

/****************************************************************
 *	3 numbers each line: r, a, Rd[r][a].
 ****/
void
WriteRd_rac(InputStruct * In_Ptr,
	    double **Rd_rac)
{
  short       ir, ia, nr = In_Ptr->nr, na = In_Ptr->na;
  double      r, a, dr = In_Ptr->dr, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

#if IBMPC
  strcpy(fname, "Rra");
#else
  strcpy(fname, "Rrac");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "%-12s\t%-12s\t%-s[J/(cm2sr)]\n",
	  "r[cm]", "a[rad]", fname);
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    for (ia = 0; ia < na; ia++) {
      a = (ia + 0.5) * da;
      fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n",
	      r, a, Rd_rac[ir][ia]);
    }
  }

  fclose(file);
}

/****
 *	2 numbers each line: r, Rd[r]
 ****/
void
WriteRd_rc(InputStruct * In_Ptr,
	   double *Rd_rc)
{
  short       ir, nr = In_Ptr->nr;
  double      dr = In_Ptr->dr;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Rrc");
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "%-12s\t%-s[J/cm2]\n", "r[cm]", fname);
  for (ir = 0; ir < nr; ir++)
    fprintf(file, "%-12.4E\t%-12.4E\n", (ir + 0.5) * dr, Rd_rc[ir]);

  fclose(file);
}

/****************************************************************
 *	3 numbers each line:r, a, Tt[r][a]. a = theta.
 ****/
void
WriteTt_rac(InputStruct * In_Ptr,
	    double **Tt_rac)
{
  short       ir, ia, nr = In_Ptr->nr, na = In_Ptr->na;
  double      r, a, dr = In_Ptr->dr, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

#if IBMPC
  strcpy(fname, "Tra");
#else
  strcpy(fname, "Trac");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "%-12s\t%-12s\t%-s[J/(cm2sr)]\n",
	  "r[cm]", "a[rad]", fname);
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    for (ia = 0; ia < na; ia++) {
      a = (ia + 0.5) * da;
      fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n",
	      r, a, Tt_rac[ir][ia]);
    }
  }

  fclose(file);
}

/****
 *	2 numbers each line: r, Tt[r].
 ****/
void
WriteTt_rc(InputStruct * In_Ptr,
	   double *Tt_rc)
{
  short       ir, nr = In_Ptr->nr;
  double      dr = In_Ptr->dr;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Trc");
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "%-12s\t%-s[J/cm2]\n", "r[cm]", fname);
  for (ir = 0; ir < nr; ir++)
    fprintf(file, "%-12.4E\t%-12.4E\n", (ir + 0.5) * dr, Tt_rc[ir]);

  fclose(file);
}

/****************************************************************
 ****/
void
BranchOutConvA(char *Cmd_Str,
	       InputStruct * In_Ptr,
	       OutStruct * Out_Ptr)
{
  switch (toupper(Cmd_Str[1])) {
    case 'R':
    if (toupper(Cmd_Str[2]) == 'Z') {	/* A_rzc. */
      if (!Out_Ptr->conved.A_rz)
	ConvA_rz(In_Ptr, Out_Ptr);
      WriteA_rzc(In_Ptr, Out_Ptr->A_rzc);
    } else
      puts("...Wrong command");
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void
BranchOutConvF(char *Cmd_Str,
	       InputStruct * In_Ptr,
	       OutStruct * Out_Ptr)
{
  switch (toupper(Cmd_Str[1])) {
    case 'R':
    if (toupper(Cmd_Str[2]) == 'Z') {	/* F_rzc. */
      if (!Out_Ptr->conved.A_rz)
	ConvA_rz(In_Ptr, Out_Ptr);
      WriteF_rzc(In_Ptr, Out_Ptr->A_rzc);
    } else
      puts("...Wrong command");
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void
BranchOutConvR(char *Cmd_Str,
	       InputStruct * In_Ptr,
	       OutStruct * Out_Ptr)
{
  char        ch;

  switch (toupper(Cmd_Str[1])) {
  case 'R':
    ch = toupper(Cmd_Str[2]);
    if (ch == '\0') {		/* Rd_rc. */
      if (!Out_Ptr->conved.Rd_r)
	ConvRd_r(In_Ptr, Out_Ptr);
      WriteRd_rc(In_Ptr, Out_Ptr->Rd_rc);
    } else if (ch == 'A') {	/* Rd_rac. */
      if (!Out_Ptr->conved.Rd_ra)
	ConvRd_ra(In_Ptr, Out_Ptr);
      WriteRd_rac(In_Ptr, Out_Ptr->Rd_rac);
    } else
      puts("...Wrong command");
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void
BranchOutConvT(char *Cmd_Str,
	       InputStruct * In_Ptr,
	       OutStruct * Out_Ptr)
{
  char        ch;

  switch (toupper(Cmd_Str[1])) {
  case 'R':
    ch = toupper(Cmd_Str[2]);
    if (ch == '\0') {		/* Tt_rc. */
      if (!Out_Ptr->conved.Tt_r)
	ConvTt_r(In_Ptr, Out_Ptr);
      WriteTt_rc(In_Ptr, Out_Ptr->Tt_rc);
    } else if (ch == 'A') {	/* Tt_rac. */
      if (!Out_Ptr->conved.Tt_ra)
	ConvTt_ra(In_Ptr, Out_Ptr);
      WriteTt_rac(In_Ptr, Out_Ptr->Tt_rac);
    } else
      puts("...Wrong command");
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 *	The global ConvVar is used.
 ****/
void
BranchOutConvCmd(char *Cmd_Str,
		 InputStruct * In_Ptr,
		 OutStruct * Out_Ptr)
{
  char        ch;

  ConvVar.in_ptr = In_Ptr;
  ConvVar.out_ptr = Out_Ptr;

  switch (toupper(Cmd_Str[0])) {
  case 'A':
    BranchOutConvA(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'F':
    BranchOutConvF(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'R':
    BranchOutConvR(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'T':
    BranchOutConvT(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'H':
    ShowOutConvMenu(In_Ptr->in_fname);
    break;
  case 'Q':
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void
OutputConvData(InputStruct * In_Ptr,
	       OutStruct * Out_Ptr)
{
  char        cmd_str[STRLEN];

  if (!Out_Ptr->allocated)
    puts("...No data to output");
  else if (In_Ptr->beam.type == pencil)
    puts("...No incident beam specified");
  else
    do {
      printf("\n> Output convolved data (h for help) => ");
      do
	gets(cmd_str);
      while (!strlen(cmd_str));	/* avoid null string. */
      BranchOutConvCmd(cmd_str, In_Ptr, Out_Ptr);
    } while (toupper(cmd_str[0]) != 'Q');
}

/****************************Contours***************************/
/****************************************************************
 *	Absorption density to fluence. A = F/mua;
 ****/
void
A2Fconv(InputStruct * In_Ptr, double **A_rz)
{
  short       nz = In_Ptr->nz, nrc = In_Ptr->nrc;
  short       ir, iz;
  double      mua;

  for (ir = 0; ir < nrc; ir++)
    for (iz = 0; iz < nz; iz++) {
      mua = In_Ptr->layerspecs[IzToLayer(iz, In_Ptr)].mua;
      if (mua > 0.0)
	A_rz[ir][iz] /= mua;
    }
}

/****************************************************************
 *	Fluence to absorption density. F = A*mua;
 ****/
void
F2Aconv(InputStruct * In_Ptr, double **A_rz)
{
  short       nz = In_Ptr->nz, nrc = In_Ptr->nrc;
  short       ir, iz;
  double      mua;

  for (ir = 0; ir < nrc; ir++)
    for (iz = 0; iz < nz; iz++) {
      mua = In_Ptr->layerspecs[IzToLayer(iz, In_Ptr)].mua;
      if (mua > 0.0)
	A_rz[ir][iz] *= mua;
    }
}

/****************************************************************
 ****/
void
ShowContConvMenu(char *in_fname)
{
  printf("A = absorption vs r & z [J/cm3]\n");
  printf("F = fluence vs r & z [J/cm2]\n");
  printf("R = diffuse reflectance vs radius and angle [J/(cm2 sr)]\n");
  printf("T = transmittance vs radius and angle [J/(cm2 sr)]\n");
  printf("Q   = Quit to main menu\n");
  printf("* input filename: %s \n", in_fname);
}

/****************************************************************
 ****/
void
BranchContConvCmd(char *Cmd_Str,
		  InputStruct * In_Ptr,
		  OutStruct * Out_Ptr)
{
  char        ch;

  ConvVar.in_ptr = In_Ptr;
  ConvVar.out_ptr = Out_Ptr;

  switch (toupper(Cmd_Str[0])) {
  case 'A':
    if (!Out_Ptr->conved.A_rz)
      ConvA_rz(In_Ptr, Out_Ptr);
    IsoPlot(Out_Ptr->A_rzc, In_Ptr->nrc - 1, In_Ptr->nz - 1,
	    In_Ptr->drc, In_Ptr->dz);
    break;
  case 'F':
    if (!Out_Ptr->conved.A_rz)
      ConvA_rz(In_Ptr, Out_Ptr);
    A2Fconv(In_Ptr, Out_Ptr->A_rzc);
    IsoPlot(Out_Ptr->A_rzc, In_Ptr->nrc - 1, In_Ptr->nz - 1,
	    In_Ptr->drc, In_Ptr->dz);
    F2Aconv(In_Ptr, Out_Ptr->A_rzc);
    break;
  case 'R':
    if (!Out_Ptr->conved.Rd_ra)
      ConvRd_ra(In_Ptr, Out_Ptr);
    IsoPlot(Out_Ptr->Rd_rac, In_Ptr->nrc - 1, In_Ptr->na - 1,
	    In_Ptr->drc, In_Ptr->da);
    break;
  case 'T':
    if (!Out_Ptr->conved.Tt_ra)
      ConvTt_ra(In_Ptr, Out_Ptr);
    IsoPlot(Out_Ptr->Tt_rac, In_Ptr->nrc - 1, In_Ptr->na - 1,
	    In_Ptr->drc, In_Ptr->da);
    break;
  case 'H':
    ShowContConvMenu(In_Ptr->in_fname);
    break;
  case 'Q':
    break;
  default:
    puts("...Wrong command");
  }
}
/****************************************************************
 ****/
void
ContourConvData(InputStruct * In_Ptr,
		OutStruct * Out_Ptr)
{
  char        cmd_str[STRLEN];

  if (!Out_Ptr->allocated)
    puts("...No data to output");
  else if (In_Ptr->beam.type == pencil)
    puts("...No incident beam specified");
  else
    do {
      printf("\n> Contour output of convolved data (h for help) => ");
      do
	gets(cmd_str);
      while (!strlen(cmd_str));	/* avoid null string. */
      BranchContConvCmd(cmd_str, In_Ptr, Out_Ptr);
    } while (toupper(cmd_str[0]) != 'Q');
}


/****************************Scanning***************************/
/****************************************************************
 ****/
void
ShowScanConvMenu(char *in_fname)
{
  printf("Ar = absorption vs r @ fixed z [J/cm3]\n");
  printf("Az = absorption vs z @ fixed r [J/cm3]\n");
  printf("Fr = fluence vs r @ fixed z [J/cm2]\n");
  printf("Fz = fluence vs z @ fixed r [J/cm2]\n");
  printf("Rr = diffuse reflectance vs r @ fixed angle [J/(cm2 sr)]\n");
  printf("Ra = diffuse reflectance vs angle @ fixed r [J/(cm2 sr)]\n");
  printf("Tr = transmittance vs r @ fixed angle [J/(cm2 sr)]\n");
  printf("Ta = transmittance vs angle @ fixed r [J/(cm2 sr)]\n");
  printf("Q  = quit\n");
  printf("* input filename: %s \n", in_fname);
}

/****************************************************************
 *	Ext is either "Ars" or "Frs".
 ****/
void
ScanConvA_r(char *Ext, InputStruct * In_Ptr, double **A_rzc)
{
  short       irc, iz, nrc = In_Ptr->nrc, nz = In_Ptr->nz;
  double      r, z, drc = In_Ptr->drc, dz = In_Ptr->dz;
  FILE       *file;

  file = GetWriteFile(Ext);
  if (file == NULL)
    return;

  printf("z grid separation is %-10.4lg cm.\n", dz);
  printf("Input fixed z index (0 - %2hd): ", nz - 1);
  iz = GetShort(0, nz - 1);
  fprintf(file, "%-12s\t%-s@z=%-9.3lg\n", "r[cm]", Ext, dz * (iz + 0.5));
  for (irc = 0; irc < nrc; irc++) {
    r = (irc + 0.5) * drc;
    fprintf(file, "%-12.4E\t%-12.4E\n", r, A_rzc[irc][iz]);
  }

  fclose(file);
}

/****************************************************************
 *	Ext is either "Azs" or "Fzs".
 ****/
void
ScanConvA_z(char *Ext, InputStruct * In_Ptr, double **A_rzc)
{
  short       irc, iz, nrc = In_Ptr->nrc, nz = In_Ptr->nz;
  double      r, z, drc = In_Ptr->drc, dz = In_Ptr->dz;
  FILE       *file;

  file = GetWriteFile(Ext);
  if (file == NULL)
    return;

  printf("r grid separation is %-10.4lg cm.\n", drc);
  printf("Input fixed r index (0 - %2hd): ", nrc - 1);
  irc = GetShort(0, nrc - 1);
  fprintf(file, "%-12s\t%-s@r=%-9.3lg\n", "z[cm]", Ext, drc * (irc + 0.5));
  for (iz = 0; iz < nz; iz++) {
    z = (iz + 0.5) * dz;
    fprintf(file, "%-12.4E\t%-12.4E\n", z, A_rzc[irc][iz]);
  }

  fclose(file);
}

/****************************************************************
 ****/
void
ScanConvRd_r(InputStruct * In_Ptr, double **Rd_rac)
{
  short       irc, ia, nrc = In_Ptr->nrc, na = In_Ptr->na;
  double      r, a, drc = In_Ptr->drc, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

#if IBMPC
  strcpy(fname, "Rrs");
#else
  strcpy(fname, "Rrsc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  printf("Angle grid separation is %-10.4lg rad.\n", da);
  printf("Input fixed angle index (0 - %2hd): ", na - 1);
  ia = GetShort(0, na - 1);

  fprintf(file, "%-12s\t%-s@a=%-9.3lg\n", "r[cm]", fname, da * (ia + 0.5));
  for (irc = 0; irc < nrc; irc++) {
    r = (irc + 0.5) * drc;
    fprintf(file, "%-12.4E\t%-12.4E\n", r, Rd_rac[irc][ia]);
  }

  fclose(file);
}

/****************************************************************
 ****/
void
ScanConvRd_a(InputStruct * In_Ptr, double **Rd_rac)
{
  short       irc, ia, nrc = In_Ptr->nrc, na = In_Ptr->na;
  double      r, a, drc = In_Ptr->drc, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

#if IBMPC
  strcpy(fname, "Ras");
#else
  strcpy(fname, "Rasc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  printf("r grid separation is %-10.4lg cm.\n", drc);
  printf("Input fixed r index (0 - %2hd): ", nrc - 1);
  irc = GetShort(0, nrc - 1);

  fprintf(file, "%-12s\t%-s@r=%-9.3lg\n", "a[rad]", fname, drc * (irc + 0.5));
  for (ia = 0; ia < na; ia++) {
    a = (ia + 0.5) * da;
    fprintf(file, "%-12.4E\t%-12.4E\n", a, Rd_rac[irc][ia]);
  }

  fclose(file);
}

/****************************************************************
 ****/
void
ScanConvTt_r(InputStruct * In_Ptr, double **Tt_rac)
{
  short       irc, ia, nrc = In_Ptr->nrc, na = In_Ptr->na;
  double      r, a, drc = In_Ptr->drc, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

#if IBMPC
  strcpy(fname, "Trs");
#else
  strcpy(fname, "Trsc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  printf("Angle grid separation is %-10.4lg rad.\n", da);
  printf("Input fixed angle index (0 - %2hd): ", na - 1);
  ia = GetShort(0, na - 1);

  fprintf(file, "%-12s\t%-s@a=%-9.3lg\n", "r[cm]", fname, da * (ia + 0.5));
  for (irc = 0; irc < nrc; irc++) {
    r = (irc + 0.5) * drc;
    fprintf(file, "%-12.4E\t%-12.4E\n", r, Tt_rac[irc][ia]);
  }

  fclose(file);
}

/****************************************************************
 ****/
void
ScanConvTt_a(InputStruct * In_Ptr, double **Tt_rac)
{
  short       irc, ia, nrc = In_Ptr->nrc, na = In_Ptr->na;
  double      r, a, drc = In_Ptr->drc, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

#if IBMPC
  strcpy(fname, "Tas");
#else
  strcpy(fname, "Tasc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  printf("r grid separation is %-10.4lg cm.\n", drc);
  printf("Input fixed r index (0 - %2hd): ", nrc - 1);
  irc = GetShort(0, nrc - 1);

  fprintf(file, "%-12s\t%-s@r=%-9.3lg\n", "a[rad]", fname, drc * (irc + 0.5));
  for (ia = 0; ia < na; ia++) {
    a = (ia + 0.5) * da;
    fprintf(file, "%-12.4E\t%-12.4E\n", a, Tt_rac[irc][ia]);
  }

  fclose(file);
}

/****************************************************************
 ****/
void
BranchScanConvA(char *Cmd_Str,
		InputStruct * In_Ptr,
		OutStruct * Out_Ptr)
{
  char        fname[STRLEN];

  if (!Out_Ptr->conved.A_rz)
    ConvA_rz(In_Ptr, Out_Ptr);

  switch (toupper(Cmd_Str[1])) {
  case 'R':
#if IBMPC
    strcpy(fname, "Ars");
    ScanConvA_r(fname, In_Ptr, Out_Ptr->A_rzc);
#else
    strcpy(fname, "Arsc");
    ScanConvA_r(fname, In_Ptr, Out_Ptr->A_rzc);
#endif
    break;
  case 'Z':
#if IBMPC
    strcpy(fname, "Azs");
    ScanConvA_z(fname, In_Ptr, Out_Ptr->A_rzc);
#else
    strcpy(fname, "Azsc");
    ScanConvA_z(fname, In_Ptr, Out_Ptr->A_rzc);
#endif
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void
BranchScanConvF(char *Cmd_Str,
		InputStruct * In_Ptr,
		OutStruct * Out_Ptr)
{
  char        fname[STRLEN];

  if (!Out_Ptr->conved.A_rz)
    ConvA_rz(In_Ptr, Out_Ptr);
  A2Fconv(In_Ptr, Out_Ptr->A_rzc);

  switch (toupper(Cmd_Str[1])) {
  case 'R':
#if IBMPC
    strcpy(fname, "Frs");
    ScanConvA_r(fname, In_Ptr, Out_Ptr->A_rzc);
#else
    strcpy(fname, "Frsc");
    ScanConvA_r(fname, In_Ptr, Out_Ptr->A_rzc);
#endif
    break;
  case 'Z':
#if IBMPC
    strcpy(fname, "Fzs");
    ScanConvA_z(fname, In_Ptr, Out_Ptr->A_rzc);
#else
    strcpy(fname, "Fzsc");
    ScanConvA_z(fname, In_Ptr, Out_Ptr->A_rzc);
#endif
    break;
  default:
    puts("...Wrong command");
  }

  F2Aconv(In_Ptr, Out_Ptr->A_rzc);
}

/****************************************************************
 ****/
void
BranchScanConvR(char *Cmd_Str,
		InputStruct * In_Ptr,
		OutStruct * Out_Ptr)
{
  if (!Out_Ptr->conved.Rd_ra)
    ConvRd_ra(In_Ptr, Out_Ptr);
  switch (toupper(Cmd_Str[1])) {
  case 'R':
    ScanConvRd_r(In_Ptr, Out_Ptr->Rd_rac);
    break;
  case 'A':
    ScanConvRd_a(In_Ptr, Out_Ptr->Rd_rac);
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void
BranchScanConvT(char *Cmd_Str,
		InputStruct * In_Ptr,
		OutStruct * Out_Ptr)
{
  if (!Out_Ptr->conved.Tt_ra)
    ConvTt_ra(In_Ptr, Out_Ptr);
  switch (toupper(Cmd_Str[1])) {
  case 'R':
    ScanConvTt_r(In_Ptr, Out_Ptr->Tt_rac);
    break;
  case 'A':
    ScanConvTt_a(In_Ptr, Out_Ptr->Tt_rac);
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void
BranchScanConvCmd(char *Cmd_Str,
		  InputStruct * In_Ptr,
		  OutStruct * Out_Ptr)
{
  char        ch;

  ConvVar.in_ptr = In_Ptr;
  ConvVar.out_ptr = Out_Ptr;

  switch (toupper(Cmd_Str[0])) {
  case 'A':
    BranchScanConvA(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'F':
    BranchScanConvF(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'R':
    BranchScanConvR(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'T':
    BranchScanConvT(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'H':
    ShowScanConvMenu(In_Ptr->in_fname);
    break;
  case 'Q':
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void
ScanConvData(InputStruct * In_Ptr,
	     OutStruct * Out_Ptr)
{
  char        cmd_str[STRLEN];

  if (!Out_Ptr->allocated)
    puts("...No data to output");
  else if (In_Ptr->beam.type == pencil)
    puts("...No incident beam specified");
  else
    do {
      printf("\n> Scans of convolved data (h for help) => ");
      do
	gets(cmd_str);
      while (!strlen(cmd_str));	/* avoid null string. */
      BranchScanConvCmd(cmd_str, In_Ptr, Out_Ptr);
    } while (toupper(cmd_str[0]) != 'Q');
}
