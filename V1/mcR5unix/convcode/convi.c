/****************************************************************
 *	Functions for reading files.
 ****/

#include "conv.h"

/****************************************************************
 *	Get the filename and open it for reading, retry until the
 *	file can be opened, or a '.' is input.
 ****/
FILE       *
GetFile(char *Fname)
{
  FILE       *file = NULL;

  do {
    printf("Input filename of mcml output(or . to quit): ");
    scanf("%s", Fname);
    if (strlen(Fname) == 1 && Fname[0] == '.')
      break;

    file = fopen(Fname, "r");
  } while (file == NULL);

  return (file);
}

/****************************************************************
 *  Allocate the arrays for the original data from mcml
 *	computation in OutStruct for each run, and they are
 *	automatically initialized to zeros.
 *
 *  Returns 1 if successful, otherwise 0.
 ****/
Boolean 
AllocOrigData(InputStruct * In_Ptr,
	      OutStruct * Out_Ptr)
{
  short       nz = In_Ptr->nz;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nl = In_Ptr->num_layers;	/* remember +2 for ambient. */

  if (nz <= 0 || nr <= 0 || na <= 0 || nl <= 0) {
    printf("Wrong grid parameters.\n");
    return (0);
  }
  /* Allocate the arrays and the matrices. */
  Out_Ptr->Rd_ra = AllocMatrix(0, nr - 1, 0, na - 1);
  Out_Ptr->Rd_r = AllocVector(0, nr - 1);
  Out_Ptr->Rd_a = AllocVector(0, na - 1);

  Out_Ptr->A_rz = AllocMatrix(0, nr - 1, 0, nz - 1);
  Out_Ptr->A_z = AllocVector(0, nz - 1);
  Out_Ptr->A_l = AllocVector(0, nl + 1);

  Out_Ptr->Tt_ra = AllocMatrix(0, nr - 1, 0, na - 1);
  Out_Ptr->Tt_r = AllocVector(0, nr - 1);
  Out_Ptr->Tt_a = AllocVector(0, na - 1);

  return (1);
}

/****************************************************************
 *  Allocate the arrays for convolution in OutStruct for each
 *	run, and they are automatically initialized to zeros.
 *
 *  Returns 1 if successful, otherwise 0.
 ****/
Boolean 
AllocConvData(InputStruct * In_Ptr,
	      OutStruct * Out_Ptr)
{
  short       nz = In_Ptr->nz;
  short       nr = In_Ptr->nrc;	/* use nrc instead of nr. */
  short       na = In_Ptr->na;

  if (nz <= 0 || nr <= 0 || na <= 0) {
    printf("Wrong grid parameters.\n");
    return (0);
  }
  Out_Ptr->Rd_rac = AllocMatrix(0, nr - 1, 0, na - 1);
  Out_Ptr->Rd_rc = AllocVector(0, nr - 1);
  Out_Ptr->A_rzc = AllocMatrix(0, nr - 1, 0, nz - 1);
  Out_Ptr->Tt_rac = AllocMatrix(0, nr - 1, 0, na - 1);
  Out_Ptr->Tt_rc = AllocVector(0, nr - 1);

  return (1);
}

/****************************************************************
 *  Free the arrays in OutStruct.
 ****/
void 
FreeOrigData(InputStruct * In_Ptr,
	     OutStruct * Out_Ptr)
{
  short       nz = In_Ptr->nz;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nl = In_Ptr->num_layers;	/* remember +2 for ambient. */

  if (nz <= 0 || nr <= 0 || na <= 0 || nl <= 0) {
    printf("Wrong grid parameters.\n");
    return;
  }
  if (Out_Ptr->Rd_ra != NULL)
    FreeMatrix(Out_Ptr->Rd_ra, 0, nr - 1, 0, na - 1);
  if (Out_Ptr->Rd_r != NULL)
    FreeVector(Out_Ptr->Rd_r, 0, nr - 1);
  if (Out_Ptr->Rd_a != NULL)
    FreeVector(Out_Ptr->Rd_a, 0, na - 1);

  if (Out_Ptr->A_rz != NULL)
    FreeMatrix(Out_Ptr->A_rz, 0, nr - 1, 0, nz - 1);
  if (Out_Ptr->A_z != NULL)
    FreeVector(Out_Ptr->A_z, 0, nz - 1);
  if (Out_Ptr->A_l != NULL)
    FreeVector(Out_Ptr->A_l, 0, nl + 1);

  if (Out_Ptr->Tt_ra != NULL)
    FreeMatrix(Out_Ptr->Tt_ra, 0, nr - 1, 0, na - 1);
  if (Out_Ptr->Tt_r != NULL)
    FreeVector(Out_Ptr->Tt_r, 0, nr - 1);
  if (Out_Ptr->Tt_a != NULL)
    FreeVector(Out_Ptr->Tt_a, 0, na - 1);
}

/****************************************************************
 *  Free the arrays in OutStruct.
 ****/
void 
FreeConvData(InputStruct * In_Ptr,
	     OutStruct * Out_Ptr)
{
  short       nz = In_Ptr->nz;
  short       nr = In_Ptr->nrc;
  short       na = In_Ptr->na;

  if (nz <= 0 || nr <= 0 || na <= 0) {
    printf("Wrong grid parameters.\n");
    return;
  }
  if (Out_Ptr->Rd_rac != NULL)
    FreeMatrix(Out_Ptr->Rd_rac, 0, nr - 1, 0, na - 1);
  if (Out_Ptr->Rd_rc != NULL)
    FreeVector(Out_Ptr->Rd_rc, 0, nr - 1);
  if (Out_Ptr->A_rzc != NULL)
    FreeMatrix(Out_Ptr->A_rzc, 0, nr - 1, 0, nz - 1);
  if (Out_Ptr->Tt_rac != NULL)
    FreeMatrix(Out_Ptr->Tt_rac, 0, nr - 1, 0, na - 1);
  if (Out_Ptr->Tt_rc != NULL)
    FreeVector(Out_Ptr->Tt_rc, 0, nr - 1);
}

/****************************************************************
 * 	Kill the ith char (counting from 0), push the following chars
 *	forward by one.
 ****/
void 
KillChar(size_t i, char *Str)
{
  size_t      sl = strlen(Str);

  for (; i < sl; i++)
    Str[i] = Str[i + 1];
}

/****************************************************************
 *	Eliminate the chars in a string which are not printing chars
 *	or spaces.
 *
 *	Space include ' ', '\f', '\t' etc.
 *
 *	Return 1 if no nonprinting chars found, otherwise return 0.
 ****/
Boolean 
CheckChar(char *Str)
{
  Boolean     found = 0;
  size_t      sl = strlen(Str);
  size_t      i = 0;

  while (i < sl)
    if (isprint(Str[i]) || isspace(Str[i]))
      i++;
    else {
      found = 1;
      KillChar(i, Str);
      sl--;
    }

  return (found);

}

/****************************************************************
 *	Return 1 if this line is comment line in which the first
 *	non-space character is "#".
 *
 *	Also return 1 if this line is space line.
 ****/
Boolean 
CommentLine(char *Buf)
{
  size_t      spn, cspn;

  spn = strspn(Buf, " \t");
  /* length spanned by space or tab chars. */
  cspn = strcspn(Buf, "#\n");
  /* length before the 1st # or return. */
  if (spn == cspn)
    return (1);
  /* comment line or space line. */
  else
    return (0);			/* the line has data. */
}

/****************************************************************
 *	Skip space or comment lines and return data line only.
 ****/
char       *
FindDataLine(FILE * File_Ptr)
{
  char        buf[STRLEN];

  do {				/* skip space or comment lines. */
    if (fgets(buf, STRLEN, File_Ptr) == NULL) {
      printf("Incomplete data.\n");
      exit(1);
    }
    CheckChar(buf);
  } while (CommentLine(buf));

  return (buf);
}

/****************************************************************
 *	Look for the key word, which is the 1st word in the line.
 ****/
Boolean 
FoundKeyWord(char *LineBuf, char *Key)
{
  char       *sub_str;
  Boolean     found = 0;

  if ((sub_str = strstr(LineBuf, Key)) != NULL)
    if (sub_str == LineBuf)
      found = 1;

  return (found);
}

/****************************************************************
 *
 ****/
void 
SeekKey(FILE * File_Ptr, char *Key)
{
  char        buf[STRLEN];

  do
    strcpy(buf, FindDataLine(File_Ptr));
  while (!FoundKeyWord(buf, Key));
}

/****************************************************************
 *
 ****/
void 
ReadLayerParm(FILE * File_Ptr,
	      short Num_Layers,
	      LayerStruct ** Layers_PP)
{
  char        buf[STRLEN];
  short       i = 0;
  double      d, n, mua, mus, g;
  double      z = 0.0;		/* z coordinate of the current layer. */

  /* layer 0 and layer Num_Layers + 1 are for ambient. */
  *Layers_PP = (LayerStruct *)
    malloc((unsigned) (Num_Layers + 2) * sizeof(LayerStruct));
  if (!(*Layers_PP))
    nrerror("allocation failure in ReadLayerParm()");

  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%lf", &n);
  (*Layers_PP)[i].n = n;
  for (i = 1; i <= Num_Layers; i++) {
    strcpy(buf, FindDataLine(File_Ptr));
    sscanf(buf, "%lf%lf%lf%lf%lf", &n, &mua, &mus, &g, &d);
    (*Layers_PP)[i].n = n;
    (*Layers_PP)[i].mua = mua;
    (*Layers_PP)[i].mus = mus;
    (*Layers_PP)[i].g = g;
    (*Layers_PP)[i].z0 = z;
    z += d;
    (*Layers_PP)[i].z1 = z;
  }
  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%lf", &n);
  (*Layers_PP)[i].n = n;
}

/****************************************************************
 *	Read in the input parameters which were used for Monte Carlo
 *	simulations.
 ****/
void 
ReadInParm(FILE * File_Ptr, InputStruct * In_Ptr)
{
  char        buf[STRLEN];

  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%c", &(In_Ptr->in_fformat));
  if (toupper(In_Ptr->in_fformat) != 'B')
    In_Ptr->in_fformat = 'A';

  /** Find the key word "InParm". */
  do
    strcpy(buf, FindDataLine(File_Ptr));
  while (!FoundKeyWord(buf, "InParm"));

  /** Escape filename & file format. */
  FindDataLine(File_Ptr);

  /** read in number of photons. **/
  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%ld", &In_Ptr->num_photons);

  /** assign in Wth (critical weight). **/
  In_Ptr->Wth = 1E-4;

  /** read in dz, dr. **/
  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%lf%lf", &In_Ptr->dz, &In_Ptr->dr);

  /** read in nz, nr, na and compute da. **/
  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%hd%hd%hd", &In_Ptr->nz,
	 &In_Ptr->nr, &In_Ptr->na);
  In_Ptr->da = 0.5 * PI / In_Ptr->na;

  /** read in number of layers. **/
  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%hd", &In_Ptr->num_layers);

  ReadLayerParm(File_Ptr, In_Ptr->num_layers,
		&In_Ptr->layerspecs);
}

/****************************************************************
 *	Read reflectance, absorbed fraction, transmittance.
 ****/
void 
ReadRAT(FILE * File_Ptr,
	OutStruct * Out_Ptr)
{
  char        buf[STRLEN];

  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%lf", &(Out_Ptr->Rsp));

  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%lf", &(Out_Ptr->Rd));

  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%lf", &(Out_Ptr->A));

  strcpy(buf, FindDataLine(File_Ptr));
  sscanf(buf, "%lf", &(Out_Ptr->Tt));
}

/****************************************************************
 ****/
void 
ReadA_layer(FILE * File_Ptr,
	    short Num_Layers,
	    OutStruct * Out_Ptr)
{
  char        buf[STRLEN];
  short       i;

  for (i = 1; i <= Num_Layers; i++) {
    strcpy(buf, FindDataLine(File_Ptr));
    sscanf(buf, "%lf", &(Out_Ptr->A_l[i]));
  }
}

/****************************************************************
 ****/
void 
ReadA_z(FILE * File_Ptr,
	short Nz,
	OutStruct * Out_Ptr)
{
  char        buf[STRLEN];
  short       i;

  for (i = 0; i < Nz; i++) {
    strcpy(buf, FindDataLine(File_Ptr));
    sscanf(buf, "%lf", &(Out_Ptr->A_z[i]));
  }
}

/****************************************************************
 ****/
void 
ReadRd_r(FILE * File_Ptr,
	 short Nr,
	 OutStruct * Out_Ptr)
{
  char        buf[STRLEN];
  short       i;

  for (i = 0; i < Nr; i++) {
    strcpy(buf, FindDataLine(File_Ptr));
    sscanf(buf, "%lf", &(Out_Ptr->Rd_r[i]));
  }
}

/****************************************************************
 ****/
void 
ReadRd_a(FILE * File_Ptr,
	 short Na,
	 OutStruct * Out_Ptr)
{
  char        buf[STRLEN];
  short       i;

  for (i = 0; i < Na; i++) {
    strcpy(buf, FindDataLine(File_Ptr));
    sscanf(buf, "%lf", &(Out_Ptr->Rd_a[i]));
  }
}

/****************************************************************
 ****/
void 
ReadTt_r(FILE * File_Ptr,
	 short Nr,
	 OutStruct * Out_Ptr)
{
  char        buf[STRLEN];
  short       i;

  for (i = 0; i < Nr; i++) {
    strcpy(buf, FindDataLine(File_Ptr));
    sscanf(buf, "%lf", &(Out_Ptr->Tt_r[i]));
  }
}

/****************************************************************
 ****/
void 
ReadTt_a(FILE * File_Ptr,
	 short Na,
	 OutStruct * Out_Ptr)
{
  char        buf[STRLEN];
  short       i;

  for (i = 0; i < Na; i++) {
    /*
     * strcpy(buf, FindDataLine(File_Ptr)); sscanf(buf, "%lf",
     * &(Out_Ptr->Tt_a[i]));
     */
    fscanf(File_Ptr, "%lf", &(Out_Ptr->Tt_a[i]));
  }
}

/****************************************************************
 ****/
void 
ReadA_rz(FILE * File_Ptr,
	 short Nr,
	 short Nz,
	 OutStruct * Out_Ptr)
{
  short       iz, ir;

  for (ir = 0; ir < Nr; ir++)
    for (iz = 0; iz < Nz; iz++)
      fscanf(File_Ptr, "%lf ", &(Out_Ptr->A_rz[ir][iz]));
}

/****************************************************************
 ****/
void 
ReadRd_ra(FILE * File_Ptr,
	  short Nr,
	  short Na,
	  OutStruct * Out_Ptr)
{
  short       ir, ia;

  for (ir = 0; ir < Nr; ir++)
    for (ia = 0; ia < Na; ia++)
      fscanf(File_Ptr, "%lf ", &(Out_Ptr->Rd_ra[ir][ia]));
}

/****************************************************************
 ****/
void 
ReadTt_ra(FILE * File_Ptr,
	  short Nr,
	  short Na,
	  OutStruct * Out_Ptr)
{
  short       ir, ia;

  for (ir = 0; ir < Nr; ir++)
    for (ia = 0; ia < Na; ia++)
      fscanf(File_Ptr, "%lf ", &(Out_Ptr->Tt_ra[ir][ia]));
}

/****************************************************************
 *	Read in the Monte Carlo output parameters.
 ****/
void 
ReadOutMC(FILE * File_Ptr,
	  InputStruct * In_Ptr,
	  OutStruct * Out_Ptr)
{
  ConvStruct  null_conved = NULLCONVSTRUCT;

  SeekKey(File_Ptr, "RAT");
  ReadRAT(File_Ptr, Out_Ptr);	/* refl.,absorption,transmission. */

  /* 1D arrays. */
  SeekKey(File_Ptr, "A_l");
  ReadA_layer(File_Ptr, In_Ptr->num_layers, Out_Ptr);

  SeekKey(File_Ptr, "A_z");
  ReadA_z(File_Ptr, In_Ptr->nz, Out_Ptr);

  SeekKey(File_Ptr, "Rd_r");
  ReadRd_r(File_Ptr, In_Ptr->nr, Out_Ptr);

  SeekKey(File_Ptr, "Rd_a");
  ReadRd_a(File_Ptr, In_Ptr->na, Out_Ptr);

  SeekKey(File_Ptr, "Tt_r");
  ReadTt_r(File_Ptr, In_Ptr->nr, Out_Ptr);

  SeekKey(File_Ptr, "Tt_a");
  ReadTt_a(File_Ptr, In_Ptr->na, Out_Ptr);

  /* 2D arrays. */
  SeekKey(File_Ptr, "A_rz");
  ReadA_rz(File_Ptr, In_Ptr->nr, In_Ptr->nz, Out_Ptr);

  SeekKey(File_Ptr, "Rd_ra");
  ReadRd_ra(File_Ptr, In_Ptr->nr, In_Ptr->na, Out_Ptr);

  SeekKey(File_Ptr, "Tt_ra");
  ReadTt_ra(File_Ptr, In_Ptr->nr, In_Ptr->na, Out_Ptr);

  Out_Ptr->conved = null_conved;
}

/****************************************************************
 *	After the Input parameters are read in, the parameters drc
 *	and nrc are initialized to dr and nr respectively.
 ****/
void 
ReadMcoFile(InputStruct * In_Ptr,
	    OutStruct * Out_Ptr)
{
  FILE       *infile;

  infile = GetFile(In_Ptr->in_fname);
  if (infile == NULL)
    return;

  if (Out_Ptr->allocated) {
    OutStruct   null_out = NULLOUTSTRUCT;
    FreeOrigData(In_Ptr, Out_Ptr);
    *Out_Ptr = null_out;
  }
  ReadInParm(infile, In_Ptr);
  In_Ptr->beam.type = pencil;
  In_Ptr->drc = In_Ptr->dr;
  In_Ptr->nrc = In_Ptr->nr;

  AllocOrigData(In_Ptr, Out_Ptr);
  AllocConvData(In_Ptr, Out_Ptr);
  Out_Ptr->allocated = 1;

  ReadOutMC(infile, In_Ptr, Out_Ptr);
  fclose(infile);
}
