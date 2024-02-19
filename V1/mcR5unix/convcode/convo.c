/****************************************************************
 *	Functions for file output.
 ****/

#include "conv.h"

/****************************************************************
 *	Center a string according to the column width.
 ****/
char       *
CenterStr(short int Wid,
	  char *InStr,
	  char *OutStr)
{
  size_t      nspaces;
  /* number of spaces to be filled before InStr. */

  nspaces = (Wid - strlen(InStr)) / 2;
  if (nspaces < 0)
    nspaces = 0;

  strcpy(OutStr, "");
  while (nspaces--)
    strcat(OutStr, " ");

  strcat(OutStr, InStr);

  return (OutStr);
}

/****************************************************************
 *	Print some messages before starting simulation.  e.g. author,
 *	address, program version, year.
 ****/
#define COLWIDTH 80
void 
ShowVersion(void)
{
  char        str[STRLEN];

  CenterStr(COLWIDTH,
	    "Convolution of mcml Simulation Data", str);
  puts("");
  puts(str);
  puts("");

  CenterStr(COLWIDTH, "Lihong Wang, Ph.D.", str);
  puts(str);

  CenterStr(COLWIDTH, "Steven L. Jacques, Ph.D.", str);
  puts(str);

  CenterStr(COLWIDTH, "Laser Biology Research Laboratory - 017", str);
  puts(str);

  CenterStr(COLWIDTH, "University of Texas M.D. Anderson Cancer Center", str);
  puts(str);

  CenterStr(COLWIDTH, "Houston, Texas 77030, USA", str);
  puts(str);

  CenterStr(COLWIDTH, "Voice: (713) 792-3664, Fax: (713) 792-3995", str);
  puts(str);

  puts("");

  CenterStr(COLWIDTH, "Version 1.1, 1994", str);
  puts(str);
  puts("\n\n");
}
#undef COLWIDTH

/****************************************************************
 *	Open a file for output with extension Ext.  If file exists,
 *	ask whether to overwrite or append or change filename.
 *
 *	Return file pointer, which could be NULL.
 *	Return the full filename as Ext.
 ****/
FILE       *
GetWriteFile(char *Ext)
{
  FILE       *file;
  char        fname[STRLEN], fmode[STRLEN];

  do {
    printf("Enter output filename with extension .%s (or . to quit): ", Ext);
    gets(fname);
    if (strlen(fname) == 1 && fname[0] == '.') {
      fmode[0] = 'q';
      break;
    } else
      fmode[0] = 'w';

    if ((file = fopen(fname, "r")) != NULL) {	/* file exists. */
      printf("File %s exists, %s",
	     fname, "w=overwrite, a=append, n=new filename, q=quit: ");
      do
	gets(fmode);
      while (!strlen(fmode));	/* avoid null line. */
      fclose(file);
    }
  } while (fmode[0] != 'w' && fmode[0] != 'a' && fmode[0] != 'q');

  if (fmode[0] != 'q')
    file = fopen(fname, fmode);
  else
    file = NULL;

  strcpy(Ext, fname);
  return (file);		/* could be NULL. */
}

/****************************************************************
 *	Return the index to the layer, where iz is in.
 *
 *	Use the center of box.
 ****/
short 
IzToLayer(short Iz,
	  InputStruct * In_Ptr)
{
  short       i = 1;		/* index to layer. */
  short       num_layers = In_Ptr->num_layers;
  double      dz = In_Ptr->dz;

  while ((Iz + 0.5) * dz >= In_Ptr->layerspecs[i].z1
	 && i < num_layers)
    i++;

  return (i);
}

/****************************************************************
 *	Write the input parameter for Monte Carlo simulation program
 *	in such a format that it can be read directly back.
 ****/
void 
WriteInParm(InputStruct * In_Ptr)
{
  short       i;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "InP");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "1.1\t# file version.\n");
  fprintf(file, "1\t# number of runs.\n\n");
  fprintf(file, "temp.out\tA\t\t#output filename.\n");

  fprintf(file, "%ld \t\t\t# No. of photons.\n",
	  In_Ptr->num_photons);

  fprintf(file, "%G\t%G\t\t# dz, dr.\n", In_Ptr->dz, In_Ptr->dr);
  fprintf(file, "%hd\t%hd\t%hd\t# No. of dz, dr, da.\n\n",
	  In_Ptr->nz, In_Ptr->nr, In_Ptr->na);

  fprintf(file, "%hd\t\t\t\t\t# Number of layers.\n",
	  In_Ptr->num_layers);
  fprintf(file,
	  "#n\tmua\tmus\tg\td\t# One line for each layer.\n");
  fprintf(file, "%G\t\t\t\t\t# n for medium above.\n",
	  In_Ptr->layerspecs[0].n);
  for (i = 1; i <= In_Ptr->num_layers; i++) {
    LayerStruct s;
    s = In_Ptr->layerspecs[i];
    fprintf(file, "%G\t%G\t%G\t%G\t%G\t# layer %hd\n",
	    s.n, s.mua, s.mus, s.g, s.z1 - s.z0, i);
  }
  fprintf(file, "%G\t\t\t\t\t# n for medium below.\n\n",
	  In_Ptr->layerspecs[i].n);

  fclose(file);
}

/****************************************************************
 *	Write reflectance, absorption, transmission.
 ****/
void 
WriteRAT(OutStruct * Out_Ptr)
{
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "RAT");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file,			/* flag. */
	  "RAT #Reflectance, absorption, transmission.\n");
  fprintf(file, "%-12.4G \t#Specular reflectance.\n",
	  Out_Ptr->Rsp);
  fprintf(file, "%-12.4G \t#Diffuse reflectance.\n",
	  Out_Ptr->Rd);
  fprintf(file, "%-12.4G \t#Absorption.\n",
	  Out_Ptr->A);
  fprintf(file, "%-12.4G \t#Transmission.\n",
	  Out_Ptr->Tt);

  fprintf(file, "\n");

  fclose(file);
}

/****************************************************************
 *	Write absorption as a function of layer.
 *	2 numbers each line: layer, A[layer].
 ****/
void 
WriteA_layer(short Num_Layers,
	     double *A_l)
{
  short       i;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Al");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "layer\t%-s[-]\n", fname);
  for (i = 1; i <= Num_Layers; i++)
    fprintf(file, "%-4hd\t%-12.4G\n", i, A_l[i]);

  fclose(file);
}

/****************************************************************
 *	2 numbers each line: z, A[z].
 ****/
void 
WriteA_z(InputStruct * In_Ptr, double *A_z)
{
  short       nz = In_Ptr->nz;
  double      dz = In_Ptr->dz;
  short       iz;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Az");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "%-12s\t%-s[1/cm]\n", "z[cm]", fname);
  for (iz = 0; iz < nz; iz++)
    fprintf(file, "%-12.4E\t%-12.4E\n", (iz + 0.5) * dz, A_z[iz]);

  fclose(file);
}

/****************************************************************
 *	3 numbers each line: r, z, A[r][z].
 ****/
void 
WriteA_rz(InputStruct * In_Ptr,
	  double **A_rz)
{
  short       ir, iz, nz = In_Ptr->nz, nr = In_Ptr->nr;
  double      r, z, dr = In_Ptr->dr, dz = In_Ptr->dz;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Arz");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "%-12s\t%-12s\t%-s[1/cm3]\n", "r[cm]", "z[cm]", fname);
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    for (iz = 0; iz < nz; iz++) {
      z = (iz + 0.5) * dz;
      fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n",
	      r, z, A_rz[ir][iz]);
    }
  }

  fclose(file);
}

/****************************************************************
 *	2 numbers each line: z, F[z].
 ****/
void 
WriteF_z(InputStruct * In_Ptr,
	 double *A_z)
{
  FILE       *file;
  short       iz, nz = In_Ptr->nz;
  double      mua, dz = In_Ptr->dz;
  char        fname[STRLEN];

  strcpy(fname, "Fz");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "%-12s\t%-s[-]\n", "z[cm]", fname);
  for (iz = 0; iz < nz; iz++) {
    mua = In_Ptr->layerspecs[IzToLayer(iz, In_Ptr)].mua;
    if (mua > 0.0)
      fprintf(file, "%-12.4E\t%-12.4E\n", (iz + 0.5) * dz,
	      A_z[iz] / mua);
    else
      fprintf(file, "%-12.4E\t%-12.4E\n", (iz + 0.5) * dz, 0.0);
  }

  fclose(file);
}

/****************************************************************
 *	3 numbers each line: r, z, F[r][z].
 ****/
void 
WriteF_rz(InputStruct * In_Ptr,
	  double **A_rz)
{
  FILE       *file;
  short       ir, iz, nz = In_Ptr->nz, nr = In_Ptr->nr;
  double      mua, r, z, dr = In_Ptr->dr, dz = In_Ptr->dz;
  char        fname[STRLEN];

  strcpy(fname, "Frz");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "%-12s\t%-12s\t%-s[1/cm2]\n",
	  "r[cm]", "z[cm]", fname);
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    for (iz = 0; iz < nz; iz++) {
      z = (iz + 0.5) * dz;
      mua = In_Ptr->layerspecs[IzToLayer(iz, In_Ptr)].mua;
      if (mua > 0.0)
	fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n",
		r, z, A_rz[ir][iz] / mua);
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
WriteRd_ra(InputStruct * In_Ptr,
	   double **Rd_ra)
{
  short       ir, ia, nr = In_Ptr->nr, na = In_Ptr->na;
  double      r, a, dr = In_Ptr->dr, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Rra");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "%-12s\t%-12s\t%-s[1/(cm2sr)]\n",
	  "r[cm]", "a[rad]", fname);
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    for (ia = 0; ia < na; ia++) {
      a = (ia + 0.5) * da;
      fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n",
	      r, a, Rd_ra[ir][ia]);
    }
  }

  fclose(file);
}

/****
 *	2 numbers each line: r, Rd[r]
 ****/
void 
WriteRd_r(InputStruct * In_Ptr,
	  double *Rd_r)
{
  short       ir, nr = In_Ptr->nr;
  double      dr = In_Ptr->dr;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Rr");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "%-12s\t%-s[1/cm2]\n", "r[cm]", fname);
  for (ir = 0; ir < nr; ir++)
    fprintf(file, "%-12.4E\t%-12.4E\n", (ir + 0.5) * dr,
	    Rd_r[ir]);

  fclose(file);
}

/****************************************************************
 *	2 numbers each line: a, Rd[a].
 ****/
void 
WriteRd_a(InputStruct * In_Ptr,
	  double *Rd_a)
{
  short       ia, na = In_Ptr->na;
  double      da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Ra");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "%-12s\t%-s[1/sr]\n", "a[rad]", fname);
  for (ia = 0; ia < na; ia++)
    fprintf(file, "%-12.4E\t%-12.4E\n", (ia + 0.5) * da,
	    Rd_a[ia]);

  fclose(file);
}

/****************************************************************
 *	3 numbers each line:r, a, Tt[r][a]. a = theta.
 ****/
void 
WriteTt_ra(InputStruct * In_Ptr,
	   double **Tt_ra)
{
  short       ir, ia, nr = In_Ptr->nr, na = In_Ptr->na;
  double      r, a, dr = In_Ptr->dr, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Tra");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "%-12s\t%-12s\t%-s[1/(cm2sr)]\n",
	  "r[cm]", "a[rad]", fname);
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    for (ia = 0; ia < na; ia++) {
      a = (ia + 0.5) * da;
      fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n",
	      r, a, Tt_ra[ir][ia]);
    }
  }

  fclose(file);
}

/****
 *	2 numbers each line: r, Tt[r].
 ****/
void 
WriteTt_r(InputStruct * In_Ptr,
	  double *Tt_r)
{
  short       ir, nr = In_Ptr->nr;
  double      dr = In_Ptr->dr;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Tr");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "%-12s\t%-s[1/cm2]\n", "r[cm]", fname);
  for (ir = 0; ir < nr; ir++)
    fprintf(file, "%-12.4E\t%-12.4E\n", (ir + 0.5) * dr,
	    Tt_r[ir]);

  fclose(file);
}

/****************************************************************
 *	2 numbers each line: theta, Tt[theta].
 ****/
void 
WriteTt_a(InputStruct * In_Ptr,
	  double *Tt_a)
{
  short       ia, na = In_Ptr->na;
  double      da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Ta");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fprintf(file, "%-12s\t%-s[1/sr]\n", "a[rad]", fname);
  for (ia = 0; ia < na; ia++)
    fprintf(file, "%-12.4E\t%-12.4E\n", (ia + 0.5) * da,
	    Tt_a[ia]);

  fclose(file);
}

/****************************************************************
 *	Write output in M. Keijzer's format so that the file can be
 *	read by the convolution program written by Keijzer in Pascal.
 ****/
void 
WriteKFormat(InputStruct * In_Ptr,
	     OutStruct * Out_Ptr)
{
  short       i, j;
  double      dz, dr;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "K");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  fputs("output.filename\n", file);
  fprintf(file, "%12hd layers\n", In_Ptr->num_layers);
  fprintf(file, "%12s %12s %12s %12s %12s %12s\n",
	  "layer", "mua", "mus", "g", "nt", "thickness");

  for (i = 1; i <= In_Ptr->num_layers; i++)
    fprintf(file,
	    "%12hd %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",
	    i, In_Ptr->layerspecs[i].mua,
	    In_Ptr->layerspecs[i].mus,
	    In_Ptr->layerspecs[i].g, In_Ptr->layerspecs[i].n,
	    In_Ptr->layerspecs[i].z1 -
	    In_Ptr->layerspecs[i].z0);
  fprintf(file, "%12.6lf index of refraction above\n",
	  In_Ptr->layerspecs[0].n);
  fprintf(file, "%12.6lf index of refraction below\n",
	  In_Ptr->layerspecs[i].n);
  fprintf(file, "\n");

  fprintf(file, "%12ld photons\n", In_Ptr->num_photons);
  fprintf(file, "%12.6lf critical weight\n", In_Ptr->Wth);
  fprintf(file, "%12.6lf depth of boxes micron\n",
	  In_Ptr->dz * 1e4);
  fprintf(file, "%12.6lf width of boxes micron\n",
	  In_Ptr->dr * 1e4);
  fprintf(file, "%12hd number of boxes in z \n", In_Ptr->nz);
  fprintf(file, "%12hd number of boxes in r \n", In_Ptr->nr);
  fprintf(file, "\n");

  fprintf(file,
	  "%12.6lf Total reflection (including direct R)\n",
	  Out_Ptr->Rsp + Out_Ptr->Rd);
  for (i = 1; i <= In_Ptr->num_layers; i++)
    fprintf(file, "%12.6lf Absorbed in layer %12hd\n",
	    Out_Ptr->A_l[i], i);
  fprintf(file, "%12.6lf Total transmission\n", Out_Ptr->Tt);
  fprintf(file, "\n");

  fprintf(file, "Reflectance and Transmission in [cm-2]\n");
  fprintf(file, "Absorption in z-layers in [cm-1]\n");
  fprintf(file, "Absorption in z/r-boxes in [cm-3]\n");

  fprintf(file, "z/r [cm] Layer");
  dr = In_Ptr->dr;
  for (i = 0; i < In_Ptr->nr; i++)
    fprintf(file, "%12.6lf ", i * dr);
  fprintf(file, "\n");

  fprintf(file, "Refl. ");
  for (i = 0; i < In_Ptr->nr; i++)
    fprintf(file, "%12.6lf ", Out_Ptr->Rd_r[i]);
  fprintf(file, "\n");

  dz = In_Ptr->dz;
  for (i = 0; i < In_Ptr->nz; i++) {
    fprintf(file, "%12.6lf %12.6lf", i * dz, Out_Ptr->A_z[i]);
    for (j = 0; j < In_Ptr->nr; j++)
      fprintf(file, "%12.6lf ", Out_Ptr->A_rz[j][i]);
    fprintf(file, "\n");
  }

  fprintf(file, "Transm.");
  for (i = 0; i < In_Ptr->nr; i++)
    fprintf(file, "%12.6lf ", Out_Ptr->Tt_r[i]);
  fprintf(file, "\n");

  fclose(file);
}

/****************************************************************
 ****/
void 
ShowOutMenu(char *in_fname)
{
  printf("I   = Input parameters of mcml\n");
  printf("3   = reflectance, absorption, and transmittance\n");

  printf("AL  = absorption vs layer [-]\n");
  printf("Az  = absorption vs z [1/cm]\n");
  printf("Arz = absorption vs r & z [1/cm3]\n");
  printf("Fz  = fluence vs z [-]\n");
  printf("Frz = fluence vs r & z [1/cm2]\n");

  printf("Rr  = diffuse reflectance vs radius r [1/cm2]\n");
  printf("Ra  = diffuse reflectance vs angle alpha [1/sr]\n");
  printf("Rra = diffuse reflectance vs radius and angle [1/(cm2 sr)]\n");

  printf("Tr  = transmittance vs radius r [1/cm2]\n");
  printf("Ta  = transmittance vs angle alpha [1/sr]\n");
  printf("Tra = transmittance vs radius and angle [1/(cm2 sr)]\n");

  printf("K   = Keijzer's format\n");
  printf("Q   = Quit to main menu\n");
  printf("* input filename: %s \n", in_fname);
}

/****************************************************************
 ****/
void 
BranchOutA(char *Cmd_Str,
	   InputStruct * In_Ptr,
	   OutStruct * Out_Ptr)
{
  switch (toupper(Cmd_Str[1])) {
    case 'L':			/* A_l. */
    WriteA_layer(In_Ptr->num_layers, Out_Ptr->A_l);
    break;
  case 'Z':
    if (toupper(Cmd_Str[2]) == '\0')	/* A_z. */
      WriteA_z(In_Ptr, Out_Ptr->A_z);
    else
      puts("...Wrong command");
    break;
  case 'R':
    if (toupper(Cmd_Str[2]) == 'Z')	/* A_rz. */
      WriteA_rz(In_Ptr, Out_Ptr->A_rz);
    else
      puts("...Wrong command");
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void 
BranchOutF(char *Cmd_Str,
	   InputStruct * In_Ptr,
	   OutStruct * Out_Ptr)
{
  switch (toupper(Cmd_Str[1])) {
    case 'Z':
    if (toupper(Cmd_Str[2]) == '\0')	/* F_z. */
      WriteF_z(In_Ptr, Out_Ptr->A_z);
    else
      puts("...Wrong command");
    break;
  case 'R':
    if (toupper(Cmd_Str[2]) == 'Z')	/* F_rz. */
      WriteF_rz(In_Ptr, Out_Ptr->A_rz);
    else
      puts("...Wrong command");
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void 
BranchOutR(char *Cmd_Str,
	   InputStruct * In_Ptr,
	   OutStruct * Out_Ptr)
{
  char        ch;

  switch (toupper(Cmd_Str[1])) {
  case 'A':			/* Rd_a. */
    WriteRd_a(In_Ptr, Out_Ptr->Rd_a);
    break;
  case 'R':
    ch = toupper(Cmd_Str[2]);
    if (ch == '\0')		/* Rd_r. */
      WriteRd_r(In_Ptr, Out_Ptr->Rd_r);
    else if (ch == 'A')		/* Rd_ra. */
      WriteRd_ra(In_Ptr, Out_Ptr->Rd_ra);
    else
      puts("...Wrong command");
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void 
BranchOutT(char *Cmd_Str,
	   InputStruct * In_Ptr,
	   OutStruct * Out_Ptr)
{
  char        ch;

  switch (toupper(Cmd_Str[1])) {
  case 'A':			/* Tt_a. */
    WriteTt_a(In_Ptr, Out_Ptr->Tt_a);
    break;
  case 'R':
    ch = toupper(Cmd_Str[2]);
    if (ch == '\0')		/* Tt_r. */
      WriteTt_r(In_Ptr, Out_Ptr->Tt_r);
    else if (ch == 'A')		/* Tt_ra. */
      WriteTt_ra(In_Ptr, Out_Ptr->Tt_ra);
    else
      puts("...Wrong command");
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void 
BranchOutCmd(char *Cmd_Str,
	     InputStruct * In_Ptr,
	     OutStruct * Out_Ptr)
{
  char        ch;

  switch (toupper(Cmd_Str[0])) {
  case 'I':
    WriteInParm(In_Ptr);
    break;
  case '3':
    WriteRAT(Out_Ptr);
    break;
  case 'K':
    WriteKFormat(In_Ptr, Out_Ptr);
    break;
  case 'A':
    BranchOutA(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'F':
    BranchOutF(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'R':
    BranchOutR(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'T':
    BranchOutT(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'H':
    ShowOutMenu(In_Ptr->in_fname);
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
OutputOrigData(InputStruct * In_Ptr,
	       OutStruct * Out_Ptr)
{
  char        cmd_str[STRLEN];

  if (!Out_Ptr->allocated)
    puts("...No data to output");
  else
    do {
      printf("\n> Output mcml data (h for help) => ");
      do
	gets(cmd_str);
      while (!strlen(cmd_str));	/* avoid null string. */
      BranchOutCmd(cmd_str, In_Ptr, Out_Ptr);
    } while (toupper(cmd_str[0]) != 'Q');
}

/****************************Contours***************************/
/****************************************************************
 ****/
void 
ShowContOrigMenu(char *in_fname)
{
  printf("A = absorption vs r & z [1/cm3]\n");
  printf("F = fluence vs r & z [1/cm2]\n");
  printf("R = diffuse reflectance vs radius and angle [1/(cm2 sr)]\n");
  printf("T = transmittance vs radius and angle [1/(cm2 sr)]\n");
  printf("Q   = Quit to main menu\n");
  printf("* input filename: %s \n", in_fname);
}

/****************************************************************
 *	Absorption density to fluence. A = F/mua;
 ****/
void 
A2F(InputStruct * In_Ptr, double **A_rz)
{
  short       nz = In_Ptr->nz, nr = In_Ptr->nr;
  double      dr = In_Ptr->dr, dz = In_Ptr->dz;
  short       ir, iz;
  double      mua;

  for (ir = 0; ir < nr; ir++)
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
F2A(InputStruct * In_Ptr, double **A_rz)
{
  short       nz = In_Ptr->nz, nr = In_Ptr->nr;
  double      dr = In_Ptr->dr, dz = In_Ptr->dz;
  short       ir, iz;
  double      mua;

  for (ir = 0; ir < nr; ir++)
    for (iz = 0; iz < nz; iz++) {
      mua = In_Ptr->layerspecs[IzToLayer(iz, In_Ptr)].mua;
      if (mua > 0.0)
	A_rz[ir][iz] *= mua;
    }
}

/****************************************************************
 ****/
void 
BranchContOrigCmd(char *Cmd_Str,
		  InputStruct * In_Ptr,
		  OutStruct * Out_Ptr)
{
  char        ch;

  switch (toupper(Cmd_Str[0])) {
  case 'A':
    IsoPlot(Out_Ptr->A_rz, In_Ptr->nr - 1, In_Ptr->nz - 1,
	    In_Ptr->dr, In_Ptr->dz);
    break;
  case 'F':
    A2F(In_Ptr, Out_Ptr->A_rz);
    IsoPlot(Out_Ptr->A_rz, In_Ptr->nr - 1, In_Ptr->nz - 1,
	    In_Ptr->dr, In_Ptr->dz);
    F2A(In_Ptr, Out_Ptr->A_rz);
    break;
  case 'R':
    IsoPlot(Out_Ptr->Rd_ra, In_Ptr->nr - 1, In_Ptr->na - 1,
	    In_Ptr->dr, In_Ptr->da);
    break;
  case 'T':
    IsoPlot(Out_Ptr->Tt_ra, In_Ptr->nr - 1, In_Ptr->na - 1,
	    In_Ptr->dr, In_Ptr->da);
    break;
  case 'H':
    ShowContOrigMenu(In_Ptr->in_fname);
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
ContourOrigData(InputStruct * In_Ptr,
		OutStruct * Out_Ptr)
{
  char        cmd_str[STRLEN];

  if (!Out_Ptr->allocated)
    puts("...No data to output");
  else
    do {
      printf("\n> Contour output of mcml data (h for help) => ");
      do
	gets(cmd_str);
      while (!strlen(cmd_str));	/* avoid null string. */
      BranchContOrigCmd(cmd_str, In_Ptr, Out_Ptr);
    } while (toupper(cmd_str[0]) != 'Q');
}

/****************************Scanning***************************/
/****************************************************************
 ****/
void 
ShowScanOrigMenu(char *in_fname)
{
  printf("Ar = absorption vs r @ fixed z [1/cm3]\n");
  printf("Az = absorption vs z @ fixed r [1/cm3]\n");
  printf("Fr = fluence vs r @ fixed z [1/cm2]\n");
  printf("Fz = fluence vs z @ fixed r [1/cm2]\n");
  printf("Rr = diffuse reflectance vs r @ fixed angle [1/(cm2 sr)]\n");
  printf("Ra = diffuse reflectance vs angle @ fixed r [1/(cm2 sr)]\n");
  printf("Tr = transmittance vs r @ fixed angle [1/(cm2 sr)]\n");
  printf("Ta = transmittance vs angle @ fixed r [1/(cm2 sr)]\n");
  printf("Q  = quit\n");
  printf("* input filename: %s \n", in_fname);
}

/****************************************************************
 *	Ext is either "Ars" or "Frs".
 ****/
void 
ScanOrigA_r(char *Ext, InputStruct * In_Ptr, double **A_rz)
{
  short       ir, iz, nr = In_Ptr->nr, nz = In_Ptr->nz;
  double      r, z, dr = In_Ptr->dr, dz = In_Ptr->dz;
  FILE       *file;

  file = GetWriteFile(Ext);
  if (file == NULL)
    return;

  printf("z grid separation is %-10.4lg cm.\n", dz);
  printf("Input fixed z index (0 - %2hd): ", nz - 1);
  iz = GetShort(0, nz - 1);
  fprintf(file, "%-12s\t%-s@z=%-9.3lg\n", "r[cm]", Ext, dz * (iz + 0.5));
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    fprintf(file, "%-12.4E\t%-12.4E\n", r, A_rz[ir][iz]);
  }

  fclose(file);
}

/****************************************************************
 *	Ext is either "Azs" or "Fzs".
 ****/
void 
ScanOrigA_z(char *Ext, InputStruct * In_Ptr, double **A_rz)
{
  short       ir, iz, nr = In_Ptr->nr, nz = In_Ptr->nz;
  double      r, z, dr = In_Ptr->dr, dz = In_Ptr->dz;
  FILE       *file;

  file = GetWriteFile(Ext);
  if (file == NULL)
    return;

  printf("r grid separation is %-10.4lg cm.\n", dr);
  printf("Input fixed r index (0 - %2hd): ", nr - 1);
  ir = GetShort(0, nr - 1);
  fprintf(file, "%-12s\t%-s@r=%-9.3lg\n", "z[cm]", Ext, dr * (ir + 0.5));
  for (iz = 0; iz < nz; iz++) {
    z = (iz + 0.5) * dz;
    fprintf(file, "%-12.4E\t%-12.4E\n", z, A_rz[ir][iz]);
  }

  fclose(file);
}

/****************************************************************
 ****/
void 
ScanOrigRd_r(InputStruct * In_Ptr, double **Rd_ra)
{
  short       ir, ia, nr = In_Ptr->nr, na = In_Ptr->na;
  double      r, a, dr = In_Ptr->dr, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Rrs");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  printf("Angle grid separation is %-10.4lg rad.\n", da);
  printf("Input fixed angle index (0 - %2hd): ", na - 1);
  ia = GetShort(0, na - 1);
  fprintf(file, "%-12s\t%-s@a=%-9.3lg\n", "r[cm]", fname, da * (ia + 0.5));
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    fprintf(file, "%-12.4E\t%-12.4E\n", r, Rd_ra[ir][ia]);
  }

  fclose(file);
}

/****************************************************************
 ****/
void 
ScanOrigRd_a(InputStruct * In_Ptr, double **Rd_ra)
{
  short       ir, ia, nr = In_Ptr->nr, na = In_Ptr->na;
  double      r, a, dr = In_Ptr->dr, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Ras");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  printf("r grid separation is %-10.4lg cm.\n", dr);
  printf("Input fixed r index (0 - %2hd): ", nr - 1);
  ir = GetShort(0, nr - 1);
  fprintf(file, "%-12s\t%-s@r=%-9.3lg\n", "a[rad]", fname, dr * (ir + 0.5));
  for (ia = 0; ia < na; ia++) {
    a = (ia + 0.5) * da;
    fprintf(file, "%-12.4E\t%-12.4E\n", a, Rd_ra[ir][ia]);
  }

  fclose(file);
}

/****************************************************************
 ****/
void 
ScanOrigTt_r(InputStruct * In_Ptr, double **Tt_ra)
{
  short       ir, ia, nr = In_Ptr->nr, na = In_Ptr->na;
  double      r, a, dr = In_Ptr->dr, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Trs");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  printf("Angle grid separation is %-10.4lg rad.\n", da);
  printf("Input fixed angle index (0 - %2hd): ", na - 1);
  ia = GetShort(0, na - 1);
  fprintf(file, "%-12s\t%-s@a=%-9.3lg\n", "r[cm]", fname, da * (ia + 0.5));
  for (ir = 0; ir < nr; ir++) {
    r = (ir + 0.5) * dr;
    fprintf(file, "%-12.4E\t%-12.4E\n", r, Tt_ra[ir][ia]);
  }

  fclose(file);
}

/****************************************************************
 ****/
void 
ScanOrigTt_a(InputStruct * In_Ptr, double **Tt_ra)
{
  short       ir, ia, nr = In_Ptr->nr, na = In_Ptr->na;
  double      r, a, dr = In_Ptr->dr, da = In_Ptr->da;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "Tas");
  file = GetWriteFile(fname);
  if (file == NULL)
    return;

  printf("r grid separation is %-10.4lg cm.\n", dr);
  printf("Input fixed r index (0 - %2hd): ", nr - 1);
  ir = GetShort(0, nr - 1);
  fprintf(file, "%-12s\t%-s@r=%-9.3lg\n", "a[rad]", fname, dr * (ir + 0.5));
  for (ia = 0; ia < na; ia++) {
    a = (ia + 0.5) * da;
    fprintf(file, "%-12.4E\t%-12.4E\n", a, Tt_ra[ir][ia]);
  }

  fclose(file);
}

/****************************************************************
 ****/
void 
BranchScanOrigA(char *Cmd_Str,
		InputStruct * In_Ptr,
		OutStruct * Out_Ptr)
{
  char        fname[STRLEN];

  switch (toupper(Cmd_Str[1])) {
  case 'R':
    strcpy(fname, "Ars");
    ScanOrigA_r(fname, In_Ptr, Out_Ptr->A_rz);
    break;
  case 'Z':
    strcpy(fname, "Azs");
    ScanOrigA_z(fname, In_Ptr, Out_Ptr->A_rz);
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void 
BranchScanOrigF(char *Cmd_Str,
		InputStruct * In_Ptr,
		OutStruct * Out_Ptr)
{
  char        fname[STRLEN];

  A2F(In_Ptr, Out_Ptr->A_rz);

  switch (toupper(Cmd_Str[1])) {
  case 'R':
    strcpy(fname, "Frs");
    ScanOrigA_r(fname, In_Ptr, Out_Ptr->A_rz);
    break;
  case 'Z':
    strcpy(fname, "Fzs");
    ScanOrigA_z(fname, In_Ptr, Out_Ptr->A_rz);
    break;
  default:
    puts("...Wrong command");
  }

  F2A(In_Ptr, Out_Ptr->A_rz);
}

/****************************************************************
 ****/
void 
BranchScanOrigR(char *Cmd_Str,
		InputStruct * In_Ptr,
		OutStruct * Out_Ptr)
{
  switch (toupper(Cmd_Str[1])) {
    case 'R':
    ScanOrigRd_r(In_Ptr, Out_Ptr->Rd_ra);
    break;
  case 'A':
    ScanOrigRd_a(In_Ptr, Out_Ptr->Rd_ra);
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void 
BranchScanOrigT(char *Cmd_Str,
		InputStruct * In_Ptr,
		OutStruct * Out_Ptr)
{
  switch (toupper(Cmd_Str[1])) {
    case 'R':
    ScanOrigTt_r(In_Ptr, Out_Ptr->Tt_ra);
    break;
  case 'A':
    ScanOrigTt_a(In_Ptr, Out_Ptr->Tt_ra);
    break;
  default:
    puts("...Wrong command");
  }
}

/****************************************************************
 ****/
void 
BranchScanOrigCmd(char *Cmd_Str,
		  InputStruct * In_Ptr,
		  OutStruct * Out_Ptr)
{
  char        ch;

  switch (toupper(Cmd_Str[0])) {
  case 'A':
    BranchScanOrigA(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'F':
    BranchScanOrigF(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'R':
    BranchScanOrigR(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'T':
    BranchScanOrigT(Cmd_Str, In_Ptr, Out_Ptr);
    break;
  case 'H':
    ShowScanOrigMenu(In_Ptr->in_fname);
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
ScanOrigData(InputStruct * In_Ptr,
	     OutStruct * Out_Ptr)
{
  char        cmd_str[STRLEN];

  if (!Out_Ptr->allocated)
    puts("...No data to output");
  else
    do {
      printf("\n> Scans of mcml data (h for help) => ");
      do
	gets(cmd_str);
      while (!strlen(cmd_str));	/* avoid null string. */
      BranchScanOrigCmd(cmd_str, In_Ptr, Out_Ptr);
    } while (toupper(cmd_str[0]) != 'Q');
}
