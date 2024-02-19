/**************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *	Program name: MCML.
 *	Input/output of data.
 ****/

#include "mcml.h"

/**************************************************************************
 *	Structure used to check against duplicated file names.
 ****/
struct NameList {
  char        name[STRLEN];
  struct NameList *next;
};

typedef struct NameList NameNode;
typedef NameNode *NameLink;


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
 ****/
double     *
AllocArray1D(short nl, short nh)
{
  double     *v;
  short       i;

  v = (double *) malloc((unsigned) (nh - nl + 1) * sizeof(double));
  if (!v)
    ErrorExit("allocation failure in AllocArray1D()");

  v = v - nl;
  for (i = nl; i <= nh; i++)
    v[i] = 0.0;			/* init. */
  return v;
}

/**************************************************************************
 *	Allocate a matrix with row index from nrl to nrh
 *	inclusive, and column index from ncl to nch
 *	inclusive.
 ****/
double    **
AllocArray2D(short nrl, short nrh,
	     short ncl, short nch)
{
  short       i, j;
  double    **m;

  m = (double **) malloc((unsigned) (nrh - nrl + 1)
			 * sizeof(double *));
  if (!m)
    ErrorExit("allocation failure 1 in AllocArray2D()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (double *) malloc((unsigned) (nch - ncl + 1) * sizeof(double));
    if (!m[i])
      ErrorExit("allocation failure 2 in AllocArray2D()");
    m[i] -= ncl;
  }

  for (i = nrl; i <= nrh; i++)
    for (j = ncl; j <= nch; j++)
      m[i][j] = 0.0;
  return m;
}

/**************************************************************************
 *	Allocate a 3D array with row index from nrl to nrh
 *	inclusive, column index from ncl to nch
 *	inclusive, and depth index from ndl to ndh
 *	inclusive.
 ****/
double   ***
AllocArray3D(short nrl, short nrh,
	     short ncl, short nch,
	     short ndl, short ndh)
{
  short       i;
  double   ***m;

  m = (double ***) malloc((unsigned) (nrh - nrl + 1)
			  * sizeof(double **));
  if (!m)
    ErrorExit("allocation failure 1 in AllocArray3D()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
    m[i] = AllocArray2D(ncl, nch, ndl, ndh);

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
 *	Center a string according to the column width.
 ****/
#define COLWIDTH 80
void
CtrPuts(char *InStr)
{
  short       nspaces;		/* number of spaces to be left-filled. */
  char        outstr[STRLEN];

  nspaces = (COLWIDTH - strlen(InStr)) / 2;
  if (nspaces < 0)
    nspaces = 0;

  strcpy(outstr, "");
  while (nspaces--)
    strcat(outstr, " ");

  strcat(outstr, InStr);

  puts(outstr);
}
#undef COLWIDTH

/**************************************************************************
 *	Print messages about MCML.
 ****/
void
AboutMCML(void)
{
  CtrPuts("MCML 2.0, Copyright (c) 1992-1996, Monte Carlo Simulation of");
  CtrPuts("Light Transport in Multi-Layered Turbid Media");

  CtrPuts(" ");
  CtrPuts("Lihong Wang, Ph.D.");
  CtrPuts("Bioengineering Program, Texas A&M University");
  CtrPuts("College Station, Texas 77843-3120, USA");

  CtrPuts("Liqiong Zheng, B.S.");
  CtrPuts("Summer student from Dept. of Computer Science,");
  CtrPuts("University of Houston, Texas, USA.");

  CtrPuts("Steven L. Jacques, Ph.D.");
  CtrPuts("Oregon Medical Laser Center, Providence/St. Vincent Hospital");
  CtrPuts("9205 SW Barnes Rd., Portland, Oregon 97225, USA");

  CtrPuts(" ");
  CtrPuts("Obtain the program thru anonymous ftp to laser.mda.uth.tmc.edu");

  CtrPuts(" ");
  CtrPuts("Please cite the following article in your publications:");
  printf("\tL.-H. Wang, S. L. Jacques, and L.-Q. Zheng, MCML - Monte \n");
  printf("\tCarlo modeling of photon transport in multi-layered\n");
  printf("\ttissues, Computer Methods and Programs in Biomedicine, 47,\n");
  printf("\t131-146 (1995)\n");
}

/**************************************************************************
 *	Kill the ith char (counting from 0), push the following
 *	chars forward by one.
 ****/
void
KillChar(size_t i, char *Str)
{
  size_t      sl = strlen(Str);

  for (; i < sl; i++)
    Str[i] = Str[i + 1];
}

/**************************************************************************
 *	Eliminate the chars in a string which are not printing
 *	chars or spaces.
 *
 *	Spaces include ' ', '\f', '\t' etc.
 *
 *	Return 1 if no nonprinting chars found, otherwise
 *	return 0.
 ****/
Boolean
CheckCharQ(char *Str)
{
  Boolean     found = 0;	/* found bad char. */
  size_t      sl = strlen(Str);
  size_t      i = 0;

  while (i < sl)
    if (Str[i] < 0 || Str[i] > 255) {
      printf("Non-ASCII file\n");
      return (0);
    } else if (isprint(Str[i]) || isspace(Str[i]))
      i++;
    else {
      found = 1;
      KillChar(i, Str);
      sl--;
    }

  return (found);
}

/**************************************************************************
 *	Return 1 if this line is a comment line in which the
 *	first non-space character is "#", or a space line.
 *	Return 0 otherwise.
 ****/
Boolean
CommentLineQ(char *Buf)
{
  size_t      spn, cspn;

  spn = strspn(Buf, " \t");	/* length spanned by space or tab chars. */

  cspn = strcspn(Buf, "#\n");	/* length before the 1st # or return. */

  if (spn == cspn)		/* comment line or space line. */
    return (1);
  else				/* the line has data. */
    return (0);
}

/**************************************************************************
 *	Skip space or comment lines and return a data line.
 ****/
char       *
FindDataLine(FILE * Fp)
{
  char        buf[STRLEN];

  buf[0] = '\0';
  do {				/* skip space or comment lines. */
    if (fgets(buf, 255, Fp) == NULL) {
      printf("Incomplete data.\n");
      buf[0] = '\0';
      break;
    } else
      CheckCharQ(buf);
  } while (CommentLineQ(buf));

  return (buf);
}

/**************************************************************************
 *	Check whether the file version is the same as Version.
 ****/
Boolean
CheckFileVersionQ(FILE * Fp, char *Version)
{
  char        buf[STRLEN];	/* line buffer. */

  buf[0] = '\0';

  /* skip comment line. */
  do {
    if (fgets(buf, 255, Fp) == NULL) {
      buf[0] = '\0';
      break;
    }
  } while (CommentLineQ(buf));

  if ((buf[0] == '\0') || (strstr(buf, Version) == NULL)) {
    puts("Wrong file version.");
    return (0);
  } else
    return (1);
}

/**************************************************************************
 *      Get a filename and open it for reading, retry until
 *      the file can be opened with a correct version or a '.' is typed.
 *	Return a NULL pointer if '.' is typed.
 ****/
FILE       *
GetFile(char *Fname, char *Version)
{
  FILE       *Fp = NULL;

  while (1) {
    /* prompt. */
    printf("Specify filename (or . to quit to main menu):");
    gets(Fname);

    /* terminate with a period. */
    if (strlen(Fname) == 1 &&  Fname[0] == '.')
      return (NULL);		/* return a NULL pointer if '.' entered. */

    /* open the file & check the version. */
    if ((Fp = fopen(Fname, "r")) == NULL)
      puts("File does not exist.");	/* cannot open the file. */
    else {			
      if (CheckFileVersionQ(Fp, Version))
	return (Fp);
      else
	fclose(Fp);
    }
  }
}

/*******************************************************************************
 *  Find number of media in the list.  At the same time, check the
 *  optical parameters.
 ****/
Boolean
FindNumMediaQ(FILE * Fp, short *NumMediaP)
{
  char        buf[STRLEN], name[STRLEN];
  short       num_media = 0;
  double      n, mua, mus, g;

  while (1) {
    strcpy(buf, FindDataLine(Fp));

    if (buf[0] == '\0') {
      printf("Missing end.\n");
      return (0);
    } else if (strstr(buf, "end") != NULL)
      break;
    else {
      num_media++;
      sscanf(buf, "%s%lf%lf%lf%lf", name, &n, &mua, &mus, &g);
      if (n <= 0 || mua < 0 || mus < 0 || g < -1 || g > 1) {
	printf("Bad optical parameters in %s\n", name);
	return (0);
      }
    }
  }

  *NumMediaP = num_media;
  return (1);
}

/*******************************************************************************
 *  Read the parameters of one medium, assumming the
 *  parameters have been checked with FindNumMediaQ().
 ****/
Boolean
ReadOneMediumQ(FILE * Fp, LayerStru * MediumP)
{
  char        buf[STRLEN];

  strcpy(buf, FindDataLine(Fp));
  if (buf[0] == '\0') {
    printf("Shouldn't happen here!");
    return (0);
  }
  sscanf(buf, "%s%lf%lf%lf%lf",
	 MediumP->medium, &MediumP->n, &MediumP->mua,
	 &MediumP->mus, &MediumP->g);

  return (1);
}

/*******************************************************************************
 *  Read the media list.
 ****/
Boolean
ReadMediumListQ(FILE * Fp, InStru * In_Ptr)
{
  long        file_pos;
  short       i;

  file_pos = ftell(Fp);
  if (!FindNumMediaQ(Fp, &(In_Ptr->num_media)))
    return (0);
  fseek(Fp, file_pos, SEEK_SET);

  /* allocate an array for the layer parameters. */
  In_Ptr->mediumlist = (LayerStru *)
    malloc((In_Ptr->num_media) * sizeof(LayerStru));
  if (!(In_Ptr->mediumlist)) {
    printf("allocation failure in ReadMediumListQ()");
    return (0);
  }
  for (i = 0; i < In_Ptr->num_media; i++)
    ReadOneMediumQ(Fp, &(In_Ptr->mediumlist[i]));
  FindDataLine(Fp);		/* skip the signal end. */

  return (1);
}

/**************************************************************************
 *	Read the file name and the file format.
 *
 *	The file format can be either A for ASCII or B for binary.
 ****/
Boolean
ReadFnameFormatQ(FILE * Fp, InStru * In_Ptr)
{
  char        buf[STRLEN];

  strcpy(buf, FindDataLine(Fp));

  if (sscanf(buf, "%s %c", In_Ptr->out_fname,
	     &(In_Ptr->out_fformat)) != 2) {
    printf("Reading file name and format.\n");
    return (0);
  }
  /* if (toupper(In_Ptr->out_fformat) != 'B') */
  In_Ptr->out_fformat = 'A';	/* now only support 'A' format. */

  return (1);
}

/**************************************************************************
 *	Read the InStru members dz, dr and dt.
 ****/
Boolean
ReadDzDrDtQ(FILE * Fp, InStru * In_Ptr)
{
  char        buf[STRLEN];

  strcpy(buf, FindDataLine(Fp));

  if (sscanf(buf, "%lf%lf%lf", &In_Ptr->dz, &In_Ptr->dr,
	     &In_Ptr->dt) != 3) {
    printf("Reading dz, dr, dt. \n");
    return (0);
  }
  if (In_Ptr->dz <= 0) {
    printf("Nonpositive dz.\n");
    return (0);
  }
  if (In_Ptr->dr <= 0) {
    printf("Nonpositive dr.\n");
    return (0);
  }
  if (In_Ptr->dt <= 0) {
    printf("Nonpositve dt. \n");
    return (0);
  }
  return (1);
}

/**************************************************************************
 *	Read the InStru members nz, nr, nt, na.
 ****/
Boolean
ReadNzNrNtNaQ(FILE * Fp, InStru * In_Ptr)
{
  char        buf[STRLEN];
  float       fz, fr, fa, ft;

  /** read in number of dz, dr, da, dt. **/
  strcpy(buf, FindDataLine(Fp));

  if (sscanf(buf, "%f%f%f%f", &fz, &fr, &ft, &fa) != 4) {
    printf("Reading number of dz, dr, dt, da's.\n");
    return (0);
  }
  if (fz <= 0) {
    printf("Nonpositive number of dz's.\n");
    return (0);
  }
  if (fr <= 0) {
    printf("Nonpositive number of dr's.\n");
    return (0);
  }
  if (fa <= 0) {
    printf("Nonpositive number of da's.\n");
    return (0);
  }
  if (ft <= 0) {
    printf("Nonpositive number of dt's.\n");
    return (0);
  }
  In_Ptr->nz = fz;
  In_Ptr->nr = fr;
  In_Ptr->nt = ft;
  In_Ptr->na = fa;
  In_Ptr->da = 0.5 * PI / In_Ptr->na;

  return (1);
}

/**************************************************************************
 *   Initialize the RecordStru.
 ****/
void
InitRecord(InStru * In_Ptr)
{
  In_Ptr->record.Rd_r = 0;
  In_Ptr->record.Rd_a = 0;
  In_Ptr->record.Rd_ra = 0;
  In_Ptr->record.Rd_t = 0;
  In_Ptr->record.Rd_rt = 0;
  In_Ptr->record.Rd_at = 0;
  In_Ptr->record.Rd_rat = 0;
  In_Ptr->record.Td_r = 0;
  In_Ptr->record.Td_a = 0;
  In_Ptr->record.Td_ra = 0;
  In_Ptr->record.Td_t = 0;
  In_Ptr->record.Td_rt = 0;
  In_Ptr->record.Td_at = 0;
  In_Ptr->record.Td_rat = 0;
  In_Ptr->record.A_z = 0;
  In_Ptr->record.A_rz = 0;
  In_Ptr->record.A_t = 0;
  In_Ptr->record.A_zt = 0;
  In_Ptr->record.A_rzt = 0;
}

/**************************************************************************
 *	Change all characters in a string to upper case.
 ****/
char *
ToUpperString(char *string)
{
  int i;

  for (i = 0; i < strlen(string); i++)
    string[i] = toupper(string[i]);

  return string;
}

/**************************************************************************
 *  Read which quantity is to be scored.
 ****/
Boolean
ReadRecordQ(FILE * Fp, InStru * In_Ptr)
{
  char        buf[STRLEN], string[STRLEN];
  char       *index;
  Boolean     error = 0;

  strcpy(buf, FindDataLine(Fp));
  index = buf;

  if (index[0] == '\0') {
    printf("Read scored quantities.\n");
    error = 1;
  }
  while (!error && index[0] != '\n' && index[0] != '#') {
    if (index[0] == '\\') {	/* use '\' to continue in next line. */
      strcpy(buf, FindDataLine(Fp));
      index = buf;
    }
    sscanf(index, "%s", string);
    index = index + strlen(string);
    while (index[0] == ' ' || index[0] == '\t')
      index++;

    if (strcmp(ToUpperString(string), "RD_R") == 0)
      In_Ptr->record.Rd_r = 1;
    else if (strcmp(ToUpperString(string), "RD_A") == 0)
      In_Ptr->record.Rd_a = 1;
    else if (strcmp(ToUpperString(string), "RD_RA") == 0)
      In_Ptr->record.Rd_ra = 1;
    else if (strcmp(ToUpperString(string), "RD_T") == 0)
      In_Ptr->record.Rd_t = 1;
    else if (strcmp(ToUpperString(string), "RD_RT") == 0)
      In_Ptr->record.Rd_rt = 1;
    else if (strcmp(ToUpperString(string), "RD_AT") == 0)
      In_Ptr->record.Rd_at = 1;
    else if (strcmp(ToUpperString(string), "RD_RAT") == 0)
      In_Ptr->record.Rd_rat = 1;
    else if (strcmp(ToUpperString(string), "TD_R") == 0)
      In_Ptr->record.Td_r = 1;
    else if (strcmp(ToUpperString(string), "TD_A") == 0)
      In_Ptr->record.Td_a = 1;
    else if (strcmp(ToUpperString(string), "TD_RA") == 0)
      In_Ptr->record.Td_ra = 1;
    else if (strcmp(ToUpperString(string), "TD_T") == 0)
      In_Ptr->record.Td_t = 1;
    else if (strcmp(ToUpperString(string), "TD_RT") == 0)
      In_Ptr->record.Td_rt = 1;
    else if (strcmp(ToUpperString(string), "TD_AT") == 0)
      In_Ptr->record.Td_at = 1;
    else if (strcmp(ToUpperString(string), "TD_RAT") == 0)
      In_Ptr->record.Td_rat = 1;
    else if (strcmp(ToUpperString(string), "A_Z") == 0)
      In_Ptr->record.A_z = 1;
    else if (strcmp(ToUpperString(string), "A_RZ") == 0)
      In_Ptr->record.A_rz = 1;
    else if (strcmp(ToUpperString(string), "A_T") == 0)
      In_Ptr->record.A_t = 1;
    else if (strcmp(ToUpperString(string), "A_ZT") == 0)
      In_Ptr->record.A_zt = 1;
    else if (strcmp(ToUpperString(string), "A_RZT") == 0)
      In_Ptr->record.A_rzt = 1;
    else {
      printf("Unknown quantity: %s\n", string);
      error = 1;
    }
  }

  return (!error);
}

/**************************************************************************
*   Filter the RecordStru.
****/
Boolean
FilterRecordQ(FILE * Fp, InStru * In_Ptr)
{
  InitRecord(In_Ptr);

  if (!ReadRecordQ(Fp, In_Ptr))
    return (0);

  if (In_Ptr->record.Rd_rat) {
    In_Ptr->record.Rd_ra = 0;
    In_Ptr->record.Rd_rt = 0;
    In_Ptr->record.Rd_at = 0;
    In_Ptr->record.Rd_r = 0;
    In_Ptr->record.Rd_a = 0;
    In_Ptr->record.Rd_t = 0;
  }
  if (In_Ptr->record.Rd_ra) {
    In_Ptr->record.Rd_r = 0;
    In_Ptr->record.Rd_a = 0;
  }
  if (In_Ptr->record.Rd_rt) {
    In_Ptr->record.Rd_r = 0;
    In_Ptr->record.Rd_t = 0;
  }
  if (In_Ptr->record.Rd_at) {
    In_Ptr->record.Rd_a = 0;
    In_Ptr->record.Rd_t = 0;
  }
  if (In_Ptr->record.Td_rat) {
    In_Ptr->record.Td_ra = 0;
    In_Ptr->record.Td_rt = 0;
    In_Ptr->record.Td_at = 0;
    In_Ptr->record.Td_r = 0;
    In_Ptr->record.Td_a = 0;
    In_Ptr->record.Td_t = 0;
  }
  if (In_Ptr->record.Td_ra) {
    In_Ptr->record.Td_r = 0;
    In_Ptr->record.Td_a = 0;
  }
  if (In_Ptr->record.Td_rt) {
    In_Ptr->record.Td_r = 0;
    In_Ptr->record.Td_t = 0;
  }
  if (In_Ptr->record.Td_at) {
    In_Ptr->record.Td_a = 0;
    In_Ptr->record.Td_t = 0;
  }
  if (In_Ptr->record.A_rzt) {
    In_Ptr->record.A_rz = 0;
    In_Ptr->record.A_zt = 0;
    In_Ptr->record.A_z = 0;
    In_Ptr->record.A_t = 0;
  }
  if (In_Ptr->record.A_rz) {
    In_Ptr->record.A_z = 0;
  }
  if (In_Ptr->record.A_zt) {
    In_Ptr->record.A_z = 0;
    In_Ptr->record.A_t = 0;
  }
  if (In_Ptr->record.A_zt) {
    In_Ptr->record.A_z = 0;
    In_Ptr->record.A_t = 0;
  }
  return (1);
}

/**************************************************************************
 *  Read the threshold weight.
 ****/
Boolean
ReadWthQ(FILE * Fp, InStru * In_Ptr)
{
  char        buf[STRLEN];

  strcpy(buf, FindDataLine(Fp));

  if (sscanf(buf, "%lf", &(In_Ptr->Wth)) != 1) {
    printf("Reading threshold weight.\n");
    return (0);
  }
  if (In_Ptr->Wth < 0 || In_Ptr->Wth >= 1.0) {
    printf("Threshold weight out of range (0-1).\n");
    return (0);
  }
  return (1);
}

/**************************************************************************
 *  Read the random number seed.
 ****/
Boolean
ReadRanSeedQ(FILE * Fp, InStru * In_Ptr)
{
  char        buf[STRLEN];

  strcpy(buf, FindDataLine(Fp));

  if (sscanf(buf, "%ld", &(In_Ptr->ran_seed)) != 1) {
    printf("Reading random number seed.\n");
    return (0);
  }
  if (In_Ptr->ran_seed < 0) {
    printf("Nonpositive random number seed.\n");
    return (0);
  }
  return (1);
}

/**************************************************************************
 *  Find number of layers.
 ****/
Boolean
FindNumLayersQ(FILE * Fp, int *NumLayerP)
{
  char        buf[STRLEN], name[STRLEN];
  short       num_layers = 0;
  double      thick;

  while (1) {
    strcpy(buf, FindDataLine(Fp));

    if (buf[0] == '\0') {
      printf("Missing end.\n");
      return (0);
    } else if (strstr(buf, "end") != NULL)
      break;
    else {
      if ((sscanf(buf, "%s %lf", name, &thick) == 2) ||
	  (sscanf(buf, "%s", name) == 1))
	num_layers++;
    }
  }

  *NumLayerP = num_layers - 2;
  return (1);
}

/**************************************************************************
 *  Check whether the medium name is in the media list.
 ****/
Boolean
ValidMediumNameQ(char *NameP, int *Index, InStru * In_Ptr)
{
  short       i;

  for (i = 0; i < In_Ptr->num_media; i++)
    if (!strcmp(NameP, In_Ptr->mediumlist[i].medium)) {
      *Index = i;
      return (1);
    }
  return (0);
}

/**************************************************************************
 *	Read the parameters of all layers.
 ****/
Boolean
ReadLayerSpecsQ(FILE * Fp, InStru * In_Ptr)
{
  char        buf[STRLEN], name[STRLEN];
  short       i = 0;
  double      z = 0.0, thick;	/* z coordinate of the current layer. */
  int         num_layers, index;
  long        file_pos;

  file_pos = ftell(Fp);
  if (!FindNumLayersQ(Fp, &num_layers))
    return (0);
  fseek(Fp, file_pos, SEEK_SET);

  In_Ptr->num_layers = num_layers;
  /* Allocate an array for the layer parameters. */
  /* layer 0 and layer Num_Layers + 1 are for ambient. */
  In_Ptr->layerspecs = (LayerStru *)
    malloc((unsigned) (num_layers + 2) * sizeof(LayerStru));
  if (!(In_Ptr->layerspecs)) {
    printf("allocation failure in ReadLayerSpecsQ()");
    return (0);
  }
  for (i = 0; i <= num_layers + 1; i++) {
    strcpy(buf, FindDataLine(Fp));
    if (i == 0 || i == num_layers + 1) {
      if (sscanf(buf, "%s", name) != 1) {
	printf("  Error in reading ambient layer name.\n");
	return (0);
      }
    } else {
      if (sscanf(buf, "%s%lf", name, &thick) != 2) {
	printf("  Error in ReadLayerSpecsQ().\n");
	return (0);
      } else if (thick <= 0.0) {
	printf("  Nonpositive layer thickness.\n");
	return (0);
      }
    }

    if (!ValidMediumNameQ(name, &index, In_Ptr)) {
      printf("  Invalid medium name. \n");
      return (0);
    }
    strcpy(In_Ptr->layerspecs[i].medium, name);
    In_Ptr->layerspecs[i].n = In_Ptr->mediumlist[index].n;
    In_Ptr->layerspecs[i].mua = In_Ptr->mediumlist[index].mua;
    In_Ptr->layerspecs[i].mus = In_Ptr->mediumlist[index].mus;
    In_Ptr->layerspecs[i].g = In_Ptr->mediumlist[index].g;

    if ((i != 0) && (i != num_layers + 1)) {
      In_Ptr->layerspecs[i].z0 = z;
      z = z + thick;
      In_Ptr->layerspecs[i].z1 = z;
    } else if (i == 0) {
      In_Ptr->layerspecs[i].z0 = 0.0;
      In_Ptr->layerspecs[i].z1 = 0.0;
    } else if (i == (num_layers + 1)) {
      In_Ptr->layerspecs[i].z0 = z;
      In_Ptr->layerspecs[i].z1 = z;
    }
  }

  return (1);
}

/**************************************************************************
 *      Read the number of photons.
 *      Read computation time limit.
 *      Type = 0, read from a .mci input file;
 *      Type = 1, read from a .mco output file.
 ****/
Boolean
ReadNumPhotonsQ(FILE * Fp, InStru * In_Ptr, char Type)
{
  char        buf[STRLEN];
  float       temp;
  int         hours, minutes;

  strcpy(buf, FindDataLine(Fp));

  if (Type == 0) {
    In_Ptr->add_num_photons = 0;
    In_Ptr->add_num_seconds = 0;
  }
  if (sscanf(buf, "%f %d:%d", &temp, &hours, &minutes) == 3) {
    if (((long) temp > 0) && (hours * 3600 + minutes * 60) >= 0) {
      if (Type == 0) {
	In_Ptr->num_photons = (long) temp;
	In_Ptr->num_seconds = hours * 3600 + minutes * 60;
      } else {
	In_Ptr->add_num_photons = (long) temp;
	In_Ptr->add_num_seconds = hours * 3600 + minutes * 60;
      }

      In_Ptr->control_bit = 3;
    } else {
      printf("Nonpositive number of photons or time limit.\n");
      return (0);
    }

  } else if (sscanf(buf, "%d:%d", &hours, &minutes) == 2) {
    if ((hours * 3600 + minutes * 60) >= 0) {
      if (Type == 0)
	In_Ptr->num_seconds = hours * 3600 + minutes * 60;
      else
	In_Ptr->add_num_seconds = hours * 3600 + minutes * 60;

      In_Ptr->control_bit = 2;
    } else {
      printf("Nonpositive time limit.\n");
      return (0);
    }

  } else if (sscanf(buf, "%f", &temp) == 1) {
    if ((long) temp > 0) {
      if (Type == 0)
	In_Ptr->num_photons = (long) temp;
      else
	In_Ptr->add_num_photons = (long) temp;
      In_Ptr->control_bit = 1;
    } else {
      printf("Nonpositive number of photons.\n");
      return (0);
    }

  } else {
    printf("Invalid number of photons or time limit.\n");
    return (0);
  }

  return (1);
}

/**************************************************************************
 *      Read the beam source type (pencil/isotropic).
 ****/
Boolean
ReadSourceTypeQ(FILE * Fp, InStru * In_Ptr)
{
  char        buf[STRLEN], b_type[STRLEN];

  strcpy(buf, FindDataLine(Fp));

  if (sscanf(buf, "%s", b_type) != 1) {
    printf("Reading photon source type. \n");
    return (0);
  }
  if (strcmp(b_type, "pencil") == 0) {
    In_Ptr->source_type = pencil;
    return (1);
  } else if (strcmp(b_type, "isotropic") == 0) {
    In_Ptr->source_type = isotropic;
    return (1);
  } else {
    printf("Unknow photon source type. \n");
    return (0);
  }
}

/**************************************************************************
 *      Compute the index to layer according to the z coordinate.
 *	If the z is on an interface between layers, the returned index 
 *	will point to the upper layer.
 *	Index 0 is the top ambient medium and index num_layers+1 is the
 *	bottom one.
 ****/
Boolean
ZToLayerQ(double z, short *index, InStru * In_Ptr)
{
  short       i = 0;		/* index to layer. */
  short       num_layers = In_Ptr->num_layers;

  if (z < 0.0) {
    printf("Nonpositive z coordinate.\n");
    return (0);
  } else if (z > In_Ptr->layerspecs[num_layers].z1) {
    printf("Source is outside of the last layer. \n");
    return (0);
  } else {
    while (z > In_Ptr->layerspecs[i].z1)
      i++;

    *index = i;
    return (1);
  }
}

/**************************************************************************
 *      Read starting position of photon source.
 ****/
Boolean
ReadStartPQ(FILE * Fp, InStru * In_Ptr)
{
  char        buf[STRLEN], smedium[STRLEN];
  double      sz;
  short       slayer;

  strcpy(buf, FindDataLine(Fp));

  if ((sscanf(buf, "%lf %s", &sz, smedium) == 2)  /* z and medium. */ 
      && (smedium[0] != '#' && smedium[0] != '\n')) {
    if (!ZToLayerQ(sz, &slayer, In_Ptr))
      return (0);

    if (strcmp(In_Ptr->layerspecs[slayer].medium, smedium) != 0) {
      if ((fabs(sz - In_Ptr->layerspecs[slayer].z1) < DBL_EPSILON) &&
	  (strcmp(In_Ptr->layerspecs[slayer + 1].medium, smedium) == 0)) {
	slayer++;
	if (slayer > In_Ptr->num_layers) {
	  puts("Source is outside of the last layer.");
	  return (0);
	}

      } else {
	printf("Medium name and z coordinate do not match.\n");
	return (0);
      }
    }

  } else if (sscanf(buf, "%lf", &(sz)) == 1) {	/* z only. */
    if (!ZToLayerQ(sz, &In_Ptr->slayer, In_Ptr))
      return (0);
    strcpy(smedium, "");

  } else {
    printf("Invalid starting position of photon source.\n");
    return (0);
  }

  if ((In_Ptr->source_type == isotropic) && (sz == 0.0)) {
    printf("Can not put isotropic source in upper ambient medium.\n");
    return 0;
  }

  In_Ptr->sz = sz;
  strcpy(In_Ptr->smedium, smedium);
  return (1);
}

/*************************************************************************
 *	Compute the critical angles for total internal
 *	reflection according to the relative refractive index
 *	of the layer.
 *	All layers are processed.
 ****/
void
CriticalAngle(short Num_Layers,
	      LayerStru ** Layerspecs_PP)
{
  short       i = 0;
  double      n1, n2;

  for (i = 1; i <= Num_Layers; i++) {
    n1 = (*Layerspecs_PP)[i].n;
    n2 = (*Layerspecs_PP)[i - 1].n;
    (*Layerspecs_PP)[i].cos_crit0 = n1 > n2 ?
      sqrt(1.0 - n2 * n2 / (n1 * n1)) : 0.0;

    n2 = (*Layerspecs_PP)[i + 1].n;
    (*Layerspecs_PP)[i].cos_crit1 = n1 > n2 ?
      sqrt(1.0 - n2 * n2 / (n1 * n1)) : 0.0;
  }
}

/**************************************************************************
 *	Read in the input parameters for one run.
 ****/
void
ReadParam(FILE * Fp, InStru * In_Ptr)
{
  if (!ReadFnameFormatQ(Fp, In_Ptr))
    exit(1);

  /* geometry. */
  if (!ReadLayerSpecsQ(Fp, In_Ptr))
    exit(1);
  FindDataLine(Fp);		/* skip the signal "end" of layers. */

  /* source. */
  if (!ReadSourceTypeQ(Fp, In_Ptr))
    exit(1);
  if (!ReadStartPQ(Fp, In_Ptr))
    exit(1);

  /* grids. */
  if (!ReadDzDrDtQ(Fp, In_Ptr))
    exit(1);
  if (!ReadNzNrNtNaQ(Fp, In_Ptr))
    exit(1);
  In_Ptr->zm = In_Ptr->dz * In_Ptr->nz;
  In_Ptr->rm = In_Ptr->dr * In_Ptr->nr;
  In_Ptr->tm = In_Ptr->dt * In_Ptr->nt;
  In_Ptr->am = In_Ptr->da * In_Ptr->na;

  /* scored data categories. */
  if (!FilterRecordQ(Fp, In_Ptr))
    exit(1);

  /* simulation control. */
  if (!ReadNumPhotonsQ(Fp, In_Ptr, 0))
    exit(1);
  if (!ReadWthQ(Fp, In_Ptr))
    exit(1);
  if (!ReadRanSeedQ(Fp, In_Ptr))
    exit(1);

  CriticalAngle(In_Ptr->num_layers, &In_Ptr->layerspecs);
}

/**************************************************************************
 *  Read the media list in interactive mode.
 ****/
void
InterReadMediumList(InStru * In_Ptr)
{
  short       i, j;
  char        string[STRLEN], medium[STRLEN];
  Boolean     name_taken;

  printf("Specify medium list. Total number of mediums: ");

  gets(string);
  while (sscanf(string, "%hd", &In_Ptr->num_media) != 1
	 || In_Ptr->num_media <= 0) {
    printf("Invalid medium number. Input again: ");
    gets(string);
  }

  /* allocate an array for the layer parameters. */
  In_Ptr->mediumlist = (LayerStru *)
    malloc((In_Ptr->num_media) * sizeof(LayerStru));
  if (!(In_Ptr->mediumlist)) {
    printf("allocation failure in ReadMediumListQ()");
    exit(1);
  }
  for (i = 0; i < In_Ptr->num_media; i++) {
    printf("Specify meidum %d: \n  Medium name: ", i + 1);
    gets(string);
    do {
      name_taken = 0;
      sscanf(string, "%s", medium);
      for (j = 0; j < i; j++)
	if (strcmp(In_Ptr->mediumlist[j].medium, medium) == 0) {
	  name_taken = 1;
	  printf("  Duplicated medium. Input again: ");
	  gets(string);
	  break;
	}
    } while (name_taken);
    strcpy(In_Ptr->mediumlist[i].medium, medium);

    printf("  Refractive index n (>= 1.0): ");
    gets(string);
    while (sscanf(string, "%lf", &In_Ptr->mediumlist[i].n) != 1
	   || In_Ptr->mediumlist[i].n < 1.0) {
      printf("  Invalid refractive index. Input again (>= 1.0): ");
      gets(string);
    }

    printf("  Absorption coefficient mua (>= 0.0 /cm): ");
    gets(string);
    while (sscanf(string, "%lf", &In_Ptr->mediumlist[i].mua) != 1
	   || In_Ptr->mediumlist[i].mua < 0.0) {
      printf("  Invalid absorption coefficient. Input again (>= 0.0): ");
      gets(string);
    }

    printf("  Scattering coefficient mus (>= 0.0 /cm): ");
    gets(string);
    while (sscanf(string, "%lf", &In_Ptr->mediumlist[i].mus) != 1
	   || In_Ptr->mediumlist[i].mus < 0.0) {
      printf("  Invalid scattering coefficient. Input again (>= 0.0): ");
      gets(string);
    }

    printf("  Anisotropy factor g (0.0 - 1.0): ");
    gets(string);
    while (sscanf(string, "%lf", &In_Ptr->mediumlist[i].g) != 1
	   || In_Ptr->mediumlist[i].g < 0.0
	   || In_Ptr->mediumlist[i].g > 1.0) {
      printf("  Invalid anisotropy factor. Input again (0.0 - 1.0): ");
      gets(string);
    }
    printf("\n");
  }
}

/**************************************************************************
 *      Read the file name and the file format interactively.
 *
 *      The file format can be either A for ASCII or B for
 *      binary.
 ****/
void
InterReadFnameFormat(InStru * In_Ptr)
{
  FILE       *file;
  char        fname[STRLEN], fmode[STRLEN];

  do {
    printf("Specify output filename with extension .mco: ");
    gets(fname);
    fmode[0] = 'w';

    if ((file = fopen(fname, "r")) != NULL) { /* file exists. */
      printf("File %s exists, %s", fname, "w=overwrite, n=new filename: ");
      do
        gets(fmode);
      while (!strlen(fmode));   /* avoid null line. */
      fclose(file);
    }
  } while (fmode[0] != 'w');

  strcpy(In_Ptr->out_fname, fname);

  /*
   * printf("Output file format (A/B): "); gets(fname); while
   * (sscanf(fname, "%c", &In_Ptr->out_fformat) != 1) { printf("Error
   * occured. Output file format (A/B): "); gets(fname); }
   * 
   * if (toupper(In_Ptr->out_fformat) != 'B')
   */
  In_Ptr->out_fformat = 'A';	/* now only support 'A' format. */

  printf("\n");
}

/**************************************************************************
 *	Read dz, dr, dt interactively.
 ****/
void
InterReadDzDrDt(InStru * In_Ptr)
{
  do {
    printf("Specify dz, dr, dt in one line\n");
    printf("(all > 0.0 cm, e.g., 0.1 0.1 0.1): ");
  } while (!ReadDzDrDtQ(stdin, In_Ptr));

  printf("\n");
}

/**************************************************************************
 *      Read the InStru members nz, nr, na interactively.
 ****/
void
InterReadNzNrNtNa(InStru * In_Ptr)
{
  do {
    printf("Specify nz, nr, nt, na in one line\n");
    printf("(all > 0, e.g., 100 100 100 100): ");
  } while (!ReadNzNrNtNaQ(stdin, In_Ptr));

  In_Ptr->da = 0.5 * PI / In_Ptr->na;
  printf("\n");
}

/**************************************************************************
 *	Read and filter the quantities to be scored interactively.
 ****/
void
InterFilterRecord(InStru * In_Ptr)
{
  do {
    printf("Select scored quantities from the following data categories:\n");
    printf("\tRd_rat\t\t\tTd_rat\t\t\tA_rzt\n");
    printf("\tRd_ra\tRd_rt\tRd_at\tTd_ra\tTd_rt\tRd_at\tA_rz\tA_zt\n");
    printf("\tRd_r\tRd_a\tRd_t\tTd_r\tTd_a\tTd_t\tA_z\tA_t\n");
  } while (!FilterRecordQ(stdin, In_Ptr));

  printf("\n");
}

/**************************************************************************
 *	Read the threshold weight interactively.
 ****/
void
InterReadWth(InStru * In_Ptr)
{
  char        string[STRLEN];

  printf("Input threshold weight (0 <= wth < 1.0, 0.0001 recommended): ");
  gets(string);
  while (sscanf(string, "%lf", &In_Ptr->Wth) != 1
	 || In_Ptr->Wth < 0 || In_Ptr->Wth >= 1) {
    printf("Invalid wth. Input again (0 <= wth < 1.0): ");
    gets(string);
  }

  printf("\n");
}

/**************************************************************************
 *	Read the random seed interactively.
 ****/
void
InterReadRanSeed(InStru * In_Ptr)
{
  char        string[STRLEN];

  printf("Input random number seed (1 <= ran_seed <= 32000): ");
  gets(string);
  while (sscanf(string, "%ld", &In_Ptr->ran_seed) != 1
	 || In_Ptr->ran_seed < 1 || In_Ptr->ran_seed > 32000) {
    printf("Invalid ran_seed. Input again (1 <= ran_seed <= 32000): ");
    gets(string);
  }

  printf("\n");
}

/**************************************************************************
 ****/
void
PrintMediumNames(InStru * In_Ptr)
{
  int         i, j;

  printf("Available medium types:\n");

  j = 1;
  for (i = 0; i < In_Ptr->num_media; i++) {
    printf("%-16s", In_Ptr->mediumlist[i].medium);
    if (j % 4 == 0)
      printf("\n");
    j++;
  }

  printf("\n");
}

/**************************************************************************
 *	Read layer specifications interactively.
 ****/
void
InterReadLayerSpecs(InStru * In_Ptr)
{
  char        string[STRLEN], name[STRLEN];
  short       i = 0;
  double      z = 0.0, thick;	/* z coordinate of the current layer. */
  int         index;
  Boolean     error;

  printf("\nSpecify layer list. ");
  PrintMediumNames(In_Ptr);
  printf("\nTotal number of layers: ");

  gets(string);
  while (sscanf(string, "%hd", &In_Ptr->num_layers) != 1
	 || In_Ptr->num_layers <= 0) {
    printf("Invalid layer number. Input again: ");
    gets(string);
  }

  /* Allocate an array for the layer parameters. */
  /* layer 0 and layer Num_Layers + 1 are for ambient. */
  In_Ptr->layerspecs = (LayerStru *)
    malloc((unsigned) (In_Ptr->num_layers + 2) * sizeof(LayerStru));
  if (!(In_Ptr->layerspecs)) {
    printf("allocation failure in ReadLayerSpecsQ()");
    exit(1);
  }
  for (i = 0; i <= In_Ptr->num_layers + 1; i++) {
    error = 1;
    while (error) {
      error = 0;
      if (i == 0) {
	printf("\n  Name of upper ambient medium: ");
	gets(string);
	sscanf(string, "%s", name);
      } else if (i == In_Ptr->num_layers + 1) {
	printf("\n  Name of lower ambient medium: ");
	gets(string);
	sscanf(string, "%s", name);
      } else {
	printf("\n  Medium name of layer %d: ", i);
	gets(string);
	sscanf(string, "%s", name);
      }

      if (!ValidMediumNameQ(name, &index, In_Ptr)) {
	printf("  Invalid medium name. Input again.");
	error = 1;
      }
    }

    strcpy(In_Ptr->layerspecs[i].medium, name);
    In_Ptr->layerspecs[i].n = In_Ptr->mediumlist[index].n;
    In_Ptr->layerspecs[i].mua = In_Ptr->mediumlist[index].mua;
    In_Ptr->layerspecs[i].mus = In_Ptr->mediumlist[index].mus;
    In_Ptr->layerspecs[i].g = In_Ptr->mediumlist[index].g;

    if ((i != 0) && (i != In_Ptr->num_layers + 1)) {
      printf("  Input the thickness of layer %d (thickness > 0.0 cm): ", i);
      gets(string);
      while (sscanf(string, "%lf", &thick) != 1 || thick <= 0) {
	printf("  Invalid thickness. Input again (thickness > 0.0 cm): ");
	gets(string);
      }
      In_Ptr->layerspecs[i].z0 = z;
      z = z + thick;
      In_Ptr->layerspecs[i].z1 = z;
    } else if (i == 0) {
      In_Ptr->layerspecs[i].z0 = 0.0;
      In_Ptr->layerspecs[i].z1 = 0.0;
    } else if (i == In_Ptr->num_layers + 1) {
      In_Ptr->layerspecs[i].z0 = z;
      In_Ptr->layerspecs[i].z1 = z;
    }
  }

  printf("\n");
}

/**************************************************************************
 *      Read the number of photons, or computation time interactively.
 ****/
void
InterReadNumPhotons(InStru * In_Ptr)
{
  printf("Specify number of photons or time in hh:mm format,\n");
  printf("or both in one line (e.g. 10000 5:30): ");

  while (!ReadNumPhotonsQ(stdin, In_Ptr, 0))
    printf("Input again: ");

  printf("\n");
}

/**************************************************************************
 *      Read the beam source type (pencil/isotropic).
 ****/
void
InterReadSourceType(InStru * In_Ptr)
{
  char        string[STRLEN], c;

  printf("Input source type (p = pencil / i = isotropic): ");
  gets(string);
  while (sscanf(string, "%c", &c) != 1 ||
	 !(toupper(c) == 'P' || toupper(c) == 'I')) {
    printf("Invalid type. Input again (p = pencil / i = isotropic): ");
    gets(string);
  }

  if (toupper(c) == 'P')
    In_Ptr->source_type = pencil;
  else
    In_Ptr->source_type = isotropic;

  printf("\n");
}

/**************************************************************************
 *      Read starting position of photon source.
 ****/
void
InterReadStartP(InStru * In_Ptr)
{
  do {
    printf("Input the z coordinate of source (0.0 - %f cm) and the medium\n",
	   In_Ptr->layerspecs[In_Ptr->num_layers].z1);
    printf("where the source is if the z is on an interface (e.g. 1.0 [air]):");
  } while (!ReadStartPQ(stdin, In_Ptr));

  printf("\n");
}

/*************************************************************************
 *     if Fp is stdout, freeze the screen and print a more
 *     message on screen every 20 lines.
 *     The Line is the line index.
 ****/
void
More(FILE * Fp, int *Line)
{
  char        c;

  if (Fp == stdout)
    if (!((*Line) % 20)) {
      printf("--More-- (Press Return to continue)");
      fflush(Fp);

      do {
	fread(&c, 1, 1, stdin);
      } while (c != '\n');
    }
}

/*************************************************************************
 *     Write medium list to the file Fp.
 *     if Fp is stdout, freeze the screen every 20 lines.
 ****/
void
PutMediumListToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  int         i;
  char        format[STRLEN];

  fprintf(Fp, "# Specify media \n");
  (*Line)++;
  fprintf(Fp, "#\tname\t\tn\tmua\tmus\tg\n");
  (*Line)++;

  for (i = 0; i < In_Ptr->num_media; i++) {
    LayerStru   s;

    More(Fp, Line);
    s = In_Ptr->mediumlist[i];
    if (strlen(s.medium) + 1 > 8)
      strcpy(format, "\t%s \t%G\t%G\t%G\t%G\n");
    else
      strcpy(format, "\t%s \t\t%G\t%G\t%G\t%G\n");

    fprintf(Fp, format, s.medium, s.n, s.mua, s.mus, s.g);
    (*Line)++;
  }
  fprintf(Fp, "end #of media\n");
  (*Line)++;
}

void
PutFnameFormatToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  More(Fp, Line);
  fprintf(Fp,
	  "%s \t%c\t\t\t# output file name, format.\n",
	  In_Ptr->out_fname, In_Ptr->out_fformat);
  (*Line)++;
}

void
PutDzDrDtToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  More(Fp, Line);
  fprintf(Fp, "%G\t%G\t%G\t\t\t# dz, dr, dt.\n",
	  In_Ptr->dz, In_Ptr->dr, In_Ptr->dt);
  (*Line)++;
}

void
PutNzNrNtNaToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  More(Fp, Line);
  fprintf(Fp, "%d\t%d\t%d\t%d\t\t# nz, nr, nt, na.\n",
	  In_Ptr->nz, In_Ptr->nr, In_Ptr->nt, In_Ptr->na);
  (*Line)++;
}

void
PutScoredToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  More(Fp, Line);
  fprintf(Fp, "# This simulation will score the following categories:\n");
  (*Line)++;

  More(Fp, Line);
  if (In_Ptr->record.Rd_r)
    fprintf(Fp, "Rd_r \t");
  if (In_Ptr->record.Rd_a)
    fprintf(Fp, "Rd_a \t");
  if (In_Ptr->record.Rd_ra)
    fprintf(Fp, "Rd_ra \t");
  if (In_Ptr->record.Rd_t)
    fprintf(Fp, "Rd_t \t");
  if (In_Ptr->record.Rd_rt)
    fprintf(Fp, "Rd_rt \t");
  if (In_Ptr->record.Rd_at)
    fprintf(Fp, "Rd_at \t");
  if (In_Ptr->record.Rd_rat)
    fprintf(Fp, "Rd_rat \t");

  if (In_Ptr->record.Td_r)
    fprintf(Fp, "Td_r \t");
  if (In_Ptr->record.Td_a)
    fprintf(Fp, "Td_a \t");
  if (In_Ptr->record.Td_ra)
    fprintf(Fp, "Td_ra \t");
  if (In_Ptr->record.Td_t)
    fprintf(Fp, "Td_t \t");
  if (In_Ptr->record.Td_rt)
    fprintf(Fp, "Td_rt \t");
  if (In_Ptr->record.Td_at)
    fprintf(Fp, "Td_at \t");
  if (In_Ptr->record.Td_rat)
    fprintf(Fp, "Td_rat \t");

  if (In_Ptr->record.A_z)
    fprintf(Fp, "A_z \t");
  if (In_Ptr->record.A_rz)
    fprintf(Fp, "A_rz \t");
  if (In_Ptr->record.A_t)
    fprintf(Fp, "A_t \t");
  if (In_Ptr->record.A_zt)
    fprintf(Fp, "A_zt \t");
  if (In_Ptr->record.A_rzt)
    fprintf(Fp, "A_rzt \t");

  fprintf(Fp, "\n");
  (*Line)++;
}

void
PutWthToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  More(Fp, Line);
  fprintf(Fp, "%G\t\t\t\t\t# threshold weight.\n", In_Ptr->Wth);
  (*Line)++;
}

void
PutRanSeedToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  More(Fp, Line);
  fprintf(Fp, "%ld\t\t\t\t\t# random number seed.\n", In_Ptr->ran_seed);
  (*Line)++;
}

void
PutLayerSpecsToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  int         i;
  char        format[STRLEN];

  More(Fp, Line);
  fprintf(Fp, "# \tmedium \t\tthickness\n");
  (*Line)++;

  for (i = 0; i <= In_Ptr->num_layers + 1; i++) {
    LayerStru   s;
    More(Fp, Line);

    s = In_Ptr->layerspecs[i];
    if (i != 0 && i != In_Ptr->num_layers + 1) {
      if (strlen(s.medium) + 1 > 8)
	strcpy(format, "\t%s \t%G\n");
      else
	strcpy(format, "\t%s \t\t%G\n");
      fprintf(Fp, format, s.medium, s.z1 - s.z0);
    } else
      fprintf(Fp, "\t%s\n", s.medium);
    (*Line)++;
  }

  More(Fp, Line);
  fprintf(Fp, "end #of layers\n");
  (*Line)++;
}

void
PutNumPhotonsToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  More(Fp, Line);

  if (In_Ptr->control_bit == 1)
    fprintf(Fp, "%ld  \t\t\t\t\t# no. of photons | time\n",
	    In_Ptr->num_photons);
  else if (In_Ptr->control_bit == 2)
    fprintf(Fp, "%ld:%ld\t\t\t\t\t# no. of photons | time\n",
	    In_Ptr->num_seconds / 3600, In_Ptr->num_seconds % 3600 / 60);
  else
    fprintf(Fp, "%ld  \t%ld:%ld\t\t\t\t# no. of photons | time\n",
	    In_Ptr->num_photons, In_Ptr->num_seconds / 3600,
	    In_Ptr->num_seconds % 3600 / 60);

  (*Line)++;
}

void
PutSourceTypeToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  More(Fp, Line);

  if (In_Ptr->source_type == pencil)
    fprintf(Fp, "pencil \t\t\t\t\t# src type: pencil/isotropic.\n");
  else
    fprintf(Fp, "isotropic \t\t\t\t# src type: pencil/isotropic.\n");

  (*Line)++;
}

void
PutStartPToFile(FILE * Fp, InStru * In_Ptr, int *Line)
{
  More(Fp, Line);

  if (strlen(In_Ptr->smedium) == 0)
    fprintf(Fp, "%G\t\t\t\t\t# starting position of source.\n",
	    In_Ptr->sz);
  else {
    if (strlen(In_Ptr->smedium) + 1 > 8)
      fprintf(Fp, "%G\t%s \t\t\t# starting position of source.\n",
	      In_Ptr->sz, In_Ptr->smedium);
    else
      fprintf(Fp, "%G\t%s \t\t\t\t# starting position of source.\n",
	      In_Ptr->sz, In_Ptr->smedium);
  }

  (*Line)++;
}

/*************************************************************************
 *     Write input parameters to the file Fp.
 *     if Fp is stdout, freeze the screen every 20 lines.
 ****/
void
PutInputToFile(FILE * Fp, InStru * In_Ptr)
{
  int         line;		/* line index. */

  fprintf(Fp, "mcmli2.0 \t\t\t# file version \n\n");
  line = 2;
  PutMediumListToFile(Fp, In_Ptr, &line);

  More(Fp, &line);
  fprintf(Fp, "\n# Specify data for run 1\n");
  line += 2;

  PutFnameFormatToFile(Fp, In_Ptr, &line);

  More(Fp, &line);	/* geometry. */
  fprintf(Fp, "\n");
  line++;
  PutLayerSpecsToFile(Fp, In_Ptr, &line);

  More(Fp, &line);	/* source. */
  fprintf(Fp, "\n");
  line++;
  PutSourceTypeToFile(Fp, In_Ptr, &line);
  PutStartPToFile(Fp, In_Ptr, &line);

  More(Fp, &line);	/* grids. */
  fprintf(Fp, "\n");
  line++;
  PutDzDrDtToFile(Fp, In_Ptr, &line);
  PutNzNrNtNaToFile(Fp, In_Ptr, &line);

  More(Fp, &line);	/* scored data categories. */
  fprintf(Fp, "\n");
  line++;
  PutScoredToFile(Fp, In_Ptr, &line);

  More(Fp, &line);	/* simulation control. */
  fprintf(Fp, "\n");
  line++;
  PutNumPhotonsToFile(Fp, In_Ptr, &line);
  PutWthToFile(Fp, In_Ptr, &line);
  PutRanSeedToFile(Fp, In_Ptr, &line);

  More(Fp, &line);
  fprintf(Fp, "end #of runs\n\n");
}

/**************************************************************************
 *      Read in the input parameters for one run
 *      in interactive mode.
 ****/
void
InterReadParam(InStru * In_Ptr)
{
  char        string[STRLEN];
  FILE       *fp;

  InterReadMediumList(In_Ptr);

  InterReadFnameFormat(In_Ptr);

  InterReadLayerSpecs(In_Ptr);

  InterReadSourceType(In_Ptr);
  InterReadStartP(In_Ptr);

  InterReadDzDrDt(In_Ptr);
  InterReadNzNrNtNa(In_Ptr);

  InterFilterRecord(In_Ptr);

  InterReadNumPhotons(In_Ptr);
  InterReadWth(In_Ptr);
  InterReadRanSeed(In_Ptr);

  printf("Do you want to save the input to a file? (Y/N)");
  gets(string);
  if (toupper(string[0]) == 'Y') {
    printf("Give the file name to save input: ( .mci): ");
    gets(string);
    if ((fp = fopen(string, "w")) == NULL)
      puts("Can not open the file to write.");
    else
      PutInputToFile(fp, In_Ptr);
  }
}

/**************************************************************************
 *      Check consistance of input parameters for one run.
 *      Such as: the consistance of medium list, layer
 *      list, souce starting position and source type.
 ****/
Boolean
CheckInputConsis(InStru * In_Ptr)
{
  int         i, index;

  for (i = 0; i <= In_Ptr->num_layers + 1; i++) {
    if (!ValidMediumNameQ(In_Ptr->layerspecs[i].medium,
			  &index, In_Ptr)) {
      printf("Invalid medium name of layer %d.\n", i);
      return 0;
    } else {
      In_Ptr->layerspecs[i].n = In_Ptr->mediumlist[index].n;
      In_Ptr->layerspecs[i].mua = In_Ptr->mediumlist[index].mua;
      In_Ptr->layerspecs[i].mus = In_Ptr->mediumlist[index].mus;
      In_Ptr->layerspecs[i].g = In_Ptr->mediumlist[index].g;
    }
  }

  if ((In_Ptr->source_type == isotropic) && (In_Ptr->sz == 0.0)) {
    printf("Can not put isotropic source in upper ambient medium.\n");
    return 0;
  }
  if (!ZToLayerQ(In_Ptr->sz, &In_Ptr->slayer, In_Ptr))
    return 0;

  if (In_Ptr->smedium[0] != '\0') {
    if (strcmp(In_Ptr->layerspecs[In_Ptr->slayer].medium,
	       In_Ptr->smedium) != 0) {
      if ((fabs(In_Ptr->sz - In_Ptr->layerspecs[In_Ptr->slayer].z1)
	   < DBL_EPSILON)
	  && (strcmp(In_Ptr->layerspecs[In_Ptr->slayer + 1].medium,
		     In_Ptr->smedium) == 0))
	In_Ptr->slayer++;
      else {
	printf("Medium name and z coordinate do not match.\n");
	return (0);
      }
    }
  }
  return 1;
}

/****************************************************************************
 *    Menu for changing input parameters.
 *****/
void
ShowChangeMenu()
{
  puts("  o = Print the input on screen.");
  puts("  m = Change media list.");
  puts("  f = Change output file name and format.");
  puts("  d = Change dz, dr, dt.");
  puts("  n = Change nz, nr, nt, na.");
  puts("  c = Change scored data categories.");
  puts("  w = Change threshold weight.");
  puts("  r = Change random number seed.");
  puts("  l = Change layer specifications.");
  puts("  p = Change photon number and computation time limit.");
  puts("  s = Change source type.");
  puts("  z = Change source starting position.");
  puts("  q = Quit from change menu and start simulation.");
  puts("  x = Exit to the main menu.");
  puts("  * Commands here are not case-sensitive");
}

/***********************************************************************
 *   Continue to change input parameters or quit.
 *****/
char
QuitOrContinue()
{
  char        string[STRLEN];

  do {
    printf("Do you want to change them? (Y/N): ");
    do {
      gets(string);
    } while (!strlen(string));
  } while (toupper(string[0]) != 'Y' && toupper(string[0]) != 'N');

  return (toupper(string[0]));
}

void
ChangeMediumList(InStru * In_Ptr)
{
  int         line = 1;

  printf("Current medium list: \n");
  PutMediumListToFile(stdout, In_Ptr, &line);
  printf("\n");

  if (QuitOrContinue() == 'Y') {
    free(In_Ptr->mediumlist);
    InterReadMediumList(In_Ptr);
  }
}

void
ChangeFnameFormat(InStru * In_Ptr)
{
  int         line = 1;

  printf("Current output file name and format: \n");
  PutFnameFormatToFile(stdout, In_Ptr, &line);
  printf("\n");

  InterReadFnameFormat(In_Ptr);
}

void
ChangeDzDrDt(InStru * In_Ptr)
{
  int         line = 1;

  printf("Current dz, dr, dt: \n");
  PutDzDrDtToFile(stdout, In_Ptr, &line);
  printf("\n");

  InterReadDzDrDt(In_Ptr);
}

void
ChangeNzNrNtNa(InStru * In_Ptr)
{
  int         line = 1;

  printf("Current nz, nr, nt, na: \n");
  PutNzNrNtNaToFile(stdout, In_Ptr, &line);
  printf("\n");

  InterReadNzNrNtNa(In_Ptr);
}

void
ChangeRecord(InStru * In_Ptr)
{
  int         line = 1;

  PutScoredToFile(stdout, In_Ptr, &line);
  printf("\n");

  if (QuitOrContinue() == 'Y')
    InterFilterRecord(In_Ptr);
}

void
ChangeWth(InStru * In_Ptr)
{
  int         line = 1;

  printf("Current threshold weight: \n");
  PutWthToFile(stdout, In_Ptr, &line);
  printf("\n");

  InterReadWth(In_Ptr);
}

void
ChangeRanSeed(InStru * In_Ptr)
{
  int         line = 1;

  printf("Current random number seed: \n");
  PutRanSeedToFile(stdout, In_Ptr, &line);
  printf("\n");

  InterReadRanSeed(In_Ptr);
}

void
ChangeLayerSpecs(InStru * In_Ptr)
{
  int         line = 1;

  printf("Current layer sepcification: \n");
  PutLayerSpecsToFile(stdout, In_Ptr, &line);
  printf("\n");

  if (QuitOrContinue() == 'Y')
    InterReadLayerSpecs(In_Ptr);
}

void
ChangeNumPhotons(InStru * In_Ptr)
{
  int         line = 1;

  printf("Current value: \n");
  PutNumPhotonsToFile(stdout, In_Ptr, &line);
  printf("\n");

  InterReadNumPhotons(In_Ptr);
}

void
ChangeSourceType(InStru * In_Ptr)
{
  int         line = 1;

  printf("Current source type: \n");
  PutSourceTypeToFile(stdout, In_Ptr, &line);
  printf("\n");

  InterReadSourceType(In_Ptr);
}

void
ChangeStartP(InStru * In_Ptr)
{
  int         line = 1;

  printf("Layer Specification: \n");
  PutLayerSpecsToFile(stdout, In_Ptr, &line);
  printf("\nCurrent starting position: \n");
  PutStartPToFile(stdout, In_Ptr, &line);
  printf("\n");

  InterReadStartP(In_Ptr);
}

/************************************************************************
 *   return 1 if string[0] = Q, quit change menu;
 *   return 2 if string[0] = X, quit to the main menu;
 *   return 0 otherwise.
 ****/
int
BranchChangeMenu(char *string, InStru * In_Ptr)
{
  switch (toupper(string[0])) {
  case 'M':
    ChangeMediumList(In_Ptr);
    break;

  case 'F':
    ChangeFnameFormat(In_Ptr);
    break;

  case 'D':
    ChangeDzDrDt(In_Ptr);
    break;

  case 'N':
    ChangeNzNrNtNa(In_Ptr);
    break;

  case 'C':
    ChangeRecord(In_Ptr);
    break;

  case 'W':
    ChangeWth(In_Ptr);
    break;

  case 'R':
    ChangeRanSeed(In_Ptr);
    break;

  case 'L':
    ChangeLayerSpecs(In_Ptr);
    break;

  case 'P':
    ChangeNumPhotons(In_Ptr);
    break;

  case 'S':
    ChangeSourceType(In_Ptr);
    break;

  case 'Z':
    ChangeStartP(In_Ptr);
    break;

  case 'O':
    PutInputToFile(stdout, In_Ptr);
    break;

  case 'H':
    ShowChangeMenu();
    break;

  case 'Q':
    return 1;

  case 'X':
    return 2;

  default:
    puts("...Unknown command");
  }

  return 0;
}

/***************************************************************************
 *   Return 1 if quit change and start simulation;
 *   return 0 if exit to main menu.
 ****/
Boolean
RunChangedInput(InStru * In_Ptr)
{
  char        string[STRLEN];
  FILE       *fp;
  int         branch;

  printf("Any changes to the input parameters? (Y/N)");
  do {
    gets(string);
  } while (!strlen(string));

  while (toupper(string[0]) == 'Y') {
    do {
      do {
	printf("\n> Change menu (h for help) => ");
	gets(string);
      } while (!strlen(string));

      if (branch = BranchChangeMenu(string, In_Ptr))
	break;			/* string[0] is 'X' or 'Q'. */
    } while (1);

    printf("Do you want to save the input to a file? (Y/N)");
    gets(string);
    if (toupper(string[0]) == 'Y') {
      printf("Give the file name to save input: ( .mci): ");
      gets(string);
      if ((fp = fopen(string, "w")) == NULL)
	puts("Can not open the file to write.");
      else
	PutInputToFile(fp, In_Ptr);
    }
    if (branch == 1) {		/* quit change menu and start simulation. */
      if (!CheckInputConsis(In_Ptr)) {
	do {
	  printf("Change input or exit to main menu (c/x): ");
	  gets(string);
	} while (!strlen(string) ||
		 toupper(string[0]) != 'X' && toupper(string[0]) != 'C');

	if (toupper(string[0]) == 'X') {
	  free(In_Ptr->mediumlist);
	  free(In_Ptr->layerspecs);
	  return 0;
	} else
	  string[0] = 'Y';	/* continue to change parameters. */
      } else
	return 1;

    } else {			/* exit to menu. */
      free(In_Ptr->mediumlist);
      free(In_Ptr->layerspecs);
      return 0;
    }
  }

  return 1;
}

/**************************************************************************
 *	Return 1, if the name in the name list.
 *	Return 0, otherwise.
 ****/
Boolean
NameInList(char *Name, NameLink List)
{
  while (List != NULL) {
    if (strcmp(Name, List->name) == 0)
      return (1);
    List = List->next;
  };
  return (0);
}

/**************************************************************************
 *	Add the name to the name list.
 ****/
void
AddNameToList(char *Name, NameLink * List_Ptr)
{
  NameLink    list = *List_Ptr;

  if (list == NULL) {		/* first node. */
    *List_Ptr = list = (NameLink) malloc(sizeof(NameNode));
    strcpy(list->name, Name);
    list->next = NULL;
  } else {			/* subsequent nodes. */
    /* Move to the last node. */
    while (list->next != NULL)
      list = list->next;

    /* Append a node to the list. */
    list->next = (NameLink) malloc(sizeof(NameNode));
    list = list->next;
    strcpy(list->name, Name);
    list->next = NULL;
  }
}

/**************************************************************************
 *	Check against duplicated file names.
 *
 *	A linked list is set up to store the file names used
 *	in this input data file.
 ****/
Boolean
FnameTaken(char *fname, NameLink * List_Ptr)
{
  if (NameInList(fname, *List_Ptr))
    return (1);
  else {
    AddNameToList(fname, List_Ptr);
    return (0);
  }
}

/**************************************************************************
 *	Free each node in the file name list.
 ****/
void
FreeFnameList(NameLink List)
{
  NameLink    next;

  while (List != NULL) {
    next = List->next;
    free(List);
    List = next;
  }
}

/*******************************************************************************
 *  Check the whether the flag end is met.
 *  The file position is restored to the current position
 *  at the end of the inquery.
 ****/
Boolean
EndOfRunsQ(FILE ** FilePP)
{
  char        buf[STRLEN];
  Boolean     found = 1;	/* found end of runs. */
  long        file_pos;

  file_pos = ftell(*FilePP);	/* record file position. */
  strcpy(buf, FindDataLine(*FilePP));
  if (buf[0] == '\0') {
    found = 0;
    printf("Missing end.\n");
  } else if (strstr(buf, "end") == NULL)
    found = 0;

  fseek(*FilePP, file_pos, SEEK_SET);	/* restore postion. */
  return (found);
}

/**************************************************************************
 *  Check the input parameters for all runs.
 *  This function will count number of runs and assign it to
 *  In_Ptr->num_runs.
 ****/
void
CheckParamFromFile(FILE * Fp, InStru * In_Ptr)
{
  short       i_run = 0;
  NameLink    head = NULL;
  Boolean     name_taken;	/* output files share the same */
  /* file name. */
  long        file_pos;

  if (!ReadMediumListQ(Fp, In_Ptr))
    exit(1);

  file_pos = ftell(Fp);
  do {
    printf("Checking input data for run %d\n", ++i_run);
    ReadParam(Fp, In_Ptr);

    name_taken = FnameTaken(In_Ptr->out_fname, &head);
    if (name_taken) {
      printf("file name %s duplicated.\n", In_Ptr->out_fname);
      exit(1);
    }
    free(In_Ptr->layerspecs);
  } while (!EndOfRunsQ(&Fp));

  In_Ptr->num_runs = i_run;
  FreeFnameList(head);
  fseek(Fp, file_pos, SEEK_SET);
}

/**************************************************************************
 *	Allocate the arrays in OutStru for one run, and
 *	array elements are automatically initialized to zeros.
 *
 *	Remember that the indices for Rd_r[], Td_r[],
 *	& A_rz[][iz] start from -1 storing the collimated
 *	responses.
 ****/
void
InitOutputData(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       nz = In_Ptr->nz;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  short       nl = In_Ptr->num_layers;
  /* remember to use nl+2 because of 2 for ambient. */

  if (nz <= 0 || nr <= 0 || na <= 0 || nl <= 0) {
    printf("Invalid grid parameters.\n");
    exit(1);
  }
  /* Init pure numbers. */
  Out_Ptr->Rsp = 0.0;
  Out_Ptr->Rb = 0.0;
  Out_Ptr->Rd = 0.0;
  Out_Ptr->Td = 0.0;
  Out_Ptr->Tb = 0.0;
  Out_Ptr->A = 0.0;

  Out_Ptr->Rbe = 0.0;
  Out_Ptr->Rde = 0.0;
  Out_Ptr->Tde = 0.0;
  Out_Ptr->Tbe = 0.0;
  Out_Ptr->Ae = 0.0;

  /* Allocate the 1D, 2D and 3D arrays. */
  Out_Ptr->Rd_rat = (In_Ptr->record.Rd_rat)
    ? AllocArray3D(0, nr - 1, 0, na - 1, 0, nt - 1) : NULL;
  Out_Ptr->Rd_ra = (In_Ptr->record.Rd_ra)
    ? AllocArray2D(0, nr - 1, 0, na - 1) : NULL;
  Out_Ptr->Rd_rt = (In_Ptr->record.Rd_rt)
    ? AllocArray2D(0, nr - 1, 0, nt - 1) : NULL;
  Out_Ptr->Rd_at = (In_Ptr->record.Rd_at)
    ? AllocArray2D(0, na - 1, 0, nt - 1) : NULL;
  Out_Ptr->Rd_r = (In_Ptr->record.Rd_r)
    ? AllocArray1D(0, nr - 1) : NULL;
  Out_Ptr->Rd_a = (In_Ptr->record.Rd_a)
    ? AllocArray1D(0, na - 1) : NULL;
  Out_Ptr->Rd_t = (In_Ptr->record.Rd_t)
    ? AllocArray1D(0, nt - 1) : NULL;

  Out_Ptr->Td_rat = (In_Ptr->record.Td_rat)
    ? AllocArray3D(0, nr - 1, 0, na - 1, 0, nt - 1) : NULL;
  Out_Ptr->Td_ra = (In_Ptr->record.Td_ra)
    ? AllocArray2D(0, nr - 1, 0, na - 1) : NULL;
  Out_Ptr->Td_rt = (In_Ptr->record.Td_rt)
    ? AllocArray2D(0, nr - 1, 0, nt - 1) : NULL;
  Out_Ptr->Td_at = (In_Ptr->record.Td_at)
    ? AllocArray2D(0, na - 1, 0, nt - 1) : NULL;
  Out_Ptr->Td_r = (In_Ptr->record.Td_r)
    ? AllocArray1D(0, nr - 1) : NULL;
  Out_Ptr->Td_a = (In_Ptr->record.Td_a)
    ? AllocArray1D(0, na - 1) : NULL;
  Out_Ptr->Td_t = (In_Ptr->record.Td_t)
    ? AllocArray1D(0, nt - 1) : NULL;

  Out_Ptr->A_rzt = (In_Ptr->record.A_rzt)
    ? AllocArray3D(0, nr - 1, 0, nz - 1, 0, nt - 1) : NULL;
  Out_Ptr->Ab_zt = (In_Ptr->record.A_rzt)
    ? AllocArray2D(0, nz - 1, 0, nt - 1) : NULL;
  Out_Ptr->A_rz = (In_Ptr->record.A_rz)
    ? AllocArray2D(0, nr - 1, 0, nz - 1) : NULL;
  Out_Ptr->Ab_z = (In_Ptr->record.A_rz)
    ? AllocArray1D(0, nz - 1) : NULL;
  Out_Ptr->A_zt = (In_Ptr->record.A_zt)
    ? AllocArray2D(0, nz - 1, 0, nt - 1) : NULL;
  Out_Ptr->A_z = (In_Ptr->record.A_z)
    ? AllocArray1D(0, nz - 1) : NULL;
  Out_Ptr->A_t = (In_Ptr->record.A_t)
    ? AllocArray1D(0, nt - 1) : NULL;
}

/**************************************************************************
 *	Undo what InitOutputData did.
 *  i.e. free the data allocations.
 ****/
void
FreeData(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       nz = In_Ptr->nz;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;

  free(In_Ptr->layerspecs);

  FreeArray3D(Out_Ptr->Rd_rat, 0, nr - 1, 0, na - 1, 0, nt - 1);
  FreeArray2D(Out_Ptr->Rd_ra, 0, nr - 1, 0, na - 1);
  FreeArray2D(Out_Ptr->Rd_rt, 0, nr - 1, 0, nt - 1);
  FreeArray2D(Out_Ptr->Rd_at, 0, na - 1, 0, nt - 1);
  FreeArray1D(Out_Ptr->Rd_r, 0, nr - 1);
  FreeArray1D(Out_Ptr->Rd_a, 0, na - 1);
  FreeArray1D(Out_Ptr->Rd_t, 0, nt - 1);

  FreeArray3D(Out_Ptr->Td_rat, 0, nr - 1, 0, na - 1, 0, nt - 1);
  FreeArray2D(Out_Ptr->Td_ra, 0, nr - 1, 0, na - 1);
  FreeArray2D(Out_Ptr->Td_rt, 0, nr - 1, 0, nt - 1);
  FreeArray2D(Out_Ptr->Td_at, 0, na - 1, 0, nt - 1);
  FreeArray1D(Out_Ptr->Td_r, 0, nr - 1);
  FreeArray1D(Out_Ptr->Td_a, 0, na - 1);
  FreeArray1D(Out_Ptr->Td_t, 0, nt - 1);

  FreeArray3D(Out_Ptr->A_rzt, 0, nr - 1, 0, nz - 1, 0, nt - 1);
  FreeArray2D(Out_Ptr->Ab_zt, 0, nz - 1, 0, nt - 1);
  FreeArray2D(Out_Ptr->A_rz, 0, nr - 1, 0, nz - 1);
  FreeArray1D(Out_Ptr->Ab_z, 0, nz - 1);
  FreeArray2D(Out_Ptr->A_zt, 0, nz - 1, 0, nt - 1);
  FreeArray1D(Out_Ptr->A_z, 0, nz - 1);
  FreeArray1D(Out_Ptr->A_t, 0, nt - 1);
}

/**************************************************************************
 *	Scale Rd and Td properly.
 *
 *	"a" stands for angle alpha.
 ****
 *	Scale Rd(r,a) and Td(r,a) by
 *      (area perpendicular to photon direction)
 *		x(solid angle)x(No. of photons).
 *	or
 *		[2*PI*r*dr*cos(a)]x[2*PI*sin(a)*da]x[No. of photons]
 *	or
 *		[2*PI*PI*dr*da*r*sin(2a)]x[No. of photons]
 ****
 *	Scale Rd(r) and Td(r) by
 *		(area on the surface)x(No. of photons).
 ****
 *	Scale Rd(a) and Td(a) by
 *		(solid angle) cos(a) x(No. of photons).
 ****
 *      Mode = 0, scale Rd and Td; Mode = 1, unscale Rd and Td.
 ****/
void
ScaleRdTd(InStru * In_Ptr, OutStru * Out_Ptr, char Mode)
{
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  short       ir, ia, it;
  double      scale1, scale2;

  scale1 = (double) In_Ptr->num_photons;
  if (Mode == 0) {
    Out_Ptr->Rde = 1 / scale1 * sqrt(Out_Ptr->Rde
				   - Out_Ptr->Rd * Out_Ptr->Rd / scale1);
    Out_Ptr->Tde = 1 / scale1 * sqrt(Out_Ptr->Tde
				   - Out_Ptr->Td * Out_Ptr->Td / scale1);
    Out_Ptr->Rbe = 1 / scale1 * sqrt(Out_Ptr->Rbe
				   - Out_Ptr->Rb * Out_Ptr->Rb / scale1);
    Out_Ptr->Tbe = 1 / scale1 * sqrt(Out_Ptr->Tbe
				   - Out_Ptr->Tb * Out_Ptr->Tb / scale1);

    Out_Ptr->Rd /= scale1;
    Out_Ptr->Td /= scale1;
    Out_Ptr->Rb = Out_Ptr->Rb / scale1 + Out_Ptr->Rsp;
    Out_Ptr->Tb /= scale1;
  } else {
    Out_Ptr->Rd *= scale1;
    Out_Ptr->Td *= scale1;
    Out_Ptr->Rb = (Out_Ptr->Rb - Out_Ptr->Rsp) * scale1;
    Out_Ptr->Tb *= scale1;

    Out_Ptr->Rde = (scale1 * Out_Ptr->Rde) * (scale1 * Out_Ptr->Rde)
      + 1 / scale1 * Out_Ptr->Rd * Out_Ptr->Rd;
    Out_Ptr->Tde = (scale1 * Out_Ptr->Tde) * (scale1 * Out_Ptr->Tde)
      + 1 / scale1 * Out_Ptr->Td * Out_Ptr->Td;
    Out_Ptr->Rbe = (scale1 * Out_Ptr->Rbe) * (scale1 * Out_Ptr->Rbe)
      + 1 / scale1 * Out_Ptr->Rb * Out_Ptr->Rb;
    Out_Ptr->Tbe = (scale1 * Out_Ptr->Tbe) * (scale1 * Out_Ptr->Tbe)
      + 1 / scale1 * Out_Ptr->Tb * Out_Ptr->Tb;
  }

  scale1 = dt * In_Ptr->num_photons;
  if (In_Ptr->record.Rd_t)
    for (it = 0; it < nt; it++)
      if (Mode == 0)		/* scale Rd_t. */
	Out_Ptr->Rd_t[it] /= scale1;
      else			/* unscale Rd_t. */
	Out_Ptr->Rd_t[it] *= scale1;

  if (In_Ptr->record.Td_t)
    for (it = 0; it < nt; it++)
      if (Mode == 0)		/* scale Td_t. */
	Out_Ptr->Td_t[it] /= scale1;
      else			/* unscale Rd_t. */
	Out_Ptr->Td_t[it] *= scale1;

  scale1 = 2.0 * PI * dr * dr * In_Ptr->num_photons;
  /* area is 2*PI*[(ir+0.5)*dr]*dr.  ir + 0.5 to be added. */

  if (In_Ptr->record.Rd_r)
    for (ir = 0; ir < nr; ir++) {
      scale2 = 1.0 / ((ir + 0.5) * scale1);
      if (Mode == 0)		/* scale Rd_r. */
	Out_Ptr->Rd_r[ir] *= scale2;
      else			/* unscale Rd_r. */
	Out_Ptr->Rd_r[ir] /= scale2;
    }

  if (In_Ptr->record.Td_r)
    for (ir = 0; ir < nr; ir++) {
      scale2 = 1.0 / ((ir + 0.5) * scale1);
      if (Mode == 0)		/* scale Td_r. */
	Out_Ptr->Td_r[ir] *= scale2;
      else			/* unscale Td_r. */
	Out_Ptr->Td_r[ir] /= scale2;
    }

  scale1 *= dt;
  if (In_Ptr->record.Rd_rt)
    for (ir = 0; ir < nr; ir++)
      for (it = 0; it < nt; it++) {
	scale2 = 1.0 / ((ir + 0.5) * scale1);
	if (Mode == 0)		/* scale Rd_rt. */
	  Out_Ptr->Rd_rt[ir][it] *= scale2;
	else			/* unscale Rd_rt. */
	  Out_Ptr->Rd_rt[ir][it] *= scale2;
      }

  if (In_Ptr->record.Td_rt)
    for (ir = 0; ir < nr; ir++)
      for (it = 0; it < nt; it++) {
	scale2 = 1.0 / ((ir + 0.5) * scale1);
	if (Mode == 0)		/* scale Td_rt. */
	  Out_Ptr->Td_rt[ir][it] *= scale2;
	else			/* unscale Td_rt. */
	  Out_Ptr->Td_rt[ir][it] /= scale2;
      }

  scale1 = PI * da * In_Ptr->num_photons;
  /* solid angle times cos(a) is PI*sin(2a)*da. sin(2a) to be added. */

  if (In_Ptr->record.Rd_a)
    for (ia = 0; ia < na; ia++) {
      scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
      if (Mode == 0)		/* scale Rd_a. */
	Out_Ptr->Rd_a[ia] *= scale2;
      else			/* unscale Rd_a. */
	Out_Ptr->Rd_a[ia] /= scale2;
    }

  if (In_Ptr->record.Td_a)
    for (ia = 0; ia < na; ia++) {
      scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
      if (Mode == 0)		/* scale Td_a. */
	Out_Ptr->Td_a[ia] *= scale2;
      else			/* unscale Td_a. */
	Out_Ptr->Td_a[ia] /= scale2;
    }

  scale1 *= dt;
  if (In_Ptr->record.Rd_at)
    for (ia = 0; ia < na; ia++)
      for (it = 0; it < nt; it++) {
	scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
	if (Mode == 0)		/* scale Rd_at. */
	  Out_Ptr->Rd_at[ia][it] *= scale2;
	else			/* unscale Rd_at. */
	  Out_Ptr->Rd_at[ia][it] /= scale2;
      }

  if (In_Ptr->record.Td_at)
    for (ia = 0; ia < na; ia++)
      for (it = 0; it < nt; it++) {
	scale2 = 1.0 / (sin(2.0 * (ia + 0.5) * da) * scale1);
	if (Mode == 0)		/* scale Td_at. */
	  Out_Ptr->Td_at[ia][it] *= scale2;
	else			/* unscale Td_at. */
	  Out_Ptr->Td_at[ia][it] /= scale2;
      }

  scale1 = 2.0 * PI * dr * dr * PI * da * In_Ptr->num_photons;
  if (In_Ptr->record.Rd_ra)
    for (ir = 0; ir < nr; ir++)
      for (ia = 0; ia < na; ia++) {
	scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
	if (Mode == 0)		/* scale Rd_ra. */
	  Out_Ptr->Rd_ra[ir][ia] *= scale2;
	else			/* unscale Rd_ra. */
	  Out_Ptr->Rd_ra[ir][ia] /= scale2;
      }

  if (In_Ptr->record.Td_ra)
    for (ir = 0; ir < nr; ir++)
      for (ia = 0; ia < na; ia++) {
	scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
	if (Mode == 0)		/* scale Td_ra. */
	  Out_Ptr->Td_ra[ir][ia] *= scale2;
	else			/* unscale Td_ra. */
	  Out_Ptr->Td_ra[ir][ia] /= scale2;
      }

  scale1 *= dt;
  if (In_Ptr->record.Rd_rat)
    for (ir = 0; ir < nr; ir++)
      for (ia = 0; ia < na; ia++)
	for (it = 0; it < nt; it++) {
	  scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
	  if (Mode == 0)	/* scale Rd_rat. */
	    Out_Ptr->Rd_rat[ir][ia][it] *= scale2;
	  else			/* unscale Rd_rat. */
	    Out_Ptr->Rd_rat[ir][ia][it] /= scale2;
	}

  if (In_Ptr->record.Td_rat)
    for (ir = 0; ir < nr; ir++)
      for (ia = 0; ia < na; ia++)
	for (it = 0; it < nt; it++) {
	  scale2 = 1.0 / ((ir + 0.5) * sin(2.0 * (ia + 0.5) * da) * scale1);
	  if (Mode == 0)	/* scale Td_rat. */
	    Out_Ptr->Td_rat[ir][ia][it] *= scale2;
	  else			/* unscale Td_rat. */
	    Out_Ptr->Td_rat[ir][ia][it] /= scale2;
	}
}

/**************************************************************************
 *	Scale absorption arrays properly.
 *      Mode = 0, scale A; Mode = 1, unscale A.
 ****/
void
ScaleA(InStru * In_Ptr, OutStru * Out_Ptr, char Mode)
{
  short       nz = In_Ptr->nz;
  short       nr = In_Ptr->nr;
  short       nt = In_Ptr->nt;
  double      dz = In_Ptr->dz;
  double      dr = In_Ptr->dr;
  double      dt = In_Ptr->dt;
  short       iz, ir, it;
  double      scale1, scale2;

  scale1 = (double) In_Ptr->num_photons;
  if (Mode == 0) {		/* scale A. */
    Out_Ptr->Ae = 1 / scale1 * sqrt(Out_Ptr->Ae
				    - Out_Ptr->A * Out_Ptr->A / scale1);
    Out_Ptr->A /= scale1;
  } else {			/* unscale A. */
    Out_Ptr->A *= scale1;
    Out_Ptr->Ae = (scale1 * Out_Ptr->Ae) * (scale1 * Out_Ptr->Ae)
      + 1 / scale1 * Out_Ptr->A * Out_Ptr->A;
  }

  scale2 = scale1 * dt;
  if (In_Ptr->record.A_t)
    for (it = 0; it < nt; it++)
      if (Mode == 0)		/* scale A_t. */
	Out_Ptr->A_t[it] /= scale2;
      else			/* unscale A_t. */
	Out_Ptr->A_t[it] *= scale2;

  scale1 *= dz;
  if (In_Ptr->record.A_z)
    for (iz = 0; iz < nz; iz++)
      if (Mode == 0)		/* scale A_z. */
	Out_Ptr->A_z[iz] /= scale1;
      else			/* unscale A_z. */
	Out_Ptr->A_z[iz] *= scale1;

  scale2 = scale1 * dt;
  if (In_Ptr->record.A_zt)
    for (iz = 0; iz < nz; iz++)
      for (it = 0; it < nt; it++)
	if (Mode == 0)		/* scale A_zt. */
	  Out_Ptr->A_zt[iz][it] /= scale2;
	else			/* unscale A_zt. */
	  Out_Ptr->A_zt[iz][it] *= scale2;

  if (In_Ptr->record.A_rz)
    for (iz = 0; iz < nz; iz++)
      if (Mode == 0)		/* scale Ab_z. */
	Out_Ptr->Ab_z[iz] /= scale1;
      else			/* unscale Ab_z. */
	Out_Ptr->Ab_z[iz] *= scale1;

  if (In_Ptr->record.A_rzt)
    for (iz = 0; iz < nz; iz++)
      for (it = 0; it < nt; it++)
	if (Mode == 0)		/* scale Ab_zt. */
	  Out_Ptr->Ab_zt[iz][it] /= scale2;
	else			/* unscale Ab_zt. */
	  Out_Ptr->Ab_zt[iz][it] *= scale2;

  scale1 = 2.0 * PI * dr * dr * dz * In_Ptr->num_photons;
  if (In_Ptr->record.A_rz)
    for (ir = 0; ir < nr; ir++)
      for (iz = 0; iz < nz; iz++)
	if (Mode == 0)		/* scale A_rz. */
	  Out_Ptr->A_rz[ir][iz] /= (ir + 0.5) * scale1;
	else			/* unscale A_rz. */
	  Out_Ptr->A_rz[ir][iz] *= (ir + 0.5) * scale1;

  scale2 = scale1 * dt;
  if (In_Ptr->record.A_rzt)
    for (ir = 0; ir < nr; ir++)
      for (iz = 0; iz < nz; iz++)
	for (it = 0; it < nt; it++)
	  if (Mode == 0)	/* scale A_rzt. */
	    Out_Ptr->A_rzt[ir][iz][it] /= (ir + 0.5) * scale2;
	  else			/* unscale A_rzt. */
	    Out_Ptr->A_rzt[ir][iz][it] *= (ir + 0.5) * scale2;
}

/**************************************************************************
 *	Scale results of current run.
 *      Mode = 0, scale result; Mode = 1, unscale result.
 ****/
void
ScaleResult(InStru * In_Ptr, OutStru * Out_Ptr, char Mode)
{
  ScaleRdTd(In_Ptr, Out_Ptr, Mode);
  ScaleA(In_Ptr, Out_Ptr, Mode);
}

/**************************************************************************
 *	Write the version number as the first string in the
 *	file.
 *	Use chars only so that they can be read as either
 *	ASCII or binary.
 ****/
void
WriteVersion(FILE * Fp, char *Version)
{
  fprintf(Fp,
	  "%s \t# Version number of the file format.\n\n",
	  Version);
  fprintf(Fp, "####\n# Data categories include: \n");
  fprintf(Fp, "# InParam, RAT, \n");
  fprintf(Fp, "# Rd_r\tRd_a\tRd_ra\tRd_t\tRd_rt\tRd_at\tRd_rat\n");
  fprintf(Fp, "# Td_r\tTd_a\tTd_ra\tTd_t\tTd_rt\tTd_at\tTd_rat\n");
  fprintf(Fp, "# A_z\tA_rz\tA_t\tA_zt\tA_rzt\n");
  fprintf(Fp, "####\n\n");
}

/***************************************************************************
 * Save the status of the random number generater to output file.
 ****/
void
SaveRandomStatus(FILE * Fp)
{
  long        status[57];
  int         i;

  RandomGen(2, 0, status);	/* get the status. */
  fprintf(Fp, "# status of the random number generator:");

  for (i = 0; i < 57; i++)
    if (i % 5)
      fprintf(Fp, "%14ld", status[i]);
    else
      fprintf(Fp, "\n%14ld ", status[i]);

  fprintf(Fp, "\n\n");
}

/***************************************************************************
 * Read and restore the status of random number generater
 * from previous output file.
 ****/
void
RestoreRandomStatus(FILE * Fp)
{
  char        buf[STRLEN];
  long        status[57];
  int         i;

  do {
    fgets(buf, STRLEN, Fp);
  } while (buf[0] != '#');

  for (i = 0; i < 57; i++)
    fscanf(Fp, "%ld", &status[i]);

  RandomGen(3, 0, status);	/* restore the status. */
}

/**************************************************************************
 *	Write reflectance, absorption, transmission.
 ****/
void
WriteRAT(FILE * Fp, OutStru * Out_Ptr)
{
  /* flag. */
  fprintf(Fp,
	  "RAT #Reflectance, absorption, transmittance.\n");
  fprintf(Fp, "# Average \tStandard Err \tRel Err\n");
  fprintf(Fp, "%-14.6G \t\t\t\t#Rsp: Specular reflectance.\n",
	  Out_Ptr->Rsp);
  fprintf(Fp, "%-14.6G \t%-14.6G %6.2f%%\t#Rb: Ballistic reflectance.\n",
	  Out_Ptr->Rb, Out_Ptr->Rbe,
	  (Out_Ptr->Rb) ? Out_Ptr->Rbe / Out_Ptr->Rb * 100 : 0);
  fprintf(Fp, "%-14.6G \t%-14.6G %6.2f%%\t#Rd: Diffuse reflectance.\n",
	  Out_Ptr->Rd, Out_Ptr->Rde,
	  (Out_Ptr->Rd) ? Out_Ptr->Rde / Out_Ptr->Rd * 100 : 0);
  fprintf(Fp, "%-14.6G \t%-14.6G %6.2f%%\t#A:  Absorbed fraction.\n",
	  Out_Ptr->A, Out_Ptr->Ae,
	  (Out_Ptr->A) ? Out_Ptr->Ae / Out_Ptr->A * 100 : 0);
  fprintf(Fp, "%-14.6G \t%-14.6G %6.2f%%\t#Tb: Ballistic transmittance.\n",
	  Out_Ptr->Tb, Out_Ptr->Tbe,
	  (Out_Ptr->Tb) ? Out_Ptr->Tbe / Out_Ptr->Tb * 100 : 0);
  fprintf(Fp, "%-14.6G \t%-14.6G %6.2f%%\t#Td: Diffuse transmittance.\n",
	  Out_Ptr->Td, Out_Ptr->Tde,
	  (Out_Ptr->Td) ? Out_Ptr->Tde / Out_Ptr->Td * 100 : 0);

  fprintf(Fp, "\n");
}

/**************************************************************************
 *	Read reflectance, absorption, transmission.
 ****/
void
ReadRAT(FILE * Fp, OutStru * Out_Ptr)
{
  char        buf[STRLEN];

  FindDataLine(Fp);		/* skip RAT line. */

  strcpy(buf, FindDataLine(Fp));
  sscanf(buf, "%lf", &(Out_Ptr->Rsp));

  strcpy(buf, FindDataLine(Fp));
  sscanf(buf, "%lf %lf", &(Out_Ptr->Rb), &(Out_Ptr->Rbe));

  strcpy(buf, FindDataLine(Fp));
  sscanf(buf, "%lf %lf", &(Out_Ptr->Rd), &(Out_Ptr->Rde));

  strcpy(buf, FindDataLine(Fp));
  sscanf(buf, "%lf %lf", &(Out_Ptr->A), &(Out_Ptr->Ae));

  strcpy(buf, FindDataLine(Fp));
  sscanf(buf, "%lf %lf", &(Out_Ptr->Tb), &(Out_Ptr->Tbe));

  strcpy(buf, FindDataLine(Fp));
  sscanf(buf, "%lf %lf", &(Out_Ptr->Td), &(Out_Ptr->Tde));
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOAb_zt(FILE * Fp,
	short Nz,
	short Nt,
	OutStru * Out_Ptr, char Mode)
{
  short       iz, it, i = 0;

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# Ab[z][t]. [1/(cm ps)]",
	    "# Ab[0][0], [0][1],..[0][nt-1]",
	    "# Ab[1][0], [1][1],..[1][nt-1]",
	    "# ...",
	    "# Ab[nz-1][0], [nz-1][1],..[nz-1][nt-1]",
	    "Ab_zt");
  else
    FindDataLine(Fp);		/* skip A_z line. */

  for (iz = 0; iz < Nz; iz++)
    for (it = 0; it < Nt; it++) {
      if (Mode == 1) {
	fprintf(Fp, "%12.4E ", Out_Ptr->Ab_zt[iz][it]);
	if (++i % 5 == 0)
	  fprintf(Fp, "\n");
      } else
	fscanf(Fp, "%lf", &(Out_Ptr->Ab_zt[iz][it]));
    }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOA_rzt(FILE * Fp,
	short Nr,
	short Nz,
	short Nt,
	OutStru * Out_Ptr, char Mode)
{
  short       iz, ir, it, i = 0;

  IOAb_zt(Fp, Nz, Nt, Out_Ptr, Mode);

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# A[r][z][t]. [1/(cm3 ps)]",
	    "# A[0][0][0], [0][0][1],..[0][0][nt-1]",
	    "# A[0][1][0], [0][1][1],..[0][1][nt-1]",
	    "# ...",
	    "# A[nr-1][nz-1][0], [nr-1][nz-1][1],..[nr-1][nz-1][nt-1]",
	    "A_rzt");
  else
    FindDataLine(Fp);

  for (ir = 0; ir < Nr; ir++)
    for (iz = 0; iz < Nz; iz++)
      for (it = 0; it < Nt; it++) {
	if (Mode == 1) {
	  fprintf(Fp, "%12.4E ", Out_Ptr->A_rzt[ir][iz][it]);
	  if (++i % 5 == 0)
	    fprintf(Fp, "\n");
	} else
	  fscanf(Fp, "%lf", &(Out_Ptr->A_rzt[ir][iz][it]));
      }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      1 number each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOAb_z(FILE * Fp,
       short Nz,
       OutStru * Out_Ptr, char Mode)
{
  short       iz;

  if (Mode == 1)
    fprintf(Fp,
	    "Ab_z #Ab[0], [1],..Ab[nz-1]. [1/cm]\n");	/* flag. */
  else
    FindDataLine(Fp);

  for (iz = 0; iz < Nz; iz++) {
    if (Mode == 1)
      fprintf(Fp, "%12.4E\n", Out_Ptr->Ab_z[iz]);
    else
      fscanf(Fp, "%lf", &(Out_Ptr->Ab_z[iz]));
  }

  if (Mode == 1)
    fprintf(Fp, "\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOA_rz(FILE * Fp,
       short Nr,
       short Nz,
       OutStru * Out_Ptr, char Mode)
{
  short       iz, ir, i = 0;

  IOAb_z(Fp, Nz, Out_Ptr, Mode);

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# A[r][z]. [1/cm3]",
	    "# A[0][0], [0][1],..[0][nz-1]",
	    "# ...",
	    "# A[nr-1][0], [nr-1][1],..[nr-1][nz-1]",
	    "A_rz");
  else
    FindDataLine(Fp);		/* skip A_rz line. */

  for (ir = 0; ir < Nr; ir++)
    for (iz = 0; iz < Nz; iz++) {
      if (Mode == 1) {
	fprintf(Fp, "%12.4E ", Out_Ptr->A_rz[ir][iz]);
	if (++i % 5 == 0)
	  fprintf(Fp, "\n");
      } else
	fscanf(Fp, "%lf", &(Out_Ptr->A_rz[ir][iz]));
    }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOA_zt(FILE * Fp,
       short Nz,
       short Nt,
       OutStru * Out_Ptr, char Mode)
{
  short       iz, it, i = 0;

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# A[z][t]. [1/(cm ps)]",
	    "# A[0][0], [0][1],..[0][nt-1]",
	    "# A[1][0], [1][1],..[1][nt-1]",
	    "# ...",
	    "# A[nz-1][0], [nz-1][1],..[nz-1][nt-1]",
	    "A_zt");
  else
    FindDataLine(Fp);		/* skip A_zt line. */

  for (iz = 0; iz < Nz; iz++)
    for (it = 0; it < Nt; it++) {
      if (Mode == 1) {
	fprintf(Fp, "%12.4E ", Out_Ptr->A_zt[iz][it]);
	if (++i % 5 == 0)
	  fprintf(Fp, "\n");
      } else
	fscanf(Fp, "%lf", &(Out_Ptr->A_zt[iz][it]));
    }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      1 number each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOA_z(FILE * Fp,
      short Nz,
      OutStru * Out_Ptr, char Mode)
{
  short       iz;

  if (Mode == 1)
    fprintf(Fp,
	    "A_z #A[0], [1],..A[nz-1]. [1/cm]\n");	/* flag. */
  else
    FindDataLine(Fp);		/* skip A_z line. */

  for (iz = 0; iz < Nz; iz++) {
    if (Mode == 1)
      fprintf(Fp, "%12.4E\n", Out_Ptr->A_z[iz]);
    else
      fscanf(Fp, "%lf", &(Out_Ptr->A_z[iz]));
  }

  if (Mode == 1)
    fprintf(Fp, "\n");
}

/**************************************************************************
 *      1 number each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOA_t(FILE * Fp,
      short Nt,
      OutStru * Out_Ptr, char Mode)
{
  short       it;

  if (Mode == 1)
    fprintf(Fp,
	    "A_t #A[0], [1],..A[nt-1]. [1/ps]\n");	/* flag. */
  else
    FindDataLine(Fp);

  for (it = 0; it < Nt; it++) {
    if (Mode == 1)
      fprintf(Fp, "%12.4E\n", Out_Ptr->A_t[it]);
    else
      fscanf(Fp, "%lf", &(Out_Ptr->A_t[it]));
  }

  if (Mode == 1)
    fprintf(Fp, "\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IORd_rat(FILE * Fp,
	 short Nr,
	 short Na,
	 short Nt,
	 OutStru * Out_Ptr, char Mode)
{
  short       ir, ia, it, i = 0;

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# Rd[r][a][t]. [1/(cm2 sr ps)]",
	    "# Rd[0][0][0], [0][0][1],..[0][0][nt-1]",
	    "# Rd[0][1][0], [0][1][1],..[0][1][nt-1]",
	    "# ...",
	    "# Rd[nr-1][na-1][0], [nr-1][na-1][1],..[nr-1][na-1][nt-1]",
	    "Rd_rat");
  else
    FindDataLine(Fp);

  for (ir = 0; ir < Nr; ir++)
    for (ia = 0; ia < Na; ia++)
      for (it = 0; it < Nt; it++) {
	if (Mode == 1) {
	  fprintf(Fp, "%12.4E ", Out_Ptr->Rd_rat[ir][ia][it]);
	  if (++i % 5 == 0)
	    fprintf(Fp, "\n");
	} else
	  fscanf(Fp, "%lf", &(Out_Ptr->Rd_rat[ir][ia][it]));
      }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IORd_ra(FILE * Fp,
	short Nr,
	short Na,
	OutStru * Out_Ptr, char Mode)
{
  short       ir, ia;

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# Rd[r][angle]. [1/(cm2 sr)].",
	    "# Rd[0][0], [0][1],..[0][na-1]",
	    "# Rd[1][0], [1][1],..[1][na-1]",
	    "# ...",
	    "# Rd[nr-1][0], [nr-1][1],..[nr-1][na-1]",
	    "Rd_ra");
  else
    FindDataLine(Fp);

  for (ir = 0; ir < Nr; ir++)
    for (ia = 0; ia < Na; ia++) {
      if (Mode == 1) {
	fprintf(Fp, "%12.4E ", Out_Ptr->Rd_ra[ir][ia]);
	if ((ir * Na + ia + 1) % 5 == 0)
	  fprintf(Fp, "\n");
      } else
	fscanf(Fp, "%lf", &(Out_Ptr->Rd_ra[ir][ia]));
    }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IORd_rt(FILE * Fp,
	short Nr,
	short Nt,
	OutStru * Out_Ptr, char Mode)
{
  short       ir, it, i = 0;

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# Rd[r][t]. [1/(cm2 ps)]",
	    "# Rd[0][0], [0][1],..[0][nt-1]",
	    "# Rd[0][0], [0][1],..[0][nt-1]",
	    "# ...",
	    "# Rd[nr-1][0], [nr-1][1],..[nr-1][nt-1]",
	    "Rd_rt");
  else
    FindDataLine(Fp);

  for (ir = 0; ir < Nr; ir++)
    for (it = 0; it < Nt; it++) {
      if (Mode == 1) {
	fprintf(Fp, "%12.4E ", Out_Ptr->Rd_rt[ir][it]);
	if (++i % 5 == 0)
	  fprintf(Fp, "\n");
      } else
	fscanf(Fp, "%lf", &(Out_Ptr->Rd_rt[ir][it]));
    }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IORd_at(FILE * Fp,
	short Na,
	short Nt,
	OutStru * Out_Ptr, char Mode)
{
  short       ia, it, i = 0;

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# Rd[a][t]. [1/(sr ps)]",
	    "# Rd[0][0], [0][1],..[0][nt-1]",
	    "# Rd[1][0], [1][1],..[1][nt-1]",
	    "# ...",
	    "# Rd[na-1][0], [na-1][1],..[na-1][nt-1]",
	    "Rd_at");
  else
    FindDataLine(Fp);

  for (ia = 0; ia < Na; ia++)
    for (it = 0; it < Nt; it++) {
      if (Mode == 1) {
	fprintf(Fp, "%12.4E ", Out_Ptr->Rd_at[ia][it]);
	if (++i % 5 == 0)
	  fprintf(Fp, "\n");
      } else
	fscanf(Fp, "%lf", &(Out_Ptr->Rd_at[ia][it]));
    }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      1 number each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IORd_r(FILE * Fp,
       short Nr,
       OutStru * Out_Ptr, char Mode)
{
  short       ir;

  if (Mode == 1)
    fprintf(Fp,
	    "Rd_r #Rd[0], [1],..Rd[nr-1]. [1/cm2]\n");	/* flag. */
  else
    FindDataLine(Fp);

  for (ir = 0; ir < Nr; ir++) {
    if (Mode == 1)
      fprintf(Fp, "%12.4E\n", Out_Ptr->Rd_r[ir]);
    else
      fscanf(Fp, "%lf", &(Out_Ptr->Rd_r[ir]));
  }

  if (Mode == 1)
    fprintf(Fp, "\n");
}

/**************************************************************************
 *      1 number each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IORd_a(FILE * Fp,
       short Na,
       OutStru * Out_Ptr, char Mode)
{
  short       ia;

  if (Mode == 1)
    fprintf(Fp,
	    "Rd_a #Rd[0], [1],..Rd[na-1]. [1/sr]\n");	/* flag. */
  else
    FindDataLine(Fp);

  for (ia = 0; ia < Na; ia++) {
    if (Mode == 1)
      fprintf(Fp, "%12.4E\n", Out_Ptr->Rd_a[ia]);
    else
      fscanf(Fp, "%lf", &(Out_Ptr->Rd_a[ia]));
  }

  if (Mode == 1)
    fprintf(Fp, "\n");
}

/**************************************************************************
 *      1 number each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IORd_t(FILE * Fp,
       short Nt,
       OutStru * Out_Ptr, char Mode)
{
  short       it;

  if (Mode == 1)
    fprintf(Fp,
	    "Rd_t #Rd[0], [1],..Rd[nt-1]. [1/ps]\n");	/* flag. */
  else
    FindDataLine(Fp);

  for (it = 0; it < Nt; it++) {
    if (Mode == 1)
      fprintf(Fp, "%12.4E\n", Out_Ptr->Rd_t[it]);
    else
      fscanf(Fp, "%lf", &(Out_Ptr->Rd_t[it]));
  }

  if (Mode == 1)
    fprintf(Fp, "\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOTd_rat(FILE * Fp,
	 short Nr,
	 short Na,
	 short Nt,
	 OutStru * Out_Ptr, char Mode)
{
  short       ir, ia, it, i = 0;

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# Td[r][a][t]. [1/(cm2 sr ps)]",
	    "# Td[0][0][0], [0][0][1],..[0][0][nt-1]",
	    "# Td[0][1][0], [0][1][1],..[0][1][nt-1]",
	    "# ...",
	    "# Td[nr-1][na-1][0], [nr-1][na-1][1],..[nr-1][na-1][nt-1]",
	    "Td_rat");
  else
    FindDataLine(Fp);

  for (ir = 0; ir < Nr; ir++)
    for (ia = 0; ia < Na; ia++)
      for (it = 0; it < Nt; it++) {
	if (Mode == 1) {
	  fprintf(Fp, "%12.4E ", Out_Ptr->Td_rat[ir][ia][it]);
	  if (++i % 5 == 0)
	    fprintf(Fp, "\n");
	} else
	  fscanf(Fp, "%lf", &(Out_Ptr->Td_rat[ir][ia][it]));
      }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOTd_ra(FILE * Fp,
	short Nr,
	short Na,
	OutStru * Out_Ptr, char Mode)
{
  short       ir, ia;

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# Td[r][angle]. [1/(cm2 sr)].",
	    "# Td[0][0], [0][1],..[0][na-1]",
	    "# Td[1][0], [1][1],..[1][na-1]",
	    "# ...",
	    "# Td[nr-1][0], [nr-1][1],..[nr-1][na-1]",
	    "Td_ra");
  else
    FindDataLine(Fp);

  for (ir = 0; ir < Nr; ir++)
    for (ia = 0; ia < Na; ia++) {
      if (Mode == 1) {
	fprintf(Fp, "%12.4E ", Out_Ptr->Td_ra[ir][ia]);
	if ((ir * Na + ia + 1) % 5 == 0)
	  fprintf(Fp, "\n");
      } else
	fscanf(Fp, "%lf", &(Out_Ptr->Td_ra[ir][ia]));
    }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOTd_rt(FILE * Fp,
	short Nr,
	short Nt,
	OutStru * Out_Ptr, char Mode)
{
  short       ir, it, i = 0;

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# Td[r][t]. [1/(cm2 ps)]",
	    "# Td[0][0], [0][1],..[0][nt-1]",
	    "# Td[0][0], [0][1],..[0][nt-1]",
	    "# ...",
	    "# Td[nr-1][0], [nr-1][1],..[nr-1][nt-1]",
	    "Td_rt");
  else
    FindDataLine(Fp);

  for (ir = 0; ir < Nr; ir++)
    for (it = 0; it < Nt; it++) {
      if (Mode == 1) {
	fprintf(Fp, "%12.4E ", Out_Ptr->Td_rt[ir][it]);
	if (++i % 5 == 0)
	  fprintf(Fp, "\n");
      } else
	fscanf(Fp, "%lf", &(Out_Ptr->Td_rt[ir][it]));
    }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      5 numbers each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOTd_at(FILE * Fp,
	short Na,
	short Nt,
	OutStru * Out_Ptr, char Mode)
{
  short       ia, it, i = 0;

  if (Mode == 1)
    fprintf(Fp,
	    "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	    "# Td[a][t]. [1/(sr ps)]",
	    "# Td[0][0], [0][1],..[0][nt-1]",
	    "# Td[1][0], [1][1],..[1][nt-1]",
	    "# ...",
	    "# Td[na-1][0], [na-1][1],..[na-1][nt-1]",
	    "Td_at");
  else
    FindDataLine(Fp);

  for (ia = 0; ia < Na; ia++)
    for (it = 0; it < Nt; it++) {
      if (Mode == 1) {
	fprintf(Fp, "%12.4E ", Out_Ptr->Td_at[ia][it]);
	if (++i % 5 == 0)
	  fprintf(Fp, "\n");
      } else
	fscanf(Fp, "%lf", &(Out_Ptr->Td_at[ia][it]));
    }

  if (Mode == 1)
    fprintf(Fp, "\n\n");
}

/**************************************************************************
 *      1 number each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOTd_r(FILE * Fp,
       short Nr,
       OutStru * Out_Ptr, char Mode)
{
  short       ir;

  if (Mode == 1)
    fprintf(Fp,
	    "Td_r #Td[0], [1],..Td[nr-1]. [1/cm2]\n");	/* flag. */
  else
    FindDataLine(Fp);

  for (ir = 0; ir < Nr; ir++) {
    if (Mode == 1)
      fprintf(Fp, "%12.4E\n", Out_Ptr->Td_r[ir]);
    else
      fscanf(Fp, "%lf", &(Out_Ptr->Td_r[ir]));
  }

  if (Mode == 1)
    fprintf(Fp, "\n");
}

/**************************************************************************
 *      1 number each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOTd_a(FILE * Fp,
       short Na,
       OutStru * Out_Ptr, char Mode)
{
  short       ia;

  if (Mode == 1)
    fprintf(Fp,
	    "Td_a #Td[0], [1],..Td[na-1]. [1/sr]\n");	/* flag. */
  else
    FindDataLine(Fp);

  for (ia = 0; ia < Na; ia++) {
    if (Mode == 1)
      fprintf(Fp, "%12.4E\n", Out_Ptr->Td_a[ia]);
    else
      fscanf(Fp, "%lf", &(Out_Ptr->Td_a[ia]));
  }

  if (Mode == 1)
    fprintf(Fp, "\n");
}

/**************************************************************************
 *      1 number each line.
 *      Mode = 0, read; Mode = 1, write.
 ****/
void
IOTd_t(FILE * Fp,
       short Nt,
       OutStru * Out_Ptr, char Mode)
{
  short       it;

  if (Mode == 1)
    fprintf(Fp,
	    "Td_t #Rd[0], [1],..Td[nt-1]. [1/ps]\n");	/* flag. */
  else
    FindDataLine(Fp);

  for (it = 0; it < Nt; it++) {
    if (Mode == 1)
      fprintf(Fp, "%12.4E\n", Out_Ptr->Td_t[it]);
    else
      fscanf(Fp, "%lf", &(Out_Ptr->Td_t[it]));
  }

  if (Mode == 1)
    fprintf(Fp, "\n");
}

/**************************************************************************
 *  Mode = 0, read result back from a output file.
 *  Mode = 1, write result to a output file;
 ****/
void
IOResult(FILE * Fp, InStru * In_Ptr,
	 OutStru * Out_Ptr,
	 char Mode)
{
  if (Mode == 1) {
    if (toupper(In_Ptr->out_fformat) == 'A')
      WriteVersion(Fp, "mcmloA2.0");
    else
      WriteVersion(Fp, "mcmloB2.0");

    PutInputToFile(Fp, In_Ptr);
    SaveRandomStatus(Fp);
    WriteRAT(Fp, Out_Ptr);
  } else {
    RestoreRandomStatus(Fp);
    ReadRAT(Fp, Out_Ptr);
  }

  /* reflectance, absorption, transmittance. */
  if (In_Ptr->record.A_rzt)
    IOA_rzt(Fp, In_Ptr->nr, In_Ptr->nz, In_Ptr->nt, Out_Ptr, Mode);
  if (In_Ptr->record.A_rz)
    IOA_rz(Fp, In_Ptr->nr, In_Ptr->nz, Out_Ptr, Mode);
  if (In_Ptr->record.A_zt)
    IOA_zt(Fp, In_Ptr->nz, In_Ptr->nt, Out_Ptr, Mode);
  if (In_Ptr->record.A_z)
    IOA_z(Fp, In_Ptr->nz, Out_Ptr, Mode);
  if (In_Ptr->record.A_t)
    IOA_t(Fp, In_Ptr->nt, Out_Ptr, Mode);

  if (In_Ptr->record.Rd_rat)
    IORd_rat(Fp, In_Ptr->nr, In_Ptr->na, In_Ptr->nt, Out_Ptr, Mode);
  if (In_Ptr->record.Rd_ra)
    IORd_ra(Fp, In_Ptr->nr, In_Ptr->na, Out_Ptr, Mode);
  if (In_Ptr->record.Rd_rt)
    IORd_rt(Fp, In_Ptr->nr, In_Ptr->nt, Out_Ptr, Mode);
  if (In_Ptr->record.Rd_at)
    IORd_at(Fp, In_Ptr->na, In_Ptr->nt, Out_Ptr, Mode);
  if (In_Ptr->record.Rd_r)
    IORd_r(Fp, In_Ptr->nr, Out_Ptr, Mode);
  if (In_Ptr->record.Rd_a)
    IORd_a(Fp, In_Ptr->na, Out_Ptr, Mode);
  if (In_Ptr->record.Rd_t)
    IORd_t(Fp, In_Ptr->nt, Out_Ptr, Mode);

  if (In_Ptr->record.Td_rat)
    IOTd_rat(Fp, In_Ptr->nr, In_Ptr->na, In_Ptr->nt, Out_Ptr, Mode);
  if (In_Ptr->record.Td_ra)
    IOTd_ra(Fp, In_Ptr->nr, In_Ptr->na, Out_Ptr, Mode);
  if (In_Ptr->record.Td_rt)
    IOTd_rt(Fp, In_Ptr->nr, In_Ptr->nt, Out_Ptr, Mode);
  if (In_Ptr->record.Td_at)
    IOTd_at(Fp, In_Ptr->na, In_Ptr->nt, Out_Ptr, Mode);
  if (In_Ptr->record.Td_r)
    IOTd_r(Fp, In_Ptr->nr, Out_Ptr, Mode);
  if (In_Ptr->record.Td_a)
    IOTd_a(Fp, In_Ptr->na, Out_Ptr, Mode);
  if (In_Ptr->record.Td_t)
    IOTd_t(Fp, In_Ptr->nt, Out_Ptr, Mode);

  fclose(Fp);
}
