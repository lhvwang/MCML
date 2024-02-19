/**************************************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *	1992.
 *
 *	Functions for reading files.
 ****/

#include "conv.h"

#define FIRSTINTERACTIONS 1

/**************************************************************************
 *      Kill the ith char (counting from 0), push the following
 *      chars forward by one.
 ****/ 
void 
KillChar(size_t i, char *Str)
{
  size_t      sl = strlen(Str);
 
  for (; i < sl; i++)
    Str[i] = Str[i + 1];
}

/**************************************************************************
 *      Eliminate the chars in a string which are not printing
 *      chars or spaces.
 *
 *      Spaces include ' ', '\f', '\t' etc.
 *
 *      Return 1 if no nonprinting chars found, otherwise
 *      return 0.
 ****/
Boolean
CheckCharQ(char *Str)
{
  Boolean     found = 0;        /* found bad char. */
  size_t      sl = strlen(Str);
  size_t      i = 0;
 
  while (i < sl)
    if (Str[i] < 0 || Str[i] > 255) {/* this condition actually happened. */
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
 *      Return 1 if this line is a comment line in which the
 *      first non-space character is "#", or a space line.
 *      Return 0 otherwise.
 ****/
Boolean
CommentLineQ(char *Buf)
{
  size_t      spn, cspn;
 
  spn = strspn(Buf, " \t");     /* length spanned by space or tab chars. */
 
  cspn = strcspn(Buf, "#\n");   /* length before the 1st # or return. */
 
  if (spn == cspn)              /* comment line or space line. */
    return (1);
  else                          /* the line has data. */
    return (0);
}

/**************************************************************************
 *	Skip space or comment lines and return data line only.
 ****/
char       *
FindDataLine(FILE * File_Ptr)
{
  static char buf[STRLEN];

  do {				/* skip space or comment lines. */
    if (fgets(buf, STRLEN, File_Ptr) == NULL) {
      printf("Incomplete data.\n");
      return (NULL);
    }
    CheckCharQ(buf);
  } while (CommentLineQ(buf));

  return (buf);
}

/**************************************************************************
 *      Get the filename and open it for reading, retry until the
 *      file can be opened, or a '.' is input.
 ****/
FILE       *
GetFile(char *Fname)
{
  FILE       *file = NULL;
 
  do {
    printf("Input filename of mcml output (or . to quit): ");
    scanf("%s", Fname);
    if (strlen(Fname) == 1 && Fname[0] == '.')
      break;
 
    file = fopen(Fname, "r");
  } while (file == NULL);
 
  return (file);
}

/**************************************************************************
 *      Check whether the file version is the same as Version.
 ****/
Boolean
CheckFileVersionQ(char *Buf, char *Version)
{
  if ((Buf[0] == '\0') || (strstr(Buf, Version) == NULL))
    return (0);
  else
    return (1);
}

/**************************************************************************
 *      Allocate the arrays in OutStru for one run, and
 *      array elements are automatically initialized to zeros.
 *
 *      Remember that the indices for Rd_r[], Td_r[],
 *      & A_rz[][iz] start from -1 storing the collimated
 *      responses.
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

  /* Allocate the 1D, 2D and 3D arrays. */
  Out_Ptr->Rd_rat = (In_Ptr->record.Rd_rat)
    ? AllocArray3D(0, nr - 1, 0, na - 1, 0, nt - 1, 1) : NULL;
  Out_Ptr->Rd_ra = (In_Ptr->record.Rd_ra)
    ? AllocArray2D(0, nr - 1, 0, na - 1, 1) : NULL;
  Out_Ptr->Rd_rt = (In_Ptr->record.Rd_rt)
    ? AllocArray2D(0, nr - 1, 0, nt - 1, 1) : NULL;
  Out_Ptr->Rd_at = (In_Ptr->record.Rd_at)
    ? AllocArray2D(0, na - 1, 0, nt - 1, 1) : NULL;
  Out_Ptr->Rd_r = (In_Ptr->record.Rd_r)
    ? AllocArray1D(0, nr - 1, 1) : NULL;
  Out_Ptr->Rd_a = (In_Ptr->record.Rd_a)
    ? AllocArray1D(0, na - 1, 1) : NULL;
  Out_Ptr->Rd_t = (In_Ptr->record.Rd_t)
    ? AllocArray1D(0, nt - 1, 1) : NULL;
 
  Out_Ptr->Td_rat = (In_Ptr->record.Td_rat)
    ? AllocArray3D(0, nr - 1, 0, na - 1, 0, nt - 1, 1) : NULL;
  Out_Ptr->Td_ra = (In_Ptr->record.Td_ra)
    ? AllocArray2D(0, nr - 1, 0, na - 1, 1) : NULL;
  Out_Ptr->Td_rt = (In_Ptr->record.Td_rt)
    ? AllocArray2D(0, nr - 1, 0, nt - 1, 1) : NULL;
  Out_Ptr->Td_at = (In_Ptr->record.Td_at)
    ? AllocArray2D(0, na - 1, 0, nt - 1, 1) : NULL;
  Out_Ptr->Td_r = (In_Ptr->record.Td_r)
    ? AllocArray1D(0, nr - 1, 1) : NULL;
  Out_Ptr->Td_a = (In_Ptr->record.Td_a)
    ? AllocArray1D(0, na - 1, 1) : NULL;
  Out_Ptr->Td_t = (In_Ptr->record.Td_t)
    ? AllocArray1D(0, nt - 1, 1) : NULL;

  Out_Ptr->A_rzt = (In_Ptr->record.A_rzt)
    ? AllocArray3D(0, nr - 1, 0, nz - 1, 0, nt - 1, 1) : NULL;
  Out_Ptr->Ab_zt = (In_Ptr->record.A_rzt)
    ? AllocArray2D(0, nz - 1, 0, nt - 1, 1) : NULL;
  Out_Ptr->A_rz = (In_Ptr->record.A_rz)
    ? AllocArray2D(0, nr - 1, 0, nz - 1, 1) : NULL;
  Out_Ptr->Ab_z = (In_Ptr->record.A_rz)
    ? AllocArray1D(0, nz - 1, 1) : NULL;
  Out_Ptr->A_zt = (In_Ptr->record.A_zt)
    ? AllocArray2D(0, nz - 1, 0, nt - 1, 1) : NULL;
  Out_Ptr->A_z = (In_Ptr->record.A_z)
    ? AllocArray1D(0, nz - 1, 1) : NULL;
  Out_Ptr->A_t = (In_Ptr->record.A_t)
    ? AllocArray1D(0, nt - 1, 1) : NULL;
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
    malloc((unsigned)(In_Ptr->num_media) * sizeof(LayerStru));
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

/*******************************************************************************
 *  Read which quantity is to be scored.
 ****/
Boolean
ReadRecordQ(FILE * Fp, InStru * In_Ptr)
{
  char        buf[STRLEN], string[STRLEN];
  char       *index;
  Boolean     error = 0;

  InitRecord(In_Ptr);
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

    if (strcmp(string, "Rd_r") == 0)
      In_Ptr->record.Rd_r = 1;
    else if (strcmp(string, "Rd_a") == 0)
      In_Ptr->record.Rd_a = 1;
    else if (strcmp(string, "Rd_ra") == 0)
      In_Ptr->record.Rd_ra = 1;
    else if (strcmp(string, "Rd_t") == 0)
      In_Ptr->record.Rd_t = 1;
    else if (strcmp(string, "Rd_rt") == 0)
      In_Ptr->record.Rd_rt = 1;
    else if (strcmp(string, "Rd_at") == 0)
      In_Ptr->record.Rd_at = 1;
    else if (strcmp(string, "Rd_rat") == 0)
      In_Ptr->record.Rd_rat = 1;
    else if (strcmp(string, "Td_r") == 0)
      In_Ptr->record.Td_r = 1;
    else if (strcmp(string, "Td_a") == 0)
      In_Ptr->record.Td_a = 1;
    else if (strcmp(string, "Td_ra") == 0)
      In_Ptr->record.Td_ra = 1;
    else if (strcmp(string, "Td_t") == 0)
      In_Ptr->record.Td_t = 1;
    else if (strcmp(string, "Td_rt") == 0)
      In_Ptr->record.Td_rt = 1;
    else if (strcmp(string, "Td_at") == 0)
      In_Ptr->record.Td_at = 1;
    else if (strcmp(string, "Td_rat") == 0)
      In_Ptr->record.Td_rat = 1;
    else if (strcmp(string, "A_z") == 0)
      In_Ptr->record.A_z = 1;
    else if (strcmp(string, "A_rz") == 0)
      In_Ptr->record.A_rz = 1;
    else if (strcmp(string, "A_t") == 0)
      In_Ptr->record.A_t = 1;
    else if (strcmp(string, "A_zt") == 0)
      In_Ptr->record.A_zt = 1;
    else if (strcmp(string, "A_rzt") == 0)
      In_Ptr->record.A_rzt = 1;
    else {
      printf("Unknown quantity: %s\n", string);
      error = 1;
    }
  }

  return (!error);
}

/*******************************************************************************
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

/*******************************************************************************
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

/*******************************************************************************
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

/*******************************************************************************
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
    malloc((unsigned) ((num_layers + 2) * sizeof(LayerStru)));
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
 *      If the z is on an interface between layers, the returned index
 *      will point to the upper layer.
 *      Index 0 is the top ambient medium and index num_layers+1 is the
 *      bottom one.
 ****/
Boolean
ZToLayerQ(double z, short *index, InStru * In_Ptr)
{
  short       i = 0;            /* index to layer. */
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
  if (!ReadRecordQ(Fp, In_Ptr))
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
 *	Initialize the ExtractStru.
 ****/
void 
InitExtractStru(ConvStru * Conv_Ptr)
{
  Conv_Ptr->extract.A_rz_t = 0;
  Conv_Ptr->extract.A_rz = 0;
  Conv_Ptr->extract.A_z_t = 0;
  Conv_Ptr->extract.A_z = 0;
  Conv_Ptr->extract.A_t_r_z = 0;
  Conv_Ptr->extract.A_t_z = 0;
  Conv_Ptr->extract.A_t = 0;
  Conv_Ptr->extract.A_l = 0;

  Conv_Ptr->extract.Rd_r_a_t = 0;                           
  Conv_Ptr->extract.Rd_r_a = 0;
  Conv_Ptr->extract.Rd_r_t = 0;
  Conv_Ptr->extract.Rd_r = 0;
 
  Conv_Ptr->extract.Rd_a_r_t = 0;
  Conv_Ptr->extract.Rd_a_r = 0;
  Conv_Ptr->extract.Rd_a_t = 0;
  Conv_Ptr->extract.Rd_a = 0;
 
  Conv_Ptr->extract.Rd_t_r_a = 0;
  Conv_Ptr->extract.Rd_t_r = 0;
  Conv_Ptr->extract.Rd_t_a = 0;
  Conv_Ptr->extract.Rd_t = 0;

  Conv_Ptr->extract.Td_r_a_t = 0;
  Conv_Ptr->extract.Td_r_a = 0;
  Conv_Ptr->extract.Td_r_t = 0;
  Conv_Ptr->extract.Td_r = 0;
 
  Conv_Ptr->extract.Td_a_r_t = 0;
  Conv_Ptr->extract.Td_a_r = 0;
  Conv_Ptr->extract.Td_a_t = 0;
  Conv_Ptr->extract.Td_a = 0;
 
  Conv_Ptr->extract.Td_t_r_a = 0;
  Conv_Ptr->extract.Td_t_r = 0;
  Conv_Ptr->extract.Td_t_a = 0;
  Conv_Ptr->extract.Td_t = 0;
}

/**************************************************************************
 *	Initialize the ExtractStru.
 ****/
void 
SetExtractStru(InStru * In_Ptr, ConvStru * Conv_Ptr)
{
  InitExtractStru(Conv_Ptr);

  if (In_Ptr->record.A_rzt) {
    Conv_Ptr->extract.A_rz_t = 1;
    Conv_Ptr->extract.A_rz = 1;
    Conv_Ptr->extract.A_z_t = 1;
    Conv_Ptr->extract.A_z = 1;
    Conv_Ptr->extract.A_t_r_z = 1;
    Conv_Ptr->extract.A_t_z = 1;
    Conv_Ptr->extract.A_t = 1;
    Conv_Ptr->extract.A_l = 1;
  }
  if (In_Ptr->record.A_rz) {  
    Conv_Ptr->extract.A_rz = 1;
    Conv_Ptr->extract.A_z = 1;
    Conv_Ptr->extract.A_l = 1;
  }
  if (In_Ptr->record.A_zt) {
    Conv_Ptr->extract.A_z_t = 1;
    Conv_Ptr->extract.A_z = 1;
    Conv_Ptr->extract.A_t_z = 1;
    Conv_Ptr->extract.A_t = 1;
    Conv_Ptr->extract.A_l = 1;
  }
  if (In_Ptr->record.A_z) { 
    Conv_Ptr->extract.A_z = 1;
    Conv_Ptr->extract.A_l = 1;
  }
  if (In_Ptr->record.A_t) 
    Conv_Ptr->extract.A_t = 1;
 
  if (In_Ptr->record.Rd_rat) {
    Conv_Ptr->extract.Rd_r_a_t = 1;
    Conv_Ptr->extract.Rd_r_a = 1;
    Conv_Ptr->extract.Rd_r_t = 1;
    Conv_Ptr->extract.Rd_r = 1;

    Conv_Ptr->extract.Rd_a_r_t = 1;
    Conv_Ptr->extract.Rd_a_r = 1;
    Conv_Ptr->extract.Rd_a_t = 1;
    Conv_Ptr->extract.Rd_a = 1;
   
    Conv_Ptr->extract.Rd_t_r_a = 1;
    Conv_Ptr->extract.Rd_t_r = 1;
    Conv_Ptr->extract.Rd_t_a = 1;
    Conv_Ptr->extract.Rd_t = 1;
  }
  if (In_Ptr->record.Rd_ra) {
    Conv_Ptr->extract.Rd_r_a = 1;
    Conv_Ptr->extract.Rd_r = 1;
    Conv_Ptr->extract.Rd_a_r = 1;
    Conv_Ptr->extract.Rd_a = 1;
  } 
  if (In_Ptr->record.Rd_rt) {
    Conv_Ptr->extract.Rd_r_t = 1;
    Conv_Ptr->extract.Rd_r = 1;
    Conv_Ptr->extract.Rd_t_r = 1;
    Conv_Ptr->extract.Rd_t = 1;
  }
  if (In_Ptr->record.Rd_at) {
    Conv_Ptr->extract.Rd_a_t = 1;
    Conv_Ptr->extract.Rd_a = 1;
    Conv_Ptr->extract.Rd_t_a = 1;
    Conv_Ptr->extract.Rd_t = 1;
  }
  if (In_Ptr->record.Rd_r) 
    Conv_Ptr->extract.Rd_r = 1;
  if (In_Ptr->record.Rd_a) 
    Conv_Ptr->extract.Rd_a = 1;
  if (In_Ptr->record.Rd_t) 
    Conv_Ptr->extract.Rd_t = 1;

  if (In_Ptr->record.Td_rat) {
    Conv_Ptr->extract.Td_r_a_t = 1;
    Conv_Ptr->extract.Td_r_a = 1;
    Conv_Ptr->extract.Td_r_t = 1;
    Conv_Ptr->extract.Td_r = 1;
 
    Conv_Ptr->extract.Td_a_r_t = 1;
    Conv_Ptr->extract.Td_a_r = 1;
    Conv_Ptr->extract.Td_a_t = 1;
    Conv_Ptr->extract.Td_a = 1;
  
    Conv_Ptr->extract.Td_t_r_a = 1;
    Conv_Ptr->extract.Td_t_r = 1;
    Conv_Ptr->extract.Td_t_a = 1;
    Conv_Ptr->extract.Td_t = 1;
  }
  if (In_Ptr->record.Td_ra) {
    Conv_Ptr->extract.Td_r_a = 1;
    Conv_Ptr->extract.Td_r = 1;
    Conv_Ptr->extract.Td_a_r = 1;
    Conv_Ptr->extract.Td_a = 1; 
  }  
  if (In_Ptr->record.Td_rt) {                              
    Conv_Ptr->extract.Td_r_t = 1;
    Conv_Ptr->extract.Td_r = 1;
    Conv_Ptr->extract.Td_t_r = 1;
    Conv_Ptr->extract.Td_t = 1;
  } 
  if (In_Ptr->record.Td_at) {
    Conv_Ptr->extract.Td_a_t = 1;
    Conv_Ptr->extract.Td_a = 1;
    Conv_Ptr->extract.Td_t_a = 1;
    Conv_Ptr->extract.Td_t = 1;
  }
  if (In_Ptr->record.Td_r)
    Conv_Ptr->extract.Td_r = 1;
  if (In_Ptr->record.Td_a)
    Conv_Ptr->extract.Td_a = 1;
  if (In_Ptr->record.Td_t)
    Conv_Ptr->extract.Td_t = 1;
} 

/**************************************************************************
 *  Mode = 0, read result back from a output file.
 *  Mode = 1, write result to a output file;
 ****/
Boolean
ReadOutMCV2(FILE * Fp, InStru * In_Ptr,
	 OutStru * Out_Ptr,
	 char Mode)
{
  FindDataLine(Fp);             /* skip the line of file version. */
  if (!ReadMediumListQ(Fp, In_Ptr))
    exit(1);
  ReadParam(Fp, In_Ptr);
  RestoreRandomStatus(Fp);
  InitOutputData(In_Ptr, Out_Ptr);
  ReadRAT(Fp, Out_Ptr);

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

  return (1);
}

/**************************************************************************
 *      Look for the key word, which is the 1st word in the line.
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
 
/**************************************************************************
 *      Return 1 if the kye is found, otherwise 0.
 ****/
Boolean
SeekKey(FILE * File_Ptr, char *Key)
{
  char       *line;
 
  do {
    if (!(line = FindDataLine(File_Ptr)))
      return (0);               /* Key not found. */
  } while (!FoundKeyWord(line, Key));
 
  return (1);                   /* Key found. */
}

/**************************************************************************
 *
 ****/
Boolean
ReadLayerParam(FILE * File_Ptr,
              short Num_Layers,
              LayerStru ** Layers_PP)
{
  char       *buf;
  short       i = 0;
  double      d, n, mua, mus, g;
  double      z = 0.0;          /* z coordinate of the current layer. */

  /* layer 0 and layer Num_Layers + 1 are for ambient. */
  *Layers_PP = (LayerStru *)
    malloc((unsigned) (Num_Layers + 2) * sizeof(LayerStru));
  if (!(*Layers_PP)) {
    printf("allocation failure in ReadLayerParam()");
    exit(1);
  }

  if (!(buf = FindDataLine(File_Ptr)))
    return (0);
  sscanf(buf, "%lf", &n);
  (*Layers_PP)[i].n = n;
  for (i = 1; i <= Num_Layers; i++) {
    if (!(buf = FindDataLine(File_Ptr)))
      return (0);
    sscanf(buf, "%lf%lf%lf%lf%lf", &n, &mua, &mus, &g, &d);
    (*Layers_PP)[i].n = n;
    (*Layers_PP)[i].mua = mua;
    (*Layers_PP)[i].mus = mus;
    (*Layers_PP)[i].g = g;
    (*Layers_PP)[i].z0 = z;
    z += d;
    (*Layers_PP)[i].z1 = z;
  }
  if (!(buf = FindDataLine(File_Ptr)))
    return (0);
  sscanf(buf, "%lf", &n);
  (*Layers_PP)[i].n = n;
 
  return (1);                   /* Data complete. */
}

/**************************************************************************
 *      Read in the input parameters which were used for Monte Carlo
 *      simulations.
 ****/
Boolean
ReadInParam(FILE * File_Ptr, InStru * In_Ptr)
{
  char       *buf;

  if (!SeekKey(File_Ptr, "InParm"))
    return (0);

  if (!ReadFnameFormatQ(File_Ptr, In_Ptr))
    return (0);
  /** read in number of photons. **/
  if (!(buf = FindDataLine(File_Ptr)))
    return (0);
  sscanf(buf, "%ld", &In_Ptr->num_photons);

  /** assign in Wth (critical weight). **/
  In_Ptr->Wth = 1E-4;

  /** read in dz, dr. **/
  if (!(buf = FindDataLine(File_Ptr)))
    return (0);
  sscanf(buf, "%lf%lf", &In_Ptr->dz, &In_Ptr->dr);

  /** read in nz, nr, na and compute da. **/
  if (!(buf = FindDataLine(File_Ptr)))
    return (0);
  sscanf(buf, "%hd%hd%hd", &In_Ptr->nz,
         &In_Ptr->nr, &In_Ptr->na);
  In_Ptr->da = 0.5 * PI / In_Ptr->na;
 
  /** read in number of layers. **/
  if (!(buf = FindDataLine(File_Ptr)))
    return (0);
  sscanf(buf, "%hd", &In_Ptr->num_layers);
 
  if (!ReadLayerParam(File_Ptr, In_Ptr->num_layers,
                     &In_Ptr->layerspecs))
    return (0);
       
  return (1);
}

/**************************************************************************
 *      Read reflectance, absorbed fraction, transmittance.
 ****/  
Boolean
ReadRATV1(FILE * File_Ptr, OutStru * Out_Ptr)
{
  char       *buf;

  if (!(buf = FindDataLine(File_Ptr)))
    return (0);
  sscanf(buf, "%lf", &(Out_Ptr->Rsp));

  if (!(buf = FindDataLine(File_Ptr)))
    return (0);
  sscanf(buf, "%lf", &(Out_Ptr->Rd));

  if (!(buf = FindDataLine(File_Ptr)))
    return (0);
  sscanf(buf, "%lf", &(Out_Ptr->A));

  if (!(buf = FindDataLine(File_Ptr)))
    return (0);
  sscanf(buf, "%lf", &(Out_Ptr->Td));

  return (1);
}

/**************************************************************************
 ****/  
Boolean
ReadA_rz(FILE * File_Ptr,
         short Ir0,             /* The lower index. */
         short Nr,
         short Nz,
         OutStru * Out_Ptr)
{
  short       iz, ir;
 
  for (iz = 0; iz < Nz; iz ++)
    Out_Ptr->Ab_z[iz] = 0.0;

  for (ir = Ir0; ir < Nr; ir++)
    for (iz = 0; iz < Nz; iz++)
      if (fscanf(File_Ptr, "%lf ",
                 &(Out_Ptr->A_rz[ir][iz])) == EOF)
        return (0);
  return (1);
}

/**************************************************************************
 ****/
Boolean
ReadRd_ra(FILE * File_Ptr,
          short Nr,
          short Na,
          OutStru * Out_Ptr)
{
  short       ir, ia;
 
  for (ir = 0; ir < Nr; ir++)
    for (ia = 0; ia < Na; ia++)
      if (fscanf(File_Ptr, "%lf ", &(Out_Ptr->Rd_ra[ir][ia])) == EOF)
        return (0);
  return (1);
}
 
/**************************************************************************
 ****/
Boolean
ReadTd_ra(FILE * File_Ptr,
          short Nr,
          short Na,
          OutStru * Out_Ptr)
{
  short       ir, ia;
 
  for (ir = 0; ir < Nr; ir++)
    for (ia = 0; ia < Na; ia++)
      if (fscanf(File_Ptr, "%lf ", &(Out_Ptr->Td_ra[ir][ia])) == EOF)
        return (0);
  return (1);
}

/**************************************************************************
 *      Read in the Monte Carlo output parameters for mco version 1.
 ****/
Boolean
ReadOutMCV1(FILE * File_Ptr,
            InStru * In_Ptr,
            OutStru * Out_Ptr)
{
  ReadInParam(File_Ptr, In_Ptr);
  InitRecord(In_Ptr);
  In_Ptr->record.A_rz = 1;
  In_Ptr->record.Rd_ra = 1;
  In_Ptr->record.Td_ra = 1;
  InitOutputData(In_Ptr, Out_Ptr);
 
  if (!SeekKey(File_Ptr, "RAT"))
    return (0);
  ReadRATV1(File_Ptr, Out_Ptr); /* refl.,absorption,transmission. */
 
  /* 2D arrays. */
  if (!SeekKey(File_Ptr, "A_rz"))
    return (0);
  ReadA_rz(File_Ptr, 0, In_Ptr->nr, In_Ptr->nz, Out_Ptr);

  if (!SeekKey(File_Ptr, "Rd_ra"))
    return (0);
  ReadRd_ra(File_Ptr, In_Ptr->nr, In_Ptr->na, Out_Ptr);
 
  if (!SeekKey(File_Ptr, "Tt_ra"))
    return (0);
  ReadTd_ra(File_Ptr, In_Ptr->nr, In_Ptr->na, Out_Ptr);
 
  puts("An mco file of version A1 has been read.");
  return (1);
}

/**************************************************************************
 *      After the Input parameters are read in, the parameters drc
 *      and nrc are initialized to dr and nr respectively.
 ****/
void
ReadMcoFile(InStru * In_Ptr,
            OutStru * Out_Ptr, ConvStru * Conv_Ptr)
{
  char        in_fname[STRLEN], buf[STRLEN];
  FILE       *infile;
  Boolean     data_complete;

  if ((infile = GetFile(in_fname)) == NULL)
    return;

  strcpy(buf, FindDataLine(infile));
  if (CheckFileVersionQ(buf, "A1")) {
    FreeData(In_Ptr, Out_Ptr);
    Conv_Ptr->datain = 0;
    Conv_Ptr->fversion = 1;
    data_complete = ReadOutMCV1(infile, In_Ptr, Out_Ptr);
  } else if (CheckFileVersionQ(buf, "mcmloA2.0")) {
    FreeData(In_Ptr, Out_Ptr);
    Conv_Ptr->datain = 0;
    Conv_Ptr->fversion = 2;
    data_complete = ReadOutMCV2(infile, In_Ptr, Out_Ptr, 0);
  } else {
    puts("...Input-file version not supported");
    return;
  }
 
  if (! data_complete) 
    FreeData(In_Ptr, Out_Ptr);
  else {
    Conv_Ptr->datain = 1;
    SetExtractStru(In_Ptr, Conv_Ptr);
    Conv_Ptr->nrc = In_Ptr->nr;
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nr = In_Ptr->nr;
    Conv_Ptr->dr = In_Ptr->dr;
  }

  fclose(infile);
}

