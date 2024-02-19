/*************************************************************************
 *  	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992.
 *
 *	Functions for file output.
 ****/

#include "conv.h"

void ConvResolution(ConvStru *);
void IsoPlot(double **, long, long, double, double);
void ConvM_rx(ConvStru *, double **);

/**************************************************************************
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
    if (Ext[0] != '\0')
      printf("Enter output filename with extension .%s (or . to quit): ", Ext);
    else
      printf("Enter output filename (or . to quit): ");
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

/**************************************************************************
 *	Return the index to the layer, where iz is in.
 *
 *	Use the center of box.
 ****/
short
IzToLayer(short Iz, InStru * In_Ptr)
{
  short       i = 0;		/* index to layer. */
  short       num_layers = In_Ptr->num_layers;
  double      dz = In_Ptr->dz;

  while ((Iz + 0.5) * dz > In_Ptr->layerspecs[i].z1
	 && i < num_layers)
    i++;

  return (i);
}

/**************************************************************************
 *      Compute A_rz from A_rzt.
 ****/
void
A_rzFromA_rzt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, iz, it;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;

  for (ir = 0; ir < nr; ir++)
    for (iz = 0; iz < nz; iz++) {
      Out_Ptr->A_rz[ir][iz] = 0.0;
      for (it = 0; it < nt; it++)
	Out_Ptr->A_rz[ir][iz] += Out_Ptr->A_rzt[ir][iz][it] * dt;
    }
}

/**************************************************************************
 *      Compute Ab_z from Ab_zt.
 ****/
void
Ab_zFromAb_zt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       iz, it;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;

  for (iz = 0; iz < nz; iz++) {
    Out_Ptr->Ab_z[iz] = 0.0;
    for (it = 0; it < nt; it++)
      Out_Ptr->Ab_z[iz] += Out_Ptr->Ab_zt[iz][it] * dt;
  }
}

#if 0
/**************************************************************************
 *      Compute A_zt from A_rzt.
 ****/
void
A_ztFromA_rzt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, iz, it;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;

  for (iz = 0; iz < nz; iz++)
    for (it = 0; it < nt; it++) {
      Out_Ptr->A_zt[iz][it] = 0.0;
      for (ir = 0; ir < nr; ir++)
	Out_Ptr->A_zt[iz][it] += Out_Ptr->A_rzt[ir][iz][it]
	  * 2.0 * PI * (ir + 0.5) * dr * dr;
      Out_Ptr->A_zt[iz][it] += Out_Ptr->Ab_zt[iz][it];
    }
}
#endif

/**************************************************************************
 *      Compute A_z from A_rzt.
 ****/
void
A_zFromA_rzt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, iz, it;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      dt = In_Ptr->dt;

  for (iz = 0; iz < nz; iz++) {
    Out_Ptr->A_z[iz] = 0.0;
    for (it = 0; it < nt; it++) {
      for (ir = 0; ir < nr; ir++)
	Out_Ptr->A_z[iz] += Out_Ptr->A_rzt[ir][iz][it] * dt
	  * 2.0 * PI * (ir + 0.5) * dr * dr;
      Out_Ptr->A_z[iz] += Out_Ptr->Ab_zt[iz][it] * dt;
    }
  }
}

/**************************************************************************
 *      Compute A_z from A_rz.
 ****/
void
A_zFromA_rz(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, iz;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  double      dr = In_Ptr->dr;

  for (iz = 0; iz < nz; iz++) {
    Out_Ptr->A_z[iz] = 0.0;
    for (ir = 0; ir < nr; ir++)
      Out_Ptr->A_z[iz] += Out_Ptr->A_rz[ir][iz]
	* 2.0 * PI * (ir + 0.5) * dr * dr;
    Out_Ptr->A_z[iz] += Out_Ptr->Ab_z[iz];
  }
}

/**************************************************************************
 *      Compute A_z from A_zt.
 ****/
void
A_zFromA_zt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       iz, it;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;

  for (iz = 0; iz < nz; iz++) {
    Out_Ptr->A_z[iz] = 0.0;
    for (it = 0; it < nt; it++)
      Out_Ptr->A_z[iz] += Out_Ptr->A_zt[iz][it] * dt;
  }
}

/**************************************************************************
 *      Compute A_t from A_rzt.
 ****/
void
A_tFromA_rzt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, iz, it;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      dz = In_Ptr->dz;

  for (it = 0; it < nt; it++) {
    Out_Ptr->A_t[it] = 0.0;
    for (iz = 0; iz < nz; iz++) {
      for (ir = 0; ir < nr; ir++)
	Out_Ptr->A_t[it] += Out_Ptr->A_rzt[ir][iz][it]
	  * 2.0 * PI * (ir + 0.5) * dr * dr * dz;
      Out_Ptr->A_t[it] += Out_Ptr->Ab_zt[iz][it] * dz;
    }
  }
}

/**************************************************************************
 *      Compute A_t from A_zt.
 ****/
void
A_tFromA_zt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       iz, it;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dz = In_Ptr->dz;

  for (it = 0; it < nt; it++) {
    Out_Ptr->A_t[it] = 0.0;
    for (iz = 0; iz < nz; iz++)
      Out_Ptr->A_t[it] += Out_Ptr->A_zt[iz][it] * dz;
  }
}

#if 0
/**************************************************************************
 *      Compute Rd_ra from Rd_rat.
 ****/
void
Rd_raFromRd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;

  for (ir = 0; ir < nr; ir++)
    for (ia = 0; ia < na; ia++) {
      Out_Ptr->Rd_ra[ir][ia] = 0.0;
      for (it = 0; it < nt; it++)
	Out_Ptr->Rd_ra[ir][ia] += Out_Ptr->Rd_rat[ir][ia][it] * dt;
    }
}

/**************************************************************************
 *      Compute Rd_rt from Rd_rat.
 ****/
void
Rd_rtFromRd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      da = In_Ptr->dt;

  for (ir = 0; ir < nr; ir++)
    for (it = 0; it < nt; it++) {
      Out_Ptr->Rd_rt[ir][it] = 0.0;
      for (ia = 0; ia < na; ia++)
	Out_Ptr->Rd_ra[ir][it] += Out_Ptr->Rd_rat[ir][ia][it]
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
    }
}

/**************************************************************************
 *      Compute Rd_at from Rd_rat.
 ****/
void
Rd_atFromRd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;

  for (ia = 0; ia < na; ia++)
    for (it = 0; it < nt; it++) {
      Out_Ptr->Rd_at[ia][it] = 0.0;
      for (ir = 0; ir < nr; ir++)
	Out_Ptr->Rd_at[ia][it] += Out_Ptr->Rd_rat[ir][ia][it]
	  * 2.0 * PI * (ir + 0.5) * dr * dr;
    }
}
#endif

/**************************************************************************
 *      Compute Rd_r from Rd_rat.
 ****/
void
Rd_rFromRd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;

  for (ir = 0; ir < nr; ir++) {
    Out_Ptr->Rd_r[ir] = 0.0;
    for (ia = 0; ia < na; ia++)
      for (it = 0; it < nt; it++)
	Out_Ptr->Rd_r[ir] += Out_Ptr->Rd_rat[ir][ia][it] * dt
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
  }
}

/**************************************************************************
 *      Compute Rd_r from Rd_ra.
 ****/
void
Rd_rFromRd_ra(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  double      da = In_Ptr->da;

  for (ir = 0; ir < nr; ir++) {
    Out_Ptr->Rd_r[ir] = 0.0;
    for (ia = 0; ia < na; ia++)
      Out_Ptr->Rd_r[ir] += Out_Ptr->Rd_ra[ir][ia]
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
  }
}

/**************************************************************************
 *      Compute Rd_r from Rd_rt.
 ****/
void
Rd_rFromRd_rt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, it;
  short       nr = In_Ptr->nr;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;

  for (ir = 0; ir < nr; ir++) {
    Out_Ptr->Rd_r[ir] = 0.0;
    for (it = 0; it < nt; it++)
      Out_Ptr->Rd_r[ir] += Out_Ptr->Rd_rt[ir][it] * dt;
  }
}

/**************************************************************************
 *      Compute Rd_a from Rd_rat.
 ****/
void
Rd_aFromRd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      dt = In_Ptr->dt;

  for (ia = 0; ia < na; ia++) {
    Out_Ptr->Rd_a[ia] = 0.0;
    for (ir = 0; ir < nr; ir++)
      for (it = 0; it < nt; it++)
	Out_Ptr->Rd_a[ia] += Out_Ptr->Rd_rat[ir][ia][it] * dt
	  * 2.0 * PI * (ir + 0.5) * dr * dr;
  }
}

/**************************************************************************
 *      Compute Rd_a from Rd_ra.
 ****/
void
Rd_aFromRd_ra(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  double      dr = In_Ptr->dr;

  for (ia = 0; ia < na; ia++) {
    Out_Ptr->Rd_a[ia] = 0.0;
    for (ir = 0; ir < nr; ir++)
      Out_Ptr->Rd_a[ia] += Out_Ptr->Rd_ra[ir][ia]
	* 2.0 * PI * (ir + 0.5) * dr * dr;
  }
}

/**************************************************************************
 *      Compute Rd_a from Rd_at.
 ****/
void
Rd_aFromRd_at(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ia, it;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;

  for (ia = 0; ia < na; ia++) {
    Out_Ptr->Rd_a[ia] = 0.0;
    for (it = 0; it < nt; it++)
      Out_Ptr->Rd_a[ia] += Out_Ptr->Rd_at[ia][it] * dt;
  }
}

/**************************************************************************
 *      Compute Rd_t from Rd_rat.
 ****/
void
Rd_tFromRd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;

  for (it = 0; it < nt; it++) {
    Out_Ptr->Rd_t[it] = 0.0;
    for (ir = 0; ir < nr; ir++)
      for (ia = 0; ia < na; ia++)
	Out_Ptr->Rd_t[it] += Out_Ptr->Rd_rat[ir][ia][it]
	  * 2.0 * PI * (ir + 0.5) * dr * dr
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
  }
}

/**************************************************************************
 *      Compute Rd_t from Rd_rt.
 ****/
void
Rd_tFromRd_rt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, it;
  short       nr = In_Ptr->nr;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;

  for (it = 0; it < nt; it++) {
    Out_Ptr->Rd_t[it] = 0.0;
    for (ir = 0; ir < nr; ir++)
      Out_Ptr->Rd_t[it] += Out_Ptr->Rd_rt[ir][it]
	* 2.0 * PI * (ir + 0.5) * dr * dr;
  }
}

/**************************************************************************
 *      Compute Rd_t from Rd_at.
 ****/
void
Rd_tFromRd_at(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ia, it;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      da = In_Ptr->da;

  for (it = 0; it < nt; it++) {
    Out_Ptr->Rd_t[it] = 0.0;
    for (ia = 0; ia < na; ia++)
      Out_Ptr->Rd_t[it] += Out_Ptr->Rd_at[ia][it]
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
  }
}

#if 0
/**************************************************************************
 *      Compute Td_ra from Td_rat.
 ****/
void 
Td_raFromTd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;

  for (ir = 0; ir < nr; ir++)
    for (ia = 0; ia < na; ia++) {
      Out_Ptr->Td_ra[ir][ia] = 0.0;
      for (it = 0; it < nt; it++)
        Out_Ptr->Td_ra[ir][ia] += Out_Ptr->Td_rat[ir][ia][it] * dt;
    }
}

/**************************************************************************
 *      Compute Td_rt from Td_rat.
 ****/
void 
Td_rtFromTd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      da = In_Ptr->dt;
 
  for (ir = 0; ir < nr; ir++)
    for (it = 0; it < nt; it++) {
      Out_Ptr->Td_rt[ir][it] = 0.0;
      for (ia = 0; ia < na; ia++)
        Out_Ptr->Td_ra[ir][it] += Out_Ptr->Td_rat[ir][ia][it]
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
    }
}

/**************************************************************************
 *      Compute Td_at from Td_rat.
 ****/
void
Td_atFromTd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
 
  for (ia = 0; ia < na; ia++)
    for (it = 0; it < nt; it++) {
      Out_Ptr->Td_at[ia][it] = 0.0;
      for (ir = 0; ir < nr; ir++)
        Out_Ptr->Td_at[ia][it] += Out_Ptr->Td_rat[ir][ia][it]
          * 2.0 * PI * (ir + 0.5) * dr * dr;
    }
}
#endif
 
/**************************************************************************
 *      Compute Td_r from Td_rat.
 ****/
void
Td_rFromTd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
 
  for (ir = 0; ir < nr; ir++) {
    Out_Ptr->Td_r[ir] = 0.0;
    for (ia = 0; ia < na; ia++)
      for (it = 0; it < nt; it++)
        Out_Ptr->Td_r[ir] += Out_Ptr->Td_rat[ir][ia][it] * dt
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
  }
}

/**************************************************************************
 *      Compute Td_r from Td_ra.
 ****/
void
Td_rFromTd_ra(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  double      da = In_Ptr->da;
 
  for (ir = 0; ir < nr; ir++) {
    Out_Ptr->Td_r[ir] = 0.0;
    for (ia = 0; ia < na; ia++)
      Out_Ptr->Td_r[ir] += Out_Ptr->Td_ra[ir][ia]
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
  }
}
 
/**************************************************************************
 *      Compute Td_r from Td_rt.
 ****/
void
Td_rFromTd_rt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, it;
  short       nr = In_Ptr->nr;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;
 
  for (ir = 0; ir < nr; ir++) {
    Out_Ptr->Td_r[ir] = 0.0;
    for (it = 0; it < nt; it++)
      Out_Ptr->Td_r[ir] += Out_Ptr->Td_rt[ir][it] * dt;
  }
}

/**************************************************************************
 *      Compute Td_a from Td_rat.
 ****/
void
Td_aFromTd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      dt = In_Ptr->dt;
 
  for (ia = 0; ia < na; ia++) {
    Out_Ptr->Td_a[ia] = 0.0;
    for (ir = 0; ir < nr; ir++)
      for (it = 0; it < nt; it++)
        Out_Ptr->Td_a[ia] += Out_Ptr->Td_rat[ir][ia][it] * dt
          * 2.0 * PI * (ir + 0.5) * dr * dr;
  }
}
 
/**************************************************************************
 *      Compute Td_a from Td_ra.
 ****/
void
Td_aFromTd_ra(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  double      dr = In_Ptr->dr;
 
  for (ia = 0; ia < na; ia++) {
    Out_Ptr->Td_a[ia] = 0.0;
    for (ir = 0; ir < nr; ir++)
      Out_Ptr->Td_a[ia] += Out_Ptr->Td_ra[ir][ia]
        * 2.0 * PI * (ir + 0.5) * dr * dr;
  }
}

/**************************************************************************
 *      Compute Td_a from Td_at.
 ****/
void 
Td_aFromTd_at(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ia, it;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;
 
  for (ia = 0; ia < na; ia++) {
    Out_Ptr->Td_a[ia] = 0.0;
    for (it = 0; it < nt; it++)
      Out_Ptr->Td_a[ia] += Out_Ptr->Td_at[ia][it] * dt;
  }
}
 
/**************************************************************************
 *      Compute Td_t from Td_rat.
 ****/
void 
Td_tFromTd_rat(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
 
  for (it = 0; it < nt; it++) {
    Out_Ptr->Td_t[it] = 0.0;
    for (ir = 0; ir < nr; ir++)
      for (ia = 0; ia < na; ia++)
        Out_Ptr->Td_t[it] += Out_Ptr->Td_rat[ir][ia][it]
          * 2.0 * PI * (ir + 0.5) * dr * dr
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
  }
}

/**************************************************************************
 *      Compute Td_t from Td_rt.
 ****/
void 
Td_tFromTd_rt(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ir, it;
  short       nr = In_Ptr->nr;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
 
  for (it = 0; it < nt; it++) {
    Out_Ptr->Td_t[it] = 0.0;
    for (ir = 0; ir < nr; ir++)
      Out_Ptr->Td_t[it] += Out_Ptr->Td_rt[ir][it]
        * 2.0 * PI * (ir + 0.5) * dr * dr;
  }
}
 
/**************************************************************************
 *      Compute Td_t from Td_at.
 ****/
void
Td_tFromTd_at(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       ia, it;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      da = In_Ptr->da;
 
  for (it = 0; it < nt; it++) {
    Out_Ptr->Td_t[it] = 0.0;
    for (ia = 0; ia < na; ia++)
      Out_Ptr->Td_t[it] += Out_Ptr->Td_at[ia][it]
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
  }
}

/**************************************************************************
 *      Extract reflectance, absorption, transmission.
 ****/
void
ExtractRAT(InStru * In_Ptr,
	   OutStru * Out_Ptr,
	   ConvStru * Conv_Ptr)
{
  FILE       *file;
  char        fname[STRLEN];

  if (!Conv_Ptr->datain) {
    puts("...No data to extract");
    return;
  }

  strcpy(fname, "RAT");
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  fprintf(file,	"# RAT -- Reflectance, absorption, transmittance.\n");

  if (Conv_Ptr->fversion == 1) {
    fprintf(file, "%-12.4G \t#Rsp: Specular reflectance.\n",
	    Out_Ptr->Rsp);
    fprintf(file, "%-12.4G \t#Rd:  Diffuse reflectance.\n",
	    Out_Ptr->Rd);
    fprintf(file, "%-12.4G \t#A:   Absorption.\n",
	    Out_Ptr->A);
    fprintf(file, "%-12.4G \t#Td:  Diffuse transmission.\n",
	    Out_Ptr->Td);
  } else {
    fprintf(file, "# %-10s\t%-12s  Rel Err\n", "Average", "Standard Err");
    fprintf(file, "%-12.4G \t\t\t\t#Rsp: Specular reflectance.\n",
	    Out_Ptr->Rsp);
    fprintf(file, "%-12.4G \t%-12.4G %6.2f%%\t#Rb: Ballistic reflectance.\n",
	    Out_Ptr->Rb, Out_Ptr->Rbe,
	    (Out_Ptr->Rb) ? Out_Ptr->Rbe / Out_Ptr->Rb * 100 : 0);
    fprintf(file, "%-12.4G \t%-12.4G %6.2f%%\t#Rd: Diffuse reflectance.\n",
	    Out_Ptr->Rd, Out_Ptr->Rde,
	    (Out_Ptr->Rd) ? Out_Ptr->Rde / Out_Ptr->Rd * 100 : 0);
    fprintf(file, "%-12.4G \t%-12.4G %6.2f%%\t#A:  Absorbed fraction.\n",
	    Out_Ptr->A, Out_Ptr->Ae,
	    (Out_Ptr->A) ? Out_Ptr->Ae / Out_Ptr->A * 100 : 0);
    fprintf(file, "%-12.4G \t%-12.4G %6.2f%%\t#Tb: Ballistic transmittance.\n",
	    Out_Ptr->Tb, Out_Ptr->Tbe,
	    (Out_Ptr->Tb) ? Out_Ptr->Tbe / Out_Ptr->Tb * 100 : 0);
    fprintf(file, "%-12.4G \t%-12.4G %6.2f%%\t#Td: Diffuse transmittance.\n",
	    Out_Ptr->Td, Out_Ptr->Tde,
	    (Out_Ptr->Td) ? Out_Ptr->Tde / Out_Ptr->Td * 100 : 0);
  }

  fclose(file);
}

/**************************************************************************
 *      Compute the mixed mua if the ith box is in two layers.
 ****/
double
IzToMua(short iz, InStru * In_Ptr) 
{
  double  tz, bz;      /* top and bottom z of ith box. */
  double  mua;
  LayerStru thelayer;  /* the layer which the ith box is in. */
  short   ilayer;      /* index to thelayer. */
   
  tz = iz * In_Ptr->dz;
  bz = (iz + 1) * In_Ptr->dz; 
  ilayer = IzToLayer(iz, In_Ptr);
  thelayer = In_Ptr->layerspecs[ilayer];  

  if (((tz < thelayer.z0) && (bz > thelayer.z0))
      || ((tz < thelayer.z1) && (bz > thelayer.z1))) {
    if ((tz < thelayer.z1) && (bz > thelayer.z1))  {
      ilayer ++;
      thelayer = In_Ptr->layerspecs[ilayer];
    }

    mua = (thelayer.z0 - tz) * In_Ptr->layerspecs[ilayer-1].mua;
    mua += (bz - thelayer.z0) * In_Ptr->layerspecs[ilayer].mua;
    mua /= In_Ptr->dz;

  } else       /* the ith box is only in one layer. */
    mua = In_Ptr->layerspecs[ilayer].mua;

  return mua;
}

/**************************************************************************
 *	Write beam parameter into a file.
 ****/  
void
SpecifyBeam(FILE * Fp, ConvStru * Conv_Ptr)
{
  if (Conv_Ptr->beam.type == flat)
    fprintf(Fp, "# Flat beam with: P = %f [J]; R = %f [cm].\n",
	    Conv_Ptr->beam.P, Conv_Ptr->beam.R);
  else if (Conv_Ptr->beam.type == gaussian)
    fprintf(Fp, "# Gaussian beam with: P = %f [J]; R = %f [cm].\n",
	    Conv_Ptr->beam.P, Conv_Ptr->beam.R);
  else if (Conv_Ptr->beam.type == arbitrary)
    fprintf(Fp, "# Arbitrary beam specified in file: %s\n", 
	   Conv_Ptr->beam.pro_fname);
}

/**************************************************************************
 *	Get a fixed t, z, r or a in range: 0.0 - Max.
 *      mode = 0, get a fixed t;  mode = 1, get a fixed z;
 *      mode = 2, get a fixed r;  mode = 3, get a fixed a;
 ****/  
double
GetFixedTZRA(double Max, char mode)
{
  double number;
  char buf[STRLEN];

  do {
    if (mode == 0)
      printf("Specify a fixed t (0 - %f ps): ", Max);
    else if (mode == 1)
      printf("Specify a fixed z (0 - %f cm): ", Max);
    else if (mode == 2)
      printf("Specify a fixed r (0 - %f cm): ", Max);
    else 
      printf("Specify a fixed a (0 - %f sr): ", Max);

    gets(buf);
  } while (sscanf(buf, "%lf", &number) != 1 || number<0.0 || number>Max);

  return number;
}

/**************************************************************************
 *      Write A_rz in 3 columns format: r, z, A_rzc.
 *      mode = 0, write A_rz@t; mode = 1, write F_rz@t.
 *      mode = 2, write A_rz; mode = 3, write F_rz.
 ****/
void
WriteA_rz3Columns(InStru *In_Ptr,
                      double **A_rz,
                      double Pt, char mode)
{
  short       ir, iz;
  double      r;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  double      dr = In_Ptr->dr;
  double      dz = In_Ptr->dz;
  FILE       *file;
  char        fname[STRLEN];
 
  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Arz@t");
  if (mode == 1) strcpy(fname, "Frz@t");
#endif
  if (mode == 2) strcpy(fname, "Arz");
  if (mode == 3) strcpy(fname, "Frz");
  if ((file = GetWriteFile(fname)) == NULL)
    return;
 
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0)      fprintf(file, "# A_rz at t = %G ps\n", Pt);
  else if (mode == 1) fprintf(file, "# F_rz at t = %G ps\n", Pt);
  else if (mode ==2)  fprintf(file, "# A_rz\n");
  else                fprintf(file, "# F_rz\n");
  if (mode == 0 || mode == 2)
    fprintf(file, "# %-10s\t%-12s\t%-12s\n", "r [cm]", "z [cm]",
           "A [1/(cm3 ps)]");
  if (mode == 1 || mode == 3)
    fprintf(file, "# %-10s\t%-12s\t%-12s\n", "r [cm]", "z [cm]",
           "F [1/(cm2 ps)]");
 
  for (ir = 0; ir < nr-1; ir ++) {
    r = ((ir + 0.5) + 1 / (12 * (ir + 0.5))) * dr;
    for (iz = 0; iz < nz-1; iz ++)
      fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n", r, (iz+0.5)*dz,
             A_rz[ir][iz]);
  }
 
  fclose(file);
}

/**************************************************************************
 *      Extract A_rz @ fix time.
 *      mode = 0, extract A_rz@t; mode = 1, extract F_rz@t.
 ****/
void
ExtractA_rz_t(InStru * In_Ptr, OutStru * Out_Ptr, char mode)
{
  short       ir, iz, it;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      dz = In_Ptr->dz;
  double      dt = In_Ptr->dt;
  double      pt, mua;
  char        buf[STRLEN];

  if ((Out_Ptr->A_rz = AllocArray2D(0, nr-1, 0, nz-1, 0)) == NULL) {
    printf("Allocating A_rz failed.\n");
    return;
  }
     
  pt = GetFixedTZRA((nt-1)*dt, 0);
  it = pt/dt;
  for (ir = 0; ir < nr; ir ++)
    for (iz = 0; iz < nz; iz ++) { 
      Out_Ptr->A_rz[ir][iz] = Out_Ptr->A_rzt[ir][iz][it]; 
      if (mode == 1) { 
	mua = IzToMua(iz, In_Ptr);
	if (mua > 0.0)
	  Out_Ptr->A_rz[ir][iz] /= mua;
      }
    }

  do {
    printf("Which output format (3 = 3 columns / c = contour)? ");
    gets(buf);
  } while ((toupper(buf[0]) != '3') && (toupper(buf[0]) != 'C'));

  if (toupper(buf[0]) == '3')
    WriteA_rz3Columns(In_Ptr, Out_Ptr->A_rz, pt, mode);
  else                    /* contour */
    IsoPlot(Out_Ptr->A_rz, nr-2, nz-2, dr, dz);

  FreeArray2D(Out_Ptr->A_rz, 0, nr-1, 0, nz-1);
}

/**************************************************************************
 *      Write convolved A_rz in 3 columns format: r, z, A_rzc.
 *      mode = 0, write A_rz@tc; mode = 1, write F_rz@tc.
 *      mode = 1, write A_rzc; mode = 3, write F_rzc.
 ****/
void
WriteConvA_rz3Columns(InStru *In_Ptr,
                      ConvStru *Conv_Ptr,
                      double **A_rz,
                      double Pt, char mode)
{
  short       ir, iz;
  double      r;
  short       nr = Conv_Ptr->nrc;
  short       nz = In_Ptr->nz;
  double      dr = Conv_Ptr->drc;
  double      dz = In_Ptr->dz;
  FILE       *file;
  char        fname[STRLEN];
 
  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0)      strcpy(fname, "Arz@tc");
  else if (mode == 1) strcpy(fname, "Frz@tc");
  else if (mode == 2) strcpy(fname, "Arzc");
  else                strcpy(fname, "Frzc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;
 
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  SpecifyBeam(file, Conv_Ptr);
  if (mode == 0)      fprintf(file, "# A_rz at t = %G ps\n", Pt);
  else if (mode == 1) fprintf(file, "# F_rz at t = %G ps\n", Pt);
  else if (mode ==2)  fprintf(file, "# A_rzc\n");
  else                fprintf(file, "# F_rzc\n");
  if (mode == 0 || mode == 2)
    fprintf(file, "# %-10s\t%-12s\t%-12s\n", "r [cm]", "z [cm]",
           "A [J/(cm3 ps)]");
  if (mode == 1 || mode == 3)
    fprintf(file, "# %-10s\t%-12s\t%-12s\n", "r [cm]", "z [cm]",
           "F [J/(cm2 ps)]");
 
  for (ir = 0; ir < nr-1; ir ++) {
    r = ((ir + 0.5) + 1 / (12 * (ir + 0.5))) * dr;
    for (iz = 0; iz < nz-1; iz ++)
      fprintf(file, "%-12.4E\t%-12.4E\t%-12.4E\n", r, (iz+0.5)*dz,
             A_rz[ir][iz]);
  }
 
  fclose(file);
}

/**************************************************************************
 *      Convolve and  extract A_rz @ a fix time.
 *      mode = 0, extract A_rz@tc; mode = 1, extract F_rz@tc.
 ****/
void
ExtractConvA_rz_t(InStru * In_Ptr, OutStru * Out_Ptr,
                  ConvStru * Conv_Ptr, char mode)
{
  short       ir, iz, it;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dz = In_Ptr->dz;
  double      dt = In_Ptr->dt;
  double      pt, mua;
  char        buf[STRLEN];
  double     **A_rzc;
 
  ConvResolution(Conv_Ptr);
  if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, nz-1, 0)) == NULL)
      || ((Conv_Ptr->Mb_x = AllocArray1D(0, nz-1, 0)) == NULL)
      || ((A_rzc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, nz-1, 0)) == NULL)) {
    printf("No enough memory to convolve A_rz@t.\n");
    return;
  }
 
  pt = GetFixedTZRA((nt-1)*dt, 0);
  it = pt/dt;
  for (iz = 0; iz < nz; iz ++) { 
    Conv_Ptr->Mb_x[iz] = Out_Ptr->Ab_zt[iz][it];
    if (mode == 1) {
      mua = IzToMua(iz, In_Ptr);
      if (mua > 0.0) 
        Conv_Ptr->Mb_x[iz] /= mua;
    }
    for (ir = 0; ir < nr; ir ++) {
      Conv_Ptr->M_rx[ir][iz] = Out_Ptr->A_rzt[ir][iz][it];
      if (mode == 1 && mua > 0.0) 
          Out_Ptr->A_rz[ir][iz] /= mua;
    }  
  }
 
  Conv_Ptr->nxc = nz;
  Conv_Ptr->dxc = dz;
  ConvM_rx(Conv_Ptr, A_rzc);
 
  do {
    printf("Which output format (3 = 3 columns / c = contour)? ");
    gets(buf);
  } while ((toupper(buf[0]) != '3') && (toupper(buf[0]) != 'C'));
 
  if (toupper(buf[0]) == '3')
    WriteConvA_rz3Columns(In_Ptr, Conv_Ptr, A_rzc, pt, mode);
  else
    IsoPlot(A_rzc, Conv_Ptr->nrc-2, nz-2, Conv_Ptr->drc, dz);

  Conv_Ptr->drc = In_Ptr->dr;
  Conv_Ptr->nrc = In_Ptr->nr;
  FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, nz-1);
  FreeArray2D(A_rzc, 0, Conv_Ptr->nrc-1, 0, nz-1);
}

/**************************************************************************
 *      Extract A_rz.
 *      mode = 0, extract A_rz; mode = 1, extract F_rz.
 ****/
void
ExtractA_rz(InStru * In_Ptr, 
	    OutStru * Out_Ptr, 
	    char mode)
{
  short       ir, iz; 
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  double      dr = In_Ptr->dr;
  double      dz = In_Ptr->dz;
  double      mua;
  char        buf[STRLEN];

  if (In_Ptr->record.A_rzt) {
    if ((Out_Ptr->A_rz = AllocArray2D(0, nr-1, 0, nz-1, 0)) == NULL) {
      printf("Allocating A_rz failed.\n");
      return;
    }
    A_rzFromA_rzt(In_Ptr, Out_Ptr);
  }

  if (mode == 1) 
    for (ir = 0; ir < nr-1; ir++) 
      for (iz = 0; iz < nz-1; iz++) {
	mua = IzToMua(iz, In_Ptr);
	if (mua > 0.0)
	  Out_Ptr->A_rz[ir][iz] /= mua;
	}

  do {
    printf("Which output format (3 = 3 columns / c = contour)? ");
    gets(buf);
  } while ((toupper(buf[0]) != '3') && (toupper(buf[0]) != 'C'));

  if (toupper(buf[0]) == '3') 
    WriteA_rz3Columns(In_Ptr, Out_Ptr->A_rz, 0, mode+2);
  else
    IsoPlot(Out_Ptr->A_rz, nr-2, nz-2, dr, dz);

  if (mode == 1)
    for (ir = 0; ir < nr-1; ir++)
      for (iz = 0; iz < nz-1; iz++) {
	mua = IzToMua(iz, In_Ptr);
	if (mua > 0.0)
	  Out_Ptr->A_rz[ir][iz] *= mua;
      }

  if (! In_Ptr->record.A_rz)
    FreeArray2D(Out_Ptr->A_rz, 0, nr-1, 0, nz-1);
}

/**************************************************************************
 *      Convolve and  extract A_rz.
 *      mode = 0, extract A_rzc; mode = 1, extract F_rzc.
 ****/
void
ExtractConvA_rz(InStru *In_Ptr, OutStru *Out_Ptr,
                ConvStru *Conv_Ptr, char mode)
{
  short       ir, iz;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  double      dz = In_Ptr->dz;
  double      mua;
  char        buf[STRLEN];
  double     **A_rzc;

  ConvResolution(Conv_Ptr);
  if ((! In_Ptr->record.A_rz &&
      ((Out_Ptr->A_rz = AllocArray2D(0, nr-1, 0, nz-1, 0)) == NULL)
      || ((Out_Ptr->Ab_z = AllocArray1D(0, nz-1, 0)) == NULL))
     || ((A_rzc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, nz-1, 0)) == NULL)) {
    printf("No enough memory to convolve A_rz.\n");
    return;
  }

  if (In_Ptr->record.A_rzt) {
    A_rzFromA_rzt(In_Ptr, Out_Ptr);
    Ab_zFromAb_zt(In_Ptr, Out_Ptr);
  }

  Conv_Ptr->M_rx = Out_Ptr->A_rz;
  Conv_Ptr->Mb_x = Out_Ptr->Ab_z;
  if (mode == 1)
    for (iz = 0; iz < nz; iz ++) {
      mua = IzToMua(iz, In_Ptr);
      if (mua > 0)
	Conv_Ptr->Mb_x[iz] /= mua;
      for (ir = 0; ir < nr; ir ++)
	if (mua > 0)
	  Conv_Ptr->M_rx[ir][iz] /= mua;
    }	

  Conv_Ptr->nxc = nz;
  Conv_Ptr->dxc = dz;
  ConvM_rx(Conv_Ptr, A_rzc);

  do {
    printf("Which output format (3 = 3 columns / c = contour)? ");
    gets(buf);
  } while ((toupper(buf[0]) != '3') && (toupper(buf[0]) != 'C'));
 
  if (toupper(buf[0]) == '3')
    WriteConvA_rz3Columns(In_Ptr, Conv_Ptr, A_rzc, 0, mode+2);
  else
    IsoPlot(A_rzc, Conv_Ptr->nrc-2, nz-2, Conv_Ptr->drc, dz);
 
  Conv_Ptr->drc = In_Ptr->dr;
  Conv_Ptr->nrc = In_Ptr->nr;
  if (In_Ptr->record.A_rzt) {
    FreeArray2D(Out_Ptr->A_rz, 0, nr-1, 0, nz-1);
    FreeArray1D(Out_Ptr->Ab_z, 0, nz-1);
  }
  FreeArray2D(A_rzc, 0, Conv_Ptr->nrc-1, 0, nz-1);
}

/**************************************************************************
 *      2 numbers each line: z, A_z@t.
 *      mode = 0, extract A_z@t; mode = 1, extract F_z@t.
 *      mode = 2, extract A_z@tc; mode = 3, extract F_z@tc.
 ****/
void
ExtractA_z_t(InStru * In_Ptr, OutStru * Out_Ptr, 
	     ConvStru * Conv_Ptr, char mode)
{
  short       ir, iz, it;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      dz = In_Ptr->dz;
  double      dt = In_Ptr->dt;
  double      pt, mua, temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0)      strcpy(fname, "Az@t");
  else if (mode == 1) strcpy(fname, "Fz@t");
  else if (mode == 2) strcpy(fname, "Az@tc");
  else                strcpy(fname, "Fz@tc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pt = GetFixedTZRA((nt-1)*dt, 0);
  it = pt / dt;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# A_z at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "z [cm]", "A [1/(cm ps)]");
  } else if (mode == 1) {
    fprintf(file, "# F_z at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "z [cm]", "F [1/ps]");
  } else if (mode == 2) { 
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# A_zc at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "z [cm]", "A [J/(cm ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# F_zc at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "z [cm]", "F [J/ps]");
  }

  if (In_Ptr->record.A_rzt) {
    for (iz = 0; iz < nz-1; iz++) {
      temp = 0.0;
      for (ir = 0; ir < nr; ir++)
	temp += Out_Ptr->A_rzt[ir][iz][it] 
	     * 2.0 * PI * (ir + 0.5) * dr * dr;
      temp += Out_Ptr->Ab_zt[iz][it];
      if (mode == 1 || mode == 3) {
	mua = IzToMua(iz, In_Ptr);
        if (mua > 0.0)
          temp /= mua;
      }
      if (mode == 2 || mode == 3)
	temp *= Conv_Ptr->beam.P;
      fprintf(file, "%-12.4E\t%-12.4E\n", (iz+0.5)*dz, temp);
    }

  } else {			/* In_Ptr->record.A_zt == 1 */
    for (iz = 0; iz < nz-1; iz++) {
      temp = Out_Ptr->A_zt[iz][it];
      if (mode == 1 || mode == 3) {
	mua = IzToMua(iz, In_Ptr);
        if (mua > 0.0)
          temp /= mua;
      }
      if (mode == 2 || mode == 3)
	temp *= Conv_Ptr->beam.P;
      fprintf(file, "%-12.4E\t %-12.4E\n", (iz+0.5)*dz, temp);
    }
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: z, A_z.
 *      mode = 0, extract A_z; mode = 1; extract F_z.
 *      mode = 2, extract A_zc; mode = 3, extract F_zc.
 ****/
void
ExtractA_z(InStru * In_Ptr, OutStru * Out_Ptr, 
	      ConvStru *Conv_Ptr, char mode)
{
  short       iz;
  short       nz = In_Ptr->nz;
  double      dz = In_Ptr->dz;
  double      mua, temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
  if (mode == 0)      strcpy(fname, "Az");
  else if (mode == 1) strcpy(fname, "Fz");
  else if (mode == 2) strcpy(fname, "Azc");
  else                strcpy(fname, "Fzc");

  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# A_z\n");
    fprintf(file, "# %-10s\t%-12s\n", "z [cm]", "A [1/cm]");
  } else if (mode == 1) {
    fprintf(file, "# F_z\n");
    fprintf(file, "# %-10s\t%-12s\n", "z [cm]", "F [-]");
  } else if (mode == 2) {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# A_zc\n");
    fprintf(file, "# %-10s\t%-12s\n", "z [cm]", "A [J/cm]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# F_zc\n");
    fprintf(file, "# %-10s\t%-12s\n", "z [cm]", "F [J]");
  }

  if (!In_Ptr->record.A_z)
    if ((Out_Ptr->A_z = AllocArray1D(0, nz - 1, 0)) == NULL) {
      printf("Allocating A_z failed.\n");
      return;
    }
  if (In_Ptr->record.A_rzt)
    A_zFromA_rzt(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.A_rz)
    A_zFromA_rz(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.A_zt)
    A_zFromA_zt(In_Ptr, Out_Ptr);

  for (iz = 0; iz < nz-1; iz++) {
    temp = Out_Ptr->A_z[iz];
    if (mode == 1 || mode == 3) {
      mua = IzToMua(iz, In_Ptr);
      if (mua > 0.0)
	temp /= mua;
    }
    if (mode == 2 || mode == 3) 
      temp *= Conv_Ptr->beam.P;
    fprintf(file, "%-12.4E\t%-12.4E\n", (iz+0.5)*dz, temp);
  }

  if (!In_Ptr->record.A_z)
    FreeArray1D(Out_Ptr->A_z, 0, nz - 1);
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, A_t@r@z.
 *      mode = 0, extract A_t@r@z; mode = 1, extract F_t@r@z.
 ****/
void
ExtractA_t_r_z(InStru * In_Ptr, 
	       OutStru * Out_Ptr, 
	       char mode)
{
  short       ir, iz, it;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      dz = In_Ptr->dz;
  double      dt = In_Ptr->dt;
  double      pr, pz, mua, temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "At@r@z");
  else           strcpy(fname, "Ft@r@z");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pr = GetFixedTZRA((nr-1)*dr, 2);
  pz = GetFixedTZRA((nz-1)*dz, 1);
  ir = pr / dr; iz = pz / dz;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# A_t at r = %G cm and z = %G cm\n", pr, pz);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "A [1/(cm3 ps)]");
  } else {
    fprintf(file, "# F_t at r = %G cm and z = %G cm\n", pr, pz);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "F [1/(cm2 ps)]");
  }

  mua = IzToMua(iz, In_Ptr);
  for (it = 0; it < nt; it++)
    if (mode == 0) 
      fprintf(file, "%-12.4E\t%-12.4E\n", 
	      (it+0.5)*dt, Out_Ptr->A_rzt[ir][iz][it]);
    else {
      temp = Out_Ptr->A_rzt[ir][iz][it];
      if (mua > 0.0)
	temp /= mua;
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, temp);
    }
 
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, A_t@r@zc.
 *      mode = 0, extract A_t@r@zc; mode = 1, extract F_t@r@zc.
 ****/
void
ExtractConvA_t_r_z(InStru * In_Ptr, OutStru * Out_Ptr,
              ConvStru *Conv_Ptr, char mode)
{
  short       ir, iz, it;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      dz = In_Ptr->dz;
  double      dt = In_Ptr->dt;
  double      pr, pz, mua;
  FILE       *file;
  char        fname[STRLEN];
  double     **A_rtc;
 
  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "At@r@zc");
  else           strcpy(fname, "Ft@r@zc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;
 
  pr = GetFixedTZRA((nr-1)*dr, 2);
  pz = GetFixedTZRA((nz-1)*dz, 1);
  ir = pr / dr; iz = pz / dz;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  SpecifyBeam(file, Conv_Ptr);
  if (mode == 0) {
    fprintf(file, "# A_tc at r = %G cm and z = %G cm\n", pr, pz);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "A [J/(cm3 ps)]");
  } else {
    fprintf(file, "# F_tc at r = %G cm and z = %G cm\n", pr, pz);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "F [J/(cm2 ps)]");
  } 

  Conv_Ptr->drc = 2 * pr;
  Conv_Ptr->nrc = 1;
  if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, nt-1, 0)) == NULL)
    || ((Conv_Ptr->Mb_x = AllocArray1D(0, nt-1, 0)) == NULL)
    || ((A_rtc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, nt-1, 0)) == NULL)) {
    printf("No enough memory to convolve A_t@r@z.\n");
    return;
  }
 
  for (it = 0; it < nt; it ++) 
    for (ir = 0; ir < nr; ir ++) 
      Conv_Ptr->M_rx[ir][it] = 0.0;
      for (iz = 0; iz < nz; iz ++)
	Conv_Ptr->M_rx[ir][it] += Out_Ptr->A_rzt[ir][iz][it] * dz;

  for (it = 0; it < nt; it ++) {
    Conv_Ptr->Mb_x[it] = 0.0;
    for (iz = 0; iz < nz; iz ++) 
      Conv_Ptr->Mb_x[it] += Out_Ptr->Ab_zt[iz][it] * dz;
  }

  mua = IzToMua(iz, In_Ptr);
  if (mode == 3)
    for (it = 0; it < nt; it ++) {
      if (mua > 0)
	Conv_Ptr->Mb_x[it] /= mua;
      for (ir = 0; ir < nr; ir ++)
	if (mua > 0)
	  Conv_Ptr->M_rx[ir][it] /= mua;
    }
 
  Conv_Ptr->dxc = dt;
  Conv_Ptr->nxc = nt;
  ConvM_rx(Conv_Ptr, A_rtc);
  for (it = 0; it < nt; it ++)
    fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, A_rtc[0][it]);

  Conv_Ptr->drc = In_Ptr->dr;
  Conv_Ptr->nrc = In_Ptr->nr;
  FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, nt-1);
  FreeArray1D(Conv_Ptr->Mb_x, 0, nt-1);
  FreeArray2D(A_rtc, 0, Conv_Ptr->nrc-1, 0, nt-1);
 
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, A_t@z.
 *      mode = 0, extract A_t@z; mode = 1, extract F_t@z.
 *      mode = 2, extract A_t@zc; mode = 3, extract F_t@zc.
 ****/
void
ExtractA_t_z(InStru * In_Ptr, OutStru * Out_Ptr, 
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, iz, it;
  short       nr = In_Ptr->nr;
  short       nz = In_Ptr->nz;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      dz = In_Ptr->dz;
  double      dt = In_Ptr->dt;
  double      pz, temp, mua;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0)      strcpy(fname, "At@z");
  else if (mode == 1) strcpy(fname, "Ft@z");
  else if (mode == 2) strcpy(fname, "At@zc");
  else                strcpy(fname, "Ft@zc");
#endif

  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pz = GetFixedTZRA((nz-1)*dz, 1);
  iz = pz / dz;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) { 
    fprintf(file, "# A_t at z = %G cm\n", pz);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "A [1/(cm ps)]");
  } else if (mode == 1) {
    fprintf(file, "# F_t at z = %G cm\n", pz);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "F [1/ps]");
  } else if (mode == 2) { 
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# A_tc at z = %G cm\n", pz);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "A [J/(cm ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# F_tc at z = %G cm\n", pz);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "F [J/ps]");
  }

  if (In_Ptr->record.A_rzt) {
    for (it = 0; it < nt; it++) {
      temp = 0.0;
      for (ir = 0; ir < nr-1; ir++)
	temp += Out_Ptr->A_rzt[ir][iz][it] 
	     * 2.0 * PI * (ir + 0.5) * dr * dr;
      temp += Out_Ptr->Ab_zt[iz][it];
      if (mode == 1 || mode == 3) {
	mua = IzToMua(iz, In_Ptr);
	if (mua > 0.0)
	  temp /= mua;
      }
      if (mode == 2 || mode == 3)
	temp *= Conv_Ptr->beam.P;
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, temp);
    }

  } else			/* In_Ptr->record.A_zt == 1 */
    for (it = 0; it < nt; it++) {
      temp = Out_Ptr->A_zt[iz][it];
      if (mode == 1 || mode == 3) {
	mua = IzToMua(iz, In_Ptr);
        if (mua > 0.0)
          temp /= mua;
      } 
      if (mode == 2 || mode == 3)
        temp *= Conv_Ptr->beam.P;
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, temp);
    }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, A_t.
 *      mode = 0, extract A_t; mode = 1, extract A_tc.
 ****/
void
ExtractA_t(InStru * In_Ptr, OutStru * Out_Ptr,
	   ConvStru *Conv_Ptr, char mode)
{
  short       it;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;
  double      temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
  if (mode == 0) strcpy(fname, "At");
  else           strcpy (fname, "Atc");
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# A_t\n");
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "A [1/ps]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# A_tc\n");
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "A [J/ps]");
  }

  if (!In_Ptr->record.A_t)
    if ((Out_Ptr->A_t = AllocArray1D(0, nt - 1, 0)) == NULL) {
      printf("Allocating A_t failed.\n");
      return;
    }
  
  if (In_Ptr->record.A_rzt)
    A_tFromA_rzt(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.A_zt)
    A_tFromA_zt(In_Ptr, Out_Ptr);

  for (it = 0; it < nt; it++) { 
    temp = Out_Ptr->A_t[it];
    if (mode == 1)
        temp *= Conv_Ptr->beam.P;
    fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, temp);
  }

  if (!In_Ptr->record.A_t)
    FreeArray1D(Out_Ptr->A_t, 0, nt - 1);
  fclose(file);
}

/**************************************************************************
 *      Write absorption as a function of layer.
 *      2 numbers each line: layer, A_l.
 ****/
void
ExtractA_l(InStru * In_Ptr, OutStru * Out_Ptr)
{
  short       iz, ilayer;
  short       nz = In_Ptr->nz;
  short       num_layers = In_Ptr->num_layers;
  double      dz = In_Ptr->dz;
  double     *A_layer;
  FILE       *file;
  char        fname[STRLEN];

  strcpy(fname, "A_l");
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  fprintf(file, "# A_l\n");
  fprintf(file, "# %-10s\t%-12s\n", "layer", "A [-]");

  if ((A_layer = AllocArray1D(1, num_layers, 0)) == NULL) {
    printf("Allocating A_layer failed.\n");
    return;
  }
  for (ilayer = 1; ilayer <= num_layers; ilayer++)
    A_layer[ilayer] = 0.0;

  if (!In_Ptr->record.A_z)
    if ((Out_Ptr->A_z = AllocArray1D(0, nz - 1, 0)) == NULL) {
      printf("Allocating A_z failed.\n");
      return;
    }
  if (In_Ptr->record.A_rzt)
    A_zFromA_rzt(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.A_rz)
    A_zFromA_rz(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.A_zt)
    A_zFromA_zt(In_Ptr, Out_Ptr);

  for (iz = 0; iz < nz; iz++) {
    ilayer = IzToLayer(iz, In_Ptr);
    A_layer[ilayer] += Out_Ptr->A_z[iz] * dz;
  }

  for (ilayer = 1; ilayer <= num_layers; ilayer++) 
    fprintf(file, "%-14d%-12.4E\n", ilayer, A_layer[ilayer]);

  if (!In_Ptr->record.A_z)
    FreeArray1D(Out_Ptr->A_z, 0, nz - 1);
  FreeArray1D(A_layer, 1, num_layers);

  fclose(file);
}

/*************************************************************************
 *      2 numbers each line: r, Rd_r@a@t.
 *      mode = 0, extract Rd_r@a@t; mode = 1, extract Rd_r@a@tc.
 ****/
void
ExtractRd_r_a_t(InStru * In_Ptr, OutStru * Out_Ptr,
              ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pa, pt;
  FILE       *file;
  char        fname[STRLEN];
  double     **Rd_rxc;
 
  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Rdr@a@t");
  else           strcpy(fname, "Rdr@a@tc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;
 
  pa = GetFixedTZRA((na-1)*da, 3);
  pt = GetFixedTZRA((nt-1)*dt, 0);
  ia = pa / da; it = pt / dt;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Rd_r at a = %G sr and t = %G ps\n", pa, pt);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Rd [1/(cm2 sr ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_rc at a = %G sr and t = %G ps\n", pa, pt);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Rd [J/(cm2 sr ps)]");
  }
 
  if (mode == 0)
    for (ir = 0; ir < nr-1; ir++)
      fprintf(file, "%-12.4E\t%-12.4E\n",
              (ir+0.5) * dr, Out_Ptr->Rd_rat[ir][ia][it]);
 
  else {           /* convolve Rd_r@a@t. */
    ConvResolution(Conv_Ptr);
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, 0, 0)) == NULL)
      || ((Rd_rxc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, 0, 0)) == NULL)) {
      printf("No enough memory to convolve Rd_r@a@t.\n");
      return;
    }
    Conv_Ptr->Mb_x = NULL;
    for (ir = 0; ir < nr-1; ir ++)
      Conv_Ptr->M_rx[ir][0] = Out_Ptr->Rd_rat[ir][ia][it];
 
    Conv_Ptr->nxc = 1;
    Conv_Ptr->dxc = dr;
    ConvM_rx(Conv_Ptr, Rd_rxc);
    for (ir = 0; ir < nr-1; ir ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Rd_rxc[0][ir]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, 0);
    FreeArray2D(Rd_rxc, 0, Conv_Ptr->nrc-1, 0, 0);
  }
 
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: r, Rd_r@a.
 *      mode = 0, extract Rd_r@a; mode = 1, extract Rd_r@ac.
 ****/
void
ExtractRd_r_a(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pa;
  FILE       *file;
  char        fname[STRLEN];
  double       **Rd_rxc;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Rdr@a");
  else strcpy(fname, "Rdr@ac");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pa = GetFixedTZRA((na-1)*da, 3);
  ia = pa / da;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Rd_r at a = %G sr\n", pa);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Rd [1/(cm2 sr)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_rc at a = %G sr\n", pa);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Rd [J/(cm2 sr)]");
  }

  if (mode == 1) ConvResolution(Conv_Ptr);

  if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, 0, 0)) == NULL)
      || ((mode == 1) && 
         (Rd_rxc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, 0, 0)) == NULL)) {
      printf("No enough memory to extract Rd_r@a@t.\n");
      return;
    }

  if (In_Ptr->record.Rd_rat) {
    for (ir = 0; ir < nr; ir++) {
      Conv_Ptr->M_rx[ir][0] = 0.0;
      for (it = 0; it < nt; it++)
	Conv_Ptr->M_rx[ir][0] += Out_Ptr->Rd_rat[ir][ia][it] * dt;
    }

  } else 			/* In_Ptr->record.Rd_ra == 1 */
    for (ir = 0; ir < nr; ir++)
      Conv_Ptr->M_rx[ir][0] = Out_Ptr->Rd_ra[ir][ia];
    
  if (mode == 0) 
    fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Conv_Ptr->M_rx[ir][0]);
  else {
    Conv_Ptr->Mb_x = NULL; 
    Conv_Ptr->nxc = 1;
    ConvM_rx(Conv_Ptr, Rd_rxc);
    for (ir = 0; ir < nr-1; ir ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Rd_rxc[0][ir]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Rd_rxc, 0, Conv_Ptr->nrc-1, 0, 0);
  }

  FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, 0);
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: r, Rd_r@t.
 *      mode = 0, extract Rd_r@t; mode = 1, extract Rd_r@tc.
 ****/
void
ExtractRd_r_t(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pt;
  FILE       *file;
  char        fname[STRLEN];
  double     **Rd_rxc;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Rdr@t");
  else strcpy(fname, "Rdr@tc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pt = GetFixedTZRA((nt-1)*dt, 0);
  it = pt / dt;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode ==0) {
    fprintf(file, "# Rd_r at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Rd [1/(cm2 ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_rc at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Rd [J/(cm2 ps)]");
  }

  if (mode == 1) ConvResolution(Conv_Ptr);

  if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, 0, 0)) == NULL)
      || ((mode == 1) &&
         (Rd_rxc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, 0, 0)) == NULL)) {
      printf("No enough memory to extract Rd_r@t.\n");
      return;
    }

  if (In_Ptr->record.Rd_rat) {
    for (ir = 0; ir < nr-1; ir++) {
      Conv_Ptr->M_rx[ir][0] = 0.0;
      for (ia = 0; ia < na; ia++)
        Conv_Ptr->M_rx[ir][0] += Out_Ptr->Rd_rat[ir][ia][it]
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
    }

  } else 			/* In_Ptr->record.Rd_rt == 1 */
    for (ir = 0; ir < nr-1; ir++)
      Conv_Ptr->M_rx[ir][0] = Out_Ptr->Rd_rt[ir][it];
 
  if (mode == 0) 
    fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Conv_Ptr->M_rx[ir][0]);
  else {
    Conv_Ptr->Mb_x = NULL;
    Conv_Ptr->nxc = 1;
    ConvM_rx(Conv_Ptr, Rd_rxc);
    for (ir = 0; ir < nr-1; ir ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Rd_rxc[0][ir]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Rd_rxc, 0, Conv_Ptr->nrc-1, 0, 0);
  }
 
  FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, 0);
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: r, Rd_r.
 *      mode = 0, extract Rd_r; mode = 1, extract Rd_rc.
 ****/
void
ExtractRd_r(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir;
  short       nr = In_Ptr->nr;
  double      dr = In_Ptr->dr;
  FILE       *file;
  char        fname[STRLEN];
  double     **Rd_rxc;

  fname[0] = '\0';
  if (mode == 0) strcpy(fname, "Rdr");
#if IBMPC == 0
  if (mode == 1) strcpy(fname, "Rdrc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) { 
    fprintf(file, "# Rd_r\n");
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Rd [1/cm2]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_rc\n");
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Rd [J/cm2]");
  }

  if (!In_Ptr->record.Rd_r)
    if ((Out_Ptr->Rd_r = AllocArray1D(0, nr - 1, 0)) == NULL) {
      printf("Allocating Rd_r failed.\n");
      return;
    }
  if (In_Ptr->record.Rd_rat)
    Rd_rFromRd_rat(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Rd_ra)
    Rd_rFromRd_ra(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Rd_rt)
    Rd_rFromRd_rt(In_Ptr, Out_Ptr);

  if (mode == 0)
    for (ir = 0; ir < nr-1; ir++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5) * dr, Out_Ptr->Rd_r[ir]);
  else {
    ConvResolution(Conv_Ptr);

    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, 0, 0)) == NULL)
      || ((Rd_rxc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, 0, 0)) == NULL)) {
      printf("No enough memory to extract Rd_r.\n");
      return;
    }

    for (ir = 0; ir < nr-1; ir ++)
      Conv_Ptr->M_rx[ir][0] = Out_Ptr->Rd_r[ir];

    Conv_Ptr->Mb_x = NULL;
    Conv_Ptr->nxc = 1;
    ConvM_rx(Conv_Ptr, Rd_rxc);
    for (ir = 0; ir < Conv_Ptr->nrc; ir ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Rd_rxc[ir][0]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Rd_rxc, 0, Conv_Ptr->nrc-1, 0, 0);
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, 0);
  }
 
  if (!In_Ptr->record.Rd_r)
    FreeArray1D(Out_Ptr->Rd_r, 0, nr - 1);
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: a, Rd_a@r@t.
 *      mode = 0, extract Rd_a@r@t; mode = 1, extract Rd_a@r@tc.
 ****/
void
ExtractRd_a_r_t(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pr, pt;
  FILE       *file;
  char        fname[STRLEN];
  double    **Rd_rac;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Rda@r@t");
  else strcpy(fname, "Rda@r@tc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pr = GetFixedTZRA((nr-1)*dr, 2);
  pt = GetFixedTZRA((nt-1)*dt, 0);
  ir = pr / dr; it = pt / dt;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Rd_a at r = %G cm and t = %G ps\n", pr, pt);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Rd [1/(cm2 sr ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_ac at r = %G cm and t = %G ps\n", pr, pt);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Rd [J/(cm2 sr ps)]");
  }

  if (mode == 0)
    for (ia = 0; ia < na; ia++)
      fprintf(file, "%-12.4E\t%-12.4E\n", 
	     (ia+0.5)*da, Out_Ptr->Rd_rat[ir][ia][it]);
  else {
    Conv_Ptr->drc = 2 * pr;
    Conv_Ptr->nrc = 1;
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, na-1, 0)) == NULL)
      || ((Rd_rac = AllocArray2D(0, Conv_Ptr->nrc-1, 0, na-1, 0)) == NULL)) {
      printf("No enough memory to convolve Rd_a@r@t.\n");
      return;
    }

    Conv_Ptr->Mb_x = NULL;
    for (ir = 0; ir < nr-1; ir ++)
      for (ia = 0; ia < na; ia ++) {
        Conv_Ptr->M_rx[ir][ia] = 0.0;
        for (it = 0; it < nt; it ++) 
          Conv_Ptr->M_rx[ir][ia] += Out_Ptr->Rd_rat[ir][ia][it] * dt; 
      }
    
    Conv_Ptr->dxc = da;
    Conv_Ptr->nxc = na;
    ConvM_rx(Conv_Ptr, Rd_rac);
    for(ia = 0; ia < na; ia ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ia+0.5)*da, Rd_rac[0][ia]);

    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, na-1);
    FreeArray2D(Rd_rac, 0, 0, 0, na-1);
  }
 
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: a, Rd_a@r.
 *      mode = 0, extract Rd_a@r; mode = 1, extract Rd_a@rc.
 ****/
void
ExtractRd_a_r(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pr, temp;
  FILE       *file;
  char        fname[STRLEN];
  double      ** Rd_rac;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Rda@r");
  else strcpy(fname, "Rda@rc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pr = GetFixedTZRA((nr-1)*dr, 2);
  ir = pr / dr;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Rd_a at r = %G cm\n", pr);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Rd [1/(cm2 sr)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_ac at r = %G cm\n", pr);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Rd [J/(cm2 sr)]");
  }

  if (mode == 0) {
    if (In_Ptr->record.Rd_rat) {
      for (ia = 0; ia < na; ia++) {
	temp = 0.0;
	for (it = 0; it < nt; it++)
	  temp += Out_Ptr->Rd_rat[ir][ia][it] * dt;
	fprintf(file, "%-12.4E\t%-12.4E\n", (ia+0.5) * da, temp);
      }

    } else {			/* In_Ptr->record.Rd_ra == 1 */
      for (ia = 0; ia < na; ia++)
	fprintf(file, "%-12.4E\t%-12.4E\n", 
		(ia+0.5) * da, Out_Ptr->Rd_ra[ir][ia]);
    }
  } else {
    Conv_Ptr->drc = 2 * pr;
    Conv_Ptr->nrc = 1;
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, na-1, 0)) == NULL)
      || ((Rd_rac = AllocArray2D(0, Conv_Ptr->nrc-1, 0, na-1, 0)) == NULL)) {
      printf("No enough memory to convolve Rd_a@r@t.\n");
      return;
    }
   
    Conv_Ptr->Mb_x = NULL;
    if (In_Ptr->record.Rd_rat) {
      for (ir = 0; ir < nr-1; ir ++)
	for (ia = 0; ia < na; ia++) {
	  Conv_Ptr->M_rx[ir][ia] = 0.0;
	  for (it = 0; it < nt; it++)
	    Conv_Ptr->M_rx[ir][ia] += Out_Ptr->Rd_rat[ir][ia][it] * dt;
	}
    } else {			/* In_Ptr->record.Rd_ra == 1 */
      for (ir = 0; ir < nr-1; ir++)
	for (ia = 0; ia < na; ia++)
	  Conv_Ptr->M_rx[ir][ia] += Out_Ptr->Rd_ra[ir][ia];
    }

    Conv_Ptr->dxc = da;
    Conv_Ptr->nxc = na;
    ConvM_rx(Conv_Ptr, Rd_rac);
    for(ia = 0; ia < na; ia ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ia+0.5)*da, Rd_rac[0][ia]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, na-1);
    FreeArray2D(Rd_rac, 0, 0, 0, na-1);
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: a, Rd_a@t.
 *      mode = 0, extract Rd_a@t; mode = 1, extract Rd_a@tc.
 ****/
void
ExtractRd_a_t(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pt, temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Rda@t");
  else           strcpy(fname, "Rda@tc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pt = GetFixedTZRA((nt-1)*dt, 0);
  it = pt / dt;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Rd_a at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Rd [1/(sr ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_ac at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Rd [J/(sr ps)]");
  }

  if (In_Ptr->record.Rd_rat) {
    for (ia = 0; ia < na; ia++) {
      temp = 0.0;
      for (ir = 0; ir < nr-1; ir++)
	temp += Out_Ptr->Rd_rat[ir][ia][it]
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
      if (mode == 1) temp *= Conv_Ptr->beam.P;
      fprintf(file, "%-12.4E\t%-12.4E\n", (ia+0.5) * da, temp);
    }

  } else {			/* In_Ptr->record.Rd_at == 1 */
    for (ia = 0; ir < na; ia++) {
      temp = Out_Ptr->Rd_at[ia][it];
      fprintf(file, "%-12.4E\t%-12.4E\n", (ia+0.5)*da, temp);
    }
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: a, Rd_a.
 *      mode = 0, extract Rd_a; mode = 1, extract Rd_ac.
 ****/
void
ExtractRd_a(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ia;
  short       na = In_Ptr->na;
  double      da = In_Ptr->da;
  double      temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
  if (mode == 0) strcpy(fname, "Rda");
#if IBMPC == 0
  if (mode == 1) strcpy(fname, "Rdac");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Rd_a\n");
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Rd [1/sr]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_ac\n");
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Rd [J/sr]");
  }

  if (!In_Ptr->record.Rd_a)
    if ((Out_Ptr->Rd_a = AllocArray1D(0, na - 1, 0)) == NULL) {
      printf("Allocating Rd_a failed.\n");
      return;
    }
  if (In_Ptr->record.Rd_rat)
    Rd_aFromRd_rat(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Rd_ra)
    Rd_aFromRd_ra(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Rd_at)
    Rd_aFromRd_at(In_Ptr, Out_Ptr);

  for (ia = 0; ia < na; ia++) {
    temp = Out_Ptr->Rd_a[ia];
    if (mode == 1) temp *= Conv_Ptr->beam.P;
    fprintf(file, "%-12.4E\t%-12.4E\n", (ia + 0.5) * da, temp);
  }

  if (!In_Ptr->record.Rd_a)
    FreeArray1D(Out_Ptr->Rd_a, 0, na - 1);
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, Rd_t@r@a.
 *      mode = 0, extract Rd_t@r@a; mode = 1, extract Rd_t@r@ac.
 ****/
void
ExtractRd_t_r_a(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pr, pa;
  FILE       *file;
  char        fname[STRLEN];
  double    **Rd_rtc;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Rdt@r@a");
  else strcpy(fname, "Rdt@r@ac");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pr = GetFixedTZRA((nr-1)*dr, 2);
  pa = GetFixedTZRA((na-1)*da, 3);
  ir = pr / dr; ia = pa / da;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Rd_a at r = %G cm and a = %G sr\n", pr, pa);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Rd [1/(cm2 sr ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_ac at r = %G cm and a = %G sr\n", pr, pa);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Rd [J/(cm2 sr ps)]");
  }

  if (mode == 0)
    for (it = 0; it < nt; it++)
      fprintf(file, "%-12.4E\t%-12.4E\n", 
	      (it+0.5)*dt, Out_Ptr->Rd_rat[ir][ia][it]);
  else {
    Conv_Ptr->drc = 2 * pr;
    Conv_Ptr->nrc = 1;
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, na-1, 0)) == NULL)
      || ((Rd_rtc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, nt-1, 0)) == NULL)) {
      printf("No enough memory to convolve Rd_a@r@t.\n");
      return;
    }  
       
    Conv_Ptr->Mb_x = NULL;
    for (ir = 0; ir < nr-1; ir ++)
      for (ia = 0; ia < na; ia ++) {
        Conv_Ptr->M_rx[ir][ia] = 0.0;
        for (it = 0; it < nt; it ++)
          Conv_Ptr->M_rx[ir][ia] += Out_Ptr->Rd_rat[ir][ia][it] * dt;
      }
       
    Conv_Ptr->dxc = dt;
    Conv_Ptr->nxc = nt;
    ConvM_rx(Conv_Ptr, Rd_rtc);
    for(it = 0; it < nt; it ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, Rd_rtc[0][it]);
       
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, nt-1);
    FreeArray2D(Rd_rtc, 0, 0, 0, nt-1);
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, Rd_t@r.
 *      mode = 0, extract Rd_t@r; mode = 1, extract Rd_t@rc.
 ****/
void
ExtractRd_t_r(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pr, temp;
  FILE       *file;
  char        fname[STRLEN];
  double    **Rd_rtc;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Rdt@r");
  else strcpy(fname, "Rdt@rc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pr = GetFixedTZRA((nr-1)*dr, 2);
  ir = pr / dr;
  fprintf(file, "# %s extracted from %s\n\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Rd_t at r = %G cm\n", pr);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Rd [1/(cm2 ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_tc at r = %G cm\n", pr);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Rd [J/(cm2 ps)]");
  }

  if (mode == 0) {
    if (In_Ptr->record.Rd_rat) {
      for (it = 0; it < nt; it++) {
	temp = 0.0;
	for (ia = 0; ia < na; ia++)
	  temp += Out_Ptr->Rd_rat[ir][ia][it]
	  * PI * sin(2.0 * (ia + 0.5) * da) * da;
	fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5) * dt, temp);
      }
    } else {			/* In_Ptr->record.Rd_rt == 1 */
      for (it = 0; it < nt; it++)
	fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, Out_Ptr->Rd_rt[ir][it]);
    }

  } else {
    Conv_Ptr->drc = 2 * pr;
    Conv_Ptr->nrc = 1;
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, nt-1, 0)) == NULL)
      || ((Rd_rtc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, nt-1, 0)) == NULL)) {
      printf("No enough memory to convolve Rd_t@rh.\n");
      return;
    }
 
    Conv_Ptr->Mb_x = NULL;
    if (In_Ptr->record.Rd_rat) {
      for (ir = 0; ir < nr-1; ir ++)
        for (it = 0; it < nt; it++) {
          Conv_Ptr->M_rx[ir][it] = 0.0;
	  for (ia = 0; ia < na; ia++)
            Conv_Ptr->M_rx[ir][ia] += Out_Ptr->Rd_rat[ir][ia][it]
	      * PI * sin(2.0 * (ia + 0.5) * da) * da;
        }
    } else {                    /* In_Ptr->record.Rd_rt == 1 */
      for (ir = 0; ir < nr-1; ir++)
        for (it = 0; it < nt; it++)
          Conv_Ptr->M_rx[ir][it] += Out_Ptr->Rd_rt[ir][it];
    }
       
    Conv_Ptr->dxc = dt;
    Conv_Ptr->nxc = nt;
    ConvM_rx(Conv_Ptr, Rd_rtc);
    for(it = 0; it < nt; it ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, Rd_rtc[0][it]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, nt-1);
    FreeArray2D(Rd_rtc, 0, 0, 0, nt-1);
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, Rd_t@a.
 *      mode = 0, extract Rd_t@a; mode = 1, extract Rd_t@ac.
 ****/
void
ExtractRd_t_a(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pa, temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Rdt@a");
  else strcpy(fname, "Rdt@ac");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pa = GetFixedTZRA((na-1)*da, 3);
  ia = pa / da;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Rd_t at a = %G sr\n", pa);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Rd [1/(sr ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_tc at a = %G sr\n", pa);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Rd [J/(sr ps)]");
  }

  if (In_Ptr->record.Rd_rat) {
    for (it = 0; it < nt; it++) {
      temp = 0.0;
      for (ir = 0; ir < nr-1; ir++)
	temp += Out_Ptr->Rd_rat[ir][ia][it]
	  * 2.0 * PI * (ir + 0.5) * dr * dr;
      if (mode == 1) temp *= Conv_Ptr->beam.P;
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5) * dt, temp);
    }

  } else {			/* In_Ptr->record.Rd_at == 1 */
    for (it = 0; it < nt; it++) {
      temp = Out_Ptr->Rd_at[ia][it];
      if (mode == 1) temp *= Conv_Ptr->beam.P;
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, Out_Ptr->Rd_at[ia][it]);
    }
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, Rd_t.
 *      mode = 0, extract Rd_t; mode = 1, extract Rd_tc.
 ****/
void
ExtractRd_t(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       it;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;
  double      temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
  if (mode == 0) strcpy(fname, "Rdt");
#if IBMPC == 0
  if (mode == 1) strcpy(fname, "Rdtc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Rd_t\n");
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Rd [1/ps]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Rd_tc\n");
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Rd [J/ps]");
  }

  if (!In_Ptr->record.Rd_t)
    if ((Out_Ptr->Rd_t = AllocArray1D(0, nt - 1, 0)) == NULL) {
      printf("Allocating Rd_t failed.\n");
      return;
    }
  if (In_Ptr->record.Rd_rat)
    Rd_tFromRd_rat(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Rd_rt)
    Rd_tFromRd_rt(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Rd_at)
    Rd_tFromRd_at(In_Ptr, Out_Ptr);

  for (it = 0; it < nt; it++) {
    temp = Out_Ptr->Rd_t[it];
    if (mode == 1)
      temp *= Conv_Ptr->beam.P;
    fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5) * dt, temp);
  }

  if (!In_Ptr->record.Rd_t)
    FreeArray1D(Out_Ptr->Rd_t, 0, nt - 1);
  fclose(file);
}

/*************************************************************************
 *      2 numbers each line: r, Td_r@a@t.
 *      mode = 0, extract Td_r@a@t; mode = 1, extract Td_r@a@tc.
 ****/
void
ExtractTd_r_a_t(InStru * In_Ptr, OutStru * Out_Ptr,
              ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pa, pt;
  FILE       *file;
  char        fname[STRLEN];
  double     **Td_rxc;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Tdr@a@t");
  else strcpy(fname, "Tdr@a@tc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pa = GetFixedTZRA((na-1)*da, 3);
  pt = GetFixedTZRA((nt-1)*dt, 0);
  ia = pa / da; it = pt / dt;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_r at a = %G sr and t = %G ps\n", pa, pt);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Td [1/(cm2 sr ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_rc at a = %G sr and t = %G ps\n", pa, pt);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Td [J/(cm2 sr ps)]");
  }
 
  if (mode == 0)
    for (ir = 0; ir < nr-1; ir++)
      fprintf(file, "%-12.4E\t%-12.4E\n",
              (ir+0.5) * dr, Out_Ptr->Td_rat[ir][ia][it]);
 
  else {           /* convolve Td_r@a@t. */
    ConvResolution(Conv_Ptr);
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, 0, 0)) == NULL)
      || ((Td_rxc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, 0, 0)) == NULL)) {
      printf("No enough memory to convolve Td_r@a@t.\n");
      return;
    }
    Conv_Ptr->Mb_x = NULL;
    for (ir = 0; ir < nr-1; ir ++)
      Conv_Ptr->M_rx[ir][0] = Out_Ptr->Td_rat[ir][ia][it];
 
    Conv_Ptr->nxc = 1;
    Conv_Ptr->dxc = dr;
    ConvM_rx(Conv_Ptr, Td_rxc);
    for (ir = 0; ir < nr-1; ir ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Td_rxc[0][ir]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, 0);
    FreeArray2D(Td_rxc, 0, Conv_Ptr->nrc-1, 0, 0);
  }
 
  fclose(file);
}
 
/**************************************************************************
 *      2 numbers each line: r, Td_r@a.
 *      mode = 0, extract Td_r@a; mode = 1, extract Td_r@ac.
 ****/
void
ExtractTd_r_a(InStru * In_Ptr, OutStru * Out_Ptr,
              ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pa;
  FILE       *file;
  char        fname[STRLEN];
  double       **Td_rxc;
 
  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Tdr@a");
  else strcpy(fname, "Tdr@ac");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;
 
  pa = GetFixedTZRA((na-1)*da, 3);
  ia = pa / da;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_r at a = %G sr\n", pa);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Td [1/(cm2 sr)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_rc at a = %G sr\n", pa);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Td [J/(cm2 sr)]");
  }
 
  if (mode == 1) ConvResolution(Conv_Ptr);
  
  if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, 0, 0)) == NULL)
      || ((mode == 1) &&
         (Td_rxc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, 0, 0)) == NULL)) {
      printf("No enough memory to extract Td_r@a@t.\n");
      return;
    }
 
  if (In_Ptr->record.Td_rat) {
    for (ir = 0; ir < nr-1; ir++) {
      Conv_Ptr->M_rx[ir][0] = 0.0;
      for (it = 0; it < nt; it++)
        Conv_Ptr->M_rx[ir][0] += Out_Ptr->Td_rat[ir][ia][it] * dt;
    }
       
  } else                        /* In_Ptr->record.Td_ra == 1 */
    for (ir = 0; ir < nr-1; ir++)
      Conv_Ptr->M_rx[ir][0] = Out_Ptr->Td_ra[ir][ia];
 
  if (mode == 0)
    fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Conv_Ptr->M_rx[ir][0]);
  else {
    Conv_Ptr->Mb_x = NULL;
    Conv_Ptr->nxc = 1;
    ConvM_rx(Conv_Ptr, Td_rxc);
    for (ir = 0; ir < nr-1; ir ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Td_rxc[0][ir]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Td_rxc, 0, Conv_Ptr->nrc-1, 0, 0);
  }
 
  FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, 0);
  fclose(file);
}
 
/**************************************************************************
 *      2 numbers each line: r, Td_r@t.
 *      mode = 0, extract Td_r@t; mode = 1, extract Td_r@tc.
 ****/
void
ExtractTd_r_t(InStru * In_Ptr, OutStru * Out_Ptr,
              ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pt;
  FILE       *file;
  char        fname[STRLEN];
  double     **Td_rxc;
 
  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Tdr@t");
  else strcpy(fname, "Tdr@tc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;
 
  pt = GetFixedTZRA((nt-1)*dt, 0);
  it = pt / dt;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode ==0) {
    fprintf(file, "# Td_r at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Td [1/(cm2 ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_rc at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Td [J/(cm2 ps)]");
  }
 
  if (mode == 1) ConvResolution(Conv_Ptr);
 
  if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, 0, 0)) == NULL)
      || ((mode == 1) &&
         (Td_rxc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, 0, 0)) == NULL)) {
      printf("No enough memory to extract Td_r@t.\n");
      return;
    }
 
  if (In_Ptr->record.Td_rat) {
    for (ir = 0; ir < nr-1; ir++) {
      Conv_Ptr->M_rx[ir][0] = 0.0;
      for (ia = 0; ia < na; ia++)
        Conv_Ptr->M_rx[ir][0] += Out_Ptr->Td_rat[ir][ia][it]
	      * PI * sin(2.0 * (ia + 0.5) * da) * da;
    }
       
  } else                        /* In_Ptr->record.Td_rt == 1 */
    for (ir = 0; ir < nr-1; ir++)
      Conv_Ptr->M_rx[ir][0] = Out_Ptr->Td_rt[ir][it];
 
  if (mode == 0)
    fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Conv_Ptr->M_rx[ir][0]);
  else {
    Conv_Ptr->Mb_x = NULL;
    Conv_Ptr->nxc = 1;
    ConvM_rx(Conv_Ptr, Td_rxc);
    for (ir = 0; ir < nr-1; ir ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Td_rxc[0][ir]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Td_rxc, 0, Conv_Ptr->nrc-1, 0, 0);
  }
 
  FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, 0);
  fclose(file);
}
 
/**************************************************************************
 *      2 numbers each line: r, Td_r.
 *      mode = 0, extract Td_r; mode = 1, extract Td_rc.
 ****/
void
ExtractTd_r(InStru * In_Ptr, OutStru * Out_Ptr,
              ConvStru *Conv_Ptr, char mode)
{
  short       ir;
  short       nr = In_Ptr->nr;
  double      dr = In_Ptr->dr;
  FILE       *file;
  char        fname[STRLEN];
  double     **Td_rxc;
 
  fname[0] = '\0';
  if (mode == 0) strcpy(fname, "Tdr");
#if IBMPC == 0
  if (mode == 1) strcpy(fname, "Tdrc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;
 
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_r\n");
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Td [1/cm2]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_rc\n");
    fprintf(file, "# %-10s\t%-12s\n", "r [cm]", "Td [J/cm2]");
  }
 
  if (!In_Ptr->record.Td_r)
    if ((Out_Ptr->Td_r = AllocArray1D(0, nr - 1, 0)) == NULL) {
      printf("Allocating Td_r failed.\n");
      return;
    }
  if (In_Ptr->record.Td_rat)
    Td_rFromTd_rat(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Td_ra)
    Td_rFromTd_ra(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Td_rt)
    Td_rFromTd_rt(In_Ptr, Out_Ptr);
 
  if (mode == 0)
    for (ir = 0; ir < nr-1; ir++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5) * dr, Out_Ptr->Td_r[ir]);
  else {
    if (mode == 1) ConvResolution(Conv_Ptr);
 
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, 0, 0)) == NULL)
      || ((Td_rxc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, 0, 0)) == NULL)) {
      printf("No enough memory to extract Td_r.\n");
      return;
    }
 
    for (ir = 0; ir < nr-1; ir ++)
      Conv_Ptr->M_rx[ir][0] = Out_Ptr->Td_r[ir];
 
    Conv_Ptr->Mb_x = NULL;
    Conv_Ptr->nxc = 1;
    ConvM_rx(Conv_Ptr, Td_rxc);
    for (ir = 0; ir < Conv_Ptr->nrc; ir ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ir+0.5)*dr, Td_rxc[ir][0]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Td_rxc, 0, Conv_Ptr->nrc-1, 0, 0);
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, 0);
  }
 
  if (!In_Ptr->record.Td_r)
    FreeArray1D(Out_Ptr->Td_r, 0, nr - 1);
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: a, Td_a@r@t.
 *      mode = 0, extract Td_a@r@t; mode = 1, extract Td_a@r@tc.
 ****/
void
ExtractTd_a_r_t(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pr, pt;
  FILE       *file;
  char        fname[STRLEN];
  double    **Td_rac;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Tda@r@t");
  else strcpy(fname, "Tda@r@tc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pr = GetFixedTZRA((nr-1)*dr, 2);
  pt = GetFixedTZRA((nt-1)*dt, 0);
  ir = pr / dr; it = pt / dt;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_a at r = %G cm and t = %G ps\n", pr, pt);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Td [1/(cm2 sr ps)]");
  } else {

    fprintf(file, "# Td_ac at r = %G cm and t = %G ps\n", pr, pt);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Td [J/(cm2 sr ps)]");
  }

  if (mode == 0)
    for (ia = 0; ia < na; ia++)
      fprintf(file, "%-12.4E\t%-12.4E\n", 
	     (ia+0.5)*da, Out_Ptr->Td_rat[ir][ia][it]);
  else {
    Conv_Ptr->drc = 2 * pr;
    Conv_Ptr->nrc = 1;
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, na-1, 0)) == NULL)
      || ((Td_rac = AllocArray2D(0, Conv_Ptr->nrc-1, 0, na-1, 0)) == NULL)) {
      printf("No enough memory to convolve Td_a@r@t.\n");
      return;
    }

    Conv_Ptr->Mb_x = NULL;
    for (ir = 0; ir < nr-1; ir ++)
      for (ia = 0; ia < na; ia ++) {
        Conv_Ptr->M_rx[ir][ia] = 0.0;
        for (it = 0; it < nt; it ++) 
          Conv_Ptr->M_rx[ir][ia] += Out_Ptr->Td_rat[ir][ia][it] * dt; 
      }
    
    Conv_Ptr->dxc = da;
    Conv_Ptr->nxc = na;
    ConvM_rx(Conv_Ptr, Td_rac);
    for(ia = 0; ia < na; ia ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ia+0.5)*da, Td_rac[0][ia]);

    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, na-1);
    FreeArray2D(Td_rac, 0, 0, 0, na-1);
  }
 
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: a, Td_a@r.
 *      mode = 0, extract Td_a@r; mode = 1, extract Td_a@rc.
 ****/
void
ExtractTd_a_r(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pr, temp;
  FILE       *file;
  char        fname[STRLEN];
  double      ** Td_rac;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Tda@r");
  else strcpy(fname, "Tda@rc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pr = GetFixedTZRA((nr-1)*dr, 2);
  ir = pr / dr;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_a at r = %G cm\n", pr);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Td [1/(cm2 sr)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_ac at r = %G cm\n", pr);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Td [J/(cm2 sr)]");
  }

  if (mode == 0) {
    if (In_Ptr->record.Td_rat) {
      for (ia = 0; ia < na; ia++) {
	temp = 0.0;
	for (it = 0; it < nt; it++)
	  temp += Out_Ptr->Td_rat[ir][ia][it] * dt;
	fprintf(file, "%-12.4E\t%-12.4E\n", (ia+0.5) * da, temp);
      }

    } else {			/* In_Ptr->record.Td_ra == 1 */
      for (ia = 0; ia < na; ia++)
	fprintf(file, "%-12.4E\t%-12.4E\n", 
		(ia+0.5) * da, Out_Ptr->Td_ra[ir][ia]);
    }
  } else {
    Conv_Ptr->drc = 2 * pr;
    Conv_Ptr->nrc = 1;
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, na-1, 0)) == NULL)
      || ((Td_rac = AllocArray2D(0, Conv_Ptr->nrc-1, 0, na-1, 0)) == NULL)) {
      printf("No enough memory to convolve Td_a@r@t.\n");
      return;
    }
   
    Conv_Ptr->Mb_x = NULL;
    if (In_Ptr->record.Td_rat) {
      for (ir = 0; ir < nr-1; ir ++)
	for (ia = 0; ia < na; ia++) {
	  Conv_Ptr->M_rx[ir][ia] = 0.0;
	  for (it = 0; it < nt; it++)
	    Conv_Ptr->M_rx[ir][ia] += Out_Ptr->Td_rat[ir][ia][it] * dt;
	}
    } else {			/* In_Ptr->record.Td_ra == 1 */
      for (ir = 0; ir < nr-1; ir++)
	for (ia = 0; ia < na; ia++)
	  Conv_Ptr->M_rx[ir][ia] += Out_Ptr->Td_ra[ir][ia];
    }

    Conv_Ptr->dxc = da;
    Conv_Ptr->nxc = na;
    ConvM_rx(Conv_Ptr, Td_rac);
    for(ia = 0; ia < na; ia ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (ia+0.5)*da, Td_rac[0][ia]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, na-1);
    FreeArray2D(Td_rac, 0, 0, 0, na-1);
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: a, Td_a@t.
 *      mode = 0, extract Td_a@t; mode = 1, extract Td_a@tc.
 ****/
void
ExtractTd_a_t(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pt, temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Tda@t");
  else strcpy(fname, "Tda@tc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pt = GetFixedTZRA((nt-1)*dt, 0);
  it = pt / dt;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_a at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Td [1/(sr ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_ac at t = %G ps\n", pt);
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Td [J/(sr ps)]");
  }

  if (In_Ptr->record.Td_rat) {
    for (ia = 0; ia < na; ia++) {
      temp = 0.0;
      for (ir = 0; ir < nr-1; ir++)
	temp += Out_Ptr->Td_rat[ir][ia][it]
	  * 2.0 * PI * (ir + 0.5) * dr * dr;
      if (mode == 1) temp *= Conv_Ptr->beam.P;
      fprintf(file, "%-12.4E\t%-12.4E\n", (ia+0.5) * da, temp);
    }

  } else {			/* In_Ptr->record.Td_at == 1 */
    for (ia = 0; ir < na; ia++) {
      temp = Out_Ptr->Td_at[ia][it];
      fprintf(file, "%-12.4E\t%-12.4E\n", (ia+0.5)*da, temp);
    }
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: a, Td_a.
 *      mode = 0, extract Td_a; mode = 1, extract Td_ac.
 ****/
void
ExtractTd_a(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ia;
  short       na = In_Ptr->na;
  double      da = In_Ptr->da;
  double      temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
  if (mode == 0) strcpy(fname, "Tda");
#if IBMPC == 0
  if (mode == 1) strcpy(fname, "Tdac");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_a\n");
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Td [1/sr]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_ac\n");
    fprintf(file, "# %-10s\t%-12s\n", "a [sr]", "Td [J/sr]");
  }

  if (!In_Ptr->record.Td_a)
    if ((Out_Ptr->Td_a = AllocArray1D(0, na - 1, 0)) == NULL) {
      printf("Allocating Td_a failed.\n");
      return;
    }
  if (In_Ptr->record.Td_rat)
    Td_aFromTd_rat(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Td_ra)
    Td_aFromTd_ra(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Td_at)
    Td_aFromTd_at(In_Ptr, Out_Ptr);

  for (ia = 0; ia < na; ia++) {
    temp = Out_Ptr->Td_a[ia];
    if (mode == 1) temp *= Conv_Ptr->beam.P;
    fprintf(file, "%-12.4E\t%-12.4E\n", (ia + 0.5) * da, temp);
  }

  if (!In_Ptr->record.Td_a)
    FreeArray1D(Out_Ptr->Td_a, 0, na - 1);
  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, Rd_t@r@a.
 *      mode = 0, extract Td_t@r@a; mode = 1, extract Td_t@r@ac.
 ****/
void
ExtractTd_t_r_a(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pr, pa;
  FILE       *file;
  char        fname[STRLEN];
  double    **Td_rtc;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Tdt@r@a");
  else strcpy(fname, "Tdt@r@ac");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pr = GetFixedTZRA((nr-1)*dr, 2);
  pa = GetFixedTZRA((na-1)*da, 3);
  ir = pr / dr; ia = pa / da;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_a at r = %G cm and a = %G sr\n", pr, pa);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Td [1/(cm2 sr ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_ac at r = %G cm and a = %G sr\n", pr, pa);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Td [J/(cm2 sr ps)]");
  }

  if (mode == 0)
    for (it = 0; it < nt; it++)
      fprintf(file, "%-12.4E\t%-12.4E\n", 
	      (it+0.5)*dt, Out_Ptr->Td_rat[ir][ia][it]);
  else {
    Conv_Ptr->drc = 2 * pr;
    Conv_Ptr->nrc = 1;
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, na-1, 0)) == NULL)
      || ((Td_rtc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, nt-1, 0)) == NULL)) {
      printf("No enough memory to convolve Td_a@r@t.\n");
      return;
    }  
       
    Conv_Ptr->Mb_x = NULL;
    for (ir = 0; ir < nr-1; ir ++)
      for (ia = 0; ia < na; ia ++) {
        Conv_Ptr->M_rx[ir][ia] = 0.0;
        for (it = 0; it < nt; it ++)
          Conv_Ptr->M_rx[ir][ia] += Out_Ptr->Td_rat[ir][ia][it] * dt;
      }
       
    Conv_Ptr->dxc = dt;
    Conv_Ptr->nxc = nt;
    ConvM_rx(Conv_Ptr, Td_rtc);
    for(it = 0; it < nt; it ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, Td_rtc[0][it]);
       
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, nt-1);
    FreeArray2D(Td_rtc, 0, 0, 0, nt-1);
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, Td_t@r.
 *      mode = 0, extract Td_t@r; mode = 1, extract Td_t@rc.
 ****/
void
ExtractTd_t_r(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pr, temp;
  FILE       *file;
  char        fname[STRLEN];
  double    **Td_rtc;

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Tdt@r");
  else strcpy(fname, "Tdt@rc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pr = GetFixedTZRA((nr-1)*dr, 2);
  ir = pr / dr;
  fprintf(file, "# %s extracted from %s\n\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_t at r = %G cm\n", pr);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Td [1/(cm2 ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_tc at r = %G cm\n", pr);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Td [J/(cm2 ps)]");
  }

  if (mode == 0) {
    if (In_Ptr->record.Td_rat) {
      for (it = 0; it < nt; it++) {
	temp = 0.0;
	for (ia = 0; ia < na; ia++)
	  temp += Out_Ptr->Td_rat[ir][ia][it]
	      * PI * sin(2.0 * (ia + 0.5) * da) * da;
	fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5) * dt, temp);
      }
    } else {			/* In_Ptr->record.Td_rt == 1 */
      for (it = 0; it < nt; it++)
	fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, Out_Ptr->Td_rt[ir][it]);
    }

  } else {
    Conv_Ptr->drc = 2 * pr;
    Conv_Ptr->nrc = 1;
    if (((Conv_Ptr->M_rx = AllocArray2D(0, nr-1, 0, nt-1, 0)) == NULL)
      || ((Td_rtc = AllocArray2D(0, Conv_Ptr->nrc-1, 0, nt-1, 0)) == NULL)) {
      printf("No enough memory to convolve Td_t@rh.\n");
      return;
    }
 
    Conv_Ptr->Mb_x = NULL;
    if (In_Ptr->record.Td_rat) {
      for (ir = 0; ir < nr-1; ir ++)
        for (it = 0; it < nt; it++) {
          Conv_Ptr->M_rx[ir][it] = 0.0;
	  for (ia = 0; ia < na; ia++)
            Conv_Ptr->M_rx[ir][ia] += Out_Ptr->Td_rat[ir][ia][it]
	      * PI * sin(2.0 * (ia + 0.5) * da) * da;
        }
    } else {                    /* In_Ptr->record.Td_rt == 1 */
      for (ir = 0; ir < nr-1; ir++)
        for (it = 0; it < nt; it++)
          Conv_Ptr->M_rx[ir][it] += Out_Ptr->Td_rt[ir][it];
    }
       
    Conv_Ptr->dxc = dt;
    Conv_Ptr->nxc = nt;
    ConvM_rx(Conv_Ptr, Td_rtc);
    for(it = 0; it < nt; it ++)
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, Td_rtc[0][it]);
 
    Conv_Ptr->drc = In_Ptr->dr;
    Conv_Ptr->nrc = In_Ptr->nr;
    FreeArray2D(Conv_Ptr->M_rx, 0, nr-1, 0, nt-1);
    FreeArray2D(Td_rtc, 0, 0, 0, nt-1);
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, Td_t@a.
 *      mode = 0, extract Td_t@a; mode = 1, extract Td_t@ac.
 ****/
void
ExtractTd_t_a(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       ir, ia, it;
  short       nr = In_Ptr->nr;
  short       na = In_Ptr->na;
  short       nt = In_Ptr->nt;
  double      dr = In_Ptr->dr;
  double      da = In_Ptr->da;
  double      dt = In_Ptr->dt;
  double      pa, temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
#if IBMPC == 0
  if (mode == 0) strcpy(fname, "Tdt@a");
  else strcpy(fname, "Tdt@ac");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  pa = GetFixedTZRA((na-1)*da, 3);
  ia = pa / da;
  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_t at a = %G sr\n", pa);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Td [1/(sr ps)]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_tc at a = %G sr\n", pa);
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Td [J/(sr ps)]");
  }

  if (In_Ptr->record.Td_rat) {
    for (it = 0; it < nt; it++) {
      temp = 0.0;
      for (ir = 0; ir < nr-1; ir++)
	temp += Out_Ptr->Td_rat[ir][ia][it]
	  * 2.0 * PI * (ir + 0.5) * dr * dr;
      if (mode == 1) temp *= Conv_Ptr->beam.P;
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5) * dt, temp);
    }

  } else {			/* In_Ptr->record.Td_at == 1 */
    for (it = 0; it < nt; it++) {
      temp = Out_Ptr->Td_at[ia][it];
      if (mode == 1) temp *= Conv_Ptr->beam.P;
      fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5)*dt, Out_Ptr->Td_at[ia][it]);
    }
  }

  fclose(file);
}

/**************************************************************************
 *      2 numbers each line: t, Td_t.
 *      mode = 0, extract Td_t; mode = 1, extract Td_tc.
 ****/
void
ExtractTd_t(InStru * In_Ptr, OutStru * Out_Ptr,
	      ConvStru *Conv_Ptr, char mode)
{
  short       it;
  short       nt = In_Ptr->nt;
  double      dt = In_Ptr->dt;
  double      temp;
  FILE       *file;
  char        fname[STRLEN];

  fname[0] = '\0';
  if (mode == 0) strcpy(fname, "Tdt");
#if IBMPC == 0
  if (mode == 1) strcpy(fname, "Tdtc");
#endif
  if ((file = GetWriteFile(fname)) == NULL)
    return;

  fprintf(file, "# %s extracted from %s\n", fname, In_Ptr->out_fname);
  if (mode == 0) {
    fprintf(file, "# Td_t\n");
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Td [1/ps]");
  } else {
    SpecifyBeam(file, Conv_Ptr);
    fprintf(file, "# Td_tc\n");
    fprintf(file, "# %-10s\t%-12s\n", "t [ps]", "Td [J/ps]");
  }

  if (!In_Ptr->record.Td_t)
    if ((Out_Ptr->Td_t = AllocArray1D(0, nt - 1, 0)) == NULL) {
      printf("Allocating Td_t failed.\n");
      return;
    }
  if (In_Ptr->record.Td_rat)
    Td_tFromTd_rat(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Td_rt)
    Td_tFromTd_rt(In_Ptr, Out_Ptr);
  else if (In_Ptr->record.Td_at)
    Td_tFromTd_at(In_Ptr, Out_Ptr);

  for (it = 0; it < nt; it++) {
    temp = Out_Ptr->Td_t[it];
    if (mode == 1)
      temp *= Conv_Ptr->beam.P;
    fprintf(file, "%-12.4E\t%-12.4E\n", (it+0.5) * dt, temp);
  }

  if (!In_Ptr->record.Td_t)
    FreeArray1D(Out_Ptr->Td_t, 0, nt - 1);
  fclose(file);
}

/**************************************************************************
 *	Print recorded quantities in MCML.
 ****/
void
PrintRecordQuan(InStru * In_Ptr)
{
  printf("Scored quantities in MCML:\n");

  if (In_Ptr->record.Rd_r) printf("Rd_r \t");
  if (In_Ptr->record.Rd_a) printf("Rd_a \t");
  if (In_Ptr->record.Rd_ra) printf("Rd_ra \t");
  if (In_Ptr->record.Rd_t) printf("Rd_t \t");
  if (In_Ptr->record.Rd_rt) printf("Rd_rt \t");
  if (In_Ptr->record.Rd_at) printf("Rd_at \t");
  if (In_Ptr->record.Rd_rat) printf("Rd_rat \t");

  if (In_Ptr->record.Td_r) printf("Td_r \t");
  if (In_Ptr->record.Td_a) printf("Td_a \t");
  if (In_Ptr->record.Td_ra) printf("Td_ra \t");
  if (In_Ptr->record.Td_t) printf("Td_t \t");
  if (In_Ptr->record.Td_rt) printf("Td_rt \t");
  if (In_Ptr->record.Td_at) printf("Td_at \t");
  if (In_Ptr->record.Td_rat) printf("Td_rat \t");

  if (In_Ptr->record.A_z) printf("A_z \t");
  if (In_Ptr->record.A_rz) printf("A_rz \t");
  if (In_Ptr->record.A_t) printf("A_t \t");
  if (In_Ptr->record.A_zt) printf("A_zt \t");
  if (In_Ptr->record.A_rzt) printf("A_rzt \t");

  printf("\n");
}

/**************************************************************************
 *	Print 6 items each line.
 ****/
PrintAndNewLine(short *Index, char *String)
{
  (*Index)++;
  if (*Index > 7) {
    printf("\n");
    *Index = 0;
  }
  printf("%-10s", String);
}

/**************************************************************************
 *	Print the extractable quantities.
 *      mode = 0, print the original quantities;
 *      mode = 1, print the convolved quantities.
 ****/
PrintExtractQuan(ConvStru * Conv_Ptr, char mode)
{
  short         index = 0;
  char          needspace = 0;

  printf("\nThe quantities available to be extracted: \n");
  if (Conv_Ptr->extract.A_t_r_z) {
    PrintAndNewLine(&index, "A_t@r@z"); needspace = 1; }
  if (Conv_Ptr->extract.A_rz_t) {
    PrintAndNewLine(&index, "A_rz@t"); needspace = 1; }
  if (Conv_Ptr->extract.A_z_t) {
    PrintAndNewLine(&index, "A_z@t"); needspace = 1; }
  if (Conv_Ptr->extract.A_t_z) {
    PrintAndNewLine(&index, "A_t@z"); needspace = 1; }
  if (Conv_Ptr->extract.A_rz) {
    PrintAndNewLine(&index, "A_rz"); needspace = 1; }
  if (Conv_Ptr->extract.A_z) {
    PrintAndNewLine(&index, "A_z"); needspace = 1; }
  if (Conv_Ptr->extract.A_t) {
    PrintAndNewLine(&index, "A_t"); needspace = 1; }
  if (Conv_Ptr->extract.A_l && mode == 0) {
    PrintAndNewLine(&index, "A_l"); needspace = 1; }

  index = 0;
  if (needspace) puts("\n");
  needspace = 0;
  if (Conv_Ptr->extract.A_t_r_z) {
    PrintAndNewLine(&index, "F_t@r@z"); needspace = 1; }
  if (Conv_Ptr->extract.A_rz_t) {
    PrintAndNewLine(&index, "F_rz@t"); needspace = 1; }
  if (Conv_Ptr->extract.A_z_t) {
    PrintAndNewLine(&index, "F_z@t"); needspace = 1; }
  if (Conv_Ptr->extract.A_t_z) {
    PrintAndNewLine(&index, "F_t@z"); needspace = 1; }
  if (Conv_Ptr->extract.A_rz) {
    PrintAndNewLine(&index, "F_rz"); needspace = 1; }
  if (Conv_Ptr->extract.A_z) {
    PrintAndNewLine(&index, "F_z"); needspace = 1; }

  index = 0;
  if (needspace) puts("\n");
  needspace = 0;
  if (Conv_Ptr->extract.Rd_r_a_t) {
    PrintAndNewLine(&index, "Rd_r@a@t"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_a_r_t) {
    PrintAndNewLine(&index, "Rd_a@r@t"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_t_r_a) {
    PrintAndNewLine(&index, "Rd_t@r@a"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_r_a) {
    PrintAndNewLine(&index, "Rd_r@a"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_r_t) {
    PrintAndNewLine(&index, "Rd_r@t"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_a_r) {
    PrintAndNewLine(&index, "Rd_a@r"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_a_t) {
    PrintAndNewLine(&index, "Rd_a@t"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_t_r) {
    PrintAndNewLine(&index, "Rd_t@r"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_t_a) {
    PrintAndNewLine(&index, "Rd_t@a"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_r) {
    PrintAndNewLine(&index, "Rd_r"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_a) {
    PrintAndNewLine(&index, "Rd_a"); needspace = 1; }
  if (Conv_Ptr->extract.Rd_t) {
    PrintAndNewLine(&index, "Rd_t"); needspace = 1; }

  index = 0;
  if (needspace) puts("\n");
  needspace = 0;
  if (Conv_Ptr->extract.Td_r_a_t) {
    PrintAndNewLine(&index, "Td_r@a@t"); needspace = 1; }
  if (Conv_Ptr->extract.Td_a_r_t) {
    PrintAndNewLine(&index, "Td_a@r@t"); needspace = 1; }
  if (Conv_Ptr->extract.Td_t_r_a) {
    PrintAndNewLine(&index, "Td_t@r@a"); needspace = 1; }
  if (Conv_Ptr->extract.Td_r_a) {
    PrintAndNewLine(&index, "Td_r@a"); needspace = 1; }
  if (Conv_Ptr->extract.Td_r_t) {
    PrintAndNewLine(&index, "Td_r@t"); needspace = 1; }
  if (Conv_Ptr->extract.Td_a_r) {
    PrintAndNewLine(&index, "Td_a@r"); needspace = 1; }
  if (Conv_Ptr->extract.Td_a_t) {
    PrintAndNewLine(&index, "Td_a@t"); needspace = 1; }
  if (Conv_Ptr->extract.Td_t_r) {
    PrintAndNewLine(&index, "Td_t@r"); needspace = 1; }
  if (Conv_Ptr->extract.Td_t_a) {
    PrintAndNewLine(&index, "Td_t@a"); needspace = 1; }
  if (Conv_Ptr->extract.Td_r) {
    PrintAndNewLine(&index, "Td_r"); needspace = 1; }
  if (Conv_Ptr->extract.Td_a) {
    PrintAndNewLine(&index, "Td_a"); needspace = 1; }
  if (Conv_Ptr->extract.Td_t) {
    PrintAndNewLine(&index, "Td_t"); needspace = 1; }

  if (needspace) puts("\n");
  puts("  A  -- absorption probability");
  puts("  Rd -- diffuse reflectance");
  puts("  Td -- diffuse transmittance");
  puts("  a  -- (alpha) light exit angle");
  puts("  r, z -- cylindrical coordinates\n");
}

/**************************************************************************
 *      Change all characters in a string to upper case.
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
 *	mode = 0, extract original data; 
 ****/
void
BranchOrigQuantity(char *Quantity,
	       InStru * In_Ptr,
	       OutStru * Out_Ptr,
	       ConvStru * Conv_Ptr)
{
  if ((strcmp(Quantity, "A_RZ@T") == 0) && Conv_Ptr->extract.A_rz_t)
    ExtractA_rz_t(In_Ptr, Out_Ptr, 0);
  else if ((strcmp(Quantity, "A_RZ") == 0) && Conv_Ptr->extract.A_rz)
    ExtractA_rz(In_Ptr, Out_Ptr, 0);
  else if ((strcmp(Quantity, "A_Z@T") == 0) && Conv_Ptr->extract.A_z_t)
    ExtractA_z_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "A_Z") == 0) && Conv_Ptr->extract.A_z)
    ExtractA_z(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "A_T@R@Z") == 0) && Conv_Ptr->extract.A_t_r_z)
    ExtractA_t_r_z(In_Ptr, Out_Ptr, 0);
  else if ((strcmp(Quantity, "A_T@Z") == 0) && Conv_Ptr->extract.A_t_z)
    ExtractA_t_z(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "A_T") == 0) && Conv_Ptr->extract.A_t)
    ExtractA_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "A_L") == 0) && Conv_Ptr->extract.A_l)
    ExtractA_l(In_Ptr, Out_Ptr);

  else if ((strcmp(Quantity, "F_RZ@T") == 0) && Conv_Ptr->extract.A_rz_t)
    ExtractA_rz_t(In_Ptr, Out_Ptr, 1);
  else if ((strcmp(Quantity, "F_RZ") == 0) && Conv_Ptr->extract.A_rz)
    ExtractA_rz(In_Ptr, Out_Ptr, 1);  
  else if ((strcmp(Quantity, "F_Z@T") == 0) && Conv_Ptr->extract.A_z_t)
    ExtractA_z_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "F_Z") == 0) && Conv_Ptr->extract.A_z)
    ExtractA_z(In_Ptr, Out_Ptr, Conv_Ptr, 1);  
  else if ((strcmp(Quantity, "F_T@R@Z") == 0) && Conv_Ptr->extract.A_t_r_z)
    ExtractA_t_r_z(In_Ptr, Out_Ptr, 1);
  else if ((strcmp(Quantity, "F_T@z") == 0) && Conv_Ptr->extract.A_t_z)
    ExtractA_t_z(In_Ptr, Out_Ptr, Conv_Ptr, 1);  

  else if ((strcmp(Quantity, "RD_R@A@T") == 0) && Conv_Ptr->extract.Rd_r_a_t)
    ExtractRd_r_a_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_R@A") == 0) && Conv_Ptr->extract.Rd_r_a)
    ExtractRd_r_a(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_R@T") == 0) && Conv_Ptr->extract.Rd_r_t)
    ExtractRd_r_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_R") == 0) && Conv_Ptr->extract.Rd_r)
    ExtractRd_r(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_A@R@T") == 0) && Conv_Ptr->extract.Rd_a_r_t)
    ExtractRd_a_r_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_A@R") == 0) && Conv_Ptr->extract.Rd_a_r)
    ExtractRd_a_r(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_A@T") == 0) && Conv_Ptr->extract.Rd_a_t)
    ExtractRd_a_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_A") == 0) && Conv_Ptr->extract.Rd_a)
    ExtractRd_a(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_T@R@A") == 0) && Conv_Ptr->extract.Rd_t_r_a)
    ExtractRd_t_r_a(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_T@R") == 0) && Conv_Ptr->extract.Rd_t_r)
    ExtractRd_t_r(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_T@A") == 0) && Conv_Ptr->extract.Rd_t_a)
    ExtractRd_t_a(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "RD_T") == 0) && Conv_Ptr->extract.Rd_t)
    ExtractRd_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);

  else if ((strcmp(Quantity, "TD_R@A@T") == 0) && Conv_Ptr->extract.Td_r_a_t)
    ExtractTd_r_a_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_R@A") == 0) && Conv_Ptr->extract.Td_r_a)
    ExtractTd_r_a(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_R@T") == 0) && Conv_Ptr->extract.Td_r_t)
    ExtractTd_r_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_R") == 0) && Conv_Ptr->extract.Td_r)
    ExtractTd_r(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_A@R@T") == 0) && Conv_Ptr->extract.Td_a_r_t)
    ExtractTd_a_r_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_A@R") == 0) && Conv_Ptr->extract.Td_a_r)
    ExtractTd_a_r(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_A@T") == 0) && Conv_Ptr->extract.Td_a_t)
    ExtractTd_a_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_A") == 0) && Conv_Ptr->extract.Td_a)
    ExtractTd_a(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_T@R@A") == 0) && Conv_Ptr->extract.Td_t_r_a)
    ExtractTd_t_r_a(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_T@R") == 0) && Conv_Ptr->extract.Td_t_r)
    ExtractTd_t_r(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_T@A") == 0) && Conv_Ptr->extract.Td_t_a)
    ExtractTd_t_a(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "TD_T") == 0) && Conv_Ptr->extract.Td_t)
    ExtractTd_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if (strcmp(Quantity, ".") == 0)
    return;
  else
    puts("Invalid quantity.");
}

/**************************************************************************
 ****/
void
ExtractOrigData(InStru * In_Ptr,
		OutStru * Out_Ptr,
		ConvStru * Conv_Ptr)
{
  char        cmd_str[STRLEN];

  if (!Conv_Ptr->datain)
    puts("...No data to extract");
  else {
    PrintRecordQuan(In_Ptr);
    PrintExtractQuan(Conv_Ptr, 1);

    printf("Specify quantity to be extracted (or . to quit): ");
    do
      gets(cmd_str);
    while (!strlen(cmd_str));     /* avoid null string. */
    BranchOrigQuantity(ToUpperString(cmd_str), In_Ptr, Out_Ptr, Conv_Ptr);
  }
}

/**************************************************************************
 ****/
BranchConvQuantity(char *Quantity,
               InStru * In_Ptr,
               OutStru * Out_Ptr,
               ConvStru * Conv_Ptr)
{
  if ((strcmp(Quantity, "A_RZ@T") == 0) && Conv_Ptr->extract.A_rz_t)
    ExtractConvA_rz_t(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "A_RZ") == 0) && Conv_Ptr->extract.A_rz)
    ExtractConvA_rz(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "A_Z@T") == 0) && Conv_Ptr->extract.A_z_t)
    ExtractA_z_t(In_Ptr, Out_Ptr, Conv_Ptr, 2);
  else if ((strcmp(Quantity, "A_Z") == 0) && Conv_Ptr->extract.A_z)
    ExtractA_z(In_Ptr, Out_Ptr, Conv_Ptr, 2);
  else if ((strcmp(Quantity, "A_T@R@Z") == 0) && Conv_Ptr->extract.A_t_r_z)
    ExtractConvA_t_r_z(In_Ptr, Out_Ptr, Conv_Ptr, 0);
  else if ((strcmp(Quantity, "A_T@Z") == 0) && Conv_Ptr->extract.A_t_z)
    ExtractA_t_z(In_Ptr, Out_Ptr, Conv_Ptr, 2);
  else if ((strcmp(Quantity, "A_T") == 0) && Conv_Ptr->extract.A_t)
    ExtractA_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
 
  else if ((strcmp(Quantity, "F_RZ@T") == 0) && Conv_Ptr->extract.A_rz_t)
    ExtractConvA_rz_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "F_RZ") == 0) && Conv_Ptr->extract.A_rz)
    ExtractConvA_rz(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "F_Z@T") == 0) && Conv_Ptr->extract.A_z_t)
    ExtractA_z_t(In_Ptr, Out_Ptr, Conv_Ptr,3);
  else if ((strcmp(Quantity, "F_Z") == 0) && Conv_Ptr->extract.A_z)
    ExtractA_z(In_Ptr, Out_Ptr, Conv_Ptr, 3);
  else if ((strcmp(Quantity, "F_T@R@Z") == 0) && Conv_Ptr->extract.A_t_r_z)
    ExtractConvA_t_r_z(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "F_T@Z") == 0) && Conv_Ptr->extract.A_t_z)
    ExtractA_t_z(In_Ptr, Out_Ptr, Conv_Ptr, 3);

  else if ((strcmp(Quantity, "RD_R@A@T") == 0) && Conv_Ptr->extract.Rd_r_a_t)
    ExtractRd_r_a_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_R@A") == 0) && Conv_Ptr->extract.Rd_r_a)
    ExtractRd_r_a(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_R@T") == 0) && Conv_Ptr->extract.Rd_r_t)
    ExtractRd_r_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_R") == 0) && Conv_Ptr->extract.Rd_r)
    ExtractRd_r(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_A@R@T") == 0) && Conv_Ptr->extract.Rd_a_r_t)
    ExtractRd_a_r_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_A@R") == 0) && Conv_Ptr->extract.Rd_a_r)
    ExtractRd_a_r(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_A@T") == 0) && Conv_Ptr->extract.Rd_a_t)
    ExtractRd_a_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_A") == 0) && Conv_Ptr->extract.Rd_a)
    ExtractRd_a(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_T@R@A") == 0) && Conv_Ptr->extract.Rd_t_r_a)
    ExtractRd_t_r_a(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_T@R") == 0) && Conv_Ptr->extract.Rd_t_r)
    ExtractRd_t_r(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_T@A") == 0) && Conv_Ptr->extract.Rd_t_a)
    ExtractRd_t_a(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "RD_T") == 0) && Conv_Ptr->extract.Rd_t)
    ExtractRd_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);

  else if ((strcmp(Quantity, "TD_R@A@T") == 0) && Conv_Ptr->extract.Td_r_a_t)
    ExtractTd_r_a_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_R@A") == 0) && Conv_Ptr->extract.Td_r_a)
    ExtractTd_r_a(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_R@T") == 0) && Conv_Ptr->extract.Td_r_t)
    ExtractTd_r_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_R") == 0) && Conv_Ptr->extract.Td_r)
    ExtractTd_r(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_A@R@T") == 0) && Conv_Ptr->extract.Td_a_r_t)
    ExtractTd_a_r_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_A@R") == 0) && Conv_Ptr->extract.Td_a_r)
    ExtractTd_a_r(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_A@T") == 0) && Conv_Ptr->extract.Td_a_t)
    ExtractTd_a_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_A") == 0) && Conv_Ptr->extract.Td_a)
    ExtractTd_a(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_T@R@A") == 0) && Conv_Ptr->extract.Td_t_r_a)
    ExtractTd_t_r_a(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_T@R") == 0) && Conv_Ptr->extract.Td_t_r)
    ExtractTd_t_r(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_T@A") == 0) && Conv_Ptr->extract.Td_t_a)
    ExtractTd_t_a(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if ((strcmp(Quantity, "TD_T") == 0) && Conv_Ptr->extract.Td_t)
    ExtractTd_t(In_Ptr, Out_Ptr, Conv_Ptr, 1);
  else if (strcmp(Quantity, ".") == 0)
    return;
  else
    puts("Invalid quantity.");
}

/**************************************************************************
 ****/
void
ExtractConvData(InStru * In_Ptr,
                OutStru * Out_Ptr,
                ConvStru * Conv_Ptr)
{
  char        cmd_str[STRLEN];

  if (!Conv_Ptr->datain)
    puts("...No data to extract");
  else if (Conv_Ptr->beam.type == original)
    puts("...No incident beam specified.");
  else {
    PrintRecordQuan(In_Ptr);
    PrintExtractQuan(Conv_Ptr, 1);

    printf("Specify quantity to be extracted (or . to quit): ");
    do
      gets(cmd_str);
    while (!strlen(cmd_str));     /* avoid null string. */
    BranchConvQuantity(ToUpperString(cmd_str), In_Ptr, Out_Ptr, Conv_Ptr);
  }
}

