/**************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *	Program name: MCML.
 *	Launch, move, and record photon weight.
 ****/

#include "mcml.h"

#define PARTIALREFLECTION 0
/* 1=split photon, 0=statistical reflection. */

#define ONE_MINUS_COSZERO 1.0E-12
/* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */

#define COS90D  1.0E-6
/* If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad. */

#define RandomNum (double)RandomGen(1, 0, NULL)

/**************************************************************************
 *	A random number generator that generates uniformly
 *	distributed random numbers between 0 and 1 inclusive.
 *	The algorithm is based on:
 *	W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *	Flannery, "Numerical Recipes in C," Cambridge University
 *	Press, 2nd edition, (1992).
 *	and
 *	D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *	of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *	When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *	When Type is 1, returns a random number.
 *	When Type is 2, gets the status of the generator.
 *	When Type is 3, restores the status of the generator.
 *
 *	The status of the generator is represented by Status[0..56].
 *
 *	Make sure you initialize the seed before you get random
 *	numbers.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double
RandomGen(char Type, long Seed, long *Status)
{
  static long i1, i2, ma[56];	/* ma[0] is not used. */
  long        mj, mk;
  short       i, ii;

  if (Type == 0) {		/* set seed. */
    mj = MSEED - (Seed < 0 ? -Seed : Seed);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (i = 1; i <= 54; i++) {
      ii = (21 * i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ)
	mk += MBIG;
      mj = ma[ii];
    }
    for (ii = 1; ii <= 4; ii++)
      for (i = 1; i <= 55; i++) {
	ma[i] -= ma[1 + (i + 30) % 55];
	if (ma[i] < MZ)
	  ma[i] += MBIG;
      }
    i1 = 0;
    i2 = 31;
  } else if (Type == 1) {	/* get a number. */
    if (++i1 == 56)
      i1 = 1;
    if (++i2 == 56)
      i2 = 1;
    mj = ma[i1] - ma[i2];
    if (mj < MZ)
      mj += MBIG;
    ma[i1] = mj;
    return (mj * FAC);
  } else if (Type == 2) {	/* get status. */
    for (i = 0; i < 55; i++)
      Status[i] = ma[i + 1];
    Status[55] = i1;
    Status[56] = i2;
  } else if (Type == 3) {	/* restore status. */
    for (i = 0; i < 55; i++)
      ma[i + 1] = Status[i];
    i1 = Status[55];
    i2 = Status[56];
  } else
    puts("Wrong parameter to RandomGen().");
  return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/**************************************************************************
 *	Compute the specular reflectance.
 *
 *	If the first layer is a turbid medium, use the Fresnel
 *	reflection from the boundary between the top abmient medium and
 *	the first layer as the specular reflectance.
 *
 *	The subroutine assumes the Layerspecs array is correctly
 *	initialized.
 ****/
double
Rspecular(LayerStru * Layerspecs_Ptr)
{
  double      r;

  r = (Layerspecs_Ptr[0].n - Layerspecs_Ptr[1].n)
    / (Layerspecs_Ptr[0].n + Layerspecs_Ptr[1].n);
  r = r * r;

  return (r);
}

/**************************************************************************
 *      Choose a new direction for photon propagation by
 *      sampling
 *	1. the polar deflection angle theta
 *	2. the azimuthal angle psi.
 *
 *      Note:
 *      theta: 0 - pi so sin(theta) is always positive
 *      feel free to use sqrt() for cos(theta).
 *
 *      psi:   0 - 2pi
 *      for 0-pi  sin(psi) is +
 *      for pi-2pi sin(psi) is -
 ****/
void
Spin(double g,
     PhotonStru * Photon_Ptr)
{
  double      cost, sint;	/* cosine and sine of theta. */
  double      cosp, sinp;	/* cosine and sine of psi. */
  double      ux = Photon_Ptr->ux;
  double      uy = Photon_Ptr->uy;
  double      uz = Photon_Ptr->uz;
  double      psi;

  /* sample theta. */
  if (g == 0.0)
    cost = 2 * RandomNum - 1;
  else {
    double      temp = (1 - g * g) / (1 - g + 2 * g * RandomNum);
    cost = (1 + g * g - temp * temp) / (2 * g);
  }
  sint = sqrt(1.0 - cost * cost);	/* sqrt() is faster than sin(). */

  /* sample psi. */
  psi = 2.0 * PI * RandomNum;
  cosp = cos(psi);
  if (psi < PI)
    sinp = sqrt(1.0 - cosp * cosp);	/* sqrt() is faster than sin(). */
  else
    sinp = -sqrt(1.0 - cosp * cosp);

  /* update directional cosines. */
  if (1 - fabs(uz) <= ONE_MINUS_COSZERO) {	/* close to perpendicular. */
    Photon_Ptr->ux = sint * cosp;
    Photon_Ptr->uy = sint * sinp;
    Photon_Ptr->uz = cost * SIGN(uz);	/* SIGN() is faster than division. */
  } else {
    double      temp = sqrt(1.0 - uz * uz);
    Photon_Ptr->ux = sint * (ux * uz * cosp - uy * sinp) / temp + ux * cost;
    Photon_Ptr->uy = sint * (uy * uz * cosp + ux * sinp) / temp + uy * cost;
    Photon_Ptr->uz = -sint * cosp * temp + uz * cost;
  }
}

/**************************************************************************
 *	Initialize a photon packet.
 *	
 *	If an isotropic source is launched inside a glass layer, we check
 *	whether the photon will be total-internally reflected.  If it 
 *	does, the photon is killed to avoid a infinite travelling inside
 *	the glass layer.
 ****/
void
LaunchPhoton(double Rsp,
	     InStru * In_Ptr,
	     OutStru * Out_Ptr,
	     PhotonStru * Photon_Ptr)
{
  Photon_Ptr->w = 1.0 - Rsp;
  Photon_Ptr->alive = 1;
  Photon_Ptr->layer = (In_Ptr->slayer != 0) ? In_Ptr->slayer : 1;
  Photon_Ptr->s = 0;
  Photon_Ptr->sleft = 0;
  Photon_Ptr->scatters = 0;
  Photon_Ptr->time = 0;

  Photon_Ptr->x = 0.0;
  Photon_Ptr->y = 0.0;
  Photon_Ptr->z = In_Ptr->sz;
  Photon_Ptr->ux = 0.0;
  Photon_Ptr->uy = 0.0;
  Photon_Ptr->uz = 1.0;

  Out_Ptr->Ai = 0.0;
  Out_Ptr->Tbi = 0.0;
  Out_Ptr->Tdi = 0.0;
  Out_Ptr->Rbi = 0.0;
  Out_Ptr->Rdi = 0.0;

  if (In_Ptr->source_type == isotropic) {
    LayerStru   lstru;

    lstru = In_Ptr->layerspecs[Photon_Ptr->layer];
    Photon_Ptr->scatters++;	/* to avoid scoring into Rb or Tb. */
    Spin(0.0, Photon_Ptr);	/* isotropically scatter the photon. */

    if (lstru.mua == 0.0 && lstru.mus == 0.0)	/* glass layer. */
      if (fabs(Photon_Ptr->uz) <= lstru.cos_crit0 &&	/* total internal */
	  fabs(Photon_Ptr->uz) <= lstru.cos_crit1)	/* reflection. */
	Photon_Ptr->alive = 0;
  }
}

/**************************************************************************
 *	Move the photon S away in the current layer of medium.
 ****/
void
Hop(PhotonStru * Photon_Ptr,
    double S,
    double n)
{
  Photon_Ptr->x += S * Photon_Ptr->ux;
  Photon_Ptr->y += S * Photon_Ptr->uy;
  Photon_Ptr->z += S * Photon_Ptr->uz;
  Photon_Ptr->time += S * n * ONEOVERC;
}

/**************************************************************************
 *	Pick a step size in dimensionless unit for a photon packet.
 *	If the member s is zero, make a new step size
 *	with: -log(rnd).
 *	Otherwise, finish the leftover in s.
 ****/
void
SetStepSize(PhotonStru * Photon_Ptr)
{
  if (Photon_Ptr->s == 0.0) {	/* make a new step. */
    double      rnd;
    while ((rnd = RandomNum) <= 0.0);	/* avoid zero. */
    Photon_Ptr->s = -log(rnd);
  }
}

/**************************************************************************
 *	Return the distance between the photon position to the
 *	boundary along the photon direction.
 ****/
double
PathToBoundary(PhotonStru * Photon_Ptr,
	       InStru * In_Ptr)
{
  double      path;		/* length to boundary. */
  short       layer = Photon_Ptr->layer;
  double      uz = Photon_Ptr->uz;

  /* Distance to the boundary. */
  if (uz > 0.0)
    path = (In_Ptr->layerspecs[layer].z1
	    - Photon_Ptr->z) / uz;	/* path>0. */
  else if (uz < 0.0)
    path = (In_Ptr->layerspecs[layer].z0
	    - Photon_Ptr->z) / uz;	/* path>0. */
  else
    path = DBL_MAX;		/* infinity. */
  return (path);
}

/**************************************************************************
 *	Drop photon weight inside the tissue (not glass).
 *
 *      The photon is assumed alive.
 *
 *	The weight drop is dw = w*mua/(mua+mus).
 *
 *	The dropped weight is assigned to the absorption array
 *	elements.
 ****/
void
Drop(InStru * In_Ptr,
     PhotonStru * Photon_Ptr,
     OutStru * Out_Ptr)
{
  double      dwa;		/* absorbed weight. */
  double      x = Photon_Ptr->x;
  double      y = Photon_Ptr->y;
  short       iz, ir, it;	/* index to z, r & t. */
  short       layer = Photon_Ptr->layer;
  double      mua, mus, temp;
  RecordStru  record;

  record = In_Ptr->record;
  /* update photon weight. */
  mua = In_Ptr->layerspecs[layer].mua;
  mus = In_Ptr->layerspecs[layer].mus;
  dwa = Photon_Ptr->w * mua / (mua + mus);
  Photon_Ptr->w -= dwa;

  /* compute array indices. */
  if (record.A_rzt || record.A_zt || record.A_z || record.A_rz) {
    if (Photon_Ptr->z >= In_Ptr->zm)
      iz = In_Ptr->nz - 1;
    else
      iz = (short) (Photon_Ptr->z / In_Ptr->dz);
  }
  if (record.A_rzt || record.A_zt || record.A_t) {
    if (Photon_Ptr->time >= In_Ptr->tm)
      it = In_Ptr->nt - 1;
    else
      it = floor(Photon_Ptr->time / In_Ptr->dt);
  }
  /* assign dwa to an absorption array element. */
  if (Photon_Ptr->scatters) {	/* scattered. */
    if (record.A_rzt || record.A_rz) {
      if ((temp = sqrt(x * x + y * y)) >= In_Ptr->rm)
	ir = In_Ptr->nr - 1;
      else
	ir = (short) (temp / In_Ptr->dr);
    }
    if (record.A_rzt)
      Out_Ptr->A_rzt[ir][iz][it] += dwa;
    if (record.A_rz)
      Out_Ptr->A_rz[ir][iz] += dwa;
  } else {			/* ballistic. */
    if (record.A_rzt)
      Out_Ptr->Ab_zt[iz][it] += dwa;
    if (record.A_rz)
      Out_Ptr->Ab_z[iz] += dwa;
  }

  if (record.A_zt)
    Out_Ptr->A_zt[iz][it] += dwa;
  if (record.A_z)
    Out_Ptr->A_z[iz] += dwa;

  if (record.A_t)
    Out_Ptr->A_t[it] += dwa;
  Out_Ptr->Ai += dwa;
}

/**************************************************************************
 *	The photon weight is small, and the photon packet tries
 *	to survive a roulette.
 ****/
void
Roulette(PhotonStru * Photon_Ptr)
{
  if (Photon_Ptr->w == 0.0)	/* already dead. */
    Photon_Ptr->alive = 0;
  else if (RandomNum < CHANCE)	/* survived the roulette. */
    Photon_Ptr->w /= CHANCE;
  else
    Photon_Ptr->alive = 0;
}

/**************************************************************************
 *	Compute the Fresnel reflectance.
 *
 *	Make sure that the cosine of the incident angle ai
 *	is positive, and the case when the angle is greater
 *	than the critical angle is ruled out.
 *
 *	cai: cosine of the incident angle ai.
 *	cat_Ptr: pointer to the cosine of the transmission
 *			 angle at.
 *
 * 	Avoid trigonometric function operations as much as
 *	possible, because they are computation-intensive.
 ****/
double
RFresnel(double ni,		/* incident refractive index. */
	 double nt,		/* transmit refractive index. */
	 double cai,		/* cosine of angle ai. +. */
	 double *cat_Ptr)
{
  double      r;

  if (ni == nt) {		/* matched boundary. */
    *cat_Ptr = cai;
    r = 0.0;
  } else if (1 - cai <= ONE_MINUS_COSZERO) {	/* normal incidence. */
    *cat_Ptr = cai;
    r = (nt - ni) / (nt + ni);
    r *= r;
  } else if (cai < COS90D) {	/* very slant incidence. */
    *cat_Ptr = 0.0;
    r = 1.0;
  } else {			/* general incidence. */
    double      sai, sat;	/* sine of the angles ai & at. */

    sai = sqrt(1 - cai * cai);
    sat = ni * sai / nt;
    if (sat >= 1.0) {		/* total internal reflection. */
      *cat_Ptr = 0.0;
      r = 1.0;
    } else {
      double      cat;		/* cosine of at. */
      double      cap, cam;	/* cosines of ai+at & ai-at. */
      double      sap, sam;	/* sines of ai+at & ai-at. */

      *cat_Ptr = cat = sqrt(1 - sat * sat);
      cap = cai * cat - sai * sat;	/* c+ = cc - ss. */
      cam = cai * cat + sai * sat;	/* c- = cc + ss. */
      sap = sai * cat + cai * sat;	/* s+ = sc + cs. */
      sam = sai * cat - cai * sat;	/* s- = sc - cs. */
      r = 0.5 * sam * sam * (cam * cam + cap * cap) / (sap * sap * cam * cam);
      /* arranged for speed. */
    }
  }
  return (r);
}

/**************************************************************************
 *	Record the photon weight exiting the first layer(uz<0),
 *	to the reflection array.
 *
 *	Update the photon weight as well.
 ****/
void
RecordR(double Refl,		/* reflectance. */
	InStru * In_Ptr,
	PhotonStru * Photon_Ptr,
	OutStru * Out_Ptr)
{
  double      x = Photon_Ptr->x;
  double      y = Photon_Ptr->y;
  short       ir, ia, it;	/* index to r & angle. */
  double      temp;
  RecordStru  record;

  record = In_Ptr->record;
  if (record.Rd_rat || record.Rd_at || record.Rd_rt || record.Rd_t) {
    if (Photon_Ptr->time >= In_Ptr->tm)
      it = In_Ptr->nt - 1;
    else
      it = floor(Photon_Ptr->time / In_Ptr->dt);
  }
  if (Photon_Ptr->scatters) {	/* scattered. */
    if (record.Rd_rat || record.Rd_rt || record.Rd_ra || record.Rd_r) {
      if ((temp = sqrt(x * x + y * y)) >= In_Ptr->rm)
	ir = In_Ptr->nr - 1;
      else
	ir = (short) (temp / In_Ptr->dr);
    }
    if (record.Rd_rat || record.Rd_at || record.Rd_ra || record.Rd_a)
      if ((ia = acos(-Photon_Ptr->uz) / In_Ptr->da) > In_Ptr->na - 1)
	ia = In_Ptr->na - 1;

    /* assign photon weight to the reflection array element. */
    if (record.Rd_rat)
      Out_Ptr->Rd_rat[ir][ia][it] += Photon_Ptr->w * (1.0 - Refl);
    if (record.Rd_ra)
      Out_Ptr->Rd_ra[ir][ia] += Photon_Ptr->w * (1.0 - Refl);

    if (record.Rd_rt)
      Out_Ptr->Rd_rt[ir][it] += Photon_Ptr->w * (1.0 - Refl);
    if (record.Rd_r)
      Out_Ptr->Rd_r[ir] += Photon_Ptr->w * (1.0 - Refl);

    if (record.Rd_at)
      Out_Ptr->Rd_at[ia][it] += Photon_Ptr->w * (1.0 - Refl);
    if (record.Rd_a)
      Out_Ptr->Rd_a[ia] += Photon_Ptr->w * (1.0 - Refl);

    if (record.Rd_t)
      Out_Ptr->Rd_t[it] += Photon_Ptr->w * (1.0 - Refl);
    Out_Ptr->Rdi += Photon_Ptr->w * (1.0 - Refl);

  } else			/* ballistic. */
    Out_Ptr->Rbi += Photon_Ptr->w * (1.0 - Refl);

  Photon_Ptr->w *= Refl;
}

/**************************************************************************
 *	Record the photon weight exiting the last layer(uz>0),
 *	no matter whether the layer is glass or not, to the
 *	transmittance array.
 *
 *	Update the photon weight as well.
 ****/
void
RecordT(double Refl,
	InStru * In_Ptr,
	PhotonStru * Photon_Ptr,
	OutStru * Out_Ptr)
{
  double      x = Photon_Ptr->x;
  double      y = Photon_Ptr->y;
  short       ir, ia, it;	/* index to r & angle. */
  double      temp;
  RecordStru  record;

  record = In_Ptr->record;
  if (record.Td_rat || record.Td_at || record.Td_rt || record.Td_t) {
    if (Photon_Ptr->time >= In_Ptr->tm)
      it = In_Ptr->nt - 1;
    else
      it = floor(Photon_Ptr->time / In_Ptr->dt);
  }
  if (Photon_Ptr->scatters) {	/* scattered. */
    if (record.Td_rat || record.Td_rt || record.Td_ra || record.Td_r) {
      if ((temp = sqrt(x * x + y * y)) >= In_Ptr->rm)
	ir = In_Ptr->nr - 1;
      else
	ir = (short) (temp / In_Ptr->dr);
    }
    if (record.Td_rat || record.Td_at || record.Td_ra || record.Td_a)
      if ((ia = acos(Photon_Ptr->uz) / In_Ptr->da) > In_Ptr->na - 1)
	ia = In_Ptr->na - 1;

    /* assign photon weight to the transmittance array element. */
    if (record.Td_rat)
      Out_Ptr->Td_rat[ir][ia][it] += Photon_Ptr->w * (1.0 - Refl);
    if (record.Td_ra)
      Out_Ptr->Td_ra[ir][ia] += Photon_Ptr->w * (1.0 - Refl);

    if (record.Td_rt)
      Out_Ptr->Td_rt[ir][it] += Photon_Ptr->w * (1.0 - Refl);
    if (record.Td_r)
      Out_Ptr->Td_r[ir] += Photon_Ptr->w * (1.0 - Refl);

    if (record.Td_at)
      Out_Ptr->Td_at[ia][it] += Photon_Ptr->w * (1.0 - Refl);
    if (record.Td_a)
      Out_Ptr->Td_a[ia] += Photon_Ptr->w * (1.0 - Refl);

    if (record.Td_t)
      Out_Ptr->Td_t[it] += Photon_Ptr->w * (1.0 - Refl);
    Out_Ptr->Tdi += Photon_Ptr->w * (1.0 - Refl);
  } else			/* Collimated. */
    Out_Ptr->Tbi += Photon_Ptr->w * (1.0 - Refl);

  Photon_Ptr->w *= Refl;
}

/**************************************************************************
 *	Decide whether the photon will be transmitted or
 *	reflected on the upper boundary (uz<0) of the current
 *	layer.
 *
 *	If "layer" is the first layer, the photon packet will
 *	be partially transmitted and partially reflected if
 *	PARTIALREFLECTION is set to 1,
 *	or the photon packet will be either transmitted or
 *	reflected determined statistically if PARTIALREFLECTION
 *	is set to 0.
 *
 *	Record the transmitted photon weight as reflection.
 *
 *	If the "layer" is not the first layer and the photon
 *	packet is transmitted, move the photon to "layer-1".
 *
 *	Update the photon parmameters.
 ****/
void
CrossUpOrNot(InStru * In_Ptr,
	     PhotonStru * Photon_Ptr,
	     OutStru * Out_Ptr)
{
  double      uz = Photon_Ptr->uz;	/* z directional cosine. */
  double      uz1;		/* cosines of transmission alpha. always
				 * +. */
  double      r;		/* reflectance */
  short       layer = Photon_Ptr->layer;
  double      ni = In_Ptr->layerspecs[layer].n;
  double      nt = In_Ptr->layerspecs[layer - 1].n;

  /* Get r. */
  if (-uz <= In_Ptr->layerspecs[layer].cos_crit0)
    r = 1.0;			/* total internal reflection. */
  else
    r = RFresnel(ni, nt, -uz, &uz1);

#if PARTIALREFLECTION
  if (layer == 1 && r < 1.0) {	/* partially transmitted. */
    Photon_Ptr->uz = -uz1;	/* escaped photon. */
    RecordR(r, In_Ptr, Photon_Ptr, Out_Ptr);
    Photon_Ptr->uz = -uz;	/* reflected photon. */
  } else if (RandomNum > r) {	/* transmitted to layer-1. */
    Photon_Ptr->layer--;
    Photon_Ptr->ux *= ni / nt;
    Photon_Ptr->uy *= ni / nt;
    Photon_Ptr->uz = -uz1;
  } else			/* reflected. */
    Photon_Ptr->uz = -uz;
#else
  if (RandomNum > r) {		/* transmitted to layer-1. */
    if (layer == 1) {
      Photon_Ptr->uz = -uz1;
      RecordR(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
      Photon_Ptr->alive = 0;	/* escaped. */
    } else {
      Photon_Ptr->layer--;
      Photon_Ptr->ux *= ni / nt;
      Photon_Ptr->uy *= ni / nt;
      Photon_Ptr->uz = -uz1;
    }
  } else			/* reflected. */
    Photon_Ptr->uz = -uz;
#endif
}

/**************************************************************************
 *	Decide whether the photon will be transmitted  or be
 *	reflected on the bottom boundary (uz>0) of the current
 *	layer.
 *
 *	If the photon is transmitted, move the photon to
 *	"layer+1". If "layer" is the last layer, record the
 *	transmitted weight as transmittance. See comments for
 *	CrossUpOrNot.
 *
 *	Update the photon parmameters.
 ****/
void
CrossDnOrNot(InStru * In_Ptr,
	     PhotonStru * Photon_Ptr,
	     OutStru * Out_Ptr)
{
  double      uz = Photon_Ptr->uz;	/* z directional cosine. */
  double      uz1;		/* cosines of transmission alpha. */
  double      r;		/* reflectance */
  short       layer = Photon_Ptr->layer;
  double      ni = In_Ptr->layerspecs[layer].n;
  double      nt = In_Ptr->layerspecs[layer + 1].n;

  /* Get r. */
  if (uz <= In_Ptr->layerspecs[layer].cos_crit1)
    r = 1.0;			/* total internal reflection. */
  else
    r = RFresnel(ni, nt, uz, &uz1);

#if PARTIALREFLECTION
  if (layer == In_Ptr->num_layers && r < 1.0) {
    Photon_Ptr->uz = uz1;
    RecordT(r, In_Ptr, Photon_Ptr, Out_Ptr);
    Photon_Ptr->uz = -uz;
  } else if (RandomNum > r) {	/* transmitted to layer+1. */
    Photon_Ptr->layer++;
    Photon_Ptr->ux *= ni / nt;
    Photon_Ptr->uy *= ni / nt;
    Photon_Ptr->uz = uz1;
  } else			/* reflected. */
    Photon_Ptr->uz = -uz;
#else
  if (RandomNum > r) {		/* transmitted to layer+1. */
    if (layer == In_Ptr->num_layers) {
      Photon_Ptr->uz = uz1;
      RecordT(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
      Photon_Ptr->alive = 0;	/* escaped. */
    } else {
      Photon_Ptr->layer++;
      Photon_Ptr->ux *= ni / nt;
      Photon_Ptr->uy *= ni / nt;
      Photon_Ptr->uz = uz1;
    }
  } else			/* reflected. */
    Photon_Ptr->uz = -uz;
#endif
}

/**************************************************************************
 *	Set a step size if the previous step has finished.
 *
 *	If the step size fits in the current layer, move the photon,
 *	drop some weight, choose a new photon direction for propagation.
 *
 *	If the step size is long enough for the photon to
 *	hit an interface, this step is divided into three steps.
 *	First, move the photon to the boundary free of
 *	absorption or scattering.
 *	Second, update the step size to the unfinished step size.
 *	Third, decide whether the photon is reflected or transmitted.
 ****/
void
HopDropSpin(InStru * In_Ptr,
	    PhotonStru * Photon_Ptr,
	    OutStru * Out_Ptr)
{
  LayerStru   layer_struct;
  double      mut;
  double      path;		/* distance between photon and boundary
				 * cm. */

  layer_struct = In_Ptr->layerspecs[Photon_Ptr->layer];
  mut = layer_struct.mua + layer_struct.mus;
  SetStepSize(Photon_Ptr);
  path = PathToBoundary(Photon_Ptr, In_Ptr);
  if (path * mut <= Photon_Ptr->s) {	/* hit boundary. */
    Hop(Photon_Ptr, path, layer_struct.n);	/* move to boundary plane. */
    Photon_Ptr->s -= path * mut;/* update s. */
    if (Photon_Ptr->uz < 0.0)
      CrossUpOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
    else
      CrossDnOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
  } else {			/* fit in layer. */
    Hop(Photon_Ptr, Photon_Ptr->s / mut, layer_struct.n);
    Photon_Ptr->s = 0;		/* update s. */
    Drop(In_Ptr, Photon_Ptr, Out_Ptr);
    Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g, Photon_Ptr);
    Photon_Ptr->scatters++;
  }
}

/**************************************************************************
 *	Trace a photon, then compute the 0D constants including A, R, T and
 *	their standard errors.
 ****/
void
TracePhoton(InStru * In_Ptr,
	    PhotonStru * Photon_Ptr,
	    OutStru * Out_Ptr)
{
  do {
    HopDropSpin(In_Ptr, Photon_Ptr, Out_Ptr);
    if (Photon_Ptr->alive && Photon_Ptr->w < In_Ptr->Wth)
      Roulette(Photon_Ptr);
  } while (Photon_Ptr->alive);

  Out_Ptr->A += Out_Ptr->Ai;
  Out_Ptr->Ae += Out_Ptr->Ai * Out_Ptr->Ai;

  Out_Ptr->Tb += Out_Ptr->Tbi;
  Out_Ptr->Tbe += Out_Ptr->Tbi * Out_Ptr->Tbi;
  Out_Ptr->Td += Out_Ptr->Tdi;
  Out_Ptr->Tde += Out_Ptr->Tdi * Out_Ptr->Tdi;

  Out_Ptr->Rb += Out_Ptr->Rbi;
  Out_Ptr->Rbe += Out_Ptr->Rbi * Out_Ptr->Rbi;
  Out_Ptr->Rd += Out_Ptr->Rdi;
  Out_Ptr->Rde += Out_Ptr->Rdi * Out_Ptr->Rdi;
}
