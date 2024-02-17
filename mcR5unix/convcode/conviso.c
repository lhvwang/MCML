/****
 *  This files plot a series of contour lines for a 2D array of double.
 *  4/14/92.
 ****/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#ifndef PI
#define PI 3.1415926
#endif

/****
 *	Return true is val is between Z[i][j] and Z[i+1][j].
 ****/
#define INTERSECTI(val, Z, i, j) \
  (Z[i][j] <= val && val <= Z[i+1][j] || \
   Z[i][j] >= val && val >= Z[i+1][j])

/****
 *	Return true is val is between Z[i][j] and Z[i][j+1].
 ****/
#define INTERSECTJ(val, Z, i, j) \
  (Z[i][j] <= val && val <= Z[i][j+1] || \
   Z[i][j] >= val && val >= Z[i][j+1])

FILE       *GetWriteFile(char *Ext);	/* in conho.c. */

struct XY {
  double      x, y;
};

struct XYpairs {
  double      x, y;
  struct XYpairs *next;
};

typedef struct XYpairs *PairList;

/* set up a que for printing the pairs. */
typedef PairList QDATA;

struct Qlinked_list {
  QDATA       d;
  struct Qlinked_list *next;
};

typedef struct Qlinked_list QELEMENT;
typedef QELEMENT *QLINK;

typedef struct {		/* a que of pairs. */
  QLINK       front, rear;
}           QUE;
/* end of que definition. */


struct Isolines {
  double      iso_val;		/* z value of the contour line. */
  PairList    pairs;
  struct Isolines *next;
};

typedef struct Isolines *IsoList;



double 
ZMin(double **Z,
     long IXmax,
     long IYmax,
     double Dx,
     double Dy,
     struct XY * Pmin_Ptr)
{
  double      zmin;
  long        i, j;

  zmin = Z[0][0];
  for (i = 0; i < IXmax; i++)
    for (j = 0; j < IYmax; j++)
      if (Z[i][j] < zmin) {
	zmin = Z[i][j];
	Pmin_Ptr->x = i * Dx;
	Pmin_Ptr->y = j * Dy;
      }
  return (zmin);
}


double 
ZMax(double **Z,
     long IXmax,
     long IYmax,
     double Dx,
     double Dy,
     struct XY * Pmax_Ptr)
{
  double      zmax;
  long        i, j;

  zmax = Z[0][0];
  for (i = 0; i < IXmax; i++)
    for (j = 0; j < IYmax; j++)
      if (Z[i][j] > zmax) {
	zmax = Z[i][j];
	Pmax_Ptr->x = i * Dx;
	Pmax_Ptr->y = j * Dy;
      }
  return (zmax);
}


/* Get the isovalues from user. */
IsoList 
GetIsoValues(double Z_Min, double Z_Max)
{
  IsoList     head;
  char        in_str[256];

  printf("Input an isovalue or . to stop: ");
  scanf("%s", in_str);
  if (strlen(in_str) == 1 && in_str[0] == '.')
    return (NULL);

  head = (IsoList) malloc(sizeof(struct Isolines));
  if (head == NULL)
    return (NULL);

  /* get the elements for the node. */
  sscanf(in_str, "%lf", &head->iso_val);
  if (head->iso_val < Z_Min)
    head->iso_val = Z_Min;
  else if (head->iso_val > Z_Max)
    head->iso_val = Z_Max;
  head->pairs = NULL;
  head->next = GetIsoValues(Z_Min, Z_Max);
  return (head);
}


/****
 *  The isoposition is between [i][j] & [i+1][j]
 *
 *  Linear interpolation is used to locate the y component of
 *  the iso position.
 ****/
void 
IsoPositionI(PairList pair,
	     double iso_val, double **Z,
	     long i, long j,
	     double Dx, double Dy)
{
  pair->y = (j + 0.5) * Dy;
  if (Z[i][j] != Z[i + 1][j])
    pair->x = (i + 0.5) * Dx + Dx * (iso_val - Z[i][j]) / (Z[i + 1][j] - Z[i][j]);
  else
    pair->x = (i + 1) * Dx;	/* take the mid point. */
}


/****
 *  The isoposition is between [i][j] & [i][j+1]
 *
 *  Linear interpolation is used to locate the y component of
 *  the iso position.
 ****/
void 
IsoPositionJ(PairList pair,
	     double iso_val, double **Z,
	     long i, long j,
	     double Dx, double Dy)
{
  pair->x = (i + 0.5) * Dx;
  if (Z[i][j + 1] != Z[i][j])
    pair->y = (j + 0.5) * Dy + Dy * (iso_val - Z[i][j]) / (Z[i][j + 1] - Z[i][j]);
  else
    pair->y = (j + 1) * Dy;	/* take the mid point. */
}

void 
IsoPosition(PairList pair,
	    double iso_val, double **Z,
	    long i, long j,
	    double Dx, double Dy)
{

  if (INTERSECTI(iso_val, Z, i, j))
    IsoPositionI(pair, iso_val, Z, i, j, Dx, Dy);
  else if (INTERSECTJ(iso_val, Z, i, j))
    IsoPositionJ(pair, iso_val, Z, i, j, Dx, Dy);
}

PairList 
AllocPair(void)
{
  PairList    pair;

  pair = (PairList) malloc(sizeof(struct XYpairs));
  if (pair == NULL)
    puts("...malloc error. Contour lines are not complete");
  return (pair);
}

void 
GetAnIsoLine(IsoList IsoNode, double **Z,
	     long IXmax, long IYmax,
	     double Dx, double Dy)
{
  long        i, j;
  double      ival = IsoNode->iso_val;
  PairList    pair_tail;

  for (j = 0; j < IYmax - 1; j++)
    for (i = 0; i < IXmax - 2; i++)
      if (INTERSECTI(ival, Z, i, j) ||
	  INTERSECTJ(ival, Z, i, j)) {	/* found a pair. */
	if (IsoNode->pairs == NULL) {	/* 1st pair. */
	  if (!(IsoNode->pairs = AllocPair()))
	    return;
	  pair_tail = IsoNode->pairs;
	} else {		/* subsequent  pairs. */
	  if (!(pair_tail->next = AllocPair()))
	    return;
	  pair_tail = pair_tail->next;
	}
	IsoPosition(pair_tail, ival, Z, i, j, Dx, Dy);
      }
  pair_tail->next = NULL;	/* end of the pair list. */
}


/****
 *	Return the quadrant.
 *	For Frz or Arz, the peak point is usually on z axis which is
 *	y here. We want start the isoline from the 4th quadrant.
 *	Therefore, we use -1 for the 4th quadrant instead of 4.
 ****/
short 
GetQuadrant(double x, double y)
{
  if (x > 0 && y >= 0)
    return (1);			/* Include +x-axis. */
  else if (x <= 0 && y > 0)
    return (2);			/* Include +y-axis. */
  else if (x < 0 && y <= 0)
    return (3);			/* Include -x-axis. */
  else if (x >= 0 && y < 0)
    return (-1);		/* Include -y-axis. */
}

/****
 *  Compare the angle wrt Pmax.  If the angle of "This" is larger
 *  than that of the "Next", return 1.  Otherwise, return 0.
 *
 *  Although atan2 returns value in the range -pi to pi, we want
 *	to avoid it because it is slow.  We compare the quadrants of
 *	the two points first.  If they are in the same quadrant, we
 *	compare the relative positions.
 ****/
char 
OutOfOrder(PairList This,	/* current pair. */
	   PairList Next,
	   struct XY Pmax)
{
  double      x0, y0, x1, y1;
  short       q0, q1;		/* Quadrants. */
  char        out_order;

  if (This == NULL || Next == NULL)
    return (0);			/* This shouldn't happen. */
  x0 = This->x - Pmax.x;
  y0 = This->y - Pmax.y;
  q0 = GetQuadrant(x0, y0);
  x1 = Next->x - Pmax.x;
  y1 = Next->y - Pmax.y;
  q1 = GetQuadrant(x1, y1);

  if (q0 < q1)
    out_order = 0;
  else if (q0 > q1)
    out_order = 1;
  else {			/* In the same quadrant. */
    if (y0 * x1 < y1 * x0)
      out_order = 0;
    else
      out_order = 1;
  }

  return (out_order);
}


/****
 *  Sort the isopositions according to the angle with respect to
 *  the position of the maximum value.
 ****/
void 
SortAnIsoLine(PairList * PairHeadPtr,
	      struct XY Pmax)
{
  char        sorted = 0;
  PairList    this, last;	/* this=the pair being compared w/ the
				 * next one. */
  PairList    sorted_head = NULL;	/* the head to the sublist of
					 * sorted pairs. */

  while (!sorted && sorted_head != *PairHeadPtr) {
    sorted = 1;			/* assume sublist is sorted. */
    last = NULL;
    this = *PairHeadPtr;	/* sublist starts at *PairHeadPtr, ends at
				 * sorted_head. */

    while (this->next != sorted_head) {	/* traverse the sublist. */
      if (OutOfOrder(this, this->next, Pmax)) {	/* swap and move to next. */
	sorted = 0;
	if (last == NULL)
	  *PairHeadPtr = this->next;
	else
	  last->next = this->next;

	last = this->next;
	this->next = last->next;
	last->next = this;
      } else {			/* move to next. */
	last = this;
	this = this->next;
      }
    }				/* end of sublist traversing. */

    sorted_head = this;
  }				/* end of big while. */
}

#if 0
/****
 *	Repeat the first pair of the isoline at the end of the list
 *	so that the isoline is a loop.
 *	Sometimes this may not be desired.
 ****/
void 
LoopAnIsoLine(PairList * PairHeadPtr)
{
  PairList    tail;

  tail = *PairHeadPtr;
  while (tail->next != NULL)
    tail = tail->next;

  tail->next = (PairList) malloc(sizeof(struct XYpairs));
  if (tail->next == NULL)
    return;
  tail = tail->next;
  tail->x = (*PairHeadPtr)->x;
  tail->y = (*PairHeadPtr)->y;
  tail->next = NULL;
}
#endif

void 
GetIsoLines(IsoList isos,
	    double **Z,
	    long IXmax,
	    long IYmax,
	    double Dx,
	    double Dy,
	    struct XY Pmax)
{
  IsoList     node = isos;

  while (node != NULL) {
    GetAnIsoLine(node, Z, IXmax, IYmax, Dx, Dy);
    SortAnIsoLine(&node->pairs, Pmax);
    /* LoopAnIsoLine(&node->pairs); */
    node = node->next;
  }
}

char 
IsEmpty(QUE q)
{
  return (q.front == NULL);
}


void 
Deque(QUE * q, QDATA * x)
{
  QLINK       temp = q->front;

  if (!IsEmpty(*q)) {
    *x = temp->d;
    q->front = temp->next;
    free(temp);
  } else
    printf("Empty que.\n");
}

void 
Enque(QUE * q, QDATA x)
{
  QLINK       temp;

  temp = (QLINK) malloc(sizeof(QELEMENT));
  temp->d = x;
  temp->next = NULL;
  if (IsEmpty(*q))
    q->front = q->rear = temp;
  else {
    q->rear->next = temp;
    q->rear = temp;
  }
}


void 
WriteIsoLines(FILE * Isofile,
	      IsoList IsoHead)
{
  QUE         qprint = {NULL, NULL};	/* que to be printed. */
  QUE         qsave = {NULL, NULL};
  QDATA       p;		/* QDATA is PairList. */
  char        all_nulls;

  while (IsoHead != NULL) {
    fprintf(Isofile, "X%-8lG\tY%-8lG\t", IsoHead->iso_val, IsoHead->iso_val);
    Enque(&qprint, IsoHead->pairs);
    IsoHead = IsoHead->next;
  }
  fprintf(Isofile, "\n");

  do {
    all_nulls = 1;
    while (!IsEmpty(qprint)) {
      Deque(&qprint, &p);
      if (p != NULL) {
	fprintf(Isofile, "%9.2lE\t%9.2lE\t", p->x, p->y);
	Enque(&qsave, p->next);	/* add the next to a new que for a new
				 * row. */
	if (all_nulls == 1)
	  all_nulls = 0;
      } else {
	/* fprintf(Isofile, "%9s\t%9s\t", " ", " "); */
	fprintf(Isofile, "\t\t");
	Enque(&qsave, p);	/* add a null to the new que. */
      }
    }
    fprintf(Isofile, "\n");
    qprint = qsave;
    qsave.front = NULL;
    qsave.rear = NULL;
  } while (!all_nulls);
}

void 
FreePairs(PairList Pairs)
{
  PairList    this_pair;

  while (Pairs != NULL) {
    this_pair = Pairs;
    Pairs = Pairs->next;
    free(this_pair);
  }
}

void 
FreeIsoLines(IsoList IsoHead)
{
  IsoList     this_iso;
  PairList    pair_head;

  while (IsoHead != NULL) {
    FreePairs(IsoHead->pairs);
    this_iso = IsoHead;
    IsoHead = IsoHead->next;
    free(this_iso);
  }
}


void 
IsoPlot(double **Z,		/* the 2D array Z[i][j]. */
	long int IXmax,
	long int IYmax,		/* the 0<=i<=IXmax, 0<=j<=IYmax. */
	double Dx,
	double Dy)
{				/* the gridline separations. */
  FILE       *isofile;		/* send isolines to this file. */
  struct XY   pmin, pmax;	/* the xy positions of the min & max. */
  double      zmin, zmax;
  IsoList     isos;
  char        fname[128] = "iso";

  /* Locate the min & max. */
  zmin = ZMin(Z, IXmax, IYmax, Dx, Dy, &pmin);
  zmax = ZMax(Z, IXmax, IYmax, Dx, Dy, &pmax);

  if (zmin == zmax || IXmax * IYmax < 4) {
    printf("...Not enough data for contour plot\n");
    return;
  }
  isofile = GetWriteFile(fname);
  if (isofile == NULL)
    return;

  printf("The range of the value is %lf to %lf.\n", zmin, zmax);
  isos = GetIsoValues(zmin, zmax);

  GetIsoLines(isos, Z, IXmax, IYmax, Dx, Dy, pmax);
  WriteIsoLines(isofile, isos);

  FreeIsoLines(isos);

  fclose(isofile);
}
