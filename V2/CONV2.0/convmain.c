/**************************************************************************
 *  	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992.
 *
 *	Convolution program for Monte Carlo simulation of photon
 *	distribution in multilayered turbid media.
 ****/

#include "conv.h"

#define EPS 0.1			/* default relative error in convolution. */

void        ReadMcoFile(InStru *, OutStru *, ConvStru *);
void        ExtractRAT(InStru *, OutStru *, ConvStru *);
void        ExtractOrigData(InStru *, OutStru *, ConvStru *);
void        ExtractConvData(InStru *, OutStru *, ConvStru *);
void        LaserBeam(BeamStru *);
void        ConvError(float *);

/**************************************************************************
 ****/
void
ShowMainMenu(void)
{
  puts("  a = About CONV.");
  puts("  i = Input a file of MCML output.");
  puts("  r = Reflectance, absorption, and transmittance.");
  puts("  o = Extract original data.\n");

  puts("  b = Specify laser beam.");
  puts("  e = Specify convolution error. ");
  puts("  c = Extract convolved data.");
  puts("  q = Quit from the program.");
  puts("  * Commands in conv are not case-sensitive");
}

/**************************************************************************
 ****/
void
QuitProgram(void)
{
  char        cmd_str[STRLEN];

  printf("Do you really want to exit CONV? (y/n): ");
  gets(cmd_str);
  if (toupper(cmd_str[0]) == 'Y') {	/* really quit. */
    exit(0);
  }
}

/**************************************************************************
 *      Center a string according to the column width.
 ****/
#define COLWIDTH 80
void
CtrPuts(char *InStr)
{
  short      nspaces;		/* number of spaces to be left-filled. */
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
 *      Print messages about conv.
 ****/
void
AboutConv(void)
{
  CtrPuts(" ");
  CtrPuts("CONV 2.0, Copyright (c) 1992-1996\n");
  /* \251 is the formal copyright sign but does not show up on some
  terminals.*/
  CtrPuts("Convolution of MCML Simulation Data\n");

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
}

/**************************************************************************
 *	main commands.
 ****/
void
BranchMainCmd(char *Cmd,	/* Cmd is command char. */
	       InStru * In_Ptr, OutStru * Out_Ptr, ConvStru * Conv_Ptr)
{
  switch (toupper(Cmd[0])) {
case 'A':
    AboutConv(); 
    break;
  case 'I':
    ReadMcoFile(In_Ptr, Out_Ptr, Conv_Ptr);
    break;
  case 'R':
    ExtractRAT(In_Ptr, Out_Ptr, Conv_Ptr);
    break;
  case 'O':
    ExtractOrigData(In_Ptr, Out_Ptr, Conv_Ptr);
    break;
  case 'B':
    LaserBeam(&Conv_Ptr->beam);
    break;
  case 'E':
    ConvError(&Conv_Ptr->eps);
    break;  
  case 'C':
    ExtractConvData(In_Ptr, Out_Ptr, Conv_Ptr);
    break; 
  case 'H':
    ShowMainMenu();
    break;  
  case 'Q':
    QuitProgram();
    break;
  default:
    puts("...Wrong command");
  }
}

/**************************************************************************
 ****/
int
main(void)
{
  InStru      in_parm;
  OutStru     out_parm;
  ConvStru    conv_parm;
  char        cmd_str[STRLEN];

  puts(" ");
  CtrPuts("CONV Version 2.0, Copyright (c) 1992-1996\n"); 
  conv_parm.eps = EPS;
  conv_parm.beam.type = original;
  conv_parm.datain = 0;         /* data is not read in yet. */

  do {
    printf("\n> Main menu (h for help) => ");
    do				/* get the command input. */
      gets(cmd_str);
    while (!strlen(cmd_str));
    BranchMainCmd(cmd_str, &in_parm, &out_parm, &conv_parm);
  } while (1);
}
