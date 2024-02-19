/**************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *	Program name: MCML.
 *	The main program for Monte Carlo simulation of light transport
 *	in multi-layered turbid media.
 ****/

/****
 *	THINKCPROFILER is defined to generate profiler calls in
 *	Think C. If 1, remember to turn on "Generate profiler
 *	calls" in the options menu.
 ****/
#define THINKCPROFILER 0

/* GNU cc does not support difftime() and CLOCKS_PER_SEC. */
#define GNUCC 0

#if THINKCPROFILER
#include <profile.h>
#include <console.h>
#endif

#include "mcml.h"

/* Declare before they are used in main(). */
void        InitOutputData(InStru *, OutStru *);
void	    CtrPuts(char *);
char       *FindDataLine(FILE *);
void        InterReadParam(InStru *);
void        LaunchPhoton(double, InStru *, OutStru *, PhotonStru *);
void        AboutMCML(void);
void        IOResult(FILE *, InStru *, OutStru *, char);
void        CheckParamFromFile(FILE *, InStru *);
void        ReadParam(FILE *, InStru *);
void        ScaleResult(InStru *, OutStru *, char);
double      Rspecular(LayerStru *);
Boolean     ReadMediumListQ(FILE *, InStru *);
void        TracePhoton(InStru *, PhotonStru *, OutStru *);
Boolean     RunChangedInput(InStru *);
FILE       *GetFile(char *, char *);
void        FreeData(InStru *, OutStru *);
Boolean     ReadNumPhotonsQ(FILE *, InStru *, char);
Boolean     CheckFileVersionQ(FILE *, char *);

#define SWAP(x, y) {long temp; temp = x; x = y; y = temp;}

/**************************************************************************
 *  If F = 0, reset the clock and return 0.
 *
 *  If F = 1, pass the real time to Msg, print Msg on the
 *  screen, return the real time elapsed since F=0.
 *
 *  If F = 2, same as F=1 except no printing.
 *
 *  Note that clock() and time() return user time and real
 *  time respectively.  User time is whatever the system allocates
 *  to the running of the program; real time is wall-clock time.
 *  In a time-shared system, they need not be the same.
 *
 *  clock() only hold 16 bit integer, which is about 32768
 *  clock ticks.  Because this fact can lead to wrap-around, it
 *  is dangerous to use clock() in Monte Carlo simulations.
 *  Therefore, we only keep track of the real time.
 *
 *  On UNIX machines, users can use time command to get the
 *  user time.  On personal computers, the user time is equal
 *  to the real time unless the system is multi-tasking during
 *  the Monte Carlo simulation, which usually is not true.
 ****/
time_t
PunchTime(char F, char *Msg, InStru * In_Ptr)
{
#if GNUCC
  return (0);
#else
  static time_t rt0;		/* real time reference. */
  double      real_time_secs;

  if (F == 0) {
    rt0 = time(NULL);
    return (0);
  }
  real_time_secs = difftime(time(NULL), rt0);

  if (F == 1) {			/* show & pass real time. */
    real_time_secs += In_Ptr->add_num_seconds;
    sprintf(Msg, "Real time for this simulation: %.1f sec = %.2f hr.",
	    real_time_secs, real_time_secs / 3600);
  }
  return (real_time_secs);
#endif
}

/**************************************************************************
 *	Print the current time and the estimated finishing time.
 *
 *	P1 is the number of computed photon packets.
 *	Pt is the total number of photon packets.
 ****/
void
PredictDoneTime(long P1, InStru * In_Ptr)
{
  time_t      now, done_time;
  struct tm  *date;
  char        s[80];

  now = time(NULL);
  date = localtime(&now);
  strftime(s, 80, "%H:%M %a %m/%d/%Y", date);
  printf("%s\t", s);

  if (In_Ptr->control_bit == 1 || In_Ptr->control_bit == 3) {
    done_time = now +
      (time_t) ((PunchTime(2, "", In_Ptr) / (double) P1)
		* (In_Ptr->num_photons - P1));
    if (In_Ptr->control_bit == 3)
      done_time = (done_time < (now + In_Ptr->num_seconds
				- PunchTime(2, "", In_Ptr)))
	? done_time : (now + In_Ptr->num_seconds - PunchTime(2, "", In_Ptr));
    date = localtime(&done_time);
    strftime(s, 80, "%H:%M %a %m/%d/%Y", date);
    printf("%s", s);
  }
  printf("\n");
}

/**************************************************************************
 *	Generate a string representing the user-specified done time.
 ****/
void
FormDateString(char *String, InStru * In_Ptr)
{
  time_t      now, done_time;
  struct tm  *date;

  now = time(NULL);
  done_time = now + In_Ptr->num_seconds - PunchTime(2, "", In_Ptr);
  date = localtime(&done_time);
  strftime(String, 80, "%H:%M %x", date);
}

/**************************************************************************
 *	Report how and when the simultion will be terminated.
 ****/
void
ReportControlInfo(short NumRunsLeft, InStru * In_Ptr)
{
  char        string[STRLEN];

  printf("%d runs left.", NumRunsLeft);
  if (In_Ptr->control_bit == 1) {
    printf("\nThe simulation will terminate after tracing %ld photons.\n\n",
	   In_Ptr->num_photons);
    printf("\tPhotons Done\tCurrent Time\t\tEstimated Done Time\n");
    printf("\t------------\t--------------------\t--------------------\n");
  } else if (In_Ptr->control_bit == 2) {
    FormDateString(string, In_Ptr);
    printf("\nThe simulation will terminate at %s.\n\n", string);
    printf("\tPhotons Done\tCurrent Time\n");
    printf("\t------------\t--------------------\n");
  } else {
    FormDateString(string, In_Ptr);
    printf("\nThe simulation will terminate after tracing %ld photons\n",
	   In_Ptr->num_photons);
    printf("or at %s, whichever comes first.\n\n", string);
    printf("\tPhotons Done\tCurrent Time\t\tEstimated Done Time\n");
    printf("\t------------\t--------------------\t--------------------\n");
  }
}

/**************************************************************************
 *	Report the estimated time, number of photons and runs left
 *	after calculating 10 photons or every 1/10 of total
 *	number of photons.
 *
 *	Pi is the number of traced photons so far.
 ****/
void
ReportStatus(long Pi, InStru * In_Ptr)
{
  if (In_Ptr->control_bit == 1) 
    printf("%11ld(%6.2f%%)\t", Pi, (float) Pi * 100 / In_Ptr->num_photons);
  else if (In_Ptr->control_bit == 2) 
    printf("\t%12ld\t", Pi);
  else 
    printf("%11ld(%6.2f%%)\t", Pi, (float) Pi * 100 / In_Ptr->num_photons);

  PredictDoneTime(Pi, In_Ptr);
}

/**************************************************************************
 *	Report time, photon number traced and write results.
 ****/
void
ReportResult(InStru * In_Ptr, OutStru * Out_Ptr)
{
  char        time_report[STRLEN];
  FILE       *fp;

  printf("\n%ld photons have been traced in this simulation.\n",
	 In_Ptr->num_photons);

  PunchTime(1, time_report, In_Ptr);
  puts(time_report);
  puts(" ");

  ScaleResult(In_Ptr, Out_Ptr, 0);

  if ((fp = fopen(In_Ptr->out_fname, "w")) == NULL) {
    printf("Can not open output file to write.\n");
    exit(1);
  }
  IOResult(fp, In_Ptr, Out_Ptr, 1);
}

/**************************************************************************
 *	Get the file name of the input data file from the
 *	argument to the command line.
 ****/
void
GetFnameFromArgv(int argc, char *argv[], char *input_filename)
{
  if (argc >= 2) {		/* filename in command line */
    strcpy(input_filename, argv[1]);
  } else
    input_filename[0] = '\0';
}

/**************************************************************************
 *	Execute Monte Carlo simulation for one independent run.
 *      Type = 0, start a new simulation;
 *      Type = 1, continue previous simulation.
 ****/
void
DoOneRun(short NumRunsLeft,
	 InStru * In_Ptr,
	 OutStru * Out_Ptr,
	 char Type)
{
  PhotonStru  photon;
  register long i_photon = 1;	/* photon number traced.  */
  Boolean     exit_switch = 0;	/* switch to terminate simulation. */
  static int  tens;

#if THINKCPROFILER
  InitProfile(200, 200);
  cecho2file("prof.rpt", 0, stdout);
#endif

  if (Type == 0) {		/* start a new simulation. */
    if (In_Ptr->slayer == 0)
      Out_Ptr->Rsp = Rspecular(In_Ptr->layerspecs);
    RandomGen(0, 1, NULL);	/* initialize the generator. */
  }
  PunchTime(0, "", In_Ptr);
  ReportControlInfo(NumRunsLeft, In_Ptr);
  tens = 10;

  do {
    LaunchPhoton(Out_Ptr->Rsp, In_Ptr, Out_Ptr, &photon);
    TracePhoton(In_Ptr, &photon, Out_Ptr);

    if (i_photon == tens) {	/* report status every ten photons. */
      tens *= 10;
      ReportStatus(i_photon, In_Ptr);
    }
    i_photon++;
    if (In_Ptr->control_bit == 1)
      exit_switch = (i_photon > In_Ptr->num_photons);
    else if (In_Ptr->control_bit == 2)
      exit_switch = (PunchTime(2, "", In_Ptr) >= In_Ptr->num_seconds);
    else
      exit_switch = (i_photon > In_Ptr->num_photons) ||
	(PunchTime(2, "", In_Ptr) >= In_Ptr->num_seconds);
  } while (!exit_switch);

#if THINKCPROFILER
  exit(0);
#endif

  In_Ptr->num_photons = In_Ptr->add_num_photons + i_photon - 1;
  In_Ptr->num_seconds = In_Ptr->add_num_seconds + PunchTime(2, "", In_Ptr);
  In_Ptr->control_bit = 3;

  ReportResult(In_Ptr, Out_Ptr);
  FreeData(In_Ptr, Out_Ptr);
}

/**************************************************************************
 *	In continuation runs, ask for additional number of photons or time.
 ****/
void
AddNumPhotons(InStru * In_Ptr)
{
  printf("\n%ld photons have been traced in the previous simulation.",
	 In_Ptr->num_photons);
  printf("\nSpecify additional photons or compution time in hh:mm format,");
  printf("\nor both in one line (e.g. 10000 5:30): ");

  while (!ReadNumPhotonsQ(stdin, In_Ptr, 1))
    printf("Input agian: ");

  printf("\n");
}

/**************************************************************************
 *	Start a new simulation non-interactively,
 *      Read input parameter from file input_filename.
 ****/
void
NonInterSimu(FILE *Fp, InStru * In_Ptr, OutStru * Out_Ptr)
{
  short         num_runs_left;

  CheckParamFromFile(Fp, In_Ptr);
  num_runs_left = In_Ptr->num_runs;
  while (num_runs_left--) {
    ReadParam(Fp, In_Ptr);
    InitOutputData(In_Ptr, Out_Ptr);
    DoOneRun(num_runs_left, In_Ptr, Out_Ptr, 0);
  }
  fclose(Fp);
  exit(0);
}

/**************************************************************************
 *	Read input parameters from a file with interactive change.
 ****/
void
FileInterSimu(InStru * In_Ptr, OutStru * Out_Ptr)
{
  char        input_filename[STRLEN];
  FILE       *input_file_ptr;

  if ((input_file_ptr = GetFile(input_filename, "mcmli2.0")) != NULL)
    if (ReadMediumListQ(input_file_ptr, In_Ptr)) {
      ReadParam(input_file_ptr, In_Ptr);
      printf("The parameters of the first run have been read in.\n");
      if (RunChangedInput(In_Ptr)) {
	InitOutputData(In_Ptr, Out_Ptr);
	DoOneRun(0, In_Ptr, Out_Ptr, 0);
	fclose(input_file_ptr);
	exit(0);
      }
      fclose(input_file_ptr);
    }
}

/**************************************************************************
 *	Continute a previous simulation.
 ****/
void
ContinueSimu(InStru * In_Ptr, OutStru * Out_Ptr)
{
  char        input_filename[STRLEN];
  FILE       *fp;

  printf("Specify the output file name of a previous simulation. \n");
  if ((fp = GetFile(input_filename, "mcmloA2.0")) == NULL)
    return;

  FindDataLine(fp);		/* skip the line of file version. */
  if (!ReadMediumListQ(fp, In_Ptr))
    exit(1);
  ReadParam(fp, In_Ptr);

  AddNumPhotons(In_Ptr);
  InitOutputData(In_Ptr, Out_Ptr);
  IOResult(fp, In_Ptr, Out_Ptr, 0);
  ScaleResult(In_Ptr, Out_Ptr, 1);

  SWAP(In_Ptr->num_photons, In_Ptr->add_num_photons);
  SWAP(In_Ptr->num_seconds, In_Ptr->add_num_seconds);
  DoOneRun(0, In_Ptr, Out_Ptr, 1);
  exit(0);
}

/**************************************************************************
****/
QuitProgram(void)
{
  char        cmd_str[STRLEN];

  printf("Do you really want to quit MCML? (y/n): ");
  gets(cmd_str);
  if (toupper(cmd_str[0]) == 'Y')      /* really quit. */
    exit(0);
}

/**************************************************************************
 ****/
void
ShowMainMenu(void)
{
  puts("  a = About MCML.");
  puts("  r = Run an input file non-interactively.");
  puts("  m = Input and modify parameters of a file (the first run only).");
  puts("  i = Input parameters interactively.");
  puts("  c = Continue a previous simulation.");
  puts("  q = Quit from the program.");
  puts("  * Commands here are not case-sensitive.");
}

/**************************************************************************
 ****/
void
BranchMainMenu(char *string, InStru * In_Ptr, OutStru * Out_Ptr)
{
  char        input_filename[STRLEN];
  FILE        *input_file_ptr;

  switch (toupper(string[0])) {
  case 'A':
    AboutMCML();
    break;

  case 'R':			/* non-interactive. */
    if ((input_file_ptr = GetFile(input_filename, "mcmli2.0")) != NULL) 
      NonInterSimu(input_file_ptr, In_Ptr, Out_Ptr);
    break;

  case 'M':			/* read a file with an interactive change. */
    FileInterSimu(In_Ptr, Out_Ptr);
    break;

  case 'I':			/* interactive. */
    InterReadParam(In_Ptr);
    if (RunChangedInput(In_Ptr)) {
      InitOutputData(In_Ptr, Out_Ptr);
      DoOneRun(0, In_Ptr, Out_Ptr, 0);
      exit(0);
    }
    break;

  case 'C':
    ContinueSimu(In_Ptr, Out_Ptr);
    break;

  case 'H':
    ShowMainMenu();
    break;

  case 'Q':
    QuitProgram();
    break;

  default:
    puts("...Unknown command");
  }
}

/**************************************************************************
 *	The argument to the command line is the input filename, if any.
 *	Macintosh does not support command line.
 ****/
int
main(int argc, char *argv[])
{
  InStru    in_param;
  OutStru   out_param;
  char        input_filename[STRLEN], string[STRLEN];
  FILE        *input_file_ptr;

  puts(" ");
  CtrPuts("MCML Version 2.0, Copyright (c) 1992-1996\n");
  /* \251 is the formal copyright sign but does not show up on some
  terminals.*/

  if (argc >= 2) {		/* non-interactive. */
    GetFnameFromArgv(argc, argv, input_filename);
    if((input_file_ptr = fopen(input_filename, "r")) != NULL) 
      if (CheckFileVersionQ(input_file_ptr, "mcmli2.0")) 
	NonInterSimu(input_file_ptr, &in_param, &out_param);
    exit(0);
  } else			/* accept commands from console. */
    while (1) {
      do {
	printf("\n> Main menu (h for help) => ");
	gets(string);
      } while (!strlen(string));

      BranchMainMenu(string, &in_param, &out_param);
    }
}
