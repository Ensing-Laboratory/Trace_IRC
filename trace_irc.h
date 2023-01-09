#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h> 

#define PRINT_TCL_BRACKET 0            /* if 1: print stuff to trace_out.tcl for gopenmol */
#define PRINT_TCL_MINIMIZE 0           /* if 1: print stuff to trace_out.tcl for gopenmol */
#define PRINT_TCL_CLEANPOINTS 1        /* if 1: print stuff to trace_out.tcl for gopenmol */
#define PRINT_TCL_SMOOTHPOINTS 1       /* if 1: print stuff to trace_out.tcl for gopenmol */
#define PRINT_MC_MOVES 1               /* if 1: print Monte Carlo points to montecarlo.out */
#define PRINT_MC_MOVES_TCL 1           /* if 1: print Monte Carlo points to trace_out.tcl */

#define INPUTFILE "trace_irc.in"
#define COLVARFILE "colvar_val"
#define DCOLVARFILE "dcolvar_val"
#define OUTPUTFILE "trace_coursegrid.out" /* irc on course grid      */
#define OUTPUTFILE2 "trace_irc.tcl"       /* bracket and optim. info */
#define OUTPUTFILE3 "umbrella.pot"        /* irc on small grid       */
#define OUTPUTFILE4 "trace_irc.scratch"   /* course vectorplot for gnuplot and */
                                          /*           integral perp2irc info  */
#define MCOUTPUTFILE "montecarlo.out"     /* Monte Carlo info  */
#define MCHISTOOUT "montecarlo.histo"     /* Monte Carlo histogram output  */
#define PRINTFORMAT "  %lf"     /* #define PRINTFORMAT "  %16.10f"  */
#define NVARMAX 13
#define MAXLINE 164
#define MAXCYCLE 200
#define NHILLSMAX 50000
#define NPOINTSMAX 100

#define TOL  0.000001
#define ETOL 0.0000001
#define DISTTOL 0.03    /* distance tolerance between point i and i-2 to decide if */
                        /* we are going back. A good value was sofar: 0.01 */

#define TUNE_MC_STEPSIZE 0
#define NADJUST 100
#define KBT 0.000949590467


typedef struct{
  int nvar;
  int nx1, nx2;      /* grid size */
  int nhills;
  double cutgauss;
  int tspherical;
  int tsym;           /* expect symmetry info in colvar files */
  int isym;           /* if isym!=-1 read only hills with this isym number */
  int tnewcode;
  double **cvpath;
  double *hllh;
  double *hllw;
  double **scale;
  double *x1pot, *x2pot, **extpot, **extdvxx, **extdvx1, **extdvx2;       /* 2D potential surface on a grid */
} HILLSPATH;

double (*energy)(double *, HILLSPATH *);

char *color[] =
{
"red",
"blue",
"green",
"yellow",
"purple",
"orange",
"pink",
"grey",
"brown",
"gold",
"silver"
};

double dseed;     /* seed for random number generator */
