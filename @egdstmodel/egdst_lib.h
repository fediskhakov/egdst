#include <stdio.h>
#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <float.h>
#include <errno.h>

#define PI 3.1415926535897932384626433832795
#define MAX(X,Y) (((X)>(Y))?(X):(Y))
#define MIN(X,Y) (((X)<(Y))?(X):(Y))
#define SIGN(X) (((X)>0)?(1.0):(-1.0))
//workaround INFINITY - really needed? Works on mac ok
#ifndef INFINITY
#define INFINITY HUGE_VAL
//#define INFINITY (-log(0))
#endif
//numerical tolerance for various routines
#ifndef TOLERANCE
#define TOLERANCE 1e-10
#endif
//artificial point on M-grid to represent zero consumption
#ifndef ZEROCONSUMPTION
#define ZEROCONSUMPTION 0.000001
#endif
//distance to separate doble points on the gird (at intersections)
#ifndef DOUBLEPOINT_DELTA
#define DOUBLEPOINT_DELTA 1e-10
#endif

#include "modelspec.h"

#ifndef EGDSTLIBguard
#define  EGDSTLIBguard

#ifndef DISTRIB
//1 for lognormal, 2 for normal
#define DISTRIB 1   
#endif

//Common globals
extern mxArray *Model;        //model itself
extern int         t0;        //main model dimentions
extern int          T;
extern int     ngridm;
extern int   ngridmax;
extern int   nthrhmax;
extern int         ny;
extern int         nd;
extern int        nnd;
extern int        nst;
extern int       nnst;
extern double    mmax;
extern double      a0;
extern int optim_UasD;
extern int optim_MUnoD;
extern int optim_UnoD;
extern int optim_TRPRnoSH;
extern int optim_stabilizer;
extern double    *stm;       //dimentions of state variables + multiplicators
extern double    *states;    //matrix of state vectors in rows
extern double *decisions;    //matrix of decision vectors in rows
extern int      byval;       //indicator for the use of values instead of indexes of states and decisions in model specific parts (used in sims)
extern char  err[300];       //error description (if it appears)

//egdst_lib.c
void parseModel();
double pdfx(double arg,double mu,double sigma);
double cdfinv(double arg,double mu,double sigma);
double expectation(PeriodVars *curr,PeriodVars *next);
double rescale(PeriodVars *curr,PeriodVars *next,double x);
void eqpgrid(double *out,int n,double mu,double sigma);
int bxsearch (double pn, double *grid,int ngp);
int optimd (double M,double *grid,double *policy,int ngp);
int bxsearch_common (double pn, double *grid,int ngp,int type);
double linter (double x,int ngrid,double *grid, double *fun);
double linter_extrap (double x,int ngrid,double *grid, double *fun, PeriodVars *period);
void printdoublearray (int dim,double *array,const char *prompt);
void printdoubletable (int n,int m,double *array,const char *prompt);
void printintarray (int dim,int *array,const char *prompt);
void tests();
double state(int ist,int colindx);
double decision(int id,int colindx);
double MatlabEval(const char *operation);
double rmax(int n, double *arr);
double cashinhandinverse(PeriodVars curr,PeriodVars next,double arg);
int intmod(int enu,int den);
double Gauss_Legendre_Integration(double a, double b,double par1,double par2, double (*f)(double,double,double));
double cdfni(double p);
void error (const char *desc);

#endif

