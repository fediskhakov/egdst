/* 
  ##########################################################################    
         This is a collection of common routines for egdst model
  ##########################################################################    
*/

#include "egdst_lib.h"

//Common globals (definitions)
mxArray *Model;        //model itself
int         t0;        //main model dimentions
int          T;
int     ngridm;
int   ngridmax;
int   nthrhmax;
int         ny;
int         nd;
int        nnd;
int        nst;
int       nnst;
double    mmax;
double      a0;
int optim_UasD;
int optim_MUnoD;
int optim_UnoD;
int optim_TRPRnoSH;
double    *stm;       //dimentions of state variables + multiplicators
double    *states;    //matrix of state vectors in rows
double *decisions;    //matrix of decision vectors in rows
int      byval;
char  err[300];

// #############################################################################################################
void parseModel()
{   /*fills out the global parameters*/
    int i;
    t0=*(mxGetPr(mxGetProperty(Model,0,"t0")));
    T=*(mxGetPr(mxGetProperty(Model,0,"T")));
    ngridm=*(mxGetPr(mxGetProperty(Model,0,"ngridm")));
    ngridmax=*(mxGetPr(mxGetProperty(Model,0,"ngridmax")));
    nthrhmax=*(mxGetPr(mxGetProperty(Model,0,"nthrhmax")));
    ny=*(mxGetPr(mxGetProperty(Model,0,"ny")));
    nd=*(mxGetPr(mxGetProperty(Model,0,"nd")));
    nnd=*(mxGetPr(mxGetProperty(Model,0,"nnd")));
    nst=*(mxGetPr(mxGetProperty(Model,0,"nst")));
    nnst=*(mxGetPr(mxGetProperty(Model,0,"nnst")));
    mmax=*(mxGetPr(mxGetProperty(Model,0,"mmax")));
    a0=*(mxGetPr(mxGetProperty(Model,0,"a0")));
    stm=(double *) mxGetPr(mxGetProperty(Model,0,"stm"));//dimentions of state variables + multiplicators
    optim_UasD=(int) mxIsLogicalScalarTrue(mxGetField(mxGetProperty(Model,0,"optim"),0,"optim_UasD"));
    optim_MUnoD=(int) mxIsLogicalScalarTrue(mxGetField(mxGetProperty(Model,0,"optim"),0,"optim_MUnoD"));
    optim_UnoD=(int) mxIsLogicalScalarTrue(mxGetField(mxGetProperty(Model,0,"optim"),0,"optim_UnoD"));
    optim_TRPRnoSH=(int) mxIsLogicalScalarTrue(mxGetField(mxGetProperty(Model,0,"optim"),0,"optim_TRPRnoSH")); 
    states=(double *) mxGetPr(mxGetProperty(Model,0,"states"));
    decisions=(double *) mxGetPr(mxGetProperty(Model,0,"decisions"));
    //initialize pointers to continuous variable grids (modelspec.c)
    loadcontinuousgrid();
    //initialize byval to default
    byval=0;
    //initialize err
    *((int*)err)=0;
}

// #############################################################################################################
// for stochastic income
#if DISTRIB == 1
double pdfx(double arg,double mupr,double sigpr)
{   /*Lognormal pdf times x*/
    return (arg/(arg*sigpr*sqrt(2*PI)))*exp(-(log(arg)-mupr)*(log(arg)-mupr)/(2*sigpr*sigpr));
}
double cdfinv(double arg,double mupr,double sigpr)
{   /*Lognormal inverse*/
    return exp(sigpr*cdfni(arg)+mupr);
}
double expectation(PeriodVars *curr,PeriodVars *next)
{
    return exp(mu_param(curr,next)+sigma_param(curr,next)*sigma_param(curr,next)/2);//it!=it+1 by convention for all two period entities
}    
double rescale(PeriodVars *curr,PeriodVars *next,double x)
{
    return exp(mu_param(curr,next)+x*sigma_param(curr,next)); 
}
#endif
#if DISTRIB == 2
double pdfx(double arg,double mupr,double sigpr)
{   /*Normal pdf times x*/
    return (arg/(sigpr*sqrt(2*PI)))*exp(-(arg-mupr)*(arg-mupr)/(2*sigpr*sigpr));
}
double cdfinv(double arg,double mupr,double sigpr)
{   /*Normal inverse*/
    return sigpr*cdfni(arg)+mupr;
}    
double expectation(PeriodVars *curr,PeriodVars *next)
{
    return mu_param(curr,next);//it!=it+1 by convention for all two period entities
}    
double rescale(PeriodVars *curr,PeriodVars *next,double x)
{
    return mu_param(curr,next)+x*sigma_param(curr,next); 
}
#endif

// #############################################################################################################
void eqpgrid(double *out,int n,double mupr,double sigpr)
{   //fills *out array with n equaprobable points
    //relies on pdfx() and cdfinv()
    //IMPORTANT: creates accending grid
    const double ismall=0.00000001;
    int i;
    double x[2];
    x[0]=ismall;
    for (i=1;i<=n;i++)
    {
        x[1]=i/(double)n;
        if (x[1]>1-ismall) x[1]=1-ismall;
        *(out+i-1)=Gauss_Legendre_Integration(cdfinv(x[0],mupr,sigpr),cdfinv(x[1],mupr,sigpr),mupr,sigpr,(double (*)(double,double,double)) pdfx)/(x[1]-x[0]);
//printf("step %d\n[x0 x1] = [%f %f]\ncdfi([x0 x1]) = [%f %f]\nquad=%f\n",i,x[0],x[1],cdfinv(x[0],mupr,sigpr),cdfinv(x[1],mupr,sigpr),*(out+i-1));
        x[0]=x[1];
    }
}

// #############################################################################################################
int bxsearch (double pn, double *grid,int ngp)
{   //finds left side index for linear interpolation
    return bxsearch_common (pn,grid,ngp,0);
}

// #############################################################################################################
int optimd (double M,double *grid,double *policy,int ngp)
{    //finds the optimal discrete policy on the grid of thresholds 
     return policy[bxsearch_common (M,grid,ngp,1)];
}

// #############################################################################################################
int bxsearch_common (double pn, double *grid,int ngp,int type)
{
    /*Return the index of grid popint so that
      pn lies between index and index+1
      pn - point
      grid[] - array of grid points
      ngp - number of grid points
      REQUIRES grid to be increasing or decreasing monotonically
      ASSUMES ngp>=2 (for type=0)
      If point outside of grid returns 0 (first point) or 
      ngp-2 (second last) when typ=0 (intended for linear extrapolation)
      ngp-1 (last) when type=1 (intended for optimal discrete choice lookup)
      Method: binary index search (hence the name)
     */
   int lim1, lim2, lim3; /*limits for the current region of indeces*/  
   if (pn<grid[1]) return 0;//first interval or below grid
   if (type==0 && pn>=grid[ngp-2]) return ngp-2;//last interval or over grid
   else if (type==1 && pn>=grid[ngp-1]) return ngp-1;//over grid
   lim1=1;
   lim2=ngp-2;
   while (lim2-lim1>1)
   {
       lim3=(int) (lim2+lim1)/2;
       if (*grid<=*(grid+ngp-1) && *(grid+lim3)>pn) 
        lim2=lim3;
       else 
        lim1=lim3;
   }
//printf("Answer: %d\n",lim1);   
   return lim1;
}

// #############################################################################################################
double linter (double x,int ngrid,double *grid, double *fun)
{   //linear interpolation
    int i;
    if (ngrid<2) {error("Error: At least two points are required for interpolation!");return -1.0;}
    i=bxsearch(x,grid,ngrid);
    //standard linear interpolation
//printf("x=%f ngrid=%d i=%d grid=[%f,%f] fun=[%f,%f] result=%f\n",x,ngrid,i,grid[i],grid[i+1],fun[i],fun[i+1],fun[i+1]*(x-grid[i])/(grid[i+1]-grid[i]) + fun[i]*(grid[i+1]-x)/(grid[i+1]-grid[i]));    
    return fun[i+1]*(x-grid[i])/(grid[i+1]-grid[i]) + fun[i]*(grid[i+1]-x)/(grid[i+1]-grid[i]);
}

// #############################################################################################################
double linter_extrap (double x,int ngrid,double *grid, double *fun, PeriodVars *prd)
{   //linear interpolation with adjusted extrapolation
    int i;
    double res;
    if (ngrid<2) {error("Error: At least two points are required for interpolation!");return -1.0;}
    i=bxsearch(x,grid,ngrid);
    //return infinities immediately without calculating the weights
    if (!isfinite(fun[i])) return fun[i];
    if (!isfinite(fun[i+1])) return fun[i+1];
//printf("w1=%f w2=%f\n",(tr(prd,x-a0)-tr(prd,grid[i]-a0))/(tr(prd,grid[i+1]-a0)-tr(prd,grid[i]-a0)),(tr(prd,grid[i+1]-a0)-tr(prd,x-a0))/(tr(prd,grid[i+1]-a0)-tr(prd,grid[i]-a0)));
    if (x>a0 && (x>grid[ngrid-1] || x<grid[0])) 
    {   //extrapolation - use extrap adjustment
        res = fun[i+1]*(tr(prd,x-a0)-tr(prd,grid[i]-a0))/(tr(prd,grid[i+1]-a0)-tr(prd,grid[i]-a0)) + fun[i]*(tr(prd,grid[i+1]-a0)-tr(prd,x-a0))/(tr(prd,grid[i+1]-a0)-tr(prd,grid[i]-a0));
//        if (isnan(res)) mexWarnMsgTxt("Extrapolation functions has returned NaN! Check specification of u.extrap");
        if (isnan(res))
        {   mexWarnMsgTxt("Extrapolation functions has returned NaN! Check specification of u.extrap");
printf("x=%f i=%d exptrap(x)=%f  grid0=%f grid1=%f diff=%e ",x,i,tr(prd,x),grid[i],grid[i+1],grid[i+1]-grid[i]);
printf("w1=%f w2=%f ",(tr(prd,x-a0)-tr(prd,grid[i]-a0))/(tr(prd,grid[i+1]-a0)-tr(prd,grid[i]-a0)),(tr(prd,grid[i+1]-a0)-tr(prd,x-a0))/(tr(prd,grid[i+1]-a0)-tr(prd,grid[i]-a0)));
printf("fun1=%f fun2=%f \n",fun[i],fun[i+1]);
        }
        return res;
    }
    else
    {   //standard linear interpolation
//printf("x=%f ngrid=%d i=%d grid=[%f,%f] fun=[%f,%f] result=%f\n",x,ngrid,i,grid[i],grid[i+1],fun[i],fun[i+1],fun[i+1]*(x-grid[i])/(grid[i+1]-grid[i]) + fun[i]*(grid[i+1]-x)/(grid[i+1]-grid[i]));    
        return fun[i+1]*(x-grid[i])/(grid[i+1]-grid[i]) + fun[i]*(grid[i+1]-x)/(grid[i+1]-grid[i]);
    }
}

// #############################################################################################################
void printdoublearray (int dim,double *array,const char *prompt)
{   //prints dim first elements from double array with prompt
    int i;
    printf("%15s ",prompt);
    for (i=0;i<dim;printf("%20.11f ",array[i++]));
    printf("\n");
}

// #############################################################################################################
void printdoubletable (int n,int m,double *array,const char *prompt)
{   //prints n by m matrix of elements from double array with prompt
    int i;
    printf("%15s :\n",prompt);
    for (i=0;i<n*m;i++)
    {
        printf("%12.11f ",array[i]);
        if (fmod((double)i+1,(double)m)<1e-10) printf("\n");
    }
    printf("\n");
}

// #############################################################################################################
void printintarray (int dim,int *array,const char *prompt)
{   //prints dim first elements from double array with prompt
    int i;
    printf("%15s ",prompt);
    for (i=0;i<dim;printf("%16d ",array[i++]));
    printf("\n");
}

// #############################################################################################################
//decode states 
double state(int ist,int colindx)
{   //returns the value of the state variable with through index ist
    //colindex should be between 1 and number of state vars in the state vector
    return states[(colindx-1)*nst+ist];
}

// #############################################################################################################
double decision(int id,int colindx)
{   //returns the value of the decision variable with through index id
    //colindex should be between 1 and number of decision vars in the state vector
    return decisions[(colindx-1)*nd+id];
}

// #############################################################################################################
double MatlabEval(const char *operation)
{
    /*wrapper for no input - single scalar output Matlab calls*/
    mxArray *out;
    mexCallMATLAB (1,&out,0,NULL,operation);
    return (double) mxGetScalar(out);
    /*Matlab clears the memory here by himself*/
}

// #############################################################################################################
double rmax(int n, double *arr)
{   //returns the max of the array
    double maybemax;
    if (n==1) return *arr;
    maybemax = rmax(n-1,arr+1);
    if (*arr>=maybemax) return *arr;
    else return maybemax;
}

// #############################################################################################################
double cashinhandinverse(PeriodVars curr,PeriodVars next,double arg)
{   //finds the value of a(t) so that cashinhand(t+1)=arg
    //ARGUMENTS BY COPY NOT POINTERS!!!
    //stops when tol<ZEROCONSUMPTION/10;
    int cnt=0, maxcnt=100;

//int print;
//if (curr.it==0) print=1;
//else print=0;
//print=1;
    next.savings=arg; //curr and next are local copies!
//if (print) printf("CashInverse: (byval=%d)\ncurr.it=%d curr.ist=%d curr.id=%d\nnext.it=%d next.ist=%d next.id=%d next.sav=%1.25f next.shock=%f\n",byval,curr.it,curr.ist,curr.id,next.it,next.ist,next.id,next.savings,next.shock);
//if (print) printf("arg=%1.25f\n",arg);
    while (fabs(cashinhand(&curr,&next)-arg)>ZEROCONSUMPTION/10) 
    {   //Newton step to find zero of cashinhand-a0
//if (print) printf(" iter %d: cash(x)=%1.10f step=%1.10f diff=%1.15f\n",cnt,cashinhand(&curr,&next),(cashinhand(&curr,&next)-arg)/cashinhand_marginal(&curr,&next),cashinhand(&curr,&next)-arg);
        next.savings-=(cashinhand(&curr,&next)-arg)/cashinhand_marginal(&curr,&next); //save iter savings in right way
//if (print) printf(" new_a=%1.15f\n",next.savings);
        if (++cnt>=maxcnt) {error("Did not manage to invert the intertemporal budget (cashinhand) after performing many-many iterations!");return -1.0;}
    }
    return next.savings;
}

// #############################################################################################################
void error (const char *desc)
{   //fills out the err[300] global
    sprintf(err,"Error:\n%s",desc);
}






////////////////////////////////////////////////////////////////////////////////
//               ROUTINES BORROWED FROM OTHERS
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// File: gauss_legendre_16pts.c                                               //
// Routines:                                                                  //
//    double Gauss_Legendre_Integration_16pts( double a, double b,            //
//                                                     double (*f)(double) )  //
//    void   Gauss_Legendre_Zeros_16pts( double zeros[] )                     //
//    void   Gauss_Legendre_Coefs_16pts( double coef[] )                      //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// The n-th Legendre polynomial is                                            //
//                 Pn(x) = 1/(2^n n!) (d/dx)^n (x^2-1)^n.                     //
// For the n point Gauss-Legendre integral approximation formula the          //
// coefficients are A[i] = 2 (1 - x[i]^2) / (n P(n-1)(x[i])^2 where x[i] is   //
// a zero of the n-th Legendre polynomial Pn(x).                              //
// Note that if x is a zero of Pn(x) then -x is also a zero of Pn(x) and the  //
// coefficients associated with x and -x are equal.                           //
////////////////////////////////////////////////////////////////////////////////
static const double x[] = {
    9.50125098376374401877e-02,    2.81603550779258913231e-01,
    4.58016777657227386350e-01,    6.17876244402643748452e-01,
    7.55404408355003033891e-01,    8.65631202387831743866e-01,
    9.44575023073232576090e-01,    9.89400934991649932601e-01
};
static const double A[] = {
    1.89450610455068496287e-01,    1.82603415044923588872e-01,
    1.69156519395002538183e-01,    1.49595988816576732080e-01,
    1.24628971255533872056e-01,    9.51585116824927848073e-02,
    6.22535239386478928628e-02,    2.71524594117540948514e-02
};
#define NUM_OF_POSITIVE_ZEROS  sizeof(x) / sizeof(double)
#define NUM_OF_ZEROS           NUM_OF_POSITIVE_ZEROS+NUM_OF_POSITIVE_ZEROS
////////////////////////////////////////////////////////////////////////////////
//  double Gauss_Legendre_Integration_16pts( double a, double b,              //
//                                                      double (*f)(double))  //
//                                                                            //
//  Description:                                                              //
//     Approximate the integral of f(x) from a to b using the 16 point Gauss- //
//     Legendre integral approximation formula.                               //
//                                                                            //
//  Arguments:                                                                //
//     double  a   Lower limit of integration.                                //
//     double  b   Upper limit of integration.                                //
//     double *f   Pointer to function of a single variable of type double.   //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f from a to b.                                         //
//                                                                            //
//  Example:                                                                  //
//     {                                                                      //
//        double f(double);                                                   //
//        double integral, lower_limit, upper_limit;                          //
//                                                                            //
//        (determine lower and upper limits of integration)                   //
//        integral = Gauss_Legendre_Integration_16pts(lower_limit,            //
//                                                          upper_limit, f);  //
//        ...                                                                 //
//     }                                                                      //
//     double f(double x) { define f }                                        //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double 
  Gauss_Legendre_Integration(double a, double b,double par1,double par2, double (*f)(double,double,double))
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *px = &x[NUM_OF_POSITIVE_ZEROS - 1];
   const double *pA = &A[NUM_OF_POSITIVE_ZEROS - 1];

   for (; px >= x; pA--, px--) {
      dum = c * *px;
      integral += *pA * ( (*f)(d - dum,par1,par2) + (*f)(d + dum,par1,par2) );
   }

   return c * integral;
}

double 
  ORIGINAL_Gauss_Legendre_Integration_16pts(double a, double b, double (*f)(double))
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *px = &x[NUM_OF_POSITIVE_ZEROS - 1];
   const double *pA = &A[NUM_OF_POSITIVE_ZEROS - 1];

   for (; px >= x; pA--, px--) {
      dum = c * *px;
      integral += *pA * ( (*f)(d - dum) + (*f)(d + dum) );
   }

   return c * integral;
}





/*
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http://www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 */
 
/* Coefficients in rational approximations. */
static const double a[] =
{
	-3.969683028665376e+01,
	 2.209460984245205e+02,
	-2.759285104469687e+02,
	 1.383577518672690e+02,
	-3.066479806614716e+01,
	 2.506628277459239e+00
};

static const double b[] =
{
	-5.447609879822406e+01,
	 1.615858368580409e+02,
	-1.556989798598866e+02,
	 6.680131188771972e+01,
	-1.328068155288572e+01
};

static const double c[] =
{
	-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
	 4.374664141464968e+00,
	 2.938163982698783e+00
};

static const double d[] =
{
	7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00
};

#define LOW 0.02425
#define HIGH 0.97575

double
cdfni(double p)
{
	double q, r;

	errno = 0;

	if (p < 0 || p > 1)
	{
		errno = EDOM;
		return 0.0;
	}
	else if (p == 0)
	{
		errno = ERANGE;
		return -HUGE_VAL /* minus "infinity" */;
	}
	else if (p == 1)
	{
		errno = ERANGE;
		return HUGE_VAL /* "infinity" */;
	}
	else if (p < LOW)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else if (p > HIGH)
	{
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else
	{
		/* Rational approximation for central region */
    		q = p - 0.5;
    		r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}

