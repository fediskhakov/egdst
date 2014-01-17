/*-------------------------------------------------------------------------------------- 
     Solves consumption savings model with an additional discrete choice and state vars
     (see @egdstmodel)
     By: Fedor Iskhakov
     INPUT: egdstmodel object containing the model
            (everything esle is done through h-file composed by matlab code, and hard 
            coded model parts in egdst_model.c)
     OUTPUT: M  = nst by (T-t0+1) 2-dimentional cell array, each cell is
                    double matrix with M, C*, A* and V as rows (number of columns <= ngridmax)
             D  = nst by (T-t0+1) 2-dimentional cell array, each cell is
                    double matrix with D*, TH as rows (number of columns <= nthrhmax)
             where M  - money-at-hand          (1 column)
                   C* - optimal consumption    (2 column)
                   A* - optimal savings        (3 column)
                   V  - value function         (4 column)
                   D* - optimal discrete decision    (1 column)
                   TH - threasholds in optimal descrete decision rules (2 column)

     Known issue: memory leaks when terminates with errors
     
     Convention about variables: t - time, st - state, d - decision
                                 n# - quantity of #
                                 i# - index of #
      in function arguments: it = index of time period t (between 0 and T-t0)
                             ist = index of state in period t
                             ist1 = index of statue in period t+1
                             id = index of decision in period t (relevant for budget at period t+1)
                             id1 = index of decision in period t+1

  --------------------------------------------------------------------------------------*/
#include "egdst_lib.h"

// MANUAL PREPROCESSOR FLAGS
#define EQPRPOINTSoff //use equeprobable representation of random distributions
#define END1off //original terminal period calculation (consumption=money-at-hand including debt)
#define END2 //restricted to A0T end of period consumption
#define MEM cache //memory allocation for envelope (DOES NOT SEEM TO WORK !!!!)

#define DEBUGOUT //additional output for debugging purposes
#define VISUAL_adrawOFF
#define VISUAL_adraw_it 14
#define VISUAL_adraw_ist 19
#define VISUAL_adraw_id 14
#define VISUAL_adraw_ncall 0



#ifdef END2
#define A0T 0.0 //credit constraint in the last period for 2 type of last period calculation
#endif

//Passed from model object preprocessor flags
//mapping from verbosity level to various diagnostics
#if VERBOSE==1
#define VERBOSE_t
#define VERBOSE_isttable
#elif VERBOSE==2
#define VERBOSE_t
#define VERBOSE_isttable
#define VERBOSE_t_results
#elif VERBOSE==3
#define VERBOSE_t
#define VERBOSE_isttable
#define VERBOSE_ist
#define VERBOSE_adraw
#define VISUAL_adraw
#elif VERBOSE==4
#define VERBOSE_t
#define VERBOSE_isttable
#define VERBOSE_ist
#define VERBOSE_adraw
#define VERBOSE_mgridvecs
#elif VERBOSE==5
#define VERBOSE_t
#define VERBOSE_isttable
#define VERBOSE_ist
#define VERBOSE_mgridvecs
#define VERBOSE_adraw
#define VERBOSE_trpr
#elif VERBOSE==6
#define VERBOSE_t
#define VERBOSE_isttable
#define VERBOSE_ist
#define VERBOSE_mgridvecs
#define VERBOSE_envelop
#define VERBOSE_thresh
#elif VERBOSE==7
#define VERBOSE_t
#define VERBOSE_isttable
#define VERBOSE_ist
#define VERBOSE_mgridvecs
#define VERBOSE_adraw
#define VERBOSE_envelop
#define VERBOSE_thresh
#define VERBOSE_trpr
#elif VERBOSE==9
#define VISUAL
#elif VERBOSE==10
#define VISUAL
#define VERBOSE_t
#define VERBOSE_ist
#define VERBOSE_envelop
#define VERBOSE_thresh
#endif


//Structure forward declarations
struct strc_main_entities;
struct strc_results;
struct aspacestruct;
//Static function declarations (in order of use)
static void solver();
static void swappointers(struct strc_main_entities **pn1,struct strc_main_entities **pn2);
static void saveoutput(int it,int ist,struct strc_results *workspace,struct strc_main_entities *currentperiod,struct strc_main_entities *nextperiod);
static void egmbellman(PeriodVars *curr,PeriodVars *next,struct strc_results *workspace,struct strc_main_entities *currentperiod,struct strc_main_entities *nextperiod);
static double valuefunc(PeriodVars *next,int ist1,struct strc_main_entities *nextperiod);
static int envelope2(double *mgridvecs,int nvd,double *evfa0, PeriodVars *curr,void *cache);
static double adraw (struct aspacestruct *aspace);
static void envelop(int it,int ist,int nd,int dim0,double *gridvecs,double *evfa0,double *outgrid,double *outfunc,double *outfunc2,double *outthrh,double *outindx,int *outn,int *outm,void *cache);
static double funcvalue(double x,int j,double***grids,double***funcs,int*curpt,double*evfa0,int it,int ist);
static int comp1(const void *a,const void *b);
static double linter2(double x, double **grid, double **fun);
static void thresholds(int nd,int *curpt,int *dims,double ***grids,double ***funcs,double ***funcs2,int it,int ist,double *evfa0,int pri,int nwi,int mode,double *outgrid,double *outfunc,double *outfunc2,double *outthrh,double *outindx,int *outi,int *outj,int *kmarks);
static void brsolve(double* br0,double* br1,double **grids,double **funcs,PeriodVars *curr,double *evfa0);
static void PeriodVarsOut(PeriodVars *v,const char *prompt);

//Own globals
static mxArray *outCells[2]; //output: pointer to output cell arrays (see above)
double *quadrature;          //column of quadrature weights + column of quadrature abscissas
#ifdef DEBUGOUT
double *dbgout,dbgoutd;
int dbgoutn,dbgouti;
#endif

#if defined(VISUAL) || defined(VISUAL_adraw)
#ifndef VISUAL_ist
#define VISUAL_ist 0
#endif
//things needed for graphing intermediate information in Matlab
char Mstring[500];
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   // The gate to be called from MatLab
    mwSize dims[2];
    int i;
    PeriodVars curr;

    if (nrhs != 1) mexErrMsgTxt("Error: wrong number of inputs!");
#ifdef DEBUGOUT
    if (nlhs != 3) mexErrMsgTxt("Error: wrong number of outputs!");
#else
    if (nlhs != 2) mexErrMsgTxt("Error: wrong number of outputs!");
#endif
    
    // Parce inputs
    Model=(mxArray *) prhs[0];
    parseModel();
    loadparameters();

    // Prepare quadrature points
    quadrature=(double *) mxGetPr(mxGetProperty(Model,0,"quadrature"));
    //recalculate quadrature abscissas into standard normally distributed points
    for (i=0;i<ny;i++) quadrature[ny+i]=cdfni(quadrature[ny+(i)]);
    
    // Create output structure
    //dimentions of the output cell arrays
    dims[0]=nst;
    dims[1]=T-t0+1;
    // 1 cell matrix for m-grid specific values
    plhs[0]=mxCreateCellArray(2,dims);
    outCells[0] = plhs[0];
    // 2 cell matrix for thresholds
    plhs[1]=mxCreateCellArray(2,dims);
    outCells[1] = plhs[1];
#ifdef DEBUGOUT
    //3 matrix for debug data
    dbgoutn = (T-t0+1)*nst*nd* 2*(T-t0+1);
    plhs[2]=mxCreateDoubleMatrix(dbgoutn,7,mxREAL);
    dbgout = mxGetPr(plhs[2]);
    dbgouti = 0;
#endif


#ifdef VERBOSE_t
    //start timer
    MatlabEval("tic");
    printf("\nRunning solver with:\n   verbosity level =%d\n   utility additively separable in d =%d\n   marginal utility independent of d =%d\n   utility independent of d =%d\n   transition probabilities independent of shocks =%d\n",VERBOSE,optim_UasD,optim_MUnoD,optim_UnoD,optim_TRPRnoSH);
#endif

#ifdef EQPRPOINTS                                                
    //if will be using eqprobable points, overwrite quadrature weights
    for(i=0;i<ny;quadrature[i++]=1/(double)ny);
#endif

#ifdef VERBOSE_isttable
//Nice table of numbers of discrete choices by time and state space
printf("\nNumber of discrete choices by it and ist\n");
printf("%4s","ist");
if (nd<10) for (curr.it=0;curr.it<=T-t0;printf("%2d",t0+curr.it++));
else for (curr.it=0;curr.it<=T-t0;printf("%3d",curr.it++));
printf(" <-- age\n");
if (nnst>1 && nd<10)  for(curr.it=0;curr.it<=T-t0;printf("--"),curr.it++);
if (nnst>1 && nd>=10) for(curr.it=0;curr.it<=T-t0;printf("---"),curr.it++);
if (nnst>1 ) printf("----\n");
for(curr.ist=0;curr.ist<nst;curr.ist++)
{
    printf("%4d",curr.ist);
    for(curr.it=0;curr.it<=T-t0;curr.it++)
    {
        i=0;
        for(curr.id=0;curr.id<nd;curr.id++)
        {
            i+=(feasible(&curr) && inchoiceset(&curr));
        }
        if (i>0 && nd<10) printf("%2d",i);
        else if (i>0 && nd>=10) printf("%3d",i);
        else if (i==0 && nd<10) printf("%2s",".");
        else if (i==0 && nd>=10) printf("%3s",".");
    }
    printf("\n");
    if (nnst>1 && (curr.ist+1)%(int)stm[nnd+1]==0 && nd<10)  for(curr.it=0;curr.it<=T-t0;printf("--"),curr.it++);
    if (nnst>1 && (curr.ist+1)%(int)stm[nnd+1]==0 && nd>=10) for(curr.it=0;curr.it<=T-t0;printf("---"),curr.it++);
    if (nnst>1 && (curr.ist+1)%(int)stm[nnd+1]==0) printf("----\n");
}
if (nnst>1)
{
    printf("%4s","ist");
    if (nd<10) for (curr.it=0;curr.it<=T-t0;printf("%2d",t0+curr.it++));
    else for (curr.it=0;curr.it<=T-t0;printf("%3d",curr.it++));
    printf("\n");
}
#endif

    // Call main part of the code (which is below)
    solver();
    if (err[0]) mexWarnMsgTxt(err);

};


// #############################################################################################################
struct strc_main_entities {
    int *egdims;     //length of endogenous grid by ist
    double **egm;    //poineters to endogenous grid points by ist
    double **egc;    //poineters to optimal consumption by ist
    double **egv;    //poineters to value functions by ist
    int *thdims;     //length of discrete policy function by ist
    double **egd;    //poineters to discrete policy function by ist
    double **egth;   //poineters to thresholds by ist
    double *evf;     //expected value at credit constraint evf(a0) by ist
};
struct strc_results {
    double *M, *C, *V, *D, *TH;     //main quantities for given ist
};

// #############################################################################################################
static void solver()
{   //implements the steps of Carroll's method
    int i,it,ist;
    PeriodVars curr,next;//all variables about current and next period
    struct strc_main_entities period[2], *currentperiod, *nextperiod; //current and next period
    struct strc_results workspace[1];  //space for newly calculated entities (max length) !!!should be private in openMP
    
    //initialize
    currentperiod=&period[0];
    nextperiod=&period[1];
    for (i=0;i<2;i++)
    {   //for both current and next periods
        period[i].egdims=(int*) calloc(nst,sizeof(int));    //dimentions of two subsequent periods endogenious grids by ist
        period[i].thdims=(int*) calloc(nst,sizeof(int));
        period[i].egm=(double **) calloc(nst,sizeof(double*));
        period[i].egc=(double **) calloc(nst,sizeof(double*));
        period[i].egv=(double **) calloc(nst,sizeof(double*));
        period[i].thdims=(int*) calloc(nst,sizeof(int));
        period[i].egd=(double **) calloc(nst,sizeof(double*));
        period[i].egth=(double **) calloc(nst,sizeof(double*));
        period[i].evf=(double *) calloc(nst,sizeof(double));
    }
    //for each it-ist
    workspace[0].M=(double*) calloc(ngridmax,sizeof(double));
    workspace[0].C=(double*) calloc(ngridmax,sizeof(double));
    workspace[0].V=(double*) calloc(ngridmax,sizeof(double));
    workspace[0].D =(double*) calloc(nthrhmax,sizeof(double));
    workspace[0].TH=(double*) calloc(nthrhmax,sizeof(double));
    
    //MAIN CYCLE OVER TIME AND SPACE
    for (it=T-t0;it>=0;it--)
    {
#ifdef VERBOSE_ist
printf("\n############ time period %d of %d, age=%d ############\n",it,T-t0,it+t0);
#endif
        /* the next for() statement should be parallelized! */
        for (ist=0;ist<nst;ist++)
        {   //for all values of st(t)
            //put it and ist into the structures
            curr.it=it;
            curr.ist=ist;
            //            
            if (feasible(&curr)==1)
            {   //consider only feasible
#ifdef VERBOSE_ist
printf("\n### ist=%d ###\n",ist);
#endif
                //call parallel subroutine with next period egm,egc,egv,egevf (but pass both periods)
                egmbellman(&curr,&next,workspace,currentperiod,nextperiod);
                //save results for given it-ist
                saveoutput(it,ist,workspace,currentperiod,nextperiod);
                //check if returned because of the error
                if (err[0]) return;
            }//feasibility
        }//ist
        /*swap pointer to the time specific elements*/
        swappointers(&currentperiod,&nextperiod);
#ifdef VERBOSE_t
printf("time period %d done at timer=%f sec.\n",it,MatlabEval("toc"));
#endif

    }//it

    //free (once, within the function just reuse space)
    for (i=0;i<2;i++)
    {   //for both current and next periods
        free(period[i].egdims);
        free(period[i].thdims);
        free(period[i].egm);
        free(period[i].egc);
        free(period[i].egv);
        free(period[i].egd);
        free(period[i].egth);
        free(period[i].evf);
    }
    //for each it-ist
    free(workspace[0].M);
    free(workspace[0].C);
    free(workspace[0].V);
    free(workspace[0].D);
    free(workspace[0].TH);
}

// #############################################################################################################
static void swappointers(struct strc_main_entities **pn1,struct strc_main_entities **pn2)
{   //interchanges two pointers to int
    struct strc_main_entities *tmp = *pn1;
    *pn1=*pn2;
    *pn2=tmp;
}

// #############################################################################################################
struct aspacestruct{    //structure to hold information used by Adraw()
    int ntogenerate;    //number of points to be generated
    int ngenerated;     //number of points generated
    int ncalls;         //number of calls to Adraw()
    double baseM;       //base point (for linear forecasting)
    double baseA;
    double lim1;        //limits for A grid (trasnformed and plain): lower
    double lim2;        //upper    
    double lim2p;
    double lim3;        //focal point
    double lim3p;
    double k3;          //number of points left of focal point
    double lastAguess;  //last generated Aguess
    double M;           //resulting M
    double M1;          //previous resulted M
    int keep;           //whether to keep the calculations from current Aguess
    PeriodVars *curr;   //corresponding current period
};

// #############################################################################################################
static void egmbellman(PeriodVars *curr,PeriodVars *next,struct strc_results *workspace,struct strc_main_entities *currentperiod,struct strc_main_entities *nextperiod)
{   //parallel block of solver algorithm called for each it-ist
    //fill out workspace.M,C,V,D,TH,egdims0 by interpolating over next period next->egm,egc,egv
    int i,j;
    double m1,m2;                   //limits for initial grids at T
    double *evfa0;                  //current period evf(a0) for each id
    double *mgridvecs;              //endogenous id specific grids (m,c,vd,id) vectors
    void *cache;                    //space for envelop: to save on calloc/free
    int nvd,nvd0;                   //id-specific and total number of endogenious grid points
    int skip=0,skip1,skipt=0,interlacing;   //id-specific indeces for vf refinement
    int idcheck;                    //check the non-emptyness of choiceset
    int iy,niy;                     //next period indexes
    double rhs,evf,c1;              //rhs of Euler and expected value function, t+1 period values
    double pr1,pr1pre,checksum;     //transition probability and check
    double *shocks;                 //next period random income representation
    struct aspacestruct aspace;

    //this allocation should be taken outside this function, but this way is easier for openMP
    evfa0=(double *) calloc(nd,sizeof(double));      //evf from A=a0 for current period for each id
    mgridvecs=(double*) calloc(4*nd*ngridmax,sizeof(double));
    cache=(void*)malloc(nd*ngridmax*sizeof(int)+3*nd*sizeof(int)+3*nd*sizeof(double**)+4*nd*ngridmax*sizeof(double*));//the space needed by envelop: to save allocate/free time
    next->it=curr->it+1; //next is about next time period

    /*at this point nextperiod[0].egdims[ist1] contains dimention of the next period
    endogenous grid for given ist1, nextperiod[0].egm[ist1] is the pointer to the
    M-grid, nextperiod[0].egc[ist1] and nextperiod[0].egv[ist1] are respectively pointers to
    optimal consumption and value function on this grid.
    nextperiod[0].thdims[ist1] contains dimention of next period discrete policy for
    given ist1, nextperiod[0].egd[ist1] and nextperiod[0].egth[ist1] contain discrete policy. 
    Also nextperiod[0].evf[ist1] contains expected value functions for exact
    calculation of previous period value function at a0
    ist1=next->ist
     */

#ifdef VISUAL
if (curr->ist==VISUAL_ist) {
 sprintf(Mstring,"point(NaN,'axes','Value functions for it=%d ist=%d');",curr->it,curr->ist);
 mexEvalString(Mstring);
}
#endif

    nvd0=0;//total number of endogenious grid point by id
    idcheck=0;//check that choiceset is non empty
    for (curr->id=0;curr->id<nd;curr->id++)
    {
        if (inchoiceset(curr)==1)
        {   //for all decisions d(t) in choiceset
#ifdef VERBOSE_ist
printf(">id=%3d\n",curr->id);
#endif

#ifdef VISUAL
if (curr->ist==VISUAL_ist) {
 sprintf(Mstring,"point(NaN,'next','id=%d');",curr->id);
 mexEvalString(Mstring);
}
#endif

            idcheck=1;
            nvd=0;//count number of points in id-specific endogenous grid
            evfa0[curr->id]=0.0;//initialize evf to distinguish from -inf
            
            //###### period T separately
            if (curr->it==T-t0)
            {   //period T
                // UNDER ASSUMPTIONS OF u'>0 everywhere !!!!!!!!!!!!!!!!!!!!!
                // 1. consumption = money-in-hand
                // 2. credit constraint holds everywhere!
#ifdef END1
                //Same credit constraint at the last period = debt if a0<0 : consumption=money-at-hand - a0
                for (i=0;i<2;i++) //set two points: for interpolation at next time iteration
                {
                    mgridvecs[4*nvd0]=mmax*(1.1+.1*i); //mgrid point ARBITRARY SCALORS
                    mgridvecs[4*nvd0+1]=mgridvecs[4*nvd0]-a0;//c point
                    mgridvecs[4*nvd0+2]=utility(curr,mgridvecs[4*nvd0+1]);//vd point
                    mgridvecs[4*nvd0+3]=(double) curr->id; //remember id 
                    evfa0[curr->id]=0; //evf(a0) no future value in last time period
                    //update indeces
                    nvd++;//number saved points for this id
                    nvd0++;//total number of A draws
                }
#endif
#ifdef END2
                //Special credit constraint at last period A0T
                //No debt to be left after last period : consumption=money-at-hand - max(a0,A0T)
                //USE LOG TRANSFORM
                m1=tr(curr,ZEROCONSUMPTION-A0T);
                m2=tr(curr,mmax-A0T);
                evfa0[curr->id]=-INFINITY; //evf(a0) no future value in last time period
                for (i=0;i<ngridm;i++) //initial grid over M
                {
                    mgridvecs[4*nvd0]=trinv(curr, m1+i*(m2-m1)/(ngridm-1) )+A0T; //mgrid point
                    mgridvecs[4*nvd0+1]=mgridvecs[4*nvd0]-A0T;//c point
                    mgridvecs[4*nvd0+2]=utility(curr,mgridvecs[4*nvd0+1]);//vd point
                    mgridvecs[4*nvd0+3]=(double) curr->id; //remember id 
#ifdef VISUAL
if (curr->ist==VISUAL_ist) {
 sprintf(Mstring,"point([%1.5f %1.5f]);",mgridvecs[4*nvd0],mgridvecs[4*nvd0+2]);
 mexEvalString(Mstring);
}
#endif
                    //update indeces
                    nvd++;//number saved points for this id
                    nvd0++;//total number of A draws
                }
#endif
            }
            else
            {   //all other periods
                // First, points for the endogenious grid construction come from "a-cycle"
                //Initialize aspace
                aspace.ncalls=0;
                aspace.ntogenerate=ngridm;
                aspace.ngenerated=0;
                aspace.M=INFINITY;
                aspace.curr=curr;//pointer only
#ifdef VERBOSE_adraw
printf("\nit=%d ist=%d id=%d\n",curr->it,curr->ist,curr->id);
#endif
                //Loop over A-guesses
                while (next->savings=adraw(&aspace),next->savings!=-INFINITY)
                {   //for each guess of A
                    if (err[0]) return;//error in adraw
                    //get ready to calculate rhs of Euler equation and expected value function
                    rhs=0;
                    if (aspace.keep==1) evf=0; //vf is only calculated when keeping the results
                    checksum=0.0; 
                    for (next->ist=0;next->ist<nst;next->ist++)
                    {   //iterate over state values at period t+1 with positive transition probability
                        if (feasible(next)==1) //it+1 because feasible is only defined within one period
                        {   //skip infeasible
                            if (optim_TRPRnoSH==1)
                            {
                                pr1pre=trpr(curr,next,1); //calculate pr1 for undefined shock because it is independent
#ifdef VERBOSE_trpr
if (pr1pre!=0.0) printf(" * trpr( it=%d ist=%d aguess=%f -> ist1=%d |id=%d,indep of shock) = %f\n",curr->it,curr->ist,next->savings,next->ist,curr->id,pr1pre);
#endif                                                
                                if (pr1pre==0.0) continue;
                            }
                            //check for degenerate distribution
                            if (sigma_param(curr,next)<=0 || ny==1) niy=1; //it!=it+1 by convention for all two period entities
                            else
                            {   //non-degenerate distribution
                                 niy=ny;
                                //populate the vector of shocks to represent income distribution at t+1
                                shocks=quadrature+ny; //by default use quadrature points
#ifdef EQPRPOINTS                                                
                                //overwrite with equapropable points if flagged
                                eqpgrid(shocks,ny,mu_param(curr,next),sigma_param(curr,next));
#endif
                            }
                            //calculate rhs of Euler equation and expected value function
                            for (iy=0;iy<niy;iy++)
                            {
                                //recalculate standard normal to needed distribution
                                if (niy==1)
                                {   //degenerate
                                    next->shock=expectation(curr,next); //use expectation in degenerate case!
                                    if (optim_TRPRnoSH!=1) pr1=trpr(curr,next,1);//it!=it+1 by convention for all two period entities
                                    else pr1=pr1pre;
#ifdef VERBOSE_trpr
if (pr1!=0.0 && optim_TRPRnoSH!=1) printf(" * trpr( it=%d ist=%d aguess=%f shock=%f -> ist1=%d |id=%d, degenerate) = %f\n",curr->it,curr->ist,next->savings,next->shock,next->ist,curr->id,pr1);
#endif                                                
                                }
                                else
                                {
#ifdef EQPRPOINTS                                                
                                    next->shock=shocks[iy];
#else
                                    next->shock=rescale(curr,next,shocks[iy]); //rescale for right parameters
#endif
                                    if (optim_TRPRnoSH!=1) pr1=trpr(curr,next,1);//it!=it+1 by convention for all two period entities
                                    else pr1=pr1pre;
#ifdef VERBOSE_trpr
if (pr1!=0.0 && optim_TRPRnoSH!=1) printf(" * trpr( it=%d ist=%d aguess=%f shock=%f -> ist1=%d |id=%d) = %f quadrature= %f\n",curr->it,curr->ist,next->savings,next->shock,next->ist,curr->id,pr1,quadrature[iy]);
#endif                                                
                                    pr1*=quadrature[iy];
                                }
                                if (pr1==0.0) continue;
                                checksum+=pr1;
                                next->cash=cashinhand(curr,next);
                                // ########## LINEAR INTERPOLATION OVER OPTIMAL CONSUMPTION NEXT PERIOD
                                c1=linter(next->cash,nextperiod[0].egdims[next->ist]+1,nextperiod[0].egm[next->ist],nextperiod[0].egc[next->ist]);//+1 for additional point for a0!
                                //if last segment in consumption function is downward, extrapolation may give negative consumption..
                                if (next->cash > *(nextperiod[0].egm[next->ist]+nextperiod[0].egdims[next->ist])) c1=MAX(c1,*(nextperiod[0].egc[next->ist]+nextperiod[0].egdims[next->ist]));//constant extrapolation
                                if (err[0]) return;
                                if (c1<=0) break; /*No point to calculate further - evf=inf for u(c)=-inf at c=0*/
                                //check if optimal discrete policy is needed
                                if (optim_MUnoD!=1 || (optim_UnoD!=1 && aspace.keep==1 && next->cash<*(nextperiod[0].egm[next->ist]+1)))
                                {   //cases: marginal utility depends on id, need id1 to calculate rhs
                                    //       utility depends on id, and need to calculate exact next period value
                                        next->id=optimd(next->cash,nextperiod[0].egth[next->ist],nextperiod[0].egd[next->ist],nextperiod[0].thdims[next->ist]);
                                }
                                else next->id=0; //fast case: pick any value
                                rhs+=pr1*utility_marginal(next,c1)*cashinhand_marginal(curr,next);
                                if (aspace.keep==1)
                                {   //expected value function
                                    evf+=pr1*valuefunc(next,next->ist,nextperiod);
                                    if (err[0]) return; 
                                    if (evf==-INFINITY) break;
                                }//aspace.keep
                            }//iy
                            if (c1<=0 || evf==-INFINITY) break; //break out of ist1 cycle too
                        }//feasible ist1
                    }//ist1
                    
                    //check returned values
                    if (c1>0 && evf>-INFINITY && fabs(checksum-1)>TOLERANCE) 
                    {   //check for the trpr to sum up
                        printf("it=%d ist=%d id=%d checksum-1=%1.5e\n",curr->it,curr->ist,curr->id,checksum-1);
                        error("Transition probabilities don't sum up! Check model specification!"); 
                        return;
                    }
                    else if (c1<=0 || evf==-INFINITY)
                    {   //emergiency break: negative consumption or -infinite value
                        if (aspace.ngenerated==0)
                        {   //Adraw is in the initial stage of searching for base point: too bad
//if (evf==-INFINITY) printf("evf==-INFINITY ");
                            printf("it=%d ist=%d id=%d ist1=%d Aguess=%f cash1=%1.3f c1=%1.3f\n",curr->it,curr->ist,curr->id,next->ist,next->savings,next->cash,c1);
                            error("Failed to find any value of savings to result in positive consumption next period! Increase mmax.");
                            return;
                        }
                        //1: save the evf
                        evfa0[curr->id]=-INFINITY;
                        //2: return M to Adraw (assuming c=0 <==> Inada)
                        aspace.M=next->cash;//mgrid point
                        if (c1<=0)
                        {   //+ZEROCONSUMPTION case
                            //1: signal to adraw()
                            aspace.M=a0-1;
                            //2: calculate the point so that next period cash in hand is above first VF point given the shock
                            //   where c1<=0 was encountered (which is shock1)
#ifdef VERBOSE_adraw
                            if (nextperiod[0].evf[next->ist] > -INFINITY)
                            {
                                aspace.lastAguess=cashinhandinverse(*curr,*next,a0)+ZEROCONSUMPTION; //ensure over a0 next time period
                                printf("Calling cashinhandinverse() to find A such that next period M=%f (case 1)\n",a0);
                            }
                            else
                            {
                                aspace.lastAguess=cashinhandinverse(*curr,*next,*(nextperiod[0].egm[next->ist]+1))+ZEROCONSUMPTION; //ensure positive C next time period
                                printf("Calling cashinhandinverse() to find A such that next period M=%f (case 2)\n",*(nextperiod[0].egm[next->ist]+1));
                            }
#else
                            if (nextperiod[0].evf[next->ist] > -INFINITY)
                            aspace.lastAguess=cashinhandinverse(*curr,*next,a0)+ZEROCONSUMPTION; //ensure over a0 next time period
                            else
                            aspace.lastAguess=cashinhandinverse(*curr,*next,*(nextperiod[0].egm[next->ist]+1))+ZEROCONSUMPTION; //ensure positive C next time period
#endif

                            if (err[0]) return;
                        }
#ifdef VERBOSE_adraw
                        else printf("Result from this Aguess is not to be saved because of evf=-inf\n");
#endif
                        //if evf=-inf just continue without saving the point to endogenous grid
                        //could do similarly to c1<=0 but since this is untypical, let grid shrink and approximation errors stay
                    }
                    else
                    {   //normal case
                        rhs*=discount(curr);
                        //pass the mgrid value to adraw() function for next Aguess calculation
                        aspace.M=next->savings+utility_marginal_inverse(curr,rhs);//mgrid point
//                        if (aspace.keep==1) //skip points corresponding to vf=-inf from resulting grid
                        if (aspace.keep==1 && isfinite(aspace.M)) //skip points corresponding to vf=-inf from resulting grid
                        {   
#ifdef VERBOSE_adraw
                            if (fabs(aspace.lastAguess-a0)<TOLERANCE && evfa0[curr->id]>-INFINITY) 
                            {
                                evfa0[curr->id]=evf;//there was no zero consumption adjustment, save evfa0
                                printf("Saved evfa0=%1.5f\n",evfa0[curr->id]);
                            }
#else
                            if (fabs(aspace.lastAguess-a0)<TOLERANCE && evfa0[curr->id]>-INFINITY) evfa0[curr->id]=evf;//there was no zero consumption adjustment, save evfa0
#endif
                            //keep the results of calculation
                            mgridvecs[4*nvd0]=aspace.M; //mgrid point first
                            mgridvecs[4*nvd0+1]=aspace.M-next->savings;//c point second
//PeriodVarsOut(curr,"within adraw");                            
                            mgridvecs[4*nvd0+2]=utility(curr,mgridvecs[4*nvd0+1])+discount(curr)*evf;//vd point third
                            mgridvecs[4*nvd0+3]=(double) curr->id; //remember id 

#ifdef VISUAL
if (curr->ist==VISUAL_ist) {
 sprintf(Mstring,"point([%1.5f %1.5f]);",mgridvecs[4*nvd0],mgridvecs[4*nvd0+2]);
 mexEvalString(Mstring);
}
#endif

                            //update indeces
                            nvd++;//number saved points for this id
                            nvd0++;//total number of A draws
                            if (nvd>=ngridmax) {error("Not enough space for endogenous grid. Increase max number of grid points for M!");return;}//not enough grid points allowed in m
                        }
                    }
                }//next Adraw

                //secondary envelope
                skip=envelope2(&mgridvecs[4*(nvd0-nvd)],nvd,&evfa0[curr->id],curr,cache);
                if (err[0]) {printf("\n\nError during period it=%d, age=%d, st=%d\n",curr->it,curr->it+t0,curr->ist);return;}
                //correct point counts
                nvd-=skip;//number of this id points
                nvd0-=skip;//total number of points
                
#ifdef VERBOSE_mgridvecs
printf("\nit=%d ist(t)=%d >>> id(t)=%d\n",curr->it,curr->ist,curr->id);
printf("mgridvecs: [M C VF id]\n");
for (i=nvd0-nvd;i<nvd0;printf("%d: %8.3f %8.3f %8.3f %8d\n",i,mgridvecs[i*4],mgridvecs[i*4+1],mgridvecs[i*4+2],(int)mgridvecs[i*4+3]),i++);
printf("Secondary upper envelope skipped %d points\n",skipt);
#endif
            }//it!=T-t0
        }// id choiceset
#ifdef VERBOSE_mgridvecs
        else
        {
printf("\nit=%d ist(t)=%d >>> id(t)=%d\nNot in choiceset\n",curr->it,curr->ist,curr->id);
        }
#endif              
    }//id

    //NOW START CALCULATION OF UPPER ENVELOPE
#ifdef VERBOSE_ist
printf(">envelope\n");
#endif
    if (idcheck==0)
    {   //empty choiceset
        printf("\nit=%d age=%d ist(t)=%d\n",curr->it,curr->it+t0,curr->ist);
        printf("state(t)=[ ");
        for (i=0;i<nnst;printf("%6.3f ",states[curr->ist+(i++)*nst]));
        printf("]\n");
        error("Empty choiceset encountered! Check model specifications!");
        return;
    }
    //special rare case when all choices produced empty grids (vf=-inf)
    if (nvd0==0)
    {   
        printf("\n\nError during period it=%d, age=%d, st=%d\n",curr->it,curr->it+t0,curr->ist);
        error("All of the choices lead to -inf value functions for all values of money-at-hand!\nThe following may help:\n > increase number of grid points\n > decrease number of time periods\n > increase discount factor\n > decrease the number of points in discrete representation of shocks\n");
        return;
    }
    //calculate upper envelope
#ifdef VISUAL
if (curr->ist==VISUAL_ist) {
 sprintf(Mstring,"point(NaN,'next','UpperEnvelope');");
 mexEvalString(Mstring);
}
#endif
#ifdef DEBUGOUT
dbgoutd=-1;    
#endif
    envelop(curr->it,curr->ist,nd,nvd0,mgridvecs,evfa0,workspace[0].M,workspace[0].V,workspace[0].C,workspace[0].TH,workspace[0].D,currentperiod[0].egdims+curr->ist,currentperiod[0].thdims+curr->ist,cache); //save dimentions to current
    if (err[0]) {printf("\n\nError during period it=%d, age=%d, st=%d\n",curr->it,curr->it+t0,curr->ist);return;}
    //check the number of points in unified grid: if d-specific value functions don't overlap we are in trouble.. 
    if (currentperiod[0].egdims[curr->ist]==0 || currentperiod[0].thdims+curr->ist==0)
    {
        printf("Upper envelope calculation failed!\n it=%d ist=%d\n",curr->it,curr->ist);
        error("Failed to compute upper envelope, most likely individual grids don't overlap!");
        return;
    }
    //save evf for exact vf calculations during next iteration
    currentperiod[0].evf[curr->ist]=evfa0[(int)workspace[0].D[0]];//evf from this period of the first optimal d
#ifdef VERBOSE_t_results 
printf("\nResulting value function and optimal consumption for it=%d ist(t)=%d\n",curr->it,curr->ist);
printdoublearray(currentperiod[0].thdims[curr->ist],workspace[0].D,"Optimal D:");
printdoublearray(currentperiod[0].thdims[curr->ist],workspace[0].TH,"Thresholds:");
printdoublearray(1,&currentperiod[0].evf[curr->ist],"Evf(m0):");
printf("%3s %10s %10s %10s\n","#","M","C","VF");
for (i=0;i<37;printf("-"),i++);
printf("\n");
printf("%3d %10.5f %10.5f %10s\n",0,a0,0.0,"NaN");
for (i=0;i<currentperiod[0].egdims[curr->ist];i++)
{
    printf("%3d %10.5f %10.5f %10.5f\n",i+1,workspace[0].M[i],workspace[0].C[i],workspace[0].V[i]);
}
for (i=0;i<37;printf("-"),i++);
printf("\n\n");
#endif

    //clean up local arrays
    free(evfa0);
    free(mgridvecs);
    free(cache);
}

// #############################################################################################################
static double valuefunc(PeriodVars *next,int ist1,struct strc_main_entities *nextperiod)
{   // Value function for next period cash
    int ist_save = next->ist;
    double res;
    next->ist = ist1;
    if (next->cash<*(nextperiod[0].egm[next->ist]+1) && nextperiod[0].evf[next->ist]>-INFINITY) //check against first point of next period grid
    {   //credit constraint binds, so calculate percise VF
        //if evf1a0[curr->id*nst+next->ist]==-INFINITY, utility which should be of ZEROCONSUMPTION is not important
        res=(utility(next,next->cash-a0)+discount(next)*nextperiod[0].evf[next->ist]);
    }
    else
    {   //credit constraint does not bind
        // ########## LINEAR INTERPOLATION OVER VALUE FUNCTION
        res=linter_extrap(next->cash,nextperiod[0].egdims[next->ist],nextperiod[0].egm[next->ist]+1,nextperiod[0].egv[next->ist]+1,next);//+1 to skip VF at a0
    }
    next->ist=ist_save;
    return res;
}


// #############################################################################################################
static int envelope2(double *mgridvecs,int nvd,double *evfa0, PeriodVars *curr,void *cache)
{   // Upper envelope for the overlapping sections of a value function (that is due to non-concave nature of the problem)
    // Work with points for this id: blocks of 4 points in mgridvecs[] between nvd0-nvd and nvd0
    //ASSUMING the Adraw is generating a-guesses in INCREASING ORDER !!!!
    //RELIES on monotonicity of savings function
    // nvd = number of points for this id

    //INPUT: mgridvec - section of mgridvec corresponding to current id
            //4*i+0 - m grid point
            //4*i+1 - consumption
            //4*i+2 - value function
            //4*i+3 - index of the discrete decision
    //       nvd - number of points for current id (generated by EGM step)
    //       evfa0[0] - evf for the analytical part of the value function
    //       curr->X - period variables (it, ist, id)
    //       cache - cache for envelop (conditional on MEM)
    
    //This function: prepares inputs to envelop(), calls it and puts the results back into mgridvecs
    //RETURN: number of points skipped in secondary envelope 

    //TODO: work on unification of envelop() function

    int i;
    double *outspace, *evf;
    int ind0, ind0i=0, dim0=0;
    int n=nvd,m;

#ifdef DEBUGOUT
dbgoutd=curr->id;    
#endif

    outspace=calloc(5*ngridmax,sizeof(double)); //space for envelop() output
    evf=calloc(10000,sizeof(double)); //space for envelop() output
//printf("\nenvelop2(): it=%d ist=%d id=%d\n",curr->it,curr->ist,curr->id);
    //split the gridvecs
    //first point
    ind0=curr->id; //start with real id index to allow for correct analytical form on the fist segment, all other indexes are never used
    ind0i=0;
    mgridvecs[3]=ind0;
    evf[ind0]=evfa0[0];
    dim0=1;
    for (i=1;i<nvd;i++)
    {
        if (mgridvecs[4*(i-1)]>mgridvecs[4*i] || mgridvecs[4*(i-1)+2]>mgridvecs[4*i+2]) //CRITERIA for split
        {   
//printf("Criteria true: %1.4f>%1.4f OR %1.4f>%1.4f\n",mgridvecs[4*(i-1)],mgridvecs[4*i],mgridvecs[4*(i-1)+2],mgridvecs[4*i+2]);
            //add point for CONSTANT extrapolation
            if (dim0>=ngridmax) {error("Not enough space for endogenous grid in envelop2()");return -1;}//not enough grid points allowed in m
            mgridvecs[4*(nvd+ind0i)+0]=1.5*mmax;
            mgridvecs[4*(nvd+ind0i)+1]=mgridvecs[4*(i-1)+1]; //?????
            mgridvecs[4*(nvd+ind0i)+2]=mgridvecs[4*(i-1)+2]; //constant extrapolation of value funciton
            mgridvecs[4*(nvd+ind0i)+3]=ind0;
            dim0++;//count added point
            //next index
            ind0++;
            ind0i++;
            if (ind0>=10000) {error("10000 is not enough in envelop2()");return -1;}
            //next evf
            evf[ind0]=-INFINITY;
        }
        mgridvecs[4*i+3]=ind0;
        dim0++;
    }
//printf("ind0i+1=%d dim0=%d\n",ind0i+1,dim0);

    if (ind0i>0)
    {   
/*
printf("\nit=%d ist(t)=%d >>> id(t)=%d\n",curr->it,curr->ist,curr->id);
printf("NEW mgridvecs [M C VF id] for envelop2\n");
for (i=0;i<dim0;printf("%d: %8.3f %8.3f %8.3f %8d\n",i,mgridvecs[i*4],mgridvecs[i*4+1],mgridvecs[i*4+2],(int)mgridvecs[i*4+3]),i++);
printdoublearray(ind0+1,evf,"evf:");
*/
#ifdef VISUAL
if (curr->ist==VISUAL_ist) {
 sprintf(Mstring,"point(NaN,'next','id=%d env2');",curr->id);
 mexEvalString(Mstring);
}
#endif

        //call envelop()
        envelop(curr->it,curr->ist,ind0+1,dim0,mgridvecs,evf,
                outspace+ngridmax*0, //m
                outspace+ngridmax*1, //vf
                outspace+ngridmax*2, //c
                outspace+ngridmax*3, //thrh
                outspace+ngridmax*4, //indx
                &n, //final number of point in the envelope
                &m, //nbr of th
                cache);
    /*    
        int it,                 //current period (for utility calculation)
        int ist,                //current state (for utility calculation)
        int nd,                 //number of functions to compare
        int dim0,               //total number of points in individual id specific grids
        double *gridvecs,       //input array with quadruples: (grid point,func value,func2 value,func index)
        double *evfa0,          //expected value corresponding to a0 by id
        double *outgrid,        //pointer to space for resulting grid
        double *outfunc,        //pointer to space for resulting func values
        double *outfunc2,       //pointer to space for resulting func2 values
        double *outthrh,        //pointer to space for thresholds
        double *outindx,        //pointer to space for max-attaining func indeces
        int *outn,              //pointer to int to hold the resulting number of points in outgrid
        int *outm,              //pointer to int to hold the resulting number of points in outthrh and outindex
        void *cache)            //preallocated space to save time on allocations
    */
//    printf("n=%d m=%d\n",n,m);

        if (n>=ngridmax) {error("Not enough space for endogenous grid in envelop2()");return -1;}//not enough grid points allowed in m

        //collect output from envelope into original mgridvecs
        for (i=0;i<n;i++)
        {
            mgridvecs[4*i+0]=outspace[ngridmax*0+i];
            mgridvecs[4*i+1]=outspace[ngridmax*2+i];
            mgridvecs[4*i+2]=outspace[ngridmax*1+i];
            mgridvecs[4*i+3]=curr->id;
        }

    }
    else
    {
        //rename back the id index
        for (i=0;i<nvd;i++)
        {
            mgridvecs[4*i+3]=curr->id;
        }
        //same number of points
        n=nvd;
    }
    
    //cleanup
    free(outspace);
    free(evf);
    
    //return the number of skipped points
    return nvd-n;
}


// #############################################################################################################
static void saveoutput(int it,int ist,struct strc_results *workspace,struct strc_main_entities *currentperiod,struct strc_main_entities *nextperiod)
{   //saves the main quantities to the output
    //and mark their locations in the memory for next period calculations
    double *pnt;
    int i,cellindex,mlength,thlength;
    cellindex=ist+it*nst;
    //First M C A V
    //Amend optimal consumption rule with the first point!
    mlength=currentperiod[0].egdims[ist]+1;
    mxSetCell(outCells[0],cellindex,mxCreateDoubleMatrix(mlength,4,mxREAL));
    pnt=mxGetData(mxGetCell(outCells[0],cellindex));
    //M
    pnt[0]=a0;
    memcpy(pnt+1,workspace[0].M,currentperiod[0].egdims[ist]*sizeof(double)); //first column is M
    currentperiod[0].egm[ist]=pnt; //remember that length is +1 !!!!!
    //C
    pnt[mlength]=0;
    memcpy(pnt+mlength+1,workspace[0].C,currentperiod[0].egdims[ist]*sizeof(double)); //second column is C
    currentperiod[0].egc[ist]=pnt+mlength; //remember that length is +1 !!!!!
    //A
    for (i=0;i<mlength;pnt[2*mlength+i]=pnt[i]-pnt[mlength+i],i++);
    //VF
    pnt[3*mlength]=currentperiod[0].evf[ist]; //skip calculation of VF at a0, but save evf!
    memcpy(pnt+3*mlength+1,workspace[0].V,currentperiod[0].egdims[ist]*sizeof(double)); //fourth column is V
    currentperiod[0].egv[ist]=pnt+3*mlength; //remember that length is +1 !!!!!
    //Second D TH
    thlength=currentperiod[0].thdims[ist];
    mxSetCell(outCells[1],cellindex,mxCreateDoubleMatrix(thlength,2,mxREAL));
    pnt=mxGetData(mxGetCell(outCells[1],cellindex));
    //D
    memcpy(pnt,workspace[0].D,thlength*sizeof(double)); //first column is D
    currentperiod[0].egd[ist]=pnt;
    //TH
    memcpy(pnt+thlength,workspace[0].TH,thlength*sizeof(double)); //second column is TH
    currentperiod[0].egth[ist]=pnt+thlength;
}

// #############################################################################################################
static double adraw (struct aspacestruct *aspace)
{   /*returns the next value in the stream of guesses of A
      terminates the sequence by sending -INFINITY
    */    
    double aa,bb,mstep,astep,upper;
    int i;
    
    aspace[0].ncalls+=1.0;
    if (aspace[0].ncalls>=ngridmax)
    {
        printf("Aspace dump:\n it=%d\n ist=%d\n id=%d\n ntogenerate=%d\n ngenerated=%d\n ncalls=%d\n baseM=%f\n baseA=%f\n lastAguess=%f\n M=%f\n keep=%d\n",
        aspace->curr->it,aspace->curr->ist,aspace->curr->id,
        aspace->ntogenerate,
        aspace->ngenerated,
        aspace->ncalls,
        aspace->baseM,
        aspace->baseA,
        aspace->lastAguess,
        aspace->M,
        aspace->keep);
        mexWarnMsgTxt("Emergiency exit from adraw, possibly infinite loop! Check model specifications!");
        aspace[0].lastAguess=-INFINITY; //stop
        return aspace[0].lastAguess;
    }
    if (aspace[0].ngenerated==0)
    {   //initial search of Aguess such that M(Aguess)<=mmax
        aspace[0].keep=0;//don't keep calculations signal
        if (aspace[0].M==INFINITY)
        {   //first step
            aspace[0].lastAguess=mmax;
#ifdef VERBOSE_adraw
printf("adraw init search: first step, new Aguess = %f\n",aspace[0].lastAguess);        
#endif
#ifdef VISUAL_adraw
if (aspace[0].curr->it==VISUAL_adraw_it && aspace[0].curr->ist==VISUAL_adraw_ist && aspace[0].curr->id==VISUAL_adraw_id) {
 sprintf(Mstring,"point(NaN,'axes','Adraw() algorithm visualization for it=%d ist=%d id=%d');",aspace[0].curr->it,aspace[0].curr->ist,aspace[0].curr->id);
 mexEvalString(Mstring);
} 
#endif
        }
        else
        {
            if (aspace[0].M<=mmax)
            {   //found base point!
                aspace[0].baseA=aspace[0].lastAguess; //abscissa (x)
                aspace[0].baseM=aspace[0].M; //ordinate (y)
                aspace[0].ngenerated=1.0; //mark change of stage
                aspace[0].keep=1; //keep calculations signal
                aspace[0].lastAguess=a0; //return a0
                aspace[0].k3=0;//init
                aspace[0].M1=a0-1;//init
#ifdef VISUAL_adraw
if (aspace[0].curr->it==VISUAL_adraw_it && aspace[0].curr->ist==VISUAL_adraw_ist && aspace[0].curr->id==VISUAL_adraw_id) {
 sprintf(Mstring,"point([%1.5f %1.5f],'next','Base point');",aspace[0].baseA,aspace[0].baseM);
 mexEvalString(Mstring);
}
#endif
            }
            else
            {   //keep looking: binary search towards a0
                if (aspace[0].lastAguess-a0<TOLERANCE)
                {
                    printf("\nLast exploratory call from adraw: A=%5.8f --> M=%5.8f, |A-a0|=%5.8e\n",aspace[0].lastAguess,aspace[0].M,fabs(aspace[0].lastAguess-a0));
                    error("Could not complete initial stage in adraw()..\nSeems like M(a0)>mmax! Increase mmax!"); 
                    return -1.0;
                }
                aspace[0].lastAguess=(aspace[0].lastAguess+a0)/2;
            }       
#ifdef VERBOSE_adraw
printf("adraw init search: returned M = %f, new Aguess = %f\n",aspace[0].M,aspace[0].lastAguess);        
#endif
        }
    }    
    else
    {   //generation of next Aguess with successive linear approximations of the point where M(Aguess)==mmax
        //from now on aspace[0].ngenerated is the number of Aguess generated

        aa=(aspace[0].M-aspace[0].baseM)/(aspace[0].lastAguess-aspace[0].baseA);
        //linear model to forcast the upper bound of A guesses to be generated
        //TEMPORARY switch off rescaling and forcasting
        if (aspace[0].ngenerated==1)
        //if (aa>0 && aspace[0].M>aspace[0].M1 && aspace[0].ngenerated<aspace[0].ntogenerate-10) //last points separatelly
        {   //update upper ONLY if new point is reasonable
            bb=aspace[0].baseM-aa*aspace[0].baseA;
            upper=MIN(mmax,(mmax-bb)/aa);//updated upper bound estimate
        }
        else
        {   //or computer aa,bb from limit saved on the initialization stage
            upper=aspace[0].lim2p;
            aa=(mmax-aspace[0].M)/(upper-aspace[0].lastAguess);
            bb=mmax-aa*upper;
        }    
        aspace[0].M1=aspace[0].M;//remember returned M

        //a0 case: initializations
        //TODO: re-design how ZEROCONSUMPTION case works here and in the main body of the program
        if (aspace[0].ngenerated==1 && aspace[0].k3==0) //k3!=0 means that limits are already defined in ZEROCONSUMPTION case
        {   //initializevariables after first 'keep' call
            aspace[0].ntogenerate=(double)ngridm; //standard number of point for the above credit constraint area
            //prepare for symmetric-log grid
            aspace[0].lim2p=MIN(mmax,(mmax-bb)/aa); //upper A (plain)
            aspace[0].lim3p=-bb/aa; //focal point (plain)
            if (a0<0 && a0<aspace[0].lim3p) //sym-log grid with focal point
               aspace[0].k3=MAX( floor(aspace[0].ntogenerate*(aspace[0].lim3p-a0)/(aspace[0].lim2p-a0)), 2.0);
            else
            { //usual log grid, redefine lim3
                aspace[0].lim3p=a0;
                aspace[0].k3=1.0;
            }
            //define transformed lim
            aspace[0].lim1=tr(aspace[0].curr,aspace[0].lim3p-a0);
            aspace[0].lim2=tr(aspace[0].curr,aspace[0].lim2p-aspace[0].lim3p);
            aspace[0].lim3=tr(aspace[0].curr,0);

#ifdef VERBOSE_adraw
if ((int)aspace[0].k3==1)
    printf("adraw limits, log grid case: lim1(lim3): tr(%f-a0)=%f, lim2: tr(%f-a0)=%f\n",a0,aspace[0].lim1,aspace[0].lim2p,aspace[0].lim2);
else
    printf("adraw limits, sym-log grid case: lim1: tr(%f)=%f, lim2: tr(%f-a0)=%f, focal point: tr(%f)=%f\n",
                                                      aspace[0].lim3p-a0,aspace[0].lim1,
                                                      aspace[0].lim2p-aspace[0].lim3p,aspace[0].lim2,
                                                      0.0,aspace[0].lim3);
#endif
        }
        //zeroconsumption adjustment or standard case or stop
        if (aspace[0].M<=a0-1+TOLERANCE)
        {   //signal of c1<=0, resend point already prepared in aspace[0].lastAguess
            aspace[0].keep=1;
            //ASSUMING inverted budget calculation (in main body of the program) found exact focal point, reset it
            //also, don't need points left of it
            //also, aa and bb are corrupt because both .lastAguess and .M are wrong
            aa=(a0-aspace[0].baseM)/(a0-aspace[0].baseA);
            bb=aspace[0].baseM-aa*aspace[0].baseA;
            aspace[0].lim2p=MIN(mmax,(mmax-bb)/aa); //upper A (plain)
            aspace[0].lim3p=aspace[0].lastAguess-ZEROCONSUMPTION; //focal point (plain)
            aspace[0].k3=1.0;
            aspace[0].lim1=tr(aspace[0].curr,aspace[0].lim3p-a0);
            aspace[0].lim2=tr(aspace[0].curr,aspace[0].lim2p-aspace[0].lim3p);
            aspace[0].lim3=tr(aspace[0].curr,0);
            
#ifdef VERBOSE_adraw
printf("adraw step %d, base=[%f,%f], returned M = %f, new Aguess = %f <resending the same point with ZEROCONSUMPTION adjustment>\n",(int)aspace[0].ngenerated,aspace[0].baseA,aspace[0].baseM,aspace[0].M,aspace[0].lastAguess);        
printf("adraw limits, log grid case: lim1(lim3): tr(0)=%f, lim2: tr(%f-%f)=%f\n",aspace[0].lim3,aspace[0].lim2p,aspace[0].lim3p,aspace[0].lim2);
#endif
        }
        else if (aspace[0].M<mmax && aspace[0].ngenerated<aspace[0].ntogenerate+0)  //allow few extra points
        {   //generate next Aguess
            aspace[0].keep=1;//keep calculations signal
#ifdef VISUAL_adraw
if (aspace[0].curr->it==VISUAL_adraw_it && aspace[0].curr->ist==VISUAL_adraw_ist && aspace[0].curr->id==VISUAL_adraw_id && aspace[0].ncalls>VISUAL_adraw_ncall) {
 sprintf(Mstring,"point([%1.5f %1.5f],'next','adraw call %d');",aspace[0].lastAguess,aspace[0].M,aspace[0].ncalls);
 mexEvalString(Mstring);
}
#endif
            //step based on log-scale over A with linear rescaling using newer approximations of max A to be generated
            if ((int)aspace[0].ngenerated<(int)aspace[0].k3-1)
                astep=-trinv(aspace[0].curr, aspace[0].lim3
                                             +(aspace[0].k3-1-aspace[0].ngenerated)
                                             *(aspace[0].lim1-aspace[0].lim3)/(aspace[0].k3-1) 
                                                              ) +aspace[0].lim3p - aspace[0].lastAguess;
            else
                astep= trinv(aspace[0].curr, aspace[0].lim3
                                             +(aspace[0].ngenerated-aspace[0].k3+1)
                                             *(aspace[0].lim2-aspace[0].lim3)/(aspace[0].ntogenerate-aspace[0].k3) 
                                                              ) +aspace[0].lim3p - aspace[0].lastAguess;
if (astep<0)
{
    printf("\nadraw call %d\n",aspace[0].ncalls);
    printf("it=%d ist=%d id=%d\n",aspace[0].curr->it,aspace[0].curr->ist,aspace[0].curr->id);
    printf("generated %d points out of %d\n",aspace[0].ngenerated,aspace[0].ntogenerate);
    printf("aa=%1.3f\n",aa);
    printf("bb=%1.3f\n",bb);
    printf("astep=%1.3f\n",astep);
    printf("upper=%1.3f\n",upper);
    printf("lim2p=%1.3f\n",aspace[0].lim2p);
    if (astep*(upper-aspace[0].lastAguess)/(aspace[0].lim2p-aspace[0].lastAguess)<0) mexWarnMsgTxt("Step < 0");
    astep=MAX(astep,1e-5);
}            

//            aspace[0].lastAguess+=astep*(upper-aspace[0].lastAguess)/(aspace[0].lim2p-aspace[0].lastAguess);
// TEMPORARILY : switch off scaling
            aspace[0].lastAguess+=astep;
#ifdef VISUAL_adraw
if (aspace[0].curr->it==VISUAL_adraw_it && aspace[0].curr->ist==VISUAL_adraw_ist && aspace[0].curr->id==VISUAL_adraw_id && aspace[0].ncalls>VISUAL_adraw_ncall) {
 sprintf(Mstring,"point([%1.5f %1.5f]);",aspace[0].lastAguess,aa*aspace[0].lastAguess+bb);
 mexEvalString(Mstring);
 sprintf(Mstring,"point([%1.5f %1.5f]);",upper,mmax);
 mexEvalString(Mstring);
}
#endif
#ifdef VERBOSE_adraw
printf("adraw step %d, base=[%f,%f], returned M = %f, new Aguess = %f\n",(int)aspace[0].ngenerated,aspace[0].baseA,aspace[0].baseM,aspace[0].M,aspace[0].lastAguess);        
#endif
            aspace[0].ngenerated+=1.0;//increment the counter
        }
        else
        {
            aspace[0].lastAguess=-INFINITY; //stop
#ifdef VERBOSE_adraw
printf("adraw step %d, base=[%f,%f], returned M = %f, stop\n",(int)aspace[0].ngenerated,aspace[0].baseA,aspace[0].baseM,aspace[0].M);        
#endif
        }    
    }
    return aspace[0].lastAguess;
}




// #############################################################################################################
static void envelop(
    int it,                 //current period (for utility calculation)
    int ist,                //current state (for utility calculation)
    int nd,                 //number of functions to compare
    int dim0,               //total number of points in individual id specific grids
    double *gridvecs,       //input array with quadruples: (grid point,func value,func2 value,func index)
    double *evfa0,          //expected value corresponding to a0 by id
    double *outgrid,        //pointer to space for resulting grid
    double *outfunc,        //pointer to space for resulting func values
    double *outfunc2,       //pointer to space for resulting func2 values
    double *outthrh,        //pointer to space for thresholds
    double *outindx,        //pointer to space for max-attaining func indeces
    int *outn,              //pointer to int to hold the resulting number of points in outgrid
    int *outm,              //pointer to int to hold the resulting number of points in outthrh and outindex
    void *cache)            //preallocated space to save time on allocations
{   //this function calculates the max of several tabulated functions,
    //merges individual grids these functions are defined at,
    //finds crossings points (thresholds) and incerts them in the resulting grid,
    //records the indexes of max-attaining functions on intervals defined by thresholds.
    int i,j,l;
    double **bigrid; //array of pointers to grid points
    int *bindx; //array indexes at each point
    int *dims; //dimentions of the specific grids
    int *curpt; //offset to current point in specific grid and function value arrays
    double ***grids, ***funcs, ***funcs2; //pointers to specific grids function values
    double grbound; //grid bounds [max of the first points, min of last points] across functions
    int outi,outj; //output indeces
    int check; //tmp indicator
    double cfunc,tmax; //tmp current function value, tmp max in trial
    int ci,cj; //tmp index
    int *kmarks; //marks for threshold function
#ifdef VISUAL
    int saveouti=0;
#endif    
#ifdef VERBOSE_envelop
printf("\nENVELOPE nd=%d dims0=%d\n",nd,dim0);    
#endif

#if MEM==calloc
    //allocate memory
    bigrid = (double**) calloc(dim0,sizeof(double*));
    bindx = (int*) calloc(dim0,sizeof(int));
    dims = (int*) calloc(nd,sizeof(int));
    curpt = (int*) calloc(nd,sizeof(int));
    grids = (double***) calloc(nd,sizeof(double**));
    funcs = (double***) calloc(nd,sizeof(double**));
    funcs2 = (double***) calloc(nd,sizeof(double**));
    for (i=0;i<nd;grids[i++]=(double**)calloc(ngridmax,sizeof(double*)));
    for (i=0;i<nd;funcs[i++]=(double**)calloc(ngridmax,sizeof(double*)));
    for (i=0;i<nd;funcs2[i++]=(double**)calloc(ngridmax,sizeof(double*)));
    kmarks = (int*) calloc(nd,sizeof(int));
    #ifdef VERBOSE_envelop
    printf("Memory allocated by calloc\n");    
    #endif
#elif MEM==cache
    //allocate cache memory (need cast to char for pointer arithmetics)
    bigrid = (double**)  cache;
    bindx  = (int*)      ((char*)cache+dim0*sizeof(double*));
    dims   = (int*)      ((char*)cache+dim0*sizeof(double*)+dim0*sizeof(int));
    curpt  = (int*)      ((char*)cache+dim0*sizeof(double*)+dim0*sizeof(int)+nd*sizeof(int));
    kmarks = (int*)      ((char*)cache+dim0*sizeof(double*)+dim0*sizeof(int)+2*nd*sizeof(int));
    grids  = (double***) ((char*)cache+dim0*sizeof(double*)+dim0*sizeof(int)+3*nd*sizeof(int));
    funcs  = (double***) ((char*)cache+dim0*sizeof(double*)+dim0*sizeof(int)+3*nd*sizeof(int)+nd*sizeof(double**));
    funcs2 = (double***) ((char*)cache+dim0*sizeof(double*)+dim0*sizeof(int)+3*nd*sizeof(int)+2*nd*sizeof(double**));
    for (i=0;i<nd;grids[i++]=(double**) ((char*)cache+dim0*sizeof(double*)+dim0*sizeof(int)+3*nd*sizeof(int)+3*nd*sizeof(double**)+i*ngridmax*sizeof(double*)));
    for (i=0;i<nd;funcs[i++]=(double**) ((char*)cache+dim0*sizeof(double*)+dim0*sizeof(int)+3*nd*sizeof(int)+3*nd*sizeof(double**)+(nd+i)*ngridmax*sizeof(double*)));
    for (i=0;i<nd;funcs2[i++]=(double**)((char*)cache+dim0*sizeof(double*)+dim0*sizeof(int)+3*nd*sizeof(int)+3*nd*sizeof(double**)+(2*nd+i)*ngridmax*sizeof(double*)));
    #ifdef VERBOSE_envelop
    printf("Memory allocated from cache\n");    
    #endif
    for(i=0;i<nd;dims[i++]=0);//initialize dims
#endif
    //Sort by point AND value AND index (see comp1):
    //for equal points higher function value is to be chosen,
    //if also equal, lower index is chosen
    qsort(gridvecs,dim0,sizeof(double)*4,&comp1);
#ifdef VERBOSE_envelop
printf("qsort completed\n");
printdoublearray(nd,evfa0,"evfa0: ");
printdoubletable(dim0,4,gridvecs,"Grid vecs [gridpoint func2 func index]");
#endif
    //Fill out specific arrays of grids, their dimentions, and full grid array
    for (i=0;i<dim0;i++)
    {
        bindx[i]=(int)gridvecs[i*4+3];                  //4th element in gridvecs = index
        bigrid[i]=gridvecs+i*4;                         //1st element in gridvecs = point
        grids[bindx[i]][dims[bindx[i]]]=gridvecs+i*4;   //1st element in gridvecs = point
        funcs[bindx[i]][dims[bindx[i]]]=gridvecs+i*4+2; //3nd element in gridvecs = functional value
        funcs2[bindx[i]][dims[bindx[i]]]=gridvecs+i*4+1; //2nd element in gridvecs = functional2 value
        dims[bindx[i]]++;
    }
#ifdef VERBOSE_envelop
printf("dims = [");
for (i=0;i<nd;printf("%d ",dims[i++]));
printf("]\n");    
#endif
    
    //initialize current points to -1, so first ++ moves to first point
    for (i=0;i<nd;curpt[i++]=-1);
    
    //upper bound for the unified grid should be min of individual last points
    grbound=INFINITY;//min of last grid points
    for (i=0;i<nd;i++)
    {
        if (dims[i]<=0) continue;
        if (*grids[i][dims[i]-1]<grbound) grbound=*grids[i][dims[i]-1];
    }

#ifdef VERBOSE_envelop
printf("bound = %f\n",grbound);    
#endif

    //initialize output indeces
    outi=0; //index for outgrid, outfunc
    outj=0; //index for outthrh, outindex


    //######## MAIN CYCLE ###########
    //Cycle through all points
    i=0;
    while (i<dim0 && *bigrid[i]<=grbound)
    {
#ifdef VERBOSE_envelop
printf("i=%3d outi=%d point=%10f indx=%3d : ",i,outi,*bigrid[i],bindx[i]);    
#endif
        if (outi>0 && outgrid[outi-1]==*bigrid[i])
        {   //SPECIAL CASE: duplicate point => skip, already calculated max here
#ifdef VERBOSE_envelop
printf("duplicate point, skipping\n");
#endif
            curpt[bindx[i]]++;//update pointer on current function
            i++;
            continue;
        }

        //current func value
        cfunc=*funcs[bindx[i]][curpt[bindx[i]]+1];//current function value

        if (outj==0)
        {   //CASE 0 : first point of the common grid
#ifdef VERBOSE_envelop
printf("initial, trying ");
#endif
            tmax = cfunc;
            ci=bindx[i];
            //find max among left ends using utility functions
            for (j=0;j<nd;j++)
            {   //function i value from utility because in first point all function are 
                if (dims[j]<=0) continue; //skip infeasible decisions   
                if (j==bindx[i]) continue;//don't compare to itself
                //func i value
                cfunc=funcvalue(*bigrid[i],j,grids,funcs,curpt,evfa0,it,ist);
#ifdef VERBOSE_envelop
printf("f(%d)=%f ",j,cfunc);
#endif
                if (cfunc>tmax) tmax=cfunc,ci=j;//update max and argmax
                if (cfunc==tmax && ci>j) ci=j;//also choose lower index among equal funcs
            }
            //save point with ci as argmax
            outthrh[outj]=a0;
            outindx[outj]=ci;
            outj++;
            if (outj>=nthrhmax) {error("Not enough space for thresholds. Increase max number of threshold points!");return;} //ERROR:max number of points in out thresholds
            if (ci==bindx[i])
            {   //only save the point if it is actually max (othervise it's below ma0 and should not be on grid)
                outgrid[outi]=*bigrid[i];//point
                outfunc[outi]=tmax;//current function value;
                outfunc2[outi]=*funcs2[bindx[i]][curpt[bindx[i]]+1];
//                if (evfa0[ci]==-INFINITY) outfunc2[outi]=ZEROCONSUMPTION;
//                else outfunc2[outi]=*bigrid[i]-a0;//consumption at credit constraint region ###
                outi++;
                if (outi>=ngridmax) {error("Not enough space for endogenous grid. Increase max number of grid points for M!");return;} //ERROR:max number of points in out grid
#ifdef VERBOSE_envelop
printf("saved >>%d [%f,%f] ",ci,outgrid[outi-1],outfunc[outi-1]);
#endif
            }
#ifdef VERBOSE_envelop
            else
            {
printf(" >>%d (%f is below ma0, so no point saved)",ci,*bigrid[i]);
            }
#endif
        }
        else if ((int)outindx[outj-1]==bindx[i])
        {   // CASE 1
#ifdef VERBOSE_envelop
printf("case 1  (prev max) : ");
#endif
            //current point is from function that was max on last step
            check=0;//check if any function is above
            for (j=0;j<nd;j++)
            {   //over all functions
                if (dims[j]<=0) continue; //skip infeasible decisions   
                if (j!=bindx[i])
                {   //don't compare to itself
                    tmax=funcvalue(*bigrid[i],j,grids,funcs,curpt,evfa0,it,ist);
                    if (cfunc<tmax)
                    {
                        check=1;
                        if (*bigrid[i]!=grbound) break; //unless last point, enough to find just one func above
                        else cfunc=tmax,ci=j;//remember max, maxarg if this is last point in the grid
                    }
                }       
            }
            if (check==0)
            {   //if no other function is above: save point
                outgrid[outi]=*bigrid[i];
                outfunc[outi]=cfunc;
                outfunc2[outi]=*funcs2[bindx[i]][curpt[bindx[i]]+1];//current function2 value
#ifdef VERBOSE_envelop
printf("saved [%f,%f] ",outgrid[outi],outfunc[outi]);
#endif
                outi++;
                if (outi==ngridmax) {error("Not enough space for endogenous grid. Increase max number of grid points for M!");return;} //ERROR:max number of points in out grid
            }
            else
            {   //some other function above: start iterative investigation
                if (*bigrid[i]!=grbound)
                {   //not last point on the grid : mark this fact and skip the point : j is func above
#ifdef VERBOSE_envelop
printf("call to thrh(mode=0) %d vs. %d : ",bindx[i],j);
#endif
                    for (l=0;l<nd;l++) kmarks[l]=(dims[l]>0?0:1);
                    //if (l=thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,bindx[i],j,0,outgrid,outfunc,outfunc2,outthrh,outindx,&outi,&outj,kmarks),l!=0)
                    //    return l; //Error in recursion
                    thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,bindx[i],j,0,outgrid,outfunc,outfunc2,outthrh,outindx,&outi,&outj,kmarks);    
                    if (err[0]) return;
                }
                else
                {   //last point on the grid : find ALL remaining thresholds : ci is max func
#ifdef VERBOSE_envelop
printf("call to thrh(mode=1) %d vs. %d : ",bindx[i],ci);
#endif
                    for (l=0;l<nd;l++) kmarks[l]=(dims[l]>0?0:1);
                    //if (l=thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,bindx[i],ci,1,outgrid,outfunc,outfunc2,outthrh,outindx,&outi,&outj,kmarks),l!=0)
                    //    return l; //Error in recursion
                    thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,bindx[i],ci,1,outgrid,outfunc,outfunc2,outthrh,outindx,&outi,&outj,kmarks);
                    if (err[0]) return;    
                    //save the last max point (remembered as ci)
                    outgrid[outi]=*bigrid[i];
                    outfunc[outi]=funcvalue(*bigrid[i],ci,grids,funcs,curpt,evfa0,it,ist);
                    if (curpt[ci]>=0) outfunc2[outi]=linter2(*bigrid[i],grids[ci]+curpt[ci],funcs2[ci]+curpt[ci]);
                    else if (evfa0[j]==-INFINITY) outfunc2[outi]=ZEROCONSUMPTION;
                    else outfunc2[outi]=*bigrid[i]-a0;
#ifdef VERBOSE_envelop
printf("saved last point [%f,%f] ",outgrid[outi],outfunc[outi]);
#endif
                    outi++;
                    if (outi>=ngridmax) {error("Not enough space for endogenous grid. Increase max number of grid points for M!");return;} //ERROR:max number of points in out grid
                } 
            }
        }
        else if ((int)outindx[outj-1]!=bindx[i])
        {   // CASE 2
#ifdef VERBOSE_envelop
printf("case 2 (!prev max) : ");
#endif
            //current point is not from function that was max on last step
            ci=(int)outindx[outj-1];//previous max function index
            //prev max function value
            tmax=funcvalue(*bigrid[i],ci,grids,funcs,curpt,evfa0,it,ist);
            //first compare to prev max
            if (tmax<cfunc)
            {   //cfunc looks like might be new max
                cj=-1;//compare to all other
                for (j=0;j<nd;j++)            
                {   //over all functions except itself and prev max
                    if (dims[j]<=0) continue; //skip infeasible decisions   
                    if (j==bindx[i]) continue;
                    if (j==ci) continue; 
                    //j func value
                    tmax=funcvalue(*bigrid[i],j,grids,funcs,curpt,evfa0,it,ist);
                    //remember using cfunc and cj
                    if ((cfunc<tmax) || (cfunc==tmax && j<cj)) cfunc=tmax,cj=j;
                }
                if (cj==-1)
                {   //if no other function is above:
                    //find all intersection points by recursive routine
#ifdef VERBOSE_envelop
printf("call to thrh(mode=1) (%d vs. %d) : ",ci,bindx[i]);
#endif
                    for (l=0;l<nd;l++) kmarks[l]=(dims[l]>0?0:1);
//                    if (l=thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,ci,bindx[i],1,outgrid,outfunc,outfunc2,outthrh,outindx,&outi,&outj,kmarks),l!=0)
//                        return l; //Error in recursion
                    thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,ci,bindx[i],1,outgrid,outfunc,outfunc2,outthrh,outindx,&outi,&outj,kmarks);    
                    if (err[0]) return;    
                    //save point
                    outgrid[outi]=*bigrid[i];
                    outfunc[outi]=cfunc;
                    outfunc2[outi]=*funcs2[bindx[i]][curpt[bindx[i]]+1];//current function value
                    outi++;
                    if (outi>=ngridmax) {error("Not enough space for endogenous grid. Increase max number of grid points for M!");return;} //ERROR:max number of points in out grid
                }
                else
                {   //even if not max, need to investigate whether this segment intersects with prev max func before leaving it
#ifdef VERBOSE_envelop
printf("call to thrh(mode=1) (%d vs. %d) : ",ci,cj);
#endif
                    for (l=0;l<nd;l++) kmarks[l]=(dims[l]>0?0:1);
//                    if (l=thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,ci,cj,1,outgrid,outfunc,outfunc2,outthrh,outindx,&outi,&outj,kmarks),l!=0)
//                        return l; //Error in recursion
                    thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,ci,cj,1,outgrid,outfunc,outfunc2,outthrh,outindx,&outi,&outj,kmarks);
                    if (err[0]) return;    
                    //if last point, save it
                    if (*bigrid[i]==grbound)
                    {
                        //last point on the grid : save the last max point (remembered as cj)
                        outgrid[outi]=*bigrid[i];
                        outfunc[outi]=funcvalue(*bigrid[i],cj,grids,funcs,curpt,evfa0,it,ist);
                        if (curpt[cj]>=0) outfunc2[outi]=linter2(*bigrid[i],grids[cj]+curpt[cj],funcs2[cj]+curpt[cj]);
                        else if (evfa0[cj]==-INFINITY) outfunc2[outi]=ZEROCONSUMPTION;
                        else outfunc2[outi]=*bigrid[i]-a0;
                        outi++;
#ifdef VERBOSE_envelop
printf("saved last point [%f,%f] ",outgrid[outi],outfunc[outi]);
#endif
                    } 
                }
            }
            else
            {   //just skip the point, it is under prev max, unless this is last point in the grid
                if (*bigrid[i]!=grbound)
                {
/* The following is good for pictures, not for solution of the model!
                    if (curpt[ci]==-1)
                    {   //save the point to resulting grid just to extend support of the curvy max function
                        outgrid[outi]=*bigrid[i];
                        outfunc[outi]=utility(curr,*bigrid[i]-a0)+discount(curr)*evfa0[ci];
                        outfunc2[outi]=*bigrid[i];
                        outi++;
                        if (outi>=ngridmax) {error("Not enough space for endogenous grid. Increase max number of grid points for M!");return;} //ERROR:max number of points in out grid
                    }
*/                                
#ifdef VERBOSE_envelop
printf(" cfunc=%f < %f (indx=prev max) : prev max above, skipping",cfunc,tmax);
#endif
                }
                else
                {   //last point on the grid: save interpolated max func (ci)
                    outgrid[outi]=*bigrid[i];
                    outfunc[outi]=funcvalue(*bigrid[i],ci,grids,funcs,curpt,evfa0,it,ist);
                    if (curpt[ci]>=0) outfunc2[outi]=linter2(*bigrid[i],grids[ci]+curpt[ci],funcs2[ci]+curpt[ci]);
                    else if (evfa0[ci]==-INFINITY) outfunc2[outi]=ZEROCONSUMPTION;
                    else outfunc2[outi]=*bigrid[i]-a0;
                    outi++;
#ifdef VERBOSE_envelop
printf("saved last point [%f,%f] ",outgrid[outi],outfunc[outi]);
#endif
                }
            }

        }
        curpt[bindx[i]]=MIN(curpt[bindx[i]]+1,dims[bindx[i]]-2);//update pointer on current function (up to next to last point to stay on grid)
        i++;//next point in bigrid
#ifdef VERBOSE_envelop
printf("\n");    
#endif

#ifdef VISUAL
if (ist==VISUAL_ist && outi>0 && saveouti!=outi)
{   /*if new point appeared*/
    sprintf(Mstring,"point([%1.5f %1.5f]);",outgrid[outi-1],outfunc[outi-1]);
    mexEvalString(Mstring);
    saveouti=outi;
}
#endif

    }
#if MEM==calloc
    //clean up
    free(bigrid);
    free(bindx);
    free(dims);
    free(curpt);
    for (i=0;i<nd;free(grids[i++]));
    for (i=0;i<nd;free(funcs[i++]));
    for (i=0;i<nd;free(funcs2[i++]));
    free(grids);
    free(funcs);
    free(funcs2);
    free(kmarks);
#endif
    //output the resulting number of points
    *outn=outi;
    *outm=outj;
}

// #############################################################################################################
static double funcvalue(double x,int j,double***grids,double***funcs,int*curpt,double*evfa0,int it,int ist)
{   //computes the value function taking into account grid, analitical interval and -inf
    // x - point of evaluation
    // j - index of the function
    // all rest are standard arguments
    PeriodVars curr;
    curr.it=it;
    curr.ist=ist;
    curr.id=j;
//printf("\nfuncvalue: it=%d ist=%d id=j=%d evf(j)=%f x=%f utility=%f df=%f anayt=%f\n",it,ist,j,evfa0[j],x,utility(&curr,x-a0),discount(&curr),utility(&curr,x-a0)+discount(&curr)*evfa0[j]);
//PeriodVarsOut(&curr,"funcvalue");
    if (curpt[j]>=0) return linter2(x,grids[j]+curpt[j],funcs[j]+curpt[j]);
    else if (evfa0[j]==-INFINITY) return -INFINITY; //value of consuming zero today and about zero from tomorrow - WRONG if utility of zero consumption is not -inf
    else return utility(&curr,x-a0)+discount(&curr)*evfa0[j];
}

// #############################################################################################################
static int comp1(const void *a,const void *b) 
{   //compares blocks of 4 doubles by 1-3desc-4
    //this comparison ensures that: 
    // 1. when grid points are identical the higher function value is taken first to speed up computation later,
    // 2. when grid points and func levels are identical, function with lower index is always prefered
    if (*(double *)a > *(double *)b) return 1;
    else if (*(double *)a < *(double *)b) return -1;
    else if ((*(double *)a == *(double *)b) && (*((double *)a+2) > *((double *)b+2))) return -1;
    else if ((*(double *)a == *(double *)b) && (*((double *)a+2) < *((double *)b+2))) return 1;
    else if ((*(double *)a == *(double *)b) && (*((double *)a+2) == *((double *)b+2)) && (*((double *)a+3) > *((double *)b+3))) return 1;
    else if ((*(double *)a == *(double *)b) && (*((double *)a+2) == *((double *)b+2)) && (*((double *)a+3) < *((double *)b+3))) return -1;
    else return 0;
}

// #############################################################################################################
static double linter2(double x, double **grid, double **fun)
{   //calculates linear interpolation of function fun at point x BETWEEN points *grid[0] and *grid[1]
    //outbehav defines how interpolation curve behaves outside of interval
    if (x==*grid[0]) return *fun[0]; //special case - for single point evaluation
    //No extrapolation allowed: this is essential for the correct envelop
    if (x<*grid[0]) return -INFINITY;
    if (x>*grid[1]) return -INFINITY;
    return *fun[1]*(x-*grid[0])/(*grid[1]-*grid[0]) + *fun[0]*(*grid[1]-x)/(*grid[1]-*grid[0]);
}

// #############################################################################################################
static void thresholds(
               int nd,          //number of functions
               int *curpt,      //pointers to current points for interpolation
               int *dims,       //dimentions of individual grids
               double ***grids, //individual grids
               double ***funcs, //individual funciton values
               double ***funcs2,//individual function2 values
               int it,          //current time period (for utility calculation)
               int ist,         //current state (for utility calculation)
               double *evfa0,    //individual constants (for utility calculation)
               int pri,         //previous function index
               int nwi,         //next function index
               int mode,        //mode for search of thresholds (0=left-most thresholds, 1=all thresholds)
               double *outgrid, //output grid
               double *outfunc, //output function values
               double *outfunc2,//output function2 values
               double *outthrh, //output thresholds
               double *outindx, //output index (discrete choice policy function)
               int *outi,       //index of next output grid point
               int *outj,       //index of next output threshold point
               int *kmarks)     //array of marks for already compared functions
{   //find thresholds and records them to output variables
    //called by envelop
    //recursive!
    int k,optk;    
    double newpoint,cmax,tmax;
    double fun2[2]; //two values for optimal policy (before and after intersection)
    double br[2]; //brackets for exact intersections
    PeriodVars curr,curr_pri,curr_nwi;

    //fill out period structure
    curr_pri.it=it;
    curr_pri.ist=ist;
    curr_pri.id=pri;
    curr_nwi.it=it;
    curr_nwi.ist=ist;
    curr_nwi.id=nwi;
    
    //mark indexes already called
    kmarks[pri]=1;
    kmarks[nwi]=1;

#ifdef VERBOSE_thresh
printf("\n    th> pri=%d vs. nwi=%d \n",pri,nwi);
printintarray(nd,curpt, "    th> curpt  : ");
printintarray(nd,kmarks,"    th> kmarks : ");
#endif

    //interpolate
    //consider cases when funcs are in exact region
    //property: if both are exact, they can not intersect (are different by a constant)
    //assume: exact function is continuous and monotone increasing (need?)
    if (curpt[pri]==-1 && curpt[nwi]!=-1)
    {   //first line is exact (not linear)
#ifdef VERBOSE_thresh
printf("    th> X-point - first line (pri) is analytical: utility(c)+beta*(%12.11f)\n",evfa0[pri]);
printf("    th> second line points are (%12.11f %12.11f) and (%12.11f %12.11f)\n",*grids[nwi][curpt[nwi]],*funcs[nwi][curpt[nwi]],*grids[nwi][curpt[nwi]+1],*funcs[nwi][curpt[nwi]+1]);
#endif
        if (evfa0[pri]==-INFINITY)
        {   //first line is vertical from -inf
            newpoint=*grids[pri][0];
        }
        else
        {   //brackets
            br[0]=*grids[nwi][curpt[nwi]];  //THESE ARE SENSITIVE TO UTILITY FUNCTION ASSUMPTIONS
            br[1]=MIN(*grids[pri][0],*grids[nwi][curpt[nwi]+1]);  //THESE ARE SENSITIVE TO UTILITY FUNCTION ASSUMPTIONS
            //calculate the exact intersection point
            brsolve(br,br+1,grids[nwi]+curpt[nwi],funcs[nwi]+curpt[nwi],&curr_pri,evfa0+pri);
            if (err[0]) return;
            newpoint=br[0];
        }
        cmax=linter2(newpoint,grids[nwi]+curpt[nwi],funcs[nwi]+curpt[nwi]);
    }
    else if (curpt[pri]!=-1 && curpt[nwi]==-1)
    {   //second line is exact (not linear)
#ifdef VERBOSE_thresh
printf("    th> X-point - second line (nwi) is analytical: utility(c)+beta*(%12.11f)\n",evfa0[nwi]);
printf("    th> first line points are (%12.11f %12.11f) and (%12.11f %12.11f)\n",*grids[pri][curpt[pri]],*funcs[pri][curpt[pri]],*grids[pri][curpt[pri]+1],*funcs[pri][curpt[pri]+1]);
#endif
        if (evfa0[nwi]==-INFINITY)
        {   //first line is vertical from -inf
            newpoint=*grids[nwi][0];
        }
        else
        {   //brackets
            br[0]=*grids[pri][curpt[pri]];  //THESE ARE SENSITIVE TO UTILITY FUNCTION ASSUMPTIONS
            br[1]=MIN(*grids[nwi][0],*grids[pri][curpt[pri]+1]);  //THESE ARE SENSITIVE TO UTILITY FUNCTION ASSUMPTIONS
            brsolve(br,br+1,grids[pri]+curpt[pri],funcs[pri]+curpt[pri],&curr_nwi,evfa0+nwi);
            if (err[0]) return;
            newpoint=br[0];
        }
        cmax=linter2(newpoint,grids[pri]+curpt[pri],funcs[pri]+curpt[pri]);
    }
    else if (curpt[pri]==-1 && curpt[nwi]==-1)
    {
        error("Fatal error in threshold module. Two analytical value functions seem to intersect. Utility is not additively separable in consumption and discrete choices.");
        return;
    }
    else
    {   //both lines are linear
#ifdef VERBOSE_thresh
printf("    th> X-point - both lines are linear\n");
printf("    th> first line points are (%12.11f %12.11f) and (%12.11f %12.11f)\n",*grids[pri][curpt[pri]],*funcs[pri][curpt[pri]],*grids[pri][curpt[pri]+1],*funcs[pri][curpt[pri]+1]);
printf("    th> secnd line points are (%12.11f %12.11f) and (%12.11f %12.11f)\n",*grids[nwi][curpt[nwi]],*funcs[nwi][curpt[nwi]],*grids[nwi][curpt[nwi]+1],*funcs[nwi][curpt[nwi]+1]);
#endif
        if (*grids[pri][curpt[pri]+1]==*grids[pri][curpt[pri]])
        {   //one line vertical
#ifdef VERBOSE_thresh
printf("    th> X-point - first line (pri) is vertical\n");
#endif
            newpoint=*grids[pri][curpt[pri]];
            cmax=   (   newpoint*
                        (   *funcs[nwi][curpt[nwi]+1]-*funcs[nwi][curpt[nwi]]
                        )/
                        (   *grids[nwi][curpt[nwi]+1]-*grids[nwi][curpt[nwi]]
                        )
                    )+
                    (   (   *funcs[nwi][curpt[nwi]]**grids[nwi][curpt[nwi]+1]-*funcs[nwi][curpt[nwi]+1]**grids[nwi][curpt[nwi]]
                        )/
                        (   *grids[nwi][curpt[nwi]+1]-*grids[nwi][curpt[nwi]]
                        )
                    );
        }
        else if (*grids[nwi][curpt[nwi]+1]==*grids[nwi][curpt[nwi]])
        {   //other line vertical
#ifdef VERBOSE_thresh
printf("    th> X-point - second line (nwi) is vertical\n");
#endif
            newpoint=*grids[nwi][curpt[nwi]];
            cmax=   (   newpoint*
                        (   *funcs[pri][curpt[pri]+1]-*funcs[pri][curpt[pri]]
                        )/
                        (   *grids[pri][curpt[pri]+1]-*grids[pri][curpt[pri]]
                        )
                    )+
                    (   (   *funcs[pri][curpt[pri]]**grids[pri][curpt[pri]+1]-*funcs[pri][curpt[pri]+1]**grids[pri][curpt[pri]]
                        )/
                        (   *grids[pri][curpt[pri]+1]-*grids[pri][curpt[pri]]
                        )
                    );
        }
        else if (   (   (   *funcs[nwi][curpt[nwi]+1]-*funcs[nwi][curpt[nwi]]
                        )/
                        (   *grids[nwi][curpt[nwi]+1]-*grids[nwi][curpt[nwi]]
                        )
                    )==
                    (   (   *funcs[pri][curpt[pri]+1]-*funcs[pri][curpt[pri]]
                        )/
                        (   *grids[pri][curpt[pri]+1]-*grids[pri][curpt[pri]]
                        )
                    )
                )
        {   //identical slopes
#ifdef VERBOSE_thresh
printf("    th> X-point - identical slopes\n");
#endif
            newpoint=   (   *grids[pri][curpt[pri]]+*grids[pri][curpt[pri]+1]+*grids[nwi][curpt[nwi]]+*grids[nwi][curpt[nwi]+1]
                        )/4;   
            cmax=   (   newpoint*
                        (   *funcs[nwi][curpt[nwi]+1]-*funcs[nwi][curpt[nwi]]
                        )/
                        (   *grids[nwi][curpt[nwi]+1]-*grids[nwi][curpt[nwi]]
                        )
                    )+
                    (   (   *funcs[nwi][curpt[nwi]]**grids[nwi][curpt[nwi]+1]-*funcs[nwi][curpt[nwi]+1]**grids[nwi][curpt[nwi]]
                        )/
                        (   *grids[nwi][curpt[nwi]+1]-*grids[nwi][curpt[nwi]]
                        )
                    );
        }
        else
        {   //normal case
#ifdef VERBOSE_thresh
printf("    th> X-point - normal case\n");
#endif
            newpoint=(  (   (   *funcs[pri][curpt[pri]]**grids[pri][curpt[pri]+1]-*funcs[pri][curpt[pri]+1]**grids[pri][curpt[pri]]
                            )/
                            (   *grids[pri][curpt[pri]+1]-*grids[pri][curpt[pri]]
                            )
                        )-
                        (   (   *funcs[nwi][curpt[nwi]]**grids[nwi][curpt[nwi]+1]-*funcs[nwi][curpt[nwi]+1]**grids[nwi][curpt[nwi]]
                            )/
                            (   *grids[nwi][curpt[nwi]+1]-*grids[nwi][curpt[nwi]]
                            )
                        )
                     )/
                     (  (   (   *funcs[nwi][curpt[nwi]+1]-*funcs[nwi][curpt[nwi]]
                            )/
                            (   *grids[nwi][curpt[nwi]+1]-*grids[nwi][curpt[nwi]]
                            )
                        )-
                        (   (   *funcs[pri][curpt[pri]+1]-*funcs[pri][curpt[pri]]
                            )/
                            (   *grids[pri][curpt[pri]+1]-*grids[pri][curpt[pri]]
                            )
                        )
                     );
            cmax=   (   newpoint*
                        (   *funcs[nwi][curpt[nwi]+1]-*funcs[nwi][curpt[nwi]]
                        )/
                        (   *grids[nwi][curpt[nwi]+1]-*grids[nwi][curpt[nwi]]
                        )
                    )+
                    (   (   *funcs[nwi][curpt[nwi]]**grids[nwi][curpt[nwi]+1]-*funcs[nwi][curpt[nwi]+1]**grids[nwi][curpt[nwi]]
                        )/
                        (   *grids[nwi][curpt[nwi]+1]-*grids[nwi][curpt[nwi]]
                        )
                    );
        }
    }
    //check max
    optk=-1;
    for (k=0;k<nd;k++)
    {   //check if max
        if (kmarks[k]==1) continue; //skip all already marked (including pri,nwi)
        if (curpt[k]>=0) tmax=linter2(newpoint,grids[k]+curpt[k],funcs[k]+curpt[k]);
        else 
        {
            curr.it=it;
            curr.ist=ist;
            curr.id=k;
//PeriodVarsOut(&curr,"th");
            tmax=utility(&curr,newpoint-a0)+discount(&curr)*evfa0[k];
        }
        if (cmax<tmax) 
        {   //found a higher value
            cmax=tmax;
            optk=k;
            if (mode==0) break;//if only looking for left-most thrh any will do
        }
    }        
    if (optk!=-1)
    {   //cmax is not max!
        //mode: 0 to only find left-most threashold
        //      1 to find all intersections
#ifdef VERBOSE_envelop
printf("recursive %d vs. %d : ",pri,optk);
#endif
        if (mode==0) thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,pri,optk,mode,outgrid,outfunc,outfunc2,outthrh,outindx,outi,outj,kmarks);
        else
        {    
            //one newpoint on the left - call first
            thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,pri,optk,mode,outgrid,outfunc,outfunc2,outthrh,outindx,outi,outj,kmarks);
            //second newpoint on the right - call second
#ifdef VERBOSE_envelop
printf("one more recursive %d vs. %d : ",optk,nwi);
#endif
            thresholds(nd,curpt,dims,grids,funcs,funcs2,it,ist,evfa0,optk,nwi,mode,outgrid,outfunc,outfunc2,outthrh,outindx,outi,outj,kmarks);
        }
    }
    else
    {
        //newpoint is valid (max attained)
        //save newpoint
        outgrid[*outi]=newpoint;
        outfunc[*outi]=cmax;
#ifdef VISUAL
if (ist==VISUAL_ist) {
 sprintf(Mstring,"point([%1.5f %1.5f]);",outgrid[*outi],outfunc[*outi]);
 mexEvalString(Mstring);
}
#endif
        //optimal policy in 2 side limits
        if (curpt[pri]>=0) fun2[0]=linter2(newpoint,grids[pri]+curpt[pri],funcs2[pri]+curpt[pri]);
        else if (evfa0[pri]==-INFINITY) fun2[0]=ZEROCONSUMPTION;
        else fun2[0]=newpoint-a0;//consumption in credit constrained region
        if (curpt[nwi]>=0) fun2[1]=linter2(newpoint,grids[nwi]+curpt[nwi],funcs2[nwi]+curpt[nwi]);
        else if (evfa0[nwi]==-INFINITY) fun2[1]=ZEROCONSUMPTION;
        else fun2[1]=newpoint-a0;//consumption in credit constrained region
        outfunc2[*outi]=(fun2[0]+fun2[1])/2;
#ifdef DEBUGOUT
if (dbgouti<dbgoutn) {
    dbgout[dbgouti+0*dbgoutn]=it;
    dbgout[dbgouti+1*dbgoutn]=ist;
    dbgout[dbgouti+2*dbgoutn]=dbgoutd;
    dbgout[dbgouti+3*dbgoutn]=newpoint;
    if (fun2[0]==ZEROCONSUMPTION) dbgout[dbgouti+4*dbgoutn]=-.9999;
    else dbgout[dbgouti+4*dbgoutn]=fun2[0];
    if (fun2[1]==ZEROCONSUMPTION) dbgout[dbgouti+5*dbgoutn]=-.9999;
    else dbgout[dbgouti+5*dbgoutn]=fun2[1];
    dbgout[dbgouti+6*dbgoutn]=fabs(fun2[1]-fun2[0]);
    dbgouti++;
}
#endif
#ifdef VERBOSE_envelop
printf("saved(%d) [%f,%f,%f] : ",*outi,outgrid[*outi],outfunc[*outi],outfunc2[*outi]);
#endif
        outthrh[*outj]=newpoint;
        outindx[*outj]=nwi;
#ifdef VERBOSE_envelop
printf("thrh(%d) %f %1.0f>>%1.0f : ",*outj,outthrh[*outj],outindx[*outj-1],outindx[*outj]);
#endif
        *outi+=1;
        *outj+=1;
        if (*outi>=ngridmax) {error("Not enough space for endogenous grid. Increase max number of grid points for M!");return;} //ERROR:max number of points in out grid
        if (*outj>=nthrhmax) {error("Not enough space for thresholds. Increase max number of threshold points!");return;} //ERROR:max number of points in out thresholds
        if (evfa0[nwi]==-INFINITY && curpt[nwi]==-1)
        {   //special case vf=-inf
            //don't have to save <- point in func2, it will be saved in main algorithm
            outfunc2[*outi-1]=fun2[0];
            outgrid[*outi-1]=outgrid[*outi-1]-TOLERANCE;
#ifdef VERBOSE_envelop
printf("special case point saved ->%f nothing<- ",outfunc2[*outi-1]);
#endif
        }
        //THRERE MUST BE A CASE HERE FOR evfa0[pri]==-INFINITY     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else if (DOUBLEPOINT_DELTA > 0)
        {    //second point on the intersection - for precise optimal policy
            outfunc2[*outi-1]=fun2[0];
            outgrid[*outi]=newpoint+DOUBLEPOINT_DELTA;
            outfunc[*outi]=cmax;
            outfunc2[*outi]=fun2[1];
#ifdef VERBOSE_envelop
printf("add point saved ->%f %f<- ",outfunc2[*outi-1],outfunc2[*outi]);
#endif
            *outi+=1;
            if (*outi>=ngridmax) {error("Not enough space for endogenous grid. Increase max number of grid points for M!");return;} //ERROR:max number of points in out grid
        }
    }
}

// #############################################################################################################
static void brsolve(double* br0,double* br1,double **grids,double **funcs,PeriodVars *curr,double *evfa0)
{   //this function implements successive bisection method to find intersection of exact and linear funcs
    //called by thresholds()
    //recursive (all args by reference to save memory)
    //ASSUMES consumption == M-a0
    double fval[2];
#ifdef VERBOSE_envelop
    static int cnt=0;
    cnt++;
    printf("br%d=[%f,%f] ",cnt,*br0,*br1);        
#endif
    //functional values
    fval[0]=utility(curr,*br0-a0)+discount(curr)**evfa0;
    fval[1]=utility(curr,*br1-a0)+discount(curr)**evfa0;

    //check for errors        
    if ( SIGN(fval[0] - linter2(*br0,grids,funcs)) ==
         SIGN(fval[1] - linter2(*br1,grids,funcs)) )
    {
            printf("\n\nbr=[%20.16f,%20.16f]\nfunc=[%20.16f,%20.16f]\nlinter=[%20.16f,%20.16f]\n",*br0,*br1,fval[0],fval[1],linter2(*br0,grids,funcs),linter2(*br1,grids,funcs));
            error("Fatal error in braketing module! Solution is outside of brackets.");
            return;
    }
    if (*br0>*br1)
    {
            error("Fatal error in braketing module! Bracket limits reversed.");
            return;
    }
    //PRESISION RULE
    if (fabs(*br0-*br1)<2*DOUBLEPOINT_DELTA || fabs(fval[0]-fval[1])<DOUBLEPOINT_DELTA)
    {   //put the answer into br0 
        *br0=(*br0+*br1)/2;
        return;
    }

    //do the bracketing
    if ( SIGN(fval[0] - linter2(*br0,grids,funcs)) ==
         SIGN(utility(curr,((*br0+*br1)/2)-a0)+discount(curr)**evfa0 - linter2((*br0+*br1)/2,grids,funcs)) )
    {
        *br0=(*br0+*br1)/2;
        return brsolve(br0,br1,grids,funcs,curr,evfa0);
        if (err[0]) return;
    }
    if ( SIGN(fval[1] - linter2(*br1,grids,funcs)) ==
         SIGN(utility(curr,((*br0+*br1)/2)-a0)+discount(curr)**evfa0 - linter2((*br0+*br1)/2,grids,funcs)) )
    {
        *br1=(*br0+*br1)/2;
        return brsolve(br0,br1,grids,funcs,curr,evfa0);
        if (err[0]) return;
    }
}

static void PeriodVarsOut(PeriodVars *v,const char *prompt)
{   //outputs SimPeriod value
    int i;
    printf("\nit=%3d ist=%3d [",v->it,v->ist);
    for (i=0;i<nnst;printf("%4.2f ",v->st[i++]));
    printf("] cash=%7.2f shock=%5.2f ",v->cash,v->shock);
    printf("id=%3d [",v->id);
    for (i=0;i<nnd;printf("%4.2f ",v->dc[i++]));
    printf("] savings=%7.2f",v->savings);
    printf(" :%s\n",prompt);
}


