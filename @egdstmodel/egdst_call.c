/*--------------------------------------------------------------------------------------   
     Calls model parts from MATLAB
     By: Fedor Iskhakov
     INPUT: egdstmodel object containing the model
            string with function to call
            array of required arguments for that function
     OUTPUT: result of the call (double)
  --------------------------------------------------------------------------------------*/
#include "egdst_lib.h"

static double vf(PeriodVars *period);

static mxArray *M;//solution to the model (cell array)
static mxArray *D;//solution to the model (cell array)


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   // The gate to be called from MatLab
    double *arg;
    PeriodVars curr,next;
    double *res;
    int i,narg,sw;
    
    if (nrhs != 3) mexWarnMsgTxt("Error in call(): wrong number of inputs!");
    if (nlhs != 1) mexWarnMsgTxt("Error in call(): wrong number of outputs!");
    
    // Parce inputs
    Model=(mxArray *) prhs[0];
    parseModel();
    loadparameters();
    //preallocate the solution to save memory (mxGetProperty actually takes a copy)
    M=(mxArray*) mxGetProperty(Model,0,"M"); //solution from the model (for all it, ist)
    D=(mxArray*) mxGetProperty(Model,0,"D"); //solution from the model (for all it, ist)
    if (M==NULL || D==NULL) mexErrMsgTxt("Error: the model has not yet been solved!");  
    sw=(int) mxGetScalar(prhs[1]); //switch

    //create output
    narg=mxGetM(prhs[2]);
    plhs[0]=mxCreateDoubleMatrix(narg,1,mxREAL);
    res=mxGetPr(plhs[0]);
    
//printf("sw=%d\n",sw);

    //allow for vector input
    for (i=0;i<narg;i++)
    {
        arg=mxGetPr(prhs[2])+i; //pointer to first argument

        //Common arguments 0=it 1=ist 2=id and their checks
        curr.it=(int)arg[0*narg]-t0;//transform to [0,T-t0]
        if (curr.it<0 || curr.it>T-t0) {mexWarnMsgTxt ("call(): it is outside of admissible range: must be in [t0,T]");sw=-1;}
        curr.ist=(int)arg[1*narg]-1;//transform to base0
        if (curr.ist<0 || curr.ist>nst) {mexWarnMsgTxt ("call(): ist is outside of admissible range: must be in [1,nst]");sw=-1;}
        if (mxGetN(prhs[2])>2 && sw!=6) 
        {
            curr.id=(int)arg[2*narg]-1;//transform to base0
            if (curr.id<0 || curr.id>nd)    {mexWarnMsgTxt ("call(): id is outside of admissible range: must be in [1,nd]");sw=-1;}
        }
//printf("it=%d ist=%d id=%d\n",curr.it,curr.ist,curr.id);        
        
        // Make the call
        switch (sw) {
        case 1:
            //utility
            //INPUT ARGS: it,ist,id,consumption
            if (mxGetN(prhs[2])!=4) {mexWarnMsgTxt ("call() Wrong number of arguments for utility!");return;}
            if (arg[3*narg]>mmax-a0) {mexWarnMsgTxt ("call(): consumption is above mmax-a0");res[i]=mxGetNaN();}
            else res[i]=utility(&curr,arg[3*narg]);
            break;
        case 2:
            //mariginal utility
            //INPUT ARGS: it,ist,id,consumption
            if (mxGetN(prhs[2])!=4) {mexWarnMsgTxt ("call(): Wrong number of arguments for utility!");return;}                                       
            if (arg[3*narg]>mmax-a0) {mexWarnMsgTxt ("call(): consumption is above mmax-a0");res[i]=mxGetNaN();}
            else res[i]=utility_marginal(&curr,arg[3*narg]);
            break;
        case 3:
            //discount
            //INPUT ARGS: it,ist
            if (mxGetN(prhs[2])!=2) {mexWarnMsgTxt ("call(): Wrong number of arguments for discount!");return;}
            res[i]=discount(&curr);
            break;
        case 4:
            //budget
            //INPUT ARGS: it,ist,id,savings,ist1,shock
            if (mxGetN(prhs[2])!=6) {mexWarnMsgTxt ("call(): Wrong number of arguments for budget!");return;}
            next.it=curr.it+1;
            next.savings=arg[3*narg];
            next.ist=(int)arg[4*narg]-1;//transform to base0
            next.shock=arg[5*narg];
            //check domain for these args
            if (next.it<0 || next.it>T-t0) {mexWarnMsgTxt ("call(): it+1 is outside of admissible range: must be in [t0,T]");res[i]=mxGetNaN();}
            else if (next.savings<a0)  {mexWarnMsgTxt ("call(): savings are below a0");res[i]=mxGetNaN();}
            else if (next.ist<0 || next.ist>nst) {mexWarnMsgTxt ("call(): ist1 is outside of admissible range: must be in [1,nst]");sw=-1;}
            else res[i]=cashinhand(&curr,&next);
            break;
        case 5:
            //marginal budget
            //INPUT ARGS: it,ist,id,savings,ist1,shock
            if (mxGetN(prhs[2])!=6) {mexWarnMsgTxt ("call(): Wrong number of arguments for budget!");return;}
            next.it=curr.it+1;
            next.savings=arg[3*narg];
            next.ist=(int)arg[4*narg]-1;//transform to base0
            next.shock=arg[5*narg];
            //check domain for these args
            if (next.it<0 || next.it>T-t0) {mexWarnMsgTxt ("call(): it+1 is outside of admissible range: must be in [t0,T]");res[i]=mxGetNaN();}
            else if (next.savings<a0)  {mexWarnMsgTxt ("call(): savings are below a0");res[i]=mxGetNaN();}
            else if (next.ist<0 || next.ist>nst) {mexWarnMsgTxt ("call(): ist1 is outside of admissible range: must be in [1,nst]");sw=-1;}
            else res[i]=cashinhand_marginal(&curr,&next);
            break;
        case 6:
            //value function
            //INPUT ARGS: it,ist,cash
            if (mxGetN(prhs[2])!=3) {mexWarnMsgTxt ("call(): Wrong number of arguments for value function!");return;}
            curr.cash=arg[2*narg];
            if (curr.cash>mmax) {mexWarnMsgTxt ("call(): cash in hand is above");res[i]=mxGetNaN();}
            else if (curr.it==T-t0) res[i]=utility(&curr,MAX(0,curr.cash));
            else res[i]=vf(&curr);
            break;
        default:
            res[i]=mxGetNaN();
        }    
    }
    
};

static double vf(PeriodVars *curr)
{   //computes value function
    mxArray* cll;
    int i,it,ist;
    int nm,nth,ith;
    double ma0,evf,*th,*cllpr,c,vfi;

    //unpack
    it=curr->it;
    ist=curr->ist;
    cll = mxGetCell(M,ist+it*nst); //solution from model object, proper cell
    if (cll==NULL) {mexWarnMsgTxt ("Solution missing for given it,ist..");return mxGetNaN();}
    cllpr=mxGetPr(cll);
    nm = (int) mxGetM(cll);//number of points in endogenous grid
    ma0=*(cllpr+0*nm+1);//m0 = second value in first column
    evf=*(cllpr+3*nm+0);//expected vf at m0 = first value in 4th column

    /*consumption in the second column of M*/
    c=linter(curr->cash,nm,cllpr,cllpr+nm);
    if (err[0]) mexWarnMsgTxt(err);
    curr->savings=curr->cash-c;

    /*optimal discrete choice*/
    cll = mxGetCell(D,ist+it*nst);
    nth=mxGetM(cll);//number of threasholds 
    th=mxGetPr(cll)+nth; //second column in D are TH
    ith=0;
    while (ith<nth && curr->cash>=th[ith]) ith++;
    curr->id=(int)*(mxGetPr(cll)+ith-1);

    /*value function for below M(a0) is exact, forth column of M*/
    if (curr->cash<ma0 && evf>-INFINITY) vfi=utility(curr,c)+discount(curr)*evf;
    else if (curr->cash<ma0 && evf==-INFINITY) vfi=-mxGetInf();
    else vfi=linter(curr->cash,nm,cllpr,cllpr+3*nm);
    if (err[0]) mexWarnMsgTxt(err);

    return vfi;
}


