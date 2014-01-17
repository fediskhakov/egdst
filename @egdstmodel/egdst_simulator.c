/*--------------------------------------------------------------------------------------   
     This function produces simulations of decision makers' behavior based on the
     solution of the egdst model
     By: Fedor Iskhakov
     
     INPUT: model
            the way to use random stream: 1 (each) = same random numbers for all simulations
                                          0 (all)  = different random numbers for different simulations (default)
            
     OUTPUT: Error code (0 for clean)
             Matrix with simulation results
                run model.smslabels' for description!
    
     Convention about variables: t - time, st - state, d - decision
                                 n# - quantity of #
                                 i# - index of #
      in function arguments: it = index of time period t (from T-t0 to 0)
                             ist = index of state in period t
                             ist1 = index of statue in period t+1
                             id = index of decision in period t (relevant for budget at period t+1)

  --------------------------------------------------------------------------------------*/
#include "egdst_lib.h"

#define VERBOSE_randidOFF

#if VERBOSE==1
#define VERBOSE_sims
#endif

//Globals for this file only
static double *sims;//output space pointer
static int nsimout;//dimention of simulated output
static mxArray *M;//solution to the model (cell array)
static mxArray *D;//solution to the model (cell array)

//Structure for sim variables
typedef struct sim_curr_variables {PeriodVars curr; double c; double vf; double mu; double sigma; double eqs[NREQ]; struct sim_curr_variables *next;} SimPeriod;
typedef struct sim_weights {double w; int ist;} SimW;

//Function declarations
static void simsoutput(int isim,SimPeriod *period);
static void policy(SimPeriod *period,int dovf);
static int simulator(int isim,int ist0,double m0,double *randstream,int *map);
static void SimPeriodOut(SimPeriod *v,int type,const char *prompt);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   // The gate to be called from MatLab
    double *init,*randstream;
    int i,ii,isim,nsim,rndtype;
    mwSize dims[3];
    int *map; //map of discrete variables in st
    
    if (nrhs != 2) mexErrMsgTxt("Error: wrong number of inputs!");
    if (nlhs != 1) mexErrMsgTxt("Error: wrong number of outputs!");
    
    // Parce inputs
    Model=(mxArray *) prhs[0];
    parseModel();
    loadparameters();
    init=(double *) mxGetPr(mxGetProperty(Model,0,"init"));
    nsim=(int) mxGetM(mxGetProperty(Model,0,"init"));
    randstream=(double *) mxGetPr(mxGetProperty(Model,0,"randstream"));
    rndtype=(int)mxGetScalar((mxArray *)prhs[1]);
    //preallocate the solution to save memory (mxGetProperty actually takes a copy)
    M=(mxArray*) mxGetProperty(Model,0,"M"); //solution from the model (for all it, ist)
    D=(mxArray*) mxGetProperty(Model,0,"D"); //solution from the model (for all it, ist)
    if (M==NULL || D==NULL) mexErrMsgTxt("Error: the model has not yet been solved!");  

    //check length of randstream (conditional on rndtype)
#define MINUNITRND 4 //per sim per period   !!! important because simulator() assumes there is enough random numbers!
    if (rndtype==1 && mxGetNumberOfElements(mxGetProperty(Model,0,"randstream"))<MINUNITRND*(T-t0+1)) 
        mexErrMsgTxt("Error: randstream is too short even to be re-used for each simulated agent!");
    else if (rndtype==0 && mxGetNumberOfElements(mxGetProperty(Model,0,"randstream"))<MINUNITRND*nsim*(T-t0+1)) 
        mexErrMsgTxt("Error: randstream is too short to be used for all requested simulations!");

    //mark discrete and continous states
    map=(int*) calloc(2*nnst+1,sizeof(int));
    for (i=0;i<nnst;i++) map[i]=(mxIsLogicalScalarTrue(mxGetField(mxGetProperty(Model,0,"s"),i,"discrete"))?1:0);
    for (i=0;i<nnst;map[nnst]+=(map[i++]==1?0:1));//additional last element contains the total number of discretes
    //fill out indexes of continous vars
    ii=1;
    for (i=0;i<nnst;i++) if (map[i]==0) map[nnst+ii++]=i;
    //map values:
    //  0..nnst-1           1 if st is discrete and 0 otherwise
    //  nnst                total number of continous variables
    //  nnst+1..nnst+nnst   indexes of continous variables in order (max number of values is nnst when all st are continuous)
    //printintarray(2*nnst+1,map,"\nmap: ");

    //discrete or continous st
    if (map[nnst]==0) byval=0; //all st vars are discrete --> use indexes in model specific functions
    else byval=1; //there are continuous vars --> use values in model specific functions

    //dimentions of the output
    nsimout=11+nnst+nnd+mxGetNumberOfElements(mxGetProperty(Model,0,"eq"));

    // output
    dims[0]=nsimout;
    dims[1]=T-t0+1;
    dims[2]=nsim;
    plhs[0]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
   
    sims=(double *) mxGetPr(plhs[0]);
    // initialize with NaNs
    for (i=1;i<dims[0]*dims[1]*dims[2];sims[i++]=mxGetNaN());
    
    //Run simulations
    //PARALLELIZE HERE
    if (rndtype==1)
        for (isim=0;isim<nsim;isim++) 
            simulator(isim,(int)init[0*nsim+isim]-1,init[1*nsim+isim],randstream,map); //base1 to base0 init indexes, same randstream
    else
        for (isim=0;isim<nsim;isim++)
            simulator(isim,(int)init[0*nsim+isim]-1,init[1*nsim+isim],randstream+MINUNITRND*(T-t0+1)*isim,map); //base1 to base0 init indexes, different randstreams
    
    free (map);//clean up
};



//use curr.cash, curr.savings and curr.shock for the same period t !!!
static void simsoutput(int isim,SimPeriod *period)
{   //records sims output as written in model.sislabels
    int i,it;
    it=period->curr.it;
    sims[isim*nsimout*(T-t0+1)+it*nsimout+0]=period->curr.cash;
    sims[isim*nsimout*(T-t0+1)+it*nsimout+1]=period->c;
    sims[isim*nsimout*(T-t0+1)+it*nsimout+2]=period->curr.savings;
    sims[isim*nsimout*(T-t0+1)+it*nsimout+3]=period->vf;
    sims[isim*nsimout*(T-t0+1)+it*nsimout+4]=(double) period->curr.id;
    sims[isim*nsimout*(T-t0+1)+it*nsimout+5]=(double) period->curr.ist;
    sims[isim*nsimout*(T-t0+1)+it*nsimout+6]=period->mu;
    sims[isim*nsimout*(T-t0+1)+it*nsimout+7]=period->sigma;
    sims[isim*nsimout*(T-t0+1)+it*nsimout+8]=period->curr.shock;
    sims[isim*nsimout*(T-t0+1)+it*nsimout+9]=utility(&(period->curr),period->c);
    sims[isim*nsimout*(T-t0+1)+it*nsimout+10]=discount(&(period->curr));
    // for (i=0;i<nnst;               sims[isim*nsimout*(T-t0+1)+it*nsimout+11+i++]=period->curr.st[i]);
    // for (i=0;i<nnd;                sims[isim*nsimout*(T-t0+1)+it*nsimout+11+nnst+i++]=period->curr.dc[i]);
    // for (i=0;i<nsimout-11-nnst-nnd;sims[isim*nsimout*(T-t0+1)+it*nsimout+11+nnst+nnd+i++]=period->eqs[i]);
    for (i=0;i<nnst;i++) sims[isim*nsimout*(T-t0+1)+it*nsimout+11+i]=period->curr.st[i];
    for (i=0;i<nnd;i++) sims[isim*nsimout*(T-t0+1)+it*nsimout+11+nnst+i]=period->curr.dc[i];
    for (i=0;i<nsimout-11-nnst-nnd;i++) sims[isim*nsimout*(T-t0+1)+it*nsimout+11+nnst+nnd+i]=period->eqs[i];
}

static void policy(SimPeriod *period,int dovf)
{   //fills in optimal period->id and period->savings using values describing the period
    //only uses id and ist (optimal policies are only defined on grids)
    // curr.it   - USED
    // curr.ist  - USED
    // curr.st[] - not used
    // curr.id   - FILLED in by this function
    // curr.dc[] - FILLED in by this function
    // curr.cash - USED
    // curr.savings - FILLED in by this function
    // curr.shock - not used
    // c         - FILLED in by this function
    // vf        - FILLED in by this function IF dovf>0
    // sigma     - not used
    // mu        - not used
    // eqs[]     - not used
    mxArray* cll;
    int i,it,ist;
    int nm,nth,ith;
    double ma0,evf,*th,*cllpr;

    //unpack
    it=period->curr.it;
    ist=period->curr.ist;
    //solution from model object, proper cell
    cll = mxGetCell(M,ist+it*nst); //solution from model object, proper cell
    if (cll==NULL) mexErrMsgTxt("Solution not found in model.M");
    cllpr=mxGetPr(cll);
    nm = (int) mxGetM(cll);//number of points in endogenous grid
    ma0=*(cllpr+0*nm+1);//m0 = second value in first column
    evf=*(cllpr+3*nm+0);//expected vf at m0 = first value in 4th column

    /*consumption in the second column of M*/
    period->c=linter(period->curr.cash,nm,cllpr,cllpr+nm);
    if (err[0]) mexErrMsgTxt(err);
    period->curr.savings=period->curr.cash-period->c;

    /*optimal discrete choice*/
    cll = mxGetCell(D,ist+it*nst);
    if (cll==NULL) mexErrMsgTxt("Solution not found in model.D");
    nth=mxGetM(cll);//number of threasholds 
    th=mxGetPr(cll)+nth; //second column in D are TH
    ith=0;
    while (ith<nth && period->curr.cash>=th[ith]) ith++;
    period->curr.id=(int)*(mxGetPr(cll)+ith-1);
    if (!dovf) return;

    //fill out decision
    for (i=0;i<nnd;i++) period->curr.dc[i]=decision(period->curr.id,i+1);
    
    /*value function for below M(a0) is exact, forth column of M*/
    if (period->curr.cash<ma0 && evf>-INFINITY) period->vf=utility(&(period->curr),period->c)+discount(&(period->curr))*evf;
    else period->vf=linter(period->curr.cash,nm,cllpr,cllpr+3*nm);
    if (err[0]) mexErrMsgTxt(err);
}



// #############################################################################################################
static int simulator(int isim,int ist0,double m0,double *randstream,int *map)
{   //simulation for single person
    int i,it,irnd=0;
    int j0,j1,ii,ist1,id1;
    SimPeriod periodspace[2],*cperiod, *nperiod;
    double rrr, rrr1, rrr2;//random vars
    double pr,th,wc,wvf;
    int nwts;
    SimW *wts;

    //checks
    if (ist0<0 || ist0>=nst) {mexWarnMsgTxt("Initial state index st(0) out of bounds! Moving to next simulation.");return 2;}
    if (m0<a0 || m0>mmax) {mexWarnMsgTxt("Initial money-at-hand out of bounds! Moving to next simulation.");return 3;}

    //make the loop and initialize pointers
    periodspace[0].next=&(periodspace[1]);
    periodspace[1].next=&(periodspace[0]);
    cperiod=&(periodspace[0]);

    //for the continuous st models
    if (map[nnst]>0)
    {   //create weights for interpolaton over st in continuous case
        nwts=(int)pow(2,map[nnst]);//number of weights for multinomial interpolation
        wts=calloc(nwts,sizeof(SimW));
    }
    
#ifdef VERBOSE_sims
printf("\nisim=%d\n",isim);
SimPeriodOut(NULL,1,"");
#endif
    //main simulation cycle - over t
    for (it=0;it<T-t0+1;it++)
    {
        if (it==0)
        {   //initial period - fill out the state vars
            //state
            cperiod->curr.it=0;
            cperiod->curr.ist=ist0;
            for (i=0;i<nnst;i++) cperiod->curr.st[i]=state(ist0,i+1);
            cperiod->curr.cash=m0;
            if (!feasible(&(cperiod->curr))) {mexWarnMsgTxt("Initial state not feasible! Moving to next simulation.");return 4;}
            //no shocks in t0
            cperiod->mu=mxGetNaN();
            cperiod->sigma=mxGetNaN();
            cperiod->curr.shock=mxGetNaN();
            //equations
            eqs_sim(&(cperiod->curr),NULL,cperiod->eqs);
        }
        else
        {   //transition from cperiod (it-1) to nperiod (it)
            nperiod=cperiod->next;
            nperiod->curr.it=it;
            nperiod->curr.savings=cperiod->curr.savings;//copy value of savings for correct interpretation in user defined functions

            //rand
            rrr=randstream[irnd++];//draw random number: trpr
            rrr1=randstream[irnd++];//draw random number: shock (same rand for all ist1)
            rrr2=randstream[irnd++];//draw random number: survival
            pr=0;

            //survival
            if (rrr2>survival(&(cperiod->curr))) break; /*death event*/

            for (nperiod->curr.ist=0;nperiod->curr.ist<nst;nperiod->curr.ist++)
            {   //iterate over state indexes
                //only look at first grid points for the continuous - ist will always have the index of the first grid point
                for (i=0;i<nnst;i++) if (map[i]==0 && (nperiod->curr.ist/(int)stm[nnst+i])%(int)stm[i]!=0) break;
                if (i<nnst) continue;
                //fill in st
                for (i=0;i<nnst;i++) if (map[i]==1) nperiod->curr.st[i]=state(nperiod->curr.ist,i+1);
                trpr_cont(&(cperiod->curr),&(nperiod->curr));//continous vars exact!
                //skip infeasible
                if (!feasible(&(nperiod->curr))) continue;
                //do the trpr simulation
                if (optim_TRPRnoSH==1) pr=trpr(&(cperiod->curr),&(nperiod->curr),0); //0 = only discrete transitions, trpr is independend of shocks
                else
                {   //trpr is dependent on shocks
                    nperiod->mu=mu_param(&(cperiod->curr),&(nperiod->curr));
                    nperiod->sigma=sigma_param(&(cperiod->curr),&(nperiod->curr));
                    if (nperiod->sigma<=0) nperiod->curr.shock=expectation(&(cperiod->curr),&(nperiod->curr)); //degenerate distribution
                    else nperiod->curr.shock=cdfinv(rrr1,nperiod->mu,nperiod->sigma); //non-degenerate distribution
                    pr=trpr(&(cperiod->curr),&(nperiod->curr),0);
                }
                //finish simulation of trpr
                rrr-=pr; 
                if (rrr<=0) break;
            }//ist next period is fixed, st filled

            if (optim_TRPRnoSH==1)
            {   //calculate shocks for independent case
                nperiod->mu=mu_param(&(cperiod->curr),&(nperiod->curr));
                nperiod->sigma=sigma_param(&(cperiod->curr),&(nperiod->curr));
                if (nperiod->sigma<=0) nperiod->curr.shock=expectation(&(cperiod->curr),&(nperiod->curr)); //degenerate distribution
                else nperiod->curr.shock=cdfinv(rrr1,nperiod->mu,nperiod->sigma); //non-degenerate distribution
            }   
            //calculate money at hand
            nperiod->curr.cash=cashinhand(&(cperiod->curr),&(nperiod->curr));
            //equations
            eqs_sim(&(cperiod->curr),&(nperiod->curr),nperiod->eqs);
            //swap periods to complete transition
            cperiod=nperiod;
        }
        //finilize the iteration

        //policy
        if (map[nnst]==0) policy(cperiod,1); //all st variables are discrete, compute vf
        else
        {   //in case of continous variables randomize consumption and discrete
            //decision using interpolation weights of the continuous states
            for (ii=0;ii<nwts;(wts+ii++)->w=1);//initialize weights
            for (ii=0;ii<nwts;(wts+ii++)->ist=cperiod->curr.ist);//initialize ist (NB all cont grids are at 0 now)
            for (i=0;i<map[nnst];i++)
            {   //over all continous states
                j0=map[nnst+i+1];//index of this state in st
                j1=bxsearch(cperiod->curr.st[j0],(double*)stgrids[j0],(int)stm[j0]);//allow extrapolation
                for (ii=0;ii<nwts;ii++)
                {   //update the weights
                    if ((ii/(int)pow(2,i))%2==0) 
                    {   //use j1
                        (wts+ii)->w*=(stgrids[j0][j1+1]-cperiod->curr.st[j0])/(stgrids[j0][j1+1]-stgrids[j0][j1]);
                        (wts+ii)->ist+=stm[nnst+j0]*j1;
                    }
                    else
                    {   //use j1+1
                        (wts+ii)->w*=(cperiod->curr.st[j0]-stgrids[j0][j1])/(stgrids[j0][j1+1]-stgrids[j0][j1]);
                        (wts+ii)->ist+=stm[nnst+j0]*(j1+1);
                    }
                }
            }

//for (ii=0;ii<nwts;ii++) printf("%d cperiod.ist=%d ist=%d cperiod.st=%f wt=%f\n",ii,cperiod->curr.ist,(wts+ii)->ist,cperiod->curr.st[j0],(wts+ii)->w);
#ifdef VERBOSE_randid
            rrr=randstream[irnd++];//draw random number: id 
#else
            rrr=.5;//no randomization: id
#endif
            //interpolate consumption and value function and randomize id-ist over calculated weights
            wc=0;
            wvf=0;
            ist1=-1;
            for (ii=0;ii<nwts;ii++)
            {   
                if ((wts+ii)->w>0)
                {   //skip zero weights
                    //update cperiod
                    cperiod->curr.ist=(wts+ii)->ist;
                    //find policy
                    policy(cperiod,0);//0 to skip vf
                    wc+=cperiod->c*(wts+ii)->w;
                    wvf+=cperiod->vf*(wts+ii)->w;
                    //randomize
                    rrr-=(wts+ii)->w;
                    if (rrr<0 && ist1==-1)
                    {   //remember the randomized 
                        ist1=(wts+ii)->ist;
                        id1=cperiod->curr.id;
                    }
                }
            }
            //update cperiod
            cperiod->c=MIN(wc,cperiod->curr.cash-a0);
            cperiod->curr.savings=cperiod->curr.cash-cperiod->c;
            // cperiod->c=wc;
            // cperiod->curr.savings=cperiod->curr.cash-wc;
            cperiod->vf=wvf;
            cperiod->curr.ist=ist1;
            cperiod->curr.id=id1;
            //fill out decision
            for (ii=0;ii<nnd;ii++) cperiod->curr.dc[ii]=decision(cperiod->curr.id,ii+1);
        } 
#ifdef VERBOSE_sims
SimPeriodOut(cperiod,1,"current period output");
#endif
        //output
        simsoutput(isim,cperiod);
    }//end of main cycle over it

    if (map[nnst]>0) free(wts);//clean up
    return 0;
}

static void SimPeriodOut(SimPeriod *v,int type,const char *prompt)
{   //outputs SimPeriod value
    int i;
    if (v==NULL)
    {   //headings
        printf("%3s%3s [","it","is");
        for (i=0;i<nnst;printf("%4s "," "),i++);
        printf("] %7s %5s %5s %5s","m","mu","sig","shk");
        if (type>0) {
            printf("%3s [","id");
            for (i=0;i<nnd;printf("%4s "," "),i++);
            printf("] %7s","a");
        } 
        else {
            printf("%3s [","id");
            for (i=0;i<nnd;printf("%4s "," "),i++);
            printf("] %7s","a");
        }
        printf(" : address : comment\n");
        for (i=0;i<90;printf("-"),i++);
        printf("\n");
    }    
    else
    {
        printf("%3d%3d [",v->curr.it,v->curr.ist);
        for (i=0;i<nnst;printf("%4.2f ",v->curr.st[i++]));
        printf("] %7.2f %5.2f %5.2f %5.2f",v->curr.cash,v->mu,v->sigma,v->curr.shock);
        if (type>0) {
            printf("%3d [",v->curr.id);
            for (i=0;i<nnd;printf("%4.2f ",v->curr.dc[i++]));
            printf("] %7.2f",v->curr.savings);
        }
        else {
            printf("%3s [",".");
            for (i=0;i<nnd;printf("%4s ","."),i++);
            printf("] %7s",".");
        }
        printf(" : %p : %s\n",v,prompt);
    }
}

