function model=model2(varargin)
%This function creates a consumption-savings model for the retirement phase

%default parameters	
%--------------------------------------------------------------------------------------------------------
p.T =                  3                       ;% terminal age
p.ngridm =             100                     ;% number of grid points
p.nquad =			   10 					   ;% number of quadrature points
p.mmax =               100                     ;% mmax
p.cc =                 0.0                     ;% credit constraint
%--------------------------------------------------------------------------------------------------------
p.df =                 1                       ;% discount factor
p.rho =                0.0                     ;% 1-crra (0=log)
p.r = 				   0.0					   ;% risk free return
p.sigma =              0.0                     ;% simga param in log-normal distributions of returns
p.duw = 			   1.0 					   ;% disutility of work
p.wage =			   5 					   ;% wage income
%--------------------------------------------------------------------------------------------------------
parseargs(p,varargin,'q'); %parse passed arguments (and get rid of p.)

%main definitions
model=egdstmodel('model2 for ZICE2014 lecture','tmp_egm_model2');
%time
model.t0=1;
model.T=T;
%grids
model.mmax=mmax;          %max cash-at-hand                    
model.ngridmax=10*ngridm; %max number of grid point in M
model.ngridm=ngridm;      %standard number of grid points in M
model.nthrhmax=ngridm;       %max number of treshold points in TH
model.ny=nquad;           %number of points in discrete representation of income shocks
model.a0=0;               %credit constraint
%states and decisions
model.s={'Labour market state',{0,'retired',1,'working'}};
model.d={'Retirement decision',{0,'Retirement',1,'Work'}};
%feasibility of states
model.feasible={'defaultfeasible',true};
%transition probabilities
model.trpr={'dc1==0',[1 0;1 0]};
model.trpr={'dc1==1',[0 1;0 1]};
%choice sets
model.choiceset={'defaultallow',true};
model.choiceset={'ist==0 && id==1','Retirement is absorbing'};
%utility
model.u={'utility','(fabs(rho)<1e-10?log(consumption):(pow(consumption,rho)-1)/rho)  - (id?duw:0.0)'};
model.u={'marginal','pow(consumption,rho-1)'};
model.u={'marginalinverse','pow(mutility,1/(rho-1))'};
model.param={'rho','1-crra parameter',rho};
model.param={'duw','scale parameter for disutility of work',duw};
%extrapolation curvature function
model.u={'extrap','pow(x,rho)'};
%intertemporal budget
model.budget={'cashinhand','savings*(1+r)*shock  + (id?wage:0.0)'};
model.budget={'marginal','(1+r)*shock'};
model.param={'r','risk free return',r};
model.discount='df';
model.param={'df','discount factor',df};
model.param={'wage','Workers wage',wage};
%credit constraint
model.a0=cc;
%shock
model.shock='lognormal';
model.shock={'sigma','sig'};
model.shock={'mu','-sigma*sigma/2'}; %to keep the expectation = 1
model.param={'sig','sigma parameter in lognormal return',sigma};

%compile
% model.cflags.VERBOSE=9; %<=8 gives text output, 9 gives graphs of d-specific values and envelopes
% model.cflags.VISUAL_ist=1; %state point where the graphs should be drawn
model.compile

% model.solve %will call solver outside

end

% end%function