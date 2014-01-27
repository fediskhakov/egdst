%This function creates a consumption-savings model for the retirement phase

function model=model1(varargin)
%default parameters	
%--------------------------------------------------------------------------------------------------------
p.T =                  30                      ;% terminal age
p.ngridm =             100                     ;% number of grid points
p.nquadrature =		   10 					   ;% number of quadrature points
p.mmax =               100                     ;% mmax
p.cc =                 0.0                     ;% credit constraint
%--------------------------------------------------------------------------------------------------------
p.df =                 .95                     ;% discount factor
p.rho =                0.0                     ;% 1-crra (0=log)
p.r = 				   0.05					   ;% risk free return
p.sigma =              0.1                     ;% simga param in log-normal distributions of returns
%--------------------------------------------------------------------------------------------------------
parseargs(p,varargin,'q'); %parse passed arguments (and get rid of p.)

model=egdstmodel('retired','tmp_model1');
%time
model.t0=1;             %initial period
model.T=p.T;            %terminal period
%grids
model.mmax=p.mmax;      %max cash-on-hand
model.ngridm=p.ngridm;  %standard number of grid points in assets
model.ny=p.nquadrature; %number of points in discrete representation of income shocks
%states and discrete choices
model.s={'Singleton state',{0,'dummy state'}};
model.trpr={'true',[1]};
model.feasible={'defaultfeasible',true};
model.d={'Dummy decision',{0,'dummy decision'}};
model.choiceset={'defaultallow',true};
%utility
model.u={'utility','(fabs(rho)<1e-10)?log(consumption):(pow(consumption,rho)-1)/rho'};
model.u={'marginal','pow(consumption,rho-1)'};
model.u={'marginalinverse','pow(mutility,1/(rho-1))'};
model.param={'rho','1-crra parameter',p.rho};
%extrapolation curvature function
model.u={'extrap','pow(x,rho)'};
%intertemporal budget
model.budget={'cashinhand','savings*(1+r)*shock'};
model.budget={'marginal','(1+r)*shock'};
model.param={'r','risk free return',p.r};
model.discount='df';
model.param={'df','discount factor',p.df};
%credit constraint
model.a0=p.cc;
%shock
model.shock='lognormal';
model.shock={'sigma','sig'};
model.shock={'mu','-sigma*sigma/2'}; %to keep the expectation = 1
model.param={'sig','sigma parameter in lognormal return',p.sigma};

%Compile and solve
model.compile;
model.solve;


end%function