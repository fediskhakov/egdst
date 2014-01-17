%Deterministic retirement model (with discontinuous consumption)

delete out.txt %diary
diary out.txt %diary

retire1=egdstmodel('retire1','tmp_retire1');
%time and space
retire1.t0=1;             %initial period
retire1.T=25;             %terminal period
%grids
retire1.mmax=10;          %max cash-on-hand
retire1.ngridmax=1000;    %max number of grid point in assets
retire1.ngridm=100;       %standard number of grid points in assets
retire1.nthrhmax=10;      %max number of treshold points in TH
retire1.ny=10;             %number of points in discrete representation of income shocks
%states and discrete choices
retire1.s={'Singleton state',{0,'dummy state'}};
retire1.trpr={'true',[1]};
retire1.feasible={'defaultfeasible',true};
retire1.d={'Labour supply',{0,'retire',1,'work'}};
retire1.choiceset={'defaultallow',true};
%utility
retire1.u={'utility','log(consumption)+duw*(id==0)'};
retire1.param={'duw','disutility of work',0.5};
retire1.u={'marginal','1/consumption'};
retire1.u={'marginalinverse','1/mutility'};
%extrapolation curvature function
retire1.u={'extrap','log(x)'};
%intertemporal budget
retire1.budget={'cashinhand','savings+wage_income*(id!=0)'};
retire1.budget={'marginal','1+interest'};
retire1.discount='1/(1+interest)';
retire1.param={'interest','return on savings',0.045};
retire1.eq={'wage_income','Realized wage income','wage*shock','next'};
retire1.param={'wage','wage (times multiplicator shock)',1.05};
%credit constraint
retire1.a0=-5;
%shock (no shock)
retire1.shock='lognormal';
retire1.shock={'sigma','0.25'};
retire1.shock={'mu','-0.5*sigma*sigma'};

%Verbosity
retire1.cflags.VERBOSE=0;

%Compile and solve
retire1.compile
retire1.solve

%Solution plots
try
    retire1.plot1('c');
    retire1.plot1('d');
    retire1.plot1('vf','it=[1:5:25 25]');
catch er
end

%Simulate and plot
retire1.sim([1 0.25]);
try
    retire1.plot2('mack-d1-q1');
catch er
end

