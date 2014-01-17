%Deterministic retirement model (with discontinuous consumption)

delete out.txt %diary
diary out.txt %diary

retire2=egdstmodel('retire2','tmp_retire2');
%time and space
retire2.t0=1;             %initial period
retire2.T=25;             %terminal period
%grids
retire2.mmax=10;          %max cash-on-hand
retire2.ngridmax=1000;    %max number of grid point in assets
retire2.ngridm=100;       %standard number of grid points in assets
retire2.nthrhmax=10;      %max number of treshold points in TH
retire2.ny=10;             %number of points in discrete representation of income shocks
%states and discrete choices
retire2.s={'Singleton state',{0,'dummy state'}};
retire2.trpr={'true',[1]};
retire2.feasible={'defaultfeasible',true};
retire2.d={'Labour supply',{0,'retire',1,'work'}};
retire2.choiceset={'defaultallow',true};
%utility
retire2.u={'utility','log(consumption)+duw*(id==0)'};
retire2.param={'duw','disutility of work',0.5};
retire2.u={'marginal','1/consumption'};
retire2.u={'marginalinverse','1/mutility'};
%extrapolation curvature function
retire2.u={'extrap','log(x)'};
%intertemporal budget
retire2.budget={'cashinhand','savings+wage_income*(id!=0)'};
retire2.budget={'marginal','1+interest'};
retire2.discount='1/(1+interest)';
retire2.param={'interest','return on savings',0.045};
retire2.eq={'wage_income','Realized wage income','wage*shock','next'};
retire2.param={'wage','wage (times multiplicator shock)',1.05};
%credit constraint
retire2.a0=-5;
%shock (no shock)
retire2.shock='lognormal';
retire2.shock={'sigma','0.25'};
retire2.shock={'mu','-0.5*sigma*sigma'};

%Verbosity
retire2.cflags.VERBOSE=0;

%Compile and solve
retire2.compile
retire2.solve

%Solution plots
try
    retire2.plot1('c');
    retire2.plot1('d');
    retire2.plot1('vf','it=[1:5:25 25]');
catch er
end

%Simulate and plot
retire2.sim([1 0.25]);
try
    retire2.plot2('mack-d1-q1');
catch er
end

