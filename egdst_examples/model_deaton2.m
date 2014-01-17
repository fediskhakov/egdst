%Deterministic deaton model with interest=discount and credit constraint a0==0

delete out.txt %diary
diary out.txt %diary

deaton1=egdstmodel('deaton1','tmp_deaton1');
%time and space
deaton1.t0=1;             %initial period
deaton1.T=25;             %terminal period
%grids
deaton1.mmax=50;          %max cash-on-hand
deaton1.ngridmax=1000;    %max number of grid point in assets
deaton1.ngridm=100;       %standard number of grid points in assets
deaton1.nthrhmax=10;      %max number of treshold points in TH
deaton1.ny=2;             %number of points in discrete representation of income shocks
%states and discrete choices
deaton1.s={'Singleton state',{0,'dummy state'}};
deaton1.trpr={'true',[1]};
deaton1.feasible={'defaultfeasible',true};
deaton1.d={'Dummy decision',{0,'dummy decision'}};
deaton1.choiceset={'defaultallow',true};
%utility
deaton1.u={'utility','log(consumption)'};
deaton1.u={'marginal','1/consumption'};
deaton1.u={'marginalinverse','1/mutility'};
%extrapolation curvature function
deaton1.u={'extrap','log(x)'};
%intertemporal budget
deaton1.budget={'cashinhand','savings*(1+interest)+income_level'};
deaton1.budget={'marginal','1+interest'};
deaton1.discount='1/(1+interest)';
deaton1.param={'interest','return on savings',0.01};
deaton1.eq={'income_level','Realized income','income*shock','next'};
deaton1.param={'income','income (times multiplicator shock)',1.25};
%credit constraint
deaton1.a0=0;
%shock (no shock)
deaton1.shock='lognormal';
deaton1.shock={'sigma','0'};
deaton1.shock={'mu','0'};

%Verbosity
deaton1.cflags.VERBOSE=0;

%Compile and solve
deaton1.compile
deaton1.solve

%Solution plots
try
    deaton1.plot1('c');
    deaton1.plot1('vf','it=[1:5:25 25]');
catch er
end

%Simulate and plot
deaton1.sim([1 8.0]);
try
    deaton1.plot2('mack');
catch er
end

