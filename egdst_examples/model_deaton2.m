%Stochastic deaton model with interest=discount and credit constraint a0<0

delete out.txt %diary
diary out.txt %diary

deaton2=egdstmodel('deaton2','tmp_deaton2');
%time and space
deaton2.t0=1;             %initial period
deaton2.T=25;             %terminal period
%grids
deaton2.mmax=100;          %max cash-on-hand
deaton2.ngridmax=1000;    %max number of grid point in assets
deaton2.ngridm=100;       %standard number of grid points in assets
deaton2.nthrhmax=10;      %max number of treshold points in TH
deaton2.ny=10;             %number of points in discrete representation of income shocks
%states and discrete choices
deaton2.s={'Singleton state',{0,'dummy state'}};
deaton2.trpr={'true',[1]};
deaton2.feasible={'defaultfeasible',true};
deaton2.d={'Dummy decision',{0,'dummy decision'}};
deaton2.choiceset={'defaultallow',true};
%utility
deaton2.u={'utility','log(consumption)'};
deaton2.u={'marginal','1/consumption'};
deaton2.u={'marginalinverse','1/mutility'};
%extrapolation curvature function
deaton2.u={'extrap','log(x)'};
%intertemporal budget
deaton2.budget={'cashinhand','savings*(1+interest)+income_level'};
deaton2.budget={'marginal','1+interest'};
deaton2.discount='1/(1+interest)';
deaton2.param={'interest','return on savings',0.01};
deaton2.eq={'income_level','Realized income','income*shock','next'};
deaton2.param={'income','income (times multiplicator shock)',1.25};
%credit constraint
deaton2.a0=-25;
%shock (no shock)
deaton2.shock='lognormal';
deaton2.shock={'sigma','0.75'};
deaton2.shock={'mu','0'};

%Verbosity
deaton2.cflags.VERBOSE=0;

%Compile and solve
deaton2.compile
deaton2.solve

%Solution plots
try
    deaton2.plot1('c');
    deaton2.plot1('vf','it=[1:5:25 25]');
catch er
end

%Simulate and plot
deaton2.sim([1 0.25]);
try
    deaton2.plot2('mack-q1');
catch er
end

