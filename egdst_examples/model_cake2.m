%Cake eating problem with discount factor = 0.75

delete out.txt %diary
diary out.txt %diary

cake2=egdstmodel('cake2','tmp_cake2');
%time and space
cake2.t0=1;             %initial period
cake2.T=25;             %terminal period
%grids
cake2.mmax=10;          %cake size
cake2.ngridmax=1000;    %max number of grid point in assets
cake2.ngridm=100;       %standard number of grid points in assets
cake2.nthrhmax=10;      %max number of treshold points in TH
cake2.ny=2;             %number of points in discrete representation of income shocks
%states and discrete choices
cake2.s={'Singleton state',{0,'dummy state'}};
cake2.trpr={'true',[1]};
cake2.feasible={'defaultfeasible',true};
cake2.d={'Dummy decision',{0,'dummy decision'}};
cake2.choiceset={'defaultallow',true};
%utility
cake2.u={'utility','log(consumption)'};
cake2.u={'marginal','1/consumption'};
cake2.u={'marginalinverse','1/mutility'};
%extrapolation curvature function
cake2.u={'extrap','log(x)'};
%intertemporal budget
cake2.budget={'cashinhand','savings'};
cake2.budget={'marginal','1'};
cake2.discount='.75';
%credit constraint
cake2.a0=0;
%shock (no shock)
cake2.shock='normal';
cake2.shock={'sigma','0'};
cake2.shock={'mu','0'};

%Verbosity
cake2.cflags.VERBOSE=0;

%Compile and solve
cake2.compile
cake2.solve

%Solution plots
try
    cake2.plot1('c');
    cake2.plot1('vf','it=[1:5:25 25]');
catch er
end

%Simulate and plot
cake2.sim([1 8.0]);
try
    cake2.plot2('mack');
catch er
end

