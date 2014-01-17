%Cake eating problem with discount factor = 0.75

delete out.txt %diary
diary out.txt %diary

cake1=egdstmodel('cake1','tmp_cake1');
%time and space
cake1.t0=1;             %initial period
cake1.T=25;             %terminal period
%grids
cake1.mmax=10;          %cake size
cake1.ngridmax=1000;    %max number of grid point in assets
cake1.ngridm=100;       %standard number of grid points in assets
cake1.nthrhmax=10;      %max number of treshold points in TH
cake1.ny=2;             %number of points in discrete representation of income shocks
%states and discrete choices
cake1.s={'Singleton state',{0,'dummy state'}};
cake1.trpr={'true',[1]};
cake1.feasible={'defaultfeasible',true};
cake1.d={'Dummy decision',{0,'dummy decision'}};
cake1.choiceset={'defaultallow',true};
%utility
cake1.u={'utility','log(consumption)'};
cake1.u={'marginal','1/consumption'};
cake1.u={'marginalinverse','1/mutility'};
%extrapolation curvature function
cake1.u={'extrap','log(x)'};
%intertemporal budget
cake1.budget={'cashinhand','savings'};
cake1.budget={'marginal','1'};
cake1.discount='.75';
%credit constraint
cake1.a0=0;
%shock (no shock)
cake1.shock='normal';
cake1.shock={'sigma','0'};
cake1.shock={'mu','0'};

%Verbosity
cake1.cflags.VERBOSE=0;

%Compile and solve
cake1.compile
cake1.solve

%Solution plots
try
    cake1.plot1('c');
    cake1.plot1('vf','it=[1:5:25 25]');
catch er
end

%Simulate and plot
cake1.sim([1 8.0]);
try
    cake1.plot2('mack');
catch er
end

