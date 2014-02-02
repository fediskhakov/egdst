%This script runs all the necessary steps during the lecture

% get ready
if ismac
    setenv('DYLD_LIBRARY_PATH',''); %to enable mex with GCC4.8 under XCode5
end
clear
clear classes
close all
addpath(['..' filesep],['..' filesep '..' filesep]); %Need ../.. on the path
clc

% HOW TO GET HELP
% help egdst
% help plot1
% help plot2
% help plot3

% PART 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch=input('Run part 1? (Matlab implementation of simple consumption model) [Y/n] ','s');
if isempty(ch) || ismember(ch(1),{'y','Y'})
	model0;
	%also draws the solution graph
end %part 1

% PART 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch=input('Run part 2? (EGDST model of retirement phase = consumption/savings only) [Y/n] ','s');
if isempty(ch) || ismember(ch(1),{'y','Y'})
	model=model1;
	model.plot1('c');
	model.plot1('vf');
end %part 2

% PART 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch=input('Run part 3? (EGDST retirement model without shocks) [Y/n] ','s');
if isempty(ch) || ismember(ch(1),{'y','Y'})
	model=model2; %compiled but not solved
	% make the solver draw graphs along the way
	%<=8 gives text output, 9 gives graphs of d-specific values and envelopes
	model.cflags.VERBOSE=9;
	%state point where the graphs should be drawn
	model.cflags.VISUAL_ist=1;
	model.compile; %have to recompile with new cflags
	%solve (and make graphs during solver run)
	model.solve;

	%optimal consumption plots (T-1 and T-2)
	model.plot1('c',sprintf('it=%d',model.T-1),'ist=2');
	set(gca,'XLim',[0 20],'YLim',[0 10]);
	model.plot1('c',sprintf('it=%d',model.T-2),'ist=2');
	set(gca,'XLim',[0 20],'YLim',[0 10]);

	%switch off the graph function along the solver
	model.cflags.VERBOSE=0;
	model.compile;

	%extend the problem to more periods
	model.T=25;
	model.solve;
	model.plot1('c',sprintf('it=[%d:%d]',model.T-6,model.T),'ist=2');
	set(gca,'XLim',[0 30],'YLim',[0 30]);
	model.plot1('d','ist=2');

end %part 3

% PART 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch=input('Run part 4? (EGDST retirement model with shocks) [Y/n] ','s');
if isempty(ch) || ismember(ch(1),{'y','Y'})
	model=model2('T',6,'ngridm',5000,'nquad',100); %large grids for nice graphs

	for v=[.10,.20,.22]
		model.setparam('sig',v);
		model.getparam;
		model.solve

		model.plot1('c','it=3','ist=2');
		set(gca,'XLim',[0 20],'YLim',[0 10]);

		model.plot1('c','it=1:3','ist=2');
		set(gca,'XLim',[0 20],'YLim',[0 10]);
	end

end %part 4

% PART 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch=input('Run part 5? (Simulations of the last model) [Y/n] ','s');
if isempty(ch) || ismember(ch(1),{'y','Y'})
	model=model2('T',25,'sigma',.25,'duw',log(5));
	model.solve

	nsims=100; %number of simulated individuals
	sigw=.5;	%sigma in initial wealth distribution
	%arbitrary way to generate initial wealth 
	initwealth=10*exp(norminv(rand(nsims,1)));
	mask=initwealth>0.9*model.mmax;
	initwealth(mask)=0.5*initwealth(mask);
	%initial conditions set according to egdstmodel specifications
	initcond=[2*ones(nsims,1) initwealth];
	%run simulator (random numbers generated and saved in model, see help egdstmodel)
	model.sim(initcond,'own_shocks'); %with individual shocks

	%plot
	model.plot2('mack-ds','sims=1:2');
	model.plot2('m','sims=1:100','combine'); 
	model.plot2('c','sims=1:100','combine');
	model.plot2('d1','sims=1:100','combine');

	%common shocks
	model.sim(initcond,'same_shocks'); %with individual shocks
	model.plot2('m','sims=1:100','combine'); 
end %part 5

