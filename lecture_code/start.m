%This script runs all the necessary steps during the lecture

% get ready
if ismac
    setenv('DYLD_LIBRARY_PATH',''); %to enable mex with GCC4.8 under XCode5
end
clear
clear classes
close all
addpath(['.' filesep],['..' filesep]);
clc

% PART 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch=input('Run part 1? (Matlab implementation of simple consumption model) [Y/n] ','s');
if isempty(ch) || ismember(ch(1),{'y','Y'})
	model0;
end %part 1

% PART 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch=input('Run part 2? (model of retirement phase = consumption/savings only) [Y/n] ','s');
if isempty(ch) || ismember(ch(1),{'y','Y'})
	model=model1;
	model.plot1('c')
	model.plot1('vf')
end %part 2

% PART 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch=input('Run part 3? (***) [Y/n] ','s');
if isempty(ch) || ismember(ch(1),{'y','Y'})


end %part 3


