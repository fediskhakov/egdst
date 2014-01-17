%Menu script to run examples in this directory

% get ready
if ismac
    setenv('DYLD_LIBRARY_PATH',''); %to enable mex with GCC4.8 under XCode5
end
clear
clear classes
close all
addpath(['.' filesep],['..' filesep]);
clc

%create randstream.mat
if numel(dir('randstream.mat'))==0
    r=rand(1000,1000);
    save ('randstream.mat','r');
    clear r;
end    

%show help of egdstmodel class
%{
ch=input('Display background info on egdstmodel class (''help egdstmodel'')? [y/N] ','s');
if ~isempty(ch) || (ischar('ch') && numel(ch)>0)
    if ismember(ch(1),{'y','Y'})
        help egdstmodel
    end
end
%}

%choice of model script
fprintf('Models in the showroom:\n');
mn={};
for d=dir'
    if ~d.isdir && strcmp(d.name(1:5),'model')
        mn{1,end+1}=d.name(7:end-2);
        mn{2,end}=d.name;
        fprintf('[%d] %12s %s',size(mn,2),mn{1,end}(1:min(end,15)),help(mn{2,end}));
    end
end
ch=input('Which one to run? ','s');
if isempty(ch)
    return
end
fprintf (['Starting ' mn{2,str2num(ch)} '\n\n\n\n\n\n']);

run(mn{2,str2num(ch)}(1:end-2));

%edit(mn{2,str2num(ch(1))}(1:end-2));

