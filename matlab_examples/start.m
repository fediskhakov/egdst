%Menu script to run examples in matlab examples directory

%choice of model script
fprintf('Matlab examples on EGM:\n');
mn={};
for d=dir'
    if ~d.isdir && strcmp(d.name(1:min(5,end)),'model')
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

