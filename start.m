%choice of model script
fprintf('Navigate to one of folders:\n');
mn={};
for d=dir'
    if d.isdir && ~ismember(d.name(1),{'.','@'})
        mn{end+1}=d.name;
        fprintf('[%d] %s\n',size(mn,2),mn{1,end});
    end
end
ch=input('Where to go? ','s');
if isempty(ch)
    return
end
try
	cd(mn{str2num(ch)});
	dn=cd();
	fprintf('Current directory: %s\n',dn);
	run('start.m');
catch err
end
clear mn d ch err
