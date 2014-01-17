function disp(m,varagin)
% overloaded display method for egdstmodel class objects

fprintf (['egdstmodel model labled "' m.label '"\n']);

%time
fprintf ('%-12s%s%1.0f','Time:','[',m.t0);
for i=m.t0:10:m.T
	fprintf ('-');
end
fprintf('%1.0f] nt=%d\n',m.T,m.nt);

%states
fprintf('%-12s','States:');
for i=m.s
	if i.discrete
		fprintf ('[%1.0f] ',numel(i.values));
	else
		fprintf ('[%1.0f''] ',numel(i.values));
	end
end
fprintf('nst=%d\n',m.nst);
for i=m.s
	if i.discrete
		fprintf ('%12s %d %s\n','',i.index,i.name);
	else
		fprintf ('%12s %d''continuous %s\n','',i.index,i.name);
	end
end

%decisions
fprintf('%-12s','Decisions:');
for i=m.d
	if i.discrete
		fprintf ('[%1.0f] ',numel(i.values));
	else
		fprintf ('[%1.0f''] ',numel(i.values));
	end
end
fprintf('nd=%d\n',m.nd);
for i=m.d
	if i.discrete
		fprintf ('%12s %d %s\n','',i.index,i.name);
	else
		fprintf ('%12s %d''continuous %s\n','',i.index,i.name);
	end
end

%wealth grid
fprintf ('%-12s%1.1f..%1.1f\n','Wealth:',m.a0,m.mmax);
fprintf ('%-12sngridm=%1.0f points in grid \n','',m.ngridm);

%shocks
fprintf('%-12s%s\n','Shocks:',m.shock.type);
if isnumeric(m.shock.mu)
	fprintf ('%-12smu=%1.3f\n','',m.shock.mu);
else
	fprintf ('%-12smu=%s\n','',m.shock.mu);
end
if isnumeric(m.shock.sigma)
	fprintf ('%-12ssigma=%1.3f\n','',m.shock.sigma);
else
	fprintf ('%-12ssigma=%s\n','',m.shock.sigma);
end

%info
fprintf('%-12s','Info:');
if m.needtocompile
	fprintf('model not compiled\n');
else
	fprintf('model compiled\n');
end
if ~isempty(m.M)
	fprintf('%-12s%s\n','','model solved');
else
	fprintf('%-12s%s\n','','model not solved');
end
if ~isempty(m.lastrun_solver)
	fprintf('%-12s%s%s\n','','last solver took ',m.ht(m.lastrun_solver));
end
if ~isempty(m.sims)
	if size(m.sims,1)>1
		fprintf('%-12s%1.0f%s\n','',size(m.sims,1),' simulations done');
	else
		fprintf('%-12s%1.0f%s\n','',size(m.sims,1),' simulation done');
	end
end
if m.quiet
	fprintf('%-12s%s\n','','messages suppressed');
end

end %function