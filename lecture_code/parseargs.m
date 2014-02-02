% This function parses parameters passed in varargin in Matlab style {..,'name',value,..}
% checking them against default values in structure p.<names> for existance and type
% The values passed through varargin are assinged in the caller's workspace

function parseargs(p,up,varargin)

%run test on static workspace
usestruct=false;
try
 assignin('caller','tmp',true);
catch err
 % Catch static workspace error
 if strcmp(err.identifier,'MATLAB:err_static_workspace_violation')
  usestruct=true;
 else
  error 'FATAL ERROR in parseargs(): failed to assign to function workspace'
 end
end
 
%assign default values form p in caller's workspace
pfields=fieldnames(p);
for i=1:numel(pfields)
    if ~usestruct
        assignin('caller',pfields{i},p.(pfields{i}));
    end    
    pp.(pfields{i})=p.(pfields{i});
end

%input from cell up
for i=1:2:numel(up)-1
    %param name and value
    nm=up{i};
    vl=up{i+1};
    %recognize param name
    if ~ismember(nm,pfields)
        error ('%s is unknown parameter!',nm);
    else
        if ~isa(vl,class(p.(nm)))
            error ('Wrong value for parameter %s, expecting class=%s!',nm,class(p.(nm)));
        else
            if usestruct
                p.(nm)=vl;
            else
                assignin('caller',nm,vl);
            end
            pp.(nm)=vl;
        end
    end
end 
if usestruct
    assignin('caller','p',p);
end

%produce output unless 'quiet' is last argument in up or passed directly in varargin
if ~( nargin>2 && ischar(varargin{1}) && ismember(varargin{1},{'q','quiet','silent'}) ) && ...
   ~( ~isempty(up) && (ischar(up{end}) && ismember(up{end},{'q','quiet','silent'})) )
   

    %printf parameters
    fprintf('----------------------------------------\n');
    fprintf('          Parameters                    \n');
    fprintf('----------------------------------------\n');
    for nm=fieldnames(pp)'
        nm=nm{1};
        switch class(pp.(nm))
        case 'char'
            fprintf(' %-18s %15s\n',nm,pp.(nm));
        case {'integer','int8','int16','int64'}
            fprintf(' %-18s %15.0d\n',nm,pp.(nm));
        case {'numeric','double','single','float'}
            if numel(pp.(nm))==0
                fprintf(' %-18s Empty\n');
            else
                if numel(pp.(nm))==1
                    fprintf(' %-18s %15.5f\n',nm,pp.(nm));
                else
                    if numel(pp.(nm))==2
                        fprintf(' %-18s [%5.3f,%5.3f]\n',nm,pp.(nm)(1),pp.(nm)(2));
                    else
                        fprintf(' %-18s [%1.1f..%5.1f]',nm,pp.(nm)(1),pp.(nm)(end));
                        fprintf(' (%dx%d)\n',size(pp.(nm)));
                    end
                end            
            end
        case 'cell'
                fprintf(' %-18s cell (%dx%d)\n',nm,size(pp.(nm)));
        case 'struct'
                fprintf(' %-18s structure (%dx%d)\n',nm,size(pp.(nm)));
        case 'function_handle'
                fprintf(' %-18s function handle (%dx%d)\n',nm,size(pp.(nm)));
        otherwise
            fprintf(' %-18s %15g\n',nm,pp.(nm));
        end
    end
    fprintf('----------------------------------------\n');

end

end
