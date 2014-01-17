function compile(model)
%creates C source files modelspec.h and modelspec.c for
%all user-defined functions and model parameters (as globals)
%checks for errors and inconsistencies
%compiles

%TODO: (double) to all id ist it age and other int in exec strings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard conversion function for executable strings + other functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function resstr=StdConvertN(execstr)
        %returns the same executable string with correct reference to standard entities
        %Adds string of arguments to equation and parameter refs for C source 
        resstr=execstr;
        %min max
        resstr=regexprep(resstr,'\<min\>','MIN');
        resstr=regexprep(resstr,'\<max\>','MAX');
        %it - alwais attribute it to current period!
        resstr=regexprep(resstr,'\<it\>','curr->it');
        %age
        resstr=regexprep(resstr,'\<age\>','(curr->it+t0)');
        %decisions
        resstr=regexprep(resstr,['\<id\>'],'curr->id');
        for iii=1:numel(model.d)
            resstr=regexprep(resstr,['\<dc' sprintf('%d',iii) '\>'],['(byval>0?curr->dc[' sprintf('%d',iii-1) ']:decisions[curr->id+' sprintf('%d',iii-1) '*nd])']);
%            resstr=regexprep(resstr,['\<dc' sprintf('%d',iii) '\>'],['decisions[curr->id+' sprintf('%d',iii-1) '*nd]']);
        end
        %states
        resstr=regexprep(resstr,['\<ist\>'],'curr->ist');
        resstr=regexprep(resstr,['\<ist1\>'],'next->ist');
        for iii=1:numel(model.s)
            resstr=regexprep(resstr,['\<st' sprintf('%d',iii) '\>'],['(byval>0?curr->st[' sprintf('%d',iii-1) ']:states[curr->ist+' sprintf('%d',iii-1) '*nst])']);
            resstr=regexprep(resstr,['\<st' sprintf('%d',iii) 'n\>'],['(byval>0?next->st[' sprintf('%d',iii-1) ']:states[next->ist+' sprintf('%d',iii-1) '*nst])']);
%            resstr=regexprep(resstr,['\<st' sprintf('%d',iii) '\>'],['states[curr->ist+' sprintf('%d',iii-1) '*nst]']);
%            resstr=regexprep(resstr,['\<st' sprintf('%d',iii) 'n\>'],['states[next->ist+' sprintf('%d',iii-1) '*nst]']);
        end
        %equations
        if numel(model.eq)>0
        for tmp=model.eq
            switch tmp.type
            case 'current'
                resstr=regexprep(resstr,['\<' tmp.ref '\>'],[tmp.ref '(curr)']);
            case 'next'
                resstr=regexprep(resstr,['\<' tmp.ref '\>'],[tmp.ref '(curr,next)']);
            end
        end
        end
        %equation arguments - as used by periods
        resstr=regexprep(resstr,'\<cash\>','curr->cash');
        resstr=regexprep(resstr,'\<savings\>','next->savings');
        resstr=regexprep(resstr,'\<shock\>','next->shock');
        %mu and sigma
        resstr=regexprep(resstr,['\<sigma\>'],['sigma_param(curr,next)']);
        resstr=regexprep(resstr,['\<mu\>'],['mu_param(curr,next)']);
        %discount
        resstr=regexprep(resstr,'\<discount\>','discount(curr)');
        %survival
        resstr=regexprep(resstr,'\<survival\>','survival(curr)');
        %check for ist
        if sum([model.s.continuous])>0 && ~isempty(regexp(execstr,'\<ist\>'))
            warning 'The model contains continuous state variables, and yet `ist` index is used in exec strings!';
        end
    end
    function res=ProhibitString(inputstring,bannedstring,where)
        %error if forbidden string is found
        if ~iscell(bannedstring)
            if ~isempty(regexpi(inputstring,['\<' bannedstring '\d*\>']))
                fprintf('String: %s\n',inputstring);
                error ('ERROR:bannedstring','Error in %s: use of `%s` is not allowed in %s!\n',where,bannedstring,where);
            else
                res=inputstring;
            end
        else
            for bstr=bannedstring
                if ~isempty(regexpi(inputstring,['\<' bstr{1} '\d*\>']))
                    fprintf('String: %s\n',inputstring);
                    error ('ERROR:bannedstring','Error in %s: use of `%s` is not allowed in %s!\n',where,bstr{1},where);
                else
                    res=inputstring;
                end
            end
        end
    end
    
    function resstr=StdConvert_DELETETHISFUNCTION(execstr)
        %returns the same executable string with correct reference to standard entities
        %variables in C (st1-stN d1-dM)
        %Adds string of arguments to equation and parameter refs for C source 
        resstr=execstr;
        %min max
        resstr=regexprep(resstr,'\<min\>','MIN');
        resstr=regexprep(resstr,'\<max\>','MAX');
        %age
        resstr=regexprep(resstr,'\<age\>','(it+t0)');
        %decisions
        for iii=1:numel(model.d)
            resstr=regexprep(resstr,['\<dc' sprintf('%d',iii) '\>'],['decisions[id+' sprintf('%d',iii-1) '*nd]']);
            resstr=regexprep(resstr,['\<dc' sprintf('%d',iii) 'n\>'],['decisions[id1+' sprintf('%d',iii-1) '*nd]']);
        end
        %states
        for iii=1:numel(model.s)
            resstr=regexprep(resstr,['\<st' sprintf('%d',iii) '\>'],['states[ist+' sprintf('%d',iii-1) '*nst]']);
            resstr=regexprep(resstr,['\<st' sprintf('%d',iii) 'n\>'],['states[ist1+' sprintf('%d',iii-1) '*nst]']);
        end
        %equations
        if numel(model.eq)>0
        for tmp=model.eq
            switch tmp.type
            case 'current'
                resstr=regexprep(resstr,['\<' tmp.ref '\>'],[tmp.ref '(it,ist,cash)']);
            case 'next'
                resstr=regexprep(resstr,['\<' tmp.ref '\>'],[tmp.ref '(it,ist,id,ist1,savings,shock)']); %not it+1 because it's done in solver itself
            end
        end
        end
        %mu and sigma
        resstr=regexprep(resstr,['\<sigma\>'],['sigma_param(it,ist,id,ist1,savings)']);
        resstr=regexprep(resstr,['\<mu\>'],['mu_param(it,ist,id,ist1,savings)']);
        %check for ist
        if sum([model.s.continuous])>0 && ~isempty(regexp(execstr,'\<ist\>'))
            warning 'The model contains continuous state variables, and yet `ist` index is used in exec strings!';
        end
    end
    function wl(nr,str)
        %writes next line to code(nr) 
        persistent linenumber;
        if nr==0
            linenumber=ones(1,10); %10 separate counters
            return;
        end
        model.code(nr).lines{linenumber(nr)}=str;
        linenumber(nr)=linenumber(nr)+1;
    end
    function WriteFile(indx)
        %writes stored C code to file
        if ~exist(model.dir,'dir')
            if ~model.quiet
                display(['Creating directory ' model.dir]);
            end
            mkdir(model.dir);
        end
        if exist(strcat(model.dir,filesep,model.code(indx).file),'file')
            if ~model.quiet
                display(['Overwriting file ' model.code(indx).file]);
            end
        end
        if isunix
         %LINUX
         fid=fopen (strcat(model.dir,filesep,model.code(indx).file),'w+');
        else
         %WINDOWS   
         fid=fopen (strcat(model.dir,filesep,model.code(indx).file),'wt');
        end
        for str=model.code(indx).lines
%            fprintf(fid,[str{1} '\n']);
            fprintf(fid,'%s\n',str{1});
        end
        fclose(fid);
    end
    function resstr=SpcConvert(execstr)
        resstr=execstr;
        %states in initial period
        for iii=1:numel(model.s)
            resstr=regexprep(resstr,['\<st' sprintf('%d',iii) '\>'],['state0[' sprintf('%d',iii-1) ']']);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C code lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~model.quiet
 fprintf('Writing .c and .h files:\n');
 tic;
end
%clear
model.code=struct('file',{},'lines',{});
wl(0,'init');
%modelspec.c code file for all specific routines
model.code(1).file='modelspec.c';
%prototypes, coefarrays as constants
model.code(2).file='modelspec.h';
%Compose
wl(1,['/*Model specific code for the model ''' model.label '''*/']);
wl(1,'#include "egdst_lib.h"');
wl(1,'');
wl(2,['/*Model specific h for the model ''' model.label '''*/']);
wl(2,'#ifndef MODELSPECguard');
wl(2,'#define  MODELSPECguard');
wl(2,'');
%structures for standard variables
wl(2,'/*Main type for withing period variables to be used for both current and next period. Not all fields are used to either interpretation.*/');
wl(2,['typedef struct curr_variables {int it; int ist; double st[' sprintf('%d',numel(model.s)) ']; int id; double dc[' sprintf('%d',numel(model.d)) ']; double cash; double savings; double shock;} PeriodVars;']);
wl(2,'');
%number of equations
wl(2,'/*Number of equations*/');
wl(2,['#define NREQ ' sprintf('%d',numel(model.eq))]);
wl(2,'');
%Coefficients
if numel(model.coef)>0
wl(1,'/*In coefs increase dimentionality for matlab style base1 indexing of the coefficient arrays*/');
for coef1=model.coef
    wl(1,['const double ' coef1.ref '[' sprintf('%d',size(coef1.array,1)+1) '][' sprintf('%d',size(coef1.array,2)+1) '] = {']);
    tmpstr=regexprep(sprintf(' %2.1f#', zeros(1,size(coef1.array,2)+1)),'#',',');
    wl(1,['{' tmpstr(2:end-1) '},']);
    for i=1:size(coef1.array,1)
        for j=1:size(coef1.array,2)
            if j==1
                tmpstr=['{0.0, ' sprintf('%18.15f',coef1.array(i,j))];
            else
                tmpstr=[tmpstr ', ' sprintf('%18.15f',coef1.array(i,j))];
            end
        end
        if i==size(coef1.array,1)
            wl(1,[tmpstr '}};']);
        else
            wl(1,[tmpstr '},']);
        end
    end
end
end
wl(1,'');
%parameters
if numel(model.param)>0                          
    for i=1:numel(model.param)
        wl(2,['extern double ' model.param(i).ref '; /*Parameter:' model.param(i).description '*/']);
        wl(1,['double ' model.param(i).ref '; /*Parameter:' model.param(i).description '*/']);
    end
end
%State grids (for continuous states)
for i=1:numel(model.s)
    if model.s(i).continuous
        wl(2,['extern double *st' sprintf('%d',i) 'grid;']);
        wl(1,['double *st' sprintf('%d',i) 'grid;']);
    end
end
wl(2,['extern double *stgrids[' sprintf('%d',numel(model.s)) ']; /*pointers to grids of continous states*/']);
wl(1,['double *stgrids[' sprintf('%d',numel(model.s)) '];']);
wl(1,'');
%FUNCTIONS START HERE
wl(1,'void loadcontinuousgrid() {');
for i=1:numel(model.s)
    if model.s(i).continuous
        %wl(2,['const double *st' sprintf('%d',i) 'grid;']);
        wl(1,['st' sprintf('%d',i) 'grid = (double *) mxGetPr(mxGetField(mxGetProperty(Model,0,"s"),' sprintf('%d',i-1) ',"grid"));']);
        wl(1,['stgrids[' sprintf('%d',i-1) '] = (double *) st' sprintf('%d',i) 'grid;']);
    else
        wl(1,['/* stgrids[' sprintf('%d',i-1) '] is never used */']);
    end
end
wl(1,'}');
wl(1,'');
wl(2,'');
wl(2,'void loadcontinuousgrid();');
%discount
if numel(model.discount)==0
    error ('egdstmodel:method','Missing .discount, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
    wl(1,'double discount(PeriodVars *curr) {');
    wl(2,'double discount(PeriodVars *curr);');
    wl(1,'static double cash=-INFINITY; /*for compatability with eq type=current*/');
    wl(1,['return ' StdConvertN(ProhibitString(model.discount,{'id','dc','cash'},'discount')) ';}']);    
    wl(1,'');
end
%survival
wl(1,'double survival(PeriodVars *curr) {');
wl(2,'double survival(PeriodVars *curr);');
wl(1,'static double cash=-INFINITY; /*for compatability with eq type=current*/');
wl(1,['return ' StdConvertN(ProhibitString(model.survival,{'id','dc','cash'},'survival')) ';}']);
wl(1,'');
%utility
if numel(model.u)==0 || isempty(model.u.utility)
    error ('egdstmodel:method','Missing .u.utility, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
    wl(1,'double utility(PeriodVars *curr,double consumption) {');
    wl(2,'double utility(PeriodVars *curr,double consumption);');
    wl(1,'static double cash=-INFINITY;/*compatability with eq type=curren t*/');
    wl(1,'if (consumption<0) {printf("it=%d ist=%d id=%d consumption=%f ",curr->it,curr->ist,curr->id,consumption);mexWarnMsgTxt ("Utility function called with negative consumption..");}');
    if ischar(model.u.utility)
        wl(1,['return ' StdConvertN(ProhibitString(model.u.utility,'cash','utility')) ';}']);
    else
        for ln=reshape(model.u.utility,1,[])
            wl(1,StdConvertN(ProhibitString(ln{1},'cash','utility')));
        end
        wl(1,'}');
    end
    wl(1,'');     
end
%marginal utility
if numel(model.u)==0 || isempty(model.u.marginal)
    error ('egdstmodel:method','Missing .u.marginal, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
    wl(1,'double utility_marginal(PeriodVars *curr,double consumption) {');
    wl(2,'double utility_marginal(PeriodVars *curr,double consumption);');
    wl(1,'static double cash=-INFINITY;/*compatability with eq type=curren t*/');
        if ischar(model.u.marginal)
            wl(1,['return ' StdConvertN(ProhibitString(model.u.marginal,'cash','marginal utility')) ';}']);
        else
            for ln=reshape(model.u.marginal,1,[])
                wl(1,StdConvertN(ProhibitString(ln{1},'cash','marginal utility')));
            end
            wl(1,'}');
        end
    wl(1,'');     
end
%marginal utility inverse
if numel(model.u)==0 || isempty(model.u.marginalinverse)
    error ('egdstmodel:method','Missing .u.marginalinverse, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
    wl(1,'double utility_marginal_inverse(PeriodVars *curr,double mutility) {');
    wl(2,'double utility_marginal_inverse(PeriodVars *curr,double mutility);');
    wl(1,'static double cash=-INFINITY;/*compatability with eq type=current*/');
        if ischar(model.u.marginalinverse)
            wl(1,['return ' StdConvertN(ProhibitString(model.u.marginalinverse,'cash','marginal utility inverse')) ';}']);
        else
            for ln=reshape(model.u.marginalinverse,1,[])
                wl(1,StdConvertN(ProhibitString(ln{1},'cash','marginal utility inverse')));
            end
            wl(1,'}');
        end
    wl(1,'');     
end
%transformation function
%extrap is always defined (has a default)
    wl(1,'double tr(PeriodVars *curr,double x) {');
    wl(2,'double tr(PeriodVars *curr,double x);');
    wl(1,['return ' StdConvertN(ProhibitString(model.transform.direct,{'id','dc','cash','savings','shock'},'extrapolation function')) ';}']);
    wl(1,'');     
%extrap is always defined (has a default)
    wl(1,'double trinv(PeriodVars *curr,double x) {');
    wl(2,'double trinv(PeriodVars *curr,double x);');
    wl(1,['return ' StdConvertN(ProhibitString(model.transform.inverse,{'id','dc','cash','savings','shock'},'extrapolation function')) ';}']);
    wl(1,'');     
%cashinhand
if numel(model.budget)==0 || isempty(model.budget.cashinhand)
    error ('egdstmodel:method','Missing .budget.cashinhand, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
    wl(1,'double cashinhand(PeriodVars *curr,PeriodVars *next) {');
    wl(2,'double cashinhand(PeriodVars *curr,PeriodVars *next);');
    wl(1,'static double cash=-INFINITY; /*compatibility*/');
        if ischar(model.budget.cashinhand)
            wl(1,['return ' StdConvertN(ProhibitString(model.budget.cashinhand,'cash','cashinhand')) ';}']);
        else
            for ln=reshape(model.budget.cashinhand,1,[])
                wl(1,StdConvertN(ProhibitString(ln{1},'cash','cashinhand')));
            end
            wl(1,'}');
        end
    wl(1,'');     
end
%cashinhand_marginal
if numel(model.budget)==0 || isempty(model.budget.marginal)
    error ('egdstmodel:method','Missing .budget.marginal, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
    wl(1,'double cashinhand_marginal(PeriodVars *curr,PeriodVars *next) {');
    wl(2,'double cashinhand_marginal(PeriodVars *curr,PeriodVars *next);');
    wl(1,'static double cash=-INFINITY; /*compatibility*/');
        if ischar(model.budget.marginal)
            wl(1,['return ' StdConvertN(ProhibitString(model.budget.marginal,'cash','cashinhand marginal')) ';}']);
        else
            for ln=reshape(model.u.budget.marginal,1,[])
                wl(1,StdConvertN(ProhibitString(ln{1},'cash','cashinhand')));
            end
            wl(1,'}');
        end
    wl(1,'');     
end
%mu_param
if numel(model.shock)==0 || isempty(model.shock.mu)
    error ('egdstmodel:method','Missing .shock.mu, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
    wl(1,'double mu_param(PeriodVars *curr,PeriodVars *next) {');
    wl(2,'double mu_param(PeriodVars *curr,PeriodVars *next);');
    wl(1,'static double shock=-INFINITY;/*compatibility with eq type=next*/');
    wl(1,'static double cash=-INFINITY; /*compatibility with eq type=current*/');
        if ischar(model.shock.mu)
            wl(1,['return ' StdConvertN(ProhibitString(model.shock.mu,'shock','mu paremeter')) ';}']);
        else
            for ln=reshape(model.shock.mu,1,[])
                wl(1,StdConvertN(ProhibitString(ln{1},'shock','mu paremeter')));
            end
            wl(1,'}');
        end
    wl(1,'');     
end
%sigma_param
if numel(model.shock)==0 || isempty(model.shock.sigma)
    error ('egdstmodel:method','Missing .shock.sigma, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
wl(1,'double sigma_param(PeriodVars *curr,PeriodVars *next) {');
wl(2,'double sigma_param(PeriodVars *curr,PeriodVars *next);');
wl(1,'static double shock=-INFINITY;/*compatibility with eq type=next*/');
wl(1,'static double cash=-INFINITY; /*compatibility with eq type=current*/');
    if ischar(model.shock.sigma)
        wl(1,['return ' StdConvertN(ProhibitString(model.shock.sigma,'shock','sigma parameter')) ';}']);
    else
        for ln=reshape(model.shock.sigma,1,[])
            wl(1,StdConvertN(ln{1},ProhibitString(model.shock.sigma,'shock','sigma parameter')) );
        end
        wl(1,'}');
    end
wl(1,'');     
end
%choiceset
if numel(model.choiceset)==0
    error ('egdstmodel:method','Missing .choiceset, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
    wl(1,'int inchoiceset(PeriodVars *curr) {');
    wl(2,'int inchoiceset(PeriodVars *curr);');
    if model.choiceset.defaultallow
        wl(1,'int res = 1; /*true*/');
    else
        wl(1,'int res = 0; /*false*/');
    end
    wl(1,'static double cash=-INFINITY; /*compatibility*/');
    for rule1=model.choiceset.rules
        if model.choiceset.defaultallow
            wl(1,['if (' StdConvertN(ProhibitString(rule1.condition,'cash','.choiceset')) ') res = 0; /*' rule1.description '*/']);
        else
            wl(1,['if (' StdConvertN(ProhibitString(rule1.condition,'cash','.choiceset')) ') res = 1; /*' rule1.description '*/']);
        end
    end    
    wl(1,'return res;}');
    wl(1,'');
end
%feasible
if numel(model.feasible)==0
    error ('egdstmodel:method','Missing .feasible, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
    wl(1,'int feasible(PeriodVars *curr) {');
    wl(2,'int feasible(PeriodVars *curr);');
    if model.feasible.defaultfeasible
        wl(1,'int res = 1; /*true*/');
    else
        wl(1,'int res = 0; /*false*/');
    end
    for rule1=model.feasible.rules
        if model.feasible.defaultfeasible
            wl(1,['if (' StdConvertN(ProhibitString(rule1.condition,{'id','dc','cash'},'.feasible')) ') res = 0; /*' rule1.description '*/']);
        else
            wl(1,['if (' StdConvertN(ProhibitString(rule1.condition,{'id','dc','cash'},'.feasible')) ') res = 1; /*' rule1.description '*/']);
        end
    end    
    wl(1,'return res;}');
    wl(1,'');
end
%equations
if numel(model.eq)>0
    for eq1=model.eq
        switch eq1.type
        case 'current'
            wl(1,['double ' eq1.ref '(PeriodVars *curr) {']);
            wl(2,['double ' eq1.ref '(PeriodVars *curr);']);
        case 'next'
            wl(1,['double ' eq1.ref '(PeriodVars *curr,PeriodVars *next) {']);
            wl(2,['double ' eq1.ref '(PeriodVars *curr,PeriodVars *next);']);
        end
        if ischar(eq1.expression)
            wl(1,['return ' StdConvertN(eq1.expression) ';}']);
        else
            for ln=reshape(eq1.expression,1,[])
                wl(1,StdConvertN(ln{1}));
            end
            wl(1,'}');
        end
        wl(1,'');
    end
end
%parameters
wl(1,'void loadparameters () {');
wl(2,'void loadparameters ();');
for i=1:numel(model.param)
    wl(1,[model.param(i).ref '=mxGetScalar(mxGetField(mxGetProperty(Model,0,"param"),' sprintf('%d',i-1) ',"value"));']);
end
wl(1,'}');
wl(1,'');
%Transition probability
if numel(model.trpr)==0
    error ('egdstmodel:method','Missing .trpr, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
else
    wl(1,'double trpr(PeriodVars *curr,PeriodVars *next,int all) {');
    wl(2,'double trpr(PeriodVars *curr,PeriodVars *next,int all);');
    wl(1,'/*if all==0 only calculated trpr for discrete variables (sim)*/');
    wl(1,'double nval,res=1.0;');
    wl(1,'int varindex, varindex1;');
    for trpr1=model.trpr
        if numel(trpr1.varindex)==0
            error ('egdstmodel:method','Missing .trpr, can not proceed with compile!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');        
        else
            wl(1,['/*' model.s(trpr1.varindex).name '*/']);
            wl(1,['varindex =(curr->ist/(int)stm[nnst+' sprintf('%d',trpr1.varindex-1) '])%(int)stm[' sprintf('%d',trpr1.varindex-1) '];']);
            wl(1,['varindex1=(next->ist/(int)stm[nnst+' sprintf('%d',trpr1.varindex-1) '])%(int)stm[' sprintf('%d',trpr1.varindex-1) '];']);
%                                           varindex=(int) fmod((int) floor(ist/(int)stm[nnst+var]),(int)stm[var]);
%            wl(1,['printf("var %d: ist=%d varindex=%d : ist1=%d varindex1=%d\n",' sprintf('%d',trpr1.varindex) ',curr->ist,varindex,next->ist,varindex1);']);
            fstflag=true;
            for case1=trpr1.cases
                probs=case1.prob;
                if fstflag
                    wl(1,['if (' StdConvertN(case1.condition) ') {']);
                    fstflag=false;
                else
                    wl(1,['else if (' StdConvertN(case1.condition) ') {']);
                end
                if model.s(trpr1.varindex).discrete
                    %discrete state variable : transition prob matrix
                    wl(1,'  switch (varindex) {');
                    for ii=1:numel(model.s(trpr1.varindex).values)
                        wl(1,['  case ' sprintf('%d',ii-1) ':']);
                        wl(1,'	    switch (varindex1) {');
                        for jj=1:numel(model.s(trpr1.varindex).values)
                            wl(1,['      case ' sprintf('%d',jj-1) ':']);
                            wl(1,['          res*=' StdConvertN(case1.prob{ii,jj}) ';']);
                            wl(1, '          break;');
                        end                
                        wl(1,'   	default:');
                        wl(1,['		    mexErrMsgTxt("Error in trpr: unknown index of the next period state variable ''' model.s(trpr1.varindex).name ''' under condition (' case1.condition ')! See modelspec.c and compile.m");']);
                        wl(1,'           break;');
                        wl(1,'	    }');
                        wl(1,'       break;');
                    end
                    wl(1,'	default:');
                    wl(1,['		mexErrMsgTxt("Error in trpr: unknown index of the current period state variable ''' model.s(trpr1.varindex).name ''' under condition (' case1.condition ')! See modelspec.c and compile.m");']);
                    wl(1,'       break;');
                    wl(1,'	}');
%wl(1,['if (res!=0.0 && res!=1.0) printf("  * st' sprintf('%d',trpr1.varindex) ' d trpr=%f\n",res);']);
                    wl(1,'}');
                else
                    %continuous state variable : deterministic motion rule
                    %return linear interpolation coefficients
                    wl(1,'if (all==1) {');
                    wl(1,'    //reuse varindex since this is continuous variable');
                    wl(1,['    nval=' StdConvertN(ProhibitString(case1.prob,'ist1','motion rules')) ';']);
                    wl(1,['    varindex = bxsearch(nval,(double*)st' sprintf('%d',trpr1.varindex) 'grid,(int)stm[' sprintf('%d',trpr1.varindex-1) ']);']);
                    wl(1,['    if (varindex==varindex1) res*=(st' sprintf('%d',trpr1.varindex) 'grid[varindex+1]-nval)/(st' sprintf('%d',trpr1.varindex) 'grid[varindex+1]-st' sprintf('%d',trpr1.varindex) 'grid[varindex]);']);
                    wl(1,['    else if (varindex+1==varindex1) res*=(nval-st' sprintf('%d',trpr1.varindex) 'grid[varindex])/(st' sprintf('%d',trpr1.varindex) 'grid[varindex+1]-st' sprintf('%d',trpr1.varindex) 'grid[varindex]);']);
                    wl(1,'    else return 0.0;');
                    wl(1,'}');
%wl(1,['if (res!=0.0 && res!=1.0) printf("  * st' sprintf('%d',trpr1.varindex) ' c trpr=%f\n",res);']);
                    wl(1,'}');
                end
            end
            wl(1,['else {']);
            wl(1,['printf("---!!!--- period it=%d ist=%d id=%d ist1=%d",curr->it,curr->ist,curr->id,next->ist);']);
            wl(1,['mexErrMsgTxt("Error in trpr: unknown combination of current state and decision (the set of cases is not complete)!");']);
            wl(1,['}']);
%            wl(1,['printf("var %d: updated res=%1.4f\n",' sprintf('%d',trpr1.varindex) ',res);']);
            wl(1,['if (res==0.0) return 0.0;//early exit if result is already zero']);
        end
    end
    wl(1,'return res;}');
    wl(1,'');
end
%Special functions for sim
wl(1,'void trpr_cont(PeriodVars *curr,PeriodVars *next) {');
wl(2,'void trpr_cont(PeriodVars *curr,PeriodVars *next);');
wl(1,'/*fills the values of continuous states into *next*/');
wl(1,'int varindex;');
for trpr1=model.trpr
    if model.s(trpr1.varindex).continuous
        wl(1,['/*' model.s(trpr1.varindex).name '*/']);
        fstflag=true;
        for case1=trpr1.cases
            if fstflag
                wl(1,['if (' StdConvertN(case1.condition) ') {']);
                fstflag=false;
            else
                wl(1,['else if (' StdConvertN(case1.condition) ') {']);
            end
            wl(1,[' next->st[' sprintf('%d',trpr1.varindex-1) ']=' StdConvertN(case1.prob) ';']);
            wl(1,['}']);
        end
    end
end
wl(1,'}');
wl(1,'');



%{
for trpr1=model.trpr
    if model.s(trpr1.varindex).discrete
        wl(1,['/*' model.s(trpr1.varindex).name '*/']);
        wl(1,['varindex =(ist/(int)stm[nnst+' sprintf('%d',trpr1.varindex-1) '])%(int)stm[' sprintf('%d',trpr1.varindex-1) '];']);
        wl(1,['varindex1=(ist1/(int)stm[nnst+' sprintf('%d',trpr1.varindex-1) '])%(int)stm[' sprintf('%d',trpr1.varindex-1) '];']);
%        wl(1,['varindex =(int) fmod((int) floor( ist/(int)stm[nnst+' sprintf('%d',trpr1.varindex-1) ']),(int)stm[' sprintf('%d',trpr1.varindex-1) ']);']);
%        wl(1,['varindex1=(int) fmod((int) floor(ist1/(int)stm[nnst+' sprintf('%d',trpr1.varindex-1) ']),(int)stm[' sprintf('%d',trpr1.varindex-1) ']);']);
%                                      varindex=(int) fmod((int) floor(ist/(int)stm[nnst+var]),(int)stm[var]);
        %wl(1,['printf("var=%d : ist=%d varindex=%d : ist1=%d varindex=%d\n",' sprintf('%d',trpr1.varindex-1) ',ist,varindex,ist1,varindex1);']);
        fstflag=true;
        for case1=trpr1.cases
            probs=case1.prob;
            if fstflag
                wl(1,['if (' StdConvert(case1.condition) ') {']);
                fstflag=false;
            else
                wl(1,['else if (' StdConvert(case1.condition) ') {']);
            end
            %discrete state variable : transition prob matrix
            wl(1,'  switch (varindex) {');
            for ii=1:numel(model.s(trpr1.varindex).values)
                wl(1,['  case ' sprintf('%d',ii-1) ':']);
                wl(1,'	    switch (varindex1) {');
                for jj=1:numel(model.s(trpr1.varindex).values)
                    wl(1,['      case ' sprintf('%d',jj-1) ':']);
                    wl(1,['          res*=' StdConvert(case1.prob{ii,jj}) ';']);
                    wl(1, '          break;');
                end                
                wl(1,'   	default:');
                wl(1,['		    mexErrMsgTxt("Error in trpr: unknown index of the next period state variable ''' model.s(trpr1.varindex).name ''' under condition (' case1.condition ')! See modelspec.c and compile.m");']);
                wl(1,'           break;');
                wl(1,'	    }');
                wl(1,'       break;');
            end
            wl(1,'	default:');
            wl(1,['		mexErrMsgTxt("Error in trpr: unknown index of the current period state variable ''' model.s(trpr1.varindex).name ''' under condition (' case1.condition ')! See modelspec.c and compile.m");']);
            wl(1,'       break;');
            wl(1,'	}');
            wl(1,'}');
        end
        wl(1,['else {']);
        wl(1,['printf("---!!!--- period it=%d ist=%d id=%d ist1=%d",it,ist,id,ist1);']);
        wl(1,['mexErrMsgTxt("Error in trpr: unknown combination of current state and decision (the set of cases is not complete)!");']);
        wl(1,['}']);
    end
end
wl(1,'return res;}');
wl(1,'');
%}

%Calculation of equations (for sims)
wl(1,'void eqs_sim(PeriodVars *curr,PeriodVars *next,double *out) {');
wl(2,'void eqs_sim(PeriodVars *curr,PeriodVars *next,double *out);');
wl(1,'int i=0;');
model.simlabels = model.simlabels_default;
for ii=1:numel(model.s)
    model.simlabels{end+1}=[sprintf('%2d ',11+ii) model.s(ii).name ' (st' sprintf('%d',ii) ')'];
end
for ii=1:numel(model.d)
    model.simlabels{end+1}=[sprintf('%2d ',11+numel(model.s)+ii) model.d(ii).name ' (dc' sprintf('%d',ii) ')'];
end
for ii=1:numel(model.eq)
    if strcmp(model.eq(ii).type,'next')
        wl(1,['if (next==NULL) out[i++]=mxGetNaN(); /*' model.eq(ii).description ', undefined in first period */' ]);
        wl(1,['else            out[i++]=' model.eq(ii).ref '(curr,next); /*' model.eq(ii).description '*/' ]);
    else
        wl(1,['out[i++]=' model.eq(ii).ref '(curr); /*' model.eq(ii).description '*/' ]);
    end
    model.simlabels{end+1}=[sprintf('%2d ',11+numel(model.s)+numel(model.d)+ii) model.eq(ii).description ' (eq' sprintf('%d',ii) ')'];
end
wl(1,'}');
wl(1,'');

%Close .h file guard
wl(2,'');
wl(2,'#endif');
wl(2,'');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteFile(1);
WriteFile(2);
if ~model.quiet
 t=toc;
 fprintf('Done in %s\n',model.ht(t));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze the optimization properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 check if marginal utility is dependent on d
if iscell(model.u.marginal)
    tmp=[model.u.marginal{:}];%combine all lines    
else
    tmp=model.u.marginal;
end
if numel(regexp(tmp,['dc[' num2str([1:model.nnd]) ']']))>0 || numel(regexp(tmp,'id'))>0
    model.optim.optim_MUnoD=false;
else
    model.optim.optim_MUnoD=true;
end
% 2 check if utility is dependent on d
if iscell(model.u.utility)
    tmp=[model.u.utility{:}];%combine all lines    
else
    tmp=model.u.utility;
end
if numel(regexp(tmp,['dc[' num2str([1:model.nnd]) ']']))>0 || numel(regexp(tmp,'id'))>0
    model.optim.optim_UnoD=false;
else
    model.optim.optim_UnoD=true;
end
% 3 check if utility is separable in d and st
if iscell(model.u.utility)
    %multiline utility
    tmpcell=model.u.utility;
else
    %one line utility function
    tmpcell={model.u.utility};
end
model.optim.optim_UasD=true;
for tmp=tmpcell
    tmp=tmp{1};%uncell
    tmp0=tmp;%save
    if ~model.optim.optim_UasD
        break
    end
    %first replace all + inside brackets with #
    while numel(regexp(tmp,'\([^+\-()]*\)'))>0
       tmp=regexprep(tmp,'\(([^+\-()]*)\)','[$1]');
    end
    while numel(regexp(tmp,'\([^()]*(+|\-)[^()]*\)'))>0
       tmp=regexprep(tmp,'(\([^()]*)(+|\-)([^()]*\))','$1#$3');
    end
    if ~strcmp(tmp,tmp0) && ~model.quiet
        fprintf ('Utility contains brackes with +/- inside. For the analysis of additive separability the original expression it was transformed to:\n%s\n',tmp);
    end
    %then split by +
    for substr=regexp(tmp,'+|-','split')
        if (numel(regexp(substr{1},['dc[' num2str([0:model.nnd-1]) ']']))>0 || ...
            numel(regexp(substr{1},'id'))>0) && ...
           numel(regexp(substr{1},'consumption'))>0
            model.optim.optim_UasD=false;
            break;
        end
    end
end
if ~model.optim.optim_UasD
    error 'Utility is not additively separable in consumption and discrete choices. This case is not yet implemented!';
end
if iscell(model.u.utility) && ~model.quiet
    warning 'The model has multi-line utility function, which is hard to analyse automatically. The best guess is that utility is additively separable in decisions and states, but check again manually! The solver is proceeding under the assumption of additive separability.';
end
% 4 trpr independent of shock
%collect all exec strings from trpr
tmp='';
for trpr1=model.trpr
 for case1=trpr1.cases
    if iscell(case1.prob)
      tmp=[tmp '#' strcat(case1.prob{:})];
    else
      tmp=[tmp '#' case1.prob];
    end
 end
end
if numel(regexp(tmp,'shock'))>0
    model.optim.optim_TRPRnoSH=false;
else
    model.optim.optim_TRPRnoSH=true;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
egdstmodel_dir=regexprep(which('egdstmodel'),'egdstmodel.m','','ignorecase');
save_dir=pwd;
cd (model.dir);%compile in model directory
switch model.shock.type
case 'lognormal'
    flags = ' -DDISTRIB=1';
case 'normal'
    flags = ' -DDISTRIB=2';
end
%preporcessor tags
for tag=fieldnames(model.cflags)'
    val0=getfield(model.cflags,tag{1});
    if isnumeric(val0) && val0==0
        val1='0';
    elseif isnumeric(val0) && val0~=0
        val1=sprintf('%1.0d',val0);
    else
        val1=val0;
    end
    flags = [flags ' -D' tag{1} '=' val1];
end
if ~model.quiet
 fprintf('Compiler flags: %s\n',flags);
end
%compile!
try
    % 1. egdst_solve
    runstr=['mex ' egdstmodel_dir 'egdst_solver.c ' egdstmodel_dir 'egdst_lib.c ' model.code(1).file ' -I./ -I' egdstmodel_dir flags];
    %fprintf('%s\n',runstr);
    if ~model.quiet
     fprintf('Compiling egdst_solver.c');
     tic;
    end
    eval (runstr);
    if ~model.quiet
     t=toc;
     fprintf(' done in %s\n',model.ht(t));
    end
    % 2. egdst_call
    runstr=['mex ' egdstmodel_dir 'egdst_call.c ' egdstmodel_dir 'egdst_lib.c ' model.code(1).file ' -I./ -I' egdstmodel_dir flags];
    %fprintf('%s\n',runstr);
    if ~model.quiet
     fprintf('Compiling egdst_call.c');
     tic;
    end
    eval (runstr);
    if ~model.quiet
     t=toc;
     fprintf(' done in %s\n',model.ht(t));
    end 
    % 3. egdst_sim
    runstr=['mex ' egdstmodel_dir 'egdst_simulator.c ' egdstmodel_dir 'egdst_lib.c ' model.code(1).file ' -I./ -I' egdstmodel_dir flags];
    %fprintf('%s\n',runstr);
    if ~model.quiet
     fprintf('Compiling egdst_simulator.c');
     tic;
    end
    eval (runstr);
    if ~model.quiet
     t=toc;
     fprintf(' done in %s\n',model.ht(t));
    end 
catch exception
    cd (save_dir);
    error 'Error(s) while compiling C code!!!'    
end
cd (save_dir);
%clear to require recompile
model.M={};
model.D={};
model.sims=[];
%mark successful compile
model.needtocompile=false;
end

