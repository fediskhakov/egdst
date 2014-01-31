classdef (ConstructOnLoad) egdstmodel < handle
%
% ==================================================================================
% Endogeneous grid method with discrete choices and additional states 
% ==================================================================================
% 'egdstmodel' is the class implementing all basic functionality for setting-up,
% solving and simulating of a discrete time dynamic programming model with 
% a continuous and discrete decisions.
% The model is solved with the endogenous grid point methods (Carroll) adopted for
% additional discrete choice and additional (discrete or discretized) state variables
% Written by: Fedor Iskhakov fediskhakov@gmail.com
% Reference on the method: IskhakovRustSchjerning 2012
%
% ==================================================================================
% RULES for defining model parts
% ==================================================================================
% Constructor of model objects
% ---------------------------------------------------------------------------------- 
%   Syntax: <model-object> = egdstmodel (['model-label']);
%   Parameter: label = model label (optional)
%
% States
% ---------------------------------------------------------------------------------- 
% Discrete states
% Syntax: state variable label + list of values with descriptions
%   <model>.s={'label',{i,'description',j,'description',..}}
% Continuous states
% Syntax: state variable label + grid limits + number of grid points
%   <model>.s={'label',[x0 x1],N}
%
% Decisions (only discrete!)
% ---------------------------------------------------------------------------------- 
% Input: decision variable label + list of values with descriptions
%   <model>.s={'label',{i,'description',j,'description',..}}
%
% General rules and notation for executable strings
% ---------------------------------------------------------------------------------- 
% Executable strings are used when generating .c and .h files by .compile
% Main rules: > executable expression in C syntax
%             > all refs must contain min 4 characters and be unique,
%             > refs must not repeat any of the constant names or reserved words
%               (run model.reswords to see the list) 
%             > multiline executable strings should be written in proper C code
%               (including ; return etc.) whereas single strings are automatically
%               formatted
% Terms for model entities for use in all executable strings:    
%           it      = time index (between 0 and T-t0)
%           age     = age at period t (age=it+t0)
%           st1-stN = state variables at period t
%           ist     = state index at period t
%                 !!! Index is accordning to the order of declaration !!!
%           dc1-dcN = decision variables for discrete decisions in period t
%           id      = index of decision at period t-1 (applying for period t)
%                 !!! Index is accordning to the order of declaration !!!
%         ref[i][j] = coefficient arrays (both indexes necessary, base1)
%               ref = parameter current values (by their own names)
%               ref = equations (by their own names)
% Where next period entities are required:
%         st1n-stNn = state variable at period t+1 (next period)
%           ist1    = state index at period t+1
%                 !!! Index is accordning to the order of declaration !!!
%         dc1n-dcNn = decision variables for discrete decisions in period t+1
%           id1     = index of decision at period t+1
%                 !!! Index is accordning to the order of declaration !!!
%
% Utility function (u)
% ---------------------------------------------------------------------------------- 
% Utility must be specified together with its derivative w.r.t. consumption
% and the inverse of this derivative (marginal and marginal-inverse).
%   Syntax: <model>.u={'utility|marginal|marginalinverse','<exec string>'}
%           <model>.u={'extrap','<exec string with x as single parameter>'}
%            Arguments:       it = period t (between 0 and T-t0)
%                            ist = period t state index
%                       st1..stN = state variables at period t
%                             id = discrete choice at period t
%                       dc1..dcN = decision variables at period t
%                    consumption = consumption at period t
%                       mutility = marginal utility (only for the inverse)
%
% Utility shape preserving transformation (transform) [optional]
% ---------------------------------------------------------------------------------- 
% To increase the accuracy of solution grids should be more dense where value
% function is "more curvy". To control grid points generating steps and also the
% extrapolation of the value function, one can define the utility-shape-preserving
% transformation and its inverse as the function of parameters, time and state,
% and x as single argument, for example for example log(x) for log utility or
% pow(x,crra_coef) for CRRA utility.
% Default transformation funciton is log(x+1) and its inverse exp(x)-1.
%   Syntax: <model>.transform={'<exec string: transform>','<exec string: inverse>'}
%            Arguments:       it = period t (between 0 and T-t0)
%                            ist = period t state index
%                       st1..stN = state variables at period t
%                              x = argument for transformation function
%
% Discount factor (discount) and survival probability (survival)
% ---------------------------------------------------------------------------------- 
% Discount factor is the function of current period (it,ist) to be multiplied
% with the expected value at next period. Executable string with standard arguments.
% Discount could depend on id, but this slows down calculation when marginal utility
% is independent of id, and complicates analysis. Use additional state variable.
% Survival probability is mainly needed for simulation model, but can be used within
% the model as well. By default survival=1.
%   Syntax: <model>.discount='<exec string>'
%   	    <model>.survival='<exec string>'
%            Arguments:       it = period t (between 0 and T-t0)
%                            ist = period t state index
%                       st1..stN = state variables at period t
%
% Intertemporal budget constraint (budget)
% ---------------------------------------------------------------------------------- 
% Budget constraint function has to express the next period cash-in-hand conditional
% on current state and decision, end of period liquid assets (savings), next period 
% state and income shock. Derivative w.r.t. savings must also be provided.
%   Syntax: <model>.budget={'cashinhand|marginal','<exec string>'};
%            Arguments:       it = period t (between 0 and T-t0)  NB: refers to t !!!
%                            ist = period t state index
%                       st1..stN = state variables at period t
%                             id = period t decision index
%                       dc1..dcN = decision variables in period t
%                           ist1 = next period state index               
%                     st1n..stNn = state variables at period t+1
%                        savings = end of period t assets (savings)
%                          shock = income shock at period t+1
%                             mu = location parameter of shock distribution
%                          sigma = shape parameter of shock distribution
% 
% Income shock (shock)
% ---------------------------------------------------------------------------------- 
% Income shock specifications include the type (normal or lognormal) and parameters
% mu and sigma. Default shock type is lognormal. Parameters mu and sigma are given 
% with executable string(s) and can depend on standard arguments (it,ist,id) 
% refering to previous period t and similarly to intertemporal budget constraint 
% on additional arguments refering to the period t+1 when shock takes place.
% It is beneficial for accuracy to specify dependence on the primitives of the
% model in parameters of income shocks, other than in budget constraint.
% When sigma==0 (as well as <0) income is treated as deterministic and in this 
% case 'shock' refers to the expectation, which is calculated analytically thus 
% saving time: all cases when shock is absent should be documented using
% this feature.
%   Syntax: <model>.shock='lognormal|normal'
%           <model>.shock={'mu|sigma','exec string'}
%            Arguments:       it = period t (between 0 and T-t0)  NB: refers to t !!!
%                            ist = period t state index
%                       st1..stN = state variables at period t
%                             id = period t decision index
%                       dc1..dcN = decision variables in period t
%                           ist1 = next period state index               
%                     st1n..stNn = state variables at period t+1
%                        savings = end of period t assets (savings)
%                             mu = location parameter of shock distribution
%                          sigma = shape parameter of shock distribution
%
% Transition probabilities (trpr)
% ---------------------------------------------------------------------------------- 
% Transition probabilities matrixes for each discrete state and deterministic 
% motion rules for continuous states. Transition probability matrixes may be both
% numeric of cell arrays of executable strings. For each variable there may be
% defined several matrices or motion rules which are invoked if the 'case' condition
% is satisfied (giving flexibility of any joint distribution). Empty strings are
% treated as 0.
%   Syntax for discrete additional state variables:
%   <model>.trpr={'condition',[square trpr matrix],..} OR
%   <model>.trpr={'condition',{square trpr matrix of exec strings},..} OR
%   <model>.trpr={varindex,'condition',[square trpr matrix],..} OR
%   <model>.trpr={varindex,'condition',{square trpr matrix of exec strings},..}
%   Syntax for continuous additional state variables:
%               NOTE: only deterministic transitions can be defined this way!
%   <model>.trpr={'condition','exec string for next period value'} OR
%   <model>.trpr={varindex,'condition','exec string for next period value'}
%            Arguments:       it = period t (between 0 and T-t0)  NB: refers to t !!!
%                            ist = period t state index
%                       st1..stN = state variables at period t
%                             id = period t decision index
%                       dc1..dcN = decision variables in period t
%                           ist1 = next period state index               
%                     st1n..stNn = state variables at period t+1
%                        savings = end of period t assets (savings)
%                             mu = location parameter of shock distribution
%                          sigma = shape parameter of shock distribution
%                          shock = income shock at period t+1
%
% Choicesets (choiceset)
% ---------------------------------------------------------------------------------- 
% Choiceset conditions help to single out cases where particular discrete decisions
% are feasible/infeasible in contrast to all other. The default behavior is given
% with model.choiceset('default',true|false), and each specific case defined with
% model.choiceset('executable string, false when ==0','description') specify when
% to apply the opposite.
%   Syntax: <model>.choiceset={'defaultallow',true|false|1|0}
%           <model>.choiceset={'condition','description'}
%            Arguments:       it = period t (between 0 and T-t0)
%                            ist = period t state index
%                       st1..stN = state variables at period t
%                             id = the index of the current discrete decision to be
%                                  checked whether it belongs to the choiceset at 
%                                  time t
%                       dc1..dcN = same, but decision values
%
% Feasibility conditions (feasible)
% ---------------------------------------------------------------------------------- 
% Similarly to Choiceset Feasibility conditions indicate whether a particular 
% (it,ist) combination is feasible, otherwise calculations for this point 
% of the state space are skipped. 
%   Syntax: <model>.feasible={'defaultfeasible',true|false|1|0}
%           <model>.feasible={'condition','description'}
%            Arguments:       it = period t (between 0 and T-t0)
%                            ist = period t state index
%                       st1..stN = state variables at period t
% 
% Equations (eq)
% ---------------------------------------------------------------------------------- 
% Equations are user defined sub-routines that return double and describe either
% current or next period. The first kind has type=current and should take standard
% input arguments (it,ist) for period t refering to the current period.
% The second kind has type=next and should take standard input arguments 
% (it,id,ist1) refering to the next period. Here time script 'it' should still
% be treated as refering to previous period (t-1). Default equation type is current.
% Chaining equations is perfectly ok, but no type=next equation can call
% type=current equation and vise versa. Only chains of the same type are allowed.
%   Syntax: <model>.eq={'ref','description','string','current|next'} OR
%           <model>.eq={'ref','description',{'string1' 'string2' ..},'current|next'}
% Available arguments: 1. for type=current equations only:
%                             it = period t (between 0 and T-t0)
%                            ist = period t state index
%                       st1..stN = state variables at period t
%                             id = period t decision index
%                       dc1..dcN = decision variables in period t
%                           cash = cash-in-hand in current period (unless
%                                  equation is refered to in trpr or budget)
%                      2. for type=next equations only:
%                             it = period t (between 0 and T-t0)  NB: refers to t !!!
%                            ist = period t state index
%                       st1..stN = state variables at period t
%                             id = period t decision index
%                       dc1..dcN = decision variables in period t
%                           ist1 = next period state index               
%                     st1n..stNn = state variables at period t+1
%                        savings = end of period t assets (savings)
%                          shock = income shock at period t+1
%                             mu = location parameter of shock distribution
%                          sigma = shape parameter of shock distribution
%
% Coefficient arrays (coef)
% ---------------------------------------------------------------------------------- 
% Coefficients are constant two dimentional arrays, it is therefore necessary to 
% index all coefficients with double subscripts [i][j]. Matlab style base1 indexing is
% preserved in C code.
%   Syntax: <model>.coef={'ref','description',[x11 x12 x13 .. ; x21 x22 x23 .. ; ..]}
%
% Parameters (param)
% ---------------------------------------------------------------------------------- 
% Parameters return current value of named parameters. Use .setparam method for 
% changing parameter values. Use .getparam method for checking parameter values.
%   Syntax: <model>.param={'ref','description',initial value}
%
%
%
% ==================================================================================
% Inputs and outputs of the class methods
% ==================================================================================
% .quietly
% ---------------------------------------------------------------------------------- 
%   Input:  none OR off
%   Output: none
% Suppresses all output OR turn back on default messages
%
% .compile
% ---------------------------------------------------------------------------------- 
%   Input:  none
%   Output: none
% Creates the object directory, .c and .h files in it, and compiles all the C code.
% Must be run before most other methods
%
% .setparam
% ---------------------------------------------------------------------------------- 
%   Input:  none or [param vector] or 'name',value,.. or index,value,..
%   Output: none
% Lists current parameter values, sets parameters by vector or by name
%
% .getparam
% ---------------------------------------------------------------------------------- 
%   Input:  none or 'name' or index
%   Output: param value
% Lists current parameter values, outputs parameter values
% 
% .solve
% ---------------------------------------------------------------------------------- 
%   Input:  none
%   Output: M = (grid,consumption,savings,value)
%                  !!! For technical reasons first value of VF is set to evf(a0)
%                      and should be replaced by -inf in interpretations !!!
%           D = (decisions,thresholds)
% Solves the model using the current values of parameters.
% Reference to method: IskhakovRustSchjerning 2011.
%
% .call
% ---------------------------------------------------------------------------------- 
%   Input:  funcname = function to call
%           funcar   = function arguments
%   Output: output from the function called
%   Supported functions and their arguments:
%  'utility'|'util'|'u' (it <t0 to T>, ist <1 to .nst>, id <1 to .nd>,consumption)
%  'mutility'|'mu'      (it <t0 to T>, ist <1 to .nst>, id <1 to .nd>,consumption)
%  'discount'|'df'      (it <t0 to T>, ist <1 to .nst>)
%  'budget'|'b' (it <t0 to T>, ist <1 to .nst>, id <1 to .nd>,savings,ist(t+1),shock)
%  'mbudget'|'mb' (it <t0 to T>, ist <1 to .nst>, id <1 to .nd>,savings,ist(t+1),shock)
%  'value'|'vf'         (it <t0 to T>, ist <1 to .nst>, cash)
% Runs internal functions of the model
% 
% .plot1 .plot2 .plot3 etc
% ---------------------------------------------------------------------------------- 
% Run <model>.plotN help
%
% .indx
% ---------------------------------------------------------------------------------- 
%   Input:  state/decision vector OR
%           set of state/decision vectors in a rows
%   Output: index/indeces
% Returns the index of the state or decision vector given the vector(s)
%
% .sim
% ---------------------------------------------------------------------------------- 
%   Input:  init AND/OR rndtype (in any combination)
%               init = [ist(t0) cash-at-hand(t0);
%			            ...	
%			            ist(t0) cash-at-hand(t0)] initial contitions in rows
%               NB: only indexes are accepted as initial conditions, i.e.
%                   for continuous state variables only grid points can 
%                   serve as initial values
%               rndtype = 'same_shocks|own_shocks' whether to use the same random 
%                          numbers for each simulation, or to use one stream of 
%                          random number for all simulations 
%   Output: sims
%           simlabels
% Simulates the model for given initial conditions and using a given stream of 
% random numbers. Init conditions should contain indexes of discrete decision
% and state describing period t=0 in rows. When many rows are given, many
% simulations are performed. 
% When model contains continuous state variables, their values are simulated 
% exactly according to the motion rules, but values of value function, optimal
% consumption and eqations are interpolated across the continues state. If optimal
% discrete decision is ambiguous, it is simulated using interpolation weights as
% probabilities. Initial values of the continuous state must belong to the grid.
% The sequence of U[0,1] random variables used for simulation should be sufficiently
% long. If randstream is not yet saved in the object, it is loaded from 
% randstream.mat in the current directory, if unsuccessful, generated with rand().
% The output is three dimentional array by people-time-simulated variables
% and the labels for those variables.
% Method .sims2panel provides easy way to generate panel of simulated data.
% Method .print_init prints initial conditions along with indexes in state space.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Properties available from outside (public)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access=public)
        %Properties that define egdst model
        label   = 'no label'        %label to mark the model
        t0      = NaN               %start period
        T       = NaN               %end period
        s     = struct('index',{}, ...
                       'name',{}, ...
                       'type',{}, ... 
                       'discrete',{}, ...
                       'continuous',{}, ...
                       'values',struct('value',{},'description',{}), ...
                       'gridlimits',{}, ...
                       'gridpoints',{}, ...
                       'grid',{})   %states variables (in structures)
        d     = struct('index',{}, ...
                       'name',{}, ...
                       'type',{}, ... 
                       'discrete',{}, ...
                       'continuous',{}, ...
                       'values',struct('value',{},'description',{}), ...
                       'gridlimits',{}, ...
                       'gridpoints',{}, ...
                       'grid',{})   %decision variables (in structures)
        mmax    = NaN               %maximum cash-at-hand
        ngridm  = 10                %standard number of endogenious grid points
        ngridmax= 100               %maximum number of endogenious grid points
        nthrhmax= 100               %maximum number of thresholds
        ny      = 1                 %number of points in discrete representation of random income
        a0      = 0                 %credit constraint
        discount= ''                %Discount factor (executable string)
        survival= '1.0'             %Survival probability (executable string)
        u     = struct('utility',{},...
                       'marginal',{},...
                       'marginalinverse',{},...
                       'extrap',{});
                                    %Utility function (executable string)
        transform =struct('direct','log(x+1)',...
                          'inverse','exp(x)-1'); 
                                    %transform function
        budget =struct('cashinhand',{},...
                       'marginal',{})
                                    %Intertemporal budget equation (executable string/strings)
        shock = struct('type','lognormal',...
                       'mu',{},...
                       'sigma',{})%Idiosyncratic shock distribution, mu, sigma parameters (executable string) 
        trpr  = struct('varindex',{}, ...
                       'cases',[struct('condition','','prob',{{}})])
                                    %Transition probabilites: matrixes of executable strings for each variable under cases   
        choiceset = struct('defaultallow',true,...
                           'rules',[struct('condition',{},'description',{})])
                                    %Choice sets D(s) (conditions are executable strings ==true iff <>0)  
        feasible = struct('defaultfeasible',true,...
                           'rules',[struct('condition',{},'description',{})])
                                    %Feasibility conditions for states (conditions are executable strings ==true iff <>0)  
        eq    = struct('ref',{},...
                       'type',{},...
                       'expression',{},...
                       'description',{}) %equations: intermediate entities (executable strings)
        coef  = struct('ref',{},...
                       'array',{},...
                       'description',{}) %coefficient arrays (matrixes)               
        param = struct('ref',{},...
                       'description',{},...
                       'value',{})  %Parameters
        cflags = struct(...         %preprocessor flags to be passed to compiler
            'TOLERANCE','1e-10',...                   
            'ZEROCONSUMPTION','1e-10',...                   
            'DOUBLEPOINT_DELTA','1e-10',...                   
            'VERBOSE','0')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Properties set inside
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess=private, GetAccess=public)
        %Properties that are calculated within the object, but readable
        nd      = 0                 %dimentionality of decision space
        nnd     = 0                 %number of decision variables (components of the decision vector)
        nst     = 0                 %dimentionality of state space
        nnst    = 0                 %number of state variables (components of the state vector)
        stm     = []                %vector of numbers of values in state vars with multiplicators
        dm      = []                %vector of numbers of values in decision vars
        code    = struct('file',{},'lines',{}) %lines of C code stored by compile method
        M       = {}                %2-dimentional cell array with M,C,A,VF (returned by .solve)
        D       = {}                %2-dimentional cell array with D,TH (returned by .solve)
        states  = []                %list of state vars
        decisions = []              %list of decision vars
        sims    = []                %3-dimentional cell array with results of the simulations
        simlabels = {}              %labels for simulated variables
        init    = []                %initial conditions for simulations
        randstream = []             %stream of random numbers
        optim = struct(...          %solution optimization properties of the model
            'optim_MUnoD',false,... %true when marginal utility does not depend on d
            'optim_UnoD',false,...  %true when utility does not depend on d
            'optim_UasD',false,...  %true when utility in separable in d
            'optim_TRPRnoSH',false ...   %true when transition probabilities are independent of shocks
                      )
        quadrature = []             %column-vector of quadrature weights and column-vector of quadrature abscissas
        quiet   = false             %suppress all output
    end 
    properties (SetAccess=private, GetAccess=public, Dependent)
        %Dependent: calculated on retrieval
        nt
    end
    properties (Access=private, Hidden)
        %Hidden: system use, not shown to user
        id                          %unique id of the model
        dir     = ''                %directory for C files
        caller  = ''                %who created the object
        needtocompile = true        %marks whether the model had been changed and needs to be compiled (for the first time or again)
        globals = { ...             %int globals
            't0','T',...
            'ngridm','ngridmax','nthrhmax','ny',...
            'nd','nnd','nst','nnst',...
            'mmax','a0'}
        reservedwords = {...        %other reserved words (inside C code)
            'age','ist','id1','id',...
            'st1','st2','st3','st4','st5','st6','st7','st8','st9','st10','st11','st12','st13','st14','st15',...
            'st1grid','st2grid','st3grid','st4grid','st5grid','st6grid','st7grid','st8grid','st9grid','st10grid','st11grid','st12grid','st13grid','st14grid','st15grid',...
            'dc1','dc2','dc3','dc4','dc5','dc6','dc7','dc8','dc9','dc10','dc11','dc12','dc13','dc14','dc15',...
            'consumption','mutility',...
            'discount','survival','utility','utility_marginal','utility_marginal_inverse',...
            'trpr','feasible','inchoiceset',...
            'cashinhand','cashinhand_marginal',...
            'shock','mu_param','sigma_param','sigma','mu',...
            'cash','savings',...
            'min','max'}
        simlabels_default = {'1  cash-in-hand (M)','2  optimal consumption (C)','3  optimal saving (A)','4  value function','5  current discrete decision index (id)',...
                     '6  current period state index (ist)','7  mu parameter of shock distribution','8  sigma parameter of shock distribution','9  income shock',...
                     '10 utility','11 discount factor'} %first labels for simulated variables (amended by .compile)
        lastrun_solver = []         %last runtime for solver
    end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Constructor and get.property methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function model = egdstmodel(varargin)
            % Constructor for class egdstmodel
            % Inputs: [label] : optional
            model.id=newid();
            model.label='<no name>';
            if (nargin>0)
                model.label=varargin{1};
            end
            model.dir=['tmp_' regexprep(model.label(1:min(10,numel(model.label))),'\W','') '_' model.id];
            model.needtocompile=true;
        end
        function delete(obj)
            %distructor of object : clean up after yourself
            try 
            	%check if Matlab is in the object directory, and try to step out before deletion
            	if strcmp(cd,obj.dir)
            		cd('..');
            	end
                rmdir(obj.dir,'s');
                if ~obj.quiet
                	display(['Deleting directory ' obj.dir]);
            	end
            catch err
            	if ~obj.quiet
                	display(['Could not delete directory ' obj.dir]);
                end
            end
        end
        function res=get.nt(object)
            %number of time periods
            res=object.T-object.t0+1;
            if res<1
                error ('egdstmodel:property','Error: t0>T!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
            end
        end
        function res=set.cflags(obj,value)
            %change in precompiler flags
            obj.cflags=value;
            obj.needtocompile=true;
        end
        function res=set.ngridm(obj,value)
            obj.ngridm=1*value(1);
            if obj.ngridm*1.5>obj.ngridmax
                obj.ngridmax=2*obj.ngridm;    
                %warning (sprintf('Automatically adjusted ngridmax to %d becasue ngridm=%d',obj.ngridmax,obj.ngridm));
            end
        end
        function res=set.ngridmax(obj,value)
            obj.ngridmax=1*value(1);
            if obj.ngridm*1.5>obj.ngridmax
                obj.ngridmax=2*obj.ngridm;    
                warning (sprintf('Automatically adjusted ngridmax to %d becasue ngridm=%d',obj.ngridmax,obj.ngridm));
            end
        end
        %model parts
        function res=set.discount(obj,value)
            if isnumeric(value(1))
                if value(1)==0
                    obj.discount='0.0';
                else
                    obj.discount=sprintf('%25.25f',value(1));
                end
            elseif ischar(value)
                obj.discount=value;
            else
                error ('egdstmodel:property','Unrecognized value for discount! Need single number of executable string.\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
            end
        end
        function res=set.survival(obj,value)
            if isnumeric(value(1))
                if value(1)==0
                    obj.survival='0.0';
                else
                    obj.survival=sprintf('%25.25f',value(1));
                end
            elseif ischar(value)
                obj.survival=value;
            else
                error ('egdstmodel:property','Unrecognized value for survival! Need single number of executable string.\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
            end
        end
        function obj=set.s(obj,value)
            %adds element to state space structure
            obj.needtocompile=true;
            %if value == [] deletes everything
            if isempty(value)
                obj.s = struct('index',{}, ...
                               'name',{}, ...
                               'type',{}, ... 
                               'discrete',{}, ...
                               'continuous',{}, ...
                               'values',struct('value',{},'description',{}), ...
                               'gridlimits',{}, ...
                               'gridpoints',{}, ...
                               'grid',{});  %states variables (in structures)
                obj.nnst=0;
                obj.trpr = struct('varindex',{}, ...
                                  'cases',[struct('condition','','prob',{{}})]); %Clear trpr too!                
            else
                if iscell(value) && numel(value)==2 && ischar(value{1}) && iscell(value{2}) && mod(numel(value{2}),2)==0
                    %name + list of values and descriptions
                    %type discrete = list of values
                    obj.nnst=obj.nnst+1;
                    obj.s(end+1)=struct(   'index',obj.nnst, ...
                                           'name',value{1}, ...
                                           'type','discrete', ... 
                                           'discrete',true, ...
                                           'continuous',false, ...
                                           'values',struct('value',{},'description',{}), ...
                                           'gridlimits',[], ...
                                           'gridpoints',NaN, ...
                                           'grid',[]);
                    for i=1:numel(value{2})/2
                        if isnumeric(value{2}{2*i-1})
                            obj.s(end).values(i)=struct('value',value{2}{2*i-1},'description',value{2}{2*i});
                        else
                            error ('egdstmodel:property','Non-numeric value of the state variable detected!..\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                        end
                    end
                    obj.stm(obj.nnst)=numel(value{2})/2;
                elseif iscell(value) && numel(value)==2 && ischar(value{1}) && isnumeric(value{2})
                    %name + list of values (numeric vector)
                    %type discrete = list of values
                    obj.nnst=obj.nnst+1;
                    obj.s(end+1)=struct(   'index',obj.nnst, ...
                                           'name',value{1}, ...
                                           'type','discrete', ... 
                                           'discrete',true, ...
                                           'continuous',false, ...
                                           'values',struct('value',{},'description',{}), ...
                                           'gridlimits',[], ...
                                           'gridpoints',NaN, ...
                                           'grid',[]);
                    for i=1:numel(value{2})
                        obj.s(end).values(i)=struct('value',value{2}(i),'description',sprintf('value %1.3f',value{2}(i)));
                    end
                    obj.stm(obj.nnst)=numel(value{2});
                elseif iscell(value) && numel(value)==3 && ischar(value{1}) && isnumeric(value{2}) && numel(value{2})==2 && isnumeric(value{3})
                    %name + grid limits + number of gridpoints
                    %type continuous = grid, but tranclate into {values} as well
                    obj.nnst=obj.nnst+1;
                    obj.s(end+1)=struct(   'index',obj.nnst, ...
                                           'name',value{1}, ...
                                           'type','continuous', ... 
                                           'discrete',false, ...
                                           'continuous',true, ...
                                           'values',struct('value',{},'description',{}), ...
                                           'gridlimits',value{2}, ...
                                           'gridpoints',value{3}, ...
                                           'grid',linspace(value{2}(1),value{2}(2),value{3}) );
                    %rewrite grid points to values for trpr compile
                    for i=1:obj.s(end).gridpoints
                        obj.s(end).values(i)=struct('value',obj.s(end).grid(i),'description','grid point');
                    end
                    obj.stm(obj.nnst)=value{3};
                else
                 error ('egdstmodel:property','Unrecognized structure for state variable!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                end
                %update dimentions
                obj.stm=[obj.stm(1:obj.nnst) stepmult(obj.stm(1:obj.nnst))];                
                obj.nst=calculate_nst(obj); %dimention of the state space
                %build states                
                obj.states=buildstates(obj,'s');
            end
        end
        function obj=set.d(obj,value)
            %adds element to state space structure
            obj.needtocompile=true;
            %if value == [] deletes everything
            if isempty(value)
                obj.d = struct('index',{}, ...
                               'name',{}, ...
                               'type',{}, ... 
                               'discrete',{}, ...
                               'continuous',{}, ...
                               'values',struct('value',{},'description',{}), ...
                               'gridlimits',{}, ...
                               'gridpoints',{}, ...
                               'grid',{});  %decision variables (in structures)
            else
                if iscell(value) && numel(value)==2 && ischar(value{1}) && iscell(value{2}) && mod(numel(value{2}),2)==0
                    %name + list of values and descriptions
                    %type discrete = list of values
                    obj.nnd=obj.nnd+1;
                    obj.d(end+1)=struct(   'index',obj.nnd, ...
                                           'name',value{1}, ...
                                           'type','discrete', ... 
                                           'discrete',true, ...
                                           'continuous',false, ...
                                           'values',struct('value',{},'description',{}), ...
                                           'gridlimits',[], ...
                                           'gridpoints',NaN, ...
                                           'grid',[]);
                    for i=1:numel(value{2})/2
                        if isnumeric(value{2}{2*i-1})
                            obj.d(end).values(i)=struct('value',value{2}{2*i-1},'description',value{2}{2*i});
                        else
                            error ('egdstmodel:property','Non-numeric value of the state variable detected!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                        end
                    end
                    obj.dm(obj.nnd)=numel(value{2})/2;
                elseif iscell(value) && numel(value)==3 && ischar(value{1}) && isnumeric(value{2}) && numel(value{2})==2 && isnumeric(value{3})
                    %name + grid limits + number of gridpoints
                    %type continuous = grid
                    obj.nnd=obj.nnd+1;
                    obj.d(end+1)=struct(   'index',obj.nnd, ...
                                           'name',value{1}, ...
                                           'type','continuous', ... 
                                           'discrete',false, ...
                                           'continuous',true, ...
                                           'values',struct('value',{},'description',{}), ...
                                           'gridlimits',value{2}, ...
                                           'gridpoints',value{3}, ...
                                           'grid',value{2}(1):(value{2}(2)-value{2}(1))/(value{3}-1):value{2}(2));
                    obj.dm(obj.nnd)=value{3};
                else
                 error ('egdstmodel:property','Unrecognized structure for decision variable!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                end
                %update dimentions
                obj.dm=[obj.dm(1:obj.nnd) stepmult(obj.dm(1:obj.nnd))];                
                obj.nd=calculate_nd(obj); %dimention of the decision space
                %build states
                obj.decisions=buildstates(obj,'d');
            end
        end
        function obj=set.u(obj,value)
           %sets utility
           obj.needtocompile=true;
           %or deletes if value == []
            if isempty(value)
                obj.u = struct('utility',{},...
                               'marginal',{},...
                               'marginalinverse',{});
            else
                if iscell(value) && numel(value)==2 && (ischar(value{2}) || iscell(value{2}))
                    switch value{1}
                    case 'utility'
                        obj.u(1).utility = value{2};
                    case 'marginal'
                        obj.u(1).marginal = value{2};
                    case 'marginalinverse'
                        obj.u(1).marginalinverse = value{2};
                    case 'extrap'
                        obj.u(1).extrap = value{2};
                    otherwise
                        error ('egdstmodel:property','Unrecognized structure for utility definition!\nNeed {''utility|marginal|marginalinverse'',''executable string''} or {''utility|marginal|marginalinverse'',{''line1'',''line3'',''line3''}}!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                    end
                end
            end            
        end
        function obj=set.transform(obj,value)
           %sets transformation function
           obj.needtocompile=true;
           %or deletes if value == []
            if isempty(value)
                obj.transform =struct('direct','log(x+1)',...
                                      'inverse','exp(x)-1'); 
            else
                if iscell(value) && numel(value)==2 && ischar(value{1}) && ischar(value{2})
                    obj.transform.direct=value{1};
                    obj.transform.inverse=value{2};
                else
                    error ('egdstmodel:property','Unrecognized structure for transformation function!\nNeed {''executable string for direct transformation function'',''executable string for inverse transformation function''}!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                end
            end            
        end
        function obj=set.budget(obj,value)
           %sets intertemporal budget constraint
           obj.needtocompile=true;
           %or deletes if value == []
            if isempty(value)
                obj.budget = struct('cashinhand',{},...
                                    'marginal',{});
            else
                if iscell(value) && numel(value)==2 && (ischar(value{2}) || iscell(value{2}))
                    switch value{1}
                    case 'cashinhand'
                        obj.budget(1).cashinhand = value{2};
                    case 'marginal'
                        obj.budget(1).marginal = value{2};
                    otherwise
                        error ('egdstmodel:property','Unrecognized structure for budget definition!\nNeed {''cashinhand|marginal'',''executable string''} or {''cashinhand|marginal'',{''line1'',''line3'',''line3''}}!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');   
                    end
                end
            end            
        end
        function obj=set.shock(obj,value)
            %sets shock structure
            obj.needtocompile=true;
            %if value == [] deletes everything
            if isempty(value)
                obj.shock = struct('type','lognormal',...
                                   'mu',{},...
                                   'sigma',{});
            else
                if ischar(value) && ismember(value,{'lognormal','normal'})
                    obj.shock(1).type=value;
                elseif iscell(value) && numel(value)==2 && ismember(value{1},{'mu','sigma'}) ...
                        && (ischar(value{2}) || iscell(value{2}))
                    obj.shock(1).(value{1})=value{2};
                elseif iscell(value) && numel(value)==2 && ismember(value{1},{'mu','sigma'}) ...
                        && (isnumeric(value{2}))
                    if value{2}==0
                        obj.shock(1).(value{1})='0.0';
                    else                    
                        obj.shock(1).(value{1})=sprintf('%25.25f',value{2});
                    end
                else
                    error ('egdstmodel:property','Unrecognized structure for shock definition!\nNeed ''lognormal|normal'' or {''mu|sigma'',''executable string''} or {''mu|sigma'',{''line1'',''lin2'',''lin3''}}!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                end
            end        
        end
        function obj=set.eq(obj,value)
            %adds element to equations
            obj.needtocompile=true;
            %if value == [] deletes everything
            %sets type=current by default
            if isempty(value)
                obj.eq = struct('ref',{},...
                                'type',{},...
                                'expression',{},...
                                'description',{});
            else
                %update eq when called with existing name
                indx=-1;
                if iscell(value) && numel(value)>1 && ischar(value{1})
                    for i=1:numel(obj.eq)
                        if strcmp(obj.eq(i).ref,value{1})
                            warning(['Replacing equation ' sprintf('%d',i) ' in the model with the updated equation!']);
                            indx=i;
                            break
                        end
                    end
                    if indx==-1
                        indx=numel(obj.eq)+1;
                    end
                end
                if iscell(value) && numel(value)==3 && ischar(value{1}) && ischar(value{2}) && iscell(value{3}) && ischar(value{3}{1})
                    %name + description + cell array of lines
                    obj.eq(indx) = struct('ref',value{1},...
                                           'type','current',...
                                           'expression',{value{3}},...
                                           'description',value{2});
                elseif iscell(value) && numel(value)==3 && ischar(value{1}) && ischar(value{2}) && ischar(value{3})
                    %name + description + line
                    obj.eq(indx) = struct('ref',value{1},...
                                           'type','current',...
                                           'expression',value(3),...
                                           'description',value{2});
                elseif iscell(value) && numel(value)==4 && ischar(value{1}) && ischar(value{2}) && iscell(value{3}) && ischar(value{3}{1}) ...
                    && ischar(value{4}) && ismember(value{4},{'current','next'})
                    %name + description + cell array of lines + type
                    obj.eq(indx) = struct('ref',value{1},...
                                           'type',value{4},...
                                           'expression',{value{3}},...
                                           'description',value{2});
                elseif iscell(value) && numel(value)==4 && ischar(value{1}) && ischar(value{2}) && ischar(value{3}) ...
                    && ischar(value{4}) && ismember(value{4},{'current','next'})
                    %name + description + line + type
                    obj.eq(indx) = struct('ref',value{1},...
                                           'type',value{4},...
                                           'expression',value(3),...
                                           'description',value{2});
                else
                 error ('egdstmodel:property','Unrecognized structure for equation definition!\nNeed \n{''ref'',''description'',''executable string''} or \n{''ref'',''description'',{''line1'',''lin2'',''lin3''}} or \n{''ref'',''description'',''executable string'',''current|next''} or \n{''ref'',''description'',{''line1'',''lin2'',''lin3'',''current|next''}}!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                end
                checkrefs(obj);
            end
        end
        function obj=set.coef(obj,value)
            %adds element to coefficients
            obj.needtocompile=true;
            %if value == [] deletes everything
            if isempty(value)
                obj.coef = struct('ref',{},...
                               'array',{},...
                               'description',{});
            else
                if iscell(value) && numel(value)==3 && ischar(value{1}) && ischar(value{2}) && isnumeric(value{3})
                    %name + description + array of coefficients
                    obj.coef(end+1) = struct('ref',value{1},...
                                             'array',value{3},...
                                             'description',value{2});
                else
                 error ('egdstmodel:property','Unrecognized structure for coefficient definition!\nNeed {''ref'',''description'',[array of numbers]}!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                end
                checkrefs(obj);
            end
        end
        function obj=set.choiceset(obj,value)
            %adds rule to choiceset
            obj.needtocompile=true;
            %if value == [] deletes everything
            %to set default allow, pass {'defaultallow',true/false} <true by default>
            if isempty(value)
                obj.choiceset = struct('defaultallow',true,...
                                       'rules',[struct('condition',{},'description',{})]);
            else
                if iscell(value) && strcmp(value{1},'defaultallow') && (isnumeric(value{2}) || islogical(value{2}))
                    %defaultallow + num/logical
                    obj.choiceset(1).defaultallow = logical(value{2});
                elseif iscell(value) && numel(value)==2 && ischar(value{1}) && ischar(value{2})
                    %rule + description
                    obj.choiceset(1).rules(end+1) = struct('condition',value{1},...
                                                         'description',value{2});
                else
                 error ('egdstmodel:property','Unrecognized structure for choiceset definition!\nNeed {''rule'',''description''} or {''defaultallow'',true/false}!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                end
            end
        end
        function obj=set.feasible(obj,value)
            %adds rule to feasibility conditions
            obj.needtocompile=true;
            %if value == [] deletes everything
            %to set default feasible, pass {'defaultfeasible',true/false} <true by default>
            if isempty(value)
                obj.feasible = struct('defaultfeasible',true,...
                                       'rules',[struct('condition',{},'description',{})]);
            else
                if iscell(value) && strcmp(value{1},'defaultfeasible') && (isnumeric(value{2}) || islogical(value{2}))
                    %defaultallow + num/logical
                    obj.feasible(1).defaultfeasible = logical(value{2});
                elseif iscell(value) && numel(value)==2 && ischar(value{1}) && ischar(value{2})
                    %rule + description
                    obj.feasible(1).rules(end+1) = struct('condition',value{1},...
                                                        'description',value{2});
                else
                 error ('egdstmodel:property','Unrecognized structure for feasibility definition!\nNeed {''rule'',''description''} or {''defaultfeasible'',true/false}!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                end
            end
        end
        function obj=set.param(obj,value)
            %adds element to parameters
            %if value == [] deletes everything
            if isempty(value)
                obj.param = struct('ref',{},...
                                   'description',{},...
                                   'value',{});
            else
                if iscell(value) && numel(value)==3 && ischar(value{1}) && ischar(value{2}) && isnumeric(value{3})
                    %name + description + startvalue
                    obj.param(end+1) = struct('ref',value{1},...
                                              'description',value{2},...
                                              'value',value{3});
                elseif isstruct(value) && numel(value)==numel(obj.param)
                 %if obj.param(i).value=X is called, value is the full structure obj.param with value changed in one place
                 for i=1:numel(obj.param)
                    for fld=fields(obj.param)'
                        if (numel(obj.param(i).(fld{1}))~=numel(value(i).(fld{1})) || ...
                            sum(sum(obj.param(i).(fld{1})~=value(i).(fld{1})))==numel(obj.param(i).(fld{1}))) ...
                           && ~strcmp(fld{1},'value')
                            error ('egdstmodel:property','Only ''model.param(i).value=X'' are allowed with model.param property\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                        end
                    end   
                 end
                 obj.param=value;
                elseif isstruct(value) && numel(value)>numel(obj.param)
                 error ('egdstmodel:property','Index out of bounds in param assignment\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                else 
                 error ('egdstmodel:property','Unrecognized structure for parameter definition! Need {''ref'',''description'',start value}!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');
                end
                checkrefs(obj);
            end
        end
        function obj=set.trpr(obj,value)
            %adds element to transition probability matrix
            obj.needtocompile=true;
            %if value == [] deletes everything
            if isempty(value)
                obj.trpr = struct('varindex',{}, ...
                                  'cases',[struct('condition','','prob',{{}})]);
            else
                if iscell(value) && numel(value)==2 ...
                   && ischar(value{1}) ...
                   && ( ...
                    isnumeric(value{2}) && size(value{2},1)==obj.stm(numel(obj.s)) && size(value{2},2)==obj.stm(numel(obj.s)) ...
                   ) || ( ...
                    iscell(value{2}) && size(value{2},1)==obj.stm(numel(obj.s)) && size(value{2},2)==obj.stm(numel(obj.s)) ...
                   ) 
                   value={numel(obj.s),value{:}}; %when varindex is skipped, add to the last added variable
                end
                if iscell(value) && numel(value)==2 ...
                   && ischar(value{1}) && ischar(value{2}) 
                   value={numel(obj.s),value{:}}; %when varindex is skipped, add to the last added variable
                end
                if iscell(value) && numel(value)==3 ...
                   && isnumeric(value{1}) && value{1}>=1 && value{1}<=obj.nnst && floor(value{1})==value{1} ...
                   && ischar(value{2}) ...
                   && isnumeric(value{3}) && size(value{3},1)==obj.stm(value{1}) && size(value{3},2)==obj.stm(value{1})
                    %varindex + condition + matrix of transition probabilities
                    obj.trpr(value{1}).varindex=value{1}; %keep varindex for checking and convenience
                    %rewrite numbers as strings
                    obj.trpr(value{1}).cases(end+1)=struct('condition',value{2},'prob',{cellfun(@(x) sprintf('%10.10f',x),num2cell(value{3}),'UniformOutput',false)});
                elseif iscell(value) && numel(value)==3 ...
                   && isnumeric(value{1}) && value{1}>=1 && value{1}<=obj.nnst && floor(value{1})==value{1} ...
                   && ischar(value{2}) ...
                   && iscell(value{3}) && size(value{3},1)==obj.stm(value{1}) && size(value{3},2)==obj.stm(value{1})
                    %varindex + condition + matrix of executable strings (square, right size)
                    obj.trpr(value{1}).varindex=value{1};
                    %replace empty strings with zeros in value{3}
                    tmp=value{3};
                    for i=1:numel(value{3})
                        if (isempty(tmp{i}))
                            tmp{i}='0.0';
                        end
                    end
                    obj.trpr(value{1}).cases(end+1)=struct('condition',value{2},'prob',{tmp});
                elseif iscell(value) && numel(value)==3 ...
                   && isnumeric(value{1}) && value{1}>=1 && value{1}<=obj.nnst && floor(value{1})==value{1} ...
                   && ischar(value{2}) && ischar(value{3}) 
                    %varindex + condition + executable string as motion rule
                    obj.trpr(value{1}).varindex=value{1};
                    obj.trpr(value{1}).cases(end+1)=struct('condition',value{2},'prob',value{3});
                else
                 error ('egdstmodel:property','Unrecognized structure for trpr definition!See <a href="matlab:help egdstmodel">help egdstmodel</a>');
                end
            end
        end
    end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Public methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access=public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %suppress all ouput
        function quietly (model,varargin)
            model.quiet=true;
            if nargin>1
                model.quiet=false;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %show reserved words
        function reswords (model)
            display ([model.globals model.reservedwords]');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %print verbal description of the model
        function print (model)
            fprintf('\n');
            fprintf('egdst model: %s\n',model.label);
            fprintf('        dir: %s\n',model.dir);
            fprintf('   n(state space): %d\n',model.nst);
            fprintf('n(desicion space): %d\n',model.nd);
            fprintf('   n(grid points): %d\n', model.ngridm);
            fprintf('    n(shock grid): %d\n', model.ny);
            %states
            fprintf('___________________________________________________________________________\n');
            fprintf('%3s%10s%37s%8s%8s%8s\n','No','Value','Description','Min','Max','N');
            fprintf('___________________________________________________________________________\n');
            for el=model.s
                if el.discrete
                    fprintf('%3d%10s%37s%8s%8s%8d\n',el.index,'',el.name(1:min(end,37)),'','',numel(el.values));
                    for el1=el.values
                        fprintf('%3s%10.3f%37s%8s%8s%8s\n','',el1.value,el1.description(1:min(end,37)),'','','');
                    end
                elseif el.continuous
                    fprintf('%3d%10s%37s%8.3f%8.3f%8.3f\n',el.index,'',el.name(1:min(end,37)),el.gridlimits(1),el.gridlimits(2),el.gridpoints);
                end
            end    
            %decisions
            fprintf('___________________________________________________________________________\n');
            for el=model.d
                if el.discrete
                    fprintf('%3d%10s%37s%8s%8s%8d\n',el.index,'',el.name(1:min(end,37)),'','',numel(el.values));
                    for el1=el.values
                        fprintf('%3s%10.3f%37s%8s%8s%8s\n','',el1.value,el1.description(1:min(end,37)),'','','');
                    end
                elseif el.continuous
                    fprintf('%3d%10s%37s%8.3f%8.3f%8.3f\n',el.index,'',el.name(1:min(end,37)),el.gridlimits(1),el.gridlimits(2),el.gridpoints);
                end
            end    
            fprintf('___________________________________________________________________________\n');
            if model.needtocompile
                fprintf('Model needs to be compiled\n');
            end
            fprintf('\n');
            %DO MORE!!!!
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %create C source files and compile
        compile(model,varargin)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %set and list parameters
        function setparam(model,varargin)
            switch nargin-1
            case 0
                %no argument is passed ==> just list params
                listparam(model)
            case 1
                %one argument ==> treat as vector of param values
                if isnumeric(varargin{1}) && numel(varargin{1})==numel(model.param)
                    for k=1:numel(model.param)
                        model.param(k).value = varargin{1}(k);
                    end
                    %model.setparam();%display
                else
                    error ('egdstmodel:method','Passed vector does not match the dimentionality of param vector in the model!\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');     
                end
            otherwise
                %many arguments ==> treat as 'ref',value pairs or index,value pairs
                if mod(numel(varargin),2)==0
                    for k=1:numel(varargin)/2
                        if isnumeric(varargin{2*k-1}) && varargin{2*k-1}<=numel(model.param) && isnumeric(varargin{2*k})
                            model.param(varargin{2*k-1}).value=varargin{2*k};
                        elseif ischar(varargin{2*k-1}) && ismember(varargin{2*k-1},{model.param.ref}) && isnumeric(varargin{2*k})
                            for j=1:numel(model.param)
                                if strcmp(model.param(j).ref,varargin{2*k-1})==1
                                    model.param(j).value=varargin{2*k};
                                end
                            end
                        else
                            error ('egdstmodel:method','Unrecognized pairs ''name'',value,.. or index,value,.. !\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');         
                        end
                    end
                    %model.setparam;%display
                else
                    error ('egdstmodel:method','Expected pairs ''name'',value,.. or index,value,.. !\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');     
                end
            end  
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function res=getparam(model,varargin)
            switch nargin-1
            case 0
                %no argument is passed ==> just list params
                listparam(model);
                res=[model.param.value];
            case 1
                val=varargin{1};
                if isnumeric(val) && val<=numel(model.param)
                    res=model.param(val).value;
                elseif ischar(val) && ismember(val,{model.param.ref})
                    for j=1:numel(model.param)
                        if strcmp(model.param(j).ref,val)==1
                            res=model.param(j).value;
                        end
                    end
                else
                    error ('egdstmodel:method','Unrecognized parameter name or parameter index out of bounds');
                end
            otherwise
                %many arguments ==> could be vector of indexes, but this is not implemented
                error ('egdstmodel:method','Please call getparam with index or param ref');
            end  
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %solve the model
        function res=solve(model)
            if model.needtocompile
                error ('egdstmodel:method','The model needs to be compiled first!\nRun <model>.compile\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');     
            end
            if model.ngridmax<=model.ngridm
                model.ngridmax=2*model.ngridm;
            end
            if ~model.quiet
                c=clock;
                if isempty(model.lastrun_solver)
                    fprintf('Running the model solver.. (started at %1.0f:%02.0f:%02.0f) ',c(4:6));
                else
                    fprintf('Running the model solver.. (started at %1.0f:%02.0f:%02.0f, last time it took %s) ',c(4:6),model.ht(model.lastrun_solver));
                end
            end
            %find quadrature points
            if model.ny>1
                [qx,qw]= quadpoints(model.ny,0,1);
                model.quadrature=[qw qx];
            end
            savedir=pwd;
        try
            cd (model.dir);
        catch err
            error (sprintf('Could not find model directory %s',model.dir));
        end
            tic;
            %[model.M model.D]=egdst_solver(model);
%debug output            
[model.M model.D dbgout]=egdst_solver(model);
%evalin('base','clear dbgout');
%assignin('base','dbgout',dbgout);
            model.lastrun_solver=toc;
            if ~model.quiet
                fprintf(' done in %s\n',model.ht(model.lastrun_solver));
            end
            cd (savedir);
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %call internal parts of the the model
        function res=call(model,func,funcargs)
			if model.needtocompile || numel(model.M)==0
                error ('egdstmodel:method','The model needs to be compiled and solved first!\nRun <model>.compile\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');     
            end
            savedir=pwd;
            cd (model.dir);
            try
                switch func %translate words to numbers
                case {'utility','util','u'}
                    res=egdst_call(model,1,funcargs);
                case {'mutility','mu'}
                    res=egdst_call(model,2,funcargs);
                case {'discount','df'}
                    res=egdst_call(model,3,funcargs);
                case {'budget','b'}
                    res=egdst_call(model,4,funcargs);
                case {'mbudget','mb'}
                    res=egdst_call(model,5,funcargs);
                case {'value','vf'}
                    res=egdst_call(model,6,funcargs);
                otherwise
                    error ('egdstmodel:method','Unknown internal model function to call! See <a href="matlab:help egdstmodel">help egdstmodel</a>');     
                end
            catch er
            end
            cd (savedir);
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % simulate the model
        function sim(model,varargin)
            if model.needtocompile || numel(model.M)==0
                error ('egdstmodel:method','The model needs to be compiled and solved first!\nRun <model>.compile and <model>.sim\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');     
            end
            if (sum(sum(cellfun(@isempty,model.M)))>0 || sum(sum(cellfun(@isempty,model.D)))>0) && numel(model.feasible.rules)==0
                error ('egdstmodel:method','The model is missing solutions for some time period in <model>.M and/or <model>.D, and needs to be solved in full before simulations can be done!\nRun <model>.sim\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');     
            end
            %argument may be (init,randtype) or init or randtype
            rndtype=0; %default is to use the random stream for all simulations (each simulation has different random numbers)
            if nargin>1
                for i=1:nargin-1
                    if isnumeric(varargin{i}) && (size(varargin{i},2)==2)
                        %looks like init
                        model.init=varargin{i};
                    elseif ischar(varargin{i}) && ismember(varargin{i},{'own_shocks' 'same_shocks'})
                        %looks like rndtype
                        rndtype = 1*strcmp(varargin{i},'same_shocks');%1 for same shocks and 0 otherwise
                    else
                        error ('egdstmodel:method','Could not recognize argumend %d!\nSee <a href="matlab:help egdstmodel">help egdstmodel for discription of inputs.</a>',i);     
                    end
                end
            end
            %checks
            %init
            if isempty(model.init) 
                model.init=[1 0];
            end
            model.init(:,2)=max(model.init(:,2),model.a0);%correct default init for case a0>0
            %randstream
            if numel(model.randstream)>0
                randstrm=model.randstream;
                if ~model.quiet
                 fprintf('Using model.randstream\n');
                end
            elseif exist('randstream.mat','file')
                %use first variable in randstream.mat
                d1=load('randstream.mat');
                d2=fieldnames(d1);
                randstrm=getfield(d1,d2{1});
                model.randstream=randstrm;
                if ~model.quiet
                 fprintf('Loaded randstream from randstream.mat into the model object.\n');
                end
            else
                randstrm=rand(max(size(model.init,1),100)*model.nt*100,1);
                model.randstream=randstrm;
                if ~model.quiet
                 fprintf('Generated new randstream and saved it into the model object (but not on disk).\n');
                end 
            end
            %run sims
            savedir=pwd;
            cd (model.dir);
            if ~model.quiet
             fprintf('Running %d model simulation(s)..',size(model.init,1));
             tic;
            end 
            model.sims=[];
            model.sims=egdst_simulator(model,rndtype);
            %go from var-time-isim to isim-time-var orger of the dimensions
            model.sims=permute(model.sims,[3 2 1]);
            if ~model.quiet
             t=toc;
             fprintf(' done in %s\n',model.ht(t));
            end
            cd (savedir);
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sims in panel style panel dataset
        function [data labels]=sims2panel(model)
            if numel(model.sims)==0
                error ('egdstmodel:method','The model needs to be simulated first!\nRun <model>.sim\nSee also <a href="matlab:help egdstmodel">help egdstmodel</a>');     
            end
            niii=size(model.sims,1);
            nttt=size(model.sims,2);
            nvvv=size(model.sims,3);
            iii=ones(nttt,1)*[1:niii];
            ttt=[1:nttt]'*ones(1,niii);
            data=cat(3,iii,ttt,permute(model.sims,[2 1 3]));
            data=reshape(data,niii*nttt,nvvv+2);
            data(isnan(data(:,3)),:)=[];%drop lines with dead people
            labels={'i' 't' model.simlabels{:}}';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculated the index/indeces of state or decision vector
        function res=indx(model,v)
            if size(v,2)==numel(model.s)
                %passed v is states, unless
                if numel(model.s)==numel(model.d)
                    warning 'As the model has equal number of states and decisions, the passed vectors are treated as states'
                end
                mask=model.states;
            elseif size(v,2)==numel(model.d)
                mask=model.decisions;
            else
                error 'Wrong dimention of the state/decision vector'
            end
            %find the index using mask
            res=[];
            for i=1:size(v,1)
                f=find(sum((mask==ones(size(mask,1),1)*v(i,:)),2)==size(mask,2));
                if isempty(f)
                    res=[res;NaN];
                else
                    res=[res;f];
                end
            end            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % print init for simulations nicely
        function res=print_init(model)
            if isempty(model.init)
                warning 'No initial conditions in the model';
                return;
            end
            %init ist indexes
            if size(model.init,1)>10
                initl=model.init(1:10,1);%first 10 inits
                initw=model.init(1:10,2);
            else
                initl=model.init(:,1);
                initw=model.init(:,2);
            end
            %matrix of init st values
            ists=[];
            for i=1:numel(initl)
                ists=[ists, model.states(initl(i),:)'];
            end
            %print
            fprintf('%-28s%12s ','State variables','sims index:');
            fprintf('%6d ',[1:size(ists,2)]);
            fprintf('\n');
            fprintf('________________________________________');
            fprintf('%s',repmat('_______',1,size(ists,2)));
            fprintf('\n');
            fprintf('%40s ','initial wealth:');
            fprintf('%6.1f ',initw);
            fprintf('\n');
            fprintf('________________________________________');
            fprintf('%s',repmat('_______',1,size(ists,2)));
            fprintf('\n');
            for i=1:numel(model.s)
                fprintf('%-40s ',model.s(i).name(1:min(end,40)));
                fprintf('%6.2f ',ists(i,:));
                fprintf('\n');
            end
            fprintf('________________________________________');
            fprintf('%s',repmat('_______',1,size(ists,2)));
            fprintf('\n');
            fprintf('%40s ','state index:');
            fprintf('%6d ',initl);
            fprintf('\n');
            fprintf('________________________________________');
            fprintf('%s',repmat('_______',1,size(ists,2)));
            fprintf('\n');
        end        
    end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access=public, Hidden)
		function str = ht(model,x)
			% This function returns elapced time in hunam readable format
			% x= time in seconds 
			names={'s','m','h',' day',' week',' year'};
			plural={'','','','s','s','s'};
			dividers=[1,60,60*60,24*60*60,7*24*60*60,52*7*24*60*60];
			str='';
			for i=numel(names):-1:1
				if dividers(i)>1
					dig=floor(x/dividers(i));
					digstr=sprintf('%1.0f',dig);		
				else
					dig=x;
					digstr=sprintf('%1.3f',dig);		
				end
				if dig>1
					str=sprintf('%s %s%s%s',str,digstr,names{i},plural{i});
				elseif dig<=1 && dig>0
					str=sprintf('%s %s%s',str,digstr,names{i});
				end
				x=x-dig*dividers(i);
			end
		end
	end


end %of classdef
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Sub-functions (not part of the class, for this file use only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res=newid()
    %Returnes the sampe dpmodel with new pseudo unique ID
    res=tempname;
    res(1:end-36)=[]; %likely unique ID for the object
    res=regexprep(res,'_','');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res=calculate_nst(model)
    %calculate n
    res=1;
    for i=1:numel(model.s)
     if model.s(i).discrete
      res=res*numel(model.s(i).values);
     elseif model.s(i).continuous
      res=res*model.s(i).gridpoints;
     end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res=calculate_nd(model)
    %calculate n
    res=1;
    for i=1:numel(model.d)
     if model.d(i).discrete
      res=res*numel(model.d(i).values);
     elseif model.d(i).continuous
      res=res*model.d(i).gridpoints;
     end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res=stepmult(a)
    %given [a1 a2 .. aN] calculates [a2*a3*..aN a3*a4*..aN  ..  aN 1]
    for i=1:numel(a)
     res(i)=prod(a(i+1:end));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function states=buildstates(o,wh)
    %build matrix of state/decision vector values
    switch wh
    case 's' 
        upp=o.nnst;
        mm=o.stm;
    case 'd' 
        upp=o.nnd;
        mm=o.dm;
    end
    for i=upp:-1:1
        %collect values
        vals=[];
        for j=1:numel(o.(wh)(i).values)
         vals(j,1)=o.(wh)(i).values(j).value;
        end
        if i==upp
            states=vals;
        else
            states=[kron(vals,ones(mm(upp+i),1)) repmat(states,mm(i),1)];         
        end
    end        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok=checkrefs(model)
    %This function sends ok=1 if no problems with refs and ok=0 otherwise 
    %MORE CHECKS ARE NEEDED
    refs=[model.globals model.reservedwords];
    k=numel(refs)+1;
    if numel(model.eq)>0
        for tmp=model.eq
            refs{k}=tmp.ref;
            k=k+1;
        end
    end
    if numel(model.coef)>0
        for tmp=model.coef
            refs{k}=tmp.ref;
            k=k+1;
        end
    end
    if numel(model.param)>0
        for tmp=model.param
            refs{k}=tmp.ref;
            k=k+1;
        end
    end
    refs=sort(refs);
    mask=strcmp(refs(1:end-1),refs(2:end));
    if sum(strcmp(refs(1:end-1),refs(2:end)))>0
       fprintf('The following refs are not unique or use the names of globals/reserved words:\n');
       display(refs(mask)');   
       error ('egdstmodel:method','Not uniques refs!\nSee <a href="matlab:help egdstmodel">help egdstmodel</a>');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function adjp(ax)
    %this adjusts the position of axes to take full width of the figure
    tightinset=get(ax,'TightInset');
    pos=get(ax,'Position');
    set(ax,'Position',[tightinset(1)+0.03 pos(2) 1-tightinset(1)-tightinset(3)-0.06 pos(4)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quadpoints: n gaussian quadrature abscissae and weights over internal [x1,x2]
% by John Rust
function [x,w]= quadpoints(n,x1,x2)
    %eps=1.0e-8; %use matlab system eps when commented out
    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for i=1:m;
      z = cos(pi*(i-.25)/(n+.5));
      z1= 2;
      while (abs(z-z1)>eps) ;
        p1 = 1.0;
        p2 = 0.0;
        for j=1:n;
          p3 = p2;
          p2 = p1;
          p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
        end;
        pp = n*(z*p1-p2)/(z*z-1.0);
        z1 = z;
        z = z1-p1/pp;
      end;
      x(i,1) = xm-xl*z;
      x(n+1-i,1) = xm+xl*z;
      w(i,1) = 2.0*xl/((1.0-z*z)*pp*pp);
      w(n+1-i,1) = w(i);
    end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listparam(model)
    if numel(model.param)>0
        fprintf('_____________________________________________________________________________________\n');
        fprintf('%3s%20s%50s%12s\n','No','Name','Description','Value');
        fprintf('_____________________________________________________________________________________\n');
        k=1;
        for parm=model.param
            fprintf('%3d%20s%50s%12.4f\n',k,parm.ref(1:min(end,19)),parm.description(1:min(end,49)),parm.value);
            k=k+1;
        end
        fprintf('_____________________________________________________________________________________\n');
    else
        fprintf('No parameters defined in the model\n');
    end    
end

