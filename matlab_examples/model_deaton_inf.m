%Simple consumption-savigns model (original EGM) - infinite horizon

function model_deaton_inf() 

%parameters
ngridm	= 	   100; %number of grid points over assets
mmax	=		50; %maximum level of assets
expn	=		10; %number of quadrature points used in calculation of expectations
tol     =     1e-8; %tolerance of conversion of optimal consumption rule
maxiter =      1e6; %maximum number of iterations
r 		=	  0.05; %interest rate
df 		=	  0.95; %discount factor
sigma   =     0.55; %sigma parameter in income shocks

%quadrature points
[quadp quadw]=quadpoints(expn,0,1);
quadstnorm=norminv(quadp,0,1);

%internal functions
function u=util(consumption) %utility
	u=log(consumption);
end %util
function mu=mutil(consumption) %marginal utility
	mu=1./consumption;
end %mutil
function cons=imutil(mutil) %inverse marginal utility
	cons=1./mutil;
end %imutil
function w=income(shock) %income in period it with given normal shock
	w=exp(1.5 + shock);
end %income

%main solver
function policy=egm_solver()
    savingsgrid=linspace(0,mmax,ngridm); %grid over savings
    %main EGM loop
    for it=1:maxiter
        if it==1
            %first iteration, assume consumption of everything is optimal rule
            policy{1}.w=[0 mmax];
            policy{1}.c=[0 mmax];
        else
            %not the terminal period
            %egm step
            policy{1}.c(1)=0;
            policy{1}.w(1)=0;
            %vectorize!!!
            %next period wealth
            wk1=ones(expn,1)*savingsgrid*(1+r)...
                +income(quadstnorm*sigma)*ones(1,ngridm);
            %next period consumption
            cons1=interp1(policy{2}.w,policy{2}.c,wk1,'linear','extrap');
            %RHS of Euler equation
            rhs=quadw'*mutil(cons1);
            %current period consumption
            cons0=imutil(df*(1+r)*rhs);
            %finilize EGM step
            policy{1}.c(2:ngridm+1)=cons0;
            policy{1}.w(2:ngridm+1)=savingsgrid+cons0;
        end % if(terminal period)
        %check for convergience
        if it>10
            measure=max(abs(policy{1}.c - interp1(policy{2}.w,policy{2}.c,policy{1}.w,'linear','extrap')));
            if measure<tol
                break;
            end
        end
        %replace next period with current
        policy{2}=policy{1};
    end %it
end %egm_solver

%Run the solver
tic;
policy=egm_solver();
t=toc;
fprintf('Solved Deaton model with\n %1.0e tolerance on consumption rule convergence,\n %d endogenous grid points and\n %d quadrature points in\n %d interations and in \n%s\n',tol,ngridm,expn,ht(t));

%export the soluton to workspace
assignin('base','deaton_inf_policy',policy{1});

%Make some graphs of the solution
figure1=figure('Color',[1 1 1],'Name','deaton_policy','NextPlot','new');
ax=axes('Parent',figure1);
title1=title(ax,'EGM solution to deaton model with infinite horizon (consumption rule)');
hold(ax,'all');
plot(ax,policy{1}.w,policy{1}.c,'DisplayName','Optimal consumption rule');
ylim=get(ax,'Ylim');
hold(ax,'all');
plot(ax,[0 mmax]',[0 mmax]','DisplayName','45''','LineStyle','--');
set(ax,'XLim',[0 mmax],'Ylim',1.2*ylim);

end %main function