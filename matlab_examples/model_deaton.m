%Simple consumption-savigns model (original Carroll's EGM method in Matlab)

function model_deaton() 

%parameters
Tbar	=	    25; %number of periods (fist period is t=1) 
ngridm	= 	   100; %number of grid points over assets
mmax	=		50; %maximum level of assets
expn	=		10; %number of quadrature points used in calculation of expectations
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
function w=income(it,shock) %income in period it with given normal shock
	%assume it=1 is age=20
	age=it+19;
	w=exp(1.5 + age*0.04 - 0.0004*age.*age + shock);
end %income

%main solver
function policy=egm_solver()
    policy=cell(Tbar,1); %output will be a cell array by time periods
    savingsgrid=linspace(0,mmax,ngridm); %grid over savings twice the size to accomodate bequest
    %main EGM loop
    for it=Tbar:-1:1
        if it==Tbar
            %terminal period
            policy{it}.w=[0 mmax];
            policy{it}.c=[0 mmax];
        else
            %not the terminal period
            %egm step
            policy{it}.c(1)=0;
            policy{it}.w(1)=0;
            %vectorize!!!
            %next period wealth
            wk1=ones(expn,1)*savingsgrid*(1+r)...
                +income(it+1,quadstnorm*sigma)*ones(1,ngridm);
            %next period consumption
            cons1=interp1(policy{it+1}.w,policy{it+1}.c,wk1,'linear','extrap');
            %RHS of Euler equation
            rhs=quadw'*mutil(cons1);
            %current period consumption
            cons0=imutil(df*(1+r)*rhs);
            %finilize EGM step
            policy{it}.c(2:ngridm+1)=cons0;
            policy{it}.w(2:ngridm+1)=savingsgrid+cons0;
        end % if(terminal period)
    end %it
end %egm_solver


%Run the solver
tic;
policy=egm_solver();
t=toc;
fprintf('Solved Deaton model with\n %d periods,\n %d endogenous grid points and\n %d quadrature points in \n%s\n',Tbar,ngridm,expn,ht(t));

%export the soluton to workspace
assignin('base','deaton_policy',policy);

%Make some graphs of the solution
%Income at all quadrature points
figure1=figure('Color',[1 1 1],'Name','deaton_income','NextPlot','new');
ax=axes('Parent',figure1);
title1=title(ax,'Incomes in deaton model');
hold(ax,'all');
for i=1:expn
	[x y]=fplot(@(it) income(it,quadstnorm(i)*sigma),[1 Tbar]);
	plot(ax,x,y,'DisplayName',sprintf('quadrature point %d normalized %1.3f',i,quadstnorm(i)));
	hold(ax,'all');
end
%Solution graph
figure1=figure('Color',[1 1 1],'Name','deaton_policy','NextPlot','new');
ax=axes('Parent',figure1);
title1=title(ax,'EGM solution to deaton model (consumption rule)');
hold(ax,'all');
for it=1:Tbar
	plot(ax,policy{it}.w,policy{it}.c,'DisplayName',sprintf('it=%d',it));
	hold(ax,'all');
end
set(ax,'XLim',[0 mmax]);

end %main function