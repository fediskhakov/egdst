%Simple consumption-savigns model + value function calculation

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
    savingsgrid=linspace(0,mmax,ngridm); %grid over savings
    %main EGM loop
    for it=Tbar:-1:1
        if it==Tbar
            %terminal period
            policy{it}.w=[0 mmax];
            policy{it}.c=[0 mmax];
            policy{it}.vf=[0 NaN]; %vf(1)=0 is important, otherwise vf is not computed in terminal period
        else
            %not the terminal period
            %egm step
            policy{it}.c(1)=0;
            policy{it}.w(1)=0;
            policy{it}.vf(1)=NaN;%filled out later as special value
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

            %now complute value function
            mask=(wk1<policy{it+1}.w(2)); %mask off credit constrained region
            if it==Tbar-1 
                mask=(ones(size(wk1))==1); %all points are credit constrained in last period
            end
            %compute expected value
            ev=NaN*ones(size(wk1));
            %for points wk1 falling into credit constrained region vf next period is computed directly
            %note vf(1) has a special value of expected future value at credit constraint (zero savings)
            ev(mask)=util(wk1(mask))+df*policy{it+1}.vf(1);
            %for points wk1 falling outside of credit constrained region vf is interpolated
            %note we drop the first point vf(1) with special value
            ev(~mask)=interp1(policy{it+1}.w(2:end),policy{it+1}.vf(2:end),wk1(~mask),'linear','extrap');
            %compute the expectation
            ev=quadw'*ev;
            %save expected value of zero savings in vf(1)
            if abs(savingsgrid(1))<eps
                policy{it}.vf(1)=ev(1);
            else
            	error 'The grid on end-of-period wealth does not start with 0 which is required!'
            end
            %finilize the computation of value function
            policy{it}.vf(2:ngridm+1)=util(cons0) + df*ev;
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
%Solution graphs
var={'c','vf'};
for i=1:numel(var)
	figure1=figure('Color',[1 1 1],'Name',sprintf('deaton_%s',var{i}),'NextPlot','new');
    ax=axes('Parent',figure1);
	title1=title(ax,sprintf('EGM solution to deaton model (%s)',var{i}));
	hold(ax,'all');
	for it=1:Tbar
		if (strcmp(var{i},'vf'))
			h=plot(ax,policy{it}.w(2:end),policy{it}.(var{i})(2:end),'DisplayName',sprintf('it=%d',it));
			hold(ax,'all');
			%analytic part of value function (credit constrained region)
			[x y]=fplot(@(x) util(x)+df*policy{it}.vf(1),[0.001 policy{it}.w(2)]);
			plot(ax,x,y,'DisplayName',sprintf('it=%d (analytic)',it),'Color',get(h,'Color'),'LineStyle','--');
			hold(ax,'all');
		else
			plot(ax,policy{it}.w,policy{it}.(var{i}),'DisplayName',sprintf('it=%d',it));
			hold(ax,'all');
		end
	end
	set(ax,'XLim',[0 mmax]);
end

end %main function