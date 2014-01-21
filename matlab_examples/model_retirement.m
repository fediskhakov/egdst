%Simple retirement model and multiple Euler equation solutions prolems

function model_retirement() 

%parameters
% Tbar	=	    70; %number of periods (fist period is t=1) 
Tbar    =        6; %number of periods (fist period is t=1) 
ngridm	= 	   100; %number of grid points over assets
mmax	=		50; %maximum level of assets
expn	=		50; %number of quadrature points used in calculation of expectations
dorefinements =  1; %whether to do the refinement of discrete choice specific value functions
doplots =        1; %whether to build plots after the solver is done
%model parameters
r 		=	  0.05; %interest rate
df 		=	  0.95; %discount factor
sigma   =     0.25; %sigma parameter in income shocks
duw     =      2.3; %disutility of work
cfloor  =    0.001; %consumption floor (safety net in retirement)

%quadrature points
[quadp quadw]=quadpoints(expn,0,1);
quadstnorm=norminv(quadp,0,1);

%internal functions
function u=util(consumption,working) %utility
	u=log(consumption)-duw*(working==1);
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
function [policy dpolicy]=egm_solver()
    policy=cell(Tbar,1); %actual solution policy, including the optimal thresholds for retirement
    dpolicy=cell(Tbar,2); %cell array of conumption rules by time periods and discrete choices
    savingsgrid=linspace(0,mmax,ngridm); %grid over savings twice the size to accomodate bequest
    if(savingsgrid(1)~=0)
        error 'The grid on end-of-period wealth does not start with 0 which is required!'
    end
    %main EGM loop for each discrete choice
    for it=Tbar:-1:1
        for id=1:2
            working=id-1;%id=1 for retired, id=2 for working
            if it==Tbar
                %terminal period: consumption independent of work
                dpolicy{it,id}.w=[0 mmax];
                dpolicy{it,id}.c=[0 mmax];
                dpolicy{it,id}.vf=[0 NaN]; %vf(1)=0 is important, otherwise vf is not computed in terminal period
            else
                %not the terminal period
                %egm step
                dpolicy{it,id}.c(1)=0;
                dpolicy{it,id}.w(1)=0;
                dpolicy{it,id}.vf(1)=NaN;%filled out later as special value
                %vectorize!!!
                %next period wealth
                wk1=ones(expn,1)*savingsgrid*(1+r)...
                    +income(it+1,quadstnorm*sigma)*ones(1,ngridm)*(working==1); %income only if working
                %next period consumption - use actual policy, whether working or not tomorrow is automatically accounted for
                cons1=interp1(policy{it+1}.w,policy{it+1}.c,wk1,'linear','extrap');
                cons1=max(cons1,cfloor);%safety net consumption floor
                %RHS of Euler equation
                rhs=quadw'*mutil(cons1);
                %current period consumption
                cons0=imutil(df*(1+r)*rhs);
                %finilize EGM step
                dpolicy{it,id}.c(2:ngridm+1)=cons0;
                dpolicy{it,id}.w(2:ngridm+1)=savingsgrid+cons0;
                %now complute value function
                mask=(wk1<policy{it+1}.w(2)); %mask off credit constrained region next period
                if it==Tbar-1 
                    mask=(ones(size(wk1))==1); %all points are credit constrained in last period
                end
                %it's safe to assume here that optimal decision in credit constrained region 
                %is to work in every period except the terminal period
                %in a more general model a mask for the discrete decisions next period
                %would be needed
                %compute expected value
                ev=NaN*ones(size(wk1));
                %for points wk1 falling into credit constrained region vf next period is computed directly
                %note vf(1) has a special value of expected future value at credit constraint (zero savings)
                %we also distinguish whether working or not since it enters utility
                if it==Tbar-1
                    ev(mask)=util(max(cfloor,wk1(mask)),0)+df*policy{it+1}.vf(1); %retired in t+1
                else
                    ev(mask)=util(max(cfloor,wk1(mask)),1)+df*policy{it+1}.vf(1); %working in t+1
                end
                %for points wk1 falling outside of credit constrained region vf is interpolated
                %note we drop the first point vf(1) with special value
                ev(~mask)=interp1(policy{it+1}.w(2:end),policy{it+1}.vf(2:end),wk1(~mask),'linear','extrap');
                %compute the expectation
                ev=quadw'*ev;
                %save expected value of zero savings in vf(1)
                %we know at credit constraint working is optimal, but save both for graphs
                if abs(savingsgrid(1))<eps
                    dpolicy{it,id}.vf(1)=ev(1);
                end
                %finilize the computation of value function
                dpolicy{it,id}.vf(2:ngridm+1)=util(cons0,working) + df*ev;

                %d-specific value function is computed
                %but because the value function in t+1 may be non-concave,
                %there can be multiple solutions to Euler equation and wiggling
                %in computed d-specific value function
                %We need to filter out these "wiggles" where grids folds back on itself
                %Justification for the search criteria: monotonicity of savings function (see paper)
                if dorefinements
                    %this part is very slow..
                    j=find(dpolicy{it,id}.w(2:end)<dpolicy{it,id}.w(1:end-1),1,'first'); %index of last ok point before fold over
                    while ~isempty(j)
                        %total number of points in fold over
                        j1=find(dpolicy{it,id}.w<dpolicy{it,id}.w(j));
                        j1=sum(j1>j);%the number of points on the fold over interval
                        %number of point to skip after the start of fold over
                        %interpolate all points in the fold over using part of grid we know is ok
                        j2=dpolicy{it,id}.vf(j+1:j+j1) < interp1(dpolicy{it,id}.w(1:j),dpolicy{it,id}.vf(1:j),dpolicy{it,id}.w(j+1:j+j1),'linear');
                        j2=sum(j2);%the number of points of the fold over lieing below the ok grid
                        %number of points to skip before the start of fold over
                        j3=find(dpolicy{it,id}.w>dpolicy{it,id}.w(j+j2+1));
                        j3=sum(j3<=j);%the number of points on the fold over interval
                        %verbosity
                        % fprintf('env2 (it=%d id=%d) point %d: w=%1.4f vf=%1.4f : %d in fold over\n%d of points to be skipped after the start of fold over\n%d of points to be skipped before the start of fold over\n',it,id,j,dpolicy{it,id}.w(j),dpolicy{it,id}.vf(j),j1,j2,j3);
                        %skip inferior points
                        dpolicy{it,id}.w(j-j3+1:j+j2)=[];
                        dpolicy{it,id}.c(j-j3+1:j+j2)=[];
                        dpolicy{it,id}.vf(j-j3+1:j+j2)=[];
                        %search for next fold over region
                        j=find(dpolicy{it,id}.w(2:end)<dpolicy{it,id}.w(1:end-1),1,'first'); %index of last ok point before fold over
                    end %while problematic fold over regions exist
                end %doing refinement of the d-specific value function
            end % if(terminal period)
        end %working
        %compare the discrete choices using d-specific value functions
        if it==Tbar
            % in the last period it's always optimal to retire (and consumpe everything)
            policy{it}.w=[0 mmax];
            policy{it}.c=[0 mmax];
            policy{it}.vf=[0 NaN]; %vf is not computed in terminal period, special value in vf(1)
            policy{it}.th=0; %retire for any level of assets --> threshold=0
        else
            %we have two different EGM grids in dpolicy{it,id} id=1,2
            %and have to compute an upper envelope of the two value functions
            %and to recover the intersection point which is a threshold
            %important: we know that there is one threshold, and vfs can not
            %intersect in the credit constraned region
            %1 points of the RETIRED grid to be on the upper envelope
            %mask off points of the RETIRED grid that fall into credit constraint of WORKING
            mask=(dpolicy{it,1}.w<dpolicy{it,2}.w(2));
            mask(1)=false; %skip the first point (with special value vf(1)
            mask1=false; %first point in mask indicating points on retirement grid where retirement is optimal
            %add indicators in credit constrained region
            mask1=[mask1 (dpolicy{it,1}.vf(mask)>util(max(cfloor,dpolicy{it,1}.w(mask)),1)+df*dpolicy{it,2}.vf(1))];
            %add indicators outside of credit constrained region
            mask=~mask;%invert the mask
            mask(1)=false;%skip the first point too (i.e. special value vf(1) 
            mask1=[mask1 (dpolicy{it,1}.vf(mask)>interp1(dpolicy{it,2}.w(2:end),dpolicy{it,2}.vf(2:end),dpolicy{it,1}.w(mask),'linear','extrap'))];
            %2 points of the WORKING grid to be on the upper envelope
            %we should do similar thing as before, but we know that credit constrained region for retirement is tiny,
            %and that retirement with little money is not optimal before the terminal period
            mask2=false; %first point in mask indicating points on working grid where working is optimal
            mask2=[mask2 (dpolicy{it,2}.vf(2:end)>interp1(dpolicy{it,1}.w(2:end),dpolicy{it,1}.vf(2:end),dpolicy{it,2}.w(2:end),'linear','extrap'))];
            %find the intersection point (threshold)
            i1=find(mask1,1,'first')-1; %left point of the retirement grid before intersection
            i2=find(mask2,1,'last');%left point of the working grid before intersection
            if i2==numel(mask2)
                i2=i2-1; %use next to last point for extrapolation
            end
            %distinguish the case when intersection point lies in the credit constraint of working
            if isempty(i1)
                %retirement is never optimal
                isect.w=[];
                isect.c=[];
                isect.vf=[];
            else
                if isempty(i2)
                    %intersection on in the cc region of working grid
                    %equation to solve
                    criterion=@(w) util(w,1)+df*dpolicy{it,2}.vf(1)-interp1(dpolicy{it,1}.w(2:end),dpolicy{it,1}.vf(2:end),w,'linear','extrap');
                    isect.w=fzero(criterion,dpolicy{it,2}.w(2));
                    isect.w=[isect.w isect.w+100*eps];
                    %value of the intersection point
                    isect.vf=[util(max(cfloor,isect.w(1)),1)+df*dpolicy{it,2}.vf(1) ... 
                              util(max(cfloor,isect.w(2)),1)+df*dpolicy{it,2}.vf(1)];
                else
                    a1=(dpolicy{it,1}.vf(i1+1)-dpolicy{it,1}.vf(i1))/(dpolicy{it,1}.w(i1+1)-dpolicy{it,1}.w(i1));
                    a2=(dpolicy{it,2}.vf(i2+1)-dpolicy{it,2}.vf(i2))/(dpolicy{it,2}.w(i2+1)-dpolicy{it,2}.w(i2));
                    b1=dpolicy{it,1}.vf(i1)-a1*dpolicy{it,1}.w(i1);
                    b2=dpolicy{it,2}.vf(i2)-a2*dpolicy{it,2}.w(i2);
                    isect.w=[(b2-b1)/(a1-a2) ((b2-b1)/(a1-a2))+100*eps];%twice the intersection point
                    %value of the intersection point
                    isect.vf=[a2*isect.w(1)+b2 a1*isect.w(2)+b1];
                end
                %consumption from two rules
                isect.c=[interp1(dpolicy{it,2}.w,dpolicy{it,2}.c,isect.w(1),'linear','extrap') ...
                         interp1(dpolicy{it,1}.w,dpolicy{it,1}.c,isect.w(2),'linear','extrap')];
            end
            %combine the points from two grids
            policy{it}.w=[0 dpolicy{it,2}.w(mask2) isect.w dpolicy{it,1}.w(mask1)];
            policy{it}.c=[0 dpolicy{it,2}.c(mask2) isect.c dpolicy{it,1}.c(mask1)];
            policy{it}.vf=[dpolicy{it,2}.vf(1) dpolicy{it,2}.vf(mask2) isect.vf dpolicy{it,1}.vf(mask1)];
        end
    end %it
end %egm_solver


%Run the solver
tic;
[policy dpolicy]=egm_solver();
t=toc;
fprintf('Solved retirement model with\n %d periods,\n %d endogenous grid points and\n %d quadrature points in \n%s\n',Tbar,ngridm,expn,ht(t));

%export the soluton to workspace
assignin('base','retirement_policy',policy);
assignin('base','retirement_dpolicy',dpolicy);

%return here if plots are not needed
if ~doplots
    return
end

%Make some graphs of the solution
%Income at all quadrature points
figure1=figure('Color',[1 1 1],'Name','retirement_income','NextPlot','new');
ax=axes('Parent',figure1);
title1=title(ax,'Incomes in deaton model');
hold(ax,'all');
for i=1:expn
	[x y]=fplot(@(it) income(it,quadstnorm(i)*sigma),[1 Tbar]);
	plot(ax,x,y,'DisplayName',sprintf('quadrature point %d normalized %1.3f',i,quadstnorm(i)));
	hold(ax,'all');
end
%Discrete chocie specific graphs
figure1=figure('Color',[1 1 1],'Name','retirement_dvf','NextPlot','new');
for j=1:6
    it=Tbar+1-j;
    if it<1
        break
    end
    ax=subplot(2,3,j,'Parent',figure1);
    title1=title(ax,sprintf('Value functions in period it=%d',it));
    hold(ax,'all');
    %working
    h=plot(ax,dpolicy{it,2}.w(2:end),dpolicy{it,2}.vf(2:end),'DisplayName',sprintf('it=%d working',it));
    hold(ax,'all');
    %analytic part of value function (credit constrained region)
    [x y]=fplot(@(x) util(x,1)+df*dpolicy{it,2}.vf(1),[0.001 dpolicy{it,2}.w(2)]);
    plot(ax,x,y,'DisplayName',sprintf('it=%d (analytic)',it),'Color',get(h,'Color'),'LineStyle','--');
    hold(ax,'all');
    %retired
    h=plot(ax,dpolicy{it,1}.w(2:end),dpolicy{it,1}.vf(2:end),'DisplayName',sprintf('it=%d retired',it));
    hold(ax,'all');
    %analytic part of value function (credit constrained region)
    [x y]=fplot(@(x) util(x,0)+df*dpolicy{it,1}.vf(1),[0.001 dpolicy{it,1}.w(2)]);
    plot(ax,x,y,'DisplayName',sprintf('it=%d (analytic)',it),'Color',get(h,'Color'),'LineStyle','--');
    hold(ax,'all');
    set(ax,'XLim',[0 mmax]);
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
			[x y]=fplot(@(x) util(x,1)+df*policy{it}.vf(1),[0.001 policy{it}.w(2)]);
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