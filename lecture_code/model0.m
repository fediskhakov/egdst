%Simple consumption-savigns model used in the lecture presentation: Matlab implementation

function policy=model0()
%parameters
TBAR	=	    25; %number of periods (fist period is t=1) 
NM	= 	    10; %number of grid points over assets
MMAX	=	   100; %maximum level of assets
EXPN	=		10; %number of quadrature points used in calculation of expectations
R 		=	  0.05; %interest rate
DF 		=	  0.95; %discount factor
SIGMA   =     0.05; %sigma parameter in returns shocks

%quadrature points
[quadp quadw]=quadpoints(EXPN,0,1);
quadstnorm=norminv(quadp,0,1);

%Run the solver
tic;%start EGM loop
savingsgrid=linspace(0,MMAX,NM);%grid over savings
policy{TBAR}.w=[0 MMAX];%terminal period
policy{TBAR}.c=[0 MMAX];%terminal period
for it=TBAR-1:-1:1
  w1=exp(quadstnorm*SIGMA)*savingsgrid*(1+R);%next period wealth
  c1=interp1(policy{it+1}.w,policy{it+1}.c,w1,'linear','extrap');%next period consumption
  rhs=quadw'*(1./c1);%RHS of Euler equation (log utility)
  policy{it}.c=1./(DF*(1+R)*rhs);%current period consumption (log utility)
  policy{it}.w=savingsgrid+policy{it}.c;
end
t=toc;%stop the timer
fprintf('Solved Deaton model with\n %d periods,\n %d endogenous grid points and\n %d quadrature points in \n%s\n',TBAR,NM,EXPN,ht(t));

%Plot solution
figure1=figure('Color',[1 1 1],'Name','deaton_policy','NextPlot','new');
ax=axes('Parent',figure1);
for it=1:TBAR
	plot(ax,policy{it}.w,policy{it}.c,'DisplayName',sprintf('it=%d',it));
	hold(ax,'all');
end
set(ax,'XLim',[0 MMAX]);
title1=title(ax,'EGM solution to simple consumption model (consumption rule)');
end %main function