%Occupational choice with 3 options

occ3=egdstmodel('Occupational choice model ','tmp_occ3');
%time and space
occ3.t0=0;
occ3.T=40;
occ3.s={'Dummy state',{0,'dummy'}};
occ3.trpr={'true',[1]};
occ3.feasible={'defaultfeasible',true};
occ3.d={'Occupational choice',{0,'public sector (lower pay, secure)',...
                               1,'private sector (hight pay, less secure)',...
                               2,'entrepreneurship'}};
occ3.choiceset={'defaultallow',true};
%utility
occ3.u={'utility','(pow(consumption,1-crra)-1)/(1-crra) - coefleisure*disutility[1][(int)dc1+1]'};
occ3.coef={'disutility','Disutility of work',[0.0 1.0 0.75]};
occ3.param={'crra','CRRA coefficient in utility',1.2};
occ3.param={'coefleisure','Weight with leisure in utility',0.2};
occ3.u={'marginal','pow(consumption,-crra)'};
occ3.u={'marginalinverse','pow(mutility,-1/crra)'};
occ3.u={'extrap','pow(x,1-crra)'};
occ3.discount='0.93';
%shocks
occ3.shock='lognormal';
occ3.shock={'sigma','sigs[1][(int)dc1+1]'};
occ3.coef={'sigs','Sigmas for different occupations',[0.15 0.35 0.75]};
occ3.shock={'mu','-0.5*sigs[1][(int)dc1+1]*sigs[1][(int)dc1+1]'};
%income equations
occ3.eq={'wage1','Realized wage in the public sector','max(ssinc,shock*0.5)','next'};
occ3.eq={'wage2','Realized wage in the private sector','max(ssinc,shock*0.5*wagegap)','next'};
occ3.param={'ssinc','Guaranteed social security income',0.01};
occ3.param={'wagegap','Wage gap between public and private sector',1.35};
occ3.eq={'entrep','Entrepreneurial income (realized)','max(ssinc,log(savings+1)*entrkap*shock)','next'};
occ3.param={'entrkap','Return on capital',0.56};
%budget
occ3.budget={'cashinhand','savings*(1+interest)+(dc1==0)*wage1+(dc1==1)*wage2+(dc1==2)*entrep'};
occ3.budget={'marginal','1+interest+(dc1==2)*max(0,entrkap*shock/(savings+1))'};
occ3.param={'interest','return on savings',0.05};
occ3.a0=0;
%grids
occ3.mmax=5;           %max cash-at-hand                    
occ3.ngridm=50;        %standard number of grid points in M
occ3.ngridmax=100;     %max number of grid point in M
occ3.nthrhmax=100;     %max number of treshold points in TH
occ3.ny=10;            %number of points in discrete representation of income shocks
%solve, plot and sim
occ3.compile

occ3.solve
occ3.plot1('c')
occ3.plot1('d')

occ3.sim([1 1.3725])
occ3.plot2('mack-d-kq');






