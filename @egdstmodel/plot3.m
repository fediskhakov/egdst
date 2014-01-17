% Plots the solution of the model along the simulated path:
% optimal consumption, savings rules, value function or discrete choices
% ---------------------------------------------------------------------------------- 
%   Input: none OR
%          '<title string>' OR    
%          'c|a|d|vf|v' OR
%          'it=NNN' OR  (it from [t0,T])
%          'sims=NNN' OR (from 1 to number of sims)
%          axes_handle OR          
%          any combination of the above
%  Output: structure containing:
%          1. vector of handles to created lineseries object (lines)
%          2. vector of handles to created lineseries object (thresholds)
%          3. vector of handles to created path line segments (path)
%          4. handle to created stairs plot (stairs)

function res=plot3(model,varargin)
hlines=[];
hthresholds=[];
hpaths=[];
if isempty(model.M)
    error 'No solution to draw! The model must first be solved.'
end
if isempty(model.sims)
    error 'No simulation to draw! The model must first be solved and simulated.'
end

%default parameter values
wh=1;%what to plot: consumption
whlabel='C*';
miny=0;%start from 0, not a0
setit=model.t0:model.T;%all time period
simiset=1:min(3,size(model.sims,1)); %max nr of figures if not asked specifically
titlestr=model.label;
ax=[];

%parse input
for i=1:nargin-1
    if numel(varargin{i})==1 && ishandle(varargin{i}) && strcmp(get(varargin{i},'Type'),'axes')
        ax=varargin{i};
        continue
    end
    if ~ischar(varargin{i})
        error 'Wrong parameter! Only strings are accepted.\nRun help egdstmodel.plot1 for help'
    end 
    if ismember(varargin{i},{'C','A','VF','c','a','vf','Vf','vF','v','V','d','id','D','ID','Id','iD'})
        %what to do
        if ismember(varargin{i},{'C','c'})
            wh=1;
            whlabel='C*';
            miny=0;
        elseif ismember(varargin{i},{'A','a'})
            wh=2;
            whlabel='A*';
            miny=model.a0;
        elseif ismember(varargin{i},{'VF','vf','Vf','vF','v','V'})
            wh=3;
            whlabel='VF';
            miny=model.a0;
        elseif ismember(varargin{i},{'d','id','D','ID','Id','iD'})
            wh=4;
            whlabel='D*';
            miny=model.a0;
        end
    elseif numel(regexp(varargin{i},'^it='))>0
        setit=str2num(regexprep(varargin{i},'^it=',''));
        setit=intersect(setit,model.t0:model.T);
        if isempty(setit) 
            warning 'Nothing to do: set of time periods to plot is empty'
            return
        end
    elseif numel(regexp(varargin{i},'^sims='))>0
        simiset=str2num(regexprep(varargin{i},'^sims=',''));
        simiset=intersect(simiset,1:size(model.sims,1));
        if isempty(simiset) 
            warning 'Nothing to do: set of simulations to plot is empty'
            return
        end
    else
        titlestr=varargin{i};
    end
end

%data
data.M=model.M;
data.D=model.D;

%cycle over simulations
for simi=simiset

    %figure and axes if no axes given, for each simi
    if isempty(ax) || simi>simiset(1)
        figure1=figure;
        set(figure1,'Name',[model.label ' sims ' sprintf('%g',simi)]);
        set(figure1,'Color',[1 1 1]);
        ax=axes('Parent',figure1);
    else
        figure1=get(ax,'Parent');
    end
    title1=title(ax,titlestr);
    hold(ax,'all');
    %colormap in vectors
    clrs=[colormap('Lines');colormap('Lines');colormap('Lines')];

    %plot
    if ismember(wh,[1 2 3]) % c a vf case

        set(ax,'xLim',[model.a0 model.mmax],'TickDir','out');
        if (wh~=3)
            set(ax,'yLim',[miny model.mmax-model.a0]);
        end
        hold(ax,'all');
        ylimmax=[];
        ylimmin=[];
        setit=sort(reshape(setit,1,numel(setit)));
        for itc=setit
            it=itc-model.t0+1;
            %find correct ist from sims
            ist=model.sims(simi,it,6)+1; %6 has ist SIMS OUTPUT SPECIFIC

            if numel(data.M{ist,it})>0
                %plot one line

                %amend vf data (and skip first line!)
                if wh==3
                    mgr=data.M{ist,it}(2:end,1);
                    fun=data.M{ist,it}(2:end,4);
                    evf=data.M{ist,it}(1,4);
                    if mgr(1)>model.a0+.01*model.mmax
                        for x=linspace(mgr(1)-.01*model.mmax,model.a0+.01*model.mmax,model.ngridm)
                            mgr=[x;mgr];
                            id=data.D{ist,it}(data.D{ist,it}(:,2)<=x,1);
                            id=id(end)+1;%optimal decision (base1)
                            consumption=interp1(data.M{ist,it}(:,1),data.M{ist,it}(:,2),x,'linear','extrap');%optimal consumption
                            fun=[model.call('utility',[itc ist id consumption])+model.call('discount',[itc ist])*evf;fun];
                        end
                    end
                else
                    mgr=data.M{ist,it}(:,1);
                    fun=data.M{ist,it}(:,wh+1);
                end
                dc=data.D{ist,it}(:,1);
                th=data.D{ist,it}(:,2);
                clr=clrs(mod(it,size(clrs,1)),1:3);
                %line
                hlines=[hlines ...
                    plot(ax,mgr,fun,'Color',clr,'DisplayName',['t=' sprintf('%0d',itc) ' st=' sprintf('%0d',ist)]) ...
                        ];
                if itc==setit(1) || itc==setit(floor(numel(setit)/2)) || itc==setit(end)
                    text(max(mgr),interp1(mgr,fun,max(mgr),'linear','extrap'),['t=' sprintf('%0d',itc)]);
                end
                %thresholds
                thp=interp1q(mgr,fun,th);
                thp(isnan(thp))=0;
                for i=0:model.nd-1
                    %each decision with its own color and text label
                    hthresholds=[hthresholds ...
                        plot(ax,th(dc==i),thp(dc==i),'Marker','square','MarkerEdgeColor',clrs(i+1,1:3),'LineStyle','none','DisplayName','');
                                ];
                end
                %path
                msim=model.sims(simi,it,1); %1 has cash-in-hand SIMS OUTPUT SPECIFIC
                fsim=interp1q(mgr,fun,msim);
                hpaths=[hpaths ...
                    line([msim msim],[0 fsim],'Color',[.2 .2 .2],'Parent',ax,'Marker','.','DisplayName',['t=' sprintf('%0d',itc) ' st=' sprintf('%0d',ist)]);
                        ];
                if itc>setit(1)
                    hpaths=[hpaths ...
                        line([lastpoint(1) msim],[lastpoint(2) 0],'Color',[.2 .2 .2],'Parent',ax,'DisplayName',['t=' sprintf('%0d',lastpoint(3)) '-->'  sprintf('%0d',itc)]);
                            ];
                end
                lastpoint=[msim fsim itc];%save for next iteration
                text(msim,fsim,sprintf('%0d',itc),'HorizontalAlignment','center','VerticalAlignment','bottom');
                text(msim,0,[sprintf('%0d',itc) ':' sprintf('%0d',ist)],'HorizontalAlignment','center','VerticalAlignment','top');

                drawnow;

            else
                fprintf ('Solution for it=%d ist=%d does not exist (feasibility?), skipping..\n',it,ist);
            end
        end
        
    elseif wh==4
        %data1 on discrete policy
        %for errorbar plot: decisions in colums, NaN to skip
        data1.x=NaN*ones(1000*model.nt,model.nd);
        data1.y=NaN*ones(1000*model.nt,model.nd);
        data1.up=NaN*ones(1000*model.nt,model.nd);
        data1.lw=NaN*ones(1000*model.nt,model.nd);
        indx=ones(model.nd);
        %collect data
        setit=sort(reshape(setit,1,numel(setit)));
        for itc=setit
            it=itc-model.t0+1;
            %find correct ist from sims
            ist=model.sims(simi,it,6)+1; %6 has ist SIMS OUTPUT SPECIFIC
            %solution for simulated ist
            if numel(model.D{ist,it})>0
                ths=[model.D{ist,it}; -1 model.mmax];
                for th=1:size(model.D{ist,it},1)
                    id=ths(th,1)+1;%decision
                    data1.x(indx(id),id)=model.t0+(it-1);
                    range=ths(th+1,2)-ths(th,2);
                    data1.y(indx(id),id)=ths(th,2)+range/2;
                    data1.up(indx(id),id)=range/2;
                    data1.lw(indx(id),id)=range/2;
                    indx(id)=indx(id)+1;
                end
            end
        end
        mask=(sum(isnan(data1.x),2)==model.nd);
        data1.x(mask,:)=[];
        data1.y(mask,:)=[];
        data1.up(mask,:)=[];
        data1.lw(mask,:)=[];
        
        %do the plot
        h=errorbar(ax,data1.x,data1.y,data1.lw,data1.up,...
            'MarkerSize',4,...
            'MarkerFaceColor',[1 1 1],...
            'LineStyle','none',...
            'LineWidth',1   ); 
        %Axes
        set(ax,'XLim',[model.t0-1 model.T+1]);
        set(ax,'YLim',[model.a0 model.mmax]);

        %line lables
        markers={'s','o','d','^','v','>','<'};
        for id=1:model.nd
            if sum(~isnan(data1.x(:,id)))>0
                nm=sprintf('%1.2f, ',model.decisions(id,:));
                nm(end-2:end)=[];
                set(h(id),'Marker',markers{mod(id,numel(markers))+1},'DisplayName',['dc=[' nm ']']);
            end
        end

        % Create legend
        legend1 = legend(ax,'show');
        set(legend1,...
            'EdgeColor',[1 1 1],...
            'Orientation','horizontal',...
            'Location','SouthOutside',...
            'YColor',[1 1 1],...
            'XColor',[1 1 1],...
            'Location','BestOutside');

        %stairs
        strs=stairs(ax,[model.t0:model.T+1]-.5,[model.sims(simi,:,1) 0]', ...
                    'DisplayName','Simulated cash-in-hand', ...
                    'Color',[.2 .2 .2]);
    end

end

%output as structure
if ismember(wh,[1 2 3])
    res=struct('lines',hlines,'dots',hthresholds,'path',hpaths);
elseif wh==4
    res=struct('erbars',h,'stairs',strs);
end
end%of function
