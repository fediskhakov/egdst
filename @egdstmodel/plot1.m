% Plots the solution of the model: optimal consumption, savings rules
% discrete choices or value function for all or selected states and time periods
% ---------------------------------------------------------------------------------- 
%   Input: none OR
%          '<title string>' OR    
%          'c|a|d|vf|v' OR
%          'it=NNN' OR  (it from [t0,T])
%          'ist=NNN' OR
%          'lines|ln|l' OR (draw vertical lines)
%          axes_handle OR          
%          any combination of the above
%  Output: structure containing:
%          1. vector of handles to created lineseries object (lines)
%          2. vector of handles to created lineseries object (thresholds)

function res=plot1(model,varargin)
hlines=[];
hthresholds=[];
if isempty(model.M)
    error 'No solution to draw! The model must first be solved.'
end

%default parameter values
wh=1;%what to plot: consumption
whlabel='C*';
miny=0;%start from 0, not a0
setit=model.t0:model.T;%all time period
setist=intersect([1 2],1:model.nst);%only first two state
titlestr=model.label;
drawlines=0;
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
    if ismember(varargin{i},{'C','A','VF','c','a','vf','Vf','vF','v','V','d','id','D','ID','Id','iD','lines','ln','l'})
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
        elseif ismember(varargin{i},{'lines','ln','l'})
            drawlines=1;
        end
    elseif numel(regexp(varargin{i},'^it='))>0
        setit=str2num(regexprep(varargin{i},'^it=',''));
        setit=intersect(setit,model.t0:model.T);
        if isempty(setit) 
            warning 'Nothing to do: set of time periods to plot is empty'
            return
        end
    elseif numel(regexp(varargin{i},'^ist='))>0
        setist=str2num(regexprep(varargin{i},'^ist=',''));
        setist=intersect(setist,1:model.nst);
        if isempty(setist) 
            warning 'Nothing to do: set of state indeces to plot is empty'
            return
        end
    else
        titlestr=varargin{i};
    end
end

%data
if ismember(wh,[1 2 3])
    data.M=model.M;
    data.D=model.D;
elseif wh==4
    %x position offset
    dx=1/(numel(setist)+1);
    %for errorbar plot: decisions in colums, NaN to skip
    data.x=NaN*ones(1000*model.nt,model.nd);
    data.y=NaN*ones(1000*model.nt,model.nd);
    data.up=NaN*ones(1000*model.nt,model.nd);
    data.lw=NaN*ones(1000*model.nt,model.nd);
    indx=ones(model.nd);
    %collect data
    for itc=reshape(setit,1,numel(setit))
        it=itc-model.t0+1;
        for isti=1:numel(setist)
            ist=setist(isti);
            if numel(model.D{ist,it})>0
                ths=[model.D{ist,it}; -1 model.mmax];
                for th=1:size(model.D{ist,it},1)
                    id=ths(th,1)+1;%decision
                    data.x(indx(id),id)=model.t0+(it-1)+(isti-1)*dx;
                    range=ths(th+1,2)-ths(th,2);
                    data.y(indx(id),id)=ths(th,2)+range/2;
                    data.up(indx(id),id)=range/2;
                    data.lw(indx(id),id)=range/2;
                    indx(id)=indx(id)+1;
                end
            end
        end
    end
    mask=(sum(isnan(data.x),2)==model.nd);
    data.x(mask,:)=[];
    data.y(mask,:)=[];
    data.up(mask,:)=[];
    data.lw(mask,:)=[];
end

%figure and axes
if isempty(ax)
    %create new figure
    ax=axes('Parent',figure);
    figure1=get(ax,'Parent');
else
    figure1=get(ax,'Parent');
end
%colormap in vectors
clrs=[colormap('Lines');colormap('Lines');colormap('Lines')];
set(figure1,'Name',[model.label ' ' whlabel]);
set(figure1,'Color',[1 1 1]);
title1=title(ax,titlestr);
hold(ax,'all');

if ismember(wh,[1 2 3])
    %plot
    set(ax,'xLim',[model.a0 model.mmax],'TickDir','out');
    if (wh~=3)
        set(ax,'yLim',[miny model.mmax-model.a0]);
    end
    hold(ax,'all');
    ylimmax=[];
    ylimmin=[];
    for itc=reshape(setit,1,numel(setit))
        it=itc-model.t0+1;
        for ist=reshape(setist,1,numel(setist))
            if numel(data.M{ist,it})>0
                %plot one line

                %amend vf data (and skip first line!)
                if wh==3
                    mgr=data.M{ist,it}(2:end,1);
                    fun=data.M{ist,it}(2:end,4);
                    evf=data.M{ist,it}(1,4);
                    if mgr(1)>model.a0+.01*model.mmax
                        for x=linspace(mgr(1),model.a0,min(model.ngridm,100))
                            if (x~=model.a0 && x~=mgr(1))
                                mgr=[x;mgr];
                                id=data.D{ist,it}(data.D{ist,it}(:,2)<=x,1);
                                id=id(end)+1;%optimal decision (base1)
                                consumption=interp1(data.M{ist,it}(:,1),data.M{ist,it}(:,2),x,'linear','extrap');%optimal consumption
                                fun=[model.call('utility',[itc ist id consumption])+model.call('discount',[itc ist])*evf;fun];
                            end
                        end
                    end
                else
                    mgr=data.M{ist,it}(:,1);
                    fun=data.M{ist,it}(:,wh+1);
                end
                dc=data.D{ist,it}(:,1);
                th=data.D{ist,it}(:,2);
                clr=clrs(mod(it,size(clrs,1)),1:3);
                hlines=[hlines ...
                    plot(ax,mgr,fun,'Color',clr,'DisplayName',['st=' sprintf('%0d',ist)]) ...
                        ];
                if itc==setit(1) || itc==setit(floor(numel(setit)/2)) || itc==setit(end)
                    mf=[mgr,fun];
                    mf=unique(mf,'rows');
                    text(max(mgr),interp1(mf(:,1),mf(:,1),max(mgr),'linear','extrap'),['t=' sprintf('%0d',itc)]);
                end
                thp=interp1q(mgr,fun,th);
                thp(isnan(thp))=0;
                for i=0:model.nd-1
                    %each decision with its own color and text label
                    hthresholds=[hthresholds ...
                        plot(ax,th(dc==i),thp(dc==i),'Marker','square','MarkerEdgeColor',clrs(i+1,1:3),'LineStyle','none','DisplayName','');
                                ];
                end
                if drawlines
                    %vertical lines through discontinuity points
                    dpd=str2num(model.cflags.DOUBLEPOINT_DELTA);
                    dscnt=[mgr [NaN;mgr(1:end-1)]];
                    dscnt=(abs(dscnt(:,1)-dscnt(:,2))<=1.5*dpd & abs(dscnt(:,1)-model.a0)>=1e-5);
                    dscnt=mgr(dscnt);%discontinous points
                    for i=1:numel(dscnt)
                        thdscnt=(abs(th-dscnt(i))<1e-6);
                        yl=get(ax,'YLim');
                        if sum(thdscnt)>0
                            %threshold
                            line([dscnt(i) dscnt(i)]',yl,'Parent',ax,'LineStyle','-','Color',[.5 .5 .5]);
                        else
                            line([dscnt(i) dscnt(i)]',yl,'Parent',ax,'LineStyle','--','Color',[.5 .5 .5]);
                        end
                    end


                end
                drawnow;

            else
                fprintf ('Solution for it=%d ist=%d does not exist (feasibility?), skipping..\n',it,ist);
            end
        end
    end

    %output as structure
    res=struct('lines',hlines,'dots',hthresholds);

elseif wh==4
    %do the plot
    h=errorbar(data.x,data.y,data.lw,data.up,...
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
        if sum(~isnan(data.x(:,id)))>0
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

    %Annotation for ist
    if numel(setist)>1
        istlabel='States shown: ';
    else
        istlabel='State shown: ';
    end
    for ist=reshape(setist,1,numel(setist))
        if numel(model.s)==1 && model.s(1).discrete
            istlabel=[istlabel model.s.values(ist).description '; '];
        else
            istlabel=[istlabel '[ ' sprintf('%1.4f ',model.states(ist,:)) '], '];
        end
    end
    istlabel(end-1:end)=[];
    text1=text(0,0,istlabel,'HorizontalAlignment','center','Units','normalized','VerticalAlignment','top');
    set(text1,'Position',[0.5 -.1])

    %output as structure
    res=struct('erbars',h);

end

end%of function
