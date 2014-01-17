% Plots the simulations of the model: any simulated time serieses
% ---------------------------------------------------------------------------------- 
%   Input: none OR
%          '<title string>' OR
%          'm|c|a|v|i|j|p|g|k|u|b|d|s|q|-' OR
%          'sims=NNN' OR
%          'combine' OR
%          'noscale' OR
%          axes_handle OR
%          any combination of the above
%  Output: structure containing:
%          1. vector of handles to created lineseries object (lines)
%
%   Codes for the lines to plot:
%               m = cash-in-hand
%               c = consumption
%               a = savings
%               v = value function
%               i = current discrete decision index
%               j = current state index
%               p = mu location parameter of shock
%               g = sigma parameter of shock
%               k = income shock in current period
%               u = utility
%               b = discount factor
%               d OR dN = all decisions OR decision N
%               s OR sN = all states OR state N
%               q OR qN = all equations OR equation N
%               - = separator to create sub-plots
%       mack-ds-q = default combination
%
%   'combine' plots all simulations on the same figure/axes
%   'noscale' plots data without scaling (which is automatically done for model equations)
%   axes_handle allows to embed plots into existing axes - ONLY if plot3 is called 
%               for a single simulation sims=N (or model.sims contains only one simulation)
%               AND no groups are defined (symbol '-' is not given)
%

function res=plot2(model,varargin)
hlines=[];
if isempty(model.sims)
    error 'No simulations to draw! The model must be first solved and simulated.'
end

%default parameter values
wh='mack-ds-q'; %what to plot: main 4 quantities
titlestr=model.label;
simiset=1:min(3,size(model.sims,1)); %max nr of figures if not asked specifically
combine=0;
ax=[];
scale={0};

%parse input
for i=1:nargin-1
    if numel(varargin{i})==1 && ishandle(varargin{i}) && strcmp(get(varargin{i},'Type'),'axes')
        ax=varargin{i};
        continue
    end
    if ~ischar(varargin{i})
        error 'Wrong parameter! Only strings are accepted.\nRun help egdstmodel.plot3 for help'
    end
    %match the control symbols
    if ~isempty(regexp(varargin{i},'^[mcavijpgkubdsq-\d]+$'))
        wh=varargin{i};
    elseif numel(regexp(varargin{i},'^sims='))>0
        simiset=str2num(regexprep(varargin{i},'^sims=',''));
        simiset=intersect(simiset,1:size(model.sims,1));
        if isempty(simiset) 
            warning 'Nothing to do: set of simulations to plot is empty'
            return
        end
    elseif strcmp(varargin{i},'combine')
        combine=1;
    elseif strcmp(varargin{i},'noscale')
        scale={};
    else
        if ~strcmp(titlestr,model.label)
            warning 'More than one argument is interpreted as title string.'
        end
        titlestr=varargin{i};
    end
end

%translate letter notation for sim values into columns of sims 
% !!!! THIS IS SPECIFIC TO THE ORDER OF OUTPUT IN SIMS
grpi=1;
lns{1}=[];
lnprop={};
lnprop_comb={};
for k=1:numel(wh)
    switch wh(k)
    case 'm'
        lns{grpi}=[lns{grpi} 1];
        lnprop{1}={'MarkerSize',8,'Marker','o','LineWidth',1.5,'Color',[0 0 .75]};
    case 'c'
        lns{grpi}=[lns{grpi} 2];
        lnprop{2}={'MarkerSize',8,'Marker','square','LineWidth',1.5,'Color',[0 .5 0]};
    case 'a'
        lns{grpi}=[lns{grpi} 3];
        lnprop{3}={'Marker','none','LineStyle','--','MarkerFaceColor','none','LineWidth',1.5,'Color',[.5 0 0]};
    case 'v'
        lns{grpi}=[lns{grpi} 4];
        lnprop{4}={};
    case 'i'
        lns{grpi}=[lns{grpi} 5];
        lnprop{5}={'MarkerSize',17,'Marker','o','LineStyle','--'};
    case 'j'
        lns{grpi}=[lns{grpi} 6];
        lnprop{6}={'MarkerSize',12,'Marker','square'};
    case 'p'
        lns{grpi}=[lns{grpi} 7];
        lnprop{7}={};
    case 'g'
        lns{grpi}=[lns{grpi} 8];
        lnprop{8}={};
    case 'k'
        lns{grpi}=[lns{grpi} 9];
        lnprop{9}={'Marker','*','LineStyle','-','LineWidth',1,'MarkerFaceColor','none','Color','red'};
    case 'u'
        lns{grpi}=[lns{grpi} 10];
        lnprop{10}={};
    case 'b'
        lns{grpi}=[lns{grpi} 11];
        lnprop{11}={};
    case 's'
        nm=regexp(wh(k+1:end),'^\d+','match','once'); %check for numericals
        if isempty(nm)
            lns{grpi}=[lns{grpi} 12:11+numel(model.s)];
            for i=12:11+numel(model.s)
                lnprop{i}={'MarkerSize',10,'Marker','square'};
            end
        else
            i=i+numel(nm);%skip numericals in the next loop
            nm=str2num(nm);
            if nm<=numel(model.s)
                lns{grpi}=[lns{grpi} 11+nm];
                lnprop{11+nm}={'MarkerSize',10,'Marker','square'};
            else
                warning 'Given index for state is too large'
            end
        end
        if numel(scale)>0
            scale={scale{:} grpi};
        end
    case 'd'
        nm=regexp(wh(k+1:end),'^\d+','match','once'); %check for numericals
        if isempty(nm)
            lns{grpi}=[lns{grpi} 12+numel(model.s):11+numel(model.s)+numel(model.d)];
            for i=12+numel(model.s):11+numel(model.s)+numel(model.d)
                lnprop{i}={'MarkerSize',15,'Marker','o','LineStyle','--'};
            end
        else
            i=i+numel(nm);%skip numericals in the next loop
            nm=str2num(nm);
            if nm<=numel(model.d)
                lns{grpi}=[lns{grpi} 11+numel(model.s)+nm];
                lnprop{11+numel(model.s)+nm}={'MarkerSize',10,'Marker','square'};
            else
                warning 'Given index for decision is too large'
            end
        end
    case 'q'
        nm=regexp(wh(k+1:end),'^\d+','match','once'); %check for numericals
        if isempty(nm)
            lns{grpi}=[lns{grpi} 12+numel(model.s)+numel(model.d):11+numel(model.s)+numel(model.d)+numel(model.eq)];
            for i=12+numel(model.s)+numel(model.d):11+numel(model.s)+numel(model.d)+numel(model.eq)
                lnprop{i}={'Marker','none','MarkerSize',6};
                lnprop_comb{i}=struct('prop','Marker','val',{'o','s','d','^','v','>','<'});%to see different sims
            end
        else
            i=i+numel(nm);%skip numericals in the next loop
            nm=str2num(nm);
            if nm<=numel(model.eq)
                lns{grpi}=[lns{grpi} 11+numel(model.s)+numel(model.d)+nm];
                lnprop{11+numel(model.s)+numel(model.d)+nm}={'MarkerSize',10,'Marker','square'};
            else
                warning 'Given index for equation is too large'
            end
        end
        if numel(scale)>0
            scale={scale{:} grpi};
        end
    case '-'
        %create next group
        grpi=grpi+1;
        lns{grpi}=[];
    otherwise
    end
end

%cycle over simulations
for simi=simiset
    %cycle over groups
    for grpi=1:numel(lns)
    
        %data
        if combine==0
            %raw data
            data=reshape(model.sims(simi,:,lns{grpi}),[],numel(lns{grpi}));
            %scale given group
            if ismember(grpi,[scale{:}])
                %scale data
                datac=data;
                datac(sum(isnan(datac),2)>0,:)=[]; %delete rows with at least one nan
                datac(:,sum(datac==0)==size(datac,1))=ones(size(datac,1),sum(sum(datac==0)==size(datac,1),2));%replace zero columns with 1
                datac(sum(datac==0,2)>0,:)=[];     %delete rows with at least one 0
                multr=mean(abs(datac),1)/max(mean(abs(datac),1));
                multr(multr>.09)=1; %dont rescale anything 10 times or less smaller
                multr(multr<.01)=3*multr(multr<.01); %represent as 1/3 of size anything that is tiny
                if min(multr)>0
                    mult=ones(size(data,1),1)*multr;
                    data=data./mult;
                end
            end
        else
            %combined sims data
            data=[];
            multr=ones(1,numel(lns{grpi}));
            for i=simiset
                data1=reshape(model.sims(i,:,lns{grpi}),[],numel(lns{grpi}));
                data=[data data1];
                %scale given group
                if ismember(grpi,[scale{:}])
                    %scale data
                    datac=data1;
                    datac(sum(isnan(datac),2)>0,:)=[]; %delete rows with at least one nan
                    datac(:,sum(datac==0)==size(datac,1))=ones(size(datac,1),sum(sum(datac==0)==size(datac,1),2));%replace zero columns with 1
                    datac(sum(datac==0,2)>0,:)=[];     %delete rows with at least one 0
                    multr1=mean(abs(datac),1)/max(mean(abs(datac),1));
                    multr1(multr1>.09)=1; %dont rescale anything 10 times or less smaller
                    multr1(multr1<.01)=3*multr1(multr1<.01); %represent as 1/3 of size anything that is tiny
                    %combine
                    if min(multr1)>0
                        multr=min(multr,multr1);
                    end
                end
            end
            if ismember(grpi,[scale{:}])
                multr=repmat(multr,1,numel(simiset));
                mult=ones(size(data,1),1)*multr;
                data=data./mult;
            end
        end

        %figure and axes
        if numel(lns)==1
            %no groups
            %create figure and axes if no axes given, for each simi
            if isempty(ax) || simi>simiset(1)
                figure1=figure;
                set(figure1,'Name',[model.label ' sims ' sprintf('%g',simi)]);
                set(figure1,'Color',[1 1 1]);
                ax{1}=axes('Parent',figure1);
            elseif ~isempty(ax) && ~iscell(ax)
                ax{1}=ax;
            end
            title1=title(ax{1},titlestr);
            hold(ax{1},'all');
        else
            if grpi==1
                %if there are many groups create new figure with subplots
                figure1=figure;
                set(figure1,'Name',[model.label ' sims ' sprintf('%g',simi)]);
                set(figure1,'Color',[1 1 1]);
            end
            ax{grpi}=subplot(numel(lns),1,grpi,'Parent',figure1);
            hold(ax{grpi},'all');
            if grpi==1
                title1=title(ax{1},titlestr);
            end
        end
        
        %plot
        if model.T-model.t0<35
            set(ax{grpi},'XTick',1:model.T-model.t0+1,'XTickLabel',model.t0:model.T,'Xlim',[0.5 model.T-model.t0+1.5],'YGrid','on');
        else
            set(ax{grpi},'XTick',1:2:model.T-model.t0+1,'XTickLabel',model.t0:2:model.T,'Xlim',[0.5 model.T-model.t0+1.5],'YGrid','on');
        end

        %plot
        hlines=plot(ax{grpi},data,'MarkerFaceColor',[1 1 1]);

        %for empty data only draw axes
        if isempty(data)
            continue;
        end

        %decorate lines
        for ln=1:numel(hlines)
            i=mod(ln-1,numel(lns{grpi}))+1;%the index of the line in lns{grdi}
            j=floor((ln-1)/numel(lns{grpi}))+1;%the index of the sim in simiset
            if ismember(grpi,[scale{:}]) && abs(multr(i)-1)>1e-10 && combine==0
                set(hlines(ln),lnprop{lns{grpi}(i)}{:},'DisplayName',[model.simlabels{lns{grpi}(i)}(4:end) ' (scaled x' sprintf('%4.3f',multr(i)) ')']);
            elseif ~(ismember(grpi,[scale{:}]) && abs(multr(i)-1)>1e-10) && combine==0
                set(hlines(ln),lnprop{lns{grpi}(i)}{:},'DisplayName',model.simlabels{lns{grpi}(i)}(4:end));
            elseif ismember(grpi,[scale{:}]) && abs(multr(i)-1)>1e-10 && combine==1
                set(hlines(ln),lnprop{lns{grpi}(i)}{:},'DisplayName',['sim' sprintf('%d',j) ' : ' model.simlabels{lns{grpi}(i)}(4:end) ' (scaled x' sprintf('%4.3f',multr(i)) ')']);
                set(hlines(ln),'Color',get(hlines(ln),'Color')/(2^(j-1)));
                if ~isempty(lnprop_comb) && ~isempty(lnprop_comb{lns{grpi}(i)})
                    ii=mod(j-1,numel(lnprop_comb{lns{grpi}(i)}))+1;
                    set(hlines(ln),lnprop_comb{lns{grpi}(i)}(ii).prop,lnprop_comb{lns{grpi}(i)}(ii).val);
                end
            elseif ~(ismember(grpi,[scale{:}]) && abs(multr(i)-1)>1e-10) && combine==1
                set(hlines(ln),lnprop{lns{grpi}(i)}{:},'DisplayName',['sim' sprintf('%d',j) ' : ' model.simlabels{lns{grpi}(i)}(4:end)]);
                set(hlines(ln),'Color',get(hlines(ln),'Color')/(2^(j-1)));
                if ~isempty(lnprop_comb) && ~isempty(lnprop_comb{lns{grpi}(i)})
                    ii=mod(j-1,numel(lnprop_comb{lns{grpi}(i)}))+1;
                    set(hlines(ln),lnprop_comb{lns{grpi}(i)}(ii).prop,lnprop_comb{lns{grpi}(i)}(ii).val);
                end
            end
            %delete some decorations for many lines plots
            if numel(hlines)>10
                set(hlines,'Marker','none','LineWidth',1);
            end
        end

        %Ylim
        ylim=get(ax{grpi},'YLim');
        ylim(1)=min(min(data));
        set(ax{grpi},'YLim',ylim);
        
        %legend
        if numel(hlines)<=10
            legend1 = legend(ax{grpi},'show');
            %set(legend1,'Location','NorthWest','Box','off','Orientation','horizontal');
            set(legend1,'Location','NorthWest','Box','off');
        end
        
        %output as structure
        if numel(lns)==1
            %no groups case
            eval(['res.lines_s' sprintf('%g',simi) '=hlines;']);
        else
            eval(['res.lines_s' sprintf('%g',simi) '_b' sprintf('%g',grpi) '=hlines;']);
        end
    end
    %exit after one iteration over simi if combined
    if combine==1 break;
end

end%of function
