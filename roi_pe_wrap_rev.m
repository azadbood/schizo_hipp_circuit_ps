%% Univariate analysis
% Reads outputs of level 2 FEAT, pool across subjects, 
% Creates a figure with NN, OO, SS condition results
% Prepares an excel file for running ststs in R

clear

label = 'both'; % indx (first session) or both
r = 10; % which ROI to tun now
ccall = 1:4;%  run all 4 both conditions together for stats

% read ROI names
roi_fnames = {'DG','CA23','CA1','v1','loc','dlpfc','mpfc','rsp','perirhinal'};

% read subject and cond info
basepath = '/data/schizo/';
fid = fopen(fullfile(basepath,'subj_info.txt'));
data = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','HeaderLines',1,'CollectOutput',1);
data = data{:};
fclose(fid);

names = data(:,1);
subjectID = data(1:end,2);
session = data(1:end,5);
group = data(:,6);
runs = data(1:end,4);
runsss =data(1:end,7:10);
bothses = data(1:end,11);
conds = {'nn','oo','ss'};

conditions = {'patient ses1','patient ses2','control ses1','control ses2'};
indx{1} = intersect(find(strcmp(session,'ses-01')),find(strcmp(group,'patient'))); %patient ses1
indx{2} = intersect(find(strcmp(session,'ses-02')),find(strcmp(group,'patient'))); %patient ses2
indx{3} = intersect(find(strcmp(session,'ses-01')),find(strcmp(group,'control'))); %control ses1
indx{4} = intersect(find(strcmp(session,'ses-02')),find(strcmp(group,'control'))); %control ses2

both{1} = intersect(indx{1},find(strcmp(bothses,'1'))); %patient ses1 who have both sessions
both{2} = intersect(indx{2},find(strcmp(bothses,'1'))); %patient ses2 who have both sessions
both{3} = intersect(indx{3},find(strcmp(bothses,'1'))); %control ses1 who have both sessions
both{4} = intersect(indx{4},find(strcmp(bothses,'1'))); %control ses2 who have both sessions

thisgroup = both;

disp(['>>>>>>>>>>>>>>' roi_fnames{r} '<<<<<<<<<<<<<<<<<'])
pe_all = {};
pp = 1;
f=figure(round(rand(1)*1000));
f.Position = [100 100 400 350];
clf; set(gcf,'Color',[1 1 1]);
for cc = ccall
    thisind = thisgroup{cc};
    for s=1:length(thisind)
        sub=thisind(s);
        
        disp([names{sub} '__' session{sub}])
        % read data
        for c = 1:length(conds)
            thiscond = conds{c};
            filename = fullfile(basepath,'glm','level2',names{sub},session{sub},[thiscond '_cope' num2str(c) '.gfeat/cope' num2str(c) '.feat'],roi_fnames{r},'report.txt');
            
            fid = fopen(filename,'r');
            pereport = textscan(fid,'%d%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
            
            if ~ exist(filename)
                disp( [names{sub},session{sub},thiscond ' -- this condition still does not exist'])
            end
        end
        
        
    end
    alldata{pp} = pedata;
    
    
    %% prepare for stats
    sz = [length(thisgroup{cc})*3 5];
    varTypes = {'double','string','string','string','double'};
    varNames = {'subnum','session','group','condition','bvalues'};
    T{pp}= table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    
    T{pp}.subnum = [names(thisgroup{cc});names(thisgroup{cc});names(thisgroup{cc})];
    T{pp}.session = [session(thisgroup{cc});session(thisgroup{cc});session(thisgroup{cc})];
    T{pp}.group = [group(thisgroup{cc});group(thisgroup{cc});group(thisgroup{cc})];
    T{pp}.condition(1:length(thisgroup{cc})) = 'nn';
    T{pp}.condition(length(thisgroup{cc})+1:length(thisgroup{cc})*2) = 'oo';
    T{pp}.condition(length(thisgroup{cc})*2+1:length(thisgroup{cc})*3) = 'ss';
    T{pp}.bvalues = [alldata{pp}(:,1);alldata{pp}(:,2);alldata{pp}(:,3)];
    
    %% plot
    subplot(1,4,pp)
    title(conditions{cc})
    hold on
    b = bar([nanmean(pedata(:,1:3)) 0]);
    b.FaceColor = 'flat';
    b.EdgeColor = 'white';
    % define colors of choice
    b.CData(1,:) = [0, 0.4470, 0.7410];
    b.CData(2,:) = [0.8500, 0.3250, 0.0980];
    b.CData(3,:) = [0.9290, 0.6940, 0.1250];
    nnn = (255 - [0, 0.4470, 0.7410]*255)*1/2;
    ooo = (255 - [0.8500, 0.3250, 0.0980]*255)*1/2;
    sss = (255 - [0.9290, 0.6940, 0.1250]*255)*1/2;
    % jitter single subjects on graph location
    for s=1:length(thisind)
        jitt = (rand(1)-0.5)/8;
        plot(1+jitt,pedata(s,1),'bo', 'Color',[0, 0.4470, 0.7410] + nnn/255 ,'MarkerSize',6);
        hold on
        plot(2+jitt,pedata(s,2),'bo', 'Color',[0.8500, 0.3250, 0.0980] + ooo/255 ,'MarkerSize',6);
        hold on
        plot(3+jitt,pedata(s,3),'bo', 'Color',[0.9290, 0.6940, 0.1250] + sss/255,'MarkerSize',6)
    end
    
    set(gca, 'XTick', [1 2 3],'TickLength', [0 0], 'XTickLabel', {'nn' 'oo' 'ss'},'FontSize',9,'FontName','Helvetica','Box','off')
    ylabel('% signal change')
    sgtitle(roi_fnames{r})
    ylim([-0.4 0.4])
    
    hold on
    plot([0 1 2 3 4 5],[0 0 0 0 0 0],'-','color', [0 0 0],'LineWidth',0.5);
    
    hold on
    nnse = nanstd(pedata(:,1))/sqrt(length(pedata(:,1)));
    oose = nanstd(pedata(:,2))/sqrt(length(pedata(:,2)));
    ssse = nanstd(pedata(:,3))/sqrt(length(pedata(:,3)));
    
    e = errorbar(1:3,nanmean(pedata(:,1:3)),[nnse, oose, ssse],'o','color', [.35 .35 .35],'MarkerSize',4,'MarkerFaceColor',[0.35 0.35 0.35]);
    e.CapSize = 0;
    e.LineWidth = 1.5;
    length(pedata)
    pe_all{cc} = pedata;
    clear b xt gca nnse oose ssse sose ymax ymin pedata%roidata
    pp=pp+1;
end
TT = [];
for t = 1:length(T)
    TT = [TT;T{t}];
end
writetable(TT,fullfile('/Users/asieh/Documents/gitprojects/schizo_pattern/results/paper/newfigures/uniparts',['roi_uni_minusSS_finall_' roi_fnames{r} '_' label '.csv']))
clear TT alldata
