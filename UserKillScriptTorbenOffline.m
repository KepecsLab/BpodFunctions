%Behavior analysis Olf2AFC
function UserKillScriptTorbenOffline()
basepath = 'C:\Data\DataPostdoc\BpodData\';
animal = 'TP50';
protocol = 'Dual2AFC';
Sessions = getDir(fullfile(basepath,animal,protocol,'Session Data'),'file',animal);
% SessionDataPath={'K:\BpodData\TP44\Dual2AFC\Session Data\TP44_Dual2AFC_Nov05_2019_Session1.mat';...
%     'K:\BpodData\TP44\Dual2AFC\Session Data\TP44_Dual2AFC_Nov06_2019_Session1.mat';...
%     'K:\BpodData\TP44\Dual2AFC\Session Data\TP44_Dual2AFC_Nov07_2019_Session1.mat';...
%     'K:\BpodData\TP44\Dual2AFC\Session Data\TP44_Dual2AFC_Nov08_2019_Session1.mat';}
if ~iscell(Sessions)
Sessions={Sessions};
end
AllData=struct('Custom',struct()); AllData.Custom.OriginalTrialIndex = []; AllData.Custom.SessionIndex=[];
for s=1:length(Sessions)
    A=load(fullfile(basepath,animal,protocol,'Session Data',Sessions{s}));
    nT=A.SessionData.nTrials-1;
    
    %filter sessions
    set = A.SessionData.Settings.GUI;
    if set.LaserTrials <0.1 || set.PercentCatch<0.001 || nT<200 || set.CatchError<1 || set.FeedbackDelayMax<7.9
        continue
    end
    
    ff=fieldnames(A.SessionData.Custom);
    for f=1:length(ff)
        if size(A.SessionData.Custom.(ff{f}),2)>=nT
            if isfield(AllData.Custom,ff{f})
            AllData.Custom.(ff{f})= [ AllData.Custom.(ff{f}),  A.SessionData.Custom.(ff{f})(1,1:nT)];
            else
                AllData.Custom.(ff{f})=  A.SessionData.Custom.(ff{f})(1,1:nT);
            end
        end
    end

    AllData.Custom.OriginalTrialIndex=[AllData.Custom.OriginalTrialIndex,1:nT];
    AllData.Custom.SessionIndex=[AllData.Custom.SessionIndex,s*ones(1,nT)];
end

AllData.Custom.ChoiceLeft(AllData.Custom.OriginalTrialIndex<50)=NaN;
AllData.nTrials=size(AllData.Custom.ChoiceLeft,2)+1;

global TaskParameters
global BpodSystem

BpodSystem.Data=AllData;
%use last session for rest
SessionData=A.SessionData;
BpodSystem.Animal=SessionData.Custom.Subject;
BpodSystem.DataPath =fullfile(basepath,animal,protocol,'Session Data',Sessions{1});

TaskParameters=SessionData.TrialSettings(end-1);


%load mail settings --> contains mail address & password & evernote e-mail address
load('MailSettings.mat');%loads MailSettings struct

%evernote mail address
MailAddress = MailSettings.EvernoteMail;

%save figure
FigureFolder = fullfile(fileparts(fileparts(BpodSystem.DataPath)),'Session Figures');
% FigureHandle = BpodSystem.GUIHandles.OutcomePlot.HandleOutcome.Parent;
% FigureString = get(FigureHandle,'Name');
% 
if ~isdir(FigureFolder)
    mkdir(FigureFolder);
end

% FigurePath = fullfile(FigureFolder,[FigureName,'.png']);
% saveas(FigureHandle,FigurePath,'png');

%Analysis
[~, FigureName] = fileparts(BpodSystem.DataPath);
try
    FigAnalysis = Analysis();
    FigurePathAnalysis = fullfile(FigureFolder,[FigureName,'Analysis.png']);
    saveas(FigAnalysis,FigurePathAnalysis,'png');
%     close(FigAnalysis);
    DidAnalysis = true;
catch
    DidAnalysis = false;
end

%send email

[~,sessionfile] = fileparts(BpodSystem.DataPath);
animal = SessionData.Custom.Subject;

Subject = strcat(sessionfile,'@',animal);
Body = sessionfile;

if DidAnalysis
    Attachment = {FigurePathAnalysis};
else
    Attachment = {};
end

sent = SendMyMail(MailSettings,MailAddress,Subject,Body,Attachment);

if sent
    fprintf('Figure "%s" sent to %s.\n','',MailAddress);
else
    fprintf('Error:SendFigureTo:Mail could not be sent to %s.\n',MailAddress);
end

%% copy data to server
% try
%     
%     %%%%%
%     %os
%     os = getenv('OS');
%     if strcmpi(os(1:min(7,length(os))),'windows')
%         servername = '\\uncertainty.cshl.edu\home'; %new uncertanity server (8/2018) works with home only
% %         user = strcat(getenv('username'));
%         user ='';
%     else
%         servername = '/media/';
%         user='torben';
%     end
%     [~,subject] = fileparts(fileparts(fileparts(fileparts(BpodSystem.DataPath))));
%     if ~isdir(fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Data'))
%         mkdir(fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Data'));
%     end
%     if ~isdir(fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Settings'))
%         mkdir(fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Settings'));
%     end
%     copyfile(BpodSystem.DataPath,fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Data'));
%     copyfile(BpodSystem.SettingsPath,fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Settings'));
% catch
%     fprintf('Error copying data to server. Files not copied!\n');
% end

end

function FigHandle = Analysis()


global TaskParameters
global BpodSystem

offline=false;
if offline
    BpodSystem.Data=SessionData;
    TaskParameters = SessionData.Settings;
    Animal ='Unknown';
else
    [~,Animal]=fileparts(fileparts(fileparts(fileparts(BpodSystem.DataPath))));
end

GracePeriodsMax = TaskParameters.GUI.FeedbackDelayGrace; %assumes same for each trial
StimTime = TaskParameters.GUI.AuditoryStimulusTime; %assumes same for each trial
MinWT = 2; %assumes same for each trial
MaxWT = 10;
AudBin = 8; %Bins for psychometric
AudBinWT = 6;%Bins for vevaiometric

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nTrials=BpodSystem.Data.nTrials;
DV = BpodSystem.Data.Custom.DV(1:nTrials-1);
ChoiceLeft = BpodSystem.Data.Custom.ChoiceLeft(1:nTrials-1);
%exclude first trials?
ChoiceLeft(1:min([nTrials,max([10,TaskParameters.GUI.StartEasyTrials])]))=NaN;
ST = BpodSystem.Data.Custom.ST(1:nTrials-1);
CatchTrial = BpodSystem.Data.Custom.CatchTrial((1:nTrials-1));
Feedback = BpodSystem.Data.Custom.Feedback(1:nTrials-1);
Correct = BpodSystem.Data.Custom.ChoiceCorrect(1:nTrials-1);
WT =  BpodSystem.Data.Custom.FeedbackTime(1:nTrials-1);


%correct WT?
% WT(WT<MinWT | WT>MaxWT)=NaN;
% meanL = nanmedian(WT(CatchTrial==1  & ChoiceLeft==1));
% meanR = nanmedian(WT(CatchTrial==1  & ChoiceLeft==0));
% grandM = nanmedian(WT(CatchTrial==1));
% WT(ChoiceLeft==1) = WT(ChoiceLeft==1) - meanL + grandM+1;
% WT(ChoiceLeft==0) = WT(ChoiceLeft==0) - meanR + grandM+1;
% WT(WT<2.5)=NaN;
% MinWT=0;
% MaxWT=12;

if isfield(BpodSystem.Data.Custom,'LaserTrial')
    LaserTrial =  BpodSystem.Data.Custom.LaserTrial(1:nTrials-1);
else
    LaserTrial=false(1,nTrials);
end
%define "completed trial"
% not that abvious for errors
%Correct vector is 1 for correct choice, 0 for incorrect choice, nan
%for no choice
%now: correct --> received feedback (reward)
%     error --> always (if error catch)
%     catch --> always (if choice happened)

CompletedTrials = (Feedback&Correct==1) | (Correct==0) | CatchTrial&~isnan(ChoiceLeft);
nTrialsCompleted = sum(CompletedTrials);

%calculate exerienced dv
ExperiencedDV=zeros(1,length(ST));

%click task
if TaskParameters.GUI.AuditoryStimulusType == 1 %click
    LeftClickTrain = BpodSystem.Data.Custom.LeftClickTrain(1:nTrials-1);
    RightClickTrain = BpodSystem.Data.Custom.RightClickTrain(1:nTrials-1);
    for t = 1 : length(ST)
        R = BpodSystem.Data.Custom.RightClickTrain{t};
        L = BpodSystem.Data.Custom.LeftClickTrain{t};
        Ri = sum(R<=ST(t));if Ri==0, Ri=1; end
        Li = sum(L<=ST(t));if Li==0, Li=1; end
        ExperiencedDV(t) = log10(Li/Ri);
        %         ExperiencedDV(t) = (Li-Ri)./(Li+Ri);
    end
elseif TaskParameters.GUI.AuditoryStimulusType == 2 %freq
    LevelsLow = 1:ceil(TaskParameters.GUI.Aud_nFreq/3);
    LevelsHigh = ceil(TaskParameters.GUI.Aud_nFreq*2/3)+1:TaskParameters.GUI.Aud_nFreq;
    AudCloud = BpodSystem.Data.Custom.AudCloud(1:nTrials-1);
    for t = 1 : length(ST)
        NLow = sum(ismember(AudCloud{t},LevelsLow)); if NLow==0, NLow=1; end
        NHigh = sum(ismember(AudCloud{t},LevelsHigh)); if NHigh==0, NHigh=1; end
        ExperiencedDV(t) = log10(NHigh/NLow);
    end
end

%caclulate  grace periods
GracePeriods=[];
GracePeriodsL=[];
GracePeriodsR=[];
% for t = 1 : length(ST)
%     GracePeriods = [GracePeriods;BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Rin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Rin_grace(:,1);BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Lin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Lin_grace(:,1)];
%     if ChoiceLeft(t) == 1
%         GracePeriodsL = [GracePeriodsL;BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Lin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Lin_grace(:,1)];
%     elseif ChoiceLeft(t)==0
%         GracePeriodsR = [GracePeriodsR;BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Rin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Rin_grace(:,1)];
%     end
% end


CompletedTrials = CompletedTrials==1;
CatchTrial = CatchTrial==1;

%laser trials?
if sum(LaserTrial)>0
    LaserCond = [false;true];
else
    LaserCond=false;
end
CondColors={[0,0,0],[.9,.1,.1]};

%%
FigHandle = figure('Position',[ 360         187        1056         598],'NumberTitle','off','Name',Animal,'Color',[1,1,1]);
% ExperiencedDV=DV;

%Psychometric
subplot(3,4,1)
hold on
for i = 1:length(LaserCond)
    CompletedTrialsCond = CompletedTrials & LaserTrial == LaserCond(i);
    AudDV = ExperiencedDV(CompletedTrialsCond);
    if ~isempty(AudDV)
        BinIdx = discretize(AudDV,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBin+1));
        PsycY = grpstats(ChoiceLeft(CompletedTrialsCond),BinIdx,'mean');
        PsycX = grpstats(ExperiencedDV(CompletedTrialsCond),BinIdx,'mean');
        plot(PsycX,PsycY,'ok','MarkerFaceColor',CondColors{i},'MarkerEdgeColor','w','MarkerSize',6)
        XFit = linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100);
        YFit = glmval(glmfit(AudDV,ChoiceLeft(CompletedTrialsCond)','binomial'),linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100),'logit');
        plot(XFit,YFit,'Color',CondColors{i});
        xlabel('DV');ylabel('p left')
        text(0.95*min(get(gca,'XLim')),0.96*max(get(gca,'YLim'))-(i-1)*.1,[num2str(round(nanmean(Correct(CompletedTrialsCond))*100)),'%,n=',num2str(sum(CompletedTrialsCond))]);
    end
end

%conditioned psychometric
subplot(3,4,2)
hold on
%low
WTmed=median(WT(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT));
AudDV = ExperiencedDV(CompletedTrials&CatchTrial&WT<=WTmed&WT>MinWT);
if ~isempty(AudDV)
    ChoiceLeftadj = ChoiceLeft(CompletedTrials&CatchTrial&WT<=WTmed&WT>MinWT);
    BinIdx = discretize(AudDV,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBin+1));
    PsycY = grpstats(ChoiceLeftadj,BinIdx,'mean');
    PsycX = grpstats(AudDV,BinIdx,'mean');
    h1=plot(PsycX,PsycY,'ok','MarkerFaceColor',[.5,.5,.5],'MarkerEdgeColor','w','MarkerSize',6);
    XFit = linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100);
    YFit = glmval(glmfit(AudDV,ChoiceLeftadj','binomial'),linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100),'logit');
    plot(XFit,YFit,'Color',[.5,.5,.5]);
    %high
    AudDV = ExperiencedDV(CompletedTrials&CatchTrial&WT>WTmed&WT<MaxWT);
    ChoiceLeftadj = ChoiceLeft(CompletedTrials&CatchTrial&WT>WTmed&WT<MaxWT);
    BinIdx = discretize(AudDV,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBin+1));
    PsycY = grpstats(ChoiceLeftadj,BinIdx,'mean');
    PsycX = grpstats(AudDV,BinIdx,'mean');
    h2=plot(PsycX,PsycY,'ok','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',6);
    XFit = linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100);
    YFit = glmval(glmfit(AudDV,ChoiceLeftadj','binomial'),linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100),'logit');
    plot(XFit,YFit,'k');
    xlabel('DV');ylabel('p left')
    legend([h2,h1],{['WT>',num2str(round(WTmed*100)/100)],['WT<',num2str(round(WTmed*100)/100)]},'Units','normalized','Position',[0.333,0.85,0.1,0.1])
end

%calibration
subplot(3,4,3)
hold on
xlabel('Waiting time (s)');ylabel('p correct')
WTBin=5;
WTEdges=linspace(min(WT(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT))-10*eps,max(WT(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT))+10*eps,WTBin+1);
ColorsCorrect = {[.1,.9,.1],[.1,.8,.6]};
ColorsError = {[.9,.1,.1],[.9,.1,.6]};

for i =1:length(LaserCond)
    WTCatch = WT(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT & LaserTrial==LaserCond(i));
    if ~isempty(WTCatch)
        BinIdx = discretize(WTCatch,WTEdges);
        WTX = grpstats(WTCatch,BinIdx,'mean');
        PerfY = grpstats(Correct(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i)),BinIdx,'mean');
        plot(WTX,PerfY,'Color',CondColors{i},'LineWidth',2);
        [r,p]=corr(WTCatch',Correct(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i))','type','Spearman');
        text(min(get(gca,'XLim'))+0.05,max(get(gca,'YLim'))-0.07*i,['r=',num2str(round(r*100)/100),', p=',num2str(round(p*100)/100)],'Color',CondColors{i});
    end
end


%Vevaiometric
subplot(3,4,4)
hold on
xlabel('DV');ylabel('Waiting time (s)')
AudDV = ExperiencedDV(CompletedTrials&CatchTrial&WT<MaxWT&WT>MinWT);
DVEdge = linspace(-max(AudDV)-10*eps,max(AudDV)+10*eps,AudBinWT+1);
Rcatch=cell(1,2);Pcatch=cell(1,2);Rerror=cell(1,2);Perror=cell(1,2);
%for confidence auc
auc = nan(AudBinWT/2,length(LaserCond),1);
auc_sem = nan(AudBinWT/2,length(LaserCond),1);
for i =1:length(LaserCond)
    
    WTCatch = WT(CompletedTrials&CatchTrial&Correct==1&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i));
    DVCatch = ExperiencedDV(CompletedTrials&CatchTrial&Correct==1&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i));
     if ~isempty(DVCatch)
        BinIdxCatch = discretize(DVCatch,DVEdge);
        if ~all(isnan(BinIdxCatch))
            WTCatchY = grpstats(WTCatch,BinIdxCatch,'mean');
            DVCatchX = grpstats(DVCatch,BinIdxCatch,'mean');
            plot(DVCatchX,WTCatchY,'Color',ColorsCorrect{i},'LineWidth',2)
        end
        WTError = WT(CompletedTrials&Correct==0&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i));
        DVError = ExperiencedDV(CompletedTrials&Correct==0&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i));
        BinIdxError = discretize(DVError,DVEdge);
        if ~all(isnan(BinIdxError))
            WTErrorY = grpstats(WTError,BinIdxError,'mean');
            DVErrorX = grpstats(DVError,BinIdxError,'mean');
            plot(DVErrorX,WTErrorY,'Color',ColorsError{i},'LineWidth',2)
        end
        
%         plot(DVCatch,WTCatch,'o','MarkerSize',2,'MarkerFaceColor',ColorsCorrect{i},'Color',ColorsCorrect{i})
%         plot(DVError,WTError,'o','MarkerSize',2,'MarkerFaceColor',ColorsError{i},'Color',ColorsError{i})
        legend('Correct Catch','Error','Location','best')
        %evaluate vevaiometric
        [Rc,Pc] = EvaluateVevaiometric(DVCatch,WTCatch);
        [Re,Pe] = EvaluateVevaiometric(DVError,WTError);
        Rcatch{i}=Rc;Pcatch{i}=Pc;Rerror{i}=Re;Perror{i}=Pe;
        %confidence auc
%         BinIdxCatch(ismember(BinIdxCatch,[1,2,3]))=NaN; %hack to only
%         look at one side
%         BinIdxError(ismember(BinIdxError,[4,5,6]))=NaN;
        for d=1:(AudBinWT/2)
            BinIdxCatch(BinIdxCatch==AudBinWT-d+1)=d;
            BinIdxError(BinIdxError==AudBinWT-d+1)=d;
        end
            
        for d=1:(AudBinWT/2)
        [auc(d,i),~,auc_sem(d,i)] = rocarea_torben(WTCatch(BinIdxCatch==d),WTError(BinIdxError==d),'bootstrap',200);
        end
        
    end
end
for i =1:length(LaserCond)
    if ~isempty(Rcatch{i}) && ~isempty(Pcatch{i}) && ~isempty(Rerror{i}) && ~isempty(Perror{i})
        unit = max(get(gca,'YLim'))-min(get(gca,'YLim'));
        text(max(get(gca,'XLim'))+0.03,max(get(gca,'YLim'))-unit*(0.1+(i-1)*.5),['r_l=',num2str(round(Rcatch{i}(1)*100)/100),' r_r=',num2str(round(Rcatch{i}(2)*100)/100)],'Color',ColorsCorrect{i});
        text(max(get(gca,'XLim'))+0.03,max(get(gca,'YLim'))-unit*(0.2+(i-1)*.5),['r=',num2str(round(Rcatch{i}(3)*100)/100),', p=',num2str(round(Pcatch{i}(3)*100)/100)],'Color',ColorsCorrect{i});
        text(max(get(gca,'XLim'))+0.03,max(get(gca,'YLim'))-unit*(0.3+(i-1)*.5),['r_l=',num2str(round(Rerror{i}(1)*100)/100),' r_r=',num2str(round(Rerror{i}(2)*100)/100)],'Color',ColorsError{i});
        text(max(get(gca,'XLim'))+0.03,max(get(gca,'YLim'))-unit*(0.4+(i-1)*.5),['r=',num2str(round(Rerror{i}(3)*100)/100),', p=',num2str(round(Perror{i}(3)*100)/100)],'Color',ColorsError{i});
    end
end


%waiting time distributions
ColorsCond = {[.5,.5,.5],[.9,.1,.1]};
WTBins = linspace(min(WT(~Feedback)),max(WT(~Feedback)),21);
if length(LaserCond)==1
    %no laser
    subplot(3,4,5)
    hold on
    xlabel('waiting time (s)'); ylabel ('n trials');
    WTnoFeedbackL = WT(~Feedback & ChoiceLeft == 1);
    WTnoFeedbackR = WT(~Feedback & ChoiceLeft == 0);
    histogram(WTnoFeedbackL,WTBins,'EdgeColor','none','FaceColor',[.2,.2,1]);
    histogram(WTnoFeedbackR,WTBins,'EdgeColor','none','FaceColor',[.8,.6,.1]);

    meanWTL = nanmean(WTnoFeedbackL);
    meanWTR = nanmean(WTnoFeedbackR);
    line([meanWTL,meanWTL],get(gca,'YLim'),'Color',[.2,.2,1]);
    line([meanWTR,meanWTR],get(gca,'YLim'),'Color',[.8,.6,.1]);
    text(meanWTL-1,1.05*(max(get(gca,'YLim'))-min(get(gca,'YLim'))),['m_l=',num2str(round(meanWTL*10)/10)],'Color',[.2,.2,1]);
    text(meanWTL-1,1.15*(max(get(gca,'YLim'))-min(get(gca,'YLim'))),['m_r=',num2str(round(meanWTR*10)/10)],'Color',[.8,.6,.1]);
    
    PshortWTL = sum(WTnoFeedbackL<MinWT)/sum(~isnan(WTnoFeedbackL));
    PshortWTR = sum(WTnoFeedbackR<MinWT)/sum(~isnan(WTnoFeedbackR));
    text(max(get(gca,'XLim'))+0.03,0.85*(max(get(gca,'YLim'))-min(get(gca,'YLim')))+min(get(gca,'YLim')),['L_{2}=',num2str(round(PshortWTL*100)/100),', R_{2}=',num2str(round(PshortWTR*100)/100)],'Color',[0,0,0]);
    
else%laser
    subplot(3,4,5)
    hold on
    xlabel('waiting time (s)'); ylabel ('p');
    subplot(3,4,6)
    hold on
    xlabel('waiting time (s)'); ylabel ('p');
    PshortWTL=cell(1,2);PshortWTR=cell(1,2);
    for i =1:length(LaserCond)
    
    WTnoFeedbackL = WT(~Feedback & ChoiceLeft == 1 & LaserTrial==LaserCond(i));
    WTnoFeedbackR = WT(~Feedback & ChoiceLeft == 0 & LaserTrial==LaserCond(i));
     meanWTL = nanmean(WTnoFeedbackL);
    meanWTR = nanmean(WTnoFeedbackR);
    subplot(3,4,5)
    histogram(WTnoFeedbackL,WTBins,'EdgeColor','none','FaceColor',ColorsCond{i},'Normalization','probability');
    line([meanWTL,meanWTL],get(gca,'YLim'),'Color',ColorsCond{i});
    text(meanWTL-1,(1.05-0.1*(i-1))*(max(get(gca,'YLim'))-min(get(gca,'YLim'))),['m_l=',num2str(round(meanWTL*10)/10)],'Color',ColorsCond{i});
    subplot(3,4,6)
    histogram(WTnoFeedbackR,WTBins,'EdgeColor','none','FaceColor',ColorsCond{i},'Normalization','probability');
    line([meanWTR,meanWTR],get(gca,'YLim'),'Color',ColorsCond{i});
    text(meanWTL-1,(1.05-0.1*(i-1))*(max(get(gca,'YLim'))-min(get(gca,'YLim'))),['m_r=',num2str(round(meanWTR*10)/10)],'Color',ColorsCond{i});

    PshortWTL{i} = sum(WTnoFeedbackL<MinWT)/sum(~isnan(WTnoFeedbackL));
    PshortWTR{i} = sum(WTnoFeedbackR<MinWT)/sum(~isnan(WTnoFeedbackR));
    
    end
    for i =1:length(LaserCond)
        text(max(get(gca,'XLim'))+0.03,(0.85/i)*(max(get(gca,'YLim'))-min(get(gca,'YLim')))+min(get(gca,'YLim')),['L_{2}=',num2str(round(PshortWTL{i}*100)/100),', R_{2}=',num2str(round(PshortWTR{i}*100)/100)],'Color',ColorsCond{i});
    end
end

%confidence index
subplot(3,4,8)
hold on
for i =1:length(LaserCond)
    errorbar(1:size(auc,1),auc(end:-1:1,i),auc_sem(end:-1:1,i),'-o','MarkerFaceColor',CondColors{i},'MarkerEdgeColor',CondColors{i},'LineWidth',2,'Color',CondColors{i})
end
xlabel('DV quantile')
ylabel('AUC')


RedoTicks(gcf);

end

function RedoTicks(h)
Chil=get(h,'Children');

for i = 1:length(Chil)
    if strcmp(Chil(i).Type,'axes')
        
        set(Chil(i),'TickDir','out','TickLength',[0.03 0.03],'box','off')
        
    end
    
    if strcmp(Chil(i).Type,'legend')
        set(Chil(i),'box','off')
    end
end

end


% sends mail from Torben's cshl gmail account
% 3 or  4 inputs: address,subject,message,cell with attachment paths
% (each as string)
function sent = SendMyMail(varargin)
sent = false;
MailSettings = varargin{1};
setpref('Internet','E_mail',MailSettings.MailFrom)
setpref('Internet','SMTP_Server','smtp.gmail.com')
setpref('Internet','SMTP_Username',MailSettings.MailFrom)
setpref('Internet','SMTP_Password',MailSettings.MailFromPassword)
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

if length(varargin)==4
    try
        sendmail(varargin{2},varargin{3},varargin{4})
        sent=true;
    catch
        display('Error:SendMyMail:E-Mail could not be sent.')
    end
elseif length(varargin)==5
    try
        %attachments need to be in full path (not ~) for linux systems
        for k =1:length(varargin{5})
            if strcmp(varargin{5}{k}(1),'~')
                varargin{5}{k} = fullfile('/home/torben',varargin{5}{k}(2:end));
            end
        end
        
        sendmail(varargin{2},varargin{3},varargin{4},varargin{5})
        sent=true;
    catch
        display('Error:SendMyMail:E-Mail could not be sent.')
    end
else
    display('Error:SendMyMail:Number of input arguments wrong.')
end

end

function [R,P] = EvaluateVevaiometric(DV,WT)
R = zeros(1,3);
P=zeros(1,3);
if sum(DV<=0)>0
    [R(1),P(1)] = corr(DV(DV<=0)',WT(DV<=0)','type','Spearman');
end
if sum(DV>0)>0
    [R(2),P(2)] = corr(DV(DV>0)',WT(DV>0)','type','Spearman');
end
if sum(~isnan(DV))>0
    [R(3),P(3)] = corr(abs(DV)',WT','type','Spearman');
end
end