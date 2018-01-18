%Behavior analysis Olf2AFC
function UserKillScriptTorben

global BpodSystem

%load mail settings --> contains mail address & password & evernote e-mail address
load('MailSettings.mat');%loads MailSettings struct

%evernote mail address
MailAddress = MailSettings.EvernoteMail;

%save figure
FigureFolder = fullfile(fileparts(fileparts(BpodSystem.DataPath)),'Session Figures');
FigureHandle = BpodSystem.GUIHandles.OutcomePlot.HandleOutcome.Parent;
FigureString = get(FigureHandle,'Name');
[~, FigureName] = fileparts(BpodSystem.DataPath);
if ~isdir(FigureFolder)
    mkdir(FigureFolder);
end

FigurePath = fullfile(FigureFolder,[FigureName,'.png']);
saveas(FigureHandle,FigurePath,'png');

%Analysis
try
    FigAnalysis = Analysis();
    FigurePathAnalysis = fullfile(FigureFolder,[FigureName,'Analysis.png']);
    saveas(FigAnalysis,FigurePathAnalysis,'png');
    close(FigAnalysis);
    DidAnalysis = true;
catch
    DidAnalysis = false;
end

%send email

[x,sessionfile] = fileparts(BpodSystem.DataPath);
[~,animal] = fileparts(fileparts(fileparts(x)));

Subject = strcat(sessionfile,'@',animal);
Body = sessionfile;

if DidAnalysis
    Attachment = {FigurePath,FigurePathAnalysis};
else
    Attachment = {FigurePath};
end

sent = SendMyMail(MailSettings,MailAddress,Subject,Body,Attachment);

if sent
    fprintf('Figure "%s" sent to %s.\n',FigureString,MailAddress);
else
    fprintf('Error:SendFigureTo:Mail could not be sent to %s.\n',MailAddress);
end

%% copy data to server
try
    
    %%%%%
    %os
    os = getenv('OS');
    if strcmpi(os(1:min(7,length(os))),'windows')
        servername = '\\kepecsdata\homes';
        user = strcat(getenv('username'));
    else
        servername = '/media/';
        user='torben';
    end
    [~,subject] = fileparts(fileparts(fileparts(fileparts(BpodSystem.DataPath))));
    if ~isdir(fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Data'))
        mkdir(fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Data'));
    end
    if ~isdir(fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Settings'))
        mkdir(fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Settings'));
    end
    copyfile(BpodSystem.DataPath,fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Data'));
    copyfile(BpodSystem.SettingsPath,fullfile(servername,user,'BpodData',subject,BpodSystem.CurrentProtocolName,'Session Settings'));
catch
    fprintf('Error copying data to server. Files not copied!\n');
end

end

function FigHandle = Analysis()

global TaskParameters
global BpodSystem

GracePeriodsMax = TaskParameters.GUI.FeedbackDelayGrace; %assumes same for each trial
StimTime = TaskParameters.GUI.AuditoryStimulusTime; %assumes same for each trial
MinWT = TaskParameters.GUI.VevaiometricMinWT; %assumes same for each trial
MaxWT = 10;
AudBin = 8; %Bins for psychometric
AudBinWT = 6;%Bins for vevaiometric
windowCTA = 150; %window for CTA (ms)

[~,Animal]=fileparts(fileparts(fileparts(fileparts(BpodSystem.DataPath))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nTrials=BpodSystem.Data.nTrials;
DV = BpodSystem.Data.Custom.DV(1:nTrials-1);
ChoiceLeft = BpodSystem.Data.Custom.ChoiceLeft(1:nTrials-1);
ST = BpodSystem.Data.Custom.ST(1:nTrials-1);
CatchTrial = BpodSystem.Data.Custom.CatchTrial((1:nTrials-1));
Feedback = BpodSystem.Data.Custom.Feedback(1:nTrials-1);
Correct = BpodSystem.Data.Custom.ChoiceCorrect(1:nTrials-1);
WT =  BpodSystem.Data.Custom.FeedbackTime(1:nTrials-1);
LaserTrial =  BpodSystem.Data.Custom.LaserTrial(1:nTrials-1);

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
        Ri = find(R>ST(t),1,'first'); if isempty(Ri), Ri=1; end
        Li = find(L>ST(t),1,'first');if isempty(Li), Li=1; end
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
for t = 1 : length(ST)
    GracePeriods = [GracePeriods;BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Rin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Rin_grace(:,1);BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Lin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Lin_grace(:,1)];
    if ChoiceLeft(t) == 1
        GracePeriodsL = [GracePeriodsL;BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Lin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Lin_grace(:,1)];
    elseif ChoiceLeft(t)==0
        GracePeriodsR = [GracePeriodsR;BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Rin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Rin_grace(:,1)];
    end
end


CompletedTrials = CompletedTrials==1;
CatchTrial = CatchTrial==1;

%%
FigHandle = figure('Position',[ 360         187        1056         431],'NumberTitle','off','Name',Animal);
% ExperiencedDV=DV;
%Psychometric
subplot(2,4,1)
hold on
AudDV = ExperiencedDV(CompletedTrials);
if ~isempty(AudDV)
    BinIdx = discretize(AudDV,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBin+1));
    PsycY = grpstats(ChoiceLeft(CompletedTrials),BinIdx,'mean');
    PsycX = grpstats(ExperiencedDV(CompletedTrials),BinIdx,'mean');
    plot(PsycX,PsycY,'ok','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',6)
    XFit = linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100);
    YFit = glmval(glmfit(AudDV,ChoiceLeft(CompletedTrials)','binomial'),linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100),'logit');
    plot(XFit,YFit,'k');
    xlabel('DV');ylabel('p left')
    text(0.95*min(get(gca,'XLim')),0.96*max(get(gca,'YLim')),[num2str(round(nanmean(Correct(CompletedTrials))*100)),'%,n=',num2str(nTrialsCompleted)]);
end

%conditioned psychometric
subplot(2,4,2)
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
subplot(2,4,3)
hold on
xlabel('Waiting time (s)');ylabel('p correct')
WTBin=5;
CondColors={[0,0,0],[.9,.1,.1]};
ColorsCorrect = {[.1,.9,.1],[.1,.8,.6]};
ColorsError = {[.9,.1,.1],[.9,.1,.6]};
if sum(LaserTrial)>0
    LaserCond = [false;true];
else
    LaserCond=false;
end
for i =1:length(LaserCond)
    WTCatch = WT(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT & LaserTrial==LaserCond(i));
    if ~isempty(WTCatch)
        BinIdx = discretize(WTCatch,linspace(min(WTCatch)-10*eps,max(WTCatch)+10*eps,WTBin+1));
        WTX = grpstats(WTCatch,BinIdx,'mean');
        PerfY = grpstats(Correct(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i)),BinIdx,'mean');
        plot(WTX,PerfY,'Color',CondColors{i},'LineWidth',2);
        [r,p]=corr(WTCatch',Correct(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i))','type','Spearman');
        text(min(get(gca,'XLim'))+0.05,max(get(gca,'YLim'))-0.07*i,['r=',num2str(round(r*100)/100),', p=',num2str(round(p*100)/100)],'Color',CondColors{i});
    end
end


%Vevaiometric
subplot(2,4,4)
hold on
xlabel('DV');ylabel('Waiting time (s)')
AudDV = ExperiencedDV(CompletedTrials&CatchTrial&WT<MaxWT&WT>MinWT);
Rcatch=cell(1,2);Pcatch=cell(1,2);Rerror=cell(1,2);Perror=cell(1,2);
for i =1:length(LaserCond)
    
    WTCatch = WT(CompletedTrials&CatchTrial&Correct==1&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i));
    DVCatch = ExperiencedDV(CompletedTrials&CatchTrial&Correct==1&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i));
    if ~isempty(DVCatch)
        BinIdx = discretize(DVCatch,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBinWT+1));
        if ~all(isnan(BinIdx))
            WTCatchY = grpstats(WTCatch,BinIdx,'mean');
            DVCatchX = grpstats(DVCatch,BinIdx,'mean');
            plot(DVCatchX,WTCatchY,'g','LineWidth',2)
        end
        WTError = WT(CompletedTrials&Correct==0&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i));
        DVError = ExperiencedDV(CompletedTrials&Correct==0&WT>MinWT&WT<MaxWT  & LaserTrial==LaserCond(i));
        BinIdx = discretize(DVError,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBinWT+1));
        if ~all(isnan(BinIdx))
            WTErrorY = grpstats(WTError,BinIdx,'mean');
            DVErrorX = grpstats(DVError,BinIdx,'mean');
            plot(DVErrorX,WTErrorY,'r','LineWidth',2)
        end
        
        plot(DVCatch,WTCatch,'o','MarkerSize',2,'MarkerFaceColor',ColorsCorrect{i},'Color',ColorsCorrect{i})
        plot(DVError,WTError,'o','MarkerSize',2,'MarkerFaceColor',ColorsError{i},'Color',ColorsError{i})
        legend('Correct Catch','Error','Location','best')
        %evaluate vevaiometric
        [Rc,Pc] = EvaluateVevaiometric(DVCatch,WTCatch);
        [Re,Pe] = EvaluateVevaiometric(DVError,WTError);
        Rcatch{i}=Rc;Pcatch{i}=Pc;Rerror{i}=Re;Perror{i}=Pe;
        
    end
end
for i =1:length(LaserCond)
    unit = max(get(gca,'YLim'))-min(get(gca,'YLim'));
    text(max(get(gca,'XLim'))+0.03,max(get(gca,'YLim'))-unit*(0.1+(i-1)*.5),['r_l=',num2str(round(Rcatch{i}(1)*100)/100),' r_r=',num2str(round(Rcatch{i}(2)*100)/100)],'Color',ColorsCorrect{i});
    text(max(get(gca,'XLim'))+0.03,max(get(gca,'YLim'))-unit*(0.2+(i-1)*.5),['r=',num2str(round(Rcatch{i}(3)*100)/100),', p=',num2str(round(Pcatch{i}(3)*100)/100)],'Color',ColorsCorrect{i});
    text(max(get(gca,'XLim'))+0.03,max(get(gca,'YLim'))-unit*(0.3+(i-1)*.5),['r_l=',num2str(round(Rerror{i}(1)*100)/100),' r_r=',num2str(round(Rerror{i}(2)*100)/100)],'Color',ColorsError{i});
    text(max(get(gca,'XLim'))+0.03,max(get(gca,'YLim'))-unit*(0.4+(i-1)*.5),['r=',num2str(round(Rerror{i}(3)*100)/100),', p=',num2str(round(Perror{i}(3)*100)/100)],'Color',ColorsError{i});
    
end



%reaction time
panel=subplot(2,4,5);
hold on
if sum(CompletedTrials)>1
    center = linspace(min(ST(CompletedTrials)),max(ST(CompletedTrials)),15);
    h=hist(ST(CompletedTrials),center);
    if ~isempty(h)
        h=h/sum(h);
        % ylabel('p')
        plot(abs(ExperiencedDV(CompletedTrials)),ST(CompletedTrials),'.k');
        xlabel('DV');ylabel('Sampling time (s)')
        ax2 = axes('Position',panel.Position);panel.Position=ax2.Position;
        plot(h,center,'r','LineWidth',2,'Parent',ax2);
        ax2.YAxis.Visible='off';ax2.XAxisLocation='top';ax2.Color='none';ax2.XAxis.FontSize = 8;ax2.XAxis.Color=[1,0,0];ax2.XLabel.String = 'p';ax2.XLabel.Position=[0.15,3.1,0];
        [r,p]=corr(abs(ExperiencedDV(CompletedTrials&~isnan(ST)))',ST(CompletedTrials&~isnan(ST))','type','Spearman');
        text(min(get(gca,'XLim'))+0.05,max(get(gca,'YLim'))-0.1,['r=',num2str(round(r*100)/100),', p=',num2str(round(p*100)/100)]);
        
    end
end

%grace periods
subplot(2,4,6)
%remove "full" grace periods
GracePeriods(GracePeriods>=GracePeriodsMax-0.001 & GracePeriods<=GracePeriodsMax+0.001 )=[];
GracePeriodsR(GracePeriodsR>=GracePeriodsMax-0.001 & GracePeriodsR<=GracePeriodsMax+0.001 )=[];
GracePeriodsL(GracePeriodsL>=GracePeriodsMax-0.001 & GracePeriodsL<=GracePeriodsMax+0.001 )=[];
center = 0:0.025:max(GracePeriods);
if ~all(isnan(GracePeriodsL)) && numel(center) > 1 && ~all(isnan(GracePeriodsR))
    g = hist(GracePeriods,center);g=g/sum(g);
    gl = hist(GracePeriodsL,center);gl=gl/sum(gl);
    gr = hist(GracePeriodsR,center);gr=gr/sum(gr);
    hold on
    plot(center,g,'k','LineWidth',2)
    plot(center,gl,'m','LineWidth',1)
    plot(center,gr,'c','LineWidth',1)
    xlabel('Grace period (s)');ylabel('p');
    text(min(get(gca,'XLim'))+0.05,max(get(gca,'YLim'))-0.05,['n=',num2str(sum(~isnan(GracePeriods))),'(',num2str(sum(~isnan(GracePeriodsL))),'/',num2str(sum(~isnan(GracePeriodsR))),')']);
end

%waiting time distributions
ColorsCond = {[.5,.5,.5],[.9,.1,.1]};
if length(LaserCond)==1
    %no laser
    subplot(2,4,7)
    hold on
    xlabel('waiting time (s)'); ylabel ('n trials');
    WTnoFeedbackL = WT(~Feedback & ChoiceLeft == 1);
    WTnoFeedbackR = WT(~Feedback & ChoiceLeft == 0);
    histogram(WTnoFeedbackL,10,'EdgeColor','none','FaceColor',[.2,.2,1]);
    histogram(WTnoFeedbackR,10,'EdgeColor','none','FaceColor',[.8,.6,.1]);

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
    subplot(2,4,7)
    hold on
    xlabel('waiting time (s)'); ylabel ('n trials');
    subplot(2,4,8)
    hold on
    xlabel('waiting time (s)'); ylabel ('n trials');
    PshortWTL=cell(1,2);PshortWTR=cell(1,2);
    for i =1:length(LaserCond)
    
    WTnoFeedbackL = WT(~Feedback & ChoiceLeft == 1 & LaserTrial==LaserCond(i));
    WTnoFeedbackR = WT(~Feedback & ChoiceLeft == 0 & LaserTrial==LaserCond(i));
     meanWTL = nanmean(WTnoFeedbackL);
    meanWTR = nanmean(WTnoFeedbackR);
    subplot(2,4,7)
    histogram(WTnoFeedbackL,10,'EdgeColor','none','FaceColor',ColorsCond{i});
    line([meanWTL,meanWTL],get(gca,'YLim'),'Color',ColorsCond{i});
    text(meanWTL-1,(1.05-0.1*(i-1))*(max(get(gca,'YLim'))-min(get(gca,'YLim'))),['m_l=',num2str(round(meanWTL*10)/10)],'Color',ColorsCond{i});
    subplot(2,4,8)
    histogram(WTnoFeedbackR,10,'EdgeColor','none','FaceColor',ColorsCond{i});
    line([meanWTR,meanWTR],get(gca,'YLim'),'Color',ColorsCond{i});
    text(meanWTL-1,(1.05-0.1*(i-1))*(max(get(gca,'YLim'))-min(get(gca,'YLim'))),['m_r=',num2str(round(meanWTR*10)/10)],'Color',ColorsCond{i});

    PshortWTL{i} = sum(WTnoFeedbackL<MinWT)/sum(~isnan(WTnoFeedbackL));
    PshortWTR{i} = sum(WTnoFeedbackR<MinWT)/sum(~isnan(WTnoFeedbackR));
    
    end
    for i =1:length(LaserCond)
        text(max(get(gca,'XLim'))+0.03,(0.85/i)*(max(get(gca,'YLim'))-min(get(gca,'YLim')))+min(get(gca,'YLim')),['L_{2}=',num2str(round(PshortWTL{i}*100)/100),', R_{2}=',num2str(round(PshortWTR{i}*100)/100)],'Color',ColorsCond{i});
    end
end

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