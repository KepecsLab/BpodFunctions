%Behavior analysis Olf2AFC
function UserKillScriptTorben

global BpodSystem

MailAddress = 'olrin.f6cfd9@m.evernote.com';

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

sent = SendMyMail(MailAddress,Subject,Body,Attachment);

if sent
    fprintf('Figure "%s" sent to %s.\n',FigureString,MailAddress);
else
    fprintf('Error:SendFigureTo:Mail could not be sent to %s.\n',MailAddress);
end

end

function FigHandle = Analysis()

global TaskParameters
global BpodSystem

GracePeriodsMax = 0.25;%only corrects for this grace period

MinWT = TaskParameters.GUI.VevaiometricMinWT;
MaxWT = 10;
AudBinWT = 6;
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
LeftClickTrain = BpodSystem.Data.Custom.LeftClickTrain(1:nTrials-1);
RightClickTrain = BpodSystem.Data.Custom.RightClickTrain(1:nTrials-1);

%define "completed trial"
% not that abvious for errors
%Correct vector is 1 for correct choice, 0 for incorrect choice, nan
%for no choice
%now: correct --> received feedback (reward)
%     error --> always (if error catch)
%     catch --> always (if choice happened)

CompletedTrials = (Feedback&Correct==1) | (Correct==0) | CatchTrial&~isnan(ChoiceLeft);

%calculate exerienced dv
ExperiencedDV=zeros(1,length(ST));
for t = 1 : length(ST)
    R = BpodSystem.Data.Custom.RightClickTrain{t};
    L = BpodSystem.Data.Custom.LeftClickTrain{t};
    Ri = find(R>ST(t),1,'first'); if isempty(Ri), Ri=1; end
    Li = find(L>ST(t),1,'first');if isempty(Li), Li=1; end
    ExperiencedDV(t) = log(Li/Ri);
    %         ExperiencedDV(t) = (Li-Ri)./(Li+Ri);
end

%caclulate  grace periods
GracePeriods=[];
for t = 1 : length(ST)
    GracePeriods = [GracePeriods;BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Rin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.rewarded_Rin_grace(:,1)];
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
AudBin = 8;
if ~isempty(AudDV)
    BinIdx = discretize(AudDV,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBin+1));
    PsycY = grpstats(ChoiceLeft(CompletedTrials),BinIdx,'mean');
    PsycX = grpstats(ExperiencedDV(CompletedTrials),BinIdx,'mean');
    plot(PsycX,PsycY,'ok','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',6)
    XFit = linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100);
    YFit = glmval(glmfit(AudDV,ChoiceLeft(CompletedTrials)','binomial'),linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100),'logit');
    plot(XFit,YFit,'k');
    xlabel('DV');ylabel('%left')
    text(0.95*min(get(gca,'XLim')),0.96*max(get(gca,'YLim')),[num2str(round(nanmean(Correct(CompletedTrials))*100)/100),'%']);
end

%conditioned psychometric
subplot(2,4,2)
hold on
%low
WTmed=median(WT(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT));
AudDV = ExperiencedDV(CompletedTrials&CatchTrial&WT<=WTmed&WT>MinWT);
if ~isempty(AudDV)
    ChoiceLeft = ChoiceLeft(CompletedTrials&CatchTrial&WT<=WTmed&WT>MinWT);
    AudBin = 8;
    BinIdx = discretize(AudDV,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBin+1));
    PsycY = grpstats(ChoiceLeft,BinIdx,'mean');
    PsycX = grpstats(AudDV,BinIdx,'mean');
    h1=plot(PsycX,PsycY,'ok','MarkerFaceColor',[.5,.5,.5],'MarkerEdgeColor','w','MarkerSize',6);
    XFit = linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100);
    YFit = glmval(glmfit(AudDV,ChoiceLeft','binomial'),linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100),'logit');
    plot(XFit,YFit,'Color',[.5,.5,.5]);
    %high
    AudDV = ExperiencedDV(CompletedTrials&CatchTrial&WT>WTmed&WT<MaxWT);
    ChoiceLeft = ChoiceLeft(CompletedTrials&CatchTrial&WT>WTmed&WT<MaxWT);
    AudBin = 8;
    BinIdx = discretize(AudDV,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBin+1));
    PsycY = grpstats(ChoiceLeft,BinIdx,'mean');
    PsycX = grpstats(AudDV,BinIdx,'mean');
    h2=plot(PsycX,PsycY,'ok','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',6);
    XFit = linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100);
    YFit = glmval(glmfit(AudDV,ChoiceLeft','binomial'),linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100),'logit');
    plot(XFit,YFit,'k');
    xlabel('DV');ylabel('%left')
    legend([h2,h1],{['WT>',num2str(round(WTmed*100)/100)],['WT<',num2str(round(WTmed*100)/100)]},'Units','normalized','Position',[0.333,0.85,0.1,0.1])
end

%calibration
subplot(2,4,3)
WTBin=5;
WTCatch = WT(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT);
if ~isempty(WTCatch)
    BinIdx = discretize(WTCatch,linspace(min(WTCatch)-10*eps,max(WTCatch)+10*eps,WTBin+1));
    WTX = grpstats(WTCatch,BinIdx,'mean');
    PerfY = grpstats(Correct(CompletedTrials&CatchTrial&WT>MinWT&WT<MaxWT),BinIdx,'mean');
    plot(WTX,PerfY,'k','LineWidth',2);
    xlabel('Waiting time (s)');ylabel('% correct')
end

%Vevaiometric
subplot(2,4,4)
hold on

WTCatch = WT(CompletedTrials&CatchTrial&Correct==1&WT>MinWT&WT<MaxWT);
DVCatch = ExperiencedDV(CompletedTrials&CatchTrial&Correct==1&WT>MinWT&WT<MaxWT);
if ~isempty(DVCatch)
    BinIdx = discretize(DVCatch,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBinWT+1));
    WTCatchY = grpstats(WTCatch,BinIdx,'mean');
    DVCatchX = grpstats(DVCatch,BinIdx,'mean');
    plot(DVCatchX,WTCatchY,'g','LineWidth',2)
    WTError = WT(CompletedTrials&Correct==0&WT>MinWT&WT<MaxWT);
    DVError = ExperiencedDV(CompletedTrials&Correct==0&WT>MinWT&WT<MaxWT);
    BinIdx = discretize(DVError,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBinWT+1));
    WTErrorY = grpstats(WTError,BinIdx,'mean');
    DVErrorX = grpstats(DVError,BinIdx,'mean');
    plot(DVErrorX,WTErrorY,'r','LineWidth',2)
    xlabel('DV');ylabel('Waiting time (s)')
    plot(DVCatch,WTCatch,'og','MarkerSize',2,'MarkerFaceColor','g')
    plot(DVError,WTError,'or','MarkerSize',2,'MarkerFaceColor','r')
    legend('Correct Catch','Error','Location','best')
end

%reaction time
subplot(2,4,5)
hold on
if sum(CompletedTrials)>1
    center = linspace(min(ST(CompletedTrials)),max(ST(CompletedTrials)),15);
    h=hist(ST(CompletedTrials),center);
    if ~isempty(h)
        h=h/sum(h);
        plot(center,h,'k','LineWidth',2)
        xlabel('Sampling time (s)')
        ylabel('p')
        subplot(2,4,6)
        plot(abs(ExperiencedDV(CompletedTrials)),ST(CompletedTrials),'.k')
        [r,p]=corrcoef(abs(ExperiencedDV(CompletedTrials&~isnan(ST))),ST(CompletedTrials&~isnan(ST)));
        xlabel('DV');ylabel('Sampling time (s)')
        text(min(get(gca,'XLim'))+0.05,0.99*max(get(gca,'YLim')),['r=',num2str(round(r(1,2)*100)/100),', p=',num2str(round(p(1,2)*100)/100)]);
    end
end

%grace periods
subplot(2,4,7)
%remove "full" grace periods
GracePeriods(GracePeriods>=GracePeriodsMax-0.001 & GracePeriods<=GracePeriodsMax+0.001 )=[];
center = 0:0.01:max(GracePeriods);
if ~all(isnan(GracePeriods))
    g = hist(GracePeriods,center);
    g=g/sum(g);
    plot(center,g,'k','LineWidth',2)
    xlabel('Grace period (s)');ylabel('p');
end

%cta
subplot(2,4,8)
hold on
Index50Fifty = abs(DV)< .1;
if sum(Index50Fifty)>1
    LData = Times2Table(LeftClickTrain(CompletedTrials&Index50Fifty),3);
    RData = Times2Table(RightClickTrain(CompletedTrials&Index50Fifty),3);
    ChoiceLeftadj = ChoiceLeft(CompletedTrials&Index50Fifty);
    window = 150;
    [cta_chosenL,~]=CTA(LData(ChoiceLeftadj==1,:),ceil(ST(CompletedTrials&Index50Fifty&ChoiceLeft==1)*1000),window);
    [cta_chosenR,~]=CTA(RData(ChoiceLeftadj==0,:),ceil(ST(CompletedTrials&Index50Fifty&ChoiceLeft==0)*1000),window);
    [cta_alternativeL,~]=CTA(LData(ChoiceLeftadj==0,:),ceil(ST(CompletedTrials&Index50Fifty&ChoiceLeft==0)*1000),window);
    [cta_alternativeR,cta_time]=CTA(RData(ChoiceLeftadj==1,:),ceil(ST(CompletedTrials&Index50Fifty&ChoiceLeft==1)*1000),window);
    
    cta_chosen = [cta_chosenL;cta_chosenR];
    cta_alternative = [cta_alternativeL;cta_alternativeR];
    
    cta_chosen_mean = nanmean(cta_chosen,1);
    cta_alternative_mean = nanmean(cta_alternative,1);
    
    nanplot(cta_time(cta_time<0)./1000,cta_chosen_mean(cta_time<0),'Color',[0,0,1],'LineWidth',2)
    nanplot(cta_time(cta_time<0)./1000,cta_alternative_mean(cta_time<0),'Color',[1,0,0],'LineWidth',2)
    xlabel('T-choice (s)')
    ylabel('Excess clicks/s')
    lowcut = quantile(ST(CompletedTrials&Index50Fifty),0.95);
    xlim([-lowcut,0])
    % ylim([-15,15])
    legend('Chosen','Alternative','Location','best')
    legend('boxoff')
    plot(get(gca,'XLim'),[.5,.5],'k--')
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

%Compute CTA
%2nd dim of E has to be "time"
function [cta2,t] = CTA(E,idx,win)


T = size(E,2);

cta = nan(size(E,1),2*T);

for i = 1:size(E,1)
    
    cta(i,T+1:T+T-idx(i)+1)=E(i,idx(i):end);
    cta(i,T-idx(i)+1:T) = E(i,1:idx(i));
    
end

%subtract mean clicks and convert to per sec
cta = (cta-repmat(mean(E,2),1,2*T))*1000;

%smooth/window
Wins = 1:win:2*T;
cta2=nan(size(E,1),length(Wins));
for i = 1:length(Wins)-1
    cta2(:,i) = mean(cta(:,Wins(i):Wins(i+1)),2);
end
t = -T:win:T-1;
end

