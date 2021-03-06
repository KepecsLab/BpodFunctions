%Behavior analysis Olf2AFC
function UserKillScriptMatching

global BpodSystem
global TaskParameters

%
% for photometry sessions, put back single-trial data files into
% Bpod struct for convenience
if TaskParameters.GUI.Photometry
    try
    for i =1:BpodSystem.Data.nTrials
        fname = fullfile(BpodSystem.DataPath(1:end-4),['NidaqData',num2str(i),'.mat']);
        tmp=load(fname);
        BpodSystem.Data.NidaqData{i}=tmp.NidaqData;
        BpodSystem.Data.Nidaq2Data{i}=tmp.Nidaq2Data;
    end
    
    %save
    SessionData = BpodSystem.Data;
    fname = [BpodSystem.DataPath(1:end-4),'_stitch.mat'];
    save(fname, 'SessionData', '-v7.3');
    
    
    %double check if file exists and has reasonable size (there were issues
    %with empty files)
    f = dir(fname);
    if ~isempty(f) && f.size/1000/1000 > 100 %larger than 100 MB
        %-->delete pre-stitch files & folder
        for i =1:BpodSystem.Data.nTrials
            fname = fullfile(BpodSystem.DataPath(1:end-4),['NidaqData',num2str(i),'.mat']);
            if exist(fname,'file')==2
                delete(fname) %only delete AFTER saving bpod file in case sth goes wrong
            end
        end
        rmdir(BpodSystem.DataPath(1:end-4))
    end
    catch
        warning('Stitching single-trial photometry data back to SessionData failed.')
    end
end

%load mail settings --> contains mail address & password & evernote e-mail address
load('MailSettings.mat');%loads MailSettings struct

%evernote mail address
MailAddress = MailSettings.EvernoteMail;

%save figure
FigureFolder = fullfile(fileparts(fileparts(BpodSystem.DataPath)),'Session Figures');
FigureHandle = BpodSystem.GUIHandles.Axes.OutcomePlot.MainHandle.Parent; %BpodSystem.GUIHandles.OutcomePlot.HandleOutcome.Parent;
FigureString = get(FigureHandle,'Name');
[~, FigureName] = fileparts(BpodSystem.DataPath);
if ~isdir(FigureFolder)
    mkdir(FigureFolder);
end

FigurePath = fullfile(FigureFolder,[FigureName,'.png']);
saveas(FigureHandle,FigurePath,'png');
try
    close(FigureHandle)
catch
end

%Analysis
try
    FigAnalysis = Analysis();
    FigurePathAnalysis = fullfile(FigureFolder,[FigureName,'Analysis.png']);
    saveas(FigAnalysis,FigurePathAnalysis,'png');
    close(FigAnalysis);
    DidAnalysis = true;
catch
    DidAnalysis = false;
    fprintf('Failed to run analysis function\n');
    try
        close(FigAnalysis)
    catch
    end
end

%check whether there are photometry figures to add to mail
try
    FigurePathPhotometry1 = fullfile(FigureFolder,[FigureName,'Photometry1.png']);
    saveas(BpodSystem.GUIHandles.PhotometryPlot1,FigurePathPhotometry1,'png');
    close(BpodSystem.GUIHandles.PhotometryPlot1)
    AddPhotometry = true;
    PhotometryPath = {FigurePathPhotometry1};
    try
        FigurePathPhotometry2 = fullfile(FigureFolder,[FigureName,'Photometry2.png']);
        saveas(BpodSystem.GUIHandles.PhotometryPlot2,FigurePathPhotometry2,'png');
        close(BpodSystem.GUIHandles.PhotometryPlot2)
        PhotometryPath = {FigurePathPhotometry1, FigurePathPhotometry2};
    catch
    end
catch
    AddPhotometry = false;
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
if AddPhotometry
    Attachment = [Attachment, PhotometryPath];
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
        [~,result]=system('ipconfig');
        if contains(result,'172') || contains(result,'wustl') %%-->WUSTL address space
            servername = '\\172.20.22.201\home';
        else %--> CSHL address space
            servername = '\\uncertainty.cshl.edu\home'; %new uncertanity server (8/2018) works with home only
        end
%         user = strcat(getenv('username'));
        user ='';
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

offline=false;
if offline
    BpodSystem.Data=SessionData;
    if isfield(SessionData,'TrialSettings')
        TaskParameters = SessionData.TrialSettings(1);
    else
        TaskParameters = SessionData.Settings;
    end
    Animal ='Unknown';
else
    [~,Animal]=fileparts(fileparts(fileparts(fileparts(BpodSystem.DataPath))));
end

% correct feedback delay. weirdly (wrongly?!) done by TG in
%% updateCustomDataFields of Matching task. TO 10/2020
nTrials=BpodSystem.Data.nTrials;
TimeChoice = nan(1,nTrials);
FeedbackDelay = nan(1,nTrials);
PokeOut = nan(1,nTrials);
PokeIn = nan(1,nTrials);
for n =1:nTrials
    statetimes = BpodSystem.Data.RawEvents.Trial{n}.States;
    if  BpodSystem.Data.Custom.ChoiceLeft(n)==1
        choicename = 'start_Lin';
        rewardedname = 'water_L';
        unrewardedname = 'unrewarded_Lin';
        earlyname = 'EarlyLout';
    elseif BpodSystem.Data.Custom.ChoiceLeft(n) == 0
        choicename = 'start_Rin';
        rewardedname = 'water_R';
        unrewardedname = 'unrewarded_Rin';
        earlyname = 'EarlyRout';
    else
        choicename = 'start_Lin'; %will be NaN
    end
    
    if BpodSystem.Data.Custom.Rewarded(n)==1
        outtime = statetimes.(rewardedname)(1,1);
    elseif BpodSystem.Data.Custom.Rewarded(n)==0 && ~isnan( BpodSystem.Data.Custom.ChoiceLeft(n))
        %for unrewarded trials, there are multiple possibilities... left as
        %an 'early trial' or trial as an 'unrewarded' trial
        if BpodSystem.Data.Custom.EarlySout(n) == 1
            outtime = statetimes.(earlyname)(1,1) - TaskParameters.GUI.Grace;
        else
            outtime = statetimes.(unrewardedname)(1,1) - TaskParameters.GUI.Grace;
        end
    elseif isnan( BpodSystem.Data.Custom.ChoiceLeft(n)) %no choice
        outtime = NaN;
    end
    
        %trial
    
    TimeChoice(n) = statetimes.(choicename)(1,1);
    FeedbackDelay(n) = outtime - TimeChoice(n);
    
    PokeIn(n) = statetimes.Cin(1);
    PokeOut(n)= statetimes.Cin(2);
    
end

%correct feedback delay field
BpodSystem.Data.Custom.FeedbackDelayCorrected = FeedbackDelay;

%%
FigHandle = figure('Position',[ 385    56   869   904],'NumberTitle','off','Name',Animal,'Color',[1,1,1]);
nTrials=BpodSystem.Data.nTrials;
ChoiceLeft = BpodSystem.Data.Custom.ChoiceLeft;
LeftHi=double(BpodSystem.Data.Custom.LeftHi);
LeftHi(LeftHi==1)=TaskParameters.GUI.pHi/100;
LeftHi(LeftHi==0)=TaskParameters.GUI.pLo/100;

%plot running choice average
subplot(4,3,[1 2 3])

if ~isempty(ChoiceLeft)
    Xdata = 1:nTrials-1;
    Ydata = LeftHi(1:nTrials-1);
    
    plot(Xdata,Ydata,'-k','Color',[.5,.5,.5],'LineWidth',2);
    hold on;
    
    smoothChoice = smooth(ChoiceLeft, 10, 'moving','omitnan');
    Ydata=smoothChoice(1:nTrials-1);
    plot(Xdata,Ydata,'-k','LineWidth',2);
    ylabel('P(Left)')
    xlabel('Trials')
    ylim([0,1])
    xlim([0,nTrials])
end
    
 %run lau/glimcher model
try
    [ Mdl, logodds ] = LauGlim( BpodSystem.Data );
    succ=true;
catch
    fprintf('Failed to run LauGlim model\n');
    succ=false;
end

if succ
    
%GLM Fit
subplot(4,3,4)
hold on

ChoiceKernelRwd_YData = Mdl.Coefficients.Estimate(7:11);
ChoiceKernelCho_YData = Mdl.Coefficients.Estimate(2:6);
intercept = Mdl.Coefficients.Estimate(1);

plot(1:length(ChoiceKernelRwd_YData), ChoiceKernelRwd_YData,'LineWidth',2,'Color',[.3,.2,.9])
plot(1:length(ChoiceKernelCho_YData), ChoiceKernelCho_YData,'LineWidth',2,'Color',[.8,.2,.8])
scatter(1,intercept,'filled','MarkerFaceColor','k')
plot([1,length(ChoiceKernelRwd_YData)],[intercept,intercept],'--k')
xlim([1,length(ChoiceKernelRwd_YData)])
ylabel('Coefficient')
xlabel('n-trials back')
l=legend('rewardKernel','choiceKernel','intercept');
l.Box='off';

%DV plot - psychometric
subplot(4,3,5)
hold on
ndxValid =~isnan(BpodSystem.Data.Custom.ChoiceLeft); ndxValid = ndxValid(:);
ChoiceLeft = BpodSystem.Data.Custom.ChoiceLeft(ndxValid);
DV = logodds(ndxValid);
dvbin=linspace(-max(abs(DV)),max(abs(DV)),10);
[x,y,e]=BinData(DV,ChoiceLeft,dvbin);
vv=~isnan(x) & ~isnan(y) & ~isnan(e);
errorbar(x(vv),y(vv),e(vv),'k','LineStyle','none','LineWidth',2,'Marker','o','MarkerFaceColor','k')

xlim([dvbin(1),dvbin(end)+eps]);
ylim([0,1])

xlabel('log odds')
ylabel('P(Left)')
%fit
mdl = fitglm(DV,ChoiceLeft(:),'Distribution','binomial');
xx=linspace(dvbin(1),dvbin(end),100);
plot(xx,predict(mdl,xx'),'-k')

%waiting time distribution plot
leave_session = TaskParameters.GUI.CatchUnrwd ==1;
if leave_session
subplot(4,3,6)
hold on

ndxBaited = (BpodSystem.Data.Custom.Baited.Left & BpodSystem.Data.Custom.ChoiceLeft==1) | (BpodSystem.Data.Custom.Baited.Right & BpodSystem.Data.Custom.ChoiceLeft==0);
ndxBaited = ndxBaited(:);
ndxValid = BpodSystem.Data.Custom.EarlyCout==0 & ~isnan(BpodSystem.Data.Custom.ChoiceLeft); ndxValid = ndxValid(:);

ti = BpodSystem.Data.Custom.FeedbackDelayCorrected(ndxValid & ~ndxBaited);
reward_delay = BpodSystem.Data.Custom.FeedbackDelayCorrected(ndxValid & BpodSystem.Data.Custom.Rewarded'==1);
cc = linspace(0,max(ti),16);
histogram(ti,cc,'Normalization','probability','FaceColor',[.6,.6,.6],'EdgeColor',[1,1,1])
hi = histcounts(reward_delay,cc,'Normalization','probability');
plot((cc(1:end-1) + cc(2:end))/2,hi,'Color',[.1,.1,1]);
xlabel('Time investment (s)')
ylabel('p')

%'calibration' plot
subplot(4,3,7)
hold on
ndxExploit = BpodSystem.Data.Custom.ChoiceLeft(:) == (logodds>0);

left = BpodSystem.Data.Custom.ChoiceLeft(ndxValid & ~ndxBaited)==1;
corr = ndxExploit(ndxValid & ~ndxBaited); %'correct'
edges = linspace(min(ti),max(ti),8);
[x,y,e]=BinData(ti,corr,edges);
vv=~isnan(x) & ~isnan(y) & ~isnan(e);
errorbar(x(vv),y(vv),e(vv),'Color','k','LineWidth',2)
xlabel('Time investment (s)')
ylabel('Percent exploit')

%plot vevaiometric    
subplot(4,3,8)
hold on

ndxBaited = (BpodSystem.Data.Custom.Baited.Left & BpodSystem.Data.Custom.ChoiceLeft==1) | (BpodSystem.Data.Custom.Baited.Right & BpodSystem.Data.Custom.ChoiceLeft==0);
ndxBaited = ndxBaited(:);
ndxValid = BpodSystem.Data.Custom.EarlyCout==0 & ~isnan(BpodSystem.Data.Custom.ChoiceLeft); ndxValid = ndxValid(:);
ndxExploit = BpodSystem.Data.Custom.ChoiceLeft(:) == (logodds>0);
ExploreScatter_XData = logodds(ndxValid & ~ndxBaited & ~ndxExploit);
ExploreScatter_YData = BpodSystem.Data.Custom.FeedbackDelayCorrected(ndxValid & ~ndxBaited & ~ndxExploit)';
ExploitScatter_XData = logodds(ndxValid & ~ndxBaited & ndxExploit);
ExploitScatter_YData = BpodSystem.Data.Custom.FeedbackDelayCorrected(ndxValid & ~ndxBaited & ndxExploit)';
[ExploreLine_XData, ExploreLine_YData] = binvevaio(ExploreScatter_XData,ExploreScatter_YData,10);
[ExploitLine_XData, ExploitLine_YData] = binvevaio(ExploitScatter_XData,ExploitScatter_YData,10);


scatter(ExploitScatter_XData, ExploitScatter_YData,'.g','MarkerFaceColor','g');
scatter(ExploreScatter_XData, ExploreScatter_YData,'.r','MarkerFaceColor','r');
h1=plot(ExploreLine_XData, ExploreLine_YData, 'r','LineWidth',3);
h2=plot(ExploitLine_XData, ExploitLine_YData,'g','LineWidth',3);
l=legend([h1,h2],{'Explore','Exploit'});
l.Box='off';
l.Location='northwest';
try
ylim([min([ExploitScatter_YData;ExploreScatter_YData]),max([ExploitScatter_YData;ExploreScatter_YData])])
xlim([-max(abs([ExploitScatter_XData;ExploreScatter_XData])),max(abs([ExploitScatter_XData;ExploreScatter_XData]))])
catch
end

xlabel('log odds')
ylabel('Time investment (s)')

%'condition psychometry'
%DV plot - psychometric
subplot(4,3,9)
Colors={[0,0,.9];[.5,.5,1]}; %high/low
hold on
ChoiceLeft = BpodSystem.Data.Custom.ChoiceLeft(ndxValid & ~ndxBaited);ChoiceLeft=ChoiceLeft(:);
DV = logodds(ndxValid & ~ndxBaited);DV=DV(:);
TI = BpodSystem.Data.Custom.FeedbackDelayCorrected(ndxValid & ~ndxBaited);
TImed=nanmedian(TI);
high = TI>TImed;
low = TI<=TImed;
dvbin=linspace(-max(abs(DV)),max(abs(DV)),10);
%high TI
[x,y,e]=BinData(DV(high),ChoiceLeft(high),dvbin);
vv=~isnan(x) & ~isnan(y) & ~isnan(e);
h1=errorbar(x(vv),y(vv),e(vv),'LineStyle','none','LineWidth',2,'Marker','o','MarkerFaceColor',Colors{1},'MarkerEdgeColor',Colors{1},'Color',Colors{1});
%fit
mdl = fitglm(DV(high),ChoiceLeft(high),'Distribution','binomial');
xx=linspace(dvbin(1),dvbin(end),100);
plot(xx,predict(mdl,xx'),'-','Color',Colors{1});

%low TI
[x,y,e]=BinData(DV(low),ChoiceLeft(low),dvbin);
vv=~isnan(x) & ~isnan(y) & ~isnan(e);
h2=errorbar(x(vv),y(vv),e(vv),'LineStyle','none','LineWidth',2,'Marker','o','MarkerFaceColor',Colors{2},'MarkerEdgeColor',Colors{2},'Color',Colors{2});
%fit
mdl = fitglm(DV(low),ChoiceLeft(low),'Distribution','binomial');
xx=linspace(dvbin(1),dvbin(end),100);
plot(xx,predict(mdl,xx'),'-','Color',Colors{2});

xlim([dvbin(1),dvbin(end)+eps]);
ylim([0,1])

xlabel('log odds')
ylabel('P(Left)')

l=legend([h1,h2],{'High TI','Low TI'});
l.Box='off';
l.Location='northwest';
end%if leave_session
end%succ

%% session diagnostic plots
subplot(4,3,10)
hold on
%caclulate  grace periods
GracePeriods=[];
GracePeriodsL=[];
GracePeriodsR=[];
for t = 1 : nTrials
    GracePeriods = [GracePeriods;BpodSystem.Data.RawEvents.Trial{t}.States.Rin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.Rin_grace(:,1);BpodSystem.Data.RawEvents.Trial{t}.States.Lin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.Lin_grace(:,1)];
    if BpodSystem.Data.Custom.ChoiceLeft(t) == 1
        GracePeriodsL = [GracePeriodsL;BpodSystem.Data.RawEvents.Trial{t}.States.Lin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.Lin_grace(:,1)];
    elseif BpodSystem.Data.Custom.ChoiceLeft(t)==0
        GracePeriodsR = [GracePeriodsR;BpodSystem.Data.RawEvents.Trial{t}.States.Rin_grace(:,2)-BpodSystem.Data.RawEvents.Trial{t}.States.Rin_grace(:,1)];
    end
end
%remove "full" grace periods
GracePeriodsMax=TaskParameters.GUI.Grace;
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

subplot(4,3,11)
hold on
%drinking time dist
DrinkingTime=[];
for t = 1 : nTrials
    if BpodSystem.Data.Custom.ChoiceLeft(t)==1
        DrinkingTime(end+1)=BpodSystem.Data.RawEvents.Trial{t}.States.DrinkingL(end,end) - BpodSystem.Data.RawEvents.Trial{t}.States.DrinkingL(1,1);
    elseif BpodSystem.Data.Custom.ChoiceLeft(t)==0
        DrinkingTime(end+1)=BpodSystem.Data.RawEvents.Trial{t}.States.DrinkingR(end,end) - BpodSystem.Data.RawEvents.Trial{t}.States.DrinkingR(1,1);
    end
end
center = 0:0.2:max(DrinkingTime);
if ~all(isnan(DrinkingTime)) && numel(center) > 1
    histogram(DrinkingTime,center,'FaceColor',[.5,.5,.5],'EdgeColor',[1,1,1]);
     xlabel('Drinking times (s)');ylabel('n');
end

subplot(4,3,12)
hold on
%actual ITI lenghts
ITI = nan(nTrials-1,1);
for t = 1 : nTrials -1
    ITI(t) = BpodSystem.Data.TrialStartTimestamp(t+1) - BpodSystem.Data.TrialStartTimestamp(t) + BpodSystem.Data.RawEvents.Trial{t+1}.States.PreITI(1,2) - BpodSystem.Data.RawEvents.Trial{t}.States.ITI(1,1);
end
cc=linspace(min(ITI),max(ITI),20);
histogram(ITI,cc,'FaceColor',[.5,.5,.5],'EdgeColor',[1,1,1])
goal = TaskParameters.GUI.ITI + TaskParameters.GUI.PreITI;
line([goal,goal],get(gca,'YLim'),'Color',[1,0,0])
xlabel('Actual ITI (s)'); ylabel('n')



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

function [ mdl, logodds ] = LauGlim( SessionData )
%LAUGLIM Statistical model to predict single trial choice behavior as in
%Lau and Glimcher's JEAB paper (2005?)

y = SessionData.Custom.ChoiceLeft(:);
C = y;
C(y==0) = -1;
R = SessionData.Custom.Rewarded(:).*C; % hopefully vectors of same length at all times
C = repmat(C,1,5);
R = repmat(R,1,5);

for j = 1:size(C,2)
    C(:,j) = circshift(C(:,j),j);
    C(1:j,j) = 0;
    R(:,j) = circshift(R(:,j),j);
    R(1:j,j) = 0;
end

X = [C, R];
X(isnan(X)) = 0;
mdl = fitglm(X,y,'distribution','binomial');
logodds = mdl.predict(X);
logodds = log(logodds) - log(1-logodds);
end

function [ newxdata, newydata ] = binvevaio( xdata, ydata, nbins )
%  UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xdata = xdata(:);
ydata = ydata(:);

if nargin < 3
    nbins = ceil(numel(xdata)/10);
end

newxdata = nan(nbins,1);
newydata = nan(nbins,1);
ndx = nan(numel(xdata),1);

for ibin = 1:nbins
%     newxdata = prctile(xdata,100*ibin/nbins);
    ndx(isnan(ndx) & xdata <= prctile(xdata,100*ibin/nbins)) = ibin;
    newxdata(ibin) = mean(xdata(ndx==ibin));
    newydata(ibin) = mean(ydata(ndx==ibin));
end
end



