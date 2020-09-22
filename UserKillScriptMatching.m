%Behavior analysis Olf2AFC
function UserKillScriptMatching

global BpodSystem

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
    TaskParameters = SessionData.Settings;
    Animal ='Unknown';
else
    [~,Animal]=fileparts(fileparts(fileparts(fileparts(BpodSystem.DataPath))));
end
    


FigHandle = figure('Position',[ 360         187        1056         598],'NumberTitle','off','Name',Animal,'Color',[1,1,1]);
nTrials=BpodSystem.Data.nTrials;
ChoiceLeft = BpodSystem.Data.Custom.ChoiceLeft;
LeftHi=double(BpodSystem.Data.Custom.LeftHi);
LeftHi(LeftHi==1)=TaskParameters.GUI.pHi/100;
LeftHi(LeftHi==0)=TaskParameters.GUI.pLo/100;

%plot running choice average
subplot(2,3,[1 2 3])

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
    
    
%plot vevaiometric    
subplot(2,3,4)
hold on

[ Mdl, logodds ] = LauGlim( BpodSystem.Data );

ndxBaited = (BpodSystem.Data.Custom.Baited.Left & BpodSystem.Data.Custom.ChoiceLeft==1) | (BpodSystem.Data.Custom.Baited.Right & BpodSystem.Data.Custom.ChoiceLeft==0);
ndxBaited = ndxBaited(:);
ndxValid = BpodSystem.Data.Custom.EarlyCout==0 & ~isnan(BpodSystem.Data.Custom.ChoiceLeft); ndxValid = ndxValid(:);
ndxExploit = BpodSystem.Data.Custom.ChoiceLeft(:) == (logodds>0);
ExploreScatter_XData = logodds(ndxValid & ~ndxBaited & ~ndxExploit);
ExploreScatter_YData = BpodSystem.Data.Custom.FeedbackDelay(ndxValid & ~ndxBaited & ~ndxExploit)';
ExploitScatter_XData = logodds(ndxValid & ~ndxBaited & ndxExploit);
ExploitScatter_YData = BpodSystem.Data.Custom.FeedbackDelay(ndxValid & ~ndxBaited & ndxExploit)';
[ExploreLine_XData, ExploreLine_YData] = binvevaio(ExploreScatter_XData,ExploreScatter_YData);
[ExploitLine_XData, ExploitLine_YData] = binvevaio(ExploitScatter_XData,ExploitScatter_YData);


scatter(ExploitScatter_XData, ExploitScatter_YData,'.g','MarkerFaceColor','g')
scatter(ExploreScatter_XData, ExploreScatter_YData,'.r','MarkerFaceColor','r')
plot(ExploreLine_XData, ExploreLine_YData, 'r','LineWidth',3)
plot(ExploitLine_XData, ExploitLine_YData,'g','LineWidth',3)
ylim([min([ExploitScatter_YData;ExploreScatter_YData]),max([ExploitScatter_YData;ExploreScatter_YData])])
xlim([-max(abs([ExploitScatter_XData;ExploreScatter_XData])),max(abs([ExploitScatter_XData;ExploreScatter_XData]))])
xlabel('log odds')
ylabel('Waiting time (s)')

%GLM Fit
subplot(2,3,5)
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
subplot(2,3,6)
hold on
ndxValid =~isnan(BpodSystem.Data.Custom.ChoiceLeft); ndxValid = ndxValid(:);
ChoiceLeft = BpodSystem.Data.Custom.ChoiceLeft(ndxValid);
DV = logodds(ndxValid);
dvbin=linspace(-max(abs(DV)),max(abs(DV)),10);
[x,y,e]=BinData(DV,ChoiceLeft,dvbin);
vv=~isnan(x) & ~isnan(y) & ~isnan(e);
errorbar(x(vv),y(vv),e(vv),'k','LineStyle','none','LineWidth',2,'Marker','o','MarkerFaceColor','k')
xlim([dvbin(1),dvbin(end)]);
ylim([0,1])
xlabel('log odds')
ylabel('P(Left)')
%fit
mdl = fitglm(DV,ChoiceLeft(:),'Distribution','binomial');
xx=linspace(dvbin(1),dvbin(end),100);
plot(xx,predict(mdl,xx'),'-k')


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
UNTITLED Summary of this function goes here
  Detailed explanation goes here

xdata = xdata(:);
ydata = ydata(:);

if nargin < 3
    nbins = ceil(numel(xdata)/10);
end

newxdata = nan(nbins,1);
newydata = nan(nbins,1);
ndx = nan(numel(xdata),1);

for ibin = 1:nbins
    newxdata = prctile(xdata,100*ibin/nbins);
    ndx(isnan(ndx) & xdata <= prctile(xdata,100*ibin/nbins)) = ibin;
    newxdata(ibin) = mean(xdata(ndx==ibin));
    newydata(ibin) = mean(ydata(ndx==ibin));
end
end



