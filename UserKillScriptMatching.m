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

[~,Animal]=fileparts(fileparts(fileparts(fileparts(BpodSystem.DataPath))));
FigHandle = figure('Position',[ 360         187        1056         598],'NumberTitle','off','Name',Animal);
subplot(2,3,[1 2 3])

nTrials=BpodSystem.Data.nTrials;
ChoiceLeft = BpodSystem.Data.Custom.ChoiceLeft;

LeftHi=double(BpodSystem.Data.Custom.LeftHi);
LeftHi(LeftHi==1)=TaskParameters.GUI.pHi/100;
LeftHi(LeftHi==0)=TaskParameters.GUI.pLo/100;
    
    if ~isempty(ChoiceLeft)
        Xdata = 1:nTrials-1;
        Ydata = LeftHi(1:nTrials-1);
        
        plot(Xdata,Ydata);
        hold on;
        
        smoothChoice = smooth(ChoiceLeft, 10, 'moving','omitnan');
        Ydata=smoothChoice(1:nTrials-1);
        plot(Xdata,Ydata);
        ylabel('pLeft')
        xlabel('trials')
    end
    
    
%vevaiometric    
subplot(2,3,4)

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

plot(ExploreLine_XData, ExploreLine_YData, 'r')
hold on;
scatter(ExploitScatter_XData, ExploitScatter_YData,'r')

plot(ExploitLine_XData, ExploitLine_YData,'g')
scatter(ExploreScatter_XData, ExploreScatter_YData,'g')



%GLM Fit
subplot(2,3,5)

ChoiceKernelRwd_YData = Mdl.Coefficients.Estimate(7:11);
ChoiceKernelCho_YData = Mdl.Coefficients.Estimate(2:6);
intercept = Mdl.Coefficients.Estimate(1);
plot(1:length(ChoiceKernelRwd_YData), ChoiceKernelRwd_YData)
hold on

plot(1:length(ChoiceKernelCho_YData), ChoiceKernelCho_YData)
scatter(1,intercept,'filled')
ylabel('coefficients')
xlabel('n-trials back')
legend('rewardKernal','choiceKernal','intercept')

%DV plot - psychometric
subplot(2,3,6)


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
%UNTITLED Summary of this function goes here
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



