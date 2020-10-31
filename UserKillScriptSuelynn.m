%Behavior analysis Olf2AFC
function UserKillScriptSuelynn

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
    fprintf('Failed to run analysis function\n');
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
            servername = '\\172.20.22.201\homes';
        else %--> CSHL address space
            servername = '\\uncertainty.cshl.edu\home'; %new uncertanity server (8/2018) works with home only
        end
%         user = strcat(getenv('username'));
        user ='suelynn';
    else
        %servername = '/media/';
        %user='torben';
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
    

%%
FigHandle = figure('Position',[  385         170        1056         762],'NumberTitle','off','Name',Animal,'Color',[1,1,1]);


nTrials=BpodSystem.Data.nTrials;
ChoiceLeft = BpodSystem.Data.Custom.ChoiceLeft(1:nTrials-1);
ST = BpodSystem.Data.Custom.ST(1:nTrials-1);   
ExperiencedDV=zeros(1,length(ST));
AudBin=8;

%define trial types
%define block transitions
logicBlock=abs(BpodSystem.Data.Custom.BlockNumberL(1:end-1)-BpodSystem.Data.Custom.BlockNumberL(2:end));
blockIndex=find(logicBlock==1); 

%for plotting - block index of completed trials (e.g. ignore nans)
a=BpodSystem.Data.Custom.BlockNumberL(~isnan(ChoiceLeft));
plotBlockIdx=a(1:end-1)-a(2:end);
plotBlockIdx=abs(plotBlockIdx);
blockIndex_omitNan=find(plotBlockIdx==1);


leftLarge=BpodSystem.Data.Custom.RewardMagnitude(:,1)>BpodSystem.Data.Custom.RewardMagnitude(:,2);
rightLarge=BpodSystem.Data.Custom.RewardMagnitude(:,1)<BpodSystem.Data.Custom.RewardMagnitude(:,2);
controlMag=BpodSystem.Data.Custom.RewardMagnitude(:,1)==BpodSystem.Data.Custom.RewardMagnitude(:,2);


%blockType %2 columns - 1 reward 2 variance
blockTable=zeros(length(unique(BpodSystem.Data.Custom.BlockNumberL)),2);
%left large = 1. right large=2, control=3
%queries block index end of each block and the last trial of the session
%for end block
blockTable(leftLarge([blockIndex, nTrials-1]),1)=1;
blockTable(rightLarge([blockIndex, nTrials-1]),1)=2;
blockTable(controlMag([blockIndex, nTrials-1]),1)=3;

%calculate variance
if sum(BpodSystem.Data.Settings.GUI.BlockTable.NoiseL ~=0) && sum(BpodSystem.Data.Settings.GUI.BlockTable.NoiseR ~=0) 
    for j=1:length(blockIndex)
        if j==1
            blockVar(1:blockIndex(j),:)=repmat(var(BpodSystem.Data.Custom.RewardMagnitude(1:blockIndex(j),:)),blockIndex(j),1);
        elseif j==length(blockIndex)
            blockVar(blockIndex(j)+1:nTrials-1,:)=repmat(var(BpodSystem.Data.Custom.RewardMagnitude(blockIndex(j)+1:nTrials-1,:)), (nTrials-1)-blockIndex(j),1);
        else
            blockVar(blockIndex(j-1)+1:blockIndex(j)+1,:)=repmat(var(BpodSystem.Data.Custom.RewardMagnitude(blockIndex(j-1)+1:blockIndex(j),:)), blockIndex(j)+1-blockIndex(j-1),1);
        end 
    end

    %var large = 1. var med=2, control=3
    xx=blockVar([blockIndex, nTrials-1],:);

    try
        blockTable(:,2)=2; %default is medium trial - replace if highest var or control var
        blockTable(sum(xx==max(max(xx)),2),2)=1;
        blockTable(xx(:,1)==xx(:,2),2)=3;
    catch
        disp('Failed to caclulate variance')   
    end
    
else
    blockTable(:,2)=nan;
    disp('no variance in reward blocks')
end


%% structure choice data into a matrix - s.t. each column is a block

if length(unique(BpodSystem.Data.Custom.BlockNumberL))>1
    ThereAreBlocks=true;
    completedTrials=ChoiceLeft(~isnan(ChoiceLeft));
    choiceLeftMat=NaN(250, length(blockIndex_omitNan)+1);
    for xi=1:length(blockIndex_omitNan)+1
        if xi==1 %first block
            choiceLeftMat(1:length(completedTrials(1:blockIndex_omitNan(xi))),xi)=completedTrials(1:blockIndex_omitNan(xi));
        elseif xi==length(blockIndex_omitNan)+1 %last block
            choiceLeftMat(1:length(completedTrials(blockIndex_omitNan(end)+1:end)),xi)=completedTrials(blockIndex_omitNan(end)+1:end);
        else
            choiceLeftMat(1:length(completedTrials(blockIndex_omitNan(xi-1)+1:blockIndex_omitNan(xi))),xi)=completedTrials(blockIndex_omitNan(xi-1)+1:blockIndex_omitNan(xi));
        end

    end

    blockTransMat=NaN(111, length(blockIndex_omitNan)+1);

    for xi=1:length(blockIndex_omitNan)
        try
            blockTransMat(1:length(completedTrials(blockIndex_omitNan(xi)-50:blockIndex_omitNan(xi)+60)),xi+1)=completedTrials(blockIndex_omitNan(xi)-50:blockIndex_omitNan(xi)+60);
        catch
            if blockIndex_omitNan(xi)+60 > length(completedTrials)
                disp(['not enough trials in block. only ' num2str(length(completedTrials)-blockIndex_omitNan(xi)) 'trials in block' num2str(xi+1) ])
                %disp(['not enough trials in block. only ' nu2str(length(completedTrials)-blockIndex_omitNan(xi) 'trials in block' num2str(xi+1) ])
            elseif blockIndex_omitNan(xi)-50 < 0
                disp('not enough trials in block. Negative trials.')
            end
        end
    end


   
else
    disp('There are no blocks')
    ThereAreBlocks=false;
end
%%            

for t = 1 : length(ST)
    R = BpodSystem.Data.Custom.RightClickTrain{t};
    L = BpodSystem.Data.Custom.LeftClickTrain{t};
    Ri = sum(R<=ST(t));if Ri==0, Ri=1; end
    Li = sum(L<=ST(t));if Li==0, Li=1; end
    %ExperiencedDV(t) = log10(Li/Ri);
    ExperiencedDV(t) = (Li-Ri)./(Li+Ri);
end

%plot running choice average
subplot(3,3,[1 2])
title('running choice avg')
plot(movmean(ChoiceLeft(~isnan(ChoiceLeft)),20, 'omitnan'),'linewidth',2);
hold on; 
plot(plotBlockIdx);


%psychometric
subplot(3,3,[4 8])
hold on


for i =1: length(unique(BpodSystem.Data.Custom.BlockNumberL))
    CompletedTrials = ~isnan(BpodSystem.Data.Custom.ChoiceLeft(1:end-1)) & BpodSystem.Data.Custom.BlockNumberL(1:length(BpodSystem.Data.Custom.ChoiceLeft(1:end-1)))==i;
    AudDV = ExperiencedDV(CompletedTrials);
    
    CondColors={[0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0, 0.4470, 0.7410]};

    if ~isempty(AudDV)
            BinIdx = discretize(AudDV,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBin+1));
            PsycY = grpstats(ChoiceLeft(CompletedTrials),BinIdx,'mean');
            PsycX = grpstats(ExperiencedDV(CompletedTrials),BinIdx,'mean');
            %plot(PsycX,PsycY,'ok','MarkerSize',6)
            XFit = linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100);
            YFit = glmval(glmfit(AudDV,ChoiceLeft(CompletedTrials)','binomial'),linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,100),'logit');
            if length(unique(BpodSystem.Data.Custom.BlockNumberL))>1
                plot(XFit,YFit, 'color',CondColors{blockTable(i,1)}, 'linewidth',2);
            else
               plot(XFit,YFit, 'linewidth',2); 
            end
            xlabel('DV');ylabel('p left')
            hold on;
            XBlock{i}=XFit;
            YBlock{i}=YFit;
            %text(0.95*min(get(gca,'XLim')),0.96*max(get(gca,'YLim')),[num2str(round(nanmean(Correct(CompletedTrialsCond))*100)),'%,n=',num2str(nTrialsCompleted)]);
    end
end


if ThereAreBlocks
    
    %psychometric avg across blocks
    subplot(3,3,[3])
    title('avg block psychometric')
    hold on
    
    try
    blockType=blockTable(:,1);
    %control block == 3
    [controlMeanX, controlMeanY]=avgPsyc(XBlock, YBlock, blockType, 3);
    %leftLarge == 1
    [LMeanX, LMeanY]=avgPsyc(XBlock, YBlock, blockType, 1);
    %rightLarge == 2
    [RMeanX, RMeanY]=avgPsyc(XBlock, YBlock, blockType, 2);
    plot(controlMeanX, controlMeanY, 'color', CondColors{3},'linewidth',2)
    plot(LMeanX, LMeanY,'color',CondColors{1}, 'linewidth',2)
    plot(RMeanX, RMeanY,'color',CondColors{2}, 'linewidth',2)
    catch
        disp('error in calculating average psychometric')
    end

    %plot avg. choice transition to control block
    title('control block transition')
    subplot(3,3,[6])
    [a,b] =size(blockTransMat);
    x = linspace(-50,60,a);
    plot(x,mean(movmean(blockTransMat(:,blockType==3),25),2,'omitnan'), 'linewidth',2);
    xline(0)

    %plot avg. choice transition to high reward block
    hold on;
    subplot(3,3,[9])
    title('high to low transition')
    [a,b] =size(blockTransMat);
    x = linspace(-50,60,a);
    plot(x,mean(movmean(blockTransMat(:,blockType==2),25),2,'omitnan'),'linewidth',2);
    xline(0)

    subplot(3,3,[9])
    [a,b] =size(blockTransMat);
    x = linspace(-50,60,a);
    plot(x,mean(movmean(blockTransMat(:,blockType==1),25),2,'omitnan'), 'linewidth',2);
    xline(0)
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
                varargin{5}{k} = fullfile(varargin{5}{k}(2:end));
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

%cleans and concat multiple SessionData files
%input: SessionData.Custom...

function [blockMeanX, blockMeanY]=avgPsyc(XBlock, YBlock, blockType, blockNum)
    %XXBlock=cell2mat(XBlock(blockType==blockNum));
    XXBlock=reshape(cell2mat(XBlock(blockType==blockNum)), 100,[]);
    blockMeanX=mean(XXBlock,2);
    
    %YYBlock=cell2mat(YBlock(blockType==blockNum));
    YYBlock=reshape(cell2mat(YBlock(blockType==blockNum)), 100,[]);
    blockMeanY=mean(YYBlock,2);
end



