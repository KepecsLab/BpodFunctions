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

function [FigHandle,blockTransRich,idxHighVar, idxLowVar] = Analysis(SessionData)

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

% trial types reward bias - control == 1, left high== 2, right high ==3 -->
% blockTable
blockTable=[];
waterTable=horzcat(BpodSystem.Data.TrialSettings(end).GUI.BlockTable.RewL,BpodSystem.Data.TrialSettings(end).GUI.BlockTable.RewR);
blockTable(waterTable(:,1)==waterTable(:,2))=1;
blockTable(waterTable(:,1)>waterTable(:,2))=2;
blockTable(waterTable(:,1)<waterTable(:,2))=3;
blockTable=blockTable';

% trial types variance - 1 == highVar, 2==medVar, 3==lowVar --> varTableIdx

if sum(BpodSystem.Data.Settings.GUI.BlockTable.NoiseL ~=0) || sum(BpodSystem.Data.Settings.GUI.BlockTable.NoiseR ~=0) 
    %var large = 1. var med=2, control=3
    varTable=horzcat(BpodSystem.Data.TrialSettings(end).GUI.BlockTable.NoiseL, BpodSystem.Data.TrialSettings(end).GUI.BlockTable.NoiseR);
    high=max(max(varTable));
    low=min(min(varTable(varTable>0)));
    
    varTableIdx=ones(length(blockTable),2);
    varTableIdx(varTableIdx==1)=2;
    varTableIdx(varTable==high)=1;
    varTableIdx(varTable==low)=3;
    varTableIdx(varTable==0)=0;
   
    
else
    varTableIdx=nan(length(blockTable),2);
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

    blockTransMat=NaN(91, length(blockIndex_omitNan)+1);

    for xi=1:length(blockIndex_omitNan)
        try
            blockTransMat(1:length(completedTrials(blockIndex_omitNan(xi)-40:blockIndex_omitNan(xi)+50)),xi+1)=completedTrials(blockIndex_omitNan(xi)-40:blockIndex_omitNan(xi)+50);
        catch
            if xi==length(blockIndex_omitNan)
                blockTransMat(1:length(completedTrials(blockIndex_omitNan(xi)-40:end)),xi+1)=completedTrials(blockIndex_omitNan(xi)-40:end);
            else
                disp('error in calculating transition matrix')
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
    
    CondColors={[0, 0.4470, 0.7410], [0.9290, 0.6940, 0.1250], [0.8500, 0.3250, 0.0980]};

    if ~isempty(AudDV)
            BinIdx = discretize(AudDV,linspace(min(AudDV)-10*eps,max(AudDV)+10*eps,AudBin+1));
            PsycY = grpstats(ChoiceLeft(CompletedTrials),BinIdx,'mean');
            PsycX = grpstats(ExperiencedDV(CompletedTrials),BinIdx,'mean');
            plot(PsycX,PsycY,'ok','MarkerSize',6, 'MarkerEdgeColor', CondColors{blockTable(i,1)})
            hold on;
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
    numBlocks=BpodSystem.Data.Custom.BlockNumberL(end);
%     if sum(BpodSystem.Data.Custom.BlockNumberL==BpodSystem.Data.Custom.BlockNumberL(end))<80
%         numBlocks=numBlocks-1;
%         disp('last block has insufficienttrials')
%     end

    blockTable=blockTable(1:numBlocks);
    varTableIdx=varTableIdx(1:numBlocks,:);
    
    %psychometric avg across blocks
    subplot(3,3,[3])
    title('avg block psychometric')
    hold on
    
    try
    %control block == 1
    [controlMeanX, controlMeanY]=avgPsyc(XBlock, YBlock, blockTable, 1);
    %leftLarge == 2
    [LMeanX, LMeanY]=avgPsyc(XBlock, YBlock, blockTable, 2);
    %rightLarge == 3
    [RMeanX, RMeanY]=avgPsyc(XBlock, YBlock, blockTable, 3);
    plot(controlMeanX, controlMeanY, 'color', CondColors{3},'linewidth',2)
    plot(LMeanX, LMeanY,'color',CondColors{1}, 'linewidth',2)
    plot(RMeanX, RMeanY,'color',CondColors{2}, 'linewidth',2)
    catch
        disp('error in calculating average psychometric')
    end

    %plot avg. choice transition to control block
     preBlock=vertcat(nan, blockTable(1:end-1));
     prepostBlock=horzcat(preBlock, blockTable);
    
     %plot reward/block conditions for this session
     subplot(3,3,[6])
     plot(1:length(BpodSystem.Data.Custom.RewardMagnitude(:,1)),BpodSystem.Data.Custom.RewardMagnitude(:,1))
     hold on;
     plot(1:length(BpodSystem.Data.Custom.RewardMagnitude(:,2)),BpodSystem.Data.Custom.RewardMagnitude(:,2))
     legend({'left','right'},'Location','northeast')
     
     
    %plot block transitions
     subplot(3,3,[9])
     %xline(50)
     hold on;
     
     idxHighVar=[find(varTableIdx(:,1)==1); find(varTableIdx(:,2)==1)];
     idxLowVar=[find(varTableIdx(:,1)==3); find(varTableIdx(:,2)==3)];
     %idx for transitions is idxLowVar+1 etc...
     completedTrialscopy=completedTrials;
     
     %flip right trials st pChoice --> pRichChoice
     %richChoice=BpodSystem.Data.Custom.RewardBase((~isnan(ChoiceLeft)),:);
     %richChoice= richChoice(:,1)>richChoice(:,2);
     %completedTrialscopy(richChoice==0)=~completedTrialscopy(richChoice==0);
     
     blockTransRich=NaN(91, length(blockIndex_omitNan)+1);

     for xi=1:length(blockIndex_omitNan)
            try
                blockTransRich(1:length(completedTrialscopy(blockIndex_omitNan(xi)-40:blockIndex_omitNan(xi)+50)),xi+1)=completedTrialscopy(blockIndex_omitNan(xi)-40:blockIndex_omitNan(xi)+50);
            catch
                if xi==length(blockIndex_omitNan)
                    blockTransRich(1:length(completedTrialscopy(blockIndex_omitNan(xi)-40:end)),xi+1)=completedTrialscopy(blockIndex_omitNan(xi)-40:end);
                else
                    disp('error in calculating transition matrix')
                end
            end
     end
        
     hold on; plot(mean(movmean(blockTransRich(:,idxHighVar+1),30,'omitnan'),2))
     hold on; plot(mean(movmean(blockTransRich(:,idxLowVar+1),30,'omitnan'),2))
     
     hold on; x=repmat(40,[10,1]); y=linspace(0,1,10); line(x,y);
     legend({'high','low',''},'Location','northeast')

          
%     title('control block transition')
%     subplot(3,3,[6])
%     [a,b] =size(blockTransMat);
%     x = linspace(-50,60,a);
%     plot(x,mean(movmean(blockTransMat(:,blockType==3),25),2,'omitnan'), 'linewidth',2);
%     xline(0)

    %plot avg. choice transition to high reward block
%     aa=histcounts(BpodSystem.Data.Custom.BlockNumberL);
%     
%     if sum(aa<60)>0
%         disp('not enough trials in last block - no transition plotted')
%     end
    
%     blockTransMat_dummy=horzcat(blockTransMat(:,1),blockTransMat); %transition matrix for postBlock
%     blockType_dummy=vertcat(0,blockTable(aa>100));
%     postBlock_dummy=vertcat(0,blockType_dummy(1:end-1));
%     
%     leftHi2rightHi= blockType_dummy==2 & postBlock_dummy==3; %if the prev. block was left high
%     subplot(3,3,[9])
%     title('high to low transition')
%     [a,b] =size(blockTransMat_dummy);
%     x = linspace(-50,60,a);
%     plot(x,mean(movmean(blockTransMat_dummy(:,leftHi2rightHi),25),2,'omitnan'),'linewidth',2, 'DisplayName','L2R');
%     xline(0)
%     hold on;
% 
%     control2rightHi=blockType_dummy==1 & postBlock_dummy==3; %control to right Hi
%     subplot(3,3,[6])
%     [a,b] =size(blockTransMat_dummy);
%     x = linspace(-50,60,a);
%     plot(x,mean(movmean(blockTransMat_dummy(:,control2rightHi),25),2,'omitnan'), 'linewidth',2, 'DisplayName','C2R');
%     xline(0)
%     hold on;
%     
%     
%     rightHi2leftHi= blockType_dummy==1 & postBlock_dummy==2; %if the prev. block was tight high
%     subplot(3,3,[9])
%     [a,b] =size(blockTransMat_dummy);
%     x = linspace(-50,60,a);
%     plot(x,mean(movmean(blockTransMat_dummy(:,rightHi2leftHi),25),2,'omitnan'), 'linewidth',2, 'DisplayName','R2L');
%     xline(0)
%     hold on;
% 
%     control2leftHi=blockType_dummy==1 & postBlock_dummy==3; %control to left Hi
%     subplot(3,3,[6])
%     title('control block transition')
%     [a,b] =size(blockTransMat_dummy);
%     x = linspace(-50,60,a);
%     plot(x,mean(movmean(blockTransMat_dummy(:,control2leftHi),25),2,'omitnan'), 'linewidth',2, 'DisplayName','C2L');
%     xline(0)
%     hold on;
%     
%         
%         
%       
%            
%     leftHi2control=blockType_dummy==3 & blockType_dummy(find(blockType_dummy==3)-1)==1; %left hi to control
%     subplot(3,3,[6])
%     [a,b] =size(blockTransMat_dummy);
%     x = linspace(-50,60,a);
%     plot(x,mean(movmean(blockTransMat_dummy(:,leftHi2control),25),2,'omitnan'), 'linewidth',2,'DisplayName','L2C');
%     xline(0)
%     hold on;
% 
%     rightHi2control=blockType_dummy==3 & blockType_dummy(find(blockType_dummy==3)-1)==2; %right hi to control
%     subplot(3,3,[6])
%     [a,b] =size(blockTransMat_dummy);
%     x = linspace(-50,60,a);
%     plot(x,mean(movmean(blockTransMat_dummy(:,rightHi2control),25),2,'omitnan'), 'linewidth',2,'DisplayName','R2C' );
%     xline(0)
%     hold on;
% 
% 
%     legend(subplot(3,3,[9]));
%     legend(subplot(3,3,[6]));

        
 
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



