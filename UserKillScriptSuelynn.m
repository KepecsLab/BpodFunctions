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

try
S=concatSession(BpodSystem.Data.Custom);
catch
end


nTrials=BpodSystem.Data.nTrials;
ChoiceLeft = BpodSystem.Data.Custom.ChoiceLeft(1:nTrials-1);
ST = BpodSystem.Data.Custom.ST(1:nTrials-1);   
ExperiencedDV=zeros(1,length(ST));
AudBin=8;

LeftClickTrain = BpodSystem.Data.Custom.LeftClickTrain(1:nTrials-1);
RightClickTrain = BpodSystem.Data.Custom.RightClickTrain(1:nTrials-1);
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
plot(movmean(S.ChoiceLeft(~isnan(S.ChoiceLeft)),20, 'omitnan'),'linewidth',2);
hold on; 
try
plot(S.plotBlockIdx);
catch
end

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
                plot(XFit,YFit, 'color',CondColors{S.blockType(i)}, 'linewidth',2);
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


if length(unique(BpodSystem.Data.Custom.BlockNumberL))>1
    %psychometric avg across blocks
    subplot(3,3,[3])
    title('avg block psychometric')
    hold on

    blockType=S.blockType(:,1);
    %control block == 3
    [controlMeanX, controlMeanY]=avgPsyc(XBlock, YBlock, blockType, 3);
    %leftLarge == 1
    [LMeanX, LMeanY]=avgPsyc(XBlock, YBlock, blockType, 1);
    %rightLarge == 2
    [RMeanX, RMeanY]=avgPsyc(XBlock, YBlock, blockType, 2);
    plot(controlMeanX, controlMeanY, 'color', CondColors{3},'linewidth',2)
    plot(LMeanX, LMeanY,'color',CondColors{1}, 'linewidth',2)
    plot(RMeanX, RMeanY,'color',CondColors{2}, 'linewidth',2)

    %plot avg. choice transition to control block
    title('control block transition')
    subplot(3,3,[6])
    [a,b] =size(S.blockTransMat);
    x = linspace(-50,60,a);
    plot(x,mean(movmean(S.blockTransMat(:,blockType==3),25),2,'omitnan'), 'linewidth',2);
    xline(0)

    %plot avg. choice transition to high reward block
    hold on;
    subplot(3,3,[9])
    title('high to low transition')
    [a,b] =size(S.blockTransMat);
    x = linspace(-50,60,a);
    plot(x,mean(movmean(S.blockTransMat(:,blockType==2),25),2,'omitnan'),'linewidth',2);
    xline(0)

    subplot(3,3,[9])
    [a,b] =size(S.blockTransMat);
    x = linspace(-50,60,a);
    plot(x,mean(movmean(S.blockTransMat(:,blockType==1),25),2,'omitnan'), 'linewidth',2);
    xline(0)
else
    disp('there are no blocks detected')
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

function S = concatSession(varargin)
    fields = fieldnames(varargin{1});
    
    %before concat. cut data length to completedTrials only
    for i = 1:length(varargin)
        completedTrials=length(varargin{i}.ChoiceLeft);
        varargin{i}.SessionNum=repmat(i,1,completedTrials);
        for ii = 1:length(fields)
            aField=fields{ii};
            if length(varargin{i}.(aField))>completedTrials
                try
                varargin{i}.(aField)=varargin{i}.(aField)(1:completedTrials, :);
                catch
                varargin{i}.(aField)=varargin{i}.(aField)(1:completedTrials);
                end
            end
        end
     


        %add blockchange idx as custom data field

        varargin{i}.blockChange=abs(varargin{i}.BlockNumberL(1:end-1)-varargin{i}.BlockNumberL(2:end));
        varargin{i}.transType=zeros(1,completedTrials); 
 
        blockIdx=find(varargin{i}.blockChange==1);
        if ~isempty(blockIdx)
            for iii=1:length(blockIdx)
                thisIdx=blockIdx(iii);
                varargin{i}.transType(thisIdx)=iii;

            end

            %added to ignore nan trials, different index
            varargin{i}.plotChoiceLeft=varargin{i}.ChoiceLeft(~isnan(varargin{i}.ChoiceLeft));
            a=varargin{i}.BlockNumberL(~isnan(varargin{i}.ChoiceLeft));
            plotBlockIdx=a(1:end-1)-a(2:end);
            varargin{i}.plotBlockIdx=abs(plotBlockIdx);
            plotIdx=find(varargin{i}.plotBlockIdx==1);
            for iii=1:length(plotIdx)
                thisIdx=plotIdx(iii);
                varargin{i}.plotTransType(thisIdx)=iii;

            end


            %add trialTypes for analysis
            %rewardMag
            varargin{i}.leftLarge=varargin{i}.RewardMagnitude(:,1)>varargin{i}.RewardMagnitude(:,2);
            varargin{i}.rightLarge=varargin{i}.RewardMagnitude(:,1)<varargin{i}.RewardMagnitude(:,2);
            varargin{i}.controlMag=varargin{i}.RewardMagnitude(:,1)==varargin{i}.RewardMagnitude(:,2);

            %index for later variance trial types
            blockIndex=find(varargin{i}.blockChange==1); 

            for j=1:length(blockIndex)
                if j==1
                    varargin{i}.variance(1:blockIndex(j),:)=repmat(var(varargin{i}.RewardMagnitude(1:blockIndex(j),:)),blockIndex(j),1);
                elseif j==length(blockIndex)
                    varargin{i}.variance(blockIndex(j)+1:completedTrials,:)=repmat(var(varargin{i}.RewardMagnitude(blockIndex(j)+1:completedTrials,:)), completedTrials-blockIndex(j),1);
                else
                    varargin{i}.variance(blockIndex(j-1)+1:blockIndex(j)+1,:)=repmat(var(varargin{i}.RewardMagnitude(blockIndex(j-1)+1:blockIndex(j),:)), blockIndex(j)+1-blockIndex(j-1),1);
                end 
            end

            %blockType %2 columns - 1 reward 2 variance
            varargin{i}.blockType=zeros(length(unique(varargin{i}.BlockNumberL)),2);
            %left large = 1. right large=2, control=3
            varargin{i}.blockType(varargin{i}.leftLarge([blockIndex, completedTrials]),1)=1;
            varargin{i}.blockType(varargin{i}.rightLarge([blockIndex, completedTrials]),1)=2;
            varargin{i}.blockType(varargin{i}.controlMag([blockIndex, completedTrials]),1)=3;

            %var large = 1. var med=2, control=3
            xx=varargin{i}.variance([blockIndex, completedTrials],:);

            if (max(max(xx))~=0)
                varargin{i}.blockType(:,2)=2; %default is medium trial - replace if highest var or control var
                varargin{i}.blockType(sum(xx==max(max(xx)),2),2)=1;
                varargin{i}.blockType(xx(:,1)==xx(:,2),2)=3;
            else
                varargin{i}.blockType(:,2)=nan;
                disp('no variance in reward blocks')
            end






            %makae choice data into matrix, sep by blocks
            blockIndex=find(varargin{i}.plotBlockIdx==1);
             varargin{i}.choiceLeftMat=zeros(250, length(blockIndex)+1);
                for xi=1:length(blockIndex)+1
                    if xi==1
                        varargin{i}.choiceLeftMat(1:length(varargin{i}.plotChoiceLeft(1:blockIndex(xi))),xi)=varargin{i}.plotChoiceLeft(1:blockIndex(xi));
                    elseif xi==length(blockIndex)+1
                        varargin{i}.choiceLeftMat(1:length(varargin{i}.plotChoiceLeft(blockIndex(end)+1:end)),xi)=varargin{i}.plotChoiceLeft(blockIndex(end)+1:end);
                    else
                        varargin{i}.choiceLeftMat(1:length(varargin{i}.plotChoiceLeft(blockIndex(xi-1)+1:blockIndex(xi))),xi)=varargin{i}.plotChoiceLeft(blockIndex(xi-1)+1:blockIndex(xi));
                    end

                end

               varargin{i}.blockTransMat=NaN(111, length(blockIndex));
               try
                for xi=1:length(blockIndex)
                  varargin{i}.blockTransMat(1:length(varargin{i}.plotChoiceLeft(blockIndex(xi)-50:blockIndex(xi)+60)),xi+1)=varargin{i}.plotChoiceLeft(blockIndex(xi)-50:blockIndex(xi)+60);
                end
               catch
                   if (blockIndex(xi)-50)<0
                       disp('not enough trials in block - starting index is negative')
                   elseif (blockIndex(xi)+60)>length(varargin{i}.plotChoiceLeft)
                       disp('not enough trials in block - end index not reached')
                       varargin{i}.blockTransMat=horzcat(varargin{i}.blockTransMat, NaN(111, 1));
                   end
               end
        else
            disp('there are no blocks detected')

        end



        % merge structs
        superStruct=cat(2, varargin{:});
        fields = fieldnames(superStruct);

        for k = 1:numel(fields)
          aField     = fields{k}; % EDIT: changed to {}
          try
          S.(aField) = cat(2, superStruct(:).(aField));
          catch
          S.(aField) = cat(1, superStruct(:).(aField));

          end

        end
end
end

function [blockMeanX, blockMeanY]=avgPsyc(XBlock, YBlock, blockType, blockNum)
    XXBlock=cell2mat(XBlock(blockType==blockNum));
    XXBlock=reshape(cell2mat(XXBlock(blockType==blockNum)), 100,[]);
    blockMeanX=mean(XXBlock,2);
    
    YYBlock=cell2mat(YBlock(blockType==blockNum));
    YYBlock=reshape(cell2mat(YYBlock(blockType==blockNum)), 100,[]);
    blockMeanY=mean(YYBlock,2);
end



