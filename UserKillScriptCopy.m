%User Kill Script
%this script, when saved in the current protocol folder, will be executed
%at Bpod termination.
%Serves for custom and protocol specific action at session end.

%% PUT USER-SPECIFIC STUFF HERE%%
% For example, save behavior figures, do additional analyses and send
% results by e-mail. For example scripts like this, see 
% https://github.com/torbencshl/BpodFunctions/

%% copy data to server
try
    %CHANGE THIS TO 'homes\YOUR-ACCOUNT-NAME' OR to 'SHARED-FOLDER-NAME'
    %examples: user = 'homes\torben'; % user name on server
    %          user = 'confidence';   % shared folder
    user = strcat('homes\',getenv('username'));
    %%%%%

    [~,subject] = fileparts(fileparts(fileparts(fileparts(BpodSystem.DataPath))));
    if ~isdir(fullfile('\\kepecsdata',user,'BpodData',BpodSystem.CurrentProtocolName,subject,'Session Data'))
        mkdir(fullfile('\\kepecsdata',user,'BpodData',BpodSystem.CurrentProtocolName,subject,'Session Data'));
    end
    if ~isdir(fullfile('\\kepecsdata',user,'BpodData',BpodSystem.CurrentProtocolName,subject,'Session Settings'))
        mkdir(fullfile('\\kepecsdata',user,'BpodData',BpodSystem.CurrentProtocolName,subject,'Session Settings'));
    end    
    copyfile(BpodSystem.DataPath,fullfile('\\kepecsdata',user,'BpodData',BpodSystem.CurrentProtocolName,subject,'Session Data'));
    copyfile(BpodSystem.SettingsPath,fullfile('\\kepecsdata',user,'BpodData',BpodSystem.CurrentProtocolName,subject,'Session Settings'));
catch
    fprintf('Error copying data to server. Files not copied!\n');
end