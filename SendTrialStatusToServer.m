function SendTrialStatusToServer(script,rig,outcome,subject,protocol)
%Calls PHP script at HTTP location specified in SCRIPT [string].
%Two required variables: RIG [string] name and OUTCOME [numerical vector of integers from 1-9].
%Two optional variables: SUBJECT [string] name and PROTOCOL [string] name.
%Torben Ott, CSHL, 2017
if nargin<5
    protocol='na';
end
if nargin<4
    subject='na';
end
if nargin<3
    error('SendTrialStatusToServer:Missing input argument.');
end
    
if length(outcome)>30
    outcome=outcome(end-29:end);
end
outcome(isnan(outcome)) = 9;
url=strcat('http://kepecsdata.cshl.edu/observer/scripts/',script,'?rig=',rig,'&outcome=',num2str(outcome,'%i'),'&subject=',subject,'&protocol=',protocol);
webread(url);
end