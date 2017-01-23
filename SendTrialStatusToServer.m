function SendTrialStatusToServer(script,rig,outcome)
%Calls PHP script at HTTP location specified in SCRIPT [string].
%Gives two variables: RIG [string] name and OUTCOME [numerical vector of integers from 1-9, max length 30].
%Torben Ott, CSHL, 2017
if length(outcome)>30
    outcome=outcome(end-29:end);
end
length(outcome)
url=strcat('http://kepecsdata.cshl.edu/scripts/',script,'?rig=',rig,'&outcome=',num2str(outcome,'%i'));
webread(url);
end