function MailAlert(persons,subject,msg)
%persons is a cell array with string containing recipient names
if iscell(persons)
    for i = 1:length(persons)
        person=persons{i};
        switch person
            case 'Torben'
                address = 'm6a1c5v4i4m5t5t0@kepecslab.slack.com';
            case 'Paul'
                address = 'n1g4g9a6v7u9s9s8@kepecslab.slack.com';
            otherwise
                fprintf('Add person to mail list on MailAlert.m\n');
                return
        end
        
        SendMyMail(address,subject,msg);
        
    end
end
end

% sends mail from Kepecslab's cshl gmail account
% 3 or  4 inputs: address,subject,message,cell with attachment paths
% (each as string)
function sent = SendMyMail(varargin)
sent = false;
setpref('Internet','E_mail','kepecslab.cshl@gmail.com')
setpref('Internet','SMTP_Server','smtp.gmail.com')
setpref('Internet','SMTP_Username','kepecslab.cshl@gmail.com')
setpref('Internet','SMTP_Password','D3cision')
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

if length(varargin)==3
    try
        sendmail(varargin{1},varargin{2},varargin{3})
        sent=true;
    catch
        display('Error:SendMyMail:E-Mail could not be sent.')
    end
elseif length(varargin)==4
    try
        sendmail(varargin{1},varargin{2},varargin{3},varargin{4})
        sent=true;
    catch
        display('Error:SendMyMail:E-Mail could not be sent.')
    end
else
    display('Error:SendMyMail:Number of input arguments wrong.')
end

end