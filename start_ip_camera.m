%starts vlc rtps stream specified in am xml file for kepecs' lab rat rigs.
%specific to default vlc installation folder on Win10 and kepecs rat rig
%xml file.
%Torben Ott, CSHL, 2017

%hostname
[~,hostname] = system('hostname');
hostname=strcat(hostname);%hack to remove spaces

%get ip address from xml file
XML = xmlread('http://kepecsdata.cshl.edu/observer/rigdata/rigs.xml');
rigs = XML.getElementsByTagName('rig');
address = [];
for k =0:rigs.getLength-1
    rig=rigs.item(k);
    name = char(rig.getElementsByTagName('name').item(0).getFirstChild.getData);
    if strcmpi(name,hostname)
        address = char(rig.getElementsByTagName('address').item(0).getFirstChild.getData);
        break
    end
end

%find vlc
if exist('C:\Program Files (x86)\VideoLAN\VLC\VLC.exe','file')==2
    path = 'C:\Program Files (x86)\VideoLAN\VLC\VLC.exe';
elseif exist('C:\Program Files\VideoLAN\VLC\VLC.exe','file')==2
    path = 'C:\Program Files\VideoLAN\VLC\VLC.exe';
% elseif exist('C:\Program Files\VideoLAN\VLC\VLC.exe','file')==2
else
    fprintf('Starting VLC from command failed.');
    return
end

%start vlc
try
    system(['start "" "',path,'" ', address]);
catch
        fprintf('Starting VLC from command failed.');
end