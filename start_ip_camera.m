%starts vlc rtps stream specified in am xml file for kepecs' lab rat rigs.
%specific to default vlc installation folder on Win10 and kepecs rat rig
%xml file.
%Torben Ott, CSHL, 2017

%hostname
[~,hostname] = system('hostname');
hostname=strcat(hostname);%hack to remove spaces

%get ip address from xml file
XML = xmlread('http://kepecsdata.cshl.edu/observer/rigdata/rigs.xml');
rigs = XML.getDocumentElement;

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

%start vlc
system(['start "" "C:\Program Files\VideoLAN\VLC\VLC.exe" ', address]);