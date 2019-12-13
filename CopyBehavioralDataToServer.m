BPODUSER = 'C:\User\torben\Documents\BpodUser\';
PROTOCOL='Dual2AFC';
SERVER = '\\uncertainty.cshl.edu\home\';

%%

%animals
animals = getDir(fullfile(BPODUSER,'Data',PROTOCOL),'folder');

for a=1:length(animals)
    an = animals{a};
    files = getDir(fullfile(BPODUSER,'Data',an,PROTOCOL),'file');
    for f =1:length(files)
        local = fullfile(BPODUSER,'Data',an,PROTOCOL,'Session Data',files{f});
        serv = fullfile(SERVER,'BpodData',an,PROTOCOL,'Session Data',files{f});
        if ~exist(serv,'file')==2
            copyfile(local,serv);
        end
    end
end