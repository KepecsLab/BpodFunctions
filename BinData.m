function [x,y,e,n,v]=BinData(datax,datay,edges,opt)


if nargin<4
    opt=struct;
end
if ~isfield(opt,'unique')
    opt.unique=false;
end
if isfield(opt,'bin_edge') && opt.bin_edge
    idx = discretize(datax,[-1/eps,edges,1/eps]);
    idx(idx==1)=2;
    idx(idx==length(edges)+1) = length(edges);
    idx = idx-1;
else


end

if  opt.unique
    idx=round(datax*1000)/1000; %resolution for unique values
    x=unique(idx);
    x=x(~isnan(x));
else
    idx = discretize(datax,edges);
    x = mean([edges(1:end-1);edges(2:end)]);

end

[y_temp,n_temp,gnames] = grpstats(datay,idx,{'mean','numel','gname'});
e_temp = grpstats(datay,idx,'sem');
if length(e_temp)<length(y_temp)
    e_temp=nan(size(y_temp));
end

if ~isempty(gnames) && ~opt.unique
    gnames = cellfun(@str2double,gnames);
    
    
    %bins with no data
    missing = setxor(1:length(edges)-1,gnames);
    missingdata = nan(length(missing),1);
    
    %fill in nan data
    namesfull = [gnames;missing(:)];
    y = [y_temp;missingdata];
    y = sortrows([y,namesfull],2); %sort
    y=y(:,1);
    
    e = [e_temp;missingdata];
    e = sortrows([e,namesfull],2); %sort
    e=e(:,1);
    
    n=[n_temp;zeros(length(missing),1)];
    n = sortrows([n,namesfull],2); %sort
    n=n(:,1);
    n=n';
    
    v = e.*sqrt(n(:));
elseif opt.unique
    %have to do it differently
    y=y_temp;
    n=length(y);
    e=e_temp;
    v = e.*sqrt(n);
else%no data at all
    y = nan(length(edges)-1,1);
    e = nan(length(edges)-1,1);
    n=zeros(length(edges)-1,1);
end

x=x(:); y=y(:);e=e(:);

end