function varargout=LoadFeatures(varargin)
% scan the given filename containing the feature vectors in the form:
% value value value ... value img_XXXX.jpg
% value value value ... value img_XX.jpg
% ...
for k=1:length(varargin)
    filename=varargin{k};
    fid = fopen(filename);
    %scan the whole file
    C=textscan(fid,'%[^\n]',25000);
    %number of observations (images)
    M=size(C{1},1);
    %remove strings 'img_' and '.jpg'
    for i=1:M
        C{1}{i}=strrep(strrep(strrep(C{1}{i},'img_',''),'.jpg',''),'.png','');
    end
    %number of individual features (+1 for image number within img_XXXX.jpg)
    N=size(str2num(C{1}{1}),2);
    %convert to the final form (double)
    varargout{k}=zeros(M,N);
    for i=1:M
        varargout{k}(i,:)=str2num(C{1}{i});
    end
    fclose(fid);
end