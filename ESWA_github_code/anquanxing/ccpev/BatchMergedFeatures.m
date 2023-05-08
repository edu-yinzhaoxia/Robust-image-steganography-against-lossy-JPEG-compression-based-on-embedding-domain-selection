function BatchMergedFeatures(DIR)
% This function extracts merged  extended DCT and Markov  features from all
% JPEG images in the given directory DIR, both original features introduced
% in [1]  (274 features)  and  Cartesian calibrated  introduced in [2] (548
% features).
%
% Comments/remarks: Jan Kodovsky, jan@kodovsky.com
%
% [1] T. Pevny and J. Fridrich.  Merging Markov and DCT features for multi-
%     class JPEG steganalysis. In E. J. Delp and P. W. Wong,  editors, Pro-
%     ceedings  SPIE,  Electronic  Imaging,  Security,  Steganography,  and
%     Watermarking of Multimedia Contents IX,  volume 6505, pages 3 1–3 14,
%     San Jose, CA, January 29–February 1, 2007.
% [2] J. Kodovsky and J. Fridrich.  Calibration revisited.  In J. Dittmann,
%     S. Craver, and J. Fridrich, editors, Proceedings of the 11th ACM Mul-
%     timedia & Security Workshop, Princeton, NJ, September 7–8, 2009.
%
% INPUT: DIR - directory with the JPEG images
% OUTPUT FILES:
%      'data/merged.fea'     - 274 features
%      'data/cc-merged.fea'  - 548 features

%create the output directory if it doesn't exist
if ~exist( './data','dir');mkdir('./data');end

%create output file
outputFileFcc='./data/cc-merged.fea';
outputFileFdiff= './data/merged.fea';

fid=fopen(outputFileFcc,'w');fclose(fid);
fid=fopen(outputFileFdiff,'w');fclose(fid);

files=dir([DIR '/*.jpg']);
for i=1:length(files)
    fprintf('processing %s: ',files(i).name);
    if files(i).isdir==0
        filename=[DIR '/' files(i).name];
        try
            F=CalculateMergedFeatures_noCalibration(filename);
            F2=CalculateMergedFeatures_onlyCalibration(filename);
            Fcc=[F;F2];
            Fdiff=F-F2;
            
            %output F diff (see publication [1])
            fid=fopen(outputFileFdiff,'a');
			for j=1:length(Fdiff)
					fprintf(fid,'%+.8e ',Fdiff(j));
			end
			fprintf(fid,'%s\n',files(i).name);
            fclose(fid);
            
            %output F cc (see publication [2])
            fid=fopen(outputFileFcc,'a');
			for j=1:length(Fcc)
					fprintf(fid,'%+.8e ',Fcc(j));
			end
			fprintf(fid,'%s\n',files(i).name);
            fclose(fid);
            fprintf('ok\n');
        catch
          errmsg = lasterr;
          fprintf('%s\n',errmsg);
        end
    end %if files(i).isdir
end %for

fprintf('End of processing.\n');

function F=CalculateMergedFeatures_noCalibration(filename)

jobj=jpeg_read(filename);
OriginalDct=PlaneToVec(jobj.coef_arrays{1});
QuantTable=jobj.quant_tables{jobj.comp_info(1).quant_tbl_no};
SpatialImage=double(Dct2Spat(OriginalDct,QuantTable));

[ExtendedDCTF]=ExtractExtendedDCTFeatures_noCalibration(OriginalDct,SpatialImage);
[OrMh,OrMv,OrMd,OrMm]=ExtractCalibratedMarkovFeatures_noCalibration(OriginalDct);

%create Merged Features
MarkovPart=0.25*(OrMh+OrMv+OrMd+OrMm);
DctPart=zeros(1,193);
DctPart(1:165)=ExtendedDCTF(1:165);
DctPart(166)=0.5*(ExtendedDCTF(166)+ExtendedDCTF(167));
DctPart(167:168)=ExtendedDCTF(168:169);
for j=0:24
    DctPart(169+j)=0.5*(ExtendedDCTF(170+j)+ExtendedDCTF(195+j));
end
F=[MarkovPart(:);DctPart(:)];

%modification in order to correspond with c++ implementation
F=F([82:274 1:81]);
F(67:165)=-F(67:165);

function F=CalculateMergedFeatures_onlyCalibration(filename)

jobj=jpeg_read(filename);
OriginalDct=PlaneToVec(jobj.coef_arrays{1});
QuantTable=jobj.quant_tables{jobj.comp_info(1).quant_tbl_no};
SpatialImage=double(Dct2Spat(OriginalDct,QuantTable));

% crop by 4x4 pixels
CroppedDct=DCTcut(SpatialImage(5:end-4,5:end-4),QuantTable);
SpatialCropped=double(Dct2Spat(CroppedDct,QuantTable));

[ExtendedDCTF]=ExtractExtendedDCTFeatures_noCalibration(CroppedDct,SpatialCropped);
[OrMh,OrMv,OrMd,OrMm]=ExtractCalibratedMarkovFeatures_noCalibration(CroppedDct);

%create Merged Features
MarkovPart=0.25*(OrMh+OrMv+OrMd+OrMm);
DctPart=zeros(1,193);
DctPart(1:165)=ExtendedDCTF(1:165);
DctPart(166)=0.5*(ExtendedDCTF(166)+ExtendedDCTF(167));
DctPart(167:168)=ExtendedDCTF(168:169);
for j=0:24
    DctPart(169+j)=0.5*(ExtendedDCTF(170+j)+ExtendedDCTF(195+j));
end
F=[MarkovPart(:);DctPart(:)];

%modification in order to correspond with c++ implementation
F=F([82:274 1:81]);
F(67:165)=-F(67:165);

function [OrMh,OrMv,OrMd,OrMm]=ExtractCalibratedMarkovFeatures_noCalibration(OriginalDct)
%[F]=ExtractMarkovFeatures(image_name);
%[OriginalDct,SpatialImage,CalibratedDct,SpatialCropped]=calibrateImage(filename);

[OrMh,OrMv,OrMd,OrMm]=calculateMatrices(vecToPlane(OriginalDct));

function [Mh,Mv,Md,Mm]=calculateMatrices(Plane)
absDctPlane=abs(Plane);

Fh=absDctPlane(:,1:end-1) - absDctPlane(:,2:end);
Fv=absDctPlane(1:end-1,:) - absDctPlane(2:end,:);
Fd=absDctPlane(1:end-1,1:end-1) - absDctPlane(2:end,2:end);
Fm=absDctPlane(2:end,1:end-1) - absDctPlane(1:end-1,2:end);


Mh=calculateMarkovField(Fh(:,1:end-1),Fh(:,2:end),4);
Mv=calculateMarkovField(Fv(1:end-1,:),Fv(2:end,:),4);
Md=calculateMarkovField(Fd(1:end-1,1:end-1),Fd(2:end,2:end),4);
Mm=calculateMarkovField(Fm(2:end,1:end-1),Fm(1:end-1,2:end),4);

function field=calculateMarkovField(reference,shifted,T)
field=zeros(2*T+1,2*T+1);
R=reference(:);
S=shifted(:);
R(find(R>T))=T;
R(find(R<-T))=-T;
S(find(S>T))=T;
S(find(S<-T))=-T;

for i=-T:T
    idx=find(R==i);
    h=hist(S(idx),-T:T);
    if length(idx)>0
        h=h/length(idx);
    end
    field(i+T+1,:)=h;
end

function plane=vecToPlane(vec)
[N M Z]=size(vec);
plane=zeros(N*8,M*8);
for i=1:N
    for j=1:M
		plane(((i-1)*8+1):i*8,((j-1)*8+1):j*8)=reshape(vec(i,j,:),8,8);
    end
end

function [F]=ExtractExtendedDCTFeatures_noCalibration(OriginalDct,SpatialImage)
% [F]=ExtractFeatures(image_name);


k=1;
F=zeros(194,1);
%----------------GLOBAL HISTOGRAM FEATURE---------------------
%firstly global histogram will be computed
R=max([max(OriginalDct(:)),5]);
L=min([min(OriginalDct(:)),-5]);

H_O=hist(OriginalDct(:),L:R);
if (sum(H_O)~=0)
  H_O=H_O/sum(H_O);       % Normalized histogram of x
end;
T=H_O;
F(k:k+10)=T(-4-L:-4-L+10);
k=k+11;

%--------------LOCAL HISTOGRAM FEATURE------------------------
%Local histograms for individual DCT modes
%statisticaly important are modes (2,1)=2 (3,1)=3 (1,2)=9, (2,2)=10 (1,3)=17
modes=[2 3 9 10 17];
for mode=modes
  Original=OriginalDct(:,:,mode);
  R_i=max([max(Original(:)),5]);
  L_i=min([min(Original(:)),-5]);
  H_oi=hist(Original(:),L_i:R_i);
  %normalization of histograms
  if (sum(H_oi)~=0)
    H_oi=H_oi/sum(H_oi);
  end;
  T=H_oi;
  F(k:k+10)=T(-4-L_i:-4-L_i+10);
  k=k+11;
end

%---------------DUAL HISTOGRAM FEATURE-------------------------
%Dual histograms in range -5 to 5
values=(-5:5);
idxs=[2,9,3,10,17,4,11,18,25];
for value=values
  h_o=zeros(1,64);
  for i=1:64
    h_o(i)=sum(sum(OriginalDct(:,:,i)==value));
  end
  if (sum(h_o)~=0)
    h_o=h_o/sum(h_o); %norm
  end;
%% Q2
  T=-h_o;
  F(k:k+8)=T(idxs);
  k=k+9;
end;

%----------------VARIATION FEATURE------------------------------

%<MOD>
%05-02-2007: modification of normalization factors in order to correspond
%            with our c++ implementation
%Horizontal Variation
Dw=abs(OriginalDct(1:end,1:end-1,:)-OriginalDct(1:end,2:end,:));
V_h=sum(Dw(:))/(size(OriginalDct,1)*(size(OriginalDct,2)-1)*64);

%vertical Variation
Dw=abs(OriginalDct(1:end-1,1:end,:)-OriginalDct(2:end,1:end,:));
V_v=sum(Dw(:))/(size(OriginalDct,2)*(size(OriginalDct,1)-1)*64);

%5-01-2007: end of modification
%</MOD>
%the feature is
F(k)=V_h;k=k+1;
F(k)=V_v;k=k+1;

%----------------BLOCKINESS FEATURE-----------------------------
%for blockiness the spatial representation of the image is used
LO=blockiness(SpatialImage);
F(k)=LO(1);k=k+1;
F(k)=LO(2);k=k+1;


%-------------CO OCCURENCE MATRIX--------------------------------
dummy=coocurrence(OriginalDct);
F(k:k+length(dummy)-1)=dummy;


% function for computing coocurence feature
% copied from jessica's routiness
function F=coocurrence(Original)
% KOD, 9-9-2008: vertical part - optimized
desiredRange=-2:2;
x=max(min(reshape(Original(1:end-1,:,2:end),1,numel(Original(1:end-1,:,2:end))),max(desiredRange)+1),min(desiredRange)-1);
y=max(min(reshape(Original(2:end,:,2:end),1,numel(Original(2:end,:,2:end))),max(desiredRange)+1),min(desiredRange)-1);
L=min([x min(desiredRange)-1]);
R=max([x max(desiredRange)+1]);
Dim=R-L+1;
B=length(x);
ind=sub2ind([Dim Dim],x-L+1,y-L+1);

Mm=ind2sub([Dim Dim],hist(ind,1:Dim*Dim))/B;
Mm=reshape(Mm,Dim,Dim);

Diff=Mm(-L+1+min(desiredRange):-L+1+max(desiredRange),-L+1+min(desiredRange):-L+1+max(desiredRange));

% horizontal part
x=max(min(reshape(Original(:,1:end-1,2:end),1,numel(Original(:,1:end-1,2:end))),max(desiredRange)+1),min(desiredRange)-1);
y=max(min(reshape(Original(:,2:end,2:end),1,numel(Original(:,2:end,2:end))),max(desiredRange)+1),min(desiredRange)-1);

ind=sub2ind([Dim Dim],x-L+1,y-L+1);
B=length(x);
Mm=ind2sub([Dim Dim],hist(ind,1:Dim*Dim))/B;

Mm=reshape(Mm,Dim,Dim);

Diffh=Mm(-L-1:-L+3,-L-1:-L+3);

F=[Diff(:);Diffh(:)];


function [L]=blockiness(SI)
[N M]=size(SI);
C=N*floor((M-1)/8)+M*floor((N-1)/8);
%horizontal lines
H_1=SI(8:8:N-1,:);
H_2=SI(9:8:N,:);
Dw_H=abs(H_1(:)-H_2(:));
%vertical lines
V_1=SI(:,8:8:M-1);
V_2=SI(:,9:8:M);
Dw_V=abs(V_1(:)-V_2(:));
%L1 blockiness
L(1)=(sum(Dw_H)+sum(Dw_V))/C;
L(2)=(sum(Dw_H.^2)+sum(Dw_V.^2))/C;

function X=jpg2raw02(QD,Qf)

% "Inverse" function to raw2jpg02.m, calculates the spatial representation for an 8x8 block
% of DCT coefficients QD quantized with matrix Qf
%
% Input:
% QD = Matrix of DCT coefficients quantized with the quantization matrix Qf (QD is 8x8 of integers)
% Qf = Quantization matrix (8x8 of quantization factors)
%
% Output:
% X  = Spatial representation of QD = DCT^(-1)(QD) rounded to integers and at 0 and 255 

% 128 is added because in calculating the DCT QD coeffs, 128 is subtracted (see raw2jpg01.m)

  X=128+round(idct2(QD.*Qf));	% Rounding to integers
  X2 = X;
  X(find(X2<0))=0;		% Truncating at 0
  X(find(X2>255))=255;		% Truncating at 255
%  X=(idct(QD',Qf'))';
  
function Mat=PlaneToVec(plane)
[Y X]=size(plane);
M=floor(X/8);	% count of columns
N=floor(Y/8);  	%count of rows
Mat=zeros(N,M,64);
for i=1:N
    for j=1:M
        Mat(i,j,:)=reshape(plane(((i-1)*8+1):i*8,((j-1)*8+1):j*8),64,1);
    end
end

function Plane=VecToPlane(Cube)
[Y X Z]=size(Cube);
M=floor(X*8);   % count of columns
N=floor(Y*8);   %count of rows
Plane=zeros([N,M]);
for i=1:Y
    for j=1:X
        Plane(((i-1)*8+1):i*8,((j-1)*8+1):j*8)=reshape(Cube(i,j,:),8,8);
    end
end

function op=Dct2Spat(X,Q)

Step1=VecToPlane(X);
op=blkproc(Step1,[8 8],@jpg2raw02,Q);
op=uint8(op); 

function D=DCTcut(X,Q)
% Calculates the DCT coefficients of the image/matrix X using quantization factor/matrix Q
% Usage: D=dctcut02(X,Q), where X could be a (0-255) grayscale matrix an image file name, such as 'len256.bmp', and Q can be either a
%       quality factor or a quantization matrix
% Note: If X is a true-color image, it will be converted to grayscale before further processing
% Note: The input images could be either JPG or non-palette BMP
% Note: Do not process palette images with this routine, the routing will screw up if you do
% Note: If the image dimensions are not multiples of 8, the image is cropped to the closest multiples of 8

s=size(Q);
if s(1)==1,Q=Qmatrix(Q);end

%  s=size(X);
%  if s(1)==1;[X,Map]=imread(X);end

%  s=size(X);
%  if length(s)==3,X=rgb2gray(X);end
X=double(X);

B=8;B2=B^2;

[M,N]=size(X);
MB=floor(M/B);
NB=floor(N/B);

D=zeros(MB,NB,B2);					% DCT coefficients of one block, D(i,j,k), (i,j) are block indices, k = 1, ..., 64

for i=1:MB
   for j=1:NB
     ib=(i-1)*B+1;
     jb=(j-1)*B+1;
     Block=X(ib:ib+7,jb:jb+7);		% B x B image block
     Dblock=raw2jpg02(Block,Q);	    % Its DCT transform
     D(i,j,:)=Dblock(:);
  end
end

function QD=raw2jpg02(Z,Qf)

% "Inverse" function to jpg2raw02.m, calculates DCT coefficients quantized with matrix Qf
% from an 8x8 spatial block Z
%
% Input:
% Z = Matrix of grayscales (QD is 8x8 of integers between 0 and 255)
% Qf = Quantization matrix (8x8 of quantization factors)
%
% Output:
% QD  = Quantized DCT transform of Z, QD = quantize(DCT(Z)) with quantization matrix Qf 

  Z=double(Z)-128;
  D=dct2(Z);  
  QD=round(D./Qf);
  %D=DCT(Z');
  %QD=round(D'./(Qf*8));  %the result of JPG DCT is 8 times larger than the normal DCT
  