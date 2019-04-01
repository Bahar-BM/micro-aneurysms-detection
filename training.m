clc;
clear all;
countertable=1;
xDoc = xmlread('annotations-consensus-ma-only.xml');
strc  = xml2struct( xDoc );
minn=zeros(1,31);
maxx=zeros(1,31);
minn(1:31)=100000;
maxx(1:31)=-100000;

for jj=0:49
d=strc.set.annotations_dash_per_dash_image(jj+1);
d=d{1,1};

if isfield(d,'annotation')
    
inputt=double(imread(strcat('training\image',int2str(jj),'_training.jpg')));
input= inputt(:, :, 2);

input_img=imread(strcat('normalization\normalization',int2str(jj),'.jpg'));
input1=double(input_img);
NormInput= input1(:, :, 2);

[R C]=size(input);
bloodVessels=mat2gray(imread(strcat('modmask\modmask',int2str(jj),'.jpg')));

[output4,Maxrr,Minrr,Avgrr]=CandidateMicroExtract(NormInput,bloodVessels,input);
[stats,co]=FeatureExtract(input1,output4,input,Maxrr,Minrr,Avgrr);

load('FeatureFile.mat');

[sz ~]=size(a);
if sz~=0        
     minn(1)=min(min(a),minn(1));
     maxx(1)=max(max(a),maxx(1));
     minn(2)=min(min(per),minn(2));
     maxx(2)=max(max(per),maxx(2)); 
     minn(3)=min(min(MrM),minn(3));
     maxx(3)=max(max(MrM),maxx(3));
     minn(4)=min(min(circularity),minn(4));
     maxx(4)=max(max(circularity),maxx(4));
     minn(5)=min(min(i_green),minn(5));
     maxx(5)=max(max(i_green),maxx(5));
     minn(6)=min(min(i_sc),minn(6));
     maxx(6)=max(max(i_sc),maxx(6));
     minn(7)=min(min(m_green),minn(7));
     maxx(7)=max(max(m_green),maxx(7));
     minn(8)=min(min(m_sc),minn(8));
     maxx(8)=max(max(m_sc),maxx(8));
     minn(9)=min(min(NI_green),minn(9));
     maxx(9)=max(max(NI_green),maxx(9));
     minn(10)=min(min(NI_sc),minn(10));
     maxx(10)=max(max(NI_sc),maxx(10));
     minn(11)=min(min(NM_green),minn(11));
     maxx(11)=max(max(NM_green),maxx(11));
     minn(12)=min(min(NM_sc),minn(12));
     maxx(12)=max(max(NM_sc),maxx(12));
     minn(13)=min(min(IDarkest),minn(13));
     maxx(13)=max(max(IDarkest),maxx(13));
     minn(14)=min(min(FiltSig1),minn(14));
     maxx(14)=max(max(FiltSig1),maxx(14));
     minn(15)=min(min(FiltSig2),minn(15));
     maxx(15)=max(max(FiltSig2),maxx(15));
     minn(16)=min(min(FiltSig4),minn(16));
     maxx(16)=max(max(FiltSig4),maxx(16));
     minn(17)=min(min(FiltSig8),minn(17));
     maxx(17)=max(max(FiltSig8),maxx(17));
     minn(18)=min(min(StdFiltSig1),minn(18));
     maxx(18)=max(max(StdFiltSig1),maxx(18));
     minn(19)=min(min(StdFiltSig2),minn(19));
     maxx(19)=max(max(StdFiltSig2),maxx(19));
     minn(20)=min(min(StdFiltSig4),minn(20));
     maxx(20)=max(max(StdFiltSig4),maxx(20));
     minn(21)=min(min(StdFiltSig8),minn(21));
     maxx(21)=max(max(StdFiltSig8),maxx(21));
     minn(22)=min(min(MaxCorr),minn(22));
     maxx(22)=max(max(MaxCorr),maxx(22));
     minn(23)=min(min(MinCorr),minn(23));
     maxx(23)=max(max(MinCorr),maxx(23));
     minn(24)=min(min(AvgCorr),minn(24));
     maxx(24)=max(max(AvgCorr),maxx(24));
     minn(26)=min(min(diffr),minn(26));
     maxx(26)=max(max(diffr),maxx(26));
     minn(27)=min(min(diffg),minn(27));
     maxx(27)=max(max(diffg),maxx(27));
     minn(28)=min(min(diffb),minn(28));
     maxx(28)=max(max(diffb),maxx(28));
     minn(29)=min(min(diffh),minn(29));
     maxx(29)=max(max(diffh),maxx(29));
     minn(30)=min(min(MajorAxis),minn(30));
     maxx(30)=max(max(MajorAxis),maxx(30));
    minn(31)= min(min(MinorAxis),minn(31));
    maxx(31)= max(max(MinorAxis),maxx(31));
end
end
end
save('minmax77file.mat','minn','maxx');