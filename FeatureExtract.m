function [stats,co]=FeatureExtract(input1,output4,input,Maxrr,Minrr,Avgrr)
[R C]=size(input);
Backgrond = medfilt2(input, [25 25]);       

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%% Feature Extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L,co] = bwlabel(output4);
stats = regionprops(L,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%% Feature 1 (a)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = [stats.Area];
%%%%%%%%%%%%%%%%%%%%%%%%%% Feature 2 (p)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
per = [stats.Perimeter];
%%%%% Feature 3 (r),4 (c),5 (i_green),7 (m_green),8 (m_sc)%%%%%%%%%%%%%%%%%
%%%%% 9 (NI_green),10 (NI_SC),11 (NM_green),12 (NM_sc),13 (I-darkest)%%%%%%
%%%%% 19 (FiltSig1),20 (FiltSig2),21 (FiltSig4),22 (FiltSig8)%%%%%%%%%%%%%%
%%%%% Feature 23 (StdFiltSig1),Feature 24 (StdFiltSig2)%%%%%%%%%%%%%%%%%%%%
%%%%% Feature 25 (StdFiltSig4),Feature 26 (StdFiltSig8)%%%%%%%%%%%%%%%%%%%%
%%%%% Feature 27 (MaxCorr),Feature 28 (MinCorr),Feature 29 (AvgCorr)%%%%%%%
MrM=zeros(1,co);
circularity=zeros(1,co);
i_green=zeros(1,co);
SC=input-Backgrond;
i_sc=zeros(1,co);
m_green=zeros(1,co);
m_sc=zeros(1,co);
NI_green=zeros(1,co);
p = Backgrond';
b = p(:)';
NI_sc=zeros(1,co);
Sig=std(b,1);
Mean=mean(b);
NM_green=zeros(1,co);
NM_sc=zeros(1,co);
IDarkest=zeros(1,co);
FiltSig1=zeros(1,co);
h = fspecial('gaussian',[7,7], 1);
gaus1= imfilter(input,h);
FiltSig2=zeros(1,co);
h = fspecial('gaussian',[13,13], 2);
gaus2= imfilter(input,h);
FiltSig4=zeros(1,co);
h = fspecial('gaussian',[25,25], 4);
gaus4= imfilter(input,h);
FiltSig8=zeros(1,co);
h = fspecial('gaussian',[33,33], 8);
gaus8= imfilter(input,h);
StdFiltSig1=zeros(1,co);
h = fspecial('gaussian',[7,7], 1);
gaus11= imfilter(input,h);
StdFiltSig2=zeros(1,co);
h = fspecial('gaussian',[13,13], 2);
gaus22= imfilter(input,h);
StdFiltSig4=zeros(1,co);
h = fspecial('gaussian',[25,25], 4);
gaus44= imfilter(input,h);
StdFiltSig8=zeros(1,co);
h = fspecial('gaussian',[33,33], 8);
gaus88= imfilter(input,h);
MaxCorr=zeros(1,co);
MinCorr=zeros(1,co);
AvgCorr=zeros(1,co);


for i=1:co
 MrM(1,i)=stats(i).MajorAxisLength/stats(i).MinorAxisLength;
 circularity(1,i)=(4*pi*a(1,i))/(per(1,i)^2);
 vl=L==i;
 i_green(i)=sum(sum(input(vl)));
 i_sc(i)=sum(sum(SC(vl)));
 m_green(1,i)=i_green(1,i)/a(1,i);
 m_sc(1,i)=i_sc(1,i)/a(1,i);
 NI_green(1,i)=(1/Sig)*(i_green(1,i)-Mean);
 NI_sc(1,i)=(1/Sig)*(i_sc(1,i)-Mean);
 NM_green(1,i)=(1/Sig)*(m_green(1,i)-Mean);
 NM_sc(1,i)=(1/Sig)*(m_sc(1,i)-Mean);
 IDarkest(i)=min(min(input(vl)));
 sums=sum(sum(vl));
 FiltSig1(i) = sum(sum(gaus1(vl)))/sums;
 FiltSig2(i) = sum(sum(gaus2(vl)))/sums;
 FiltSig4(i) = sum(sum(gaus4(vl)))/sums;
 FiltSig8(i) = sum(sum(gaus8(vl)))/sums;
 %[r0 c0]=size(vl);
 n1 = gaus11(vl);
 non0 = nonzeros(n1);
 StdFiltSig1(1,i)=std(non0,1);
 n1 = gaus22(vl);
 non0 = nonzeros(n1);
 StdFiltSig2(1,i)=std(non0,1);
 n1 = gaus44(vl);
 non0 = nonzeros(n1);
 StdFiltSig4(1,i)=std(non0,1);
 n1 = gaus88(vl);
 non0 = nonzeros(n1);
 StdFiltSig8(1,i)=std(non0,1);
 MaxCorr(1,i)=sum(sum(Maxrr(vl)))/sums;
 MinCorr(1,i)=sum(sum(Minrr(vl)))/sums;
 AvgCorr(1,i)=sum(sum(Avgrr(vl)))/sums;     
end

%%%%%%%%%%%%%%%%%%%%%%% Feature 15,16,17,18 (diff)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsv = rgb2hsv(input1);
Hueinput=hsv(:,:,1);
Redinput= input1(:, :, 1);
Greeninput= input1(:, :, 2);
Blueinput= input1(:, :, 3);

mm_green=zeros(1,co);
diffr=zeros(1,co);
diffg=zeros(1,co);
diffb=zeros(1,co);
diffh=zeros(1,co);
for i=1:co
    num=stats(i).FilledImage;
    [r c]=size(num);
    s=sum(sum(num));
    mm_green(1,i)=i_green(1,i)/s;
    
    temp2=stats(i).SubarrayIdx;
    temp3=temp2{1,1};
    temp4=temp2{1,2};
    ceny=round(r/2);
    cenx=round(c/2);
    r1=temp3(1,ceny)-6;
    rm=temp3(1,ceny)+6;
    c1=temp4(1,cenx)-6;
    cm=temp4(1,cenx)+6;
    if c1<=0
        c1=1;
    end
    if r1<=0
        r1=1;
    end
    rrr=0;
    gg=0;
    bb=0;
    hh=0;
    if rm>R
        rm=R;
    end
    if cm>C
        cm=C;
    end
    for p=r1:rm
        for q=c1:cm
           rrr=rrr+Redinput(p,q);
           gg=gg+Greeninput(p,q);
           bb=bb+Blueinput(p,q);
           hh=hh+Hueinput(p,q);
        end
    end
    meanr=rrr/(r1*c1);
    meang=gg/(r1*c1);
    meanb=bb/(r1*c1);
    meanh=hh/(r1*c1);
    diffr(1,i)=meanr-mm_green(1,i);
    diffg(1,i)=meang-mm_green(1,i);
    diffb(1,i)=meanb-mm_green(1,i);
    diffh(1,i)=meanh-mm_green(1,i);
end
%%%%%%%%%%%%%%%%%%%%% Feature 30 (MajorAxis)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 MajorAxis = [stats.MajorAxisLength];
 
%%%%%%%%%%%%%%%%%%%%% Feature 31 (MinorAxis)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 MinorAxis = [stats.MinorAxisLength];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('FeatureFile.mat');