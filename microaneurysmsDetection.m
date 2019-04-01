% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%Vessle Map(1)%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clc;clear;close all 
% warning off
% global BlkMsk TetaStp;
% 
% %Parameters
% % 	n=13;								%######### Subimage Size
% % 	stp=3;								%######### Overlap Control 
%     n=13;								%######### Subimage Size
%     LinVldThr = 3.3;					%######### Line Validation
% 	stp=3;								%######### Overlap Control 
% 	TetaStp = 4;						%########### Teta Steps in Rad. Trans.
% % 	LinVldThr = 3.3;					%######### Line Validation
% %End Parameters
% 
% % Circular Mask
% BlkMsk = ones(n,n);
% nh = fix(n/2)+1;
% rad2 = nh^2;
% for r = 1:n
% 	for c = 1:n
% 		if (r-nh)^2+(c-nh)^2 > rad2
% 			BlkMsk(r,c) = 0;
% 		end
% 	end
% end
% 
% 
% for ii=0:0
% f=imread(strcat('images\image',int2str(ii),'_test.jpg'));
% %f=imresize(f,[480 640]);
% [R,C,D]=size(f);
% msk=zeros(R, C, 'uint8');
% vesselWidth=zeros(R, C, n, 'uint8');
% for r=1:fix(n/stp):R-n+1
% 	for c=1:fix(n/stp):C-n+1
% 		fl=f(r:r+n-1,c:c+n-1,:);
% 		[mx LineStrt LineEnd LineAngle LineVld]=LocRadFun1(fl);
% 		if LineVld>(LinVldThr/n)
% 				ml=zeros(n,n); mz=zeros(n,n);
% 				ml(:,LineStrt:LineEnd)=1;
% 				ml=[mz ml mz; mz ml mz; mz ml mz];
% 				ml=imrotate(ml,LineAngle,'crop');
% 				mskl=uint8(ml(n+1:2*n,n+1:2*n));
% 				f2=fl(:,:,2);	f2=f2.*mskl;		t1=mean(f2(f2~=0));
% 				f2=fl(:,:,2);	f2=f2.*(1-mskl);	t2=mean(f2(f2~=0));
% 				f2=fl(:,:,2);	f2=uint8(im2bw(im2double(f2), (t1+t2)/2/255));
% 				mskl=mskl.*(1-f2);
% 				msk(r:r+n-1,c:c+n-1)=msk(r:r+n-1,c:c+n-1)|mskl;
% 		end
% 	end
% end
% 
% 		
% g1=255-f(:,:,2);
% imwrite(g1, strcat('graytest\graytest',int2str(ii),'.jpg'));
% 	%figure(2); imshow(g1);
% 
% 	%figure(3); imshow(msk*255);
% imwrite(msk*255, strcat('masktest\masktest',int2str(ii),'.jpg'));
% msk=imclose(msk,strel('disk', 2));
% % 			figure(4); imshow(msk*255);
% imwrite(msk*255,strcat('modmasktest\modmasktest',int2str(ii),'.jpg'));
% 	
% em=edge(msk);	em=1-em;
% imwrite(em,strcat('edgtest\edgtest',int2str(ii),'.jpg'));
% % 	figure(5); imshow(em);
% 	
% 	%f(:,:,1)=f(:,:,1).*em;
% f(:,:,1)=f(:,:,1).*uint8(em);	f(:,:,3)=f(:,:,3).*uint8(em);
% % 	figure(6); imshow(f);
% imwrite(f,strcat('resulttest\resulttest',int2str(ii),'.jpg'));
% end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%Vessle Map(2)%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % clc;
% % % input=imread('image13_training.jpg');
% % % %input= input(:, :, 2);
% % % figure;
% % % imshow(input);
% % % [R C]=size(input);
% % % 
% % % 
% % % %Extract Blood Vessels
% % % Threshold = 7;
% % % bloodVessels = VesselExtract(input, Threshold);
% % % 
% % % %Output Blood Vessels image
% % % figure;
% % % imshow(bloodVessels);
% % % 
% % % se = strel('disk',1);
% % % dilate = imdilate(bloodVessels,se);
% % % figure, imshow(dilate)
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%% Phase1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
countertable=1;
xDoc = xmlread('annotations-consensus-ma-only.xml');
strc  = xml2struct( xDoc );

%MicroIndx=zeros(1,co);
countindx=1;

for jj=1:1d=strc.set.annotations_dash_per_dash_image(jj+1);
d=d{1,1};
if isfield(d,'annotation')

input1=double(imread(strcat('training\image',int2str(jj),'_training.jpg')));
input= input1(:, :, 2);
[R C]=size(input);
bloodVessels=imread(strcat('modmask\modmask',int2str(jj),'.jpg'));
bloodVessels=mat2gray(bloodVessels);

gussian1=zeros(9,9);
gussian2=zeros(9,9);
gussian3=zeros(9,9);
gussian4=zeros(11,11);
gussian5=zeros(11,11);

sigma1=1.1;
sigma2=1.2;
sigma3=1.3;
sigma4=1.4;
sigma5=1.5;

for i=1:9
  for j=1:9
      gussian1(i,j)=(1/(2*pi*(sigma1^2)))*exp(-((i-5)^2+(j-5)^2)/(2*(sigma1^2)));
  end
end

for i=1:9
  for j=1:9
      gussian2(i,j)=(1/(2*pi*(sigma2^2)))*exp(-((i-5)^2+(j-5)^2)/(2*(sigma2^2)));
  end
end

for i=1:9
  for j=1:9
      gussian3(i,j)=(1/(2*pi*(sigma3^2)))*exp(-((i-5)^2+(j-5)^2)/(2*(sigma3^2)));
  end
end
%%%
for i=1:11
  for j=1:11
      gussian4(i,j)=(1/(2*pi*(sigma4^2)))*exp(-((i-6)^2+(j-6)^2)/(2*(sigma4^2)));
  end
end
%%%
for i=1:11
  for j=1:11
      gussian5(i,j)=(1/(2*pi*(sigma5^2)))*exp(-((i-6)^2+(j-6)^2)/(2*(sigma5^2)));
  end
end

gussian1=imcomplement((mat2gray(gussian1)));
gussian2=imcomplement((mat2gray(gussian2)));
gussian3=imcomplement((mat2gray(gussian3)));
gussian4=imcomplement((mat2gray(gussian4)));
gussian5=imcomplement((mat2gray(gussian5)));

Maxrr=zeros(R,C);
Minrr=zeros(R,C);
Avgrr=zeros(R,C);
 
r1=CorrAB(input,gussian1);
r2=CorrAB(input,gussian2);
r3=CorrAB(input,gussian3);
r4=CorrAB(input,gussian4);
r5=CorrAB(input,gussian5);   

for i=5:R-6
  for j=5:C-6
    r=i-4;
    c=j-4;
    A=[r1(r,c),r2(r,c),r3(r,c),r4(r,c),r5(r,c)];
    Maxrr(i,j)=max(A);
    Minrr(i,j)=min(A);
    Avgrr(i,j)=mean(A);
  end
end

imwrite(Maxrr, strcat('Out1\Out1',int2str(jj),'.jpg'));

%threshold1
output=zeros(R,C);
for r=1:R
    for c=1:C
        if Maxrr(r,c)>0.7
            output(r,c)=1;
        end
    end
end

imwrite(output, strcat('Out2\Out2',int2str(jj),'.jpg'));

%Removing any candidates on the vessels
Stage1=output-bloodVessels;
for r=1:R
    for c=1:C
        if Stage1(r,c)<0
            Stage1(r,c)=0;
        end
    end
end
imwrite(Stage1, strcat('Out3\Out3',int2str(jj),'.jpg'));
output4=Stage1;
Backgrond = medfilt2(input, [25 25]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%% Region Growing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx=0.5;
mask1 = false(size(input));
[L1,num1]=bwlabel(output4);
stats1 = regionprops(L1,'all');

for i=1:num1
    vlcent=stats1(i).Centroid;
    jm=vlcent(1,1);
    im=vlcent(1,2);
    rms=2;
    [x,y]=meshgrid((im-rms):(im+rms),(jm-rms):(jm+rms));
    c_mask=(((x-im).^2+(y-jm).^2)<=rms^2);
    [r1 c1]=size(c_mask);
    for cnt1=1:r1
      for cnt2=1:c1
        if c_mask(cnt1,cnt2)==1
         x1=round(x(cnt1,cnt2));
         y1=round(y(cnt1,cnt2));
         if x1>0 && x1<R && y1>0 && y1<C
         mask1(x1,y1)=1;  
         end
        end
      end
    end
end
mask = false(size(input));
[L1,num]=bwlabel(mask1);
for i=1:num
    vl = L1==i;
    idark = min(min(input(vl)));
    [rr,cc]=find(input==idark & vl,1,'first');
    thd = idark-xx*(idark-Backgrond(rr,cc));
    region = (input<=(thd)) & vl;
    [L2,num2]=bwlabel(region,8);
    region(L2~=L2(rr,cc))=0;
    mask(region)=true;
end

%figure(3);imshow(mask);
output4=mask;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%% Feature Extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
 [r0 c0]=size(vl);
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

%%%%%%%%%%%%%%%%%%%%%%% Feature 14 (v)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%baporsaaaaaaaaaam
v=zeros(1,co);
for i=1:co
    temp=stats(i).FilledImage;
    [r c]=size(temp);
    up=temp(1,:);
    down=temp(r,:);
    k=1;
    for j=1:c
        if up(1,j)==1
            dd(1,k)=sqrt((r/2-1)^2+(c/2-j)^2);
            k=k+1;
        end
        
        if down(1,j)==1
            dd(1,k)=sqrt((r/2-r)^2+(c/2-j)^2);
            k=k+1;
        end
    end
    left=temp(:,1);
    right=temp(:,c);
    
    for j=1:r
        if left(j,1)==1
            dd(1,k)=sqrt((r/2-j)^2+(c/2-1)^2);
            k=k+1;
        end
        
        if right(j,1)==1
            dd(1,k)=sqrt((r/2-j)^2+(c/2-c)^2);
            k=k+1;
        end
    end
    dm=mean(dd);
    dd=dd-dm;
    for l=1:k-1
        dd(1,l)=dd(1,l)^2;
    end
    v(i)=sqrt(sum(dd)/(k-1));
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
 
load('minmaxfile.mat','minn','maxx');
ind1=(a>=minn(1) & a<=maxx(1)); 
ind2=(per>=minn(2) & per<=maxx(2));
ind3=(MrM>=minn(3) & MrM<=maxx(3));
ind4=(circularity>=minn(4) & circularity<=maxx(4));
ind5=(i_green>=minn(5) & i_green<=maxx(5));
ind6=(i_sc>=minn(6) & i_sc<=maxx(6));
ind7=(m_green>=minn(7) & m_green<=maxx(7));
ind8=(m_sc>=minn(8) & m_sc<=maxx(8));
ind9=(NI_green>=minn(9) & NI_green<=maxx(9));
ind10=(NI_sc>=minn(10) & NI_sc<=maxx(10));
ind11=(NM_green>=minn(11) & NM_green<=maxx(11));
ind12=(NM_sc>=minn(12) & NM_sc<=maxx(12));
ind13=(IDarkest>=minn(13) & IDarkest<=maxx(13));
ind14=(FiltSig1>=minn(14) & FiltSig1<=maxx(14));
ind15=(FiltSig2>=minn(15) & FiltSig2<=maxx(15));
ind16=(FiltSig4>=minn(16) & FiltSig4<=maxx(16));
ind17=(FiltSig8>=minn(17) & FiltSig8<=maxx(17));
ind18=(StdFiltSig1>=minn(18) & StdFiltSig1<=maxx(18));
ind19=(StdFiltSig2>=minn(19) & StdFiltSig2<=maxx(19));
ind20=(StdFiltSig4>=minn(20) & StdFiltSig4<=maxx(20));
ind21=(StdFiltSig8>=minn(21) & StdFiltSig8<=maxx(21));
ind22=(MaxCorr>=minn(22) & MaxCorr<=maxx(22));
ind23=(MinCorr>=minn(23) & MinCorr<=maxx(23));
ind24=(AvgCorr>=minn(24) & AvgCorr<=maxx(24));
ind25=(v>=minn(25) & v<=maxx(25));
ind26=(diffr>=minn(26) & diffr<=maxx(26));
ind27=(diffg>=minn(27) & diffg<=maxx(27));
ind28=(diffb>=minn(28) & diffb<=maxx(28));
ind29=(diffh>=minn(29) & diffh<=maxx(29));
ind30=(MajorAxis>=minn(30) & MajorAxis<=maxx(30));
ind31=(MinorAxis>=minn(31) & MinorAxis<=maxx(31));

indd=ind1&ind2&ind3&ind4&ind5&ind6&ind7&ind8&ind9&ind10&ind11&ind12&ind13&ind14&ind15&ind16&ind17;
ind=indd&ind18&ind19&ind20&ind21&ind22&ind23&ind24&ind25&ind26&ind27&ind28&ind29&ind30&ind31;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

finalmask = false(size(input));

for i=1:co
  if ind(i)==1
    vlcent=stats1(i).Centroid;
    jm=vlcent(1,1);
    im=vlcent(1,2);
    rms=10;
    [x,y]=meshgrid((im-rms):(im+rms),(jm-rms):(jm+rms));
    c_mask=(((x-im).^2+(y-jm).^2)<=rms^2);
    [r1 c1]=size(c_mask);
    for cnt1=1:r1
      for cnt2=1:c1
        if c_mask(cnt1,cnt2)==1
         x1=round(x(cnt1,cnt2));
         y1=round(y(cnt1,cnt2));
         if x1>0 && x1<R && y1>0 && y1<C
            finalmask(x1,y1)=1;  
         end
        end
      end
    end
 end    
end
imshow(finalmask);

em=edge(finalmask);	em=1-em;
imwrite(em,strcat('ResultEdgTest\ResultEdgTest',int2str(jj),'.jpg'));
input_img(:,:,1)=input_img(:,:,1).*uint8(em);	input_img(:,:,3)=input_img(:,:,3).*uint8(em);
imwrite(input_img,strcat('FinalResultTest\FinalResultTest',int2str(jj),'.jpg'));
end
