function [Stage1,Maxrr,Minrr,Avgrr]=CandidateMicroExtract(counter,NormInput,bloodVessels,input)
[R C]=size(input);

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
 
r1=CorrAB(NormInput,gussian1);
r2=CorrAB(NormInput,gussian2);
r3=CorrAB(NormInput,gussian3);
r4=CorrAB(NormInput,gussian4);
r5=CorrAB(NormInput,gussian5);   

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
%imwrite(Maxrr,strcat('Maxrr\maxrr',int2str(counter),'.jpg'));

%threshold1
output=(Maxrr>0.5);
%imwrite(output,strcat('Threshold\Thr',int2str(counter),'.jpg'));
%Removing any candidates on the vessels
Stage1=output-bloodVessels;
for r=1:R
    for c=1:C
        if Stage1(r,c)<0
            Stage1(r,c)=0;
        end
    end
end

%imwrite(Stage1,strcat('Stage1\Stage',int2str(counter),'.jpg'));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Region Growing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask1 = false(size(input));
[L1,num1]=bwlabel(Stage1);
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
Stage1=mask1-bloodVessels;
imwrite(Stage1,strcat('RegionGrowing\MsImani\Thresh05\RegionGrowing',int2str(counter),'.jpg'));