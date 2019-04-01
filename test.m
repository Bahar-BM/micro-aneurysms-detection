jj=0;
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

figure;
imshow(Maxrr);
imwrite(Maxrr, strcat('Out1\Out1test',int2str(jj),'.jpg'));

%threshold1
output=zeros(R,C);
for r=1:R
    for c=1:C
        if Maxrr(r,c)>0.4
            output(r,c)=1;
        end
    end
end
figure;
imshow(output);
imwrite(output, strcat('Out2\Outtest2',int2str(jj),'.jpg'));

%Removing any candidates on the vessels
Stage1=output-bloodVessels;
figure;
imshow(Stage1);
imwrite(Stage1, strcat('Out3\Out3test',int2str(jj),'.jpg'));
% 
% %threshold2
% output4=zeros(R,C);
% for r=1:R
%     for c=1:C
%         if Stage1(r,c)>0.7
%             output4(r,c)=1;
%         end
%     end
% end
% figure;
% imshow(output4);
% imwrite(output4, strcat('Out4\Out4test',int2str(jj),'.jpg'));
% Backgrond = medfilt2(input, [25 25]);