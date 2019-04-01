clc;
clear all;
close all;
path=input('Enter path of the image :','s');
% rgb=imread(path);
% grimg=rgb2gray(rgb);
% siz=size(grimg);
% row=siz(1);
% col=siz(2);
% fin=zeros(size(grimg));
% dup=fin;
% imshow(grimg);
% impixelinfo;
% grimg=double(grimg);
% for a=1:2
% stx=input('Enter Row no :');
% sty=input('Enter column no :');
% th=input('Enter thershold :');
% count=0;
% while length(stx)~=0
% count=count+1;
% display(['Pixels added :',num2str(count)]);
% ls=length(stx);
% sr=stx(ls);
% sc=sty(ls);
% svalue=grimg(sr,sc);
% stx=stx(1:ls-1);
% sty=sty(1:ls-1);
% if ls==1
% stx=[];
% sty=[];
% end
% for i=sr-1:sr+1
%     if i==0 || i==row+1
%         continue
%     end
%     for j=sc-1:sc+1;
%         if j==0 || j==col+1
%             continue
%         end
%         dif=grimg(i,j)- svalue;
%         if abs(dif)<=th
%             if dup(i,j)~=1
%             stx=[stx i];
%             sty=[sty j];
%             fin(i,j)=grimg(i,j);
%             dup(i,j)=1;
%             end
%         end
%     end
% end
% end
% figure;
% subplot(121);
% imshow(uint8(grimg));title('Original image');
% impixelinfo;
% subplot(122);
% imshow(uint8(fin));hold on;
% title('Grown image');
% impixelinfo;
% end
% imhist(uint8(fin));
