%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%Vessle Map(1)%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clc;clear;close all 
warning off
global BlkMsk TetaStp;

%Parameters
% 	n=13;								%######### Subimage Size
% 	stp=3;								%######### Overlap Control 
    n=13;								%######### Subimage Size
    LinVldThr = 3.3;					%######### Line Validation
	stp=3;								%######### Overlap Control 
	TetaStp = 4;						%########### Teta Steps in Rad. Trans.
% 	LinVldThr = 3.3;					%######### Line Validation
%End Parameters

% Circular Mask
BlkMsk = ones(n,n);
nh = fix(n/2)+1;
rad2 = nh^2;
for r = 1:n
	for c = 1:n
		if (r-nh)^2+(c-nh)^2 > rad2
			BlkMsk(r,c) = 0;
		end
	end
end


for ii=0:49
f=imread(strcat('images\image',int2str(ii),'_test.jpg'));
%f=imresize(f,[480 640]);
[R,C,D]=size(f);
msk=zeros(R, C, 'uint8');
vesselWidth=zeros(R, C, n, 'uint8');
for r=1:fix(n/stp):R-n+1
	for c=1:fix(n/stp):C-n+1
		fl=f(r:r+n-1,c:c+n-1,:);
		[mx LineStrt LineEnd LineAngle LineVld]=LocRadFun1(fl);
		if LineVld>(LinVldThr/n)
				ml=zeros(n,n); mz=zeros(n,n);
				ml(:,LineStrt:LineEnd)=1;
				ml=[mz ml mz; mz ml mz; mz ml mz];
				ml=imrotate(ml,LineAngle,'crop');
				mskl=uint8(ml(n+1:2*n,n+1:2*n));
				f2=fl(:,:,2);	f2=f2.*mskl;		t1=mean(f2(f2~=0));
				f2=fl(:,:,2);	f2=f2.*(1-mskl);	t2=mean(f2(f2~=0));
				f2=fl(:,:,2);	f2=uint8(im2bw(im2double(f2), (t1+t2)/2/255));
				mskl=mskl.*(1-f2);
				msk(r:r+n-1,c:c+n-1)=msk(r:r+n-1,c:c+n-1)|mskl;
		end
	end
end

		
g1=255-f(:,:,2);
imwrite(g1, strcat('graytest\graytest',int2str(ii),'.jpg'));
	%figure(2); imshow(g1);

	%figure(3); imshow(msk*255);
imwrite(msk*255, strcat('masktest\masktest',int2str(ii),'.jpg'));
msk=imclose(msk,strel('disk', 2));
% 			figure(4); imshow(msk*255);
imwrite(msk*255,strcat('modmasktest\modmasktest',int2str(ii),'.jpg'));
	
em=edge(msk);	em=1-em;
imwrite(em,strcat('edgtest\edgtest',int2str(ii),'.jpg'));
% 	figure(5); imshow(em);
	
	%f(:,:,1)=f(:,:,1).*em;
f(:,:,1)=f(:,:,1).*uint8(em);	f(:,:,3)=f(:,:,3).*uint8(em);
% 	figure(6); imshow(f);
imwrite(f,strcat('resulttest\resulttest',int2str(ii),'.jpg'));
end