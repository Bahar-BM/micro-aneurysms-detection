function [mx LStrt LEnd LineAngle LineVld]=locradfun1(f)
global BlkMsk TetaStp;
%figure(1); imshow(f(:,:,2));
Level=0.7;                          %###########Line Width Control

[R C D]=size(f);
f=double(f);
g1=255-f(:,:,2);
g1=g1/255;	g1=imadjust(g1);
immean = mean2(g1);
g1=g1.*BlkMsk;
%figure(9); imshow(g1);

theta = 0:TetaStp:180;
[Rd,xp] = radon(g1,theta);

h=fspecial('average',3);
Rd=imfilter(Rd,h);

[mv mi]=max(Rd);
[mx mc]=max(mv);
mr=mi(mc);

LineAngle=mc*TetaStp;

if LineAngle > 90
	OrtAngle=fix((LineAngle-90)/TetaStp);
else
	OrtAngle=fix((LineAngle+90)/TetaStp);
end

maxcul=Rd(:,mc);
maxcul=maxcul-R*immean;		maxcul=maxcul.*((sign(maxcul)+1)/2);
%figure(2); plot(maxcul);
[mx LinePos]=max(maxcul);

if mx == 0
	LineVld=0;
	LStrt=0;
	LEnd=0;
else
	%LinePos=mr;
	LineVld=mx/R;
	[RR RC]=size(Rd);

	LStrt=LinePos;
	LEnd=LinePos;
	for i=LinePos:-1:1
	    LStrt=i;
	    if maxcul(i)<Level*mx
	        break
		end
	end
	for i=LinePos:1:RC
	    LEnd=i;
	    if maxcul(i)<Level*mx
	        break
		end
	end

	LStrt=LStrt-fix(RR/2)+fix(R/2);
	LEnd=LEnd-fix(RR/2)+fix(R/2);
	%[LStrt LEnd]
	if LStrt<1
	    LStrt=1;
	end
	if LEnd>R
	    LEnd=R;
	end
	%[LStrt LEnd]
end
