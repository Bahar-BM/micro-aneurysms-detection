function bloodVessels = VesselExtract(inImg, threshold)
[r c]=size(inImg);
bloodVessels=zeros(r,c);
%Kirsch's Templates
h1=[5 -3 -3;
    5  0 -3;
    5 -3 -3]/15;
h2=[-3 -3 5;
    -3  0 5;
    -3 -3 5]/15;
h3=[-3 -3 -3;
     5  0 -3;
     5  5 -3]/15;
h4=[-3  5  5;
    -3  0  5;
    -3 -3 -3]/15;
h5=[-3 -3 -3;
    -3  0 -3;
     5  5  5]/15;
h6=[ 5  5  5;
    -3  0 -3;
    -3 -3 -3]/15;
h7=[-3 -3 -3;
    -3  0  5;
    -3  5  5]/15;
h8=[ 5  5 -3;
     5  0 -3;
    -3 -3 -3]/15;

%Spatial Filtering by Kirsch's Templates
t1=filter2(h1,inImg);
t2=filter2(h2,inImg);
t3=filter2(h3,inImg);
t4=filter2(h4,inImg);
t5=filter2(h5,inImg);
t6=filter2(h6,inImg);
t7=filter2(h7,inImg);
t8=filter2(h8,inImg);

%s=size(inImg);

%temp=zeros(1,8);
%
bV = max(t1,t2); 
bV = max(bV,t3); 
bV = max(bV,t4); 
bV = max(bV,t5); 
bV = max(bV,t6); 
bV = max(bV,t7); 
bV = max(bV,t8);


for i=1:r
    for j=1:c
        if(bV(i,j)>threshold)
            bloodVessels(i,j)=bV(i,j);
        end
    end
end
