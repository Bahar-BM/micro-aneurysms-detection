for jj=44:44
image=imread(strcat('training\image',int2str(jj),'_training.jpg'));
thresh=0;
window=20;
mask=MaskGeneration(image,thresh);
img=ColorNormalization(image,mask,window);
imwrite(img,strcat('normalization\normalization',int2str(jj),'.jpg'));
%imshow(img);
end