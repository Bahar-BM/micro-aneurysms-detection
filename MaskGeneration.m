function [mask]=MaskGeneration(image,thresh)
%--------------------------------------------------------------------------
% generating mask of retinal image
%--------------------------------------------------------------------------
% Input:
%   image: Input  image
%   thresh: threshold
%          
% Output:
%   mask: mask of retinal image
%--------------------------------------------------------------------------
% Programmer: Elaheh Imani 1392/12/19
%--------------------------------------------------------------------------
    mask=image(:,:,2)>thresh;
    se=strel('disk',5);
    mask=imopen(mask,se);
    mask=imclose(mask,se);
    clear 'se';
end