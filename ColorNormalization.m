function [img]=ColorNormalization(image,mask,window)
%--------------------------------------------------------------------------
% color normalization of retinal images
%--------------------------------------------------------------------------
% Input:
%   window: window size
%          
% Output:
%   img: normalized retinal image
%--------------------------------------------------------------------------
% Programmer: Elaheh Imani 1393/2/10
%--------------------------------------------------------------------------
%     global image;
%     global mask;
    b=load('color_normalization');
    Mu=b.Mu;
    STD=b.STD;
    clear b
    min_value=0;max_value=255;
    R=image(:,:,1);
    G=image(:,:,2);
    B=image(:,:,3);
    window_size=round(size(image,1)/window);
    R_filtered=double(medfilt2(R,[window_size window_size]));
    G_filtered=double(medfilt2(G,[window_size window_size]));
    B_filtered=double(medfilt2(B,[window_size window_size]));
    
    R_filtered =(double(R)-R_filtered);
    G_filtered = (double(G)-G_filtered);
    B_filtered = (double(B)-B_filtered);
    clear R G B
    m_r=min(R_filtered(mask(:)));
    m_g=min(G_filtered(mask(:)));
    m_b=min(B_filtered(mask(:)));
    
    R_filtered=min_value+((R_filtered-m_r)./(max(R_filtered(mask(:)))-m_r))*max_value;
    G_filtered=min_value+((G_filtered-m_g)./(max(G_filtered(mask(:)))-m_g))*max_value;
    B_filtered=min_value+((B_filtered-m_b)./(max(B_filtered(mask(:)))-m_b))*max_value;
    clear m_r m_g m_b min_value max_value
    Mu_r=mean(R_filtered(mask(:)));
    Mu_g=mean(G_filtered(mask(:)));
    Mu_b=mean(B_filtered(mask(:)));
    
    STD_r=std(R_filtered(mask(:)));
    STD_g=std(G_filtered(mask(:)));
    STD_b=std(B_filtered(mask(:)));
    
    R_filtered=(R_filtered-Mu_r)./STD_r;
    G_filtered=(G_filtered-Mu_g)./STD_g;
    B_filtered=(B_filtered-Mu_b)./STD_b;
    clear  STD_r STD_g STD_b Mu_r Mu_g Mu_b
    new_r=(R_filtered*STD(1))+Mu(1);
    new_g=(G_filtered*STD(2))+Mu(2);
    new_b=(B_filtered*STD(3))+Mu(3);
    clear R_filtered G_filtered B_filtered STD Mu
    new_r(~mask)=0;
    new_g(~mask)=0;
    new_b(~mask)=0;
    img=zeros(size(image),'uint8');
    img(:,:,1)=new_r;
    img(:,:,2)=new_g;
    img(:,:,3)=new_b;
    clear new_r new_g new_b
end