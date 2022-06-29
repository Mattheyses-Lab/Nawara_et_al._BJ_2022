clear;close all;
im=imread('rice.png');
fim=mat2gray(im);
lnfim=localnormalize(fim,4,4);
lnfim=mat2gray(lnfim);
imshow(fim);
figure,imshow(lnfim);

%%%% that might be th way to go

% Good for finding dim puncta next to brigth plaques   
close all
    [Gmag,~] = imgradient(Pseudo_Cell_488_001);
    %imshowpair(Gmag, Gdir, 'montage'); imagesc(Gmag);  axis equal
    %G = Gmag > 200; 
    figure; imagesc(G)
    %BW_Gmag = imdilate(G, ones(1,1));
    %BW_Gmag_filled = imfill(G,'holes'); %figure; imshow(BW_Gmag_filled); 
    
 

%%%
cell = BB; 
fim=mat2gray(double(cell));
lnfim=localnormalize(fim,4,4);
lnfim=mat2gray(lnfim);
%imshow(fim);
%figure;imshowpair(BB, lnfim);
figure;imshow(lnfim);


BW_cell = imbinarize(lnfim); %figure; imshow(BW_cell)
%%%%%%%%% Add dialtion at the end

[Gmag,~] = imgradient(cell*1000);
G = Gmag > 100;

BW_G = imfill(~G,'holes'); figure; imshow(BW_G)
BW_G_exlud = bwpropfilt(BW_G,'Area',[10 10*10^7]); figure; imshow(BW_G_exlud)
BW_G_dial = imdilate(BW_G_exlud, ones(10,10)); figure; imshow(BW_G_dial)

BW_both = BW_cell .* BW_G_dial;

C = imfuse(cell*150,BW_both,'falsecolor','Scaling','independent','ColorChannels',[2 1 2]);
figure; imshow(C); title('488 mask on 488 image');

%%%




BW_adapt = imbinarize(lnfim, 'adaptive', 'Sensitivity', 0.4); figure; imshow(BW_adapt)
% BW_dial = imdilate(BW_adapt, ones(3,3)); figure; imshow(BW_dial)
% BW_eroded = imerode(BW_dial, ones(3,3)); figure; imshow(BW_eroded)
BW_holes = imfill(BW_adapt,'holes'); figure; imshow(BW_holes)
BW_size_exlud = bwpropfilt(BW_holes,'Area',[3 400]); figure; imshow(BW_size_exlud)


BW_holes = imfill(BW_cell,'holes'); figure; imshow(BW_holes)
BW_dial = imdilate(BW_cell, ones(3,3)); figure; imshow(BW_dial)

BW_size_exlud = bwpropfilt(BW_cell,'Area',[3 400]); figure; imshow(BW_size_exlud)

BW_dial = imdilate(BW_cell, ones(3,3)); figure; imshow(BW_dial)
BW_eroded = imerode(BW_cell, ones(3,3)); figure; imshow(BW_eroded)


BW_adapt = imbinarize(Gmag, 'adaptive', 'Sensitivity', 0.1); figure; imshow(BW_adapt) 


    
    % Good for finding dim puncta next to brigth plaques     
    [Gmag,~] = imgradient(Pseudo_Cell_488_001);
    %imshowpair(Gmag, Gdir, 'montage'); imagesc(Gmag);  axis equal
    %G = Gmag > 450; 
    figure; imagesc(Gmag)
    %BW_Gmag = imdilate(G, ones(1,1));
    BW_Gmag_filled = imfill(G,'holes'); %figure; imshow(BW_Gmag_filled); 
    
    zzz = imgaussfilt(bbb,2) .* BW_cor;
    bb(isnan(bb)) = 0;
    zzzz = zzz + bb;
    
    figure; imagesc(zzzz)
    
    
    bb+
    

    
    
    
    
    
    