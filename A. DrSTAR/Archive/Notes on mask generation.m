%this one icludes some puncta on dark background
%BW = imbinarize(AA{k, 1}, 'adaptive', 'Sensitivity', 0.6); figure; imshow(BW);
BW_adapt = imbinarize(AA{k, 1}, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity', 0.6); figure; imshow(BW); title('BW adaptive')

%to exlude objects bigger than (remove brigtPM then repalce with sobel?)
BW_size_exlud = bwpropfilt(BW_adapt,'Area',[0 100]); figure; imshow(BW_size_exlud); title('BW adaptive size exluded 100pix')
BW_dial = imdilate(BW_size_exlud, ones(3,3)); figure; imshow(BW_dial); title('BW adaptive size exluded 100pix dialted')

%I am aftid this sensitivity will be an isue...

%good for finding brigts spots on bright background
BW_Sobel = edge(AA{k, 1}, 'sobel'); figure;imshow(BW_Sobel); title('BW sobel')
BW_Sobel_dil = imdilate(BW_Sobel, ones(3,3));
BW_Sobel_filled = imfill(BW_Sobel_dil,'holes'); figure; imshow(BW_Sobel_filled); title('BW sobel dilated and filled')

BW_both = BW_dial + BW_Sobel_filled;
BW_both(BW_both == 2) = 1; imshow(BW_both); title('Combined mask')

picresized = imresize(AA{k, 1},5);

f = 11;
sigma = 10;
I488 = imgaussfilt(AA{k, 1}, sigma);
BW_Harris_xy = detectHarrisFeatures(I488);
figure; imshow(imadjust(AA{k, 1})); hold on; plot(BW_Harris_xy); title('Detection'); hold off;


%this one icludes some puncta on dark background
%BW = imbinarize(AA{k, 1}, 'adaptive', 'Sensitivity', 0.6); figure; imshow(BW);
BW_adapt = imbinarize(AA{k, 1}, 'adaptive'); figure; imshow(BW); title('BW adaptive')

%to exlude objects bigger than (remove brigtPM then repalce with sobel?)
BW_size_exlud = bwpropfilt(BW_adapt,'Area',[0 450]); figure; imshow(BW_size_exlud); title('BW adaptive size exluded 100pix')
 BW_dial = imdilate(BW_size_exlud, ones(1,1)); figure; imshow(BW_dial); title('BW adaptive size exluded 100pix dialted')
 BW_dial = imerode(BW_size_exlud, ones(1,1));
%BW_dial = BW_size_exlud;

%good for finding brigts spots on bright background
BW_Sobel = edge(AA{k, 1}, 'sobel'); figure;imshow(BW_Sobel); title('BW sobel')
BW_Sobel_dil = imdilate(BW_Sobel, ones(3,3));
BW_Sobel_filled = imfill(BW_Sobel_dil,'holes'); figure; imshow(BW_Sobel_filled); title('BW sobel dilated and filled')

clims = [0 800];
p=FastPeakFind(AA{k, 1}, 0, ones(3,3));
figure; imagesc(AA{k, 1}, clims); axis equal; hold on
plot(p(1:2:end),p(2:2:end),'r+')

[imageSizeY,imageSizeX] = size(AA{k, 1});

[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
circlePixels = false(imageSizeY,imageSizeX);
final_mask = false(imageSizeY,imageSizeX);

for ii = 1:2:length(p)
    centerX = p(ii);
    centerY = p(ii+1);
    radius = 4;
    circlePixels = (rowsInImage - centerY).^2 ...
        + (columnsInImage - centerX).^2 <= radius.^2;
    final_mask = final_mask + circlePixels;
end
imshow(final_mask) ;

BW_both = BW_dial + final_mask + BW_Sobel_filled;
BW_both(BW_both > 1) = 1; %imshow(BW_both); title('Combined mask')
%BW_both = bwpropfilt(logical(BW_both),'Area',[0 900]);

C = imfuse(AA{k, 1}*150,BW_both,'falsecolor','Scaling','independent','ColorChannels',[2 1 2]);
figure; imshow(C)

imwrite(BW_both,'myImage.tiff','TIFF');