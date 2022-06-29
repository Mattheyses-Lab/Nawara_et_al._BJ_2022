function PM_only_pic = gap_filling_v2(AA)

%Filll gaps based on mean (chage to zero to fill based on meddian)
mean_filt = 1;
smothing = 1; 

close all


% sigma = 5;
% AA =
% bb = imgaussfilt(AA, sigma);
% % imagesc(bb)
% % imshowpair(AA,bb,'montage')
% cell_mat2tiff('C:\Users\tnawara\Desktop\test.tif', bb)
%
% bb = edge(AA, 'sobel');
%
% se = offsetstrel('ball',3,3);
% bb = imdilate(AA,se);
% figure(1), imagesc(AA), figure(2), imagesc(bb)
%
% level = graythresh(AA)
% BW = imbinarize(AA,level);
% imshowpair(AA,BW,'montage').

%mask to detect puncta

PM_only_pic = cell(size(AA));
BW_triple = cell(size(AA));
for ii = 1:length(AA)
    PM_only_pic{ii,1} = zeros(size(AA{1,1}));
    BW_triple{ii,1} = zeros(size(AA{1,1}));
end

parfor k = 1:length(AA) %Generate blured 488 to match 647 intrploation
    
    %     BW = imbinarize(AA{k, 1}, 'adaptive', 'ForegroundPolarity', 'bright'); figure; imshow(BW);
    %     BW_dial = imdilate(BW, ones(4,4));
    %     BW_eroded = imerode(BW_dial, ones(2,2));
    %     mask = double(~BW_eroded);
    
    
    %this one icludes some puncta on dark background
    %BW = imbinarize(AA{k, 1}, 'adaptive', 'Sensitivity', 0.6); figure; imshow(BW);
    BW_adapt = imbinarize(AA{k, 1}, 'adaptive', 'Sensitivity', 0.1); %figure; imshow(BW_adapt); title('BW adaptive')
    
    %to exlude objects bigger than (remove brigtPM then repalce with sobel?)
    BW_size_exlud = bwpropfilt(BW_adapt,'Area',[10 900]); %figure; imshow(BW_size_exlud); title('BW adaptive size exluded 100pix')
    BW_dial = imdilate(BW_size_exlud, ones(4,4)); %figure; imshow(BW_dial); title('BW adaptive size exluded 100pix dialted')
    BW_dial = imerode(BW_size_exlud, ones(1,1));
    %BW_dial = BW_size_exlud;
    %figure; imshow(BW_dial)
    
    %good for finding brigts spots on bright background
    BW_Sobel = edge(AA{k, 1}, 'canny'); %figure;imshow(BW_Sobel); title('BW sobel') %%% used to be sobel method
    BW_Sobel_dil = imdilate(BW_Sobel, ones(3,3));
    BW_Sobel_filled = imfill(BW_Sobel_dil,'holes'); %figure; imshow(BW_Sobel_filled); title('BW sobel dilated and filled')
    
    clims = [0 800];
    p = FastPeakFind(AA{k, 1}, 0, ones(3,3));
    %     figure; imagesc(AA{k, 1}, clims); axis equal; hold on
    %     plot(p(1:2:end),p(2:2:end),'r+')
    
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
    %imshow(final_mask) ;
    
    BW_triple_x = BW_dial + final_mask + BW_Sobel_filled;
    BW_triple_x(BW_triple_x > 1) = 1;
    BW_triple{k,1} = BW_triple_x;
    
    
    %imshow(BW_triple) ;
    %      C = imfuse(AA{k, 1}*150,BW_triple,'falsecolor','Scaling','independent','ColorChannels',[2 1 2]);
    %      figure; imshow(C)

end

if smothing  == 1
    matrix_3d = zeros(size(BW_triple{1, 1},1), size(BW_triple{1, 1},2), size(BW_triple,1));
    matrix_3d_smooth = matrix_3d;
    averaginng_window = 10;
    
    for k = 1:size(BW_triple, 1)
        matrix_3d(:,:,k) = BW_triple{k, 1};
    end
    
    padding_for_matrix_3d = nan(size(BW_triple{1, 1},1), size(BW_triple{1, 1},2), averaginng_window); %in agremwnt wityh movmean padding of a window
    
    matrix_3d = cat(3,matrix_3d, padding_for_matrix_3d);
    
    for k = 1:size(matrix_3d, 3)
        if k > size(matrix_3d_smooth,3)
            break
        end
        matrix_3d_smooth(:,:,k) = mean(matrix_3d(:,:,k:k+(averaginng_window-1)),3,'omitnan');
    end
    
    for k = 1:size(matrix_3d_smooth, 3)
        matrix_3d_smooth(:,:,k) = ceil(matrix_3d_smooth(:,:,k));
    end
    
    for k = 1:size(matrix_3d_smooth,3)
        BW_triple{k, 1} = matrix_3d_smooth(:,:,k);
    end
end
    

parfor k = 1:length(AA)
    
    mask = double(~BW_triple{k, 1});
    
    mask(mask == 0) = NaN;
    bb = double(AA{k, 1}) .* mask;
    
    frame = 5;
    frame4median = frame - 1;
    bb_frame = double(padarray(bb,[frame frame],0,'both'));
    bb_frame(1:frame,1:size(bb_frame,2)) = NaN;
    bb_frame(size(bb_frame,1)-(frame-1):size(bb_frame,1) ,1:size(bb_frame,2)) = NaN;
    bb_frame(1:size(bb_frame,1),1:frame) = NaN;
    bb_frame(1:size(bb_frame,1),size(bb_frame,2)-(frame-1):size(bb_frame,2)) = NaN;
    bb_frame_corrected = bb_frame;
    
    while isnan(sum(bb_frame_corrected((frame+1):end-(frame),(frame+1):end-(frame)), 'all'))
        for col = frame+1:size(bb_frame,2)-(frame)
            for row = frame+1:size(bb_frame,1)-(frame)
                if isnan(bb_frame(row,col))
                    matrix = bb_frame(row-frame4median:row+frame4median, col-frame4median:col+frame4median);
                    if mean_filt == 1
                        bb_frame_corrected(row, col) = round(mean(matrix,'all', 'omitnan'));
                    else
                        bb_frame_corrected(row, col) = round(median(matrix,'all', 'omitnan'));
                    end
                end
            end
        end
        bb_frame = bb_frame_corrected;
        %imagesc(bb_frame_corrected)
    end
    
    bbb = bb_frame((frame+1):end-(frame),(frame+1):end-(frame));
    
    %bbb_smooth = imgaussfilt(bbb, 2); %imagesc(bbb_smooth);
    PM_only_pic{k, 1} = imgaussfilt(bbb, 2); %imagesc(bbb_smooth);
end

%%%%






% average_window = 2; % this controls the size of the moving average filter - i.e. the filter will be [average_window x average_window]
% average_filter = fspecial('average',[average_window average_window]);
% matrix_3d_smooth = imfilter(matrix_3d,average_filter);



% clims = [0 800];
% subplot(3,1,1); imagesc(AA{1, 1}, clims); title('Oryginal'); axis equal; axis([0 size(AA,2) 0 size(AA,1)]);
% subplot(3,1,2); imagesc(double(BW_eroded)); title('Puncta mask');  axis equal; axis([0 size(AA,2) 0 size(AA,1)]);
% subplot(3,1,3); imagesc(PM_only_pic{1, 1},clims); title('Puncta removed'); axis equal;  axis([0 size(AA,2) 0 size(AA,1)]);

% subplot(2,1,1); imagesc(AA, clims); title('Oryginal'); axis equal; axis([0 size(AA,2) 0 size(AA,1)]);
% subplot(2,1,2); imagesc(bbb_smooth,clims); title('Puncta removed'); axis equal;  axis([0 size(AA,2) 0 size(AA,1)]);
%
%
cell_mat2tiff('C:\Users\Nikon\Desktop\Optmization of gapfilt\test_full_smoothing_of_mask_ceil_sens0.1.tif', PM_only_pic)

%%%% itrs still notr perfect work on the bliking big poimts next maybe
%%%% anther mask just for big objects

end

