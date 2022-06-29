function [PM_only_pic, Cell_mask] = Signal_removal_V2(AA, lambda, Cell_mask_488)

%Filll gaps based on mean (chage to zero to fill based on meddian)
smothing = 1;
averaginng_window = 10; % make this always even as it will be devided by 2

if ~iscell(AA) %incase its only one pictre draged into matlab
    BB = AA;
    clear AA
    AA{1,1} = BB;
    clear BB
    
elseif sum(size(AA) == [1 4]) == 2 %incase pictre was opend usinf the bfopen
    BB = AA{1,1};
    clear AA
    AA = BB;
    clear BB  
end

%close all

%mask to detect puncta
PM_only_pic = cell(size(AA));
Mask_only = cell(size(AA));
for ii = 1:length(AA)
    PM_only_pic{ii,1} = zeros(size(AA{1,1}));
    Mask_only{ii,1} = zeros(size(AA{1,1}));
end

if isequal(lambda, 488)   
    % a way better cell mask
    [Gmag,~] = imgradient(AA{1, 1}*5000);
    G = Gmag > 10;
    BW_cell = imfill(~G,'holes');  %figure; imshow(BW_cell)
    BW_cell = imerode(BW_cell, ones(10,10));
    BW_cell_exlud = bwpropfilt(BW_cell,'Area',[20 10*10^7]); %figure; imshow(BW_G_exlud)
    Cell_mask = imdilate(BW_cell_exlud, ones(10,10)); %figure; imshow(Cell_mask)
    
elseif isequal(lambda, 647)
    Cell_mask = Cell_mask_488;
    %     C = imfuse(AA{1, 1}*150,Cell_mask,'falsecolor','Scaling','independent','ColorChannels',[2 1 2]);
    %     figure; imshow(C); title('488 mask on 647 image');
end

BB = uint16(zeros([size(AA{1,1}), size(AA,1)]));
prealoc_mat = BB;

for dd = 1:1000 
    BB(:,:,dd) = AA{dd,1};
end
  
%%% a two mask system one for all and one for remaing signal that is 3>SD
%%% of remaining signal
%     Gaus_1 = imgaussfilt(AA{k, 1}, 1, 'FilterSize',5); %figure; imagesc(Gaus_1, [0 1000])
%     Gaus_2 = imgaussfilt(AA{k, 1}, 3, 'FilterSize',5); %figure; imagesc(Gaus_2, [0 1000])
    
    Gaus_1 = imgaussfilt(BB, 2, 'FilterSize',7); %figure; imagesc(Gaus_1, [0 1000])
    Gaus_2 = imgaussfilt(BB, 4, 'FilterSize',7); %figure; imagesc(Gaus_2, [0 1000])
    Gaus_dif = (Gaus_1 - Gaus_2) .* uint16(Cell_mask); %figure; imagesc(Gaus_dif, [0 5])
    Gaus_dif_smooth = imgaussfilt(Gaus_dif,1); %figure; imagesc(Gaus_dif_smooth, [0 5])
    
    Gaus_dif_NaN = double(Gaus_dif_smooth);
    Gaus_dif_NaN(Gaus_dif_NaN < 0) = 0;
    Gaus_dif_NaN(Gaus_dif_NaN == 0) = NaN;
    
    mode_Gaus_dif = zeros(1,1,size(AA,1));
    prealoc_vect = mode_Gaus_dif;
    for dd = 1: size(AA,1)
        vec_Gaus_dif_NaN = Gaus_dif_NaN(:,:,dd);
        vec_Gaus_dif = vec_Gaus_dif_NaN(~isnan(vec_Gaus_dif_NaN));
        mode_Gaus_dif(:,:,dd) = mode(vec_Gaus_dif);
    end
    
    BW = Gaus_dif_NaN > mode_Gaus_dif ; %figure; imshow(BW);
    BW2 = logical(prealoc_mat);
    for dd = 1: size(AA,1)
        BW2(:,:,dd) = bwpropfilt(BW(:,:,dd), 'Area', [4 2000]); %figure; imshow(BW2);
    end
    BW3 = imdilate(BW2, ones(3,3)); %figure; imshow(BW3);
    
    remaining = imgaussfilt(double(BB),3) .* ~BW3 .* Cell_mask; %figure; imagesc(remaining, [0 5])
    remaining(remaining == 0) = NaN;
    
    std_rem = prealoc_vect;
    m_rem = prealoc_vect;
    
    for dd = 1: size(AA,1)
        std_rem(:,:,dd) = std(remaining(:,:,dd),[], 'all', 'omitnan');
        m_rem(:,:,dd) = mean(remaining(:,:,dd), 'all', 'omitnan');
    end
    
    %for anything left that is super brigt outside the mask
    mask_rem = remaining > (m_rem + 3*std_rem); %figure; imshow(mask_rem);
    %mask_rem = imdilate(mask_rem, ones(3,3));
    
    BW_cor = BW3 + mask_rem;
    BW_cor(BW_cor > 1) = 1;
       
    BW_cor = BW_cor .* Cell_mask; % figure; imshow(BW_cor);
    
%      C = imfuse(AA{k, 1}*150, BW_cor,'falsecolor','Scaling','independent','ColorChannels',[2 1 2]);
%      figure; imshow(C)

    mask = double(~BW_cor);
    %Mask_only{k,1} = BW_cor; 
    
    mask(mask == 0) = NaN;
    bb = double(BB) .* mask; %figure; imagesc(bb, [0 1000]);
    
    bbb = double(prealoc_mat);
    
    parfor dd = 1: size(AA,1)
        bbb(:,:,dd) = inpaint_nans(bb(:,:,dd), 2);  %figure; imagesc(bbb, [0 1000]);
    end
    
    
      %figure; imagesc(bbb, [0 1000]);   
%     buffer_blured_mask = imdilate(double(BW_cor), ones(3,3)); figure; imshow(buffer_blured_mask);
%     buffer_blured_crop =  bbb .* ~buffer_blured_mask; figure; imagesc(buffer_blured_crop, [0 1000]);
%     buffer_blured_signal =  imgaussfilt(bbb, 5) .* buffer_blured_mask; figure; imagesc(buffer_blured_signal, [0 1000]);
%     PM_only_pic{k, 1} = buffer_blured_crop + buffer_blured_signal; figure; imagesc(PM_only_pic{k, 1}, [0 1000]);
    
    BBBB = double(prealoc_mat);
    BBBB = imgaussfilt(bbb, 3); %figure; imagesc(PM_only_pic{k, 1}, [0 1000]);
    
    parfor dd = 1:1: size(AA,1)
        PM_only_pic{dd, 1} = BBBB(:,:,dd);
    end
    
    %PM_only_pic{k, 1} = medfilt2(bbb,[5 5]);   figure; imagesc(PM_only_pic{k, 1}, [0 1000]);


%%%% 

if smothing  == 1
    
    if size(PM_only_pic,1) > averaginng_window
        matrix_3d = zeros(size(PM_only_pic{1, 1},1), size(PM_only_pic{1, 1},2), size(PM_only_pic,1));
        matrix_3d_smooth = matrix_3d;
        %averaginng_window = 10; % make this always even as it will be devided by 2
        
        for k = 1:size(PM_only_pic, 1)
            matrix_3d(:,:,k) = PM_only_pic{k, 1};
        end
        
        padding_for_matrix_3d = nan(size(PM_only_pic{1, 1},1), size(PM_only_pic{1, 1},2), averaginng_window/2); %in agremwnt wityh movmean padding of a window
        
        matrix_3d = cat(3, padding_for_matrix_3d, matrix_3d, padding_for_matrix_3d);
        
        for k = (averaginng_window/2+1):size(matrix_3d, 3)-(averaginng_window/2)
            
            w = gausswin(averaginng_window + 1);
            gausian_mean = matrix_3d(:,:,k-(averaginng_window/2):k+(averaginng_window/2));
            sum_of_gaus_pic = zeros(size(matrix_3d(:,:,k)));
            sum_of_gaus_weight = 0;
            
            parfor dd = 1:size(w,1)
                if ~isnan(sum(gausian_mean(:,:,dd),'all'))
                    sum_of_gaus_pic = sum_of_gaus_pic + (gausian_mean(:,:,dd) .*w(dd));
                    sum_of_gaus_weight = sum_of_gaus_weight + w(dd);
                end
            end
            matrix_3d_smooth(:,:,k-(averaginng_window/2)) = sum_of_gaus_pic / sum_of_gaus_weight;
        end
        
        for k = 1:size(matrix_3d_smooth,3)
            PM_only_pic{k, 1} = matrix_3d_smooth(:,:,k);
        end
    end
end

% clims = [0 800];
% subplot(3,1,1); imagesc(AA{1, 1}, clims); title('Oryginal'); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]);
% subplot(3,1,2); imagesc(double(BW_eroded)); title('Puncta mask');  axis equal; axis([0 size(AA,2) 0 size(AA,1)]);
% subplot(3,1,3); imagesc(PM_only_pic{1, 1},clims); title('Puncta removed'); axis equal;  axis([0 size(AA,2) 0 size(AA,1)]);

% subplot(2,1,1); imagesc(AA, clims); title('Oryginal'); axis equal; axis([0 size(AA,2) 0 size(AA,1)]);
% subplot(2,1,2); imagesc(bbb_smooth,clims); title('Puncta removed'); axis equal;  axis([0 size(AA,2) 0 size(AA,1)]);
%
%
%cell_mat2tiff('C:\Users\tnawara\Desktop\Test of gap filing\Fake cells 2\PM_only_488_0010.tif', PM_only_pic)

%%%% itrs still notr perfect work on the bliking big poimts next maybe
%%%% anther mask just for big objects
end

% close all
% clims = [0 800];
% figure; imagesc(AA{1, 1}, clims); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]); set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);
% figure; imshow(Cell_mask); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]);
% figure; imagesc(Gaus_1, clims); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]); set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);
% figure; imagesc(Gaus_2, clims); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]); set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);
% figure; imagesc(Gaus_dif, [0 100]); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]); set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);
% figure; imshow(BW3); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]);
% figure; imagesc(remaining, clims); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]); set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);
% figure; imshow(mask_rem); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]);
% figure; imshow(BW_cor); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]);
% figure; imagesc(bb, clims); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]); set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);
% figure; imagesc(bbb, clims); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]); set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);
% figure; imagesc(PM_only_pic{1,1}, clims); axis equal; axis([0 size(AA{1,1},2) 0 size(AA{1,1},1)]); set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);