function [error_pic_PM_masked_pixels] = error_percent_meas(Pseudo_cell, Pseudo_PM, PM_only)
%close all
%Calculates error as ((PM_only - Pseudo_PM) / Pseudo_Cell) * 100
% erode = 5;

%mask based on resolved pictre so edge is deicarded
[Gmag,~] = imgradient(Pseudo_cell*1000);
    G = Gmag > 100;
    BW_cell = imfill(~G,'holes'); % figure; imshow(BW_cell)
    BW_cell = imerode(BW_cell, ones(10,10));
    BW_cell_exlud = bwpropfilt(BW_cell,'Area',[20 10*10^7]); %figure; imshow(BW_G_exlud)
    Cell_mask = imdilate(BW_cell_exlud, ones(10,10)); %figure; imshow(Cell_mask)


% BW_cell = imerode(Cell_mask, ones(erode+5,erode+5)); %figure; imshow(BW_cell);
% BW_cell(1:erode,1:size(BW_cell,2)) = 0;
% BW_cell(size(BW_cell,1)-(erode-1):size(BW_cell,1) ,1:size(BW_cell,2)) = 0;
% BW_cell(1:size(BW_cell,1),1:erode) = 0;
% BW_cell(1:size(BW_cell,1),size(BW_cell,2)-(erode-1):size(BW_cell,2)) = 0;
    
%  C = imfuse(PM_only,BW_cell,'falsecolor','Scaling','independent','ColorChannels',[2 1 2]);
%  figure; imshow(C); title('PM only and cell mask overlay');
BW_cell = double(Cell_mask);
BW_cell(BW_cell == 0) = NaN;
error_pic_PM = ((double(PM_only) - double(Pseudo_PM)) ./ double(Pseudo_PM)) * 100; %figure; imagesc(error_pic_PM, [-100 100]);
error_pic_PM_masked = error_pic_PM .* BW_cell; %figure; imagesc(error_pic_PM_masked, [-100 100]);

error_pic_PM_masked_pixels = error_pic_PM_masked(:);
error_pic_PM_masked_pixels = error_pic_PM_masked_pixels(~isnan(error_pic_PM_masked_pixels));
mean_error = mean(error_pic_PM_masked_pixels);
std_error = std(error_pic_PM_masked_pixels);

high_STD_zone_neg = error_pic_PM_masked < (std_error * -1);
high_STD_zone_poz = error_pic_PM_masked > std_error;
high_STD_zone = high_STD_zone_neg + high_STD_zone_poz;
high_STD_zone(high_STD_zone > 1) = 1; %figure; imagesc(high_STD_zone);

high_2STD_zone_neg = error_pic_PM_masked < (std_error * -2);
high_2STD_zone_poz = error_pic_PM_masked > (std_error * 2);
high_2STD_zone = high_2STD_zone_neg + high_2STD_zone_poz;
high_2STD_zone(high_2STD_zone > 1) = 1; %figure; imagesc(high_2STD_zone);

% C = imfuse(Pseudo_cell*150,high_STD_zone,'falsecolor','Scaling','independent','ColorChannels',[2 1 2]);
% figure; imshow(C); title('Areas below and above SD');
figure
all_together = cat(3, high_2STD_zone*100000, Pseudo_cell*100, high_STD_zone*100000); imshow(all_together); 
title(sprintf('Blue error > %.2f%%, Magenta error > %.2f%%', std_error, std_error*2));
    
end
% 
% errr1 = error_percent_meas(Pseudo_Cell_488_001, Pseudo_PM_488_001, PM_only_488_001);
% errr2 = error_percent_meas(Pseudo_Cell_488_002, Pseudo_PM_488_002, PM_only_488_002);
% errr3 = error_percent_meas(Pseudo_Cell_488_003, Pseudo_PM_488_003, PM_only_488_003);
% errr4 = error_percent_meas(Pseudo_Cell_488_004, Pseudo_PM_488_004, PM_only_488_004);
% errr5 = error_percent_meas(Pseudo_Cell_488_005, Pseudo_PM_488_005, PM_only_488_005);
% errr6 = error_percent_meas(Pseudo_Cell_488_006, Pseudo_PM_488_006, PM_only_488_006);
% errr7 = error_percent_meas(Pseudo_Cell_488_007, Pseudo_PM_488_007, PM_only_488_007);
% errr8 = error_percent_meas(Pseudo_Cell_488_008, Pseudo_PM_488_008, PM_only_488_008);
% errr9 = error_percent_meas(Pseudo_Cell_488_009, Pseudo_PM_488_009, PM_only_488_009);
% errr10 = error_percent_meas(Pseudo_Cell_488_0010, Pseudo_PM_488_0010, PM_only_488_0010);