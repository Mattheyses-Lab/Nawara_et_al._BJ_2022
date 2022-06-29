function modeling_for_gap_fill
clear; close all;

background = 200;
diffraction_signal = 300;
bigger_signal = 1500;

mat_size = [300, 600];
n_diffraction_limited = 500;
n_bigger = 200;

modeled_background_no_signal = zeros(mat_size);
modeled_background_no_signal(:,:) = background; %generates the background picture
signal_only_mask_difraction = zeros(mat_size);
signal_only_mask_bigger_structres = zeros(mat_size);
spread_diffraction_signal_noisy = zeros(mat_size);
spread_bigger_signal_noisy = zeros(mat_size);

modeled_poisson_pictre = poissrnd(modeled_background_no_signal); %generates the poisson noise added background picture

clims = [min(modeled_poisson_pictre, [],'all'), max(modeled_poisson_pictre, [],'all')];
plot_axis = [1 size(modeled_background_no_signal,2) 1 size(modeled_background_no_signal,1)];

x = rand(n_diffraction_limited,1) * mat_size(2); %geenrates random x coordiantes acros the size of image for diffraction limited structres
y = rand(n_diffraction_limited,1) * mat_size(1); %geenrates random x coordiantes acros the size of image for diffraction limited structres

xy = [x y];

%figure; imagesc(modeled_background_no_signal, clims); axis equal; axis(plot_axis); hold on; scatter(xy(:,1), xy(:,2), 'r*'); hold off;
%figure; imagesc(modeled_poisson_pictre, clims); axis equal; axis(plot_axis); hold on; scatter(xy(:,1), xy(:,2), 'r*'); hold off;

[columnsInImage, rowsInImage] = meshgrid(1:mat_size(2), 1:mat_size(1));
circlePixels = false(mat_size(1),mat_size(2));
final_mask = false(mat_size(1),mat_size(2));

for radius = 1:3
    for ii = 1:size(xy,1)
        centerX = xy(ii,1);
        centerY = xy(ii,2);
        %radius = 2;
        circlePixels = (rowsInImage - centerY).^2 ...
            + (columnsInImage - centerX).^2 <= radius.^2;
        signal_only_mask_difraction = signal_only_mask_difraction + circlePixels;
    end
    
    signal_only_mask_difraction(signal_only_mask_difraction > 1) = 1; %figure; imagesc(signal_only_mask_difraction); axis equal; axis(plot_axis);
    diffraction_signal_pic = signal_only_mask_difraction * (diffraction_signal/radius); %figure; imagesc(diffraction_signal_pic)
    diffraction_signal_noisy = poissrnd(diffraction_signal_pic); %figure; imagesc(diffraction_signal_noisy); axis equal; axis(plot_axis);
    spread_diffraction_signal_noisy = spread_diffraction_signal_noisy + diffraction_signal_noisy;
end

noisy_back_plus_signal = modeled_poisson_pictre + spread_diffraction_signal_noisy; %figure; imagesc(noisy_back_plus_signal);

x_big = rand(n_bigger,1) * mat_size(2); %geenrates random x coordiantes acros the size of image for bigger structres
y_big = rand(n_bigger,1) * mat_size(1); %geenrates random x coordiantes acros the size of image for bigger structres

xy_big = [x_big y_big];

for radius = 2:5
    for ii = 1:size(xy_big,1)
        centerX = xy_big(ii,1);
        centerY = xy_big(ii,2);
        %radius = 5;
        circlePixels = (rowsInImage - centerY).^2 ...
            + (columnsInImage - centerX).^2 <= radius.^2;
        signal_only_mask_bigger_structres = signal_only_mask_bigger_structres + circlePixels;
    end
    
    signal_only_mask_bigger_structres(signal_only_mask_bigger_structres > 1) = 1; %figure; imagesc(signal_only_mask_bigger_structres);; axis equal; axis(plot_axis);
    bigger_signal_pic = signal_only_mask_bigger_structres * (bigger_signal/(radius)); %figure; imagesc(bigger_signal_pic)
    bigger_signal_noisy = poissrnd(bigger_signal_pic); %figure; imagesc(bigger_signal_noisy); axis equal; axis(plot_axis);
    spread_bigger_signal_noisy = spread_bigger_signal_noisy + bigger_signal_noisy;
end

noisy_back_plus_signal = noisy_back_plus_signal + spread_bigger_signal_noisy; 

%figure; imagesc(noisy_back_plus_signal); axis equal; axis(plot_axis);

%generates floppy PM matrix
 N = mat_size; % size in pixels of image
 F = 5;        % frequency-filter width
 [X,Y] = ndgrid(1:N(1),1:N(2));
 i = min(X-1,N(1)-X+1);
 j = min(Y-1,N(2)-Y+1);
 H = exp(-.5*(i.^2+j.^2)/F^2);
 Z = real(ifft2(H.*fft2(randn(N))));
 Z = rescale(Z); 
 %surf(X,Y,Z,'edgecolor','none'); light;

%frame = 50;
%pseudo_cell = imgaussfilt(noisy_back_plus_signal, 1);
% make a shape of the cell
noisy_back_plus_signal = noisy_back_plus_signal .* Z;
figure; imagesc(noisy_back_plus_signal, clims); axis equal; axis(plot_axis);
h = drawfreehand('Smoothing', 0);
ROI_mask = createMask(h);
pseudo_cell = noisy_back_plus_signal .* ROI_mask;
pseudo_cell = round(imgaussfilt(pseudo_cell, 1));
%pseudo_cell = padarray(pseudo_cell,[frame frame],0,'both');
%chage zeros to sometning between 0-5
id = pseudo_cell == 0;
r = rand(1,nnz(id))*5;
pseudo_cell(id) = r;
clims = [0 max(pseudo_cell, [], 'all')];
figure; imagesc(pseudo_cell, clims); axis equal; axis(plot_axis);

modeled_poisson_pictre = modeled_poisson_pictre .* Z;
pseudo_PM = modeled_poisson_pictre .* ROI_mask;
pseudo_PM = round(imgaussfilt(pseudo_PM, 1));

%pseudo_PM = padarray(pseudo_PM,[frame frame],0,'both');
%chage zeros to sometning between 0-5 uses the same values from above
%id = pseudo_PM == 0;
pseudo_PM(id) = r;
figure; imagesc(pseudo_PM, clims); axis equal; axis(plot_axis);


end

% AA = cell(1,1);
%AA{1,1} = uint16(pseudo_cell);
% gap_filling_testdata(AA, 488, [])
% set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);
% title('Simulated cell')
% 
% figure; imagesc(pseudo_cell, clims); axis equal; axis([1 cell_pic_size(2) 1 cell_pic_size(1)]);
%  set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);
%  title('Simulated cell')
% 
% 
% figure; imagesc(PM_only_pic{1, 1}, [0 1900]); axis equal; axis([1 150 1 200]);
% set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);
% title('PM only')