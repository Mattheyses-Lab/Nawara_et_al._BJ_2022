function modeling_for_gap_fill_V4(n_cell, PSF_infocus)
%PSF_infocus can be generated using FiJi PSF generator, its the middle
%slice (NA 1.48, ref of media 1.515, pix size 110nm, born and wolf 3d,
%lambda = 488)

close all;
PSF_infocus = double(PSF_infocus);

generate_sim = 1;

while generate_sim == 1
    
    background = 300; % put aboout 1.5 times more than what you want as it will be scaled down
    diffraction_signal = [50, 200, 350];
    bigger_signal = [200 400 600];
    
    mat_size = [300, 600];
    n_diffraction_limited = 1000;
    n_bigger = 200;
    
    modeled_background_no_signal = zeros(mat_size);
    modeled_background_no_signal(:,:) = background; %generates the background picture
    signal_only_mask_difraction = zeros(mat_size);
    signal_only_mask_bigger_structres = zeros(mat_size);
    spread_diffraction_signal = zeros(mat_size);
    %spread_diffraction_signal_noisy = zeros(mat_size);
    spread_bigger_signal = zeros(mat_size);
    %spread_bigger_signal_noisy = zeros(mat_size);
    
    %modeled_poisson_pictre = poissrnd(modeled_background_no_signal); %generates the poisson noise added background picture
    
    %clims = [min(modeled_poisson_pictre, [],'all'), max(modeled_poisson_pictre, [],'all')];
    plot_axis = [1 size(modeled_background_no_signal,2) 1 size(modeled_background_no_signal,1)];
    
    x = rand(n_diffraction_limited,1) * mat_size(2); %geenrates random x coordiantes acros the size of image for diffraction limited structres
    y = rand(n_diffraction_limited,1) * mat_size(1); %geenrates random y coordiantes acros the size of image for diffraction limited structres
    
    xy = [x y];
    
    %figure; imagesc(modeled_background_no_signal, clims); axis equal; axis(plot_axis); hold on; scatter(xy(:,1), xy(:,2), 'r*'); hold off;
    %figure; imagesc(modeled_poisson_pictre, clims); axis equal; axis(plot_axis); hold on; scatter(xy(:,1), xy(:,2), 'r*'); hold off;
    
    [columnsInImage, rowsInImage] = meshgrid(1:mat_size(2), 1:mat_size(1));
    circlePixels = false(mat_size(1),mat_size(2));
    final_mask = false(mat_size(1),mat_size(2));
    
    for nn = 1:3
        %for radius = 1:3
            for ii = nn:3:size(xy,1)
                centerX = xy(ii,1);
                centerY = xy(ii,2);
                radius = 2;
                circlePixels = (rowsInImage - centerY).^2 ...
                    + (columnsInImage - centerX).^2 <= radius.^2;
                signal_only_mask_difraction = signal_only_mask_difraction + circlePixels;
            end
            
            signal_only_mask_difraction(signal_only_mask_difraction > 1) = 1; %figure; imagesc(signal_only_mask_difraction); axis equal; axis(plot_axis);
            diffraction_signal_pic = signal_only_mask_difraction * (diffraction_signal(nn)); %figure; imagesc(diffraction_signal_pic)
            spread_diffraction_signal = spread_diffraction_signal + diffraction_signal_pic; %figure; imagesc(spread_diffraction_signal)
            %         diffraction_signal_noisy = poissrnd(diffraction_signal_pic); %figure; imagesc(diffraction_signal_noisy); axis equal; axis(plot_axis);
            %         spread_diffraction_signal_noisy = spread_diffraction_signal_noisy + diffraction_signal_noisy; %figure; imagesc(spread_diffraction_signal_noisy)
            signal_only_mask_difraction = zeros(mat_size);
        %end
    end
    
    %figure; imagesc(spread_diffraction_signal)
    %spread_diffraction_signal_PSFed = conv2(spread_diffraction_signal, PSF_infocus, 'same');
    %figure; imagesc(spread_diffraction_signal_PSFed)
    
    x_big = rand(n_bigger,1) * mat_size(2); %geenrates random x coordiantes acros the size of image for bigger structres
    y_big = rand(n_bigger,1) * mat_size(1); %geenrates random x coordiantes acros the size of image for bigger structres
    
    xy_big = [x_big y_big];
    
    for nn = 1:3
        %for radius = (1+nn):(4+nn)
            for ii = nn:3:size(xy_big,1)
                centerX = xy_big(ii,1);
                centerY = xy_big(ii,2);
                radius = 5;
                circlePixels = (rowsInImage - centerY).^2 ...
                    + (columnsInImage - centerX).^2 <= radius.^2;
                signal_only_mask_bigger_structres = signal_only_mask_bigger_structres + circlePixels;
            end
            
            signal_only_mask_bigger_structres(signal_only_mask_bigger_structres > 1) = 1; %figure; imagesc(signal_only_mask_bigger_structres);; axis equal; axis(plot_axis);
            bigger_signal_pic = signal_only_mask_bigger_structres * (bigger_signal(nn)); %figure; imagesc(bigger_signal_pic)
            spread_bigger_signal = spread_bigger_signal + bigger_signal_pic; %figure; imagesc(spread_bigger_signal); axis equal; axis(plot_axis);
            %     bigger_signal_noisy = poissrnd(bigger_signal_pic); %figure; imagesc(bigger_signal_noisy); axis equal; axis(plot_axis);
            %     spread_bigger_signal_noisy = spread_bigger_signal_noisy + bigger_signal_noisy;
            signal_only_mask_bigger_structres = zeros(mat_size); %figure; imagesc(signal_only_mask_bigger_structres); axis equal; axis(plot_axis);
        %end
    end
    
    %figure; imagesc(spread_bigger_signal)
    %spread_bigger_signal_PSFed = conv2(spread_bigger_signal, PSF_infocus, 'same');
    %figure; imagesc(spread_bigger_signal_PSFed)
    
    %generates Intensity variation mask
    N = mat_size; % size in pixels of image
    F = 40;        % frequency-filter width will chage the cut of bigger signal
    [X,Y] = ndgrid(1:N(1),1:N(2));
    i = min(X-1,N(1)-X+1);
    j = min(Y-1,N(2)-Y+1);
    H = exp(-.8*(i.^2+j.^2)/F^2);
    Z = real(ifft2(H.*fft2(randn(N))));
    Z = rescale(Z, -1, 2);
    Z(Z < 0) = 0; %surf(X,Y,Z,'edgecolor','none'); light;
    
    spread_bigger_signal_irregular = spread_bigger_signal .* Z; %figure; imagesc(spread_bigger_signal_irregular); axis equal; axis(plot_axis);
    
    %figure; imagesc(noisy_back_plus_signal); axis equal; axis(plot_axis);
    
    %generates floppy PM matrix
    N = mat_size; % size in pixels of image
    F = 7;        % frequency-filter width
    [X,Y] = ndgrid(1:N(1),1:N(2));
    i = min(X-1,N(1)-X+1);
    j = min(Y-1,N(2)-Y+1);
    H = exp(-1*(i.^2+j.^2)/F^2);
    Z = real(ifft2(H.*fft2(randn(N))));
    Z = rescale(Z);
    Z(Z < 0.2) = 0.2;
    Z(Z > 0.8) = 0.8;
    %Z = imgaussfilt(Z, 2);
    % figure; surf(X,Y,Z,'edgecolor','none'); light;
    
    spread_bigger_signal_irregular(Z < 0.6) = 0; % rmoves big brigts from dim PM areas (not prexent in real data)
    spread_diffraction_signal(Z < 0.35) = 0; % rmoves small brigts from dim PM areas (not prexent in real data)
    
    %%%% add them all back otheher also add noise
    %back_plus_signal = (modeled_background_no_signal + spread_diffraction_signal + spread_bigger_signal_irregular);
    %figure; imagesc(noisy_back_plus_signal, [0 1000]); axis equal; axis(plot_axis);
    
    %frame = 50;
    %pseudo_cell = imgaussfilt(noisy_back_plus_signal, 1);
    % make a shape of the cell
    clims = [0 3000];
    modeled_background = modeled_background_no_signal .* Z;
    modeled_background_PSFed = conv2(modeled_background, PSF_infocus, 'same');
    modeled_signal = (spread_diffraction_signal + spread_bigger_signal_irregular) .* Z;
    modeled_signal_PSFed = conv2(modeled_signal, PSF_infocus, 'same');
    back_plus_signal = modeled_background_PSFed + modeled_signal_PSFed;
    imagesc(back_plus_signal, clims); axis equal; axis(plot_axis);
    
    answer = questdlg('Accept simulation?', ...
        'Simulation conformation', ...
        'Yes','No', 'Just STOP...', 'Yes');
    % Handle response
    switch answer
        case 'Yes'
            generate_sim = 0;
        case 'No'
            generate_sim = 1;
        case 'Just STOP...'
           close all
           return
    end
end
h = drawfreehand('Smoothing', 3);
ROI_mask = createMask(h);close all;

ROI_cell_edge = imdilate(ROI_mask, ones(30,30));
ROI_cell_edge = ROI_cell_edge .* ~ROI_mask; %figure; imshow(ROI_cell_edge)
ROI_cell_edge = double(~ROI_cell_edge);
ROI_cell_edge(ROI_cell_edge == 0) = NaN;
%  ww = window2(mat_size(1),mat_size(2),@gausswin);
%  ww(ROI_mask == 1) = 1; h = surf(ww); set(h,'LineStyle','none')
%  ww = ww .*;

pseudo_PM_cell = modeled_background_PSFed .* ROI_mask;
%pseudo_cell = imgaussfilt(pseudo_cell, 1);
%pseudo_cell = padarray(pseudo_cell,[frame frame],0,'both');
%chage zeros to sometning between 0-5
id = pseudo_PM_cell == 0;
r = rand(1,nnz(id))*20;
pseudo_PM_cell(id) = r;

pseudo_PM = pseudo_PM_cell .* ROI_cell_edge;
pseudo_PM = inpaint_nans(pseudo_PM, 1);
pseudo_PM_noisy = poissrnd(abs(pseudo_PM));
pseudo_PM_noisy = imgaussfilt(pseudo_PM_noisy, 1); %figure; imagesc(pseudo_PM_noisy, clims); axis equal; axis(plot_axis);

pseudo_signal_raw = modeled_signal_PSFed .* ROI_mask; %figure; imagesc(pseudo_signal_raw, clims); axis equal; axis(plot_axis);
pseudo_signal_noisy = poissrnd(abs(pseudo_signal_raw));
pseudo_signal_noisy = imgaussfilt(pseudo_signal_noisy, 1); %figure; imagesc(pseudo_signal_noisy, clims); axis equal; axis(plot_axis);


pseudo_cell = pseudo_PM_noisy + pseudo_signal_noisy;


clims = [0 max(pseudo_cell, [], 'all')];
figure; imagesc(pseudo_cell, clims); axis equal; axis(plot_axis);
figure; imagesc(pseudo_PM_noisy, clims); axis equal; axis(plot_axis);

% modeled_background_no_signal = modeled_background_no_signal .* Z;
% pseudo_PM = modeled_background_no_signal .* ROI_mask;
% pseudo_PM = imgaussfilt(pseudo_PM, 1);
% 
% %pseudo_PM = padarray(pseudo_PM,[frame frame],0,'both');
% %chage zeros to sometning between 0-5 uses the same values from above
% %id = pseudo_PM == 0;
% pseudo_PM(id) = r;
%figure; imshowpair(pseudo_cell,pseudo_PM);% axis equal; axis(plot_axis);



cell_mat2tiff(['C:\Users\tnawara\Desktop\Pseudo_PM_488.00', num2str(n_cell), '.tif'], pseudo_PM_noisy)
cell_mat2tiff(['C:\Users\tnawara\Desktop\Pseudo_Cell_488.00', num2str(n_cell), '.tif'], pseudo_cell)

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