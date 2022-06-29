function TIRF_cell_modeling(n_cell)
close all;

%% PSF generated using PSF_infocus can be generated using FiJi PSF generator, its the middle
%slice (NA 1.48, ref of media 1.515, pix size 110nm, born and wolf 3d,
%lambda = 488, size 21)
PSF_infocus = [2.22232956730295e-05,7.74849104345776e-05,0.000108448111859616,6.28387351753190e-05,4.43680874013808e-05,4.04352686018683e-05,9.18487203307450e-05,0.000167843216331676,0.000223630195250735,0.000257747800787911,0.000269232696155086,0.000257747800787911,0.000223630195250735,0.000167843216331676,9.18487203307450e-05,4.04352686018683e-05,4.43680874013808e-05,6.28387351753190e-05,0.000108448111859616,7.74849104345776e-05,2.22232956730295e-05;...
    7.74849104345776e-05,0.000102491030702367,5.02926995977759e-05,4.22435077780392e-05,8.11816353234462e-05,0.000201157876290381,0.000247673364356160,0.000196033171960153,0.000157906426466070,0.000134490313939750,0.000126590020954609,0.000134490313939750,0.000157906426466070,0.000196033171960153,0.000247673364356160,0.000201157876290381,8.11816353234462e-05,4.22435077780392e-05,5.02926995977759e-05,0.000102491030702367,7.74849104345776e-05;...
    0.000108448111859616,5.02926995977759e-05,4.15244576288387e-05,0.000124128870083950,0.000269232696155086,0.000188494246685877,0.000124249709188007,0.000107440362626221,9.49345048866235e-05,8.72093369252980e-05,8.45948088681325e-05,8.72093369252980e-05,9.49345048866235e-05,0.000107440362626221,0.000124249709188007,0.000188494246685877,0.000269232696155086,0.000124128870083950,4.15244576288387e-05,5.02926995977759e-05,0.000108448111859616;...
    6.28387351753190e-05,4.22435077780392e-05,0.000124128870083950,0.000254896382102743,0.000157906426466070,0.000109889588202350,8.72093369252980e-05,0.000351401889929548,0.000584486348088831,0.000729645136743784,0.000778994639404118,0.000729645136743784,0.000584486348088831,0.000351401889929548,8.72093369252980e-05,0.000109889588202350,0.000157906426466070,0.000254896382102743,0.000124128870083950,4.22435077780392e-05,6.28387351753190e-05;...
    4.43680874013808e-05,8.11816353234462e-05,0.000269232696155086,0.000157906426466070,0.000104974307760131,0.000216357395402156,0.000632405048236251,0.000689819746185094,0.000572574092075229,0.000498680572491139,0.000473387772217393,0.000498680572491139,0.000572574092075229,0.000689819746185094,0.000632405048236251,0.000216357395402156,0.000104974307760131,0.000157906426466070,0.000269232696155086,8.11816353234462e-05,4.43680874013808e-05;...
    4.04352686018683e-05,0.000201157876290381,0.000188494246685877,0.000109889588202350,0.000216357395402156,0.000729645136743784,0.000596585276070982,0.000443090946646407,0.000363196944817901,0.000311913958285004,0.000294167664833367,0.000311913958285004,0.000363196944817901,0.000443090946646407,0.000596585276070982,0.000729645136743784,0.000216357395402156,0.000109889588202350,0.000188494246685877,0.000201157876290381,4.04352686018683e-05;...
    9.18487203307450e-05,0.000247673364356160,0.000124249709188007,8.72093369252980e-05,0.000632405048236251,0.000596585276070982,0.000411889166571200,0.000294167664833367,0.00233249855227768,0.00368026853539050,0.00415563723072410,0.00368026853539050,0.00233249855227768,0.000294167664833367,0.000411889166571200,0.000596585276070982,0.000632405048236251,8.72093369252980e-05,0.000124249709188007,0.000247673364356160,9.18487203307450e-05;...
    0.000167843216331676,0.000196033171960153,0.000107440362626221,0.000351401889929548,0.000689819746185094,0.000443090946646407,0.000294167664833367,0.00321868737228215,0.00409260857850313,0.00402177870273590,0.00399584881961346,0.00402177870273590,0.00409260857850313,0.00321868737228215,0.000294167664833367,0.000443090946646407,0.000689819746185094,0.000351401889929548,0.000107440362626221,0.000196033171960153,0.000167843216331676;...
    0.000223630195250735,0.000157906426466070,9.49345048866235e-05,0.000584486348088831,0.000572574092075229,0.000363196944817901,0.00233249855227768,0.00409260857850313,0.00412771664559841,0.00458299228921533,0.00476442975923419,0.00458299228921533,0.00412771664559841,0.00409260857850313,0.00233249855227768,0.000363196944817901,0.000572574092075229,0.000584486348088831,9.49345048866235e-05,0.000157906426466070,0.000223630195250735;...
    0.000257747800787911,0.000134490313939750,8.72093369252980e-05,0.000729645136743784,0.000498680572491139,0.000311913958285004,0.00368026853539050,0.00402177870273590,0.00458299228921533,0.171250596642494,0.288974076509476,0.171250596642494,0.00458299228921533,0.00402177870273590,0.00368026853539050,0.000311913958285004,0.000498680572491139,0.000729645136743784,8.72093369252980e-05,0.000134490313939750,0.000257747800787911;...
    0.000269232696155086,0.000126590020954609,8.45948088681325e-05,0.000778994639404118,0.000473387772217393,0.000294167664833367,0.00415563723072410,0.00399584881961346,0.00476442975923419,0.288974076509476,1,0.288974076509476,0.00476442975923419,0.00399584881961346,0.00415563723072410,0.000294167664833367,0.000473387772217393,0.000778994639404118,8.45948088681325e-05,0.000126590020954609,0.000269232696155086;...
    0.000257747800787911,0.000134490313939750,8.72093369252980e-05,0.000729645136743784,0.000498680572491139,0.000311913958285004,0.00368026853539050,0.00402177870273590,0.00458299228921533,0.171250596642494,0.288974076509476,0.171250596642494,0.00458299228921533,0.00402177870273590,0.00368026853539050,0.000311913958285004,0.000498680572491139,0.000729645136743784,8.72093369252980e-05,0.000134490313939750,0.000257747800787911;...
    0.000223630195250735,0.000157906426466070,9.49345048866235e-05,0.000584486348088831,0.000572574092075229,0.000363196944817901,0.00233249855227768,0.00409260857850313,0.00412771664559841,0.00458299228921533,0.00476442975923419,0.00458299228921533,0.00412771664559841,0.00409260857850313,0.00233249855227768,0.000363196944817901,0.000572574092075229,0.000584486348088831,9.49345048866235e-05,0.000157906426466070,0.000223630195250735;...
    0.000167843216331676,0.000196033171960153,0.000107440362626221,0.000351401889929548,0.000689819746185094,0.000443090946646407,0.000294167664833367,0.00321868737228215,0.00409260857850313,0.00402177870273590,0.00399584881961346,0.00402177870273590,0.00409260857850313,0.00321868737228215,0.000294167664833367,0.000443090946646407,0.000689819746185094,0.000351401889929548,0.000107440362626221,0.000196033171960153,0.000167843216331676;...
    9.18487203307450e-05,0.000247673364356160,0.000124249709188007,8.72093369252980e-05,0.000632405048236251,0.000596585276070982,0.000411889166571200,0.000294167664833367,0.00233249855227768,0.00368026853539050,0.00415563723072410,0.00368026853539050,0.00233249855227768,0.000294167664833367,0.000411889166571200,0.000596585276070982,0.000632405048236251,8.72093369252980e-05,0.000124249709188007,0.000247673364356160,9.18487203307450e-05;...
    4.04352686018683e-05,0.000201157876290381,0.000188494246685877,0.000109889588202350,0.000216357395402156,0.000729645136743784,0.000596585276070982,0.000443090946646407,0.000363196944817901,0.000311913958285004,0.000294167664833367,0.000311913958285004,0.000363196944817901,0.000443090946646407,0.000596585276070982,0.000729645136743784,0.000216357395402156,0.000109889588202350,0.000188494246685877,0.000201157876290381,4.04352686018683e-05;...
    4.43680874013808e-05,8.11816353234462e-05,0.000269232696155086,0.000157906426466070,0.000104974307760131,0.000216357395402156,0.000632405048236251,0.000689819746185094,0.000572574092075229,0.000498680572491139,0.000473387772217393,0.000498680572491139,0.000572574092075229,0.000689819746185094,0.000632405048236251,0.000216357395402156,0.000104974307760131,0.000157906426466070,0.000269232696155086,8.11816353234462e-05,4.43680874013808e-05;...
    6.28387351753190e-05,4.22435077780392e-05,0.000124128870083950,0.000254896382102743,0.000157906426466070,0.000109889588202350,8.72093369252980e-05,0.000351401889929548,0.000584486348088831,0.000729645136743784,0.000778994639404118,0.000729645136743784,0.000584486348088831,0.000351401889929548,8.72093369252980e-05,0.000109889588202350,0.000157906426466070,0.000254896382102743,0.000124128870083950,4.22435077780392e-05,6.28387351753190e-05;...
    0.000108448111859616,5.02926995977759e-05,4.15244576288387e-05,0.000124128870083950,0.000269232696155086,0.000188494246685877,0.000124249709188007,0.000107440362626221,9.49345048866235e-05,8.72093369252980e-05,8.45948088681325e-05,8.72093369252980e-05,9.49345048866235e-05,0.000107440362626221,0.000124249709188007,0.000188494246685877,0.000269232696155086,0.000124128870083950,4.15244576288387e-05,5.02926995977759e-05,0.000108448111859616;...
    7.74849104345776e-05,0.000102491030702367,5.02926995977759e-05,4.22435077780392e-05,8.11816353234462e-05,0.000201157876290381,0.000247673364356160,0.000196033171960153,0.000157906426466070,0.000134490313939750,0.000126590020954609,0.000134490313939750,0.000157906426466070,0.000196033171960153,0.000247673364356160,0.000201157876290381,8.11816353234462e-05,4.22435077780392e-05,5.02926995977759e-05,0.000102491030702367,7.74849104345776e-05;...
    2.22232956730295e-05,7.74849104345776e-05,0.000108448111859616,6.28387351753190e-05,4.43680874013808e-05,4.04352686018683e-05,9.18487203307450e-05,0.000167843216331676,0.000223630195250735,0.000257747800787911,0.000269232696155086,0.000257747800787911,0.000223630195250735,0.000167843216331676,9.18487203307450e-05,4.04352686018683e-05,4.43680874013808e-05,6.28387351753190e-05,0.000108448111859616,7.74849104345776e-05,2.22232956730295e-05];
%figure; imagesc(PSF_infocus, [0.00000001, 0.01])
%%
generate_sim = 1;

while generate_sim == 1
    
    background = 300; % put aboout 1.5 times more than what you want as it will be scaled down
    diffraction_signal = [200, 300, 400];
    bigger_signal = [300, 600, 900];
    
    mat_size = [300, 600]; % 300 600
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
    y = rand(n_diffraction_limited,1) * mat_size(1); %geenrates random x coordiantes acros the size of image for diffraction limited structres
    
    xy = [x y];
    
    %figure; imagesc(modeled_background_no_signal, clims); axis equal; axis(plot_axis); hold on; scatter(xy(:,1), xy(:,2), 'r*'); hold off;
    %figure; imagesc(modeled_poisson_pictre, clims); axis equal; axis(plot_axis); hold on; scatter(xy(:,1), xy(:,2), 'r*'); hold off;
    
    [columnsInImage, rowsInImage] = meshgrid(1:mat_size(2), 1:mat_size(1));
    circlePixels = false(mat_size(1),mat_size(2));
    final_mask = false(mat_size(1),mat_size(2));
    
    for nn = 1:3
%         for radius = 1:3
            for ii = nn:3:size(xy,1)
                centerX = xy(ii,1);
                centerY = xy(ii,2);
                radius = randi(2);
                circlePixels = (rowsInImage - centerY).^2 ...
                    + (columnsInImage - centerX).^2 <= radius.^2;
                signal_only_mask_difraction = signal_only_mask_difraction + circlePixels;
            end
            
            signal_only_mask_difraction(signal_only_mask_difraction > 1) = 1; %figure; imagesc(signal_only_mask_difraction); axis equal; axis(plot_axis);
            diffraction_signal_pic = signal_only_mask_difraction * diffraction_signal(nn); %figure; imagesc(diffraction_signal_pic)
            spread_diffraction_signal = spread_diffraction_signal + diffraction_signal_pic; %figure; imagesc(spread_diffraction_signal, [0 1000]); ; axis equal; axis(plot_axis); set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
            %         diffraction_signal_noisy = poissrnd(diffraction_signal_pic); %figure; imagesc(diffraction_signal_noisy); axis equal; axis(plot_axis);
            %         spread_diffraction_signal_noisy = spread_diffraction_signal_noisy + diffraction_signal_noisy; %figure; imagesc(spread_diffraction_signal_noisy)
            signal_only_mask_difraction = zeros(mat_size);
%         end
    end
    
    %noisy_back_plus_signal = modeled_poisson_pictre + spread_diffraction_signal_noisy; %figure; imagesc(noisy_back_plus_signal);
    
    x_big = rand(n_bigger,1) * mat_size(2); %geenrates random x coordiantes acros the size of image for bigger structres
    y_big = rand(n_bigger,1) * mat_size(1); %geenrates random x coordiantes acros the size of image for bigger structres
    
    xy_big = [x_big y_big];
    
    for nn = 1:3
%         for radius = (1+nn):(4+nn)
            for ii = nn:3:size(xy_big,1)
                centerX = xy_big(ii,1);
                centerY = xy_big(ii,2);
                radius = randi([4 6]);
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
%         end
    end
    
    %generates Intensity variation mask
    N = mat_size; % size in pixels of image
    F = 40;        % frequency-filter width will chage the cut of bigger signal
    [X,Y] = ndgrid(1:N(1),1:N(2));
    i = min(X-1,N(1)-X+1);
    j = min(Y-1,N(2)-Y+1);
    H = exp(-.8*(i.^2+j.^2)/F^2);
    Z = real(ifft2(H.*fft2(randn(N))));
    Z = rescale(Z, -1, 2);
    Z(Z < 0) = 0; %surf(X,Y,Z,'edgecolor','none');
    
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
    BW_big = imbinarize(spread_bigger_signal_irregular);
    BW_big_err = imerode(BW_big, ones(3,3));
    spread_bigger_signal_irregular = spread_bigger_signal_irregular .* double(BW_big_err);
    spread_diffraction_signal(Z < 0.35) = 0; % rmoves small brigts from dim PM areas (not prexent in real data)
    
    %%%% add them all back otheher also add noise
    %back_plus_signal = (modeled_background_no_signal + spread_diffraction_signal + spread_bigger_signal_irregular);
    %figure; imagesc(noisy_back_plus_signal, [0 1000]); axis equal; axis(plot_axis);
    
    %frame = 50;
    %pseudo_cell = imgaussfilt(noisy_back_plus_signal, 1);
    % make a shape of the cell
    clims = [0 3000];
    modeled_background = modeled_background_no_signal .* Z;
    modeled_signal = (spread_diffraction_signal + spread_bigger_signal_irregular) .* Z;
    modeled_signal = conv2(modeled_signal, PSF_infocus, 'same'); %adds PSF blur to data
    modeled_signal = round(modeled_signal) - 1; %clears up the PM area
    modeled_signal(modeled_signal < 1) = 0; %clears up the PM area
    back_plus_signal =  modeled_background + modeled_signal;
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

pseudo_PM_cell = modeled_background .* ROI_mask;
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

pseudo_signal_raw = modeled_signal .* ROI_mask; %figure; imagesc(pseudo_signal_raw, clims); axis equal; axis(plot_axis);
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



cell_mat2tiff(['C:\Users\tnawara\Desktop\Test of gap filing\Fake cells 2\Pseudo_PM_488.00', num2str(n_cell), '.tif'], pseudo_PM_noisy)
cell_mat2tiff(['C:\Users\tnawara\Desktop\Test of gap filing\Fake cells 2\Pseudo_Cell_488.00', num2str(n_cell), '.tif'], pseudo_cell)

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