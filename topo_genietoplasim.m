function [] = topo_genietoplasim(fname, tx)
% topo_genietoplasim
%
%   *********************************************************
%   *** RE_GRID GENIE TOPOGRAPHY TO PLASIM ***
%   *********************************************************
%
% (1) read in '.k1' file
% (2) extract core data plus imax, jmax
% (3) extract topo
% (4) return '_surf_0172.sra' 32x64 topo file
%
%   *********************************************************
%
%   ***********************************************************************
%%
% *********************************************************************** %
% *** INITIALIZE ******************************************************** %
% *********************************************************************** %
%
disp(['       * Initializing grid ...']);
switch tx
    case 'T21'
        % create PLASIM grid
        imax = 32;
        jmax = 64;
        %lon-lat grid
        axislatce = [-88.5731, -83.0813, -77.5570, -72.0255, -66.4911, -60.9555, ...
        -55.4190, -49.8821, -44.3450, -38.8076, -33.2701, -27.7324, -22.1947, ...
        -16.6570, -11.1192, -5.5814, -0.0436, 5.4942, 11.0320, 16.5697, 22.1074, ...
        27.6451, 33.1826, 38.7200, 44.2571, 49.7940, 55.3305, 60.8661, 66.4005, ...
        71.9320, 77.4563, 82.9481, 90.0000];
        axislatce = transpose(axislatce);
        axislonce = [-2.8125, 2.8125, 8.4375, 14.0625, 19.6875, 25.3125, 30.9375 ...
        36.5625, 42.1875, 47.8125, 53.4375, 59.0625, 64.6875, 70.3125, 75.9375, 81.5625 ...
        87.1875, 92.8125, 98.4375, 104.0625, 109.6875, 115.3125, 120.9375, 126.5625 ..., 
        132.1875, 137.8125, 143.4375, 149.0625, 154.6875, 160.3125, 165.9375, 171.5625 ...
        177.1875, 182.8125, 188.4375, 194.0625, 199.6875, 205.3125, 210.9375, 216.5625 ...
        222.1875, 227.8125, 233.4375, 239.0625, 244.6875, 250.3125, 255.9375, 261.5625 ...
        267.1875, 272.8125, 278.4375, 284.0625, 289.6875, 295.3125, 300.9375, 306.5625 ...
        312.1875, 317.8125, 323.4375, 329.0625, 334.6875, 340.3125, 345.9375, 351.5625, 357.1875];
        axislonce = transpose(axislonce);
    case 'T42'
        % create PLASIM grid
        imax = 64;
        jmax = 128;
        %lon-lat grid
        axislatce = [-89.2701, -86.5028, -83.7192, -80.9319, -78.1432, ...
        -75.3538, -72.5640, -69.7740, -66.9839, -64.1936, -61.4033, -58.6129 ... 
        -55.8224, -53.0320, -50.2415, -47.4510, -44.6604, -41.8699, -39.0793 ...
        -36.2888, -33.4982, -30.7076, -27.9170, -25.1264, -22.3358, -19.5452 ...
        -16.7546, -13.9640, -11.1734, -8.3828, -5.5922, -2.8016, -0.0109 ...
        2.7797, 5.5703, 8.3609, 11.1515, 13.9421, 16.7327, 19.5233, 22.3139 ...
        25.1045, 27.8951, 30.6857, 33.4763, 36.2668, 39.0574, 41.8479 ...
        44.6385, 47.4290, 50.2195, 53.0099, 55.8004, 58.5908, 61.3811 ...
        64.1714, 66.9615, 69.7515, 72.5413, 75.3307, 78.1194, 80.9067 ...
        83.6903, 86.4576, 90.0000];
        axislatce = transpose(axislatce);
        axislonce = [-1.4062, 1.4062, 4.2188, 7.0312, 9.8438, 12.6562 ...
        15.4688, 18.2812, 21.0938, 23.9062, 26.7188, 29.5312, 32.3438 ...
        35.1562, 37.9688, 40.7812, 43.5938, 46.4062, 49.2188, 52.0312 ...
        54.8438, 57.6562, 60.4688, 63.2812, 66.0938, 68.9062, 71.7188 ...
        74.5312, 77.3438, 80.1562, 82.9688, 85.7812, 88.5938, 91.4062 ...
        94.2188, 97.0312, 99.8438, 102.6562, 105.4688, 108.2812, 111.0938 ...
        113.9062, 116.7188, 119.5312, 122.3438, 125.1562, 127.9688, 130.7812 ...
        133.5938, 136.4062, 139.2188, 142.0312, 144.8438, 147.6562, 150.4688 ...
        153.2812, 156.0938, 158.9062, 161.7188, 164.5312, 167.3438, 170.1562 ...
        172.9688, 175.7812, 178.5938, 181.4062, 184.2188, 187.0312, 189.8438 ...
        192.6562, 195.4688, 198.2812, 201.0938, 203.9062, 206.7188, 209.5312 ...
        212.3438, 215.1562, 217.9688, 220.7812, 223.5938, 226.4062, 229.2188 ...
        232.0312, 234.8438, 237.6562, 240.4688, 243.2812, 246.0938, 248.9062 ...
        251.7188, 254.5312, 257.3438, 260.1562, 262.9688, 265.7812, 268.5938 ...
        271.4062, 274.2188, 277.0312, 279.8438, 282.6562, 285.4688, 288.2812 ...
        291.0938, 293.9062, 296.7188, 299.5312, 302.3438, 305.1562, 307.9688 ...
        310.7812, 313.5938, 316.4062, 319.2188, 322.0312, 324.8438, 327.6562 ...
        330.4688, 333.2812, 336.0938, 338.9062, 341.7188, 344.5312, 347.3438 ...
        350.1562, 352.9688, 355.7812, 358.5938];
        axislonce = transpose(axislonce);
end
% determine output grid size (remember: [rows columns])
%[jmax, imax] = size(gomask); 
[go_lonm,go_lone,go_latm,go_late,go_dm,go_de] = make_genie_grid(36,36,16,5000.0,-260.0,true,0);
%[axis_imid,axis_iedge,axis_jmid,axis_jedge,axis_kmid,axis_kedge] = make_genie_grid(imax,jmax,kmax,par_max_D,par_lon_off,opt_equalarea,par_add_Dk)
%1x36, 1x37, 1x36, 1x37, 1x16, 1x17
% create GENIE u and v edge axes
% latitude
golatue = go_late;
golatve = [go_latm 90.0]; 
% longitude
golonue = [go_lonm go_lonm(end)+360.0/36];
golonve = go_lone;

% *********************************************************************** %
% *** load input data *************************************************** %
% *********************************************************************** %

disp(['       * Loading input data ...']);

loc_k1 = load([fname '.k1']);
%disp(size(loc_k1))
%process k1 file into mask
loc_k1(loc_k1 < 17) = 0;
loc_k1(loc_k1 > 90) = 1;
k1 = zeros(37,37);

k1(2:37,:) = loc_k1(2:37,2:38);
k1(1,:) = loc_k1(2,2:38);
%disp(size(k1))

% *********************************************************************** %

% *********************************************************************** %
% *** RE-GRID *********************************************************** %
% *********************************************************************** %
%
disp(['       * Regridding topography ...']);

switch tx
    case 'T21'
        %interpolate k1 onto 73x73 grid 
        [X,Y] = meshgrid(1:size(k1,2), 1:size(k1,1));
        [U,V] = meshgrid(1:0.5:size(k1,2), 1:0.5:size(k1,1));
        %k1_into = zeros(3)
        k1_interp = interp2(X,Y,k1,U,V,'cubic');
        k1_interpsort = k1_interp;
        %threshold 1
        thr1 = 0.9;
        k1_interpsort(k1_interpsort < thr1) = 0;
        %maxint = max(k1_interpsort, [], 'all');
        %thr2 = maxint - (0.07*maxint);
        k1_interpsort(k1_interpsort >= thr1) = 1;
        %interpolate lat, lon coords to match 
        golon_q = 1.5:36.5;
        golon_itp = interp1(golonue,golon_q);
        golongrid = zeros(1,73);
        for i=0:36
            golongrid(:,((2*i)+1)) = golonue(:,i+1);
        end
        for i=1:36
            golongrid(:,(2*i)) = golon_itp(:,i);
        end
        
        golat_q = 1.5:36.5;
        golat_itp = interp1(golatue,golat_q);
        golatgrid = zeros(1,73);
        for i=0:36
            golatgrid(:,((2*i)+1)) = golatue(:,i+1);
        end
        for i=1:36
            golatgrid(:,(2*i)) = golat_itp(:,i);
        end

        k1_into = zeros(72,72);
        k1_into = k1_interpsort(2:73,1:72);

        [sra_int, dscd] = make_regrid_2d(golongrid,golatgrid,k1_into',axislonce,axislatce,false);
        sra_int = sra_int';
        sra_int_cl = sra_int;
        %adjust binary threshold;
        thr2 = 0.2; 
        sra_int_cl(sra_int_cl > thr2) = 1;
        sra_int_cl(sra_int_cl <= thr2) = 0;
        %sra_int_cl = flipud(sra_int_cl);
        %figure;
        %plot_2dgridded(topobase,0.9E19,'','sra','');
    case 'T42'
        %interpolate k1 onto 144x144 grid 
        [X,Y] = meshgrid(1:size(k1,2), 1:size(k1,1));
        [U,V] = meshgrid(1:0.25:size(k1,2), 1:0.25:size(k1,1));
        k1_interp = interp2(X,Y,k1,U,V,'cubic');
        k1_interpsort = k1_interp;
        k1_thr = k1_interpsort;
        thr1 = 0.6;
        k1_thr(k1_thr < thr1) = 0;
        k1_thr(k1_thr >= thr1) = 1;
        %interpolate lat, lon coords to match 
        golon_q = 1.25:0.25:36.75;
        golon_itp = interp1(golonue,golon_q);
        golongrid = zeros(1,145);
        golongrid(:,1) = golonue(:,1);
        golongrid(:,145) = golonue(:,37);
        golongrid(:,2:144) = golon_itp;

        golat_q = 1.25:0.25:36.75;
        golat_itp = interp1(golatue,golat_q);
        golatgrid = zeros(1,145);
        golatgrid(:,1) = golatue(:,1);
        golatgrid(:,145) = golatue(:,37);
        golatgrid(:,2:144) = golat_itp;

        k1_into = zeros(144,144);
        k1_into = k1_thr(2:145,1:144);
        [sra_int, dscd] = make_regrid_2d(golongrid,golatgrid,k1_into',axislonce,axislatce,false);
        sra_int = sra_int';
        sra_int_cl = sra_int;
        %adjust binary threshold;
        thr2 = 0.3;
        sra_int_cl(sra_int_cl > thr2) = 1;
        sra_int_cl(sra_int_cl <= thr2) = 0;
end

disp(['       * Regridding from low-res to high-res requires manual cleaning. Initializing user editing of mask ...']);
%
%
% *********************************************************************** %
% *** USER EDITING OF MASK ********************************************** %
% *********************************************************************** %
%
[jmax imax] = size(sra_int_cl);
grid_fmask = zeros(jmax,imax) + 1;

gm_ex = sra_int_cl;
gm_ex = [gm_ex gm_ex(:,end)];
gm_ex = [gm_ex(1,:); gm_ex];
gfm_ex = grid_fmask;
gfm_ex = [gfm_ex gfm_ex(:,end)];
gfm_ex = [gfm_ex(1,:); gfm_ex];
% Do the topography alteration
fprintf('       * Mask alteration procedure :\n');
fprintf('         (1) left button to turn grid cell to ocean\n');
fprintf('         (2) right button to turn grid cell to land\n');
fprintf('         (3) click outside the grid to finish\n');
cmap = [0.5 0.5 1; 0.5 1 0.5];
%
flag = 0;
while flag == 0
    figure(1); clf
    colormap(cmap);
    pcolor (flipud(gm_ex)); axis image;
    caxis ([-0.5 1.5]); h = colorbar ('horiz');
    %set(h,'XTick',0:1,'XTickLabel',cbartxt);
    set(h,'XTick',0:1);
    title ('Ocean land/sea mask');
    [x,y,button] = ginput(1);
    ix = floor(x); iy = floor(y);
    if ix < 1 | ix > imax | iy < 1 | iy > jmax
        fprintf('       - Out of grid range\n');
        flag = 1;
        fprintf('       * Mask alteration complete\n');
    elseif button == 2
        flag = 1;
        fprintf('       * Mask alteration complete\n');
    else
        loc_j = jmax-iy+2;
        if (gfm_ex(loc_j,ix) == 0.0)
            fprintf('       - WARNING: there is no ocean depth information available at cell (%d, %d)\n',ix,iy);
        end
        if button == 1
            fprintf('         -> Cell at (%d, %d) now ocean\n', ix, iy);
            gm_ex(loc_j,ix) = 0;
        elseif button == 3
            fprintf('         -> Cell at (%d, %d) now land\n', ix, iy);
            gm_ex(loc_j,ix) = 1;
        else
            fprintf('       - Cannot switch mask value of cell at (%d, %d)\n',ix,iy);
        end
    end
    clf
end
close(1)
% return new mask
sra_mask = gm_ex(2:end,1:end-1);

% *********************************************************************** %
% *** EXPORT DATA ******************************************************* %
% *********************************************************************** %
%
% *********************************************************************** %
disp(['       * Regridding complete']);
disp(['       * Cleaning and exporting topo file ...']);

%output file
outname = [fname '_surf_0172.sra'];

sra_out = fix(sra_mask);
header = [172 0 20170927 0 64 32 0 0];

%write to output file

%** uncomment these lines if you don't have WriteMode **
dlmwrite(outname,header,'delimiter','\t', 'precision','%1.0f','roffset',0)
dlmwrite(outname,sra_out,'-append','delimiter','\t')
% *********************************************************************** %

%** uncomment these two lines if you have WriteMode**
%writematrix(header, outname, 'FileType', 'text', 'Delimiter', 'tab')
%writematrix(sra_out, outname, 'FileType', 'text', 'WriteMode','append', 'Delimiter', 'tab')
% *********************************************************************** %

%figure;
%plot_2dgridded(sra_out,0.9E19,'','sra','');
%sra_out = flipud(sra_mask);







