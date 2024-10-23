function [] = topo_plasimtogenie(fname)
% topo_plasimtogenie
%
%   *********************************************************
%   *** RE_GRID PLASIM TOPOGRAPHY TO GENIE ***
%   *********************************************************
%
% (1) read in '_surf_0172.sra' file
% (2) extract core data plus imax, jmax
% (3) extract topo
% (4) return .k1 topo file
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
% create GENIE grid
%imax = par_max_i; %defined in input file
%jmax = par_max_j;
imax = 36;
jmax = 36;
% determine output grid size (remember: [rows columns])
%[jmax, imax] = size(gomask); 
[go_lonm,go_lone,go_latm,go_late,go_dm,go_de] = make_genie_grid(36,36,16,5000.0,-260.0,true,0);
%[axis_imid,axis_iedge,axis_jmid,axis_jedge,axis_kmid,axis_kedge] = make_genie_grid(imax,jmax,kmax,par_max_D,par_lon_off,opt_equalarea,par_add_Dk)
%1x36, 1x37, 1x36, 1x37, 1x16, 1x17

% *********************************************************************** %
% *** load input data *************************************************** %
% *********************************************************************** %

disp(['       * Loading input data ...']);
sra = readmatrix([fname '_surf_0172.sra'], 'FileType', 'text', 'NumHeaderLines', 1);

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

% *********************************************************************** %
% *** RE-GRID *********************************************************** %
% *********************************************************************** %
%
% *** set up GOLDSTEIN grid ********************************************* %
%
% Need to determine wind stress at u and v points of the Arakawa C
% grid in GENIE ...
%
%   ----v----
%   |       |
%   |   c   u
%   |       |
%   ---------
% 
% remember: [rows columns] == [j i]
% grid boundaries are as follows:
% GENIE c-grid: (golate,golone)   == jmax+1 x imax+1
% GENIE u-grid: (golatue,golonue) == jmax   x imax+1
% GENIE v-grid: (golatve,golonve) == jmax+1 x imax
%
% create GENIE u and v edge axes
% latitude
golatue = go_late;
golatve = [go_latm 90.0]; 
% longitude
golonue = [go_lonm go_lonm(end)+360.0/imax];
golonve = go_lone;

% *** Put on a GOLDSTEIN grid ******************************************* %


disp(['       * Regridding topography ...']);
[topobase,wtopobase] = make_regrid_2d(axislonce,axislatce,sra',golonue,golatue,false);
topobase = topobase';
figure; imagesc(topobase); colorbar

%figure;
%plot_2dgridded(topobase,0.9E19,'','sra','');

%make binary
topobase(topobase >= 0.3) = 1;
topobase(topobase < 0.3) = 0;
%topobase = flipud(topobase);

disp(['       * Initializing user editing of mask ...']);

% *********************************************************************** %
% *** USER EDITING OF MASK ********************************************** %
% *********************************************************************** %
%
[jmax imax] = size(topobase);
grid_fmask = zeros(jmax,imax) + 1;

gm_ex = topobase;
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
k1_mask = gm_ex(2:end,1:end-1);
figure; imagesc(k1_mask); colorbar
%k1_mask = flipud(k1_mask);

% *********************************************************************** %
% *** EXPORT DATA ******************************************************* %
% *********************************************************************** %

% *** Clean and export .k1 file ***************************************** %

disp(['       * Regridding complete']);

k1_mask(k1_mask > 0) = 91;
k1_mask(k1_mask == 0) = 1;

%correct values
k1 = zeros(jmax+2,imax+2);
k1(1,:) = 17;
k1(38,:) = 18;
k1(2:37,2:37) = k1_mask;
k1(2:37,1) = k1_mask(:,end);
k1(2:37,38) = k1_mask(:,1);

str = struct('gcm', {}, 'exp', {}, 'path', {}, 'dir', {}, 'nc', {});
str = setfield(str, {1}, 'gcm', 'k1');

k1rn = make_grid_runoff_rnd(k1,str,false);
k1rn(k1rn == 88) = 17;
k1rn(k1rn == 89) = 18;
k1out = fix(k1rn);
figure; imagesc(k1out); colorbar

disp(['       * Exporting .k1 file ...']);

%figure;
%plot_2dgridded(k1out,0.9E19,'','sra','');

%output file
outname = [fname '.k1'];
fid = fopen(outname, 'w');
for i = 1:38
    for j = 1:38
        if k1out(i,j) < 10
            fprintf(fid,'  %d',k1out(i,j));
        else
            fprintf(fid,' %d',k1out(i,j));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);



