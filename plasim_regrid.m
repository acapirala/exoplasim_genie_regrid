function [] = plasim_regrid(fname)
% make_grid_winds_plasim
%
%   *********************************************************
%   *** RE_GRID WIND PRODUCTS                             ***
%   *********************************************************
%
%
%   wstress is organised :   level 1   tau_x at u point
%                            level 2   tau_x at v point
%                            level 3   tau_y at u point
%                            level 4   tau_y at v point
%
%   wvelocity is organised : level 1   x velocity at grid point
%                            level 2   y velocity at grid point
%
%   wspeed is organised :    level 1   speed at grid point
%
%   str input KEY:
% str(1).nc == par_nc_axes_name
% str(2).nc == par_nc_topo_name
% str(3).nc == par_nc_atmos_name
% str(4).nc == par_nc_mask_name
% str(5).nc == par_nc_coupl_name
% str(6).nc == par_nc_ocean_name
%
%   *********************************************************
%
%   ***********************************************************************
%%

% *********************************************************************** %
% *** INITIALIZE ******************************************************** %
% *********************************************************************** %
%
% create GENIE grid
%imax = par_max_i; %defined in input file
%jmax = par_max_j;
imax = 36;
jmax = 36;
% determine output grid size (remember: [rows columns])
%[jmax, imax] = size(gomask); 
% *********************************************************************** %

% *********************************************************************** %
% *** load input data *************************************************** %
% *********************************************************************** %
%
% NOTE: wind stress in units of: N/m^2
%
% *** load and process data -- wind stress ****************************** %
%
ncid = netcdf.open([fname '.nc'],'nowrite');

%varid  = netcdf.inqVarID(ncid,'puma_u-stress');
%loctaux(:,:) = netcdf.getVar(ncid,varid);
%varid  = netcdf.inqVarID(ncid,'puma_y-stress');
%loctauy(:,:,:) = netcdf.getVar(ncid,varid);

varid  = netcdf.inqVarID(ncid,'monthly_ustress');
loctaux(:,:,:) = netcdf.getVar(ncid,varid);
varid  = netcdf.inqVarID(ncid,'monthly_vstress');
loctauy(:,:,:) = netcdf.getVar(ncid,varid);

%set up input arrays
giwtauu = zeros(size(loctaux(:,:,1))); %64x32 db
giwtauv = zeros(size(loctaux(:,:,1)));
giwvelu = zeros(size(loctaux(:,:,1)));
giwvelv = zeros(size(loctaux(:,:,1)));


% LOOP
for t=1:12 
    giwtauu = giwtauu + loctaux(:,:,t)/12.0; 
    giwtauv = giwtauv + loctauy(:,:,t)/12.0; 
end

% re-orientate
giwtauu = flipud(giwtauu'); %32x64 db
giwtauv = flipud(giwtauv');
%
[go_lonm,go_lone,go_latm,go_late,go_dm,go_de] = make_genie_grid(36,36,16,5000.0,-260.0,true,0);
%[axis_imid,axis_iedge,axis_jmid,axis_jedge,axis_kmid,axis_kedge] = make_genie_grid(imax,jmax,kmax,par_max_D,par_lon_off,opt_equalarea,par_add_Dk)
%1x36, 1x37, 1x36, 1x37, 1x16, 1x17

varid = netcdf.inqVarID(ncid,'latitude');
axislatcm = netcdf.getVar(ncid,varid);
axislatcm = flipud(axislatcm);
axislatce = [axislatcm-180/length(axislatcm)/2.0; 90.0];

varid  = netcdf.inqVarID(ncid,'longitude');
axisloncm = netcdf.getVar(ncid,varid);
axislonce = [axisloncm-360/length(axisloncm)/2; 360.0-360/length(axisloncm)/2];
gilatae = axislatce; 
gilonae = axislonce;

%loc_k1 = load('shnm_gwy.k1'); %38x38
%gk1 = loc_k1(2:end-1,2:end-1); %36x36
%gm(intersect(find(gk1<=89),find(gk1>=1))) = 1;
%go_mask = gm %36x36 binary mask 
%varid     = netcdf.inqVarID(ncid,'latitude');
%axislatcm = netcdf.getVar(ncid,varid); %32x1 single
%axislatcm = flipud(axislatcm);
%axislatce = [axislatcm-180/length(axislatcm)/2.0; 90.0] %convert to
%c grid, 33x1
%varid  = netcdf.inqVarID(ncid,'longitude');
%axisloncm = netcdf.getVar(ncid,varid);
%axislonce = [axisloncm-360/length(axisloncm)/2; 360.0-360/length(axisloncm)/2]; %65x1
%gilatae = axislatce;
%gilonae = axislonce;
%
% *** load data -- wind velocity **************************************** %
%
% LOOP
seas = ['DJF'; 'MAM'; 'JJA'; 'SON'];
% open netCDF file
%ncid = netcdf.open([str(1).path '/' str_month(t,:) str(1).exp '.' str(3).nc '.nc'],'nowrite');
for t=1:4
    % U COMPONENT OF SURFACE AIR WIND
    varid  = netcdf.inqVarID(ncid, [seas(t,:) '_uspeed']);
    locwvelu = netcdf.getVar(ncid,varid);
    % V COMPONENT OF SURFACE AIR WIND
    varid  = netcdf.inqVarID(ncid,[seas(t,:) '_vspeed']);
    locwvelv = netcdf.getVar(ncid,varid);
    % create annual averages
    giwvelu  = giwvelu  + locwvelu(:,:,10)/4.0;
    giwvelv  = giwvelv  + locwvelv(:,:,10)/4.0;
end
% re-orientate
giwvelu = flipud(giwvelu'); %
giwvelv = flipud(giwvelv');
%

% *** set wind stress mask ********************************************** %
%
% copy wind stress mask from ocean mask
%gimaskp = gimask;
% % re-orientate
%?% gimaskp = flipud(gimaskp');
% derive mask values (ocean == 1, land == NaN)
%gimaskp(find(gimaskp~=1.0))  = NaN;
%
% *********************************************************************** %

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
% create GENIE NaN mask
%gm = go_mask;
%gm(find(gm == 0)) = NaN;
%
% *** Put on a GOLDSTEIN grid ******************************************* %
%
% NOTE: don't forget to flip the re-gridded orientation back around again
%
% plot raw wind stress
%if (optplots), plot_2dgridded(flipud(giwtauu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_u.IN'],['wind stress in -- u']); end
%if (optplots), plot_2dgridded(flipud(giwtauv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_v.IN'],['wind stress in -- v']); end
% 
% -> wind stress @ u point
% NOTE: GENIE u-grid: (golatue,golonue) == jmax   x imax+1
% apply p-grid mask to input wind stress grid
% replace NaNs with zeros
% apply GENIE mask to output wind stress
% 
[gowtauuu,gofwtauuu] = make_regrid_2d(gilonae,gilatae,giwtauu',golonue,golatue,false);
gowtauuu(find(isnan(gowtauuu))) = 0.0;
gowtauuu = gowtauuu';

%if (optplots), plot_2dgridded(gm.*gowtauuu,999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_xATu.out'],['wind stress out -- x @ u']); end
% v
%fprintf('       - Regridding wind stress (y) to GOLDSTEIN u-grid\n');
[gowtauvu,gofwtauvu] = make_regrid_2d(gilonae,gilatae,giwtauv',golonue,golatue,false);
gowtauvu(find(isnan(gowtauvu))) = 0.0;
gowtauvu = gowtauvu';

%if (optplots), plot_2dgridded(flipud(gm.*gowtauvu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_yATu.out'],['wind stress out -- y @ u']); end
% 
% -> wind stress @ v point
% NOTE: GENIE v-grid: (golatve,golonve) == jmax+1 x imax
% apply p-grid mask to input wind stress grid
% replace NaNs with zeros
%fprintf('       - Regridding wind stress (x) to GOLDSTEIN v-grid\n');

[gowtauuv,gofwtauuv] = make_regrid_2d(gilonae,gilatae,giwtauu',golonve,golatve,false);
gowtauuv(find(isnan(gowtauuv))) = 0.0;
gowtauuv = gowtauuv';

%if (optplots), plot_2dgridded(flipud(gm.*gowtauuv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_xATv.out'],['wind stress out -- x @ v']); end
% v
%fprintf('       - Regridding wind stress (y) to GOLDSTEIN v-grid\n');

[gowtauvv,gofwtauvv] = make_regrid_2d(gilonae,gilatae,giwtauv',golonve,golatve,false);
gowtauvv(find(isnan(gowtauvv))) = 0.0;
gowtauvv = gowtauvv';

%if (optplots), plot_2dgridded(flipud(gm.*gowtauvv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_yATv.out'],['wind stress out -- y @ v']); end
% 
% -> wind velocity
% NOTE: GENIE c-grid: (golate,golone) == jmax+1 x imax+1
% NOTE: DO NOT apply a mask (as the output is used globally by the EMBM)
% NOTE: input grid is the FOAM atmosphere grid
% replace NaNs with zeros
% u
%fprintf('       - Regridding wind velocity (x) to GOLDSTEIN c-grid\n');
%if (optplots), plot_2dgridded(flipud(giwvelu),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_x.IN'],['wind velocity in -- x']); end

[gowvelu,gofwvelu] = make_regrid_2d(gilonae,gilatae,giwvelu',go_lone,go_late,false);
gowvelu(find(isnan(gowvelu))) = 0.0;
gowvelu = gowvelu'; 
gofwvelu = gofwvelu';

%if (optplots), plot_2dgridded(flipud(gowvelu),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_x.OUT'],['wind velocity out -- x']); end
% 
%fprintf('       - Regridding wind velocity (y) to GOLDSTEIN c-grid\n');
%if (optplots), plot_2dgridded(flipud(giwvelv),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_y.IN'],['wind velocity in -- y']); end
[gowvelv,gofwvelv] = make_regrid_2d(gilonae,gilatae,giwvelv',go_lone,go_late,false);
gowvelv(find(isnan(gowvelv))) = 0.0;
gowvelv = gowvelv';
gofwvelv = gofwvelv';

%if (optplots), plot_2dgridded(flipud(gowvelv),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_y.OUT'],['wind velocity out -- y']); end

%
% *** Copy to output arrays ********************************************* %
%
% -> wind stress
gowtauuu = flipud(gowtauuu); %g_taux_u
gowtauuv = flipud(gowtauuv); %g_taux_v
gowtauvu = flipud(gowtauvu); %g_tauy_u
gowtauvv = flipud(gowtauvv); %g_tauy_v

% -> wind velocity
gowvelu = flipud(gowvelu);
gowvelv = flipud(gowvelv);

%
% *********************************************************************** %

% *********************************************************************** %
% *** SAVE DATA ********************************************************* %
% *********************************************************************** %
%
% Save regridded data to file
% Taux at u point (g_taux_u == gowtauuu)

simname = fname;

%
outname = [simname '.taux_u.dat'];
c = gowtauuu; b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');

%fprintf('       - Written tau u (u point) data to %s\n',outname);
% Taux at v point (g_taux_v == gowtauuv)
outname = [simname '.taux_v.dat'];
c = gowtauuv; b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');

%fprintf('       - Written tau u (v point) data to %s\n',outname);
% Tauy at u point (g_tauy_u == gowtauvu)
outname = [simname '.tauy_u.dat'];
c = gowtauvu; b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
%fprintf('       - Written tau v (u point) data to %s\n',outname);

% Tauy at v point (g_tauy_v == gowtauvv)
outname = [simname '.tauy_v.dat'];
c = gowtauvv; b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
%fprintf('       - Written tau v (v point) data to %s\n',outname);

% U wind velocity 
outname = [simname '.wvelx.dat'];
c = gowvelu; b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
%fprintf('       - Written u wind speed data to %s\n',outname);
% V wind velocity
outname = [simname '.wvely.dat'];
c = gowvelv; b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
%fprintf('       - Written v wind speed data to %s\n',outname);

%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%
% *********************************************************************** %
