function [  ] = load_and_grid_BEHR( start_date, end_date )
%LOAD_AND_GRID_OMNO2 Read in and grid level 2 BEHR data.
%   This function is a special case for gridding BEHR data and should not
%   be used normally. The point of this is to grid BEHR pixels in exactly
%   the same way as worldwide OMNO2 and DOMINO data for an apples-to-apples
%   comparison.

E = JLLErrors;
DEBUG_LEVEL = 2;

global onCluster
if isempty(onCluster)
    onCluster = false;
end

if onCluster
    behr_path = '/global/home/users/laughner/myscratch/SAT/BEHR/BEHR_Files';
    save_path = '/global/home/users/laughner/myscratch/MATLAB/Data/OMI/BEHR/0.25x0.25-avg';
else
    behr_path = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_2014';
    save_path = '/Volumes/share2/USERS/LaughnerJ/DOMINO-OMNO2_comparision/BEHR/0.25x0.25';
end

% Existing files must be more recent than this to be considered complete
min_date = '09-Jan-2016 12:00:00';

xres = 0.25; yres = 0.25;
[loncorn, latcorn] = grid_corners(xres, yres);
loncorn = loncorn';
lon = loncorn(1:end-1,1:end-1) + xres/2;
latcorn = latcorn';
lat = latcorn(1:end-1,1:end-1) + yres/2;

parfor d=datenum(start_date):datenum(end_date)
    t = getCurrentTask();
%    t.ID = 0;
    if DEBUG_LEVEL > 0; fprintf('W%d: Loading/gridding data for %s\n',t.ID,datestr(d)); end

    % Skip this day if the file has already been created and is new enough
    save_name = sprintf('OMI_BEHRL2b_%04d%02d%02d.mat',year(d),month(d),day(d));
    if exist(fullfile(save_path, save_name),'file')
        F = dir(fullfile(save_path, save_name));
        if F.datenum > datenum(min_date)
            fprintf('W%d: Skipping %s, already complete\n',t.ID,save_name);
            continue
        end
    end

    if DEBUG_LEVEL > 2; fprintf('W%d: Findings files for %s\n',t.ID,datestr(d)); end
    % BEHR data should all be in one folder
    file_pat = sprintf('OMI_BEHR_%04d%02d%02d.mat',year(d),month(d),day(d));
    F = dir(fullfile(behr_path, file_pat));

    
    if DEBUG_LEVEL > 2; fprintf('W%d: Initializing GC structure\n', t.ID); end
    omi_mat = nan(size(loncorn)-1);
    omi_cell = cell(size(loncorn)-1);
    GC = struct('BEHRAMFTrop', omi_mat, 'CloudFraction', omi_mat, 'CloudPressure', omi_mat, 'MODISCloud', omi_mat,...
        'MODISAlbedo', omi_mat, 'GLOBETerpres', omi_mat, 'BEHRColumnAmountNO2Trop', omi_mat,...
        'Areaweight', omi_mat, 'vcdQualityFlags', {omi_cell}, 'XTrackQualityFlags', {omi_cell});

    
    %lonmin = -180;  lonmax = 180;
    %latmin = -90;   latmax = 90;
    %reslat = 2; reslon = 2.5;
    D = load(fullfile(behr_path, F.name));
    GC = repmat(GC,1,numel(D.Data));
    for s=1:numel(D.Data)
        if DEBUG_LEVEL > 1; fprintf('\t W%d: Handling swath %d of %d\n',t.ID,s,numel(D.Data)); end
        %OMI(s) = add2grid_DOMINO(Data(s),OMI(s),reslat,reslon,[lonmin, lonmax],[latmin, latmax]);
        GC(s) = grid_to_gc(D.Data(s), GC(s),loncorn,latcorn,DEBUG_LEVEL);
    end

    % Average over the days' files to get a single day's data
    if DEBUG_LEVEL > 2; fprintf('W%d: Averaging GC structure\n',t.ID); end
    GC_avg = struct(GC(1));
    GC_avg = rmfield(GC_avg,'Areaweight');
    fns = fieldnames(GC_avg);
    aw = cat(3,GC.Areaweight);
    if DEBUG_LEVEL > 3; fprintf('W%d: Percent nans areaweight = %.4f\n', t.ID, sum(isnan(aw(:)))/numel(aw)*100); end
    total_aw = nansum2(aw,3);
    if DEBUG_LEVEL > 3; fprintf('W%d: Percent nans total areaweight = %.4f\n', t.ID, sum(isnan(total_aw(:)))/numel(total_aw)*100); end
    for f=1:numel(fns)
        if DEBUG_LEVEL > 3; fprintf('\tW%d: Avg. field %s\n',t.ID,fns{f}); end
        if ~iscell(GC_avg.(fns{f}))
            GC_avg.(fns{f}) = nansum2(cat(3, GC.(fns{f})) .* aw, 3) ./ total_aw;
        else
            GC_avg = rmfield(GC_avg, fns{f});
        end
    end

    GC_avg.Longitude = lon;
    GC_avg.Latitude = lat;
  

    save_name = sprintf('OMI_BEHRL2b_%04d%02d%02d.mat',year(d),month(d),day(d));
    if DEBUG_LEVEL > 0; fprintf('W%d: Saving file %s\n', t.ID, fullfile(save_path, save_name)); end
    saveData(fullfile(save_path,save_name),GC_avg);

end


end

function saveData(filename,GC_avg)
    save(filename,'GC_avg','-v7.3')
end

function [loncorn, latcorn] = grid_corners(lonres, latres)
    loncorn = -125:lonres:-65;
    latcorn = 25:latres:50;
    % Add a little give on the last corner since we test >= the lower index
    % corner but < the higher index one; this will prevent issues with pixels
    % at exactly 180 E or 90 N
    loncorn(end)=loncorn(end) + 0.1;
    latcorn(end)=latcorn(end) + 0.1;
    [loncorn, latcorn] = meshgrid(loncorn, latcorn);
end

function GC = grid_to_gc(Data, GC, gloncorn, glatcorn, DEBUG_LEVEL)
%[gloncorn, glatcorn] = geos_chem_corners();

sz = size(gloncorn)-1;

Data = calc_pixel_area(Data);
TotalAreaweight = zeros(size(GC.BEHRColumnAmountNO2Trop));

% Filter for clouds, row anomaly and albedo (DOMINO recommends only using
% albedo < 0.3 as snow/ice can interfere in cloud retrieval)
row_anom = Data.XTrackQualityFlags > 0;
vcd_qual = Data.vcdQualityFlags > 0;
clds = Data.CloudFraction > 0.3;
alb = Data.MODISAlbedo > 0.3;
negvcds = Data.BEHRColumnAmountNO2Trop < 0;
extreme_vcd = Data.BEHRColumnAmountNO2Trop > 1e17;

fns = fieldnames(GC);
for c=1:numel(fns)
    if DEBUG_LEVEL > 2; fprintf('File %s: Data.%s is %.4f %% nans before fill removal\n',Data.Filename(19:32), fns{c}, sum(isnan(Data.(fns{c})(:)))/numel(Data.(fns{c}))); end
    Data.(fns{c})(row_anom | vcd_qual | clds | alb | negvcds | extreme_vcd) = nan;
    if DEBUG_LEVEL > 2; fprintf('File %s: Data.%s is %.4f %% nans after fill removal\n',Data.Filename(19:32), fns{c}, sum(isnan(Data.(fns{c})(:)))/numel(Data.(fns{c}))); end
end

glon_vec = gloncorn(:,1);
glat_vec = glatcorn(1,:);

fprintf('Starting loop\n');
tval = tic;
for a=1:numel(Data.Longitude)
    %if mod(a,100) == 1
     %   fprintf('Pixel %d: time elapsed = %f\n',a,toc(tval));
    %end
    if Data.Longitude(a) < -200 || Data.Latitude(a) < -200 || isnan(Data.Longitude(a)) || isnan(Data.Latitude(a))
        % Fill values are ~1x10^-30, this will skip any pixels with
        % fill values for coordinates
        if DEBUG_LEVEL > 3; fprintf('Skipping pixel in %s b/c lat/lon is fill\n',Data.Filename); end
        continue
    end
    xx = Data.Longitude(a) >= glon_vec(1:end-1) & Data.Longitude(a) < glon_vec(2:end);
    yy = Data.Latitude(a) >= glat_vec(1:end-1) & Data.Latitude(a) < glat_vec(2:end);
    if sum(xx) > 1 || sum(yy) > 1
        error('load_and_grid_OMNO2:pixel_assignment_error', 'Pixel %d in %s fell into more than 1 grid cell', a, Data.Date);
    elseif sum(xx) < 1 || sum(yy) < 1 
        % Because we only really care about BEHR (for now) between -125 and
        % -65 W and 25 to 50 N, but make an effort to keep rows in OMI
        % intact, there will be pixels outside the domain that we need to
        % skip.
        if DEBUG_LEVEL > 2; fprintf('Pixel %d in %s could not be assigned to a grid cell', a, Data.Date); end
        continue
    end
    TotalAreaweight(xx,yy) = nansum2([TotalAreaweight(xx,yy), Data.Areaweight(a)]);
    for c=1:numel(fns)
        if iscell(GC.(fns{c}))
           GC.(fns{c}){xx,yy} = cat(1,GC.(fns{c}){xx,yy},Data.(fns{c})(a));
        else
            % The areaweight will be divided out at the end of the loop
            GC.(fns{c})(xx,yy) = nansum2([GC.(fns{c})(xx,yy), Data.(fns{c})(a) .* Data.Areaweight(a)]);
        end 
    end
end

for c=1:numel(fns)
    if ~iscell(GC.(fns{c}))
        GC.(fns{c}) = GC.(fns{c}) ./ TotalAreaweight;
    end
end
fprintf('\t Time for one file = %f\n',toc(tval));

%for a=1:sz(1)
%    if mod(a,25) == 1
%        fprintf('Now on %d of %d\n',a,sz(1));
%    end
%    tval_a = tic;
%    for b=1:sz(2)
%        tval=tic;
%        xx = Data.Longitude >= gloncorn(a,b) & Data.Longitude < gloncorn(a+1,b+1);
%        yy = Data.Latitude >= glatcorn(a,b) & Data.Latitude < glatcorn(a+1,b+1);
%        fprintf('\t Timer: xx yy = %f\n',toc(tval));
%        for c=1:numel(fns)
%            if sum(xx(:)&yy(:)) > 0
%                if iscell(GC.(fns{c}))
%                    tval=tic;
%                    GC.(fns{c}){a,b} = Data.(fns{c})(xx&yy);
%                    fprintf('\t Timer: binning cell arrays = %f\n',toc(tval));
%                else
%                    % Calculate an area-weighted mean
%                    tval=tic;
%                    GC.(fns{c})(a,b) = nansum2(Data.(fns{c})(xx&yy) .* Data.Areaweight(xx&yy)) ./ nansum2(Data.Areaweight(xx&yy));
%                    fprintf('\t Timer: averaging = %f\n',toc(tval));
%                end
%            end
%        end
%    end
%    fprintf('  Timer: one loop over b = %f\n',toc(tval_a));
%end

end

function Data = calc_pixel_area(Data)
Lon1 = squeeze(Data.Loncorn(1,:,:));
Lon2 = squeeze(Data.Loncorn(2,:,:));
Lon3 = squeeze(Data.Loncorn(3,:,:));
Lon4 = squeeze(Data.Loncorn(4,:,:));

Lat1 = squeeze(Data.Latcorn(1,:,:));
Lat2 = squeeze(Data.Latcorn(2,:,:));
Lat3 = squeeze(Data.Latcorn(3,:,:));
Lat4 = squeeze(Data.Latcorn(4,:,:));

Data.Areaweight = nan(size(Data.Longitude));

for x=1:size(Lon1,1)
    for y=1:size(Lon1,2)
        pixelarea = (m_lldist([Lon1(x,y)-180 Lon2(x,y)-180],[Lat1(x,y) Lat2(x,y)]))*(m_lldist([Lon1(x,y)-180, Lon4(x,y)-180],[Lat1(x,y), Lat4(x,y)]));
        Data.Areaweight(x,y) = 1/pixelarea;
    end
end
end

function Data = read_fields(Data, hi, gnum, field_names)
% Reads in specified fields from the .he5 file given. hi should be the
% result of h5info.
for f=1:size(field_names,1)
    fillval = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}), 'MissingValue'));
    offset = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}), 'Offset'));
    scalefactor = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1}), 'ScaleFactor'));
    
    
    Data.(field_names{f,2}) = double(Data.(field_names{f,2}));
    Data.(field_names{f,2}) = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,gnum ,field_names{f,1})));
    Data.(field_names{f,2})(Data.(field_names{f,2}) == fillval) = nan;
    Data.(field_names{f,2}) = (Data.(field_names{f,2}) - offset) * scalefactor;
end
end

