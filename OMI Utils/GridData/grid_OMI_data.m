% grid_GOME_data - grids GOME data to 0.05 degree grid cells

date_start = '9/01/2014';
date_end = '9/24/2014';

latbdy = [25 50];
lonbdy = [-125 -65];
reslat = 0.05;
reslon = 0.05;

load_dir = '/Volumes/share-sat/SAT/OMI/Bare_SP_Files/';
save_dir = '/Volumes/share-sat/SAT/OMI/Gridded_SP_Files/';

load_prefix = 'OMI_SP_';
save_prefix = 'OMI_SP_griddedKingFire_';

DEBUG_LEVEL = 2;

datenums = datenum(date_start):datenum(date_end);
for d=1:numel(datenums);
    curr_date = datestr(datenums(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
    fname = sprintf('%s%s%s%s.mat',load_prefix,year,month,day);
    load_path = fullfile(load_dir,fname);
    if ~exist(load_path,'file')
        if DEBUG_LEVEL > 0; fprintf('%s not found in %s\n',fname, load_path); end
    else
        if DEBUG_LEVEL > 0; fprintf('%s found\n',fname); end
        load(load_path);
        nd = numel(Data); % Data is loaded from the GOME mat file
        E=0;
        for D=1:nd
            if DEBUG_LEVEL > 0; fprintf('    Swath %d\n',D); end
            inbdy = Data(D).Latitude(:) > min(latbdy) & Data(D).Latitude(:) < max(latbdy)...
                & Data(D).Longitude(:) > min(lonbdy) & Data(D).Longitude(:) < max(lonbdy);
            if sum(inbdy) == 0; 
                if DEBUG_LEVEL > 1; fprintf('\t No point in boundaries\n'); end
                continue; 
            end
            E=E+1;
            OMI(E) = add2grid_OMI(Data(D),reslat,reslon,lonbdy,latbdy);
        end
        
        save_name = sprintf('%s%s%s%s.mat',save_prefix,year,month,day);
        save_path = fullfile(save_dir,save_name);
        save(save_path,'Data','OMI','-v7.3');
        
        clear('Data'); clear('OMI');
    end
end