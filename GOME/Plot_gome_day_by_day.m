date_start = '6/21/2012';
date_end = '6/21/2012';

latbdy = [25 50];
lonbdy = [-125 -65];
reslat = 0.05;
reslon = 0.05;

load_dir = '/Volumes/share-sat/SAT/GOME-2/Matlab Files/Gridded';

load_prefix = 'GOME2_gridded_';

screensize = get(0,'ScreenSize');

datenums = datenum(date_start):datenum(date_end);
for d=1:numel(datenums);
    curr_date = datestr(datenums(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
        
    fname = sprintf('%s%s%s%s.mat',load_prefix,year,month,day);
    load_path = fullfile(load_dir,year,fname);
    if ~exist(load_path,'file')
        fprintf('No data available for %s\n',curr_date);
    else
        load(load_path);
        for s=1:numel(OMI)
            l(s) = figure;
            pos = get(l(s),'Position');
            set(l(s),'Position',[(s-1)*pos(3), pos(2:4)]);
            pcolor(OMI(s).SurfacePressure); shading flat;
            colorbar; %caxis([0 10]);
        end
        pause;
        close(l);
    end
end