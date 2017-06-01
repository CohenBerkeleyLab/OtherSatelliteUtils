function [ Data ] = read_hdf5( filename, geofields, datafields )
%READ_HDF5 Reads geographic and data variables from HDF5 file, handles fills, scaling, and offset
%   DATA = READ_HDF5( FILENAME, GEOFIELDS, DATAFIELDS ) Reads the HDF5 file
%   specified by FILENAME, reading the specified GEOFIELDS and DATAFIELDS,
%   given as cell arrays. Returns as structure DATA with each dataset as a
%   field. All data has fill values removed (tolerance 0.1%), then scale
%   factor and offset applied, in that order.

give_offset_warn = true;

hi = h5info(filename);

for a=1:numel(geofields)
    Data.(geofields{a}) = read_dataset(2, geofields{a});
end
for a=1:numel(datafields)
    Data.(datafields{a}) = read_dataset(1, datafields{a});
end

    function vals = read_dataset(grp, dsetname)
        fillval = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,grp,dsetname), '_FillValue'));
        scalefac = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,grp,dsetname), 'ScaleFactor'));
        offset = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,grp,dsetname), 'Offset'));
        if offset ~= 0 && give_offset_warn
            warning('Offset in %s not 0, offset is applied after scaling but this order is unconfirmed.\nFurther offset warning will be suppressed',dsetname)
            give_offset_warn = false;
        end
        
        vals = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,grp,dsetname)));
        fills = abs((vals - fillval)/fillval) < 1e-3;
        vals(fills) = nan;
        vals = (vals * scalefac) + offset;
    end

end

