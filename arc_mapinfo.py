import math

from pyresample.geometry import AreaDefinition
from pyresample import save_quicklook, SwathDefinition, image
from pyresample.kd_tree import resample_nearest
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset


class ArcMapInfo:
    def __init__(self, area_id):
        if area_id is None:  # DEFAULT
            area_id = 'npsn_olci'
        self.area_def = self.get_area_definition(area_id)
        self.print_area_def_info()
        self.olci_l2_bands = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75, 753.75, 778.75]


    def create_nc_filegrid(self, ofname, createMask):

        # lons, lats = self.area_def.get_lonlats()

        try:
            OFILE = Dataset(ofname, 'w', format='NETCDF4')
        except PermissionError:
            # print('Permission denied: ', ofname)
            return False

        print(f'Starting file:{ofname}')

        OFILE.createDimension('lat', self.area_def.height)  # ny
        OFILE.createDimension('lon', self.area_def.width)  # nx

        # latitude
        satellite_latitude = OFILE.createVariable('latitude', 'f4', ('lat', 'lon'), fill_value=-999, zlib=True,
                                                  complevel=6)
        # satellite_latitude[:, :] = [np.transpose(lats[:, :])]
        satellite_latitude.units = "degrees_north"
        satellite_latitude.long_name = "latitude"
        satellite_latitude.standard_name = "latitude"

        # longitude
        satellite_longitude = OFILE.createVariable('longitude', 'f4', ('lat', 'lon'), fill_value=-999, zlib=True,
                                                   complevel=6)
        # satellite_longitude[:, :] = [np.transpose(lons[:, :])]
        satellite_longitude.units = "degrees_east"
        satellite_longitude.long_name = "longitude"
        satellite_longitude.standard_name = "longitude"

        # mask
        if createMask:
            min_lat = self.get_lat_min_spherical()
            satellite_mask = OFILE.createVariable('smask', 'i2', ('lat', 'lon'), fill_value=-999, zlib=True,
                                                  complevel=6)
            # satellite_mask[:, :] = [np.transpose(maskma[:, :])]
            satellite_mask.long_name = f'coordinates_mask'
            satellite_mask.standard_name = f'coordinates_mask'
            satellite_mask.description = f'Pixels with latitude lower than: {min_lat} are masked'

        tileY = 5000
        tileX = 5000

        for y in range(0, self.area_def.height, tileY):
            # if self.verbose and (y == 0 or ((y % tileY) == 0)):
            #     print(f'[INFO] Processing line {y}/{ny}')
            print(f'[INFO] Line-> {y}')
            for x in range(0, self.area_def.width, tileX):
                yini = y
                yend = y + tileY
                if yend > self.area_def.height:
                    yend = self.area_def.height
                xini = x
                xend = x + tileX
                if xend > self.area_def.width:
                    xend = self.area_def.width
                width_here = xend - xini
                height_here = yend - yini
                xval = np.zeros((height_here, width_here))
                yval = np.zeros((height_here, width_here))
                for yvalhere in range(height_here):
                    yval[yvalhere, :] = np.linspace(xini, xend, width_here)
                for xvalhere in range(width_here):
                    xval[:, xvalhere] = np.linspace(yini, yend, height_here)
                lons, lats = self.area_def.get_lonlat_from_array_coordinates(xval, yval)
                print(lons.shape)
                satellite_latitude[yini:yend, xini:xend] = [lats[:, :]]
                satellite_longitude[yini:yend, xini:xend] = [lons[:, :]]
                if createMask:
                    mask = np.zeros(lats.shape, dtype=bool)
                    mask[lats < min_lat] = True
                    maskma = ma.MaskedArray(mask, mask=mask)
                    satellite_mask[yini:yend, xini:xend] = [maskma[:, :]]

        OFILE.close()

        return True

    def create_nc_reflectance_file(self, ofname):
        print('Copy file')
        datasetout = self.copy_nc_base(ofname)
        ##create rrs variables
        for wl in self.olci_l2_bands:
            wlstr = str(wl).replace('.', '_')
            bandname = f'RRS{wlstr}'
            var = datasetout.createVariable(bandname, 'f4', ('lat', 'lon'), fill_value=-999, zlib=True,complevel=6)
            var.wavelength = wl

        return datasetout



    def copy_nc_base(self, ofile):
        ifile = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/ArcGridOlci.nc'
        with Dataset(ifile) as src:
            dst = Dataset(ofile, 'w', format='NETCDF4')
            # copy global attributes all at once via dictionary
            dst.setncatts(src.__dict__)

            # copy dimensions
            for name, dimension in src.dimensions.items():
                dst.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))

            # copy all file data except for the excluded
            for name, variable in src.variables.items():
                dst.createVariable(name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)

                dst[name][:] = src[name][:]
        return dst

    def save_quick_look_impl(self, fileout, data):
        save_quicklook(fileout, self.area_def, data)

    def save_quick_look_fgrid(self, fileout, fgrid):
        dataset = Dataset(fgrid)
        data = np.array(dataset.variables['smask'][:, :])

        datavis = np.ma.zeros(data.shape)
        for y in range(datavis.shape[0]):
            datavis[y, :] = y
        datavis[data == 1] = -999
        datavis = np.ma.masked_values(datavis, -999)

        print('NO MASKED: ', datavis.count())
        print('MASKED: ', np.ma.count_masked(datavis))
        self.save_quick_look_impl(fileout, datavis)

    def save_quick_look_fdata(self,fileout,fdata):
        dataset = Dataset(fdata)
        data = np.ma.array(dataset.variables['RRS560'][:, :])
        data = np.ma.masked_values(data,0)
        self.save_quick_look_impl(fileout, data)

    def get_lat_min_spherical(self):
        xmid = math.floor(self.area_def.width / 2)
        lon, lat = self.area_def.get_lonlat_from_array_coordinates(xmid, 0)
        lat_min = np.floor(lat)
        return lat_min

    def print_area_def_info(self):
        print('GRID SIZE:')
        print('Width: ', self.area_def.width)
        print('Height: ', self.area_def.height)

        print('RESOLUTION')
        print('Pixel Size X: ', self.area_def.pixel_size_x)
        print('Pixel Size Y: ', self.area_def.pixel_size_y)

        print('OUTER BOUNDARY CORNERS:')
        ob = self.area_def.outer_boundary_corners
        print('Upper Left: ', np.degrees(ob[0].lat), np.degrees(ob[0].lon))
        print('Upper Rigth: ', np.degrees(ob[1].lat), np.degrees(ob[1].lon))
        print('Lower Right: ', np.degrees(ob[2].lat), np.degrees(ob[2].lon))
        print('Lower Left ', np.degrees(ob[3].lat), np.degrees(ob[3].lon))

        print('MIDDLE POINTS: ')
        xmid = math.floor(self.area_def.width / 2)
        ymid = math.floor(self.area_def.height / 2)
        xend = self.area_def.width - 1
        yend = self.area_def.height - 1
        lon1, lat = self.area_def.get_lonlat_from_array_coordinates(xmid, 0)
        print('Upper Middle: ', lat, lon1)
        lon2, lat = self.area_def.get_lonlat_from_array_coordinates(xend, ymid)
        print('Right Middle: ', lat, lon2)
        lon3, lat = self.area_def.get_lonlat_from_array_coordinates(xmid, yend)
        print('Lower Middle: ', lat, lon3)
        lon4, lat = self.area_def.get_lonlat_from_array_coordinates(0, ymid)
        print('Left Middle: ', lat, lon4)

        print('SPHERIC LIMITS:')
        lons = [np.degrees(ob[0].lon), np.degrees(ob[1].lon), np.degrees(ob[2].lon), np.degrees(ob[3].lon), lon1, lon2,
                lon3, lon4]
        lonmin = np.min(lons)
        lonmax = np.max(lons)
        if lonmin == -180 or lonmax == 180:
            lonmin = -180
            lonmax = 180
        latmin = lat
        lon, latmax = self.area_def.get_lonlat_from_array_coordinates(xmid, ymid)
        print('Lat min: ', latmin)
        print('Lat max: ', latmax)
        print('Lon min: ', lonmin)
        print('Lon max: ', lonmax)

    def get_area_definition(self, area_id):
        area_info = self.get_area_info(area_id)
        proj_id = area_info['proj_id']
        projection = self.get_projection_info(proj_id)
        extent = area_info['extent']
        area_def = AreaDefinition(area_id, area_info['description'], proj_id, projection,
                                  area_info['width'], area_info['height'], extent)
        return area_def

    def get_projection_info(self, proj_id):
        projections = {
            'npsn': '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

        }
        return projections[proj_id]

    def get_area_info(self, area_id):
        # extent: xmin, ymin, xmax, ymax
        areas = {
            'npsn': {
                'description': 'WGS 84 / NSIDC Sea Ice Polar Stereographic North',
                'proj_id': 'npsn',
                'width': 7600,
                'height': 11200,
                'extent': (-3850000, -5350000, 3750000, 5850000)
            },
            'npsn_a': {
                'description': 'WGS 84 / NSIDC Sea Ice Polar Stereographic North',
                'proj_id': 'npsn',
                'width': 6625,
                'height': 6625,
                'extent': (-3314693.24, -3314693.24, 3314693.24, 3314693.24)
            },
            'npsn_olci': {
                'description': 'WGS 84 / NSIDC Sea Ice Polar Stereographic North',
                'proj_id': 'npsn',
                'width': 22098,
                'height': 22098,
                'extent': (-3314693.24, -3314693.24, 3314693.24, 3314693.24)
            },
        }
        return areas[area_id]

    def get_subarea_def(self, lats, lons):

        xcoords, ycoords = self.area_def.get_array_coordinates_from_lonlat(lons, lats)
        print(xcoords.shape)
        print(ycoords.shape)
        xmin = np.floor(np.min(xcoords))
        xmax = np.ceil(np.max(xcoords)) + 1
        ymin = np.floor(np.min(ycoords))
        ymax = np.ceil(np.max(ycoords)) + 1
        print(xmin)
        print(xmax)
        print(ymin)
        print(ymax)
        width = (xmax - xmin)
        height = (ymax - ymin)
        print(width, height)
        projection = self.get_projection_info(self.area_def.proj_id)
        print(projection)

        xvalues = [xmin, xmax, xmax, xmin]
        yvalues = [ymin, ymin, ymax, ymax]
        xcoords, ycoords = self.area_def.get_projection_coordinates_from_array_coordinates(xvalues, yvalues)
        extent = (np.min(xcoords), np.min(ycoords), np.max(xcoords), np.max(ycoords))
        print(extent)

        area_def = AreaDefinition('SubArea', 'SubArea', self.area_def.proj_id, projection,
                                  width, height, extent)

        print(area_def.pixel_size_x, area_def.pixel_size_y)
        print(area_def.width, area_def.height)

        limits = [xmin, xmax, ymin, ymax]
        return limits, area_def

    def make_resample_impl(self, olimage, fileout):
        lats, lons = olimage.get_lat_long_arrays()
        limits, sub_area_def = self.get_subarea_def(lats, lons)
        swath_def = SwathDefinition(lons=lons, lats=lats)

        datasetout = self.create_nc_reflectance_file(fileout)
        for var in datasetout.variables:
            if var.startswith('RRS'):
                var_rrs = datasetout.variables[var]
                wlref = var_rrs.wavelength
                print(var,wlref)
                if wlref>0:
                    rrsarray = olimage.get_reflectance_band_array(wlref)
                    print('doing resample...')
                    result = resample_nearest(swath_def, rrsarray, sub_area_def, radius_of_influence=2000)
                    ymin = int(limits[2])
                    ymax = int(limits[3])
                    xmin = int(limits[0])
                    xmax = int(limits[1])
                    ancho = xmax-xmin
                    alto = ymax-ymin
                    print(ymin,ymax,xmin,xmax)
                    print(alto,ancho)
                    print(result.shape)
                    var_rrs[ymin:ymax,xmin:xmax]=[result[:,:]]

        datasetout.close()
        # refband = olimage.get_reflectance_band_array(560)
        # print('doing resample...')
        # result = resample_nearest(swath_def, refband, sub_area_def, radius_of_influence=2000)
        # print(result.shape)

        # refband = olimage.get_reflectance_band_array(560)
        # print(type(lats), lats.shape)
        # print(type(lons), lons.shape)
        # print(type(refband), refband.shape)
        # print('defining image container...')
        # lats = lats[0:500, 0:500]
        # lons = lons[0:500, 0:500]
        # refband = refband[0:500, 0:500]
        # swath_def = SwathDefinition(lons=lons, lats=lats)
        # # swath_con = image.ImageContainerNearest(refband, swath_def, radius_of_influence=2000)
        # print('doing resample...')
        # # area_con = swath_con.resample(self.area_def)
        #
        # result = resample_nearest(swath_def, refband, self.area_def, radius_of_influence=50000, epsilon=0.5)
        #
        # print(type(result), result.shape)

        # result = kd_tree.resample_nearest(swath_def, data,
        #                                   area_def, radius_of_influence=50000, epsilon=0.5)
