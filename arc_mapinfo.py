import math, os
import configparser
import shutil

import matplotlib.pyplot as plt
from pyresample.geometry import AreaDefinition
from pyresample import save_quicklook, SwathDefinition, utils, image

from pyresample.kd_tree import resample_nearest
from pyresample.kd_tree import get_sample_from_neighbour_info, get_neighbour_info
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import warnings

warnings.resetwarnings()
warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', DeprecationWarning)
warnings.simplefilter('ignore', ResourceWarning)


class ArcMapInfo:
    def __init__(self, arc_options, verbose):
        # DEFAULTS
        area_id = 'polar_stereographic'
        ifile_base_default = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/GRID_FILES/ArcGrid_65_90_300mNOGEO.nc'
        # ifile_base_default = '/store/COP2-OC-TAC/arc/code/ArcGrid_65_90_300m.nc'
        self.ifile_base = ifile_base_default

        if arc_options is not None:
            area_id = arc_options.get_value_param('GENERAL', 'area_id', 'polar_stereographic', 'str')
            self.ifile_base = arc_options.get_value_param('GENERAL', 'grid_base', ifile_base_default, 'str')

        self.area_def = self.get_area_definition(area_id)
        self.verbose = verbose
        if self.verbose:
            self.print_area_def_info()
            print('---------------------------------------------------------------------------------------------')

        self.olci_l2_bands = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75, 753.75, 778.75]
        self.climatology_path = None

    # option config_file for adding globat attributes
    def create_nc_filegrid(self, ofname, createMask, createLatLong):

        try:
            OFILE = Dataset(ofname, 'w', format='NETCDF4')
        except PermissionError:
            return False

        print(f'Starting file:{ofname}')

        OFILE.createDimension('y', self.area_def.height)  # ny
        OFILE.createDimension('x', self.area_def.width)  # nx
        OFILE.createDimension('time', 1)

        # variable stereographic
        stereographic = OFILE.createVariable('stereographic', 'i4')
        stereographic.grid_mapping_name = "polar_stereographic"
        stereographic.latitude_of_projection_origin = 90.0
        stereographic.longitude_of_projection_origin = -45.0
        stereographic.standard_parallel = 70.0
        stereographic.straight_vertical_longitude_from_pole = -45.0
        stereographic.false_easting = 0.0
        stereographic.false_northing = 0.0
        stereographic.earth_radius = 6378273
        stereographic.proj4 = '+proj=stere +lon_0=-45 +lat_0=90 +k=1 +R=6378273 +no_defs'

        # y
        satellite_y = OFILE.createVariable('y', 'f4', ('y',), zlib=True, shuffle=True, complevel=4)
        satellite_y.units = "metres"
        satellite_y.standard_name = "projection_y_coordinate"
        satellite_y.axis = "Y"
        ymin = self.area_def.area_extent[3]
        ymax = self.area_def.area_extent[1]
        array = np.linspace(ymin, ymax, self.area_def.height)
        satellite_y[:] = [array[:]]

        # x
        satellite_x = OFILE.createVariable('x', 'f4', ('x',), zlib=True, shuffle=True, complevel=4)
        satellite_x.units = "metres"
        satellite_x.standard_name = "projection_x_coordinate"
        satellite_x.axis = "X"
        xmin = self.area_def.area_extent[0]
        xmax = self.area_def.area_extent[2]
        array = np.linspace(xmin, xmax, self.area_def.width)
        satellite_x[:] = [array[:]]

        # time
        time = OFILE.createVariable('time', 'i4', ('time',), zlib=True, shuffle=True, complevel=4)
        time.long_name = 'reference time'
        time.standard_name = 'time'
        time.axis = 'T'
        time.calendar = 'Gregorian'
        time.units = 'seconds since 1981-01-01 00:00:00'
        time[0] = 0

        if createLatLong:
            # latitude
            satellite_latitude = OFILE.createVariable('lat', 'f4', ('y', 'x'), zlib=True, shuffle=True,
                                                      complevel=4, least_significant_digit=3)
            satellite_latitude.units = "degrees_north"
            satellite_latitude.standard_name = "latitude"
            satellite_latitude.long_name = "latitude"
            # satellite_latitude.valid_min = self.get_lat_min()
            # satellite_latitude.valid_max = 90.0
            satellite_latitude.comment = "Spherical latidude from 65 to 90 degrees north"

            # longitude
            satellite_longitude = OFILE.createVariable('lon', 'f4', ('y', 'x'), zlib=True, shuffle=True,
                                                       complevel=4, least_significant_digit=3)
            satellite_longitude.units = "degrees_east"
            satellite_longitude.standard_name = "longitude"
            satellite_longitude.long_name = "longitude"
            # satellite_longitude.valid_min = -180.0
            # satellite_longitude.valid_max = 180.0
            satellite_longitude.comment = "Spherical longitude from -180 to 180 degrees east"

        # mask
        if createMask:
            min_lat = self.get_lat_min_spherical()
            satellite_mask = OFILE.createVariable('SENSORMASK', 'i2', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                                  shuffle=True, complevel=4)
            satellite_mask.long_name = f'Sensor Mask'
            satellite_mask.comment = 'Each SENSORMASK pixel is the sum of all available sensor values. Pixels with ' \
                                     'latitude lower than 65 degrees are masked '
            satellite_mask.grid_mapping = "stereographic"
            satellite_mask.coordinates = "longitude latitude"

        tileY = 5000
        tileX = 5000

        for y in range(0, self.area_def.height, tileY):
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
                    yval[yvalhere, :] = yini + yvalhere
                for xvalhere in range(width_here):
                    xval[:, xvalhere] = xini + xvalhere
                lons, lats = self.area_def.get_lonlat_from_array_coordinates(xval, yval)

                if createLatLong:
                    satellite_latitude[yini:yend, xini:xend] = [lats[:, :]]
                    satellite_longitude[yini:yend, xini:xend] = [lons[:, :]]

                if createMask:
                    mask = np.zeros(lats.shape, dtype=int)
                    mask[lats < min_lat] = -999
                    # mask[np.where((lats >= 70) & (lats <= 75))] = 10
                    # mask[np.where((lons >= 10) & (lons <= 20))] = mask[np.where((lons >= 10) & (lons <= 20))] + 10
                    satellite_mask[0, yini:yend, xini:xend] = [mask[:, :]]

        OFILE.close()

        return True

    def create_nc_file_resample_base(self, olimage, fileout, arc_opt):
        if self.verbose:
            print('--------------------')
            print(f'[INFO] Starting creating resampling file base using: {olimage.path_source}')
        if self.verbose:
            print(f'[INFO] Creating output file from file grid...')
        all_band_names = self.get_all_data_bands_names(arc_opt)
        for name in all_band_names:
            print(f'[INFO] Band name -> {name}')
        datasetout = self.create_nc_file_resampled(fileout, olimage, arc_opt)

        if datasetout is None:
            print(f'[ERROR] Output file {fileout} could not be created')
        datasetout.close()

    def create_nc_file_resampled(self, ofname, olimage, arc_opt):
        section = 'RESAMPLE'
        file_base = arc_opt.get_value_param(section, 'file_base', None, 'str')
        if file_base is not None:
            if os.path.exists(file_base):
                if self.verbose:
                    print(f'[INFO] Check variables in file base: {file_base}')
                check_variables = True
                dataset = Dataset(file_base)
                all_band_names = self.get_all_data_bands_names(arc_opt)
                for name in all_band_names:
                    if not name in dataset.variables:
                        if self.verbose:
                            print(f'[INFO] Variable: {name} is not present in file base.')
                        check_variables = False
                dataset.close()
                if check_variables:
                    if self.verbose:
                        print(f'[INFO] All the variables are present. Working with file base.')
                    self.ifile_base = file_base
                    datasetout = self.copy_nc_base(ofname)
                    return datasetout

        datasetout = self.copy_nc_base(ofname)

        if datasetout is None:
            return datasetout

        ##create rrs variables
        rrs_bands = arc_opt.get_value_param(section, 'rrs_bands', self.olci_l2_bands, 'floatlist')
        for wl in rrs_bands:
            wlstr = str(wl).replace('.', '_')
            bandname = f'RRS{wlstr}'
            if bandname.endswith('_0'):
                bandname = bandname[:6]
            if self.verbose:
                print(f'[INFO] Creating band: {bandname}')
            var = datasetout.createVariable(bandname, 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
            var[:] = -999
            var.wavelength = wl
            var.standard_name = f'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air'
            var.long_name = f'Sea surface reflectance defined as the ratio of water-leaving radiance to surface irradiance at {wl} nm.'
            var.grid_mapping = "stereographic"
            var.coordinates = "longitude latitude"
            var.units = "sr^-1"
            var.source = "Sentinel-3a, Sentinel-3b, OLCIA-L3, OLCIB-L3"
            var.valid_min = 0.0
            var.valid_max = 1.0

        ##create angle variable
        name_angles = arc_opt.get_value_param(section, 'angle_bands', [], 'strlist')
        for angle in name_angles:
            if self.verbose:
                print(f'[INFO] Creating band: {angle}')
            info = olimage.get_angle_info(angle)
            var = datasetout.createVariable(angle, 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
            var.grid_mapping = "stereographic"
            var.coordinates = "longitude latitude"
            for at in info.keys():
                if at == '_FillValue':
                    continue
                if at == 'scale_factor':
                    continue
                var.setncattr(at, info[at])

        # create transparence bands
        name_transp = arc_opt.get_value_param(section, 'transp_bands', [], 'strlist')
        for transp in name_transp:
            if self.verbose:
                print(f'[INFO] Creating band: {transp}')
            var = datasetout.createVariable(transp, 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
            var.grid_mapping = "stereographic"
            var.coordinates = "longitude latitude"
            var.standard_name = "volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water"
            var.long_name = "OLCI Diffuse Attenuation Coefficient at 490nm"
            var.units = "m^-1"

        # create mask variable
        satellite_mask = datasetout.createVariable('mask', 'i2', ('y', 'x'), fill_value=-999, zlib=True,
                                                   shuffle=True, complevel=4)
        satellite_mask.standard_name = f'mask'
        satellite_mask.description = f'Masked pixels: 0; Valid pixels: 1'
        satellite_mask.grid_mapping = "stereographic"
        satellite_mask.coordinates = "longitude latitude"

        return datasetout

    def get_all_data_bands_names(self, arc_opt):
        section = 'RESAMPLE'
        rrs_bands = arc_opt.get_value_param(section, 'rrs_bands', self.olci_l2_bands, 'rrslist')
        angles_bands = arc_opt.get_value_param(section, 'angle_bands', [], 'strlist')
        transp_bands = arc_opt.get_value_param(section, 'transp_bands', [], 'strlist')
        all_bands = rrs_bands + angles_bands + transp_bands
        return all_bands

    def create_nc_file_resampled_from_pmldataset(self, ofname, ncpml):
        datasetout = self.copy_nc_base(ofname)
        if datasetout is None:
            return datasetout

        ##create variables
        for var in ncpml.variables:

            if var.lower().startswith('lat') or var.lower().startswith('lon') or var.lower().startswith('time'):
                continue
            if self.verbose:
                print(f'[INFO] [PML] Adding variable: {var}')
            datasetout.createVariable(var, 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=4, shuffle=True)

        return datasetout

    def copy_nc_base(self, ofile):
        dst = None
        if not os.path.exists(self.ifile_base):
            print(f'[ERROR] Grid file base: {self.ifile_base} does not exist')
            return dst
        if self.verbose:
            print(f'[INFO] Copying file grid base: {self.ifile_base}...')

        # shutil.copy(self.ifile_base,ofile)

        cmd = f'cp -a "{self.ifile_base}" "{ofile}"'
        if self.verbose:
            print(f'[INFO] cmd: {cmd}')
        import subprocess
        import time
        subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        originalSize = os.path.getsize(self.ifile_base)
        historicalSize = -1
        while historicalSize != originalSize:
            if os.path.exists(ofile):
                historicalSize = os.path.getsize(ofile)
                if self.verbose:
                    porc = (historicalSize / originalSize) * 100
                    print(f'[INFO] Copying {porc:.2f} %')
            time.sleep(1)
        if self.verbose:
            print('[INFO] Copy completed')

        dst = Dataset(ofile, 'a', format='NETCDF4')
        return dst

    def show_map_chl_impl(self, ax, data, lats, lons):
        limits, subareadef = self.get_subarea_def(lats, lons)
        crs = subareadef.to_cartopy_crs()
        from matplotlib.colors import LogNorm
        img = ax.imshow(data, transform=crs, extent=crs.bounds, origin='upper', norm=LogNorm(vmin=0.001, vmax=100))
        return img

    def save_full_map_impl(self, fileout, data, lats, lons):
        import matplotlib.pyplot as plt
        import cartopy
        if lats is not None and lons is not None:
            limits, subareadef = self.get_subarea_def(lats, lons)
            crs = subareadef.to_cartopy_crs()
        else:
            crs = self.area_def.to_cartopy_crs()

        ax = plt.axes(projection=crs)
        coastlines = ax.coastlines()
        img = plt.imshow(data, transform=crs, extent=crs.bounds, origin='upper')
        cbar = plt.colorbar()
        plt.savefig(fileout)

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

    def save_quick_look_fdata(self, fileout, fdata, name_var):
        dataset = Dataset(fdata)
        if dataset.variables[name_var].ndim == 3:
            data = np.ma.array(dataset.variables[name_var][0, :, :])
        if dataset.variables[name_var].ndim == 2:
            data = np.ma.array(dataset.variables[name_var][:, :])
        # data = np.ma.masked_values(data, 0)
        # print(type(data))
        # import numpy.ma as ma
        # data = ma.masked_greater(data, 0.5)
        self.save_quick_look_impl(fileout, data)

    def save_full_fdata(self, fileout, fdata, name_var):
        if self.verbose:
            print(f'[INFO] Saving Quick Look Image from file: {fdata}')
            print(f'[INFO] Variable: {name_var}')

        dataset = Dataset(fdata)

        dateherestr = None
        if 'time' in dataset.variables:
            from datetime import datetime as dt
            from datetime import timedelta
            datehere = dt(1981, 1, 1) + timedelta(seconds=float(dataset.variables['time'][0]))
            var_count = f'{name_var}_count' ##count variable in the file indicates that this is a monthly file
            if var_count in dataset.variables:
                dateherestr =  datehere.strftime('%B %Y')
            else:
                dateherestr = datehere.strftime('%Y-%m-%d')

        if dateherestr is None and name_var == 'MEDIAN':
            try:
                from datetime import datetime as dt
                datehere_n = os.path.basename(fileout).split('_')[0]
                datehere = dt.strptime(datehere_n,'%Y%m%d')
                dateherestr = datehere.strftime('%d %B')
            except:
                pass

        if self.verbose:
            print(f'[INFO] Date: {dateherestr}')
            print(f'[INFO] Starting figure and axes')

        fig, ax = self.start_full_figure()


        ##plotting images
        xstep = int(np.ceil(self.area_def.width / 2))
        ystep = int(np.ceil(self.area_def.height / 2))
        nsteps = 4
        iprogress = 1

        ##plotting images
        # xstep = int(np.ceil(self.area_def.width / 4))
        # ystep = int(np.ceil(self.area_def.height / 4))
        # nsteps = 8
        # iprogress = 1


        dataset_grid = None
        if 'lat' not in dataset.variables or 'lon' not in dataset.variables:
            dataset_grid = Dataset(self.ifile_base)



        for y in range(0, self.area_def.height, ystep):
            for x in range(0, self.area_def.width, xstep):
                if self.verbose:
                    print(f'[INFO] Plotting image({iprogress}/{nsteps})...')
                iprogress = iprogress + 1
                limits = self.get_limits(y, x, ystep, xstep, self.area_def.height, self.area_def.width)
                yini = limits[0]
                yfin = limits[1]
                xini = limits[2]
                xfin = limits[3]
                if len(dataset.variables[name_var].shape)==3:
                    data = np.ma.array(dataset.variables[name_var][0, yini:yfin, xini:xfin])
                elif len(dataset.variables[name_var].shape) == 2:
                    data = np.ma.array(dataset.variables[name_var][yini:yfin, xini:xfin])

                if name_var=='MEDIAN':
                    data = np.ma.power(10,data)

                if dataset_grid is not None:
                    lats = np.array(dataset_grid.variables['lat'][yini:yfin, xini:xfin])
                    lons = np.array(dataset_grid.variables['lon'][yini:yfin, xini:xfin])
                else:
                    lats = np.array(dataset.variables['lat'][yini:yfin, xini:xfin])
                    lons = np.array(dataset.variables['lon'][yini:yfin, xini:xfin])

                if name_var == 'CHL' or name_var == 'MEDIAN':
                    img = self.show_map_chl_impl(ax, data, lats, lons)

        if dataset_grid is not None:
            dataset_grid.close()

        if self.verbose:
            print(f'[INFO] Setting colorbar and title...')

        ##color bar and title
        if name_var == 'CHL' or name_var == 'MEDIAN':
            cbar = fig.colorbar(img, format="$%.2f$", anchor=(0.1, 0.5))
            cbar.ax.tick_params(labelsize=20)
            units = r'mg m$^-$$^3$'
            cbar.set_label(label=f'CHL ({units})', size=20)
            title = f'Chlorophyll a concentration ({units})'
            if dateherestr is not None:
                title = f'{title} - {dateherestr}'
            ax.set_title(title, fontsize=25, pad=36)



        fig.savefig(fileout, dpi=150, bbox_inches='tight')
        self.close_figure(fig)
        dataset.close()

        if self.verbose:
            print(f'[INFO] Completed')

    ##CREATE FIGURE AND AXES
    def close_figure(self, fig):
        from matplotlib import pyplot as plt
        plt.close(fig)

    def start_full_figure(self):
        from matplotlib import pyplot as plt
        import matplotlib.path as mpath
        import matplotlib.ticker as mticker


        # start figure and axes
        crs = self.area_def.to_cartopy_crs()
        fig, ax = plt.subplots(subplot_kw=dict(projection=crs))
        fig.set_figwidth(15)
        fig.set_figheight(15)

        # coastlines
        ax.coastlines(resolution='50m')

        ##NOT WORKING WITH PROBLEMS WITH THE BOUNDARIES OF THE POLYGONS
        ##ax.add_feature(cartopy.feature.LAND, edgecolor='black', linewidth=0.5, transform = cartopy.crs.PlateCarree())


        # Prep circular boundary
        r_extent = self.area_def.area_extent[1]
        r_extent *= 1.005
        circle_path = mpath.Path.unit_circle()
        circle_path = mpath.Path(circle_path.vertices.copy() * r_extent, circle_path.codes.copy())
        ax.set_boundary(circle_path)
        ax.set_frame_on(False)

        # gri lines
        gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=True)
        gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135, 180])
        gl.ylocator = mticker.FixedLocator([65, 70, 75, 80, 85])
        gl.rigth_labels = True
        gl.left_labels = True
        gl.bottom_labels = True
        gl.top_labels = True
        gl.xlabel_style = {'size': 15}
        gl.ylabel_style = {'size': 15}
        plt.draw()
        for ea in gl.label_artists:
            txt = ea.get_text()
            pos = ea.get_position()
            if txt == '135°W' and pos[0] < (-2000000):
                ea.set_visible(True)
            if txt == '45°E' and pos[0] > 2000000:
                ea.set_visible(True)
            if pos[0] == 90:
                ea.set_visible(True)
                ea.set_position([135, pos[1]])
            if pos[0] == -90:
                ea.set_visible(True)
                ea.set_position([-45, pos[1]])
            if pos[1] == 65:
                ea.set_visible(False)

        return fig, ax

    def get_limits(self, y, x, ystep, xstep, ny, nx):
        yini = y
        xini = x
        yfin = y + ystep
        xfin = x + xstep
        if yfin > ny:
            yfin = ny
        if xfin > nx:
            xfin = nx

        limits = [yini, yfin, xini, xfin]
        return limits

    def get_lat_min_spherical(self):
        xmid = math.floor(self.area_def.width / 2)
        lon, lat = self.area_def.get_lonlat_from_array_coordinates(xmid, 0)
        lat_min = np.round(lat)
        return lat_min

    def get_lat_min(self):
        lon, lat_min = self.area_def.get_lonlat_from_array_coordinates(0, 0)
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
        lon5, lat5 = self.area_def.get_lonlat_from_array_coordinates(xmid, ymid)
        print('Middle: ', lat5, lon5)

        print('SPHERIC LIMITS:')
        lons = [np.degrees(ob[0].lon), np.degrees(ob[1].lon), np.degrees(ob[2].lon), np.degrees(ob[3].lon), lon1, lon2,
                lon3, lon4]
        lonmin = round(np.min(lons))
        lonmax = round(np.max(lons))
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

    def set_area_definition(self, area_id):
        self.area_def = self.get_area_definition(area_id)
        ifile_base_default = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/GRID_FILES/ArcGrid_65_90_300mNOGEO.nc'
        print(self.ifile_base)
        print(ifile_base_default)
        if area_id == 'polar_stereographic_4km' and self.ifile_base==ifile_base_default:
            self.ifile_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/GRID_FILES/ArcGrid_65_90_4KM_GridBase.nc'
        if self.verbose:
            print('--------------------------------------------------------------------------------------------')
            print('UPDATE AREA DEFINITION')
            self.print_area_def_info()
            print('---------------------------------------------------------------------------------------------')

    def get_projection_info(self, proj_id):
        projections = {
            'npsn': '+proj=stere +lat_0=90 +lon_0=-45 +k_0=1 +x_0=0 +y_0=0 +ellps=WGS84',
            'polarn': '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k_0=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs +type=crs',
            'polarn_orig': '+proj=stere +lat_0=90 +lon_0=-45 +k_0=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs'
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
            'polar_stereographic': {
                'description': 'WGS 84 / Polar Stereographic North',
                'proj_id': 'polarn',
                'width': 18345,
                'height': 18345,
                'extent': (-2751778, 2751778, 2751778, -2751778)
            },
            'polar_stereographic_orig': {
                'description': 'WGS 84 / Polar Stereographic North',
                'proj_id': 'polarn_orig',
                'width': 18915,
                'height': 18915,
                'extent': (-2837300, -2837300, 2837200, 2837200)
            },
            'polar_stereographic_3km': {
                'description': 'WGS 84 / Polar Stereographic North',
                'proj_id': 'polarn',
                'width': 1892,
                'height': 1892,
                'extent': (-2837300, -2837300, 2837200, 2837200)
            },
            'polar_stereographic_4km': {
                'description': 'WGS 84 / Polar Stereographic North',
                'proj_id': 'polarn',
                'width': 1375,
                'height': 1375,
                'extent': (-2752500, 2747500, 2747500, -2752500)
            },
            'polar_stereographic_600m': {
                'description': 'WGS 84 / Polar Stereographic North',
                'proj_id': 'polarn',
                'width': 6000,
                'height': 6000,
                'extent': (-2837300, -2837300, 2837200, 2837200)
            },
            'polar_stereographic_nerc': {
                'description': 'WGS 84 / Polar Stereographic North',
                'proj_id': 'npsn',
                'width': 609,
                'height': 881,
                'extent': (-3800000, 5500000, 3800000, -5500000)
            }
        }
        return areas[area_id]

    def get_subarea_def_from_array_coordinates(self, ymin, ymax, xmin, xmax):
        width = (xmax - xmin)
        height = (ymax - ymin)
        projection = self.get_projection_info(self.area_def.proj_id)
        xvalues = [xmin, xmax, xmax, xmin]
        yvalues = [ymin, ymin, ymax, ymax]
        xcoords, ycoords = self.area_def.get_projection_coordinates_from_array_coordinates(xvalues, yvalues)
        extent = (np.min(xcoords), np.max(ycoords), np.max(xcoords), np.min(ycoords))

        area_def = AreaDefinition('SubArea', 'SubArea', self.area_def.proj_id, projection,
                                  width, height, extent)

        lons, lats = self.area_def.get_lonlat_from_array_coordinates(xvalues, yvalues)
        min_lon = np.min(lons[:])
        max_lon = np.max(lons[:])
        min_lat = np.min(lats[:])
        max_lat = np.max(lats[:])
        geo_limits = [min_lat, max_lat, min_lon, max_lon]

        return geo_limits, area_def

    def get_subarea_def(self, lats, lons):

        xcoords, ycoords = self.area_def.get_array_coordinates_from_lonlat(lons, lats)

        xmin = np.floor(np.min(xcoords))
        if xmin < 0:
            xmin = 0
        xmax = np.ceil(np.max(xcoords)) + 1
        if xmax > self.area_def.width:
            xmax = self.area_def.width
        ymin = np.floor(np.min(ycoords))
        if ymin < 0:
            ymin = 0
        ymax = np.ceil(np.max(ycoords)) + 1
        if ymax > self.area_def.height:
            ymax = self.area_def.height

        width = (xmax - xmin)
        height = (ymax - ymin)
        # print(width, height)
        projection = self.get_projection_info(self.area_def.proj_id)
        # print(projection)

        xvalues = [xmin, xmax, xmax, xmin]
        yvalues = [ymin, ymin, ymax, ymax]
        # xvalues = [xmin_new, xmax_new, xmax_new, xmin_new]
        # yvalues = [ymin_new, ymin_new, ymax_new, ymax_new]
        xcoords, ycoords = self.area_def.get_projection_coordinates_from_array_coordinates(xvalues, yvalues)
        # lontal,lattal = self.area_def.get_lonlat_from_array_coordinates(xvalues,yvalues)
        # print('******************')
        # print('yvalues:',yvalues)
        # print('xvalues:',xvalues)
        # print('lontal:',lontal)
        # print('lattal:',lattal)

        extent = (np.min(xcoords), np.max(ycoords), np.max(xcoords), np.min(ycoords))

        # print(extent)
        # print('*********************')

        area_def = AreaDefinition('SubArea', 'SubArea', self.area_def.proj_id, projection,
                                  width, height, extent)

        # print(area_def.pixel_size_x, area_def.pixel_size_y)
        # print(area_def.width, area_def.height)

        limits = [xmin, xmax, ymin, ymax]

        # dif = xmax-xmin

        return limits, area_def

    def check_make_resample_impl(self, olimage, arc_opt):
        if self.verbose:
            print('--------------------')
            print(f'[INFO] Checking resampling for granule: {olimage.path_source}')
        section = 'RESAMPLE'
        # Source;Width;Height;NTotal;NFlagged;NWater1;NWater2;NValid;PValid'
        rrs_bands = arc_opt.get_value_param(section, 'rrs_bands', self.olci_l2_bands, 'floatlist')
        olimage.set_reflectance_bands_mask(rrs_bands)
        mask, res_original, line_original = olimage.get_mask_default()
        th_nvalid = arc_opt.get_value_param(section, 'th_nvalid', -1, 'int')

        if res_original[7] <= th_nvalid and th_nvalid >= 0:
            print('[WARNING] No valid pixels were found. Skipping granule...')
            return None
        params_mask = [mask, res_original, line_original]
        return params_mask

    def make_resample_impl(self, olimage, fileout, granuleindex, orbitindex, arc_opt, params_mask):
        if self.verbose:
            print('--------------------')
            print(f'[INFO] Starting resampling for granule: {olimage.path_source}')
        section = 'RESAMPLE'
        rrs_bands = arc_opt.get_value_param(section, 'rrs_bands', self.olci_l2_bands, 'floatlist')
        olimage.set_reflectance_bands_mask(rrs_bands)
        start_date = olimage.get_start_date()
        start_date_str = start_date.strftime('%Y%m%dT%H%M%S')
        th_nvalid = arc_opt.get_value_param(section, 'th_nvalid', -1, 'int')

        if params_mask is None:
            mask, res_original, line_original = olimage.get_mask_default()
        else:
            mask = params_mask[0]
            res_original = params_mask[1]
            line_original = params_mask[2]

        ##CHECKING MASK
        # Source;Width;Height;NTotal;NFlagged;NWater1;NWater2;NValid;PValid'
        if res_original[7] <= th_nvalid and th_nvalid >= 0:
            res_resampled = [-999] * 9
            line_resampled = ';'.join(str(l) for l in res_resampled)
            line_output = f'{olimage.name_source};{start_date_str};{olimage.get_rel_pass()};{granuleindex};-999;{line_original};{line_resampled}'
            print('[WARNING] No valid pixels were found. Skipping granule...')
            return line_output
        ##END CHECKING

        if self.verbose:
            print(f'[INFO][{granuleindex}] Defining subarea for the granule...')
        lats, lons = olimage.get_lat_long_arrays()

        if self.verbose:
            print(
                f'[INFO][{granuleindex}] Granule geographic limits. Lat: {np.min(lats)} - {np.max(lats)}  Long: {np.min(lons)} - {np.max(lons)}')

        limits, sub_area_def = self.get_subarea_def(lats, lons)
        swath_def = SwathDefinition(lons=lons, lats=lats)
        if self.verbose:
            print(f'[INFO][{granuleindex}] Array limits after resampling: {limits}')
            print(f'[INFO][{granuleindex}] Getting neighbour info...')
        # row_indices, col_indices = utils.generate_nearest_neighbour_linesample_arrays(swath_def, sub_area_def, 2000)
        valid_input_index, valid_output_index, index_array, distance_array = get_neighbour_info(swath_def, sub_area_def,
                                                                                                500, neighbours=1)

        output_shape = (sub_area_def.height, sub_area_def.width)
        ymin = int(limits[2])
        ymax = int(limits[3])
        xmin = int(limits[0])
        xmax = int(limits[1])

        if self.verbose:
            print(f'[INFO][{granuleindex}] Creating output file from file grid...')
        datasetout = self.create_nc_file_resampled(fileout, olimage, arc_opt)

        if datasetout is None:
            print(f'[ERROR][{granuleindex}] Output file {fileout} could not be created')

        if self.verbose:
            print(f'[INFO][{granuleindex}] Resampling variables...')
            print(f'[INFO][{granuleindex}] --->Mask')
        # mask (using fill_value=10 becuase of the problems resampling with fill_value=-999)
        result = get_sample_from_neighbour_info('nn', output_shape, mask, valid_input_index,
                                                valid_output_index, index_array,
                                                distance_array=distance_array, fill_value=10)
        result_m = np.zeros(result.shape, dtype=np.int)
        result_m[result == 10] = -999
        result_m[result == 0] = 1
        resampled_n_total = np.count_nonzero(result_m >= 0)
        if resampled_n_total == 0:
            print(
                f'[WARNING][{granuleindex}] No valid pixels were resampled for granule {olimage.path_source}. Removing file and skipping')
            datasetout.close()
            os.remove(fileout)
            return
        resampled_n_valid = np.count_nonzero(result_m == 1)
        var = datasetout.variables['mask']
        var[ymin:ymax, xmin:xmax] = [result_m[:, :]]

        ##REMAINING VARIABLES
        all_variables = self.get_all_data_bands_names(arc_opt)
        angles_bands = arc_opt.get_value_param(section, 'angle_bands', [], 'strlist')
        for var_name in all_variables:
            if self.verbose:
                print(f'[INFO][{granuleindex}] --->{var_name}')
            var = datasetout.variables[var_name]
            fvalue = var.getncattr('_FillValue')
            if var_name.startswith('RRS'):
                wlref = var.wavelength
                array = olimage.get_reflectance_band_array(wlref, fvalue)
                array[mask == 1] = -999
            elif var_name in angles_bands:
                array = olimage.get_angle_array(var_name)
            else:
                array = olimage.get_general_array(var_name, fvalue)

            result = get_sample_from_neighbour_info('nn', output_shape, array, valid_input_index,
                                                    valid_output_index, index_array,
                                                    distance_array=distance_array, fill_value=fvalue)
            var[ymin:ymax, xmin:xmax] = [result[:, :]]

        # sensor mask
        if self.verbose:
            print(f'[INFO][{granuleindex}] --->Sensor mask')
        var = datasetout.variables['sensor_mask']
        array = np.array(datasetout.variables['sensor_mask'][ymin:ymax, xmin:xmax])
        # sensorflag = granuleindex
        # if orbitindex >= 0:
        #     sensorflag = math.pow(2, orbitindex / 100)
        array[result_m != -999] = array[result_m != -999] + orbitindex
        var[ymin:ymax, xmin:xmax] = [array[:, :]]

        datasetout.granule_index = int(granuleindex)
        datasetout.relative_orbit = int(olimage.get_rel_pass())
        datasetout.orbit_index = int(orbitindex)
        datasetout.platform = olimage.get_platform()

        datasetout.start_date = start_date_str
        lc = []
        for c in olimage.coords_image:
            lc.append(str(c))
        str_coords = ';'.join(lc)
        datasetout.source = olimage.name_source
        datasetout.geo_polygon = str_coords
        datasetout.resampled_ymin = ymin
        datasetout.resampled_ymax = ymax
        datasetout.resampled_xmin = xmin
        datasetout.resampled_xmax = xmax
        datasetout.resampled_width = sub_area_def.width
        datasetout.resampled_height = sub_area_def.height
        datasetout.resampled_n_total = resampled_n_total
        datasetout.resampled_n_valid = resampled_n_valid
        resampled_p_valid = (resampled_n_valid / resampled_n_total) * 100
        datasetout.resampled_p_valid = resampled_p_valid

        # Source;Width;Height;NTotal;NFlagged;NWater1;NWater2;NValid;PValid'
        datasetout.original_width = res_original[1]
        datasetout.original_height = res_original[2]
        datasetout.original_total = res_original[3]
        datasetout.original_n_flagged = res_original[4]
        datasetout.original_n_water1 = res_original[5]
        datasetout.original_n_water2 = res_original[6]
        datasetout.original_n_valid = res_original[7]
        datasetout.original_p_valid = res_original[8]

        datasetout.close()

        # line_original = ';'.join(str(l) for l in res_original)
        res_resampled = [ymin, ymax, xmin, xmax, sub_area_def.width, sub_area_def.height, resampled_n_total,
                         resampled_n_valid, resampled_p_valid]
        line_resampled = ';'.join(str(l) for l in res_resampled)
        line_original = line_original[line_original.find(';'):]  # remove source from line original
        line_output = f'{olimage.name_source};{start_date_str};{olimage.get_rel_pass()};{granuleindex};{orbitindex};{line_original};{line_resampled}'

        if self.verbose:
            print(f'[INFO][{granuleindex}] Completed. Output file: {fileout}')

        return line_output

    def make_resample_from_file_orig_multi_dep(self, input_file, output_file, date_here):
        if not os.path.exists(self.ifile_base):
            print(f'[ERROR] File grid base: {self.ifile_base} does not exist')
            return
        if self.verbose:
            print('[INFO] Starting resample...')
            print(f'[INFO] Output file: {output_file}')

        datasetout = self.copy_nc_base(output_file)

        if datasetout is None:
            return datasetout

        if self.verbose:
            print(f'[INFO] Adding variables from file orig: {input_file}')

        ncdataset = Dataset(input_file, 'r')
        # vartime = False
        for name in ncdataset.variables:
            if name.lower().startswith('lat'):
                variable_lat = name
            if name.lower().startswith('lon'):
                variable_lon = name
            # if name.lower().startswith('time'):
            #     vartime = True

        lat_array = ncdataset.variables[variable_lat][0:600]
        lon_array = ncdataset.variables[variable_lon][:]
        # nlon = len(lon_array)
        # xstep = 600
        step = 600

        for name, variable in ncdataset.variables.items():
            if name.lower().startswith('lat') or name.lower().startswith('lon') or name.lower().startswith('time'):
                continue
            if name.startswith('ZSD') or name.startswith('SPM') or name == 'flags' or name.endswith('uncertainty'):
                continue

            if self.verbose:
                print(f'[INFO]  Adding variable: {name}')

            var_output = datasetout.createVariable(name, 'f4', ('time', 'y', 'x'), fill_value=-999.0, zlib=True,
                                                   complevel=4,
                                                   shuffle=True)

            var_output.stardard_name = ncdataset[name].standard_name
            var_output.type = ncdataset[name].type
            var_output.coordinates = 'time lon lat'
            var_output.grid_mapping = 'stereographic'

            valid_min_at = None
            valid_max_at = None
            missing_value_at = -999.0
            if 'valid_min' in ncdataset[name].ncattrs():
                valid_min_at = ncdataset[name].valid_min
            elif 'min_val' in ncdataset[name].ncattrs():
                valid_min_at = ncdataset[name].min_val

            if 'valid_max' in ncdataset[name].ncattrs():
                valid_max_at = ncdataset[name].valid_min
            elif 'max_val' in ncdataset[name].ncattrs():
                valid_max_at = ncdataset[name].max_val

            if 'missing_value' in ncdataset[name].ncattrs():
                missing_value_at = ncdataset[name].missing_value
            elif '_FillValue' in ncdataset[name].ncattrs():
                missing_value_at = ncdataset[name]._FillValue

            var_output.valid_min = valid_min_at
            var_output.valid_max = valid_max_at
            var_output.missing_value = missing_value_at

            if name.startswith('RRS'):
                var_output.long_name = f'Remote Sensing Reflectance at {name.lower()}'
                var_output.units = 'sr^-1'

            if name == 'KD490':
                var_output.long_name = 'Diffuse Attenuation Coefficient at 490nm'
                var_output.units = 'm^-1'

            for y in range(0, self.area_def.width, step):
                for x in range(0, self.area_def.height, step):
                    ymax = y + step
                    xmax = x + step
                    if ymax > self.area_def.height:
                        ymax = self.area_def.height
                    if xmax > self.area_def.width:
                        xmax = self.area_def.width
                    geo_limits, area_def = self.get_subarea_def_from_array_coordinates(y, ymax, x, xmax)

                    lon_diff = geo_limits[3] - geo_limits[2]
                    print(y, ymax, x, xmax, geo_limits, lon_diff)

    def make_resample_from_file_orig_multi(self, input_file, output_file, date_here):

        if not os.path.exists(self.ifile_base):
            print(f'[ERROR] File grid base: {self.ifile_base} does not exist')
            return
        if self.verbose:
            print('[INFO] Starting resample...')
            print(f'[INFO] Output file: {output_file}')

        datasetout = self.copy_nc_base(output_file)

        if datasetout is None:
            return datasetout

        # lat_big = np.array(datasetout.variables['lat'])
        # lon_big = np.array(datasetout.variables['lon'])
        addrrsqibands = False
        if self.climatology_path is not None and os.path.exists(self.climatology_path):
            addrrsqibands = True

        if self.verbose:
            print(f'[INFO] Adding variables from file orig: {input_file}')

        ncdataset = Dataset(input_file, 'r')
        # vartime = False
        for name in ncdataset.variables:
            if name.lower().startswith('lat'):
                variable_lat = name
            if name.lower().startswith('lon'):
                variable_lon = name
            # if name.lower().startswith('time'):
            #     vartime = True

        lat_array = ncdataset.variables[variable_lat][0:600]
        lon_array = ncdataset.variables[variable_lon][:]
        nlon = len(lon_array)
        xstep = 2160

        for name, variable in ncdataset.variables.items():
            if name.lower().startswith('lat') or name.lower().startswith('lon') or name.lower().startswith('time'):
                continue
            if name.startswith('ZSD') or name.startswith('SPM') or name == 'flags' or name.endswith('uncertainty'):
                continue

            if self.verbose:
                print(f'[INFO]  Adding variable: {name}')

            var_output = datasetout.createVariable(name, 'f4', ('time', 'y', 'x'), fill_value=-999.0, zlib=True,
                                                   complevel=4,
                                                   shuffle=True)

            var_output.stardard_name = ncdataset[name].standard_name
            var_output.type = ncdataset[name].type
            var_output.coordinates = 'time lon lat'
            var_output.grid_mapping = 'stereographic'

            valid_min_at = None
            valid_max_at = None
            missing_value_at = -999.0
            if 'valid_min' in ncdataset[name].ncattrs():
                valid_min_at = ncdataset[name].valid_min
            elif 'min_val' in ncdataset[name].ncattrs():
                valid_min_at = ncdataset[name].min_val

            if 'valid_max' in ncdataset[name].ncattrs():
                valid_max_at = ncdataset[name].valid_min
            elif 'max_val' in ncdataset[name].ncattrs():
                valid_max_at = ncdataset[name].max_val

            if 'missing_value' in ncdataset[name].ncattrs():
                missing_value_at = ncdataset[name].missing_value
            elif '_FillValue' in ncdataset[name].ncattrs():
                missing_value_at = ncdataset[name]._FillValue

            var_output.valid_min = valid_min_at
            var_output.valid_max = valid_max_at
            var_output.missing_value = missing_value_at

            if name.startswith('RRS'):
                var_output.long_name = f'Remote Sensing Reflectance at {name.lower()}'
                var_output.units = 'sr^-1'

            if name == 'KD490':
                var_output.long_name = 'Diffuse Attenuation Coefficient at 490nm'
                var_output.units = 'm^-1'

            for x in range(0, nlon, xstep):
                print(f'[INFO] -> {x}')
                yini = 0
                yfin = 600
                xini = x
                xfin = x + xstep
                if xfin > nlon:
                    xfin = nlon

                lon_here = lon_array[xini:xfin]
                lat_here = lat_array[yini:yfin]

                lat2D, lon2D = self.get_2D_lat_lon_arrays(lat_here, lon_here)
                limits, sub_area_def = self.get_subarea_def(lat2D, lon2D)

                swath_def = SwathDefinition(lons=lon2D, lats=lat2D)
                valid_input_index, valid_output_index, index_array, distance_array = get_neighbour_info(swath_def,
                                                                                                        sub_area_def,
                                                                                                        6800,
                                                                                                        neighbours=1)

                rrsarray = np.array(variable[0, yini:yfin, xini:xfin])
                result = get_sample_from_neighbour_info('nn', sub_area_def.shape, rrsarray, valid_input_index,
                                                        valid_output_index, index_array,
                                                        distance_array=distance_array,
                                                        fill_value=-999)
                ymin = int(limits[2])
                ymax = int(limits[3])
                xmin = int(limits[0])
                xmax = int(limits[1])

                var_output_array = np.array(var_output[0, ymin:ymax, xmin:xmax])
                var_output_array[result != -999] = result[result != -999]
                var_output[0, ymin:ymax, xmin:xmax] = var_output_array[:, :]

            if name.startswith('RRS') and addrrsqibands:
                nameqi = f'QI_{name}'
                var_output_qi = datasetout.createVariable(nameqi, 'f4', ('time', 'y', 'x'), fill_value=-999.0,
                                                          zlib=True,
                                                          complevel=4,
                                                          shuffle=True)
                rrs = name[-3:]
                var_output_qi.long_name = f'Quality Index for  remote sensing reflectance at {rrs} nm'
                var_output_qi.comment = 'QI=(DailyData-ClimatologyMedian)/ClimatologyStandardDeviation'
                var_output_qi.type = 'surface'
                var_output_qi.units = '1'
                var_output_qi.missing_value = -999.0
                var_output_qi.valid_min = -5.0
                var_output_qi.valid_max = 5.0
                file_clima = self.get_file_climatology(name, date_here)
                if file_clima is not None:
                    if self.verbose:
                        print(f'[INFO] Computing  climatology uisng file: {file_clima}')
                    output_array = self.compute_climatology(file_clima, var_output)
                    var_output_qi[0, :, :] = output_array[:, :]

        from datetime import datetime as dt
        if date_here is not None:
            nsec = (date_here - dt(1981, 1, 1)).total_seconds()
            datasetout.variables['time'][0] = [np.int(nsec)]
            datasetout.start_date = date_here.strftime('%Y-%m-%d')
            datasetout.stop_date = date_here.strftime('%Y-%m-%d')
        cdate = dt.utcnow()
        datasetout.creation_date = cdate.strftime('%Y-%m-%d')
        datasetout.creation_time = cdate.strftime('%H:%M:%S UTC')

        ncdataset.close()
        datasetout.close()

    def get_file_climatology(self, variable, date_work):
        path = os.path.join(self.climatology_path, variable.lower())
        if not os.path.exists(path):
            return None
        if date_work.month == 2 and date_work.day == 29:
            date_work = date_work.replace(day=28)
        ddmm = date_work.strftime('%m%d')
        if variable.lower() == 'kd490':
            variable = 'transp'
        fname = f'1998{ddmm}_2022{ddmm}_{variable.lower()}_arc_multi_clima.nc'
        file_clima = os.path.join(path, fname)
        if os.path.isfile(file_clima):
            return file_clima
        else:
            return None

    def compute_climatology(self, file_clima, input_variable):
        dclima = Dataset(file_clima)
        central = np.array(dclima.variables['MEDIAN'])
        dispersion = np.array(dclima.variables['N_STD'])
        dclima.close()
        input_array = np.array(input_variable)
        output_array = np.zeros((input_array.shape[1], input_array.shape[2]))
        output_array[:, :] = input_array[0, :, :]
        output_array[output_array != -999] = (output_array[output_array != -999] - central[output_array != -999]) / \
                                             dispersion[output_array != -999]

        return output_array

    def make_resample_pml(self, fpml, fileout):
        if self.verbose:
            print('--------------------')
            print(f'[INFO] Starting resampling for PML product: {fpml}')

        ncpml = Dataset(fpml)
        datasetout = self.create_nc_file_resampled_from_pmldataset(fileout, ncpml)

        # print(ncpml.variables)
        lat_array = ncpml.variables['latitude'][:]
        lon_array = ncpml.variables['longitude'][:]
        nlat = len(lat_array)
        nlon = len(lon_array)
        ystep = 2400
        xstep = 2400

        for y in range(0, nlat, ystep):
            # CHECKING NON MASKED VALUES

            yini = y
            yfin = y + ystep
            if yfin > nlat:
                yfin = nlat
            xini = 0
            xfin = nlon
            height = (yfin - yini)
            width = (xfin - xini)
            shape = (height, width)

            nnonmasked = 0
            for var in ncpml.variables:
                if var.lower().startswith('lat') or var.lower().startswith('lon') or var.lower().startswith('time'):
                    continue
                rrsarray = np.ma.zeros(shape)
                rrsarray[:, :] = np.ma.array(ncpml.variables[var][0, yini:yfin, xini:xfin])
                # print('Non masked: ',rrsarray.count())
                nnonmasked = nnonmasked + rrsarray.count()

            if self.verbose:
                print(f'[INFO] [PML] Counting non-masked values: {yini} to {yfin} / {nlat}: {nnonmasked}')
            if nnonmasked == 0:
                continue

            for x in range(0, nlon, xstep):
                xini = x
                xfin = x + xstep
                if xfin > nlon:
                    xfin = nlon
                height = (yfin - yini)
                width = (xfin - xini)
                shape = (height, width)

                nnonmasked = 0
                for var in ncpml.variables:
                    if var.lower().startswith('lat') or var.lower().startswith('lon') or var.lower().startswith('time'):
                        continue
                    rrsarray = np.ma.zeros(shape)
                    rrsarray[:, :] = np.ma.array(ncpml.variables[var][0, yini:yfin, xini:xfin])
                    # print('Non masked: ',rrsarray.count())
                    nnonmasked = nnonmasked + rrsarray.count()

                if self.verbose:
                    print(f'[INFO] [PML] Non-masked pixels in block x: {xini}-{xfin} y: {yini}-{yfin} -> {nnonmasked}')

                if nnonmasked == 0:
                    continue

                lon_here = lon_array[xini:xfin]
                lat_here = lat_array[yini:yfin]

                lat2D, lon2D = self.get_2D_lat_lon_arrays(lat_here, lon_here)
                limits, sub_area_def = self.get_subarea_def(lat2D, lon2D)
                swath_def = SwathDefinition(lons=lon2D, lats=lat2D)

                if self.verbose:
                    print(f'[INFO] [PML RESAMPLE]     Getting neighbour info...')
                # row_indices, col_indices = utils.generate_nearest_neighbour_linesample_arrays(swath_def, sub_area_def, 2000)

                valid_input_index, valid_output_index, index_array, distance_array = get_neighbour_info(swath_def,
                                                                                                        sub_area_def,
                                                                                                        500,
                                                                                                        neighbours=1)

                for var in ncpml.variables:
                    if var.lower().startswith('lat') or var.lower().startswith('lon') or var.lower().startswith('time'):
                        continue
                    if var.lower().endswith('count') or var.lower().endswith('error'):
                        continue
                    if self.verbose:
                        print(f'[INFO] [PML RESAMPLE]     Resampling variable: {var}...')
                    var_rrs = datasetout.variables[var]

                    # rrsarray = np.ma.zeros(shape)
                    rrsarray = np.array(ncpml.variables[var][0, yini:yfin, xini:xfin])

                    # print('Non masked before: ',rrsarray.count())
                    # print(rrsarray.shape)
                    # output_shape = (sub_area_def.height, sub_area_def.width)
                    # print('Before:', rrsarray.shape, np.max(rrsarray))
                    result = get_sample_from_neighbour_info('nn', sub_area_def.shape, rrsarray, valid_input_index,
                                                            valid_output_index, index_array,
                                                            distance_array=distance_array,
                                                            fill_value=-999)
                    # print('Afger:', result.shape, np.max(result), type(result))
                    # print('Non masked later: ',result.count())
                    ymin = int(limits[2])
                    ymax = int(limits[3])
                    xmin = int(limits[0])
                    xmax = int(limits[1])
                    var_rrs_array = np.array(var_rrs[ymin:ymax, xmin:xmax])
                    var_rrs_array[result != -999] = result[result != -999]
                    var_rrs[ymin:ymax, xmin:xmax] = [var_rrs_array[:, :]]

                    for yp in range(ymin, ymax):
                        for xp in range(xmin, xmax):
                            val = result[yp - ymin, xp - xmin]
                            if val != -999:
                                # print('-->',yp,xp)
                                var_rrs[yp, xp] = [val]
                    # var_rrs[ymin:ymax, xmin:xmax] = [result[:, :]]

        datasetout.close()
        ncpml.close()
        print('DONE')
        # print('Dim',nlat,nlon)
        # for x in range(nlon):
        #     print(lon_array[x])

    def get_2D_lat_lon_arrays(self, lat1D, lon1D):
        height = len(lat1D)
        width = len(lon1D)
        lat2D = np.zeros((height, width), lat1D.dtype)
        lon2D = np.zeros((height, width), lat1D.dtype)
        for y in range(height):
            lat2D[y, :] = lat1D[y]
        for x in range(width):
            lon2D[:, x] = lon1D[x]
        return lat2D, lon2D

    def deprecated(self):
        print('deprecated')
# for var in datasetout.variables:
#     if var.startswith('RRS'):
#         if self.verbose:
#             print(f'[INFO] --->{var}')
#         var_rrs = datasetout.variables[var]
#
#         wlref = var_rrs.wavelength
#         if wlref > 0:
#             rrsarray = olimage.get_reflectance_band_array(wlref)
#             rrsarray[mask == 1] = -999
#
#             result = get_sample_from_neighbour_info('nn', output_shape, rrsarray, valid_input_index,
#                                                     valid_output_index, index_array,
#                                                     distance_array=distance_array, fill_value=-999)
#             var_rrs[ymin:ymax, xmin:xmax] = [result[:, :]]
# name_angles = ['OZA', 'OAA', 'SZA', 'SAA']
# for angle in name_angles:
#     if self.verbose:
#         print(f'[INFO] --->{angle}')
#     var = datasetout.variables[angle]
#     fvalue = var.getncattr('_FillValue')
#
#     array = olimage.get_angle_array(angle)
#     result = get_sample_from_neighbour_info('nn', output_shape, array, valid_input_index, valid_output_index,
#                                             index_array, distance_array=distance_array, fill_value=fvalue)
#     # print(np.min(array), np.max(array), np.min(result), np.max(result))
#     var[ymin:ymax, xmin:xmax] = [result[:, :]]
