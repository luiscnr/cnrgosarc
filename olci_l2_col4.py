import os, math
import numpy as np
from netCDF4 import Dataset
import Class_Flags_OLCI as flag
from check_geo import CHECK_GEO
from datetime import datetime as dt

class OLCI_L2_COL4:

    def __init__(self, path_source, verbose):
        self.verbose = verbose
        self.path_source = path_source
        self.name_source = self.path_source.split('/')[-1]
        self.wl_list = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75, 753.75, 778.75, 865, 885,1020.5]
        self.reflectance_bands, self.nbands = self.get_reflectance_bands_info()
        self.reflectance_bands_mask = self.reflectance_bands
        self.max_dif_wl = 1
        self.width = -1
        self.height = -1
        self.params = {}
        self.coords_image = None

    def get_reflectance_bands_info(self):
        reflectance_bands = {}
        bands_ints = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18, 21]

        nbands = len(bands_ints)
        for v,wl in zip(bands_ints,self.wl_list):
            # v = bands_ints[index]
            # wl = self.wl_list[index]
            band_name = f'Oa{v:02d}_reflectance'
            band_path = f'{band_name}.nc'
            file_path = os.path.join(self.path_source, band_path)
            reflectance_bands[band_name] = {
                'wavelength': wl,
                'file_path': file_path
            }
        return reflectance_bands, nbands

    def check_granule(self):
        nfiles = 0
        valid = True
        for name in os.listdir(self.path_source):
            if name!='iop_lsd.nc':
                nfiles = nfiles + 1
            if name.endswith('.nc'):
                try:
                    dataset = Dataset(os.path.join(self.path_source,name))
                    dataset.close()
                except Exception as ex:
                    print(f'[ERROR] Error reading {name}:  {ex}')
                    valid = False
                    break
        if nfiles<32:
            valid = False
        return valid

    def get_geo_and_params(self):
        self.get_dimensions()
        cgeo = CHECK_GEO()
        cgeo.start_polygon_from_prod_manifest_file(self.path_source)
        self.params = cgeo.params
        self.coords_image = cgeo.coords_image



    def get_dimensions(self):
        coordinates_filename = 'geo_coordinates.nc'
        filepah = os.path.join(self.path_source, coordinates_filename)
        nc_sat = Dataset(filepah, 'r')
        lat = nc_sat.variables['latitude']
        shape = lat.shape
        if lat.dimensions[0] == 'rows' and lat.dimensions[1] == 'columns':
            self.height = shape[0]
            self.width = shape[1]
        else:
            self.height = shape[0]
            self.width = shape[1]

    def get_platform(self):
        if self.path_source.startswith('S3A'):
            platform = 'A'
        elif self.path_source.startswith('S3B'):
            platform = 'B'
        else:
            self.get_geo_and_params()
            platform = self.params['number']
        return platform

    def get_rel_pass(self):
        if 'relativeOrbitNumber' in self.params.keys():
            return int(self.params['relativeOrbitNumber'])
        else:
            name_list = self.path_source.split('/')[-1][:-3].split('_')
            try:
                rel_pass = int(name_list[12])
            except:
                rel_pass = -1
            return rel_pass

    def get_start_date(self):
        if 'startTime' in self.params.keys():
            stime = self.params['startTime']
            try:
                dtime = dt.strptime(stime, '%Y-%m-%dT%H:%M:%S.%fZ')
            except:
                dtime = None
        else:
            name_list = self.path_source.split('/')[-1][:-3].split('_')
            stime = name_list[7]
            try:
                dtime = dt.strptime(stime, '%Y%m%dT%H%M%S')
            except:
                dtime = None

        return dtime

    def get_other_bands_info(self):
        other_bands = {
            'kd_490': {
                'file_path': os.path.join(self.path_source, 'trsp.nc')
            }
        }
        return other_bands

    def get_general_array(self, band_name, fill_value):
        other_bands = self.get_other_bands_info()
        if band_name in other_bands:
            file_path = other_bands[band_name]['file_path']
            nc_sat = Dataset(file_path, 'r')
            array = np.ma.array(nc_sat.variables[band_name][:, :])
            array = np.ma.filled(array, fill_value=fill_value)
            return array
        return None

    def get_lat_long_arrays(self):
        coordinates_filename = 'geo_coordinates.nc'
        filepath = os.path.join(self.path_source, coordinates_filename)
        nc_sat = Dataset(filepath, 'r')
        lat = nc_sat.variables['latitude'][:, :]
        lon = nc_sat.variables['longitude'][:, :]

        return lat, lon

    def get_reflectance_band_array(self, wlref, fill_value):
        for band_name in self.reflectance_bands:
            dif = abs(wlref - self.reflectance_bands[band_name]['wavelength'])
            if dif < self.max_dif_wl:
                nc_sat = Dataset(self.reflectance_bands[band_name]['file_path'], 'r')
                array_reflectance = np.ma.array(nc_sat.variables[band_name][:, :])
                #array_reflectance = array_reflectance / np.pi NOT REQUIRED, DATA ARE GIVEN IN RRS
                array_reflectance = np.ma.filled(array_reflectance, fill_value=fill_value)
                return array_reflectance
        return None

    def get_reflectance_band_name(self, wlref):
        for band_name in self.reflectance_bands:
            dif = abs(wlref - self.reflectance_bands[band_name]['wavelength'])
            if dif < self.max_dif_wl:
                return band_name
        return None

    def set_reflectance_bands_mask(self, wlvalues):
        self.reflectance_bands_mask = {}
        for wl in wlvalues:
            band_name = self.get_reflectance_band_name(wl)
            if band_name is not None:
                self.reflectance_bands_mask[band_name] = self.reflectance_bands[band_name]


    def get_val_from_tie_point_grid(self, yPoint, xPoint, ySubsampling, xSubsampling, dataset):
        grid_height = dataset.shape[0]
        grid_width = dataset.shape[1]
        fi = (xPoint + 0.5) / xSubsampling
        fj = (yPoint + 0.5) / ySubsampling
        i0 = self.floor_and_crop(fi, 0, grid_width - 2)
        j0 = self.floor_and_crop(fj, 0, grid_height - 2)
        i1 = i0 + 1
        j1 = j0 + 1
        wi = fi - i0
        wj = fj - j0
        x00 = dataset[j0, i0]
        x10 = dataset[j0, i1]
        x01 = dataset[j1, i0]
        x11 = dataset[j1, i1]
        val = x00 + (wi * (x10 - x00)) + (wj * (x01 - x00)) + (wi * wj * (x11 + x00 - x01 - x10))
        return val

    def floor_and_crop(self, v, minV, maxV):
        rv = math.floor(v)
        if rv < minV:
            return minV
        if rv > maxV:
            return maxV
        return rv

    # name_band: OZA,OAA, SZA, SAA
    def get_angle_array(self, name_band):
        file_path = os.path.join(self.path_source, 'tie_geometries.nc')
        nc_sat = Dataset(file_path, 'r')
        xsubsampling = nc_sat.getncattr('ac_subsampling_factor')
        ysubsampling = nc_sat.getncattr('al_subsampling_factor')

        dataset = nc_sat.variables[name_band][:]

        if self.height == -1 and self.width == -1:
            self.get_dimensions()
        shape = (self.height, self.width)

        array = np.zeros(shape, dtype=np.float64)
        # scale_factor = nc_sat.variables[name_band].scale_factor
        for y in range(0, self.height, ysubsampling):
            if self.verbose and (y == 0 or (y % 1000) == 0):
                print(f'[INFO] Creating angle array: {y}/{self.height}')
            for x in range(0, self.width, xsubsampling):
                xini = x
                xfin = x + xsubsampling
                vali = self.get_val_from_tie_point_grid(y, xini, ysubsampling, xsubsampling, dataset)
                valf = self.get_val_from_tie_point_grid(y, xfin, ysubsampling, xsubsampling, dataset)
                increm = (valf - vali) / xsubsampling
                for ix in range(0, xsubsampling, 1):
                    xp = x + ix
                    if xp < self.width:
                        val = vali + (ix * increm)
                        array[y, xp] = float(val)

        nc_sat.close()

        return array

    def get_angle_info(self, name_band):
        file_path = os.path.join(self.path_source, 'tie_geometries.nc')
        nc_sat = Dataset(file_path, 'r')
        info = {}
        for at in nc_sat.variables[name_band].ncattrs():
            info[at] = nc_sat.variables[name_band].getncattr(at)
        nc_sat.close()
        return info

    def get_oza_angle_array(self):
        array = self.get_angle_array('OZA')
        return array

    def get_oaa_angle_array(self):
        array = self.get_angle_array('OAA')
        return array

    def get_sza_angle_array(self):
        array = self.get_angle_array('SZA')
        return array

    def get_saa_angle_array(self):
        array = self.get_angle_array('SAA')
        return array

    def get_oza_angle_info(self):
        info = self.get_angle_info('OZA')
        return info

    def get_oaa_angle_info(self):
        info = self.get_angle_info('OAA')
        return info

    def get_sza_angle_info(self):
        info = self.get_angle_info('SZA')
        return info

    def get_saa_angle_info(self):
        info = self.get_angle_info('SAA')
        return info

    def get_mask_default(self):
        print(f'[INFO] Creating default mask...')
        file_path = os.path.join(self.path_source, 'wqsf.nc')
        nc_sat = Dataset(file_path, 'r')
        satellite_flag = nc_sat.variables['WQSF']
        flagging = flag.Class_Flags_OLCI(satellite_flag.flag_masks, satellite_flag.flag_meanings)
        flag_list = 'LAND,COASTLINE,CLOUD,CLOUD_AMBIGUOUS,CLOUD_MARGIN,INVALID,COSMETIC,SATURATED,SUSPECT,HISOLZEN,HIGHGLINT,SNOW_ICE,AC_FAIL,WHITECAPS,RWNEG_O2,RWNEG_O3,RWNEG_O4,RWNEG_O5,RWNEG_O6,RWNEG_O7,RWNEG_O8'
        flag_list = flag_list.replace(" ", "")
        flag_list = str.split(flag_list, ',')

        mask_array = np.array(satellite_flag)
        flag_mask = flagging.Mask(mask_array, flag_list)
        flag_mask[np.where(flag_mask != 0)] = 1
        self.width = flag_mask.shape[1]
        self.height = flag_mask.shape[0]
        ntotal = self.width * self.height
        nflagged = np.count_nonzero(flag_mask)
        nvalid = ntotal - nflagged
        # pvalid = (nvalid / ntotal) * 100
        nc_sat.close()
        print(f'[INFO] Number of non-masked pixels: {nvalid}')

        nvalidrrs = ntotal
        for band_name in self.reflectance_bands_mask:
            nc_sat = Dataset(self.reflectance_bands_mask[band_name]['file_path'], 'r')
            array_reflectance = np.ma.array(nc_sat.variables[band_name][:, :])
            nvalid = ntotal - np.ma.count_masked(array_reflectance)
            if nvalid < nvalidrrs:
                nvalidrrs = nvalid
            print(f'[INFO] Valid RRS values for {band_name}: {nvalid}')
            flag_mask[array_reflectance.mask] = 1
            nc_sat.close()

        ##WATER MASKS
        flist = ['WATER']
        water_mask = flagging.Mask(mask_array, flist)
        nwater1 = np.count_nonzero(water_mask)
        print('[INFO] #WATER WITH WATER FLAG: ', nwater1)
        flist = ['INVALID', 'LAND', 'COASTLINE', 'SNOW_ICE']
        land_mask = flagging.Mask(mask_array, flist)
        nwater2 = ntotal - np.count_nonzero(land_mask)
        print('[INFO] #WATER WITH NOT (INVALID + LAND + COASTLINE + SNOW_ICE): ', nwater2)

        nmasked = np.count_nonzero(flag_mask)
        nvalid = ntotal - nmasked
        if nwater2 == 0:
            pvalid = 0
            print(f'[INFO] Number of non-masked pixels: {nvalid} ({pvalid:.2f}%)')
            print(f'----------------------------------------------> {file_path}')
        else:
            pvalid = (nvalid / nwater2) * 100
            print(f'[INFO] Number of non-masked pixels: {nvalid} ({pvalid:.2f}%)')

        # lines = ['Source;Width;Height;NTotal;NFlagged;NWater1;NWater2;NValid;PValid']
        res = [self.name_source, self.width, self.height, ntotal, nflagged, nwater1, nwater2, nvalid, pvalid]
        res_line = f'{self.name_source};{self.width};{self.height};{ntotal};{nflagged};{nwater1};{nwater2};{nvalid};{pvalid}'

        return flag_mask, res, res_line