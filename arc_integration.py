import json
import math

import numpy as np
import os
from arc_mapinfo import ArcMapInfo
from netCDF4 import Dataset
from datetime import datetime as dt
import simplekml


class ArcIntegration():

    def __init__(self, arc_opt, verbose, dir_input, output_type):

        self.verbose = verbose
        self.dir_input = dir_input
        self.ami = ArcMapInfo(arc_opt, False)
        self.width = self.ami.area_def.width
        self.height = self.ami.area_def.height
        self.olci_l2_bands = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75]
        self.olci_l2_min_values = [-0.0063 ,-0.0058 ,-0.0046,-0.0029,-0.0024 ,-0.0017 ,-0.0012,-0.00083,-0.000794,-0.00071,-0.00065]

        self.info = {}
        if output_type is None:
            output_type = 'RRS'
        self.output_type = output_type
        self.time_min = -1
        self.time_max = -1

        self.rrs_variables_all = []
        for wl in self.olci_l2_bands:
            wls = str(wl)
            wls = wls.replace('.', '_')
            bname = f'RRS{wls}'
            self.rrs_variables_all.append(bname)

        self.transp_variables_all = ['KD490_M07']

        if arc_opt is None:  ##only for creating base file
            return

        section = 'INTEGRATE'
        self.arc_integration_method = arc_opt.get_value_param(section, 'method', 'average', 'str')
        self.th_nvalid = arc_opt.get_value_param(section, 'th_nvalid', -1, 'int')
        self.mask_negatives = arc_opt.get_value_param(section, 'mask_negatives', False, 'boolean')
        self.ystep = arc_opt.get_value_param(section, 'ystep', 6500, 'int')
        self.xstep = arc_opt.get_value_param(section, 'xstep', 6500, 'int')
        self.platform = arc_opt.get_value_param(section, 'platform', 'S3', 'str')
        if self.output_type == 'RRS':
            self.average_variables = self.rrs_variables_all
        elif self.output_type == 'TRANSP':
            self.average_variables = self.transp_variables_all
        else:
            rrs_variables = arc_opt.get_value_param(section, 'rrs_bands', self.rrs_variables_all, 'rrslist')
            transp_variables = arc_opt.get_value_param(section, 'transp_bands', self.transp_variables_all, 'strlist')
            self.average_variables = rrs_variables + transp_variables

        if self.verbose:
            print(f'[INFO] Integration method: {self.arc_integration_method}')
            print(f'[INFO] YStep: {self.ystep} XStep: {self.xstep}')
            print(f'[INFO] Platform: {self.platform}')
            print(f'[INFO] Variables: {self.average_variables}')

    def create_nc_file_out(self, ofname):
        if self.verbose:
            print(f'[INFO] Copying file base to start output file...')
        datasetout = self.ami.copy_nc_base(ofname)
        if datasetout is None:
            return datasetout

        ##create rrs variables
        if self.output_type == 'RRS':
            for idx in range(len(self.olci_l2_bands)):
                wl = self.olci_l2_bands[idx]
                wlstr = str(wl).replace('.', '_')
                bandname = f'RRS{wlstr}'
                if bandname in datasetout.variables:
                    continue
                if self.verbose:
                    print(f'[INFO] Creating RRS band: {bandname}')
                var = datasetout.createVariable(bandname, 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
                var[:] = 0
                #var.wavelength = wl
                var.long_name = f'Remote Sensing Reflectance at {bandname.lower()}'
                var.standard_name = f'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air'
                var.units = 'sr^-1'
                var.grid_mapping = 'stereographic'
                var.coordinates = 'longitude latitude'
                var.valid_min = self.olci_l2_min_values[idx]
                var.valid_max = 1.0
                var.type = 'surface'
                var.applied_flags = '(WATER, INLAND_WATER) and not (CLOUD, CLOUD_AMBIGUOUS, CLOUD_MARGIN, INVALID, ' \
                                    'COSMETIC, SATURATED, SUSPECT, HISOLZEN, HIGHGLINT, SNOW_ICE, AC_FAIL, WHITECAPS, ' \
                                    'ADJAC, RWNEG_O2, RWNEG_O3, RWNEG_O4, RWNEG_O5, RWNEG_O6, RWNEG_O7, RWNEG_O8) '
                var.source = 'OLCI - Level2'



        if self.verbose:
            print(f'[INFO] Creating other bands...')

        ##create sum_weights variable
        var = datasetout.createVariable('sum_weights', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = 0

        # create n_granules
        var = datasetout.createVariable('n_granules', 'i4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = 0

        if self.output_type == 'RRS' or self.output_type == 'TRANSP':
            return datasetout

        # time_dif
        var = datasetout.createVariable('time_dif', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = -999

        # time_min
        var = datasetout.createVariable('time_min', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = self.time_max + 1

        # time_max
        var = datasetout.createVariable('time_max', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = self.time_min - 1

        # oza_min
        var = datasetout.createVariable('oza_min', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = 360

        # oza_max
        var = datasetout.createVariable('oza_max', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = 0

        return datasetout

    ##TYPE: RRS, KD, CONFIG
    def make_integration(self, file_out):
        if self.verbose:
            print('[INFO] Checking bands: ')
        for band in self.average_variables:
            if band not in self.average_variables_all:
                print(f'[ERROR] RRS variable {band} is not available. Exiting...')
                return
        if self.arc_integration_method == 'average':
            self.make_integration_avg(file_out)

    # TYPE: RRS, KD, CONFIG
    def make_integration_avg(self, file_out):
        if self.verbose:
            print('[INFO] Retrieving info from granules...')
        self.get_info()
        ngranules = len(self.info)
        if ngranules == 0:
            print(f'[WARNING] No valid granules were found. Check date and platform values. Skiping...')
            return
        if self.verbose:
            print(f'[INFO] Creating ouptput file: {file_out}')
        datasetout = self.create_nc_file_out(file_out)

        var_sensor_mask = datasetout.variables['sensor_mask']
        var_n_granules = datasetout.variables['n_granules']
        var_weighted_mask = datasetout.variables['sum_weights']

        # var_time_dif = datasetout.variables['time_dif']
        # var_time_min = datasetout.variables['time_min']
        # var_time_max = datasetout.variables['time_min']
        # var_min_oza = datasetout.variables['oza_min']
        # var_max_oza = datasetout.variables['oza_max']

        for name in self.info:
            if self.th_nvalid >= 0:
                nvalid = self.info[name]['n_valid']
                if nvalid <= self.th_nvalid:
                    print(
                        f'[INFO] Number of valid pixels {nvalid} must be greater than {self.th_nvalid}. Granule {name} is skipped.')
                    continue

            if self.verbose:
                print(f'[INFO]Working with granule: {name}')
            file = os.path.join(self.dir_input, name)
            yini = self.info[name]['y_min']
            yfin = self.info[name]['y_max']
            xini = self.info[name]['x_min']
            xfin = self.info[name]['x_max']

            dataset = Dataset(file)
            sensor_mask_granule = np.array(dataset.variables['sensor_mask'][yini:yfin, xini:xfin])
            weigthed_mask_granule = np.array(dataset.variables['mask'][yini:yfin, xini:xfin])

            if self.mask_negatives:
                for idx in range(1, 7):
                    var_rrs = self.average_variables_all[idx]
                    var_rrs_array = np.array(dataset.variables[var_rrs][yini:yfin, xini:xfin])
                    weigthed_mask_granule = np.where(var_rrs_array > 0, weigthed_mask_granule, 0)

            sensor_mask = np.array(var_sensor_mask[yini:yfin, xini:xfin])
            ngranules = np.array(var_n_granules[yini:yfin, xini:xfin])
            weigthed_mask = np.array(var_weighted_mask[yini:yfin, xini:xfin])

            indices = np.where(sensor_mask != -999)
            sensor_mask[indices] = sensor_mask[indices] + sensor_mask_granule[indices]

            indices = np.where(weigthed_mask_granule >= 0)
            weigthed_mask[indices] = weigthed_mask[indices] + weigthed_mask_granule[indices]

            ngranules[indices] = ngranules[indices] + 1

            var_sensor_mask[yini:yfin, xini:xfin] = [sensor_mask[:, :]]
            var_n_granules[yini:yfin, xini:xfin] = [ngranules[:, :]]
            var_weighted_mask[yini:yfin, xini:xfin] = [weigthed_mask[:, :]]

            if self.info[name]['n_valid'] == 0:
                if self.verbose:
                    print('[INFO] Granule without valid pixels...')
                dataset.close()
                continue

            for var_avg_name in self.average_variables:
                if self.verbose:
                    print(f'[INFO]--> {var_avg_name}')
                var_avg = datasetout.variables[var_avg_name]
                avg_array = var_avg[yini:yfin, xini:xfin]
                avg_granule = np.array(dataset.variables[var_avg_name][yini:yfin, xini:xfin])
                indices = np.where(weigthed_mask_granule > 0)
                # avg_array[weigthed_mask_granule > 0] = avg_array[weigthed_mask_granule > 0] + (
                #         avg_granule[weigthed_mask_granule > 0] * weigthed_mask_granule[weigthed_mask_granule > 0])
                avg_array[indices] = avg_array[indices] + (avg_granule[indices] * weigthed_mask_granule[indices])
                var_avg[yini:yfin, xini:xfin] = [avg_array[:, :]]
            dataset.close()

        if self.verbose:
            print('[INFO] Computing average...')

        for y in range(0, self.height, self.ystep):
            if self.verbose:
                print(f'[INFO] -> {y}')
            for x in range(0, self.width, self.xstep):
                limits = self.get_limits(y, x, self.ystep, self.xstep, self.height, self.width)
                weigthed_mask = np.array(var_weighted_mask[limits[0]:limits[1], limits[2]:limits[3]])

                if np.max(weigthed_mask[:]) == 0:
                    continue
                indices_good = np.where(weigthed_mask > 0)
                indices_mask = np.where(weigthed_mask == 0)
                for var_avg_name in self.average_variables:
                    var_avg = datasetout.variables[var_avg_name]
                    avg_array = var_avg[limits[0]:limits[1], limits[2]:limits[3]]
                    avg_array[indices_good] = avg_array[indices_good] / weigthed_mask[indices_good]
                    avg_array[indices_mask] = -999
                    var_avg[limits[0]:limits[1], limits[2]:limits[3]] = [avg_array[:, :]]

        datasetout.close()

    def get_overlapping_images(self, name_ref, yini, yfin, xini, xfin, onlyvalid):
        info_over = {}
        for name in self.info:
            if name_ref == name:
                continue
            inside_y = False
            inside_x = False
            if yini <= self.info[name]['y_min'] <= yfin:
                inside_y = True
            if yini <= self.info[name]['y_max'] <= yfin:
                inside_y = True
            if xini <= self.info[name]['x_min'] <= xfin:
                inside_x = True
            if xini <= self.info[name]['x_max'] <= xfin:
                inside_x = True
            if inside_x and inside_y:
                if onlyvalid:
                    if self.info[name]['n_valid'] > 0:
                        info_over[name] = self.info[name]
                else:
                    info_over[name] = self.info[name]

        return info_over

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

    def get_info(self):
        self.info = {}
        for name in os.listdir(self.dir_input):
            if not name.endswith('.nc'):
                continue
            if not name.startswith(self.platform):
                continue

            finput = os.path.join(self.dir_input, name)
            try:
                dataset = Dataset(finput)
            except:
                print(f'[WARNING] File: {finput} is not a valid NetCDF4 dataset. Skipping...')
                continue
            self.info[name] = {
                'rel_orbit': dataset.relative_orbit,
                'y_min': dataset.resampled_ymin,
                'y_max': dataset.resampled_ymax,
                'x_min': dataset.resampled_xmin,
                'x_max': dataset.resampled_xmax,
                'n_total': dataset.resampled_n_total,
                'n_valid': dataset.resampled_n_valid,
                'granule_index': dataset.granule_index,
                'sensor_flag': dataset.sensorflag,
                'start_date': dataset.start_date
            }
            time_here = float(dt.strptime(dataset.start_date, '%Y%m%dT%H%M%S').timestamp())
            if self.time_min == -1 and self.time_max == -1:
                self.time_min = time_here
                self.time_max = time_here
            else:
                if time_here <= self.time_min:
                    self.time_min = time_here
                if time_here >= self.time_max:
                    self.time_max = time_here
            dataset.close()

        # file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/infoImages.json'
        # json.dumps(file_out)
        # kml.save(file_out)

    def compute_nvalid_array(self):
        array = np.zeros((self.height, self.height))
        for name in os.listdir(self.dir_input):
            if not name.endswith('.nc'):
                continue
            print(name)
            finput = os.path.join(self.dir_input, name)
            try:
                dataset = Dataset(finput)
            except:
                print(f'[WARNING] File: {finput} is not a valid NetCDF4 dataset. Skipping...')
                continue
            if not 'mask' in dataset.variables:
                print(f'[WARNING] Dataset {name} does not contain mask variable. Skipping...')
                continue
            mask = np.array(dataset.variables['mask'])
            array = array + mask
            dataset.close()
        return array

    def compute_sensor_mask_array(self, limits):
        if limits is None:
            yini = 0
            yfin = self.height
            xini = 0
            xfin = self.width
        else:
            yini = limits[0]
            yfin = limits[1]
            xini = limits[2]
            xfin = limits[3]
        h = yfin - yini
        w = xfin - xini
        smask = np.zeros((h, w), dtype=np.int)
        sensor_id = 1
        for name in os.listdir(self.dir_input):
            if not name.endswith('.nc'):
                continue
            finput = os.path.join(self.dir_input, name)
            try:
                dataset = Dataset(finput)
            except:
                print(f'[WARNING] File: {finput} is not a valid NetCDF4 dataset. Skipping...')
                continue
            if sensor_id > 10:
                continue
            print(name, sensor_id)

            oza = np.array(dataset.variables['OZA'][yini:yfin, xini:xfin])
            nid = np.count_nonzero(np.where(oza != -999))
            # print(nid)
            if nid == 0:
                sensor_id = sensor_id + 1
                continue
            # smask[np.where(smask > 0) and np.where(oza != -999)] = 1000
            # smask[oza != -999] = smask[oza != -999] + sensor_id
            smask[oza != -999] = sensor_id
            sensor_id = sensor_id + 1

            dataset.close()

        return smask
