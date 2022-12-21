import json
import math

import numpy as np
import os
from arc_mapinfo import ArcMapInfo
from netCDF4 import Dataset
from datetime import datetime as dt
import simplekml


class ArcIntegration():

    def __init__(self, fconfig, verbose, dir_input):

        if fconfig is None:
            fconfig = 'arc_config.ini'
        self.verbose = verbose
        self.dir_input = dir_input
        self.ami = ArcMapInfo(fconfig, False)
        self.width = self.ami.area_def.width
        self.height = self.ami.area_def.height
        self.ystep = 6500
        self.xstep = 6500
        self.olci_l2_bands = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75, 753.75, 778.75]
        self.info = {}
        self.time_min = -1
        self.time_max = -1
        self.average_variables = ['RRS510']

    def create_nc_file_out(self, ofname):
        if self.verbose:
            print(f'[INFO] Copying file base to start output file {ofname} ...')
        datasetout = self.ami.copy_nc_base(ofname)
        if datasetout is None:
            return datasetout

        ##create rrs variables
        for wl in self.olci_l2_bands:
            wlstr = str(wl).replace('.', '_')
            bandname = f'RRS{wlstr}'
            #self.average_variables.append(bandname)
            var = datasetout.createVariable(bandname, 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
            var[:] = 0
            var.wavelength = wl

        ##create sum_weights variable
        var = datasetout.createVariable('sum_weights', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = 0

        # create n_granules
        var = datasetout.createVariable('n_granules', 'i4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = 0

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

    def make_integration_avg(self, file_out):
        if self.verbose:
            print('[INFO] Retrieving info from granules...')
        self.get_info()
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
            if self.verbose:
                print(f'Working with granule: {name}')
            file = os.path.join(self.dir_input, name)
            yini = self.info[name]['y_min']
            yfin = self.info[name]['y_max']
            xini = self.info[name]['x_min']
            xfin = self.info[name]['x_max']
            # sensor_flag = self.info[name]['sensor_flag']
            # h = yfin - yini
            # w = xfin - xini
            # info_overlap = self.get_overlapping_images(name, yini, yfin, xini, xfin)
            dataset = Dataset(file)
            sensor_mask_granule = np.array(dataset.variables['sensor_mask'][yini:yfin, xini:xfin])
            weigthed_mask_granule = np.array(dataset.variables['mask'][yini:yfin, xini:xfin])

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
                dataset.close()
                continue

            for var_avg_name in self.average_variables:
                if self.verbose:
                    print(f'--> {var_avg_name}')
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
            print(f'[INFO] Line-> {y}')
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

        for name in os.listdir(self.dir_input):
            if not name.endswith('.nc'):
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
