import json
import math

import numpy as np
import os
from arc_mapinfo import ArcMapInfo
from netCDF4 import Dataset
from datetime import datetime as dt


class ArcIntegration():

    def __init__(self, arc_opt, verbose, dir_input, output_type, file_attributes):

        self.verbose = verbose
        self.dir_input = dir_input
        self.ami = ArcMapInfo(arc_opt, False)
        self.width = self.ami.area_def.width
        self.height = self.ami.area_def.height
        self.olci_l2_bands = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75]
        self.olci_l2_min_values = [-0.0063, -0.0058, -0.0046, -0.0029, -0.0024, -0.0017, -0.0012, -0.00083, -0.000794,
                                   -0.00071, -0.00065]
        self.file_attributes = file_attributes

        self.apply_pool = 0

        self.info = {}
        if output_type is None:
            output_type = 'RRS'
        self.output_type = output_type
        self.time_min = -1
        self.time_max = -1

        self.rrs_variables_all = {}
        self.rrs_variables_rneg_flag = []
        for idx in range(len(self.olci_l2_bands)):
            wl = self.olci_l2_bands[idx]
            wls = str(wl)
            wls = wls.replace('.', '_')
            bname = f'RRS{wls}'
            self.rrs_variables_all[bname] = {
                'wavelength': wl,
                'min_value': self.olci_l2_min_values[idx],
                'max_value': 1.0,
                'mask': None
            }
            if 410 < wl < 670:
                self.rrs_variables_rneg_flag.append(bname)

        self.transp_variables_all = {
            'KD490': {
                'min_value': 0.0,
                'max_value': 10.0,
                'mask': None
            }
        }

        self.granule_variables = {
            'KD490': 'KD490_M07'
        }

        if self.output_type == 'RRS':
            self.average_variables = list(self.rrs_variables_all.keys())
        elif self.output_type == 'TRANSP':
            self.average_variables = list(self.transp_variables_all)
        elif self.output_type == 'OPERATIVE':
            self.average_variables = list(self.rrs_variables_all.keys()) + list(self.transp_variables_all)

        if arc_opt is None:  ##only for creating base file
            return

        section = 'INTEGRATE'
        self.arc_integration_method = arc_opt.get_value_param(section, 'method', 'average', 'str')
        self.th_nvalid = arc_opt.get_value_param(section, 'th_nvalid', -1, 'int')
        self.mask_negatives = arc_opt.get_value_param(section, 'mask_negatives', False, 'boolean')
        self.ystep = arc_opt.get_value_param(section, 'ystep', 6500, 'int')
        self.xstep = arc_opt.get_value_param(section, 'xstep', 6500, 'int')
        self.platform = arc_opt.get_value_param(section, 'platform', 'S3', 'str')

        if self.output_type == 'TEST':
            rrs_variables = arc_opt.get_value_param(section, 'rrs_bands', list(self.rrs_variables_all.keys()),
                                                    'rrslist')
            transp_variables = arc_opt.get_value_param(section, 'transp_bands', list(self.transp_variables_all.keys()),
                                                       'strlist')
            self.average_variables = rrs_variables + transp_variables

        if self.verbose:
            print(f'[INFO] Integration method: {self.arc_integration_method}')
            print(f'[INFO] YStep: {self.ystep} XStep: {self.xstep}')
            print(f'[INFO] Platform: {self.platform}')
            print(f'[INFO] Variables: {self.average_variables}')

    def get_global_attributes(self, timeliness):
        if self.file_attributes is None:
            return None
        if not os.path.exists(self.file_attributes):
            return None
        import configparser
        try:
            options = configparser.ConfigParser()
            options.read(self.file_attributes)
        except:
            return None
        if not options.has_section('GLOBAL_ATTRIBUTES'):
            return None
        at_dict = dict(options['GLOBAL_ATTRIBUTES'])
        if timeliness is None:
            return at_dict

        sensor = 'Ocean and Land Colour Instrument'
        if 'sensor' in at_dict.keys():
            sensor = at_dict['sensor']

        if self.output_type == 'RRS':
            at_dict['parameter'] = 'Remote Sensing Reflectance'
            at_dict['parameter_code'] = 'RRS'
            if timeliness == 'NR':
                at_dict['timeliness'] = 'NR'
                at_dict['cmems_product_id'] = 'OCEANCOLOUR_ARC_BGC_L3_NRT_009_121'
                if sensor == 'Ocean and Land Colour Instrument':
                    at_dict['title'] = 'cmems_obs-oc_arc_bgc-reflectance_nrt_l3-olci-300m_P1D'
                if sensor == 'ESA Ocean Colour Climate Initiative v6':
                    at_dict['title'] = 'cmems_obs-oc_arc_bgc-reflectance_nrt_l3-multi-4km_P1D'
            if timeliness == 'NT':
                at_dict['timeliness'] = 'NT'
                at_dict['cmems_product_id'] = 'OCEANCOLOUR_ARC_BGC_L3_MY_009_123'
                if sensor == 'Ocean and Land Colour Instrument':
                    at_dict['title'] = 'cmems_obs-oc_arc_bgc-reflectance_my_l3-olci-300m_P1D'
                if sensor == 'ESA Ocean Colour Climate Initiative v6':
                    at_dict['title'] = 'cmems_obs-oc_arc_bgc-reflectance_my_l3-multi-4km_P1D'

        if self.output_type == 'TRANSP':
            at_dict['parameter'] = 'Diffuse attenuation coefficient at 490nm'
            at_dict['parameter_code'] = 'KD490'
            if timeliness == 'NR':
                at_dict['timeliness'] = 'NR'
                at_dict['cmems_product_id'] = 'OCEANCOLOUR_ARC_BGC_L3_NRT_009_121'
                if sensor == 'Ocean and Land Colour Instrument':
                    at_dict['title'] = 'cmems_obs-oc_arc_bgc-transp_nrt_l3-olci-300m_P1D'
                if sensor == 'ESA Ocean Colour Climate Initiative v6':
                    at_dict['title'] = 'cmems_obs-oc_arc_bgc-transp_nrt_l3-multi-4km_P1D'
            if timeliness == 'NT':
                at_dict['timeliness'] = 'NT'
                at_dict['cmems_product_id'] = 'OCEANCOLOUR_ARC_BGC_L3_MY_009_123'
                if sensor == 'Ocean and Land Colour Instrument':
                    at_dict['title'] = 'cmems_obs-oc_arc_bgc-transp_my_l3-olci-300m_P1D'
                if sensor == 'ESA Ocean Colour Climate Initiative v6':
                    at_dict['title'] = 'cmems_obs-oc_arc_bgc-transp_my_l3-multi-4km_P1D'

        return at_dict

    def create_nc_file_out_avg(self, ofname):
        if self.verbose:
            print(f'[INFO] Copying file base {self.ami.ifile_base} to start output file...')
        datasetout = self.ami.copy_nc_base(ofname)
        if datasetout is None:
            return datasetout
        ##create sum_weights variable
        if 'sum_weights' not in datasetout.variables:
            var = datasetout.createVariable('sum_weights', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
            var[:] = 0

        ##create average variable
        if 'average' not in datasetout.variables:
            var = datasetout.createVariable('average', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
            var[:] = 0

        # create n_granules
        if 'n_granules' not in datasetout.variables:
            var = datasetout.createVariable('n_granules', 'i4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
            var[:] = 0

        return datasetout

    def create_nc_file_out(self, ofname, timeliness):
        if self.verbose:
            print(f'[INFO] Copying file base to start output file...')
        datasetout = self.ami.copy_nc_base(ofname)
        if datasetout is None:
            return datasetout

        create_variables = True
        atribs = self.get_global_attributes(timeliness)
        if atribs is not None:  ##atrib is None, atribs ard defined in file base
            for at in atribs:
                if at == 'conventions':
                    datasetout.setncattr('Conventions', atribs[at])
                else:
                    datasetout.setncattr(at, atribs[at])
                if at == 'sensor' and atribs[at] == 'ESA Ocean Colour Climate Initiative v6':
                    create_variables = False

        if not create_variables:
            return datasetout

        ##create rrs variables
        if self.output_type == 'RRS' or self.output_type == 'TEST' or self.output_type == 'OPERATIVE':
            for idx in range(len(self.olci_l2_bands)):
                wl = self.olci_l2_bands[idx]
                wlstr = str(wl).replace('.', '_')
                bandname = f'RRS{wlstr}'
                if bandname in datasetout.variables:
                    continue
                if not bandname in self.average_variables:
                    continue
                if self.verbose:
                    print(f'[INFO] Creating RRS band: {bandname}')
                var = datasetout.createVariable(bandname, 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                                complevel=6)
                # var[:] = 0
                # var.wavelength = wl
                var.long_name = f'Remote Sensing Reflectance at {bandname.lower()}'
                var.standard_name = f'surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air'
                var.units = 'sr^-1'
                var.grid_mapping = 'stereographic'
                var.coordinates = 'lon lat'
                var.valid_min = self.rrs_variables_all[bandname]['min_value']
                var.valid_max = self.rrs_variables_all[bandname]['max_value']
                var.type = 'surface'
                var.missing_value = -999.0
                var.applied_flags = '(WATER, INLAND_WATER) and not (CLOUD, CLOUD_AMBIGUOUS, CLOUD_MARGIN, INVALID, ' \
                                    'COSMETIC, SATURATED, SUSPECT, HISOLZEN, HIGHGLINT, SNOW_ICE, AC_FAIL, WHITECAPS, ' \
                                    'ADJAC, RWNEG_O2, RWNEG_O3, RWNEG_O4, RWNEG_O5, RWNEG_O6, RWNEG_O7, RWNEG_O8) '
                var.source = 'OLCI - Level2'

        bandname = 'KD490'
        addtransp = False
        if bandname in self.average_variables and not bandname in datasetout.variables:
            addtransp = True
        if self.output_type == 'RRS':
            addtransp = False
        if addtransp:
            var = datasetout.createVariable(bandname, 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True, complevel=6)
            # var[:] = 0
            var.band_name = 'OLCI band name KD490_M07'
            var.long_name = 'OLCI Diffuse Attenuation Coefficient at 490nm'
            var.standard_name = 'volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water'
            var.units = 'm^-1'
            var.grid_mapping = 'stereographic'
            var.coordinates = 'lon lat'
            var.type = 'surface'
            var.valid_min = self.transp_variables_all[bandname]['min_value']
            var.valid_max = self.transp_variables_all[bandname]['max_value']
            var.missing_value = -999.0
            var.applied_flags = '(WATER, INLAND_WATER) and not (CLOUD, CLOUD_AMBIGUOUS, CLOUD_MARGIN, INVALID, ' \
                                'COSMETIC, SATURATED, SUSPECT, HISOLZEN, HIGHGLINT, SNOW_ICE, AC_FAIL, WHITECAPS, ' \
                                'ADJAC, RWNEG_O2, RWNEG_O3, RWNEG_O4, RWNEG_O5, RWNEG_O6, RWNEG_O7, RWNEG_O8) '
            var.source = 'OLCI - Level2'

        if self.output_type == 'RRS' or self.output_type == 'TRANSP' or self.output_type == 'OPERATIVE' or self.output_type == 'NONE':
            return datasetout

        if self.verbose:
            print(f'[INFO] Creating other bands...')

        ##create sum_weights variable
        if 'sum_weights' not in datasetout.variables:
            var = datasetout.createVariable('sum_weights', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
            var[:] = 0

        # create n_granules
        if 'n_granules' not in datasetout.variables:
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

    ##TYPE: RRS, TRANSP, CONFIG
    def make_integration(self, output_path):
        if self.verbose:
            print('[INFO] Checking bands: ')
        for band in self.average_variables:
            if band not in self.rrs_variables_all and band not in self.transp_variables_all:
                print(f'[ERROR] Variable {band} is not available. Exiting...')
                return

        if self.arc_integration_method == 'average':
            self.make_integration_avg(output_path)

    def make_integration_avg_deprecate(self, file_out, date_run, timeliness):

        ##DEPRECATED 1, METHOD BASED IN GRANULE GROUPS
        if self.verbose:
            print('[INFO] Retrieving info from granules...')
        self.get_info()
        ngranules = len(self.info)
        if ngranules == 0:
            print(f'[WARNING] No valid granules were found. Check date and platform values. Skipping...')
            return
        if self.verbose:
            print(f'[INFO] Creating ouptput file: {file_out}')
        datasetout = self.create_nc_file_out(file_out, timeliness)

        if date_run is not None:
            datasetout.start_date = date_run.strftime('%Y-%m-%d')
            datasetout.stop_date = date_run.strftime('%Y-%m-%d')
        if timeliness is not None:
            datasetout.timeliness = timeliness
        cdate = dt.utcnow()
        datasetout.creation_date = cdate.strftime('%Y-%m-%d')
        datasetout.creation_time = cdate.strftime('%H:%M:%S UTC')

        var_sensor_mask = datasetout.variables['SENSORMASK']
        var_n_granules = datasetout.variables['n_granules']
        var_weighted_mask = datasetout.variables['sum_weights']

        infog = self.group_granules()
        ngroups = len(infog)
        for group in infog:
            granules = infog[group]['granules_list']
            limits = [infog[group]['y_min'], infog[group]['y_max'], infog[group]['x_min'], infog[group]['x_max']]

            # Step 0: getting general arrays
            sensor_mask = np.array(var_sensor_mask[limits[0]:limits[1], limits[2]:limits[3]])
            ngranules = np.array(var_n_granules[limits[0]:limits[1], limits[2]:limits[3]])
            weigthed_mask = np.array(var_weighted_mask[limits[0]:limits[1], limits[2]:limits[3]])

            # step 1, sum values for weighted mask, ngranules and for each average band. Setting sensor mask
            ngranules_rec = len(granules)
            igranule_rec = 1
            for name in granules:
                if self.verbose:
                    print(f'[INFO] Working with granule: {name} ({igranule_rec}/{ngranules_rec})')
                    igranule_rec = igranule_rec + 1
                file = os.path.join(self.dir_input, name)
                dataset_granule = Dataset(file)
                weigthed_mask_granule = np.array(
                    dataset_granule.variables['mask'][limits[0]:limits[1], limits[2]:limits[3]])

                # assuring that pixels lower than 65 degress are masked
                weigthed_mask_granule[sensor_mask == -999] = -999

                # assuring that pixels lower than valid_min in central bands are masked
                for var_rrs_name in self.rrs_variables_rneg_flag:
                    if not var_rrs_name in self.average_variables:
                        continue
                    min_value = self.rrs_variables_all[var_rrs_name]['min_value']
                    max_value = self.rrs_variables_all[var_rrs_name]['max_value']
                    var_rrs_array_granule = np.array(
                        dataset_granule.variables[var_rrs_name][limits[0]:limits[1], limits[2]:limits[3]])
                    weigthed_mask_granule = np.where(
                        np.logical_and(var_rrs_array_granule >= min_value, var_rrs_array_granule <= max_value),
                        weigthed_mask_granule, 0)

                if self.mask_negatives:  ##Option to mask all the negative reflectances (not used in operational processing)
                    for var_rrs_name in self.rrs_variables_rneg_flag:
                        if not var_rrs_name in self.average_variables:
                            continue
                        var_rrs_array_granule = np.array(
                            dataset_granule.variables[var_rrs_name][limits[0]:limits[1], limits[2]:limits[3]])
                        weigthed_mask_granule = np.where(var_rrs_array_granule > 0, weigthed_mask_granule, 0)

                # ngranules, only for testing
                ngranules[weigthed_mask_granule >= 0] = ngranules[weigthed_mask_granule >= 0] + 1

                # weighted_mask (sum of n valid pixels)
                weigthed_mask[weigthed_mask_granule > 0] = weigthed_mask[weigthed_mask_granule > 0] + \
                                                           weigthed_mask_granule[weigthed_mask_granule > 0]

                # sensor_mask: 1 ,2, 3
                if name.startswith('S3A'):
                    sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 0)] = 1
                    sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 2)] = 3
                if name.startswith('S3B'):
                    sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 0)] = 2
                    sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 1)] = 3

                # assigning general masks
                var_sensor_mask[limits[0]:limits[1], limits[2]:limits[3]] = [sensor_mask[:, :]]
                var_n_granules[limits[0]:limits[1], limits[2]:limits[3]] = [ngranules[:, :]]
                var_weighted_mask[limits[0]:limits[1], limits[2]:limits[3]] = [weigthed_mask[:, :]]

                for var_avg_name in self.average_variables:
                    if self.verbose:
                        print(f'[INFO]--> {var_avg_name}')

                    # destination var
                    var_avg = datasetout.variables[var_avg_name]
                    avg_array = np.array(var_avg[limits[0]:limits[1], limits[2]:limits[3]])

                    # origin var
                    var_granule = var_avg_name
                    if var_avg_name in self.granule_variables.keys():
                        var_granule = self.granule_variables[var_avg_name]
                    avg_granule = np.array(
                        dataset_granule.variables[var_granule][limits[0]:limits[1], limits[2]:limits[3]])
                    if var_granule == 'KD490_M07':
                        avg_granule[avg_granule != -999] = np.power(10, avg_granule[avg_granule != -999])

                    mask_granule = weigthed_mask_granule
                    # for varibles out of the central band, independent mask
                    if var_avg_name not in self.rrs_variables_rneg_flag:
                        if var_avg_name in self.rrs_variables_all:
                            min_value = self.rrs_variables_all[var_avg_name]['min_value']
                            max_value = self.rrs_variables_all[var_avg_name]['max_value']
                            mask_all = self.rrs_variables_all[var_avg_name]['mask']
                        if var_avg_name in self.transp_variables_all:
                            min_value = self.transp_variables_all[var_avg_name]['min_value']
                            max_value = self.transp_variables_all[var_avg_name]['max_value']
                            mask_all = self.transp_variables_all[var_avg_name]['mask']
                        mask_granule = np.where(np.logical_and(avg_granule >= min_value, avg_granule <= max_value),
                                                mask_granule, 0)

                        if mask_all is None:
                            mask_all = np.zeros((self.height, self.width))
                        mask_all_here = mask_all[limits[0]:limits[1], limits[2]:limits[3]]
                        mask_all_here[sensor_mask == -999] = -999
                        mask_all_here[mask_granule > 0] = mask_all_here[mask_granule > 0] + mask_granule[
                            mask_granule > 0]
                        mask_all[limits[0]:limits[1], limits[2]:limits[3]] = mask_all_here[:, :]

                        if var_avg_name in self.rrs_variables_all:
                            self.rrs_variables_all[var_avg_name]['mask'] = mask_all
                        if var_avg_name in self.transp_variables_all:
                            self.transp_variables_all[var_avg_name]['mask'] = mask_all

                    indices = np.where(mask_granule > 0)
                    # making sum in destination var.
                    avg_array = np.where(np.logical_and(mask_granule > 0, avg_array == -999), 0, avg_array)
                    avg_array[indices] = avg_array[indices] + (avg_granule[indices] * mask_granule[indices])
                    var_avg[limits[0]:limits[1], limits[2]:limits[3]] = [avg_array[:, :]]

                dataset_granule.close()

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

                for var_avg_name in self.average_variables:
                    var_avg = datasetout.variables[var_avg_name]
                    avg_array = np.array(var_avg[limits[0]:limits[1], limits[2]:limits[3]])
                    mask = weigthed_mask
                    if var_avg_name not in self.rrs_variables_rneg_flag:
                        if var_avg_name in self.rrs_variables_all:
                            mask_o = self.rrs_variables_all[var_avg_name]['mask']
                            mask = mask_o[limits[0]:limits[1], limits[2]:limits[3]]
                        if var_avg_name in self.transp_variables_all:
                            mask_o = self.transp_variables_all[var_avg_name]['mask']
                            mask = mask_o[limits[0]:limits[1], limits[2]:limits[3]]
                    indices_good = np.where(mask > 0)
                    indices_mask = np.where(mask <= 0)
                    avg_array[indices_good] = avg_array[indices_good] / weigthed_mask[indices_good]
                    avg_array[indices_mask] = -999

                    var_avg[limits[0]:limits[1], limits[2]:limits[3]] = [avg_array[:, :]]

        datasetout.close()

        ##DEPRECATED 2, METHOD BASED IN AREAS
        # for y in range(0, self.height, self.ystep):
        #     for x in range(0, self.width, self.xstep):
        #         if self.verbose:
        #             print(f'[INFO] -> {y} <-> {x} --------------------------------------------------------------------')
        #         limits = self.get_limits(y, x, self.ystep, self.xstep, self.height, self.width)
        #
        #         # step 1, getting general arrays
        #         sensor_mask = np.array(var_sensor_mask[limits[0]:limits[1], limits[2]:limits[3]])
        #         ngranules = np.array(var_n_granules[limits[0]:limits[1], limits[2]:limits[3]])
        #         weigthed_mask = np.array(var_weighted_mask[limits[0]:limits[1], limits[2]:limits[3]])
        #
        #         # step 1, sum values for weighted mask, ngranules and for each average band. Setting sensor mask
        #         info_rec = self.get_rec_info(limits)
        #         ngranules_rec = len(info_rec)
        #         igranule_rec = 1
        #         for name in info_rec:
        #             if self.verbose:
        #                 print(f'[INFO] Working with granule: {name} ({igranule_rec}/{ngranules_rec})')
        #                 igranule_rec = igranule_rec + 1
        #             file = os.path.join(self.dir_input, name)
        #             dataset = Dataset(file)
        #             weigthed_mask_granule = np.array(dataset.variables['mask'][limits[0]:limits[1], limits[2]:limits[3]])
        #             # assuring that pixels lower than 65 degress are masked
        #             weigthed_mask_granule[sensor_mask == -999] = -999
        #             if self.mask_negatives:  ##Option to mask all the negative reflectances (not used in operational processing)
        #                 for idx in range(1, 7):
        #                     var_rrs = self.rrs_variables_all[idx]
        #                     if var_rrs in self.average_variables:
        #                         var_rrs_array = np.array(dataset.variables[var_rrs][limits[0]:limits[1], limits[2]:limits[3]])
        #                         weigthed_mask_granule = np.where(var_rrs_array > 0, weigthed_mask_granule, 0)
        #
        #             indices = np.where(weigthed_mask_granule >= 0)
        #             ngranules[indices] = ngranules[indices] + 1
        #
        #             indices = np.where(weigthed_mask_granule > 0)
        #
        #             if len(indices[0])==0:
        #                 continue
        #
        #             weigthed_mask[indices] = weigthed_mask[indices] + weigthed_mask_granule[indices]
        #             if name.startswith('S3A'):
        #                 sensor_mask[indices] = sensor_mask[indices] + 1
        #             if name.startswith('S3B'):
        #                 sensor_mask[indices] = sensor_mask[indices] + 2
        #
        #             for var_avg_name in self.average_variables:
        #                 if self.verbose:
        #                     print(f'[INFO]--> {var_avg_name}')
        #                 # destination var
        #                 var_avg = datasetout.variables[var_avg_name]
        #                 avg_array = np.array(var_avg[limits[0]:limits[1], limits[2]:limits[3]])
        #                 # origin var
        #                 var_granule = var_avg_name
        #                 if var_avg_name in self.granule_variables.keys():
        #                     var_granule = self.granule_variables[var_avg_name]
        #                 avg_granule = np.array(dataset.variables[var_granule][limits[0]:limits[1], limits[2]:limits[3]])
        #                 if var_granule == 'KD490_M07':
        #                     avg_granule[avg_granule != -999] = np.power(10, avg_granule[avg_granule != -999])
        #
        #                 # making sum in destination var
        #                 avg_array = np.where(np.logical_and(weigthed_mask_granule > 0, avg_array == -999), 0, avg_array)
        #                 avg_array[indices] = avg_array[indices] + (
        #                             avg_granule[indices] * weigthed_mask_granule[indices])
        #                 var_avg[limits[0]:limits[1], limits[2]:limits[3]] = [avg_array[:, :]]
        #
        #             dataset.close()
        #
        #
        #
        #         #STEP 2: Averaging
        #         if np.max(weigthed_mask[:]) == 0: #nothing to average
        #             continue
        #
        #         mask_modified = False
        #         if self.verbose:
        #             print('[INFO] Computing average...')
        #         for var_avg_name in self.average_variables:
        #             if self.verbose:
        #                 print(f'[INFO]--> {var_avg_name}')
        #             ##As mask it's modified according to min/max valid values, indices_good and indices_mask are obtained each time
        #             indices_good = np.where(weigthed_mask > 0)
        #             indices_mask = np.where(weigthed_mask == 0)
        #             var_avg = datasetout.variables[var_avg_name]
        #             avg_array = np.array(var_avg[limits[0]:limits[1], limits[2]:limits[3]])
        #             avg_array[indices_good] = avg_array[indices_good] / weigthed_mask[indices_good]
        #             avg_array[indices_mask] = -999
        #             # checking again min/max valid values
        #             valid_min = var_avg.valid_min
        #             valid_max = var_avg.valid_max
        #             indices_novalid_min = np.where(np.logical_and(avg_array < valid_min, weigthed_mask > 0))
        #             if len(indices_novalid_min[0]) > 0:
        #                 avg_array[indices_novalid_min] = -999
        #                 weigthed_mask[indices_novalid_min] = 0
        #                 mask_modified = True
        #             indices_novalid_max = np.where(avg_array > valid_max)
        #             if len(indices_novalid_max[0]) > 0:
        #                 avg_array[indices_novalid_max] = -999
        #                 weigthed_mask[indices_novalid_max] = 0
        #                 mask_modified = True
        #             var_avg[limits[0]:limits[1], limits[2]:limits[3]] = [avg_array[:, :]]
        #
        #         if mask_modified:
        #             if self.verbose:
        #                 print('Modifying mask...')
        #             for var_avg_name in self.average_variables:
        #                 var_avg = datasetout.variables[var_avg_name]
        #                 avg_array = np.array(var_avg[limits[0]:limits[1], limits[2]:limits[3]])
        #                 avg_array[weigthed_mask == 0] = -999
        #                 var_avg[limits[0]:limits[1], limits[2]:limits[3]] = [avg_array[:, :]]
        #
        #         var_sensor_mask[limits[0]:limits[1], limits[2]:limits[3]] = [sensor_mask[:, :]]
        #         var_n_granules[limits[0]:limits[1], limits[2]:limits[3]] = [ngranules[:, :]]
        #         var_weighted_mask[limits[0]:limits[1], limits[2]:limits[3]] = [weigthed_mask[:, :]]
        #
        # datasetout.close()

        # var_time_dif = datasetout.variables['time_dif']
        # var_time_min = datasetout.variables['time_min']
        # var_time_max = datasetout.variables['time_min']
        # var_min_oza = datasetout.variables['oza_min']
        # var_max_oza = datasetout.variables['oza_max']

    def make_integration_avg(self, output_path):
        if self.verbose:
            print('[INFO] Retrieving info from granules...')
        self.get_info()
        nvalidgranules = len(self.info)
        if nvalidgranules == 0:
            print(f'[WARNING] No valid granules were found. Check date and platform values. Skipping...')
            return
        infof = self.filter_granules()
        nvalidgranules = len(infof)
        if nvalidgranules == 0:
            print(f'[WARNING] No valid granules were found. Check date and platform values. Skipping...')
            return
        if self.verbose:
            print(f'[INFO] Number of valid granules to be averaged: {nvalidgranules}')

        file_mask = os.path.join(output_path, 'Mask.nc')

        ##AVERAGE MASK
        self.create_basic_mask(file_mask, infof, nvalidgranules)

        # CHECKING MASK
        dataset_mask = Dataset(file_mask, 'r+')
        var_mask = dataset_mask.variables['sum_weights']
        var_mask_array = np.array(var_mask)
        ngood = np.count_nonzero(var_mask_array > 0)
        if self.verbose:
            print(f'[INFO] # Number of good pixels (flag-based mask) First check: {ngood}')
        dataset_mask.close()

        if ngood == 0:
            print(f'[WARNING] No valid pixels were found. All the pixels are masked')
            return

        # AVERAGE VARIABLES
        if self.apply_pool == 0:
            for var_avg_name in self.average_variables:
                file_var = os.path.join(output_path, f'{var_avg_name}.nc')
                self.create_single_average(var_avg_name, file_var, infof, nvalidgranules)
        else:
            from multiprocessing import Pool
            params_list = []
            if self.apply_pool < 0:
                poolhere = Pool()
            else:
                # POOL THRESHOLD
                if ngood > 2000000:
                    self.apply_pool = 7
                poolhere = Pool(self.apply_pool)

            for var_avg_name in self.average_variables:
                file_var = os.path.join(output_path, f'{var_avg_name}.nc')

                params = [var_avg_name, file_var, infof, nvalidgranules]
                params_list.append(params)
                if len(params_list) == self.apply_pool:
                    poolhere.map(self.create_single_average_parallel, params_list)
                    params_list = []
            if len(params_list) > 0:
                poolhere.map(self.create_single_average_parallel, params_list)

        # CHECKING MASK
        dataset_mask = Dataset(file_mask, 'r+')
        var_mask = dataset_mask.variables['sum_weights']
        var_mask_array = np.array(var_mask)
        ngood = np.count_nonzero(var_mask_array > 0)
        if self.verbose:
            print(f'[INFO] # Number of good pixels (flag-based mask): {ngood}')
        for var_avg_name in self.average_variables:
            if var_avg_name in self.rrs_variables_rneg_flag:
                if self.verbose:
                    print(f'[INFO] Checking range in variable: {var_avg_name}')
                file_var = os.path.join(output_path, f'{var_avg_name}.nc')
                dataset_var = Dataset(file_var)
                var_avg_array = np.array(dataset_var.variables['average'])
                var_mask_array[var_avg_array < 0] = 0
                dataset_var.close()
        # var_mask_array[var_mask_array==0] = -999
        var_mask[:] = [var_mask_array[:]]
        ngood = np.count_nonzero(var_mask_array > 0)
        if self.verbose:
            print(f'[INFO] # Number of good pixels (flag-based mask + range): {ngood}')
        dataset_mask.close()

    def start_rrs_or_transp_file(self, output_path, file_out, date_run, timeliness):
        if self.verbose:
            print(f'[INFO] Creating ouptput file: {file_out}')
        datasetout = self.create_nc_file_out(file_out, timeliness)

        if date_run is not None:
            datasetout.start_date = date_run.strftime('%Y-%m-%d')
            datasetout.stop_date = date_run.strftime('%Y-%m-%d')
            if 'time' in datasetout.variables:
                timeseconds = (date_run - dt(1981, 1, 1, 0, 0, 0)).total_seconds()
                datasetout.variables['time'][0] = [np.int32(timeseconds)]
        if timeliness is not None:
            datasetout.timeliness = timeliness
        cdate = dt.utcnow()
        datasetout.creation_date = cdate.strftime('%Y-%m-%d')
        datasetout.creation_time = cdate.strftime('%H:%M:%S UTC')

        ##SENSOR MASK
        var_sensor_mask = datasetout.variables['SENSORMASK']
        file_mask = os.path.join(output_path, 'Mask.nc')

        # Mask.nc only exist if there are data valid
        wmask = None
        if os.path.exists(file_mask):
            dataset_mask = Dataset(file_mask)
            smask = np.array(dataset_mask.variables['SENSORMASK'])
            wmask = np.array(dataset_mask.variables['sum_weights'])
            var_sensor_mask[:] = [smask[:]]
            dataset_mask.close()

        return datasetout, wmask

    def create_rrs_file(self, output_path, file_out, date_run, timeliness, correctRrs):
        datasetout, wmask = self.start_rrs_or_transp_file(output_path, file_out, date_run, timeliness)

        # RRS BANDS
        for avg_name in self.rrs_variables_all:
            file_avg = os.path.join(output_path, f'{avg_name}.nc')
            if os.path.exists(file_avg) and avg_name in datasetout.variables:
                if self.verbose:
                    print(f'[INFO] Assigning variable: {avg_name}')
                variable = datasetout.variables[avg_name]
                dataset_var = Dataset(file_avg)
                var_array = np.array(dataset_var.variables['average'])
                if avg_name in self.rrs_variables_rneg_flag:
                    if self.verbose:
                        print(f'[INFO] Adapting mask for variable {avg_name}')
                    var_array[wmask <= 0] = -999.0
                if correctRrs:
                    var_array[var_array != -999] = var_array[var_array != -999] / np.pi
                variable[0, :, :] = [var_array[:, :]]
                dataset_var.close()

        datasetout.close()
        if self.verbose:
            print(f'[INFO] Completed')

    def create_transp_file(self, output_path, file_out, date_run, timeliness):
        datasetout, wmask = self.start_rrs_or_transp_file(output_path, file_out, date_run, timeliness)

        # TRANSP BANDS
        for avg_name in self.transp_variables_all:
            file_avg = os.path.join(output_path, f'{avg_name}.nc')
            if os.path.exists(file_avg) and avg_name in datasetout.variables:
                if self.verbose:
                    print(f'[INFO] Assigning variable: {avg_name}')
                variable = datasetout.variables[avg_name]
                dataset_var = Dataset(file_avg)
                var_array = np.array(dataset_var.variables['average'])
                variable[0, :, :] = [var_array[:, :]]
                dataset_var.close()
        datasetout.close()
        if self.verbose:
            print(f'[INFO] Completed')

    def create_single_average_parallel(self, params):
        self.create_single_average(params[0], params[1], params[2], params[3])

    def create_single_average(self, var_avg_name, file_var, infof, nvalidgranules):
        if self.verbose:
            print(f'[INFO] Creating basic average file: {file_var}')

        if os.path.exists(file_var):
            os.remove(file_var)

        min_value = 0.0
        max_value = 1.0
        if var_avg_name in self.rrs_variables_all:
            min_value = self.rrs_variables_all[var_avg_name]['min_value']
            max_value = self.rrs_variables_all[var_avg_name]['max_value']
        if var_avg_name in self.transp_variables_all:
            min_value = self.transp_variables_all[var_avg_name]['min_value']
            max_value = self.transp_variables_all[var_avg_name]['max_value']

        if self.verbose:
            print(f'[INFO] Minimal value: {min_value}')
            print(f'[INFO] Maximum value: {max_value}')

        datasetout = self.create_nc_file_out_avg(file_var)
        var_sensor_mask = datasetout.variables['SENSORMASK']
        var_sum = datasetout.variables['average']
        var_num = datasetout.variables['sum_weights']

        igranule = 1
        for name in infof:
            if self.verbose:
                print(f'[INFO] Working with granule: {name} ({igranule}/{nvalidgranules})')
                igranule = igranule + 1

            file = os.path.join(self.dir_input, name)
            dataset_granule = Dataset(file)
            yini = self.info[name]['y_min']
            yfin = self.info[name]['y_max']
            xini = self.info[name]['x_min']
            xfin = self.info[name]['x_max']

            # destination var
            sum_array = np.array(var_sum[yini:yfin, xini:xfin])
            num_array = np.array(var_num[yini:yfin, xini:xfin])

            # sensor mask
            sensor_mask_overall = np.array(var_sensor_mask[yini:yfin, xini:xfin])

            # origin var
            var_granule = var_avg_name
            if var_avg_name in self.granule_variables.keys():
                var_granule = self.granule_variables[var_avg_name]
            avg_granule = np.array(dataset_granule.variables[var_granule][yini:yfin, xini:xfin])
            if var_granule == 'KD490_M07':
                avg_granule[avg_granule != -999] = np.power(10, avg_granule[avg_granule != -999])

            # origin mask
            weigthed_mask_granule = np.array(dataset_granule.variables['mask'][yini:yfin, xini:xfin])

            # assuring that pixels lower than 65 degress are masked
            weigthed_mask_granule[sensor_mask_overall == -999] = -999

            # average is computed only for valid min/max values
            weigthed_mask_granule = np.where(np.logical_and(avg_granule >= min_value, avg_granule <= max_value),
                                             weigthed_mask_granule, 0)

            # computing sum
            indices = np.where(weigthed_mask_granule > 0)
            sum_array[indices] = sum_array[indices] + (avg_granule[indices] * weigthed_mask_granule[indices])
            num_array[indices] = num_array[indices] + weigthed_mask_granule[indices]
            var_sum[yini:yfin, xini:xfin] = [sum_array]
            var_num[yini:yfin, xini:xfin] = [num_array]

            dataset_granule.close()

        # print('ystep: ', self.ystep, ' xstep: ', self.xstep)
        for y in range(0, self.height, self.ystep):
            if self.verbose:
                print(f'[INFO] -> {y}')
            for x in range(0, self.width, self.xstep):
                limits = self.get_limits(y, x, self.ystep, self.xstep, self.height, self.width)
                sum_array = np.array(var_sum[limits[0]:limits[1], limits[2]:limits[3]])
                if np.max(sum_array[:]) == 0:
                    sum_array[:] = -999
                    var_sum[limits[0]:limits[1], limits[2]:limits[3]] = [sum_array[:, :]]
                    continue
                num_array = np.array(var_num[limits[0]:limits[1], limits[2]:limits[3]])
                indices_good = np.where(num_array > 0)
                indices_mask = np.where(num_array <= 0)
                sum_array[indices_good] = sum_array[indices_good] / num_array[indices_good]
                sum_array[indices_mask] = -999
                sum_array[np.logical_or(sum_array < min_value, sum_array > max_value)] = -999

                var_sum[limits[0]:limits[1], limits[2]:limits[3]] = [sum_array[:, :]]

        datasetout.close()

    def create_basic_mask(self, file_mask, infof, nvalidgranules):
        if self.verbose:
            print(f'[INFO] Creating basic mask file: {file_mask}')

        if os.path.exists(file_mask):
            os.remove(file_mask)

        datasetout = self.create_nc_file_out_avg(file_mask)

        var_sensor_mask = datasetout.variables['SENSORMASK']
        var_n_granules = datasetout.variables['n_granules']
        var_weighted_mask = datasetout.variables['sum_weights']
        igranule = 1
        for name in infof:
            if self.verbose:
                print(f'[INFO] Working with granule: {name} ({igranule}/{nvalidgranules})')
                igranule = igranule + 1

            file = os.path.join(self.dir_input, name)
            yini = self.info[name]['y_min']
            yfin = self.info[name]['y_max']
            xini = self.info[name]['x_min']
            xfin = self.info[name]['x_max']

            # general variables
            sensor_mask = np.array(var_sensor_mask[yini:yfin, xini:xfin])
            ngranules = np.array(var_n_granules[yini:yfin, xini:xfin])
            weigthed_mask = np.array(var_weighted_mask[yini:yfin, xini:xfin])

            # weighted mask granule
            dataset_granule = Dataset(file)
            weigthed_mask_granule = np.array(dataset_granule.variables['mask'][yini:yfin, xini:xfin])

            # assuring that pixels lower than 65 degress are masked
            weigthed_mask_granule[sensor_mask == -999] = -999

            ngranules[weigthed_mask_granule >= 0] = ngranules[weigthed_mask_granule >= 0] + 1

            # indices = np.where(weigthed_mask_granule > 0)

            weigthed_mask[weigthed_mask_granule > 0] = weigthed_mask[weigthed_mask_granule > 0] + weigthed_mask_granule[
                weigthed_mask_granule > 0]
            if name.startswith('S3A'):
                sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 0)] = 1
                sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 2)] = 3
            if name.startswith('S3B'):
                sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 0)] = 2
                sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 1)] = 3

            var_sensor_mask[yini:yfin, xini:xfin] = [sensor_mask[:, :]]
            var_n_granules[yini:yfin, xini:xfin] = [ngranules[:, :]]
            var_weighted_mask[yini:yfin, xini:xfin] = [weigthed_mask[:, :]]

        datasetout.close()

    # TYPE: RRS, TRANSP, OPERATIVE, TEST
    def make_integration_avg_v1(self, file_out, date_run, timeliness):
        if self.verbose:
            print('[INFO] Retrieving info from granules...')
        self.get_info()
        nvalidgranules = len(self.info)
        if nvalidgranules == 0:
            print(f'[WARNING] No valid granules were found. Check date and platform values. Skipping...')
            return
        infof = self.filter_granules()
        nvalidgranules = len(infof)
        if nvalidgranules == 0:
            print(f'[WARNING] No valid granules were found. Check date and platform values. Skipping...')
            return
        if self.verbose:
            print(f'[INFO] Number of valid granules to be averaged: {nvalidgranules}')

        if self.verbose:
            print(f'[INFO] Creating ouptput file: {file_out}')
        datasetout = self.create_nc_file_out(file_out, timeliness)

        if date_run is not None:
            datasetout.start_date = date_run.strftime('%Y-%m-%d')
            datasetout.stop_date = date_run.strftime('%Y-%m-%d')
        if timeliness is not None:
            datasetout.timeliness = timeliness
        cdate = dt.utcnow()
        datasetout.creation_date = cdate.strftime('%Y-%m-%d')
        datasetout.creation_time = cdate.strftime('%H:%M:%S UTC')

        var_sensor_mask = datasetout.variables['SENSORMASK']
        var_n_granules = datasetout.variables['n_granules']
        var_weighted_mask = datasetout.variables['sum_weights']

        igranule = 1
        for name in infof:
            if self.verbose:
                print(f'[INFO] Working with granule: {name} ({igranule}/{nvalidgranules})')
                igranule = igranule + 1

            file = os.path.join(self.dir_input, name)
            yini = self.info[name]['y_min']
            yfin = self.info[name]['y_max']
            xini = self.info[name]['x_min']
            xfin = self.info[name]['x_max']

            # general variables
            sensor_mask = np.array(var_sensor_mask[yini:yfin, xini:xfin])
            ngranules = np.array(var_n_granules[yini:yfin, xini:xfin])
            weigthed_mask = np.array(var_weighted_mask[yini:yfin, xini:xfin])

            # weighted mask granule
            dataset_granule = Dataset(file)
            weigthed_mask_granule = np.array(dataset_granule.variables['mask'][yini:yfin, xini:xfin])

            # assuring that pixels lower than 65 degress are masked
            weigthed_mask_granule[sensor_mask == -999] = -999

            # assuring that pixels lower than valid_min in central bands are masked
            for var_rrs_name in self.rrs_variables_rneg_flag:
                if not var_rrs_name in self.average_variables:
                    continue
                min_value = self.rrs_variables_all[var_rrs_name]['min_value']
                max_value = self.rrs_variables_all[var_rrs_name]['max_value']
                var_rrs_array_granule = np.array(dataset_granule.variables[var_rrs_name][yini:yfin, xini:xfin])
                weigthed_mask_granule = np.where(
                    np.logical_and(var_rrs_array_granule >= min_value, var_rrs_array_granule <= max_value),
                    weigthed_mask_granule, 0)

            if self.mask_negatives:  ##Option to mask all the negative reflectances (not used in operational processing)
                for var_rrs_name in self.rrs_variables_rneg_flag:
                    if not var_rrs_name in self.average_variables:
                        continue
                    min_value = self.rrs_variables_all[var_rrs_name]['min_value']
                    var_rrs_array_granule = np.array(dataset_granule.variables[var_rrs_name][yini:yfin, xini:xfin])
                    weigthed_mask_granule = np.where(var_rrs_array_granule > 0, weigthed_mask_granule, 0)

            # ngranules, only for testing
            ngranules[weigthed_mask_granule >= 0] = ngranules[weigthed_mask_granule >= 0] + 1

            # indices = np.where(weigthed_mask_granule > 0)

            weigthed_mask[weigthed_mask_granule > 0] = weigthed_mask[weigthed_mask_granule > 0] + weigthed_mask_granule[
                weigthed_mask_granule > 0]
            if name.startswith('S3A'):
                sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 0)] = 1
                sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 2)] = 3
            if name.startswith('S3B'):
                sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 0)] = 2
                sensor_mask[np.logical_and(weigthed_mask_granule > 0, sensor_mask == 1)] = 3

            var_sensor_mask[yini:yfin, xini:xfin] = [sensor_mask[:, :]]
            var_n_granules[yini:yfin, xini:xfin] = [ngranules[:, :]]
            var_weighted_mask[yini:yfin, xini:xfin] = [weigthed_mask[:, :]]

            for var_avg_name in self.average_variables:
                if self.verbose:
                    print(f'[INFO]--> {var_avg_name}')

                # destination var
                var_avg = datasetout.variables[var_avg_name]
                avg_array = np.array(var_avg[yini:yfin, xini:xfin])

                # origin var
                var_granule = var_avg_name
                if var_avg_name in self.granule_variables.keys():
                    var_granule = self.granule_variables[var_avg_name]
                avg_granule = np.array(dataset_granule.variables[var_granule][yini:yfin, xini:xfin])
                if var_granule == 'KD490_M07':
                    avg_granule[avg_granule != -999] = np.power(10, avg_granule[avg_granule != -999])

                mask_granule = weigthed_mask_granule
                # for varibles out of the central band, independent mask
                if var_avg_name not in self.rrs_variables_rneg_flag:
                    if var_avg_name in self.rrs_variables_all:
                        min_value = self.rrs_variables_all[var_avg_name]['min_value']
                        max_value = self.rrs_variables_all[var_avg_name]['max_value']
                        mask_all = self.rrs_variables_all[var_avg_name]['mask']
                    if var_avg_name in self.transp_variables_all:
                        min_value = self.transp_variables_all[var_avg_name]['min_value']
                        max_value = self.transp_variables_all[var_avg_name]['max_value']
                        mask_all = self.transp_variables_all[var_avg_name]['mask']
                    mask_granule = np.where(np.logical_and(avg_granule >= min_value, avg_granule <= max_value),
                                            mask_granule, 0)

                    if mask_all is None:
                        mask_all = np.zeros((self.height, self.width))
                    mask_all_here = mask_all[yini:yfin, xini:xfin]
                    mask_all_here[sensor_mask == -999] = -999
                    mask_all_here[mask_granule > 0] = mask_all_here[mask_granule > 0] + mask_granule[mask_granule > 0]
                    mask_all[yini:yfin, xini:xfin] = mask_all_here[:, :]

                    if var_avg_name in self.rrs_variables_all:
                        self.rrs_variables_all[var_avg_name]['mask'] = mask_all
                    if var_avg_name in self.transp_variables_all:
                        self.transp_variables_all[var_avg_name]['mask'] = mask_all

                indices = np.where(mask_granule > 0)
                # making sum in destination var.
                avg_array = np.where(np.logical_and(mask_granule > 0, avg_array == -999), 0, avg_array)
                avg_array[indices] = avg_array[indices] + (avg_granule[indices] * mask_granule[indices])
                var_avg[yini:yfin, xini:xfin] = [avg_array[:, :]]

            dataset_granule.close()

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

                for var_avg_name in self.average_variables:
                    var_avg = datasetout.variables[var_avg_name]
                    avg_array = np.array(var_avg[limits[0]:limits[1], limits[2]:limits[3]])
                    mask = weigthed_mask
                    if var_avg_name not in self.rrs_variables_rneg_flag:
                        if var_avg_name in self.rrs_variables_all:
                            mask_o = self.rrs_variables_all[var_avg_name]['mask']
                            mask = mask_o[limits[0]:limits[1], limits[2]:limits[3]]
                        if var_avg_name in self.transp_variables_all:
                            mask_o = self.transp_variables_all[var_avg_name]['mask']
                            mask = mask_o[limits[0]:limits[1], limits[2]:limits[3]]
                    indices_good = np.where(mask > 0)
                    indices_mask = np.where(mask <= 0)
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

    def filter_granules(self):
        info_filtered = []
        if len(self.info) == 0:
            return info_filtered
        for name in self.info:
            if self.th_nvalid >= 0:
                nvalid = self.info[name]['n_valid']
                if nvalid <= self.th_nvalid:
                    print(
                        f'[INFO] Number of valid pixels {nvalid} must be greater than {self.th_nvalid}. Granule {name} is skipped.')
                    continue
            info_filtered.append(name)

        return info_filtered

    def group_granules(self):
        info_agrup = {}
        ngranules = len(self.info)
        if ngranules == 0:
            return info_agrup
        granules_included = [False] * ngranules
        names = list(self.info.keys())
        igroup = 1
        ngranules_included = 0
        for idx in range(ngranules):
            if granules_included[idx]:
                continue
            name = names[idx]
            if self.th_nvalid >= 0:
                nvalid = self.info[name]['n_valid']
                if nvalid <= self.th_nvalid:
                    granules_included[idx] = True
                    ngranules_included = ngranules_included + 1
                    continue

            limits = [self.info[name]['y_min'], self.info[name]['y_max'], self.info[name]['x_min'],
                      self.info[name]['x_max']]
            granules_included[idx] = True
            granules_list = [name]
            extended_limits, granules_list, granules_included = self.get_extended_limits(limits, granules_list,
                                                                                         granules_included)
            if self.verbose:
                print(f'[INFO] Group of granules {igroup} with {len(granules_list)} granules')
            ngranules_included = ngranules_included + len(granules_list)
            info_agrup[str(igroup)] = {
                'granules_list': granules_list,
                'y_min': extended_limits[0],
                'y_max': extended_limits[1],
                'x_min': extended_limits[2],
                'x_max': extended_limits[3]
            }
            igroup = igroup + 1

        if ngranules == ngranules_included and self.verbose:
            print(f'[INFO] Grouping granules compleded with {igroup - 1} {len(info_agrup)} groups')

        return info_agrup

    def check_dataset_attrs(self, dataset):
        at_list = ['relative_orbit', 'resampled_ymin', 'resampled_ymax', 'resampled_xmin', 'resampled_xmax',
                   'resampled_n_total', 'resampled_n_valid', 'granule_index', 'start_date']
        at_list_dataset = dataset.ncattrs()
        for at in at_list:
            if at not in at_list_dataset:
                print(f'[ERROR] Attribute {at} is not available in dataset')
                return False
        return True

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

            b = self.check_dataset_attrs(dataset)
            if not b:
                print(f'[WARNING] File: {finput} is not a valid resampled dataset. Skipping...')
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

    # def get_masks_bands(self):
    #     mask_bands = []
    #     for band in self.average_variables:
    #         if self.verbose:
    #             print(f'[INFO] Starting mask for band: {band}')
    #         mask = np.zeros((self.height,self.width))
    #         mask_bands.append(mask)
    #     return mask_bands
    def get_extended_limits(self, limits, granules_list, granules_included):
        names = list(self.info.keys())
        extended_limits = limits
        for idx in range(len(self.info)):
            if len(granules_list) == 10:
                break
            if granules_included[idx]:
                continue
            name = names[idx]
            if self.th_nvalid >= 0:
                nvalid = self.info[name]['n_valid']
                if nvalid <= self.th_nvalid:
                    # granules_included[idx] = True
                    continue

            ymin_here = self.info[name]['y_min']
            ymax_here = self.info[name]['y_max']
            xmin_here = self.info[name]['x_min']
            xmax_here = self.info[name]['x_max']
            add_granule = self.check_add_granule(ymin_here, ymax_here, xmin_here, xmax_here, limits)
            if add_granule:
                granules_included[idx] = True
                granules_list.append(name)
                if ymin_here < extended_limits[0]:
                    extended_limits[0] = ymin_here
                if ymax_here > extended_limits[1]:
                    extended_limits[1] = ymax_here
                if xmin_here < extended_limits[2]:
                    extended_limits[2] = xmin_here
                if xmax_here > extended_limits[3]:
                    extended_limits[3] = xmax_here

        return extended_limits, granules_list, granules_included

    def get_rec_info(self, limits):
        info_rec = {}
        for name in self.info:
            if self.th_nvalid >= 0:
                nvalid = self.info[name]['n_valid']
                if nvalid <= self.th_nvalid:
                    print(
                        f'[INFO] Number of valid pixels {nvalid} must be greater than {self.th_nvalid}. Granule {name} is skipped.')
                    continue
            ymin_here = self.info[name]['y_min']
            ymax_here = self.info[name]['y_max']
            xmin_here = self.info[name]['x_min']
            xmax_here = self.info[name]['x_max']
            add_granule = self.check_add_granule(ymin_here, ymax_here, xmin_here, xmax_here, limits)

            if add_granule:
                info_rec[name] = self.info[name]
        return info_rec

    def check_add_granule(self, ymin_here, ymax_here, xmin_here, xmax_here, limits):
        add_granule = False
        if limits[0] <= ymin_here <= limits[1] and limits[2] <= xmin_here <= limits[3]:
            add_granule = True
        if limits[0] <= ymin_here <= limits[1] and limits[2] <= xmax_here <= limits[3]:
            add_granule = True
        if limits[0] <= ymax_here <= limits[1] and limits[2] <= xmin_here <= limits[3]:
            add_granule = True
        if limits[0] <= ymax_here <= limits[1] and limits[2] <= xmax_here <= limits[3]:
            add_granule = True
        if limits[0] <= ymin_here <= limits[1] and xmin_here < limits[2] and xmax_here > limits[3]:
            add_granule = True
        if limits[0] <= ymax_here <= limits[1] and xmin_here < limits[2] and xmax_here > limits[3]:
            add_granule = True
        if limits[2] <= xmin_here <= limits[3] and ymin_here < limits[0] and ymax_here > limits[1]:
            add_granule = True
        if limits[2] <= xmax_here <= limits[3] and ymin_here < limits[0] and ymax_here > limits[1]:
            add_granule = True
        if ymin_here < limits[0] and ymax_here > limits[1] and xmin_here < limits[2] and xmax_here > limits[3]:
            add_granule = True
        return add_granule

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
