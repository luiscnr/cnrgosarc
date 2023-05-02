import os.path

from arc_mapinfo import ArcMapInfo
from arc_gpr_model import ARC_GPR_MODEL
from netCDF4 import Dataset
from datetime import datetime as dt
import numpy as np
import warnings

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', RuntimeWarning)


class ArcProcessing:

    def __init__(self, arc_opt, verbose, output_type, file_attributes):
        self.verbose = verbose
        self.arc_opt = arc_opt
        self.ami = ArcMapInfo(self.arc_opt, False)
        self.width = self.ami.area_def.width
        self.height = self.ami.area_def.height
        self.chla_model = None
        if output_type is None:
            output_type = 'CHLA'
        self.output_type = output_type
        self.file_attributes = file_attributes

        # FILE DEFAULTS
        file_model_default = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/SeSARC/ArcModel.json'
        file_grid_default = self.ami.ifile_base

        if arc_opt is None:  # only for creating the file base
            self.file_grid = file_grid_default
            self.file_model = file_model_default
        else:
            ##GETTING GENERAL PARAMETERS
            section = 'PROCESSING'
            self.file_model = self.arc_opt.get_value_param(section, 'file_model', file_model_default, 'str')
            self.file_grid = arc_opt.get_value_param(section, 'file_base', file_grid_default, 'str')
            self.ystep = arc_opt.get_value_param(section, 'ystep', 6500, 'int')
            self.xstep = arc_opt.get_value_param(section, 'xstep', 6500, 'int')

        if os.path.exists(self.file_model):
            self.chla_model = ARC_GPR_MODEL(self.file_model)

    def create_source_list(self,list_files,file_out):
        if os.path.exists(file_out):
            os.remove(file_out)
        f1 = open(file_out,'w')
        for file in list_files:
            f1.write(file)
            f1.write('\n')
        f1.close()


    def compute_month(self, year, month, timeliness, list_files, file_out):
        section = 'PROCESSING'
        file_base = self.arc_opt.get_value_param(section, 'file_base', None, 'str')
        if file_base is None:
            file_base = self.file_grid
        if not os.path.exists(file_base):
            print(f'[ERROR] File base {file_base} could not be found.')
            return

        ##file_grid is used for getting SENSORMARK
        file_grid = self.arc_opt.get_value_param(section, 'file_grid', None, 'str')
        if file_grid is None:
            file_grid = self.file_grid
        if not os.path.exists(file_grid):
            print(f'[ERROR] File grid {file_grid} could not be found.')
            return

        # note that output_type (TRANSP or CHLA) is defined when the class is started
        datasetout = self.create_nc_file_out_month(file_out, file_base, timeliness)
        if datasetout is None:
            print('[ERROR] Output dataset could not be started. Exiting.')
            return
        datasetgrid = Dataset(file_grid)

        from calendar import monthrange
        last_day = monthrange(year, month)[1]
        datasetout.start_date = dt(year, month, 1).strftime('%Y-%m-%d')
        datasetout.stop_date = dt(year, month, last_day).strftime('%Y-%m-%d')
        cdate = dt.utcnow()
        datasetout.creation_date = cdate.strftime('%Y-%m-%d')
        datasetout.creation_time = cdate.strftime('%H:%M:%S UTC')

        timeseconds = (dt(year, month, 1) - dt(1981, 1, 1, 0, 0, 0)).total_seconds()
        datasetout.variables['time'][0] = [np.int32(timeseconds)]

        if self.output_type == 'CHLA':
            var_avg_name = 'CHL'
            var_avg_count_name = 'CHL_count'
            var_avg_error_name = 'CHL_error'
        if self.output_type == 'TRANSP':
            var_avg_name = 'KD490'
            var_avg_count_name = 'KD490_count'
            var_avg_error_name = 'KD490_error'

        var_avg = datasetout.variables[var_avg_name]
        var_avg_count = datasetout.variables[var_avg_count_name]
        var_avg_error = datasetout.variables[var_avg_error_name]
        var_smask = datasetgrid.variables['SENSORMASK']

        self.ystep = self.arc_opt.get_value_param(section, 'ystep', 6500, 'int')
        self.xstep = self.arc_opt.get_value_param(section, 'xstep', 6500, 'int')

        iprogress = 1
        iprogress_end = np.ceil((self.height / self.ystep) * (self.width / self.xstep))

        for y in range(0, self.height, self.ystep):
            for x in range(0, self.width, self.xstep):
                if self.verbose:
                    print(f'[INFO] AVERAGING {iprogress}/{iprogress_end}')
                iprogress = iprogress + 1
                limits = self.get_limits(y, x, self.ystep, self.xstep, self.height, self.width)

                array = np.array(var_avg[0, limits[0]:limits[1], limits[2]:limits[3]])
                count_array = np.array(var_avg_count[0, limits[0]:limits[1], limits[2]:limits[3]])
                error_array = np.array(var_avg_error[0, limits[0]:limits[1], limits[2]:limits[3]])
                # print(array.shape)
                mask_array = np.array(var_smask[0, limits[0]:limits[1], limits[2]:limits[3]])
                # print(var_smask.shape)
                # print(mask_array.shape)

                x2array = np.zeros(array.shape)
                xarray = np.zeros(array.shape)

                for file_here in list_files:
                    nc_sat = Dataset(file_here)
                    array_here = np.array(nc_sat.variables[var_avg_name][0, limits[0]:limits[1], limits[2]:limits[3]])
                    nc_sat.close()
                    indices = np.where(array_here > 0)
                    if len(indices[0]) == 0:
                        continue

                    array = np.where(np.logical_and(array_here > 0, array == -999.0), 0, array)
                    array[indices] = array[indices] + array_here[indices]
                    count_array[indices] = count_array[indices] + 1

                    xarray[indices] = xarray[indices] + array_here[indices]
                    x2array[indices] = x2array[indices] + np.power(array_here[indices], 2)

                indices_good = np.where(array > 0)
                array[indices_good] = array[indices_good] / count_array[indices_good]

                xarray[indices_good] = np.power(xarray[indices_good], 2) / count_array[indices_good]
                coef_array = np.zeros(count_array.shape)
                coef_array[indices_good] = 1 / (count_array[indices_good] - 1)
                error_array[indices_good] = coef_array[indices_good] * (x2array[indices_good] - xarray[indices_good])

                array[mask_array == -999.0] = -999.0
                count_array[mask_array == -999.0] = -999.0
                error_array[mask_array == -999.0] = -999.0

                var_avg[0, limits[0]:limits[1], limits[2]:limits[3]] = [array[:, :]]
                var_avg_count[0, limits[0]:limits[1], limits[2]:limits[3]] = [count_array[:, :]]
                var_avg_error[0, limits[0]:limits[1], limits[2]:limits[3]] = [error_array[:, :]]

        datasetout.close()
        datasetgrid.close()

    def compute_chla_month(self, fileout, timeliness):
        section = 'PROCESSING'
        if self.chla_model is None:
            print('[ERROR] Chla model could not be initiated. Please review file_model option in PROCESSING section.')
            return
        file_base = self.arc_opt.get_value_param(section, 'file_base', None, 'str')
        if file_base is None:
            file_base = self.file_grid
        if not os.path.exists(file_base):
            print(f'[ERROR] File base {file_base} could not be found.')
            return

        datasetout = self.create_nc_file_out(fileout, file_base, timeliness)
        if datasetout is None:
            print('[ERROR] Output dataset could not be started. Exiting.')
            return

        datasetout.start_date = "2019-07-01"
        datasetout.stop_date = "2019-07-31"
        datasetout.timeliness = timeliness
        cdate = dt.utcnow()
        datasetout.creation_date = cdate.strftime('%Y-%m-%d')
        datasetout.creation_time = cdate.strftime('%H:%M:%S UTC')

        datasetout.cmems_product_id = "OCEANCOLOUR_ARC_BGC_L4_MY_009_124"
        datasetout.title = "cmems_obs-oc_arc_bgc-plankton_my_l4-olci-300m_P1M"

        var = datasetout.variables['CHL']
        var.cell_methods = "time: mean (interval: 1 month  comment: sampled instantaneously)"

        var = datasetout.createVariable('CHL_count', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = 0
        var.grid_mapping = 'stereographic'
        var.coordinates = 'longitude latitude'
        var.long_name = "OLCI Number Of Observations Of Monthly Chlorophyll a concentration"
        var.type = "surface"
        var.units = "1"
        var.missing_value = -999.0
        var.valid_min = 0.0
        var.valid_max = 31.0

        var = datasetout.createVariable('CHL_error', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = 0
        var.grid_mapping = 'stereographic'
        var.coordinates = 'longitude latitude'
        var.long_name = "OLCI Standard Deviation Of Monthly Chlorophyll a concentration"
        var.type = "surface"
        var.units = "milligram m^-3"
        var.missing_value = -999.0
        var.valid_min = 0.003
        var.valid_max = 100.0

        file1 = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/INTEGRATED/2019/207/O2019207_plankton-arc-fr.nc'
        file2 = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/INTEGRATED/2019/175/O2019175_plankton-arc-fr.nc'
        files = [file2, file1]
        ifile = [0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1]
        iweight = [0.1, 1, 0.5, 1.2, 1, 1, 0, 0, 0.25, 0.60, 0.75, 0.60, 0.50, 0.4, 0.8, 1, 1, 1.2, 1.3, 0.9, 0.1, 0.4,
                   0.6, 0.3, 0.2, 0.8, 1, 0, 0.9, 1.2, 1.0]
        var_chl = datasetout.variables['CHL']
        var_chl_count = datasetout.variables['CHL_count']
        var_chl_error = datasetout.variables['CHL_error']
        var_smask = datasetout.variables['SENSORMASK']
        self.xstep = 6500
        self.ystep = 6500
        iprogress = 1
        iprogress_end = np.ceil((self.height / self.ystep) * (self.width / self.xstep))

        for y in range(0, self.height, self.ystep):
            for x in range(0, self.width, self.xstep):
                print(f'AVERAGING {iprogress}/{iprogress_end}')
                iprogress = iprogress + 1
                limits = self.get_limits(y, x, self.ystep, self.xstep, self.height, self.width)

                chl_array = np.array(var_chl[limits[0]:limits[1], limits[2]:limits[3]])
                chl_count_array = np.array(var_chl_count[limits[0]:limits[1], limits[2]:limits[3]])
                chl_error_array = np.array(var_chl_error[limits[0]:limits[1], limits[2]:limits[3]])
                mask_array = np.array(var_smask[limits[0]:limits[1], limits[2]:limits[3]])

                x2array = np.zeros(chl_array.shape)
                xarray = np.zeros(chl_array.shape)

                for idx in range(len(ifile)):
                    file_here = files[ifile[idx]]
                    nc_sat = Dataset(file_here)
                    chl_here = np.array((nc_sat.variables['CHL'][0, limits[0]:limits[1], limits[2]:limits[3]]))
                    nc_sat.close()
                    indices = np.where(chl_here > 0)
                    chl_here[indices] = np.multiply(chl_here[indices], np.float(iweight[idx]))
                    indices = np.where(chl_here > 0)
                    if len(indices[0]) == 0:
                        continue

                    chl_array = np.where(np.logical_and(chl_here > 0, chl_array == -999.0), 0, chl_array)
                    # print('ATENCION: Max de chla array deberia ser al menos zero: ', np.max(chl_array))

                    # chl_errr_array = np.where(np.logical_and(chl_here > 0, chl_errr_array == -999), 0, chl_errr_array)

                    chl_array[indices] = chl_array[indices] + chl_here[indices]
                    chl_count_array[indices] = chl_count_array[indices] + 1

                    xarray[indices] = xarray[indices] + chl_here[indices]
                    x2array[indices] = x2array[indices] + np.power(chl_here[indices], 2)

                indices_good = np.where(chl_array > 0)
                chl_array[indices_good] = chl_array[indices_good] / chl_count_array[indices_good]

                xarray[indices_good] = np.power(xarray[indices_good], 2) / chl_count_array[indices_good]
                coef_array = np.zeros(chl_count_array.shape)
                coef_array[indices_good] = 1 / (chl_count_array[indices_good] - 1)
                chl_error_array[indices_good] = coef_array[indices_good] * (
                        x2array[indices_good] - xarray[indices_good])

                chl_array[mask_array == -999.0] = -999.0
                chl_count_array[mask_array == -999.0] = -999.0
                chl_error_array[mask_array == -999.0] = -999.0

                var_chl[limits[0]:limits[1], limits[2]:limits[3]] = [chl_array[:, :]]
                var_chl_count[limits[0]:limits[1], limits[2]:limits[3]] = [chl_count_array[:, :]]
                var_chl_error[limits[0]:limits[1], limits[2]:limits[3]] = [chl_error_array[:, :]]

        datasetout.close()

    def compute_chla_image(self, filein, fileout, timeliness):
        section = 'PROCESSING'
        if self.chla_model is None:
            print('[ERROR] Chla model could not be initiated. Please review file_model option in PROCESSING section.')
            return
        file_base = self.arc_opt.get_value_param(section, 'file_base', None, 'str')
        if file_base is None:
            file_base = self.file_grid
        if not os.path.exists(file_base):
            print(f'[ERROR] File base {file_base} could not be found.')
            return

        ncsat = Dataset(filein)
        rrs_bands = ['RRS442_5', 'RRS490', 'RRS560', 'RRS665']
        if self.verbose:
            print('[INFO] Checking RRS bands')
        for band in rrs_bands:
            if band not in ncsat.variables:
                print(f'[ERROR] RRS variable {band} is not available. Exiting...')
                return
        var443 = ncsat.variables[rrs_bands[0]]
        var490 = ncsat.variables[rrs_bands[1]]
        var560 = ncsat.variables[rrs_bands[2]]
        var665 = ncsat.variables[rrs_bands[3]]

        if self.verbose:
            print('[INFO] Checking longitude... ')
        if 'lon' in ncsat.variables:
            varLong = ncsat.variables['lon']
            ncgrid = None
        else:
            ncgrid = Dataset(file_base)
            if 'lon' in ncgrid.variables:
                varLong = ncgrid.variables['lon']
            else:
                print(f'[ERROR] Longitude variable is not available in {file_base}. Exiting.')
                return
        if self.verbose:
            print('[INFO] Checking date...')
        sat_date = None
        if 'start_date' in ncsat.ncattrs():
            try:
                sat_date = dt.strptime(ncsat.start_date, '%Y-%m-%d')
            except:
                sat_date = None
        else:
            fname = filein.split('/')[-1]
            if fname.startswith('O'):
                try:
                    sat_date = dt.strptime(fname[1:8], '%Y%j')
                except:
                    sat_date = None

        if sat_date is None:
            print(
                f'[ERROR] Satellite date could not be set from file name and start_date attribute is not available/valid')
            return

        jday = int(sat_date.strftime('%j'))

        if self.verbose:
            print(f'[INFO] Creating ouptput file: {fileout}')
        datasetout = self.create_nc_file_out(fileout, file_base, timeliness)
        if datasetout is None:
            print('[ERROR] Output dataset could not be started. Exiting.')
            return

        if self.verbose:
            print(f'[INFO] Setting times...')
        if sat_date is not None:
            datasetout.start_date = sat_date.strftime('%Y-%m-%d')
            datasetout.stop_date = sat_date.strftime('%Y-%m-%d')
        if timeliness is not None:
            datasetout.timeliness = timeliness
        cdate = dt.utcnow()
        datasetout.creation_date = cdate.strftime('%Y-%m-%d')
        datasetout.creation_time = cdate.strftime('%H:%M:%S UTC')

        timeseconds = (sat_date - dt(1981, 1, 1, 0, 0, 0)).total_seconds()
        datasetout.variables['time'][0] = [np.int32(timeseconds)]


        if self.verbose:
            print(f'[INFO] Getting sensor mask...')
        sensor_mask_array = np.array(ncsat.variables['SENSORMASK'])
        datasetout.variables['SENSORMASK'] = [sensor_mask_array]

        if self.verbose:
            print(f'[INFO] Getting chl variable...')
        var_chla = datasetout.variables['CHL']
        min_value = var_chla.valid_min
        max_value = var_chla.valid_max

        if self.output_type == 'COMPARISON':
            var_chla_prev = datasetout.variables['chla']
            var_diff = datasetout.variables['DIFF']

        if self.verbose:
            print(f'[INFO] Computing chla. YStep: {self.ystep} XStep: {self.xstep}')

        iprogress = 1
        iprogress_end = np.ceil((self.height / self.ystep) * (self.width / self.xstep))
        for y in range(0, self.height, self.ystep):
            for x in range(0, self.width, self.xstep):
                iprogress = iprogress + 1
                limits = self.get_limits(y, x, self.ystep, self.xstep, self.height, self.width)
                array_443 = np.array(var443[0, limits[0]:limits[1], limits[2]:limits[3]])
                array_490 = np.array(var490[0, limits[0]:limits[1], limits[2]:limits[3]])
                array_560 = np.array(var560[0, limits[0]:limits[1], limits[2]:limits[3]])
                array_665 = np.array(var665[0, limits[0]:limits[1], limits[2]:limits[3]])
                nvalid = self.chla_model.check_chla_valid(array_443, array_490, array_560, array_665)
                if self.verbose:
                    print(f'[INFO] -> {self.ystep} {self.xstep} ({iprogress} / {iprogress_end}) -> {nvalid}')
                if nvalid > 0:
                    array_long = np.array(varLong[limits[0]:limits[1], limits[2]:limits[3]])

                    array_chla = self.chla_model.compute_chla_from_2d_arrays(array_long, jday, array_443, array_490,
                                                                             array_560, array_665)

                    array_chla[array_chla < min_value] = -999.0
                    array_chla[array_chla > max_value] = -999.0
                    var_chla[0, limits[0]:limits[1], limits[2]:limits[3]] = [array_chla[:, :]]
                    if self.output_type == 'COMPARISON':
                        array_chla_prev = np.array(var_chla_prev[limits[0]:limits[1], limits[2]:limits[3]])
                        array_diff = array_chla_prev / array_chla
                        var_diff[limits[0]:limits[1], limits[2]:limits[3]] = [array_diff]
        ncsat.close()
        if ncgrid is not None:
            ncgrid.close()
        datasetout.close()
        if self.verbose:
            print('[INFO] Chla computation completed. ')

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
        if self.output_type == 'CHLA':
            at_dict['parameter'] = 'Chlorophyll-a concentration'
            at_dict['parameter_code'] = 'PLANKTON'
            if timeliness == 'NR':
                at_dict['timeliness'] = 'NR'
                at_dict['cmems_product_id'] = 'OCEANCOLOUR_ARC_BGC_L3_NRT_009_121'
                at_dict['title'] = 'cmems_obs-oc_arc_bgc-plankton_nrt_l3-olci-300m_P1D'
            if timeliness == 'NT':
                at_dict['timeliness'] = 'NT'
                at_dict['cmems_product_id'] = 'OCEANCOLOUR_ARC_BGC_L3_MY_009_123'
                at_dict['title'] = 'cmems_obs-oc_arc_bgc-plankton_my_l3-olci-300m_P1D'

        if self.output_type == 'TRANSP':  # ONLY NRT MONTHLY
            at_dict['parameter'] = 'Diffuse attenuation coefficient at 490nm'
            at_dict['parameter_code'] = 'KD490'
            at_dict['timeliness'] = timeliness

        return at_dict

    def create_nc_file_out(self, ofname, file_base, timeliness):
        if self.verbose:
            print(f'[INFO] Copying file base {file_base} to start output file {ofname}...')
        self.ami.ifile_base = file_base
        self.ami.verbose = self.verbose
        datasetout = self.ami.copy_nc_base(ofname)
        if datasetout is None:
            return datasetout

        ##global attributes
        if self.verbose:
            print(f'[INFO] Setting global attributes...')
        atribs = self.get_global_attributes(timeliness)
        if atribs is not None:  ##atrib could be defined in file base
            for at in atribs:
                if at == 'conventions':
                    datasetout.setncattr('Conventions', atribs[at])
                else:
                    datasetout.setncattr(at, atribs[at])

        ##CREATE CHL VARIABLE
        if self.output_type == 'CHLA' or self.output_type == 'COMPARISON':
            if 'CHL' not in datasetout.variables:
                if self.verbose:
                    print('[INFO] Creating CHL variable...')
                var = datasetout.createVariable('CHL', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                                complevel=6)
                var[:] = -999
                var.grid_mapping = 'stereographic'
                var.coordinates = 'lon lat'
                var.long_name = "Chlorophyll a concentration"
                var.standard_name = "mass_concentration_of_chlorophyll_a_in_sea_water"
                var.type = "surface"
                var.units = "milligram m^-3"
                var.missing_value = -999.0
                var.valid_min = 0.003
                var.valid_max = 100.0
                var.comment = "OLCI - WFR STANDARD PROCESSOR - Gaussian Processor Regressor (GPR) Algorithm"

        if self.output_type == 'COMPARISON':
            if self.verbose:
                print('[INFO] Creating DIFF variable')
            if not 'DIFF' in datasetout.variables:
                var = datasetout.createVariable('DIFF', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
                var[:] = -999

        return datasetout

    def create_nc_file_out_month(self, ofname, file_base, timeliness):
        if self.verbose:
            print(f'[INFO] Copying file base {file_base} to start output file {ofname}...')
        self.ami.ifile_base = file_base
        self.ami.verbose = self.verbose
        datasetout = self.ami.copy_nc_base(ofname)
        if datasetout is None:
            return datasetout

        ##global attributes
        if self.verbose:
            print(f'[INFO] Setting global attributes...')
        atribs = self.get_global_attributes(timeliness)
        if atribs is not None:  ##atrib could be defined in file base
            for at in atribs:
                if at == 'conventions':
                    datasetout.setncattr('Conventions', atribs[at])
                else:
                    datasetout.setncattr(at, atribs[at])

        datasetout.product_level = 'L4'
        if self.output_type == 'CHLA':
            if timeliness == 'NR':
                datasetout.cmems_product_id = 'OCEANCOLOUR_ARC_BGC_L4_NRT_009_122'
                datasetout.title = 'cmems_obs-oc_arc_bgc-plankton_nrt_l4-olci-300m_P1M'
            if timeliness == 'NT':
                datasetout.cmems_product_id = 'OCEANCOLOUR_ARC_BGC_L4_MY_009_124'
                datasetout.title = 'cmems_obs-oc_arc_bgc-plankton_my_l4-olci-300m_P1M'

        if self.output_type == 'TRANSP':
            if timeliness == 'NR':
                datasetout.cmems_product_id = 'OCEANCOLOUR_ARC_BGC_L4_NRT_009_122'
                datasetout.title = 'cmems_obs-oc_arc_bgc-transp_nrt_l4-olci-300m_P1M'

        # CHLA
        if self.output_type == 'CHLA':
            if 'CHL' not in datasetout.variables:
                if self.verbose:
                    print('[INFO] Creating CHL variable...')
                var = datasetout.createVariable('CHL', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                                complevel=6)
                var[:] = -999
                var.grid_mapping = 'stereographic'
                var.coordinates = 'lon lat'
                var.long_name = "Chlorophyll a concentration"
                var.standard_name = "mass_concentration_of_chlorophyll_a_in_sea_water"
                var.type = "surface"
                var.units = "milligram m^-3"
                var.missing_value = -999.0
                var.valid_min = 0.003
                var.valid_max = 100.0
                var.comment = "OLCI - WFR STANDARD PROCESSOR - Gaussian Processor Regressor (GPR) Algorithm"
                var.cell_methods = "time: mean (interval: 1 month  comment: sampled instantaneously)"

            if 'CHL_count' not in datasetout.variables:
                var = datasetout.createVariable('CHL_count', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                                complevel=6)
                var[:] = 0
                var.grid_mapping = 'stereographic'
                var.coordinates = 'lon lat'
                var.long_name = "OLCI Number Of Observations Of Monthly Chlorophyll a concentration"
                var.type = "surface"
                var.units = "1"
                var.missing_value = -999.0
                var.valid_min = 0.0
                var.valid_max = 31.0

            if 'CHL_error' not in datasetout.variables:
                var = datasetout.createVariable('CHL_error', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                                complevel=6)
                var[:] = 0
                var.grid_mapping = 'stereographic'
                var.coordinates = 'lon lat'
                var.long_name = "OLCI Standard Deviation Of Monthly Chlorophyll a concentration"
                var.type = "surface"
                var.units = "milligram m^-3"
                var.missing_value = -999.0
                var.valid_min = 0.003
                var.valid_max = 100.0

        # TRANSP
        if self.output_type == 'TRANSP':
            if 'KD490' not in datasetout.variables:
                if self.verbose:
                    print('[INFO] Creating KD490 variable...')
                var = datasetout.createVariable('KD490', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                                complevel=6)
                var[:] = -999
                var.grid_mapping = 'stereographic'
                var.coordinates = 'lon lat'
                var.long_name = "OLCI Diffuse Attenuation Coefficient at 490nm"
                var.standard_name = "volume_attenuation_coefficient_of_downwelling_radiative_flux_in_sea_water"
                var.type = "surface"
                var.units = "m^-1"
                var.missing_value = -999.0
                var.valid_min = 0.0
                var.valid_max = 10.0
                var.comment = "OLCI - WFR STANDARD PROCESSOR"
                var.cell_methods = "time: mean (interval: 1 month  comment: sampled instantaneously)"

            if 'KD490_count' not in datasetout.variables:
                var = datasetout.createVariable('KD490_count', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                                complevel=6)
                var[:] = 0
                var.grid_mapping = 'stereographic'
                var.coordinates = 'lon lat'
                var.long_name = "OLCI Number Of Observations Of Monthly Diffuse Attenuation Coefficient at 490nm"
                var.type = "surface"
                var.units = "1"
                var.missing_value = -999.0
                var.valid_min = 0.0
                var.valid_max = 31.0

            if 'KD490_error' not in datasetout.variables:
                var = datasetout.createVariable('KD490_error', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                                complevel=6)
                var[:] = 0
                var.grid_mapping = 'stereographic'
                var.coordinates = 'lon lat'
                var.long_name = "OLCI Standard Deviation Of Monthly Diffuse Attenuation Coefficient at 490nm"
                var.type = "surface"
                var.units = "m^-1"
                var.missing_value = -999.0
                var.valid_min = 0.0
                var.valid_max = 10.0

        return datasetout
