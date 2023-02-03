import os.path

from arc_mapinfo import ArcMapInfo
from arc_gpr_model import ARC_GPR_MODEL
from netCDF4 import Dataset
from datetime import datetime as dt
import numpy as np


class ArcProcessing:

    def __init__(self, arc_opt, verbose):
        self.verbose = verbose
        self.arc_opt = arc_opt
        self.ami = ArcMapInfo(self.arc_opt, False)
        self.width = self.ami.area_def.width
        self.height = self.ami.area_def.height
        self.chla_model = None

        # FILE DEFAULTS
        file_model_default = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/SeSARC/ArcModel.json'
        file_grid_default = self.ami.ifile_base

        ##GETTING GENERAL PARAMETERS
        section = 'PROCESSING'
        self.file_model = self.arc_opt.get_value_param(section, 'file_model', file_model_default, 'str')
        if os.path.exists(self.file_model):
            self.chla_model = ARC_GPR_MODEL(self.file_model)
        self.file_grid = arc_opt.get_value_param(section, 'file_base', file_grid_default, 'str')
        self.ystep = arc_opt.get_value_param(section, 'ystep', 6500, 'int')
        self.xstep = arc_opt.get_value_param(section, 'xstep', 6500, 'int')

    def compute_chla_image(self, filein, fileout):
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
            print('[INFO] Checking RRS bands... ')
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
        if 'longitude' in ncsat.variables:
            varLong = ncsat.variables['longitude']
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

        jday = int(sat_date.strftime('%j'))

        if self.verbose:
            print(f'[INFO] Creating ouptput file: {fileout}')
        datasetout = self.create_nc_file_out_chl(fileout, file_base)
        if datasetout is None:
            print('[ERROR] Output dataset could not be started. Exiting.')
            return
        var_chla = datasetout.variables['CHL']

        var_chla_prev = datasetout.variables['chla']
        var_diff = datasetout.variables['DIFF']

        if self.verbose:
            print(f'[INFO] Computing chla. YStep: {self.ystep} XStep: {self.xstep}')
        valid_tiles = []

        for y in range(0, self.height, self.ystep):
            if self.verbose:
                print(f'[INFO] -> {y}')
            for x in range(0, self.width, self.xstep):
                if self.verbose:
                    print(f'[INFO] -> {y} {x}')
                limits = self.get_limits(y, x, self.ystep, self.xstep, self.height, self.width)
                #array_long = np.array(varLong[limits[0]:limits[1], limits[2]:limits[3]])
                array_443 = np.array(var443[limits[0]:limits[1], limits[2]:limits[3]])
                array_490 = np.array(var490[limits[0]:limits[1], limits[2]:limits[3]])
                array_560 = np.array(var560[limits[0]:limits[1], limits[2]:limits[3]])
                array_665 = np.array(var665[limits[0]:limits[1], limits[2]:limits[3]])
                nvalid = self.chla_model.check_chla_valid(array_443,array_490,array_560,array_665)
                if nvalid>0:
                    valid_tiles.append(limits)

                # array_chla_prev = np.array(var_chla_prev[limits[0]:limits[1], limits[2]:limits[3]])
                #
                # array_chla = self.chla_model.compute_chla_from_2d_arrays(array_long, jday, array_443, array_490,
                #                                                          array_560, array_665)
                # array_diff = array_chla_prev / array_chla
                #
                # var_chla[limits[0]:limits[1], limits[2]:limits[3]] = [array_chla[:, :]]
                # var_diff[limits[0]:limits[1], limits[2]:limits[3]] = [array_diff]

        nvalidtiles = len(valid_tiles)
        itile = 0
        print('NValidTiles: ',nvalidtiles)
        for limits in valid_tiles:

            if self.verbose:
                itile = itile + 1
                print(f'[INFO] Computing chla for tile: {itile} / {nvalidtiles} ')
            array_long = np.array(varLong[limits[0]:limits[1], limits[2]:limits[3]])
            array_443 = np.array(var443[limits[0]:limits[1], limits[2]:limits[3]])
            array_490 = np.array(var490[limits[0]:limits[1], limits[2]:limits[3]])
            array_560 = np.array(var560[limits[0]:limits[1], limits[2]:limits[3]])
            array_665 = np.array(var665[limits[0]:limits[1], limits[2]:limits[3]])
            array_chla_prev = np.array(var_chla_prev[limits[0]:limits[1], limits[2]:limits[3]])

            array_chla = self.chla_model.compute_chla_from_2d_arrays(array_long, jday, array_443, array_490,
                                                                     array_560, array_665)
            array_diff = array_chla_prev / array_chla

            var_chla[limits[0]:limits[1], limits[2]:limits[3]] = [array_chla[:, :]]
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

    def create_nc_file_out_chl(self, ofname, file_base):
        if self.verbose:
            print(f'[INFO] Copying file base {file_base} to start output file {ofname}...')
        self.ami.ifile_base = file_base
        datasetout = self.ami.copy_nc_base(ofname)
        if datasetout is None:
            return datasetout

        ##CREATE CHL VARIABLE
        if self.verbose:
            print('[INFO] Creating CHL variable...')
        var = datasetout.createVariable('CHL', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = -999
        var.coordinates = "lat lon"
        var.long_name = "Chlorophyll a concentration"
        var.standard_name = "mass_concentration_of_chlorophyll_a_in_sea_water"
        var.type = "surface"
        var.units = "milligram m^-3"
        var.missing_value = -999.0
        var.valid_min = np.float(0.01)
        var.valid_max = np.float(300.0)
        var.comment = ""
        var.source = "OLCI - WFR STANDARD PROCESSOR - GPR CHL-A ALGORITHM"

        if self.verbose:
            print('[INFO] Creating DIFF variable')
        var = datasetout.createVariable('DIFF', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
        var[:] = -999

        return datasetout
