import argparse
import os.path
from datetime import datetime as dt
from datetime import timedelta
from netCDF4 import Dataset

import numpy as np
import warnings

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', RuntimeWarning)

parser = argparse.ArgumentParser(description="Arctic climatology")
parser.add_argument("-m", "--mode", help="Mode", choices=["MULTI"], required=True)
parser.add_argument("-c", "--config_file", help="Config file", required=True)
# parser.add_argument('-i', "--inputpath", help="Input directory")
parser.add_argument('-sd', "--start_date", help="Start date (yyyy-mm-dd)", required=True)
parser.add_argument('-ed', "--end_date", help="End Date (yyyy-mm-dd)", required=True)
# parser.add_argument('-o', "--output", help="Output directory")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
args = parser.parse_args()


def do_test():
    file_nc = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_MULTI_CLIMATOLOGY/TEMP_07_19/TempResults_800_1000_800_1000.nc'
    from netCDF4 import Dataset
    dataset = Dataset(file_nc)
    dc = np.array(dataset.variables['DATA_COMPUTATION'])
    dw = np.array(dataset.variables['DATA_WEIGHT'])

    avg_var = np.array(dataset.variables['W_AVG'])
    avg_var_t = np.zeros(dc.shape)
    avg_var_t[:, :, :] = avg_var[:, :]
    ndata = np.array(dataset.variables['NOFOBS'])
    num_wstd = np.power((dc - avg_var_t), 2)
    num_wstd = np.multiply(dw, num_wstd)
    el = num_wstd[:, 14, 173]
    for e in el:
        print('---> el cuadrado', e)
    lasuma = np.nansum(num_wstd, axis=0)
    print('la suma: ', lasuma[14, 173])
    # num_wstd = np.multiply(num_wstd, dw)
    factor = (ndata - 1) / ndata
    print('Ndata ', ndata[14, 173], 'factor: ', factor[14, 173])
    weighted_std = np.nansum(num_wstd, axis=0) / (factor * np.nansum(dw, axis=0))
    weighted_std = np.sqrt(weighted_std)

    data = dc[:, 14, 173]
    weights = dw[:, 14, 173]
    sum_data = 0
    sum_weights = 0
    for d, w in zip(data, weights):
        # print(d,w)
        if not np.isnan(d):
            sum_data = sum_data + (d * w)
            sum_weights = sum_weights + w
    avg_w = sum_data / sum_weights

    n_obs = 0
    num = 0
    denom = 0
    for d, w in zip(data, weights):
        # print(d,w)
        if not np.isnan(d):
            num = num + (w * ((d - avg_w) * (d - avg_w)))
            denom = denom + w
            n_obs = n_obs + 1
    print('la suma tal: ', num)
    factor = (n_obs - 1) / n_obs
    print('nobs: ', n_obs, 'factor: ', factor)
    denom = ((n_obs - 1) / n_obs) * denom
    w_std = np.sqrt(num / denom)
    # avg = np.nanmean(dc,axis=0)
    # print('-----------------------')
    # print(np.nanmean(tal))
    # print(avg[14,173])
    print('------------')
    # avg_var = np.array(dataset.variables['W_AVG'])
    # print(avg_var[14,173])
    # print(avg_w)
    # std_nor = np.array(dataset.variables['N_STD'])
    std_var = np.array(dataset.variables['W_STD'])
    # print(std_nor[14,173])
    # print('NDATA',ndata[14,173],n_obs)
    print(weighted_std[14, 173])
    print(w_std)
    print(std_var[14, 173])

    dataset.close()

    return True


def do_test2():
    file_nc = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_MULTI_CLIMATOLOGY/TEMP_07_19/Temp_200_400_400_600_2022.nc'
    from netCDF4 import Dataset
    dataset = Dataset(file_nc)
    w = np.array(dataset.variables['WEIGHTS'])
    w[w == -999] = np.nan
    sw = np.nansum(w, axis=0)
    print(np.min(sw), np.max(sw))
    dataset.close()

    return True


def main():
    # if do_test2():
    #     return
    print('[INFO] STARTED CLIMATOLOGY')
    if not os.path.exists(args.config_file):
        print(f'[ERROR] Config file {args.config_file} does not exist. Exiting...')
        return
    try:
        import configparser
        options = configparser.ConfigParser()
        options.read(args.config_file)
    except:
        print(f'[ERROR] Config file {args.config_file} could not be read. Exiting...')

    ref = 'CLIMATOLOGY'
    if not options.has_section(ref):
        print(f'[ERROR] Config file {args.config_file} must containt a section named CLIMATOLOGY')
        return

    if not options.has_option(ref, 'input_path'):
        print(f'[ERROR] input_path option is required in config file {args.config_file}')
        return

    input_path = options[ref]['input_path']
    if not os.path.exists(input_path):
        print(f'[ERROR] Input path {input_path} does not exist')
        return
    if options.has_option(ref, 'output_path'):
        output_path = options[ref]['output_path']
    else:
        output_path = input_path
    if not os.path.exists(output_path):
        try:
            os.mkdir(output_path)
        except:
            print(f'[ERROR] Output path {output_path} does not exist and could not be created')
            return
    # defaults
    options_clim = {
        'temporal_window': 11,
        'spatial_window': 3,
        'nywindow': 200,
        'nxwindow': 200,
        'use_weights': True,
        'variable': 'CHL',
        'ny': -1,
        'nx': -1,
        'file_ref': '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/GRID_FILES/ArcGrid_65_90_4KM_GridBase.nc',
        'year_ini': 1998,
        'year_end': 2022,
        'yref': -1,
        'xref': -1,
        'name_out': '$DATE$_chl_arc_multi_clima.nc',
        'applyPool': 0,
        'input_name': 'C$DATE$_chl-arc-4km.nc',
        'log_scale': False
    }

    # setting options
    for option in options_clim:
        if options.has_option(ref, option):
            if option == 'file_ref' or option == 'variable' or option == 'name_out' or option == 'input_name':
                options_clim[option] = options[ref][option].strip()
            elif option == 'use_weights' or option == 'log_scale':
                if options[ref][option].strip().lower() == 'true' or options[ref][option].strip() == '1':
                    options_clim[option] = True
                else:
                    options_clim[option] = False
            else:
                options_clim[option] = int(options[ref][option].strip())

    ## getting date
    date_here_str = args.start_date
    try:
        date_here = dt.strptime(date_here_str, '%Y-%m-%d')
    except:
        print(f'[ERROR] Start Date for climatology is not in the correct format YYYY-mm-dd')
        return

    if args.end_date:
        date_here_end_str = args.end_date
        try:
            date_here_end = dt.strptime(date_here_end_str, '%Y-%m-%d')
        except:
            print(f'[WARNING] Start Date for climatology is not in the correct format YYYY-mm-dd.')
            return
    else:
        date_here_end = date_here

    while date_here<=date_here_end:
        if args.mode == 'MULTI':
            year_ini = options_clim['year_ini']
            year_end = options_clim['year_end']
            apply_pool = options_clim['applyPool']
            if apply_pool != 0:
                params_list = []
                from multiprocessing import Pool
                if args.verbose:
                    print('[INFO] Starting parallel processing')
            else:
                if args.verbose:
                    print('[INFO] Starting sequential processing')

            for year in range(year_ini, year_end + 1):
                date_here_y = date_here.replace(year=year)
                options_clim_here = options_clim.copy()
                options_clim_here['year_ini'] = year
                options_clim_here['year_end'] = year
                warnings.simplefilter('ignore', UserWarning)
                warnings.simplefilter('ignore', RuntimeWarning)
                if apply_pool == 0:  ##sequential
                    make_clim_multi_day_extracts(input_path, output_path, date_here_y, options_clim_here)
                else:
                    params_here = [input_path, output_path, date_here_y, options_clim_here]
                    params_list.append(params_here)
                    if len(params_list) == options_clim['applyPool']:
                        if args.verbose:
                            print(f'[INFO] Running parallel processes: {apply_pool}')
                        poolhere = Pool(apply_pool)
                        poolhere.map(make_clim_multi_day_extracts_parallel, params_list)
                        params_list = []
            if apply_pool != 0 and len(params_list) > 0:
                if args.verbose:
                    print(f'[INFO] Running parallel processes: {apply_pool} -> {len(params_list)}')
                poolhere.map(make_clim_multi_day_extracts_parallel, params_list)
            make_clime_multi_day_computation(output_path, date_here, options_clim)
            make_clim_together(output_path, date_here, options_clim)

        date_here = date_here + timedelta(hours=24)



def make_clim_multi_day_extracts_parallel(params):
    make_clim_multi_day_extracts(params[0], params[1], params[2], params[3])


def make_clim_multi_day_extracts(input_path, output_path, date_here, options):
    temporal_window = options['temporal_window']  # 11
    spatial_window = options['spatial_window']  # 3
    nywindow = options['nywindow']
    nxwindow = options['nxwindow']
    use_weights = options['use_weights']
    variable = options['variable']
    ny = options['ny']
    nx = options['nx']
    input_name = options['input_name']
    log_scale = options['log_scale']

    date_list = get_input_file_list_day(date_here, temporal_window, input_path, date_here.year, date_here.year,
                                        input_name)
    if ny == -1 and nx == -1:
        for date_str in date_list:
            input_file = date_list[date_str]['input_file']
            if input_file is not None:
                ny, nx = get_ny_nx(input_file, variable)

    ndata = temporal_window * spatial_window * spatial_window

    data_computation = np.zeros((ndata, ny, nx))
    data_weights = np.ones((ndata, ny, nx))
    data_computation[:] = -999.0
    if use_weights:
        if args.verbose:
            print('[INFO] Retrieving weights...')
        spatial_weights = get_spatial_weight_filter_3x3(None, 1, 1)
        temporal_weights = get_temporal_weight_filter(None, temporal_window, 1, 1)
        index = 0
        for itemporal in range(len(temporal_weights)):
            for ispatial in range(len(spatial_weights)):
                data_weights[index, :, :] = spatial_weights[ispatial] * temporal_weights[itemporal]
                index = index + 1

    itime = 0
    idate = 0
    for date_str in date_list:
        # date_h = dt.strptime(date_str, '%Y-%m-%d')
        input_file = date_list[date_str]['input_file']
        if args.verbose:
            print(
                f'[INFO] Index date: {idate} Index all: {itime} Date: {date_str} Input file: {input_file}')
        if input_file is None:
            itime = itime + (spatial_window * spatial_window)
            continue

        if nywindow >= ny and nxwindow >= nx:
            indices = [itime, 0, ny, 0, nx, ny, nx]
            data_computation, weighted_avg = set_data_window_3x3(data_computation, input_file, variable, indices, None,
                                                                 None)
        else:
            for y in range(0, ny, nywindow):
                for x in range(0, nx, nxwindow):
                    limits = get_limits(y, x, nywindow, nxwindow, ny, nx)
                    nywindow_here = limits[1] - limits[0]
                    nxwindow_here = limits[3] - limits[2]
                    indices = [itime, limits[0], nywindow_here, limits[2], nxwindow_here, ny, nx]
                    data_computation, weighted_avg = set_data_window_3x3(data_computation, input_file, variable,
                                                                         indices, None, None)
        idate = idate + 1
        itime = itime + (spatial_window * spatial_window)

    data_weights[data_computation == -999] = -999

    mstr = date_here.strftime('%m')
    dstr = date_here.strftime('%d')
    output_temp = os.path.join(output_path, f'TEMP_{mstr}_{dstr}')
    if not os.path.isdir(output_temp):
        os.mkdir(output_temp)
    for y in range(0, ny, nywindow):
        for x in range(0, nx, nxwindow):
            limits = get_limits(y, x, nywindow, nxwindow, ny, nx)
            nywindow_here = limits[1] - limits[0]
            nxwindow_here = limits[3] - limits[2]
            output_file = os.path.join(output_temp,
                                       f'Temp_{limits[0]}_{limits[1]}_{limits[2]}_{limits[3]}_{date_here.year}.nc')
            dout = start_data_temporal_file(output_file, nywindow_here, nxwindow_here, ndata)
            array = data_computation[:, limits[0]:limits[1], limits[2]:limits[3]]
            if log_scale:
                array[array!=-999] = np.log10(array[array!=-999])
            dout.variables['DATA'][:, :, :] = array
            dout.variables['WEIGHTS'][:, :, :] = data_weights[:, limits[0]:limits[1], limits[2]:limits[3]]
            dout.close()

    print('[INFO] Completed')


def make_clime_multi_day_computation(output_path, date_here, options):
    temporal_window = options['temporal_window']  # 11
    spatial_window = options['spatial_window']  # 3
    nywindow = options['nywindow']
    nxwindow = options['nxwindow']
    ny = options['ny']
    nx = options['nx']
    year_ini = options['year_ini']  # 1998
    year_end = options['year_end']  # 2022
    apply_pool = options['applyPool']
    if apply_pool != 0:
        params_list = []
        from multiprocessing import Pool
        if args.verbose:
            print('[INFO] Starting parallel processing')

    mstr = date_here.strftime('%m')
    dstr = date_here.strftime('%d')
    output_temp = os.path.join(output_path, f'TEMP_{mstr}_{dstr}')

    nyears = (year_end - year_ini) + 1
    nbyyear = temporal_window * spatial_window * spatial_window
    ndata = nyears * temporal_window * spatial_window * spatial_window

    for y in range(0, ny, nywindow):
        for x in range(0, nx, nxwindow):
            limits = get_limits(y, x, nywindow, nxwindow, ny, nx)
            if apply_pool == 0:
                make_clime_multi_day_computation_impl(ndata, nbyyear, nyears, limits, year_ini, year_end, output_temp)
            else:
                params_here = [ndata, nbyyear, nyears, limits, year_ini, year_end, output_temp]
                params_list.append(params_here)
                if len(params_list) == apply_pool:
                    if args.verbose:
                        print(f'[INFO] Running parallel processes: {apply_pool}')
                    poolhere = Pool(apply_pool)
                    poolhere.map(make_clime_multi_day_computation_impl_parallel, params_list)
                    params_list = []
    if apply_pool != 0 and len(params_list) > 0:
        if args.verbose:
            print(f'[INFO] Running parallel processes: {apply_pool} -> {len(params_list)}')
        poolhere = Pool(apply_pool)
        poolhere.map(make_clime_multi_day_computation_impl_parallel, params_list)


def make_clime_multi_day_computation_impl_parallel(params):
    make_clime_multi_day_computation_impl(params[0], params[1], params[2], params[3], params[4], params[5], params[6])


def make_clime_multi_day_computation_impl(ndata, nbyyear, nyears, limits, year_ini, year_end, output_temp):
    nywindow_here = limits[1] - limits[0]
    nxwindow_here = limits[3] - limits[2]
    data_computation = np.zeros((ndata, nywindow_here, nxwindow_here))
    data_weights = np.ones((ndata, nywindow_here, nxwindow_here))
    output_file = os.path.join(output_temp, f'TempResults_{limits[0]}_{limits[1]}_{limits[2]}_{limits[3]}.nc')
    idata = 0
    for year in range(year_ini, year_end + 1):
        input_file = os.path.join(output_temp, f'Temp_{limits[0]}_{limits[1]}_{limits[2]}_{limits[3]}_{year}.nc')
        dset = Dataset(input_file)
        data_computation[idata:idata + nbyyear, :, :] = np.array(dset.variables['DATA'][:, :, :])
        data_weights[idata:idata + nbyyear, :, :] = np.array(dset.variables['WEIGHTS'][:, :, :])

        dset.close()
        idata = idata + nbyyear

    data_weights[data_computation == -999] = np.nan
    data_computation[data_computation == -999] = np.nan

    print(f'[INFO] Computing for {limits[0]} - {limits[1]} {limits[2]} - {limits[3]}')

    avg = np.nanmean(data_computation, axis=0)
    median = np.nanmedian(data_computation, axis=0)
    std = np.nanstd(data_computation, axis=0)
    min_array = np.nanmin(data_computation, axis=0)
    max_array = np.nanmax(data_computation, axis=0)
    nofobs = np.count_nonzero(~np.isnan(data_weights), axis=0)
    weigthed_prod = np.multiply(data_computation, data_weights)
    sum_weights = np.nansum(data_weights, axis=0)
    weighted_avg = np.nansum(weigthed_prod, axis=0) / sum_weights
    diff_n_w = np.abs(avg - weighted_avg)

    # weighted str
    weighted_avg_t = np.zeros(data_computation.shape)
    weighted_avg_t[:, :, :] = weighted_avg[:, :]
    num_wstd = np.power((data_computation - weighted_avg_t), 2)
    num_wstd = np.multiply(num_wstd, data_weights)
    factor = (nofobs - 1) / nofobs
    weighted_std = np.nansum(num_wstd, axis=0) / (factor * sum_weights)
    weighted_std = np.sqrt(weighted_std)

    pw = (sum_weights / (24 * nyears)) * 100

    # ntimes = data_weights.shape[0]
    dout = start_results_temporal_file(output_file, nywindow_here, nxwindow_here)
    dout.variables['W_AVG'][:, :] = weighted_avg[:, :]
    dout.variables['N_AVG'][:, :] = avg[:, :]
    dout.variables['W_STD'][:, :] = weighted_std[:, :]
    dout.variables['N_STD'][:, :] = std[:, :]
    dout.variables['MEDIAN'][:, :] = median[:, :]
    dout.variables['MIN'][:, :] = min_array[:, :]
    dout.variables['MAX'][:, :] = max_array[:, :]
    dout.variables['NOFOBS'][:, :] = nofobs[:, :]
    dout.variables['SUM_W'][:, :] = sum_weights[:, :]
    dout.variables['PORC_W'][:, :] = pw[:, :]
    dout.variables['DIFF_N_W'][:, :] = diff_n_w[:, :]

    # dout.variables['DATA_WEIGHT'][:, :, :] = data_weights[:, :, :]
    # dout.variables['DATA_COMPUTATION'][:, :, :] = data_computation[:, :, :]

    dout.close()


def make_clim_together(output_path, date_here, options):
    nywindow = options['nywindow']
    nxwindow = options['nxwindow']
    ny = options['ny']
    nx = options['nx']
    year_ini = options['year_ini']  # 1998
    year_end = options['year_end']  # 2022
    name_out = options['name_out']

    mstr = date_here.strftime('%m')
    dstr = date_here.strftime('%d')
    output_temp = os.path.join(output_path, f'TEMP_{mstr}_{dstr}')
    date_here_str = f'{year_ini}{mstr}{dstr}_{year_end}{mstr}{dstr}'
    name_out = name_out.replace('$DATE$', date_here_str)
    file_out = os.path.join(output_path, name_out)
    dout = start_results_temporal_file(file_out, ny, nx)

    for y in range(0, ny, nywindow):
        for x in range(0, nx, nxwindow):
            limits = get_limits(y, x, nywindow, nxwindow, ny, nx)
            # nywindow_here = limits[1] - limits[0]
            # nxwindow_here = limits[3] - limits[2]
            file_results = os.path.join(output_temp, f'TempResults_{limits[0]}_{limits[1]}_{limits[2]}_{limits[3]}.nc')
            din = Dataset(file_results)
            for var in din.variables:
                dout.variables[var][limits[0]:limits[1], limits[2]:limits[3]] = din.variables[var][:, :]
            din.close()

    dout.close()


def get_input_file_list_day(date_here, temporal_window, input_path, year_ini, year_end, input_name):
    if args.verbose:
        print(f'[INFO] Getting file list...')
    ndaysw = int(np.floor(temporal_window / 2))
    date_list = {}
    if year_ini == -1 and year_end == -1:
        year_ini = date_here.year
        year_end = date_here.year
    for year in range(year_ini, year_end + 1):
        date_here_min = date_here.replace(year=year) - timedelta(days=ndaysw)
        date_here_max = date_here.replace(year=year) + timedelta(days=ndaysw)
        date_here_ref = date_here_min
        while date_here_ref <= date_here_max:
            date_here_str = date_here_ref.strftime('%Y-%m-%d')
            year_str = date_here_ref.strftime('%Y')
            jjj_str = date_here_ref.strftime('%j')
            date_str = f'{year_str}{jjj_str}'
            name_in = input_name.replace('$DATE$', date_str)
            # input_file = os.path.join(input_path, year_str, jjj_str, f'C{year_str}{jjj_str}_chl-arc-4km.nc')
            input_file = os.path.join(input_path, year_str, jjj_str, name_in)
            if os.path.exists(input_file):
                date_list[date_here_str] = {'input_file': input_file, 'input_array': None}
            else:
                date_list[date_here_str] = {'input_file': None, 'input_array': None}
            if args.verbose:
                print(f'[INFO]  Date: {date_here_str} File: {input_file}')
            date_here_ref = date_here_ref + timedelta(days=1)
    if args.verbose:
        print(f'[INFO] Getting file list: Completed')
    return date_list


def start_data_temporal_file(output_file, ny, nx, ntimes):
    try:
        dout = Dataset(output_file, 'w', format='NETCDF4')
    except PermissionError:
        return None

    print(f'Starting file:{output_file}')

    dout.createDimension('y', ny)  # ny
    dout.createDimension('x', nx)  # nx
    dout.createDimension('time', ntimes)
    dout.createVariable('DATA', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('WEIGHTS', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True, complevel=6)
    return dout


def start_results_temporal_file(output_file, ny, nx):
    try:
        dout = Dataset(output_file, 'w', format='NETCDF4')
    except PermissionError:
        return None

    print(f'Starting file:{output_file}')

    # dout.createDimension('time', ntimes)
    dout.createDimension('y', ny)  # ny
    dout.createDimension('x', nx)  # nx
    dout.createVariable('W_AVG', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('N_AVG', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('DIFF_N_W', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('N_STD', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('W_STD', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('SUM_W', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('PORC_W', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('MEDIAN', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('MIN', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('MAX', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)
    dout.createVariable('NOFOBS', 'f4', ('y', 'x'), fill_value=-999, zlib=True, complevel=6)

    # dout.createVariable('DATA_WEIGHT', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True, complevel=6)
    # dout.createVariable('DATA_COMPUTATION', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True, complevel=6)

    return dout


def start_climatology_temporal_file(output_file, file_ref):
    from arc_mapinfo import ArcMapInfo
    ami = ArcMapInfo(None, False)
    if file_ref is not None:
        ami.ifile_base = file_ref
    datasetout = ami.copy_nc_base(output_file)
    if 'W_AVG' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating W_AVG variable...')
        var = datasetout.createVariable('W_AVG', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0
    if 'N_AVG' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating N_AVG variable...')
        var = datasetout.createVariable('N_AVG', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0
    if 'DIFF_N_W' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating DIFF_N_W variable...')
        var = datasetout.createVariable('DIFF_N_W', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0

    if 'SUM_W' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating SUM_W variable...')
        var = datasetout.createVariable('SUM_W', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0

    if 'STD' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating STD variable...')
        var = datasetout.createVariable('STD', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0

    if 'MEDIAN' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating MEDIAN variable...')
        var = datasetout.createVariable('MEDIAN', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0

    if 'MIN' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating MIN variable...')
        var = datasetout.createVariable('MIN', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0

    if 'MAX' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating MAX variable...')
        var = datasetout.createVariable('MAX', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0

    if 'PORC_W' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating PORC_W variable...')
        var = datasetout.createVariable('PORC_W', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0

    if 'N_YEARS' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating N_YEARS variable...')
        var = datasetout.createVariable('N_YEARS', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0

    if 'N_OBS' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating N_OBS variable...')
        var = datasetout.createVariable('N_OBS', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0

    if 'N_DAYS' not in datasetout.variables:
        if args.verbose:
            print('[INFO] Creating N_DAYS variable...')
        var = datasetout.createVariable('N_DAYS', 'f4', ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                        complevel=6)
        var[:] = -999.0

    return datasetout


# first dimension is time
def get_ny_nx(input_file, variable):
    dataset = Dataset(input_file)
    variable = dataset.variables[variable]
    ny = variable.shape[1]
    nx = variable.shape[2]
    return ny, nx


def get_limits(y, x, ystep, xstep, ny, nx):
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


def test_compute_avg(input_file, variable):
    # array_avg, array_avg_w = test_compute_avg(
    #     '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_MULTI_CLIMATOLOGY/2000/195/C2000195_chl-arc-4km.nc', 'CHL')
    #
    # print('==============', avg[615, 1364])
    # print('==============', weighted_avg[615, 1364])
    # print('++++++++++++++', array_avg[615, 1364])
    # print('++++++++++++++', array_avg_w[615, 1364])
    #
    # check_avg = np.divide(avg, array_avg)
    # check_avg_w = np.divide(weighted_avg, array_avg_w)
    # print('CHECKING: ')
    # print(np.nanmean(check_avg))
    # print(np.nanmean(check_avg_w))

    dataset = Dataset(input_file)
    input_data = np.array(dataset.variables[variable])
    dataset.close()
    ny = 1375
    nx = 1375
    array_avg = np.zeros((1375, 1375))
    array_avg_w = np.zeros((1375, 1375))
    mini_weigths = np.array([[0.25, 0.5, 0.25], [0.5, 1, 0.5], [0.25, 0.5, 0.25]])
    for y in range(ny):
        # print(y)
        for x in range(nx):
            yini = y - 1
            yfin = y + 2
            xini = x - 1
            xfin = x + 2
            yiniw = 0
            yfinw = 3
            xiniw = 0
            xfinw = 3
            if y == 0:
                yini = 0
                yiniw = 1
            if x == 0:
                xini = 0
                xiniw = 1
            if y == (ny - 1):
                yfin = ny
                yfinw = 2
            if x == (nx - 1):
                xfin = nx
                xfinw = 2

            mini_array = input_data[0, yini:yfin, xini:xfin]
            mini_array[mini_array == -999] = np.nan

            mini_weigths_w = mini_weigths[yiniw:yfinw, xiniw:xfinw]

            mini_weigths_w[mini_array == -999] = np.nan

            mini_array_w = np.multiply(mini_array, mini_weigths_w)
            array_avg_w[y, x] = np.nansum(mini_array_w) / np.nansum(mini_weigths_w)

            array_avg[y, x] = np.nanmean(mini_array)

    return array_avg, array_avg_w


def set_data_window_3x3(data_computation, input_file, variable, indices, spatial_weights, input_data):
    time_min = indices[0]
    yref = indices[1]
    ny = indices[2]
    xref = indices[3]
    nx = indices[4]
    ny_abs = indices[5]
    nx_abs = indices[6]

    if input_data is None:
        dataset = Dataset(input_file)
        input_data = np.array(dataset.variables[variable])
        dataset.close()

    limit_output = False
    if data_computation.shape[1] == ny and data_computation.shape[2] == nx and ny_abs != ny and nx_abs != nx:
        limit_output = True

    # 0,0
    xoutput = xref
    youtput = yref
    nyhere = ny
    nxhere = nx
    xinput = xoutput - 1
    yinput = youtput - 1
    if xoutput == 0:
        xoutput = 1
        xinput = 0
        nxhere = nxhere - 1
    if youtput == 0:
        youtput = 1
        yinput = 0
        nyhere = nyhere - 1
    if limit_output:
        if youtput > 1:
            youtput = 0
        if xoutput > 1:
            xoutput = 0
    data_computation[time_min + 0, youtput:youtput + nyhere, xoutput:xoutput + nxhere] = input_data[0,
                                                                                         yinput:yinput + nyhere,
                                                                                         xinput:xinput + nxhere]
    # 0,1
    xoutput = xref
    youtput = yref
    nxhere = nx
    nyhere = ny
    xinput = xoutput
    yinput = youtput - 1
    if youtput == 0:
        youtput = 1
        yinput = 0
        nyhere = nyhere - 1
    if limit_output:
        xoutput = 0
        if youtput > 1:
            youtput = 0
    data_computation[time_min + 1, youtput:youtput + nyhere, xoutput:xoutput + nxhere] = input_data[0,
                                                                                         yinput:yinput + nyhere,
                                                                                         xinput:xinput + nxhere]
    # 0,2
    xoutput = xref
    youtput = yref
    nxhere = nx
    nyhere = ny
    xinput = xoutput + 1
    yinput = youtput - 1
    if youtput == 0:
        youtput = 1
        yinput = 0
        nyhere = nyhere - 1
    if (xoutput + nx) == nx_abs:
        nxhere = nxhere - 1
    if limit_output:
        if youtput > 1:
            youtput = 0
        xoutput = 0
    data_computation[time_min + 2, youtput:youtput + nyhere, xoutput:xoutput + nxhere] = input_data[0,
                                                                                         yinput:yinput + nyhere,
                                                                                         xinput:xinput + nxhere]
    # 1,0
    xoutput = xref
    youtput = yref
    nxhere = nx
    nyhere = ny
    xinput = xoutput - 1
    yinput = youtput
    if xoutput == 0:
        xoutput = 1
        xinput = 0
        nxhere = nxhere - 1
    if limit_output:
        if xoutput > 1:
            xoutput = 0
        youtput = 0
    data_computation[time_min + 3, youtput:youtput + nyhere, xoutput:xoutput + nxhere] = input_data[0,
                                                                                         yinput:yinput + nyhere,
                                                                                         xinput:xinput + nxhere]
    # 1,1
    xoutput = xref
    youtput = yref
    nxhere = nx
    nyhere = ny
    xinput = xoutput
    yinput = youtput
    if limit_output:
        xoutput = 0
        youtput = 0
    data_computation[time_min + 4, youtput:youtput + nyhere, xoutput:xoutput + nxhere] = input_data[0,
                                                                                         yinput:yinput + nyhere,
                                                                                         xinput:xinput + nxhere]
    # 1,2
    xoutput = xref
    youtput = yref
    nxhere = nx
    nyhere = ny
    xinput = xoutput + 1
    yinput = youtput
    if (xoutput + nx) == nx_abs:
        nxhere = nxhere - 1
    if limit_output:
        xoutput = 0
        youtput = 0
    data_computation[time_min + 5, youtput:youtput + nyhere, xoutput:xoutput + nxhere] = input_data[0,
                                                                                         yinput:yinput + nyhere,
                                                                                         xinput:xinput + nxhere]
    # 2,0
    xoutput = xref
    youtput = yref
    nxhere = nx
    nyhere = ny
    xinput = xoutput - 1
    yinput = youtput + 1
    if xoutput == 0:
        xoutput = 1
        xinput = 0
        nxhere = nxhere - 1
    if (youtput + ny) == ny_abs:
        nyhere = nyhere - 1
    if limit_output:
        if xoutput > 1:
            xoutput = 0
        youtput = 0
    data_computation[time_min + 6, youtput:youtput + nyhere, xoutput:xoutput + nxhere] = input_data[0,
                                                                                         yinput:yinput + nyhere,
                                                                                         xinput:xinput + nxhere]
    # 2,1
    xoutput = xref
    youtput = yref
    nxhere = nx
    nyhere = ny
    xinput = xoutput
    yinput = youtput + 1
    if (youtput + ny) == ny_abs:
        nyhere = nyhere - 1
    if limit_output:
        xoutput = 0
        youtput = 0
    data_computation[time_min + 7, youtput:youtput + nyhere, xoutput:xoutput + nxhere] = input_data[0,
                                                                                         yinput:yinput + nyhere,
                                                                                         xinput:xinput + nxhere]
    # 2,2
    xoutput = xref
    youtput = yref
    nxhere = nx
    nyhere = ny
    xinput = xoutput + 1
    yinput = youtput + 1
    if (xoutput + nx) == nx_abs:
        nxhere = nxhere - 1
    if (youtput + ny) == ny_abs:
        nyhere = nyhere - 1
    if limit_output:
        xoutput = 0
        youtput = 0
    data_computation[time_min + 8, youtput:youtput + nyhere, xoutput:xoutput + nxhere] = input_data[0,
                                                                                         yinput:yinput + nyhere,
                                                                                         xinput:xinput + nxhere]
    weighted_avg = None
    if spatial_weights is not None:
        data_avg = data_computation[time_min:time_min + 9, yref:yref + ny, xref:xref + nx]
        data_avg[data_avg == -999] = np.nan
        spatial_weights[np.isnan(data_avg)] = np.nan
        weighted_data = np.multiply(data_avg, spatial_weights)
        weighted_data_sum = np.nansum(weighted_data, axis=0)
        spatial_weights_sum = np.nansum(spatial_weights, axis=0)
        weighted_avg = weighted_data_sum / spatial_weights_sum

    return data_computation, weighted_avg


def get_temporal_weight_filter(options, temporal_window, ny, nx):
    temporal_weights = np.ones((temporal_window, ny, nx))
    increm = int(temporal_window / 2)
    wsum = 0
    for idx, ir in zip(range(temporal_window), range(increm * (-1), increm + 1)):
        weight = ((increm + 1) - np.abs(ir)) / (increm + 1)
        temporal_weights[idx, :, :] = weight
        wsum = wsum + weight

    if ny == 1 and nx == 1:
        temporal_weights = temporal_weights.flatten()

    return temporal_weights


def get_spatial_weight_filter_3x3(options, ny, nx):
    spatial_weights = np.ones((9, ny, nx))
    if options is None:
        spatial_weights[0, :, :] = 0.25
        spatial_weights[1, :, :] = 0.5
        spatial_weights[2, :, :] = 0.25
        spatial_weights[3, :, :] = 0.5
        spatial_weights[4, :, :] = 1.0
        spatial_weights[5, :, :] = 0.5
        spatial_weights[6, :, :] = 0.25
        spatial_weights[7, :, :] = 0.5
        spatial_weights[8, :, :] = 0.25
    if ny == 1 and nx == 1:
        spatial_weights = spatial_weights.flatten()
    return spatial_weights


def mainV2():
    # ##version 2 DEPRECATED
    # date_list = get_input_file_list_day(date_here, temporal_window, input_path, year_ini, year_end)
    # for y in range(0, ny, nywindow):
    #     for x in range(0, nx, nxwindow):
    #         limits = get_limits(y, x, nywindow, nxwindow, ny, nx)
    #         nywindow_here = limits[1] - limits[0]
    #         nxwindow_here = limits[3] - limits[2]
    #         options_clim_here = options_clim.copy()
    #         options_clim_here['nywindow'] = nywindow_here
    #         options_clim_here['nxwindow'] = nxwindow_here
    #         options_clim_here['yref'] = limits[0]
    #         options_clim_here['xref'] = limits[2]
    #
    #         make_clim_multi_day_v2(input_path, output_path, date_here, options_clim_here, date_list, ny, nx)
    #
    # date_here_str = options[ref]['date']
    # date_here = dt.strptime(date_here_str, '%Y-%m-%d')
    return
def make_clim_multi_day_v2(input_path, output_path, date_here, options, date_list, ny, nx):
    if options is None:
        temporal_window = 11
        spatial_window = 3
        xref = 0
        yref = 0
        nywindow = 100
        nxwindow = 100
        use_weights = True
        variable = 'CHL'
        ny = -1
        nx = -1
        file_ref = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/GRID_FILES/ArcGrid_65_90_4KM_GridBase.nc'
        year_ini = 1998
        year_end = 2022
    else:
        temporal_window = options['temporal_window']  # 11
        spatial_window = options['spatial_window']  # 3
        nywindow = options['nywindow']
        nxwindow = options['nxwindow']
        use_weights = options['use_weights']
        variable = options['variable']
        ny = options['ny']
        nx = options['nx']
        file_ref = options['file_ref']
        year_ini = options['year_ini']  # 1998
        year_end = options['year_end']  # 2022
        yref = options['yref']
        xref = options['xref']
    input_name = options['input_name']
    if args.verbose:
        print(f'[INFO]**********************************************************')
        print(f'[INFO] YRef: {yref} XRef: {xref}  NYWindow: {nywindow} NXWindow: {nxwindow}')
        print(f'[INFO] Temporal window: {temporal_window} Spatial window: {spatial_window}')
    if date_list is None:
        date_list = get_input_file_list_day(date_here, temporal_window, input_path, year_ini, year_end, input_name)
    if ny == -1 and nx == -1:
        for date_str in date_list:
            input_file = date_list[date_str]['input_file']
            if input_file is not None:
                ny, nx = get_ny_nx(input_file, variable)
    nyears = (year_end - year_ini) + 1
    ndata = nyears * temporal_window * spatial_window * spatial_window

    data_computation = np.zeros((ndata, nywindow, nxwindow))
    data_weights = np.ones((ndata, nywindow, nxwindow))
    if use_weights:
        if args.verbose:
            print('[INFO] Retrieving weights...')
        spatial_weights = get_spatial_weight_filter_3x3(None, 1, 1)
        temporal_weights = get_temporal_weight_filter(None, temporal_window, 1, 1)
        index = 0
        for iyear in range(nyears):
            for itemporal in range(len(temporal_weights)):
                for ispatial in range(len(spatial_weights)):
                    data_weights[index, :, :] = spatial_weights[ispatial] * temporal_weights[itemporal]
                    index = index + 1

    if args.verbose:
        print('[INFO] Getting data...')
    data_computation[:] = -999.0
    itime = 0
    idate = 0
    year_ref = -1
    for date_str in date_list:
        date_h = dt.strptime(date_str, '%Y-%m-%d')
        year = date_h.year
        iyear = year - year_ini
        input_file = date_list[date_str]['input_file']
        input_data = None
        # input_data = date_list[date_str]['input_array']
        # if input_data is None:
        #     dataset = Dataset(input_file)
        #     input_data = np.array(dataset.variables[variable])
        #     dataset.close()
        if args.verbose and year != year_ref:
            year_ref = year
            print(
                f'[INFO] Date: {date_str} Index all: {itime}  Index date: {idate}  Index year: {iyear} Input file: {input_file}')
        if input_file is None:
            itime = itime + (spatial_window * spatial_window)
            idate = idate + 1
            continue

        # print(itime,yref,nywindow,xref,nxwindow,ny,nx)
        indices = [itime, yref, nywindow, xref, nxwindow, ny, nx]
        data_computation, weighted_avg = set_data_window_3x3(data_computation, input_file, variable, indices, None,
                                                             input_data)
        itime = itime + (spatial_window * spatial_window)
        idate = idate + 1

    # data_computation[data_computation == -999.0] = np.nan
    # data_weights[data_computation == -999.0] = np.nan
    #
    #
    # # r, c = np.unravel_index(np.nanargmax(std), std.shape)
    # # print(std[r, c])
    # # print(data_computation[:, r, c])
    #
    # # avg = np.nanmean(data_computation[:, :, :], axis=0)
    # std = np.std(data_computation[:, :, :], axis=0)
    # median = np.median(data_computation[:, :, :], axis=0)
    #
    # prod = data_computation * data_weights
    # weighted_avg = np.sum(prod) / np.sum(data_weights)
    # nofobs = np.count_nonzero(~np.isnan(data_computation), axis=0)
    #
    # # min_array = np.nanmin(data_computation[:, :, :], axis=0)
    # # max_array = np.nanmax(data_computation[:, :, :], axis=0)

    if args.verbose:
        print('[INFO] COMPLETED')



if __name__ == '__main__':
    main()
