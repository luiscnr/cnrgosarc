import argparse
import os
from datetime import timedelta

import numpy as np

from arc_multi_sources import ARC_MULTI_SOURCES
import warnings

warnings.filterwarnings(action='ignore', category=ResourceWarning)

parser = argparse.ArgumentParser(description="Artic resampler")
parser.add_argument("-m", "--mode", help="Mode",
                    choices=["CHECKPY", "CHECK_ALGORITHM", "GRID", "RESAMPLE", "RESAMPLEFILE", "INTEGRATE", "CHLA", "KD490", "QL",
                             "MONTHLY_CHLA", "MONTHLY_KD490", "MONTHLY_RRS_TEST","QI","CORRECT_TIMESTAMP","CHECK_TIMESTAMP"],
                    required=True)
parser.add_argument("-p", "--product", help="Input product (testing)")
parser.add_argument('-i', "--inputpath", help="Input directory")
parser.add_argument('-o', "--outputpath", help="Output file")
parser.add_argument('-tp', "--temp_path", help="Temporary directory")
parser.add_argument('-sd', "--start_date", help="Start date (yyyy-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date (yyyy-mm-dd")
parser.add_argument('-c', "--config_file", help="Configuration file (Default: arc_config.ini)")
parser.add_argument('-bf', "--base_file", help="Create base file for the following mode: RESAMPLE,INTEGRATE.")
parser.add_argument('-var', "--variable_plot", help="Variable to be plot using QL mode. Default: CHL",
                    choices=["CHL", "KD490"])
parser.add_argument("-chla_algo", "--chla_algorithm",help="Chl-a algorithm",choices=["SeaSARC","CIAO"],default="CIAO")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
args = parser.parse_args()


def only_test():
    path = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/S3A_OL_2_WFR____20230715T232353_20230715T232653_20230717T115842_0179_087_301_1800_MAR_O_NT_003.SEN3'
    from netCDF4 import Dataset
    import numpy as np
    from kd_algorithm import KD_ALGORITHMS
    path_kd = os.path.join(path, 'trsp.nc')
    dataset_kd = Dataset(path_kd)
    yini = 2100
    yfin = 2800
    xini = 4000
    xfin = 4700
    kd_array_olci = np.array(dataset_kd.variables['KD490_M07'][yini:yfin, xini:xfin])
    scale_kd = dataset_kd.variables['KD490_M07'].scale_factor
    offset_kd = dataset_kd.variables['KD490_M07'].add_offset
    dataset_kd.close()

    path_443 = os.path.join(path, 'Oa03_reflectance.nc')
    dataset_443 = Dataset(path_443)
    array_443 = np.array(dataset_443.variables['Oa03_reflectance'][yini:yfin, xini:xfin])

    path_490 = os.path.join(path, 'Oa04_reflectance.nc')
    dataset_490 = Dataset(path_490)
    array_490 = np.array(dataset_490.variables['Oa04_reflectance'][yini:yfin, xini:xfin])

    path_510 = os.path.join(path, 'Oa05_reflectance.nc')
    dataset_510 = Dataset(path_510)
    array_510 = np.array(dataset_510.variables['Oa05_reflectance'][yini:yfin, xini:xfin])

    path_560 = os.path.join(path, 'Oa06_reflectance.nc')
    dataset_560 = Dataset(path_560)
    array_560 = np.array(dataset_560.variables['Oa06_reflectance'][yini:yfin, xini:xfin])

    path_chla = os.path.join(path, 'chl_oc4me.nc')
    dataset_chla = Dataset(path_chla)
    scale_chl = dataset_chla.variables['CHL_OC4ME'].scale_factor
    offset_chl = dataset_chla.variables['CHL_OC4ME'].add_offset
    chla_array_olci = np.array(dataset_chla.variables['CHL_OC4ME'][yini:yfin, xini:xfin])

    kda = KD_ALGORITHMS('OK2-560')
    kd_array_new = kda.compute_kd490_ok2_560(array_490, array_560, chla_array_olci)
    chl_array_new = kda.compute_chla_ocme4(array_443, array_490, array_510, array_560, chla_array_olci)

    ##KD TRANSFORMATION FOR COMPARISON
    # invert
    kd_array_new[kd_array_olci != 255] = np.log10(kd_array_new[kd_array_olci != 255])
    kd_array_new[kd_array_olci != 255] = (kd_array_new[kd_array_olci != 255] - offset_kd) / scale_kd
    # round
    kd_array_new[kd_array_olci != 255] = np.ceil(kd_array_new[kd_array_olci != 255])
    # values shoud be between 0 and 254 (as 255 is used as invalid value)
    kd_array_new[kd_array_new < 0] = 0.0
    kd_array_new[kd_array_new > 254] = 254.0
    kd_array_new[kd_array_olci == 255] = 255  ##invalid value, from kd_olci
    kd_array_new = kd_array_new.astype(np.uint8).astype(np.float32)
    ##reconvert
    kd_array_new[kd_array_olci != 255] = (kd_array_new[kd_array_olci != 255] * scale_kd) + offset_kd
    kd_array_new[kd_array_olci != 255] = np.power(10, kd_array_new[kd_array_olci != 255])
    ##masks
    kd_array_new[kd_array_olci == 255] = -999.0
    kd_array_new[chla_array_olci == 255] = -999.0

    ##CHL TRANSFORMATION FOR COMPARISON
    # invert
    chl_array_new[chla_array_olci != 255] = np.log10(chl_array_new[chla_array_olci != 255])
    chl_array_new[chla_array_olci != 255] = (chl_array_new[chla_array_olci != 255] - offset_chl) / scale_chl
    # round
    chl_array_new[chla_array_olci != 255] = np.ceil(chl_array_new[chla_array_olci != 255])
    # values shoud be between 0 and 254 (as 255 is used as invalid value)
    chl_array_new[chl_array_new < 0] = 0.0
    chl_array_new[chl_array_new > 254] = 254.0
    chl_array_new[chla_array_olci == 255] = 255  ##invalid value, from chl_olci
    chl_array_new = chl_array_new.astype(np.uint8).astype(np.float32)
    ##reconvert
    chl_array_new[chla_array_olci != 255] = (chl_array_new[chla_array_olci != 255] * scale_chl) + offset_chl
    chl_array_new[chla_array_olci != 255] = np.power(10, chl_array_new[chla_array_olci != 255])
    ##masks
    chl_array_new[chla_array_olci == 255] = -999.0

    ##NON LOG10 TRANSFORMED VALUES FOR KD Y CHLA FROM OLCI
    kd_array_olci[kd_array_olci != 255] = np.power(10, kd_array_olci[kd_array_olci != 255])
    chla_array_olci[chla_array_olci != 255] = np.power(10, chla_array_olci[chla_array_olci != 255])

    file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/chlacomparison_1.csv'
    f1 = open(file_out, 'w')
    f1.write('Y;X;KD_OLCI;KD_NEW;CHLA_OLCI;CHLA_NEW;RATIO_KD;RATIO_CHLA')
    for y in range(kd_array_olci.shape[0]):
        for x in range(kd_array_olci.shape[1]):
            ypoint = y + yini
            xpoint = x + xini
            if kd_array_new[y, x] != -999.0:
                val_kd_olci = kd_array_olci[y, x]
                val_kd_new = kd_array_new[y, x]
                val_chla_olci = chla_array_olci[y, x]
                val_chla_new = chl_array_new[y, x]
                ron_kd = kd_array_olci[y, x] / kd_array_new[y, x]
                ron_chla = chla_array_olci[y, x] / chl_array_new[y, x]
                # if y==691 and x==4:
                #     print(f'OLCI: {val_old:.4f};COMPUTED: {val_new:.4f}; CHLA: {val_chla:.4f}; RON: {ron:.4f}')
                line = f'{ypoint};{xpoint};{val_kd_olci};{val_kd_new};{val_chla_olci};{val_chla_new};{ron_kd};{ron_chla}'
                f1.write('\n')
                f1.write(line)
    f1.close()

    dataset_443.close()
    dataset_490.close()
    dataset_510.close()
    dataset_560.close()
    dataset_chla.close()

    return True


def do_test():
    file_check = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_COMPARISON_OLCI_MULTI/ALGORITHMS/C2016122_rrs-arc-4km.nc'
    from netCDF4 import Dataset
    from kd_algorithm import KD_ALGORITHMS
    dcheck = Dataset(file_check)
    var_names = ['RRS443', 'RRS490', 'RRS510', 'RRS560', 'RRS665']

    val443 = dcheck.variables['RRS443'][0, 357, 1136]
    val490 = dcheck.variables['RRS490'][0, 357, 1136]
    val510 = dcheck.variables['RRS510'][0, 357, 1136]
    val560 = dcheck.variables['RRS560'][0, 357, 1136]
    # val665 = dcheck.variables['RRS560'][0, 357, 1136]

    dcheck.close()
    kda = KD_ALGORITHMS('OK2-560')
    res = kda.compute_kd_param(val443, val490, val510, val560)
    print(res[1])
    file_output = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_COMPARISON_OLCI_MULTI/ALGORITHMS/C2016122_kd490-arc-4km.nc'
    dout = Dataset(file_output)
    print(dout.variables['KD490'][0, 357, 1136])
    dout.close()

    return True


def do_global_grid_monthly():
    input_file = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/SOURCES/2019/07/20190701_c3s_obs-oc_glo_bgc-reflectance_my_l3-multi-4km_P1D.nc'
    output_file = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/GRID_FILES/GlobalGrid_Base.nc'

    from netCDF4 import Dataset
    input_dataset = Dataset(input_file)
    ncout = Dataset(output_file, 'w')

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    variable = input_dataset.variables['latitude']
    ncout.createVariable('lat', variable.datatype, variable.dimensions, zlib=True, shuffle=True, complevel=6)
    ncout['lat'][:] = input_dataset['latitude'][:]
    variable = input_dataset.variables['longitude']
    ncout.createVariable('lon', variable.datatype, variable.dimensions, zlib=True, shuffle=True, complevel=6)
    ncout['lon'][:] = input_dataset['longitude'][:]
    variable = input_dataset.variables['time']
    ncout.createVariable('time', variable.datatype, variable.dimensions, zlib=True, shuffle=True, complevel=6)
    ncout['time'][:] = input_dataset['time'][:]

    variable = input_dataset.variables['RRS510']
    ncout.createVariable('RRS510', variable.datatype, variable.dimensions, fill_value=-999.0, zlib=True, shuffle=True,
                         complevel=6)
    ncout.createVariable('RRS510_count', variable.datatype, variable.dimensions, fill_value=-999.0, zlib=True,
                         shuffle=True,
                         complevel=6)
    ncout.createVariable('RRS510_error', variable.datatype, variable.dimensions, fill_value=-999.0, zlib=True,
                         shuffle=True,
                         complevel=6)

    input_dataset.close()
    ncout.close()

    return True


def check_point():
    from datetime import datetime as dt
    from netCDF4 import Dataset
    dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/OUTPUT/2019'
    dir_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/SOURCES/2019/07'
    for day in range(1, 32):
        date_here = dt(2019, 7, day)
        yyyy = '2019'
        dd = date_here.strftime('%d')
        mm = date_here.strftime('%m')
        jjj = date_here.strftime('%j')
        fout = os.path.join(dir_out, jjj, f'C{yyyy}{jjj}_rrs-arc-4km.nc')
        fout = os.path.join(dir_out, f'{yyyy}{mm}{dd}_c3s_obs-oc_glo_bgc-reflectance_my_l3-multi-4km_P1D.nc')
        dataset = Dataset(fout)
        array = np.array(dataset.variables['RRS510'])
        # print(day,'-->',array[0,1070,826])
        print(day, '-->', array[0, 359, 7080])
        dataset.close()

    return True


def check_sources():
    from datetime import datetime as dt
    file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/missing_sources.csv'
    fout = open(file_out,'w')
    input_path = '/mnt/c/DATA_LUIS'
    moi_user = 'lgonzalezvilas'
    moi_pass = 'MegaRoma17!'
    ams = ARC_MULTI_SOURCES(input_path, None, moi_user, moi_pass, False, True)

    for year in range(1997,2023):
        for month in range(1,13):
            if year==1997 and month<9:
                continue
            print('-->',year,month)
            missing_dates = ams.check_month_files_ftp(year,month)
            for date in missing_dates:
                fout.write(date.strftime('%Y-%m-%d'))
                fout.write('\n')
    # date_ref = dt(1997, 9, 4)
    # date_end = dt(2023, 12, 31)
    # while date_ref <= date_end:
    #     if date_ref.day == 1:
    #         print('-->', date_ref)
    #     b = ams.check_file_ftp(date_ref)
    #     if not b:
    #         print('NOT AVAILABLE')
    #         fout.write('\n')
    #         fout.write(date_ref.strftime('%Y-%m-%d'))
    #     date_ref = date_ref + timedelta(hours=24)

    fout.close()

    return True

def check_coverage():
    file_grid = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/GRID_FILES/ArcGrid_65_90_4KM_GridBase.nc'
    from netCDF4 import Dataset
    import numpy as np
    dataset = Dataset(file_grid)

    lat_array = np.array(dataset.variables['lat'])
    lon_array = np.array(dataset.variables['lon'])
    from global_land_mask import globe

    ncoverage = 0
    for y in range(0,1375):
        print(y, 'de 1375')
        for x in range(0,1375):
            lat_here = lat_array[y,x]
            lon_here = lon_array[y,x]
            if lat_here<65.0:
                continue
            if globe.is_land(lat_here,lon_here):
                continue
            ncoverage = ncoverage + 1


    print(ncoverage)

    return True

def do_test_array():
    file_csv = '/mnt/c/Users/LuisGonzalez/OneDrive - NOLOGIN OCEANIC WEATHER SYSTEMS S.L.U/CNR/OCTAC_WORK/ARC_WORK/CIAO/CIAO_Test_Dataset.csv'
    import pandas as pd
    df = pd.read_csv(file_csv,sep=',')
    ndata = len(df.index)
    matrix = np.ones((ndata,6))
    col_names = ['','DOY','Rrs_443','Rrs_490','Rrs_510','Rrs_560']
    for idx in range(1,6):
        if idx>1:
            matrix[:, idx] = np.log10(np.array(df[col_names[idx]]))
        else:
            matrix[:,idx] = np.array(df[col_names[idx]])

    from arc_gpr_model import ARC_GPR_MODEL
    fModel = '/home/lois/PycharmProjects/cnrgosarc/CIAO_Algorithm.json'
    gModel = ARC_GPR_MODEL(fModel,True)
    model_chl = gModel.compute_chla_ciao_from_matrix(matrix)

    for idx in range(ndata):
        val_fin = gModel.compute_chla_impl(matrix[idx,1:])
        computed = 10 ** val_fin
        print(computed,'??',model_chl[idx],'-------------------->',f'{computed/model_chl[idx]:.8f}')


    return True

def make_sbatch():
    from datetime import datetime as dt
    from datetime import timedelta
    ncores = 10
    index_job = 1
    ifile = 1
    work_date = dt(2018,1,1)
    end_date = dt(2018,12,31)
    base_file = '/home/gosuser/Processing/gos-oc-processingchains_v202411/arcProcessing/config/cnrarc_config_hpc01_operational_DT_'
    pr = f" | awk '{{print $NF}}')"

    print('#!/bin/bash')
    print('')
    print(f'tfile=/home/gosuser/Processing/gos-oc-processingchains_v202411/arcProcessing/s3olciProcessing/mail_2018.txt')
    print(f'subject="Multiple OLCI ARC processing {work_date.strftime("%Y-%m-%d")} {end_date.strftime("%Y-%m-%d")}"')
    print(f'mailrcpt="luis.gonzalezvilas@artov.ismar.cnr.it"')
    print('')
    print(f'echo "Launching multiple OLCI ARC processing from {work_date.strftime("%Y-%m-%d")} to {end_date.strftime("%Y-%m-%d")}">$tfile')
    print(f'echo "">>$tfile')
    print('')
    while work_date <= end_date:

        work_date_str = work_date.strftime('%Y-%m-%d')
        if index_job<=ncores:
            line=f'job{index_job}=$(sbatch make_processing.slurm {base_file}{ifile}.ini {work_date_str} {work_date_str})'
            print(line)
            line=f'job{index_job}id=$(echo "$job{index_job}"{pr}'
            print(line)
        else:
            index_job_wait = index_job-ncores
            line = f'job{index_job}=$(sbatch --dependency=afterany:$job{index_job_wait}id make_processing.slurm {base_file}{ifile}.ini {work_date_str} {work_date_str})'
            print(line)
            line = f'job{index_job}id=$(echo "$job{index_job}"{pr}'
            print(line)

        print(f'echo "Date: {work_date_str} ID: $job{index_job}id">>$tfile')

        index_job = index_job + 1
        ifile = ifile + 1
        if ifile>ncores:
            ifile = 1

        work_date = work_date + timedelta(days=1)

    print(f'cat $tfile | mail -s "$subject" "$mailrcpt"')
    return True

def run_ql(arc_opt, start_date, end_date):
    options = arc_opt.get_ql_options()
    if options is None:
        return
    output_type = arc_opt.get_value_param('QL', 'output_type', 'CHL', 'str')
    # if not output_type == 'CHL':
    #     return

    print(f'[INFO] Output type for quick looks: {output_type}')

    name_file_format_default = None
    name_file_date_format_default = '%Y%j'
    if output_type == 'CHL':
        name_file_format_default = 'C$DATE$_chl-arc-4km.nc'
        output_var = 'CHL'
    elif output_type=='CHL_CLIMA':
        name_file_format_default = '1998$DATE$_2022$DATE$_chl_arc_multi_clima.nc'
        name_file_date_format_default = '%m%d'
        output_var = 'MEDIAN'
    else:
        output_var = output_type

    name_file_format = arc_opt.get_value_param('QL', 'name_file_format', name_file_format_default, 'str')
    name_file_date_format = arc_opt.get_value_param('QL', 'name_file_date_format_default',
                                                    name_file_date_format_default, 'str')

    if start_date is None or end_date is None:
        start_date = options['start_date']
        end_date = options['end_date']
    date_run = start_date

    from arc_mapinfo import ArcMapInfo
    ami = ArcMapInfo(arc_opt, args.verbose)
    ami.set_area_definition('polar_stereographic_4km')

    while date_run <= end_date:
        if args.verbose:
            print('*****************************')
            print(f'[INFO] Date: {date_run}')
        make_ql = True
        input_path = arc_opt.get_folder_date(options['input_path'], options['input_path_organization'], date_run, False)
        date_file_str = date_run.strftime(name_file_date_format)
        name_file = name_file_format.replace('$DATE$', date_file_str)
        input_file = os.path.join(input_path, name_file)
        if not os.path.exists(input_file):
            print(f'[WARNING] Input file {input_file} for date {date_run} is not available. Skiping...')
            make_ql = False
        output_path = arc_opt.get_folder_date(options['output_path'], options['output_path_organization'], date_run,
                                              True)
        if output_path is None:
            print(f'[WARNING] Output path {input_path} for date {date_run} is not available. Skiping...')
            make_ql = False

        output_name = f'{name_file[:-3]}_{output_type}.png'
        output_file = os.path.join(output_path, output_name)
        if os.path.exists(output_file):
            print(f'[INFO] Output file {output_file} already exists. Skipping...')
            make_ql = False

        if make_ql:
            if args.verbose:
                print(f'[INFO] Input file: {input_file}')
                print(f'[INFO] Output file: {output_file}')
            #ami.save_quick_look_fdata(output_file, input_file, output_type)

            ami.save_full_fdata(output_file, input_file, output_var)
        date_run = date_run + timedelta(hours=24)

def main():
    # if do_global_grid_monthly():
    #     return
    # if do_test():
    #     return
    # if only_test():
    #     return
    # if check_point():
    #     return
    # if check_sources():
    #     return

    # if check_coverage():
    #     return
    # if do_test_array():
    #     return
    # if make_sbatch():
    #     return

    print('[INFO] Started Artic Processing Tool [MULTI 4 KM]')
    if args.mode == "CHECKPY":
        check_py()
        return

    if args.mode == 'CHECK_ALGORITHM':
        return

    if args.mode == 'GRID' and args.outputpath:  ##creating single file grid
        from arc_mapinfo import ArcMapInfo
        if args.verbose:
            print('[INFO] Create grid file')
        ami = ArcMapInfo(None, args.verbose)
        ami.set_area_definition('polar_stereographic_4km')
        file_out = args.outputpath
        if not os.path.exists(file_out):
            ami.create_nc_filegrid(file_out, False, True)
        else:
            print(f'[WARNING] Grid file already exists. Skipping...')
        return

    ##BASE FILES (WITHOUT VARIABLES)
    if args.mode == 'INTEGRATE' and args.base_file and args.config_file and args.outputpath:
        output_file = args.outputpath
        if not os.path.isdir(os.path.dirname(output_file)) or not output_file.endswith('.nc'):
            print(f'[ERROR] Output file: {output_file} is not valid')
            return
        # file_at = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/CONFIG_FILES/global_attributes.ini'
        file_at = args.config_file
        if not os.path.isfile(file_at):
            print(f'[ERROR] Attributes files: {file_at} is not valid')
        file_base = args.base_file
        if not os.path.isfile(file_at):
            print(f'[ERROR] File base: {file_base} is not valid')

        from arc_integration import ArcIntegration
        if args.verbose:
            print('[INFO] Started creation of base files for integration...')
        params = ['RRS', 'TRANSP']
        param = None
        for p in params:
            if output_file.find(p) > 0:
                param = p
        modes = ['NR', 'NT']
        mode = None
        for m in modes:
            if output_file.find(m) > 0:
                mode = m
        if param is None:
            print(f'[ERROR] Param. is not defined in output file name. It should be RRS, TRANSP')
            return
        if mode is None:
            print(f'[ERROR] Mode is not defined in output file name. It should be NR or NT')
            return
        arcInt = ArcIntegration(None, args.verbose, None, param, file_at)
        arcInt.ami.set_area_definition('polar_stereographic_4km')
        arcInt.ami.ifile_base = file_base
        dout = arcInt.create_nc_file_out(output_file, mode)
        dout.close()

    # BASE FILE WITH VARIABLES
    if args.mode == 'CHLA' and args.base_file and args.config_file and args.outputpath:
        output_file = args.outputpath
        if not os.path.isdir(os.path.dirname(output_file)) or not output_file.endswith('.nc'):
            print(f'[ERROR] Output file: {output_file} is not valid. It should be a nc file in an existing directory')
            return
        # file_at = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/CONFIG_FILES/global_attributes.ini'
        file_at = args.config_file
        if not os.path.isfile(file_at):
            print(f'[ERROR] Attributes files: {file_at} does not exist or is not a valid file')

        file_base = args.base_file
        if not os.path.isfile(file_base):
            print(f'[ERROR] File base: {file_base} does not exist of is not a valid file')

        mode = 'NR' if output_file.find('NR')>0 else 'NT' if output_file.find('NT')>0 else None
        if mode is None:
            print(f'[ERROR] Mode is not defined in output file name. It should be NR or NT')
            return

        from arc_processing import ArcProcessing
        output_type  = f'CHLA_{args.chla_algorithm}'
        arcProc = ArcProcessing(None, args.verbose, output_type, file_at)

        arcProc.ami.set_area_definition('polar_stereographic_4km')
        dout = arcProc.create_nc_file_out(output_file, file_base, mode)

        if dout is not None:
            dout.close()

        return

    if args.mode == 'MONTHLY_CHLA' and args.base_file and args.config_file and args.outputpath:
        output_file = args.outputpath
        if not os.path.isdir(os.path.dirname(output_file)) or not output_file.endswith('.nc'):
            print(f'[ERROR] Output file: {output_file} is not valid')
            return
        # file_at = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/CONFIG_FILES/global_attributes.ini'
        file_at = args.config_file
        if not os.path.isfile(file_at):
            print(f'[ERROR] Attributes files: {file_at} is not valid')
        file_base = args.base_file
        if not os.path.isfile(file_at):
            print(f'[ERROR] File base: {file_base} is not valid')
        modes = ['NR', 'NT']
        mode = None
        for m in modes:
            if output_file.find(m) > 0:
                mode = m
        if mode is None:
            print(f'[ERROR] Mode is not defined in output file name. It should be NR or NT')
            return
        from arc_processing import ArcProcessing
        arcProc = ArcProcessing(None, args.verbose, 'CHLA', file_at)
        arcProc.create_nc_file_out_month(output_file, file_base, mode)

    # RESAMPLE FILES
    if args.mode == 'RESAMPLEFILE' and args.product and args.outputpath:
        input_file = args.product
        if not os.path.isfile(input_file) or not input_file.endswith('.nc'):
            print(f'[ERROR] {args.input_file} is not a valid input file')
        if input_file.find('reflectance') > 0:
            base_file = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/GRID_FILES/ArcGrid_65_90_4KM_RRS_NT_Base.nc'
        if input_file.find('transp') > 0:
            base_file = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/GRID_FILES/ArcGrid_65_90_4KM_TRANSP_NT_Base.nc'
        from datetime import datetime as dt
        from netCDF4 import Dataset
        date_here_str = input_file.split('/')[-1][:8]
        date_here = dt.strptime(date_here_str, '%Y%m%d')
        dbase = Dataset(base_file)
        dname = dbase.title
        dbase.close()
        if os.path.isdir(args.outputpath):
            name_output = f'{date_here_str}_{dname}.nc'
            output_file = os.path.join(args.outputpath, name_output)
        elif os.path.isdir(os.path.dirname(args.outputpath)) and args.outputpath.endswith('nc'):
            output_file = args.outputpath
        else:
            print(f'[ERROR] {args.outputpath} is not a valid output file or directory')
            return

        from arc_mapinfo import ArcMapInfo
        ami = ArcMapInfo(None, args.verbose)
        ami.set_area_definition('polar_stereographic_4km')
        ami.ifile_base = base_file

        ami.make_resample_from_file_orig_multi(input_file, output_file, date_here)

        return

    if args.mode == 'CORRECT_TIMESTAMP' and args.start_date and args.inputpath and args.outputpath:
        start_date, end_date = get_dates_from_arg()
        if start_date is None or end_date is None:
             return
        if not os.path.isdir(args.inputpath):
            print(f'[ERROR] {args.inputpath} is not a valid directory')
            return
        if not os.path.isdir(args.outputpath):
            try:
                os.mkdir(args.outputpath)
            except:
                print(f'[ERROR] {args.outputpath} is not a valid directory')
                return
        run_correct_time_stamp(start_date,end_date,args.inputpath,args.outputpath)
        return

    if args.mode == 'CHECK_TIMESTAMP' and  args.inputpath:
        start_date, end_date = get_dates_from_arg()
        if start_date is None or end_date is None:
            return
        if not os.path.isdir(args.inputpath):
            print(f'[ERROR] {args.inputpath} is not a valid directory')
            return
        run_correct_time_stamp(start_date,end_date,args.inputpath,None)
        return

    ##FROM HERE, ALL THE MODES REQUIRE CONFIGURATION FILES, BUT DATES COULD BE ALSO PASSED AS ARGS
    if not args.config_file:
        print(f'[ERROR] Config file or input product should be defined for {args.mode} option. Exiting...')
        return
    if not os.path.exists(args.config_file):
        print(f'[ERROR] Config file {args.config_file} does not exist. Exiting...')
        return
    try:
        import configparser
        options = configparser.ConfigParser()
        options.read(args.config_file)
    except:
        print(f'[ERROR] Config file {args.config_file} could not be read. Exiting...')

    from arc_options import ARC_OPTIONS
    arc_opt = ARC_OPTIONS(options)

    # check if dates from args should be used
    start_date = None
    end_date = None
    if args.start_date:
        start_date, end_date = get_dates_from_arg()
        if start_date is None or end_date is None:
            return

    if args.mode == 'RESAMPLE':
        run_resample(arc_opt, start_date, end_date)
        return

    if args.mode == 'CHLA':
        run_chla(arc_opt, start_date, end_date)
        return

    if args.mode == 'KD490':
        run_kd490(arc_opt, start_date, end_date)
        return

    if args.mode == 'MONTHLY_CHLA':
        # operative_mode = arc_opt.get_value_param('PROCESSING', 'operative_mode', False, 'boolean')
        # if operative_mode:
        #     date = check_monthly_operative_mode(arc_opt)
        #     if date is not None:
        #         run_month(arc_opt,'CHLA',date,date)
        # else:
        run_month(arc_opt, 'CHLA', start_date, end_date)
        return

    if args.mode == 'MONTHLY_RRS_TEST':
        run_month(arc_opt, 'RRS510', start_date, end_date)
        return

    if args.mode == 'QL':
        # print('to be done')
        run_ql(arc_opt, start_date, end_date)



def run_correct_time_stamp(start_date,end_date,input_path,output_path):
    if start_date is None or end_date is None:
        return
    date_ref = start_date
    while date_ref <= end_date:
        print(f'[INFO]******************************************************************************->{date_ref}')
        input_path_date = os.path.join(input_path,date_ref.strftime('%Y'),date_ref.strftime('%j'))
        dateyj = date_ref.strftime('%Y%j')
        name = f'C{dateyj}_rrs-arc-4km.nc'
        file_input = os.path.join(input_path_date,name)
        if os.path.exists(file_input):
            if output_path is None: ##only info
                from netCDF4 import Dataset
                from datetime import datetime as dt
                dataset = Dataset(file_input)
                ts = float(dataset.variables['time'][0])
                date_file = dt(1981,1,1)+timedelta(seconds=ts)
                btime = date_file.strftime('%Y-%m-%d')==date_ref.strftime('%Y-%m-%d')
                print(f'[INFO] {date_ref.strftime("%Y-%m-%d")}->{date_file.strftime("%Y-%m-%d")}->{btime}')
                dataset.close()
            else:
                output_path_year = os.path.join(output_path,date_ref.strftime('%Y'))
                if not os.path.exists(output_path_year): os.mkdir(output_path_year)
                output_path_date = os.path.join(output_path_year,date_ref.strftime('%j'))
                if not os.path.exists(output_path_date): os.mkdir(output_path_date)
                file_output = os.path.join(output_path_date, name)
                create_copy_correcting_time_stamp(file_input,file_output,date_ref)

        date_ref = date_ref + timedelta(hours=24)

def create_copy_correcting_time_stamp(file_input,file_output,date_here):
    from datetime import datetime as dt
    from netCDF4 import Dataset
    time_ref = dt(1981,1,1,0,0,0)
    new_ts = (date_here.replace(hour=0,minute=0,second=0,microsecond=0)-time_ref).total_seconds()
    new_ts = np.int32(new_ts)
    input_dataset = Dataset(file_input)
    ncout = Dataset(file_output, 'w', format='NETCDF4')

    ##global attributes
    ncout.setncatts(input_dataset.__dict__)

    # copy dimensions
    for name, dimension in input_dataset.dimensions.items():
        ncout.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    for name, variable in input_dataset.variables.items():
        fill_value = None
        if '_FillValue' in list(variable.ncattrs()):
            fill_value = variable._FillValue

        ncout.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value, zlib=True, complevel=6)

        # copy variable attributes all at once via dictionary
        ncout[name].setncatts(input_dataset[name].__dict__)

        # copy data
        if name == 'time':
            ncout[name][0] = new_ts
        else:
            ncout[name][:] = input_dataset[name][:]

    ncout.close()
    input_dataset.close()

def run_resample(arc_opt, start_date, end_date):
    options = arc_opt.get_resample_options()
    if options is None:
        print(f'[ERROR] Error getting the RESAMPLE options. Please review the config file: {args.config_file}')
        return
    if args.verbose:
        print('[INFO] Resample options -------------------------------------------------------')
        for option in options:
            print(f'[INFO]  {option} -> {options[option]}')
        print('[INFO] ------------------------------------------------------------------------')
    if start_date is not None and end_date is not None:
        date_ref = start_date
        date_fin = end_date
    else:
        date_ref = options['start_date']
        date_fin = options['end_date']
    options = arc_opt.add_moi_credentials(options)
    options['file_base'] = arc_opt.get_value_param('RESAMPLE', 'file_base', None, 'file')
    overwrite = arc_opt.get_value_param('RESAMPLE', 'overwrite', False, 'boolean')
    if options['file_base'] is None:
        print(f'[ERROR] Option file_base is not available in RESAMPLE section, or file does not exist')
        return
    from netCDF4 import Dataset
    base_file = options['file_base']
    climatology_path = options['climatology_path']
    # dbase = Dataset(base_file)
    # dname = dbase.title
    # dbase.close()

    use_myint_sources = options['use_myint_sources']
    ams = ARC_MULTI_SOURCES(options['input_path'], options['input_path_organization'], options['moi_user'],
                            options['moi_pass'], use_myint_sources, args.verbose)

    while date_ref <= date_fin:
        print(f'[INFO]******************************************************************************->{date_ref}')

        date_here_str = date_ref.strftime('%Y%m%d')
        output_path = ams.get_folder_date(date_ref, True, options['output_path'], options['output_path_organization'])
        if output_path is None:
            date_ref = date_ref + timedelta(hours=24)
            print(f'[ERROR] Output path: {output_path} is not available and could not be created')
            continue

        dateyj = date_ref.strftime('%Y%j')
        name_output = f'C{dateyj}_rrs-arc-4km.nc'
        # name_output = f'{date_here_str}_{dname}.nc'
        file_output = os.path.join(output_path, name_output)

        if os.path.exists(file_output) and not overwrite:
            date_ref = date_ref + timedelta(hours=24)
            print(f'[WARNING] Output file: {output_path} already exists. Skipping...')
            continue

        ##implements 3 attemps to download the file
        nattemps = 0
        file_date = None
        while nattemps < 3:
            file_date = ams.get_file_date(date_ref, True)
            if os.path.exists(file_date):
                break
            nattemps = nattemps + 1

        if file_date is None:
            date_ref = date_ref + timedelta(hours=24)
            print(f'[ERROR] File for date: {date_ref} is not available and could not be downloaded')
            continue
        if file_date == 'NOFILE':
            date_ref = date_ref + timedelta(hours=24)
            print(f'[ERROR] File for date: {date_ref} is not available')
            continue
        if file_date == 'INVALIDFILE':
            date_ref = date_ref + timedelta(hours=24)
            print(f'[ERROR] File for date: {date_ref} is not valid, maybe it was not downloaded correctly')
            continue

        from arc_mapinfo import ArcMapInfo
        ami = ArcMapInfo(None, args.verbose)
        ami.set_area_definition('polar_stereographic_4km')
        if climatology_path is not None and os.path.exists(climatology_path):
            ami.climatology_path = climatology_path
        ami.ifile_base = base_file
        if args.verbose:
            print(f'[INFO] Make resample from file {file_date} to file: {file_output}')
        ami.make_resample_from_file_orig_multi(file_date, file_output, date_ref)

        date_ref = date_ref + timedelta(hours=24)


def run_chla(arc_opt, start_date, end_date):
    section = 'CHLA'
    from datetime import timedelta
    options = arc_opt.get_basic_options(section)
    if options is None:
        return
    overwrite = arc_opt.get_value_param(section, 'overwrite', False, 'boolean')
    file_base = arc_opt.get_value_param(section, 'file_base', None, 'file')
    file_att = arc_opt.get_value_param(section, 'file_att', None, 'file')
    chla_algo = args.chla_algorithm if args.chla_algorithm else arc_opt.get_value_param(section,'chla_algorithm','CIAO','str')
    if chla_algo not in ["SeaSARC","CIAO"]:
        print(f'[ERROR] {chla_algo} is not a valid chl-a algorithm. Possible values are: SeaSARC, CIAO')
        return

    if file_base is None:
        print(f'[ERROR] file_base does not exist or is not available in the config file. Please review it')
        return

    timeliness = 'NR' if os.path.basename(file_base).find('NR')>0 else 'NT' if os.path.basename(file_base).find('NT')>0 else None
    if timeliness is None:
        print(f'[ERROR] Timeliness is not defined in the file base. It should be NR or NT')
        return

    from arc_processing import ArcProcessing
    if args.verbose:
        print('[INFO] CHLA PROCESSING OPTIONS:')
        for opt in options:
            print(f'[INFO]  {opt}->{options[opt]}')
        print(f'[INFO]  file_base:->{file_base}')
        print(f'[INFO]  file_att:->{file_att}')
        print(f'[INFO]  overwrite:->{overwrite}')
        print(f'[INFO]  timeliness:->{timeliness}')
        print(f'[INFO]  chl-a algorithm:->{chla_algo}')

    climatology_path = options['climatology_path']

    ##WORKING WITH SINGLE GRANULE, ONLY CHLA
    input_name = arc_opt.get_value_param('PROCESSING', 'name_input', None, 'str')
    if input_name is not None:
        input_file = os.path.join(options['input_path'], input_name)
        if os.path.exists(input_file):
            if args.verbose:
                print(f'[INFO] Working with the single file: {input_file}')
        else:
            print(f'[ERROR] File {input_file} does not exist')
        output_name = arc_opt.get_value_param('PROCESSING', 'name_output', None, 'str')
        if output_name is None:
            output_name = 'SingleOutputChla.nc'
        output_file = os.path.join(options['output_path'], output_name)
        if args.verbose:
            print(f'[INFO] Output file: {output_file}')
        # defining arc_proc, last parameters (file_at) is none because it's defined a file base with attributes
        arc_proc = ArcProcessing(arc_opt, args.verbose, 'CHLA', file_att)
        arc_proc.climatology_path = climatology_path

        arc_proc.compute_chla_image(input_file, output_file, timeliness)
        return

    ##WORKING WITH DATES
    if start_date is None or end_date is None:
        start_date = options['start_date']
        end_date = options['end_date']
    date_run = start_date

    output_type = f'CHLA_{chla_algo}'

    arc_proc = ArcProcessing(arc_opt, args.verbose, output_type, None)
    arc_proc.climatology_path = climatology_path
    while date_run <= end_date:
        if args.verbose:
            print('*****************************')
            print(f'[INFO] Date: {date_run}')
        make_processing = True
        input_path = arc_opt.get_folder_date(options['input_path'], options['input_path_organization'], date_run, False)
        dateyj = date_run.strftime('%Y%j')
        name_rrs = f'C{dateyj}_rrs-arc-4km.nc'
        input_file = os.path.join(input_path, name_rrs)
        if not os.path.exists(input_file):
            print(f'[WARNING] Input file {input_file} for date {date_run} is not available. Skiping...')
            make_processing = False

        output_path = arc_opt.get_folder_date(options['output_path'], options['output_path_organization'], date_run,
                                              True)
        if output_path is None:
            print(f'[WARNING] Output path {input_path} for date {date_run} is not available. Skiping...')
            make_processing = False

        output_name = f'C{dateyj}_chl-arc-4km.nc'
        output_file = os.path.join(output_path, output_name)
        if os.path.exists(output_file) and not overwrite:
            print(f'[INFO] Output file {output_file} already exists. Skipping...')
            make_processing = False

        if make_processing:
            if args.verbose:
                print(f'[INFO] Input file: {input_file}')
                print(f'[INFO] Output file: {output_file}')
            arc_proc.compute_chla_image(input_file, output_file, timeliness)
        date_run = date_run + timedelta(hours=24)


def run_kd490(arc_opt, start_date, end_date):
    section = 'KD490'
    from datetime import timedelta
    options = arc_opt.get_basic_options(section)
    if options is None:
        return
    overwrite = arc_opt.get_value_param(section, 'overwrite', False, 'boolean')
    file_base = arc_opt.get_value_param(section, 'file_base', None, 'file')
    file_att = arc_opt.get_value_param(section, 'file_att', None, 'file')
    if file_base is None:
        print(f'[ERROR] file_base does not exist or is not available in the config file. Please review it')
        return

    timelinesses = ['NR', 'NT']
    for t in timelinesses:
        if file_base.split('/')[-1].find(t) > 0:
            timeliness = t
    if timeliness is None:
        print(f'[ERROR] Timeliness is not defined in the file base. It should be NR or NT')
        return

    from arc_processing import ArcProcessing
    if args.verbose:
        print('[INFO] KD400 PROCESSING OPTIONS:')
        for opt in options:
            print(f'[INFO]  {opt}->{options[opt]}')
        print(f'[INFO]  file_base:->{file_base}')
        print(f'[INFO]  file_att:->{file_att}')
        print(f'[INFO]  overwrite:->{overwrite}')
        print(f'[INFO]  timeliness:->{timeliness}')

    climatology_path = options['climatology_path']

    ##WORKING WITH SINGLE GRANULE
    input_name = arc_opt.get_value_param('KD490', 'name_input', None, 'str')
    if input_name is not None:
        input_file = os.path.join(options['input_path'], input_name)
        if os.path.exists(input_file):
            if args.verbose:
                print(f'[INFO] Working with the single file: {input_file}')
        else:
            print(f'[ERROR] File {input_file} does not exist')
        output_name = arc_opt.get_value_param('KD490', 'name_output', None, 'str')
        if output_name is None:
            output_name = 'SingleOutputkd490.nc'
        output_file = os.path.join(options['output_path'], output_name)
        if args.verbose:
            print(f'[INFO] Output file: {output_file}')
        # defining arc_proc, last parameters (file_at) is none because it's defined a file base with attributes
        arc_proc = ArcProcessing(arc_opt, args.verbose, 'KD490', file_att)
        arc_proc.climatology_path = climatology_path
        arc_proc.compute_kd490_image(input_file, output_file, timeliness)
        return

    ##WORKING WITH DATES
    if start_date is None or end_date is None:
        start_date = options['start_date']
        end_date = options['end_date']
    date_run = start_date

    arc_proc = ArcProcessing(arc_opt, args.verbose, 'KD490', file_att)
    arc_proc.climatology_path = climatology_path
    while date_run <= end_date:
        if args.verbose:
            print('*****************************')
            print(f'[INFO] Date: {date_run}')
        make_processing = True
        input_path = arc_opt.get_folder_date(options['input_path'], options['input_path_organization'], date_run, False)
        dateyj = date_run.strftime('%Y%j')
        name_rrs = f'C{dateyj}_rrs-arc-4km.nc'
        input_file = os.path.join(input_path, name_rrs)
        if not os.path.exists(input_file):
            print(f'[WARNING] Input file {input_file} for date {date_run} is not available. Skiping...')
            make_processing = False

        output_path = arc_opt.get_folder_date(options['output_path'], options['output_path_organization'], date_run,
                                              True)
        if output_path is None:
            print(f'[WARNING] Output path {input_path} for date {date_run} is not available. Skiping...')
            make_processing = False

        output_name = f'C{dateyj}_kd490-arc-4km.nc'
        output_file = os.path.join(output_path, output_name)
        if os.path.exists(output_file) and not overwrite:
            print(f'[INFO] Output file {output_file} already exists. Skipping...')
            make_processing = False

        if make_processing:
            if args.verbose:
                print(f'[INFO] Input file: {input_file}')
                print(f'[INFO] Output file: {output_file}')
            arc_proc.compute_kd490_image(input_file, output_file, timeliness)
        date_run = date_run + timedelta(hours=24)


# mode: CHLA or TRANSP
def run_month(arc_opt, mode, start_date, end_date):
    options = arc_opt.get_processing_options()
    if mode is None:
        output_type = arc_opt.get_value_param('PROCESSING', 'output_type', 'CHLA', 'str')
    else:
        output_type = mode

    from arc_processing import ArcProcessing
    from calendar import monthrange

    if args.verbose:
        print('[INFO] PROCESSING OPTIONS:')
        for opt in options:
            print(f'[INFO]  {opt}->{options[opt]}')

    file_base, timeliness = get_monthly_timeliness(arc_opt)

    if args.verbose:
        print(f'[INFO] File base: {file_base}')
        print(f'[INFO] Timeliness: {timeliness}')
    arc_proc = ArcProcessing(arc_opt, args.verbose, output_type, None)

    if start_date is None or end_date is None:
        start_date = options['start_date']
        end_date = options['end_date']
    start_date = start_date.replace(day=15)
    end_date = end_date.replace(day=15)
    date_run = start_date

    # from datetime import datetime as dt
    # date_run = dt(2023,10,11)
    # date_run.replace(day=15)

    if output_type == 'CHLA':
        file_date = 'CDATE_chl-arc-4km.nc'
        param_name = 'plankton'
    if output_type == 'TRANSP':
        file_date = 'CDATE_kd490-arc-4km.nc'
        param_name = 'transp'

    file_date_format = '%Y%j'

    if output_type == 'RRS510':
        file_date = 'DATE_c3s_obs-oc_glo_bgc-reflectance_my_l3-multi-4km_P1D.nc'
        param_name = 'rrs510'
        file_date_format = '%Y%m%d'

    while date_run <= end_date:
        if args.verbose:
            print('*****************************')
            print(f'[INFO] Date: {date_run.month}/{date_run.year}')

        output_path = arc_opt.get_folder_year(options['output_path'], date_run, True)
        if output_path is None:
            print(f'[WARNING] Output path {output_path} for date {date_run} is not available. Skiping...')
            date_run = date_run + timedelta(days=30)
            date_run = date_run.replace(day=15)
            continue

        nfiles_month = monthrange(date_run.year, date_run.month)[1]

        input_files, input_files_timeliness = arc_opt.get_list_files_month(options['input_path'],
                                                                           options['input_path_organization'],
                                                                           date_run.year,
                                                                           date_run.month, file_date, file_date_format,
                                                                           timeliness)

        nfiles_available = len(input_files)
        ntimeliness = len(input_files_timeliness)
        if nfiles_available == 0:
            print(
                f'[ERROR] No files avaiable for computing the {output_type} average for {date_run.year}/{date_run.month}. Skipping...')
            date_run = date_run + timedelta(days=30)
            date_run = date_run.replace(day=15)
            continue
        if nfiles_available < nfiles_month:
            print(f'[WARNING] Only {nfiles_available} of {nfiles_month} are available for computing the average')
        if ntimeliness < nfiles_available:
            print(f'[WARNING] Only {ntimeliness} of {nfiles_available} are in the correct timeliness {timeliness}')

        sdate = date_run.replace(day=1).strftime('%Y%j')
        edate = date_run.replace(day=nfiles_month).strftime('%j')

        name_out_end = f'C{sdate}{edate}-{param_name}_monthly-arc-4km.nc'
        file_out = os.path.join(output_path, name_out_end)
        if args.verbose:
            print(f'[INFO] Output file: {file_out}')

        name_source = f'C{sdate}{edate}-{param_name}_monthly-arc-4km.sources'
        file_source = os.path.join(output_path, name_source)
        name_source_timeliness = f'C{sdate}{edate}-{param_name}_monthly-arc-4km_{timeliness}.sources'
        file_source_timeliness = os.path.join(output_path, name_source_timeliness)

        compute_month = True
        # CHECK IF FILE MUST BE RUN AGAIN
        if os.path.exists(file_out) and os.path.exists(file_source_timeliness):
            nsources_prev = 0
            f1 = open(file_source_timeliness, 'r')
            for line in f1:
                if len(line.strip()) > 0:
                    nsources_prev = nsources_prev + 1
            f1.close()
            if nsources_prev == ntimeliness:
                compute_month = False
                print(f'[INFO] Month file has already been computed. Skipping...')

        if compute_month:
            arc_proc.create_source_list(input_files, file_source)
            arc_proc.create_source_list(input_files_timeliness, file_source_timeliness)
            arc_proc.compute_month(date_run.year, date_run.month, timeliness, input_files, file_out)

        date_run = date_run + timedelta(days=30)
        date_run = date_run.replace(day=15)


def get_monthly_timeliness(arc_opt):
    file_base = arc_opt.get_value_param('PROCESSING', 'file_base', None, 'str')
    timeliness = arc_opt.get_value_param('PROCESSING', 'timeliness', None, 'str')
    if file_base is not None:
        if os.path.exists(file_base):
            if timeliness is None:
                if file_base.find('_NR_') > 0:
                    timeliness = 'NR'
                if file_base.find('_NT_') > 0:
                    timeliness = 'NT'
    return file_base, timeliness


def check_py():
    print('[INFO] Checking py imports...')
    check = True
    try:
        import simplekml
    except:
        print('[ERROR] simplekml was not found...')
        check = False
    try:
        import netCDF4
    except:
        print('[ERROR] netCDF4 was not found...')
        check = False
    try:
        import numpy as np
    except:
        print('[ERROR] numpy was not found')
        check = False
    try:
        import zipfile as zp
    except:
        print('[ERROR] zipfile was not found')
        check = False
    try:
        import json
    except:
        print('[ERROR] json was not found')
        check = False
    try:
        import pyresample
    except:
        print('[ERROR] pyresample was not found')
        check = False
    try:
        import pyproj
    except:
        print('[ERROR] pyproj was not found')
        check = False
    try:
        import shapely
    except:
        print('[ERROR] shapely was not found')
        check = False
    try:
        import multiprocessing
        print(f'[INFO] Multiprocessing is allowed with {os.cpu_count()} CPUs')
    except:
        print('[ERROR] multiprocessing was not found')
        check = False
    if not check:
        print('[ERROR] Some packages are not available...')
    else:
        print('[INFO] All the packages are available')


def get_dates_from_arg():
    from datetime import datetime as dt
    from datetime import timedelta
    start_date = None
    end_date = None
    if args.start_date:
        try:
            start_date = dt.strptime(args.start_date, '%Y-%m-%d')
        except:
            try:
                tdelta = int(args.start_date)
                start_date = dt.now() + timedelta(days=tdelta)
                start_date = start_date.replace(hour=12, minute=0, second=0, microsecond=0)
            except:
                print(f'[ERROR] Start date {args.start_date} is not in the correct format: YYYY-mm-dd or integer')
    if args.end_date:
        try:
            end_date = dt.strptime(args.end_date, '%Y-%m-%d')
        except:
            try:
                tdelta = int(args.end_date)
                end_date = dt.now() + timedelta(days=tdelta)
                end_date = end_date.replace(hour=12, minute=0, second=0, microsecond=0)
            except:
                print(f'[ERROR] End date {args.end_date} is not in the correct format: YYYY-mm-dd or integer')
    if args.start_date and not args.end_date:
        end_date = start_date

    return start_date, end_date


if __name__ == '__main__':
    main()
