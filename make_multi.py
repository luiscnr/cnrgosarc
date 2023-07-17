import argparse
import os
from datetime import timedelta
from arc_multi_sources import ARC_MULTI_SOURCES
import warnings

warnings.filterwarnings(action='ignore', category=ResourceWarning)

parser = argparse.ArgumentParser(description="Artic resampler")
parser.add_argument("-m", "--mode", help="Mode",
                    choices=["CHECKPY", "CHECK", "GRID", "RESAMPLE", "RESAMPLEFILE", "INTEGRATE", "CHLA", "KD490", "QL",
                             "MONTHLY_CHLA", "MONTHLY_KD490"],
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
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
args = parser.parse_args()


def only_test():
    path = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/S3A_OL_2_WFR____20220715T232353_20220715T232653_20220717T115842_0179_087_301_1800_MAR_O_NT_003.SEN3'
    from netCDF4 import Dataset
    import numpy as np
    from kd_algorithm import KD_ALGORITHMS
    path_kd = os.path.join(path, 'trsp.nc')
    dataset_kd = Dataset(path_kd)
    kd_array_olci = np.array(dataset_kd.variables['KD490_M07'][2100:2800, 4000:4700])
    dataset_kd.close()

    path_490 = os.path.join(path, 'Oa04_reflectance.nc')
    dataset_490 = Dataset(path_490)
    array_490 = np.array(dataset_490.variables['Oa04_reflectance'][2100:2800, 4000:4700])

    path_560 = os.path.join(path, 'Oa06_reflectance.nc')
    dataset_560 = Dataset(path_560)
    array_560 = np.array(dataset_560.variables['Oa06_reflectance'][2100:2800, 4000:4700])

    path_chla = os.path.join(path, 'chl_oc4me.nc')
    dataset_chla = Dataset(path_chla)
    array_chla = np.array(dataset_chla.variables['CHL_OC4ME'][2100:2800, 4000:4700])


    kda = KD_ALGORITHMS('OK2-560')
    kda.qratio = 0.97941
    kd_array_new = kda.compute_kd490_ok2_560(array_490, array_560)


    kd_array_olci[kd_array_olci != 255] = np.power(10, kd_array_olci[kd_array_olci != 255])
    array_chla[array_chla != 255] = np.power(10,array_chla[array_chla != 255])


    kd_array_new[kd_array_olci == 255] = -999.0
    kd_array_new[kd_array_olci == 10.0] = -999.0
    kd_array_new[array_chla == 255] = -999.0

    file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/MULTI/kdcomparison.csv'
    f1 = open(file_out,'w')
    f1.write('Y;X;OLCI;COMPUTED;CHLA;RATIO')
    for y in range(kd_array_olci.shape[0]):
        for x in range(kd_array_olci.shape[1]):
            ypoint = y + 2100
            xpoint = x + 4000
            if kd_array_new[y,x]!=-999.0:
                val_old = kd_array_olci[y, x]
                val_new = kd_array_new[y, x]
                val_chla = array_chla[y,x]
                ron = kd_array_olci[y, x]/kd_array_new[y, x]
                if y==691 and x==4:
                    print(f'OLCI: {val_old:.4f};COMPUTED: {val_new:.4f}; CHLA: {val_chla:.4f}; RON: {ron:.4f}')
                line = f'{ypoint};{xpoint};{val_old};{val_new};{val_chla};{ron}'
                f1.write('\n')
                f1.write(line)
    f1.close()

    dataset_490.close()
    dataset_560.close()

    return True


def main():
    # if only_test():
    #     return
    print('[INFO] Started Artic Processing Tool [MULTI 4 KM]')
    if args.mode == "CHECKPY":
        check_py()
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
            print(f'[ERROR] Output file: {output_file} is not valid')
            return
        # file_at = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/CONFIG_FILES/global_attributes.ini'
        file_at = args.config_file
        if not os.path.isfile(file_at):
            print(f'[ERROR] Attributes files: {file_at} is not valid')
        file_base = args.base_file
        if not os.path.isfile(file_base):
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

    ##FROM HERE, ALL THE MODES REQUIRE CONFIGURATION MODEL. DATES COULD BE ALSO PASSED AS ARGS
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
    # dbase = Dataset(base_file)
    # dname = dbase.title
    # dbase.close()

    ams = ARC_MULTI_SOURCES(options['input_path'], options['input_path_organization'], options['moi_user'],
                            options['moi_pass'], args.verbose)
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
        print('[INFO] CHLA PROCESSING OPTIONS:')
        for opt in options:
            print(f'[INFO]  {opt}->{options[opt]}')
        print(f'[INFO]  file_base:->{file_base}')
        print(f'[INFO]  file_att:->{file_att}')
        print(f'[INFO]  overwrite:->{overwrite}')
        print(f'[INFO]  timeliness:->{timeliness}')

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

        arc_proc.compute_chla_image(input_file, output_file, timeliness)
        return

    ##WORKING WITH DATES
    if start_date is None or end_date is None:
        start_date = options['start_date']
        end_date = options['end_date']
    date_run = start_date

    arc_proc = ArcProcessing(arc_opt, args.verbose, 'CHLA', None)
    while date_run <= end_date:
        if args.verbose:
            print('*****************************')
            print(f'[INFO] Date: {date_run}')
        make_processing = True
        input_path = arc_opt.get_folder_date(options['input_path'], options['input_path_organization'], date_run, False)
        dateyj = date_run.strftime('%Y%j')
        # dateymd = date_run.strftime('%Y%m%d')
        name_rrs = f'C{dateyj}_rrs-arc-4km.nc'
        # name_rrs = f'{dateymd}_cmems_obs-oc_arc_bgc-reflectance_my_l3-multi-4km_P1D.nc'
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

    ##WORKING WITH SINGLE GRANULE
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
            output_name = 'SingleOutputkd490.nc'
        output_file = os.path.join(options['output_path'], output_name)
        if args.verbose:
            print(f'[INFO] Output file: {output_file}')
        # defining arc_proc, last parameters (file_at) is none because it's defined a file base with attributes
        arc_proc = ArcProcessing(arc_opt, args.verbose, 'KD490', file_att)
        arc_proc.compute_kd490_image(input_file, output_file, timeliness)
        return

    ##WORKING WITH DATES
    if start_date is None or end_date is None:
        start_date = options['start_date']
        end_date = options['end_date']
    date_run = start_date

    arc_proc = ArcProcessing(arc_opt, args.verbose, 'KD490', file_att)
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
            arc_proc.compute_kd490_image(input_file,output_file,timeliness)
        date_run = date_run + timedelta(hours=24)


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
