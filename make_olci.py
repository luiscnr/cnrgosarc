
import __init__,math,os,argparse,configparser
from datetime import timedelta
# import netCDF4
# import numpy as np
# import zipfile as zp
# from arc_mapinfo import ArcMapInfo
# from arc_integration import ArcIntegration
# from olci_l2 import OLCI_L2
# import simplekml
import numpy as np
from arc_mapinfo import ArcMapInfo
from olci_l2 import OLCI_L2
from arc_processing import ArcProcessing
from arc_gpr_model import ARC_GPR_MODEL

parser = argparse.ArgumentParser(description="Arctic OLCI processor")
parser.add_argument("-m", "--mode", help="Mode",
                    choices=["CHECKPY", "CHECKMODEL", "GRID", "RESAMPLE", "RESAMPLEPML", "INTEGRATE", "CHLA", "QL",
                             "MONTHLY_CHLA", "MONTHLY_KD490"],
                    required=True)
parser.add_argument("-p", "--product", help="Input product (testing)")
parser.add_argument('-i', "--inputpath", help="Input directory")
parser.add_argument('-o', "--outputpath", help="Output file")
parser.add_argument('-tp', "--temp_path", help="Temporary directory")
parser.add_argument('-sd', "--start_date", help="Start date (yyyy-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date (yyyy-mm-dd")
parser.add_argument('-c', "--config_file", help="Configuration file (Default: arc_config.ini)")
parser.add_argument('-bf', "--base_file", help="Create base file for the following mode: RESAMPLE,INTEGRATE.",
                    action="store_true")
parser.add_argument('-var', "--variable_plot", help="Variable to be plot using QL mode. Default: CHL",
                    choices=["CHL", "KD490"])
parser.add_argument("-chla_algo", "--chla_algorithm",help="Chl-a algorithm",choices=["SeaSARC","CIAO"],default="CIAO")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
args = parser.parse_args()


def main():
    print('[INFO] Started Artic Processing Tool')

    if args.mode == "CHECKPY":
        check_py()
        return

    if args.mode == "CHECKMODEL":
        check_model()
        return

    if args.mode == 'GRID' and args.outputpath:  ##creating single file grid
        if args.verbose:
            print('Create grid file')
        ami = ArcMapInfo(None, args.verbose)
        file_out = args.outputpath
        ami.create_nc_filegrid(file_out, True, True)
        return

    if args.mode == "RESAMPLEPML" and args.product and args.outputpath:  ##testing, resampling of a single PML file
        ami = ArcMapInfo(None, args.verbose)
        file_out = args.outputpath
        fpml = args.product
        ami.make_resample_pml(fpml, file_out)
        return

    if args.mode == 'QL' and args.product and args.outputpath:  ##Quick Look Generation
        ami = ArcMapInfo(None, args.verbose)
        file_out = args.outputpath
        fdataset = args.product
        variable = 'CHL'
        if args.variable_plot:
            variable = args.variable_plot
        ami.save_full_fdata(file_out, fdataset, variable)
        return

    if args.mode == 'INTEGRATE' and args.base_file and args.outputpath:
        output_path = args.outputpath
        if not os.path.exists(output_path) and not os.path.isdir(output_path):
            print(f'[ERROR] Output path {output_path} does not exist or nor it is a directory')
            return
        if args.verbose:
            print('[INFO] Started creation of base files for integration...')
        from arc_integration import ArcIntegration
        file_at = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/CONFIG_FILES/global_attributes.ini'

        arcInt = ArcIntegration(None, args.verbose, None, 'TRANSP', file_at)
        arcInt.ami.ifile_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/GRID_FILES/ArcGrid_65_90_300m_GridBase.nc'
        fout = os.path.join(output_path, 'ArcGrid_69_90_300m_TRANSP_NR_Base.nc')
        arcInt.create_nc_file_out(fout, 'NR')
        fout = os.path.join(output_path, 'ArcGrid_69_90_300m_TRANSP_NT_Base.nc')
        arcInt.create_nc_file_out(fout, 'NT')

        return

    if args.mode == 'CHLA' and args.base_file and args.outputpath:
        output_path = args.outputpath
        if not os.path.exists(output_path) and not os.path.isdir(output_path):
            print(f'[ERROR] Output path {output_path} does not exist or nor it is a directory')
            return
        if args.verbose:
            print(f'[INFO] Started creation of base files for processing CHLA')
        file_at = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/CONFIG_FILES/global_attributes.ini'
        from arc_processing import ArcProcessing
        arcProc = ArcProcessing(None, args.verbose, 'CHLA', file_at)
        file_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/GRID_FILES/ArcGrid_65_90_300m_GridBase.nc'

        fout = os.path.join(output_path, 'ArcGrid_65_90_300m_PLANKTON_NR_Base.nc')
        arcProc.create_nc_file_out(fout, file_base, 'NR')
        fout = os.path.join(output_path, 'ArcGrid_65_90_300m_PLANKTON_NT_Base.nc')
        arcProc.create_nc_file_out(fout, file_base, 'NT')
        return

    if args.mode == 'MONTHLY_CHLA' and args.base_file and args.outputpath:
        output_path = args.outputpath
        if not os.path.exists(output_path) and not os.path.isdir(output_path):
            print(f'[ERROR] Output path {output_path} does not exist or nor it is a directory')
            return
        if args.verbose:
            print(f'[INFO] Started creation of base files for processing monthly CHLA')
        file_at = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/CONFIG_FILES/global_attributes.ini'
        from arc_processing import ArcProcessing
        arcProc = ArcProcessing(None, args.verbose, 'CHLA', file_at)
        file_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/GRID_FILES/ArcGrid_65_90_300m_GridBase_NOSENSORMASK.nc'

        fout = os.path.join(output_path, 'ArcGrid_65_90_300m_PLANKTON_NR_MONTHLY_Base.nc')
        arcProc.create_nc_file_out_month(fout, file_base, 'NR')
        fout = os.path.join(output_path, 'ArcGrid_65_90_300m_PLANKTON_NT_MONTHLY_Base.nc')
        arcProc.create_nc_file_out_month(fout, file_base, 'NT')
        return

    if args.mode == 'MONTHLY_KD490' and args.base_file and args.outputpath:
        output_path = args.outputpath
        if not os.path.exists(output_path) and not os.path.isdir(output_path):
            print(f'[ERROR] Output path {output_path} does not exist or nor it is a directory')
            return
        if args.verbose:
            print(f'[INFO] Started creation of base files for processing monthly KD490')
        file_at = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/CONFIG_FILES/global_attributes.ini'
        from arc_processing import ArcProcessing
        arcProc = ArcProcessing(None, args.verbose, 'TRANSP', file_at)
        file_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/GRID_FILES/ArcGrid_65_90_300m_GridBase_NOSENSORMASK.nc'

        fout = os.path.join(output_path, 'ArcGrid_65_90_300m_TRANSP_NR_MONTHLY_Base.nc')
        arcProc.create_nc_file_out_month(fout, file_base, 'NR')
        # fout = os.path.join(output_path, 'ArcGrid_65_90_300m_PLANKTON_NT_MONTHLY_Base.nc')
        # arcProc.create_nc_file_out_month(fout, file_base, 'NT')
        return

    ##FROM HERE, ALL THE MODES REQUIRE CONFIGURATION MODEL. DATES COULD BE ALSO PASSED AS ARGS
    if not args.config_file:
        print(f'[ERROR] Config file or input product should be defined for {args.mode} option. Exiting...')
        return
    if not os.path.exists(args.config_file):
        print(f'[ERROR] Config file {args.config_file} does not exist. Exiting...')
        return
    try:
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
            ##error definint start or end dates
            return

    if args.mode == 'RESAMPLE' and args.base_file and args.product:  ##options are required
        folci = args.product
        olimage = OLCI_L2(folci, args.verbose)
        olimage.get_geo_and_params()
        ami = ArcMapInfo(None, args.verbose)
        file_out = args.outputpath
        ami.create_nc_file_resample_base(olimage, file_out, arc_opt)

    if args.mode == "RESAMPLE" and args.product:  ##resampling of a single granule. Options are also required
        ami = ArcMapInfo(None, args.verbose)
        file_out = args.outputpath
        folci = args.product
        olimage = OLCI_L2(folci, args.verbose)
        olimage.get_geo_and_params()
        line_output = ami.make_resample_impl(olimage, file_out, 1, -1, arc_opt, None)
        if args.verbose:
            print(f'[INFO] Output line: {line_output}')
        return

    if args.mode == 'RESAMPLE':
        run_resample(arc_opt, start_date, end_date)
        return

    if args.mode == 'INTEGRATE':
        run_integration(arc_opt, start_date, end_date)
        return

    if args.mode == 'CHLA':
        run_chla(arc_opt, start_date, end_date)
        return

    if args.mode == 'MONTHLY_CHLA':
        operative_mode = arc_opt.get_value_param('PROCESSING', 'operative_mode', False, 'boolean')
        if operative_mode:
            date = check_monthly_operative_mode(arc_opt)
            if date is not None:
                run_month(arc_opt, 'CHLA', date, date)
        else:
            run_month(arc_opt, 'CHLA', start_date, end_date)
        return

    if args.mode == 'MONTHLY_KD490':
        operative_mode = arc_opt.get_value_param('PROCESSING', 'operative_mode', False, 'boolean')
        if operative_mode:
            date = check_monthly_operative_mode(arc_opt)
            if date is not None:
                run_month(arc_opt, 'TRANSP', date, date)
        else:
            run_month(arc_opt, 'TRANSP', start_date, end_date)
        return

    if args.mode == 'QL':
        run_ql(arc_opt, start_date, end_date)


def check_monthly_operative_mode(arc_opt):
    from datetime import datetime as dt
    file_base, timeliness = get_monthly_timeliness(arc_opt)
    day_today = dt.utcnow().day
    if timeliness == 'NR' and day_today >= 8:
        print('[WARNING] Month operative NR files are only processed between days 1 and 8 of the month')
        return None
    if timeliness == 'NR' and day_today >= 8:
        print('[WARNING] Month operative NT files are only processed after day 8 of the month')
        return None
    month_today = dt.utcnow().month
    year_today = dt.utcnow().year
    if month_today == 1:
        month_processing = 12
        year_processing = year_today - 1
    else:
        month_processing = month_today - 1
        year_processing = year_today
    date = dt(year_processing, month_processing, 15)
    return date


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
    # date_run = dt(2022,10,11)
    # date_run.replace(day=15)

    if output_type == 'CHLA':
        file_date = 'ODATE_plankton-arc-fr.nc'
        param_name = 'plankton'
    if output_type == 'TRANSP':
        file_date = 'ODATE_transp-arc-fr.nc'
        param_name = 'transp'

    file_date_format = '%Y%j'

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

        name_out_end = f'O{sdate}{edate}-{param_name}_monthly-arc-fr.nc'
        file_out = os.path.join(output_path, name_out_end)
        if args.verbose:
            print(f'[INFO] Output file: {file_out}')

        name_source = f'O{sdate}{edate}-{param_name}_monthly-arc-fr.sources'
        file_source = os.path.join(output_path, name_source)
        name_source_timeliness = f'O{sdate}{edate}-{param_name}_monthly-arc-fr_{timeliness}.sources'
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


def compute_month_chl(arc_opt):
    options = arc_opt.get_processing_options()
    if options is None:
        return
    ##ONLY CHLA MODE IS IMPLEMENTED
    output_type = arc_opt.get_value_param('PROCESSING', 'output_type', 'CHLA', 'str')
    if not output_type == 'CHLA':
        return
    from arc_processing import ArcProcessing
    if args.verbose:
        print('[INFO] PROCESSING OPTIONS:')
        for opt in options:
            print(f'[INFO]  {opt}->{options[opt]}')

    file_base = arc_opt.get_value_param('PROCESSING', 'file_base', None, 'str')
    timeliness = arc_opt.get_value_param('PROCESSING', 'timeliness', None, 'str')
    if file_base is not None:
        if os.path.exists(file_base):
            if timeliness is None:
                if file_base.find('NR') > 0:
                    timeliness = 'NR'
                if file_base.find('NT') > 0:
                    timeliness = 'NT'
    if args.verbose:
        print(f'[INFO] File base: {file_base}')
        print(f'[INFO] Timeliness: {timeliness}')
    arc_proc = ArcProcessing(arc_opt, args.verbose, output_type, None)
    fileout = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/INTEGRATED/O201907_plankton-arc-fr.nc'
    arc_proc.compute_chla_month(fileout, timeliness)




def adding_time():
    from datetime import datetime as dt
    from datetime import timedelta
    dir_base = '/store/COP2-OC-TAC/arc/integrated'
    date_here = dt(2019, 6, 10)
    date_end = dt(2019, 6, 23)
    while date_here <= date_end:
        yyyy = date_here.strftime('%Y')
        jjj = date_here.strftime('%j')
        file_in = os.path.join(dir_base, yyyy, jjj, f'O{yyyy}{jjj}_rrs-arc-fr_NOTIME.nc')
        file_out = os.path.join(dir_base, yyyy, jjj, f'O{yyyy}{jjj}_rrs-arc-fr.nc')
        print(file_in, file_out)
        # os.rename(file_in,file_out)
        copy_nc_adding_time_variable(file_in, file_out)
        date_here = date_here + timedelta(hours=24)

    # file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/175/O2019175_rrs-arc-fr_NOTIME.nc'
    # file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/175/O2019175_rrs-arc-fr.nc'
    # copy_nc_adding_time_variable(file_in,file_out)
    # file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/175/O2019175_plankton-arc-fr_NOTIME.nc'
    # file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/175/O2019175_plankton-arc-fr2.nc'
    # copy_nc_adding_time_variable(file_in, file_out)
    #
    # file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/207/O2019207_rrs-arc-fr_NOTIME.nc'
    # file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/207/O2019207_rrs-arc-fr.nc'
    # copy_nc_adding_time_variable(file_in, file_out)
    # file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/207/O2019207_plankton-arc-fr_NOTIME.nc'
    # file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/207/O2019207_plankton-arc-fr.nc'
    # copy_nc_adding_time_variable(file_in, file_out)

    # file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2016/183/O2016183_transp-arc-fr_NOTIME.nc'
    # file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2016/183/O2016183_transp-arc-fr.nc'
    # copy_nc_adding_time_variable(file_in, file_out)
    # file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2016/197/O2016197_transp-arc-fr_NOTIME.nc'
    # file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2016/197/O2016197_transp-arc-fr.nc'
    # copy_nc_adding_time_variable(file_in, file_out)

    # file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/O201906_plankton-arc-fr_NOTIME.nc'
    # file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/O201906_plankton-arc-fr_NOSENSOR.nc'
    # copy_nc_excluding_variables(file_in,file_out,['SENSORMASK'])
    # file_out_end = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/O201906_plankton-arc-fr.nc'
    # copy_nc_adding_time_variable(file_out, file_out_end)

    # file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/O201907_plankton-arc-fr_NOTIME.nc'
    # file_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/O201907_plankton-arc-fr_NOSENSOR.nc'
    # copy_nc_excluding_variables(file_in,file_out,['SENSORMASK'])
    # file_out_end = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/O201907_plankton-arc-fr.nc'
    # copy_nc_adding_time_variable(file_out, file_out_end)

    ##TO SET ATRIBUTES
    # file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/integrated/2019/175/O2019175_plankton-arc-fr_NOTIME.nc'
    # from netCDF4 import Dataset
    # ncsat = Dataset(file_in, 'a')
    # ncsat.cmems_product_id = "OCEANCOLOUR_ARC_BGC_L3_MY_009_123"
    # ncsat.title = "cmems_obs-oc_arc_bgc-plankton_my_l3-olci-300m_P1D"
    # ncsat.timeliness = "NT"

    # import configparser
    # options = configparser.ConfigParser()
    # file_at = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/CONFIG_FILES/global_attributes.ini'
    # options.read(file_at)
    # ats = list(ncsat.ncattrs())
    # at_dict = dict(options['GLOBAL_ATTRIBUTES'])
    # at_dict['parameter'] = 'Chlorophyll-a concentration'
    # at_dict['parameter_code'] = 'PLANKTON'
    # at_dict['product_level'] = 'L4'
    # from datetime import datetime as dt
    # date = dt(2019,7,1)
    # date_fin = dt(2019,7,31)
    # names = []
    # while date<=date_fin:
    #     print(date)
    #     date_str = date.strftime('%Y%j')
    #     name = f'O{date_str}_plankton-arc-fr.nc'
    #     names.append(name)
    #     date  = date+timedelta(hours=24)
    # namestr = ' '.join(names)
    # at_dict['source_files'] = namestr
    # for at in at_dict:
    #     if at not in ats:
    #         print(at, at_dict[at])
    #         ncsat.setncattr(at, at_dict[at])
    # ncsat.close()


def modify_chunksizes():
    file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/INTEGRATED/2016/183/O2016183_transp-arc-fr.nc'
    file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/INTEGRATED/2019/O2019152181-plankton_monthly-arc-fr.nc'
    chunk_sizes = [1223, 3669]  # [1223, 3669, 6115, 18345]
    from datetime import datetime as dt

    for idx in range(len(chunk_sizes)):
        chunk_size = chunk_sizes[idx]
        # iday = idx +1
        # date_here = dt(2023,2,iday)
        imonth = idx + 1
        date_here = dt(2023, imonth, 1)
        if imonth == 1:
            date_here_end = dt(2023, imonth, 31)
        if imonth == 2:
            date_here_end = dt(2023, imonth, 28)

        year_str = date_here.strftime('%Y')
        jjj_str = date_here.strftime(('%j'))
        jjj_end = date_here_end.strftime(('%j'))
        # name_file = f'O{year_str}{jjj_str}_transp-arc-fr.nc'
        # file_out = os.path.join('/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/INTEGRATED/',year_str,jjj_str,name_file)
        name_file = f'O{year_str}{jjj_str}{jjj_end}-transp_monthly-arc-fr.nc'
        file_out = os.path.join('/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/INTEGRATED/', year_str, name_file)
        print(file_out, chunk_size)
        copy_nc_with_chunksize(file_in, file_out, chunk_size, date_here, date_here_end)


def copy_nc_with_chunksize(ifile, ofile, chunk_size, date_here, date_here_end):
    # csizes = (1,chunk_size,chunk_size)
    csizes = None
    from netCDF4 import Dataset
    from datetime import datetime as dt
    with Dataset(ifile) as src:
        dst = Dataset(ofile, 'w', format='NETCDF4')
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # modifying
        dst.start_date = date_here.strftime('%Y-%m-%d')
        if date_here_end is None:
            date_here_end = date_here
        dst.stop_date = date_here_end.strftime('%Y-%m-%d')
        dst.timeliness = 'NR'
        cdate = dt.utcnow()
        dst.creation_date = cdate.strftime('%Y-%m-%d')
        dst.creation_time = cdate.strftime('%H:%M:%S UTC')
        dst.cmems_product_id = 'OCEANCOLOUR_ARC_BGC_L4_NRT_009_122'
        dst.title = 'cmems_obs-oc_arc_bgc-transp_nrt_l4-olci-300m_P1M'
        timeseconds = (date_here - dt(1981, 1, 1, 0, 0, 0)).total_seconds()
        # copy dimensions
        for name, dimension in src.dimensions.items():
            if args.verbose:
                print(f'[INFO] -> Copying dimension: {name}')
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        # copy var
        for name, variable in src.variables.items():
            if args.verbose:
                print(f'[INFO] -> Copying variables: {name}')
            if name == 'time' or name == 'y' or name == 'x' or name == 'lat' or name == 'lon' or name == 'stereographic':
                dst.createVariable(name, variable.datatype, variable.dimensions, fill_value=-999, zlib=True,
                                   shuffle=True,
                                   complevel=6)
            else:
                dst.createVariable(name, variable.datatype, variable.dimensions, fill_value=-999, zlib=True,
                                   shuffle=True,
                                   complevel=6, chunksizes=csizes)
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)
            # copy variable data
            dst[name][:] = src[name][:]
        # modifying time
        dst.variables['time'][0] = [np.int32(timeseconds)]

    if args.verbose:
        print('COMPLETED')


def correcting_time_variable_in_plankton_files(start_date, end_date):
    dir_base = '/store/COP2-OC-TAC/arc/integrated'
    dir_output = '/store/COP2-OC-TAC/arc/daily'
    import stat
    if not os.path.exists(dir_output):
        os.mkdir(dir_output)
        # os.chmod(dir_output,0o666)
    run_date = start_date
    while run_date <= end_date:
        print(f'DATE-------------------------------------------------> {run_date}')
        yearstr = run_date.strftime('%Y')
        jjjstr = run_date.strftime('%j')
        name_file = f'O{yearstr}{jjjstr}_plankton-arc-fr.nc'
        input_path = os.path.join(dir_base, yearstr, jjjstr, name_file)
        if os.path.exists(input_path):
            output_year = os.path.join(dir_output, yearstr)
            if not os.path.exists(output_year):
                os.mkdir(output_year)
                # os.chmod(output_year,0o666)
            output_jday = os.path.join(output_year, jjjstr)
            if not os.path.exists(output_jday):
                os.mkdir(output_jday)
                # os.chmod(output_jday,0o666)
            output_path = os.path.join(output_jday, name_file)
            copy_nc_setting_time_variable(input_path, output_path)

        run_date = run_date + timedelta(hours=24)


def copy_nc_setting_time_variable(ifile_base, ofile):
    dst = None
    if not os.path.exists(ifile_base):
        print(f'[ERROR] File base: {ifile_base} does not exist')
        return dst
    if args.verbose:
        print(f'[INFO] Copying file grid: {ifile_base}...')

    # shutil.copy(self.ifile_base,ofile)

    cmd = f'cp -a {ifile_base} {ofile}'
    if args.verbose:
        print(f'[INFO] cmd: {cmd}')
    import subprocess
    import time
    subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    originalSize = os.path.getsize(ifile_base)
    historicalSize = -1
    while historicalSize != originalSize:
        if os.path.exists(ofile):
            historicalSize = os.path.getsize(ofile)
            if args.verbose:
                porc = (historicalSize / originalSize) * 100
                print(f'[INFO] Copying {porc:.2f} %')
        time.sleep(1)
    if args.verbose:
        print('[INFO] Copy completed')

    from netCDF4 import Dataset
    dst = Dataset(ofile, 'a', format='NETCDF4')
    from datetime import datetime as dt
    date_here = dt.strptime(dst.start_date, '%Y-%m-%d')
    print(f'[INFO] Setting time: {dst.start_date}')
    timeseconds = (date_here - dt(1981, 1, 1, 0, 0, 0)).total_seconds()
    var_time = dst.variables['time']
    var_time[0] = [np.int32(timeseconds)]

    dst.close()


def copy_nc_adding_time_variable(ifile, ofile):
    from netCDF4 import Dataset

    with Dataset(ifile) as src:
        dst = Dataset(ofile, 'w', format='NETCDF4')

        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)

        # copy dimensions
        for name, dimension in src.dimensions.items():
            if args.verbose:
                print(f'[INFO] -> Copying dimension: {name}')
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))

        # add time dimension
        dst.createDimension('time', 1)

        from datetime import datetime as dt
        date_here = dt.strptime(src.start_date, '%Y-%m-%d')
        print(date_here)
        timeseconds = (date_here - dt(1981, 1, 1, 0, 0, 0)).total_seconds()

        # add time variable
        var_time = dst.createVariable('time', 'i4', ('time',), fill_value=-999, zlib=True, complevel=6)
        var_time.long_name = "reference time"
        var_time.standard_name = "time"
        var_time.axis = "T"
        var_time.calendar = "Gregorian"
        var_time.units = "seconds since 1981-01-01 00:00:00"
        var_time[0] = [np.int32(timeseconds)]

        # for the remaining variables, creating as time,y,x
        for name, variable in src.variables.items():
            if args.verbose:
                print(f'[INFO] -> Copying variable: {name}')
            if name == 'KD490':
                print('Skipping...')
                continue
            if name == 'time' or name == 'y' or name == 'x' or name == 'lat' or name == 'lon' or name == 'stereographic':
                dst.createVariable(name, variable.datatype, variable.dimensions, fill_value=-999, zlib=True,
                                   shuffle=True,
                                   complevel=6)
            else:
                dst.createVariable(name, variable.datatype, ('time', 'y', 'x'), fill_value=-999, zlib=True,
                                   shuffle=True,
                                   complevel=6)
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)
            # copy data
            if name == 'time' or name == 'y' or name == 'x' or name == 'lat' or name == 'lon' or name == 'stereographic':
                dst[name][:] = src[name][:]
            else:
                # print('variable: ', name)
                dst[name][0, :, :] = src[name][:, :]

        dst.close()


def copy_nc_excluding_variables(ifile, ofile, excluded_variables):
    from netCDF4 import Dataset
    with Dataset(ifile) as src:
        dst = Dataset(ofile, 'w', format='NETCDF4')

        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)

        # copy dimensions
        for name, dimension in src.dimensions.items():
            if args.verbose:
                print(f'[INFO] -> Copying dimension: {name}')
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))

        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if name in excluded_variables:
                continue
            if args.verbose:
                print(f'[INFO] -> Copying variable: {name}')
            dst.createVariable(name, variable.datatype, variable.dimensions, fill_value=-999, zlib=True, shuffle=True,
                               complevel=6)
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)
            dst[name][:] = src[name][:]


def check_chla():
    file = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/INTEGRATED/2019/175/O2019175_AVERAGE-arc-fr.nc'
    from netCDF4 import Dataset
    ncsat = Dataset(file)
    wmask = np.array(ncsat.variables['sum_weights'])
    nvalid = np.count_nonzero(wmask[wmask > 0])
    print('N Valid mask: ', nvalid)
    bands = ['RRS400', 'RRS412_5', 'RRS442_5', 'RRS490', 'RRS510', 'RRS560', 'RRS620', 'RRS665', 'RRS673_75',
             'RRS681_25', 'RRS708_75']
    for band in bands:
        array = np.array(ncsat.variables[band])
        nvalid = np.count_nonzero(array[array > -999.0])
        print(band, nvalid)

    ncsat.close()

    return True


def check_model():
    file_model = None
    if args.chla_algorithm == 'CIAO':
        file_model = os.path.join(os.path.dirname(__init__.__file__), 'CIAO_Algorithm.json')
    elif args.chla_algorithm == 'SeaSARC':
        file_model = os.path.join(os.path.dirname(__init__.__file__), 'SeaSARC_Algorithm.json')

    if not os.path.exists(file_model):
        print(f'[ERROR] File model for algorithm: {args.chla_algorithm} is not available')
    else:
        print(f'[INFO] File model for algorithm {args.chla_algorithm} is avaiable. ')


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


def run_resampling_info():
    lines = ['Source;Width;Height;NTotal;NWater1;NWater2;NValid;PValid;NValidNew;PValidNew;NErrors']

    ##SINGLE IMAGE: LOCAL
    import zipfile as zp
    from olci_l2 import OLCI_L2
    import os
    # path = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/S3A_OL_2_WFR____20190624T004648_20190624T004948_20211129T080410_0180_046_145______MAR_R_NT_003.SEN3'
    path = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/S3B_OL_2_WFR____20190624T164944_20190624T165103_20210813T114041_0078_027_012______MAR_R_NT_003.SEN3'
    oimage = OLCI_L2(path, True)
    flag_mask, line = oimage.get_mask_default()
    lines.append(line)
    path_out = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/inforesampling_new.csv'
    f1 = open(path_out, 'w')
    for line in lines:
        f1.write(line)
        f1.write('\n')
    f1.close()
    ##VM(SINGLE PATH)
    # base_path = '/store/COP2-OC-TAC/arc/sources/20190624'
    # run_resampling_info_dir(base_path)


def run_resampling_info_dir(base_path):
    import zipfile as zp
    from olci_l2 import OLCI_L2
    import os
    lines = ['Source;Width;Height;NTotal;NWater1;NWater2;NValid;PValid;NValidNew;PValidNew;NErrors']
    path_out = os.path.join(base_path, 'INFO_RESAMPLING.csv')
    unzip_path = '/store/COP2-OC-TAC/arc/'
    for name in os.listdir(base_path):
        if not name.endswith('.zip'):
            continue
        prod_path = os.path.join(base_path, name)
        if zp.is_zipfile(prod_path):
            do_zip_here = True
            path_prod_u = prod_path.split('/')[-1][0:-4]
            path_prod_u = os.path.join(unzip_path, path_prod_u)
        if do_zip_here:
            with zp.ZipFile(prod_path, 'r') as zprod:
                if args.verbose:
                    print(f'[INFO] Unziping {name} to {unzip_path}')
                zprod.extractall(path=unzip_path)
        olimage = OLCI_L2(path_prod_u, args.verbose)
        flag_mask, line = olimage.get_mask_default()
        lines.append(line)
        if do_zip_here:
            # os.remove(prod_path)
            for fn in os.listdir(path_prod_u):
                os.remove(os.path.join(path_prod_u, fn))

    f1 = open(path_out, 'w')
    for line in lines:
        f1.write(line)
        f1.write('\n')
    f1.close()


def run_resample(arc_opt, start_date, end_date):
    options = arc_opt.get_resample_options()
    if options is None:
        print('[ERROR] Error getting the RESAMPLE options. Please review the config file')
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
    while date_ref <= date_fin:
        print(f'[INFO]******************************************************************************->{date_ref}')
        input_dir = arc_opt.get_folder_date_o(options, 'input_path', 'input_path_organization', date_ref, False)
        print(f'[INFO] Input dir: {input_dir}')
        if not os.path.exists(input_dir):
            print(f'[WARNING] Input directory {input_dir} is not available. Skiping...')
            date_ref = date_ref + timedelta(hours=24)
            continue
        output_dir_day = arc_opt.get_folder_date_o(options, 'output_path', 'output_path_organization', date_ref, True)
        # date_ref_str = date_ref.strftime('%Y%m%d')
        # output_name = f'{date_ref_str}_cmems_cnr_arc_rrs_resampled.nc'
        # file_output = os.path.join(output_dir_day, output_name)
        # if os.path.exists(file_output):
        #     print(f'[WARNING] Output file {file_output} already exist. Skiping...')
        #     date_ref = date_ref + timedelta(hours=24)
        #     continue
        unzip_path = options['unzip_path']
        make_resample_dir(input_dir, output_dir_day, unzip_path, arc_opt)
        date_ref = date_ref + timedelta(hours=24)


def run_integration(arc_opt, start_date, end_date):
    options = arc_opt.get_integrate_options()
    if options is None:
        return
    if args.verbose:
        print('[INFO] INTEGRATE OPTIONS:')
        for opt in options:
            print(f'[INFO]  {opt}->{options[opt]}')

    if start_date is None or end_date is None:
        start_date = options['start_date']
        end_date = options['end_date']

    platform = options['platform']
    date_run = start_date

    while date_run <= end_date:
        if args.verbose:
            print('*****************************')
            print(f'[INFO] Date: {date_run}')
        make_integration = True
        input_path = arc_opt.get_folder_date(options['input_path'], options['input_path_organization'], date_run, False)
        if not os.path.exists(input_path):
            print(f'[WARNING] Input path {input_path} for date {date_run} is not available. Skiping...')
            make_integration = False
        output_path = arc_opt.get_folder_date(options['output_path'], options['output_path_organization'], date_run,
                                              True)
        if output_path is None:
            print(f'[WARNING] Output path {input_path} for date {date_run} is not available. Skiping...')
            make_integration = False

        output_type = arc_opt.get_value_param('INTEGRATE', 'output_type', 'OPERATIVE', 'str')
        if output_type == 'CORRECT_RRS':
            alternative_path = arc_opt.get_folder_date(options['alternative_path'],
                                                       options['alternative_path_organization'], date_run, False)
            if not os.path.exists(alternative_path):
                print(f'[WARNING] Alternative path {alternative_path} for date {date_run} is not available. Skiping...')
                make_integration = False

        if make_integration:
            datestr = date_run.strftime('%Y%j')
            pl = platform[-1]
            if pl == '3':
                pl = ''
            from arc_integration import ArcIntegration
            dir_base = arc_opt.get_value_param('INTEGRATE', 'file_base', None, 'str')
            if args.verbose:
                print(f'[INFO] Output type: {output_type}')
            arc_integration = ArcIntegration(arc_opt, args.verbose, input_path, output_type, None)
            timeliness = arc_opt.get_value_param('INTEGRATE', 'timeliness', 'NT', 'str')
            if dir_base is not None:
                file_base = os.path.join(dir_base, 'ArcGrid_65_90_300m_AVERAGE_Base.nc')
                if os.path.exists(file_base):
                    arc_integration.ami.ifile_base = file_base
                    if args.verbose:
                        print(f'[INFO] File base: {file_base}')
            arc_integration.apply_pool = arc_opt.get_value_param('INTEGRATE', 'apply_pool', 0, 'int')
            arc_integration.timeliness = timeliness
            if args.verbose:
                print(f'[INFO] Timeliness: {timeliness}')
                print(f'[INFO] Input path: {input_path}')
                print(f'[INFO] Output file: {output_path}')
                print(f'[INFO] Apply pool: {arc_integration.apply_pool}')

            if output_type != 'CORRECT_RRS':
                arc_integration.make_integration(output_path)

            if output_type == 'RRS' or output_type == 'OPERATIVE' or output_type == 'CORRECT_RRS':
                file_base = os.path.join(dir_base, f'ArcGrid_65_90_300m_RRS_{timeliness}_Base.nc')
                arc_integration.ami.ifile_base = file_base
                name_out_end = f'O{pl}{datestr}_rrs-arc-fr.nc'
                file_out = os.path.join(output_path, name_out_end)
                if os.path.exists(file_out):
                    os.remove(file_out)
                arc_integration.output_type = 'RRS'
                if output_type == 'CORRECT_RRS':  ##COPY FILES FROM ALTERNATIVE PATH
                    if os.path.exists(alternative_path):
                        nfiles = len(os.listdir(alternative_path))
                        if nfiles == 3 or nfiles == 4:
                            for name in os.listdir(alternative_path):
                                fname = os.path.join(alternative_path, name)
                                fcopy = os.path.join(output_path, name)
                                copy_file(fname, fcopy)
                        elif nfiles == 16:
                            for name in os.listdir(alternative_path):
                                if name.find('_rrs-arc-fr.nc') > 0:
                                    continue
                                if name.find('_plankton-arc-fr.nc') > 0:
                                    continue
                                fname = os.path.join(alternative_path, name)
                                fcopy = os.path.join(output_path, name)
                                copy_file(fname, fcopy)
                            arc_integration.create_rrs_file(output_path, file_out, date_run, timeliness, True)

                else:
                    arc_integration.create_rrs_file(output_path, file_out, date_run, timeliness, False)

            if output_type == 'TRANSP' or output_type == 'OPERATIVE':
                file_base = os.path.join(dir_base, f'ArcGrid_65_90_300m_TRANSP_{timeliness}_Base.nc')
                arc_integration.ami.ifile_base = file_base
                name_out_end = f'O{pl}{datestr}_transp-arc-fr.nc'
                file_out = os.path.join(output_path, name_out_end)
                if os.path.exists(file_out):
                    os.remove(file_out)
                arc_integration.output_type = 'TRANSP'
                arc_integration.create_transp_file(output_path, file_out, date_run, timeliness)

        date_run = date_run + timedelta(hours=24)


def copy_file(ifile, ofile):
    if args.verbose:
        print(f'[INFO] Copying file : {ifile}...')

    cmd = f'cp -a {ifile} {ofile}'
    if args.verbose:
        print(f'[INFO] cmd: {cmd}')
    import subprocess
    import time
    subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    originalSize = os.path.getsize(ifile)
    historicalSize = -1
    while historicalSize != originalSize:
        if os.path.exists(ofile):
            historicalSize = os.path.getsize(ofile)
            if args.verbose:
                porc = (historicalSize / originalSize) * 100
                print(f'[INFO] Copying {porc:.2f} %')
        time.sleep(1)
    if args.verbose:
        print('[INFO] Copy completed')


def run_chla(arc_opt, start_date, end_date):
    options = arc_opt.get_processing_options()
    if options is None:
        return

    ##ONLY CHLA MODE IS IMPLEMENTED
    output_type = arc_opt.get_value_param('PROCESSING', 'output_type', 'CHLA', 'str')
    overwrite = arc_opt.get_value_param('PROCESSING', 'overwrite', False, 'boolean')
    if not output_type == 'CHLA':
        return

    if args.verbose:
        print('[INFO] PROCESSING OPTIONS:')
        for opt in options:
            print(f'[INFO]  {opt}->{options[opt]}')

    file_base = arc_opt.get_value_param('PROCESSING', 'file_base', None, 'str')
    timeliness = arc_opt.get_value_param('PROCESSING', 'timeliness', None, 'str')
    if file_base is not None:
        if os.path.exists(file_base):
            if timeliness is None:
                if file_base.find('NR') > 0:
                    timeliness = 'NR'
                if file_base.find('NT') > 0:
                    timeliness = 'NT'
    if args.verbose:
        print(f'[INFO] File base: {file_base}')
        print(f'[INFO] Timeliness: {timeliness}')

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
        arc_proc = ArcProcessing(arc_opt, args.verbose, output_type, None)
        arc_proc.compute_chla_image(input_file, output_file, timeliness)
        return

    ##WORKING WITH DATES
    if start_date is None or end_date is None:
        start_date = options['start_date']
        end_date = options['end_date']
    if start_date is None or end_date is None:
        return

    date_run = start_date
    if output_type=='CHLA':
        output_type = f'CHLA_{args.chla_algorithm}'
    arc_proc = ArcProcessing(arc_opt, args.verbose, output_type, None)

    while date_run <= end_date:
        if args.verbose:
            print('*****************************')
            print(f'[INFO] Date: {date_run}')
        make_processing = True
        input_path = arc_opt.get_folder_date(options['input_path'], options['input_path_organization'], date_run, False)
        dateyj = date_run.strftime('%Y%j')
        name_rrs = f'O{dateyj}_rrs-arc-fr.nc'
        input_file = os.path.join(input_path, name_rrs)
        if not os.path.exists(input_file):
            print(f'[WARNING] Input file {input_file} for date {date_run} is not available. Skiping...')
            make_processing = True
        output_path = arc_opt.get_folder_date(options['output_path'], options['output_path_organization'], date_run,
                                              True)
        if output_path is None and not overwrite:
            print(f'[WARNING] Output path {input_path} for date {date_run} is not available. Skipping...')
            make_processing = False

        output_name = f'O{dateyj}_plankton-arc-fr.nc'
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


def run_ql(arc_opt, start_date, end_date):
    options = arc_opt.get_ql_options()
    if options is None:
        return
    output_type = arc_opt.get_value_param('QL', 'output_type', 'CHL', 'str')
    if not output_type == 'CHL':
        return

    if args.verbose:
        print('[INFO] QL OPTIONS:')
        for opt in options:
            print(f'[INFO]  {opt}->{options[opt]}')

    from arc_mapinfo import ArcMapInfo
    ami = ArcMapInfo(None, args.verbose)

    ##WORKING WITH SINGLE GRANULE, ONLY CHLA
    input_name = arc_opt.get_value_param('QL', 'name_input', None, 'str')
    if input_name is not None:
        input_file = os.path.join(options['input_path'], input_name)
        if os.path.exists(input_file):
            if args.verbose:
                print(f'[INFO] Working with the single file: {input_file}')
        else:
            print(f'[ERROR] File {input_file} does not exist')
        output_name = arc_opt.get_value_param('QL', 'name_output', None, 'str')
        if output_name is None:
            output_name = 'QuickLookChla.nc'
        output_file = os.path.join(options['output_path'], output_name)
        if args.verbose:
            print(f'[INFO] Output file: {output_file}')
        # ami.save_quick_look_fdata(output_file, input_file, output_type)
        ami.save_full_fdata(output_file, input_file, output_file)
        return

    ##WORKING WITH DATES
    name_file_format_default = None
    name_file_date_format_default = '%Y%j'
    if output_type == 'CHL':
        name_file_format_default = 'O$DATE$_plankton-arc-fr.nc'
    name_file_format = arc_opt.get_value_param('QL', 'name_file_format', name_file_format_default, 'str')
    name_file_date_format = arc_opt.get_value_param('QL', 'name_file_date_format_default',
                                                    name_file_date_format_default, 'str')

    if start_date is None or end_date is None:
        start_date = options['start_date']
        end_date = options['end_date']
    date_run = start_date
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
            # ami.save_quick_look_fdata(output_file, input_file, output_type)
            ami.save_full_fdata(output_file, input_file, output_type)
        date_run = date_run + timedelta(hours=24)

    # file_out = args.outputpath
    # fdataset = args.product
    # # ami.save_quick_look_fdata(file_out, fdataset, 'sensor_mask')
    # ami.save_quick_look_fdata(file_out, fdataset, 'chla')


def compute_statistics(variable):
    # print(variable)
    width = variable.shape[1]
    height = variable.shape[2]
    ystep = 1000
    xstep = 1000
    import numpy.ma as ma
    nvalid_all = 0
    for y in range(0, height, ystep):
        # print(y)
        for x in range(0, width, xstep):
            try:
                limits = get_limits(y, x, ystep, xstep, height, width)
                array_lim = ma.array(variable[0, limits[0]:limits[1], limits[2]:limits[3]])
                nvalid = ma.count(array_lim)
                nvalid_all = nvalid_all + nvalid
            except:
                return -1

    return nvalid_all


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



# ami.make_resample_impl(olimage, file_out, granule_index, orbit_index, arc_opt)
# params: 0: ami; 1 olimage; 2: file_out; 3: granule_index; 4: orbit_index; 5: arc_opt; 6:param_mask
def make_resample_dir_parallel(params):
    ami = params[0]
    ami.make_resample_impl(params[1], params[2], params[3], params[4], params[5], params[6])
    path_prod_u = params[7]
    if path_prod_u is not None:
        for fn in os.listdir(path_prod_u):
            os.remove(os.path.join(path_prod_u, fn))


def make_resample_dir(dirorig, dirdest, unzip_path, arc_opt):
    section = 'RESAMPLE'
    doresample = arc_opt.get_value_param(section, 'do_resample', True, 'boolean')
    dokml = arc_opt.get_value_param(section, 'do_kml', False, 'boolean')
    timeliness = arc_opt.get_value_param(section, 'timeliness', 'NT', 'str')
    if dirdest is None:
        dirdest = dirorig
    apply_pool = arc_opt.get_value_param(section, 'apply_pool', 0, 'int')
    import zipfile as zp
    from arc_mapinfo import ArcMapInfo
    from olci_l2 import OLCI_L2

    ami = ArcMapInfo(arc_opt, args.verbose)
    if apply_pool != 0:
        params_list = []
        from multiprocessing import Pool
        if args.verbose:
            print('[INFO] Starting parallel processing')
    else:
        if args.verbose:
            print('[INFO] Starting sequencial processing')

    first_line = ['Source', 'StartDate', 'RelOrbit', 'GranuleIndex', 'OrbitIndex', 'OrigWidth', 'OrigHeight',
                  'OrigNTotal', 'OrigNFlagged', 'OrigNWater1', 'OrigNWatet2', 'OrigNValid', 'OrigPValid', 'YMin',
                  'YMax',
                  'XMin', 'XMax', 'Width', 'Height', 'NTotal', 'NValid', 'PValid']
    lines_out = [';'.join(first_line)]
    rel_pass_dict = {}
    idref = -1  # it's used for defining orbit index in rel_pass_dict as pow(2,idref)
    granule_index = 0  # 1,2,3...

    foutkml = None
    kml = None
    red = [166, 31, 178, 51, 251, 227, 253, 256, 202, 106, 255, 177, 0]
    green = [206, 120, 223, 160, 154, 26, 191, 127, 178, 61, 255, 89, 0]
    blue = [227, 180, 138, 44, 153, 28, 111, 0, 214, 154, 153, 40, 0]
    idcolor = 0

    idx = 1
    nfiles = len(os.listdir(dirorig))
    for name in os.listdir(dirorig):

        if args.verbose:
            print('------------------------------------------------------------------------')
            print(f'[INFO] File: {name} ({idx}/{nfiles})')
            idx = idx + 1
        pt = f'_{timeliness}_'
        if name.find(pt) < 0:
            continue
        prod_path = os.path.join(dirorig, name)
        do_zip_here = False
        if name.endswith('.SEN3'):
            path_prod_u = prod_path
        elif zp.is_zipfile(prod_path):
            do_zip_here = True
            path_prod_u = prod_path.split('/')[-1][0:-4]
            path_prod_u = os.path.join(unzip_path, path_prod_u)
        else:
            continue

        if path_prod_u.endswith('.SEN3'):
            output_name = path_prod_u.split('/')[-1][0:-5]
        else:
            output_name = path_prod_u.split('/')[-1]

        resample_done = False
        if doresample:
            file_out = os.path.join(dirdest, f'{output_name}_resampled.nc')
            if os.path.exists(file_out):
                print(f'[INFO] File {file_out} already exist. Skipping resample...')
                resample_done = True

        if do_zip_here:
            with zp.ZipFile(prod_path, 'r') as zprod:
                if args.verbose:
                    print(f'[INFO] Unziping {name} to {unzip_path}')
                zprod.extractall(path=unzip_path)

        olimage = OLCI_L2(path_prod_u, args.verbose)

        if olimage.check_granule():
            olimage.get_geo_and_params()
            granule_index = granule_index + 1
            # A new orbit_index is assigned to each platform + relative pass
            rel_pass = str(olimage.get_rel_pass())
            platform = olimage.get_platform()
            rel_pass = f'{platform}_{rel_pass}'
            if rel_pass not in rel_pass_dict.keys():
                idref = idref + 1
                rel_pass_dict[rel_pass] = math.pow(2, idref)
                if dokml:
                    import simplekml
                    if foutkml is not None and kml is not None:
                        kml.save(foutkml)
                        idcolor = idcolor + 1
                        if idcolor >= len(red):
                            idcolor = 0
                    foutkml = os.path.join(dirdest, f'Passes_RelativeOrbit_{rel_pass}.kml')
                    kml = simplekml.Kml()
            orbit_index = rel_pass_dict[rel_pass]
        else:
            print(f'[WARNING] Granule {output_name} is corrupted. Skipping granule...')
            dokml = False
            doresample = False
            ##Granule must be cancelled
            if zp.is_zipfile(prod_path):
                for fn in os.listdir(path_prod_u):
                    os.remove(os.path.join(path_prod_u, fn))
                os.remove(prod_path)
            else:
                for fn in os.listdir(prod_path):
                    os.remove(os.path.join(prod_path,fn))
                os.rename(prod_path,os.path.join(unzip_path,f'{output_name}.SEN3'))
            ##Granule name is save to a file of corrupted files
            f_corrupted_list = os.path.join(dirorig,'CorruptedFiles.txt')
            already_existed = os.path.exists(f_corrupted_list)
            fout = open(f_corrupted_list,'a')
            if already_existed:
                fout.write('\n')
            fout.write(output_name)
            fout.close()


        if doresample and not resample_done:
            if apply_pool == 0:
                line_out = ami.make_resample_impl(olimage, file_out, granule_index, orbit_index, arc_opt, None)
                if line_out is not None:
                    lines_out.append(line_out)
                if zp.is_zipfile(prod_path):
                    for fn in os.listdir(path_prod_u):
                        os.remove(os.path.join(path_prod_u, fn))
            else:
                params_mask = ami.check_make_resample_impl(olimage, arc_opt)
                if params_mask is None:  ##no resampling, deleting path_prod_u if working with zip file
                    if zp.is_zipfile(prod_path):
                        for fn in os.listdir(path_prod_u):
                            os.remove(os.path.join(path_prod_u, fn))
                else:
                    if args.verbose:
                        print('[INFO] Granule added to the parallel process list')
                    if apply_pool < 0:
                        params_granule = [ami, olimage, file_out, granule_index, orbit_index, arc_opt, None, None]
                    else:
                        params_granule = [ami, olimage, file_out, granule_index, orbit_index, arc_opt, params_mask,
                                          None]
                    if zp.is_zipfile(prod_path):
                        params_granule[7] = path_prod_u
                    params_list.append(params_granule)
                    if len(params_list) == apply_pool:
                        if args.verbose:
                            print(f'[INFO] Running parallel processes: {apply_pool}')
                        poolhere = Pool(apply_pool)
                        poolhere.map(make_resample_dir_parallel, params_list)
                        params_list = []

        # if resample is done, we must delete unzipped file
        if not dokml and doresample and resample_done:
            if zp.is_zipfile(prod_path):
                for fn in os.listdir(path_prod_u):
                    os.remove(os.path.join(path_prod_u, fn))

        if dokml:
            start_date = olimage.get_start_date()
            start_date_s = 'UNKNOWN START DATE'
            if start_date is not None:
                start_date_s = start_date.strftime('%Y%m%dT%H%M:%S')
            lin = kml.newlinestring(name=start_date_s, description=name, coords=olimage.coords_image)
            lin.style.linestyle.color = simplekml.Color.rgb(red[idcolor], green[idcolor], blue[idcolor], 255)
            lin.style.linestyle.width = 3
            # name = name_list[7], description = name, coords = coordinates

        if dokml and not doresample:
            if zp.is_zipfile(prod_path):
                for fn in os.listdir(path_prod_u):
                    os.remove(os.path.join(path_prod_u, fn))

    if doresample and apply_pool != 0 and len(params_list) > 0:

        if args.verbose:
            print(f'[INFO]******************************************************************************************')
            print(f'[INFO] Starting parallel processing. Number of products: {len(params_list)}')
            print(f'[INFO] CPUs: {os.cpu_count()}')
            print(f'[INFO] Parallel processes: {apply_pool}')
            print(f'[INFO]******************************************************************************************')

        if apply_pool < 0:
            poolhere = Pool()
        else:
            poolhere = Pool(apply_pool)

        poolhere.map(make_resample_dir_parallel, params_list)

    if doresample and apply_pool == 0:
        fcsvout = os.path.join(dirdest, 'ResampleInfo.csv')
        if args.verbose:
            print(f'[INFO] Creating info file:{fcsvout}')
        fw = open(fcsvout, 'w')
        for line in lines_out:
            fw.write(line)
            fw.write('\n')
        fw.close()


def get_dates_from_arg():
    from datetime import datetime as dt
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
