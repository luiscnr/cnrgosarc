import argparse
import warnings
from datetime import datetime as dt
from datetime import timedelta
import os

import numpy as np

warnings.filterwarnings(action='ignore', category=ResourceWarning)

parser = argparse.ArgumentParser(description="Artic resampler")
parser.add_argument('-i', "--inputpath", help="Input directory", required=True)
parser.add_argument('-o', "--outputpath", help="Output dirctory", required=True)
parser.add_argument('-sd', "--start_date", help="Start date (yyyy-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date (yyyy-mm-dd")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
args = parser.parse_args()


def main():
    print('[INFO] Started correct reflectance')
    input_path = args.inputpath
    output_path = args.outputpath
    start_date, end_date = get_dates_from_arg()
    if not os.path.isdir(input_path):
        print(f'[ERROR] {input_path} does not exist or is not a directory')
        return
    if not os.path.isdir(output_path):
        try:
            os.mkdir(output_path)
        except:
            print(f'[ERROR] Output path: {output_path} does not exist and could not be created')
            return
    if start_date is None or end_date is None:
        return

    date_here = start_date
    while date_here <= end_date:
        if args.verbose:
            print(f'[INFO] Date: {date_here}')

        input_dirdate = get_folder_date(input_path, date_here, False)
        output_dirdate = get_folder_date(output_path, date_here, True)
        if input_dirdate is None:
            print(f'[ERROR] Input dir date does not exist. Skipping...')
            return
        if output_dirdate is None:
            print(f'[ERROR] Input dir date: {output_dirdate} is not available. Skipping...')
            return
        yyyy = date_here.strftime('%Y')
        jjj = date_here.strftime('%j')
        input_file = os.path.join(input_dirdate, f'O{yyyy}{jjj}_rrs-arc-fr.nc')
        if not os.path.exists(input_file):
            print(f'[ERROR] Input file: {input_file} does not exist. Skipping...')
            return
        output_file = os.path.join(output_dirdate, f'O{yyyy}{jjj}_rrs-arc-fr.nc')
        correct_rrs_impl(input_file, output_file)

        date_here = date_here + timedelta(hours=24)


def correct_rrs_impl(input_file, output_file):
    from netCDF4 import Dataset

    copy_nc_base(input_file, output_file)
    height = 18345
    width = 18345
    ystep = 6500
    xstep = 6500
    dst = Dataset(output_file, 'a', format='NETCDF4')
    for var in dst.variables:

        if var.startswith('RRS'):
            if args.verbose:
                print(f'[INFO] Updating variable: {var}')

            variable = dst.variables[var]
            for y in range(0, height, ystep):
                for x in range(0, width, xstep):
                    limits = get_limits(y, x, ystep, xstep,height, width)
                    array = np.array(variable[0, limits[0]:limits[1], limits[2]:limits[3]])
                    array[array != -999] = array[array!=-999]/np.pi
                    variable[0, limits[0]:limits[1], limits[2]:limits[3]] = [array[:, :]]
            #dst.variables[var] = variable
    dst.close()
    if args.verbose:
        print('COMPLETED')

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

def copy_nc_base(ifile, ofile):
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


def get_folder_date(dirbase, datehere, create):
    try:
        yyyy = datehere.strftime('%Y')
        jjj = datehere.strftime('%j')
        dirdate = os.path.join(dirbase, yyyy, jjj)

        if os.path.isdir(dirdate):
            return dirdate
        else:
            if create:
                try:
                    diryear = os.path.join(dirbase, yyyy)
                    if not os.path.isdir(diryear):
                        os.mkdir(diryear)
                    os.mkdir(dirdate)
                    return dirdate
                except:
                    return None
            else:
                return None
    except:
        return None


def get_dates_from_arg():
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
