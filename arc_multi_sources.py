import calendar
import os
from ftplib import FTP
from datetime import datetime as dt

class ARC_MULTI_SOURCES:
    def __init__(self, input_dir, input_dir_org,moi_user,moi_pass,use_myint_sources, verbose):
        self.verbose = verbose
        self.input_dir = input_dir
        if input_dir_org is None:
            self.input_dir_org = 'YYYY/jjj'
        self.input_dir_org = input_dir_org
        self.name_product = 'OCEANCOLOUR_GLO_BGC_L3_MY_009_107'
        self.name_dataset = 'c3s_obs-oc_glo_bgc-reflectance_my_l3-multi-4km_P1D'
        self.use_myint_sources = use_myint_sources
        self.server = 'my.cmems-du.eu'
        self.variable_check = ['latitude','longitude','time','RRS412','RRS443','RRS490','RRS510','RRS560','RRS665']

        if moi_user is not None and moi_pass is not None:
            self.moi_user = moi_user
            self.moi_pass = moi_pass

    # def get_output_file(self,date_here,datasetname,output_path,output_path_organization):
    #   print('fdf')

    def get_file_date(self, date_here, make_download):
        dirdate = self.get_input_folder_date(date_here, True)
        if dirdate is None:
            return None
        date_here_str = date_here.strftime('%Y%m%d')
        name = f'{date_here_str}_{self.name_dataset}.nc'

        fdate = os.path.join(dirdate, name)

        if not os.path.exists(fdate):
            name_myint = name.replace('my','myint')
            fdate_myint = os.path.join(dirdate,name_myint)
            if os.path.exists(fdate_myint):
                fdate = os.path.join(dirdate,name_myint)

        if not os.path.exists(fdate) and make_download:
            self.download_file(date_here, dirdate, name)

        if os.path.exists(fdate):
            from netCDF4 import Dataset
            try:
                dataset = Dataset(fdate)
                for var in self.variable_check:
                    if not var in dataset.variables:
                        print(f'[ERROR] Variable {var} is not available in file: {fdate}')
                        dataset.close()
                        return 'INVALIDFILE'
                dataset.close()
                return fdate
            except:
                print(f'[ERROR] File {fdate} is not a valid NetCDF file')
                return 'INVALIDFILE'

        else:
            return 'NOFILE'


    def check_file_ftp(self,date_here):

        date_here_str = date_here.strftime('%Y%m%d')
        name = f'{date_here_str}_{self.name_dataset}.nc'
        fwork = FTPWork(self.server, self.moi_user, self.moi_pass)
        fwork.path_dataset = f'/Core/{self.name_product}/{self.name_dataset}'
        return fwork.check_daily_file(date_here,name)

    def check_month_files_ftp(self,year,month):
        fwork = FTPWork(self.server, self.moi_user, self.moi_pass)
        fwork.path_dataset = f'/Core/{self.name_product}/{self.name_dataset}'
        missing_dates = fwork.check_daily_files_month(year,month,self.name_dataset)
        return missing_dates

    def download_file(self, date_here, dirdate, name):
        if self.verbose:
            print(f'[INFO] Downloading file {name} in folder: {dirdate}')

        fwork = FTPWork(self.server,self.moi_user,self.moi_pass)
        fwork.path_dataset = f'/Core/{self.name_product}/{self.name_dataset}'
        fwork.donwload_daily_file(date_here,dirdate,name)
        fwork.close()


    def get_input_folder_date(self, date_here, create):
        return self.get_folder_date(date_here,create,self.input_dir,self.input_dir_org)

    def get_folder_date(self,date_here,create, input_dir,input_dir_org):
        replaces = {
            'YYYY': date_here.strftime('%Y'),
            'jjj': date_here.strftime('%j'),
            'mm': date_here.strftime('%m'),
            'dd': date_here.strftime('%d')
        }

        dirs = input_dir_org.split('/')
        dirdate = input_dir
        for d in dirs:
            dhere = d
            for r in replaces:
                if d.find(r) >= 0:
                    dhere = dhere.replace(r, replaces[r])
            dirdate = os.path.join(dirdate, dhere)
            if not os.path.exists(dirdate) and create:
                try:
                    os.mkdir(dirdate)
                except:
                    pass

        if not os.path.exists(dirdate):
            if create:
                print(f'[ERROR] Folder {dirdate} could not be created')
            else:
                print(f'[WARNING] Folder {dirdate} does not exist')
            return None

        return dirdate

class FTPWork():

    def __init__(self, server,uname,passwd):
        self.server = server
        self.ftpdu = FTP(server, uname, passwd)
        self.path_dataset = None

    def go_subdir(self, rpath):
        # print('Changing directory to: ', rpath)
        self.ftpdu.cwd(rpath)

    def go_year_subdir(self, rpath, year):
        if rpath is None:
            rpath = self.path_dataset
        if rpath is None:
            return False
        dateref = dt(year, 1, 1)
        yearstr = dateref.strftime('%Y')
        self.ftpdu.cwd(rpath)
        if yearstr in self.ftpdu.nlst():
            self.ftpdu.cwd(yearstr)
            return True
        else:
            return False

    def go_month_subdir(self, rpath, year, month):
        if rpath is None:
            rpath = self.path_dataset
        if rpath is None:
            return False
        dateref = dt(year, month, 1)
        yearstr = dateref.strftime('%Y')  # dt.strptime(str(year), '%Y')
        monthstr = dateref.strftime('%m')  # dt.strptime(month, '%m')
        self.ftpdu.cwd(rpath)
        if yearstr in self.ftpdu.nlst():
            self.ftpdu.cwd(yearstr)
            if monthstr in self.ftpdu.nlst():
                self.ftpdu.cwd(monthstr)
                return True

        return False

    def check_daily_file(self,date_here,name_file):
        b = self.go_month_subdir(None,date_here.year,date_here.month)
        if not b:
            return False

        if name_file in self.ftpdu.nlst():
            return True
        else:
            return False

    def check_daily_files_month(self,year,month,name_ref):

        res = calendar.monthrange(year,month)
        last_day_month = res[1]
        b = self.go_month_subdir(None, year, month)
        if not b:
            missing_dates = []
            for day in range(1, last_day_month + 1):
                date_here = dt(year, month, day)
                missing_dates.append(date_here)
            return missing_dates
        list = self.ftpdu.nlst()
        missing_dates = []
        for day in range(1,last_day_month+1):
            date_here = dt(year,month,day)
            date_here_str = date_here.strftime('%Y%m%d')
            name_file = f'{date_here_str}_{name_ref}.nc'
            if name_file not in list:
                missing_dates.append(date_here)

        return missing_dates


    def donwload_daily_file(self,date_here,dir_out,name_file):
        b = self.check_daily_file(date_here,name_file)
        if not b:
            print(f'[ERROR] {self.path_dataset}/{name_file} is not available for download from {self.server}')
            return False
        # Write file in binary mode
        file_out = os.path.join(dir_out,name_file)
        with open(file_out, "wb") as file:
            # Command for Downloading the file "RETR filename"
            self.ftpdu.retrbinary(f"RETR {name_file}", file.write)

        return True

    def close(self):
        self.ftpdu.close()