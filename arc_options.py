import os.path
from datetime import datetime as dt
from datetime import timedelta

from fontTools.misc.plistlib import end_date


def create_folder(folder):
    if os.path.exists(folder):
        return True
    try:
        os.mkdir(folder)
    except:
        print(f'[ERROR] Path  {folder} could not be created')
        return False
    return True


def get_dates_from_stroptions(start_date_str, end_date_str):
    from datetime import timedelta
    start_date = None
    end_date = None
    try:
        start_date = dt.strptime(start_date_str, '%Y-%m-%d')
    except:
        try:
            tdelta = int(start_date_str)
            start_date = dt.now() + timedelta(days=tdelta)
            start_date = start_date.replace(hour=12, minute=0, second=0, microsecond=0)
        except:
            print(f'[ERROR] Start date {start_date_str} is not valid')
    try:
        end_date = dt.strptime(end_date_str, '%Y-%m-%d')
    except:
        try:
            tdelta = int(end_date_str)
            end_date = dt.now() + timedelta(days=tdelta)
            end_date = end_date.replace(hour=12, minute=0, second=0, microsecond=0)
        except:
            print(f'[ERROR] End date {end_date_str} is not valid')
    return start_date, end_date


def check_folder_organization(organization):
    valid_options = ['YYYYmmdd', 'YYYY/mm/dd', 'YYYY/jjj','YYYY','YYYY/mm']
    if organization in valid_options:
        return True
    else:
        print(f'[ERROR] Folder organization {organization} is not implemented')
        print(f'[ERROR] Valid options: {valid_options}')
        return False


class ARC_OPTIONS:
    def __init__(self, options):
        # import configparser
        # options = configparser.ConfigParser()
        # options.read(config_file)
        self.options = options

    def add_moi_credentials(self,options):
        if options is None:
            return options
        section = 'GENERAL'

        options['moi_user'] = self.get_value_param(section,'moi_user','lgonzalezvilas','str')
        options['moi_pass'] = self.get_value_param(section,'moi_pass','MegaRoma17!','str')

        return options

    def get_resample_options(self):
        # input path, output path, path organizations, platform, start date, end date
        section = 'RESAMPLE'
        options = self.get_basic_options(section)
        if options is None:
            return options
        if options['platform']=='MULTI':
            return options
        unzip_path = self.get_path(section, 'unzip_path', False)
        if unzip_path is None:
            unzip_path = options['input_path']
            print(f'[WARNING] Indepentent unzip_path was not defined, using input path as unzip_path')
        options['unzip_path'] = unzip_path



        return options

    def get_integrate_options(self):
        return self.get_basic_options('INTEGRATE')

    def get_processing_options(self):
        return self.get_basic_options('PROCESSING')


    def get_ql_options(self):
        return self.get_basic_options('QL')

    def get_basic_options(self, section):
        # section = 'INTEGRATE'
        if not self.options.has_section(section):
            print(f'[ERROR] {section} section is not included in config file')
            return None
        #compulsory_options = ['input_path', 'output_path', 'start_date', 'end_date']
        compulsory_options = ['input_path', 'output_path']
        if not self.check_compulsory_options(section, compulsory_options):
            return None
        input_path = self.get_path(section, 'input_path', False)
        output_path = self.get_path(section, 'output_path', True)
        if input_path is None or output_path is None:
            return None
        if self.options.has_option(section,'start_date'):
            start_date_str = self.options[section]['start_date']
        else:
            start_date_str = (dt.now()-timedelta(days=1)).strftime('%Y-%m-%d')
        if self.options.has_option(section, 'end_date'):
            end_date_str = self.options[section]['end_date']
        else:
            end_date_str = start_date_str
        start_date, end_date = get_dates_from_stroptions(start_date_str, end_date_str)
        if start_date is None or end_date is None:
            return None

        climatology_path = self.get_path(section,'climatology_path',False)
        use_myint_sources = self.get_value_param(section, 'use_myint_sources', False, 'boolean')

        input_path_organization = self.get_path_organization(section, 'input_path_organization')
        output_path_organization = self.get_path_organization(section, 'output_path_organization')
        if input_path_organization == 'INVALID' or output_path_organization == 'INVALID':
            return None

        platform = self.get_platform(section)

        alternative_path = self.get_path(section, 'alternative_path', False)
        alternative_path_organization = self.get_path_organization(section, 'alternative_path_organization')

        options_out = {
            'input_path': input_path,
            'output_path': output_path,
            'input_path_organization': input_path_organization,
            'output_path_organization': output_path_organization,
            'start_date': start_date,
            'end_date': end_date,
            'platform': platform,
            'alternative_path': alternative_path,
            'alternative_path_organization': alternative_path_organization,
            'climatology_path': climatology_path,
            'use_myint_sources': use_myint_sources
        }
        return options_out

    def get_path_organization(self, section, key):
        if self.options.has_option(section, key):
            path_organization = self.options[section][key]
            if not check_folder_organization(path_organization):
                return 'INVALID'
            return path_organization
        else:
            return None

    def get_path(self, section, key, create):
        if not self.options.has_option(section,key):
            print(f'[WARNING] {key} is not available in {section}')
            return None
        path = self.options[section][key]
        if not os.path.exists(path):
            if not create:
                print(f'[ERROR] {key} {path} does not exist')
                return None
            else:
                try:
                    os.mkdir(path)
                except:
                    print(f'[ERROR] {key} {path} does not exist and could not be created')
                    return None
        return path

    def get_platform(self, section):
        platform = 'S3'
        if self.options.has_option(section, 'platform'):
            platform = self.options[section]['platform']
            valid_values = ['S3', 'S3A', 'S3B','MULTI']
            if not platform in valid_values:
                print(f'[ERROR] Platform value {platform} is not valid (valid options: S3, S3A, S3B')
                return None
        return platform

    def check_compulsory_options(self, section, compulsory_options):
        for coption in compulsory_options:
            if not self.options.has_option(section, coption):
                print(f'[ERROR] Option {section}/{coption} is compulsory')
                return False
        return True

    def get_value(self, section, key):
        value = None
        if self.options.has_option(section, key):
            value = self.options[section][key].strip()
        return value

    def get_value_param(self, section, key, default, type):
        value = self.get_value(section, key)
        if value is None:
            return default
        if type == 'str':
            return value
        if type == 'int':
            return int(value)
        if type == 'file':
            if not os.path.exists(value):
                return default
            else:
                return value
        if type == 'boolean':
            if value == '1' or value.upper() == 'TRUE':
                return True
            elif value == '0' or value.upper() == 'FALSE':
                return False
            else:
                return True
        if type == 'rrslist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip().replace('.', '_')
                list.append(f'RRS{vals}')
            return list
        if type == 'strlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                list.append(vals.strip())
            return list
        if type == 'floatlist':
            list_str = value.split(',')
            list = []
            for vals in list_str:
                vals = vals.strip()
                list.append(float(vals))
            return list

    def get_folder_date_o(self, options, key_path, key_organization, date_here, create):
        path_base = options[key_path]
        org = options[key_organization]
        return self.get_folder_date(path_base, org, date_here, create)

    def get_list_files_month(self, path_base, org, year, month, file_date, file_date_format, timeliness):
        from calendar import monthrange
        from netCDF4 import Dataset
        ndays = monthrange(year, month)[1]
        list_files = []
        list_files_timeliness = []
        #ntimeliness = 0
        for iday in range(1, ndays + 1):
            date_here = dt(year, month, iday)
            folder_date = self.get_folder_date(path_base, org, date_here, False)
            date_here_str = date_here.strftime(file_date_format)
            name_file = file_date.replace('DATE', date_here_str)
            file_out = os.path.join(folder_date, name_file)
            if os.path.exists(file_out):
                list_files.append(file_out)
                dataset = Dataset(file_out)
                if 'timeliness' in dataset.ncattrs():
                    if timeliness == dataset.timeliness:
                        #ntimeliness = ntimeliness + 1
                        list_files_timeliness.append(file_out)
                        print(f'[INFO] File {file_out} is added for computing the average..')
                    else:
                        print(f'[WARNING] File {file_out} is {timeliness} but should be {timeliness}. Added anyway...')
                dataset.close()
            else:
                print(f'[WARNING] File {file_out} is not available for computing the average.')
        return list_files, list_files_timeliness

    def get_folder_year(self,path_base,date_here,create):
        year_str = date_here.strftime('%Y')
        path_year = os.path.join(path_base, year_str)
        if create:
            if not create_folder(path_year):
                return None
        return path_year

    def get_folder_date(self, path_base, org, date_here, create):
        if org is None:
            return path_base
        if org == 'YYYYmmdd':
            date_here_str = date_here.strftime('%Y%m%d')
            path_date = os.path.join(path_base, date_here_str)
            if create:
                if not create_folder(path_date):
                    return None

        if org == 'YYYY/mm/dd':
            year_str = date_here.strftime('%Y')
            month_str = date_here.strftime('%m')
            day_str = date_here.strftime('%d')
            path_year = os.path.join(path_base, year_str)
            path_month = os.path.join(path_year, month_str)
            path_date = os.path.join(path_month, day_str)
            if create:
                if not create_folder(path_year):
                    return None
                if not create_folder(path_month):
                    return None
                if not create_folder(path_date):
                    return None
        if org == 'YYYY/mm':
            year_str = date_here.strftime('%Y')
            month_str = date_here.strftime('%m')
            path_year = os.path.join(path_base, year_str)
            path_date = os.path.join(path_year, month_str)
            if create:
                if not create_folder(path_year):
                    return None
                if not create_folder(path_date):
                    return None

        if org == 'YYYY/jjj':
            year_str = date_here.strftime('%Y')
            jjj_str = date_here.strftime('%j')
            path_year = os.path.join(path_base, year_str)
            path_date = os.path.join(path_year, jjj_str)
            if create:
                if not create_folder(path_year):
                    return None
                if not create_folder(path_date):
                    return None

        if org == 'YYYY':
            year_str = date_here.strftime('%Y')
            path_date = os.path.join(path_base, year_str)
            if create:
                if not create_folder(path_date):
                    return None

        return path_date
