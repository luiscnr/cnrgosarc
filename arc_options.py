import os.path
from datetime import datetime as dt


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
    start_date = None
    end_date = None
    try:
        start_date = dt.strptime(start_date_str, '%Y-%m-%d')
    except:
        print(f'[ERROR] Start date {start_date_str} is not valid')
    try:
        end_date = dt.strptime(end_date_str, '%Y-%m-%d')
    except:
        print(f'[ERROR] End date {end_date_str} is not valid')
    return start_date, end_date


def check_folder_organization(organization):
    valid_options = ['YYYYmmdd', 'YYYY/mm/dd', 'YYYY/jjj']
    if organization in valid_options:
        return True
    else:
        print(f'[ERROR] Folder organization {organization} is not implemented')
        print(f'[ERROR] Valid optons: {valid_options}')
        return False


class ARC_OPTIONS:
    def __init__(self, options):
        # import configparser
        # options = configparser.ConfigParser()
        # options.read(config_file)
        self.options = options

    def get_integrate_options(self):
        if not self.options.has_section('INTEGRATE'):
            print(f'[ERROR] Integrate section is not included in config file')
            return None
        compulsory_options = ['input_path', 'output_path', 'start_date', 'end_date']
        for coption in compulsory_options:
            if not self.options.has_option('INTEGRATE', coption):
                print(f'[ERROR] Option INTEGRATE/{coption} is compulsory')
                return None
        input_path = self.options['INTEGRATE']['input_path']
        output_path = self.options['INTEGRATE']['output_path']
        start_date_str = self.options['INTEGRATE']['start_date']
        end_date_str = self.options['INTEGRATE']['end_date']
        if not os.path.exists(input_path):
            print(f'[ERROR] Input path {input_path} does not exist')
            return None
        if not os.path.exists(output_path):
            try:
                os.mkdir(output_path)
            except:
                print(f'[ERROR] Output path {output_path} does not exist and coul not be created')
        start_date, end_date = get_dates_from_stroptions(start_date_str, end_date_str)
        if start_date is None or end_date is None:
            return None

        input_path_organization = None
        if self.options.has_option('INTEGRATE', 'input_path_organization'):
            input_path_organization = self.options['INTEGRATE']['input_path_organization']
            if not check_folder_organization(input_path_organization):
                return None

        output_path_organization = None
        if self.options.has_option('INTEGRATE', 'output_path_organization'):
            output_path_organization = self.options['INTEGRATE']['output_path_organization']
            if not check_folder_organization(output_path_organization):
                return None

        platform = 'S3'
        if self.options.has_option('INTEGRATE', 'platform'):
            platform = self.options['INTEGRATE']['platform']
            valid_values = ['S3', 'S3A', 'S3B']
            if not platform in valid_values:
                print(f'[ERROR] Platform value {platform} is not valid (valid options: S3, S3A, S3B')
                return None

        options_out = {
            'input_path': input_path,
            'output_path': output_path,
            'input_path_organization': input_path_organization,
            'output_path_organization': output_path_organization,
            'start_date': start_date,
            'end_date': end_date,
            'platform': platform
        }
        return options_out

    def get_value(self, section, key):
        value = None
        if self.options.has_option(section, key):
            value = self.options[section][key]
        return value

    def get_value_param(self, section, key, default, type):
        value = self.get_value(section, key)
        if value is None:
            return default
        if type == 'str':
            return value
        if type == 'int':
            return int(value)
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

        return path_date
