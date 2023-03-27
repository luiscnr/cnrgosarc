import argparse

parser = argparse.ArgumentParser(description="Artic config options")
parser.add_argument("-m", "--mode", help="Mode", choices=["GET", "SET", "LOGFOLDER"], required=True)
# parser.add_argument("-p", "--product", help="Input product (testing)")
parser.add_argument('-i', "--inputpath", help="Input directory")
parser.add_argument('-pr', "--prefix",help="Prefix for getting log folder")
# parser.add_argument('-o', "--outputpath", help="Output file")
# parser.add_argument('-tp', "--temp_path", help="Temporary directory")
parser.add_argument('-d', "--date", help="Start date (yyyy-mm-dd)")
# parser.add_argument('-ed', "--end_date", help="End date (yyyy-mm-dd")
parser.add_argument('-c', "--config_file", help="Configuration file (Default: arc_config.ini)")
parser.add_argument('-s', "--section", help="Config file section")
parser.add_argument('-k', "--key", help="Config file key")
parser.add_argument('-val', "--value", help="Config file key value to be set")
# parser.add_argument('-bf', "--base_file", help="Create base file for the following mode: RESAMPLE,INTEGRATE.",
#                     action="store_true")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
args = parser.parse_args()


def main():
    #print('[INFO] Started config')
    if args.mode=='SET':
        if not args.config_file or not args.section or not args.key or not args.value:
            return
        import configparser
        options = configparser.ConfigParser()
        options.read(args.config_file)

        options[args.section][args.key] = args.value.strip()

        with open(args.config_file,'w') as configfile:
            options.write(configfile)
        if args.verbose:
            print(f'[INFO] {args.section}/{args.key} set to {args.value} in config file: {args.config_file}')

    if args.mode=='GET':
        if not args.config_file or not args.section or not args.key:
            return
        import configparser
        options = configparser.ConfigParser()
        options.read(args.config_file)
        print(options[args.section][args.key])

    if args.mode=='LOGFOLDER':
        if not args.inputpath or not args.date:
            return
        from datetime import datetime as dt
        import os
        date = dt.strptime(args.date,'%Y-%m-%d')
        datestr = date.strftime('%Y%m%d')
        prefix = ''
        if args.prefix:
            prefix = args.prefix

        name_base = f'{prefix}{datestr}_'
        logfolder = None
        exist = True
        index = 1
        while exist:
            fhere = os.path.join(args.inputpath,f'{name_base}{index}')
            if os.path.exists(fhere) and os.path.isdir(fhere):
                index = index+1
            else:
                try:
                    os.mkdir(fhere)
                    logfolder = fhere
                    exist = False
                except:
                    exist = False
        if logfolder is not None:
            print(logfolder)
            if args.config_file and os.path.exists(args.config_file):
                import configparser
                options = configparser.ConfigParser()
                options.read(args.config_file)
                options['GENERAL']['log_folder'] = logfolder
                with open(args.config_file, 'w') as configfile:
                    options.write(configfile)
            return logfolder



if __name__ == '__main__':
    main()
