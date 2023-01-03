import math
import os
import argparse

# import netCDF4
# import numpy as np
# import zipfile as zp
# from arc_mapinfo import ArcMapInfo
# from arc_integration import ArcIntegration
# from olci_l2 import OLCI_L2
# import simplekml

parser = argparse.ArgumentParser(description="Artic resampler")
parser.add_argument("-m", "--mode", help="Mode", choices=["CHECKPY", "CHECK", "GRID", "RESAMPLE", "RESAMPLEPML", "QL"],
                    required=True)
parser.add_argument("-p", "--product", help="Input product (testing)")
parser.add_argument('-i', "--inputpath", help="Input directory")
parser.add_argument('-o', "--outputpath", help="Output file")
parser.add_argument('-tp', "--temp_path", help="Temporary directory")
parser.add_argument('-sd', "--start_date", help="Start date (yyyy-mm-dd)")
parser.add_argument('-ed', "--end_date", help="End date (yyyy-mm-dd")
parser.add_argument('-c', "--config_file", help="Configuration file (Default: arc_config.ini)")
parser.add_argument("-v", "--verbose", help="Verbose mode.", action="store_true")
args = parser.parse_args()


def main():
    print('[INFO] Started Artic Resampler')
    if args.mode == "CHECKPY":
        check_py()
        return

    from arc_mapinfo import ArcMapInfo

    from olci_l2 import OLCI_L2

    if args.mode == "CHECK":
        do_resampled_vm_christmas()
        # do_check6()
        # dirorig = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/15072018'
        # unzip_path = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/temp'
        # make_resample_dir(dirorig, dirorig,unzip_path, True, False)
        return

    fconfig = None
    if args.config_file:
        fconfig = args.config_file

    ami = ArcMapInfo(fconfig, args.verbose)

    if args.mode == "CHECK":
        return

    file_out = args.outputpath

    if args.mode == 'GRID':  ##creating grid
        print('create grid')
        # lon,lat = ami.area_def.get_lonlat_from_projection_coordinates(0,0)
        # file = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/ArcGrid_65_90_300mORIGts70.nc'
        # from netCDF4 import Dataset
        # ncsat = Dataset(file,'w')
        # ncsat.variables['sensor_mask'].grid_mapping = 'latitude_longitude'
        # ncsat.close()
        ami.create_nc_filegrid(file_out, True, False)

    if args.mode == "RESAMPLE" and not args.inputpath and args.product:  ##testing, resampling of a single granule
        folci = args.product
        olimage = OLCI_L2(folci, args.verbose)
        olimage.get_geo_and_params()
        line_output = ami.make_resample_impl(olimage, file_out, 1, -1)
        print(line_output)

    if args.mode == "RESAMPLEPML" and args.product:  ##testing, resampling of a PML file
        fpml = args.product
        ami.make_resample_pml(fpml, file_out)

    if args.mode == 'QL' and args.product:
        fdataset = args.product
        ami.save_quick_look_fdata(file_out, fdataset, 'sensor_mask')

        # olimage = OLCI_L2(folci, args.verbose)
        # ami.make_resample_impl(olimage, file_out, 1)

    # fgrid = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/ArcGrid_65_90_300m.nc'
    # fout = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/ArcGridOlciQuickLook.png'
    # ami.create_nc_filegrid(fgrid,True)
    # ami.save_quick_look_fgrid(fout,fgrid)

    # folci = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/S3A_OL_2_WFR____20180715T110812_20180715T111112_20211121T193318_0180_033_251______MAR_R_NT_003.SEN3'
    # fout = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/test.nc'
    # olimage = OLCI_L2(folci)
    # array, info = olimage.get_observation_angle()
    # print(info)

    # ami.make_resample_impl(olimage,fout)

    # fdata = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/test.nc'
    # fout = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/test.png'
    # ami.save_quick_look_fdata(fout,fdata)


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

    if not check:
        print('[ERROR] Some packages are not available...')
    else:
        print('[INFO] All the packages are available')


def do_check6():
    import numpy as np

    file_in = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/NERC_EXAMPLES/202207_mm-metno-MODEL-topaz4-ARC-fv02.0.nc'
    from netCDF4 import Dataset
    dataset = Dataset(file_in)
    xcoords = dataset.variables['x']
    ycoords = dataset.variables['y']
    xmin = xcoords[0] * 100000
    xmax = xcoords[-1] * 100000
    ymin = ycoords[0] * 100000
    ymax = ycoords[-1] * 100000
    print(xmin, xmax, ymin, ymax)
    latitude = dataset.variables['latitude']
    longitude = dataset.variables['longitude']
    print(np.min(latitude[:]), np.max(latitude[:]), np.min(longitude[:]), np.max(longitude[:]))
    from arc_mapinfo import ArcMapInfo
    ami = ArcMapInfo(None, True)

    x0 = xcoords[0] * 100000
    y0 = ycoords[0] * 100000
    xn = xcoords[608] * 100000
    yn = ycoords[880] * 100000
    lon1, lat1 = ami.area_def.get_lonlat_from_projection_coordinates(x0, y0)
    lon2, lat2 = ami.area_def.get_lonlat_from_projection_coordinates(xn, y0)
    lon3, lat3 = ami.area_def.get_lonlat_from_projection_coordinates(x0, yn)
    lon4, lat4 = ami.area_def.get_lonlat_from_projection_coordinates(xn, yn)
    print(lon1, lat1)
    print(lon2, lat2)
    print(lon3, lat3)
    print(lon4, lat4)

    lon1, lat1 = ami.area_def.get_lonlat_from_array_coordinates(0, 0)
    lon2, lat2 = ami.area_def.get_lonlat_from_array_coordinates(608, 0)
    lon3, lat3 = ami.area_def.get_lonlat_from_array_coordinates(0, 880)
    lon4, lat4 = ami.area_def.get_lonlat_from_array_coordinates(608, 880)
    print('..........................................')
    print(lon1, lat1)
    print(lon2, lat2)
    print(lon3, lat3)
    print(lon4, lat4)
    print('##################################################################')

    for ypix in range(881):
        xc, yc = ami.area_def.get_lonlat_from_array_coordinates(0, ypix)
        print(ypix, '->', yc)

    # lats = np.zeros((881,609))
    # lons = np.zeros((881,609))
    # for y in range(881):
    #     print(y)
    #     for x in range(609):
    #         xc = xcoords[x]*100000
    #         yc = ycoords[y]*100000
    #         lon, lat = ami.area_def.get_lonlat_from_projection_coordinates(xc,yc)
    #         lats[y,x] = lat
    #         lons[y,x] = lon
    #
    # latdif = lats-latitude
    # londif = lons-longitude
    # print(np.min(latdif),np.max(latdif),np.min(londif),np.max(londif))

    # print(np.min(latdif),np.max(latdif),np.min(londif),np.max(londif))


def do_resampled_vm_christmas():
    input_dir_base = '/store/COP2-OC-TAC/arc/'
    output_dir_base = '/store/COP2-OC-TAC/arc/resampled'
    unzip_path = '/store/COP2-OC-TAC/arc/unzip'
    from datetime import datetime as dt
    from datetime import timedelta
    date_ref = dt(2016, 5, 1)
    date_fin = dt(2016, 5, 31)
    while date_ref <= date_fin:
        print(f'[INFO]******************************************************************************->{date_ref}')
        input_dir = os.path.join(input_dir_base, date_ref.strftime('%Y%m%d'))
        do_resample = True
        if not os.path.exists(input_dir):
            print(f'[WARNING] Input directory {input_dir} is not available. Skiping...')
            do_resample = False
        output_dir_year = os.path.join(output_dir_base, date_ref.strftime('%Y'))
        output_dir_month = os.path.join(output_dir_year, date_ref.strftime('%m'))
        output_dir_day = os.path.join(output_dir_month,date_ref.strftime('%d'))
        date_ref_str = date_ref.strftime('%Y%m%d')
        output_name = f'{date_ref_str}_cmems_cnr_arc_rrs_resampled.nc'
        file_output = os.path.join(output_dir_day, output_name)
        if os.path.exists(file_output):
            print(f'[WARNING] Output file {file_output} already exist. Skiping...')
            do_resample = False
        if not do_resample:
            date_ref = date_ref + timedelta(hours=24)
            continue
        if not os.path.exists(output_dir_year):
            os.mkdir(output_dir_year)
        if not os.path.exists(output_dir_month):
            os.mkdir(output_dir_month)
        if not os.path.exists(output_dir_day):
            os.mkdir(output_dir_day)

        make_resample_dir(input_dir, output_dir_day, unzip_path, True, False)
        date_ref = date_ref + timedelta(hours=24)


def do_check5():
    print('do check 5')
    input_dir = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/15072018'
    fout = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/20180715_CNR_ARC_Resampled_v3.nc'
    from arc_integration import ArcIntegration
    arc_integration = ArcIntegration(None, args.verbose, input_dir)
    arc_integration.make_integration_avg(fout)
    # arc_integration.get_info()

    # for name in arc_integration.info:
    #     print(name,'-------------------------------------------')
    #     yini = arc_integration.info[name]['y_min']
    #     yfin = arc_integration.info[name]['y_max']
    #     xini = arc_integration.info[name]['x_min']
    #     xfin = arc_integration.info[name]['x_max']
    #     info_over = arc_integration.get_overlapping_images(name,yini,yfin,xini,xfin,False)
    #     print(name, '------------------------------------------->',len(info_over))
    #     # for nameo in info_over:
    #     #     print('        ',nameo)
    # # print(arc_integration.info)

    # print(arc_integration.time_min,arc_integration.time_max)
    # from datetime import datetime as dt
    # print(dt.fromtimestamp(arc_integration.time_min))
    # print(dt.fromtimestamp(arc_integration.time_max))
    # print(dt.fromtimestamp(arc_integration.time_min-1))
    # print((dt.fromtimestamp(arc_integration.time_max)-dt.fromtimestamp(arc_integration.time_min)).total_seconds()/3600)


def do_check44():
    file = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/inforesampling.txt'
    f1 = open(file, 'r')
    lines = []
    for line in f1:
        line = line.strip()
        if line.startswith('[INFO] Number of non-masked pixels:'):
            istart = line.find(':') + 1
            iend = line.find('(')
            nhere = line[istart:iend].strip()
            lines.append(nhere)
    f1.close()

    for idx in range(0, len(lines), 2):
        print(lines[idx])


def make_resample_dir(dirorig, dirdest, unzip_path, doresample, dokml):
    import simplekml
    # import netCDF4
    # import numpy as np
    import zipfile as zp
    # dirorig = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/15072018'
    # unzip_path = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/temp'
    if dirdest is None:
        dirdest = dirorig

    fconfig = None
    if args.config_file:
        fconfig = args.config_file
    from arc_mapinfo import ArcMapInfo
    ami = ArcMapInfo(fconfig, args.verbose)

    first_line = ['Source', 'StartDate', 'RelOrbit', 'GranuleId', 'SensorFlag', 'OrigWidth', 'OrigHeight', 'OrigNTotal',
                  'OrigNRRS', 'OrigNValid', 'OrigPValid', 'YMin', 'YMax', 'XMin', 'XMax', 'Width', 'Height', 'NTotal',
                  'NValid', 'PValid']
    lines_out = [';'.join(first_line)]

    rel_pass_dict = {}
    idref = -1
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
            print(f'File: {name} ({idx}/{nfiles})')
            idx = idx + 1
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

        if doresample:
            file_out = os.path.join(dirdest, f'{output_name}_resampled.nc')
            if os.path.exists(file_out):
                print(f'[INFO] File {file_out} already exist. Skyping resample...')
                if not dokml:
                    continue
            if do_zip_here:
                with zp.ZipFile(prod_path, 'r') as zprod:
                    if args.verbose:
                        print(f'[INFO] Unziping {name} to {unzip_path}')
                    zprod.extractall(path=unzip_path)

        from olci_l2 import OLCI_L2
        olimage = OLCI_L2(path_prod_u, args.verbose)
        olimage.get_geo_and_params()

        rel_pass = str(olimage.get_rel_pass())
        if rel_pass not in rel_pass_dict.keys():
            if idref == -1:
                idref = 0
            else:
                idref = idref + 100
            rel_pass_dict[rel_pass] = 0
            if dokml:
                if foutkml is not None and kml is not None:
                    kml.save(foutkml)
                    idcolor = idcolor + 1
                    if idcolor >= len(red):
                        idcolor = 0
                foutkml = os.path.join(dirdest, f'Passes_RelativeOrbit_{rel_pass}.kml')
                kml = simplekml.Kml()

        else:
            rel_pass_dict[rel_pass] = rel_pass_dict[rel_pass] + 1
        sensor_id = math.pow(2, rel_pass_dict[rel_pass]) + idref

        if doresample:
            # file_out = os.path.join(dirdest, f'{output_name}_resampled.nc')
            # if os.path.exists(file_out):
            #     print(f'[INFO] File {file_out} already exist. Skyping...')
            #     continue
            line_out = ami.make_resample_impl(olimage, file_out, sensor_id, idref)
            lines_out.append(line_out)
            if zp.is_zipfile(prod_path):
                #os.remove(prod_path)
                for fn in os.listdir(path_prod_u):
                    os.remove(os.path.join(path_prod_u, fn))
               

        if dokml:
            start_date = olimage.get_start_date()
            start_date_s = 'UKNOWN START DATE'
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

    if doresample:
        fcsvout = os.path.join(dirdest, 'ResampleInfo.csv')
        if args.verbose:
            print(f'[INFO] Creating info file:{fcsvout}')
        fw = open(fcsvout, 'w')
        for line in lines_out:
            fw.write(line)
            fw.write('\n')
        fw.close()


def do_check3():
    import netCDF4

    file = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/20220713_cmems_obs-oc_arc_bgc-reflectance_nrt_l3-olci-300m_P1D.nc'
    dataset = netCDF4.Dataset(file)
    import numpy as np
    var = dataset.variables['RRS510']
    filters = var.filters()
    print(filters)


def do_check2():
    import netCDF4
    # import numpy as np
    # import zipfile as zp
    print('-----------------------')
    print('ORIGINAL')
    file = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/PML_EXAMPLES/20220713_cmems_obs-oc_arc_bgc-reflectance_nrt_l3-olci-300m_P1D.nc'
    dataset = netCDF4.Dataset(file)
    import numpy as np
    # for namevar in dataset.variables:
    #     if namevar.startswith('RRS'):
    #         var = dataset.variables[namevar]
    #         ny = var.shape[1]
    #         nx = var.shape[2]
    #         total = ny*nx
    # nvalid = 0
    # ystep = 1200
    # xstep = nx
    #
    # for y in range(0, ny, ystep):
    #     for x in range(0,nx,xstep):
    #         #print(y,x)
    #         yini = y
    #         yfin = y + ystep
    #         if yfin > ny:
    #             yfin = ny
    #         xini = 0
    #         xfin = nx
    #         rrsarray = np.ma.array(var[0, yini:yfin, xini:xfin])
    #         nvalid = nvalid + rrsarray.count()

    # porc = (nvalid/total)*100
    # line = f'{namevar};{var.shape[1]};{var.shape[2]};{total};{nvalid};{porc}'
    # print(line)
    dataset.close()

    print('----------------------')
    print('RESAMPLED')
    file = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/PML_EXAMPLES/20220713_cmems_obs-oc_arc_bgc-reflectance_nrt_l3-olci-300m_P1D_resampled.nc'
    dataset = netCDF4.Dataset(file)
    for namevar in dataset.variables:
        if namevar.startswith('RRS'):
            var = dataset.variables[namevar]
            ny = var.shape[0]
            nx = var.shape[1]
            total = ny * nx
            nvalid = 0
            ystep = 5000
            xstep = nx

            for y in range(0, ny, ystep):
                for x in range(0, nx, xstep):
                    yini = y
                    yfin = y + ystep
                    if yfin > ny:
                        yfin = ny
                    xini = 0
                    xfin = nx
                    rrsarray = np.ma.array(var[yini:yfin, xini:xfin])
                    nvalid = nvalid + rrsarray.count()

            porc = (nvalid / total) * 100
            line = f'{namevar};{var.shape[0]};{var.shape[1]};{total};{nvalid};{porc}'
            print(line)

    dataset.close()

    # lat_array = np.array(dataset.variables['latitude'])
    # lon_array = np.array(dataset.variables['longitude'])
    # from geopy import distance
    # lines = []
    # for y in range(len(lat_array)-2,1,-50):
    #     print(y)
    #     #dist = []
    #     for x in range(0,len(lon_array)-1,1000):
    #         coords_1 = (lat_array[y],lon_array[x])
    #         coords_2 = (lat_array[y],lon_array[x+1])
    #         d1 = distance.distance(coords_1,coords_2).meters
    #         coords_1 = (lat_array[y], lon_array[x])
    #         coords_2 = (lat_array[y], lon_array[x - 1])
    #         d2 = distance.distance(coords_1, coords_2).meters
    #         coords_1 = (lat_array[y-1], lon_array[x])
    #         coords_2 = (lat_array[y], lon_array[x])
    #         d3 = distance.distance(coords_1, coords_2).meters
    #         coords_1 = (lat_array[y+1], lon_array[x])
    #         coords_2 = (lat_array[y], lon_array[x])
    #         d4 = distance.distance(coords_1, coords_2).meters
    #         l = [d1,d2,d3,d4]
    #         ly = [d3,d4]
    #         lx = [d1,d2]
    #         m = np.mean(np.array(l))
    #         my = np.mean(np.array(ly))
    #         mx = np.mean(np.array(lx))
    #         #dist.append(m)
    #         line = f'{y};{x};{lat_array[y]};{lon_array[x]};{my};{mx};{m}'
    #         lines.append(line)
    # fout = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/GridPML.csv'
    # with open(fout, 'w') as f:
    #     for line in lines:
    #         f.write(line)
    #         f.write('\n')
    # y = 7573
    # x = 3600
    # coords_1 = (lat_array[y], lon_array[x])
    # coords_2 = (lat_array[y], lon_array[x + 1])
    # d1 = distance.distance(coords_1, coords_2).meters
    # coords_1 = (lat_array[y], lon_array[x])
    # coords_2 = (lat_array[y], lon_array[x - 1])
    # d2 = distance.distance(coords_1, coords_2).meters
    # coords_1 = (lat_array[y - 1], lon_array[x])
    # coords_2 = (lat_array[y], lon_array[x])
    # d3 = distance.distance(coords_1, coords_2).meters
    # coords_1 = (lat_array[y+1], lon_array[x])
    # coords_2 = (lat_array[y], lon_array[x])
    # d4 = distance.distance(coords_1, coords_2).meters
    # l = [d1,d2,d3,d4]
    # ly = [d3,d4]
    # lx = [d1,d2]
    # m = np.mean(np.array(l))
    # my = np.mean(np.array(ly))
    # mx = np.mean(np.array(lx))
    # line = f'{y};{x};{lat_array[y]};{lon_array[x]};{my};{mx};{m}'
    # print(line)


def do_check():
    from olci_l2 import OLCI_L2
    folci = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/S3A_OL_2_WFR____20180715T110812_20180715T111112_20211121T193318_0180_033_251______MAR_R_NT_003.SEN3'
    olimage = OLCI_L2(folci, args.verbose)
    array = olimage.get_observation_angle_array()
    print(array.shape)


if __name__ == '__main__':
    main()
