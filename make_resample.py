from arc_mapinfo import ArcMapInfo
from olci_l2 import OLCI_L2
def main():
    print('Started make_resample')
    ami = ArcMapInfo(None)
    #fgrid = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/ArcGridOlci.nc'
    #fout = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/ArcGridOlciQuickLook.png'
    #ami.create_nc_filegrid('/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/ArcGridOlci.nc',True)
    #ami.save_quick_look_fgrid(fout,fgrid)

    folci = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/S3A_OL_2_WFR____20180715T110812_20180715T111112_20211121T193318_0180_033_251______MAR_R_NT_003.SEN3'
    fout = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/test.nc'
    olimage = OLCI_L2(folci)
    ami.make_resample_impl(olimage,fout)

    # fdata = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/test.nc'
    # fout = '/mnt/c/DATA_LUIS/OCTAC_WORK/ARC_TEST/test.png'
    # ami.save_quick_look_fdata(fout,fdata)

if __name__ == '__main__':
    main()
