import os
from netCDF4 import Dataset


class OLCI_L2():
    def __init__(self, path_source):
        self.path_source = path_source
        self.wl_list = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75, 753.75, 778.75, 865, 885,
                        1020.5]
        self.reflectance_bands, self.nbands = self.get_reflectance_bands_info()
        self.max_dif_wl = 1

    def get_reflectance_bands_info(self):
        reflectance_bands = {}
        bands_ints = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18, 21]

        nbands = len(bands_ints)
        for index in range(nbands):
            v = bands_ints[index]
            wl = self.wl_list[index]
            band_name = f'Oa{v:02d}_reflectance'
            band_path = f'{band_name}.nc'
            file_path = os.path.join(self.path_source, band_path)
            reflectance_bands[band_name] = {
                'wavelenght': wl,
                'file_path': file_path
            }
        return reflectance_bands, nbands

    def get_lat_long_arrays(self):
        coordinates_filename = 'geo_coordinates.nc'
        filepah = os.path.join(self.path_source, coordinates_filename)
        nc_sat = Dataset(filepah, 'r')
        lat = nc_sat.variables['latitude'][:, :]
        lon = nc_sat.variables['longitude'][:, :]

        return lat, lon

    def get_reflectance_band_array(self,wlref):
        for band_name in self.reflectance_bands:
            dif = abs(wlref-self.reflectance_bands[band_name]['wavelenght'])
            if dif<self.max_dif_wl:
                nc_sat = Dataset(self.reflectance_bands[band_name]['file_path'],'r')
                array_reflectance = nc_sat.variables[band_name][:,:]
                return array_reflectance
        return None


