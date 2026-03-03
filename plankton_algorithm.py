import os,configparser,__init__
import numpy as np
from arc_options import ARC_OPTIONS


class PLANKTON_ALGORITHMS:

    def __init__(self):
        self.functions = {
            'PICO':{
                'type':'poly3',
                'a': 0.0196,
                'b': 0.0074,
                'c': -0.1900,
                'd': 0.1800
            },
            'NANO':{
                'type': 'exp_af',
                'a': 0.1540,
                'b': -1.2400,
                'c': 0.3980,
                'd': 0.3650,
                'e': 0.0266,
                'f': 0.7190
            },
            'MICRO':{
                'type': 'MICRO',
                'value': "1-NANO-PICO"
            },
            'DINO':{
                'type': 'poly2',
                'a': -0.0215,
                'b': -0.0145,
                'c': 0.0528,
            },
            'DIATO':{
                'type': 'DIATO',
                'value': 'MICRO - DINO',
            },
            'CRYPTO':{
                'type': 'exp_af',
                'a': 0.0780,
                'b': -4.0000,
                'c': 4.0000,
                'd': 0.0908,
                'e': -0.0709,
                'f': 0.5260
            },
            'GREEN':{
                'type': 'exp_ac',
                'a': -0.6380,
                'b': 1.8900,
                'c': 9.1400,
            },
            'PROKAR':{
                'type':'poly3',
                'a': -0.0033,
                'b': -0.0172,
                'c': -0.0451,
                'd': 0.0304
            },
            'HAPT':{
                'type': 'HAPT',
                'value': '1 - MICRO - CRYPTO - GREEN - PROKAR'
            }
        }
        self.min_chla = 0.010
        self.max_chla = 4.580
        self.file_attrs = os.path.join(os.path.dirname(__init__.__file__),'PSC_PFT_Attrs.ini')

    def get_attrs(self):
        if not os.path.exists(self.file_attrs):
            print(f'[WARNING] Attribute file for PSC/PFT {self.file_attrs} is not available')
            return None
        try:
            options = configparser.ConfigParser()
            options.read(self.file_attrs)
        except:
            print(f'[ERROR] Config file {self.file_attrs} could not be read. Exiting...')
            return

        arc_opt = ARC_OPTIONS(options)
        attrs = {x:arc_opt.get_str_options_from_section_as_dict(x) for x in self.functions}
        return attrs

    def compute_poly3(self,array,a,b,c,d,function=''):
        try:
            res = (a*np.power(array,3))+(b*np.power(array,2))+(c*array)+d
        except Exception as ex:
            print(f'[ERROR] Error computing function {function}. Exception: {ex}')
            return None
        return res
    def compute_poly3_dict(self,array,dvalues,function=''):
        return self.compute_poly3(array,dvalues['a'],dvalues['b'],dvalues['c'],dvalues['d'],function=function)

    def compute_exp_af(self,array,a,b,c,d,e,f,function=''):
        try:
            res = (a * (np.exp((-1)*np.power((array-b)/c,2)))) +  (d * (np.exp((-1)*np.power((array-e)/f,2))))
        except Exception as ex:
            print(f'[ERROR] Error computing function {function}. Exception: {ex}')
            return None
        return res
    def compute_exp_af_dict(self, array,dvalues,function=''):
        return self.compute_exp_af(array,dvalues['a'],dvalues['b'],dvalues['c'],dvalues['d'],dvalues['e'],dvalues['f'],function=function)

    def compute_poly2(self,array,a,b,c,function=''):
        try:
            res = (a*np.power(array,2))+(b*array)+(c)
        except Exception as ex:
            print(f'[ERROR] Error computing function {function}. Exception: {ex}')
            return None
        return res
    def compute_poly2_dict(self,array,dvalues,function=''):
        return self.compute_poly2(array,dvalues['a'],dvalues['b'],dvalues['c'],function=function)

    #[exp(a x + b) + c x]⁻¹
    def compute_exp_ac(self,array,a,b,c,function=''):
        try:
            res = np.power(np.exp((a*array)+b)+(c * array),-1)
        except Exception as ex:
            print(f'[ERROR] Error computing function {function}. Exception: {ex}')
            return None
        return res
    def compute_exp_ac_dict(self,array,dvalues,function=''):
        return self.compute_exp_ac(array,dvalues['a'],dvalues['b'],dvalues['c'],function=function)

    def compute_micro(self,nano,pico):
        try:
            return 1 - nano - pico
        except Exception as ex:
            print(f'[ERROR] Error computing MICRO. Exception: {ex}')
            return None
    def compute_micro_dict(self,res):
        if res is None:
            print(f'[ERROR] MICRO can not be computed. Argument is None, but it should be a dictionary including NANO and PICO arrays')
            return None
        if not isinstance(res,dict):
            print(f'[ERROR] MICRO can not be computed. Argument should be a dictionary including NANO and PICO arrays, but it is a {type(res)} object')
            return None
        if not 'NANO' in res:
            print(f'NANO is not available in the dictionary argument for computing MICRO')
            return None
        if not 'PICO' in res:
            print(f'PICO is not available in the dictionary argument for computing MICRO')
            return None
        return self.compute_micro(res['NANO'],res['PICO'])

    def compute_diato(self,micro,dino):
        try:
            return micro-dino
        except Exception as ex:
            print(f'[ERROR] Error computing DIATO. Exception: {ex}')
            return None

    def compute_diato_dict(self,res):
        if res is None:
            print(f'[ERROR] DIATO can not be computed. Argument is None, but it should be a dictionary including MICRO and DINO arrays')
            return None
        if not isinstance(res,dict):
            print(f'[ERROR] DIATO can not be computed. Argument should be a dictionary including MICRO and DINO arrays, but it is a {type(res)} object')
            return None
        if not 'MICRO' in res:
            print(f'MICRO is not available in the dictionary argument for computing DIATO')
            return None
        if not 'DINO' in res:
            print(f'DINO is not available in the dictionary argument for computing DIATO')
            return None
        return self.compute_diato(res['MICRO'],res['DINO'])

    def compute_hapt(self,micro,crypto,green,prokar):
        try:
            return 1 - micro - crypto - green - prokar
        except Exception as ex:
            print(f'[ERROR] Error computing HAPT. Exception: {ex}')
            return None

    def compute_hapt_dict(self,res):
        if res is None:
            print(f'[ERROR] HAPT can not be computed. Argument is None, but is should be a dictionary including MICRO, CRYPTO, GREEN and PROKAR arrays')
            return None
        if not isinstance(res,dict):
            print(f'[ERROR] HAPT can not be computed. Argument should be a dictionary including MICRO, CRYPTO, GREEN and PROKAR arrays, but it is a {type(res)} object')
            return None
        if not 'MICRO' in res:
            print(f'MICRO is not available in the dictionary argument for computing HAPT')
            return None
        if not 'CRYPTO' in res:
            print(f'CRYPTO is not available in the dictionary argument for computing HAPT')
            return None
        if not 'GREEN' in res:
            print(f'GREEN is not available in the dictionary argument for computing HAPT')
            return None
        if not 'PROKAR' in res:
            print(f'PROKAR is not available in the dictionary argument for computing HAPT')
            return None
        return self.compute_hapt(res['MICRO'],res['CRYPTO'],res['GREEN'],res['PROKAR'])

    ##result is a dictionary to save results, if it is None, array is directly returned. Return ONLY FRACTIONS
    def compute_function(self,function,chl,make_log_chl,result):
        if make_log_chl:
            chl = np.log10(chl)

        if not function in self.functions:
            return None
        type_f = self.functions[function]['type']

        if type_f=='poly3':
            array  = self.compute_poly3_dict(chl,self.functions[function],function=function)
        elif type_f=='poly2':
            array = self.compute_poly2_dict(chl, self.functions[function],function=function)
        elif type_f=='exp_af':
            array = self.compute_exp_af_dict(chl, self.functions[function],function=function)
        elif type_f=='exp_ac':
            array = self.compute_exp_ac_dict(chl,self.functions[function],function=function)
        elif type_f=='MICRO':
            array = self.compute_micro_dict(result)
        elif type_f=='DIATO':
            array = self.compute_diato_dict(result)
        elif type_f=='HAPT':
            array = self.compute_hapt_dict(result)
        else:
            print(f'[ERROR] Function type {funtion} is not defined.')
            array = None

        if result is None:
            return array
        elif isinstance(result,dict):
            result[function] = array
            return result
        else:
            print(f'[WARNING] result should be a dict or None, returning the array')
            return array

    ##implementation for all the function based on an array of log10 chl values.
    def compute_all_functions_impl(self,array_log_chl,output_as_conc,transform_array,array_no_log_chl = None):
        if array_no_log_chl is None and output_as_conc:
            array_no_log_chl = np.power(10,array_log_chl)
        result = {f:None for f in self.functions}

        ##compute functions
        for f in self.functions:
            print(f'[INFO] Computing function {f} using a chl array with shape {array_log_chl.shape}')
            result = self.compute_function(f,array_log_chl,False,result)

        ##output as concentration and shape transformation
        for f in self.functions:
            str_info = f'[INFO] {f}: '
            if output_as_conc:
                array_c = np.multiply(result[f],array_no_log_chl)
                result[f] = array_c
                str_info = f'{str_info}Fractions converted to concentrations. '
            if transform_array and transform_array['shape_orig']!=result[f].shape:
                array_t = np.ma.masked_all(transform_array['shape_orig'])
                array_t[transform_array['indices_valid']] = result[f]
                result[f] = array_t
                str_info = f'{str_info}Array transformed to original shape:{transform_array["shape_orig"]}. '
            print(str_info)



        return result

    ##compute all function for masked array of any dimension
    def compute_all_functions(self,chl,log_input_chl,filter_min_max,output_as_conc):
        chl_valid, chl_valid_no_log, transform_array = self.get_input_valid_array(chl,log_input_chl,filter_min_max,output_as_conc)
        if chl_valid is None:
            return
        result = self.compute_all_functions_impl(chl_valid,output_as_conc,transform_array,array_no_log_chl=chl_valid_no_log)
        return result

    ##retrieve input valid arrays
    def get_input_valid_array(self,chl,log_input_chl,filter_min_max,output_as_conc):
        if not isinstance(chl, np.ndarray):
            chl = np.ma.array([chl])
        chl = np.ma.array(chl)  ##work with masked arrays
        shape_orig = chl.shape

        if not filter_min_max:
            indices_valid = np.where(np.logical_and(chl.mask == False, chl > 0))
        else:
            indices_valid = np.where(
                np.logical_and(chl.mask == False, np.logical_and(chl > self.min_chla, chl < self.max_chla)))

        if len(indices_valid[0])==0:
            print(f'[WARNING] No valid chl-a data is available within the applicability range ({self.min_chla}-{self.max_chla}) to compute PSC-PFT')
            return [None]*3

        chl_valid = chl[indices_valid]
        if len(indices_valid[0]) == 1 and np.ma.isarray(chl_valid) == False:  ##make sure chl_valid is a masked array
            chl_valid = np.ma.array([chl_valid])

        chl_valid_no_log = None
        if log_input_chl:
            if output_as_conc:
                chl_valid_no_log = chl_valid.copy()
            chl_valid = np.ma.log10(chl_valid)
        else:
            if output_as_conc:
                chl_valid_no_log = np.power(10, chl_valid)

        if len(indices_valid[0])==chl.size:
            transform_array = None
        else:
            transform_array = {
                'shape_orig':shape_orig,
                'indices_valid':indices_valid
            }

        return chl_valid,chl_valid_no_log,transform_array




