

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

    def compute_poly3(self,array,a,b,c,d):
        res = (a*np.power(array,3))+(b*np.power(array,2))+(c*array)+d
        return res
    def compute_poly3_dict(self,array,dvalues):
        return self.compute_poly3(array,dvalues['a'],dvalues['b'],dvalues['c'],dvalues['d'])

    def compute_exp_af(self,array,a,b,c,d,e,f):
        res = (a * (np.exp((-1)*np.power((array-b)/c,2)))) +  (d * (np.exp((-1)*np.power((array-e)/f,2))))
        return res

    def compute_exp_af_dict(self, array,dvalues):
        return self.compute_exp_af(array,dvalues['a'],dvalues['b'],dvalues['c'],dvalues['d'],dvalues['e'],dvalues['f'])
