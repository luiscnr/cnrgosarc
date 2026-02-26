import numpy as np
from scipy.interpolate import interp1d


class KD_ALGORITHMS:

    def __init__(self, kdalgorithm):
        self.kdalgorithm = kdalgorithm

        self.fillValue = -999.0

        self.q0ratio_default = 1.0

        q0_chla = np.array([0.03, 0.1, 0.3, 1.0, 3.0, 10.0])
        q0_443 = np.array([3.29125, 3.3857, 3.52968, 3.74683, 3.99138, 4.31326])
        q0_490 = np.array([3.24564, 3.40843, 3.61341, 3.86893, 4.11003, 4.39395])
        q0_510 = np.array([3.1765, 3.35968, 3.60166, 3.89409, 4.14114, 4.40546])
        q0_560 = np.array([3.13802, 3.33689, 3.60695, 3.94203, 4.18865, 4.36813])

        self.interp443 = interp1d(q0_chla, q0_443, kind='linear', fill_value=(0.03, 10.0), bounds_error=False)
        self.interp490 = interp1d(q0_chla, q0_490, kind='linear', fill_value=(0.03, 10.0), bounds_error=False)
        self.interp510 = interp1d(q0_chla, q0_510, kind='linear', fill_value=(0.03, 10.0), bounds_error=False)
        self.interp560 = interp1d(q0_chla, q0_560, kind='linear', fill_value=(0.03, 10.0), bounds_error=False)
        # self.chla_range = [0.03, 10.0]

    def get_q0_ratio(self, chla, type):
        if type == '443/560':
            q0_ratio = self.interp443(chla) / self.interp560(chla)
        if type == '490/560':
            q0_ratio = self.interp490(chla) / self.interp560(chla)
        if type == '510/560':
            q0_ratio = self.interp510(chla) / self.interp560(chla)


        return q0_ratio

    def compute_kd_param(self,rrs443, rrs490, rrs510, rrs560):
        array443 = np.array([rrs443])
        array490 = np.array([rrs490])
        array510 = np.array([rrs510])
        array560 = np.array([rrs560])



        chla = self.compute_chla_ocme4(array443,array490,array510,array560,None)
        if chla is None:
            return -999.0, -999.0, -999.0
        chla2 = self.compute_chla_ocme4(array443, array490, array510, array560, chla)
        q0,kd = self.compute_kd490_ok2_560(array490,array560,chla2)
        return q0[0],kd[0], chla[0], chla2[0]

    def compute_kd(self, *args):
        if self.kdalgorithm == 'OK2-560':
            if len(args) == 2:
                return self.compute_kd490_ok2_560(args[0], args[1], None)
            elif len(args) == 3:
                return self.compute_kd490_ok2_560(args[0], args[1], args[2])
        return None

    def check_kd(self, *args):
        if self.kdalgorithm == 'OK2-560':
            return self.check_kd490_ok2_560(self, args[0], args[1])

    def check_kd490_ok2_560(self, rrs490, rrs560):
        rrs490[rrs490 < 0] = self.fillValue
        rrs560[rrs560 < 0] = self.fillValue
        valid = np.logical_and(rrs490 != self.fillValue, rrs560 != self.fillValue)
        return np.count_nonzero(valid)

    def compute_chla_ocme4(self, rrs443, rrs490, rrs510, rrs560, chla):
        #print(rrs443.shape)

        rrs443[rrs443 <= 0] = self.fillValue
        rrs490[rrs490 <= 0] = self.fillValue
        rrs510[rrs510 <= 0] = self.fillValue
        rrs560[rrs560 <= 0] = self.fillValue
        #print(rrs443.shape)
        # if len(rrs443[rrs443==self.fillValue])==len(rrs443):
        #     return None
        # if len(rrs490[rrs490==self.fillValue])==len(rrs490):
        #     return None
        # if len(rrs510[rrs510==self.fillValue])==len(rrs510):
        #     return None
        # if len(rrs560[rrs560==self.fillValue])==len(rrs560):
        #     return None

        sh = rrs490.shape
        # valid1 = np.logical_and(rrs443 != self.fillValue, rrs560 != self.fillValue)
        # valid2 = np.logical_and(rrs490 != self.fillValue, rrs510 != self.fillValue)
        # valid = np.logical_and(valid1==True, valid2==True)

        # 443-560
        valid = np.logical_and(rrs443 != self.fillValue, rrs560 != self.fillValue)
        log_ratio1 = np.zeros(sh)
        q0_ratio = np.zeros(sh)
        log_ratio1[:] = self.fillValue
        q0_ratio[:] = self.q0ratio_default
        if chla is not None:
            q0_ratio[valid] = self.get_q0_ratio(chla[valid], '443/560')
        log_ratio1[valid] = rrs443[valid] / rrs560[valid]
        log_ratio1[valid] = log_ratio1[valid] / q0_ratio[valid]

        # 490-560
        valid = np.logical_and(rrs490 != self.fillValue, rrs560 != self.fillValue)
        log_ratio2 = np.zeros(sh)
        q0_ratio = np.zeros(sh)
        log_ratio2[:] = self.fillValue
        q0_ratio[:] = self.q0ratio_default
        if chla is not None:
            q0_ratio[valid] = self.get_q0_ratio(chla[valid], '490/560')
        log_ratio2[valid] = rrs490[valid] / rrs560[valid]
        log_ratio2[valid] = log_ratio2[valid] / q0_ratio[valid]

        # 510-560
        valid = np.logical_and(rrs510 != self.fillValue, rrs560 != self.fillValue)
        log_ratio3 = np.zeros(sh)
        q0_ratio = np.zeros(sh)
        log_ratio3[:] = self.fillValue
        q0_ratio[:] = self.q0ratio_default
        if chla is not None:
            q0_ratio[valid] = self.get_q0_ratio(chla[valid], '510/560')
        log_ratio3[valid] = rrs510[valid] / rrs560[valid]
        log_ratio3[valid] = log_ratio3[valid] / q0_ratio[valid]

        log_ratio = np.maximum(log_ratio1, np.maximum(log_ratio2, log_ratio3))
        log_ratio[log_ratio != self.fillValue] = np.log10(log_ratio[log_ratio != self.fillValue])

        coeffs = [0.4502748, -3.259491, 3.522731, -3.359422, 0.949586]
        res = np.zeros(sh)
        res[log_ratio == self.fillValue] = self.fillValue
        for x in range(len(coeffs)):
            res[valid] = res[valid] + (coeffs[x] * np.power(log_ratio[valid], x))

        res[valid] = np.power(10, res[valid])

        return res

    def compute_kd490_ok2_560(self,rrs490, rrs560, chla):
        ##Morel et al., 2017 Examining the consistency of products derived from various ocean color sensors in open ocean (Case 1) waters in the perspective of a multi-sensor approach
        ##Remote Sens. Environ 111(1), 69-88

        rrs490[rrs490 <= 0] = self.fillValue
        rrs560[rrs560 <= 0] = self.fillValue
        valid = np.logical_and(rrs490 != self.fillValue, rrs560 != self.fillValue)


        # if chla is not None:
        #     valid = np.logical_and(valid == True, chla >= self.chla_range[0])

        coeffs = [-0.8278866, -1.642189, 0.90261, -1.626853, 0.0885039]

        log_ratio = np.zeros(rrs490.shape)
        q0_ratio = np.zeros(rrs490.shape)
        log_ratio[:] = self.fillValue
        q0_ratio[:] = self.q0ratio_default
        if chla is not None:
            q0_ratio[valid] = self.get_q0_ratio(chla[valid], '490/560')
        #print('We are here: ',log_ratio.shape,valid.shape,)
        log_ratio[valid] = np.divide(rrs490[valid],rrs560[valid])
        #print('---->', rrs490[valid],rrs560[valid],log_ratio[valid], q0_ratio[valid])
        log_ratio[valid] = log_ratio[valid] / q0_ratio[valid]
        log_ratio[valid] = np.log10(log_ratio[valid])
        #log_ratio[valid] = log_ratio[valid]

        res = np.zeros(rrs490.shape)
        res[valid == False] = self.fillValue
        for x in range(len(coeffs)):
            res[valid] = res[valid] + (coeffs[x] * np.power(log_ratio[valid], x))
            #print(x,log_ratio[valid],res[valid])

        #print(res[valid])

        res[valid] = np.power(10, res[valid]) + 0.0166



        return q0_ratio,res
