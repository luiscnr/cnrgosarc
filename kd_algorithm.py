import numpy as np


class KD_ALGORITHMS():

    def __init__(self, kdalgorithm):
        self.kdalgorithm = kdalgorithm

        self.fillValue = -999.0

        self.qratio = 1.0

    def compute_kd(self, *args):
        if self.kdalgorithm == 'OK2-560':
            return self.compute_kd490_ok2_560(args[0], args[1])

        return None

    def check_kd(self, *args):
        if self.kdalgorithm == 'OK2-560':
            return self.check_kd490_ok2_560(self, args[0], args[1])

    def check_kd490_ok2_560(self, rrs490, rrs560):
        rrs490[rrs490 < 0] = self.fillValue
        rrs560[rrs560 < 0] = self.fillValue
        valid = np.logical_and(rrs490 != self.fillValue, rrs560 != self.fillValue)
        return np.count_nonzero(valid)

    def compute_kd490_ok2_560(self, rrs490, rrs560):
        ##Morel et al., 2017 Examining the consistency of products derived from various ocean color sensors in open ocean (Case 1) waters in the perspective of a multi-sensor approach
        ##Remote Sens. Environ 111(1), 69-88

        rrs490[rrs490 < 0] = self.fillValue
        rrs560[rrs560 < 0] = self.fillValue
        valid = np.logical_and(rrs490 != self.fillValue, rrs560 != self.fillValue)

        coeffs = [-0.8278866, -1.642189, 0.90261, -1.626853, 0.0885039]


        log_ratio = np.zeros(rrs490.shape)
        log_ratio[:] = self.fillValue

        log_ratio[valid] = rrs490[valid]/ rrs560[valid]
        log_ratio[valid] = log_ratio[valid]/self.qratio
        log_ratio[valid] = np.log10(log_ratio[valid])
        log_ratio[valid] = log_ratio[valid]


        res = np.zeros(rrs490.shape)
        res[valid == False] = self.fillValue
        for x in range(len(coeffs)):
            res[valid] = res[valid] + (coeffs[x] * np.power(log_ratio[valid], x))

        res[valid] = np.power(10,res[valid])+0.0166

        return res
