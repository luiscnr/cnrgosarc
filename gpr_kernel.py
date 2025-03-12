import numpy as np


class GPR_KERNEL:
    def __init__(self, kernel_dict):
        self.kernel_type = kernel_dict['Name']
        self.alpharq = None
        self.sigmaf = None
        self.sigmal = None
        self.length_scale = None
        if self.kernel_type == 'ARDRationalQuadratic':
            self.start_kernel_ARDRationalQuadratic(kernel_dict)
        if self.kernel_type == 'RationalQuadratic':
            self.start_kernel_RationalQuadratic(kernel_dict)

        # print('Kernel alpha:',self.alpharq)
        # print('Kernel sigma F :',self.sigmaf)
        # print('Kernel lscale :', self.length_scale)

    def start_kernel_RationalQuadratic(self, kernel_dict):
        self.alpharq = kernel_dict['AlphaRQ']
        self.sigmaf = kernel_dict['SigmaF']
        self.sigmal = kernel_dict['SigmaL']

    def start_kernel_ARDRationalQuadratic(self, kernel_dict):
        self.alpharq = kernel_dict['AlphaRQ']
        self.sigmaf = kernel_dict['SigmaF']
        predictors_added = False
        index_predictor = 1
        self.length_scale = []
        while not predictors_added:
            key = f'LengthScale{index_predictor}'
            if key in kernel_dict:
                self.length_scale.append(kernel_dict[key])
            else:
                predictors_added = True
            index_predictor = index_predictor + 1

    def compute_kernel(self, vector1, vector2):
        if self.kernel_type == 'ARDRationalQuadratic':
            return self.compute_ard_rational_quadratic(vector1, vector2)
        if self.kernel_type == 'RationalQuadratic':
            return self.compute_rational_quadratic(vector1, vector2)

    def compute_rational_quadratic(self,vector1,vector2):
        if len(vector1) != len(vector2):
            print('[WARNING] Kernel can not be computed, both vectors must have the same length')
            return np.NaN
        result = (self.sigmaf ** 2) * ((1 + (np.sum((vector1-vector2)**2)/(2*self.alpharq*(self.sigmal**2))))**-self.alpharq)
        return result




    def compute_ard_rational_quadratic(self, vector1, vector2):
        if len(vector1) != len(vector2):
            print('[WARNING] Kernel can not be computed, both vectors must have the same length')
            return np.NaN
        if len(vector1) != len(self.length_scale):
            print(
                f'[WARNING] ARD Rational Quadratic Kernel can not be computed, input vectors must have a length = {len(self.length_scale)}')
            return np.NaN
        # IN TWO LINES
        sum = np.sum(np.power((vector1 - vector2), 2) / np.power(self.length_scale, 2))
        result = (self.sigmaf ** 2) * ((1 + ((1 / (2 * self.alpharq)) * sum)) ** -self.alpharq)

        ##to check:
        # https://it.mathworks.com/help/stats/kernel-covariance-function-options.html
        # STEP BY STEP
        # sum = 0
        # nvec = len(vector1)
        # for idx in range(nvec):
        #     sum = sum + ((vector1[idx] - vector2[idx]) ** 2) / (self.length_scale[idx] ** 2)
        # coef = 1 / (2 * self.alpharq)
        # val = 1 + (coef * sum)
        # valp = val ** -self.alpharq
        # result = (self.sigmaf ** 2) * valp

        return result
