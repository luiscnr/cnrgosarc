import numpy as np
import os

class ARC_GPR_MODEL():

    def __init__(self, fmodel):
        self.use_ciao = False
        if os.path.basename(fmodel).startswith('CIAO_'):
            print('[INFO] Started GPR model - CIAO Algorithm')
            self.use_ciao = True
        else:
            print('[INFO] Started GPR model - SeaSARC Algorithm')
        import json
        f = open(fmodel)
        model_dict = json.load(f)
        f.close()
        kernel_dict = model_dict['Kernel']
        from gpr_kernel import GPR_KERNEL
        self.kernel = GPR_KERNEL(kernel_dict)
        self.active_set_vectors = np.array(model_dict['ActiveSetVectors'])
        self.alpha = np.array(model_dict['Alpha'])
        self.nactive_set_vectors = len(self.active_set_vectors)
        self.npredictors = len(self.active_set_vectors[0])

        self.sigma = model_dict['Sigma']
        self.beta = np.array(model_dict['Beta'])

        # print(model_dict.keys())
        # print(kernel_dict)
        # print(type(self.active_set_vectors),self.active_set_vectors.shape)

        # for idx in range(self.nactive_set_vectors):
        #     active_vector = self.active_set_vectors[idx]
        #     print(active_vector)
        # print(model_dict['BasisFunction'])
        # print(model_dict['Beta'])
        # print(model_dict['Sigma'])

    ##Basic implmentation, with log-transformed RRS data, valid for SeaSARC and CIAO
    def compute_chla_impl(self, feature_vector):
        # term 1->linear equation
        X = np.concatenate(([1], feature_vector))

        H = X @ self.beta


        # term 2->kernel
        kresult = np.zeros(self.nactive_set_vectors)
        for idx in range(self.nactive_set_vectors):
            active_vector = self.active_set_vectors[idx]
            kresult[idx] = self.kernel.compute_kernel(active_vector, feature_vector)

        K = kresult @ self.alpha

        val_fin = H + K

        return val_fin

    ##Basic implementation valid for SeaSARC and CIAO. RRS could log-transormed or no transformed.
    def compute_chla(self, feature_vector, transform_input):
        if len(feature_vector) != self.npredictors:
            return np.nan
        if transform_input:
            if self.use_ciao:
                feature_vector[1:] = np.log10(feature_vector[1:])
            else:
                feature_vector[2:] = np.log10(feature_vector[2:])
            # for idx in range(2, self.npredictors):
            #     feature_vector[idx] = np.log10(feature_vector[idx])
        val_fin = self.compute_chla_impl(feature_vector)
        chla = 10 ** val_fin
        return chla
    ##SeaSARC#######################################################################################################
    ##compute chla from params (no log-transformed rrs) for SeaSARC algorithm
    def compute_chla_from_param(self, long_value, day, val_443, val_490, val_560, val_665):
        feature_vector = [long_value, day, val_443, val_490, val_560, val_665]
        return self.compute_chla(feature_vector, True)

    ##fast implementantion, with transformed data, for SeaSARC algorithm
    def compute_chla_from_matrix(self, matrix):
        if len(matrix.shape) != 2 or matrix.shape[1] != (self.npredictors + 1):
            return None
        H = self.beta @ matrix.transpose()
        nobs = matrix.shape[0]
        KResults = np.zeros((nobs, self.nactive_set_vectors))
        KTemp = np.zeros((nobs, self.npredictors))
        for idv in range(self.nactive_set_vectors):
            active_vector = self.active_set_vectors[idv]
            KTemp[:, 0] = np.power((active_vector[0] - matrix[:, 1]), 2) / np.power(self.kernel.length_scale[0], 2)
            KTemp[:, 1] = np.power((active_vector[1] - matrix[:, 2]), 2) / np.power(self.kernel.length_scale[1], 2)
            KTemp[:, 2] = np.power((active_vector[2] - matrix[:, 3]), 2) / np.power(self.kernel.length_scale[2], 2)
            KTemp[:, 3] = np.power((active_vector[3] - matrix[:, 4]), 2) / np.power(self.kernel.length_scale[3], 2)
            KTemp[:, 4] = np.power((active_vector[4] - matrix[:, 5]), 2) / np.power(self.kernel.length_scale[4], 2)
            KTemp[:, 5] = np.power((active_vector[5] - matrix[:, 6]), 2) / np.power(self.kernel.length_scale[5], 2)
            KSum = np.sum(KTemp, 1)
            KResults[:, idv] = (self.kernel.sigmaf ** 2) * (
                    (1 + ((1 / (2 * self.kernel.alpharq)) * KSum)) ** -self.kernel.alpharq)
            KResults[:, idv] = KResults[:, idv] * self.alpha[idv]
        K = np.sum(KResults, 1)
        result = H + K
        result = np.power(10, result)
        return result

    ##FAST CHECK OF INPUT VALID
    def check_chla_valid(self, array_443, array_490, array_560, array_665):
        indices = np.where(
            np.logical_and(np.logical_and(array_443 > 0, array_490 > 0), np.logical_and(array_560 > 0, array_665 > 0)))
        nvalid = len(indices[0])
        return nvalid

    # EFFECTIVE METHOD FOR RETREIVING CHLA FROM ARRAYS, CHECK THE VALID INDICES AND USE FAST IMPLEMENTATION -SeaSARC
    def compute_chla_from_2d_arrays(self, array_long, day, array_443, array_490, array_560, array_665):
        chla_array = np.zeros(array_long.shape)
        chla_array[:] = -999
        indices = np.where(
            np.logical_and(np.logical_and(array_443 > 0, array_490 > 0), np.logical_and(array_560 > 0, array_665 > 0)))
        nvalid = len(indices[0])
        if nvalid == 0:
            return chla_array
        input_matrix = np.ones((nvalid, self.npredictors + 1))
        input_matrix[:, 1] = array_long[indices]
        input_matrix[:, 2] = day
        input_matrix[:, 3] = np.log10(array_443[indices])
        input_matrix[:, 4] = np.log10(array_490[indices])
        input_matrix[:, 5] = np.log10(array_560[indices])
        input_matrix[:, 6] = np.log10(array_665[indices])
        chla_1d = self.compute_chla_from_matrix(input_matrix)
        chla_array[indices] = chla_1d
        return chla_array

    #CIAO###############################################################################################################
    ##compute chla from params (no log-transformed rrs) for CIAO algorithm
    def compute_ciao_from_param(self, day, val_443, val_490, val_510, val_560):
        feature_vector = [day, val_443, val_490, val_510, val_560]
        return self.compute_chla(feature_vector, True)