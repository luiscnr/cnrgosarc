import numpy as np


class ARC_GPR_MODEL():

    def __init__(self, fmodel):
        print('STARTED MODEL')
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

    def compute_chla(self, feature_vector, transform_input):
        if len(feature_vector) != self.npredictors:
            return np.nan
        if transform_input:
            for idx in range(2, self.npredictors):
                feature_vector[idx] = np.log10(feature_vector[idx])
        val_fin = self.compute_chla_impl(feature_vector)
        chla = 10 ** val_fin
        return chla

    def compute_chla_from_param(self, long_value, day, val_443, val_490, val_560, val_665):
        feature_vector = [long_value, day, val_443, val_490, val_560, val_665]
        return self.compute_chla(feature_vector, True)

    def compute_chla_from_matrix(self, matrix, transform_input):
        matrix = np.zeros((16, 6))
        if len(matrix.shape) != 2 or matrix.shape[1] != self.npredictors:
            return None
        nobs = matrix.shape[0]
        result = np.zeros(nobs)
        for idx in range(nobs):
            result[idx] = self.compute_chla(matrix[idx], transform_input)
        return result

    def compute_chla_from_2d_arrays(self, array_long, day, array_443, array_490, array_560, array_665):

        chla_array = np.zeros(array_long.shape)
        chla_array[:] = -999

        indices = np.where(
            np.logical_and(np.logical_and(array_443 > 0, array_490 > 0), np.logical_and(array_560 > 0, array_665 > 0)))
        nvalid = len(indices[0])
        if nvalid == 0:
            return chla_array

        input_matrix = np.zeros(nvalid, self.npredictors)
        input_matrix[:, 0] = array_long[indices]
        input_matrix[:, 1] = day
        input_matrix[:, 2] = array_443[indices]
        input_matrix[:, 3] = array_490[indices]
        input_matrix[:, 4] = array_560[indices]
        input_matrix[:, 5] = array_665[indices]

        chla_1d = self.compute_chla_from_matrix(input_matrix, True)

        chla_array[indices] = chla_1d

        return chla_array
