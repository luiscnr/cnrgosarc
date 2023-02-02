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

        self.sigma = model_dict['Sigma']
        print('Model sigma es: ', self.sigma)
        self.beta = np.array(model_dict['Beta'])
        print('Beta data es: ', self.beta)

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

        kresult = np.zeros(self.nactive_set_vectors)
        for idx in range(self.nactive_set_vectors):
            active_vector = self.active_set_vectors[idx]
            kresult[idx] = self.kernel.compute_kernel(active_vector, feature_vector)
        K = kresult @ self.alpha

        print('FINAL: ', val)

        X = np.concatenate(([1],feature_vector))
        H = X @ self.beta
        # valnew = self.beta[0]
        # for idx in range(1, 7):
        #     valnew = valnew + (self.beta[idx] * feature_vector[idx - 1])

        print('CHECK: ', H)

        val_fin = H + val

        print('VALFIN: ', val_fin)

        return val
