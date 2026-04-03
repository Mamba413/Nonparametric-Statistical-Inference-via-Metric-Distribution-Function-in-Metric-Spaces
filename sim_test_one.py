import sys
import numpy as np
from geomstats.geometry.pre_shape import PreShapeSpace, KendallShapeMetric
from geomstats.geometry.euclidean import EuclideanMetric
from geomstats.geometry.hypersphere import HypersphereMetric
from geomstats.geometry.spd_matrices import SPDMetricAffine
from geomstats.geometry.positive_lower_triangular_matrices import CholeskyMetric
from frechetKSampleTest import frechetKSampleTest
import pyreadr

R_data = pyreadr.read_r('hypothesis_test_{0}_{1}_{2}.RData'.format(sys.argv[1], sys.argv[2], sys.argv[3]))
# ## For debuging:
# import os
# os.chdir('./code/')
## read data:
# R_data = pyreadr.read_r('hypothesis_test_{0}_{1}_{2}.RData'.format("2-2", "1", "10"))

data_type = R_data["data_type"].to_numpy()
data_type = data_type[0][0]
distance_type = R_data["distance_type"].to_numpy()
distance_type = distance_type[0][0]
testing_data = R_data["testing_data"].to_numpy()

## convert to python layout:
if data_type == "shape" or data_type == "spd":
    k_landmarks, m_ambient, num = testing_data.shape
    testing_data_py = np.zeros([num, k_landmarks, m_ambient])
    for i in range(num):
        testing_data_py[i, :, :] = testing_data[:, :, i]
        pass
elif data_type == "multivariate" or data_type == "sphere":
    num, dim = testing_data.shape
    pass    

if data_type == "shape":
    if distance_type == "Euclidean":
        test_metric = EuclideanMetric(k_landmarks*m_ambient)
        testing_data_py_new = np.zeros([num, k_landmarks*m_ambient])
        for i in range(num):
            testing_data_py_new[i, :] = testing_data_py[i, :, :].flatten()
            pass
        
        testing_data_py = testing_data_py_new
    else:
        preshape = PreShapeSpace(m_ambient=m_ambient, k_landmarks=k_landmarks)
        testing_data_py = preshape.projection(testing_data_py)
        testing_data_py = preshape.align(point=testing_data_py, base_point=testing_data_py[0])
        test_metric = KendallShapeMetric(m_ambient=m_ambient, k_landmarks=k_landmarks)
elif data_type == "spd":
    if distance_type == "Riemanne":
        test_metric = SPDMetricAffine(testing_data_py.shape[1])
    else:
        for i in range(num):
            testing_data_py[i, :, :] = np.linalg.cholesky(testing_data_py[i, :, :])
            pass
        test_metric = CholeskyMetric(testing_data_py.shape[1])
    pass
elif data_type == "multivariate":
    test_metric = EuclideanMetric(dim)
    testing_data_py = np.copy(testing_data)
elif data_type == "sphere":
    test_metric = HypersphereMetric(dim-1)
    testing_data_py = np.copy(testing_data)

group_size = R_data["group_size"]
group_size = group_size.to_numpy().flatten().astype(int)
test_sample = np.split(testing_data_py, np.cumsum(group_size)[:-1], axis=0)

import time
t1 = time.time()
try:
    res = frechetKSampleTest(test_sample, test_metric)
except:
    res = np.array([0.0, 1.0])
# print("stats: {}, pvalue: {}, runtime: {:2f}".format(res[0], res[1], time.time() - t1))

np.savetxt("frechet_test_result_{0}_{1}_{2}.txt".format(sys.argv[1], sys.argv[2], sys.argv[3]), res)
