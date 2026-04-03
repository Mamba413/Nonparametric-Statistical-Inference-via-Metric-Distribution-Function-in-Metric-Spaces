from geomstats.learning.frechet_mean import FrechetMean, variance
from geomstats.geometry.euclidean import EuclideanMetric
from geomstats.geometry.hypersphere import HypersphereMetric
import numpy as np

def frechetKSampleTest(samples, metric, num_rounds=199, seed=0, max_iter=32, epsilon=0.0001):
    pool_frechet_var = computePoolFrechetVariance(samples, metric, max_iter, epsilon)
    statistics = frechetKSampleStatistic(pool_frechet_var, samples, metric, max_iter, epsilon)    
    sample_size = np.array([x.shape[0] for x in samples])
    pool_samples = np.concatenate(samples, axis=0)
    permuted_statistics_vec = []
    for r in range(num_rounds):
        np.random.seed(seed)
        np.random.shuffle(pool_samples)
        new_samples = np.split(pool_samples, np.cumsum(sample_size)[:-1])
        permuted_statistics = frechetKSampleStatistic(pool_frechet_var, new_samples, metric, max_iter, epsilon)
        permuted_statistics_vec.append(permuted_statistics)
        pass
    permuted_statistics_vec = np.array(permuted_statistics_vec)

    p_value = np.append(permuted_statistics_vec >= statistics, [True]).mean()
    return np.array([statistics, p_value])

def computePoolFrechetVariance(samples, metric, max_iter=32, epsilon=0.0001):
    # init_pool_mean = np.concatenate(init_pool_mean, axis=0)
    # mean.fit(init_pool_mean)
    # init_pool_mean = mean.estimate_
    pool_samples = np.concatenate(samples, axis=0)
    if metric.__class__ is EuclideanMetric or metric.__class__ is HypersphereMetric:
        point_type='vector'
    else:
        point_type='matrix'
    mean = FrechetMean(metric=metric, max_iter=max_iter, epsilon=epsilon, point_type=point_type)
    mean.fit(pool_samples)
    pool_f_mean = mean.estimate_
    pool_f_var = variance(pool_samples, base_point=pool_f_mean, metric=metric, point_type=point_type)
    return pool_f_var

def frechetKSampleStatistic(pool_frechet_variance, samples, metric, max_iter=32, epsilon=0.0001):
    '''Frechet test statistic T_n in Reference 
    The k-sample Anderson-Darling test is a modification of the
    one-sample Anderson-Darling test. It tests the null hypothesis
    that k-samples are drawn from the same population without having
    to specify the distribution function of that population. The
    critical values depend on the number of samples.
    Parameters
    ----------
    samples : sequence of array_like
        Array of sample data in arrays.
    metric : function
        metric function
    Returns
    -------
    statistic : float
        Normalized k-sample Anderson-Darling test statistic.
    p-value : float
        p-value of hypothesis test
    References
    ----------
    .. [1] Fréchet analysis of variance for random objects, Biometrika.
    Examples
    --------
    '''
    sample_size = np.array([x.shape[0] for x in samples])
    pool_sample_size = np.sum(sample_size)
    sample_prop = sample_size / pool_sample_size
    f_var_vec = []
    f_var_est_vec = []
    is_euclidean = metric.__class__ is EuclideanMetric or metric.__class__ is HypersphereMetric
    
    if is_euclidean:
        point_type='vector'
    else:
        point_type='matrix'
    mean = FrechetMean(metric=metric, max_iter=max_iter, epsilon=epsilon, point_type=point_type)
    # init_pool_mean = []
    for one_samples in samples:
        ## sample Frechet mean:
        mean.fit(one_samples)
        f_mean = mean.estimate_
        # init_pool_mean.append([f_mean])
        
        ## sample Frechet variance:
        f_var = variance(one_samples, base_point=f_mean, metric=metric, point_type=point_type)
        f_var_vec.append(f_var)
        
        ## variance estimation of Frechet variance:
        fourth_degree_dist = np.square(metric.squared_dist(one_samples, f_mean))
        f_var_est = np.mean(fourth_degree_dist) - np.square(f_var)
        f_var_est_vec.append(f_var_est)
        pass
    f_var_vec = np.array(f_var_vec)
    f_var_est_vec = np.array(f_var_est_vec)

    F_n = pool_frechet_variance - np.sum(sample_prop * f_var_vec)
    # print("F_n (shall be a.s. nonnegative): {}".format(F_n))
    U_n_term1 = f_var_vec.reshape(-1, 1) - f_var_vec.transpose()
    U_n_term2 = sample_prop.reshape(-1, 1) * sample_prop.transpose()
    U_n_term3 = f_var_est_vec.reshape(-1, 1) * f_var_est_vec.transpose()
    U_n_mat = np.square(U_n_term1) * U_n_term2 / U_n_term3
    U_n_mat = np.triu(U_n_mat, k=1).flatten()
    U_n = U_n_mat[np.nonzero(U_n_mat)]
    U_n = np.mean(U_n)

    U_n_weight = pool_sample_size / np.sum(sample_prop / f_var_est_vec)
    F_n_weight = pool_sample_size / np.sum(np.square(sample_prop) * f_var_est_vec)
    T_n = U_n * U_n_weight + np.square(F_n) * F_n_weight

    statistics = T_n
    return statistics
