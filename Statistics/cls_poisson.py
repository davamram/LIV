import numpy as np
from scipy.stats import poisson

def calculate_likelihood(data, lambdas):
    likelihood = 1.0
    for i, x in enumerate(data):
        likelihood *= poisson.pmf(round(x), round(lambdas[i]))
        print("x : ", x)
        print("lambda : ", lambdas[i])
        print("pmf : ", poisson.pmf(round(x), round(lambdas[i])))
        print("likelihood : ", likelihood)
        sep='-'*10
        print(sep)
    return likelihood

def calculate_test_statistic(data, lambda_h1, lambda_h0):
    likelihood_h1 = calculate_likelihood(data, lambda_h1)
    likelihood_h0 = calculate_likelihood(data, lambda_h0)
    a = likelihood_h1 / likelihood_h0
    return -2 * np.log(likelihood_h1 / likelihood_h0)

def calculate_pvalue(data, lambda_h1, lambda_h0):
    observed_statistic = calculate_test_statistic(data, lambda_h1, lambda_h0)
    p_value = 1 - poisson.sf(observed_statistic, 1)
    return p_value

def calculate_cls(data, lambda_h1, lambda_h0):
    p_value_h1 = calculate_pvalue(data, lambda_h1, lambda_h0)
    p_value_h0 = calculate_pvalue(data, lambda_h0, lambda_h0)
    cls = p_value_h1 / p_value_h0
    return cls

# Note that those values are not right. They will change after having perform the good simulation
data = [2.284340e+00, 1.067280e+00, 5.424560e-01, 2.317930e-01, 8.713560e-02, 3.699980e-02, 1.789660e-02, 8.427990e-03, 3.592930e-03, 1.498080e-03, 6.316090e-04, 2.499870e-04, 7.636430e-05, 1.399750e-05, 1.098470e-06, 2.185840e-07]
liv_sample = [9.776820e-09, 1.000000e-10, 9.776820e-09, 1.000000e-10, 1.000000e-10, 4.888410e-09, 9.776822e-09, 6.983444e-09, 1.527629e-08, 9.776821e-09, 8.310298e-08, 2.102017e-07, 1.371443e-05, 8.595660e-06, 8.823582e-07, 6.306050e-08]
sm_sample = [1.132219e+00, 6.113484e-01, 3.190417e-01, 1.421000e-01, 5.346908e-02, 2.384420e-02, 1.092171e-02, 5.495240e-03, 2.086177e-03, 8.215948e-04, 3.935730e-04, 1.067489e-04, 4.660717e-05, 8.345068e-06, 2.044124e-09, 1.000000e-10]

# Energy-bin size
Energy_bin = [0.25e+02, 0.25e+02, 0.25e+02, 0.5e+02, 0.5e+02, 0.5e+02, 0.5e+02, 0.7e+02, 0.8e+02, 1e+02, 1e+02, 1.5e+02, 2e+02, 4e+02, 5e+02, 5e+02]
luminosity = 36e03

# This is done to go from d(sigma)/dE to Nevents
data = [x * Energy_bin[i] * luminosity for i, x in enumerate(data)]
liv_sample = [x * Energy_bin[i] * luminosity for i, x in enumerate(liv_sample)]
sm_sample = [x * Energy_bin[i] * luminosity for i, x in enumerate(sm_sample)]

cls = calculate_cls(data, liv_sample, sm_sample)
print("CLs:", cls)