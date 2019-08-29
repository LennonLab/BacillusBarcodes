from __future__ import division
import math
import numpy as np
#from statsmodels.base.model import GenericLikelihoodModel
import scipy.optimize as op
import scipy.io as sio

'''
This code modifies the approach introduced in Li et al. to estimate fitness
of genotypes using temporal barcoded sequencing. We examine only pair-wise
competition between a mutant and ancestor, we're only interested in the
fitness of the mutant. We assume two timepoints.

Li et al., 2018, Cell Systems 7, 521–525
https://doi.org/10.1016/j.cels.2018.09.004
'''

R = 2

num_observed_genotypes = [2, 2]
read_counts_mut = [100000, 120000]
read_counts_anc = [100300, 90000]
times = [0, 6]

k = 2.5


# multiply get_initial_fitness by frequency to get mean fitness of population





def estimate_fitness():
    num_observed_genotypes = [2, 2]
    read_counts_mut = [100000, 120000]
    read_counts_anc = [100300, 90000]
    times = [0, 6]


    def get_fitness(t0, t1, f_mut_t0, f_mut_t1):
        return np.exp((np.log(f_mut_t1) - np.log(f_mut_t0)) / (t1-t0) ) - 1


    def lik_fitness(paramters):
        x_mut = paramters[0]
        x_mean = paramters[1]
        L = 0
        for i in np.arange(0, len(times)):
            r_mut_i = read_counts_mut[i]
            if i == 0:
                r_hat = r_mut_i
            else:

                r_hat1 = ( ((1+x_mut) ** (times[1]-times[0])) * r_mut_i ) / (1 + x_mean)

                r_hat2 = read_counts_mut[i-1] / (2**(times[i]-times[i-1]))
                r_hat = max(r_hat1, r_hat2)

            L -= np.log(  np.sqrt(np.sqrt(r_hat) / ( 4*math.pi*k*(r_mut_i**(3/2)))) * np.exp(-1* ((np.sqrt(r_mut_i)-np.sqrt(r_hat))**2)/k) )
        return L

    f_mut_t0_init = read_counts_mut[0] / (read_counts_anc[0] + read_counts_mut[0])
    f_mut_t1_init = read_counts_mut[1] / (read_counts_anc[1] + read_counts_mut[1])
    x_mut_init = get_fitness(times[0], times[1], f_mut_t0_init, f_mut_t1_init) + 0.01

    f_anc_t0 = read_counts_anc[0] / (read_counts_anc[0] + read_counts_mut[0])
    f_anc_t1 = read_counts_anc[1] / (read_counts_anc[1] + read_counts_mut[1])

    f_mut_t0 = read_counts_mut[0] / (read_counts_anc[0] + read_counts_mut[0])
    f_mut_t1 = read_counts_mut[1] / (read_counts_anc[1] + read_counts_mut[1])

    x_anc_init = get_fitness(times[0], times[1], f_anc_t0, f_anc_t1)

    x_mean_init = (x_mut_init*f_mut_t1) + (x_anc_init*f_anc_t1)

    lik_model = op.minimize(lik_fitness, np.array([x_mut_init, x_mean_init ]), method='L-BFGS-B')

    print(lik_model)

    print(x_mut_init, lik_model['x'][0])



#estimate_fitness()


mat_contents = sio.loadmat('/Users/WRShoemaker/GitHub/BacillusBarcodes/Python/FitSeq.m')




#def lik(parameters):
#    m = parameters[0]
#    b = parameters[1]
#    sigma = parameters[2]
#    for i in np.arange(0, len(x)):
#        y_exp = m * x + b
#    L = (len(x)/2 * np.log(2 * np.pi) + len(x)/2 * np.log(sigma ** 2) + 1 /
#         (2 * sigma ** 2) * sum((y - y_exp) ** 2))

#    print(L)
#    return L



#x = np.array([1,2,3,4,5])
#y = np.array([2,5,8,11,14])
#lik_model = op.minimize(lik, np.array([1]), method='L-BFGS-B')
#plt.scatter(x,y)
#plt.plot(x, lik_model['x'][0] * x + lik_model['x'][1])
#plt.show()
