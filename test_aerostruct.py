# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 15:23:40 2016

@author: samfriedman
"""
from __future__ import print_function
# import math
import numpy as np
import time
# from scipy.spatial.distance import squareform, pdist
# from OpenAeroStruct import coupled
import gibbs_sampler
import cProfile


if __name__ == "__main__":
    doProfile = True
    nsamp = 100
    maxiter = 1
    burn = 10  # burn-in iterations
    tol = 1e-6  # convergence tolerance
    npts = np.array([3,5])
    discL = np.array([100,20])
    discM = np.array([1e-8,20])
    print('Run gibbs sampler with {0} samples, {1} iterations, and npts={2}'.format(nsamp,maxiter,npts))

    if doProfile:
        cProfile.run('gibbs_sampler.gibbs_sampler(discL[0], discL[1], discM[0], discM[1], nsamp, burn, tol, maxiter, npts[0], npts[1])','gibbs_sampler.prof')
        # cProfile.run('gibbs_sampler_python(discL, discM, nsamp, npts, maxiter)')
    else:
        # print('Time Python version...')
        # tic = time.clock()
        # gibbs_sampler_python(discL, discM, nsamp, npts, maxiter)
        # toc_python = time.clock() - tic
        # print('Python time = {0} sec ({1} min)'.format(toc_python,toc_python/60))

        # print('Time Cython version...')
        tic = time.clock()
        gibbs_sampler.gibbs_sampler(discL[0], discL[1], discM[0], discM[1], nsamp, burn, tol, maxiter, npts[0], npts[1])
        toc_cython = time.clock() - tic
        print('Cython time = {0} sec ({1} min)'.format(toc_cython,toc_cython/60))
