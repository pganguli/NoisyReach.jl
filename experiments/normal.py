#!/bin/env python3

import sys

import numpy as np

#lambdas = np.random.normal(loc=0.0, scale=0.533, size=1_000_000)
#lambdas = np.random.normal(loc=0.0, scale=0.264, size=1_000_000)
lambdas = np.random.normal(loc=0.0, scale=float(sys.argv[1]), size=1_000_000)

print(np.mean(np.abs(lambdas)))
