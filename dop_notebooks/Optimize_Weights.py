#   This is a code to optimize the weights for the radii of gyration
#   It calls the method weight_optimizer.py

import sys
import os
import numpy as np
from array import array

import BSTCalcMethods
import BSTFileMethods
import BSTPlotMethods
import BSTUtilities

#   -----------------------------------------------------------------

initial_weights, adjusted_weights = BSTCalcMethods.weight_optimizer()

# print('')
# print('Initial Weights: ',initial_weights)
# print('')



