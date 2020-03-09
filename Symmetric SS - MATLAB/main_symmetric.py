import math
import numpy as np
import control
import scipy.io
from ss_symmetric import ss_sym

mat = scipy.io.loadmat('flight1_clean.mat')
array1 = mat['clean_data']

sys = ss_sym()

print(sys)