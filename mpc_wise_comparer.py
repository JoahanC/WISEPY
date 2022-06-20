from visualizer_functions import visualizer
from mpc_wise_functions import *
import sys

new_epochs = data_comparer(sys.argv[1], sys.argv[2], True)
visualizer(sys.argv[1], sys.argv[2], sys.argv[3])