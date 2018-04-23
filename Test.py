import numpy as np
import json
import sys
import FiniteVolume as fv
import PlotUtil as pu
import TestRiemann as tr

with open('config.json', 'r') as f:
    config = json.load(f)

test_inputs = config['TEST']

domain_length = test_inputs['DOMAIN_LENGTH']
cell_count = test_inputs['CELL_COUNT']

test_case = test_inputs['TEST_CASE']
IC_l = test_inputs[test_case]['IC']['LEFT']
IC_r = test_inputs[test_case]['IC']['RIGHT']
time_eval = test_inputs[test_case]['TIME_EVAL']

plot_factor = config['CONSTANTS']['PLOT_FACTOR']

bc = (0, 0)

if cell_count%2 != 0:
    print('Please use even number of cells')
    sys.exit()


ic = np.zeros([2, 5])

ic[0][0] = 1
ic[0][1] = cell_count/2
ic[0][2] = IC_l['DENS']
ic[0][3] = IC_l['VELO']
ic[0][4] = IC_l['PRES']

ic[1][0] = cell_count/2+1
ic[1][1] = cell_count
ic[1][2] = IC_r['DENS']
ic[1][3] = IC_r['VELO']
ic[1][4] = IC_r['PRES']

if __name__ == '__main__':

    cells = fv.run(domain_length, cell_count, ic, bc, time_eval)

    x, d, u, p = tr.test(domain_length, cell_count*plot_factor, IC_l, IC_r, time_eval)

    pu.plot(cells, domain_length/cell_count, plot_factor, x, d, u, p)
