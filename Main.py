import numpy as np
import json
import FiniteVolume as fv
import PlotUtil as pu


# Extracting configuration properties

with open('config.json', 'r') as f:
    config = json.load(f)

test_inputs = config['MAIN']

domain_length = test_inputs['DOMAIN_LENGTH']
cell_count = test_inputs['CELL_COUNT']

bc = (config['MAIN']['BC']['LEFT'], config['MAIN']['BC']['RIGHT'])

ic = np.zeros([np.size(config['MAIN']['IC']), 5])
for i in range(0, np.size(config['MAIN']['IC'])):
    ic[i][0] = config['MAIN']['IC'][i]['START']
    ic[i][1] = config['MAIN']['IC'][i]['END']
    ic[i][2] = config['MAIN']['IC'][i]['DENS']
    ic[i][3] = config['MAIN']['IC'][i]['VELO']
    ic[i][4] = config['MAIN']['IC'][i]['PRES']

time_eval = config['MAIN']['TIME_EVAL']

plot_factor = config['CONSTANTS']['PLOT_FACTOR']

cells = fv.run(domain_length, cell_count, ic, bc, time_eval)

pu.plot(cells, domain_length/cell_count, plot_factor)



