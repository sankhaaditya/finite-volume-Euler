import numpy as np
import json
import FiniteVolume as fv
import PlotUtil as pu


# Extracting configuration properties

with open('config.json', 'r') as f:
    config = json.load(f)

length = config['DOMAIN']['LENGTH'] # Length of Domain

cell_count = config['DOMAIN']['CELL_COUNT'] # Cell count

bc = (config['BC']['LEFT'], config['BC']['RIGHT'])  # Boundary Conditions

ic = np.zeros([np.size(config['IC']), 5])   # Initial Conditions
for i in range(0, np.size(config['IC'])):
    ic[i][0] = config['IC'][i]['START']
    ic[i][1] = config['IC'][i]['END']
    ic[i][2] = config['IC'][i]['DENS']
    ic[i][3] = config['IC'][i]['VELO']
    ic[i][4] = config['IC'][i]['PRES']

time_eval = config['TIME_EVAL']


# Initialize list of cells with values

cells, interfaces, shadowCells = fv.initCellsInterfaces(length, cell_count, ic, bc)

dx = length/cell_count


# Starting time marching

t = 0.0

done = False

while True:

    # Solve Riemann problem at each interface

    interfaces, S_max = fv.solveRiemann(cells, interfaces)

    # Calculating time step

    time_step = fv.calcTimeStep(dx, S_max)

    if t+time_step >= time_eval:

        time_step = time_eval-t

        done = True

    # Update cells to next time

    cells = fv.update(cells, interfaces, shadowCells, dx, time_step)

    t += time_step

    if done:

        break

pu.plot(cells, dx)
