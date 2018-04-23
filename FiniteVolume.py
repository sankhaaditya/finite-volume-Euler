import numpy as np
import time
import multiprocessing as mp
import itertools as it
import VarHolders as vh
import Config as cf
import Riemann as R


def initCellsInterfaces(length, cell_count, ic, bc):

    cells = []

    index = 1
    
    for i in range(0, len(ic)):
        for j in range(int(ic[i][0]), int(ic[i][1])+1):
            
            cell = vh.Cell(ic[i][2], ic[i][3], ic[i][4], None)
            
            cell.setInterfaces(index-1, index)
            
            cells.append(cell)

            index += 1

    if bc[0] == 0:  # Transparent boundary condition

        d, u, p = cells[0].getPrimitiveVars()
        
        cells.insert(0, vh.Cell(d, u, p, 1)) # Shadow cell at beginning

    if bc[1] == 0:

        d, u, p = cells[len(cells)-1].getPrimitiveVars()
        
        cells.insert(len(cells), vh.Cell(d, u, p, len(cells)-1))   # Shadow cell at end

    interfaces = []

    for i in range(0, cell_count+1):

        interface = vh.Interface(i, i+1)

        interfaces.append(interface)

    return cells, interfaces, (0, len(cells)-1)


def solveRiemann(interface, cells):

    index_l, index_r = interface.getCells()

    cell_l = cells[index_l]

    cell_r = cells[index_r]

    d, u, p, S = R.exact(cell_l, cell_r)

    interface.setPrimitiveVars(d, u, p)

    interface.setS(S)

    return interface


def calcTimeStep(dx, S_max):

    return cf.cfl_coefficient*dx/S_max


def update(cell, interfaces, dx, time_step):

    if cell.shadow != None:
            
        return cell

    index_l, index_r = cell.getInterfaces()

    F_1l, F_2l, F_3l = interfaces[index_l].getFluxes()

    F_1r, F_2r, F_3r = interfaces[index_r].getFluxes()

    CV_1, CV_2, CV_3 = cell.getConservedVars()

    CV_1 += time_step/dx*(F_1l-F_1r)

    CV_2 += time_step/dx*(F_2l-F_2r)

    CV_3 += time_step/dx*(F_3l-F_3r)

    cell.setConservedVars(CV_1, CV_2, CV_3)

    return cell


def run(domain_length, cell_count, ic, bc, time_eval):
    
    
    # Initialize list of cells with values

    cells, interfaces, shadowCells = initCellsInterfaces(domain_length, cell_count, ic, bc)

    dx = domain_length/cell_count


    # Starting time marching

    t = 0.0
    
    cpuTimeRiemann = 0.0
    
    cpuTimeUpdate = 0.0

    if cf.parallel == 'On':

        pool = mp.Pool(processes=2)
    
    done = False

    while True:
        
            
        # Solve Riemann problem at each interface
        
        t0 = time.time()

        if cf.parallel == 'On':
        
            interfaces = pool.starmap(solveRiemann, zip(interfaces, it.repeat(cells)))

        else:

            for i in range(0, len(interfaces)):

                interfaces[i] = solveRiemann(interfaces[i], cells)
                
        t1 = time.time()

        cpuTimeRiemann += (t1-t0)

        S_max = 0.0

        for i in range(0, len(interfaces)):
            
            if abs(interfaces[i].getS()) > S_max:

                S_max = interfaces[i].getS()
                

        # Calculating time step

        time_step = calcTimeStep(dx, S_max)

        if t+time_step >= time_eval:

            time_step = time_eval-t

            done = True
            

        # Update cells to next time
        
        t0 = time.time()
        
        if cf.parallel == 'On':
        
            cells = pool.starmap(update, zip(cells, it.repeat(interfaces), it.repeat(dx), it.repeat(time_step)))

        else:

            for i in range(0, len(cells)):

                cells[i] = update(cells[i], interfaces, dx, time_step)
                
        t1 = time.time()

        cpuTimeUpdate += (t1-t0)
        
        for i in range(0, len(shadowCells)):    # CHECK LATER

            cells[shadowCells[i]].d = cells[cells[shadowCells[i]].shadow].d

            cells[shadowCells[i]].u = cells[cells[shadowCells[i]].shadow].u

            cells[shadowCells[i]].p = cells[cells[shadowCells[i]].shadow].p

        t += time_step

        if done:

            break
        
    print(cpuTimeRiemann, cpuTimeUpdate)

    if cf.parallel == 'On':

        pool.close()

        pool.join()

    return cells
