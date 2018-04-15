import numpy as np
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


def solveRiemann(cells, interfaces):

    S = np.zeros([len(interfaces)])

    for i in range(0, len(interfaces)):

        index_l, index_r = interfaces[i].getCells()

        cell_l = cells[index_l]

        cell_r = cells[index_r]

        d, u, p, S[i] = R.exact(cell_l, cell_r)

        interfaces[i].setPrimitiveVars(d, u, p)

    S_max = max(S)

    return interfaces, S_max


def calcTimeStep(dx, S_max):

    return cf.cfl_coefficient*dx/S_max


def update(cells, interfaces, shadowCells, dx, time_step):
    
    for i in range(0, len(cells)):

        if cells[i].shadow != None:
            
            continue

        index_l, index_r = cells[i].getInterfaces()

        F_1l, F_2l, F_3l = interfaces[index_l].getFluxes()

        F_1r, F_2r, F_3r = interfaces[index_r].getFluxes()

        CV_1, CV_2, CV_3 = cells[i].getConservedVars()

        CV_1 += time_step/dx*(F_1l-F_1r)

        CV_2 += time_step/dx*(F_2l-F_2r)

        CV_3 += time_step/dx*(F_3l-F_3r)

        cells[i].setConservedVars(CV_1, CV_2, CV_3)

    for i in range(0, len(shadowCells)):    # CHECK LATER

        cells[shadowCells[i]].d = cells[cells[shadowCells[i]].shadow].d

        cells[shadowCells[i]].u = cells[cells[shadowCells[i]].shadow].u

        cells[shadowCells[i]].p = cells[cells[shadowCells[i]].shadow].p

    return cells
