import numpy as np
import matplotlib.pyplot as plt
import json
import VarHolders as vh
import Riemann as R


def test():

    with open('config.json', 'r') as f:
        config = json.load(f)

    test_case = config['TEST_RIEMANN_INPUT']['TEST_CASE']
    d_l = config['TEST_RIEMANN_INPUT'][test_case]['DENS_L']
    u_l = config['TEST_RIEMANN_INPUT'][test_case]['VELO_L']
    p_l = config['TEST_RIEMANN_INPUT'][test_case]['PRES_L']
    d_r = config['TEST_RIEMANN_INPUT'][test_case]['DENS_R']
    u_r = config['TEST_RIEMANN_INPUT'][test_case]['VELO_R']
    p_r = config['TEST_RIEMANN_INPUT'][test_case]['PRES_R']
    domain_length = config['TEST_RIEMANN_INPUT'][test_case]['DOMAIN_LENGTH']
    divisions = config['TEST_RIEMANN_INPUT'][test_case]['DIVISIONS']
    time_eval = config['TEST_RIEMANN_INPUT'][test_case]['TIME_EVAL']
    
    # Initializing variables

    cell_l = vh.Cell(d_l, u_l, p_l, False)
    cell_r = vh.Cell(d_r, u_r, p_r, False)

    x = np.zeros([divisions+1])

    d = np.zeros([divisions+1])

    u = np.zeros([divisions+1])

    p = np.zeros([divisions+1])


    # Calculating exact Riemann solutions

    d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r = R.extractValues(cell_l, cell_r)

    p_0, u_0 = R.starRegionEstimate(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r)

    S_max = R.calcMaxSlope(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r, p_0, u_0)


    # Evaluating at t = time_eval

    for i in range(0, divisions+1):

        S = (i-int(divisions/2))*domain_length/time_eval/divisions

        x[i] = i*domain_length/divisions

        d[i], u[i], p[i] = R.sample(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r, p_0, u_0, S)

    return x, d, u, p

# Plotting

#plt.plot(x, d, label = 'Density')

#plt.plot(x, u, label = 'Velocity')

#plt.plot(x, p, label = 'Pressure')

#plt.xlabel('Position')

#plt.ylabel('Quantity')

#plt.legend()

#plt.show()
