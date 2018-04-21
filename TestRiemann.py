import numpy as np
import VarHolders as vh
import Riemann as R


def test(domain_length, points, IC_l, IC_r, time_eval):
    
    d_l = IC_l['DENS']
    u_l = IC_l['VELO']
    p_l = IC_l['PRES']
    d_r = IC_r['DENS']
    u_r = IC_r['VELO']
    p_r = IC_r['PRES']

    
    # Initializing variables

    cell_l = vh.Cell(d_l, u_l, p_l, False)
    cell_r = vh.Cell(d_r, u_r, p_r, False)

    x = np.zeros([points])

    d = np.zeros([points])

    u = np.zeros([points])

    p = np.zeros([points])


    # Calculating exact Riemann solutions

    d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r = R.extractValues(cell_l, cell_r)

    p_0, u_0 = R.starRegionEstimate(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r)

    S_max = R.calcMaxSlope(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r, p_0, u_0)


    # Evaluating at t = time_eval

    for i in range(0, points):

        S = (i-int(points/2))*domain_length/time_eval/points

        x[i] = i*domain_length/points

        d[i], u[i], p[i] = R.sample(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r, p_0, u_0, S)

    return x, d, u, p
