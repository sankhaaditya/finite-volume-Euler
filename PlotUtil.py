import numpy as np
import matplotlib.pyplot as plt
import TestRiemann as tr


factor = 10

def plot(cells, dx):    # Assuming 1st and last cells are shadow, CHECK LATER
    
    x = np.zeros([(len(cells)-2)*factor])

    d = np.zeros([(len(cells)-2)*factor])

    u = np.zeros([(len(cells)-2)*factor])

    p = np.zeros([(len(cells)-2)*factor])

    pointer = 0.0

    for i in range(0, len(cells)-2):

        for j in range(0, factor):

            x[i*factor+j] = pointer

            d[i*factor+j] = cells[i+1].d

            u[i*factor+j] = cells[i+1].u

            p[i*factor+j] = cells[i+1].p

            pointer += dx/factor

    x_exact, d_exact, u_exact, p_exact = tr.test()
      
    plt.subplot(311)

    plt.plot(x_exact, d_exact, label='Analytical')

    plt.plot(x, d, label='FV: '+str(len(cells)-2)+' cells')

    plt.ylabel('Density')

    plt.subplot(312)

    plt.plot(x_exact, u_exact, label='Analytical')

    plt.plot(x, u, label='FV: '+str(len(cells)-2)+' cells')

    plt.ylabel('Velocity')

    plt.subplot(313)

    plt.plot(x_exact, p_exact, label='Analytical')

    plt.plot(x, p, label='FV: '+str(len(cells)-2)+' cells')

    plt.ylabel('Pressure')

    plt.xlabel('Position')

    plt.legend()

    plt.show()
