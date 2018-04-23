import numpy as np
import matplotlib.pyplot as plt
import TestRiemann as tr


def extract(cells, dx, plot_factor):

    x = np.zeros([(len(cells)-2)*plot_factor])

    d = np.zeros([(len(cells)-2)*plot_factor])

    u = np.zeros([(len(cells)-2)*plot_factor])

    p = np.zeros([(len(cells)-2)*plot_factor])

    pointer = 0.0

    for i in range(0, len(cells)-2):

        for j in range(0, plot_factor):

            x[i*plot_factor+j] = pointer

            d[i*plot_factor+j] = cells[i+1].d

            u[i*plot_factor+j] = cells[i+1].u

            p[i*plot_factor+j] = cells[i+1].p

            pointer += dx/plot_factor

    return x, d, u, p


def plot(cells, dx, plot_factor, x_exact=None, d_exact=None, u_exact=None, p_exact=None):    # Assuming 1st and last cells are shadow, CHECK LATER
    
    x, d, u, p = extract(cells, dx, plot_factor)
      
    plt.subplot(311)

    if d_exact != None:

        plt.plot(x_exact, d_exact, label='Analytical')

    plt.plot(x, d, label='FV: '+str(len(cells)-2)+' cells')

    plt.ylabel('Density')

    plt.subplot(312)

    if u_exact != None:

        plt.plot(x_exact, u_exact, label='Analytical')

    plt.plot(x, u, label='FV: '+str(len(cells)-2)+' cells')

    plt.ylabel('Velocity')

    plt.subplot(313)

    if p_exact != None:

        plt.plot(x_exact, p_exact, label='Analytical')

    plt.plot(x, p, label='FV: '+str(len(cells)-2)+' cells')

    plt.ylabel('Pressure')

    plt.xlabel('Position')

    plt.legend()

    plt.show()
