import json
import Config as cf


with open('config.json', 'r') as f:
    config = json.load(f)

pressureApprox = config['CONSTANTS']['RIEMANN']['PRESSURE_APPROX'] # Pressure Approximation Method
TOL = config['CONSTANTS']['RIEMANN']['TOL'] # Tolerence for error in pressure guessing



def exact(cell_l, cell_r):
    

    #Extracting values from cells

    d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r = extractValues(cell_l, cell_r)


    # Escaping with error if vacuum state CHECK LATER

    if (u_r-u_l) >= 2*(a_l+a_r)/(cf.gama-1):

        print('VACUUM STATE! ABORT!')

        return None


    # Calculating pressure and velocity estimate in Star Region
    
    p_0, u_0 = starRegionEstimate(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r)


    # Sampling values at S = x/t = 0

    d, u, p = sample(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r, p_0, u_0, 0)


    # Calculating the maximum slope

    S_max = calcMaxSlope(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r, p_0, u_0)


    return d, u, p, S_max
    
    

def extractValues(cell_l, cell_r):
    

    d_l = cell_l.d
    d_r = cell_r.d

    u_l = cell_l.u
    u_r = cell_r.u

    p_l = cell_l.p
    p_r = cell_r.p

    a_l = cell_l.getSoundSpeed()
    a_r = cell_r.getSoundSpeed()
    

    return d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r



def starRegionEstimate(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r):
    

    # Initial pressure guessing functions

    def TRApprox():

        pTR = ((a_l+a_r-0.5*(cf.gama-1)*(u_r-u_l))/(a_l/p_l**((cf.gama-1)/2/cf.gama)+a_r/p_r**((cf.gama-1)/2/cf.gama)))**(2*cf.gama/(cf.gama-1))

        return pTR

    def PVApprox():

        pPV = max(TOL, 0.5*(p_l+p_r)-0.125*(u_r-u_l)*(d_l+d_r)*(a_l+a_r))

        return pPV
    

    # Guessing pressure in Star Region

    if pressureApprox == 'TR':
        
        p_0 = TRApprox()
        
    elif pressureApprox == 'PV':
        
        p_0 = PVApprox()


    # Pressure function

    def f_k(p, k, deriv):

        if k == 'l':

            d_k = d_l
            u_k = u_l
            p_k = p_l
            a_k = a_l

        if k == 'r':

            d_k = d_r
            u_k = u_r
            p_k = p_r
            a_k = a_r

        A_k = 2/(cf.gama+1)/d_k
        B_k = (cf.gama-1)/(cf.gama+1)*p_k

        if deriv == 0:  # No derivative

            if p > p_k:

                return (p-p_k)*(A_k/(p+B_k))**0.5

            else:

                return 2*a_k/(cf.gama-1)*((p/p_k)**((cf.gama-1)/2/cf.gama)-1)

        if deriv == 1:  # First derivative

            if p > p_k:

                return (A_k/(B_k+p))**0.5*(1-(p-p_k)/2/(B_k+p))

            else:

                return 1/d_k/a_k*(p/p_k)**(-(cf.gama+1)/2/cf.gama)

    def f(p):

        return f_k(p,'l',0)+f_k(p,'r',0)+(u_r-u_l)

    def fDeriv(p):

        return f_k(p,'l',1)+f_k(p,'r',1)
    

    # Newton-Raphson iteration

    while True:

        p_next = p_0-f(p_0)/fDeriv(p_0)

        CHA = abs(p_next-p_0)/0.5/(p_next+p_0)

        p_0 = p_next

        if CHA < TOL:
            
            break


    # Velocity in Star Region

    u_0 = 0.5*(u_l+u_r)+0.5*(f_k(p_0,'r',0)-f_k(p_0,'l',0))

    return p_0, u_0



def sample(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r, p_0, u_0, S):

    if S < u_0: # To left of contact discontinuity

        if p_l < p_0:   # Left shock

            S_l = u_l-a_l*((cf.gama+1)*p_0/2/cf.gama/p_l+(cf.gama-1)/2/cf.gama)**0.5

            if S < S_l:

                return d_l, u_l, p_l

            else:

                d_0l = d_l*((p_0/p_l+(cf.gama-1)/(cf.gama+1))/((cf.gama-1)/(cf.gama+1)*p_0/p_l+1))

                return d_0l, u_0, p_0

        else:   # Left fan

            a_0l = a_l*(p_0/p_l)**((cf.gama-1)/2/cf.gama)

            S_hl = u_l-a_l

            S_tl = u_0-a_0l

            if S < S_hl:

                return d_l, u_l, p_l

            else:

                if S > S_tl:

                    d_0l = d_l*(p_0/p_l)**(1/cf.gama)

                    return d_0l, u_0, p_0

                else:

                    d_lfan = d_l*(2/(cf.gama+1)+(cf.gama-1)/(cf.gama+1)/a_l*(u_l-S))**(2/(cf.gama-1))

                    u_lfan = 2/(cf.gama+1)*(a_l+(cf.gama-1)/2*u_l+S)

                    p_lfan = p_l*(2/(cf.gama+1)+(cf.gama-1)/(cf.gama+1)/a_l*(u_l-S))**(2*cf.gama/(cf.gama-1))

                    return d_lfan, u_lfan, p_lfan

    else:   # To right of contact discontinuity

        if p_r < p_0:   # Right shock

            S_r = u_r+a_r*((cf.gama+1)*p_0/2/cf.gama/p_r+(cf.gama-1)/2/cf.gama)**0.5

            if S > S_r:

                return d_r, u_r, p_r

            else:

                d_0r = d_r*((p_0/p_r+(cf.gama-1)/(cf.gama+1))/((cf.gama-1)/(cf.gama+1)*p_0/p_r+1))

                return d_0r, u_0, p_0

        else:   # Right fan

            a_0r = a_r*(p_0/p_r)**((cf.gama-1)/2/cf.gama)

            S_hr = u_r+a_r

            S_tr = u_0+a_0r

            if S > S_hr:

                return d_r, u_r, p_r

            else:

                if S < S_tr:

                    d_0r = d_r*(p_0/p_r)**(1/cf.gama)

                    return d_0r, u_0, p_0

                else:

                    d_rfan = d_r*(2/(cf.gama+1)-(cf.gama-1)/(cf.gama+1)/a_r*(u_r-S))**(2/(cf.gama-1))

                    u_rfan = 2/(cf.gama+1)*(-a_r+(cf.gama-1)/2*u_r+S)

                    p_rfan = p_r*(2/(cf.gama+1)-(cf.gama-1)/(cf.gama+1)/a_r*(u_r-S))**(2*cf.gama/(cf.gama-1))

                    return d_rfan, u_rfan, p_rfan

                    

def calcMaxSlope(d_l, d_r, u_l, u_r, p_l, p_r, a_l, a_r, p_0, u_0):
    

    if p_l < p_0:   # Left shock

        S_l = u_l-a_l*((cf.gama+1)*p_0/2/cf.gama/p_l+(cf.gama-1)/2/cf.gama)**0.5

        S_lmax = S_l

    else:       # Left fan

        a_0l = a_l*(p_0/p_l)**((cf.gama-1)/2/cf.gama)

        S_hl = u_l-a_l

        S_lmax = S_hl

    if p_r < p_0:   # Right shock

        S_r = u_r+a_r*((cf.gama+1)*p_0/2/cf.gama/p_r+(cf.gama-1)/2/cf.gama)**0.5

        S_rmax = S_r

    else:       # Right fan

        a_0r = a_r*(p_0/p_r)**((cf.gama-1)/2/cf.gama)

        S_hr = u_r+a_r

        S_rmax = S_hr


    return max(abs(S_lmax), abs(S_rmax))
