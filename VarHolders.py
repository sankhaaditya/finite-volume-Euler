import math as ma
import Config as cf



class VarHolder:

    
    def __init__(self, d, u, p):
        
        self.setPrimitiveVars(d, u, p)


    def setPrimitiveVars(self, d, u, p):

        self.d = d
        self.u = u
        self.p = p


    def getPrimitiveVars(self):

        return self.d, self.u, self.p


    def getInternalEnergy(self):    # Using ideal gas EOS

        return self.p/(cf.gama-1)/self.d


    def getSoundSpeed(self):

        return ma.sqrt(cf.gama*self.p/self.d)


    def getConservedVars(self):

        return self.d, self.d*self.u, self.d*(0.5*self.u**2+self.getInternalEnergy())


    def setConservedVars(self, CV_1, CV_2, CV_3):

        self.d = CV_1
        self.u = CV_2/self.d
        self.p = (CV_3/self.d-0.5*self.u**2)*(cf.gama-1)*self.d

    

class Cell(VarHolder):
    

    def __init__(self, d, u, p, shadow):
        
        self.setPrimitiveVars(d, u, p)
        self.shadow = shadow
        

    def setInterfaces(self, index_l, index_r):
        
        self.index_l = index_l
        self.index_r = index_r


    def getInterfaces(self):

        return self.index_l, self.index_r 
    


class Interface(VarHolder):
    

    def __init__(self, index_l, index_r):

        self.index_l = index_l
        self.index_r = index_r


    def getCells(self):

        return self.index_l, self.index_r


    def getFluxes(self):

        F_1 = self.d*self.u

        F_2 = self.d*self.u**2+self.p

        F_3 = self.u*(self.d*(0.5*self.u**2+self.getInternalEnergy())+self.p)

        return F_1, F_2, F_3
