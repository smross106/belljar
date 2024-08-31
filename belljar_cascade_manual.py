import belljar_lockstep as bjls
import numpy as np
from CoolProp.CoolProp import PropsSI

fluid="Nitrogen"

class Volume(object):
    def __init__(self):
        self.ID = ""

        self.rho  = 1
        self.rhoV = 1
        self.rhoE = 1

        self.volume = 1e-1
        self.areaIn = 1e-1
        self.areaOut = 1e-1

        self.fluidProp = "Simple"
        self.fluid = "Nitrogen"

        self.rhoFlux  = 0
        self.rhoVFlux = 0
        self.rhoEFlux = 0

        self.drhodt  = 0
        self.drhoVdt = 0
        self.drhoEdt = 0

        self.rho_nm1  = 0
        self.rhoV_nm1 = 0
        self.rhoE_nm1 = 0

        self.rho_np1  = 0
        self.rhoV_np1 = 0
        self.rhoE_np1 = 0

        self.rhoHistory  = []
        self.rhoVHistory = []
        self.rhoEHistory = []

        self.drhodtHistory  = []
        self.drhoVdtHistory = []
        self.drhoEdtHistory = []

        self.intermediaterhoHistory  = []
        self.intermediaterhoVHistory = []
        self.intermediaterhoEHistory = []
        self.intermediatedrhodtHistory  = []
        self.intermediatedrhoEdtHistory = []
        self.intermediatedrhoVdtHistory = []

    def computedVardt(self, time, dt):
        self.drhodt  = self.rhoFlux   / self.volume
        self.drhoVdt = self.rhoVFlux  / self.volume
        self.drhoEdt = self.rhoEFlux  / self.volume

        self.rhoFlux  = 0
        self.rhoVFlux = 0
        self.rhoEFlux = 0
    
    def setNextVar(self, dt):
        self.rhoHistory.append(self.rho_nm1)
        self.drhodtHistory.append((self.rho_np1 - self.rho_nm1)/dt)
        self.rho  = self.rho_np1  * 1.0
        self.rho_np1  = 0
        self.intermediaterhoHistory  = []
        self.intermediatedrhodtHistory  = []

        self.rhoVHistory.append(self.rhoV_nm1)
        self.drhoVdtHistory.append((self.rhoV_np1 - self.rhoV_nm1)/dt)
        self.rhoV = self.rhoV_np1 * 1.0
        self.rhoV_np1 = 0
        self.intermediaterhoVHistory = []
        self.intermediatedrhoVdtHistory = []

        self.rhoEHistory.append(self.rhoE_nm1)
        self.drhoEdtHistory.append((self.rhoE_np1 - self.rhoE_nm1)/dt)
        self.rhoE = self.rhoE_np1 * 1.0
        self.rhoE_np1 = 0
        self.intermediaterhoEHistory = []
        self.intermediatedrhoEdtHistory = []

        self.setOtherProperties()
    
    def resetSolver(self):
        self.rho  = self.rho_nm1
        self.rho_np1  = 0
        self.intermediaterhoHistory  = []
        self.intermediatedrhodtHistory  = []

        self.rhoV = self.rho_nm1
        self.rhoV_np1 = 0
        self.intermediaterhoVHistory = []
        self.intermediatedrhoVdtHistory = []

        self.rhoE = self.rhoE_nm1
        self.rhoE_np1 = 0
        self.intermediaterhoEHistory = []
        self.intermediatedrhoEdtHistory = []
    
    def setOtherProperties(self):
        self.V = self.rhoV / self.rho
        E = self.rhoE / self.rho
        self.mDotIn = self.rhoV * self.areaIn
        self.mDotOut = self.rhoV * self.areaOut

        U = E - (0.5 * self.rho * np.power(self.V, 2))

        if self.fluidProp == "CoolProp":
            self.T = PropsSI("T", "Umass", U, "Dmass", self.rho, self.fluid)
            self.P = PropsSI("P", "Umass", U, "Dmass", self.rho, self.fluid)
            self.s = PropsSI("Smass", "Umass", U, "Dmass", self.rho, self.fluid)
            A = PropsSI("A", "Umass", U, "Dmass", self.rho, self.fluid)
            self.M = self.V / A
            self.h0 = E + (self.P / self.rho)
        else:
            self.T = 298 - (0.5*np.power(self.V, 2)/1005)
            self.P = self.rho * 287 * self.T
            self.h0 = (self.rhoE + self.P) / self.rho
            A = np.sqrt(1.4 * 287 * self.T)
            self.M = self.V / A

        
    def setState(self, M=None, T=None, mDot=None, P=None):
        if M != None and T != None and P != None:
            self.M = M
            self.T = T
            self.P = P
            A = PropsSI("A", "T", self.T, "P", self.P, self.fluid)
            self.V = M * A
            U = PropsSI("Umass", "T", self.T, "P", self.P, self.fluid)
            self.rho = PropsSI("Dmass", "T", self.T, "P", self.P, self.fluid)
            E = U + (0.5 * self.rho * self.V**2)
            self.rhoE = E * self.rho
            self.rhoV = self.V * self.rho
            self.h0 = E + (self.P / self.rho)

volumes = [Volume(), Volume(), Volume(), Volume()]
for v in volumes:
    v.setState(M=0.1, T=298, P=1e5)


pUp = 1.01e5
pDown = 1e5

volumes[0].setState(M=0.1, T=298, P=pUp)

t = 0
dt = 0.0001#0.0001

for v in volumes:
    print(v.rho, v.rhoV, v.rhoE)

while t < dt*100:
    for v in range(1, len(volumes)-1):
        volumes[v].rhoFlux += volumes[v-1].rho * volumes[v-1].V * volumes[v-1].areaOut
        volumes[v].rhoFlux -= volumes[v].rho * volumes[v].V * volumes[v].areaOut

        volumes[v].rhoVFlux += volumes[v-1].rho * volumes[v-1].V * volumes[v-1].V * volumes[v-1].areaOut
        volumes[v].rhoVFlux -= volumes[v].rho * volumes[v].V * volumes[v].V * volumes[v].areaOut
        volumes[v].rhoVFlux += volumes[v-1].areaOut * (volumes[v-1].P)
        volumes[v].rhoVFlux -= volumes[v].areaOut * (volumes[v].P)

        volumes[v].rhoEFlux += volumes[v-1].rho * volumes[v-1].V * volumes[v-1].h0 * volumes[v-1].areaOut
        volumes[v].rhoEFlux -= volumes[v].rho * volumes[v].V * volumes[v].h0 * volumes[v].areaOut


    for v in range(1, len(volumes)-1):
        #print(volumes[v].rho, volumes[v].rhoV, volumes[v].rhoE)
        #print(volumes[v].rhoFlux, volumes[v].rhoVFlux, volumes[v].rhoEFlux)
        volumes[v].computedVardt(t, dt)

        volumes[v].rho += volumes[v].drhodt * dt
        volumes[v].rhoV += volumes[v].drhoVdt * dt
        volumes[v].rhoE += volumes[v].drhoEdt * dt
        volumes[v].setOtherProperties()
        #print(volumes[v].rho, volumes[v].rhoV, volumes[v].rhoE)
        #print("")
    
    t += dt
    if t % 0.1 < dt:
        print(t)
print("")

for v in volumes:
    print(v.rho, v.rhoV, v.rhoE)

for v in volumes:
    print(v.P, v.T, v.M)