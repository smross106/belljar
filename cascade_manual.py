import numpy as np

T0_in = 298
P_in = 1e5
P_out = 0.85e5
cp = 1005
R = 287
gamma = 1.4
cv = cp/gamma

rho_in = (P_in) / (T0_in * R)
T_out = T0_in * np.power(P_out/P_in, (gamma-1)/gamma)

class Volume(object):
    def __init__(self):
        self.rho = 1
        self.rhoV = 1
        self.rhoE = 1

        self.rhoFlux = 0
        self.rhoVFlux = 0
        self.rhoEFlux = 0

        self.areaIn = 1e-1
        self.areaOut = 1e-1
        self.volume = 1e-2

        self.calculateProperties()

        self.BC = 0
        self.time = []
        self.pHistory = []
        self.VHistory = []
        self.THistory = []
        self.rhoHistory = []
        self.rhoVHistory = []
        self.rhoEHistory = []

    def calculateProperties(self):
        self.V = self.rhoV / self.rho
        v2 = self.V**2
        Tstat = T0_in - 0.5*v2/cp
        self.P = Tstat * self.rho * R
        self.h0 = (self.rhoE + self.P) / self.rho
        self.T = Tstat
        A = np.sqrt(gamma * R * Tstat)
        self.M = self.V / A
    
    def setValues(self, P, T, M):
        self.P = P
        self.T = T
        A = np.sqrt(gamma * R * T)
        self.M = M
        self.V = M * A
        self.rho = P / (T * R)
        self.rhoV = self.rho * self.V
        self.rhoE = self.rho * (cv * T + 0.5*self.V*self.V)
        self.h0 = (self.rhoE + self.P) / self.rho
    
    def logStates(self, time):
        self.time.append(time)
        self.pHistory.append(self.P)
        self.VHistory.append(self.V)
        self.THistory.append(self.T)
        self.rhoHistory.append(self.rho)
        self.rhoVHistory.append(self.rhoV)
        self.rhoEHistory.append(self.rhoE)

    
volumes = [Volume(), Volume(), Volume(), Volume()]

for v_i, v in enumerate(volumes):
    pLocal = P_in + (P_out - P_in) * (v_i/(len(volumes)-1))
    v.setValues(pLocal, T0_in, 0.1)

volumes[0].setValues(P_in, T0_in, 0.5)


volumes[-1].P = P_out

for v in volumes:
    print(v.P, v.T, v.V, v.M, "\t", v.rho, v.rhoV)
print("")

t = 0
dt = 0.25 * 1e-2 / 500


while t < dt*2000:
    # Apply inlet boundary condition
    volumes[0].rho = (0.25 * volumes[0].rho) + (0.75 * rho_in)
    volumes[0].rho = min(volumes[0].rho, rho_in)
    tstat = T0_in * np.power(volumes[0].rho/rho_in, (gamma-1)/gamma)
    vel = np.sqrt(2*cp*(T0_in-tstat))
    eke = cv*tstat + 0.5*vel**2
    volumes[0].P = volumes[0].rho*R*tstat
    volumes[0].rhoV = volumes[0].rho * vel
    volumes[0].rhoE = volumes[0].rho * eke
    volumes[0].calculateProperties()

    volumes[-1].P = P_out

    for v_i, v in enumerate(volumes[1:]):
        mFluxIn = volumes[v_i-1].areaOut * volumes[v_i-1].rhoV
        mFluxOut = volumes[v_i].areaOut * volumes[v_i].rhoV
        volumes[v_i].rhoFlux += mFluxIn - mFluxOut

        volumes[v_i].rhoVFlux += (mFluxIn * volumes[v_i-1].V) + (volumes[v_i-1].P * volumes[v_i-1].areaOut)
        volumes[v_i].rhoVFlux -= (mFluxOut * volumes[v_i].V) + (volumes[v_i].P * volumes[v_i].areaOut)

        volumes[v_i].rhoEFlux += (mFluxIn * volumes[v_i-1].h0)
        volumes[v_i].rhoEFlux -= (mFluxOut * volumes[v_i].h0)
    
    for v in volumes[1:3]:
        #print(v.P, v.T, v.rho, v.rhoFlux, v.rhoV, v.rhoVFlux)
        v.rho += v.rhoFlux * dt / v.volume
        v.rhoV += v.rhoVFlux * dt / v.volume
        v.rhoE += v.rhoEFlux * dt / v.volume

        v.calculateProperties()

        v.rhoFlux = 0
        v.rhoVFlux = 0
        v.rhoEFlux = 0
    
    for v in volumes:
        v.logStates(t)
    
    dt = 0.2 * 1e-2 / max([v.V for v in volumes])


    t += dt 

import matplotlib.pyplot as plt
for v_i, v in enumerate(volumes):
    plt.plot(v.time, v.pHistory, label=str(v_i))

plt.legend()
#plt.ylim(-100, 1000)

plt.show()