import belljar_lockstep as bjls
import numpy as np
from CoolProp.CoolProp import PropsSI

class CascadeFluidState(object):
    def __init__(self):
        self.ID = ""

        self.rho  = 1
        self.rhoV = 1
        self.rhoE = 1

        self.volume = 1e-1
        self.areaIn = 1e-1
        self.areaOut = 1e-1

        self.fluidProp = "simple"
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
        self.drhodt  = self.rhoFlux  / self.volume
        self.drhoVdt = self.rhoVFlux / self.volume
        self.drhoEdt = self.rhoEFlux / self.volume

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
            if self.T < 0:
                print(self.V, self.rho, self.rhoV, self.rhoE, self.T, self.P)
                raise ValueError
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

class CascadeFluidStation(object):
    def __init__(self):
        self.state = None

        self.nVarHistory = 0
        self.nDerivHistory = 0

        self.countErrors = True

        self.ID = ""
    
    def computeDerivatives(self, time, dt):
        self.state.computedVardt(time, dt)

        self.state.intermediatedrhodtHistory.append( self.state.drhodt)
        self.state.intermediatedrhoVdtHistory.append(self.state.drhoVdt)
        self.state.intermediatedrhoEdtHistory.append(self.state.drhoEdt)
    
    def backupVariables(self):
        self.state.rho_nm1 =  self.state.rho * 1.0
        self.state.rhoV_nm1 = self.state.rhoV * 1.0
        self.state.rhoE_nm1 = self.state.rhoE * 1.0
    
    def updateVariables(self, variableHistoryValues={}, intermediateDerivativeHistoryValues={}, 
                        derivativeHistoryValues={},
                        updatePlusOne=False, returnValue=False, prevValue=1):
        nextrho =  self.state.rho_nm1  * prevValue
        nextrhoV = self.state.rhoV_nm1 * prevValue
        nextrhoE = self.state.rhoE_nm1 * prevValue

        for key in variableHistoryValues:
            nextrho  += self.state.rhoHistory[key]  * variableHistoryValues[key]
            nextrhoV += self.state.rhoVHistory[key] * variableHistoryValues[key]
            nextrhoE += self.state.rhoEHistory[key] * variableHistoryValues[key]
        for key in intermediateDerivativeHistoryValues:
            nextrho  += self.state.intermediatedrhodtHistory[key]  * intermediateDerivativeHistoryValues[key]
            nextrhoV += self.state.intermediatedrhoVdtHistory[key] * intermediateDerivativeHistoryValues[key]
            nextrhoE += self.state.intermediatedrhoEdtHistory[key] * intermediateDerivativeHistoryValues[key]

        for key in derivativeHistoryValues:
            nextrho  += self.state.drhodtHistory[key]  * derivativeHistoryValues[key]
            nextrhoV += self.state.drhoVdtHistory[key] * derivativeHistoryValues[key]
            nextrhoE += self.state.drhoEdtHistory[key] * derivativeHistoryValues[key]
        
        if returnValue:
            return([nextrho, nextrhoV, nextrhoE])
        else:
            if updatePlusOne:
                self.state.rho_np1  = nextrho  * 1.0
                self.state.rhoV_np1 = nextrhoV * 1.0
                self.state.rhoE_np1 = nextrhoE * 1.0
            else:
                self.state.rho  = nextrho * 1.0
                self.state.rhoV = nextrhoV * 1.0
                self.state.rhoE = nextrhoE * 1.0
                self.state.setOtherProperties()
    
    def setNextVar(self, dt):
        self.state.setNextVar(dt)

        self.nVarHistory += 1
        self.nDerivHistory += 1
    
    def resetNextVar(self):
        self.state.resetSolver()
    
    def calculateErrorValueInput(self, correctedValues):
        errors = []
        errors.append(abs(self.state.rho_np1  - correctedValues[0])/self.state.rho_nm1 )
        errors.append(abs(self.state.rhoV_np1 - correctedValues[1])/self.state.rhoV_nm1)
        errors.append(abs(self.state.rhoE_np1 - correctedValues[2])/self.state.rhoE_nm1)
        if self.countErrors:
            return(errors)
        else:
            return(0)

    def calculateErrorPreviousValues(self):
        errors = []
        errors.append(abs(self.state.rho_np1  - self.state.rho_nm1 )/self.state.rho_nm1 )
        errors.append(abs(self.state.rhoV_np1 - self.state.rhoV_nm1)/self.state.rhoV_nm1)
        errors.append(abs(self.state.rhoE_np1 - self.state.rhoE_nm1)/self.state.rhoE_nm1)
        if self.countErrors:
            return(errors)
        else:
            return(0)

class BoundaryConditionUpstream(object):
    def __init__(self, M, T, P):
        self.M = M
        self.T = T
        self.P = P

        self.fluidOutlet = CascadeFluidStation()
        self.linkedStations = [self.fluidOutlet]

        self.refState = CascadeFluidState()
        self.refState.setState(M = self.M, T = self.T, P = self.P)
    
    def compute(self, time):
        self.fluidOutlet.state.rhoFlux  += self.refState.V * self.refState.rho

        self.fluidOutlet.state.rhoVFlux += -0.75*self.fluidOutlet.state.rhoV + 0.75*self.refState.V

        #self.fluidOutlet.state.rhoVFlux += self.refState.V * self.refState.rho * self.refState.V
        #self.fluidOutlet.state.rhoVFlux += self.fluidOutlet.state.areaIn * (self.refState.P - self.fluidOutlet.state.P) * 0.25
        
        self.fluidOutlet.state.rhoEFlux += self.refState.V * self.refState.rho * self.refState.h0

class BoundaryConditionDownstream(object):
    def __init__(self, M, T, P):
        self.M = M
        self.T = T
        self.P = P

        self.fluidInlet = CascadeFluidStation()
        self.linkedStations = [self.fluidInlet]

        self.refState = CascadeFluidState()
        self.refState.setState(M = self.M, T = self.T, P = self.P)
    
    def compute(self, time):
        self.fluidInlet.state.rhoFlux  -= self.fluidInlet.state.V * self.fluidInlet.state.rho

        self.fluidInlet.state.rhoVFlux -= self.fluidInlet.state.V * self.fluidInlet.state.rho * self.fluidInlet.state.V
        self.fluidInlet.state.rhoVFlux += self.fluidInlet.state.areaOut * (self.fluidInlet.state.P - self.refState.P)
        
        self.fluidInlet.state.rhoEFlux -= self.fluidInlet.state.V * self.fluidInlet.state.rho * self.fluidInlet.state.h0


class FlowObject(object):
    def __init__(self):
        self.fluidInlet = CascadeFluidStation()
        self.fluidOutlet = CascadeFluidStation()

        self.linkedStations = [self.fluidInlet, self.fluidOutlet]
    
    def compute(self, time):
        self.fluidInlet.state.rhoFlux   -= self.fluidOutlet.state.V * self.fluidOutlet.state.rho
        self.fluidOutlet.state.rhoFlux  += self.fluidInlet.state.V  * self.fluidInlet.state.rho

        self.fluidInlet.state.rhoVFlux  -= self.fluidOutlet.state.V * self.fluidOutlet.state.rho * self.fluidOutlet.state.V
        self.fluidOutlet.state.rhoVFlux += self.fluidInlet.state.V  * self.fluidInlet.state.rho  * self.fluidInlet.state.V

        self.fluidInlet.state.rhoVFlux  -= self.fluidInlet.state.areaOut * (self.fluidInlet.state.P - self.fluidOutlet.state.P)
        self.fluidOutlet.state.rhoVFlux += self.fluidOutlet.state.areaIn * (self.fluidOutlet.state.P - self.fluidOutlet.state.P)
        
        self.fluidInlet.state.rhoEFlux  -= self.fluidOutlet.state.V * self.fluidOutlet.state.rho * self.fluidOutlet.state.h0
        self.fluidOutlet.state.rhoEFlux += self.fluidInlet.state.V  * self.fluidInlet.state.rho  * self.fluidInlet.state.h0



class CascadeTimestepSolver(bjls.TimestepSolver):
    def connectStationsFluidCascade(self, component1, station1, component2, station2):
        if getattr(getattr(self.components[component1], station1), "state") != None:
            fs = getattr(getattr(self.components[component1], station1), "state")
        elif getattr(getattr(self.components[component2], station2), "state") != None:
            fs = getattr(getattr(self.components[component2], station2), "state")
        else:
            fs = CascadeFluidState()

        setattr(getattr(self.components[component1], station1), "state", fs)
        setattr(getattr(self.components[component2], station2), "state", fs)
    
    def initialiseSolverCascade(self):
        for component_key in self.components.keys():
            for station in self.components[component_key].linkedStations:
                if station.state == None:
                    pass
                else:
                    if type(station.state) == bjls.SolidState:
                        station.state.stateInitialise(298)
                        if station.state.ID == "":
                            station.state.ID = component_key
                    else:
                        station.state.setState(M=0.1, T=298, P=1e5)
                        if station.state.ID == "":
                            station.state.ID = component_key

                    if station.state not in self.states:
                        self.stationsUnique.append(station)
                        if type(station) == bjls.SolidStation:
                            self.stations.append([station, station.state, None])
                        else:
                            self.stations.append([station, station.state, None])
                        self.states.append(station.state)

solver = CascadeTimestepSolver(precision=5)
solver.mode = "RK4ERROR"

solver.components["source"] = BoundaryConditionUpstream(M=0.1, T=298, P=1.01e5)
solver.components["f1"] = FlowObject()
solver.components["f2"] = FlowObject()
solver.components["sink"] = BoundaryConditionDownstream(M=0.1, T=298, P=1e5)

solver.connectStationsFluidCascade("source", "fluidOutlet", "f1", "fluidInlet")
solver.connectStationsFluidCascade("f1", "fluidOutlet", "f2", "fluidInlet")
solver.connectStationsFluidCascade("f2", "fluidOutlet", "sink", "fluidInlet")

solver.initialiseSolverCascade()

solver.runSolver(0, 1, 0.001)

print(solver.components["source"].refState.P, solver.components["source"].refState.M)
print(solver.components["f1"].fluidInlet.state.P, solver.components["f1"].fluidInlet.state.M)
print(solver.components["f1"].fluidOutlet.state.P, solver.components["f1"].fluidOutlet.state.M)
print(solver.components["f2"].fluidInlet.state.P, solver.components["f2"].fluidInlet.state.M)
print(solver.components["f2"].fluidOutlet.state.P, solver.components["f2"].fluidOutlet.state.M)

import matplotlib.pyplot as plt

plt.plot(solver.time, [solver.components["f1"].fluidInlet.state.rhoVHistory[t]/solver.components["f1"].fluidInlet.state.rhoHistory[t] for t in range(len(solver.time))])
plt.plot(solver.time, [solver.components["f2"].fluidInlet.state.rhoVHistory[t]/solver.components["f2"].fluidInlet.state.rhoHistory[t] for t in range(len(solver.time))])

plt.show()

# x1 = CascadeFluidState()
# x1.rhoV = 330
# x1.rhoE = 330**2 + 1000*200

# x1.setState(M=0.1, T=2000, P=1e5)
# x1.setOtherProperties()
# print("")
# x1.rhoV += 100
# x1.setOtherProperties()

