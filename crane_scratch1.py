import numpy as np
import time as pytime
from CoolProp.CoolProp import PropsSI


class SolidState(object):
    def __init__(self, mass=1, specificHeatCapacity=1, TRef=298):
        self._T = TRef
        self.mass = mass
        self.specificHeatCapacity = specificHeatCapacity
        self.heatCapacity = self.mass * self.specificHeatCapacity

        self.boundaryConditionT = False

        self.storedE = 0
        self.dTdt = 0
        self.THistory = []
        self.dTHistory = []

        self.intermediateTHistory = []
        self.intermediatedTdtHistory = []
        self.Tnp1 = 0
        self.Tnm1_temp = 0
    
    @property
    def T(self):
        return self._T
    
    @T.setter
    def T(self, value):
        if not self.boundaryConditionT:
            self._T = value
        else:
            self._T = self.boundaryConditionT

        
    def computedTdt(self, time):
        self.dTdt = self.storedE/(self.heatCapacity)
        self.storedE = 0

    def setNextT(self, dt):
        self.THistory.append(self.Tnm1_temp)
        self.dTHistory.append((self.Tnp1 - self.Tnm1_temp)/dt)
        if not self.boundaryConditionT:
            self._T = self.Tnp1 * 1.0
        self.Tnp1 = 0
        self.intermediateTHistory = []
        self.intermediatedTdtHistory = []
    
    def stateInitialise(self, TRef):
        self.THistory = []
        self.intermediateTHistory = []
        self.intermediatedTdtHistory = []
        self.Tnm1_temp = 0
        self.Tnp1 = 0 
        self._T = TRef
        self.storedE = 0

class SolidStation(object):
    def __init__(self):
        self.state = None

class SolidComponent(object):
    def __init__(self):
        self.linkedStations = []
    
    def compute(self):
        pass

class ThermalConductance(SolidComponent):
    def __init__(self, C=1):
        super().__init__()
        self.C = C
        # Thermal conductance, W/K
        
        self.inlet = SolidStation()
        self.outlet = SolidStation()

        self.linkedStations = [self.inlet, self.outlet]
    
    def compute(self, time):
        dT = self.inlet.state.T - self.outlet.state.T

        self.inlet.state.storedE -= dT * self.C
        self.outlet.state.storedE += dT * self.C

class ConstantTemperatureBoundary(SolidComponent):
    def __init__(self, TSet):
        super().__init__()

        self.TSet = TSet

        self.outlet = SolidStation()

        self.linkedStations = [self.outlet]
    
    def compute(self, time):
        if callable(self.TSet):
            TSet = self.TSet(time)
        else:
            TSet = self.TSet * 1.0
        self.outlet.state.boundaryConditionT = TSet
        dT = TSet - self.outlet.state.T

        self.outlet.state.storedE += dT * self.outlet.state.heatCapacity

class HeatInputBoundary(SolidComponent):
    def __init__(self, QIn):
        super().__init__()

        self.QIn = QIn

        self.outlet = SolidStation()

        self.linkedStations = [self.outlet]
    
    def compute(self, time):
        if callable(self.QIn):
            QIn = self.QIn(time)
        else:
            QIn = self.QIn * 1.0

        self.outlet.state.storedE += QIn

class SolverSolid(object):
    def __init__(self, precision=4):
        self.components = {}
        self.stations = []
        self.states = []

        self.precision = precision
        self.errors = []
        self.time = []
        self.allowTimeAdjust = True

        self.mode = "RK4GILL"

    def connectStations(self, component1, station1, component2, station2, mass=1, specificHeatCapacity=1, TRef=298):
            
        fs = SolidState(mass=mass, specificHeatCapacity=specificHeatCapacity, TRef=TRef)

        self.states.append(fs)
        setattr(getattr(self.components[component1], station1), "state", fs)
        setattr(getattr(self.components[component2], station2), "state", fs)
    
    def initialiseSolver(self, Tref):
        for component_key in self.components.keys():
            for station in self.components[component_key].linkedStations:
                station.state.stateInitialise(Tref)
                if station not in self.stations:
                    self.stations.append(station)

        self.errors = []
        self.time = []
    
    def calculateDerivatives(self, time):
        for component_key in self.components.keys():
            self.components[component_key].compute(time)
        
        for state in self.states:
            state.computedTdt(time)
    
    def changeTimestep(self, dt, errorFailed):
        if self.allowTimeAdjust:
            if errorFailed:
                return(dt * 0.625)
            else:
                return(dt * 2)
        else:
            return(dt)

    def runSolver(self, t0, t1, dtTrial):
        time = t0
        dt = dtTrial
        minRunsBetweenTimestepIncrease = 3
        runsSinceTimestepIncrease = 0
        numIters = 0
        runsSinceLastError = 0

        while time < t1:
            #print(time, dt)
            if self.mode == "RK4ERROR":
                success, error = self.RK4OneStepError(time, dt)
                if success:
                    for state in self.states:
                        state.setNextT(dt)
                    time += dt
                    self.time.append(time)
                    self.errors.append(error)

                    runsSinceTimestepIncrease += 1
                    if runsSinceTimestepIncrease >= minRunsBetweenTimestepIncrease:
                        dt = self.changeTimestep(dt, False)
                        runsSinceTimestepIncrease = 0
                else:
                    dt = self.changeTimestep(dt, True)
            elif self.mode == "RK4GILL":

                success, error = self.RK4OneStepGill(time, dt)
                for state in self.states:
                    state.setNextT(dt)
                time += dt
                self.time.append(time)
                self.errors.append(np.nan)

            elif self.mode == "CRANE":
                pass
                #success, error = self.PredictorCorrectorCraneError(time, dt)

            elif self.mode == "RK4CRANE":
                # Crane algorithm, falling back to RK4 with the previous timestep
                if runsSinceLastError < 3:
                    # We don't have enough previous dT/dt values to use the Crane method - fall back to RK4
                    success, error = self.RK4OneStepGill(time, dt)
                    if success:
                        runsSinceLastError += 1
                        for state in self.states:
                            state.setNextT(dt)
                        time += dt
                        self.time.append(time)
                        self.errors.append(error)

                        runsSinceTimestepIncrease += 1
                else:
                    # We can happily use the Crane algorithm
                    success, error = self.PredictorCorrectorCraneError(time, dt)
                    if success:
                        runsSinceLastError += 1
                        for state in self.states:
                            state.setNextT(dt)
                        time += dt
                        self.time.append(time)
                        self.errors.append(error)

                        runsSinceTimestepIncrease += 1

                        if runsSinceTimestepIncrease >= minRunsBetweenTimestepIncrease:
                            dt = self.changeTimestep(dt, False)
                            runsSinceTimestepIncrease = 0
                    else:
                        dt = self.changeTimestep(dt, True)

            else:
                break
            #print(time, success, error, dt)
            numIters += 1
            if numIters > 100000:
                break
    
    def RK4OneStepError(self, time, dt):
        # Back up previous temperatures to allow for RK4 to occur
        for state in self.states:
            state.Tnm1_temp = state.T * 1.0
        
        # Step 0 - generates k0
        self.calculateDerivatives(time)
        for state in self.states:
            state.intermediatedTdtHistory.append(state.dTdt)
            # Feed in T_k1
            state.T = state.Tnm1_temp + state.intermediatedTdtHistory[0]*dt/2

        # Step 1 - generates k1
        self.calculateDerivatives(time+dt/2)
        for state in self.states:
            state.intermediatedTdtHistory.append(state.dTdt)
            # Feed in T_k2
            state.T = state.Tnm1_temp + (state.intermediatedTdtHistory[0] + state.intermediatedTdtHistory[1])*dt/4
        
        # Step 2 - generates k2
        self.calculateDerivatives(time+dt/2)
        for state in self.states:
            state.intermediatedTdtHistory.append(state.dTdt)
            # Feed in T_k3
            state.T = state.Tnm1_temp + (2*state.intermediatedTdtHistory[2] - state.intermediatedTdtHistory[1])*dt
        
        # Step 3 - generates k3
        self.calculateDerivatives(time+dt)
        for state in self.states:
            state.intermediatedTdtHistory.append(state.dTdt)
            # Feed in T_k4
            state.T = state.Tnm1_temp + (
                7*state.intermediatedTdtHistory[0] + 10*state.intermediatedTdtHistory[1] + state.intermediatedTdtHistory[3])*dt/27
        
        # Step 4 - generates k4
        self.calculateDerivatives(time+2*dt/3)
        for state in self.states:
            state.intermediatedTdtHistory.append(state.dTdt)
            # Feed in T_k5
            state.T = state.Tnm1_temp + (
                28*state.intermediatedTdtHistory[0] - 125*state.intermediatedTdtHistory[1] + 546*state.intermediatedTdtHistory[2]
                +54*state.intermediatedTdtHistory[3] - 378*state.intermediatedTdtHistory[4])*dt/625
        

        maxError = 0
        # Step 5 - generates k5
        self.calculateDerivatives(time+dt/5)
        for state in self.states:
            state.intermediatedTdtHistory.append(state.dTdt)

            # Calculate error
            state.Tnp1 = state.Tnm1_temp + (dt/6)*(
                state.intermediatedTdtHistory[0] + 4*state.intermediatedTdtHistory[2] + state.intermediatedTdtHistory[3])

            Tnp1_bar = state.Tnm1_temp + (dt/336)*(
                14*state.intermediatedTdtHistory[0] + 35*state.intermediatedTdtHistory[3] + 162*state.intermediatedTdtHistory[4]
                + 125*state.intermediatedTdtHistory[5])
            
            #state.T = state.Tnp1
            
            maxError = max(abs(state.Tnp1- Tnp1_bar)/state.Tnm1_temp, maxError)
            
            #print(state.Tnp1, Tnp1_bar, abs(state.Tnp1- Tnp1_bar)/state.THistory[-1])
        print(maxError, dt)
        if -np.log10(maxError) > self.precision:
            return(True, maxError)
        else:
            return(False, maxError)

    def RK4OneStepGill(self, time, dt):
        # Back up previous temperatures to allow for RK4 to occur
        for state in self.states:
            state.Tnm1_temp = state.T * 1.0
        
        # Step 0 - generates k0
        self.calculateDerivatives(time)
        for state in self.states:
            state.intermediatedTdtHistory.append(state.dTdt)
            # Feed in T_k1
            state.T = state.Tnm1_temp + state.intermediatedTdtHistory[0]*dt/2

        # Step 1 - generates k1
        self.calculateDerivatives(time+dt/2)
        for state in self.states:
            state.intermediatedTdtHistory.append(state.dTdt)
            # Feed in T_k2
            state.T = state.Tnm1_temp + ((state.intermediatedTdtHistory[0]*(-1+np.sqrt(2)) + 
                                          state.intermediatedTdtHistory[1]*(2-np.sqrt(2))))*dt/2
        
        # Step 2 - generates k2
        self.calculateDerivatives(time+dt/2)
        for state in self.states:
            state.intermediatedTdtHistory.append(state.dTdt)
            # Feed in T_k3
            state.T = state.Tnm1_temp + ((state.intermediatedTdtHistory[1]*(-np.sqrt(2)) + 
                                          state.intermediatedTdtHistory[2]*(1+np.sqrt(2))))*dt/2
        
        # Step 3 - generates k3
        self.calculateDerivatives(time+dt)
        for state in self.states:
            state.intermediatedTdtHistory.append(state.dTdt)
            # Calculate T_np1
            state.Tnp1 = state.Tnm1_temp + (dt/6)*(state.intermediatedTdtHistory[0] +
                                                    state.intermediatedTdtHistory[1] * (2 - np.sqrt(2)) +
                                                    state.intermediatedTdtHistory[2] * (2 + np.sqrt(2)) +
                                                    state.intermediatedTdtHistory[3]
            )
        return(True, 0)
        
    def PredictorCorrectorCraneError(self, time, dt):
        # Back up previous temperatures to allow for RK4 to occur
        for state in self.states:
            if len(state.dTHistory) < 3 or len(state.THistory) < 3:
                return(False, 1000)
            state.Tnm1_temp = state.T * 1.0
        
        # Carry out predictor step
        self.calculateDerivatives(time)
        for state in self.states:
            # This is y'_n
            state.intermediatedTdtHistory.append(state.dTdt)

            # This is p_n+1
            state.T = ( 1.54765200 * state.Tnm1_temp + 
                       -1.86750300 * state.THistory[-1] + 
                        2.01720400 * state.THistory[-2] + 
                       -0.69735300 * state.THistory[-3]) + dt*(
                            2.00224700 * state.intermediatedTdtHistory[-1] +
                           -2.03169000 * state.dTHistory[-1] +
                            1.81860900 * state.dTHistory[-2] +
                           -0.71432000 * state.dTHistory[-3]
                       )
        
        self.calculateDerivatives(time+dt)
        maxError = 0
        for state in self.states:
            # This is p'_n+1
            state.intermediatedTdtHistory.append(state.dTdt)

            state.Tnp1 = state.Tnm1_temp + dt * (
                 0.375000000 * state.intermediatedTdtHistory[-1] +
                 0.791666667 * state.intermediatedTdtHistory[-2] +
                -0.208333333 * state.dTHistory[-1] +
                 0.041666667 * state.dTHistory[-2]
            )
            state.T = state.Tnm1_temp

            maxError = max(abs(state.Tnp1- state.T)/(16.21966*state.Tnm1_temp), maxError)

        if -np.log10(maxError) > self.precision:
            return(True, maxError)
        else:
            return(False, maxError)

cs = SolverSolid()

"""cs.components["s1"] = ConstantTemperatureBoundary(302)
cs.components["s2"] = ConstantTemperatureBoundary(296)

cs.components["l1"] = ThermalConductance()
cs.components["l2"] = ThermalConductance()

cs.connectStations("s1", "outlet", "l1", "inlet", specificHeatCapacity=10)
cs.connectStations("l1", "outlet", "l2", "inlet", specificHeatCapacity=10)
cs.connectStations("l2", "outlet", "s2", "outlet", specificHeatCapacity=10)

cs.initialiseSolver()

#print(cs.RK4OneStepError(0, 0.01))
cs.runSolver(0, 10, 0.01)
"""

def hat(time):
    return(200+10*np.sin(time/5))

cs.components["s1"] = HeatInputBoundary(500)
cs.components["s2"] = ConstantTemperatureBoundary(hat)
cs.components["s3"] = ConstantTemperatureBoundary(300)
cs.components["l1"] = ThermalConductance(C=10)
cs.components["l2"] = ThermalConductance(C=20)
cs.components["l3"] = ThermalConductance(C=10)
cs.components["l4"] = ThermalConductance(C=10)


cs.connectStations("s1", "outlet", "l1", "inlet", specificHeatCapacity=10)
cs.connectStations("l1", "outlet", "l4", "inlet", specificHeatCapacity=10)
cs.connectStations("l4", "outlet", "s3", "outlet", specificHeatCapacity=10)
cs.connectStations("s2", "outlet", "l2", "inlet", specificHeatCapacity=10)
cs.connectStations("l2", "outlet", "l3", "inlet", specificHeatCapacity=10)
cs.connectStations("l1", "outlet", "l3", "outlet", specificHeatCapacity=10)

#cs.connectStations(cs.components["s1"].outlet, cs.components["l1"].inlet, 1, 100, 300)
#cs.connectStations(cs.components["l1"].outlet, cs.components["l4"].inlet, 1, 100, 300)
#cs.connectStations(cs.components["l4"].outlet, cs.components["s3"].outlet, 1, 100, 300)
#cs.connectStations(cs.components["s2"].outlet, cs.components["l2"].inlet, 1, 100, 300)
#cs.connectStations(cs.components["l2"].outlet, cs.components["l3"].inlet, 1, 100, 300)
#cs.connectStations(cs.components["l1"].outlet, cs.components["l3"].outlet, 1, 100, 300)

cs.initialiseSolver(298)
cs.precision = 4



print("crane")
cs.initialiseSolver(298)
print("crane")
cs.initialiseSolver(298)
cs.mode="RK4CRANE"

cs.runSolver(0, 1, 0.01)

import matplotlib.pyplot as plt

plt.plot(cs.time, cs.components["l1"].outlet.state.THistory, label="Outlet")
plt.plot(cs.time, cs.components["l2"].inlet.state.THistory, label="BC")
plt.plot(cs.time, cs.components["l2"].outlet.state.THistory, label="Inter")
plt.legend()

print("#")
print(cs.components["l2"].inlet.state.T)
print(cs.components["l3"].inlet.state.T)
print(cs.components["l3"].outlet.state.T)
print(cs.components["l4"].outlet.state.T)
print(cs.components["l1"].inlet.state.T)


plt.show()