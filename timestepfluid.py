from CoolProp.CoolProp import PropsSI
import numpy as np
from crane_scratch2 import *


class FluidState(object):
    def __init__(self, fluid):
        self.fluid = fluid
        self.components = [self.fluid]
        self.ratios = [1]

        self.updated = False

        self.boundaryConditionT = False
        self.boundaryConditionh = False
        self.boundaryConditionp = False
        self.boundaryConditionrho = False

        self.storedE = 0
        self.dhdt = 0
        self.hHistory = []
        self.dhHistory = []

        self.intermediatehHistory = []
        self.intermediatedhdtHistory = []
        self.hnp1 = 0
        self.hnm1_temp = 0

        self._T = 298
        self._p = 1
        self._h = 0
        self._rho = 0
        self.update_Tp(self.T, self.p)

    @property
    def T(self):
        return self._T
    
    @T.setter
    def T(self, value):
        if not self.boundaryConditionT:
            self._T = value
        else:
            self._T = self.boundaryConditionT
    
    @property
    def h(self):
        return self._h
    
    @h.setter
    def h(self, value):
        if not self.boundaryConditionh:
            self._h = value
        else:
            self._h = self.boundaryConditionh
    
    @property
    def p(self):
        return self._p
    
    @p.setter
    def p(self, value):
        if not self.boundaryConditionp:
            self._p = value
        else:
            self._p = self.boundaryConditionp
    
    @property
    def rho(self):
        return self._rho
    
    @p.setter
    def rho(self, value):
        if not self.boundaryConditionrho:
            self._rho = value
        else:
            self._rho = self.boundaryConditionrho
    
    def stateInitialise(self, TRef, pRef):
        self.pHistory = []
        self.intermediatepHistory = []
        self.intermediatedpdtHistory = []
        self.pnm1_temp = 0
        self.pnp1 = 0 

        self.hHistory = []
        self.intermediatehHistory = []
        self.intermediatedhdtHistory = []
        self.hnm1_temp = 0
        self.hnp1 = 0 

        self._T = TRef
        self._p = pRef
        self.storedE = 0
    
    def setNexth(self, dt):
        self.hHistory.append(self.hnm1_temp)
        self.dhHistory.append((self.hnp1 - self.hnm1_temp)/dt)
        if not self.boundaryConditionh:
            self._h = self.hnp1 * 1.0
        self.update_ph(self.p, self._h)
        self.hnp1 = 0
        self.intermediatehHistory = []
        self.intermediatedhdtHistory = []
    
    def stateInitialise(self, hRef):
        self.hHistory = []
        self.intermediatehHistory = []
        self.intermediatedhdtHistory = []
        self.hnm1_temp = 0
        self.hnp1 = 0 
        self._h = hRef
        self.storedE = 0
    
    def update_Tp(self, T_k, p_bar):
        self.T = T_k
        self.p = p_bar

        # Check for saturation if T < T_triple
        var1 = "T"
        var2 = "P"
        if self.T < PropsSI("Tcrit", self.fluid):
            p_sat = PropsSI("P", "T", self.T, "Q", 0, self.fluid)/1e5
            if abs(p_sat/p_bar - 1) < 1e-6:
                var1 = "T|liquid"

        try:
            self.h = PropsSI("Hmass", var1, self.T, var2, self.p*1e5, self.fluid)
            self.s = PropsSI("Smass", var1, self.T, var2, self.p*1e5, self.fluid)
            self.rho = PropsSI("Dmass", var1, self.T, var2, self.p*1e5, self.fluid)
            self.Q = PropsSI("Q", var1, self.T, var2, self.p*1e5, self.fluid)
            self.b = self.h - (self.s * 298)
        except:
            raise ValueError("T {:.1f}K p{:.2f}bar fluid ".format(self.T, self.p)+self.fluid)
    
    def update_ps(self, p_bar, s):
        self.p = p_bar
        self.s = s
        
        try:
            self.h = PropsSI("Hmass", "Smass", self.s, "P", self.p*1e5, self.fluid)
            self.T = PropsSI("T", "Smass", self.s, "P", self.p*1e5, self.fluid)
            self.rho = PropsSI("Dmass", "Smass", self.s, "P", self.p*1e5, self.fluid)
            self.Q = PropsSI("Q", "Smass", self.s, "P", self.p*1e5, self.fluid)
            self.b = self.h - (self.s * 298)
        except:
            raise ValueError("s {:.1f}J/kgK p{:.2f}bar fluid ".format(self.s, self.p)+self.fluid)

    def update_hs(self, h, s):
        self.h = h
        self.s = s
        
        try:
            self.p = PropsSI("P", "Smass", self.s, "Hmass", self.h, self.fluid)/1e5
            self.T = PropsSI("T", "Smass", self.s, "Hmass", self.h, self.fluid)
            self.rho = PropsSI("Dmass", "Smass", self.s, "Hmass", self.h, self.fluid)
            self.Q = PropsSI("Q", "Smass", self.s, "Hmass", self.h, self.fluid)
            self.b = self.h - (self.s * 298)
        except:
            raise ValueError("s {:.1f}J/kgK h{:.1f}kJ/kg fluid ".format(self.s, self.h/1e3)+self.fluid)
    
    def update_ph(self, p_bar, h):
        self.p = p_bar
        self.h = h
        
        try:
            self.s = PropsSI("Smass", "Hmass", self.h, "P", self.p*1e5, self.fluid)
            self.T = PropsSI("T", "Hmass", self.h, "P", self.p*1e5, self.fluid)
            self.rho = PropsSI("Dmass", "Hmass", self.h, "P", self.p*1e5, self.fluid)
            self.Q = PropsSI("Q", "Hmass", self.h, "P", self.p*1e5, self.fluid)
            self.b = self.h - (self.s * 298)
        except:
            pass#raise ValueError("h {:.0f}J/kg p{:.2f}bar fluid ".format(self.s, self.p)+self.fluid)

    def update_Trho(self, T_k, rho):
        self.T = T_k
        self.rho = rho

        try:
            self.h = PropsSI("Hmass", "T", self.T, "Dmass", self.rho, self.fluid)
            self.s = PropsSI("Smass", "T", self.T, "Dmass", self.rho, self.fluid)
            self.p = PropsSI("P", "T", self.T, "Dmass", self.rho, self.fluid)/1e5
            self.Q = PropsSI("Q", "T", self.T, "Dmass", self.rho, self.fluid)
            self.b = self.h - (self.s * 298)
        except:
            raise ValueError("T {:.1f}K rho{:.2f}kg/m3 fluid ".format(self.T, self.rho)+self.fluid)

    def update_pQ(self, p_bar, Q):
        self.p = p_bar
        self.Q = Q
        
        try:
            self.s = PropsSI("Smass", "Q", self.Q, "P", self.p*1e5, self.fluid)
            self.T = PropsSI("T", "Q", self.Q, "P", self.p*1e5, self.fluid)
            self.rho = PropsSI("Dmass", "Q", self.Q, "P", self.p*1e5, self.fluid)
            self.h = PropsSI("Hmass", "Q", self.Q, "P", self.p*1e5, self.fluid)
            self.b = self.h - (self.s * 298)
        except:
            raise ValueError("Q {:.2f} p{:.2f}bar fluid ".format(self.Q, self.p)+self.fluid)

class FluidFlow(object):
    def __init__(self, volume):
        self.mDot = 0

        self.volume = volume

        self._m = 0
        self.dmIn = 0
        self.dmOut = 0
        self.stored_dm = 0
        self.dmdt = 0
        self.mHistory = []
        self.dmHistory = []

        self.boundaryConditionm = False

        self.intermediatemHistory = []
        self.intermediatedmdtHistory = []
        self.m_np1 = 0
        self.m_nm1_temp = 0
    
    @property
    def m(self):
        return self._m
    
    @m.setter
    def m(self, value):
        if not self.boundaryConditionm:
            self._m = value
        else:
            self._m = self.boundaryConditionm
    
    def flowInitialise(self, mDotRef):
        self.mHistory = []
        self.dmHistory = []
        self.intermediatemHistory = []
        self.intermediatedmdtHistory = []
        self.m_np1 = 0
        self.m_nm1_temp = 0
        self._mDot = mDotRef
        self.dmIn = 0
        self.dmOut = 0
        self.stored_dm = 0

class Station(object):
    def __init__(self, volume):
        self.state = None
        self.flow = None
    
    def __repr__(self):
        txt = " Fluid:{}, mDot {:.2f} P {:.1f}bar T{:.0f}K".format(self.state.fluid, self.flow.mDot, self.state.p, self.state.T)
        return(txt)
    
class CycleComponent(object):
    def __init__(self):
        self.compType = "None"
    
    def compute(self):
        pass

class FluidSource_TP(CycleComponent):
    def __init__(self, mDot, p, T, W=None, RH=None):
        super().__init__()
        self.compType = "Source"

        self.mDot = mDot
        self.p = p
        self.T = T

        self.outlet = Station()

        self.stations = [self.outlet]
    
    def compute(self, time):
        if callable(self.mDot):
            mDot = self.mDot(time)
        else:
            mDot = self.mDot * 1.0

        if callable(self.mDot):
            T = self.T(time)
        else:
            T = self.T * 1.0

        if callable(self.p):
            p = self.p(time)
        else:
            p = self.p * 1.0
        
        self.outlet.flow.dmIn += mDot

        self.outlet.state.update_Tp(T, p)
        self.outlet.state.storedE = self.outlet.state.h - self.outlet.state.hnm1_temp
    
    def __repr__(self):
        return(self.compType+ str(self.outlet))

class VolumeLimitedOutflow(CycleComponent):
    def __init__(self, maxOutflow):
        super().__init__()
        self.compType = "Volume"

        self.maxOutflow = maxOutflow

        self.inlet = Station()
        self.outlet = Station()
        
        self.stations = [self.inlet, self.outlet]
    
    def compute(self, time):
        if callable(self.maxOutflow):
            maxOutflow = self.maxOutflow(time)
        else:
            maxOutflow = self.maxOutflow

        self.outlet.state.update_Tp(self.inlet.state.T, self.inlet.state.p)
        self.outlet.state.storedE = self.outlet.state.h - self.outlet.state.hnm1_temp

        self.outlet.flow.dmIn += min(maxOutflow, self.inlet.flow.mDot)

    def __repr__(self):
        return(self.compType+ str(self.inlet))

class FluidSink(CycleComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Sink"

        self.inlet = Station()

        self.stations = [self.inlet]
    
    def __repr__(self):
        return(self.compType+ str(self.inlet))


class SolverMixedConstantPressure(object):
    def __init__(self, precision=4, fluid="Nitrogen"):
        self.fluid = fluid
        self.components = {}
        self.stations = []
        self.states = []
        self.flows = []

        self.precision = precision
        self.errors = []
        self.time = []
        self.allowTimeAdjust = True

        self.mode = "RK4GILL"

    def connectStations(self, component1, station1, component2, station2, solid=False, mass=1, specificHeatCapacity=1, TRef=298, volume=1):
        
        if solid:
            fs = SolidState(mass=mass, specificHeatCapacity=specificHeatCapacity, TRef=TRef)
        else:
            fs = FluidState(self.fluid)
            fl = FluidFlow(volume)

        self.states.append(fs)
        setattr(getattr(self.components[component1], station1), "state", fs)
        setattr(getattr(self.components[component2], station2), "state", fs)

        if not self.solid:
            self.flows.append(fl)
            setattr(getattr(self.components[component1], station1), "flow", fl)
            setattr(getattr(self.components[component2], station2), "flow", fl)

    
    def initialiseSolver(self, TRef, pRef, mDotRef):
        for component_key in self.components.keys():
            for station in self.components[component_key].linkedStations:
                station.state.stateInitialise(TRef, pRef)
                if hasattr(station, "flow"):
                    station.flow.flowInitialise(mDotRef)
                if station not in self.stations:
                    self.stations.append(station)

        self.errors = []
        self.time = []
    
    def calculateDerivatives(self, time):
        for component_key in self.components.keys():
            self.components[component_key].compute(time)
        
        for state in self.states:
            if type(state) == SolidState:
                state.computedTdt(time)
            elif type(state)
    
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
            if type(state) == FluidState:
                state.hnm1_temp = state.h * 1.0
            elif type(state) == SolidState:
                state.Tnm1_temp = state.T * 1.0
        for flow in self.flows:
            flow.m_nm1_temp = flow.m * 1.0    
    
        # Step 0 - generates k0
        self.calculateDerivatives(time)
        for state in self.states:
            if type(state) == SolidState:
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
            
            state.T = state.Tnp1
            
            maxError = max(abs(state.Tnp1- Tnp1_bar)/state.Tnm1_temp, maxError)
            
            #print(state.Tnp1, Tnp1_bar, abs(state.Tnp1- Tnp1_bar)/state.THistory[-1])
        #print(maxError)
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


