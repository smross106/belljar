from CoolProp.CoolProp import PropsSI
import numpy as np
import crane_scratch2 as cs2


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
    
    def computedHdt(self, time):
        self.dhdt = self.storedE/(self.heatCapacity)
        self.storedE = 0

    def setNextH(self, dt):
        self.hHistory.append(self.hnm1_temp)
        self.dhHistory.append((self.hnp1 - self.hnm1_temp)/dt)
        if not self.boundaryConditionh:
            self.update_ph(self.p, self.hnp1)
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
    
    def computedmdt(self, time):
        self.dmdt = self.dmIn - self.dmOut
        self.stored_dm = 0
        self.dmIn = 0
        self.dmOut = 0
    
    def setNextm(self, dt):
        self.mHistory.append(self.m_nm1_temp)
        self.dmHistory.append((self.m_np1-self.m_nm1_temp)/dt)
        if not self.boundaryConditionm:
            self._m = self.m_np1 * 1.0
        self.m_np1 = 0
        
        self.intermediatemHistory = []
        self.intermediatedmdtHistory = []
        

class FluidStation(object):
    def __init__(self, volume=0):
        self.state = None
        self.flow = None
    
        self.nVarHistory = 0
        self.nDerivHistory = 0
        self.volume = volume
    
    def computeDerivatives(self, time):
        self.state.computedHdt(time)
        self.state.intermediatedHdtHistory.append(self.state.dHdt)

        self.flow.computedmdt(time)
        self.flow.intermediatedmdtHistory.append(self.flow.dmdt)
    
    def backupVariables(self):
        self.state.Tnm1_temp = self.state.T * 1.0
        self.flow.m_nm1_temp = self.flow.m * 1.0
    
    def updateVariables(self, variableHistoryValues={}, intermediateDerivativeHistoryValues={}, 
                        derivativeHistoryValues={},
                        updatePlusOne=False, returnValue=False, prevValue=1):
        nextH = self.state.Tnm1_temp * prevValue
        nextm = self.flow.m_nm1_temp * prevValue

        for key in variableHistoryValues:
            nextH += self.state.THistory[key] * variableHistoryValues[key]
        for key in intermediateDerivativeHistoryValues:
            nextH += self.state.intermediatedTdtHistory[key] * intermediateDerivativeHistoryValues[key]
        for key in derivativeHistoryValues:
            nextH += self.state.dTHistory[key] * derivativeHistoryValues[key]
        
        if returnValue:
            return([nextH, nextm])
        else:
            if updatePlusOne:
                self.state.h_np1 = nextH * 1.0
                self.flow.m_np1 = nextm * 1.0
            else:
                self.state.update_ph(self.state.p, nextH)
                self.flow.m = nextm * 1.0
    
    def setNextVar(self, dt):
        self.state.setNexth(dt)
        self.flow.setNextm(dt)

        self.nVarHistory += 1
        self.nDerivHistory += 1
    
    def calculateError(self, correctedValues):
        errors = []
        errors.append(abs(self.state.hnp1 - correctedValues[0])/self.state.hnm1_temp)
        errors.append(abs(self.flow.m_np1 - correctedValues[1])/self.flow.m_nm1_temp)
        return(errors)
    
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

        self.outlet = FluidStation()

        self.linkedStations = [self.outlet]
    
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

        self.inlet = FluidStation()
        self.outlet = FluidStation()
        
        self.linkedStations = [self.inlet, self.outlet]
    
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

        self.inlet = FluidStation()

        self.linkedStations = [self.inlet]
    
    def __repr__(self):
        return(self.compType+ str(self.inlet))


cs = cs2.SolverSolid()

def connectStationsFluid(self, component1, station1, component2, station2, fluid, volume):
    fs = FluidState(fluid)
    fl = FluidFlow(volume)

    setattr(getattr(self.components[component1], station1), "state", fs)
    setattr(getattr(self.components[component2], station2), "state", fs)

    setattr(getattr(self.components[component1], station1), "flow", fl)
    setattr(getattr(self.components[component2], station2), "flow", fl)

cs.connectStationsFluid = connectStationsFluid.__get__(cs)

def flowvar(time):
    return(0.5 + 0.5*np.sin(time))

cs.components["so"] = FluidSource_TP(flowvar, 1, 298)
cs.components["weir"] = VolumeLimitedOutflow(0.5)
cs.components["si"] = FluidSink()

cs.connectStationsFluid("so", "outlet", "weir", "inlet", "Nitrogen", 1)
cs.connectStationsFluid("weir", "outlet", "si", "inlet", "Nitrogen", 1)

cs.initialiseSolver(298)

for station in cs.stations:
    print(station)

