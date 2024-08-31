from CoolProp.CoolProp import PropsSI


class FluidState(object):
    def __init__(self, fluid):
        self.fluid = fluid
        self.components = [self.fluid]
        self.ratios = [1]

        self.updated = False
    
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
    def __init__(self):
        self.mDot = 0

class Station(object):
    def __init__(self):
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

class TwoStationComponent(CycleComponent):
    def __init__(self):
        self.inlet = Station()
        self.outlet = Station()

        self.stations = [self.inlet, self.outlet]
    
    def __repr__(self):
        return(self.compType+ str(self.inlet))

class FluidSource_TP(CycleComponent):
    def __init__(self, mDot, p, T, W=None, RH=None):
        super().__init__()
        self.compType = "Source"

        self.mDot = mDot
        self.p = p
        self.T = T

        self.outlet = Station()

        self.stations = [self.outlet]
    
    def compute(self):
        self.outlet.flow.mDot = self.mDot

        self.outlet.state.update_Tp(self.T, self.p)
    
    def __repr__(self):
        return(self.compType+ str(self.outlet))

class FluidSink(CycleComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Sink"

        self.inlet = Station()

        self.stations = [self.inlet]
    
    def __repr__(self):
        return(self.compType+ str(self.inlet))

class Compressor(TwoStationComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Compressor"

        self.WIn = 0
    
    def compute(self, pOut, eta=1):
        #pOut = self.inlet.state.p * pressureRatio
        if pOut < self.inlet.state.p:
            raise ValueError("Compressor outlet pressure greater than inlet pressure", pOut, self.inlet.state.p)
        

        self.outlet.state.update_ps(pOut, self.inlet.state.s)

        dHIdeal = self.outlet.state.h - self.inlet.state.h 
        dHActual = dHIdeal / eta

        self.outlet.state.update_ph(pOut, self.inlet.state.h + dHActual)

        self.WIn = (self.outlet.state.h - self.inlet.state.h) * self.outlet.flow.mDot

        self.outlet.flow.mDot = self.inlet.flow.mDot

class Turbine(TwoStationComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Turbine"

        self.WIn = 0
    
    def compute(self, pOut, eta=1):
        if pOut > self.inlet.state.p:
            raise ValueError("Turbine outlet pressure greater than inlet pressure")
        
        self.outlet.state.update_ps(pOut, self.inlet.state.s)

        dHIdeal = self.outlet.state.h - self.inlet.state.h 
        dHActual = dHIdeal * eta
        
        self.outlet.state.update_ph(pOut, self.inlet.state.h + dHActual)

        self.WIn = dHActual * self.inlet.flow.mDot

        self.outlet.flow.mDot = self.inlet.flow.mDot

class Heater(TwoStationComponent):
    def __init__(self, mode):
        super().__init__()
        self.compType = "Heater"

        self.mode = mode

        self.QIn = 0

    def compute(self, TOut=None, QIn=None):
        if self.mode == "T":
            if TOut != None:
                set_TOut = max(self.inlet.state.T, TOut)

                self.outlet.state.update_Tp(set_TOut, self.inlet.state.p)
                
                self.QIn = (self.outlet.state.h - self.inlet.state.h) * self.inlet.flow.mDot
            else:
                raise ValueError("No TOut applied to heater in T mode")
        elif self.mode == "Q":
            if QIn != None:
                self.QIn = QIn
                self.outlet.state.update_ph(self.inlet.state.p, self.inlet.state.h + (QIn/self.inlet.flow.mDot))
            else:
                raise ValueError("No QIn applied to heater in Q mode")
            
        self.outlet.flow.mDot = self.inlet.flow.mDot

class Cooler(TwoStationComponent):
    def __init__(self, mode):
        super().__init__()
        self.compType = "Cooler"

        self.mode = mode

        self.QIn = 0

    def compute(self, TOut=None, QOut=None):
        if self.mode == "T":
            if TOut != None:
                set_TOut = min(self.inlet.state.T, TOut)
                self.outlet.state.update_Tp(set_TOut, self.inlet.state.p)

                self.QIn = (self.outlet.state.h - self.inlet.state.h) * self.inlet.flow.mDot
            else:
                raise ValueError("No TOut applied to cooler in T mode")
        elif self.mode == "Q":
            if QOut != None:
                self.QIn = -QOut

                self.outlet.state.update_ph(self.inlet.state.p, self.inlet.state.h - (QOut/self.inlet.flow.mDot))
            else:
                raise ValueError("No QOut applied to cooler in Q mode")
            
        self.outlet.flow.mDot = self.inlet.flow.mDot

class Recuperator(CycleComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Recuperator"
        
        self.inlet1 = Station()
        self.outlet1 = Station()
        self.inlet2 = Station()
        self.outlet2 = Station()

        self.stations = [self.inlet1, self.outlet1, self.inlet2, self.outlet2]

        self.QDot = 0
    
    def compute(self, eta=1):
        self.outlet1.state.update_Tp(self.inlet2.state.T, self.inlet1.state.p)
        self.outlet2.state.update_Tp(self.inlet1.state.T, self.inlet2.state.p)
        Q1Max = self.inlet1.flow.mDot * (self.inlet1.state.h - self.outlet1.state.h)
        Q2Max = self.inlet2.flow.mDot * (self.inlet2.state.h - self.outlet2.state.h)

        if (abs(Q1Max) > abs(Q2Max)):
            self.QDot = Q2Max * eta
        else:
            self.QDot = -Q1Max * eta
    
    def computeStream1(self):
        m1Dot = self.inlet1.flow.mDot
        self.outlet1.flow.mDot = m1Dot
        self.outlet1.state.update_ph(self.inlet1.state.p, self.inlet1.state.h + self.QDot / m1Dot)
    
    def computeStream2(self):
        m2Dot = self.inlet2.flow.mDot
        self.outlet2.flow.mDot = m2Dot
        self.outlet2.state.update_ph(self.inlet2.state.p, self.inlet2.state.h - self.QDot / m2Dot)
    
    def __repr__(self):
        return(self.compType+ str(self.inlet1)+ str(self.inlet2))

class IdealThrottle(TwoStationComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Throttle"
    
    def compute(self, pOut):
        if pOut > self.inlet.state.p:
            raise ValueError("Throttle outlet pressure greater than inlet pressure")

        self.outlet.state.update_ph(pOut, self.inlet.state.h)

        self.outlet.flow.mDot = self.inlet.flow.mDot

class PhaseSeperator(CycleComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Phase Seperator"

        self.inlet = Station()
        self.liquidOutlet = Station()
        self.gasOutlet = Station()
    
    def compute(self):
        if (0 < self.inlet.state.Q < 1):
            fq = self.inlet.state.Q
        else:
            fq = 1.0
		
        self.liquidOutlet.state.update_pQ(self.inlet.state.p, 0)
        self.liquidOutlet.flow.mDot = (1 - fq) * self.inlet.flow.mDot

        self.gasOutlet.state.update_pQ(self.inlet.state.p, 1)
        self.gasOutlet.flow.mDot = fq * self.inlet.flow.mDot
    
    def __repr__(self):
        return(self.compType + str(self.inlet))

class FlowJunction(CycleComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Flow Junction"
        
        self.inlet1 = Station()
        self.inlet2 = Station()
        self.outlet = Station()

        self.stations = [self.inlet1, self.inlet2, self.outlet]
    
    def compute(self):
        
        mDotIn = self.inlet1.flow.mDot + self.inlet2.flow.mDot

        hIn = (self.inlet1.state.h * self.inlet1.flow.mDot) + (self.inlet2.state.h * self.inlet2.flow.mDot)
        hIn /= mDotIn

        sIn = (self.inlet1.state.s * self.inlet1.flow.mDot) + (self.inlet2.state.s * self.inlet2.flow.mDot)
        sIn /= mDotIn

        TIn = (self.inlet1.state.T * self.inlet1.flow.mDot) + (self.inlet2.state.T * self.inlet2.flow.mDot)
        TIn /= mDotIn


        if abs(self.inlet1.state.p/self.inlet2.state.p - 1) > 1e-3:
            print("Junction with mismatched pressures", self.inlet1.flow.mDot, self.inlet1.state.p, self.inlet2.state.p)

        self.outlet.state.update_ph(self.inlet1.state.p, hIn)

        self.outlet.flow.mDot = mDotIn
    
    def __repr__(self):
        return(self.compType+ str(self.outlet))

class FlowSplitter(CycleComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Flow Splitter"
        
        self.inlet = Station()
        self.outlet1 = Station()
        self.outlet2 = Station()

        self.stations = [self.inlet, self.outlet1, self.outlet2]
    
    def compute(self, flowFraction2):
        self.outlet1.state.update_Trho(self.inlet.state.T, self.inlet.state.rho)
        self.outlet2.state.update_Trho(self.inlet.state.T, self.inlet.state.rho)

        self.outlet1.flow.mDot = self.inlet.flow.mDot * (1 - flowFraction2)
        self.outlet2.flow.mDot = self.inlet.flow.mDot * flowFraction2
    
    def __repr__(self):
        return(self.compType+ str(self.inlet))


class Evaporator(TwoStationComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Evaporator"

        self.QIn = 0
    
    def compute(self):
        TSat = PropsSI("T", "P", self.inlet.state.p*1e5, "Q", 0, self.inlet.state.fluid)
        if self.inlet.state.T > TSat:
            self.outlet.state.update_Tp(self.inlet.state.T, self.inlet.state.p)
            self.QIn = 0
        else:
            self.outlet.state.update_pQ(self.inlet.state.p, 1)
            self.QIn = (self.outlet.state.h - self.inlet.state.h) * self.inlet.flow.mDot

        self.outlet.flow.mDot = self.inlet.flow.mDot


class CycleSolver(object):
    def __init__(self, fluid, constant_fluid=True):
        self.fluid = fluid
        self.constant_fluid = constant_fluid

        self.components = {}
        self.fluidStates = []
        self.fluidFlows = []

        self.oldComponentStates = []
        self.residualHistory = []
    
    def connectStations(self, component1, station1, component2, station2):
            
        fs = FluidState(self.fluid)

        fl = FluidFlow()

        self.fluidStates.append(fs)
        setattr(getattr(self.components[component1], station1), "state", fs)
        setattr(getattr(self.components[component2], station2), "state", fs)

        self.fluidFlows.append(fl)
        setattr(getattr(self.components[component1], station1), "flow", fl)
        setattr(getattr(self.components[component2], station2), "flow", fl)
    
    def printStation(self, stationNumber, component, station):
        stat = getattr(self.components[component], station)

        print("{} \t {:.2f}K \t {:.2f}bar \t {:.2f}kJ/kg \t {:.2f}kJ/kg-K \t {:.2f} \t {:.2f}kg/min".format(
            stationNumber, stat.state.T, stat.state.p, stat.state.h/1000, stat.state.s/1000, stat.state.Q, stat.flow.mDot*60
        ))

    
    def initialiseSolver(self, T, p, mDot):
        
        for flow in self.fluidFlows:
            flow.mDot = mDot
        for state in self.fluidStates:
            state.update_Tp(T, p)
    
    def computeCycle(self):
        pass
