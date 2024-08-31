from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI
from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL
from thermo.interaction_parameters import IPDB
import matplotlib.pyplot as plt
import numpy as np
import copy


portNames = ["inlet", "outlet", 
             "inlet1", "inlet2", "outlet1", "outlet2", 
             "liquidOutlet", "gasOutlet", "airOutlet", "waterOutlet"]

class SolidState(object):
    def __init__(self, mass, specificHeatCapacity, TRef=0):
        self.mass = mass
        self.specificHeatCapacity = specificHeatCapacity
        self.TRef = 0

        self.heatCapacity = self.mass * self.specificHeatCapacity

        self.T = self.TRef
        self.h = 0

        self.accumulateH = 0
    
    def applyAccumulation(self, dt):
        self.accumulateH *= dt
        changeInT = self.accumulateH / self.heatCapacity
        self.T += changeInT
        self.h += self.accumulateH
    
    def resetAccumulation(self):
        self.accumulateH = 0

    def update_T(self, T):
        self.T = T
        self.h = self.heatCapacity * (self.T - self.TRef)
    
    def accumulate_T(self, TIn):
        self.accumulateH += (TIn - self.T) * self.heatCapacity
    
    def accumulate_H(self, HIn):
        self.accumulateH += HIn

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
        except:
            raise ValueError("Q {:.2f} p{:.2f}bar fluid ".format(self.Q, self.p)+self.fluid)

class MixtureState(FluidState):
    def __init__(self, fluid, components, ratios):
        self.fluid = fluid
        self.components = components
        self.ratios = ratios

        constants, properties = ChemicalConstantsPackage.from_IDs(['methane', 'carbon monoxide'])
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}

        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
        self.flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)

    def update_Tp(self, T_k, p_bar):
        self.T = T_k
        self.p = p_bar

        try:
            state = self.flasher.flash(T=self.T, P=self.p*1e5, zs=self.ratios)
            self.h = state.H_mass()
            self.s = state.S_mass()
            self.rho = state.rho_mass()
            self.Q = state.VF

        except:
            raise ValueError("T {:.1f}K p{:.2f}bar fluid ".format(self.T, self.p)+self.fluid)

class HumidAirState(FluidState):
    def __init__(self, fluid, components, ratios):
        super().__init__(fluid, components, ratios)
    
    def update_TpW(self, T_k, p_bar, W):
        self.T = T_k
        self.p = p_bar
        self.W = W

        # Check if temperature is below dew point
        try:
            self.T_dew = HAPropsSI("T_dp", "T", self.T, "P", self.p*1e5, "W", W)
        except:
            raise ValueError("Condensation occurs T {:.1f}K p{:.2f}bar {:.5f}kg/kg".format(self.T, self.p, W))

        try:
            self.h = HAPropsSI("Hha", "T", self.T, "P", self.p*1e5, "W", self.W)
            self.s = HAPropsSI("Sha", "T", self.T, "P", self.p*1e5, "W", self.W)
            self.RH = HAPropsSI("RH", "T", self.T, "P", self.p*1e5, "W", self.W)
        except:
            raise ValueError("T {:.1f}K p{:.2f}bar {:.2f} HumidAir".format(self.T, self.p, W))
    
    def update_TpRH(self, T_k, p_bar, RH):
        if not (0 <= RH <= 1):
            raise ValueError("Humidity outside 0-1 entered")
        self.T = T_k
        self.p = p_bar
        self.RH = RH

        # Check if temperature is below dew point
        try:
            self.T_dew = HAPropsSI("T_dp", "T", self.T, "P", self.p*1e5, "RH", RH)
        except:
            raise ValueError("Condensation occurs T {:.1f}K p{:.2f}bar {:.2f}%".format(self.T, self.p, RH))

        try:
            self.h = HAPropsSI("Hha", "T", self.T, "P", self.p*1e5, "RH", self.RH)
            self.s = HAPropsSI("Sha", "T", self.T, "P", self.p*1e5, "RH", self.RH)
            self.W = HAPropsSI("W", "T", self.T, "P", self.p*1e5, "RH", self.RH)
        except:
            raise ValueError("T {:.1f}K p{:.2f}bar {:.2f}% HumidAir".format(self.T, self.p, RH))
    
    def update_psW(self, p_bar, s, W):
        self.p = p_bar
        self.s = s
        self.W = W

        # Check if temperature is below dew point
        try:
            self.T_dew = HAPropsSI("T_dp", "Sha", self.s, "P", self.p*1e5, "W", W)
        except:
            raise ValueError("Condensation occurs p{:.1f}K s{:.0f}J/kgK {:.5f}kg/kg".format(self.p, self.s, W))

        try:
            self.h = HAPropsSI("Hha", "Sha", self.s, "P", self.p*1e5, "W", self.W)
            self.T = HAPropsSI("T", "Sha", self.s, "P", self.p*1e5, "W", self.W)
            self.RH = HAPropsSI("RH", "Sha", self.s, "P", self.p*1e5, "W", self.W)
        except:
            raise ValueError("p{:.2f}bar s{:.0f}J/kgK {:.2f}kg/kg HumidAir".format(self.p, self.s, W))
        
    def update_phW(self, p_bar, h, W):
        self.p = p_bar
        self.h = h
        self.W = W

        # Check if temperature is below dew point
        try:
            self.T_dew = HAPropsSI("T_dp", "Hha", self.h, "P", self.p*1e5, "W", W)
        except:
            raise ValueError("Condensation occurs p{:.1f}K h{:.0f}J/kgK {:.5f}kg/kg".format(self.p, self.h, W))

        try:
            self.s = HAPropsSI("Sha", "Hha", self.h, "P", self.p*1e5, "W", self.W)
            self.T = HAPropsSI("T", "Hha", self.h, "P", self.p*1e5, "W", self.W)
            self.RH = HAPropsSI("RH", "Hha", self.h, "P", self.p*1e5, "W", self.W)
        except:
            raise ValueError("p{:.2f}bar h{:.0f}J/kgK {:.2f}kg/kg HumidAir".format(self.p, self.h, W))

class ConstantFluidState(FluidState):
    def __init__(self, fluid, components, ratios):
        super().__init__(fluid, components, ratios)
    
    def initialise_Tp(self, T_k, p_bar, force_liquid=False, force_gas = False):
        self.T = T_k
        self.p = p_bar

        if not force_liquid and not force_gas:
            self.update_Tp(self.T, self.p)
        else:
            T_text = "T"
            if force_liquid:
                T_text = "T|liquid"
            elif force_gas:
                T_text = "T|gas"
            
            self.h = PropsSI("Hmass", T_text, self.T, "P", self.p*1e5, self.fluid)
            self.s = PropsSI("Smass", T_text, self.T, "P", self.p*1e5, self.fluid)
            self.rho = PropsSI("Dmass", T_text, self.T, "P", self.p*1e5, self.fluid)
            self.Q = PropsSI("Q", T_text, self.T, "P", self.p*1e5, self.fluid)
    
    def update_Tp(self, T_k, p_bar):
        pass

    def update_ps(self, p_bar, s):
        pass

    def update_ph(self, p_bar, h):
        pass

    def update_pQ(self, p_bar, Q):
        pass

    def update_Trho(self, T_k, rho):
        pass

class FluidFlow(object):
    def __init__(self):
        self.mDot = 0

        self.updated = False

class Station(object):
    def __init__(self):
        self.state = None
        self.flow = None
    
    def __repr__(self):
        if type(self.state) != SolidState:
            txt = " Fluid:{}, mDot {:.2f} P {:.1f}bar T{:.0f}K".format(self.state.fluid, self.flow.mDot, self.state.p, self.state.T)
        else:
            txt = " Solid, T{:.0f}K".format(self.state.T)
        return(txt)

class SolidStation(object):
    def __init__(self):
        self.state = None        

class CycleComponent(object):
    def __init__(self):
        self.compType = "None"
    
    def compute(self):
        pass

class SolidComponent(object):
    def __init__(self):
        self.compType = "None"
        self.linkedStations = []
    
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

        self.W = W
        self.RH = RH
        # W takes priority if both are filled in

        self.outlet = Station()

        self.stations = [self.outlet]
    
    def compute(self):
        self.outlet.flow.mDot = self.mDot

        if type(self.outlet.state) == HumidAirState:
            if self.W != None:
                self.outlet.state.update_TpW(self.T, self.p, self.W)
            elif self.RH != None:
                self.outlet.state.update_TpRH(self.T, self.p, self.W)
            else:
                raise ValueError("Humid air used but no humidity entered")
        else:
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
        
        if type(self.inlet.state) == FluidState:
            if self.inlet.state.Q != -1:
                self.outlet.state.update_ps(self.inlet.state.p+0.5, self.inlet.state.s)
                self.outlet.flow.mDot = self.inlet.flow.mDot
                return
        
        if type(self.inlet.state) == HumidAirState:
            self.outlet.state.update_psW(pOut, self.inlet.state.s, self.inlet.state.W)
        else:
            self.outlet.state.update_ps(pOut, self.inlet.state.s)

        dHIdeal = self.outlet.state.h - self.inlet.state.h 
        dHActual = dHIdeal / eta

        if type(self.inlet.state) == HumidAirState:
            self.outlet.state.update_phW(pOut, self.inlet.state.h + dHActual, self.inlet.state.W)
        else:
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
        #pOut = self.inlet.state.p / pressureRatio
        if type(self.inlet.state) == HumidAirState:
            self.outlet.state.update_psW(pOut, self.inlet.state.s, self.inlet.state.W)
        else:
            self.outlet.state.update_ps(pOut, self.inlet.state.s)

        dHIdeal = self.outlet.state.h - self.inlet.state.h 
        dHActual = dHIdeal * eta
        
        if type(self.inlet.state) == HumidAirState:
            self.outlet.state.update_phW(pOut, self.inlet.state.h + dHActual, self.inlet.state.W)
        else:
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

                if type(self.inlet.state) == HumidAirState:
                    self.outlet.state.update_TpW(set_TOut, self.inlet.state.p, self.inlet.state.W)
                else:
                    self.outlet.state.update_Tp(set_TOut, self.inlet.state.p)
                
                self.QIn = (self.outlet.state.h - self.inlet.state.h) * self.inlet.flow.mDot
            else:
                raise ValueError("No TOut applied to heater in T mode")
        elif self.mode == "Q":
            if QIn != None:
                self.QIn = QIn
                if type(self.inlet.state) == HumidAirState:
                    self.outlet.state.update_phW(self.inlet.state.p, self.inlet.state.h + (QIn/self.inlet.flow.mDot), 
                                                 self.inlet.state.W)
                else:
                    self.outlet.state.update_ph(self.inlet.state.p, self.inlet.state.h + (QIn/self.inlet.flow.mDot))
            else:
                raise ValueError("No QIn applied to heater in Q mode")
            
        self.outlet.flow.mDot = self.inlet.flow.mDot

class Cooler(TwoStationComponent):
    def __init__(self, mode):
        super().__init__()
        self.compType = "Cooler"

        self.mode = mode

        self.QOut = 0

    def compute(self, TOut=None, QOut=None):
        if self.mode == "T":
            if TOut != None:
                set_TOut = min(self.inlet.state.T, TOut)
                if type(self.inlet.state) == HumidAirState:
                    self.outlet.state.update_TpW(set_TOut, self.inlet.state.p, self.inlet.state.W)
                else:
                    self.outlet.state.update_Tp(set_TOut, self.inlet.state.p)

                self.QOut = (self.outlet.state.h - self.inlet.state.h) * self.inlet.flow.mDot
            else:
                raise ValueError("No TOut applied to cooler in T mode")
        elif self.mode == "Q":
            if QOut != None:
                self.QOut = QOut
                if type(self.inlet.state) == HumidAirState:
                    self.outlet.state.update_phW(self.inlet.state.p, self.inlet.state.h - (QOut/self.inlet.flow.mDot), 
                                                 self.inlet.state.W)
                else:
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

        self.Qdot = 0
    
    def computeQ(self, eta=1):
        old_outlet1_T = self.outlet1.state.T
        old_outlet1_p = self.outlet1.state.p
        old_outlet2_T = self.outlet2.state.T
        old_outlet2_p = self.outlet2.state.p

        
        if type(self.inlet1.state) == HumidAirState:
            self.outlet1.state.update_TpW(self.inlet2.state.T, self.inlet1.state.p, self.inlet1.state.W)
        else:
            self.outlet1.state.update_Tp(self.inlet2.state.T, self.inlet1.state.p)

        if type(self.inlet2.state) == HumidAirState:
            self.outlet2.state.update_TpW(self.inlet1.state.T, self.inlet2.state.p)
        else:
            self.outlet2.state.update_Tp(self.inlet1.state.T, self.inlet2.state.p)

        QdotFlow1 = self.inlet1.flow.mDot * (self.inlet1.state.h - self.outlet1.state.h)
        QdotFlow2 = self.inlet2.flow.mDot * (self.inlet2.state.h - self.outlet2.state.h)

        if abs(QdotFlow1) > abs(QdotFlow2):
            self.Qdot = QdotFlow2 * eta
        else:
            self.Qdot = -QdotFlow1 * eta
        
        if type(self.inlet1.state) == HumidAirState:
            self.outlet1.state.update_TpW(old_outlet1_T, old_outlet1_p, self.inlet1.state.W)
        else:
            self.outlet1.state.update_Tp(old_outlet1_T, old_outlet1_p)

        if type(self.inlet2.state) == HumidAirState:
            self.outlet2.state.update_TpW(old_outlet2_T, old_outlet2_p)
        else:
            self.outlet2.state.update_Tp(old_outlet2_T, old_outlet2_p)
        
    
    def computeStream1(self):
        mDot = self.inlet1.flow.mDot
        if type(self.inlet1.state) == HumidAirState:
            self.outlet1.state.update_phW(self.inlet1.state.p, self.inlet1.state.h + (self.Qdot / mDot), 
                                            self.inlet1.state.W)
        else:
            self.outlet1.state.update_ph(self.inlet1.state.p, self.inlet1.state.h + (self.Qdot / mDot))

        self.outlet1.flow.mDot = self.inlet1.flow.mDot
    
    def computeStream2(self):
        if self.inlet2.flow.mDot > 0:
            if type(self.inlet2.state) == HumidAirState:
                self.outlet2.state.update_phW(self.inlet2.state.p, self.inlet2.state.h - (self.Qdot / self.inlet2.flow.mDot), 
                                                self.inlet2.state.W)
            else:
                self.outlet2.state.update_ph(self.inlet2.state.p, self.inlet2.state.h + (self.Qdot / self.inlet2.flow.mDot))
        else:
            if type(self.inlet2.state) == HumidAirState:
                self.outlet2.state.update_TpW(self.inlet2.state.T, self.inlet2.state.p, self.inlet2.state.W)
            else:
                self.outlet2.state.update_Tp(self.inlet2.state.T, self.inlet2.state.p)
        
        self.outlet2.flow.mDot = self.inlet2.flow.mDot
    
    def __repr__(self):
        return(self.compType+ str(self.inlet1)+ str(self.inlet2))

class Evaporator(TwoStationComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Evaporator"

        self.QIn = 0
    
    def compute(self):
        if type(self.inlet.state) == HumidAirState:
            raise ValueError("Evaporator cannot work on humid air")
        TSat = PropsSI("T", "P", self.inlet.state.p*1e5, "Q", 0, self.inlet.state.fluid)
        if self.inlet.state.T > TSat:
            self.outlet.state.update_Tp(self.inlet.state.T, self.inlet.state.p)
            self.QIn = 0
        else:
            self.outlet.state.update_pQ(self.inlet.state.p, 1)
            self.QIn = (self.outlet.state.h - self.inlet.state.h) * self.inlet.flow.mDot

        self.outlet.flow.mDot = self.inlet.flow.mDot

class Condenser(TwoStationComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Condenser"

        self.QIn = 0
    
    def compute(self):
        if type(self.inlet.state) == HumidAirState:
            raise ValueError("Condenser cannot work on humid air")
        TSat = PropsSI("T", "P", self.inlet.state.p*1e5, "Q", 1, self.inlet.state.fluid)
        if self.inlet.state.T < TSat:
            self.outlet.state.update_Tp(self.inlet.state.T, self.inlet.state.p)
            self.QIn = 0
        else:
            self.outlet.state.update_pQ(self.inlet.state.p, 0)
            self.QIn = (self.outlet.state.h - self.inlet.state.h) * self.inlet.flow.mDot

        self.outlet.flow.mDot = self.inlet.flow.mDot

class FlowSplitter(CycleComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Flow Splitter"
        
        self.inlet = Station()
        self.outlet1 = Station()
        self.outlet2 = Station()

        self.stations = [self.inlet, self.outlet1, self.outlet2]
    
    def compute(self, flowFraction2):
        if type(self.inlet.state) == HumidAirState:
            self.outlet1.state.update_TpW(self.inlet.state.T, self.inlet.state.p, self.inlet.state.W)
            self.outlet2.state.update_TpW(self.inlet.state.T, self.inlet.state.p, self.inlet.state.W)
        else:
            self.outlet1.state.update_Tp(self.inlet.state.T, self.inlet.state.p)
            self.outlet2.state.update_Tp(self.inlet.state.T, self.inlet.state.p)


        self.outlet1.flow.mDot = self.inlet.flow.mDot * (1 - flowFraction2)
        self.outlet2.flow.mDot = self.inlet.flow.mDot * flowFraction2
    
    def __repr__(self):
        return(self.compType+ str(self.inlet))



class FlowJunction(CycleComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Flow Junction"
        
        self.inlet1 = Station()
        self.inlet2 = Station()
        self.outlet = Station()

        self.stations = [self.inlet1, self.inlet2, self.outlet]
    
    def compute(self):
        if type(self.inlet1.state) != type(self.inlet2.state):
            raise ValueError("Trying to mix two unlike fluids")
        
        mDotIn = self.inlet1.flow.mDot + self.inlet2.flow.mDot

        hIn = (self.inlet1.state.h * self.inlet1.flow.mDot) + (self.inlet2.state.h * self.inlet2.flow.mDot)
        hIn /= mDotIn

        sIn = (self.inlet1.state.s * self.inlet1.flow.mDot) + (self.inlet2.state.s * self.inlet2.flow.mDot)
        sIn /= mDotIn

        TIn = (self.inlet1.state.T * self.inlet1.flow.mDot) + (self.inlet2.state.T * self.inlet2.flow.mDot)
        TIn /= mDotIn

        if type(self.inlet1.state) == HumidAirState:
            WIn = (self.inlet1.state.W * self.inlet1.flow.mDot) + (self.inlet2.state.W * self.inlet2.flow.mDot)
            WIn /= mDotIn

        if abs(self.inlet1.state.p/self.inlet2.state.p - 1) > 1e-3:
            print("Junction with mismatched pressures", self.inlet1.flow.mDot, self.inlet1.state.p, self.inlet2.state.p)

        if type(self.inlet1.state) == HumidAirState:
            self.outlet.state.update_phW(self.inlet1.state.p, hIn, WIn)
        else:
            #self.outlet.state.update_hs(hIn, sIn)
            self.outlet.state.update_ps(self.inlet1.state.p, sIn)
            #self.outlet.state.update_Tp(TIn, self.inlet1.state.p)
        #print(self.outlet.flow.mDot, mDotIn, self.inlet1.flow.mDot, self.inlet2.flow.mDot)
        self.outlet.flow.mDot = self.inlet1.flow.mDot + self.inlet2.flow.mDot
    
    def __repr__(self):
        return(self.compType+ str(self.outlet))

class IdealThrottle(TwoStationComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Throttle"
    
    def compute(self, POut):
        if POut > self.inlet.state.p:
            raise ValueError("Throttle outlet pressure greater than inlet pressure")
        if type(self.inlet.state) == HumidAirState:
            self.outlet.state.update_ph(POut, self.inlet.state.h, self.inlet.state.W)
        else:
            self.outlet.state.update_ph(POut, self.inlet.state.h)
        self.outlet.flow.mDot = self.inlet.flow.mDot

class PhaseSeperator(CycleComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Phase Seperator"

        self.inlet = Station()
        self.liquidOutlet = Station()
        self.gasOutlet = Station()

        self.stations = [self.inlet, self.liquidOutlet, self.gasOutlet]

        self.prevQ = 1.0
    
    def compute(self, relaxLiquidRate=0.4):
        if type(self.inlet.state) == HumidAirState:
            raise ValueError("Phase seperator cannot work on humid air")

        nonCondense = False
        if 0 < self.inlet.state.Q < 1:
            iQ = self.inlet.state.Q
            # Adjust iQ to account for the previous cycle being less/more liquid - basically add a drag to this changing too rapidly
            
            self.prevQ = iQ
        else:
            iQ = 1.0

            """
            self.prevQ = -1.0
            #print("phase sep out of range")
            if PropsSI("Tcrit", self.inlet.state.fluid) > self.inlet.state.T:
                # Fluid is below critical temperature
                if PropsSI("P", "T", self.inlet.state.T, "Q", 0, self.inlet.state.fluid) > self.inlet.state.p*1e5:
                    iQ = 1
                else:
                    iQ = 0
            else:
                # Supercritical phase - pass to gas outlet
                iQ = 1.0
                self.liquidOutlet.flow.mDot = 0
        
        if self.prevQ > 0 and iQ > 0:
            iQ = self.prevQ + (iQ - self.prevQ)*relaxLiquidRate
        elif self.prevQ < 0 and iQ > 0:
            iQ = 1 + (iQ - 1)*relaxLiquidRate
        """
        #if self.prevQ < 0:
        #    iQ = self.prevQ + (iQ - 1)*relaxLiquidRate
        #else:
        #    iQ = self.prevQ + (iQ - self.prevQ)*relaxLiquidRate
        
        if False:
            self.gasOutlet.state.update_Tp(self.inlet.state.T, self.inlet.state.p)
            self.liquidOutlet.state.update_Tp(self.inlet.state.T, self.inlet.state.p)
        else:
            self.liquidOutlet.state.update_pQ(self.inlet.state.p, 0)
            self.gasOutlet.state.update_pQ(self.inlet.state.p, 1)

        self.liquidOutlet.flow.mDot = self.inlet.flow.mDot * (1 - iQ)
        #self.gasOutlet.state.update_pQ(self.inlet.state.p, 1)
        self.gasOutlet.flow.mDot = self.inlet.flow.mDot * iQ
    
    def __repr__(self):
        return(self.compType + str(self.inlet))

class PressureDrop(TwoStationComponent):
    def __init__(self):
        super().__init__()
        self.compType = "PressureDrop"
    
    def compute(self, POut):
        if POut > self.inlet.state.p:
            raise ValueError("Outlet pressure greater than inlet pressure")
        if type(self.inlet.state) == HumidAirState:
            self.outlet.state.update_ph(POut, self.inlet.state.h, self.inlet.state.W)
        else:
            self.outlet.state.update_ph(POut, self.inlet.state.h)
        self.outlet.flow.mDot = self.inlet.flow.mDot
    
class Dehumidifier(CycleComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Dehumidifier"

        self.inlet = Station()
        self.waterOutlet = Station()
        self.airOutlet = Station()

        self.stations = [self.inlet, self.waterOutlet, self.airOutlet]
    
    def compute(self, outletT=None, outletRH=None, outletW=None):
        if type(self.inlet.state) != HumidAirState:
            raise ValueError("Dehumidifier cannot work on non-humid air")

        if outletT != None:
            saturatedW = HAPropsSI("W", "T", outletT, "P", self.inlet.state.p, "RH", 1.0)
            if saturatedW < self.inlet.state.W:
                # Saturated air would have less water than the water has now - no water comes out of the air
                self.airOutlet.state.update_TpW(outletT, self.inlet.state.p, self.inlet.state.W)
                self.airOutlet.flow.mDot = self.inlet.flow.mDot
                self.waterOutlet.flow.mDot = 0
            else:
                self.airOutlet.state.update_TpRH(outletT, self.inlet.state.p, 1.0)
                mDotAir = self.inlet.flow.mDot * (1 / (1 + self.inlet.state.W ))
                mDotWaterIn = self.inlet.flow.mDot * (self.inlet.state.W / (1 + self.inlet.state.W ))
                mDotWaterOutInAir = mDotAir * self.airOutlet.state.W
                mDotWaterLoss = mDotWaterIn - mDotWaterOutInAir

                self.waterOutlet.state.update_Tp()

        elif outletRH != None:
            pass

        elif outletW != None:
            pass
        else:
            raise ValueError("Dehumidifier has no inputs")



        if 0 < self.inlet.state.Q < 1:
            self.liquidOutlet.state.update_pQ(self.inlet.state.p, 0)
            self.liquidOutlet.flow.mDot = self.inlet.flow.mDot * (1 - self.inlet.state.Q)

            self.gasOutlet.state.update_pQ(self.inlet.state.p, 1)
            self.gasOutlet.flow.mDot = self.inlet.flow.mDot * self.inlet.state.Q
        else:
            if PropsSI("Tcrit", self.inlet.state.fluid) > self.inlet.state.T:
                # Fluid is below critical temperature
                if PropsSI("P", "T", self.inlet.state.T, "Q", 0, self.inlet.state.fluid) > self.inlet.state.p:
                    # Fluid is above saturation pressure - all vapour
                    self.liquidOutlet.flow.mDot = 0

                    self.gasOutlet.state.update_Tp(self.inlet.state.T, self.inlet.state.p)
                    self.gasOutlet.flow.mDot = self.inlet.flow.mDot
                else:
                    self.liquidOutlet.state.update_Tp(self.inlet.state.T, self.inlet.state.p)
                    self.liquidOutlet.flow.mDot = self.inlet.flow.mDot

                    self.gasOutlet.flow.mDot = 0
            else:
                # Supercritical phase - pass to gas outlet
                self.liquidOutlet.flow.mDot = 0

                self.gasOutlet.state.update_Tp(self.inlet.state.T, self.inlet.state.p)
                self.gasOutlet.flow.mDot = self.inlet.flow.mDot

class SolidTemperatureBoundary(SolidComponent):
    def __init__(self):
        super().__init__()
        self.compType = "Solid temperature boundary"
    
    def compute(self, TIn):
        for station in self.linkedStations:
            station.accumulate_T(TIn)

class SolidThermalConductivity(SolidComponent):
    def __init__(self, thermalConductance=None):
        self.linkedStations = [SolidStation(), SolidStation]

        self.thermalConductance = thermalConductance
    
    def compute(self, thermalConductance=None):
        if self.thermalConductance != None:
            tC = self.thermalConductance
        elif thermalConductance != None:
            tC = thermalConductance
        else:
            raise ValueError("No thermal conductance set")
        
        dT = self.self.linkedStations[1].state.T - self.self.linkedStations[0].state.T
        Q = tC * dT
        self.inlet.state.accumulate_H(-Q)
        self.outlet.state.accumulate_H(Q)

class CycleSolver(object):
    def __init__(self, fluid, constant_fluid=True):
        self.fluid = fluid
        self.constant_fluid = constant_fluid

        self.components = {}
        self.fluidStates = []
        self.fluidFlows = []

        self.oldComponentStates = []
        self.residualHistory = []
    
    def connectStations(self, station1, station2, stateType = FluidState, fluidOverride = None):
        if fluidOverride == None:
            fs = stateType(self.fluid)
        else:
            fs = stateType(fluidOverride)

        fl = FluidFlow()

        self.fluidStates.append(fs)
        station1.state = fs
        station2.state = fs

        self.fluidFlows.append(fl)
        station1.flow = fl
        station2.flow = fl
    
    def initialiseSolver(self, T, p, mDot, W=0):
        
        for flow in self.fluidFlows:
            flow.mDot = mDot
        for state in self.fluidStates:
            if type(state) == HumidAirState:
                state.update_TpW(T, p, W)
            else:
                state.update_Tp(T, p)
    
    def computeCycle(self):
        pass

    def solveCycleRelaxationFactor(self, relaxationFactor=1, maxIterations=100, residualTarget=1e-6, limitAcceleration=None, limitVelocity=0.1):
        
        residual = 1e10
        previousResidual = 1e10

        for iteration in range(maxIterations):
            self.oldComponentStates.append(copy.deepcopy(self.components))

            self.computeCycle()

            previousResidual = copy.copy(residual)
            residual = 0

            #print(iteration, "Q2 {:.3f} Q {:.2f}, mLiquid {:.3f} mGas {:.3f}".format(
            #    self.components["th2"].outlet.state.Q, self.components["ps1"].inlet.state.Q, 
            #    self.components["ps1"].liquidOutlet.flow.mDot, self.components["ps1"].gasOutlet.flow.mDot))

            for c in self.components.keys():
                stationStates_n = []
                stationStates_nm1 = []
                stationStates_nm2 = []

                for s_i, s in enumerate(self.components[c].stations):
                    if self.components[c].stations[s_i].state.updated == True:
                        pass
                    else:
                        stationStates_n.append(self.components[c].stations[s_i])
                        stationStates_nm1.append(self.oldComponentStates[-1][c].stations[s_i])
                        if iteration >= 1:
                            stationStates_nm2.append(self.oldComponentStates[-2][c].stations[s_i])  
                        self.components[c].stations[s_i].state.updated = True
                
                #print(stationOldStates)
                #print(stationUpdatedStates) 

                if len(stationStates_n)>0:
                    for s_i, s in enumerate(stationStates_n):
                        
                        H_nm1 = stationStates_nm1[s_i].state.h 
                        H_n = copy.copy(stationStates_n[s_i].state.h)
                        dH_n = (H_n - H_nm1)

                        P_nm1 = stationStates_nm1[s_i].state.p
                        P_n = copy.copy(stationStates_n[s_i].state.p)

                        newT = copy.copy(stationStates_n[s_i].state.T)

                        mDot_nm1 = stationStates_nm1[s_i].flow.mDot
                        mDot_n = stationStates_n[s_i].flow.mDot
                        dMDot_n = (mDot_n - mDot_nm1)


                        if limitAcceleration!=None and iteration>=1:
                            dH_nm1 = (stationStates_nm1[s_i].state.h - stationStates_nm2[s_i].state.h)
                            if abs(dH_n/H_nm1) > (abs(dH_nm1/stationStates_nm2[s_i].state.h)*limitAcceleration):
                                dH_n = np.sign(dH_n) * dH_nm1 * limitAcceleration
                        
                        if limitVelocity!=None:
                            if abs(dH_n/H_nm1) > limitVelocity:
                                dH_n = np.sign(dH_n) * H_nm1 * limitVelocity
                            if mDot_nm1 != 0:
                                if abs(dMDot_n/mDot_nm1) > limitVelocity:
                                    dMDot_n = np.sign(dMDot_n) * mDot_nm1 * limitVelocity
                        
                        finalH = H_nm1 + (dH_n * relaxationFactor)
                        finalP = P_nm1 + ((P_n - P_nm1) * relaxationFactor)
                        finalMDot = mDot_nm1 + (dMDot_n * relaxationFactor)
                        
                        s.state.update_ph(finalP, finalH)
                        s.flow.mDot = finalMDot

                        residual += s.state.h

                        #print(stationOldStates[s_i].state.T, newT, s.state.T)


            for c in self.components.keys():
                for s in self.components[c].stations:
                    s.flow.updated = False
                    s.state.updated = False

            residualChange = abs(residual/previousResidual-1)
            self.residualHistory.append([iteration, residual, residualChange])
            #print(iteration, residual, previousResidual, abs(residual/previousResidual-1))

            #print(iteration, "Q {:.2f}, mLiquid {:.3f} mGas {:.3f}".format(
            #    self.components["ps1"].inlet.state.Q, self.components["ps1"].liquidOutlet.flow.mDot, self.components["ps1"].gasOutlet.flow.mDot))

        
                
                


    
    def cycleDiagramSimple_hs(self, componentOrder):
        hs_pairs = []

        for componentID in componentOrder:
            if componentID[1] == 0:
                hs_pairs.append([
                    self.components[componentID[0]].outlet.state.s/1000,
                    self.components[componentID[0]].outlet.state.h/1000
                ])
            elif componentID[1] == 1:
                hs_pairs.append([
                    self.components[componentID[0]].outlet1.state.s/1000,
                    self.components[componentID[0]].outlet1.state.h/1000
                ])
            elif componentID[1] == 2:
                hs_pairs.append([
                    self.components[componentID[0]].outlet2.state.s/1000,
                    self.components[componentID[0]].outlet2.state.h/1000
                ])
        

        fig, ax = plt.subplots()
        for index in range(len(hs_pairs)-1):
            pr_i = hs_pairs[index]
            pr_i1 = hs_pairs[index+1]
            halfway = [0.5*(pr_i[0]+pr_i1[0]), 0.5*(pr_i[1]+pr_i1[1])]
            ax.plot([pr_i[0], pr_i1[0]], [pr_i[1], pr_i1[1]], c="k")
            #ax.arrow(ph_i[0], ph_i[1], 0.5*(ph_i1[0]-ph_i[0]),  0.5*(ph_i1[1] - ph_i[1]), color="k", head_width=0.1, shape='full')
            ax.annotate("", xy=(halfway), xytext=pr_i, arrowprops=dict(arrowstyle="->"), size=20)
        ax.scatter([i[0] for i in hs_pairs], [i[1] for i in hs_pairs], c="k")
        
        ax.set_xlabel("Entropy (kJ/kgK)")
        ax.set_ylabel("Enthalpy (kJ/kg)")
        ax.set_title("Enthalpy-entropy diagram")
        
        plt.show()
    
    def cycleDiagramSimple_ph(self, componentOrder):
        ph_pairs = []

        for componentID in componentOrder:
            if componentID[1] == 0:
                ph_pairs.append([
                    self.components[componentID[0]].outlet.state.h/1000,
                    self.components[componentID[0]].outlet.state.p
                ])
            elif componentID[1] == 1:
                ph_pairs.append([
                    self.components[componentID[0]].outlet1.state.h/1000,
                    self.components[componentID[0]].outlet1.state.p
                ])
            elif componentID[1] == 2:
                ph_pairs.append([
                    self.components[componentID[0]].outlet2.state.h/1000,
                    self.components[componentID[0]].outlet2.state.p
                ])
        

        fig, ax = plt.subplots()
        for index in range(len(ph_pairs)-1):
            pr_i = ph_pairs[index]
            pr_i1 = ph_pairs[index+1]
            halfway = [0.5*(pr_i[0]+pr_i1[0]), np.exp(0.5*(np.log(pr_i[1])+np.log(pr_i1[1])))]
            ax.plot([pr_i[0], pr_i1[0]], [pr_i[1], pr_i1[1]], c="k")
            #ax.arrow(ph_i[0], ph_i[1], 0.5*(ph_i1[0]-ph_i[0]),  0.5*(ph_i1[1] - ph_i[1]), color="k", head_width=0.1, shape='full')
            ax.annotate("", xy=(halfway), xytext=pr_i, arrowprops=dict(arrowstyle="->"), size=20)
        ax.scatter([i[0] for i in ph_pairs], [i[1] for i in ph_pairs], c="k")
        
        ax.set_xlabel("Enthalpy (kJ/kg)")
        ax.set_ylabel("Pressure (bar)")
        ax.set_yscale("log")
        ax.set_title("Pressure-enthalpy diagram")
        
        plt.show()

class RK4SolverSolid(object):
    def __init__(self):
        self.components = {}
        self.stations = []

    
    def connectStations(self, station1, station2, mass, specificHeatCapacity, TRef):
        if station1.state != None:
            state = station1.state

            station2.state = state
        if station2.state != None:
            state = station2.state

            station1.state = state
        else:
            state = SolidState(mass, specificHeatCapacity, TRef)

            station1.state = state
            station2.state = state
    
    def initialiseSolver(self):
        for component_key in self.components.keys():
            for station in self.components[component_key].linkedStations:
                station.resetAccumulation()
                if station not in self.stations:
                    self.stations.append(station)
    
    def calculateDerivatives(self):
        for component_key in self.components.keys():
            self.components[component_key].calculate()
    
    def craneSolverMetman(self, dT, FLAG_INITL):
        pass

        
    


if __name__ == "__main":


    rk = CycleSolver("Water")

    rk.components["source"] = FluidSource_TP(0.001, 1, 298)
    rk.components["recup"] = Recuperator()
    rk.components["compr"] = Compressor()
    rk.components["combust"] = Heater("T")
    rk.components["turb"] = Turbine()
    rk.components["sink"] = FluidSink()

    def computeCycle(self):
        #rk.components["source"].compute()
        self.components["compr"].compute(10, 0.8)
        self.components["recup"].computeQ(eta=0.8)
        self.components["recup"].computeStream1()
        self.components["combust"].compute(TOut=1000)
        self.components["turb"].compute(1, 0.8)
        self.components["recup"].computeStream2()
        
        
    rk.computeCycle = computeCycle.__get__(rk)

    rk.connectStations(rk.components["source"].outlet, rk.components["compr"].inlet)
    rk.connectStations(rk.components["compr"].outlet, rk.components["recup"].inlet1)
    rk.connectStations(rk.components["recup"].outlet1, rk.components["combust"].inlet)
    rk.connectStations(rk.components["combust"].outlet, rk.components["turb"].inlet)
    rk.connectStations(rk.components["turb"].outlet, rk.components["recup"].inlet2)
    rk.connectStations(rk.components["recup"].outlet2, rk.components["sink"].inlet)

    """rk.connectStations(rk.components["source"].outlet, rk.components["compr"].inlet)
    rk.connectStations(rk.components["compr"].outlet, rk.components["combust"].inlet)
    rk.connectStations(rk.components["combust"].outlet, rk.components["turb"].inlet)
    rk.connectStations(rk.components["turb"].outlet, rk.components["sink"].inlet)"""

    rk.initialiseSolver()

    for i in range(5):
        rk.computeCycle()
        for s in rk.fluidStates:
            print("P {:.2f}bar, T {:.1f}K, Q{:.3f} Fluid:".format(s.p, s.T, s.Q)+s.fluid)
        print("")

    work = (rk.components["turb"].WIn) - (rk.components["compr"].WIn)
    heat = (rk.components["combust"].QIn)
    print("Work output {:.1f}W".format(work))
    print("Heat input {:.1f}W".format(heat))
    print("Engine efficiency {:.1f}%".format(work*100/heat))

    #rk.cycleDiagramSimple_ph([("source", 0), ("compr", 0), ("recup", 1), ("combust", 0), ("turb", 0), ("recup", 2)])
    rk.cycleDiagramSimple_ph([("source", 0), ("compr", 0), ("combust", 0), ("turb", 0)])