import numpy as np
from CoolProp.CoolProp import PropsSI

fluidPropOptions = ["simple", "coolprop"]
simpleGasProperties = {
    "Air": {"cp": 1005, "R": 287},
    "Water": {"cp": 4100, "R": None, "rho": 1000}
}

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

        self.ID = ""
    
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
    
    def resetSolver(self):
        self.T = self.Tnm1_temp
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

        self.nVarHistory = 0
        self.nDerivHistory = 0

        self.ID = ""
    
    def computeDerivatives(self, time, dt):
        self.state.computedTdt(time)
        self.state.intermediatedTdtHistory.append(self.state.dTdt)
    
    def backupVariables(self):
        self.state.Tnm1_temp = self.state.T * 1.0
    
    def updateVariables(self, variableHistoryValues={}, intermediateDerivativeHistoryValues={}, 
                        derivativeHistoryValues={},
                        updatePlusOne=False, returnValue=False, prevValue=1):
        nextT = self.state.Tnm1_temp * prevValue

        for key in variableHistoryValues:
            nextT += self.state.THistory[key] * variableHistoryValues[key]
        for key in intermediateDerivativeHistoryValues:
            nextT += self.state.intermediatedTdtHistory[key] * intermediateDerivativeHistoryValues[key]
        for key in derivativeHistoryValues:
            nextT += self.state.dTHistory[key] * derivativeHistoryValues[key]
        
        if returnValue:
            return([nextT])
        else:
            if updatePlusOne:
                self.state.Tnp1 = nextT * 1.0
            else:
                self.state.T = nextT * 1.0
    
    def setNextVar(self, dt):
        self.state.setNextT(dt)

        self.nVarHistory += 1
        self.nDerivHistory += 1
    
    def resetNextVar(self):
        self.state.resetSolver()
    
    def calculateErrorValueInput(self, correctedValues):
        errors = []
        errors.append(abs(self.state.Tnp1 - correctedValues[0])/self.state.Tnm1_temp)
        return(errors)

    def calculateErrorPreviousValues(self):
        errors = []
        errors.append(abs(self.state.Tnp1 - self.state.Tnm1_temp)/self.state.Tnm1_temp)
        return(errors)

class FluidState(object):
    def __init__(self, TRef=298, fluid="Nitrogen", mDot=1, p_bar=1, dilutesConcentrations={}, volume=1e-3, fluidProp="simple"):
        self._T = TRef
        self.fluid = fluid
        self.mDot = mDot
        self.p = p_bar
        self.volume = volume

        self.fluidProp = fluidProp

        if self.fluidProp == "simple":
            if fluid not in simpleGasProperties.keys():
                raise ValueError("Fluid is in simple mode but properties are not available")

        self.dilutesConcentration = dilutesConcentrations
        self.dilutes = {i:self.dilutesConcentration[i] * self.mass for i in self.dilutesConcentration}
        
        self.dilutesInflows = {i:0 for i in self.dilutes}
        self.dDilutes = {i:0 for i in self.dilutes}
        self.dilutesHistory = {i:[] for i in self.dilutes}
        self.dilutesConcentrationHistory = {i:[] for i in self.dilutes}

        self.boundaryConditionT = False
        self.boundaryConditionh = False

        self.T = self._T
        self.mass = self.volume * self.rho

        self.storedE = 0
        self.dhdt = 0
        self.hHistory = []
        self.dhHistory = []
        self.THistory = []
        self.dTHistory = []

        self.intermediatehHistory = []
        self.intermediatedhdtHistory = []
        self.hnp1 = 0
        self.hnm1_temp = 0

        self.ID = ""
    
    @property
    def T(self):
        return self._T
    
    @T.setter
    def T(self, value):
        if not self.boundaryConditionT:
            self._T = value * 1.0
        else:
            self._T = self.boundaryConditionT
        
        if self.fluidProp == "simple":
            self._h = simpleGasProperties[self.fluid]["cp"]*self._T
            self.cp = simpleGasProperties[self.fluid]["cp"]
            if simpleGasProperties[self.fluid]["R"] == None:
                self.rho = simpleGasProperties[self.fluid]["rho"]
            else:
                self.rho = self.p*1e5 / (self._T * simpleGasProperties[self.fluid]["R"])
        elif self.fluidProp == "coolprop":
            self._h = PropsSI("Hmass", "T", self._T, "P", self.p*1e5, self.fluid)
            self.rho = PropsSI("Dmass", "T", self._T, "P", self.p*1e5, self.fluid)
            self.rho = PropsSI("CP0MASS", "T", self._T, "P", self.p*1e5, self.fluid)
        self.mass = self.volume * self.rho
    
    @property
    def h(self):
        return self._h
    
    @h.setter
    def h(self, value):
        if not self.boundaryConditionh:
            self._h = value * 1.0
        else:
            self._h = self.boundaryConditionh
        try:
            if self.fluidProp == "simple":
                self.T = self._h / simpleGasProperties[self.fluid]["cp"]
            elif self.fluidProp == "coolprop":
                self.T = PropsSI("T", "Hmass", self.h, "P", self.p*1e5, self.fluid)
            self.mass = self.volume * self.rho
        except:
            print(self.ID, ">", self.h, self.hnm1_temp, self.dhHistory, self.intermediatedhdtHistory, self.T, self.mDot)
            exit()

        
    def computedhdt(self, time, dt):
        self.dhdt = self.storedE / self.mass
        self.storedE = 0

        for dilute in self.dilutesInflows:
            self.dDilutes[dilute] = self.dilutesInflows[dilute] * 1.0
            self.dilutesInflows[dilute] = 0
        
            
    def setNexth(self, dt):
        self.hHistory.append(self.hnm1_temp)
        self.dhHistory.append((self.hnp1 - self.hnm1_temp)/dt)
        self.THistory.append(self.T)
        self.h = self.hnp1
        self.dTHistory.append((self.T - self.THistory[-1])/dt)

        for dilute in self.dDilutes:
            self.dilutesHistory[dilute].append(self.dilutes[dilute])
            self.dilutesConcentrationHistory[dilute].append(self.dilutesConcentration[dilute])

            self.dilutes[dilute] = max(0, self.dilutes[dilute] + self.dDilutes[dilute]*dt)

            #self.dDilutes[dilute] = 0
            self.dilutesInflows[dilute] = 0

            self.dilutesConcentration[dilute] = self.dilutes[dilute] / self.mass

        self.hnp1 = 0
        self.intermediatehHistory = []
        self.intermediatedhdtHistory = []

        self.dDilutes = {i:0 for i in self.dilutes}
        self.dilutesInflows = {i:0 for i in self.dilutes}
    
    def resetSolver(self):
        self.h = self.hnm1_temp
        self.hnp1 = 0
        self.intermediatehHistory = []
        self.intermediatedhdtHistory = []

        self.dDilutes = {i:0 for i in self.dilutes}
        self.dilutesInflows = {i:0 for i in self.dilutes}

    
    def stateInitialise(self, TRef, mDot, p):
        self.THistory = []
        self.intermediateTHistory = []
        self.intermediatedTdtHistory = []
        self.Tnm1_temp = 0
        self.Tnp1 = 0 
        self.T = TRef
        if mDot != None:
            self.mDot = mDot
        self.p = p
        self.storedE = 0

    def dilutesInitialise(self, dilutes):
        self.dilutes = dilutes
        self.dilutesConcentration = {i:dilutes[i]/self.mass for i in dilutes}
        self.dilutesInflows = {i:0 for i in dilutes}
        self.dDilutes = {i:0 for i in dilutes}
        self.dilutesHistory = {i:[] for i in dilutes}
        self.dilutesConcentrationHistory = {i:[] for i in dilutes}

class FluidStation(object):
    def __init__(self):
        self.state = None

        self.nVarHistory = 0
        self.nDerivHistory = 0

        self.countErrors = True

        self.ID = ""
    
    def computeDerivatives(self, time, dt):
        self.state.computedhdt(time, dt)
        self.state.intermediatedhdtHistory.append(self.state.dhdt)
    
    def backupVariables(self):
        self.state.hnm1_temp = self.state.h * 1.0
    
    def updateVariables(self, variableHistoryValues={}, intermediateDerivativeHistoryValues={}, 
                        derivativeHistoryValues={},
                        updatePlusOne=False, returnValue=False, prevValue=1):
        nexth = self.state.hnm1_temp * prevValue

        for key in variableHistoryValues:
            nexth += self.state.hHistory[key] * variableHistoryValues[key]
        for key in intermediateDerivativeHistoryValues:
            nexth += self.state.intermediatedhdtHistory[key] * intermediateDerivativeHistoryValues[key]
        for key in derivativeHistoryValues:
            nexth += self.state.dhHistory[key] * derivativeHistoryValues[key]
        
        if returnValue:
            return([nexth])
        else:
            if updatePlusOne:
                self.state.hnp1 = nexth * 1.0
            else:
                self.state.h = nexth * 1.0
    
    def setNextVar(self, dt):
        self.state.setNexth(dt)

        self.nVarHistory += 1
        self.nDerivHistory += 1
    
    def resetNextVar(self):
        self.state.resetSolver()
    
    def calculateErrorValueInput(self, correctedValues):
        errors = []
        errors.append(abs(self.state.hnp1 - correctedValues[0])/self.state.hnm1_temp)
        if self.countErrors:
            return(errors)
        else:
            return(0)

    def calculateErrorPreviousValues(self):
        errors = []
        errors.append(abs(self.state.hnp1 - self.state.hnm1_temp)/self.state.hnm1_temp)
        if self.countErrors:
            return(errors)
        else:
            return(0)

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

        self.conductionHeatTransfer = 0

        self.linkedStations = [self.inlet, self.outlet]
    
    def compute(self, time):
        dT = self.inlet.state.T - self.outlet.state.T
        dE = dT * self.C

        self.inlet.state.storedE -= dE
        self.outlet.state.storedE += dE
        self.conductionHeatTransfer = dE

class ConstantTemperatureBoundary(SolidComponent):
    def __init__(self, TSet):
        super().__init__()

        self.TSet = TSet

        self.Q = 0

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

        self.Q = dT * self.outlet.state.heatCapacity

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

class CoolantPseudoConduction(SolidComponent):
    def __init__(self, SHC=1, mDot=1):
        super().__init__()
        self.SHC = SHC
        self.mDot = mDot
        # Thermal conductance, W/K
        
        self.inlet = SolidStation()
        self.outlet = SolidStation()

        self.linkedStations = [self.inlet, self.outlet]
    
    def compute(self, time):
        if callable(self.mDot):
            mDot = self.mDot(time)
        else:
            mDot = self.mDot
        conductivity = mDot * self.SHC
        dT = self.inlet.state.T - self.outlet.state.T

        self.inlet.state.storedE -= dT * conductivity
        self.outlet.state.storedE += dT * conductivity

class RadiationTransfer(SolidComponent):
    def __init__(self, areaInlet=1, areaOutlet=1, 
                 emittanceInlet=1, emittanceOutlet=1, 
                 absorbanceInlet=None, absorbanceOutlet=None,
                 viewFactor=1):
        super().__init__()
        self.areaInlet = areaInlet
        self.areaOutlet = areaOutlet

        self.emittanceInlet = emittanceInlet
        self.emittanceOutlet = emittanceOutlet

        self.viewFactor = viewFactor

        if absorbanceInlet != None: 
            self.absorbanceInlet = absorbanceInlet
        else:
            self.absorbanceInlet = self.emittanceInlet

        if absorbanceOutlet != None:
            self.absorbanceOutlet = absorbanceOutlet
        else:
            self.absorbanceOutlet = self.emittanceOutlet
        
        self.inlet = SolidStation()
        self.outlet = SolidStation()

        self.radiationHeatTransfer = 0

        self.linkedStations = [self.inlet, self.outlet]
    
    def compute(self, time):
        energy = 5.670374419e-8 * self.areaInlet * self.viewFactor * (
            (self.emittanceInlet * self.absorbanceOutlet * np.power(self.inlet.state.T, 4)) -
            (self.emittanceOutlet * self.absorbanceInlet * np.power(self.outlet.state.T, 4))
        ) 

        self.inlet.state.storedE -= energy
        self.outlet.state.storedE += energy
        self.radiationHeatTransfer = energy

class FluidSource(SolidComponent):
    def __init__(self, mDot=0, fluid="Nitrogen", p=1, T=298, dilutesFlow={}):
        super().__init__()
        self.mDot = mDot
        self.T = T
        self.p = p
        self.dilutesFlow = dilutesFlow
        self.fluid = fluid
        
        self.fluidOutlet = FluidStation()

        self.linkedStations = [self.fluidOutlet]
    
    def compute(self, time):
        if callable(self.mDot):
            mDot = self.mDot(time)
        else:
            mDot = self.mDot
        if callable(self.T):
            T = self.T(time)
        else:
            T = self.T

        self.inputH = 1005*T # PropsSI("Hmass", "T", T, "P", self.p*1e5, self.fluid)
        
        self.fluidOutlet.state.boundaryConditionT = T
        self.fluidOutlet.state.mDot =   mDot
        self.fluidOutlet.state.p =      self.p

        self.fluidOutlet.state.storedE += self.inputH  * mDot

        for dilute in self.dilutesFlow:
            diluteFlow = self.dilutesFlow[dilute]
            self.fluidOutlet.state.dilutesInflows[dilute] += diluteFlow * 1.0
        
class FluidSink(SolidComponent):
    def __init__(self, mDot=1,):
        super().__init__()
        self.mDot = mDot

        self.fluidInlet = FluidStation()

        self.linkedStations = [self.fluidInlet]
    
    def compute(self, time):
        self.fluidInlet.state.storedE -= self.fluidInlet.state.mDot * self.fluidInlet.state.h

        for dilute in self.fluidInlet.state.dilutesInflows:
            self.fluidInlet.state.dilutesInflows[dilute] -= self.fluidInlet.state.dilutesConcentration[dilute] * self.fluidInlet.state.mDot
        # if callable(self.mDot):
        #     mDot = self.mDot(time)
        # else:
        #     mDot = self.mDot
        # if callable(self.T):
        #     T = self.T(time)
        # else:
        #     T = self.T
        
        # self.outlet.state.boundaryConditionT = T
        # self.fluidOutlet.state.mDot =   self.fluidInlet.state.mDot
        # self.fluidOutlet.state.p =      self.fluidInlet.state.p

        # for dilute in self.dilutesFlow:
        #     diluteFlow = self.dilutesFlow[dilute]
        #     if dilute in self.fluidOutlet.state.dilutesInflows:
        #         self.fluidOutlet.state.dilutesInflows[dilute] += diluteFlow
        #     else:
        #         self.fluidOutlet.state.dilutesInflows[dilute] =  diluteFlow

class FluidLink(SolidComponent):
    def __init__(self, mDot=0):
        super().__init__()
        self.mDot = 0
        self.flowQ = 0
        
        self.fluidInlet = FluidStation()
        self.fluidOutlet = FluidStation()

        self.linkedStations = [self.fluidInlet, self.fluidOutlet]
    
    def compute(self, time):
        if callable(self.fluidInlet.state.mDot):
            self.mDot = self.fluidInlet.state.mDot(time)
        else:
            self.mDot = self.fluidInlet.state.mDot

        self.fluidOutlet.state.mDot =   self.fluidInlet.state.mDot
        self.fluidOutlet.state.p =      self.fluidInlet.state.p

        self.fluidInlet.state.storedE  -= self.fluidInlet.state.mDot * self.fluidInlet.state.h
        self.fluidOutlet.state.storedE += self.fluidInlet.state.mDot * self.fluidInlet.state.h
        self.flowQ = self.fluidInlet.state.mDot * self.fluidInlet.state.h

        for dilute in self.fluidInlet.state.dilutes:
            diluteFlow = self.mDot * self.fluidInlet.state.dilutesConcentration[dilute]
            diluteFlow = min(diluteFlow, self.fluidInlet.state.dilutes[dilute])
            
            self.fluidOutlet.state.dilutesInflows[dilute] += diluteFlow
            
            self.fluidInlet.state.dilutesInflows[dilute] -=  diluteFlow

class FluidHeatTransfer(SolidComponent):
    def __init__(self, thermalConductance):
        super().__init__()
        self.mDot = 0
        self.solidHeatTransfer = 0
        self.thermalConductance = thermalConductance
        self.heatTransferArea = 0
        self.heatTransferCoefficient = 0
        self.useHeatTransferCoefficient = False
        
        self.fluidInlet = FluidStation()
        self.fluidOutlet = FluidStation()
        self.solidLink = SolidStation()

        self.linkedStations = [self.fluidInlet, self.fluidOutlet, self.solidLink]
    
    def compute(self, time):
        if callable(self.fluidInlet.state.mDot):
            self.mDot = self.fluidInlet.state.mDot(time)
        else:
            self.mDot = self.fluidInlet.state.mDot

        self.fluidOutlet.state.mDot =   self.fluidInlet.state.mDot
        self.fluidOutlet.state.p =      self.fluidInlet.state.p

        if self.useHeatTransferCoefficient:
            self.thermalConductance = self.area * self.heatTransferCoefficient

        # Theoretical maximum conduction heat
        conductionHeat = self.thermalConductance * (
            self.fluidInlet.state.T - self.solidLink.state.T)
        # Test how much heat could be carried out
        fluidOutletOldT = self.fluidOutlet.state.T * 1.0
        self.fluidOutlet.state.T = self.solidLink.state.T
        maxConvectionHeat = self.mDot * (self.fluidOutlet.state.h - self.fluidInlet.state.h)

        conductionHeat = min(abs(conductionHeat), abs(maxConvectionHeat)) * np.sign(conductionHeat)

        self.fluidOutlet.state.T = fluidOutletOldT

        self.solidLink.state.storedE += conductionHeat

        self.fluidInlet.state.storedE -= self.fluidInlet.state.mDot * self.fluidInlet.state.h
        self.fluidOutlet.state.storedE += self.fluidInlet.state.mDot * self.fluidInlet.state.h - (
            conductionHeat)
        
        self.solidHeatTransfer = conductionHeat


        for dilute in self.fluidInlet.state.dilutes:
            diluteFlow = self.mDot * self.fluidInlet.state.dilutesConcentration[dilute]
            diluteFlow = min(diluteFlow, self.fluidInlet.state.dilutes[dilute])

            self.fluidOutlet.state.dilutesInflows[dilute] += diluteFlow
            self.fluidInlet.state.dilutesInflows[dilute] -= diluteFlow

class DiluteSourceSink(SolidComponent):
    def __init__(self, diluteRate):
        super().__init__()
        self.diluteRate = diluteRate
        
        self.fluidOutlet = FluidStation()

        self.linkedStations = [self.fluidOutlet]
    
    def compute(self, time):
        for dilute in self.diluteRate:
            if callable(self.diluteRate[dilute]):
                rate = self.diluteRate[dilute](time)
            else:
                rate = self.diluteRate[dilute]
                
            self.fluidOutlet.state.dilutesInflows[dilute] += rate

class GeneralFluidLink(SolidComponent):
    def __init__(self, 
                 flowModel=True,
                 pressureDropModel={"active": False,}, 
                 convectionModel={"active": False, "useLMTD": False}, 
                 dilutesModel={"active": False}):
        super().__init__()
        self.mDot = 0
        self.convectionHeatTransfer = 0
        self.evaporationHeatTransfer = 0

        self.flowModel = flowModel
        self.pressureDropModel = pressureDropModel
        self.convectionModel = convectionModel
        self.dilutesModel = dilutesModel
        
        self.fluidInlet = FluidStation()
        self.fluidOutlet = FluidStation()
        self.solidLink = SolidStation()

        self.linkedStations = [self.fluidInlet, self.fluidOutlet, self.solidLink]
    
    def compute(self, time):
        if callable(self.fluidInlet.state.mDot):
            self.mDot = self.fluidInlet.state.mDot(time)
        else:
            self.mDot = self.fluidInlet.state.mDot
        

        if self.flowModel:
            self.fluidOutlet.state.mDot =   self.mDot
            self.fluidOutlet.state.storedE += self.fluidInlet.state.h * self.mDot
            self.fluidInlet.state.storedE -= self.fluidInlet.state.h * self.mDot


        if self.pressureDropModel["active"]:
            pass
        else:
            self.fluidOutlet.state.p =      self.fluidInlet.state.p

        if self.convectionModel["active"]:
            if self.convectionModel["mode"] == "HTC":
                thermalConductance = self.convectionModel["area"] * self.convectionModel["HTC"]
            elif self.convectionModel["mode"] == "simple":
                thermalConductance = self.convectionModel["conductance"]
            elif self.convectionModel["mode"] == "NTU":
                thermalConductance = (1 - np.exp(-self.convectionModel["NTU"])) * self.fluidInlet.state.mDot * self.fluidInlet.state.cp

            if self.convectionModel["useLMTD"]:
                # Use the log-mean temperature difference
                # Assume the solid link has a constant temperature 
                deltaTInlet = (self.fluidInlet.state.T - self.solidLink.state.T)
                deltaTOutlet = (self.fluidInlet.state.T - self.solidLink.state.T)
                try:
                    if abs(deltaTInlet-deltaTOutlet) < 1e-5:
                        raise ValueError
                    LMTD = (deltaTInlet - deltaTOutlet) / np.log(deltaTInlet/deltaTOutlet)
                    if LMTD > max(deltaTInlet, deltaTOutlet) or LMTD < min(deltaTInlet, deltaTOutlet):
                        raise ValueError
                    else:
                        deltaT = LMTD
                except:
                    deltaT = np.mean((deltaTInlet, deltaTOutlet))
            else:
                # Use the absolute temperature difference
                deltaT = (self.fluidInlet.state.T - self.solidLink.state.T)

            # Theoretical maximum conduction heat
            conductionHeat = thermalConductance * deltaT
            # Test how much heat could be carried out
            fluidOutletOldT = self.fluidOutlet.state.T * 1.0
            self.fluidOutlet.state.T = self.solidLink.state.T
            maxConvectionHeat = self.mDot * (self.fluidOutlet.state.h - self.fluidInlet.state.h)
            
            convectionHeat = min(abs(conductionHeat), abs(maxConvectionHeat)) * np.sign(conductionHeat)

            self.fluidOutlet.state.T = fluidOutletOldT

            self.solidLink.state.storedE += convectionHeat
            self.fluidOutlet.state.storedE -= convectionHeat

            self.convectionHeatTransfer = convectionHeat
        else:
            self.convectionHeatTransfer = 0
        
        self.evaporationHeatTransfer = 0
        for dilute in self.fluidInlet.state.dilutes:
            diluteInletFlow = 0
            if self.flowModel:
                diluteInletFlow = self.mDot * self.fluidInlet.state.dilutesConcentration[dilute]
                diluteInletFlow = min(diluteInletFlow, self.fluidInlet.state.dilutes[dilute])

            diluteOutletFlow = diluteInletFlow

            if self.dilutesModel["active"]:
                if self.dilutesModel[dilute]["usemDot"]:
                    if callable(self.dilutesModel[dilute]["mDot"]):
                        evapFlow = self.dilutesModel[dilute]["mDot"](time)
                    else:
                        evapFlow = self.dilutesModel[dilute]["mDot"]
                    evapHeat = evapFlow * self.dilutesModel[dilute]["enthalpy"]
                else:
                    # Use energy
                    if callable(self.dilutesModel[dilute]["Q"]):
                        evapHeat = self.dilutesModel[dilute]["Q"](time)
                    else:
                        evapHeat = self.dilutesModel[dilute]["Q"]
                    evapFlow = abs(evapHeat) / self.dilutesModel[dilute]["enthalpy"]

                diluteOutletFlow += evapFlow

                self.evaporationHeatTransfer += evapHeat


            self.fluidOutlet.state.dilutesInflows[dilute] += diluteOutletFlow
            self.fluidInlet.state.dilutesInflows[dilute] -= diluteInletFlow

        if self.solidLink.state != None:
            self.solidLink.state.storedE -= self.evaporationHeatTransfer

class FluidJunction(SolidComponent):
    def __init__(self, nInputs=1, nOutputs=1, outputSplitRatio=[1]):
        self.nInputs = nInputs
        self.nOutputs = nOutputs
        self.linkedStations = []
        self.fluidInlets = []
        self.fluidOutlets = []
        for nI in range(nInputs):
            setattr(self, "fluidInlet"+str(nI), FluidStation())
            self.linkedStations.append(getattr(self, "fluidInlet"+str(nI)))
            self.fluidInlets.append(getattr(self, "fluidInlet"+str(nI)))
        
        for nO in range(nOutputs):
            setattr(self, "fluidOutlet"+str(nO), FluidStation())
            self.linkedStations.append(getattr(self, "fluidOutlet"+str(nO)))
            self.fluidOutlets.append(getattr(self, "fluidOutlet"+str(nO)))
        
        self.outputSplitRatio = outputSplitRatio
        if len(self.outputSplitRatio) != nOutputs:
            raise ValueError("Wrong number of split ratios in junction")
    
    def compute(self, time):
        linkedFluids = [i.state.fluid for i in self.linkedStations]
        if len(set(linkedFluids)) > 1:
            raise ValueError("Multiple fluids in the same system")
        
        mDotIn = 0
        hIn = 0
        dilutesIn = {}

        for i in self.fluidInlets:
            inletmDot = 0
            if callable(i.state.mDot):
                inletmDot = i.state.mDot(time)
            else:
                inletmDot = i.state.mDot
            
            mDotIn += inletmDot
            hIn += inletmDot * i.state.h
            i.state.storedE -= inletmDot * i.state.h

            for dilute in i.state.dilutes:
                if dilute not in dilutesIn.keys():
                    dilutesIn[dilute] = inletmDot * i.state.dilutesConcentration[dilute]
                else:
                    dilutesIn[dilute] += inletmDot * i.state.dilutesConcentration[dilute]
                
                i.state.dilutesInflows[dilute] -=  inletmDot * i.state.dilutesConcentration[dilute]
        
        #print(mDotIn, hIn)
        #print(i.state.h, i.state.storedE)

        for o_index, o in enumerate(self.fluidOutlets):
            o.state.mDot = mDotIn * self.outputSplitRatio[o_index]
            o.state.storedE += hIn * self.outputSplitRatio[o_index]

            for dilute in dilutesIn:
                o.state.dilutesInflows[dilute] += dilutesIn[dilute] * self.outputSplitRatio[o_index]
    
class TimestepSolver(object):
    def __init__(self, precision=4, maxTimestep=1, fluidProp="simple"):
        self.components = {}
        self.stations = []
        self.stationsUnique = []
        self.stationsDuplicates = []
        self.states = []

        self.precision = precision
        self.maxTimestep = maxTimestep
        self.minRunsBetweenTimestepIncrease = 3
        self.errors = []
        self.time = []
        self.timesteps = []
        self.allowTimeAdjust = True

        self.mode = "RK4ERROR"

        if fluidProp not in fluidPropOptions:
            raise ValueError("fluidProp not a valid option")
        self.fluidProp = fluidProp

        self.controlParameters = {}

    def connectStationsSolid(self, component1, station1, component2, station2, 
                             mass=1, specificHeatCapacity=1, TRef=298):
        
        if getattr(getattr(self.components[component1], station1), "state") != None:
            fs = getattr(getattr(self.components[component1], station1), "state")
        elif getattr(getattr(self.components[component2], station2), "state") != None:
            fs = getattr(getattr(self.components[component2], station2), "state")
        else:
            fs = SolidState(mass=mass, specificHeatCapacity=specificHeatCapacity, TRef=TRef)

        setattr(getattr(self.components[component1], station1), "state", fs)
        setattr(getattr(self.components[component2], station2), "state", fs)

    def connectStationsFluid(self, component1, station1, component2, station2, 
                             fluid="Nitrogen", p_bar=1, mDot=1, TRef=298, dilutes={}, volume=1e-3):

        if getattr(getattr(self.components[component1], station1), "state") != None:
            fs = getattr(getattr(self.components[component1], station1), "state")
        elif getattr(getattr(self.components[component2], station2), "state") != None:
            fs = getattr(getattr(self.components[component2], station2), "state")
        else:
            fs = FluidState(fluid=fluid, p_bar=p_bar, mDot=mDot, TRef=TRef, dilutesConcentrations=dilutes, volume=volume)

        setattr(getattr(self.components[component1], station1), "state", fs)
        setattr(getattr(self.components[component2], station2), "state", fs)
    
    def initialiseSolver(self, Tref, mDot=1, p_bar=1, dilutes={}):
        for component_key in self.components.keys():
            for station in self.components[component_key].linkedStations:
                if station.state == None:
                    pass
                else:
                    if type(station.state) == SolidState:
                        station.state.stateInitialise(Tref)
                        if station.state.ID == "":
                            station.state.ID = component_key
                    else:
                        station.state.stateInitialise(Tref, mDot, p_bar)
                        station.state.dilutesInitialise(dilutes)
                        if station.state.ID == "":
                            station.state.ID = component_key

                    if station.state not in self.states:
                        self.stationsUnique.append(station)
                        if type(station) == SolidStation:
                            self.stations.append([station, station.state, None])
                        else:
                            self.stations.append([station, station.state, None])
                        self.states.append(station.state)

        self.errors = []
        self.time = []
    
    def calculateDerivatives(self, time, dt):
        for component_key in self.components.keys():
            self.components[component_key].compute(time)
        
        for station in self.stations:
            station[0].computeDerivatives(time, dt)
    
    def changeTimestep(self, dt, errorFailed):
        if self.allowTimeAdjust:
            if errorFailed:
                return(dt * 0.625)
            else:
                return(min(self.maxTimestep, dt * 2))
        else:
            return(dt)

    def solverStepSuccess(self, time, dt, error):
        for station in self.stations:
            station[0].setNextVar(dt)
        self.time.append(time)
        self.errors.append(error)
        self.timesteps.append(dt)
        
        self.systemController(time, dt)

        time += dt

        return(time)

    def runSolver(self, t0, t1, dtTrial):
        time = t0
        dt = dtTrial
        runsSinceTimestepIncrease = 0
        numIters = 0
        runsSinceLastError = 0

        milestones = np.linspace(t0, t1, 26)

        while time < t1:
            #print(time, dt)
            if self.mode == "RK4ERROR":
                success, error = self.RK4OneStepError(time, dt)
                if success:
                    time = self.solverStepSuccess(time, dt, error)

                    runsSinceTimestepIncrease += 1
                    if runsSinceTimestepIncrease >= self.minRunsBetweenTimestepIncrease:
                        dt = self.changeTimestep(dt, False)
                        runsSinceTimestepIncrease = 0
                else:
                    for station in self.stations:
                        station[0].resetNextVar()
                    dt = self.changeTimestep(dt, True)

            elif self.mode == "RK4GILL":

                success, error = self.RK4OneStepGill(time, dt)
                time = self.solverStepSuccess(time, dt, np.nan)

            elif self.mode == "CRANE":
                pass
                #success, error = self.PredictorCorrectorCraneError(time, dt)

            elif self.mode == "CRANERK4":
                # Crane algorithm, falling back to RK4 with the previous timestep
                if runsSinceLastError < 3:
                    # We don't have enough previous dT/dt values to use the Crane method - fall back to RK4
                    success, error = self.RK4OneStepGill(time, dt)
                    if success:
                        runsSinceLastError += 1
                        time = self.solverStepSuccess(time, dt, error)

                        runsSinceTimestepIncrease += 1
                else:
                    # We can happily use the Crane algorithm
                    success, error = self.PredictorCorrectorCraneError(time, dt)
                    if success:
                        runsSinceLastError += 1
                        time = self.solverStepSuccess(time, dt, error)

                        runsSinceTimestepIncrease += 1

                        if runsSinceTimestepIncrease >= self.minRunsBetweenTimestepIncrease:
                            dt = self.changeTimestep(dt, False)
                            runsSinceTimestepIncrease = 0
                    else:
                        for station in self.stations:
                            station[0].resetNextVar()
                        dt = self.changeTimestep(dt, True)

            else:
                break

            if time > milestones[0]:
                print("Progress: {:.1f}%".format(100*time/(t1-t0)))
                milestones = np.delete(milestones, 0)
            numIters += 1
            if numIters > 100000:
                break
    
    def RK4OneStepError(self, time, dt):
        """
        A. S. Chai. 1968. 
        Error estimate of a fourth-order Runge-Kutta method with only one initial derivative evaluation. 
        In Proceedings of the April 30--May 2, 1968, spring joint computer conference (AFIPS '68 (Spring)). 
        Association for Computing Machinery, New York, NY, USA, 467–471. 
        https://doi.org/10.1145/1468075.1468144
        
        """

        # Back up previous temperatures to allow for RK4 to occur
        for station in self.stations:
            station[0].backupVariables()
        
        
        # Step 0 - generates k0
        self.calculateDerivatives(time, dt)
        for station in self.stations:
            # Feed in T_k1
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0: dt/2
                })

        # Step 1 - generates k1
        self.calculateDerivatives(time+dt/2, dt)
        for station in self.stations:
            # Feed in T_k2
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0: dt/4,
                    1: dt/4
                })

        # Step 2 - generates k2
        self.calculateDerivatives(time+dt/2, dt)
        for station in self.stations:
            # Feed in T_k3
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    1: 2*dt,
                    2:  -dt
                })

        # Step 3 - generates k3
        self.calculateDerivatives(time+dt, dt)
        for station in self.stations:
            # Feed in T_k4
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0:  7*dt/27,
                    1: 10*dt/27,
                    3:  1*dt/27
                })

        # Step 4 - generates k4
        self.calculateDerivatives(time+2*dt/3, dt)
        for station in self.stations:
            # Feed in T_k5
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0:   28*dt/625,
                    1: -125*dt/625,
                    2:  546*dt/625,
                    3:   54*dt/625,
                    4: -378*dt/625
                })

        maxError = 0
        worstComponent = ""
        # Step 5 - generates k5
        self.calculateDerivatives(time+dt/5, dt)
        for station in self.stations:

            # Calculate error
            station[0].updateVariables(
                variableHistoryValues = {}, updatePlusOne=True,
                intermediateDerivativeHistoryValues={
                    0:   dt/6,
                    2: 4*dt/6,
                    3:   dt/6,
                })

            var_np1_bar = station[0].updateVariables(
                variableHistoryValues = {}, updatePlusOne=False, returnValue=True,
                intermediateDerivativeHistoryValues={
                    0:  14*dt/336,
                    3:  35*dt/336,
                    4: 162*dt/336,
                    5: 125*dt/336
                })
            
            error = max(station[0].calculateErrorValueInput(var_np1_bar))
            if error > maxError:
                maxError = error * 1.0
                #worstComponent = station[0].state.ID

        #print(time, dt, maxError, worstComponent)
        if -np.log10(maxError) > self.precision:
            return(True, maxError)
        else:
            return(False, maxError)

    def RK4OneStepGill(self, time, dt):
        """
        S. Gill 
        A process for the step-by-step integration of differential equations in an automatic digital computing machine. 
        Mathematical Proceedings of the Cambridge Philosophical Society. 
        1951;47(1):96-108. 
        doi:10.1017/S0305004100026414
        
        """
        # Back up previous temperatures to allow for RK4 to occur
        for station in self.stations:
            station[0].backupVariables()
        
        # Step 0 - generates k0
        self.calculateDerivatives(time, dt)
        for station in self.stations:
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0: dt/2
                })

        # Step 1 - generates k1
        self.calculateDerivatives(time+dt/2, dt)
        for station in self.stations:
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0: (-1 + np.sqrt(2)) * dt/2,
                    1: ( 2 - np.sqrt(2)) * dt/2
                })
        
        # Step 2 - generates k2
        self.calculateDerivatives(time+dt/2, dt)
        for station in self.stations:
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    1: (   -np.sqrt(2)) * dt/2,
                    2: (1 + np.sqrt(2)) * dt/2
                })
        
        # Step 3 - generates k3
        self.calculateDerivatives(time+dt, dt)
        for station in self.stations:
            station[0].updateVariables(
                variableHistoryValues = {}, updatePlusOne=True, 
                intermediateDerivativeHistoryValues={
                    0:                    dt/6,
                    1: (2 - np.sqrt(2)) * dt/6,
                    2: (2 + np.sqrt(2)) * dt/6,
                    3:                    dt/6,
                })
        return(True, 0)
        
    def PredictorCorrectorCraneError(self, time, dt):
        """
        R. L. Crane and R. W. Klopfenstein. 1965. 
        A Predictor-Corrector Algorithm with an Increased Range of Absolute Stability. 
        J. ACM 12, 2 (April 1965), 227–241. 
        https://doi.org/10.1145/321264.321272
        """
        
        # Back up previous temperatures to allow for RK4 to occur
        for station in self.stations:
            if station[0].nVarHistory < 3 or station[0].nDerivHistory < 3:
                return(False, 1000)
            station[0].backupVariables()
        
        # Carry out predictor step
        self.calculateDerivatives(time, dt)
        for station in self.stations:
            # This is y'_n

            station[0].updateVariables(
                prevValue = 1.54765200,
                variableHistoryValues = {
                    -1: -1.86750300,
                    -2:  2.01720400,
                    -3: -0.69735300
                }, 
                intermediateDerivativeHistoryValues={
                    -1:  2.00224700 * dt
                },
                derivativeHistoryValues={
                    -1: -2.03169000 * dt,
                    -2:  1.81860900 * dt,
                    -3: -0.71432000 * dt 
                })
        
        self.calculateDerivatives(time+dt, dt)
        maxError = 0
        for station in self.stations:
            # This is p'_n+1

            station[0].updateVariables(
                variableHistoryValues = {}, updatePlusOne=True, returnValue=False, 
                intermediateDerivativeHistoryValues={
                    -1:  0.375000000 * dt,
                    -2:  0.791666667 * dt
                },
                derivativeHistoryValues={
                    -1: -0.208333333 * dt,
                    -2:  0.041666667 * dt,
               })
            maxError = max(max(station[0].calculateErrorPreviousValues()), maxError)

        #print(time, dt, maxError)
        if -np.log10(maxError) > self.precision:
            return(True, maxError)
        else:
            return(False, maxError)

    def systemController(self, time, dt):
        pass
