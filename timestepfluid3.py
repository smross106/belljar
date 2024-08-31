import crane_scratch2 as cs2
from CoolProp.CoolProp import PropsSI
import numpy as np

class FluidState(object):
    def __init__(self, TRef=298, fluid="Nitrogen", mDot=1, p_bar=1):
        self._T = TRef
        self.fluid = fluid
        self.mDot = mDot
        self.p = p_bar
        self._h = PropsSI("Hmass", "T", self.T, "P", self.p*1e5, self.fluid)

        self.boundaryConditionT = False
        self.boundaryConditionh = False

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
    
    @property
    def T(self):
        return self._T
    
    @T.setter
    def T(self, value):
        if not self.boundaryConditionT:
            self._T = value * 1.0
        else:
            self._T = self.boundaryConditionT
        
        self._h = PropsSI("Hmass", "T", self.T, "P", self.p*1e5, self.fluid)
    
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
            self.T = PropsSI("T", "Hmass", self.h, "P", self.p*1e5, self.fluid)
        except:
            print(self.h, self.hnm1_temp, self.dhHistory, self.intermediatedhdtHistory)
            exit()

        
    def computedhdt(self, time):
        self.dhdt = self.storedE / self.mDot
        self.storedE = 0

    def setNexth(self, dt):
        self.hHistory.append(self.hnm1_temp)
        self.dhHistory.append((self.hnp1 - self.hnm1_temp)/dt)
        self.THistory.append(self.T)
        self.h = self.hnp1
        self.dTHistory.append((self.T - self.THistory[-1])/dt)
        self.hnp1 = 0
        self.intermediatehHistory = []
        self.intermediatedhdtHistory = []
    
    def resetSolver(self):
        self.h = self.hnm1_temp
        self.hnp1 = 0
        self.intermediatehHistory = []
        self.intermediatedhdtHistory = []

    
    def stateInitialise(self, TRef, mDot, p):
        self.THistory = []
        self.intermediateTHistory = []
        self.intermediatedTdtHistory = []
        self.Tnm1_temp = 0
        self.Tnp1 = 0 
        self.T = TRef
        self.mDot = mDot
        self.p = p
        self.storedE = 0

class FluidStation(object):
    def __init__(self):
        self.state = None

        self.nVarHistory = 0
        self.nDerivHistory = 0
    
    def computeDerivatives(self, time):
        self.state.computedhdt(time)
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
        return(errors)

    def calculateErrorPreviousValues(self):
        errors = []
        errors.append(abs(self.state.hnp1 - self.state.hnm1_temp)/self.state.hnm1_temp)
        return(errors)

class FluidLink(cs2.SolidComponent):
    def __init__(self, volume=0, mDot=0):
        super().__init__()
        self.volume = volume
        self.mDot = 0
        
        self.fluidInlet = FluidStation()
        self.fluidOutlet = FluidStation()

        self.linkedStations = [self.fluidInlet, self.fluidOutlet]
    
    def compute(self, time):
        self.mDot = self.fluidInlet.state.mDot

        self.fluidOutlet.state.mDot =   self.fluidInlet.state.mDot
        self.fluidOutlet.state.p =      self.fluidInlet.state.p

        self.fluidInlet.state.storedE -= self.fluidInlet.state.mDot * self.fluidInlet.state.h
        self.fluidOutlet.state.storedE += self.fluidInlet.state.mDot * self.fluidInlet.state.h

class FluidHeatTransfer(cs2.SolidComponent):
    def __init__(self, heatTransferCoefficient):
        super().__init__()
        self.mDot = 0
        self.solidHeatTransfer = 0
        self.heatTransferCoefficient = heatTransferCoefficient
        
        self.fluidInlet = FluidStation()
        self.fluidOutlet = FluidStation()
        self.solidLink = cs2.SolidStation()

        self.linkedStations = [self.fluidInlet, self.fluidOutlet, self.solidLink]
    
    def compute(self, time):
        self.mDot = self.fluidInlet.state.mDot

        self.fluidOutlet.state.mDot =   self.fluidInlet.state.mDot
        self.fluidOutlet.state.p =      self.fluidInlet.state.p

        # Theoretical maximum conduction heat
        conductionHeat = self.heatTransferCoefficient * (
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

def connectStationsFluid(self, component1, station1, component2, station2, fluid="Nitrogen", p_bar=1, mDot=1, TRef=298):
            
        fs = FluidState(fluid=fluid, p_bar=p_bar, mDot=mDot, TRef=TRef)

        #self.states.append(fs)
        setattr(getattr(self.components[component1], station1), "state", fs)
        setattr(getattr(self.components[component2], station2), "state", fs)

def systemController(self, time, dt):
    # Aim for a target temperature at the heater output by varying the cooler power
    heatError = self.components["heater"].fluidInlet.state.T - self.controlParameters["Target T"]

    self.controlParameters["Integral error"] += heatError * dt
    
    prop = heatError * -self.controlParameters["Proportional"] + (
        self.controlParameters["Derivative"] * self.components["heater"].fluidInlet.state.dTHistory[-1]) + (
            self.controlParameters["Integral error"] * -self.controlParameters["Integral"]
        )

    if abs(prop/dt) > self.controlParameters["Max cooler Q ramp"]:
        prop = self.controlParameters["Max cooler Q ramp"]*dt*np.sign(prop)

    targetQ = self.components["cool source"].QIn + prop
    if targetQ>self.controlParameters["Max cooler Q"]:
        targetQ = self.controlParameters["Max cooler Q"]
    if targetQ<self.controlParameters["Min cooler Q"]:
        targetQ = self.controlParameters["Min cooler Q"]

    self.components["cool source"].QIn = targetQ




loop = cs2.SolverSolid(precision=3)


loop.connectStationsFluid = connectStationsFluid.__get__(loop)
loop.systemController = systemController.__get__(loop)

loop.components["heater"] = FluidHeatTransfer(1)
loop.components["hx1 hot"] = FluidHeatTransfer(0.5)
loop.components["hx1 wall"] = cs2.ThermalConductance(5)
loop.components["hx1 cool"] = FluidHeatTransfer(0.5)
loop.components["cooler"] = FluidHeatTransfer(1)
loop.components["hx to cooler"] = FluidLink()
loop.components["hx to heater"] = FluidLink()
loop.components["heat source"] = cs2.ConstantTemperatureBoundary(320)
loop.components["cool source"] = cs2.HeatInputBoundary(-40)
#loop.components["cool source"] = cs2.ConstantTemperatureBoundary(280)

loop.connectStationsFluid("heater", "fluidOutlet", "hx1 hot", "fluidInlet")
loop.connectStationsFluid("hx1 hot", "fluidOutlet", "hx to cooler", "fluidInlet")
loop.connectStationsFluid("hx to cooler", "fluidOutlet", "cooler", "fluidInlet")
loop.connectStationsFluid("cooler", "fluidOutlet", "hx1 cool", "fluidInlet")
loop.connectStationsFluid("hx1 cool", "fluidOutlet", "hx to heater", "fluidInlet")
loop.connectStationsFluid("hx to heater", "fluidOutlet", "heater", "fluidInlet")
loop.connectStationsSolid("hx1 hot", "solidLink", "hx1 wall", "inlet", specificHeatCapacity=1)
loop.connectStationsSolid("hx1 wall", "outlet", "hx1 cool", "solidLink", specificHeatCapacity=1)
loop.connectStationsSolid("heater", "solidLink", "heat source", "outlet", specificHeatCapacity=10)
loop.connectStationsSolid("cooler", "solidLink", "cool source", "outlet", specificHeatCapacity=10)

loop.initialiseSolver(298.0, mDot=0.001)
loop.mode="RK4ERROR"

loop.controlParameters = {"Target T":270, "Max cooler Q ramp":1, "Max cooler Q":0, "Min cooler Q":-100, 
                          "Proportional":0.005*0.6, "Derivative":0*0.01*3*60/40, "Integral":0.005*0.54/60., "Integral error":0}
loop.runSolver(0, 250, 0.005)

print("T heater in", loop.components["heater"].fluidInlet.state.T)
print("T heater out", loop.components["hx1 hot"].fluidInlet.state.T)
print("T heater wall", loop.components["heater"].solidLink.state.T)
print("T HX hot out", loop.components["hx1 hot"].fluidOutlet.state.T)
print("T HX hot wall", loop.components["hx1 hot"].solidLink.state.T)
print("T cooler in", loop.components["cooler"].fluidInlet.state.T)
print("T cooler out", loop.components["cooler"].fluidOutlet.state.T)
print("T HX cold out", loop.components["hx1 cool"].fluidOutlet.state.T)
print("T HX cold wall", loop.components["hx1 cool"].solidLink.state.T)

print("")
print("Q cooler solid", loop.components["cooler"].solidHeatTransfer)
print("Q cooler fluid", loop.components["cooler"].mDot * (loop.components["cooler"].fluidOutlet.state.h - loop.components["cooler"].fluidInlet.state.h))

print("Q heater solid", loop.components["heater"].solidHeatTransfer)
print("Q heater fluid", loop.components["heater"].mDot * (loop.components["heater"].fluidOutlet.state.h - loop.components["heater"].fluidInlet.state.h))

print("")
print("Q hx1 hot solid", loop.components["hx1 hot"].solidHeatTransfer)
print("Q hx1 hot fluid", loop.components["hx1 hot"].mDot * (loop.components["hx1 hot"].fluidOutlet.state.h - loop.components["hx1 hot"].fluidInlet.state.h))

print("Q hx1 cool solid", loop.components["hx1 cool"].solidHeatTransfer)
print("Q hx1 hot fluid", loop.components["hx1 cool"].mDot * (loop.components["hx1 cool"].fluidOutlet.state.h - loop.components["hx1 cool"].fluidInlet.state.h))
        
import matplotlib.pyplot as plt

plt.plot(loop.time, loop.components["heater"].fluidInlet.state.THistory, label="Heater inlet", linestyle="--")
plt.plot(loop.time, loop.components["hx1 hot"].fluidInlet.state.THistory, label="Hot")
plt.plot(loop.time, loop.components["hx1 hot"].solidLink.state.THistory, label="Wall hot")
plt.plot(loop.time, loop.components["hx1 cool"].solidLink.state.THistory, label="Wall cold")
plt.plot(loop.time, loop.components["hx1 cool"].fluidInlet.state.THistory, label="Cold")
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("Temperature evolution in heat exchanger")

print("Hot out", loop.components["hx1 hot"].fluidInlet.state.T)
print("Cold hot", loop.components["hx1 cool"].fluidInlet.state.T)

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(loop.time, loop.timesteps, 'g-')
ax2.plot(loop.time, [i*100 for i in loop.errors], 'b-')

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Timestep (s)', color='g')
ax2.set_ylabel('Relative error in step, %', color='b')
ax1.set_title("Solver progress")
ax2.set_ylim(0, 0.2)
ax1.set_ylim(0, 0.5)

plt.show()