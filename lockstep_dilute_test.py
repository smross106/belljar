import belljar_lockstep as bjls


def bell(time):
    return(295+time)

loop = bjls.TimestepSolver(precision=4, maxTimestep=0.5)

loop.components["inlet"] = bjls.FluidSource(mDot=0.01, fluid="Nitrogen", T=bell, p=1, dilutesFlow={"H2O":0.01*1e-4})
loop.components["vent in"] = bjls.FluidLink()
loop.components["vent out"] = bjls.FluidLink(mDot=1)
loop.components["outlet"] = bjls.FluidSink()

loop.connectStationsFluid("inlet", "fluidOutlet", "vent in", "fluidInlet",fluid="Nitrogen",volume=1e-3)
loop.connectStationsFluid("vent in", "fluidOutlet", "vent out", "fluidInlet", fluid="Nitrogen",volume=0.1)
loop.connectStationsFluid("vent out", "fluidOutlet", "outlet", "fluidInlet", fluid="Nitrogen",volume=1e-3)

loop.initialiseSolver(305, mDot=0.01, dilutes={"H2O":0})

#loop.components["inlet"].fluidOutlet.state.T = 295

for state in loop.states:
    print(state.T, state.mDot, state.mass)

loop.mode = "RK4ERROR"
loop.runSolver(0, 50, 0.01)

for state in loop.states:
    print(state.T, state.mDot, state.dDilutes, state.dilutesHistory["H2O"][-7:-1])

print(loop.components["inlet"].dilutesFlow)

import matplotlib.pyplot as plt

plt.plot(loop.time, loop.components["vent in"].fluidInlet.state.THistory)
plt.plot(loop.time, loop.components["vent in"].fluidOutlet.state.THistory)
#plt.plot(loop.time, loop.components["vent out"].fluidOutlet.state.THistory)

plt.show()


exit()

loop.components["heater"] = bjls.FluidHeatTransfer(1)
loop.components["hx1 hot"] = bjls.FluidHeatTransfer(0.5)
loop.components["hx1 wall"] = bjls.ThermalConductance(5)
loop.components["hx1 cool"] = bjls.FluidHeatTransfer(0.5)
loop.components["cooler"] = bjls.FluidHeatTransfer(1)
loop.components["hx to cooler"] = bjls.FluidLink()
loop.components["hx to heater"] = bjls.FluidLink()
loop.components["heat source"] = bjls.ConstantTemperatureBoundary(320)
loop.components["cool source"] = bjls.HeatInputBoundary(-40)
loop.components["CO2 source"] = bjls.DiluteSourceSink({"CO2":1e-5})
#loop.components["CO2 scrub"] = bjls.DiluteSourceSink({"CO2":bell})
#loop.components["cool source"] = cs2.ConstantTemperatureBoundary(280)

loop.connectStationsFluid("heater", "fluidOutlet", "hx1 hot", "fluidInlet", dilutes={"CO2":1e-5, "H2O":0})
loop.connectStationsFluid("hx1 hot", "fluidOutlet", "hx to cooler", "fluidInlet", dilutes={"CO2":1e-5, "H2O":0})
loop.connectStationsFluid("hx to cooler", "fluidOutlet", "cooler", "fluidInlet", dilutes={"CO2":1e-5, "H2O":0})
loop.connectStationsFluid("cooler", "fluidOutlet", "hx1 cool", "fluidInlet", dilutes={"CO2":1e-5, "H2O":0})
loop.connectStationsFluid("hx1 cool", "fluidOutlet", "hx to heater", "fluidInlet", dilutes={"CO2":1e-5, "H2O":0})
loop.connectStationsFluid("hx to heater", "fluidOutlet", "heater", "fluidInlet", dilutes={"CO2":1e-5, "H2O":0})

#loop.connectStationsFluid("heater", "fluidInlet", "CO2 source", "fluidOutlet", dilutes={"CO2":1e-5, "H2O":0})
loop.connectStationsFluid("CO2 source", "fluidOutlet", "heater", "fluidInlet", dilutes={"CO2":1e-5, "H2O":0})

loop.connectStationsSolid("hx1 hot", "solidLink", "hx1 wall", "inlet", specificHeatCapacity=1)
loop.connectStationsSolid("hx1 wall", "outlet", "hx1 cool", "solidLink", specificHeatCapacity=1)
loop.connectStationsSolid("heater", "solidLink", "heat source", "outlet", specificHeatCapacity=10)
loop.connectStationsSolid("cooler", "solidLink", "cool source", "outlet", specificHeatCapacity=10)

loop.initialiseSolver(298.0, mDot=0.01)

print(len([i for i in loop.states if type(i) == bjls.FluidState]))

co2SumBefore = 0
for s in loop.states:
    if type(s) == bjls.FluidState:
        #print(c, loop.components[c].fluidOutlet.state.dilutes, loop.components[c].fluidOutlet.state.dilutesHistory)
        co2SumBefore += s.dilutes["CO2"] 
print("before", co2SumBefore)

loop.mode="RK4GILL"
loop.runSolver(0, 100, 0.1)

#print(loop.components["hx1 hot"].fluidInlet.state.dilutesHistory)
#print(loop.components["hx1 hot"].fluidOutlet.state.dilutesHistory)


import matplotlib.pyplot as plt

co2Sum = 0

for t in range(2):
    tSum = 0
    for c in loop.components:
        if hasattr(loop.components[c], "fluidOutlet"): 
            tSum += loop.components[c].fluidOutlet.state.dilutesHistory["CO2"][t]
    print(t, tSum)

for s in loop.states:
    if type(s) == bjls.FluidState:
        co2Sum += s.dilutes["CO2"] 
        #plt.plot(loop.time, loop.components[c].fluidOutlet.state.THistory, label=c)
        plt.plot(loop.time, s.dilutesConcentrationHistory["CO2"])
#plt.plot(loop.time, loop.components["hx1 hot"].fluidOutlet.state.dilutesHistory["CO2"], label="co2")
#plt.plot(loop.time, loop.components["hx1 hot"].fluidInlet.state.dilutesHistory["H2O"], label="h2o")
plt.legend()
print("after", co2Sum)
print("delta", co2Sum - co2SumBefore, (co2Sum-co2SumBefore)/(100*1e-5))

plt.show()