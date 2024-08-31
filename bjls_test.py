import belljar_lockstep as bjls
import numpy as np


solver = bjls.TimestepSolver(precision=2, maxTimestep=10)

solver.components["Source"] = bjls.FluidSource(mDot=0.01, p=1, T=298, fluid="Air")

solver.components["P1"] = bjls.GeneralFluidLink(flowModel=True)


solver.components["J1"] = bjls.FluidJunction(1, 2, [0.5, 0.5])

solver.components["P2"] = bjls.GeneralFluidLink(flowModel=True)
solver.components["H1"] = bjls.GeneralFluidLink(flowModel=False, 
                                                convectionModel={"active":True, "conductance":1, "useHTC":False, "useLMTD":False})
solver.components["HS1"] = bjls.ConstantTemperatureBoundary(318)

solver.components["P3"] = bjls.GeneralFluidLink(flowModel=True)
solver.components["C1"] = bjls.GeneralFluidLink(flowModel=False, 
                                                convectionModel={"active":True, "conductance":1, "useHTC":False, "useLMTD":False})
solver.components["CS1"] = bjls.ConstantTemperatureBoundary(278)

solver.components["J2"] = bjls.FluidJunction(2, 1)
solver.components["P4"] = bjls.GeneralFluidLink(flowModel=True)
solver.components["Sink"] = bjls.FluidSink()
#solver.components["Sink2"] = bjls.FluidSink()


solver.connectStationsFluid("Source", "fluidOutlet", "P1", "fluidInlet", fluid="Air")
solver.connectStationsFluid("P1", "fluidOutlet", "J1", "fluidInlet0", fluid="Air")

solver.connectStationsFluid("J1", "fluidOutlet0", "P2", "fluidInlet", fluid="Air")

solver.connectStationsFluid("H1", "fluidInlet", "P2", "fluidInlet", fluid="Air")
solver.connectStationsFluid("H1", "fluidOutlet", "P2", "fluidOutlet", fluid="Air")
solver.connectStationsSolid("H1", "solidLink", "HS1", "outlet")

solver.connectStationsFluid("J1", "fluidOutlet1", "P3", "fluidInlet", fluid="Air")

solver.connectStationsFluid("C1", "fluidInlet", "P3", "fluidInlet", fluid="Air")
solver.connectStationsFluid("C1", "fluidOutlet", "P3", "fluidOutlet", fluid="Air")
solver.connectStationsSolid("C1", "solidLink", "CS1", "outlet")

solver.connectStationsFluid("P2", "fluidOutlet", "J2", "fluidInlet0", fluid="Air")
solver.connectStationsFluid("P3", "fluidOutlet", "J2", "fluidInlet1", fluid="Air")
solver.connectStationsFluid("J2", "fluidOutlet0", "P4", "fluidInlet", fluid="Air")
solver.connectStationsFluid("P4", "fluidOutlet", "Sink", "fluidInlet", fluid="Air")

solver.components["J1"].fluidInlet0.state.ID = "J1 in"
solver.components["J1"].fluidOutlet0.state.ID = "J1 out 1"
solver.components["J1"].fluidOutlet1.state.ID = "J1 out 2"


solver.initialiseSolver(Tref=298, mDot=.01, p_bar=1)
#for o_i, o in enumerate(solver.components["J1"].fluidOutlets):
#    o.state.mDot *= solver.components["J1"].outputSplitRatio[o_i]

for state in solver.states:
    print(state.T)

solver.runSolver(0, 20, 0.1)

for c in solver.components:
    if "H" not in c:
        try:
            print(c, solver.components[c].fluidOutlet.state.T)
        except: pass

print("J2", solver.components["J2"].fluidOutlet0.state.T)
print(solver.components["H1"].convectionHeatTransfer)
print(solver.components["C1"].convectionHeatTransfer)
