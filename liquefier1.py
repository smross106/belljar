from components import *


lh = CycleSolver("Nitrogen")

lh.components["so1"] = FluidSource_TP(0.1, 1, 15+273)
lh.components["jn1"] = FlowJunction()
lh.components["cp1"] = Compressor()
lh.components["co1"] = Cooler(mode="T")
lh.components["cp2"] = Compressor()
lh.components["hx1"] = Cooler(mode="T")
lh.components["hx2"] = Recuperator()
lh.components["th1"] = IdealThrottle()
lh.components["ps1"] = PhaseSeperator()
lh.components["th2"] = IdealThrottle()
lh.components["si1"] = FluidSink()

lh.connectStations(lh.components["so1"].outlet, lh.components["jn1"].inlet1)
lh.connectStations(lh.components["jn1"].outlet, lh.components["cp1"].inlet)
lh.connectStations(lh.components["cp1"].outlet, lh.components["co1"].inlet)
lh.connectStations(lh.components["co1"].outlet, lh.components["cp2"].inlet)
lh.connectStations(lh.components["cp2"].outlet, lh.components["hx1"].inlet)
lh.connectStations(lh.components["hx1"].outlet, lh.components["hx2"].inlet1)
lh.connectStations(lh.components["hx2"].outlet1, lh.components["th1"].inlet)
lh.connectStations(lh.components["th1"].outlet, lh.components["ps1"].inlet)
lh.connectStations(lh.components["ps1"].liquidOutlet, lh.components["si1"].inlet)
lh.connectStations(lh.components["ps1"].gasOutlet, lh.components["hx2"].inlet2)
lh.connectStations(lh.components["hx2"].outlet2, lh.components["th2"].inlet)
lh.connectStations(lh.components["th2"].outlet, lh.components["jn1"].inlet2)


def computeCycle(self):
    self.components["so1"].compute()
    self.components["jn1"].compute()
    self.components["cp1"].compute(17, 0.5)
    self.components["co1"].compute(TOut=273+15)
    self.components["cp2"].compute(300, 0.5)

    self.components["hx1"].compute(TOut=273+15)
    self.components["hx2"].computeQ(eta=0.9)
    self.components["hx2"].computeStream1()

    self.components["th1"].compute(POut=1.2)
    self.components["ps1"].compute(relaxLiquidRate=1.0)
    self.components["hx2"].computeStream2()
    self.components["th2"].compute(POut=1)

    self.components["jn1"].compute()
    
    
lh.computeCycle = computeCycle.__get__(lh)

lh.initialiseSolver(288, 1, 0.2)
#lh.components["so1"].outlet.state.T = 288
#lh.components["so1"].outlet.state.p = 1
#lh.components["so1"].outlet.flow.mDot = 0.1


cycle = []
ps_temps = []
liq_flows = []
gas_flows = []

for i in range(100):
    lh.computeCycle()
    #print(i, lh.components["jn1"].inlet1.flow.mDot, lh.components["jn1"].inlet2.flow.mDot, lh.components["jn1"].outlet.flow.mDot)
    #print(i, lh.components["ps1"].inlet.flow.mDot, lh.components["ps1"].gasOutlet.flow.mDot, lh.components["ps1"].liquidOutlet.flow.mDot)
    cycle.append(i)
    ps_temps.append(lh.components["ps1"].inlet.state.T)
    liq_flows.append(lh.components["ps1"].liquidOutlet.flow.mDot)
    gas_flows.append(lh.components["ps1"].gasOutlet.flow.mDot)

print("")
#print("W "+str(lh.components["cp1"].WIn/1000 + lh.components["cp2"].WIn/1000)+"kW")
#print("Q "+str(-lh.components["hx1"].QOut/1000)+"kW")
#print("E "+str((lh.components["cp1"].WIn+lh.components["cp2"].WIn)/lh.components["ps1"].liquidOutlet.flow.mDot / 1000)+"kJ/kg")

for c in lh.components.keys():
    print(lh.components[c])


fig, ax1 = plt.subplots(figsize=(9, 6))

# Instantiate a second axes that shares the same x-axis
ax2 = ax1.twinx()  
ax1.plot(cycle, liq_flows, linestyle="-", c="b", label="Liquid")
ax1.plot(cycle, gas_flows, linestyle="--", c="b", label="Gas")
ax1.set_ylabel("Flow rates (kg/s)", color="b")
ax1.legend()

ax2.plot(cycle, ps_temps, c="r")
ax2.set_ylabel("Outlet temperature (K)", color="r")
ax1.set_title("Linde-Hampson Cycle (Nitrogen) Solver Progress")
ax1.set_xlabel("Solver cycle")

plt.show()
