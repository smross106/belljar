from components import *

cl = CycleSolver("Nitrogen")

cl.components["so1"] = FluidSource_TP(1, 5, 273+15)
cl.components["jn1"] = FlowJunction()
cl.components["cp1"] = Compressor()
cl.components["co1"] = Cooler(mode="T")
cl.components["hx1"] = Recuperator()

cl.components["sp1"] = FlowSplitter()
cl.components["ex1"] = Turbine()
cl.components["jn2"] = FlowJunction()

cl.components["hx2"] = Recuperator()
cl.components["hx3"] = Recuperator()

cl.components["th1"] = IdealThrottle()
cl.components["ps1"] = PhaseSeperator()
cl.components["th2"] = IdealThrottle()
cl.components["si1"] = FluidSink()

cl.connectStations(cl.components["so1"].outlet, cl.components["jn1"].inlet1)
cl.connectStations(cl.components["jn1"].outlet, cl.components["cp1"].inlet)
cl.connectStations(cl.components["cp1"].outlet, cl.components["co1"].inlet)
cl.connectStations(cl.components["co1"].outlet, cl.components["hx1"].inlet1)

cl.connectStations(cl.components["hx1"].outlet1, cl.components["sp1"].inlet)

cl.connectStations(cl.components["sp1"].outlet2, cl.components["ex1"].inlet)
cl.connectStations(cl.components["ex1"].outlet, cl.components["jn2"].inlet2)

cl.connectStations(cl.components["sp1"].outlet1, cl.components["hx2"].inlet1)
cl.connectStations(cl.components["hx2"].outlet1, cl.components["hx3"].inlet1)

cl.connectStations(cl.components["hx3"].outlet1, cl.components["th1"].inlet)
cl.connectStations(cl.components["th1"].outlet, cl.components["ps1"].inlet)
cl.connectStations(cl.components["ps1"].liquidOutlet, cl.components["si1"].inlet)
cl.connectStations(cl.components["ps1"].gasOutlet, cl.components["th2"].inlet)
cl.connectStations(cl.components["th2"].outlet, cl.components["hx3"].inlet2)
cl.connectStations(cl.components["hx3"].outlet2, cl.components["jn2"].inlet1)

cl.connectStations(cl.components["jn2"].outlet, cl.components["hx2"].inlet2)
cl.connectStations(cl.components["hx2"].outlet2, cl.components["hx1"].inlet2)
cl.connectStations(cl.components["hx1"].outlet2, cl.components["jn1"].inlet2)


def computeCycle(self):
    #self.components["so1"].compute()
    
    self.components["cp1"].compute(15, 0.6)
    self.components["co1"].compute(TOut=273+100)

    self.components["hx1"].computeQ(eta=0.8)
    self.components["hx1"].computeStream1()
    self.components["sp1"].compute(flowFraction2=0.6)
    
    

    self.components["hx2"].computeQ(eta=0.9)
    self.components["hx2"].computeStream1()
    self.components["hx3"].computeQ(eta=0.9)
    self.components["hx3"].computeStream1()

    self.components["th1"].compute(POut=5.01)
    self.components["ps1"].compute()
    
    
    self.components["si1"].compute()

    self.components["th2"].compute(POut=5)
    self.components["hx3"].computeStream2()
    self.components["ex1"].compute(pOut=5.0, eta=1.0)
    self.components["jn2"].compute()
    self.components["hx2"].computeStream2()
    self.components["hx1"].computeStream2()
    

    self.components["jn1"].compute()

    #self.components["ps1"].compute()
    

cl.computeCycle = computeCycle.__get__(cl)

cl.initialiseSolver(273+15, 6, 2)
cl.components["so1"].compute()
cl.components["ps1"].liquidOutlet.flow.mDot = 0
cycle = []
ps_temps = []
liq_flows = []
gas_flows = []


cl.solveCycleRelaxationFactor(relaxationFactor=1.0, maxIterations=100, limitAcceleration=None, limitVelocity=None)
#cl.solveCycleRelaxationFactor(relaxationFactor=1.0, maxIterations=100, limitAcceleration=None, limitVelocity=None)
cl.residualHistory = np.array(cl.residualHistory)
#cl.components["ps1"].compute(relaxLiquidRate=1.0)

print("phase sep")
print(cl.components["ps1"].inlet.state.Q)
print(cl.components["ps1"].inlet.flow.mDot)
print(cl.components["ps1"].liquidOutlet.flow.mDot)
print(cl.components["ps1"].gasOutlet.flow.mDot)
print("after phase sep throttle")
print(cl.components["th2"].outlet.state.Q)
print("hx3")
print(cl.components["hx3"].Qdot)
print(cl.components["hx3"].inlet1.state.T, cl.components["hx3"].inlet1.state.h)
print(cl.components["hx3"].outlet1.state.T, cl.components["hx3"].outlet1.state.h)
print(cl.components["hx3"].inlet2.state.Q, cl.components["hx3"].inlet2.state.h)
print(cl.components["hx3"].outlet2.state.Q, cl.components["hx3"].outlet2.state.h)

print("expander")
print(cl.components["ex1"].inlet.state.Q)
print(cl.components["ex1"].outlet.state.Q)
print("junction")
print(cl.components["jn1"].outlet.state.Q)

fig, ax1 = plt.subplots(figsize=(9, 6))
ax2 = ax1.twinx()  

ax1.plot(cl.residualHistory[:,0], cl.residualHistory[:,2], linestyle="dotted")
ax1.set_yscale("log")

#ax2.plot(cl.residualHistory[:,0], [i["ps1"].liquidOutlet.flow.mDot for i in cl.oldComponentStates], c="r", linestyle="dotted")
#ax2.plot(cl.residualHistory[:,0], [i["ps1"].gasOutlet.flow.mDot for i in cl.oldComponentStates], c="r", linestyle="dashed")
ax2.plot(cl.residualHistory[:,0], [i["ps1"].liquidOutlet.flow.mDot for i in cl.oldComponentStates],c="r")
ax2.plot(cl.residualHistory[:,0], [i["ps1"].gasOutlet.flow.mDot for i in cl.oldComponentStates],c="r", linestyle=":")
ax2.plot(cl.residualHistory[:,0], [i["ps1"].inlet.flow.mDot for i in cl.oldComponentStates],c="r", linestyle="--")
ax2.plot(cl.residualHistory[:,0], [i["ps1"].inlet.state.Q for i in cl.oldComponentStates], c="g")

ax2.plot(cl.residualHistory[:,0], [i["cp1"].inlet.state.Q for i in cl.oldComponentStates], c="k")
ax2.plot(cl.residualHistory[:,0], [i["ex1"].outlet.state.Q for i in cl.oldComponentStates], c="k", linestyle="--")

#ax2.plot(cl.residualHistory[:,0], [i["ps1"].inlet.state.Q for i in cl.oldComponentStates], c="r", linestyle="dashed")
plt.show()
#print("\n", len(cl.fluidStates))




exit()
for i in range(151):
    print("cycle", i)
    print(cl.components["ps1"].inlet.state.T, cl.components["ps1"].liquidOutlet.flow.mDot)
    print(cl.components["cp1"].outlet.state.T)
    cl.computeCycle(relaxComp=min(i/100,1))

    cycle.append(i)
    ps_temps.append(cl.components["ps1"].inlet.state.T)
    liq_flows.append(cl.components["ps1"].liquidOutlet.flow.mDot)
    gas_flows.append(cl.components["ps1"].gasOutlet.flow.mDot)

for c in cl.components.keys():
    try:
        print(c, cl.components[c].outlet.state.Q, cl.components[c].outlet.state.T)
    except:
        try:
            print(c, cl.components[c].outlet2.state.Q, cl.components[c].outlet2.state.T)
        except:
            pass
    print(cl.components[c])


print(cl.components["cp1"].inlet.state.Q)
print(cl.components["cp1"].outlet)
print(cl.components["cp1"].inlet.state.h, cl.components["cp1"].outlet.state.h, cl.components["cp1"].outlet.flow.mDot, cl.components["cp1"].WIn, cl.components["ps1"].liquidOutlet.flow.mDot)



try:
    workIn = cl.components["cp1"].WIn + cl.components["ex1"].WIn
    specificWork = workIn / cl.components["ps1"].liquidOutlet.flow.mDot

    Hfg = (PropsSI("Hmass", "P", 5.2e5, "T", 288, cl.fluid) - PropsSI("Hmass", "P", 5.2e5, "Q", 0, cl.fluid))
    print(workIn/1e3, specificWork/1e3, cl.components["co1"].QOut/(1e3*cl.components["co1"].inlet.flow.mDot), Hfg/1e3, Hfg/specificWork)
except:
    pass

fig, ax1 = plt.subplots(figsize=(9, 6))

# Instantiate a second axes that shares the same x-axis
ax2 = ax1.twinx()  
ax1.plot(cycle, liq_flows, linestyle="-", c="b", label="Liquid")
ax1.plot(cycle, gas_flows, linestyle="--", c="b", label="Gas")
ax1.set_ylabel("Flow rates (kg/s)", color="b")
ax1.legend()

ax2.plot(cycle, ps_temps, c="r")
ax2.set_ylabel("Outlet temperature (K)", color="r")
ax1.set_title("Zapitza Cycle (Methane, 5 bar inlet) Solver Progress")
ax1.set_xlabel("Solver cycle")

plt.show()