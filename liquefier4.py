from components import *

lh = CycleSolver("Argon")

lh.components["so1"] = FluidSource_TP(1.0, 1, 300)
lh.components["jn1"] = FlowJunction()
lh.components["cp1"] = Compressor()
lh.components["hx1"] = Cooler(mode="T")
lh.components["hx2"] = Recuperator()
lh.components["th1"] = IdealThrottle()
lh.components["ps1"] = PhaseSeperator()
lh.components["th2"] = IdealThrottle()
lh.components["si1"] = FluidSink()

lh.connectStations(lh.components["so1"].outlet, lh.components["jn1"].inlet1)
lh.connectStations(lh.components["jn1"].outlet, lh.components["cp1"].inlet)
lh.connectStations(lh.components["cp1"].outlet, lh.components["hx1"].inlet)
lh.connectStations(lh.components["hx1"].outlet, lh.components["hx2"].inlet1)
lh.connectStations(lh.components["hx2"].outlet1, lh.components["th1"].inlet)
lh.connectStations(lh.components["th1"].outlet, lh.components["ps1"].inlet)
lh.connectStations(lh.components["ps1"].liquidOutlet, lh.components["si1"].inlet)
lh.connectStations(lh.components["ps1"].gasOutlet, lh.components["th2"].inlet)
lh.connectStations(lh.components["th2"].outlet, lh.components["hx2"].inlet2)
lh.connectStations(lh.components["hx2"].outlet2, lh.components["jn1"].inlet2)


def computeCycle(self):
    #self.components["so1"].compute()
    self.components["cp1"].compute(200, 1.0)
    self.components["hx1"].compute(TOut=350)

    self.components["hx2"].computeQ(eta=0.99)
    self.components["hx2"].computeStream1()

    self.components["th1"].compute(POut=1.2)
    self.components["ps1"].compute(relaxLiquidRate=0.9)
    self.components["th2"].compute(POut=1)
    self.components["hx2"].computeStream2()
    

    self.components["jn1"].compute()

lh.computeCycle = computeCycle.__get__(lh)

lh.initialiseSolver(300, 1, 1.0/0.5)

#for i in range(100):
#    lh.computeCycle()

lh.solveCycleRelaxationFactor(relaxationFactor=1, maxIterations=500, limitAcceleration=None, limitVelocity=None)

print("phase sep")
print(lh.components["ps1"].inlet.state.Q)
print(lh.components["ps1"].inlet.flow.mDot)
print(lh.components["ps1"].liquidOutlet.flow.mDot)
print(lh.components["ps1"].gasOutlet.flow.mDot)

lh.components["ps1"].compute(relaxLiquidRate=1.0)

print("phase sep")
print(lh.components["ps1"].inlet.state.Q)
print(lh.components["ps1"].inlet.flow.mDot)
print(lh.components["ps1"].liquidOutlet.flow.mDot)
print(lh.components["ps1"].gasOutlet.flow.mDot)

fig, ax1 = plt.subplots(figsize=(9, 6))
ax2 = ax1.twinx()  

ax1.plot(lh.residualHistory[:,0], lh.residualHistory[:,2], linestyle="dotted")
ax1.set_yscale("log")

ax2.plot(lh.residualHistory[:,0], [i["ps1"].liquidOutlet.flow.mDot for i in lh.oldComponentStates],c="r")
ax2.plot(lh.residualHistory[:,0], [i["ps1"].gasOutlet.flow.mDot for i in lh.oldComponentStates],c="r", linestyle=":")
ax2.plot(lh.residualHistory[:,0], [i["ps1"].inlet.flow.mDot for i in lh.oldComponentStates],c="r", linestyle="--")
ax2.plot(lh.residualHistory[:,0], [i["ps1"].inlet.state.Q for i in lh.oldComponentStates], c="g")

plt.show()