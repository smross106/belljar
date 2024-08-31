from components import *

cl = CycleSolver("Nitrogen")

cl.components["so1"] = FluidSource_TP(1, 1, 300)
cl.components["jn1"] = FlowJunction()
cl.components["cp1"] = Compressor()
cl.components["co1"] = Cooler(mode="T")
cl.components["ps1"] = PhaseSeperator()
cl.components["siL"] = FluidSink()
cl.components["th1"] = IdealThrottle()
cl.components["hx1"] = Recuperator()


cl.connectStations(cl.components["so1"].outlet, cl.components["jn1"].inlet1)
cl.connectStations(cl.components["jn1"].outlet, cl.components["cp1"].inlet)
cl.connectStations(cl.components["cp1"].outlet, cl.components["co1"].inlet)
cl.connectStations(cl.components["co1"].outlet, cl.components["hx1"].inlet1)
cl.connectStations(cl.components["hx1"].outlet1, cl.components["th1"].inlet)
cl.connectStations(cl.components["th1"].outlet, cl.components["ps1"].inlet)
cl.connectStations(cl.components["ps1"].gasOutlet, cl.components["hx1"].inlet2)
cl.connectStations(cl.components["ps1"].liquidOutlet, cl.components["siL"].inlet)
cl.connectStations(cl.components["hx1"].outlet2, cl.components["jn1"].inlet2)


def computeCycle(self):
    self.components["so1"].compute()
    self.components["cp1"].compute(200, 1)
    self.components["co1"].compute(TOut=350)
    self.components["hx1"].computeQ(eta=1)
    self.components["hx1"].computeStream1()
    self.components["th1"].compute(POut=1)
    self.components["ps1"].compute()
    self.components["siL"].compute()
    self.components["hx1"].computeStream2()
    self.components["jn1"].compute()

cl.computeCycle = computeCycle.__get__(cl)

cl.initialiseSolver(273, 1, 2)

for i in range(100):
    cl.computeCycle()
    print(i, cl.components["ps1"].inlet.flow.mDot, cl.components["ps1"].liquidOutlet.flow.mDot, cl.components["ps1"].gasOutlet.flow.mDot, cl.components["ps1"].prevQ)
    print("\t", cl.components["th1"].inlet.state.T, cl.components["th1"].outlet.state.T, cl.components["th1"].outlet.state.Q)

for c in cl.components.keys():
    print(cl.components[c])
