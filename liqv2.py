from components2 import *

cyc = CycleSolver("Argon")

# Linde-Hampson Cycle Argon Example

cyc.components["jn1"] = FlowJunction()
cyc.components["comp"] = Compressor()
cyc.components["cool"] = Cooler("T")
cyc.components["hx1"] = Recuperator()
cyc.components["th1"] = IdealThrottle()
cyc.components["phase"] = PhaseSeperator()
cyc.components["th2"] = IdealThrottle()
cyc.components["sink"] = FluidSink()
cyc.components["src"] = FluidSource_TP(1/60, 1, 288.15)

cyc.connectStations("src", "outlet", "jn1", "inlet1")
cyc.connectStations("jn1", "outlet", "comp", "inlet")
cyc.connectStations("comp", "outlet", "cool", "inlet")
cyc.connectStations("cool", "outlet", "hx1", "inlet1")
cyc.connectStations("hx1", "outlet1", "th1", "inlet")
cyc.connectStations("th1", "outlet", "phase", "inlet")
cyc.connectStations("phase", "liquidOutlet", "sink", "inlet")
cyc.connectStations("phase", "gasOutlet", "th2", "inlet")
cyc.connectStations("th2", "outlet", "hx1", "inlet2")
cyc.connectStations("hx1", "outlet2", "jn1", "inlet2")

def computeCycle(self):
    self.components["comp"].compute(pOut=200, eta=0.8)
    self.components["cool"].compute(TOut=288.15)
    self.components["hx1"].compute(eta=1.0)
    self.components["hx1"].computeStream1()
    self.components["th1"].compute(pOut=1.2)
    self.components["phase"].compute()
    self.components["th2"].compute(pOut=1.0)
    self.components["hx1"].computeStream2()
    self.components["jn1"].compute()

cyc.computeCycle = computeCycle.__get__(cyc)
cyc.initialiseSolver(288.15, 1, 2/60)
cyc.components["src"].compute()


for i in range(50):
    print(i, "\t", cyc.components["th1"].outlet.state.Q)
    #for c in cyc.components:
    #    print(cyc.components[c])

    cyc.computeCycle()

cyc.printStation(1, "src", "outlet")
cyc.printStation(2, "comp", "inlet")
cyc.printStation(3, "comp", "outlet")
cyc.printStation(4, "cool", "outlet")
cyc.printStation(5, "hx1", "outlet1")
cyc.printStation(6, "th1", "outlet")
cyc.printStation(7, "phase", "gasOutlet")
cyc.printStation(8, "th2", "outlet")
cyc.printStation(9, "hx1", "outlet2")

cyc.printStation(10, "phase", "liquidOutlet")

liqEnth = PropsSI("Hmass", "T", 273+15, "P", 1.0e5, "Argon") - PropsSI("Hmass", "Q", 0, "P", 1.0e5, "Argon")
dH = cyc.components["comp"].outlet.state.h - cyc.components["comp"].inlet.state.h
dE = cyc.components["comp"].WIn / cyc.components["comp"].inlet.flow.mDot


print("E", cyc.components["comp"].outlet.state.h - cyc.components["comp"].inlet.state.h)
print("W", (cyc.components["comp"].outlet.state.h - cyc.components["comp"].inlet.state.h)*cyc.components["comp"].inlet.flow.mDot)
print("W2", cyc.components["comp"].WIn)
print("dH", PropsSI("Hmass", "T", 298, "P", 1.0e5, "Argon") - PropsSI("Hmass", "Q", 0, "P", 1.0e5, "Argon"))
print("eta", liqEnth/dH)
print("eta", liqEnth/dE)