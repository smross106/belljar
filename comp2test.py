from components2 import *

cycle = CycleSolver("R134A")

cycle.components["co"] = Compressor()
cycle.components["cn"] = Cooler("T")
cycle.components["hx"] = Recuperator()
cycle.components["th"] = IdealThrottle()
cycle.components["ev"] = Evaporator()

cycle.connectStations("co", "outlet", "cn", "inlet")
cycle.connectStations("cn", "outlet", "hx", "inlet1")
cycle.connectStations("hx", "outlet1", "th", "inlet")
cycle.connectStations("th", "outlet", "ev", "inlet")
cycle.connectStations("ev", "outlet", "hx", "inlet2")
cycle.connectStations("hx", "outlet2", "co", "inlet")

def computeCycle(self):
    self.components["hx"].compute(eta=1)
    self.components["hx"].computeStream2()
    self.components["co"].compute(pOut=10, eta=0.75)
    self.components["cn"].compute(TOut=36+273)
    self.components["hx"].computeStream1()
    self.components["th"].compute(pOut=1.333)
    self.components["ev"].compute()

cycle.computeCycle = computeCycle.__get__(cycle)

cycle.initialiseSolver(298, 1.333, 1/60)

for i in range(40):
    
    cycle.computeCycle()

states = [["ev", "outlet"], ["co", "inlet"], ["co", "outlet"], ["cn", "outlet"],["hx", "outlet1"],["th", "outlet"]]

for s in states:
    val = getattr(cycle.components[s[0]], s[1]).state
    print(s[0], "{:.2f}K {:.2f}bar \t {:.2f}kJ/kg {:.2f}kJ/kg-K {:.2f}".format(val.T, val.p, val.h/1000, val.s/1000, val.Q))