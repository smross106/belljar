import belljar_fruitloop as bjfl

cyc = bjfl.CycleSolver(fluid="Oxygen")

cyc.components["ox inlet"] = bjfl.FluidSource_TP(mDot=1, T=80, p=2)
cyc.components["cooler"] = bjfl.Cooler(mode="Q")
cyc.components["orifice"] = bjfl.CompressibleOrifice()
cyc.components["pump"] = bjfl.Compressor()
cyc.components["ox outlet"] = bjfl.FluidSink()

cyc.connectStations("ox inlet", "outlet", "cooler", "inlet")
cyc.connectStations("cooler", "outlet", "orifice", "inlet")
cyc.connectStations("orifice", "outlet", "pump", "inlet")
cyc.connectStations("pump", "outlet", "ox outlet", "inlet")

def computeCycle(self, pOutlet, CdA, coolQ, pumpPR, forward=True):
    self.components["ox outlet"].inlet.state.update_Tp(self.components["ox outlet"].inlet.state.T, pOutlet)

    self.components["cooler"].compute(QOut=coolQ, forward=False)
    self.components["orifice"].compute(CdA, forward=forward)
    self.components["pump"].compute(pOut=None, pressureRatio=pumpPR, forward=True)

cyc.computeCycle = computeCycle.__get__(cyc)

cyc.initialiseSolver(80, 1, 1)
cyc.components["ox inlet"].compute()

for i in range(10):
    print(i,cyc.components["cooler"].inlet.flow.mDot, cyc.components["cooler"].outlet.flow.mDot, cyc.components["pump"].inlet.flow.mDot, cyc.components["pump"].outlet.flow.mDot)
    cyc.computeCycle(1, 1e-5, 1e4, 5, forward=i%2==1)

for c in ["cooler", "orifice", "pump"]:
    print(cyc.components[c])