from components2 import *
from scipy import optimize
import numpy as np


def makeCycle(fluid):

    cyc = CycleSolver(fluid)

    # Claude Cycle Nitrogen Example

    cyc.components["jn1"] = FlowJunction()
    cyc.components["comp"] = Compressor()
    cyc.components["cool"] = Cooler("T")
    cyc.components["hx1"] = Recuperator()

    cyc.components["sp1"] = FlowSplitter()
    cyc.components["ex1"] = Turbine()
    cyc.components["jn2"] = FlowJunction()

    cyc.components["hx2"] = Recuperator()
    cyc.components["hx3"] = Recuperator()

    cyc.components["th1"] = IdealThrottle()
    cyc.components["phase"] = PhaseSeperator()
    cyc.components["th2"] = IdealThrottle()

    cyc.components["sink"] = FluidSink()
    cyc.components["src"] = FluidSource_TP(11.28/60, 1, 288.15)

    cyc.connectStations("src", "outlet", "jn1", "inlet1")
    cyc.connectStations("jn1", "outlet", "comp", "inlet")
    cyc.connectStations("comp", "outlet", "cool", "inlet")
    cyc.connectStations("cool", "outlet", "hx1", "inlet1")
    cyc.connectStations("hx1", "outlet1", "sp1", "inlet")

    cyc.connectStations("sp1", "outlet1", "hx2", "inlet1")
    cyc.connectStations("hx2", "outlet1", "hx3", "inlet1")

    cyc.connectStations("hx3", "outlet1", "th1", "inlet")
    cyc.connectStations("th1", "outlet", "phase", "inlet")
    cyc.connectStations("phase", "liquidOutlet", "sink", "inlet")
    cyc.connectStations("phase", "gasOutlet", "th2", "inlet")
    cyc.connectStations("th2", "outlet", "hx3", "inlet2")
    cyc.connectStations("hx3", "outlet2", "jn2", "inlet1")
    cyc.connectStations("jn2", "outlet", "hx2", "inlet2")
    cyc.connectStations("hx2", "outlet2", "hx1", "inlet2")
    cyc.connectStations("hx1", "outlet2", "jn1", "inlet2")

    cyc.connectStations("sp1", "outlet2", "ex1", "inlet")
    cyc.connectStations("ex1", "outlet", "jn2", "inlet2")

    return(cyc)

def computeCycle(self, pHigh, pLiq, pIn, expanderFrac, TCool, HX3eff):

    PR_comp = pHigh/pIn
    log_PR = np.log(PR_comp) / np.log(10)
    eta_poly = 0.7578 - 0.1954 * log_PR
    compr_eff = ((np.power(PR_comp, 0.4/1.4) - 1) / 
                (np.power(PR_comp, 0.4/(1.4*eta_poly)) - 1))
    #compr_eff = 0.8
    expand_eff = 1 - (1 - compr_eff)*0.8
    
    self.components["comp"].compute(pOut=pHigh, eta=compr_eff)
    self.components["cool"].compute(TOut=TCool)

    self.components["hx1"].compute(eta=0.99)
    self.components["hx1"].computeStream1()

    self.components["sp1"].compute(flowFraction2=expanderFrac)

    self.components["hx2"].compute(eta=1.0)
    self.components["hx2"].computeStream1()

    self.components["hx3"].compute(eta=HX3eff)
    self.components["hx3"].computeStream1()

    self.components["th1"].compute(pOut=pLiq)
    self.components["phase"].compute()
    self.components["th2"].compute(pOut=pIn)

    self.components["hx3"].computeStream2()
    self.components["ex1"].compute(pOut=pIn, eta=expand_eff)
    self.components["jn2"].compute()
    self.components["hx2"].computeStream2()
    self.components["hx1"].computeStream2()

    self.components["jn1"].compute()

pIn = 5.0
pHigh = pIn*30
pLiq = pIn + 0.1
TCool = 273+30
expanderFrac=0.5
TIn = 273
mDotIn = 1


def ESM_radiator(T_reject):
    # kg/W
    delta_T = T_reject - 235
    if delta_T < 0:
        return(np.nan)
    one_over_ESM = (5.858E-6 * np.power(delta_T, 2)) - (1.004E-4 * delta_T) + (1.498E-2)
    ESM = 1/one_over_ESM

    return(ESM/1000)

def ESM_FCHX(T_reject):
    # kg/W
    delta_T = T_reject - 235
    if delta_T < 0:
        return(np.nan)
    one_over_ESM = 2.853E-2 * np.log(delta_T) + 7.468E-2
    #one_over_ESM = (1.7225E-8 * np.power(T_reject, 3)) - (2.141E-5 * np.power(T_reject, 2)) + (8.863E-3 * T_reject) - 1.0003E0
    ESM = 1/one_over_ESM

    return(ESM/1000)

def compressor_mass(power_kw):
    mass = 31.83 * np.power(power_kw, 0.476)

    return(mass)

def getCyclePerformance(input_args, *args, output=False):
    (compPressureRatio, liquidPressureRatio, TCool, expanderFrac, HX3eff) = input_args
    (pIn, TIn, mDotIn, ESM_vector, fluid) = args

    pHigh = pIn * compPressureRatio
    pLiq = pIn * liquidPressureRatio

    cyc = makeCycle(fluid)

    cyc.fluid = fluid
    cyc.computeCycle = computeCycle.__get__(cyc)
    cyc.initialiseSolver(273, pIn, mDotIn*2)
    cyc.components["src"].p = pIn
    cyc.components["src"].T = TIn
    cyc.components["src"].mDot = mDotIn
    cyc.components["src"].compute()


    prevQ = -1
    nWithGoodConv = 0
    for i in range(50):
        if cyc.components["th1"].outlet.state.Q > 0:
            if abs(cyc.components["th1"].outlet.state.Q / prevQ - 1) < 0.001:
                nWithGoodConv += 1
        prevQ = cyc.components["th1"].outlet.state.Q
        if nWithGoodConv > 5:
            break
        #print(i, "\t", cyc.components["th1"].outlet.state.Q, nWithGoodConv)
        #for c in cyc.components:
        #    print(cyc.components[c])

        cyc.computeCycle(pHigh, pLiq, pIn, expanderFrac, TCool, HX3eff)

        
    
    if nWithGoodConv < 3 or cyc.components["phase"].liquidOutlet.flow.mDot<=0:
        return(1e10)

    liqEnth = PropsSI("Hmass", "T", TIn, "P", pIn*1e5, fluid) - PropsSI("Hmass", "Q", 0, "P", pIn*1e5, fluid)
    W = cyc.components["comp"].WIn + cyc.components["ex1"].WIn
    Q = -cyc.components["cool"].QIn
    qHX1 = abs(cyc.components["hx1"].QDot)
    qHX2 = abs(cyc.components["hx2"].QDot)
    qHX3 = abs(cyc.components["hx3"].QDot)

    massHX = 2.72/1000 * (qHX1 + qHX2 + qHX3 + Q)
    massComp = compressor_mass(cyc.components["comp"].WIn/1000)
    massExpa = compressor_mass(-cyc.components["ex1"].WIn/1000)
    massW = ESM_vector[1] * W 
    massQ = ESM_vector[2](TCool) * Q

    spW = W/cyc.components["phase"].liquidOutlet.flow.mDot

    print("pHigh {:.1f}bar, exp frac {:.2f} eta {:.1f}% M {:.0f}kg".format(pHigh, expanderFrac, 100*liqEnth/spW, massComp+massExpa+massHX+massQ+massW))
    #print(pHigh, expanderFrac,"eta",liqEnth/spW)
    if output:
        print(pHigh, pLiq, pIn, expanderFrac, TCool, HX3eff)
        print(spW, liqEnth)
        print(W, Q)
        print(cyc.components["phase"].liquidOutlet.flow.mDot)
        print((massComp+massExpa+massHX+massQ+massW))

    return((massComp+massExpa+massHX+massQ+massW) / cyc.components["phase"].liquidOutlet.flow.mDot)

ESM_vec = [1, 150/1000, ESM_radiator]
otherargs = (5, 273, 1/60, ESM_vec, "Oxygen")

x0 = np.array([30, 1.01, 273+30, 0.5, 0.5])
x0_upper = np.array([100, 3, 373, 1, 1])
N = len(x0)
zdelt = 0.00025
sim = np.empty((N + 1, N), dtype=x0.dtype)
sim[0] = x0
for k in range(N):
    y = 0.8*np.array(x0, copy=True)
    if y[k] != 0:
        y[k] = x0_upper[k]
        y[k] = zdelt
    sim[k + 1] = y

properties = optimize.minimize(getCyclePerformance, 
        x0, 
        args=otherargs,
        bounds=[(1.1, 100), (1.01, 3), (262, 373), (0, 1), (0, 1)],
        method="Nelder-Mead", options={"maxiter": 100, "disp":False, "adaptive":False, "return_all":False, "xatol":1e-4, "initial_simplex":sim}
        #method='SLSQP'
        #constraints=object_cons, method="trust-constr", jac="3-point", hess=optimize.BFGS(), options={"maxiter":2000}
        )
xOpt = properties.x
print(xOpt)

xOptMethane5bar273K = [2.17122696, 1.01000000, 263.967033, 0.931878221, 4.14349595e-04]
xOptMethane5bar150K = [2.14897084,   1.01,       262.,           0.9324098,    0.        ]
print("initial")
getCyclePerformance(x0, otherargs[0], otherargs[1], otherargs[2], otherargs[3], otherargs[4], output=True)
print("optimal")
getCyclePerformance(xOpt, otherargs[0], otherargs[1], otherargs[2], otherargs[3], otherargs[4], output=True)

#print(getCyclePerformance((30, 1.01, 273+30, 0.5, 0.5), 5, 273, 1/60, ESM_vec, "Methane"))

"""
print("#")
cyc.printStation(1, "src", "outlet")
cyc.printStation(2, "comp", "inlet")
cyc.printStation(3, "comp", "outlet")
cyc.printStation(4, "cool", "outlet")
cyc.printStation(5, "hx1", "outlet1")

cyc.printStation(6, "hx2", "inlet1")
cyc.printStation(7, "hx3", "inlet1")
cyc.printStation(8, "th1", "inlet")
cyc.printStation(9, "th1", "outlet")

cyc.printStation(10, "th2", "inlet")
cyc.printStation(11, "th2", "outlet")

cyc.printStation(12, "hx3", "outlet2")
cyc.printStation(13, "hx2", "inlet2")

cyc.printStation(14, "hx1", "inlet2")
cyc.printStation(15, "hx1", "outlet2")

cyc.printStation(16, "ex1", "inlet")
cyc.printStation(17, "ex1", "outlet")

cyc.printStation(18, "phase", "liquidOutlet")"""


exit()
liqEnth = PropsSI("Hmass", "T", 273+15, "P", 1.0e5, "Nitrogen") - PropsSI("Hmass", "Q", 0, "P", 1.0e5, "Nitrogen")
dB = -cyc.components["src"].outlet.state.b + cyc.components["sink"].inlet.state.b
dH = cyc.components["comp"].outlet.state.h - cyc.components["comp"].inlet.state.h
dE = cyc.components["comp"].WIn / cyc.components["comp"].inlet.flow.mDot + cyc.components["ex1"].WIn/cyc.components["ex1"].inlet.flow.mDot

W = cyc.components["comp"].WIn + cyc.components["ex1"].WIn
print("W", cyc.components["comp"].WIn, cyc.components["ex1"].WIn)
spW = W/cyc.components["phase"].liquidOutlet.flow.mDot
print(W)
print(spW)
print(liqEnth)
print(liqEnth/spW)
print(dB)
print(dB/spW)