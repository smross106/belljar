import numpy as np

class SolidState(object):
  def __init__(self, mass, specificHeatCapacity, TRef=0):
    self.mass = mass
    self.specificHeatCapacity = specificHeatCapacity
    self.TRef = TRef

    self.T = TRef
    self.H = 0
    self.heatCapacity = self.mass * self.specificHeatCapacity

    self.accumulateE = 0
    self.previousdT = []
    self.previousAccumulateE = []

  def accumulateEnergy(self, E):
    self.accumulateE += E

  def accumulateTemperature(self, T):
    dE = self.heatCapacity * (T - self.T)
    self.accumulateEnergy(dE)

  def applyAccumulatedEnergy(self, dt):
    dE = self.accumulateE * dt
    dT = dE / self.heatCapacity
    self.H += dE

    self.previousdT.insert(0, dT*1.0)
    self.previousAccumulateE.insert(0, self.accumulateE*1.0)

    self.accumulateE = 0

class Station(object):
  def __init__(self):
    self.flow = None
    self.state = None

class SolidComponent(object):
  def __init__(self):
    pass
  def compute(self):
    pass

class SolidTwoStationComponent(object):
  def __init__(self):
    self.inlet = Station()
    self.outlet = Station()

  def compute(self):
    pass

class SolidConstantTemperatureBoundary(SolidComponent):
  def __init__(self):
    self.outlet = Station()

  def compute(self, T):
    self.outlet.state.accumulateTemperature(T)

class SolidThermalConductance(SolidTwoStationComponent):
  def __init__(self):
    self.inlet = Station()
    self.outlet = Station()

  def compute(self, thermalConductance):
    dT = self.outlet.state.T - self.inlet.state.T
    dE = dT * thermalConductance
    self.inlet.state.accumulateEnergy(dE)
    self.outlet.state.accumulateEnergy(-dE)



class CraneSolver(object):
    def __init__(self):
        self.components = {}
        self.states = []
        self.nStates = 0

        self.flags ={
            "INITL": 0,
            "MINUT": 3,
            "PREC": 2.5,
            "PRECN": 0,
            "NB": [],
            "DONT": False,
            "OKPC": False,
            "IS": 0,
            "KS": 0,
            "KT": 0,
            "KBC": 0,
            "dT":1,
            "E": np.zeros(10),
            "VNU": 0.14011612E-07,
            "BETA": 0,
            "BOUNDA": 0,
            "BOUNDB": 0,
            "BOUND": 0,
            "K": 0,
            "ENDONX": True,
            "XHAT": 10,
            "H1": 0,
            "BETA": 0,
            "XN": 0,
            "RX": 0,
            "AK": 0,
            "S": 0,
            "T": 0,
            "PMC": 0
        }
        self.stateValues = {
            "Y0": [],
            "Y1": [],
            "Y2": [],
            "Y3": [],
            "Y4": [],
            "Y5": [],
            "dY0": [],
            "dY1": [],
            "dY2": [],
            "dY3": [],
            "dY4": [],
            "dY5": []
        }

        self.Y = []
        self.dY = []
        self.Q = []

        self.X = 0  # equal to time
        self.H = .05  # equal to delta time

    def connectStations(self, station1, station2, mass, specificHeatCapacity, TRef):
        if station1.state != None:
            station2.state = station1.state
        elif station2.state != None:
            station1.state = station1.state
        else:
            st = SolidState(mass, specificHeatCapacity, TRef)
            station1.state = st
            station2.state = st

    def initialiseSolver(self):
        for comp_key in self.components.keys():
            comp = self.components[comp_key]
            try:
                if comp.inlet.state not in self.states:
                    self.states.append(comp.inlet.state)
            except:
                pass
            try:
                if comp.outlet.state not in self.states:
                    self.states.append(comp.outlet.state)
            except:
                pass

        for state in self.states:
            self.stateValues["Y0"].append(state.T)
            self.nStates += 1

        self.Y = np.zeros(self.nStates)
        self.dY = np.zeros(self.nStates)
        self.Q = np.zeros(self.nStates)

        for key in self.stateValues:
            self.stateValues[key] = np.zeros(self.nStates)
        
        self.flags["NB"] = np.zeros(self.nStates)
    
    def deriv(self):
        for i in range(self.nStates):
            self.dY[i] = 0.01 * np.sign(self.Y[i] - i)  * np.power(i - self.Y[i], 2)
    
    def b0(self):
        if self.flags["INITL"] < 0:
            self.b1()
        elif self.flags["INITL"] == 0:
           self.b2()
        else:
           self.b3()
    
    def b1(self):
       self.flags["INITL"] = 0
       self.flags["E"][self.flags["MINUT"]] = 10

       self.b2()
    
    def b2(self):
       # Original - check if the number of nodes is correct
       self.b5()
    
    def b5(self):
        print("B5")
        if self.flags["MINUT"] < 3:
            self.flags["MINUT"] = 3
        elif self.flags["MINUT"] > 10:
            self.flags["MINUT"] = 10
        
        if (self.flags["ENDONX"] and (self.flags["XHAT"]-self.X)/self.H < 0):
            self.b2040()
        
        self.flags["BETA"] = -np.log10(self.flags["VNU"]) * 0.9

        self.b70()

    def b70(self):
        if self.flags["PREC"] - self.flags["BETA"] < 0:
            self.b11()
        elif self.flags["PREC"] - self.flags["BETA"] == 0:
            self.b13()
        else:
            self.b12()
    
    def b12(self):
        self.flags["PREC"] = self.flags["BETA"] 
        self.b11()
    
    def b11(self):
        if self.flags["PREC"] < 2:
           self.flags["PREC"] = 4
        self.b13()
    
    def b13(self):
        self.flags["PRECN"] = self.flags["PREC"]*1.0
        self.flags["T"] = 10 **(-self.flags["PREC"])
        self.flags["BOUNDA"] = self.flags["T"] * 16.21966
        self.flags["BOUNDB"] = self.flags["T"] * 0.05
        self.flags["BOUND"] = self.flags["BOUNDA"] * self.flags["MINUT"]

        if self.flags["INITL"] <= 0:
           self.b1155()
        else:
           self.b16()
    
    def b1155(self):
        self.deriv()
        self.flags["KT"] = 18
        self.b17()
    
    def b17(self):
        self.flags["K"] = 0
        if (self.flags["ENDONX"] and (
           (self.flags["XHAT"] - self.X)/(4*self.H) < 1
        )):
            self.H = (self.flags["XHAT"]-self.X)/4
        
        self.flags["IS"] = 55
        self.flags["H1"] = 0.5 * self.H
        self.flags["KBC"] = 3
        self.flags["OKPC"] = False

        if self.flags["KT"] == 18:
            self.b18()
        elif self.flags["KT"] == 20:
            self.b20()
        else:
           raise ValueError("KT is wrong!")
    
    def b15(self):
        self.deriv()
        self.b18()

    def b18(self):
        self.flags["INITL"] += 1
        if self.flags["ENDONX"] and abs((self.flags["XHAT"] - self.X)/self.flags["XHAT"]) < self.flags["VNU"]:
           self.flags["ENDONX"] = False
        self.b19()
    
    def b19(self):
        if (self.flags["OKPC"]==1 or self.flags["INITL"]==1):
            #self.b999()
            self.b14()
            pass
        else:
            self.b14()
    
    def b999(self):
        raise ValueError("Block 999")
        #exit()

    def b14(self):
        self.b3()
    
    def b3(self):
        print("B3")
        if abs(self.flags["PREC"] - self.flags["PRECN"] > 0.1e-7):
           self.b70()
        else:
            self.b16()
    
    def b16(self):
        print("B16")
        if self.flags["DONT"] and self.flags["OKPC"]:
            self.b213()
        else:
            if self.flags["OKPC"]:
                self.b212()
            else:
                self.b20()
    
    def b20(self):
        print("B20")
        self.flags["K"] = 4
        self.b213()
    
    def b213(self):
        print("B213")
        # Includes b22
        for i in range(self.nStates):
            self.stateValues["dY3"][i] = self.stateValues["dY2"][i]
            self.stateValues["dY2"][i] = self.stateValues["dY1"][i]
            self.stateValues["dY1"][i] = self.stateValues["dY0"][i]
            self.stateValues["dY0"][i] = self.dY[i]

            self.stateValues["Y3"][i] = self.stateValues["Y2"][i]
            self.stateValues["Y2"][i] = self.stateValues["Y1"][i]
            self.stateValues["Y1"][i] = self.stateValues["Y0"][i]
            self.stateValues["Y0"][i] = self.Y[i]
        
        if self.flags["KBC"] == 0 or self.flags["DONT"]:
            self.b23()
        else:
            self.b24()
    
    def b24(self):
        self.flags["KBC"] -= 1
        self.flags["XN"] = self.X + self.H
        self.X = self.X + self.flags["H1"]
        self.flags["RX"] = 0.292893219
        self.flags["KS"] = 33

        self.b37()
    
    def b37(self):
        print("B37", self.H, self.flags["INITL"], self.flags["K"], self.X, self.flags["PREC"], self.flags["PRECN"])
        # Start of one Runge-Kutta step
        for i in range(self.nStates):
            #print(self.flags["K"], i)
            self.flags["AK"] = self.H * self.dY[i]
            if self.flags["KS"] == 33:
                self.b33(i)
            elif self.flags["KS"] == 34:
                self.b34(i)
            elif self.flags["KS"] == 30:
                self.b30(i)

        if self.flags["K"] == 1:
            self.b1551()
        elif self.flags["K"] == 2:
            self.b224()
        elif self.flags["K"] == 3:
            self.b219()
        elif self.flags["K"] == 4:
            self.b251()
        else:
            self.b219()
    
    def b33(self, i):
        self.Y[i] += 0.5*self.flags["AK"]
        self.Q[i] = self.flags["AK"]
        
        # self.b29()
        # Not needed as this just continues the loop of b37()
    
    def b30(self, i):
        self.Y[i] += (-self.Q[i] - self.Q[i] + self.flags["AK"])/6

        # self.b29()
        # Not needed as this just continues the loop of b37()
    
    def b34(self, i):
        R = self.flags["RX"] * (self.flags["AK"] - self.Q[i])
        self.Q[i] += ((3*R) - self.flags["RX"]*self.flags["AK"])
        self.Y[i] += R

        # self.b29()
        # Not needed as this just continues the loop of b37()
        
    def b219(self):
        self.flags["RX"] = 1.70710678
        self.b1551()

    def b251(self):
        self.flags["KS"] = 34
        self.b1551()
    
    def b224(self):
        self.flags["KS"] = 34
        self.X = self.flags["XN"]
        self.b1551()
    
    def b1551(self):
        self.deriv()
        self.b36()
    
    def b36(self):
        self.flags["K"] -= 1
        if self.flags["K"] <= 0:
            self.b18()
        else:
            self.b37()
    
    def b212(self):
        print("B212")
        if self.flags["ENDONX"] and ((self.flags["XHAT"]-self.X)/(4*self.H)) < 1:
            self.b17()
        elif self.flags["E"][self.flags["MINUT"]] > self.flags["BOUNDA"]:
            self.b28()
        else:
            self.b100()
    
    def b100(self):
        for i in range(1, self.flags["MINUT"]):
            self.flags["E"][0] +=  self.flags["E"][i]
        
        if self.flags["E"][0] < self.flags["BOUND"]:
            self.b38()
        else:
            self.b28()
    
    def b28(self):
        for i in range(self.nStates):
            self.stateValues["dY5"][i] = self.stateValues["dY4"][i]
            self.stateValues["dY4"][i] = self.stateValues["dY3"][i]
            self.stateValues["dY3"][i] = self.stateValues["dY2"][i]
            self.stateValues["dY2"][i] = self.stateValues["dY1"][i]
            self.stateValues["dY1"][i] = self.stateValues["dY0"][i]
            self.stateValues["dY0"][i] = self.dY[i]

            self.stateValues["Y5"][i] = self.stateValues["Y4"][i]
            self.stateValues["Y4"][i] = self.stateValues["Y3"][i]
            self.stateValues["Y3"][i] = self.stateValues["Y2"][i]
            self.stateValues["Y2"][i] = self.stateValues["Y1"][i]
            self.stateValues["Y1"][i] = self.stateValues["Y0"][i]
            self.stateValues["Y0"][i] = self.Y[i]

            # Block 39
            self.Y[i] = ((1.547625*self.Y[i]) + 
                         (2.017204*self.stateValues["Y2"][i]) + 
                         (-1.867503*self.stateValues["Y1"][i]) +
                         (0.697353*self.stateValues["Y3"][i]) + self.flags["S"]*(
                             (0.985508124*self.dY[i]) +
                             (0.895121303*self.stateValues["dY2"][i]) + 
                             (-0.351589071*self.stateValues["dY3"][i]) + 
                             (self.stateValues["dY1"][i]))
                         )
        
        self.b40()
    
    def b40(self):
        self.X += self.H
        self.deriv()

        self.b42()
    
    def b42(self):
        # Block 10
        for i in range(1, self.flags["MINUT"]):
            self.flags["E"][i-1] = self.flags["E"][i]
        self.flags["E"][self.flags["MINUT"]] = 0

        # Block 10, end of loop is Block 43
        for i in range(self.nStates):
            self.flags["T"] = self.Y[i] + self.flags["AK"]*(
                9*self.dY[i] + 
                19*self.stateValues["dY0"][i] + 
                self.stateValues["dY2"][i] - 
                5*self.stateValues["dY1"][i]
            )
            self.flags["PMC"] = self.Y[i] - self.flags["T"]

            if self.flags["DONT"]:
                # Go to b43
                pass
            else:
                self.b72(i)
                
            # block 43
            self.Y[i] = self.flags["T"]
        
        print("LOOP", self.X, self.H)
        print(self.Y, "\n")
        #print(self.X, self.H, self.dY)
        #print(self.flags["E"][self.flags["MINUT"]], self.flags["BOUNDB"], self.flags["AK"], self.flags["PMC"])
        #print(self.flags["S"])
        if self.X>=100 or self.H<1e-5:
            return

        if (self.flags["E"][self.flags["MINUT"]] - self.flags["BOUNDB"]) > 0:
            self.b50()
        else:
            self.b51()

    def b72(self, i):
        if self.flags["NB"][i] < 0:
            return
            # Return to block 43
        elif self.flags["NB"][i] == 0:
            self.b44(i)
        else:
            self.b47(i)
    
    def b44(self, i):
        if self.Y[i] == 0:
            return
            # Return to block 43
        else:
            self.b46(i)
    
    def b46(self, i):
        self.flags["PMC"] = self.flags["PMC"]/self.Y[i]

        self.b47(i)
    
    def b47(self, i):
        if abs(self.flags["PMC"]) > self.flags["E"][self.flags["MINUT"]]:
            self.b48(i)
        else:
            return
            # Return to block 43

    def b48(self, i):
        self.flags["E"][self.flags["MINUT"]] = abs(self.flags["PMC"])
        
    def b50(self):
        if self.flags["IS"] == 15:
            self.b15()
        elif self.flags["IS"] == 55:
            self.b55()
        else:
            raise ValueError("IS not good!")
    
    def b51(self):
        if self.flags["OKPC"] == False:
            self.b59()
        else:
            self.b52()
    
    def b52(self):
        print("#### B52")
        self.X -= self.H

        for i in range(self.nStates):
            self.Y[i] = self.stateValues["Y0"][i]
            self.dY[i] = self.stateValues["dY0"][i]

        # Block 54
        self.flags["KT"] = 20

        self.b53()
    
    def b53(self):
        if self.flags["DONT"]:
            self.b1000()
        else:
            self.H *= 0.625
            
            self.b17()
    
    def b55(self):
        print("B55")
        self.flags["OKPC"] = True
        self.flags["IS"] = 15
        self.b15()
    
    def b59(self):
        print("RESET RK4")
        self.X -= 4*self.H

        for i in range(self.nStates):
            self.Y[i] = self.stateValues["Y3"][i]
            self.dY[i] = self.stateValues["dY3"][i]
        
        self.b61()
    
    def b61(self):
        self.flags["KT"] = 18
        self.b53()
    
    def b38(self):
        self.H *= 2
        for i in range(self.nStates):
            self.stateValues["dY2"][i] = self.stateValues["dY3"][i]
            self.stateValues["dY3"][i] = self.stateValues["dY5"][i]
            self.stateValues["dY0"][i] = self.dY[i]

            self.stateValues["Y2"][i] = self.stateValues["Y3"][i]
            self.stateValues["Y3"][i] = self.stateValues["Y5"][i]
            self.stateValues["Y0"][i] = self.Y[i]  

        self.b23()          

    def b23(self):
        self.flags["AK"] = self.H/24
        self.flags["E"][self.flags["MINUT"]] = 10
        self.flags["S"] = self.H * 2.03169

        for i in range(self.nStates):
            self.Y[i] = ((1.547625*self.Y[i]) + 
                         (2.017204*self.stateValues["Y2"][i]) + 
                         (-1.867503*self.stateValues["Y1"][i]) +
                         (0.697353*self.stateValues["Y3"][i]) + self.flags["S"]*(
                             (0.985508124*self.dY[i]) +
                             (0.895121303*self.stateValues["dY2"][i]) + 
                             (-0.351589071*self.stateValues["dY3"][i]) + 
                             (self.stateValues["dY1"][i]))
                         )
        
        self.b40()




def tracefunc(frame, event, arg, indent=[0]):
      ignore = ["getstate", "decode", "__init__", "daemon", "<genexpr>", "call", "recurser", "extend", "_", "max", "<", "dictcomp", "err"]
      for i in ignore:
        if frame.f_code.co_name in i or i in frame.f_code.co_name:
            return
      if event == "call":
          indent[0] += 2
          print("-" * indent[0] + "> call function", frame.f_code.co_name)
      elif event == "return":
          print("<" + "-" * indent[0], "exit function", frame.f_code.co_name)
          indent[0] -= 2
      return tracefunc

import sys
#sys.setprofile(tracefunc)
sys.setrecursionlimit(5000)
        
cs = CraneSolver()

cs.components["s1"] = SolidConstantTemperatureBoundary()
cs.components["s2"] = SolidConstantTemperatureBoundary()
cs.components["s3"] = SolidConstantTemperatureBoundary()
cs.components["l1"] = SolidThermalConductance()
cs.components["l2"] = SolidThermalConductance()
cs.components["l3"] = SolidThermalConductance()
cs.components["l4"] = SolidThermalConductance()

cs.connectStations(cs.components["s1"].outlet, cs.components["l1"].inlet, 1, 100, 300)
cs.connectStations(cs.components["l1"].outlet, cs.components["l4"].inlet, 1, 100, 300)
cs.connectStations(cs.components["l4"].outlet, cs.components["s3"].outlet, 1, 100, 300)
cs.connectStations(cs.components["s2"].outlet, cs.components["l2"].inlet, 1, 100, 300)
cs.connectStations(cs.components["l2"].outlet, cs.components["l3"].inlet, 1, 100, 300)
cs.connectStations(cs.components["l1"].outlet, cs.components["l3"].outlet, 1, 100, 300)

cs.initialiseSolver()

cs.b0()

#print(cs.X)
#print(cs.H)
#print(cs.flags["E"][cs.flags["MINUT"]], cs.flags["BOUNDB"])
#print(cs.flags["PMC"])