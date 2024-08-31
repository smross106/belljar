import numpy as np
class SolidState(object):
    def __init__(self, mass=1, specificHeatCapacity=1, TRef=298):
        self._T = TRef
        self.mass = mass
        self.specificHeatCapacity = specificHeatCapacity
        self.heatCapacity = self.mass * self.specificHeatCapacity

        self.boundaryConditionT = False

        self.storedE = 0
        self.dTdt = 0
        self.THistory = []
        self.dTHistory = []

        self.intermediateTHistory = []
        self.intermediatedTdtHistory = []
        self.Tnp1 = 0
        self.Tnm1_temp = 0
    
    @property
    def T(self):
        return self._T
    
    @T.setter
    def T(self, value):
        if not self.boundaryConditionT:
            self._T = value
        else:
            self._T = self.boundaryConditionT

        
    def computedTdt(self, time):
        self.dTdt = self.storedE/(self.heatCapacity)
        self.storedE = 0

    def setNextT(self, dt):
        self.THistory.append(self.Tnm1_temp)
        self.dTHistory.append((self.Tnp1 - self.Tnm1_temp)/dt)
        if not self.boundaryConditionT:
            self._T = self.Tnp1 * 1.0
        self.Tnp1 = 0
        self.intermediateTHistory = []
        self.intermediatedTdtHistory = []
    
    def resetSolver(self):
        self.T = self.Tnm1_temp
        self.Tnp1 = 0
        self.intermediateTHistory = []
        self.intermediatedTdtHistory = []
    
    def stateInitialise(self, TRef):
        self.THistory = []
        self.intermediateTHistory = []
        self.intermediatedTdtHistory = []
        self.Tnm1_temp = 0
        self.Tnp1 = 0 
        self._T = TRef
        self.storedE = 0

class SolidStation(object):
    def __init__(self):
        self.state = None

        self.nVarHistory = 0
        self.nDerivHistory = 0
    
    def computeDerivatives(self, time):
        self.state.computedTdt(time)
        self.state.intermediatedTdtHistory.append(self.state.dTdt)
    
    def backupVariables(self):
        self.state.Tnm1_temp = self.state.T * 1.0
    
    def updateVariables(self, variableHistoryValues={}, intermediateDerivativeHistoryValues={}, 
                        derivativeHistoryValues={},
                        updatePlusOne=False, returnValue=False, prevValue=1):
        nextT = self.state.Tnm1_temp * prevValue

        for key in variableHistoryValues:
            nextT += self.state.THistory[key] * variableHistoryValues[key]
        for key in intermediateDerivativeHistoryValues:
            nextT += self.state.intermediatedTdtHistory[key] * intermediateDerivativeHistoryValues[key]
        for key in derivativeHistoryValues:
            nextT += self.state.dTHistory[key] * derivativeHistoryValues[key]
        
        if returnValue:
            return([nextT])
        else:
            if updatePlusOne:
                self.state.Tnp1 = nextT * 1.0
            else:
                self.state.T = nextT * 1.0
    
    def setNextVar(self, dt):
        self.state.setNextT(dt)

        self.nVarHistory += 1
        self.nDerivHistory += 1
    
    def resetNextVar(self):
        self.state.resetSolver()
    
    def calculateErrorValueInput(self, correctedValues):
        errors = []
        errors.append(abs(self.state.Tnp1 - correctedValues[0])/self.state.Tnm1_temp)
        return(errors)

    def calculateErrorPreviousValues(self):
        errors = []
        errors.append(abs(self.state.Tnp1 - self.state.Tnm1_temp)/self.state.Tnm1_temp)
        return(errors)

class SolidComponent(object):
    def __init__(self):
        self.linkedStations = []
    
    def compute(self):
        pass

class ThermalConductance(SolidComponent):
    def __init__(self, C=1):
        super().__init__()
        self.C = C
        # Thermal conductance, W/K
        
        self.inlet = SolidStation()
        self.outlet = SolidStation()

        self.linkedStations = [self.inlet, self.outlet]
    
    def compute(self, time):
        dT = self.inlet.state.T - self.outlet.state.T

        self.inlet.state.storedE -= dT * self.C
        self.outlet.state.storedE += dT * self.C

class ConstantTemperatureBoundary(SolidComponent):
    def __init__(self, TSet):
        super().__init__()

        self.TSet = TSet

        self.outlet = SolidStation()

        self.linkedStations = [self.outlet]
    
    def compute(self, time):
        if callable(self.TSet):
            TSet = self.TSet(time)
        else:
            TSet = self.TSet * 1.0
        self.outlet.state.boundaryConditionT = TSet
        dT = TSet - self.outlet.state.T

        self.outlet.state.storedE += dT * self.outlet.state.heatCapacity

class HeatInputBoundary(SolidComponent):
    def __init__(self, QIn):
        super().__init__()

        self.QIn = QIn

        self.outlet = SolidStation()

        self.linkedStations = [self.outlet]
    
    def compute(self, time):
        if callable(self.QIn):
            QIn = self.QIn(time)
        else:
            QIn = self.QIn * 1.0

        self.outlet.state.storedE += QIn

class SolverSolid(object):
    def __init__(self, precision=4):
        self.components = {}
        self.stations = []
        self.stationsUnique = []
        self.stationsDuplicates = []
        self.states = []

        self.precision = precision
        self.errors = []
        self.time = []
        self.timesteps = []
        self.allowTimeAdjust = True

        self.mode = "RK4ERROR"

        self.controlParameters = {}

    def connectStationsSolid(self, component1, station1, component2, station2, mass=1, specificHeatCapacity=1, TRef=298):
            
        fs = SolidState(mass=mass, specificHeatCapacity=specificHeatCapacity, TRef=TRef)

        setattr(getattr(self.components[component1], station1), "state", fs)
        setattr(getattr(self.components[component2], station2), "state", fs)

    
    def initialiseSolver(self, Tref, mDot=1, p_bar=1):
        for component_key in self.components.keys():
            for station in self.components[component_key].linkedStations:

                if type(station.state) == SolidState:
                    station.state.stateInitialise(Tref)
                else:
                    station.state.stateInitialise(Tref, mDot, p_bar)

                if station.state not in self.states:
                    self.stationsUnique.append(station)
                    if type(station) == SolidStation:
                        self.stations.append([station, station.state, None])
                    else:
                        self.stations.append([station, station.state, None])
                    self.states.append(station.state)

        self.errors = []
        self.time = []
    
    def calculateDerivatives(self, time):
        for component_key in self.components.keys():
            self.components[component_key].compute(time)
        
        for station in self.stations:
            station[0].computeDerivatives(time)
    
    def changeTimestep(self, dt, errorFailed):
        if self.allowTimeAdjust:
            if errorFailed:
                return(dt * 0.625)
            else:
                return(dt * 2)
        else:
            return(dt)

    def solverStepSuccess(self, time, dt, error):
        for station in self.stations:
            station[0].setNextVar(dt)
        time += dt
        self.time.append(time)
        self.errors.append(error)
        self.timesteps.append(dt)

        self.systemController(time, dt)

        return(time)

    def runSolver(self, t0, t1, dtTrial):
        time = t0
        dt = dtTrial
        minRunsBetweenTimestepIncrease = 3
        runsSinceTimestepIncrease = 0
        numIters = 0
        runsSinceLastError = 0

        milestones = np.linspace(t0, t1, 25)

        while time < t1:
            #print(time, dt)
            if self.mode == "RK4ERROR":
                success, error = self.RK4OneStepError(time, dt)
                if success:
                    time = self.solverStepSuccess(time, dt, error)

                    runsSinceTimestepIncrease += 1
                    if runsSinceTimestepIncrease >= minRunsBetweenTimestepIncrease:
                        dt = self.changeTimestep(dt, False)
                        runsSinceTimestepIncrease = 0
                else:
                    for station in self.stations:
                        station[0].resetNextVar()
                    dt = self.changeTimestep(dt, True)

            elif self.mode == "RK4GILL":

                success, error = self.RK4OneStepGill(time, dt)
                time = self.solverStepSuccess(time, dt, np.nan)

            elif self.mode == "CRANE":
                pass
                #success, error = self.PredictorCorrectorCraneError(time, dt)

            elif self.mode == "CRANERK4":
                # Crane algorithm, falling back to RK4 with the previous timestep
                if runsSinceLastError < 3:
                    # We don't have enough previous dT/dt values to use the Crane method - fall back to RK4
                    success, error = self.RK4OneStepGill(time, dt)
                    if success:
                        runsSinceLastError += 1
                        time = self.solverStepSuccess(time, dt, error)

                        runsSinceTimestepIncrease += 1
                else:
                    # We can happily use the Crane algorithm
                    success, error = self.PredictorCorrectorCraneError(time, dt)
                    if success:
                        runsSinceLastError += 1
                        time = self.solverStepSuccess(time, dt, error)

                        runsSinceTimestepIncrease += 1

                        if runsSinceTimestepIncrease >= minRunsBetweenTimestepIncrease:
                            dt = self.changeTimestep(dt, False)
                            runsSinceTimestepIncrease = 0
                    else:
                        for station in self.stations:
                            station[0].resetNextVar()
                        dt = self.changeTimestep(dt, True)

            else:
                break

            if time > milestones[0]:
                print("Progress: {:.0f}%".format(100*time/(t1-t0)))
                milestones = np.delete(milestones, 0)
            numIters += 1
            if numIters > 100000:
                break
    
    def RK4OneStepError(self, time, dt):
        # Back up previous temperatures to allow for RK4 to occur
        for station in self.stations:
            station[0].backupVariables()
        
        
        # Step 0 - generates k0
        self.calculateDerivatives(time)
        for station in self.stations:
            # Feed in T_k1
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0: dt/2
                })

        # Step 1 - generates k1
        self.calculateDerivatives(time+dt/2)
        for station in self.stations:
            # Feed in T_k2
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0: dt/4,
                    1: dt/4
                })

        # Step 2 - generates k2
        self.calculateDerivatives(time+dt/2)
        for station in self.stations:
            # Feed in T_k3
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    1: 2*dt,
                    2:  -dt
                })

        # Step 3 - generates k3
        self.calculateDerivatives(time+dt)
        for station in self.stations:
            # Feed in T_k4
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0:  7*dt/27,
                    1: 10*dt/27,
                    3:  1*dt/27
                })

        # Step 4 - generates k4
        self.calculateDerivatives(time+2*dt/3)
        for station in self.stations:
            # Feed in T_k5
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0:   28*dt/625,
                    1: -125*dt/625,
                    2:  546*dt/625,
                    3:   54*dt/625,
                    4: -378*dt/625
                })

        maxError = 0
        # Step 5 - generates k5
        self.calculateDerivatives(time+dt/5)
        for station in self.stations:

            # Calculate error
            station[0].updateVariables(
                variableHistoryValues = {}, updatePlusOne=True,
                intermediateDerivativeHistoryValues={
                    0:   dt/6,
                    2: 4*dt/6,
                    3:   dt/6,
                })

            var_np1_bar = station[0].updateVariables(
                variableHistoryValues = {}, updatePlusOne=False, returnValue=True,
                intermediateDerivativeHistoryValues={
                    0:  14*dt/336,
                    3:  35*dt/336,
                    4: 162*dt/336,
                    5: 125*dt/336
                })
            
            maxError = max(max(station[0].calculateErrorValueInput(var_np1_bar)), maxError)

        if -np.log10(maxError) > self.precision:
            return(True, maxError)
        else:
            return(False, maxError)

    def RK4OneStepGill(self, time, dt):
        # Back up previous temperatures to allow for RK4 to occur
        for station in self.stations:
            station[0].backupVariables()
        
        # Step 0 - generates k0
        self.calculateDerivatives(time)
        for station in self.stations:
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0: dt/2
                })

        # Step 1 - generates k1
        self.calculateDerivatives(time+dt/2)
        for station in self.stations:
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    0: (-1 + np.sqrt(2)) * dt/2,
                    1: ( 2 - np.sqrt(2)) * dt/2
                })
        
        # Step 2 - generates k2
        self.calculateDerivatives(time+dt/2)
        for station in self.stations:
            station[0].updateVariables(
                variableHistoryValues = {}, intermediateDerivativeHistoryValues={
                    1: (   -np.sqrt(2)) * dt/2,
                    2: (1 + np.sqrt(2)) * dt/2
                })
        
        # Step 3 - generates k3
        self.calculateDerivatives(time+dt)
        for station in self.stations:
            station[0].updateVariables(
                variableHistoryValues = {}, updatePlusOne=True, 
                intermediateDerivativeHistoryValues={
                    0:                    dt/6,
                    1: (2 - np.sqrt(2)) * dt/6,
                    2: (2 + np.sqrt(2)) * dt/6,
                    3:                    dt/6,
                })
        return(True, 0)
        
    def PredictorCorrectorCraneError(self, time, dt):
        # Back up previous temperatures to allow for RK4 to occur
        for station in self.stations:
            if station[0].nVarHistory < 3 or station[0].nDerivHistory < 3:
                return(False, 1000)
            station[0].backupVariables()
        
        # Carry out predictor step
        self.calculateDerivatives(time)
        for station in self.stations:
            # This is y'_n

            station[0].updateVariables(
                prevValue = 1.54765200,
                variableHistoryValues = {
                    -1: -1.86750300,
                    -2:  2.01720400,
                    -3: -0.69735300
                }, 
                intermediateDerivativeHistoryValues={
                    -1:  2.00224700 * dt
                },
                derivativeHistoryValues={
                    -1: -2.03169000 * dt,
                    -2:  1.81860900 * dt,
                    -3: -0.71432000 * dt 
                })
        
        self.calculateDerivatives(time+dt)
        maxError = 0
        for station in self.stations:
            # This is p'_n+1

            station[0].updateVariables(
                variableHistoryValues = {}, updatePlusOne=True, returnValue=False, 
                intermediateDerivativeHistoryValues={
                    -1:  0.375000000 * dt,
                    -2:  0.791666667 * dt
                },
                derivativeHistoryValues={
                    -1: -0.208333333 * dt,
                    -2:  0.041666667 * dt,
               })
            
            maxError = max(max(station[0].calculateErrorPreviousValues()), maxError)

        #print(time, dt, maxError)
        if -np.log10(maxError) > self.precision:
            return(True, maxError)
        else:
            return(False, maxError)

    def systemController(self, time, dt):
        pass


if False:

    cs = SolverSolid()


    def hat(time):
        return(200+10*np.sin(time/5))

    cs.components["s1"] = HeatInputBoundary(500)
    cs.components["s2"] = ConstantTemperatureBoundary(hat)
    cs.components["s3"] = ConstantTemperatureBoundary(300)
    cs.components["l1"] = ThermalConductance(C=10)
    cs.components["l2"] = ThermalConductance(C=20)
    cs.components["l3"] = ThermalConductance(C=10)
    cs.components["l4"] = ThermalConductance(C=10)


    cs.connectStationsSolid("s1", "outlet", "l1", "inlet", specificHeatCapacity=10)
    cs.connectStationsSolid("l1", "outlet", "l4", "inlet", specificHeatCapacity=10)
    cs.connectStationsSolid("l4", "outlet", "s3", "outlet", specificHeatCapacity=10)
    cs.connectStationsSolid("s2", "outlet", "l2", "inlet", specificHeatCapacity=10)
    cs.connectStationsSolid("l2", "outlet", "l3", "inlet", specificHeatCapacity=10)
    cs.connectStationsSolid("l1", "outlet", "l3", "outlet", specificHeatCapacity=10)

    #cs.connectStations(cs.components["s1"].outlet, cs.components["l1"].inlet, 1, 100, 300)
    #cs.connectStations(cs.components["l1"].outlet, cs.components["l4"].inlet, 1, 100, 300)
    #cs.connectStations(cs.components["l4"].outlet, cs.components["s3"].outlet, 1, 100, 300)
    #cs.connectStations(cs.components["s2"].outlet, cs.components["l2"].inlet, 1, 100, 300)
    #cs.connectStations(cs.components["l2"].outlet, cs.components["l3"].inlet, 1, 100, 300)
    #cs.connectStations(cs.components["l1"].outlet, cs.components["l3"].outlet, 1, 100, 300)

    cs.initialiseSolver(298)
    cs.precision = 4



    print("crane")
    cs.initialiseSolver(298)
    cs.mode="CRANERK4"

    cs.runSolver(0, 1, 0.01)


    import matplotlib.pyplot as plt

    plt.plot(cs.time, cs.components["l1"].outlet.state.THistory, label="Outlet")
    plt.plot(cs.time, cs.components["l2"].inlet.state.THistory, label="BC")
    plt.plot(cs.time, cs.components["l2"].outlet.state.THistory, label="Inter")
    plt.legend()

    print("#")
    print(cs.components["l2"].inlet.state.T)
    print(cs.components["l3"].inlet.state.T)
    print(cs.components["l3"].outlet.state.T)
    print(cs.components["l4"].outlet.state.T)
    print(cs.components["l1"].inlet.state.T)


    plt.show()