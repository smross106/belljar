import belljar_lockstep as bjls
import numpy as np
from copy import copy
from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI
import time

zones = ["Head", "Trunk", "Arm L", "Arm R", "Hand L", "Hand R", "Leg L", "Leg R", "Foot L", "Foot R"]
thermalLayers = ["Core", "Muscle", "Fat", "Skin"]
calculationLayers = ["Core skeleton", "Core viscera", "Muscle", "Fat", "Skin"]


valueDict = {
        "Head":   0, "Trunk":  0,
        "Arm L":  0, "Arm R":  0,
        "Hand L": 0, "Hand R": 0,
        "Leg L":  0, "Leg R":  0,
        "Foot L": 0, "Foot R": 0
    }



def makeSolverBodyParams(height_cm, weight_kg, varyMetabolism=True, useNavyBodyFat=True, 
                         male=True, waistCircum=50, neckCircum=20, hipCircum=60):
    
    # Du Bois formula
    totalSurfaceArea = 0.007184 * np.power(height_cm, 0.725) * np.power(weight_kg, 0.45)

    zoneAreas = {
        "Head": 0.070 , "Trunk": 0.3602,
        "Arm L": 0.06075, "Arm R": 0.06075,
        "Hand L": 0.025, "Hand R": 0.025,
        "Leg L": 0.1587, "Leg R": 0.1587,
        "Foot L": 0.0343, "Foot R": 0.0343
    }
    for zone in zoneAreas:
        zoneAreas[zone] *= totalSurfaceArea


    if useNavyBodyFat:
        # US Navy body fat estimation
        # https://www.omnicalculator.com/health/navy-body-fat#what-formula-does-the-us-navy-body-fat-calculator-use
        if male:
            fatWeight = 495 / (1.0324
                                - 0.19077 * np.log10(waistCircum - neckCircum)
                                + 0.15456 * np.log10(height_cm)) - 450
        else:
            fatWeight = 495 / (1.29579 
                                - 0.35004 * np.log10(waistCircum + hipCircum - neckCircum)
                                + 0.22100 * np.log10(height_cm)) - 450
        fatWeight *= weight_kg/100
    else:
        # Original METMAN formula
        fatWeight = (0.8 * np.power(height_cm, 0.242)/np.power(weight_kg*1000, 0.1)) + 0.162
        fatWeight = max((5.548/fatWeight) - 5.044, 0)
        fatWeight *= weight_kg

    leanWeight = weight_kg - fatWeight

    calculationNodeValues = {}
    thermalNodeValues = {}
    parameters = {"Weight":0, "Length":0, "Radius":0 , "Radius outer": 0, "Volume":0, 
                  "Capacitance":0, "Basal metabolism":0, "Basal bloodflow":0,
                  "Conductance out":0, "Conductance out to":""}
    for zone in zones:
        for layer in calculationLayers:
            calculationNodeValues[zone+"/"+layer] = copy(parameters)
        for layer in thermalLayers:
            thermalNodeValues[zone+"/"+layer] = copy(parameters)
    thermalNodeValues["Blood"] = copy(parameters)
    
    weightDistribution = {"Core skeleton":{
            "Head":   0.0192920, "Trunk":  0.0470400,
            "Arm L":  0.0118780, "Arm R":  0.0118780,
            "Hand L": 0.0018240, "Hand R": 0.0018240,
            "Leg L":  0.0395840, "Leg R":  0.0395840,
            "Foot L": 0.0029242, "Foot R": 0.0029242
        },
        "Core viscera":{
            "Head":   0.0282320, "Trunk":  0.1870400,
            "Arm L":  0.0058820, "Arm R":  0.0058820,
            "Hand L": 0.0002352, "Hand R": 0.0002352,
            "Leg L":  0.0151760, "Leg R":  0.0151760,
            "Foot L": 0.0004706, "Foot R": 0.0004706
        },
        "Muscle":{
            "Head":   0.005880, "Trunk":  0.283400,
            "Arm L":  0.026640, "Arm R":  0.026640,
            "Hand L": 0.000594, "Hand R": 0.000594,
            "Leg L":  0.080500, "Leg R":  0.080500,
            "Foot L": 0.000594, "Foot R": 0.000594
        },
        "Fat": {
            "Head":   0.033300, "Trunk":  0.633300,
            "Arm L":  0.043350, "Arm R":  0.043350,
            "Hand L": 0.006665, "Hand R": 0.006665,
            "Leg L":  0.106650, "Leg R":  0.106650,
            "Foot L": 0.010000, "Foot R": 0.010000
        },
        "Skin":{
            "Head":   0.00423, "Trunk":  0.02130,
            "Arm L":  0.00382, "Arm R":  0.00382,
            "Hand L": 0.00147, "Hand R": 0.00147,
            "Leg L":  0.00947, "Leg R":  0.00947,
            "Foot L": 0.00188, "Foot R": 0.00188
        }
    }
    # Specific heat capacities in J/kg
    heatCapacity = {
        "Skeleton": 2092,
        "Viscera": 3765.6,
        "Muscle": 3765.6,
        "Fat": 2510.4,
        "Skin": 3765.6
    }
    zoneThermalConductivity = {
        "Core":{
        "Head":   0.80749, "Trunk":  0.08041 ,
        "Arm L":  0.22854, "Arm R":   0.22854,
        "Hand L": 0.20284, "Hand R": 0.20284 ,
        "Leg L":  0.18255, "Leg R":  0.18255,
        "Foot L": 0.15586, "Foot R": 0.15586
    },
    "Muscle":{
        "Head":   0.03771, "Trunk":  0.06950,
        "Arm L":  0.15078, "Arm R":  0.15078,
        "Hand L": 0.04367, "Hand R": 0.04367 ,
        "Leg L":  0.21164, "Leg R":  0.21164,
        "Foot L": 0.05060, "Foot R": 0.05060
    },
    "Fat":{
        "Head":   0.06480, "Trunk":  0.03874,
        "Arm L":  0.08103, "Arm R":  0.08103,
        "Hand L": 0.07041, "Hand R": 0.07041,
        "Leg L":  0.05181, "Leg R":  0.05181,
        "Foot L": 0.06403, "Foot R": 0.06403
    },
    "Skin":{
        "Head":   0.06480, "Trunk":  0.03874,
        "Arm L":  0.08103, "Arm R":  0.08103,
        "Hand L": 0.07041, "Hand R": 0.07041,
        "Leg L":  0.05181, "Leg R":  0.05181,
        "Foot L": 0.06403, "Foot R": 0.06403
    }}
    # Thermal conductivities of each layer and zone in W/(m-K)
    # Unit converted from Metman manual data, appendix D.2 table D.1 

    layerWeight = {
        "Core":   leanWeight * 0.4354,
        "Muscle": leanWeight * 0.5058,
        "Fat":    fatWeight,
        "Skin":   leanWeight * 0.0588
    }

    # Conversion of per-area scaling from kcal/hr-m2 to W/m2, and kcal/hr-kg to W/kg
    metabolismAreaFactor = 45.02    # W/m2
    metabolismWeightFactor = 0.3487 # W/kg
    basalMetabolism = {
        "Core":   metabolismAreaFactor * totalSurfaceArea * 0.10,
        "Muscle": metabolismAreaFactor * totalSurfaceArea * 0.18,
        "Fat":    metabolismWeightFactor * layerWeight["Fat"],
        "Skin":   metabolismWeightFactor * layerWeight["Skin"]
    }

    # Skin bloodflow all converted to kg/s/kg-skin (1.060kg/L). Original values in L/hr/kg-skin
    skinBloodflowFactor = {
        "Head":   5.34, "Trunk":  1.56,
        "Arm L":  1.35, "Arm R":  1.35,
        "Hand L": 5.58, "Hand R": 5.58,
        "Leg L":  0.84, "Leg R":  0.84,
        "Foot L": 2.34, "Foot R": 2.34
    }
    for zone in skinBloodflowFactor:
        skinBloodflowFactor[zone] *= 1.060 / 3600

    for zone in zones:
        # Weights and volumes in each calculation zone (including skeleton core and viscera core)
        sumVolume = 0
        for layer in calculationLayers:
            if "Fat" in layer:
                calculationNodeValues[zone+"/"+layer]["Weight"] = (fatWeight * weightDistribution[layer][zone])
            else:
                calculationNodeValues[zone+"/"+layer]["Weight"] = (leanWeight * weightDistribution[layer][zone])

            calculationNodeValues[zone+"/"+layer]["Volume"] = (calculationNodeValues[zone+"/"+layer]["Weight"]/1000) + sumVolume
            sumVolume = calculationNodeValues[zone+"/"+layer]["Volume"]
        
        # Weights and thermal conductances of each thermal zone (unified core)
        thermalNodeValues[zone+"/Core"]["Weight"] = (
            calculationNodeValues[zone+"/Core skeleton"]["Weight"] + calculationNodeValues[zone+"/Core viscera"]["Weight"])
        thermalNodeValues[zone+"/Core"]["Capacitance"] = (
            calculationNodeValues[zone+"/Core skeleton"]["Weight"]*heatCapacity["Skeleton"] + 
            calculationNodeValues[zone+"/Core viscera"]["Weight"]*heatCapacity["Viscera"])
        
        thermalNodeValues[zone+"/Muscle"]["Weight"] = calculationNodeValues[zone+"/Muscle"]["Weight"]
        thermalNodeValues[zone+"/Muscle"]["Capacitance"] = calculationNodeValues[zone+"/Muscle"]["Weight"] * heatCapacity["Muscle"]

        thermalNodeValues[zone+"/Fat"]["Weight"] = calculationNodeValues[zone+"/Fat"]["Weight"]
        thermalNodeValues[zone+"/Fat"]["Capacitance"] = calculationNodeValues[zone+"/Fat"]["Weight"] * heatCapacity["Fat"]

        thermalNodeValues[zone+"/Skin"]["Weight"] = calculationNodeValues[zone+"/Skin"]["Weight"]
        thermalNodeValues[zone+"/Skin"]["Capacitance"] = calculationNodeValues[zone+"/Skin"]["Weight"] * heatCapacity["Skin"]

        # Lengths and radii of each calculation zone
        zoneArea = zoneAreas[zone]
        zoneLength = np.power(zoneArea, 2) / (4*np.pi * sumVolume)

        for layer_i, layer in enumerate(calculationLayers):
            calculationNodeValues[zone+"/"+layer]["Length"] = zoneLength

            calculationNodeValues[zone+"/"+layer]["Radius outer"] = np.power(
                calculationNodeValues[zone+"/"+layer]["Volume"]/(zoneLength*np.pi), 0.5)
            if layer_i == 0:
                radius_nm1 = 0
            else:
                radius_nm1 = calculationNodeValues[zone+"/"+calculationLayers[layer_i-1]]["Radius outer"]
            calculationNodeValues[zone+"/"+layer]["Radius"] = (radius_nm1 + calculationNodeValues[zone+"/"+layer]["Radius outer"]) / 2

            if layer_i > 1:
                thermalNodeValues[zone+"/"+layer]["Radius outer"] = calculationNodeValues[zone+"/"+layer]["Radius outer"]
                thermalNodeValues[zone+"/"+layer]["Radius"] = calculationNodeValues[zone+"/"+layer]["Radius"]
                thermalNodeValues[zone+"/"+layer]["Length"] = zoneLength
                #thermalNodeValues[zone+"/"+layer]["Volume"] = calculationNodeValues[zone+"/"+layer]["Volume"]
            elif layer_i == 1:
                thermalNodeValues[zone+"/Core"]["Radius outer"] = calculationNodeValues[zone+"/Core viscera"]["Radius outer"]
                thermalNodeValues[zone+"/Core"]["Radius"] = thermalNodeValues[zone+"/Core"]["Radius outer"]/2
                thermalNodeValues[zone+"/Core"]["Length"] = zoneLength
                #thermalNodeValues[zone+"/Core"]["Volume"] = calculationNodeValues[zone+"/Core viscera"]["Volume"]
        

        # Thermal conductance between zones
        # Calculate thermal conductance up to the next zone
        for layer_i, layer in enumerate(thermalLayers[:-1]):
            layerUp = thermalLayers[layer_i+1]
            firstThermalResistance = (
                np.log(thermalNodeValues[zone+"/"+layer]["Radius outer"]/thermalNodeValues[zone+"/"+layer]["Radius"])/
                (2*np.pi*thermalNodeValues[zone+"/"+layer]["Length"]*zoneThermalConductivity[layer][zone]))
            secondThermalResistance = (
                np.log(thermalNodeValues[zone+"/"+layerUp]["Radius"]/thermalNodeValues[zone+"/"+layer]["Radius outer"])/
                (2*np.pi*thermalNodeValues[zone+"/"+layerUp]["Length"]*zoneThermalConductivity[layerUp][zone]))
            combinedThermalConductance = 1/(firstThermalResistance+secondThermalResistance)

            thermalNodeValues[zone+"/"+layer]["Conductance out"] = combinedThermalConductance
            thermalNodeValues[zone+"/"+layer]["Conductance out to"] = zone+"/"+layerUp
        
        # Basal metabolism and bloodflow in each zone
        for layer in thermalLayers:
            heat = basalMetabolism[layer] * thermalNodeValues[zone+"/"+layer]["Weight"] / layerWeight[layer]
            if zone == "Head" and layer == "Core":
                heat += metabolismAreaFactor * totalSurfaceArea * 0.16
            elif zone == "Trunk" and layer == "Core":
                heat += metabolismAreaFactor * totalSurfaceArea * 0.58
            
            if layer == "Skin":
                bloodflow = skinBloodflowFactor[zone] * thermalNodeValues[zone+"/"+layer]["Weight"]
            elif layer == "Core" and zone == "Head":
                bloodflow = 45 * 1.060/3600
            elif layer == "Core" and zone == "Trunk":
                bloodflow = 210 * 1.060/3600
            else:
                bloodflow = heat * 1.2 * 1.060/4184

            thermalNodeValues[zone+"/"+layer]["Basal metabolism"] = heat
            thermalNodeValues[zone+"/"+layer]["Basal bloodflow"] = bloodflow


    thermalNodeValues["Trunk/Core"]["Weight"] -= 2.5 * 1.060
    thermalNodeValues["Blood"]["Weight"] += 2.5 * 1.060
    thermalNodeValues["Trunk/Core"]["Capacitance"] -= 2.5 * 1.060 * 3617 
    thermalNodeValues["Blood"]["Capacitance"] += 2.5 * 1.060 * 3617 
    #https://itis.swiss/virtual-population/tissue-properties/database/heat-capacity/ 

    for zone in zones:
        thermalNodeValues[zone+"/Skin"]["Area"] = zoneAreas[zone]
        thermalNodeValues[zone+"/Skin"]["Sunlit area vertical"] = np.power(thermalNodeValues[zone+"/Skin"]["Radius outer"], 2)
        thermalNodeValues[zone+"/Skin"]["Sunlit area horizontal"] =  2 * (thermalNodeValues[zone+"/Skin"]["Radius outer"] * 
                                                                          thermalNodeValues[zone+"/Skin"]["Length"])


    return(thermalNodeValues)

def stepWork(time):
    if time < 10*60*60:
        return(50)
    elif time < 12*60*60:
        return(500)
    else:
        return(150)

def makeBodyThermalConnections(solver, skinPlaceholderC=None):
    
    # Add the metabolic heat input to each layer
    for node in solver.thermalNodeValues:
        if "Blood" not in node:
            solver.components["Metabolism @ "+node] = bjls.HeatInputBoundary(QIn = solver.thermalNodeValues[node]["Basal metabolism"])
            solver.components["Metabolism @ "+node].QMin = solver.thermalNodeValues[node]["Basal metabolism"]    

    # Make the conduction links between each layer
    for node in solver.thermalNodeValues:
        if solver.thermalNodeValues[node]["Conductance out to"] != "":
            solver.components[node + " -> " + solver.thermalNodeValues[node]["Conductance out to"]] = bjls.ThermalConductance(
                C = solver.thermalNodeValues[node]["Conductance out"])
        elif "Skin" in node and skinPlaceholderC!=None:
            solver.components[node + " -> " + "Env"] = bjls.ThermalConductance(C = skinPlaceholderC)
        
        if "Blood" not in node:
            solver.components[node + " -> " + "Blood"] = bjls.CoolantPseudoConduction(
                SHC = 3617, mDot=solver.thermalNodeValues[node]["Basal bloodflow"])
            solver.components[node + " -> " + "Blood"].mDotMin= solver.thermalNodeValues[node]["Basal bloodflow"]
    
    if skinPlaceholderC != None:
        solver.components["@ Env"] = bjls.ConstantTemperatureBoundary(298)
    
    for zone in zones:
        solver.components["Metabolism @ "+zone+"/Skin"].area = solver.thermalNodeValues[zone+"/Skin"]["Area"]
    
    # Connect the nodes between conduction layers
    for c in solver.components:
        if "->" in c and "Blood" not in c:
            # Component is a thermal conduction type
            [inNode, outNode] = c.split(" -> ")
            try:
                for c2 in solver.components:
                    if c != c2 and "@ "+outNode in c2:
                        capacitance = solver.thermalNodeValues[outNode]["Capacitance"]
                        weight = solver.thermalNodeValues[outNode]["Weight"]
                        solver.connectStationsSolid(c, "outlet", c2, "outlet", mass=weight, specificHeatCapacity=capacitance/weight)
                        solver.components[c].outlet.ID = outNode
                        solver.components[c].outlet.state.ID = outNode
            except:
                if "Env" in c:
                    capacitance = solver.thermalNodeValues[inNode]["Capacitance"]
                    weight = solver.thermalNodeValues[inNode]["Weight"]
                    solver.connectStationsSolid(c, "outlet", "@ Env", "outlet", mass=1e5, specificHeatCapacity=1)
                else:
                    print("FAIL", c)
                    pass
        elif " @ " in c and "Blood" not in c:
            # Component is a metabolic heat input type
            outNode = c.split(" @ ")[1]
            try:
                for c2 in solver.components:
                    if c != c2 and outNode+" -> " in c2:
                        capacitance = solver.thermalNodeValues[outNode]["Capacitance"]
                        weight = solver.thermalNodeValues[outNode]["Weight"]
                        #print(c, c2, capacitance)
                        solver.connectStationsSolid(c, "outlet", c2, "inlet", mass=weight, specificHeatCapacity=capacitance/weight)
                        solver.components[c].outlet.ID = outNode
                        solver.components[c].outlet.state.ID = outNode
            except:
                print(outNode)
                return
        
        elif " -> Blood" in c and "Head/Core" not in c:
            bloodWeight = solver.thermalNodeValues["Blood"]["Weight"]
            bloodSHC = solver.thermalNodeValues["Blood"]["Capacitance"] / bloodWeight

            solver.connectStationsSolid(c, "outlet", "Head/Core -> Blood", "outlet", 
                                        mass=bloodWeight, specificHeatCapacity=bloodSHC)
    
    solver.connectStationsSolid("Head/Core -> Blood", "outlet", "Trunk/Core -> Blood", "outlet", 
                                        mass=bloodWeight, specificHeatCapacity=bloodSHC)
    
    if skinPlaceholderC != None:
        solver.components["@ Env"].outlet.ID = "Environment"
        solver.components["@ Env"].outlet.state.ID = "Environment"
    solver.components["Head/Core -> Blood"].outlet.ID = "Blood"
    solver.components["Head/Core -> Blood"].outlet.state.ID = "Blood"

def makeShirtsleveesThermalConnections(solver, 
                                       cabinPressurebar, cabinAirflowms, gravityRatio, TCabin, approxTGarment,
                                       cabinWallTemperature, cabinAirflowkgs, cabinAirflowInletHumidity, cabinAirflowCO2, cabinWallArea,
                                       cabinVolume,
                                       undergarmentEmissivity, undergarmentConductivities, undergarmentMasses):

    inletWaterDensity = HAPropsSI("W", "T", TCabin, "P", cabinPressurebar*1e5, "R", cabinAirflowInletHumidity)
    waterInletMass = inletWaterDensity * cabinAirflowkgs
    CO2InletMass = cabinAirflowCO2 * cabinAirflowkgs

    solver.components["Cabin air inlet"] = bjls.FluidSource(mDot=cabinAirflowkgs, fluid="Air", 
                                                            p=cabinPressurebar, T=TCabin,
                                                            dilutesFlow={"H2O":waterInletMass,
                                                                         "CO2": CO2InletMass})
    solver.components["Cabin air outlet"] = bjls.FluidSink(mDot=cabinAirflowkgs)


    solver.components["Cabin air flow"] = bjls.GeneralFluidLink(flowModel=True)
    solver.connectStationsFluid("Cabin air inlet", "fluidOutlet", "Cabin air flow", "fluidInlet",
                                fluid="Air", p_bar=cabinPressurebar, mDot=cabinAirflowkgs, volume=cabinVolume/2)
    solver.connectStationsFluid("Cabin air flow", "fluidOutlet", "Cabin air outlet", "fluidInlet",
                                fluid="Air", p_bar=cabinPressurebar, mDot=cabinAirflowkgs, volume=cabinVolume/2)
    

    enthalpyEvap = PropsSI("Hmass", "T", approxTGarment, "Q", 1, "Water") - PropsSI("Hmass", "T", approxTGarment, "Q", 0, "Water")

    solver.components["Cabin wall convection"] = bjls.GeneralFluidLink(flowModel=False, 
                                                            convectionModel={"active":True, "mode":"HTC",
                                                                             "area":cabinWallArea, "HTC":1, "useLMTD": False})
    solver.components["Cabin wall BC"] = bjls.ConstantTemperatureBoundary(cabinWallTemperature)
    solver.connectStationsFluid("Cabin air inlet", "fluidOutlet", "Cabin wall convection", "fluidInlet")
    solver.connectStationsFluid("Cabin air outlet", "fluidInlet", "Cabin wall convection", "fluidOutlet")
    solver.connectStationsSolid("Cabin wall convection", "solidLink", "Cabin wall BC", "outlet")

    solver.conductionComponents = []
    solver.convectionComponents = []
    solver.radiationComponents = []
    solver.LCGComponents = []

    for zone in zones:
        # Create the connection from garment node to airflow

        if undergarmentConductivities[zone] > 0:
            
            # Create the undergarment node with thermal conductivity
            undergarmSHC = 1300
            undergarmMass = undergarmentMasses[zone]*solver.thermalNodeValues[zone+"/Skin"]["Area"]
            solver.components[zone+"/Skin -> Garment"] = bjls.ThermalConductance(
                C=undergarmentConductivities[zone]*solver.thermalNodeValues[zone+"/Skin"]["Area"])
            solver.connectStationsSolid("Metabolism @ "+zone+"/Skin", "outlet", zone+"/Skin -> Garment", "inlet",
                                        mass=undergarmMass/2, specificHeatCapacity=undergarmSHC)
            
            solver.conductionComponents.append(zone+"/Skin -> Garment")

            solver.components[zone+"/Garment -> Cabin"] = bjls.GeneralFluidLink(flowModel=False,
                convectionModel={"active":True, "mode":"HTC", "HTC":1, "area":solver.thermalNodeValues[zone+"/Skin"]["Area"], "useLMTD": False}
            )
            solver.connectStationsSolid(zone+"/Skin -> Garment", "outlet", zone+"/Garment -> Cabin", "solidLink",
                                       mass=undergarmMass/2, specificHeatCapacity=undergarmSHC)
            solver.connectStationsFluid("Cabin air inlet", "fluidOutlet", zone+"/Garment -> Cabin", "fluidInlet")
            solver.connectStationsFluid(zone+"/Garment -> Cabin", "fluidOutlet", "Cabin air outlet", "fluidInlet")

            # Create the radiation from garment to cabin wall
            solver.components[zone+"/Garment | Cabin"] = bjls.RadiationTransfer(areaInlet=0.795*solver.thermalNodeValues[zone+"/Skin"]["Area"],
                                                                                emittanceInlet=undergarmentEmissivity, emittanceOutlet=1,
                                                                                absorbanceInlet=undergarmentEmissivity, absorbanceOutlet=1)
            solver.connectStationsSolid(zone+"/Garment | Cabin", "inlet", zone+"/Skin -> Garment", "outlet")
            solver.connectStationsSolid(zone+"/Garment | Cabin", "outlet", "Cabin wall BC", "outlet")
            
        else:
            # Create the convection from skin node to airflow
            solver.components[zone+"/Skin -> Cabin"] = bjls.GeneralFluidLink(flowModel=False,
                convectionModel={"active":True, "mode":"HTC", "HTC":1, "area":solver.thermalNodeValues[zone+"/Skin"]["Area"], "useLMTD": False}
            )
            solver.connectStationsSolid("Metabolism @ "+zone+"/Skin", "outlet", zone+"/Skin -> Cabin", "solidLink",
                                       mass=1, 
                                       specificHeatCapacity=1000)
            solver.connectStationsFluid(zone+"/Skin -> Cabin", "fluidInlet", "Cabin air inlet", "fluidOutlet")
            solver.connectStationsFluid(zone+"/Skin -> Cabin", "fluidOutlet", "Cabin air outlet", "fluidInlet")

            solver.convectionComponents.append(zone+"/Skin -> Cabin")

            # Create the radiation from skin node to cabin wall
            solver.components[zone+"/Skin | Cabin"] = bjls.RadiationTransfer(areaInlet=0.795*solver.thermalNodeValues[zone+"/Skin"]["Area"],
                                                                                emittanceInlet=undergarmentEmissivity, emittanceOutlet=1,
                                                                                absorbanceInlet=undergarmentEmissivity, absorbanceOutlet=1)
            solver.connectStationsSolid(zone+"/Skin | Cabin", "inlet", "Metabolism @ "+zone+"/Skin", "outlet")
            solver.connectStationsSolid(zone+"/Skin | Cabin", "outlet", "Cabin wall BC", "outlet")

            solver.radiationComponents.append(zone+"/Skin | Cabin")
            
        
        # Create the nodes for sweat
        solver.components[zone+" ~ Sweat"] = bjls.GeneralFluidLink(flowModel=False,
                dilutesModel={"active":True, "H2O":{"MTC":1, "usemDot":True, "mDot":0, "enthalpy":enthalpyEvap, 
                                                    "EMax": 1, "EMax without pressure": 1},
                              "CO2":{"usemDot":True, "mDot":0, "enthalpy":0}})

        solver.connectStationsSolid("Metabolism @ "+zone+"/Skin", "outlet", zone+" ~ Sweat", "solidLink", 
                                    mass=1, specificHeatCapacity=1)
        
        solver.connectStationsFluid(zone+" ~ Sweat", "fluidInlet", "Cabin air inlet", "fluidOutlet")
        solver.connectStationsFluid(zone+" ~ Sweat", "fluidOutlet", "Cabin air outlet", "fluidInlet")
    
    # Create the connections for respiratory heat transfer
    respiratoryNodes = ["Head/Core", "Head/Muscle", "Head/Fat", "Trunk/Core", "Trunk/Muscle"]
    for node in respiratoryNodes:
        solver.components[node+" ~ Breath"] = bjls.GeneralFluidLink(flowModel=False, 
                                convectionModel={"active":True, "mode":"simple", "conductance":0, "useLMTD": False},
                                dilutesModel={"active":True,
                                              "H2O": {"usemDot":False, "Q":0, "enthalpy": enthalpyEvap},
                                              "CO2": {"usemDot":True, "mDot":0, "enthalpy":0}})
        solver.connectStationsFluid(node+" ~ Breath", "fluidInlet", "Cabin air inlet", "fluidOutlet")
        solver.connectStationsFluid(node+" ~ Breath", "fluidOutlet", "Cabin air outlet", "fluidInlet")
        solver.connectStationsSolid(node+" ~ Breath", "solidLink", "Metabolism @ "+node, "outlet")
        
    
    solver.components["Cabin air inlet"].fluidOutlet.state.ID = "Cabin air in"
    solver.components["Cabin air outlet"].fluidInlet.state.ID = "Cabin air out"

    calculateShirtsleevesHTCs(solver, cabinPressurebar, cabinAirflowms, gravityRatio, TCabin, 
                              approxTGarment, cabinWallTemperature,  undergarmentConductivities)
    calculateShirtsleevesMTCs(solver, cabinPressurebar, cabinAirflowms, gravityRatio, TCabin, approxTGarment,
                                    cabinAirflowInletHumidity)

def calculateShirtsleevesHTCs(solver, 
                                       cabinPressurebar, cabinAirflowms, gravityRatio, TCabin, approxTGarment,
                                       cabinWallTemperature, undergarmentConductivities):
    # Metman gives 0.0212 * sqrt(P * v) * A * delta T * K_1
    # Following through appendix G.1 with SI units gets constant goes to 0.0110759
    heatTransferCoefficientForcedConvection = 0.0110759 * np.sqrt(cabinPressurebar*1e5 * cabinAirflowms)

    # Metman gives 0.06 * (G * P^2 * delta T)^0.25 * A * delta T * K_2
    # Following through Appendix G.1 with SI units gets constant = 0.043793
    approxDeltaT = approxTGarment - TCabin
    heatTransferCoefficientNaturalConvection = (
        0.0043793 * 1e2 * np.power(gravityRatio * (np.power(cabinPressurebar*10, 2) * approxDeltaT), 0.25))

    cabinWallForcedConvection = heatTransferCoefficientForcedConvection
    cabinWallNaturalConvection = (
        0.0043793 * 1e2 * np.power(gravityRatio * (np.power(cabinPressurebar*10, 2) * (TCabin - cabinWallTemperature)), 0.25))
    
    cabinWallHTC = max(cabinWallNaturalConvection, cabinWallForcedConvection)

    solver.components["Cabin wall convection"].convectionModel["HTC"] = cabinWallHTC
    
    for zone in zones:
        # Create the connection from garment node to airflow
        if heatTransferCoefficientForcedConvection > heatTransferCoefficientNaturalConvection:
            HTC = heatTransferCoefficientForcedConvection
            HTC *= np.sqrt(1/(2*solver.thermalNodeValues[zone+"/Skin"]["Radius outer"]))
        else:
            HTC = heatTransferCoefficientNaturalConvection
            HTC *= np.power(solver.thermalNodeValues[zone+"/Skin"]["Length"], -0.25)

        if undergarmentConductivities[zone] > 0:

            solver.components[zone+"/Garment -> Cabin"].convectionModel["HTC"] = HTC
        else:
            solver.components[zone+"/Skin -> Cabin"].convectionModel["HTC"] = HTC
    
def calculateShirtsleevesMTCs(solver, 
                                       cabinPressurebar, cabinAirflowms, gravityRatio, TCabin, approxTGarment,
                                    cabinAirflowInletHumidity):
    # Create the basic mass transfer terms for sweating
    evaporationCoefficientForcedConvection = (3.071e-3 * 
                                              np.power(cabinAirflowms/(cabinPressurebar*1e5), 0.5) * 
                                              np.power(TCabin, 1.04))
    cabinPressurePa = cabinPressurebar * 1e5
    cabinPressurepsi = cabinPressurePa * 1.45e-4
    cabinTempR = TCabin * 9/5
    inletWaterPressureSaturationpsi = PropsSI("P", "T", TCabin, "Q", 1, "Water") * 1.45e-4
    undergarmentWaterPressureSaturationpsi = PropsSI("P", "T", approxTGarment, "Q", 1, "Water") * 1.45e-4
    deltaTR = (approxTGarment - TCabin) * 9/5

    evaporationCoefficientNaturalConvection = (3.175e-5 * (cabinTempR/cabinPressurepsi) * np.power(
        cabinPressurePa * gravityRatio * (0.005 * cabinPressurepsi * deltaTR + 
                                          1.02 * (undergarmentWaterPressureSaturationpsi - inletWaterPressureSaturationpsi)), 0.25)
    )
    enthalpyEvap = PropsSI("Hmass", "T", approxTGarment, "Q", 1, "Water") - PropsSI("Hmass", "T", approxTGarment, "Q", 0, "Water")

    for zone in zones:
        if evaporationCoefficientForcedConvection > evaporationCoefficientNaturalConvection:
            # Use forced convection mass transfer
            MTC = evaporationCoefficientForcedConvection
            MTC *= np.sqrt(1/(2*solver.thermalNodeValues[zone+"/Skin"]["Radius outer"]))
        else:
            MTC = evaporationCoefficientNaturalConvection
            MTC *= np.power(3.5/(solver.thermalNodeValues[zone+"/Skin"]["Length"]/0.3), 0.25)
        
        EmaxWithoutPressure = MTC * enthalpyEvap * (solver.thermalNodeValues[zone+"/Skin"]["Area"] / (287 * TCabin))

        print(zone, MTC)
        
        EMax = EmaxWithoutPressure * (
            PropsSI("P", "T", approxTGarment, "Q", 1, "Water") - 
            HAPropsSI("Tdp", "T", TCabin, "P", cabinPressurePa, "RH", cabinAirflowInletHumidity)
        )

        solver.components[zone+" ~ Sweat"].dilutesModel["H2O"]["MTC"] = MTC
        solver.components[zone+" ~ Sweat"].dilutesModel["H2O"]["EMax"] = EMax
        solver.components[zone+" ~ Sweat"].dilutesModel["H2O"]["EMax without pressure"] = EmaxWithoutPressure
    
    exit()

def calculateSpacesuitInternalThermalConnections(solver, 
                                       suitPressureBar, TInlet, approxTGarment, inletHumidity, inletCO2, 
                                       undergarmentEmissivity, suitInnerEmissivity, suitOuterEmissivity, suitOuterAbsorbivity, suitSHC, 
                                       undergarmentConductivities, undergarmentMasses,
                                       suitMasses, suitConductivities, 
                                       suitAreas,
                                       useLCG, LCGTInlet, LCGFlowRates, 
                                       zoneAirflows, zoneAirConnectivity, zoneAirVolumes, useFixedSource=True,
                                       externalIsVacuum=True, externalT=250, externalEmissivity=0.8):
    
    inletFlow = zoneAirflows["Air inlet"]
    inletWaterDensity = HAPropsSI("W", "T", TInlet, "P", suitPressureBar*1e5, "R", inletHumidity)
    waterInletMass = inletWaterDensity * inletFlow
    CO2InletMass = inletCO2 * inletFlow
    enthalpyEvap = PropsSI("Hmass", "T", approxTGarment, "Q", 1, "Water") - PropsSI("Hmass", "T", approxTGarment, "Q", 0, "Water")

    LCGFlowsNumbers = [LCGFlowRates[i] for i in LCGFlowRates.keys()]
    LCGFlowTotal = sum(LCGFlowsNumbers)
    LCGNumLoops = sum([1 for i in LCGFlowsNumbers if i>0])
    LCGFlowRatio = [i/LCGFlowTotal for i in LCGFlowsNumbers if i>0]

    solver.conductionComponents = []
    solver.convectionComponents = []
    solver.radiationComponents = []
    solver.LCGComponents = []

    solver.components["Environment"] = bjls.ConstantTemperatureBoundary(externalT)
    if not externalIsVacuum:
        pass
        # TODO Add nodes and conductivity for suit-inside-cabin
        

    if useFixedSource:
        solver.components["In Air"] = bjls.FluidSource(mDot=inletFlow, fluid="Air", 
                                                            p=suitPressureBar, T=TInlet,
                                                            dilutesFlow={"H2O":waterInletMass,
                                                                         "CO2": CO2InletMass})

        solver.components["Out Air"] = bjls.FluidSink()

        if LCGNumLoops > 0 and useLCG:
            solver.components["LCG In"] = bjls.FluidSource(mDot=LCGFlowTotal, fluid="Water",
                                                        T=LCGTInlet)
            solver.components["LCG Inlet Manifold"] = bjls.FluidJunction(1, LCGNumLoops, LCGFlowRatio)
            solver.components["LCG Outlet Manifold"] = bjls.FluidJunction(LCGNumLoops, 1, [1])
            solver.components["LCG Out"] = bjls.FluidSink()
    else:
        # TODO Leave the inlet/outlet hanging, so it can be attached to the PLSS
        pass

    if useLCG:
        # Set up the LCG
        if LCGNumLoops > 0:
            solver.connectStationsFluid("LCG In", "fluidOutlet", "LCG Inlet Manifold", "fluidInlet0",
                                        TRef=LCGTInlet, fluid="Water")
            solver.connectStationsFluid("LCG Outlet Manifold", "fluidOutlet0", "LCG Out", "fluidInlet",
                                        TRef=LCGTInlet, fluid="Water")
            
            zoneCounter = 0
            for zone in zones:
                if LCGFlowRates[zone] > 0:
                    solver.components[zone+"/LCG"] = bjls.GeneralFluidLink(convectionModel={"active":True, "useLMTD":False,
                                                                                            "mode":"NTU", "NTU": 1})
                    solver.connectStationsFluid("LCG Inlet Manifold", "fluidOutlet"+str(zoneCounter), zone+"/LCG", "fluidInlet",
                                                fluid="Water", TRef=LCGTInlet)
                    solver.connectStationsFluid(zone+"/LCG", "fluidOutlet", "LCG Outlet Manifold", "fluidInlet"+str(zoneCounter),
                                                fluid="Water", TRef=LCGTInlet)
                    solver.connectStationsSolid(zone+"/LCG", "solidLink", "Metabolism @ "+zone+"/Skin", "outlet")
                    zoneCounter += 1

                    solver.LCGComponents.append(zone+"/LCG")

    
    outputsFromZone = copy(valueDict)
    inputsToZone = copy(valueDict)
    for z in outputsFromZone.keys():
        outputsFromZone[z] = []
        inputsToZone[z] = []
    outputsFromZone["In"] = []
    inputsToZone["In"] = []
    inputsToZone["Out"] = []

    zonesPlusOut = zones + ["Out"]
    # Create the airflow nodes to allow for hooking them together in order
    for zone in zones:
        solver.components[zone+" Air"] = bjls.GeneralFluidLink(flowModel=True)

        inputsToZone[zone].append([zoneAirConnectivity[zone]["From"], zoneAirConnectivity[zone]["Fraction"]])
        outputsFromZone[zoneAirConnectivity[zone]["From"]].append([zone, zoneAirConnectivity[zone]["Fraction"]])
    
    outputsFromZone["Out"] = []
    
    for zone in zonesPlusOut:
        if outputsFromZone[zone] == [] and "Out" not in zone:
            inputsToZone["Out"].append([zone, 1])
            outputsFromZone[zone].append(["Out", 1])
        if len(outputsFromZone[zone]) > 1:
            # Create a splitter junction
            solver.components[zone+" Air Outlet Junction"] = bjls.FluidJunction(1, len(outputsFromZone[zone]),
                                                                                [0 for i in outputsFromZone[zone]])
        if len(inputsToZone[zone]) > 1:
            # Create a splitter junction
            solver.components[zone+" Air Inlet Junction"] = bjls.FluidJunction(len(inputsToZone[zone]), 1)

    # Connections use the terminology
    # (outletZone, inletZone, outletIndex, inletIndex, outletFraction)
    # outletIndex and inletIndex is added to the string for connection - eg fluidOutlet0
    # If outletIndex or inletIndex are None, they are not added to the string
    connections = []
    for zone in zonesPlusOut:
        # Case 1 - goes from 1 zone to another zone
        # Case 2 - goes from 1 zone to outlet junction, outlet junction to zone
        # Case 3 - goes from 1 zone to inlet junction, inlet junction to zone
        # Case 4 - goes from 1 zone to outlet junction, outlet junction to inlet junction, inlet junction to zone
        if len(inputsToZone[zone]) == 1:
            # Only one input to the zone - Case 1 or 2
            if len(outputsFromZone[inputsToZone[zone][0][0]]) == 1:
                # The inputting zone only outputs to one zone - Case 1
                connections.append((inputsToZone[zone][0][0]+" Air", zone+" Air", None, None, 1, zoneAirflows[zone], zoneAirVolumes[zone]))
            
            elif len(outputsFromZone[inputsToZone[zone][0][0]]) > 1:
                # The inputting zone outputs to more than 1 zone - Case 2
                outputtingZonesFromInput = [i[0] for i in outputsFromZone[inputsToZone[zone][0][0]]]
                connections.append((inputsToZone[zone][0][0]+" Air", inputsToZone[zone][0][0]+" Air Outlet Junction", None, 0, 1, 
                                    zoneAirflows[inputsToZone[zone][0][0]], zoneAirVolumes[inputsToZone[zone][0][0]]))
                connections.append((inputsToZone[zone][0][0]+" Air Outlet Junction", zone+" Air", outputtingZonesFromInput.index(zone), None, 
                                    inputsToZone[zone][0][1], zoneAirflows[inputsToZone[zone][0][0]] * inputsToZone[zone][0][1], 1e-2))
        
        elif len(inputsToZone[zone]) > 1:
            # More than one input to the zone - Case 3 or 4
            inputtingZonesToOutput = [i[0] for i in inputsToZone[zone]]
            totalMassFlow = 0
            for inputtingZone_i, inputtingZone in enumerate(inputsToZone[zone]):
                if len(outputsFromZone[inputtingZone[0]]) == 1:
                    # The inputting zone only has a single output - Case 3
                    connections.append((inputtingZone[0]+" Air", zone+" Air Inlet Junction", None, 
                                        inputtingZonesToOutput.index(inputtingZone[0]), inputtingZone[1], 
                                        zoneAirflows[inputtingZone[0]], zoneAirVolumes[inputtingZone[0]]))
                    totalMassFlow += zoneAirflows[inputtingZone[0]]
                else:
                    # Case 4
                    raise ValueError("Connectivity case not handled yet")
            connections.append((zone+" Air Inlet Junction", zone+" Air", 0, None, 1, totalMassFlow, 1e-2))
    
    connections = list(set(connections))
    
    for c in connections:
        print(c)
        connection1Name = "fluidOutlet"
        connection2Name = "fluidInlet"
        if c[2] != None:
            connection1Name += str(c[2])
        if c[3] != None:
            connection2Name += str(c[3])

        solver.connectStationsFluid(c[0], connection1Name, c[1], connection2Name, mDot=c[5], volume=c[6], 
                                    fluid="Air", p_bar=suitPressureBar, TRef=TInlet)
        if "Outlet Junction" in c[0]:
            solver.components[c[0]].outputSplitRatio[c[2]] = c[4]
            for i in solver.components[c[0]].linkedStations: i.countErrors = False
    
    for zone in zones:
        # # Create the suit shell
        solver.components[zone+"/Suit"] = bjls.ThermalConductance(suitConductivities[zone])
        solver.components[zone+"/Suit | Environment"] = bjls.RadiationTransfer(areaInlet=suitAreas[zone],
                                                                               areaOutlet=1, 
                                                                               emittanceInlet=suitOuterEmissivity,
                                                                               emittanceOutlet=externalEmissivity, 
                                                                               absorbanceInlet=suitOuterAbsorbivity,
                                                                               viewFactor=1)
        solver.components[zone+"/Suit | Environment"] = bjls.ThermalConductance()
        solver.connectStationsSolid(zone+"/Suit", "outlet", zone+"/Suit | Environment", "inlet", 
                                    mass=suitMasses[zone]/2, specificHeatCapacity = suitSHC)
        solver.connectStationsSolid(zone+"/Suit | Environment", "outlet", "Environment", "outlet",
                                    mass=1e10, specificHeatCapacity = suitSHC)

        
        # Create the heat transfer from the air to the suit inner
        solver.components[zone+"/Suit -> Air"] = bjls.GeneralFluidLink(flowModel=False,
                convectionModel={"active":True, "mode":"HTC", "HTC":0, "area":solver.thermalNodeValues[zone+"/Skin"]["Area"], "useLMTD": False}
            )
        solver.connectStationsFluid(zone+" Air", "fluidInlet", zone+"/Suit -> Air", "fluidInlet")
        solver.connectStationsFluid(zone+"/Suit -> Air", "fluidOutlet", zone+" Air", "fluidOutlet")
        solver.connectStationsSolid(zone+"/Suit -> Air", "solidLink", zone+"/Suit", "inlet", 
                                    mass=suitMasses[zone]/2, specificHeatCapacity = suitSHC)
        

        if undergarmentConductivities[zone] > 0:
            
            # Create the undergarment node with thermal conductivity
            undergarmSHC = 1300
            undergarmMass = undergarmentMasses[zone]*solver.thermalNodeValues[zone+"/Skin"]["Area"]
            solver.components[zone+"/Skin -> Garment"] = bjls.ThermalConductance(
                C=undergarmentConductivities[zone]*solver.thermalNodeValues[zone+"/Skin"]["Area"])
            
            solver.connectStationsSolid("Metabolism @ "+zone+"/Skin", "outlet", zone+"/Skin -> Garment", "inlet",
                                        mass=undergarmMass/2, specificHeatCapacity=undergarmSHC)
            
            solver.conductionComponents.append(zone+"/Skin -> Garment")

            solver.components[zone+"/Garment -> Air"] = bjls.GeneralFluidLink(flowModel=False,
                convectionModel={"active":True, "mode":"HTC", "HTC":0, "area":solver.thermalNodeValues[zone+"/Skin"]["Area"], "useLMTD": False}
            )
            solver.connectStationsSolid(zone+"/Skin -> Garment", "outlet", zone+"/Garment -> Air", "solidLink",
                                       mass=undergarmMass/2, specificHeatCapacity=undergarmSHC)
            
            solver.connectStationsFluid(zone+" Air", "fluidInlet", zone+"/Garment -> Air", "fluidInlet")
            solver.connectStationsFluid(zone+"/Garment -> Air", "fluidOutlet", zone+" Air", "fluidOutlet")

            # Create the radiation from garment to suit inner surface
            solver.components[zone+"/Garment | Suit"] = bjls.RadiationTransfer(areaInlet=0.795*solver.thermalNodeValues[zone+"/Skin"]["Area"],
                                                                                emittanceInlet=undergarmentEmissivity, emittanceOutlet=1,
                                                                                absorbanceInlet=undergarmentEmissivity, absorbanceOutlet=1)
            solver.connectStationsSolid(zone+"/Garment | Suit", "inlet", zone+"/Skin -> Garment", "outlet")
            solver.connectStationsSolid(zone+"/Garment | Suit", "outlet", zone+"/Suit", "inlet")
            
        else:
            # Create the convection from skin node to airflow
            solver.components[zone+"/Skin -> Air"] = bjls.GeneralFluidLink(flowModel=False,
                convectionModel={"active":True, "mode":"HTC", "HTC":0, "area":solver.thermalNodeValues[zone+"/Skin"]["Area"], "useLMTD": False}
            )
            solver.connectStationsSolid("Metabolism @ "+zone+"/Skin", "outlet", zone+"/Skin -> Air", "solidLink",
                                       mass=1, 
                                       specificHeatCapacity=1000)
            
            solver.connectStationsFluid(zone+" Air", "fluidInlet", zone+"/Skin -> Air", "fluidInlet")
            solver.connectStationsFluid(zone+"/Skin -> Air", "fluidOutlet", zone+" Air", "fluidOutlet")

            solver.convectionComponents.append(zone+"/Skin -> Air")

            # Create the radiation from skin node to cabin wall
            solver.components[zone+"/Skin | Suit"] = bjls.RadiationTransfer(areaInlet=0.795*solver.thermalNodeValues[zone+"/Skin"]["Area"],
                                                                                emittanceInlet=undergarmentEmissivity, emittanceOutlet=suitInnerEmissivity,
                                                                                absorbanceInlet=undergarmentEmissivity, absorbanceOutlet=suitInnerEmissivity)

            solver.connectStationsSolid(zone+"/Skin | Suit", "inlet", zone+"/Skin -> Air", "solidLink")
            solver.connectStationsSolid(zone+"/Skin | Suit", "outlet", zone+"/Suit", "inlet")

            solver.radiationComponents.append(zone+"/Skin | Suit")
            
        
        # Create the nodes for sweat
        solver.components[zone+" ~ Sweat"] = bjls.GeneralFluidLink(flowModel=False,
                dilutesModel={"active":True, "H2O":{"MTC":0, "usemDot":True, "mDot":0, "enthalpy":enthalpyEvap, 
                                                    "EMax": 1, "EMax without pressure": 1},
                              "CO2":{"usemDot":True, "mDot":0, "enthalpy":0}})

        solver.connectStationsSolid("Metabolism @ "+zone+"/Skin", "outlet", zone+" ~ Sweat", "solidLink", 
                                    mass=1, specificHeatCapacity=1)
        
        solver.connectStationsFluid(zone+" ~ Sweat", "fluidInlet", zone+" Air", "fluidInlet")
        solver.connectStationsFluid(zone+" ~ Sweat", "fluidOutlet", zone+" Air", "fluidOutlet")
    
    # Create the connections for respiratory heat transfer
    respiratoryNodes = ["Head/Core", "Head/Muscle", "Head/Fat", "Trunk/Core", "Trunk/Muscle"]
    for node in respiratoryNodes:
        solver.components[node+" ~ Breath"] = bjls.GeneralFluidLink(flowModel=False, 
                                convectionModel={"active":True, "mode":"simple", "conductance":0, "useLMTD": False},
                                dilutesModel={"active":True,
                                              "H2O": {"usemDot":False, "Q":0, "enthalpy": enthalpyEvap},
                                              "CO2": {"usemDot":True, "mDot":0, "enthalpy":0}})
        solver.connectStationsFluid(node+" ~ Breath", "fluidInlet", "Head Air", "fluidOutlet")
        solver.connectStationsFluid(node+" ~ Breath", "fluidOutlet", "Head Air", "fluidInlet")
        solver.connectStationsSolid(node+" ~ Breath", "solidLink", "Metabolism @ "+node, "outlet")
    
    # Connect all zones to the outlet that don't currently have a connection on outflow
    for zone in outputsFromZone:
        if zone != "Out":
            solver.components[zone+" Air"].fluidOutlet.state.ID = zone+" Air"
    
    solver.components["Out Air"].fluidInlet.state.ID = "Out Air"

    calculateSpacesuitInternalHTCs(solver, zoneAirflows, undergarmentConductivities)
    calculateSpacesuitInternalMTCs(solver, approxTGarment, suitPressureBar)
    calculateLCGPerformance(solver, 70)

def calculateSpacesuitExternalThermalConnections(solver,
                                                 suitAreas, suitMasses, suitSHC,
                                                 suitOuterEmissivity, suitOuterAbsorbivity,
                                                 groundTemp, groundViewFactor,
                                                 skyTemp,
                                                 sunAltitude, sunIntensity,
                                                 isVacuum, airTemp):
    
    
    solver.components["Ground"] = bjls.ConstantTemperatureBoundary(groundTemp)
    solver.components["Sky"] = bjls.ConstantTemperatureBoundary(skyTemp)

    for zone in zones:
        # Create radiative transfer to ground
        solver.components[zone+"/Suit | Ground"] = bjls.RadiationTransfer(
            areaInlet=suitAreas[zone], emittanceInlet=solver.controlParameters["Suit outer emissivity"], 
            absorbanceInlet=solver.controlParameters["Suit outer absorbivity"],
            viewFactor=groundViewFactor
        )
        solver.connectStationsSolid(zone+"/Suit", "outlet", zone+"/Suit | Ground", "inlet", 
                                    mass=suitMasses[zone]/2, specificHeatCapacity = suitSHC)
        solver.connectStationsSolid(zone+"/Suit | Ground", "outlet", "Ground", "outlet",
                                    mass=1e10, specificHeatCapacity = suitSHC)
        
        # Create radiative transfer to sky

        solver.components[zone+"/Suit | Sky"] = bjls.RadiationTransfer(
            areaInlet=suitAreas[zone], emittanceInlet=solver.controlParameters["Suit outer emissivity"], 
            absorbanceInlet=solver.controlParameters["Suit outer absorbivity"],
            viewFactor=1-groundViewFactor
        )
        solver.connectStationsSolid(zone+"/Suit", "outlet", zone+"/Suit | Sky", "inlet", 
                                    mass=suitMasses[zone]/2, specificHeatCapacity = suitSHC)
        solver.connectStationsSolid(zone+"/Suit | Sky", "outlet", "Sky", "outlet",
                                    mass=1e10, specificHeatCapacity = suitSHC)
        
        areaIlluminated = ((solver.thermalNodeValues[zone+"/Skin"]["Sunlit area vertical"] * np.sin(np.deg2rad(sunAltitude))) + 
                           (solver.thermalNodeValues[zone+"/Skin"]["Sunlit area horizontal"] * np.cos(np.deg2rad(sunAltitude))))
        
        solarPower = areaIlluminated * sunIntensity * solver.controlParameters["Suit outer absorbivity"]

        solver.components["Sunlight "+zone] = bjls.HeatInputBoundary(solarPower)
        solver.connectStationsSolid("Sunlight "+zone, "outlet", zone+"/Suit", "outlet")

def spacesuitGetExternalHeatLoads(solver, 
                                  groundTemp,
                                    skyTemp,
                                    sunAltitude, sunIntensity,
                                    isVacuum, airTemp):
    solver.components["Ground"].TSet = groundTemp
    solver.components["Sky"].TSet = skyTemp
    for zone in zones:
        solver.components[zone+"/Suit | Ground"].viewFactor = groundTemp
        solver.components[zone+"/Suit | Sky"].viewFactor = 1 - groundTemp

        areaIlluminated = ((solver.thermalNodeValues[zone+"/Skin"]["Sunlit area vertical"] * np.sin(np.deg2rad(sunAltitude))) + 
                           (solver.thermalNodeValues[zone+"/Skin"]["Sunlit area horizontal"] * np.cos(np.deg2rad(sunAltitude))))
        
        solarPower = areaIlluminated * sunIntensity * solver.controlParameters["Suit outer absorbivity"]
        solver.components["Sunlight "+zone].QIn = solarPower

def calculateSpacesuitInternalHTCs(solver,
                                   zoneAirflows, undergarmentConductivities):
    for zone in zones:
        heatTransferCoefficientForcedConvection = 14.99 * np.power(zoneAirflows[zone], 1/3) * np.power(1.067 / solver.thermalNodeValues[zone+"/Skin"]["Length"], 1/3)

        if undergarmentConductivities[zone] > 0:
            solver.components[zone+"/Garment -> Air"].convectionModel["HTC"] = heatTransferCoefficientForcedConvection
        else:
            solver.components[zone+"/Skin -> Air"].convectionModel["HTC"] = heatTransferCoefficientForcedConvection
        
        solver.components[zone+"/Suit -> Air"].convectionModel["HTC"] = heatTransferCoefficientForcedConvection

def calculateSpacesuitInternalMTCs(solver, 
                                       approxTGarment, suitPressureBar, useRealHumidity=False):
    
    enthalpyEvap = PropsSI("Hmass", "T", approxTGarment, "Q", 1, "Water") - PropsSI("Hmass", "T", approxTGarment, "Q", 0, "Water")
    
    for zone in zones:
        mDot = solver.components[zone+" ~ Sweat"].fluidInlet.state.mDot
        rho = solver.components[zone+" ~ Sweat"].fluidInlet.state.rho
        TAir = solver.components[zone+" ~ Sweat"].fluidInlet.state.T
        TSkin = approxTGarment
        L = solver.thermalNodeValues[zone+"/Skin"]["Length"]

        MTC = 0.01957 * np.power(mDot * np.power(TSkin, 3.62) / (rho * np.power(suitPressureBar*1e5, 2)), 1/3)
        MTC *= np.power(1/L, 1/3)

        EMaxWithoutPressure = MTC * enthalpyEvap * (solver.thermalNodeValues[zone+"/Skin"]["Area"] / (287 * TAir))

        if useRealHumidity:
            EMax = EMaxWithoutPressure * (
                PropsSI("P", "T", approxTGarment, "Q", 1, "Water") - 
                HAPropsSI("Tdp", "T", TAir, "P", suitPressureBar*1e5, "W", solver.components[zone+" ~ Sweat"].fluidInlet.state.dilutesConcentration["H2O"])
            )
        else:
            EMax = EMaxWithoutPressure * (
                PropsSI("P", "T", approxTGarment, "Q", 1, "Water") - 
                HAPropsSI("Tdp", "T", TAir, "P", suitPressureBar*1e5, "RH", 0)
            )

        solver.components[zone+" ~ Sweat"].dilutesModel["H2O"]["MTC"] = MTC
        solver.components[zone+" ~ Sweat"].dilutesModel["H2O"]["EMax"] = EMax
        solver.components[zone+" ~ Sweat"].dilutesModel["H2O"]["EMax without pressure"] = EMaxWithoutPressure

        print(zone, MTC, EMax)

def calculateLCGPerformance(solver, metabolicRate):
    UACorrelationMetabolicRate = [0, 14.655, 29.31, 43.965, 58.62, 73.275, 87.93, 102.585, 205.17, 219.825, 234.48, 249.135, 263.79, 586.2]
    UACorrelationUA = [0, 5.275, 15.03375, 19.78125, 22.41875, 24.00125, 25.05625, 26.11125, 34.2875, 35.3425, 36.13375, 36.66125, 36.925, 36.925]

    LCGUA = np.interp(metabolicRate, UACorrelationMetabolicRate, UACorrelationUA)

    for component in solver.LCGComponents:
        solver.components[component].convectionModel["NTU"] = LCGUA / (
            solver.components[component].fluidInlet.state.mDot * solver.components[component].fluidInlet.state.cp)

def systemController(self, time, dt):
    # For each node, calculate the difference to the setpoint
    temperatureErrors = {
        "Core": copy(valueDict),
        "Muscle": copy(valueDict),
        "Fat": copy(valueDict),
        "Skin": copy(valueDict)
    }
    energyStorage = 0
    for layer, layerDict in temperatureErrors.items():
        for zone in layerDict:
            Tnode = self.components["Metabolism @ "+zone+"/"+layer].outlet.state.T
            error = self.controlParameters["Set points"][layer][zone] - Tnode
            temperatureErrors[layer][zone] = error
            energyStorage += -error * self.thermalNodeValues[zone+"/"+layer]["Capacitance"]
    
    # Calculate the control variable for shivering
    shiveringHeat = self.controlParameters["Shivering coefficient"] * (
        temperatureErrors["Core"]["Head"]) * sum(
            max(0, temperatureErrors["Skin"][i]*self.controlParameters["Skin compartment distribution"][i]) for i in
            temperatureErrors["Skin"]
        )

    # Calculate heat generation in each muscle zone
    if callable(self.workHeat):
        workHeat = self.workHeat(time)
    else:
        workHeat = self.workHeat

    muscleWorks = copy(valueDict)
    for zone in muscleWorks:
        QBasal = self.components["Metabolism @ "+zone+"/Muscle"].QMin
        QWork = self.controlParameters["Work distribution"][zone] * workHeat
        QShiv = self.controlParameters["Shivering distribution"][zone] * shiveringHeat
        muscleWorks[zone] = QBasal + QWork + QShiv
    
    totalMetabolicHeat = sum(muscleWorks.values())


    # Calculate the blood flow targets
    muscleBloodflows = copy(valueDict)
    skinBloodFlows = copy(valueDict)
    
    for zone in muscleBloodflows:
        muscleBloodflows[zone] = (self.controlParameters["Bloodflow per W muscles"] * muscleWorks[zone])
        
    dilationMassFlow = self.controlParameters["Vasodilation coefficient"] * (
        self.components["Metabolism @ Head/Core"].outlet.state.T - self.controlParameters["Set points"]["Core"]["Head"])
    dilationMassFlow = max(dilationMassFlow, 0)

    constrictionControl = self.controlParameters["Vasoconstriction coefficient"] * (
        temperatureErrors["Core"]["Head"]) * sum(
            max(0, temperatureErrors["Skin"][i]*self.controlParameters["Skin compartment distribution"][i]) for i in
            temperatureErrors["Skin"]
        )

    for zone in skinBloodFlows:
        mDotBasal = self.components[zone+"/Skin -> Blood"].mDotMin
        skinBloodFlows[zone] = ((mDotBasal + self.controlParameters["Vasodilation distribution"][zone]*dilationMassFlow) / 
                                (1 + self.controlParameters["Vasoconstriction distribution"][zone]*constrictionControl))
        
    # Calculate the control value for sweating
    sweatControl = sum(max(0, -temperatureErrors["Skin"][i]*self.controlParameters["Skin compartment distribution"][i]) for i in
            temperatureErrors["Skin"])
    sweatControl = 1591.2 + 132.12 * sweatControl
    sweatControl *= max(0, -temperatureErrors["Core"]["Head"])
    targetSweatRates = copy(valueDict)
    for zone in targetSweatRates:
        targetSweatRates[zone] = (sweatControl * self.controlParameters["Sweat distribution"][zone] * 
                                  np.power(2, -temperatureErrors["Skin"][zone]/4))
    
    # Recalculate maximum values of sweating
    if self.controlParameters["Recalculate max sweat rate"]:
        pass
        # TODO - find the new EMax by doing EMax without pressure * (Psat(T_skin) - P_dew(air))
    else:
        pass

    # Calculate the quantites for diffusion
    diffusionRates =  copy(valueDict)
    for zone in diffusionRates:
        diffusionRate = 0.003047 * self.components["Metabolism @ "+zone+"/Skin"].area * (
            PropsSI("P", "T", self.components["Metabolism @ "+zone+"/Skin"].outlet.state.T, "Q", 1, "Water") - 
            (1e5 * self.components[zone+" ~ Sweat"].fluidInlet.state.p * 
             self.components[zone+" ~ Sweat"].fluidInlet.state.dilutesConcentration["H2O"])
        )
        diffusionRates[zone] = max(diffusionRate, 0)
    
    # Calculate the average skin temperature
    TResp = 0
    for node in self.controlParameters["Respiratory temperature distribution"]:
        TResp += (self.controlParameters["Respiratory temperature distribution"][node] * 
                  self.components["Metabolism @ "+node].outlet.state.T)
    vapourPressureResp = PropsSI("P", "T", TResp, "Q", 1, "Water")
    vapourPressureRespInlet = (self.components[zone+" ~ Sweat"].fluidInlet.state.p * 1e5 * 
                               self.components[zone+" ~ Sweat"].fluidInlet.state.dilutesConcentration["H2O"])

    convectionConductanceRespiratory = copy(self.controlParameters["Respiratory temperature distribution"])
    evaporationPowerRespiratory = copy(self.controlParameters["Respiratory temperature distribution"])
    QConvectionResp = 0
    for node in convectionConductanceRespiratory:
        convectionConductanceRespiratory[node] = (14.5308 * self.components[node+" ~ Breath"].fluidInlet.state.p / 
                                                  (48.3 * 9/5 * self.components[zone+" ~ Sweat"].fluidInlet.state.T))
        evaporationPowerRespiratory[node] = convectionConductanceRespiratory[node] * 1.0

        convectionConductanceRespiratory[node] *= totalMetabolicHeat * 0.24 * 4200 * 0.2931 * 5/9
        convectionConductanceRespiratory[node] *= self.controlParameters["Respiratory heat ratios"][node]
        QConvectionResp += abs(convectionConductanceRespiratory[node])

        evaporationPowerRespiratory[node] *= 18 * 1040 / 32
        evaporationPowerRespiratory[node] *= (vapourPressureResp - 0.8 * vapourPressureRespInlet) / (
            self.components[zone+" ~ Sweat"].fluidInlet.state.p * 1e5)
        evaporationPowerRespiratory[node] *= self.controlParameters["Respiratory heat ratios"][node]
    
    # If the LCG is active - calculate the coefficients
    if len(self.LCGComponents) > 0:
        calculateLCGPerformance(self, totalMetabolicHeat)
    
    # Calculate rates of CO2 consumption, oxygen production
    mDotO2 = totalMetabolicHeat * 4.299e-4 * (2.0265e-4 -  4.5055e-5 * self.controlParameters["Respiratory quotient"])
    mDotCO2 = mDotO2 * 44/32

    # Calculate cabin dew point
    cabinDewPoint = HAPropsSI("Tdp", "T", self.components["Head/Core ~ Breath"].fluidOutlet.state.T, 
                              "P", self.components["Head/Core ~ Breath"].fluidOutlet.state.p*1e5,
                              "W", self.components["Head/Core ~ Breath"].fluidOutlet.state.dilutesConcentration["H2O"])

    # Calculate average temperatures for storage
    TSkinAv = 0
    TMuscleAv = 0
    TGarmAv = 0
    for zone in zones:
        TSkinAv += self.components["Metabolism @ "+zone+"/Skin"].outlet.state.T * self.controlParameters["Skin temperature distribution"][zone]
        TMuscleAv += self.components["Metabolism @ "+zone+"/Muscle"].outlet.state.T * self.controlParameters["Muscle temperature distribution"][zone]
    for zone in self.controlParameters["Garment temperature distribution"]:
        TGarmAv += self.components[zone+"/Skin -> Garment"].outlet.state.T * self.controlParameters["Garment temperature distribution"][zone]

    # Calculate total heat flows for storage
    QConvection = 0
    QRadiation = 0
    QConductionGarment = 0
    QLCG = 0

    for c in self.conductionComponents:
        QConductionGarment += self.components[c].conductionHeatTransfer
    for c in self.convectionComponents:
        QConvection += -self.components[c].convectionHeatTransfer
    for c in self.radiationComponents:
        QRadiation += self.components[c].radiationHeatTransfer
    for c in self.LCGComponents:
        QLCG += -self.components[c].convectionHeatTransfer

    # Enact heating
    for zone in muscleWorks:
        self.components["Metabolism @ "+zone+"/Muscle"].QIn = muscleWorks[zone]
    
    # Enact muscle bloodflow
    for zone in muscleBloodflows:
        self.components[zone+"/Muscle -> Blood"].mDot = muscleBloodflows[zone]
    
    # Enact skin bloodflow
    for zone in skinBloodFlows:
        self.components[zone+"/Skin -> Blood"].mDot = skinBloodFlows[zone]

    # Enact sweat rate, respecting max value
    QSweat = 0
    QDiffusion = 0
    for zone in targetSweatRates:
        EMaxZone = self.components[zone+" ~ Sweat"].dilutesModel["H2O"]["EMax"]
        
        sweatPower = min(targetSweatRates[zone], EMaxZone)
        QSweat += abs(sweatPower)
        sweatPower +=  diffusionRates[zone]
        QDiffusion += abs(diffusionRates[zone])
        
        sweatRate = sweatPower / self.components[zone+" ~ Sweat"].dilutesModel["H2O"]["enthalpy"]
        sweatRate = max(0, sweatRate)
        self.components[zone+" ~ Sweat"].dilutesModel["H2O"]["mDot"] = sweatRate
    
    # Enact conductance and flow in respiratory nodes
    QLatentBreath = 0
    for node in self.controlParameters["Respiratory heat ratios"]:
        self.components[node+" ~ Breath"].convectionModel["Conductance"] = convectionConductanceRespiratory[node]
        self.components[node+" ~ Breath"].dilutesModel["H2O"]["Q"] = max(evaporationPowerRespiratory[node], 0)
        QLatentBreath += abs(max(evaporationPowerRespiratory[node], 0))

    self.components["Head/Core ~ Breath"].dilutesModel["CO2"]["mDot"] = mDotCO2

    energyStorageRate = 0
    if len(self.outputHistory["Total storage (J)"]) > 0:
        energyStorageRate = (energyStorage - self.outputHistory["Total storage (J)"][-1])/dt

    # Record data in the outputHistory
    self.outputHistory["Time (min)"].append(time/60)
    self.outputHistory["Temp head core (K)"].append(self.components["Metabolism @ Head/Core"].outlet.state.T)
    self.outputHistory["Temp skin average (K)"].append(TSkinAv)
    self.outputHistory["Temp muscle average (K)"].append(TMuscleAv)
    self.outputHistory["Temp undergarment average (K)"].append(TGarmAv)
    self.outputHistory["Q sensible (W)"].append(QRadiation + QConvection + QConductionGarment)
    self.outputHistory["Q evap (W)"].append(QSweat)
    self.outputHistory["Q latent (W)"].append(QSweat + QDiffusion + QLatentBreath)
    self.outputHistory["Q LCG (W)"].append(QLCG)
    self.outputHistory["Q shiver (W)"].append(shiveringHeat)
    self.outputHistory["Heat storage rate (W)"].append(energyStorageRate)
    self.outputHistory["Total storage (J)"].append(energyStorage)
    self.outputHistory["Total storage (BTU)"].append(energyStorage * 0.000947817)

def initialiseBody(solver):

    for layer in solver.controlParameters["Set points"]:
        for zone in solver.controlParameters["Set points"][layer]:
            solver.components["Metabolism @ "+zone+"/"+layer].outlet.state.T = solver.controlParameters["Set points"][layer][zone]

def initialiseSpacesuit(solver, TEnvironment=250, TAir=290, TLCG=290):
    THalf = (TEnvironment + TAir)/2
    TInner = min(TAir - 20, THalf)
    TOuter = max(TEnvironment+20, THalf)
    if TLCG != TAir:
        TSkin = min(TLCG, 309)
    else:
        TSkin = min(TAir+5, 309)
    TLCGOut = TSkin*0.8 + TLCG*0.2
    for zone in zones:
        solver.components[zone+"/Suit"].inlet.state.T = TInner
        solver.components[zone+"/Suit"].outlet.state.T = TOuter
        solver.components["Metabolism @ "+zone+"/Skin"].outlet.state.T = TSkin
        try:
            solver.components[zone+"/LCG"].fluidOutlet.state.T = TLCGOut
        except:
            pass



solver = bjls.TimestepSolver(maxTimestep=60, precision=10)
solver.workHeat = stepWork
solver.referenceValues = {}
solver.outputHistory = {"Time (min)": [], "Temp head core (K)": [], 
                        "Temp skin average (K)": [], "Temp muscle average (K)": [],
                        "Temp undergarment average (K)": [],
                        "Q sensible (W)": [], "Q evap (W)": [], "Q latent (W)": [], "Q LCG (W)":[],
                        "Heat storage rate (W)": [], "Q shiver (W)": [], 
                        "Total storage (J)": [], "Total storage (BTU)": []}
solver.systemController = systemController.__get__(solver)
solver.controlParameters = {
    "Bloodflow per W muscles": 4.299e-4 * 0.554, #kg/s of additional bloodflow per W of mechanical work in muscles
    "Vasodilation coefficient": 2.268e-4 * 183.5, #kg/(s-K) of additional bloodflow per K error in head core temp
    "Vasoconstriction coefficient": 9.99, # 1/K in vasoconstriction term
    "Shivering coefficient": 11.604,    # W/K^2, coefficient in front of shivering control sum
    "Work distribution": {
        "Head":   0.000, "Trunk":  0.300,
        "Arm L":  0.040, "Arm R":  0.040,
        "Hand L": 0.005, "Hand R": 0.005,
        "Leg L":  0.300, "Leg R":  0.300,
        "Foot L": 0.005, "Foot R": 0.005
    },
    "Vasodilation distribution": {
        "Head":   0.1320, "Trunk":  0.3220,
        "Arm L":  0.0475, "Arm R":  0.0475,
        "Hand L": 0.0610, "Hand R": 0.0610,
        "Leg L":  0.1150, "Leg R":  0.1150,
        "Foot L": 0.0500, "Foot R": 0.0500
    },
    "Vasoconstriction distribution": {
        "Head":   0.050, "Trunk":  0.150,
        "Arm L":  0.025, "Arm R":  0.025,
        "Hand L": 0.175, "Hand R": 0.175,
        "Leg L":  0.025, "Leg R":  0.025,
        "Foot L": 0.175, "Foot R": 0.175
    },
    "Skin compartment distribution": {
        "Head":   0.08270, "Trunk":  0.58700,
        "Arm L":  0.04110, "Arm R":  0.04110,
        "Hand L": 0.01110, "Hand R": 0.01110,
        "Leg L":  0.09300, "Leg R":  0.09300,
        "Foot L": 0.01995, "Foot R": 0.01995
    },
    "Shivering distribution": {
        "Head":   0.02300, "Trunk":  0.94800,
        "Arm L":  0.00265, "Arm R":  0.00254,
        "Hand L": 0.00115, "Hand R": 0.00115,
        "Leg L":  0.00950, "Leg R":  0.00950,
        "Foot L": 0.01200, "Foot R": 0.01200
    },
    "Recalculate max sweat rate": False,
    "Respiratory quotient": 0.82,
    "Sweat distribution": {
        "Head":   0.0810, "Trunk":  0.4820,
        "Arm L":  0.0765, "Arm R":  0.0765,
        "Hand L": 0.0155, "Hand R": 0.0155,
        "Leg L":  0.1090, "Leg R":  0.1090,
        "Foot L": 0.0175, "Foot R": 0.0175
    },
    "Skin temperature distribution": {
        "Head":   0.07000, "Trunk":  0.36020,
        "Arm L":  0.06075, "Arm R":  0.06075,
        "Hand L": 0.02500, "Hand R": 0.02500,
        "Leg L":  0.15870, "Leg R":  0.15870,
        "Foot L": 0.03430, "Foot R": 0.03430
    },
    "Muscle temperature distribution": {
        "Head":   0.02325, "Trunk":  0.54900,
        "Arm L":  0.05270, "Arm R":  0.05270,
        "Hand L": 0.00115, "Hand R": 0.00115,
        "Leg L":  0.15920, "Leg R":  0.15920,
        "Foot L": 0.00115, "Foot R": 0.00115
    },
    "Garment temperature distribution": {
        "Trunk":  0.33170,
        "Arm L":  0.10400, "Arm R":  0.10400,
        "Leg L":  0.23015, "Leg R":  0.23015,
    },
    "Respiratory heat ratios":{
        "Head/Core":0.3855, "Head/Muscle":0.0860, 
        "Head/Fat": 0.0287, "Trunk/Core": 0.0238, "Trunk/Muscle":0.2615
    },
    "Respiratory temperature distribution":{
        "Head/Core":0.3850, "Head/Muscle":0.0860, 
        "Head/Fat": 0.0287, "Trunk/Core": 0.2380, "Trunk/Muscle":0.2615
    },
    "Suit outer absorbivity": 0.5,
    "Suit outer emissivity": 0.5,
    "Set points": {
        "Core": {
            "Head":   310.15, "Trunk":  310.26,
            "Arm L":  308.76, "Arm R":  308.76,
            "Hand L": 308.65, "Hand R": 308.65,
            "Leg L":  309.59, "Leg R":  309.59,
            "Foot L": 308.59, "Foot R": 308.59
        },
        "Muscle": {
            "Head":   309.59, "Trunk":  309.98,
            "Arm L":  308.21, "Arm R":  308.21,
            "Hand L": 308.54, "Hand R": 308.54,
            "Leg L":  308.98, "Leg R":  308.98,
            "Foot L": 308.43, "Foot R": 308.43
        },
        "Fat": {
            "Head":   309.26, "Trunk":  308.65,
            "Arm L":  307.65, "Arm R":  307.65,
            "Hand L": 308.48, "Hand R": 308.48,
            "Leg L":  308.21, "Leg R":  308.21,
            "Foot L": 308.54, "Foot R": 308.54
        },
        "Skin": {
            "Head":   309.04, "Trunk":  307.82,
            "Arm L":  307.43, "Arm R":  307.43,
            "Hand L": 308.43, "Hand R": 308.43,
            "Leg L":  307.87, "Leg R":  307.87,
            "Foot L": 307.87, "Foot R": 307.87
        }
    
    } 
}
#makeSolverBodyParams(165, 48, useNavyBodyFat=True, male=False, waistCircum=60, neckCircum=30, hipCircum=93)
thermalNodeValues = makeSolverBodyParams(175, 80)

solver.thermalNodeValues = thermalNodeValues

makeBodyThermalConnections(solver)

"""makeShirtsleveesThermalConnections(solver, 
                                    cabinPressurebar=1, cabinAirflowms=1.0, gravityRatio=1, TCabin=290, approxTGarment=300,
                                    cabinWallTemperature=285, cabinAirflowkgs=0.05, cabinAirflowInletHumidity=0.9, cabinAirflowCO2=400e-6, 
                                    cabinWallArea=5, cabinVolume=10, undergarmentEmissivity=0.8,
                                   undergarmentConductivities={
                                        "Head":   0, "Trunk":  0.05,
                                        "Arm L":  0.05, "Arm R":  0.05,
                                        "Hand L": 0, "Hand R": 0,
                                        "Leg L":  0.05, "Leg R":  0.05,
                                        "Foot L": 0, "Foot R": 0
                                    }, 
                                    undergarmentMasses={
                                        "Head":   0, "Trunk":  0.1,
                                        "Arm L":  0.1, "Arm R":  0.1,
                                        "Hand L": 0, "Hand R": 0,
                                        "Leg L":  0.1, "Leg R":  0.1,
                                        "Foot L": 0, "Foot R": 0
                                    },)"""
calculateSpacesuitInternalThermalConnections(solver,
                                     suitPressureBar=1, TInlet=290, approxTGarment=295, 
                                     inletHumidity=0.2, inletCO2=400e-6, 
                                     undergarmentEmissivity=0.8, suitInnerEmissivity=0.8, suitOuterEmissivity=0.5, 
                                     suitOuterAbsorbivity=0.5, suitSHC=1000, 
                                     LCGTInlet=290, useLCG=True,
                                     undergarmentConductivities={
                                        "Head":   0, "Trunk":  0.05,
                                        "Arm L":  0.05, "Arm R":  0.05,
                                        "Hand L": 0, "Hand R": 0,
                                        "Leg L":  0.05, "Leg R":  0.05,
                                        "Foot L": 0, "Foot R": 0}, 
                                    undergarmentMasses={
                                        "Head":   0, "Trunk":  0.1,
                                        "Arm L":  0.1, "Arm R":  0.1,
                                        "Hand L": 0, "Hand R": 0,
                                        "Leg L":  0.1, "Leg R":  0.1,
                                        "Foot L": 0, "Foot R": 0},
                                    LCGFlowRates={
                                        "Head":   0, "Trunk":  0.01,
                                        "Arm L":  0.003, "Arm R":  0.003,
                                        "Hand L": 0, "Hand R": 0,
                                        "Leg L":  0.005, "Leg R":  0.005,
                                        "Foot L": 0, "Foot R": 0},
                                    suitMasses={
                                        "Head":   4, "Trunk":  100,
                                        "Arm L":  3, "Arm R":  3,
                                        "Hand L": 1, "Hand R": 1,
                                        "Leg L":  4, "Leg R":  4,
                                        "Foot L": 1, "Foot R": 1}, 
                                    suitConductivities={
                                        "Head":   0.05, "Trunk":  0.03,
                                        "Arm L":  0.04, "Arm R":  0.04,
                                        "Hand L": 0.1, "Hand R": 0.1,
                                        "Leg L":  0.03, "Leg R":  0.03,
                                        "Foot L": 0.1, "Foot R": 0.1}, 
                                    suitAreas={
                                        "Head":   0.04, "Trunk":  0.1,
                                        "Arm L":  0.02, "Arm R":  0.03,
                                        "Hand L": 0.01, "Hand R": 0.01,
                                        "Leg L":  0.05, "Leg R":  0.05,
                                        "Foot L": 0.01, "Foot R": 0.01},
                                    zoneAirflows={
                                        "Air inlet": 0.04,
                                        "Head":   0.04, "Trunk":  0.04,
                                        "Arm L":  0.01, "Arm R":  0.01,
                                        "Hand L": 0.01, "Hand R": 0.01,
                                        "Leg L":  0.01, "Leg R":  0.01,
                                        "Foot L": 0.01, "Foot R": 0.01}, 
                                    zoneAirConnectivity={
                                        "Head":   {"From": "In", "Fraction": 1}, "Trunk":  {"From": "Head", "Fraction": 1},
                                        "Arm L":  {"From": "Trunk", "Fraction": 0.25}, "Arm R":  {"From": "Trunk", "Fraction": 0.25},
                                        "Hand L": {"From": "Arm L", "Fraction": 1}, "Hand R": {"From": "Arm R", "Fraction": 1},
                                        "Leg L":  {"From": "Trunk", "Fraction": 0.25}, "Leg R":  {"From": "Trunk", "Fraction": 0.25},
                                        "Foot L": {"From": "Leg L", "Fraction": 1}, "Foot R": {"From": "Leg R", "Fraction": 1}}, 
                                    zoneAirVolumes={
                                        "Head":   0.08, "Trunk":  0.1,
                                        "Arm L":  0.01, "Arm R":  0.01,
                                        "Hand L": 0.01, "Hand R": 0.01,
                                        "Leg L":  0.01, "Leg R":  0.01,
                                        "Foot L": 0.01, "Foot R": 0.01},
                                    useFixedSource=True
                                        )


solver.initialiseSolver(Tref=295, mDot=None, dilutes={"H2O":0, "CO2":400e-6})

initialiseBody(solver)
initialiseSpacesuit(solver)


print(len(solver.components.values()))

for c in solver.components:
    if "Air" in c:
        try:print(c, "...", solver.components[c].fluidInlet.state.ID, "\t", solver.components[c].fluidOutlet.state.ID, solver.components[c].fluidOutlet.state.mDot)
        except:
            try:
                print(c, "...", solver.components[c].fluidOutlet0.state.ID, "\t", solver.components[c].fluidOutlet0.state.ID)
            except:pass

solver.mode = "RK4ERROR"
start = time.time()
solver.runSolver(0, 1*60*60, 1)
print(time.time() - start)

import matplotlib.pyplot as plt


fig, ax1 = plt.subplots()

ax2 = ax1.twinx()


#plt.plot(solver.time, [i-273.15 for i in solver.components["Metabolism @ Head/Core"].outlet.state.THistory], label="Head core", c="red")
#plt.plot(solver.time, [i-273.15 for i in solver.components["Metabolism @ Leg L/Core"].outlet.state.THistory], label="Leg core", c="blue")
#plt.plot(solver.time, [i-273.15 for i in solver.components["Metabolism @ Head/Skin"].outlet.state.THistory], label="Head skin", c="red", linestyle="--")
#plt.plot(solver.time, [i-273.15 for i in solver.components["Metabolism @ Leg L/Skin"].outlet.state.THistory], label="Leg skin", c="blue", linestyle="--")
#plt.plot(solver.time, [i-273 for i in solver.components["Cabin air inlet"].fluidOutlet.state.THistory], label="Cabin air in")
#plt.plot(solver.time, [i-273 for i in solver.components["Cabin air outlet"].fluidInlet.state.THistory], label="Cabin air out")
#ax1.plot(solver.time, solver.components["Cabin air inlet"].fluidOutlet.state.dilutesConcentrationHistory["H2O"], label="inlet", c="blue", linestyle=":")
#ax1.plot(solver.time, solver.components["Cabin air outlet"].fluidInlet.state.dilutesConcentrationHistory["H2O"], label="outlet", c="red", linestyle=":")
#ax2.plot(solver.time, solver.components["Cabin air inlet"].fluidOutlet.state.THistory, label="inlet", c="blue")
#ax2.plot(solver.time, solver.components["Cabin air outlet"].fluidInlet.state.THistory, label="outlet", c="red")
#plt.plot(solver.time, solver.components["Cabin air outlet"].fluidInlet.state.THistory)
#ax2.set_ylim(280, 300)

ax1.plot(solver.outputHistory["Time (min)"], solver.outputHistory["Q sensible (W)"], label="Q sensible")
ax1.plot(solver.outputHistory["Time (min)"], solver.outputHistory["Q LCG (W)"], label="Q LCG", linestyle="--")
ax1.plot(solver.outputHistory["Time (min)"], solver.outputHistory["Heat storage rate (W)"], label="Q stor", linestyle=":")
ax1.plot(solver.outputHistory["Time (min)"], solver.outputHistory["Q latent (W)"], label="Q latent")

#ax1.plot(solver.time, [i for i in solver.components["Metabolism @ Trunk/Core"].outlet.state.THistory], label="Trunk core", c="red")
#ax1.plot(solver.time, [i for i in solver.components["Metabolism @ Trunk/Muscle"].outlet.state.THistory], label="Trunk Muscle", c="blue")
#ax1.plot(solver.time, [i for i in solver.components["Metabolism @ Trunk/Skin"].outlet.state.THistory], label="Trunk Skin", c="red")
#ax1.plot(solver.time, [i for i in solver.components["Trunk/Suit"].inlet.state.THistory], label="Trunk suit inner", c="red", linestyle=":")
#ax1.plot(solver.time, [i for i in solver.components["Trunk/Suit"].outlet.state.THistory], label="Trunk suit outer", c="red", linestyle="--")

#ax1.plot(solver.time, [i for i in solver.components["Trunk/LCG"].fluidOutlet.state.THistory], label="Trunk LCG")
#ax1.plot(solver.time, [i for i in solver.components["Arm L/LCG"].fluidOutlet.state.THistory], label="Arm LCG")
#ax1.hlines([37+273.15], [0], [solver.time[-1]], linestyles="dashed", colors=["k"], )


#ax2.plot(solver.outputHistory["Time (min)"], solver.components["Cabin air outlet"].fluidInlet.state.dilutesConcentrationHistory["H2O"], linestyle=":")
#ax2.plot(solver.outputHistory["Time (min)"], [i for i in solver.outputHistory["Q sensible (W)"]], label="Q sensible", linestyle=":")
#ax2.plot(solver.outputHistory["Time (min)"], [i for i in solver.outputHistory["Q latent (W)"]], label="Q latent", linestyle=":")
ax2.plot(solver.outputHistory["Time (min)"], solver.outputHistory["Total storage (J)"], label="E store", linestyle=":", c="k")
#ax2.plot(solver.time, solver.timesteps, linestyle=":")

print("Stored E: ", str(solver.outputHistory["Total storage (J)"][-1]))
print("skin temp: ", str(solver.outputHistory["Temp skin average (K)"][-1]))
print("latent heat: ", str(solver.outputHistory["Q latent (W)"][-1]))
print("sensible heat: ", str(solver.outputHistory["Q sensible (W)"][-1]))
print("lcg outlet: ", str(solver.components["LCG Out"].fluidInlet.state.T))


ax1.legend()


plt.xlabel("Time (s)")
#plt.ylabel("Node temperature (C)")

plt.title("Head layer temperatures, bloodflow and thermoreg")

#plt.yticks(np.arange(310.15, 330.15), [str(int(i)) for i in np.arange(310.15-273.15, 330.15-273.15)])

plt.show()