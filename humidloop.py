from components import *

greenhouse = CycleSolver("HumidAir")

greenhouse.components["Plants"] = FluidSource_TP(1/3600, 1E5, 300, RH=1.0)
greenhouse.components["V1"] = PressureDrop()
greenhouse.components["V2"] = Heater(mode="T")
greenhouse.components["V3"] = FlowJunction()
greenhouse.components["F1"] = Compressor()
greenhouse.components["HX1"] = Dehumidifier()
greenhouse.components["HXP"] = PressureDrop()
greenhouse.components["Drain"] = FluidSink()

greenhouse.connectPorts(greenhouse.components["HXP"].outlet, greenhouse.components["V1"].inlet, 
                        stateType=HumidAirState)
greenhouse.connectPorts(greenhouse.components["V1"].outlet, greenhouse.components["V2"].inlet, 
                        stateType=HumidAirState)
greenhouse.connectPorts(greenhouse.components["V2"].outlet, greenhouse.components["V3"].inlet1, 
                        stateType=HumidAirState)
greenhouse.connectPorts(greenhouse.components["Plants"].outlet, greenhouse.components["V3"].inlet2, 
                        stateType=HumidAirState)
greenhouse.connectPorts(greenhouse.components["V3"].outlet, greenhouse.components["F1"].inlet, 
                        stateType=HumidAirState)
greenhouse.connectPorts(greenhouse.components["F1"].outlet, greenhouse.components["HX1"].inlet, 
                        stateType=HumidAirState)
greenhouse.connectPorts(greenhouse.components["HX1"].airOutlet, greenhouse.components["HXP"].inlet, 
                        stateType=HumidAirState)
greenhouse.connectPorts(greenhouse.components["F1"].outlet, greenhouse.components["HX1"].inlet, 
                        stateType=HumidAirState)
greenhouse.connectPorts(greenhouse.components["HX1"].waterOutlet, greenhouse.components["Drain"].inlet, 
                        stateType=ConstantFluidState, fluidOverride="Water")


