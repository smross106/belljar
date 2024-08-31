from components import *
import CoolProp.CoolProp as CP
from matplotlib import ticker
from scipy import optimize
import pyforfluids as pff
from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL, FlashVLN
from thermo.interaction_parameters import IPDB

from CoolProp.HumidAirProp import HAPropsSI

"""
T 288
P 1e5
Hmass 290078.2882959461
Smass 6852.412494872328
Q -1.0
Dmass 1.205301577686857

T, p: all variables
p, h: none
p, s
"""


print(HAPropsSI("Hha", "P", 1e5, "T", 298, "RH", 0.7))
print(HAPropsSI("Hha", "P", 1e5, "T", 298, "RH", 1.0))

print(HAPropsSI("Sha", "P", 1e5, "T", 298, "RH", 0.7))
print(HAPropsSI("Sha", "P", 1e5, "T", 298, "RH", 1.0))

print(HAPropsSI("W", "P", 1e5, "T", 298, "RH", 0.7))
print(HAPropsSI("W", "P", 1e5, "T", 298, "RH", 1.0))

print(HAPropsSI("RH", "Hha", 59916, "P", 1e5, "W", 0.014))

exit()

CP.apply_simple_mixing_rule("Oxygen","Water","linear")
args = ["T", "P", "Hmass", "Smass", "Dmass"]
vals = [288, 1e5, 470687.37620*2, 6694.12*1.09, 1.2053]
fluid = "Oxygen[0.5]&Water[0.5]"
comp = {"Oxygen":0.5, "CarbonMonoxide":0.5}
#fluids = ["Oxygen", "CarbonMonoxide"]
#proportions = [0.5, 0.5]

#print(PropsSI("P", "Hmass", 290078, "Smass", 6852, fluid))
#print(PropsSI("Smass", "P", 1e5, "T", 288, fluid))

def resid_hs(input_args, *args):
    (P, T) = input_args
    (fluid, HTarget, STarget) = args
    H = PropsSI("Hmass", "P", P, "T", T, fluid)
    S = PropsSI("Smass", "P", P, "T", T, fluid)
    HResid = 100 * ((H/HTarget)-1)
    SResid = 100 * ((S/STarget)-1)
    resid = HResid**2 + SResid**2
    return(resid)

def resid_ps(input_args, *args):
    (T) = input_args
    (fluid, P, STarget) = args
    try:
        S = PropsSI("Smass", "P", P, "T", T, fluid)
        SResid = 100 * ((S/STarget)-1)
        resid = SResid**2
    except:
        print("!")
        resid = 1e9 + T
    return(resid)

def resid_ph(input_args, *args):
    (T) = input_args
    (fluid, P, HTarget) = args
    try:
        H = PropsSI("Hmass", "P", P, "T", T, fluid)
        HResid = 100 * ((H/HTarget)-1)
        resid = HResid**2
    except:
        print("!")
        resid = 1e9 + T
    return(resid)


def resid_pQ(input_args, *args):
    (T) = input_args
    (fluid, P, QTarget) = args
    try:
        Q = PropsSI("Q", "P", P, "T", T, fluid)
        QResid = 100 * abs(Q-QTarget)
        resid = QResid**2
    except:
        print("!")
        resid = 1e9 + T
    return(resid)


constants, properties = ChemicalConstantsPackage.from_IDs(['methane', 'carbon monoxide'])
kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}

gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)

zs = [0.5, 0.5]

PT = flasher.flash(VF=0.5, P=5e5, zs=zs)

vals = [PT.T, PT.P, PT.H_mass(), PT.S_mass(), PT.rho_mass(), PT.VF]

print(vals)

exit()
"""out_hs = optimize.minimize(resid_hs, 
        [1e5, 288], 
        args=(fluid, vals[2], vals[3]),
        bounds=[(1, 500e5), (100, 1000)],
        method="Nelder-Mead")"""
T = []
Q = []
for q in range(0, 21):
    print(q/20)
#print(PropsSI("T", "P", 1e5, "Q", 0.08, fluid))
    Q.append(q/20)
    #T.append(PropsSI("T", "P", 1e5, "Q", q/20, fluid))
    out_pQ = optimize.minimize(resid_pQ, 
        [250], 
        args=(fluid, vals[2], q/20),
        bounds=[(100, 1000)],
        method="Nelder-Mead")
    T.append(out_pQ.x[0])

plt.plot(Q, T)
plt.show()

exit()

out = optimize.minimize(resid_ps,
                        [288],
                        args = (fluid, vals[1], vals[3]),
                        bounds=[(100, 1000)],
                        method="Nelder-Mead")
P=1e5
T = out.x[0]
print(P, T)
print(PropsSI("Smass", "P", P, "T", T, fluid), vals[3])
#print(PropsSI("Smass", "P", P, "T", T, fluid), vals[3])