import numpy as np
import copy 

ni = 60
nj = 20

xlow = np.linspace(0, 10, ni)
xhigh = np.linspace(0, 10, ni)
ylow = np.zeros(ni)
ylow = [0.02 * x for x in xlow[0:30]]+[0.02 * x for x in xlow[0:30]][::-1]
yhigh = [3-0.02 * x for x in xlow[0:30]]+[3-0.02 * x for x in xlow[0:30]][::-1]


rgas = 287
gamma = 1.4
pstagin = 100000.   
tstagin = 300.   
alpha1 = 00.0   
pdown = 85000.
cfl = 0.5
smooth_fac_in = 1
nsteps = 3000
conlim_in = 0.5

emax       = 1000000.0
eavg       = emax
cp         = rgas*gamma/(gamma-1.0)
cv         = cp/gamma
fga        = (gamma - 1.0)/gamma
smooth_fac = smooth_fac_in*cfl
conlim     = conlim_in*cfl
alpha1     = alpha1*3.14159/180.0

ncells    =  ni * nj
jmid      =  (1 + nj)/2
roin      =  pstagin/rgas/tstagin
ref_ro    =  (pstagin-pdown)/rgas/tstagin
ref_t     =  tstagin*(pdown/pstagin)**fga
ref_v     =  np.sqrt(2*cp*(tstagin-ref_t))
ref_rovx  =  roin*ref_v
ref_rovy  =  roin*ref_v
ref_roe   =  roin*cv*(tstagin-ref_t)

x = np.zeros((ni, nj))
y = np.zeros((ni, nj))
area = np.zeros((ni, nj))
dli = np.zeros((ni, nj))
dlj = np.zeros((ni, nj))
dljx = np.zeros((ni, nj)) 
dljy = np.zeros((ni, nj)) 
dlix = np.zeros((ni, nj))
dliy = np.zeros((ni, nj))
ro_inc = np.zeros((ni, nj))
rovx_inc = np.zeros((ni, nj))
rovy_inc = np.zeros((ni, nj))
roe_inc = np.zeros((ni, nj))
ro_start = np.zeros((ni, nj))
rovx_start = np.zeros((ni, nj))
rovy_start = np.zeros((ni, nj))
roe_start = np.zeros((ni, nj))
ro = np.zeros((ni, nj))
rovx = np.zeros((ni, nj))
rovy = np.zeros((ni, nj))
roe = np.zeros((ni, nj))
vx = np.zeros((ni, nj))
vy = np.zeros((ni, nj))
tstatic = np.zeros((ni, nj))
p = np.zeros((ni, nj))
hstag = np.zeros((ni, nj))

fluxi_mass = np.zeros((ni, nj))
fluxj_mass = np.zeros((ni, nj))
fluxi_xmom = np.zeros((ni, nj))
fluxj_xmom = np.zeros((ni, nj))
fluxi_ymom = np.zeros((ni, nj))
fluxj_ymom = np.zeros((ni, nj))
fluxi_enth = np.zeros((ni, nj))
fluxj_enth = np.zeros((ni, nj))
previous = np.zeros((ni, nj))
delro = np.zeros((ni, nj))
delrovx = np.zeros((ni, nj))
delrovy = np.zeros((ni, nj))
delroe = np.zeros((ni, nj))
store = np.zeros((ni, nj))

pin = np.zeros(nj)
roinlet = np.zeros(nj)
flow = np.zeros(ni)

dmin = 9999999

# Generate grid

for i in range(ni):
    for j in range(nj):
        x[i, j] = xlow[i]+(j/nj)*(xhigh[i] - xlow[i])
        y[i, j] = ylow[i]+(j/nj)*(yhigh[i] - ylow[i])

for i in range(ni-1):
    for j in range(nj-1):
        A1 = x[i, j+1]-x[i+1, j]
        B1 = y[i, j+1]-y[i+1, j]
        A2 = x[i+1, j+1]-x[i, j]
        B2 = y[i+1, j+1]-y[i, j]
        area[i, j] = 0.5 * abs((A2*B1) - (A1*B2))

for i in range(ni):
    for j in range(nj-2):
        dlix[i, j] = y[i, j+1] - y[i, j]
        dliy[i, j] = - (x[i, j+1] - x[i, j])
        ilength = np.sqrt(dlix[i, j]**2 + dliy[i, j]**2)
        dmin = min(ilength, dmin)

for i in range(ni-2):
    for j in range(nj):
        dljx[i, j] = -(y[i+1, j] - y[i, j])
        dljy[i, j] = x[i+1, j] - x[i, j]
        jlength = np.sqrt(dljx[i, j]**2 + dljy[i, j]**2)
        dmin = min(jlength, dmin)


print(np.max(area), np.min(area))


# Guess the flow
jmid   = int(nj/2)
tdown  = tstagin*(pdown/pstagin)**fga
vdown  = np.sqrt(2*cp*(tstagin - tdown))
rodown = pdown/rgas/tdown
pinlet = 55000.
tinlet  = tstagin*(pinlet/pstagin)**fga
vinlet  = np.sqrt(2*cp*(tstagin - tinlet))
roin    = pinlet/rgas/tinlet

for j in range(nj):
    for i in range(ni-1):
        dx = 1e-7 + x[i+1, jmid] - x[i, jmid]
        dy = 1e-7 + y[i+1, jmid] - y[i, jmid]
        ds = np.sqrt(dx*dx + dy*dy)
        vlocal = vinlet  + (vdown-vinlet) * (i/ni)
        rolocal = roin + (rodown - roin) * (i/ni)
        tlocal = tinlet + (tdown - tinlet) * (i/ni)

        xvel = vlocal * (dx/ds)
        yvel = vlocal * (dy/ds)

        rovx[i, j] = rolocal * xvel 
        rovy[i, j] = rolocal * yvel
        ro[i, j] = rolocal
        roe[i, j] = rolocal * (cv*tlocal + 0.5*vlocal*vlocal)

    rovx[ni-1, j] = rovx[ni-2, j]
    rovy[ni-1, j] = rovy[ni-2, j]
    ro[ni-1, j] = ro[ni-2, j]
    roe[ni-1, j] = roe[ni-2, j]

# set timestep
astag  = np.sqrt(gamma*rgas*tstagin)
umax   = astag
deltat = cfl*dmin/(umax+astag)

def sum_fluxes(iflux, jflux, delprop, prop_inc):
    previous = copy.copy(delprop)

    for i in range(ni-1):
        for j in range(nj-1):
            delprop[i,j] = (deltat/area[i,j])*(iflux[i,j]-iflux[i+1,j]+jflux[i,j]-jflux[i,j+1])

    for i in range(1, ni-2):
        for j in range(1, nj-2):
            prop_inc[i, j] = 0.25 * (delprop[i, j] + delprop[i-1, j] + delprop[i, j-1] + delprop[i-1, j-1])
    
    for i in range(2, ni-1):
        prop_inc[i, 0] = 0.5 * (delprop[i, 0] + delprop[i-1, 0])
        prop_inc[i, nj-1] = 0.5 * (delprop[i, nj-1] + delprop[i-1, nj-1])
    
    for j in range(2, nj-1):
        prop_inc[ni-1, j] = 0.5 * (delprop[ni-2, j-1] + delprop[ni-2, j])
        prop_inc[0, j] = 0.5 * (delprop[0, j-1] + delprop[0, j])
    
    prop_inc[0, 0] = delprop[0, 0]
    prop_inc[0, nj-1] = delprop[0, nj-1]
    prop_inc[ni-1, 0] = delprop[ni-1, 0]
    prop_inc[ni-1, nj-1] = delprop[ni-1, nj-1]

def smooth(prop, maxprop):
    sf = smooth_fac
    sfm1 = 1 - sf

    for i in range(ni):
        ip1 = min(i+1, ni-1)
        im1 = max(i-1, 0)
        for j in range(1, nj-1):
            avg = 0.25 * (prop[ip1, j]+prop[im1, j]+prop[i, j-1] + prop[i, j+1])
        
            store[i,j] = min(sfm1*prop[i,j] + sf * avg, maxprop)
        
        avg0 = (prop[im1, 0]+prop[ip1, 0]+2*prop[i, 1]-prop[i, 2])/3
        avgnj = (prop[im1, nj-1]+prop[ip1,nj-1]+prop[i,nj-2])/3
        store[i,0] = sfm1*prop[i,0] + sf*avg0
        store[i, nj-1] = sfm1*prop[i,nj-1] + sf*avgnj
    
    prop = copy.copy(store)


for nstep in range(100):

    for i in range(ni): 
        for j in range(nj):
            ro_start[i, j] = ro[i, j]
            rovx_start[i, j] = rovx[i, j]
            rovy_start[i, j] = rovy[i, j]
            roe_start[i, j] = roe[i, j]
    
    vx = rovx/ro
    vy = rovy/ro
    tstatic = tstagin - 0.5 * np.sqrt(vx**2 + vy**2)/cp
    p = ro * rgas * tstatic
    hstag = (roe + p) / ro

    rfin     = 0.25
    rfin1    = 1.0-rfin
    rostagin = pstagin/(rgas*tstagin)
    gm1      = gamma - 1.0

    for j in range(nj):
        if nstep == 0:
            roinlet[j] = ro[0, j]
        else:
            roinlet[j] = rfin*ro[0,j] + rfin1*roinlet[j]
        
        if roinlet[j] > 0.9999 * rostagin:
            roinlet[j] = 0.9999 * rostagin
        
        tstat = tstagin * ((roinlet[j]/rostagin)**(gm1))
        vel = np.sqrt(2*cp*(tstagin-tstat))
        eke = cv*tstat + 0.5*vel**2


        vx[0, j] = vel * np.cos(alpha1)
        vy[0, j] = vel * np.sin(alpha1)
        p[0, j] = roinlet[j]*rgas*tstat 
        rovx[0, j] = roinlet[j]*vx[0, j]
        rovy[0, j] = roinlet[j]*vy[0, j]
        roe[0, j] = roinlet[j] * eke

        hstag[0, j] = (roe[0, j] + p[0, j])/roinlet[j]

        p[ni-1, j] = pdown

    # Set fluxes
    for i in range(ni):
        flow[i] = 0
        for j in range(nj-1):
            fluxi_mass[i, j] = 0.5* (((rovx[i, j]+rovx[i, j+1])*dlix[i,j]) + 
                                     ((rovy[i, j]+rovy[i, j+1])*dliy[i,j]))
            flow[i] += fluxi_mass[i, j]
    
    for i in range(ni-1):
        for j in range(1, nj-2):
            fluxj_mass[i, j] = 0.5* (((rovx[i, j]+rovx[i, j+1])*dljx[i,j]) + 
                                     ((rovy[i, j]+rovy[i, j+1])*dljy[i,j]))
            flow[j] += fluxj_mass[i, j]
    
    for i in range(ni-1):
        fluxj_mass[i, 0] = 0
        fluxj_mass[i, nj-1] = 0
    
    for i in range(ni):
        for j in range(nj-1):
            fluxi_xmom[i, j] = 0.5 * ( fluxi_mass[i, j] * (vx[i, j] + vx[i, j+1]) + 
                                      (p[i, j] + p[i, j+1]) * dlix[i, j])
    for i in range(ni-1):
        for j in range(nj):
            fluxj_xmom[i, j] = 0.5 * ( fluxj_mass[i, j] * (vx[i, j] + vx[i+1, j]) + 
                                      (p[i, j] + p[i+1, j]) * dljx[i, j])
    
    for i in range(ni):
        for j in range(nj-1):
            fluxi_ymom[i, j] = 0.5 * ( fluxi_mass[i, j] * (vy[i, j] + vy[i, j+1]) + 
                                      (p[i, j] + p[i, j+1]) * dliy[i, j])
    for i in range(ni-1):
        for j in range(nj):
            fluxj_ymom[i, j] = 0.5 * ( fluxj_mass[i, j] * (vy[i, j] + vy[i+1, j]) + 
                                      (p[i, j] + p[i+1, j]) * dljy[i, j])


    for i in range(ni):
        for j in range(nj-1):
            fluxi_enth[i, j] = 0.5 * (fluxi_mass[i, j] * (hstag[i, j] + hstag[i, j+1]))
    for i in range(ni-1):
        for j in range(nj):
            fluxj_enth[i, j] = 0.5 * (fluxj_mass[i, j] * (hstag[i, j] + hstag[i+1, j]))
    

    sum_fluxes(fluxi_mass, fluxj_mass, delro, ro_inc)
    sum_fluxes(fluxi_enth, fluxj_enth, delroe, roe_inc)
    sum_fluxes(fluxi_xmom, fluxj_xmom, delrovx, rovx_inc)
    sum_fluxes(fluxi_ymom, fluxj_ymom, delrovy, rovy_inc)

    ro = ro_start + ro_inc
    rovx = rovx_start + rovx_inc
    rovy = rovy_start + rovy_inc
    roe = roe_start + roe_inc

    smooth(ro, 50)
    smooth(rovx, 5*500)
    smooth(rovy, 5*500)
    smooth(roe, 50*50000)

    print(nstep)


import matplotlib.pyplot as plt

plt.contourf(x[:,0], y[0,:], np.transpose(tstatic))
plt.plot(xlow, ylow)
plt.plot(xhigh, yhigh)
plt.colorbar()
plt.show()

for i in range(ni):
        for j in range(nj):
            pass