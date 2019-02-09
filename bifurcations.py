import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from PyDSTool import *
from matplotlib import rc

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# ODE system
dX1  = '(s+mu)*I-(2.*mu+rho+s)*X1-(gam + uT)*X1';
dP01 = 'rho*(1.-h*(1-uC))*X1*(1.-X1/Xe)-(s+phi*h*(1-uC)+2.*mu)*P01+(gam + uT)*(I-X1-2.*P01)'
dI   = 'rho*h*(1-uC)*X1*(1.-X1/Xe)+phi*h*(1-uC)*P01-mu*I-(gam + uT)*I'

# setup solver
DSargs = args(name='ODE')
DSargs.tdomain = [0.0, 20.0]

DSargs.pars = {
'rho': 5.0,
's': 2.0,
'mu': 1./9.,
'phi': 52.,
'N': 1.e6,
##-- tricho
#'h': 0.115,
#'gam': 0.727,
##-- gono
#'h': 0.348,
#'gam': 1.538,
#-- chlamy
#'h': 0.129,
#'gam': 0.855,
#-- HPV
'h': 0.073,
'gam': 0.5,
'uT': 0.0,
'uC': 0.0}

dummie = DSargs.pars.copy()

dummie2 = dummie.update({
'Xe': DSargs.pars['N']*(DSargs.pars['s']+2.*DSargs.pars['mu'])/(DSargs.pars['s']+2.*DSargs.pars['mu']+DSargs.pars['rho']),
'Pe':
DSargs.pars['N']*DSargs.pars['rho']/(2.*(DSargs.pars['s']+2.*DSargs.pars['mu']+DSargs.pars['rho']))})

DSargs.pars = dummie;

print DSargs.pars['Xe']
print DSargs.pars['Pe']

DSargs.varspecs = {'X1': dX1, 'P01': dP01, 'I': dI}
DSargs.ics = {'X1': 1.e5, 'P01': 3.e5, 'I': 7.e5}
testODE = Generator.Vode_ODEsystem(DSargs)

# solve ODE
print 'Integrating...'
start = clock()
testtraj = testODE.compute('test')
print '  ... finished in %.3f seconds.\n' % (clock()-start)

# plot ODE
plt.figure()
plotData = testtraj.sample(dt=0.001)  # Pointset
solu1=plt.plot(plotData['t'], plotData['X1'])
solu1=plt.plot(plotData['t'], plotData['P01'])
solu1=plt.plot(plotData['t'], plotData['I'])
plt.legend(['X1','P01','I'])

print plotData['X1'][-1]
print plotData['P01'][-1]
print plotData['I'][-1]

#%%
# bifurcation

# initial condition near a steady-state
testODE.set(ics = {
'X1':  plotData['X1'][-1], 
'P01': plotData['P01'][-1], 
'I':   plotData['I'][-1]
})

PyCont = ContClass(testODE)

# curve 1
PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['uC']
PCargs.MaxNumPoints = 1000
PCargs.StepSize = 1e3
PCargs.MaxStepSize = 1e3
PCargs.MinStepSize = 1e2
PCargs.LocBifPoints = 'all'

PyCont.newCurve(PCargs)
PyCont['EQ1'].backward()

fig = plt.figure(figsize=(4,3))
PyCont.computeEigen()
PyCont['EQ1'].display(
stability=True,
linewidth=2)

# curve 2
PCargs = args(name='EQ2', type='EP-C')
PCargs.initpoint = {
'X1': 0, 
'P01': 0, 
'I': 0
}
PCargs.freepars = ['uC']
PCargs.MaxNumPoints = 200
PCargs.MaxStepSize = 1e-2
PCargs.LocBifPoints = 'all'

PyCont.newCurve(PCargs)
PyCont['EQ2'].forward()
PyCont.computeEigen()

plt.figure()

PyCont['EQ1'].display(
stability=True,
linewidth=2)

PyCont['EQ2'].display(
stability=True,
linewidth=2)

PyCont.plot.togglePoints(visible='off', bylabel='P1')
PyCont.plot.toggleLabels(visible='off', bylabel='P1')
PyCont.plot.togglePoints(visible='off', bylabel='P2')
PyCont.plot.toggleLabels(visible='off', bylabel='P2')
PyCont.plot.setLabels('$BP1$', bylabel='BP1')

plt.title('')
plt.xlabel('$u_C^0$',fontsize=18)
plt.ylabel('Infected population $I^*$',fontsize=18)
plt.tick_params(axis='both', labelsize=12)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.xlim([0,1.0])
plt.ylim([-1e5,1e6])
plt.tight_layout()
PyCont.plot.refresh()
