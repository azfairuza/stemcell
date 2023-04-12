import sys
sys.path.append('/mnt/c/mydata/code/stemcell')
import src.physica as psc
from cell import Cell, Cells
from datetime import datetime
from forces import total_force_1, total_force_2
from input_reader import readFile, getValue
from integrin import Integrin
from ligand import Ligand
from misc import filter_by_dist
from nanopattern import Nanopattern
from plotter import showAll, contourPlot, init_figure, build_GIF
from save_procedure import save, save_Input
from simlog import init_simulation
from warnings import filterwarnings

# for debuging purpose
cond = 'debug'
filterwarnings("ignore", category=FutureWarning)

# initiate the code
time, log = init_simulation()

# open file SIMCON file
simcon = readFile('SIMCON')
metadata = getValue(simcon, 'METADATA')
print(f'SYSTEM: simulation run by \t\t: {metadata["username"]}')
print(f'SYSTEM: title of simulation\t\t: {metadata["title"]}')

# get simulation configuration value from SIMCON file
n_iteration = int(getValue(simcon, 'iteration'))
savefig = getValue(simcon, 'savefig')
centerofmass = getValue(simcon, 'centerofmass')
cellarea = getValue(simcon, 'cellarea')
alphaValue = getValue(simcon, 'alpha')
cellmaping = getValue(simcon, 'cellmaping')
patternmaping = getValue(simcon, 'patternmaping')
showintegrin = getValue(simcon, 'showintegrin')
savegap = getValue(simcon, 'savegap')
gif = getValue(simcon, 'gif')
forcearrow = getValue(simcon, 'forcearrow')
getcontour = getValue(simcon, 'contour')

# get the physical configuration value from SIMCON file
epsilon = getValue(simcon, 'epsilon')
spring_constant = getValue(simcon, 'springconstant')
damping_coefficient = getValue(simcon, 'dampingcoeff')
min_force = getValue(simcon, 'minForce')
viscosity = getValue(simcon, 'viscosity')
dt = getValue(simcon, 'timestep')

# reset all the simulation dependent variables
Ligand.resetCount()
Cell.resetCount()
Integrin.resetCount()

# create substrate and cells 
substrate = Nanopattern()
cells = Cells()
near_dist = psc.nearest_dist_LJ(epsilon, cells.integrin_size, min_force)

# change value into boolean or default value
if patternmaping == 1: save(substrate, time)
if n_iteration is None: n_iteration = 1
savefig = True if savefig == 1 else False
showintegrin = True if showintegrin == 1 else False
centerofmass = True if centerofmass == 1 else False
cellarea = True if cellarea == 1 else False
gif = True if gif == 1 else False
cellmaping = True if cellmaping == 1  else False
forcearrow = True if forcearrow == 1 else False
contour = True if getcontour == 1 else False

# update cell condition after creation
for cell in cells.members:
    cell.updateAlphaShape(alphaValue=alphaValue)
    for integrin in cell.integrins:
        integrin.updateTargetBound(cells, substrate)
for cell in cells.members:
    for integrin in cell.integrins:
        integrin.bonding()

# initiate the figure for plot
fig_cell = init_figure()
fig_contour = init_figure()

# show result of creation
showAll(fig_cell, cells, substrate, time, dt,
        show_substrate=True,
        save=savefig,
        folder='newsimulate',
        number=0,
        forcearrow=forcearrow,
        showintegrin=showintegrin)

# get contour plot if necessary
if getcontour is True:
    contourPlot(fig_contour, cells, substrate, time, dt,
                number=0, 
                folder='newsimulate')

# save the cell energy and cell map
save(cells, time, type="CELLEN")
save(cells, time, type="CELMAP")

# get the cell area data if necessary
if cellarea is True: save(cells, time, type="CELLAR")
if centerofmass is True: save(cells, time, type="CELLCM")

# region <simulation>
iter_simulation = 0
while(iter_simulation <= n_iteration):
    iter_simulation += 1
    print(f'SYSTEM: iteration number {iter_simulation}')
    for cell in cells.members:
        if cells.many:
            surface_integrin = cells.surfaceIntegrinsTarget(cell)
        for integrin in cell.integrins:
            if integrin.bound is False:
                # create the equation of motion (EOM)
                if integrin.isSurface and cells.many:
                    # in this case the force acting on the integrin are:
                    # 1. nearest another surface integrin
                    # 2. neighboring integrin in the form of spring potential
                    nearest_surface_integrin = filter_by_dist(surface_integrin, near_dist, integrin.position)                   
                    eom = lambda x,v: total_force_1(x, v, 
                                                    nearest_surface_integrin, 
                                                    integrin.neighbors, 
                                                    cell.normal_length,
                                                    spring_constant,
                                                    damping_coefficient,
                                                    epsilon,
                                                    integrin.size)
                    integrin.temp_position, integrin.temp_speed = psc.eom_rungekutta(integrin.position, 
                                                                                     integrin.speed, 
                                                                                     eom, 
                                                                                     integrin.mass, 
                                                                                     dt)
                else:
                    # if only one cell exist, every integrin attracted to nearest ligand
                    # also if the integrin is not a surface it would be attracted to nearest ligand
                    nearest_ligands = substrate.nearest(integrin.x, integrin.y, near_dist)
                    eom = lambda x,v: total_force_2(x, v, 
                                                    nearest_ligands, 
                                                    integrin.neighbors, 
                                                    cell.normal_length,
                                                    spring_constant,
                                                    damping_coefficient,
                                                    epsilon,
                                                    integrin.size)
                    integrin.temp_position, integrin.temp_speed = psc.eom_rungekutta(integrin.position, 
                                                                                     integrin.speed, 
                                                                                     eom, 
                                                                                     integrin.mass, 
                                                                                     dt)    
    # Update all the cell   
    for cell in cells.members:
        for integrin in cell.integrins:
            integrin.update()
            cell.updatePosition()
    for cell in cells.members:
        for integrin in cell.integrins:
            integrin.updateTargetBound(cells, substrate)
    for cell in cells.members:
        for integrin in cell.integrins:
            integrin.bonding()
    
    # save the data
    if iter_simulation%savegap == 0 or iter_simulation > n_iteration:
        # to get the cell shape
        for cell in cells.members:
            cell.updateAlphaShape(alphaValue=alphaValue)
        # generate the image
        showAll(fig_cell, 
                cells, 
                substrate, 
                time, 
                dt, 
                show_substrate=True,
                save=savefig,
                folder='newsimulate',
                number=iter_simulation,
                showintegrin=showintegrin)
        # generate contour
        if getcontour is True: contourPlot(fig_contour, 
                                           cells, 
                                           substrate, 
                                           time, 
                                           dt, 
                                           number=iter_simulation, 
                                           folder='newsimulate')
        # save area
        if cellarea is True: save(cells, time, type="CELLAR")
        # save centter of mass
        if centerofmass is True: save(cells, time, type="CELLCM")
        # saave cell mapping
        if cellmaping is True: save(cells, time, type="CELLMAP")   
    
    
# endregion
    # save energy
    save(cells, time, type="CELLEN")

if gif is True:
    build_GIF(time)
    print(f'SYSTEM: GIF created!')

save_Input(time)
print(f'SYSTEM: simulation done!')
elapse_time = datetime.now() - time
print(f'SYSTEM: execution time: {elapse_time}')
log.close()