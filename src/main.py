"""The main module of the simulation"""

# built-in import
import sys
from datetime import datetime
from warnings import filterwarnings

# local import
import cell as cel
import forces
import inputfile as ifile
import integrin as ign
import ligand as lig
import misc
import nanopattern as npt
import physica as psc
import plotter
import save
import simlog

# for debuging purpose
COND = "debug"
filterwarnings("ignore", category=FutureWarning)

# initiate the code
ori = sys.stdout
time, log = simlog.init_simulation()

# open file SIMCON file
simcon = ifile.Read("SIMCON")
metadata = simcon.get("METADATA")
print(f'SYSTEM: simulation run by \t\t: {metadata["username"]}')
print(f'SYSTEM: title of simulation\t\t: {metadata["title"]}')

# get simulation configuration value from SIMCON file
N_ITERATION = int(simcon.get("iteration"))
SAVE_GAP = simcon.get("savegap")
SAVE_FIG = bool(simcon.get("savefig") == 1)
ALPHAVALUE = simcon.get("alpha")
SAVE_PATTERN_MAP = bool(simcon.get("patternmaping") == 1)
SHOW_INTEGRIN = bool(simcon.get("showintegrin") == 1)
SAVE_CENTER_OF_MASS = bool(simcon.get("centerofmass") == 1)
SAVE_CELL_AREA = bool(simcon.get("cellarea") == 1)
SAVE_GIF = bool(simcon.get("gif") == 1)
SAVE_CELL_MAP = bool(simcon.get("cellmaping") == 1)
FORCE_ARROW = bool(simcon.get("forcearrow") == 1)
GET_CONTOUR = bool(simcon.get("contour") == 1)
SHOW_PROGRESS = bool(simcon.get("showprogress") == 1)

# get the physical configuration value from SIMCON file
EPSILON = simcon.get("epsilon")
SPRING_CONSTANT = simcon.get("springconstant")
DAMPING_COEFFICIENT = simcon.get("dampingcoeff")
MIN_FORCE = simcon.get("minForce")
VISCOSITY = simcon.get("viscosity")
TIMESTEP = simcon.get("timestep")


# reset all the simulation dependent variables
lig.Ligand.reset_count()
cel.Cell.reset_count()
ign.Integrin.reset_count()

# create substrate and cells
substrate = npt.Nanopattern()
cells = cel.Cells()
NEAR_DIST = psc.force.nearest_dist_LJ(EPSILON, cells.integrin_size, MIN_FORCE)


# change value into boolean or default value
if SAVE_PATTERN_MAP:
    save.save(substrate, time)

# update cell condition after creation
for cell in cells.members:
    cell.update_alphashape(alpha_value=ALPHAVALUE)
    for integrin_ in cell.integrins:
        integrin_.update_target_bound(cells, substrate)
for cell in cells.members:
    for integrin_ in cell.integrins:
        integrin_.bonding()

# initiate the figure for plot
fig_cell = plotter.init_figure()
fig_contour = plotter.init_figure()

# show result of creation
plotter.show_all(
    fig_cell,
    cells,
    substrate,
    time,
    TIMESTEP,
    show_substrate=True,
    save=SAVE_FIG,
    folder="newsimulate",
    number=0,
    forcearrow=FORCE_ARROW,
    showintegrin=SHOW_INTEGRIN,
)

# get contour plot if necessary
if GET_CONTOUR:
    plotter.contour_plot(
        fig_contour, 
        cells, 
        substrate, 
        time, 
        TIMESTEP, 

        number=0, 
        folder="newsimulate"
    )

# save the cell energy and cell map
save.save(cells, time, timestep=TIMESTEP,data_type="CELLEN")
save.save(cells, time, timestep=TIMESTEP, data_type="CELMAP")
save.save(cells, time, timestep=TIMESTEP, data_type="CELNBR")

# get the cell area data if necessary
if SAVE_CELL_AREA:
    save.save(cells, time, timestep=TIMESTEP, data_type="CELLAR")
if SAVE_CENTER_OF_MASS:
    save.save(cells, time, timestep=TIMESTEP, data_type="CELLCM")

# region <simulation>
iter_simulation = 0
while iter_simulation <= N_ITERATION:
    percent_progress = round(iter_simulation*100/N_ITERATION,3)
    iter_simulation += 1
    print(f"SYSTEM: iteration number {iter_simulation}")
    for cell in cells.members:
        if cells.many:
            surface_integrin = cells.surface_integrins_target(cell, NEAR_DIST)
        for integrin_ in cell.integrins:
            if integrin_.bound is False:
                # create the equation of motion (EOM)
                if integrin_.issurface and cells.many:
                    # in this case the force acting on the integrin are:
                    # 1. nearest another surface integrin
                    # 2. neighboring integrin in the form of spring potential
                    nearest_surface_integrin = misc.filter_by_dist(
                        surface_integrin, NEAR_DIST, integrin_.position
                    )
                    integrin_._nearest = nearest_surface_integrin
                    integrin_._radar_radius = NEAR_DIST
                    eom = lambda x, v: forces.total_force_1(
                        x,
                        v,
                        nearest_surface_integrin,
                        integrin_.neighbors,
                        cell.normal_length,
                        SPRING_CONSTANT,
                        DAMPING_COEFFICIENT,
                        VISCOSITY,
                        EPSILON,
                        integrin_.size,
                        dim=2
                    )
                    integrin_.force = eom(integrin_.position, integrin_.velocity)
                    integrin_.temp_position, integrin_.temp_velocity = psc.integration.eom_rungekutta(
                        integrin_.position,
                        integrin_.velocity,
                        eom,
                        integrin_.mass,
                        TIMESTEP,
                    )
                else:
                    # if only one cell exist, every integrin attracted to nearest ligand
                    # also if the integrin is not a surface it would be attracted to nearest ligand
                    nearest_ligands = substrate.nearest(
                        integrin_.x_position, integrin_.y_position, NEAR_DIST
                    )
                    integrin_._nearest = nearest_ligands
                    integrin_._radar_radius = NEAR_DIST
                    eom = lambda x, v: forces.total_force_2(
                        x,
                        v,
                        nearest_ligands,
                        integrin_.neighbors,
                        cell.normal_length,
                        SPRING_CONSTANT,
                        DAMPING_COEFFICIENT,
                        VISCOSITY,
                        EPSILON,
                        integrin_.size,
                        dim=2
                    )
                    integrin_.force = eom(integrin_.position, integrin_.velocity)
                    integrin_.temp_position, integrin_.temp_velocity = psc.integration.eom_rungekutta(
                        integrin_.position,
                        integrin_.velocity,
                        eom,
                        integrin_.mass,
                        TIMESTEP,
                    )
    # Update all the cell
    for cell in cells.members:
        for integrin_ in cell.integrins:
            integrin_.update()
            cell.update_position()
    for cell in cells.members:
        for integrin_ in cell.integrins:
            integrin_.update_target_bound(cells, substrate)
    for cell in cells.members:
        for integrin_ in cell.integrins:
            integrin_.bonding()
    
    # Calculate potential energy
    for cell in cells.members:
        if cells.many:
            surface_integrin = cells.surface_integrins_target(cell, NEAR_DIST)
        for integrin_ in cell.integrins:
            if integrin_.issurface and cells.many:
                nearest_surface_integrin = misc.filter_by_dist(
                    surface_integrin, NEAR_DIST, integrin_.position
                )
                integrin_.calc_potential(nearest_surface_integrin,
                                        cell.normal_length,
                                        SPRING_CONSTANT,
                                        EPSILON,
                                        integrin_.size
                                        )
            else:
                nearest_ligands = substrate.nearest(
                    integrin_.x_position, integrin_.y_position, NEAR_DIST
                )
                integrin_.calc_potential(nearest_ligands,
                                        cell.normal_length,
                                        SPRING_CONSTANT,
                                        EPSILON,
                                        integrin_.size
                                        )                

    # save the data
    if iter_simulation % SAVE_GAP == 0 or iter_simulation > N_ITERATION:
        # to get the cell shape
        for cell in cells.members:
            cell.update_alphashape(alpha_value=ALPHAVALUE)
        # generate the image
        plotter.show_all(
            fig_cell,
            cells,
            substrate,
            time,
            TIMESTEP,
            show_substrate=True,
            save=SAVE_FIG,
            folder="newsimulate",
            number=iter_simulation,
            showintegrin=SHOW_INTEGRIN,
            forcearrow=FORCE_ARROW,
        )
        # generate contour
        if GET_CONTOUR:
            plotter.contour_plot(
                fig_contour,
                cells,
                substrate,
                time,
                TIMESTEP,
                number=iter_simulation,
                folder="newsimulate",
            )
        # save area
        if SAVE_CELL_AREA is True:
            save.save(cells, time, iter_simulation, timestep=TIMESTEP, data_type="CELLAR")
        # save centter of mass
        if SAVE_CENTER_OF_MASS is True:
            save.save(cells, time, iter_simulation, timestep=TIMESTEP, data_type="CELLCM")
        # save cell mapping
        if SAVE_CELL_MAP is True:
            save.save(cells, time, iter_simulation, timestep=TIMESTEP, data_type="CELMAP",)
        # save nanopattern
        if SAVE_PATTERN_MAP is True:
            save.save(substrate, time, iter_simulation, timestep=TIMESTEP, data_type="PATMAP")

    # endregion
    # save energy
    save.save(cells, time, iter_simulation, timestep=TIMESTEP, data_type="CELLEN")
    if SHOW_PROGRESS:
        sys.stdout = ori
        print(f"\rprogress: [{percent_progress}%]", end="")
        sys.stdout = log

if SAVE_GIF:
    plotter.build_GIF(time)
    print("SYSTEM: GIF created!")

save.save("Input", time)
print("SYSTEM: simulation done!")
elapse_time = datetime.now() - time
print(f"SYSTEM: execution time: {elapse_time}")
log.close()
sys.stdout = ori
print(f"\n{psc.time_format(time)}")
