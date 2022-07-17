import stemcell as sc
import warnings
import sys
from pathlib import Path


warnings.filterwarnings("ignore", category=FutureWarning)
current_time = sc.datetime.now()
namefolder = f'./output/{sc.getTime(current_time)}-output/file'
Path(namefolder).mkdir(parents=True, exist_ok=True)
namefile = f'{namefolder}/SIMLOG.txt'
log = open(namefile, 'w')
sys.stdout = log

print('==================================================================')
print(f'Program made by\t: Achmad Zacky Fairuza')
print(f'email\t\t\t: fairuza.zacky1@gmail.com')
print(f'this program is still under development')
print('==================================================================')

#current running time
print(f'SYSTEM: Simulation is start \t: {current_time}')

# read SIMCON.txt

# open file
simcon = sc.readFile('SIMCON')
metadata = sc.getValue(simcon, 'METADATA')
print(f'SYSTEM: simulation run by \t\t: {metadata["username"]}')
print(f'SYSTEM: title of simulation\t\t: {metadata["title"]}')

# get value
iter_simulation = sc.getValue(simcon, 'iteration')
dstlimit = sc.getValue(simcon, 'dstlimit')
savefig = bool(sc.getValue(simcon, 'savefig'))
movespeed = sc.getValue(simcon, 'movespeed')
alphaValue = sc.getValue(simcon, 'alpha')
gif = sc.getValue(simcon, 'gif')
savegap = sc.getValue(simcon, 'savegap')


# Read PATCON.txt
# open file
patcon = sc.readFile('PATCON')

# get value
dot_number = sc.getValue(patcon, 'dotnumber')
grid_size = sc.getValue(patcon, 'unitsize')
substrate_size = sc.getValue(patcon, 'subsize')
dot_size = sc.getValue(patcon, 'dotsize')
ligand_position = sc.getValue(patcon,'LIGAND') 

sc.Ligand.resetNumber()
substrate = sc.Nanopattern(substrate_size[0], substrate_size[1], 
    grid_size[0], grid_size[1], ligand_position, dot_size)

#Read CELCON.txt

# open file
celcon = sc.readFile('CELCON')

# get value
cell_number = sc.getValue(celcon, 'cellnumber')
integrin_size = sc.getValue(celcon, 'ressize')
integrin_mass = sc.getValue(celcon, 'resmass')
cell_properties = sc.getValue(celcon, 'CELL')

# reset the cell and integrin counter 
sc.Integrin.resetNumber()
sc.Cell.resetNumber()

# build the cell
cells = []
for obj in cell_properties:
    build_cell = sc.Cell(obj[0], obj[1], obj[2], obj[3], obj[4], integrin_size)
    cells.append(build_cell)
    # build_cell.show(substrate)

# Simulate

sc.simulate2(cells, substrate, current_time)
substrate.show(current_time, True, folder='nanopatern')
sc.showAll(cells, substrate, current_time, alphaValue=alphaValue, 
    show_substrate=True, save=savefig, folder='all', line=True)
sc.saveInput(current_time)

if gif == 1:
    sc.buildGIF(current_time)
    print(f'SYSTEM: GIF created!')
sc.filterPicture(int(savegap), current_time)

print(f'SYSTEM: simulation done!')
log.close()
