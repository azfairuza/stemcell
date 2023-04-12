import src.physica as psc
import shutil
from cell import Cells
from datetime import datetime
from pathlib import Path
from .nanopattern import Nanopattern

def save(save_obj, time: datetime, num_iteration: int=0, type: str=None):
    """function to save the data into file

    The objects that can be converted into text data are as follows:
    - Nanopattern object
        - PATMAP (Map)
    - Cells object
        - CELLEN (Energy)
        - CELLMAP (Map)
        - CELLAR (Area)
        - CELLCM (Center of Mass)
    
    Parameter
    ---------
    save_obj
       the object that want to be saved as data text
    time: datetime
        time of simulation
    num_iteration: int, default=0
        the current number of iteration
    type: str, default=None
        type of of generated file, it need to be specified if the
        object is `Cells`
    """
    # PATMAP
    if isinstance(save_obj, Nanopattern):
        namefolder = f'./output/{psc.timeFormat(time)}-output/file'
        # build the folder
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        # determine the name of the file
        namefile = f'{namefolder}/PATMAP.txt'
        head_text = 'ligand_id\tx_pos\ty_pos\n'
        with open(namefile, 'w') as output:
            output.write(head_text)
        for ligand in save_obj.ligands:
            content = f'{ligand.id}\t{ligand.x}\t{ligand.y}\n'
            with open(namefile, 'a') as output:
                output.write(content)
    #CELLS
    elif isinstance(save_obj, Cells):
        
        # CELLEN
        if type == "energy" or type == "E" or type == "EN" or type == "Energy" or type == "CELLEN":
            namefolder = f'./output/{psc.timeFormat(time)}-output/file'
            # build the folder
            Path(namefolder).mkdir(parents=True, exist_ok=True)
            # determine the name of the file
            namefile = f'{namefolder}/CELLEN.txt'
            if num_iteration <= 0:
                head_text = 't\t'
                for cell in save_obj.members:
                    head_cell = f'E{cell.id}\t'
                    head_text += head_cell
                head_text += '\n'
                # save the data
                with open(namefile, 'w') as output:
                    output.write(head_text)
            output_text = f'{num_iteration}\t'
            for cell in save_obj.members:
                cell_output = f'{cell.kineticEnergy}\t'
                output_text += cell_output
            output_text += '\n'
            # save the data
            with open(namefile, 'a') as output:
                output.write(output_text)
        
        # CELMAP
        elif type == "MAP" or type == "CELMAP" or type == "map":
            namefolder = f'./output/{psc.timeFormat(time)}-output/file/CELMAP'
            # build the folder
            Path(namefolder).mkdir(parents=True, exist_ok=True)
            # determine the name of the file
            namefile = f'{namefolder}/CELMAP{num_iteration:06}.txt'
            head_text = 'cell_id\tintegrin_id\tx_pos\ty_pos\n'
            # save the data
            with open(namefile, 'w') as output:
                output.write(head_text)
            for cell in save_obj.members:
                for integrin in cell.integrins:
                    content = f'{cell.id}\t{integrin.id}\t{integrin.x}\t{integrin.y}\n'
                    with open(namefile, 'a') as output:
                        output.write(content)
        
        # CELLAR
        elif type == "area" or type == "Area" or type == "CELLAR":
            namefolder = f'./output/{psc.timeFormat(time)}-output/file'
            # build the folder
            Path(namefolder).mkdir(parents=True, exist_ok=True)
            # determine the name of the file
            namefile = f'{namefolder}/CELLAR.txt'
            if num_iteration <= 0:
                head_text = 't\t'
                for cell in save_obj.members:
                    head_cell = f'A{cell.id}\t'
                    head_text += head_cell
                head_text += '\n'
                # save the data
                with open(namefile, 'w') as output:
                    output.write(head_text)
            output_text = f'{num_iteration}\t'
            for cell in save_obj.members:
                area = round(cell.alpha_shape.area, 2)
                cell_output = f'{area}\t'
                output_text += cell_output
            output_text += '\n'
            # save the data
            with open(namefile, 'a') as output:
                output.write(output_text)
            print(f'SYSTEM: CELLAR updated on {namefolder}')
        
        #CELLCM
        elif type == "CELLCM" or type == "CM" or type == "COM":
            namefolder = f'./output/{psc.timeFormat(time)}-output/file'
            # build the folder
            Path(namefolder).mkdir(parents=True, exist_ok=True)
            # determine the name of the file
            namefile = f'{namefolder}/CELLCM.txt'
            if num_iteration <= 0:
                head_text = 't\t'
                for cell in save_obj.members:
                    head_cell = f'x{cell.id}\ty{cell.id}\tn{cell.id}\t'
                    head_text += head_cell
                head_text += '\n'
                # save the data
                with open(namefile, 'w') as output:
                    output.write(head_text)
            output_text = f'{num_iteration}\t'
            for cell in save_obj.members:
                cell_output = f'{cell.x}\t{cell.y}\t{cell.total_bound}\t'
                output_text += cell_output
            output_text += '\n'
            # save the data
            with open(namefile, 'a') as output:
                output.write(output_text)
            print(f'SYSTEM: CELLCM updated on {namefolder}')

def save_Input(time: datetime):
    """procedure to save the input file
    
    Parameter
    ---------
    time: datetime
        the start time of simulation, for generating a unique folder
    """
    namefolder = f'./output/{psc.timeFormat(time)}-output/input'
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    shutil.copy2('./PATCON.txt', namefolder)
    shutil.copy2('./CELCON.txt', namefolder)
    shutil.copy2('./SIMCON.txt', namefolder)
    print(f'SYSTEM: input file has been copied on {namefolder}')
        