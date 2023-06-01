"""The saving data procedure functions"""

# built-in import
import shutil
from datetime import datetime
from pathlib import Path

# local import
import cell as cel
import nanopattern as npt
import physica as psc


def save(save_obj, time: datetime, num_iteration: int = 0, timestep=1.0, data_type: str = None):
    """function to save the data into file

    The objects that can be converted into text data are as follows:
    - Nanopattern object
        - PATMAP (Map)
    - Cells object
        - CELLEN (Energy)
        - CELLMAP (Map)
        - CELLAR (Area)
        - CELLCM (Center of Mass)
    - Input files

    Parameter
    ---------
    save_obj
       the object that want to be saved as data text
    time: datetime
        time of simulation
    num_iteration: int, default=0
        the current number of iteration
    data_type: str, default=None
        type of of generated file, it need to be specified if the
        object is `Cells`
    """
    # PATMAP
    if isinstance(save_obj, npt.Nanopattern):
        namefolder = f"./output/{psc.time_format(time)}-output/file/PATMAP"
        # build the folder
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        # determine the name of the file
        namefile = f"{namefolder}/PATMAP{num_iteration:06}.txt"
        head_text = "ligand_id\tbound\t"
        head_text += "x_grid\ty_grid\t"
        head_text += "x_pos\ty_pos\t"
        head_text += "cell_target\tint_target\t"
        head_text += "\n"
        with open(namefile, "w", encoding="utf-8") as output:
            output.write(f"size\t{save_obj.ligand_size}\n\n")
            output.write(head_text)
        for i in range(len(save_obj._grid)):
            for j in range(len(save_obj._grid[i])):
                for member in save_obj._grid[i][j]:
                    content = f"{member.id_}\t{int(member.bound)}\t"
                    content += f"{j}\t{i}\t"
                    content += f"{member.x_position}\t{member.y_position}\t"
                    if int(member.bound) == 1:
                        content += f"{member.target_cell_id}\t{member.target_integrin_id}\t"
                    else:
                        content += f"{None}\t{None}\t"
                    content += "\n"
                    with open(namefile, "a", encoding="utf-8") as output:
                        output.write(content)


    # CELLS
    elif isinstance(save_obj, cel.Cells):
        # CELLEN
        if data_type in ("energy", "E", "EN", "Energy", "CELLEN"):
            namefolder = f"./output/{psc.time_format(time)}-output/file"
            # build the folder
            Path(namefolder).mkdir(parents=True, exist_ok=True)
            # determine the name of the file
            namefile = f"{namefolder}/CELLEN.txt"
            if num_iteration <= 0:
                head_text = "t\t"
                for cell in save_obj.members:
                    head_cell = f"EK{cell.id_}\tEP{cell.id_}\tEB{cell.id_}"
                    head_text += head_cell
                head_text += "\n"
                # save the data
                with open(namefile, "w", encoding="utf-8") as output:
                    output.write(head_text)
            output_text = f"{round(num_iteration*timestep, 3)}\t"
            for cell in save_obj.members:
                cell_output = f"{cell.kinetic_energy}\t{cell.potential_energy}\t{cell.bonding_energy}\t"
                output_text += cell_output
            output_text += "\n"
            # save the data
            with open(namefile, "a", encoding="utf-8") as output:
                output.write(output_text)

        # CELMAP
        elif data_type in ("MAP", "CELMAP", "map"):
            namefolder = f"./output/{psc.time_format(time)}-output/file/CELMAP"
            # build the folder
            Path(namefolder).mkdir(parents=True, exist_ok=True)
            # determine the name of the file
            namefile = f"{namefolder}/CELMAP{num_iteration:06}.txt"
            head_text = "cell_id\tintegrin_id\tbound\t"
            head_text += "x_pos\ty_pos\t"
            head_text += "vx_pos\tvy_pos\t"
            head_text += "ax_pos\tay_pos\t"
            head_text += "\n"
            # save the data
            with open(namefile, "w", encoding="utf-8") as output:
                output.write(head_text)
            for cell in save_obj.members:
                for integrin in cell.integrins:
                    content = f"{cell.id_}\t{integrin.id_}\t{int(integrin.bound)}\t"
                    content += f"{integrin.x_position}\t{integrin.y_position}\t"
                    content += f"{integrin.x_velocity}\t{integrin.y_velocity}\t"
                    content += f"{integrin.x_acceleration}\t{integrin.y_acceleration}"
                    content += f"{integrin.x_force}\t{integrin.y_force}"
                    content += "\n"
                    with open(namefile, "a", encoding="utf-8") as output:
                        output.write(content)
        
        #CELNBR
        elif data_type in ("NBR", "CELNBR", "Neighbors"):
            namefolder = f"./output/{psc.time_format(time)}-output/file"
            # build the folder
            Path(namefolder).mkdir(parents=True, exist_ok=True)
            # determine the name of the file
            namefile = f"{namefolder}/CELNBR.txt"
            head_text = "cell_id\tintegrin_id\tneighbour_id\n"
            # save the data
            with open(namefile, "w", encoding="utf-8") as output:
                output.write(head_text)
            for cell in save_obj.members:
                for integrin in cell.integrins:
                    content = f"{cell.id_}\t{integrin.id_}"
                    for id_ in integrin.neighbors_id:
                        content += f"\t{id_}"
                    content += "\n"
                    with open(namefile, "a", encoding="utf-8") as output:
                        output.write(content)

        # CELLAR
        elif data_type in ("area", "Area", "CELLAR"):
            namefolder = f"./output/{psc.time_format(time)}-output/file"
            # build the folder
            Path(namefolder).mkdir(parents=True, exist_ok=True)
            # determine the name of the file
            namefile = f"{namefolder}/CELLAR.txt"
            if num_iteration <= 0:
                head_text = "t\t"
                for cell in save_obj.members:
                    head_cell = f"A{cell.id_}\t"
                    head_text += head_cell
                head_text += "\n"
                # save the data
                with open(namefile, "w", encoding="utf-8") as output:
                    output.write(head_text)
            output_text = f"{round(num_iteration*timestep, 3)}\t"
            for cell in save_obj.members:
                area = round(cell.alpha_shape.area, 2)
                cell_output = f"{area}\t"
                output_text += cell_output
            output_text += "\n"
            # save the data
            with open(namefile, "a", encoding="utf-8") as output:
                output.write(output_text)
            print(f"SYSTEM: CELLAR updated on {namefolder}")

        # CELLCM
        elif data_type in ("CELLCM", "CM", "COM"):
            namefolder = f"./output/{psc.time_format(time)}-output/file"
            # build the folder
            Path(namefolder).mkdir(parents=True, exist_ok=True)
            # determine the name of the file
            namefile = f"{namefolder}/CELLCM.txt"
            if num_iteration <= 0:
                head_text = "t\t"
                for cell in save_obj.members:
                    head_cell = f"x{cell.id_}\ty{cell.id_}\tn{cell.id_}\t"
                    head_text += head_cell
                head_text += "\n"
                # save the data
                with open(namefile, "w", encoding="utf-8") as output:
                    output.write(head_text)
            output_text = f"{round(num_iteration*timestep, 3)}\t"
            for cell in save_obj.members:
                cell_output = (
                    f"{cell.x_position}\t{cell.y_position}\t{cell.total_bound}\t"
                )
                output_text += cell_output
            output_text += "\n"
            # save the data
            with open(namefile, "a", encoding="utf-8") as output:
                output.write(output_text)
            print(f"SYSTEM: CELLCM updated on {namefolder}")
    
    # INPUT
    elif isinstance(save_obj, str):
        if save_obj in ("input", "INPUT", "Input"):
            namefolder = f"./output/{psc.time_format(time)}-output/input"
            Path(namefolder).mkdir(parents=True, exist_ok=True)
            shutil.copy2("./input/PATCON.txt", namefolder)
            shutil.copy2("./input/CELCON.txt", namefolder)
            shutil.copy2("./input/SIMCON.txt", namefolder)
            print(f"SYSTEM: input file has been copied on {namefolder}")
    
