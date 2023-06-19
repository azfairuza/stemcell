"""module to plot the data"""

# built-in imports
import os
from pathlib import Path
from datetime import datetime
import itertools

# third party imports
import imageio
import matplotlib.pyplot as plt
import numpy as np
import cv2 as cv
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from PIL import Image

# local import
import physica as psc
import cell as cel
import nanopattern as npt
import inputfile as ifile


def init_figure():
    """Procedure to create the figure"""
    fig = plt.figure(figsize=(20, 20))
    fig.dpi = 100
    return fig


def show_all(
    fig: Figure,
    cells: cel.Cells,
    substrate: npt.Nanopattern,
    time: datetime,
    timestep,
    show_substrate: bool = False,
    save: bool = False,
    number: int = 0,
    folder: str = None,
    showintegrin: bool = True,
    forcearrow: bool = False,
    line: bool = False,
    showID: bool = False,
    tracking: bool = False
):
    """Procedure to show all element of simulation including cells and
    nanopattern.

    Parameter
    ---------
    cells: :obj: `Cells`
        The collection of cell
    substrate: :obj: `Nanopattern`
        The nanopatterned substrate which is consisted of ligands
    time: datetime
        The simulation start time
    show_substrate: default=False
        The plot is not printing the substrate when `False`
    save: default=False
        Whether the generated plot be saved or not
    number: int, default=1
        The number of the file
    folder: str, default=None
        The folder to place the saved plot image
    showintegrin: default=True
        Whether the plot will include the integrin or not
    forcearrow: default=False
        Whether the plot will include arrow which indicates the force
        vector or not
    """
    # determine the folder's name
    if folder is None:
        namefolder = f"./output/{psc.time_format(time)}-output/figure"
    else:
        namefolder = f"./output/{psc.time_format(time)}-output/figure/{folder}"

    # build the folder
    Path(namefolder).mkdir(parents=True, exist_ok=True)

    # determine the file name
    namefile = f"{namefolder}/{number:06}.jpg"
    
    fig.clf()
    axis = fig.gca()
    axis.set_aspect('equal')
    
    # print the cells
    for cell in cells.members:
        if show_substrate is True:
            # draw the free ligands
            psc.circles(
                axis,
                substrate.x_pos_list[0],
                substrate.y_pos_list[0],
                substrate.ligand_size,
                "green",
                alpha=1,
                ec="none",
            )
            # draw the bound ligands
            psc.circles(
                axis,
                substrate.x_pos_list[1],
                substrate.y_pos_list[1],
                substrate.ligand_size,
                "yellow",
                alpha=1,
                ec="none",
            )
        # draw center of mass
        cm_position = (cell.x_position, cell.y_position)
        cm_patch = Circle(cm_position, 2, color="yellow", alpha=0.2)
        axis.add_patch(cm_patch)
        # draw the outer layer
        print(cell.alpha_shape)
        psc.polygon_patch(axis, cell.alpha_shape, alpha=0.2)
        if showintegrin is True:
            # draw the free integrins
            psc.circles(
                axis,
                cell.x_pos_list[0],
                cell.y_pos_list[0],
                cell.integrin_size,
                "red",
                alpha=1,
                ec="none",
            )
            # draw the bound integrins
            psc.circles(
                axis,
                cell.x_pos_list[1],
                cell.y_pos_list[1],
                cell.integrin_size,
                "yellow",
                alpha=1,
                ec="none",
            )
            # draw the id
            if showID is True:
                for cell_obj in cells.members:
                    for integrin_obj in cell_obj.integrins:
                        axis.text(integrin_obj.x_position, integrin_obj.y_position, integrin_obj.id_)

        # draw line between integrins
        for integrin_ in cell.integrins:
            if line is True:
                for neighbor in integrin_.neighbors:
                    x_line = [integrin_.x_position, neighbor.x_position]
                    y_line = [integrin_.y_position, neighbor.y_position]
                    axis.plot(x_line, y_line, color="black", linewidth="1", alpha=0.1)
            if forcearrow is True and integrin_.bound is False:
                axis.arrow(
                    integrin_.x_position,
                    integrin_.y_position,
                    integrin_.x_force,
                    integrin_.y_force,
                    color="blue",
                    width=0.3,
                )
            if tracking:
                if integrin_.id_ == 145:
                    x_neighbor_pos = []
                    y_neighbor_pos = []
                    monitor = Circle(integrin_.position, integrin_.size, color="black")
                    radar = Circle(integrin_.position, integrin_._radar_radius, color="blue", alpha=0.2)
                    for neighbor in integrin_._nearest:
                        x_neighbor_pos.append(neighbor.x_position)
                        y_neighbor_pos.append(neighbor.y_position)
                    psc.circles(axis, x_neighbor_pos, y_neighbor_pos, integrin_.size, color="blue")
                    axis.add_patch(monitor)
                    axis.add_patch(radar)
    axis.set_xlim(0, substrate.width)
    axis.set_ylim(0, substrate.height)
    axis.set_title(f"time = {number*timestep}", fontsize=22)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)

    # save if necessary
    if save is True:
        print(f"SYSTEM: figure {number:06}.jpg saved on {namefolder}")
        fig.savefig(namefile, bbox_inches="tight", dpi=100)


def contour_plot(
    fig: Figure,
    cells: cel.Cells,
    substrate: npt.Nanopattern,
    time: datetime,
    timestep,
    number=0,
    folder=None,
):
    """Procedure to show contour plot of integrin in cells"""
    # determine the folder's name
    if folder is None:
        namefolder = f"./output/{psc.time_format(time)}-output/figure"
    else:
        namefolder = f"./output/{psc.time_format(time)}-output/figure/{folder}"

    # build the folder
    Path(namefolder).mkdir(parents=True, exist_ok=True)

    # determine the file name
    namefile = f"{namefolder}/C{number:06}.jpg"

    # draw the figure
    fig.clf()
    axis = fig.gca()
    axis.set_aspect('equal')

    # concatenate the list as we treated equally between bounded
    # integrins and unbounded integrins
    x_pos_list = list(itertools.chain.from_iterable(cells.x_list))
    y_pos_list = list(itertools.chain.from_iterable(cells.y_list))

    # get the limit
    xmin = 0
    xmax = substrate.width
    ymin = 0
    ymax = substrate.height

    xx_grid, yy_grid = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    fig_matrix = np.zeros((100, 100))
    x_gap = substrate.width / 100
    y_gap = substrate.height / 100
    for i in range(len(x_pos_list)):
        x_index = int(x_pos_list[i] // x_gap)
        y_index = int(y_pos_list[i] // y_gap)
        fig_matrix[x_index][y_index] += 1

    fig_matrix_blur = cv.GaussianBlur(fig_matrix, (7, 7), 0)

    # position = np.vstack([xx.ravel(), yy.ravel()])
    # values = np.vstack([x, y])
    # kernel = st.gaussian_kde(values)
    # f = np.reshape(kernel(position).T, xx.shape)

    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    cfset = axis.contourf(
        xx_grid, yy_grid, fig_matrix_blur, cmap="turbo", vmin=0, vmax=0.5
    )
    # cset = fig.gca().contour(xx, yy, f, colors='k')
    # fig.gca().clabel(cset, inline=1, fontsize=10)
    cbar = fig.colorbar(cfset, aspect=5, shrink=0.5)
    cbar.set_label(
        "density (integrin/${\mu}m^2$)", rotation=270, fontsize=22, labelpad=40
    )
    cbar.ax.set_label("density (integrin/${\mu}m^2$)")
    cbar.ax.tick_params(labelsize=22)
    axis.set_xlabel("X (${\mu}m$)", fontsize=22)
    axis.set_ylabel("Y (${\mu}m$)", fontsize=22)
    axis.set_title(f"time = {number*timestep}", fontsize=22)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)

    # save if necessary
    fig.savefig(namefile, bbox_inches="tight", dpi=100)


def build_GIF(time: datetime):
    """function to generate GIF using python

    parameter
    ---------
    time: datetime
        start time of simulation, to generate a specific folder
    """
    image_dir = f"./output/{psc.time_format(time)}-output/figure/newsimulate"
    gif_dir = f"./output/{psc.time_format(time)}-output/figure/gif"
    Path(gif_dir).mkdir(parents=True, exist_ok=True)
    images = []
    for namefile in sorted(os.listdir(image_dir)):
        if namefile.endswith(".jpg"):
            pathfile = os.path.join(image_dir, namefile)
            images.append(
                Image.fromarray(imageio.imread(pathfile)).resize((1000, 1000))
            )
    pathgif = os.path.join(gif_dir, "simulation.gif")
    imageio.mimsave(pathgif, images, fps=5)

def plot_energy():
    """method to plot the energy from energy file"""
    folder_code = input("folder time code: ")
    file_dir = f"./output/{folder_code}-output/file/CELLEN.txt"
    print(f"open: {file_dir}")
    with open(file_dir, "r", encoding="utf-8") as data_file:
        lst_strng = data_file.readlines()
    stripped_strng = ifile.filter_item(lst_strng)
    counter = []
    kinetic_energy = []
    potential_energy = []
    bonding_energy = []
    total_energy = []
    for i in range(1, len(stripped_strng)):
        dummy = stripped_strng[i].split()
        counter.append(float(dummy[0]))
        kinetic_energy.append(float(dummy[1]))
        potential_energy.append(float(dummy[2]))
        bonding_energy.append(float(dummy[3]))
        total_energy.append(float(dummy[1]) + float(dummy[2]) + float(dummy[3]))
        # total_energy.append(float(dummy[1]) + float(dummy[2]))
    fig = plt.figure()
    plt.plot(counter, kinetic_energy, label="EK")
    plt.plot(counter, potential_energy, label="EP")
    plt.plot(counter, bonding_energy, label="EB")
    plt.plot(counter, total_energy, label="ET")
    plt.legend()
    plt.show()    

def _plot_map(fig: Figure, map_type: str, folder_code: str, file_timestamp: str):
    """method to plot map"""
    if map_type in ["patmap", "PATMAP", "nanopattern"]:
        file_dir = f"./output/{folder_code}-output/file/PATMAP/PATMAP{int(file_timestamp):06}.txt"
        print(f"open: {file_dir}")
        with open(file_dir, "r", encoding="utf-8") as data_file:
            lst_strng = data_file.readlines()
        stripped_strng = ifile.filter_item(lst_strng)
        x_pos_list = [[],[]]
        y_pos_list = [[],[]]
        dummy = stripped_strng[0].split()
        ligand_size = float(dummy[1])
        for i in range(2, len(stripped_strng)):
            dummy = stripped_strng[i].split()
            x_pos_list[int(dummy[1])].append(float(dummy[4]))
            y_pos_list[int(dummy[1])].append(float(dummy[5]))
        axis = fig.gca()
        axis.set_aspect("equal")
        # draw the free ligands
        psc.circles(
            axis,
            x_pos_list[0],
            y_pos_list[0],
            ligand_size,
            "green",
            alpha=1,
            ec="none",
        )
        # draw the bound ligands
        psc.circles(
            axis,
            x_pos_list[1],
            y_pos_list[1],
            ligand_size,
            "yellow",
            alpha=1,
            ec="none",
        )
    elif map_type in ["celmap", "CELMAP", "cell", "cells"]:
        file_dir = f"./output/{folder_code}-output/file/CELMAP/CELMAP{int(file_timestamp):06}.txt"
        print(f"open: {file_dir}")
        with open(file_dir, "r", encoding="utf-8") as data_file:
            lst_strng = data_file.readlines()
        stripped_strng = ifile.filter_item(lst_strng)
        x_pos_list = [[],[]]
        y_pos_list = [[],[]]
        dummy = stripped_strng[0].split()
        ligand_size = float(dummy[1])
        for i in range(2, len(stripped_strng)):
            dummy = stripped_strng[i].split()
            x_pos_list[int(dummy[2])].append(float(dummy[3]))
            y_pos_list[int(dummy[2])].append(float(dummy[4]))
        axis = fig.gca()
        axis.set_aspect("equal")
        # draw the free integrins
        psc.circles(
            axis,
            x_pos_list[0],
            y_pos_list[0],
            ligand_size,
            "green",
            alpha=1,
            ec="none",
        )
        # draw the bound integrins
        psc.circles(
            axis,
            x_pos_list[1],
            y_pos_list[1],
            ligand_size,
            "yellow",
            alpha=1,
            ec="none",
        )


def plot_patmap():
    """method to plot the nanopattern"""
    folder_code = input("folder time code: ")
    file_timestamp = input("file timestamp: ")
    fig = plt.figure()
    _plot_map(fig, "PATMAP", folder_code, file_timestamp)
    plt.show() 

def plot_celmap():
    """method to plot the cell"""
    folder_code = input("folder time code: ")
    file_timestamp = input("file timestamp: ")
    fig = plt.figure()
    _plot_map(fig, "cELMAP", folder_code, file_timestamp)
    plt.show() 

    



