import itertools
import imageio
import matplotlib.pyplot as plt
import numpy as np
import os
import src.physica as psc
from cell import Cells
from cv2 import GaussianBlur
from datetime import datetime
from descartes import PolygonPatch
from nanopattern import Nanopattern
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from pathlib import Path
from PIL import Image



def init_figure():
    """Procedure to create the figure"""
    fig = plt.figure(figsize=(20,20))
    fig.dpi = 100
    fig.add_axes(aspect="equal")
    return fig

def showAll(
    fig: Figure,
    cells: Cells,
    substrate: Nanopattern,
    time: datetime,
    dt,
    show_substrate: bool=False,
    save: bool=False,
    number: int=0,
    folder: str=None,
    showintegrin: bool=True,
    forcearrow: bool=False,
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
        namefolder = f"./output/{psc.timeFormat(time)}-output/figure"
    else:
        namefolder = f"./output/{psc.timeFormat(time)}-output/figure/{folder}"
    
    # build the folder
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    
    # determine the file name
    namefile = f"{namefolder}/{number:06}.jpg"
    
    axis = fig.gca()
    axis.collections.clear()

    # print the cells
    for cell in cells.members:
        if show_substrate is True:
            # draw the free ligands
            psc.circles(
                substrate.x_pos_list[0],
                substrate.y_pos_list[0],
                substrate.ligand_size,
                "green",
                alpha=1,
                ec="none",
            )
            # draw the bound ligands
            psc.circles(
                substrate.x_pos_list[1],
                substrate.y_pos_list[1],
                substrate.ligand_size,
                "yellow",
                alpha=1,
                ec="none",
            )
        # draw center of mass
        cm_position = (cell.x, cell.y)
        cm_patch = Circle(cm_position, 2, color="yellow", alpha=0.2)
        axis.add_patch(cm_patch)
        # draw the outer layer
        outerLayer = PolygonPatch(cell.alpha_shape, alpha=0.2)
        axis.add_patch(outerLayer)
        if showintegrin is True:
            # draw the free integrins
            psc.circles(
                cell.x_pos_list[0],
                cell.y_pos_list[0],
                cell.integrin_size,
                "red",
                alpha=1,
                ec="none",
            )
            # draw the bound integrins
            psc.circles(
                cell.x_pos_list[1],
                cell.y_pos_list[1],
                cell.integrin_size,
                "yellow",
                alpha=1,
                ec="none",
            )
            
    # draw line between integrins
        for integrin in cell.integrins:
            line = False
            if line is True:
                for neighbor in integrin.neighbors:
                    x_line = [integrin.x, neighbor.x]
                    y_line = [integrin.y, neighbor.y]
                    line_plot = plt.plot(x_line, y_line, color="black", linewidth="1")
            if forcearrow is True and integrin.bound is False:
                arrow_plot = plt.arrow(
                    integrin.x,
                    integrin.y,
                    integrin.fx,
                    integrin.fy,
                    color="blue",
                    width=0.3,
                )
    
    axis.set_xlim(0, substrate.width)
    axis.set_ylim(0, substrate.height)
    axis.set_title(f"time = {number*dt}", fontsize=22)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    
    # save if necessary
    if save is True:
        print(f"SYSTEM: figure {number:06}.jpg saved on {namefolder}")
        fig.savefig(namefile, bbox_inches="tight", dpi=100)
    
def contourPlot(
    fig: Figure,
    cells: Cells, 
    substrate: Nanopattern, 
    time: datetime, 
    dt, 
    number=0, 
    folder=None
):
    """Procedure to show contour plot of integrin in cells
    
    """
    # determine the folder's name
    if folder is None:
        namefolder = f'./output/{psc.timeFormat(time)}-output/figure'
    else:
        namefolder = f'./output/{psc.timeFormat(time)}-output/figure/{folder}'
    
    # build the folder
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    
    # determine the file name        
    namefile = f'{namefolder}/C{number:06}.jpg'
    
    # draw the figure
    axis = fig.gca()
    axis.collections.clear()
    
    # concatenate the list as we treated equally between bounded
    # integrins and unbounded integrins
    x = list(itertools.chain.from_iterable(cells.x_list))
    y = list(itertools.chain.from_iterable(cells.y_list))
    
    # get the limit
    xmin = 0
    xmax = substrate.width
    ymin = 0
    ymax = substrate.height

    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    f = np.zeros((100,100))
    x_gap = substrate.width/100
    y_gap = substrate.height/100
    for i in range(len(x)):
        x_index = int(x[i]//x_gap)
        y_index = int(y[i]//y_gap)
        f[x_index][y_index] += 1

    fBlur = GaussianBlur(f, (9,9), 0)
    
    # position = np.vstack([xx.ravel(), yy.ravel()])
    # values = np.vstack([x, y])
    # kernel = st.gaussian_kde(values)
    # f = np.reshape(kernel(position).T, xx.shape)

    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    cfset = axis.contourf(xx, yy, fBlur, cmap='inferno', vmin=0, vmax=1)
    # cset = fig.gca().contour(xx, yy, f, colors='k')
    # fig.gca().clabel(cset, inline=1, fontsize=10)
    cbar = fig.colorbar(cfset, aspect=5, shrink=0.5)
    cbar.set_label('density (integrin/${\mu}m^2$)', rotation=270, fontsize=22, labelpad=40)
    cbar.ax.set_label('density (integrin/${\mu}m^2$)')
    cbar.ax.tick_params(labelsize=22)
    axis.set_xlabel('X (${\mu}m$)', fontsize=22)
    axis.set_ylabel('Y (${\mu}m$)', fontsize=22)
    axis.set_title(f"time = {number*dt}", fontsize=22)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)

    
    # save if necessary
    fig.savefig(namefile, bbox_inches='tight', dpi=100)

def build_GIF(time: datetime):
    """function to generate GIF using python

    parameter
    ---------
    time: datetime
        start time of simulation, to generate a specific folder
    """
    image_dir = f'./output/{psc.timeFormat(time)}-output/figure/newsimulate'
    gif_dir = f'./output/{psc.timeFormat(time)}-output/figure/gif'
    Path(gif_dir).mkdir(parents=True, exist_ok=True)
    images = []
    for namefile in sorted(os.listdir(image_dir)):
        if namefile.endswith('.jpg'):
            pathfile = os.path.join(image_dir, namefile)
            images.append(Image.fromarray(imageio.imread(pathfile)).resize((1000,1000)))
    pathgif = os.path.join(gif_dir, 'simulation.gif')
    imageio.mimsave(pathgif, images, fps=5)