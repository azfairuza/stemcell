"""Module to init the simulation"""

# built-in import
import sys
from datetime import datetime
from pathlib import Path

# local import
import physica as psc


def init_simulation(cond=None, logfile_name="SIMLOG.txt"):
    """Initiate the simulation log

    Parameter
    ---------
    cond: str, default=None
        It give the special condition for example 'debug' log
    logfile_name: str, default="SIMLOG.txt"
        the name of log file.

    Return
    ------
    time: datetime
        the start time of the simulation
    log: TextIOWrapper
        the log object
    """
    # simulation time
    start_time = datetime.now()
    if cond in ("debug", -1):
        time = "debug"
    elif cond is None:
        time = start_time
    else:
        print("wrong condition")
        raise ValueError()

    namefolder = f"./output/{psc.time_format(time)}-output/file"
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    namefile = f"{namefolder}/{logfile_name}"
    log = open(namefile, "w", encoding="utf-8")
    sys.stdout = log

    print("==================================================================")
    print("Program made by\t: Achmad Zacky Fairuza")
    print("email\t\t\t: fairuza.zacky1@gmail.com")
    print("this program is still under development")
    print("==================================================================")

    # current running time
    print(f"SYSTEM: Simulation is start \t: {time}")

    return time, log
