import physica as psc
import sys
from datetime import datetime
from pathlib import Path

def init_simulation(cond=None, logfile_name='SIMLOG.txt'):
    # simulation time
    start_time =datetime.now()
    if cond == 'debug' or cond == -1:
        time = 'debug'
    elif cond == None:
        time = start_time
    else:
        print('wrong condition')
        raise ValueError()

    namefolder = f'./output/{psc.timeFormat(time)}-output/file'
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    namefile = f'{namefolder}/{logfile_name}'
    log = open(namefile, 'w')
    sys.stdout = log

    print('==================================================================')
    print(f'Program made by\t: Achmad Zacky Fairuza')
    print(f'email\t\t\t: fairuza.zacky1@gmail.com')
    print(f'this program is still under development')
    print('==================================================================')

    #current running time
    print(f'SYSTEM: Simulation is start \t: {time}')

    return time, log