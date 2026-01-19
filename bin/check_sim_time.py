import argparse
from argparse import ArgumentTypeError as err
import os
import logging
import sys
import glob
import time


# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# Add sys.sdout as handler
logger.addHandler(logging.StreamHandler(sys.stdout))


class PathType(object):
    def __init__(self, exists=True, type='file', dash_ok=True):
        '''exists:
                True: a path that does exist
                False: a path that does not exist, in a valid parent directory
                None: don't care
           type: file, dir, symlink, None, or a function returning True for valid paths
                None: don't care
           dash_ok: whether to allow "-" as stdin/stdout'''

        assert exists in (True, False, None)
        assert type in ('file','dir','symlink',None) or hasattr(type,'__call__')

        self._exists = exists
        self._type = type
        self._dash_ok = dash_ok

    def __call__(self, string):
        if string=='-':
            # the special argument "-" means sys.std{in,out}
            if self._type == 'dir':
                raise err('standard input/output (-) not allowed as directory path')
            elif self._type == 'symlink':
                raise err('standard input/output (-) not allowed as symlink path')
            elif not self._dash_ok:
                raise err('standard input/output (-) not allowed')
        else:
            e = os.path.exists(string)
            if self._exists==True:
                if not e:
                    raise err("path does not exist: '%s'" % string)

                if self._type is None:
                    pass
                elif self._type=='file':
                    if not os.path.isfile(string):
                        raise err("path is not a file: '%s'" % string)
                elif self._type=='symlink':
                    if not os.path.symlink(string):
                        raise err("path is not a symlink: '%s'" % string)
                elif self._type=='dir':
                    if not os.path.isdir(string):
                        raise err("path is not a directory: '%s'" % string)
                elif not self._type(string):
                    raise err("path not valid: '%s'" % string)
            else:
                if self._exists==False and e:
                    raise err("path exists: '%s'" % string)

                p = os.path.dirname(os.path.normpath(string)) or '.'
                if not os.path.isdir(p):
                    raise err("parent path is not a directory: '%s'" % p)
                elif not os.path.exists(p):
                    raise err("parent directory does not exist: '%s'" % p)

        return string

def parser_input():

    # Parse arguments :
    parser = argparse.ArgumentParser(
        description='Check remaining simulation time.')
    #parser.add_argument('-f', action="store", dest="f",
    #                    help='Input folders', type=list, required=True)

    parser.add_argument('-f', type=PathType(exists=True, type=None),
                        help='Input folders', required=True,
                        nargs='+', dest="f")
    parser.add_argument('-t', type=float,
                        help='Expected simulation time (us), default=1us',
                        dest="t", default=1.0)
    parser.add_argument('-dt', type=float,
                        help='Integration time step (fs), default=4fs',
                        dest="dt", default=4.0)

    return parser

def get_last_line(file):

    with open(file, 'rb') as f:
        try:  # catch OSError in case of a one line file 
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        last_line = f.readline().decode('utf-8')

    return last_line

def get_first_csv_line(file):

    with open(file, 'r') as f:
        f.readline()
        first_line = f.readline()

    return first_line

if __name__ == "__main__":

    my_parser = parser_input()
    args = my_parser.parse_args()
    # logger.info(args)

    chekpoint_frame_num = 100000
    dt = args.dt * 1e-3

    # print("args.f", args.f)

    max_str_size = 0
    for f in args.f:
        files = glob.glob(f + "/*_ST*.csv") + glob.glob(f + "/*_sst2*.csv") + glob.glob(f + "/*_prod*csv")
        # print(files)
        if len(files) > 0:
            # print("f", f)
            for file in files:
                if file.find("full") < 0:
                    if len(f) > max_str_size:
                        max_str_size = len(f)

    for f in args.f:

        # list files with "*_ST*.csv" and "_sst2*.csv"
        files = glob.glob(f + "/*_ST*.csv") + glob.glob(f + "/*_sst2*.csv") + glob.glob(f + "/*_prod*csv")
        csv_list = []


        if len(files) > 0:
            # print("f", f)
            for file in files:
                if file.find("full") < 0:
                    csv_list.append(file)
        
            csv_list.sort()
            last_step = 0
            gap = 0
            running = False

            for csv in csv_list:
                
                time_mod = os.path.getmtime(csv)
                if time.time() -  time_mod < 60:
                    running = True

                #print(csv, time_mod)
                #print(time.time() -  time_mod)
                first_line = get_first_csv_line(csv)
                first_step = int(first_line.split(',')[0])

                if first_step < last_step - chekpoint_frame_num:
                    # print("Timestep was reset")
                    gap += last_step - chekpoint_frame_num

                last_line = get_last_line(csv)
                last_step = int(last_line.split(',')[0])
                ns_a_day = float(last_line.split(',')[-1])

            sim_time = (last_step + gap) * dt / 1e6
            run_perc = sim_time / args.t * 100

            if running:

                us_to_compute = args.t - sim_time
                expected_time = us_to_compute / (ns_a_day * 1e-3)
                expected_days = expected_time
                expected_hours = (expected_days % 1) * 24
                expected_minutes = (expected_hours % 1) * 60
                print(f"{f:{max_str_size}} {sim_time:7.3f} us {run_perc:5.1f}% Complete in {int(expected_days):2d}:{int(expected_hours):02d}:{int(expected_minutes):02d} d:h:min ({ns_a_day:6.1f} ns/day)")
        
            else:
                print(f"{f:{max_str_size}} {sim_time:7.3f} us {run_perc:5.1f}%")
