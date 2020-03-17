#!/usr/bin/python3

import glob, os, sys, shutil
import multiprocessing
import subprocess
import argparse
import yaml
import itertools
import re
import time

def create_argparser():
    parser = argparse.ArgumentParser(description='runs gams processes in parallel\n')
    parser.add_argument('-s', '--sim', default='sim.gms', dest='sim',
                                        help='defines the yaml file')
    parser.add_argument('-y', '--yaml', default='csets.yml', dest='yml',
                                        help='defines the yaml file')
    parser.add_argument('-p', '--procs', type=int, default=multiprocessing.cpu_count()-1, dest='max_procs',
                                        help='how many processors should be used,\n\t default and maximum is number of processors - 1')
    parser.add_argument('--clean', dest='clean', action='store_true',
                                        help='')
    parser.add_argument('--resume', dest='resume', action='store_true',
                                        help='')
    parser.add_argument('--run', dest='run', action='store_true',
                                        help='')

    return parser

def create_tuples(ydata):
    comblines = list()
    permlines = list()
    
    print(ydata)

    if 'comb' in ydata:
        ks = ydata['comb']['keys']
        for t in ydata['comb']['values']:
            line = ''
            for i,v in enumerate(t):
                line+="".join(("--",ks[i], "=", str(v)))
                line+=";"
            line = line.strip()
            comblines.append(line)

    if 'perm' in ydata:
        values = list(ydata['perm'].values())
        permvals = list(itertools.product(*values))

        for v in permvals:
            line = ""
            for index, key in enumerate(ydata['perm']):
                line+="".join(("--",key, "=", str(v[index])))
                line+=";"
            permlines.append(line)
    if not comblines:
        comblines = ['']
    if not permlines:
        permlines = ['']
    tuples = list(itertools.product(comblines, permlines))

    result = list()
    for tup in tuples:
        tup = tup[0]+tup[1]
        tup = tup.split(';')[:-1]
        result.append(tup)
    return result

def cleandir(path, remove):
    if remove:
        try:
            shutil.rmtree(path)
        except:
            pass
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError:
            print ("Creation of the directory %s failed" % path)
    return 0

def execgams_parallel(params):
    return execgams(*params)

def execgams(gamsbin, sim, tup, simdir, resdir, perc):
    tupdir = ''.join(tup)
    tupdir = tupdir.replace('--','-')
    tupdir = tupdir.replace('=','.')
    tupdir = tupdir[1:]
    pdir = os.path.join(simdir, tupdir)
    
    try:
        os.mkdir(pdir)
    except:
        pass

    # make list with commands and parameters
    # each parameter has to be a separate item in the list
    cmd = []
    cmd.append(gamsbin)
    cmd.append(str(os.path.join(os.getcwd(),sim)))
    cmd.extend(tup)
    cmd.append(str('PutDir='+pdir))
    cmd.append('LogOption=2')
    #print(cmd)
    print("%s%%: Processing %s %s" % (perc, sim, tup))
    sp = subprocess.run(cmd)

    return sp.returncode

def run_simpool(gamsbin, sim, tuples, simdir, resdir, mprocs, resume):

    if resume:
        for tup in tuples:
            pdir = os.path.join(simdir,'-'.join(re.findall(r'\d+', tup)))
            donefile = os.path.join(pdir,"gams.done")
            if os.path.isfile(donefile):
                tuples.remove(tup)

    pool = multiprocessing.Pool(processes=mprocs)
    totaliterations=len(tuples)

    print("Starting a total of %s iterations using %s of %s CPU Threads" % (totaliterations,mprocs,multiprocessing.cpu_count()))

    for index, tup in enumerate(tuples):
        perc = int(index/totaliterations*100)
        params = [(gamsbin, sim, tup, simdir, resdir, perc)]
        res = pool.map_async(execgams_parallel, params)
        
    pool.close()
    pool.join()
    
    while True:
        if res.ready():
            break
        else:
            time.sleep(10)
    
    for index, tup in enumerate(tuples):
        move_results(tup, simdir, resdir)

    if totaliterations == 0:
        sys.exit("Nothing to do")

    return 0

def move_results(tup, simdir, resdir):
    tupdir = ''.join(tup)
    tupdir = tupdir.replace('--','-')
    tupdir = tupdir.replace('=','.')
    tupdir = tupdir[1:]
    pdir = os.path.join(simdir, tupdir)
    
    dstfile = os.path.join(resdir, tupdir+'.csv')
    donefile = os.path.join(pdir,"gams.done")
    
    open(donefile, 'x').close()

    resfile = os.path.join(pdir,"result.csv")
    #print(resfile, '->', dstfile)

    if os.path.isfile(resfile):
        shutil.move(resfile, dstfile)

def main():

    parser = create_argparser()
    args = parser.parse_args()

    gamsbin = "gams"
    if sys.platform == 'win32':
        gamsbin = "gams.exe"

    if args.run != True:
        print('Not running with args:',args)
        sys.exit()

    with open(args.yml, 'r') as stream:
        ydata = yaml.load(stream, Loader=yaml.BaseLoader)

    try:
        simdir = ydata['config']['simdir']
    except:
        simdir = os.path.join(os.getcwd(), 'sim')

    try:
        resdir = ydata['config']['resdir']
    except:
        resdir = os.path.join(os.getcwd(), 'res')

    if args.clean == True:
        cleandir(simdir, 1)
        cleandir(resdir, 1)
    else:
        cleandir(simdir, 0)
        cleandir(resdir, 0)

    tuples = create_tuples(ydata['tuples'])
    run_simpool(gamsbin, args.sim, tuples, simdir, resdir, args.max_procs, args.resume)
    

if __name__ == "__main__":
        main()
