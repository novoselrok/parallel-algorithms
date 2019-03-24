import os
import sys
import subprocess
from collections import defaultdict
import json
import numpy as np
import matplotlib.pyplot as plt

executables = [
    # ('chapel', './main --filename={file} --k={k} --maxIter={maxiter} --numPoints={npoints}'),
    # ('c', './main {file} {k} {maxiter} {npoints}'),
    ('julia', 'julia main.jl {file} {k} {maxiter}')
]

datafiles = [
    '../data/test_64000'
]

threads = [
    1,
    2,
    4,
    8,
    16
]

nclusters = [
    256,
    # 512
]

prefixes = [
    'sequential',
    # 'sm'
]

N_REPEATS = 3

def get_filename(lang, filename):
    if lang == 'chapel':
        return filename + '.chpl.txt'
    return filename + '.txt'


def set_sm_env(env, lang, threads):
    if lang == 'julia':
        env['JULIA_NUM_THREADS'] = str(threads)
    elif lang == 'chapel':
        env['CHPL_RT_NUM_THREADS_PER_LOCALE'] = str(threads)
    elif lang == 'c':
        env['OMP_NUM_THREADS'] = str(threads)
    return env


def run_cmd(cmd, env=None):
    result = subprocess.run(cmd.split(' '), stdout=subprocess.PIPE, env=env)
    return float(result.stdout.decode('utf-8'))


def main():
    results = {}

    for datafile in datafiles:
        results[datafile] = {}

        for prefix in prefixes:
            print(prefix)

            results[datafile][prefix] = {}

            if prefix in ('sequential', 'sm'):
                for lang, cmd_template in executables:
                    print(lang)

                    results[datafile][prefix][lang] = {}
                    filename = get_filename(lang, datafile)
                    os.chdir(os.path.join('..', prefix + '-' + lang))

                    for k in nclusters:
                        print(k)
                        _results = None
                        cmd = cmd_template.format(file=filename, k=k, maxiter=1, npoints=100000)

                        if prefix == 'sequential':
                            _results = []
                            for _ in range(N_REPEATS):
                                t = run_cmd(cmd)
                                print(t)
                                _results.append(t)
                        else:
                            _results = defaultdict(list)
                            for nthreads in threads:
                                print(nthreads)
                                my_env = set_sm_env(os.environ.copy(), lang, nthreads)
                                for _ in range(N_REPEATS):
                                    t = run_cmd(cmd, env=my_env)
                                    print(t)
                                    _results[nthreads].append(t)
                        
                        results[datafile][prefix][lang][k] = _results

    os.chdir(os.path.join('..', 'benchmark'))
    json.dump(results, open('results.json', 'w', encoding='utf-8'))

def plot():
    results = json.load(open('results.json', encoding='utf-8'))
    pass

if __name__ == '__main__':
    main()
