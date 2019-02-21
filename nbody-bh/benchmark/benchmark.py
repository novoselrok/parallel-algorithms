import os
import sys
import subprocess
from collections import defaultdict
import json
import numpy as np
import matplotlib.pyplot as plt

executables = [
    ('chapel', './main --filename={file} --nBodies={nbodies} --iterations={iterations}'),
    ('c', './main {file} {nbodies} {iterations}'),
    ('julia', 'julia -O3 --check-bounds=no main.jl {file} {iterations}')
]

threads = [
    1,
    2,
    4,
    8,
    16
]

nbodies = [
    10,
    20,
    40,
    80
]

prefixes = [
    'sequential',
    'sm'
]

N_REPEATS = 3

def get_filename(lang, n):
    file = None
    if lang == 'chapel':
        file = 'test' + str(n) + 'k.chpl.txt'
    else:
        file = 'test' + str(n) + 'k.txt'
    return '../../nbody/data/' + file

def set_sm_env(env, lang, threads):
    if lang == 'julia':
        env['JULIA_NUM_THREADS'] = str(threads)
        env['JL_NRETRIES'] = '5'
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

    for n in nbodies:
        print(n)
        results[n] = {}

        for prefix in prefixes:
            print(prefix)

            results[n][prefix] = {}

            if prefix in ('sequential', 'sm'):
                for lang, cmd_template in executables:
                    print(lang)

                    results[n][prefix][lang] = {}
                    os.chdir(os.path.join('..', prefix + '-' + lang))

                    _results = None
                    file = get_filename(lang, n)
                    cmd = cmd_template.format(file=file, nbodies=n * 1000, iterations=1)

                    if prefix == 'sequential':
                        my_env = set_sm_env(os.environ.copy(), lang, 1)
                        _results = []
                        for _ in range(N_REPEATS):
                            t = run_cmd(cmd, env=my_env)
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
                        
                    results[n][prefix][lang] = _results

    os.chdir(os.path.join('..', 'benchmark'))
    json.dump(results, open('results.json', 'w', encoding='utf-8'))

if __name__ == '__main__':
    main()
