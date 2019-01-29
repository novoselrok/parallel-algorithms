import os
import sys
import subprocess
from collections import defaultdict
import json
import numpy as np
import matplotlib.pyplot as plt

executables = [
    ('chapel', './main --filename={file} --nkeys={nkeys}'),
    ('c', './main {file} {nkeys}'),
    ('julia', 'julia main.jl {file}')
]

threads = [
    1,
    2,
    4,
    8,
    16
]

nkeys = [
    # int(1e5),
    # int(1e6),
    int(1e7)
]

prefixes = [
    'sequential',
    'sm'
]

N_REPEATS = 3

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

    for n in nkeys:
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
                    file = '../data/arr{}.txt'.format(n)
                    cmd = cmd_template.format(file=file, nkeys=n)

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
