import os
import sys
import subprocess
from collections import defaultdict
import json
import numpy as np
import matplotlib.pyplot as plt

executables = [
    ('chapel', './main --filename={} --n_bodies={}'),
    ('openmp', './main {} {} 1 out.txt'),
    ('julia', 'julia main.jl {} 1')
]

data = [10000, 20000, 40000, 80000, 160000]

def main():
    results = {}
    for n_bodies in data:
        results[n_bodies] = defaultdict(list)
        thousands = n_bodies // 1000
        fname = '../data/test{}k.txt'.format(thousands)
        fname_chpl = '../data/test{}k.chpl.txt'.format(thousands)
        print(n_bodies)
        
        for folder, cmd in executables:
            os.chdir(os.path.join('..', folder))
            my_env = os.environ.copy()
            if folder == 'julia':
                my_env['JULIA_NUM_THREADS'] = '16'
            elif folder == 'chapel':
                my_env['CHPL_RT_NUM_THREADS_PER_LOCALE'] = '16'
            f = fname_chpl if folder == 'chapel' else fname
            print(folder)
            for _ in range(10):
                result = subprocess.run(cmd.format(f, n_bodies).split(' '), stdout=subprocess.PIPE, env=my_env)
                t = float(result.stdout.decode('utf-8'))
                print(t)
                results[n_bodies][folder].append(t)
    os.chdir(os.path.join('..', 'benchmark'))
    json.dump(results, open('results.json', 'w', encoding='utf-8'))


def plot():
    results = json.load(open('results.json', encoding='utf-8'))
    xlabels = list(map(int, results.keys()))
    langs = ['julia', 'openmp', 'chapel']

    for lang in langs:
        y = []
        yerr = []
        for n_bodies, times in results.items():
            lang_times = np.array(times[lang])
            mean = lang_times.mean()
            minimum = lang_times.min()
            maximum = lang_times.max()
            y.append(mean)
            yerr.append([mean - minimum, maximum - mean])
        plt.errorbar(xlabels, y, yerr=np.array(yerr).T, label=lang)
    
    plt.xlabel('number of bodies')
    plt.ylabel('time of 1 iteration (seconds)')
    plt.legend()
    plt.savefig('nbodies.png')
    # plt.show()


if __name__ == '__main__':
    if sys.argv[1] == 'plot':
        plot()
    else:
        main()
