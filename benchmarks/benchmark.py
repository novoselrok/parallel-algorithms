import sys
import os
import subprocess
import itertools
import json
import time
import socket
import numpy as np

BASE_DIR = os.path.dirname(os.path.realpath(__file__))
N_REPEATS = 3
NUM_WORKERS = [1, 2, 4, 8, 16]
INPUT_FILE_TEMPLATE = 'test_{size}{postfix}.txt'
SLEEP_TIME = 1
IS_LOCAL = socket.gethostname() == 'box'

workers = {
    'sequential': [1],
    'sm': [1, 2, 4, 8, 16],
    'dm': [1, 2, 4]
}

env_vars = {
    'c': {
        'nthreads': 'OMP_NUM_THREADS'
    },
    'chapel': {
        'nthreads': 'CHPL_RT_NUM_THREADS_PER_LOCALE'
    },
    'julia': {
        'nthreads': 'JULIA_NUM_THREADS'
    }
}

inputs = {
    'kmeans': {
        'sizes': [64000, 128000, 256000, 512000],
        'consts': {
            'k': 256,
            'maxiter': 10
        },
        'executables': {
            'c': '{cmdoptions} {exepath}/main {inputfile} {k} {maxiter} {inputsize}',
            'julia': '{cmdoptions} -O3 --check-bounds=no --math-mode=fast {exepath}/main.jl {inputfile} {k} {maxiter}',
            'chapel': '{exepath}/main {cmdoptions} --filename={inputfile} --k={k} --maxIter={maxiter} --numPoints={inputsize}'
        }
    },
    'sample-sort': {
        'sizes': [10000000, 20000000, 40000000, 80000000],
        'consts': {},
        'executables': {
            'c': '{cmdoptions} {exepath}/main {inputsize}',
            'julia': '{cmdoptions} -O3 --check-bounds=no --math-mode=fast {exepath}/main.jl {inputsize}',
            'chapel': '{exepath}/main {cmdoptions} --nkeys={inputsize}'
        }
    },
    'nbody-bh': {
        'sizes': [10000, 20000, 40000, 80000],
        'consts': {
            'iterations': 10
        },
        'executables': {
            'c': '{cmdoptions} {exepath}/main {inputfile} {inputsize} {iterations}',
            'julia': '{cmdoptions} -O3 --check-bounds=no --math-mode=fast {exepath}/main.jl {inputfile} {iterations}',
            'chapel': '{exepath}/main {cmdoptions} --filename={inputfile} --nBodies={inputsize} --iterations={iterations}'
        }
    }
}

with open('hostfile', encoding='utf-8') as f:
    GASNET_SSH_SERVERS = ' '.join(f.read().split('\n'))

def get_input_file_postfix(language, problem):
    if language == 'chapel' and problem != 'sample-sort':
        return '.chpl'
    return ''

def get_env(language, problem, implementation, nworkers):
    language_env_vars = env_vars[language]
    env = {}
    if implementation != 'dm':
        env[language_env_vars['nthreads']] = str(nworkers)

    if language == 'julia':
        env['JL_NRETRIES'] = str(3) if problem != 'sample-sort' else str(1)
    
    if language == 'chapel' and implementation == 'dm':
        env['CHPL_TARGET_ARCH'] = 'native'
        env['GASNET_SPAWNFN'] = 'S'
        env['GASNET_SSH_CMD'] = 'ssh'
        env['GASNET_SSH_OPTIONS'] = '-x'
        env['GASNET_SSH_SERVERS'] = GASNET_SSH_SERVERS

    return env

def get_cmd_options(language, problem, implementation, nworkers):
    if language == 'julia':
        if not IS_LOCAL and implementation == 'dm':
            return '/home/guest/roknovosel/julia-1.1.0/bin/julia --machine-file hostfile.{}'.format(nworkers)
        else:
            return '/home/rok/julia-1.1.0/bin/julia'
    elif language == 'c' and implementation == 'dm':
        return 'mpirun.mpich -np {} --hostfile hostfile'.format(nworkers)
    elif language == 'chapel' and implementation == 'dm':
        return '-nl {}'.format(nworkers)

    return ''

def run_cmd(cmd, language, problem, implementation, env=None):
    cmd_list = [arg for arg in cmd.strip().split(' ') if arg]
    result = subprocess.run(cmd_list, stdout=subprocess.PIPE, env=env)
    output = result.stdout.decode('utf-8').strip()
    if language == 'julia' and implementation == 'dm':
        return float(output.split('\n')[-1])
    else:
        return float(output)

def main(args):
    print(BASE_DIR)
    problems = args[0].split(',') # kmeans,sample-sort,nbody-bh
    languages = args[1].split(',') # c,chapel,julia
    implementations = args[2].split(',') # sequential,sm,dm
    
    results = []
    for (problem, language, implementation) in itertools.product(problems, languages, implementations):
        print((problem, language, implementation))

        path_to_implementation = os.path.join(problem, implementation + '-' + language)

        problem_inputs = inputs[problem]

        for size in problem_inputs['sizes']:
            input_file = os.path.join(problem, 'data', INPUT_FILE_TEMPLATE.format(size=size, postfix=get_input_file_postfix(language, problem)))

            for nworkers in workers[implementation]:
                problem_input = {
                    'cmdoptions': get_cmd_options(language, problem, implementation, nworkers),
                    'exepath': path_to_implementation,
                    'inputfile': input_file,
                    'inputsize': size,
                    **problem_inputs['consts']
                }
                cmd = problem_inputs['executables'][language].format(**problem_input)
                print(cmd, "nworkers=", nworkers)

                problem_times = []
                for _ in range(N_REPEATS):
                    try:
                        if language == 'chapel' and problem == 'sample-sort' and implementation != 'sequential':
                            times = []
                            for initial_seed in range(100):
                                new_cmd = cmd + ' --initialSeed=' + str(initial_seed)
                                t = run_cmd(new_cmd, language, problem, implementation, get_env(language, problem, implementation, nworkers))
                                times.append(t)
                            problem_time = np.array(times).mean()
                        else:
                            problem_time = run_cmd(cmd, language, problem, implementation, get_env(language, problem, implementation, nworkers))

                        print(problem_time)
                        problem_times.append(problem_time)
                        time.sleep(SLEEP_TIME)
                    except Exception as e:
                        print("error: " + str(e))

                results.append({
                    'problem': problem,
                    'language': language,
                    'implementation': implementation,
                    'size': size,
                    'nworkers': nworkers,
                    'times': problem_times
                })

    output_filename = '-'.join(['results'] + args + [str(int(time.time()))]) + '.json'
    with open(output_filename, 'w', encoding='utf-8') as f:
        json.dump(results, f)


if __name__ == "__main__":
    main(sys.argv[1:])

