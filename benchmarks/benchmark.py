import sys
import os
import subprocess
import itertools
import json
import time

BASE_DIR = os.path.dirname(os.path.realpath(__file__))
N_REPEATS = 5
NUM_WORKERS = [1, 2, 4, 8, 16]
INPUT_FILE_TEMPLATE = 'test_{size}{postfix}.txt'
SLEEP_TIME = 1

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
        'sizes': [64_000, 128_000, 256_000, 512_000],
        'consts': {
            'k': 256,
            'maxiter': 10
        },
        'executables': {
            'c': '{exepath}/main {inputfile} {k} {maxiter} {inputsize}',
            'julia': '/home/rok/julia-1.1.0/bin/julia -O3 --check-bounds=no --math-mode=fast {exepath}/main.jl {inputfile} {k} {maxiter}',
            'chapel': '{exepath}/main --filename={inputfile} --k={k} --maxIter={maxiter} --numPoints={inputsize}'
        }
    },
    'sample-sort': {
        'sizes': [10_000_000, 20_000_000, 40_000_000, 80_000_000],
        'consts': {},
        'executables': {
            'c': '{exepath}/main {inputfile} {inputsize}',
            'julia': '/home/rok/julia-1.1.0/bin/julia -O3 --check-bounds=no --math-mode=fast {exepath}/main.jl {inputfile}',
            'chapel': '{exepath}/main --filename={inputfile} --nkeys={inputsize}'
        }
    },
    'nbody-bh': {
        'sizes': [10_000, 20_000, 40_000, 80_000],
        'consts': {
            'iterations': 10
        },
        'executables': {
            'c': '{exepath}/main {inputfile} {inputsize} {iterations}',
            'julia': '/home/rok/julia-1.1.0/bin/julia -O3 --check-bounds=no --math-mode=fast {exepath}/main.jl {inputfile} {iterations}',
            'chapel': '{exepath}/main --filename={inputfile} --nBodies={inputsize} --iterations={iterations}'
        }
    }
}

def get_input_file_postfix(language, problem):
    if language == 'chapel' and problem != 'sample-sort':
        return '.chpl'
    return ''

def get_env(language, implementation, nworkers):
    language_env_vars = env_vars[language]
    env = {}
    env[language_env_vars['nthreads']] = str(nworkers)

    if language == 'julia':
        env['JL_NRETRIES'] = str(3)

    return env

def run_cmd(cmd, env=None):
    result = subprocess.run(cmd.split(' '), stdout=subprocess.PIPE, env=env)
    return float(result.stdout.decode('utf-8'))

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
            problem_input = {
                'exepath': path_to_implementation,
                'inputfile': input_file,
                'inputsize': size,
                **problem_inputs['consts']
            }
            cmd = problem_inputs['executables'][language].format(**problem_input)
            print(cmd)

            for nworkers in workers[implementation]:
                print(nworkers)

                problem_times = []
                for _ in range(N_REPEATS):
                    try:
                        problem_time = run_cmd(cmd, get_env(language, implementation, nworkers))
                        print(problem_time)
                        problem_times.append(problem_time)
                        time.sleep(SLEEP_TIME)
                    except e:
                        print(e)
            
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

