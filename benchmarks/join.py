import os
import json
import re

import numpy as np
import pandas as pd

# results-(.*)-(\d*).json
def main():
    json_files = [f for f in os.listdir('.') if f.endswith('.json')]
    print(json_files)
    
    benchmarks = []
    for json_file in json_files:
        with open(json_file, encoding='utf-8') as f:
            benchmarks.extend(json.load(f))
    
    for benchmark in benchmarks:
        benchmark['time'] = np.array(benchmark['times']).min()
        del benchmark['times']

    pd.DataFrame(benchmarks).to_csv('benchmarks.csv', index=False)


if __name__ == "__main__":
    main()

