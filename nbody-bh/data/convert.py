import numpy as np
import sys

def format_tuple(tpl):
    return "({})".format(", ".join(["{:e}".format(num) for num in tpl]))

if __name__ == "__main__":
    arg = sys.argv[1]
    name = arg.split('.')[0]
    a = np.loadtxt(arg)

    with open(name + '.chpl.txt', 'w', encoding='utf-8') as f:
        for row in a:
            f.write(format_tuple(row[:3]))
            f.write(' ')
            f.write(format_tuple(row[3:6]))
            f.write(' ')
            f.write(str(row[-1]))
            f.write('\n')
