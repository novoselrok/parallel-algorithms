import sys

import numpy as np
from sklearn.datasets import make_blobs


def format_tuple(tpl):
    return "({})".format(", ".join(["{:e}".format(num) for num in tpl]))


# How to run: python3 gen.py <filename> <rows> <cols> <k>
def main():
    rows = int(sys.argv[2])
    cols = int(sys.argv[3])
    k = int(sys.argv[4])
    name = sys.argv[1] + '_' + str(rows) + '_' + str(cols) + '_' + str(k)
    fname = name + '.txt'
    fname_chpl = name + '.chpl.txt'
    
    (X, y) = make_blobs(rows, cols, centers=k)

    # Save labels
    np.savetxt(name + '_labels.txt', y, fmt='%d')
    # Save regular matrix
    np.savetxt(fname, X, fmt='%E')
    
    with open(fname_chpl, 'w', encoding='utf-8') as f_chpl:
        for point in X:
            f_chpl.write(format_tuple(point) + '\n')

if __name__ == '__main__':
    main()
