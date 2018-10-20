import sys

import numpy as np
from sklearn.datasets import make_blobs


def main():
    rows = int(sys.argv[1])
    cols = int(sys.argv[2])
    k = int(sys.argv[3])
    fname = sys.argv[4]

    (X, y) = make_blobs(rows, cols, centers=k)

    np.savetxt(fname + '.txt', X, fmt='%E')
    np.savetxt(fname + '_labels.txt', y, fmt='%d')


if __name__ == '__main__':
    main()
