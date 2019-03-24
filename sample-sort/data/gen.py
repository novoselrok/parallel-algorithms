import numpy as np
import sys

if __name__ == "__main__":
    arr = np.arange(int(sys.argv[1]), dtype=np.int)
    np.random.shuffle(arr)
    np.savetxt('test_' + str(arr.shape[0]) + '.txt', arr, delimiter=' ', fmt='%d')
