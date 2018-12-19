import numpy as np

if __name__ == "__main__":
    arr = np.arange(5e7, dtype=np.int)
    np.random.shuffle(arr)
    np.savetxt('arr' + str(arr.shape[0]) + '.txt', arr, delimiter=' ', fmt='%d')
