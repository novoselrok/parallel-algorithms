import numpy as np
from sklearn.cluster import KMeans

import time

def main():
    kmeans = KMeans(n_clusters=100, init='random', n_init=1, max_iter=10)

    X = np.loadtxt('../data/test_R256k_C16_k100.txt')
    print('Start simulation', X.shape)
    start = time.time()
    kmeans.fit(X)
    print(time.time() - start)
    print(kmeans.labels_)

if __name__ == '__main__':
    main()