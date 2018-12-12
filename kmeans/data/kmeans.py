import sys

from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import numpy as np

# Sanity check
def main():
    X = np.loadtxt(sys.argv[1])
    labels = np.loadtxt(sys.argv[2], dtype=int)

    kmeans = KMeans(n_clusters=int(sys.argv[3]), init='random', n_init=1, max_iter=int(sys.argv[4]))
    kmeans.fit(X)
    print(adjusted_rand_score(labels, kmeans.labels_))


if __name__ == "__main__":
    main()