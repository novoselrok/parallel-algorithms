import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
import numpy as np

np.random.seed(0)

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12
plt.style.use('ggplot')

colors = ['#4EACC5', '#FF9C34', '#4E9A06']

if __name__ == "__main__":
    (X, y) = make_blobs(1000, 2, centers=3, cluster_std=2, random_state=0)

    fig = plt.figure(figsize=(8, 3))
    for idx, mi in enumerate([1, 5, 10]):
        kmeans = KMeans(init='random', n_clusters=3, n_init=1, random_state=0, max_iter=mi)
        kmeans.fit(X)

        ax = fig.add_subplot(1, 3, idx + 1)
        for k, col in zip(range(3), colors):
            my_members = kmeans.labels_ == k
            cluster_center = kmeans.cluster_centers_[k]
            ax.plot(X[my_members, 0], X[my_members, 1], 'w', markerfacecolor=col, marker='.', linestyle='None')
            ax.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col, markeredgecolor='k', markersize=6, linestyle='None')
        ax.set_title('Iteration: ' + str(mi), fontsize=9)
        ax.set_xticks(())
        ax.set_yticks(())

    fig.suptitle('k-means iterations')
    fig.tight_layout()
    fig.subplots_adjust(top=0.85)
    # plt.show()
    plt.savefig('figures/kmeans-iterations-eng.pdf', transparent=False, dpi=80)
    plt.clf()
