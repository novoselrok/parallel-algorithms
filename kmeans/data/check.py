import sys

from sklearn.metrics import adjusted_rand_score
import numpy as np

def main():
    true_labels = np.loadtxt(sys.argv[1], dtype=int)
    pred_labels = np.loadtxt(sys.argv[2], dtype=int)

    print(adjusted_rand_score(true_labels, pred_labels))


if __name__ == "__main__":
    main()
