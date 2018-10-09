import sys
import random

MIN_POSITION = -6E05
MAX_POSITION = 6E05
MIN_VELOCITY = -3E03
MAX_VELOCITY = 3E03
MIN_MASS = 1E19
MAX_MASS = 1E21

# How to run: python3 gen.py <output filename> <# bodies> -> python3 gen.py out1000.txt 1000
def main():
    fname = sys.argv[1]
    nbodies = int(sys.argv[2])

    with open(fname, 'w', encoding='utf-8') as f:
        f.write(str(nbodies) + '\n')
        for _ in range(nbodies):
            body = [
                random.uniform(MIN_POSITION, MAX_POSITION),
                random.uniform(MIN_POSITION, MAX_POSITION),
                random.uniform(MIN_POSITION, MAX_POSITION),

                random.uniform(MIN_VELOCITY, MAX_VELOCITY),
                random.uniform(MIN_VELOCITY, MAX_VELOCITY),
                random.uniform(MIN_VELOCITY, MAX_VELOCITY),

                random.uniform(MIN_MASS, MAX_MASS)
            ]

            f.write(" ".join(map(lambda x: '{0:E}'.format(x), body)) + '\n')

if __name__ == '__main__':
    main()
