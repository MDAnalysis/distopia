import numpy as np
import sys


def main(ncoords):
    box = np.array([10, 20, 30], dtype=np.float32)

    coords = np.random.random((ncoords, 3))

    with open('data.txt', 'w') as f:
        f.write('{}\n'.format(ncoords))
        f.write('{} {} {}\n'.format(box[0], box[1], box[2]))
        for c in coords:
            f.write('{} {} {}\n'.format(c[0], c[1], c[2]))

if __name__ == '__main__':
    try:
        natoms = int(sys.argv[1])
    except:
        print("Usage: generate_coords.py <n>")
        exit(1)

    main(natoms)
