import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    nx = ny = int(sys.argv[1])
    filt = bool(int(sys.argv[2])) if len(sys.argv) > 2 else False

    x = range(nx,     0, -1)
    y = range(nx+ny, nx, -1)

    xy = [{xi, yi} for xi in x for yi in y]

    prod = []
    for i in range(len(xy)):
        prod.append([])
        for j in range(len(xy)):
            if not filt or j < i:
                value = 1 if (len(xy[i].union(xy[j])) == 4) else 0
            else:
                value = 0
            prod[-1].append(value)
            # print(f"{1 if len(i1.union(i2)) == 4 else 0:3d}", end=" ")
        # print("")
    prod = np.array(prod)

    soma = int(prod.sum())
    if not filt:
        soma /= 2

    print(f"{soma}")

    plt.imshow(prod, cmap='binary')
    plt.show()


if __name__ == "__main__":
    main()