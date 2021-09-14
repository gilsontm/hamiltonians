import numpy as np
import matplotlib.pyplot as plt


def number_of_bits_for_x(N):
    floor = np.floor(np.sqrt(N))
    floor -= 1 - (floor % 2)
    nx = np.ceil(np.log2(floor)) - 1
    return nx

def number_of_bits_for_y(N):
    ny = np.ceil(np.log2(np.floor(N / 3))) - 1
    return ny

def main():
    interval_max = 10 ** 6
    Ns = np.array(range(10, interval_max))
    Xs = number_of_bits_for_x(Ns)
    Ys = number_of_bits_for_y(Ns)
    XsYs = Xs + Ys
    Ls = np.log2(Ns)

    fig, ax = plt.subplots()
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.plot(Ns, XsYs, color="y")
    ax.plot(Ns, Xs, color="b")
    ax.plot(Ns, Ys, color="r")
    ax.plot(Ns, Ls, color="g")
    ax.set_xlabel("N")
    ax.set_ylabel("number of bits")
    ax.legend(["total bits", "bits for x", "bits for y", "log2(N)"])
    plt.show()



if __name__ == "__main__":
    main()