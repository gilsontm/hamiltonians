from pysat.formula import WCNF
# from pysat.examples.fm import FM
from pysat.examples.rc2 import RC2

def main():
    wcnf = WCNF()
    wcnf.append([-1,-3], weight=(+16))
    wcnf.append([-1,-3,-4], weight=(+64))
    wcnf.append([-1,-3,-5], weight=(+128))
    wcnf.append([-1,-4], weight=(+64))
    wcnf.append([-1,-4,-5], weight=(+256))
    wcnf.append([-1,-5], weight=(+256))
    wcnf.append([-1,-2,-3], weight=(+64))
    wcnf.append([-1,-2,-3,-4], weight=(+256))
    wcnf.append([-1,-2,-3,-5], weight=(+512))
    wcnf.append([-1,-2,-4], weight=(+256))
    wcnf.append([-1,-2,-4,-5], weight=(+1024))
    wcnf.append([-1,-2,-5], weight=(+1024))
    wcnf.append([-2,-3], weight=(+64))
    wcnf.append([-2,-3,-4], weight=(+256))
    wcnf.append([-2,-3,-5], weight=(+512))
    wcnf.append([-2,-4], weight=(+256))
    wcnf.append([-2,-4,-5], weight=(+1024))
    wcnf.append([-2,-5], weight=(+1024))
    wcnf.append([-1,-3], weight=(+16))
    wcnf.append([-1,-4], weight=(+32))
    wcnf.append([-1,-5], weight=(+64))
    wcnf.append([-1,-2,-3], weight=(+64))
    wcnf.append([-1,-2,-4], weight=(+128))
    wcnf.append([-1,-2,-5], weight=(+256))
    wcnf.append([-1,-3], weight=(+16))
    wcnf.append([-1,-3,-4], weight=(+64))
    wcnf.append([-1,-3,-5], weight=(+128))
    wcnf.append([-1,-4], weight=(+64))
    wcnf.append([-1,-4,-5], weight=(+256))
    wcnf.append([-1,-5], weight=(+256))
    wcnf.append([-2,-3], weight=(+64))
    wcnf.append([-2,-4], weight=(+128))
    wcnf.append([-2,-5], weight=(+256))
    wcnf.append([-2,-3], weight=(+32))
    wcnf.append([-2,-3,-4], weight=(+128))
    wcnf.append([-2,-3,-5], weight=(+256))
    wcnf.append([-2,-4], weight=(+128))
    wcnf.append([-2,-4,-5], weight=(+512))
    wcnf.append([-2,-5], weight=(+512))
    wcnf.append([-1,-2], weight=(+16))
    wcnf.append([-3,-4], weight=(+16))
    wcnf.append([-3,-5], weight=(+32))
    wcnf.append([-4,-5], weight=(+64))
    wcnf.append([+1,-3], weight=(+248))
    wcnf.append([+1,-4], weight=(+496))
    wcnf.append([+1,-5], weight=(+992))
    wcnf.append([+2,-3], weight=(+496))
    wcnf.append([+2,-4], weight=(+992))
    wcnf.append([+2,-5], weight=(+1984))
    wcnf.append([-1], weight=(+4))
    wcnf.append([-2], weight=(+16))
    wcnf.append([-3], weight=(+4))
    wcnf.append([+3], weight=(+248))
    wcnf.append([+3], weight=(+496))
    wcnf.append([-4], weight=(+16))
    wcnf.append([+4], weight=(+496))
    wcnf.append([+4], weight=(+992))
    wcnf.append([-5], weight=(+64))
    wcnf.append([+5], weight=(+992))
    wcnf.append([+5], weight=(+1984))
    wcnf.append([+1], weight=(+128))
    wcnf.append([+2], weight=(+256))
    wcnf.append([+3], weight=(+128))
    wcnf.append([+4], weight=(+256))
    wcnf.append([+5], weight=(+512))

    # wcnf.append([-1,-3], weight=(-248))
    # wcnf.append([-1,-4], weight=(-496))
    # wcnf.append([-1,-5], weight=(-992))
    # wcnf.append([-2,-3], weight=(-496))
    # wcnf.append([-2,-4], weight=(-992))
    # wcnf.append([-2,-5], weight=(-1984))
    # wcnf.append([-1], weight=(-128))
    # wcnf.append([-2], weight=(-256))
    # wcnf.append([-3], weight=(-128))
    # wcnf.append([-4], weight=(-256))
    # wcnf.append([-5], weight=(-512))
    # wcnf.normalize_negatives([
    #     [[-1,-3], -248],
    #     [[-1,-4], -496],
    #     [[-1,-5], -992],
    #     [[-2,-3], -496],
    #     [[-2,-4], -992],
    #     [[-2,-5], -1984],
    #     [[-1], -128],
    #     [[-2], -256],
    #     [[-3], -128],
    #     [[-4], -256],
    #     [[-5], -512],])

    with RC2(wcnf) as solver:
        for model in solver.enumerate():
            literals = [f"{literal:+}" for literal in model]
            print(f"model [{','.join(literals)}] has cost {solver.cost:5d} (actual {solver.cost - 5464:5d})")



if __name__ == "__main__":
    main()