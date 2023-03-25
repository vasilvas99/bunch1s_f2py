import numpy as np
from typing import Dict
import pandas as pd


def l_stat(Y: np.ndarray, tk: float, bdef: float) -> Dict:
    w = h = stime = dyi = avtd = tw = ntd = 0.0  #saving space!
    M = len(Y)

    mind = bdef
    mintd = M
    maxbd = 0.0
    maxtd = 0.0

    dy = np.zeros(M)

    # ...initializations...

    # Compute the interstep distances, the minimal 1 and
    # some other quantities that characterize the bunching:
    for i in range(1, M):
        dyi = Y[i] - Y[i - 1]
        dy[i] = dyi

        if dyi <= bdef:
            # counting the number of distances in bunches ...
            h += 1

            # ...and the total bunch width
            w += dyi

            if dyi < mind:
                mind = dyi

            if dyi > maxbd:
                maxbd = dyi
        else:
            tw += dyi
            ntd += 1

            if dyi < mintd:
                mintd = dyi

            if dyi > maxtd:
                maxtd = dyi

    dyi = Y[0] - Y[M - 1] + M  # take care of the PBC
    dy[0] = dyi

    if dyi <= bdef:
        w += dyi
        h += 1

        if dyi < mind:
            mind = dyi

        if dyi > maxbd:
            maxbd = dyi
    else:
        tw += dyi
        ntd += 1

        if dyi < mintd:
            mintd = dyi

        if dyi > maxtd:
            maxtd = dyi

    avtd = tw / ntd

    if maxbd > bdef:
        print('Alert:maxbd!')
        return 0

    if mintd < bdef:
        print('Alert:mintd!')
        return 0

    # Counting the terraces and their width:
    Tno = 0

    if dy[0] > bdef:
        Tno = 1

    for i in range(1, M):
        dyi = dy[i]

        if dyi > bdef and dy[i - 1] <= bdef:
            Tno += 1

    # check if the terrace in the beginning spannes the boundaries
    if dy[0] > bdef and dy[M - 1] > bdef:
        Tno -= 1

    if Tno > 0:
        stime = float(M) / float(Tno)
        tw = tw / float(Tno)
        h = h / float(Tno)
        w = w / float(Tno)

        statistics = {
            "dy": dy,
            "M": M,
            "bdef": bdef,
            "tk": tk,
            "h": h,
            "stime": stime,
            "w": w,
            "tw": tw,
            "avtd": avtd,
            "Tno": Tno,
            "mind": mind,
            "maxbd": maxbd,
            "mintd": mintd,
            "maxtd": maxtd
        }

        return statistics
    else:
        print("terrace alert")

# Heavy statistics, returns pandas dataframe with results
# TODO: add better colum names with final statistics
def h_stat(dy, bdef, max_bunch_size=300) -> pd.DataFrame:
    M = len(dy)

    bszi = np.zeros(max_bunch_size)
    absmin = np.zeros(max_bunch_size)
    fdmin = np.zeros(max_bunch_size)
    bwmin = np.zeros(max_bunch_size)
    bdmin = np.zeros(max_bunch_size)
    avbw = np.zeros((max_bunch_size, 2))
    avmn = np.zeros(max_bunch_size)
    avbd = np.zeros(max_bunch_size)
    avf = np.zeros(max_bunch_size)
    avl = np.zeros(max_bunch_size)

    bw = np.zeros(M // 2)
    bd = np.zeros(M // 2)
    mindi = np.zeros(M // 2)
    lbd = np.zeros(M // 2)
    fbd = np.zeros(M // 2)
    bsz = np.zeros(M // 2, dtype=np.int64)

    tbw = 0.0
    tmin = 0.0
    tla = 0.0

    bno = 0
    inc = 0
    dyi = dy[0]
    if dyi <= bdef:
        tla = dyi
        tmin = dyi
        tbw = dyi
        bno = 1
        fbd[0] = dyi
        inc = 2
        for Ip in range(1, M):
            dyip = dy[Ip]
            if dyip <= bdef:
                inc += 1
                if dyip < tmin:
                    tmin = dyip
                tbw += dyip
                tla = dyip
            else:
                break

    if bno == 1:
        bd[0] = tbw / float(inc - 1)
        bw[0] = tbw
        mindi[0] = tmin
        bsz[0] = inc
        lbd[0] = tla

    for i in range(1, M):
        dyi = dy[i]
        if dyi <= bdef:
            if dy[i - 1] > bdef:
                tla = dyi
                tmin = dyi
                tbw = dyi
                bno += 1
                fbd[bno - 1] = dyi
                inc = 2
                for Ip in range(i + 1, M):
                    dyip = dy[Ip]
                    if dyip <= bdef:
                        inc += 1
                        if dyip < tmin:
                            tmin = dyip
                        tbw += dyip
                        tla = dyip
                    else:
                        break
                bd[bno - 1] = tbw / float(inc - 1)
                bw[bno - 1] = tbw
                mindi[bno - 1] = tmin
                bsz[bno - 1] = inc
                lbd[bno - 1] = tla


#                  To account for the situation when the bunch crosses the
#                   boundary conditions:
    if dy[0] <= bdef and dy[M - 1] <= bdef:
        ibw = bsz[0] + bsz[bno - 1] - 1
        bw[0] += bw[bno - 1]
        bd[0] = bw[0] / float(ibw - 1)
        mindi[0] = min(mindi[0], mindi[bno - 1])
        bsz[0] = ibw
        fbd[0] = fbd[bno - 1]
        bsz[bno - 1] = 0
        bd[bno - 1] = 0.0
        bw[bno - 1] = 0.0
        mindi[bno - 1] = 0.0
        lbd[bno - 1] = 0.0
        bno -= 1

    if bno > 0:
        for i in range(0, bno + 1):
            bszi = bsz[i - 1]
            if bszi > max_bunch_size:
                print(f'Alert: Bunch Size exceeds {max_bunch_size=}')
                continue
            tmin = mindi[i - 1]
            avmn[bszi - 1] += tmin
            if tmin < absmin[bszi - 1]:
                absmin[bszi - 1] = tmin
            tmin = fbd[i - 1]
            avf[bszi - 1] += tmin
            if tmin < fdmin[bszi - 1]:
                fdmin[bszi - 1] = tmin
            tmin = bd[i - 1]
            avbd[bszi - 1] += tmin
            if tmin < bdmin[bszi - 1]:
                bdmin[bszi - 1] = tmin
            tmin = bw[i - 1]
            avbw[bszi - 1, 0] += tmin
            avbw[bszi - 1, 1] += 1
            if tmin < bwmin[bszi - 1]:
                bwmin[bszi - 1] = tmin
            avl[bszi - 1] += lbd[i - 1]

        # # LOOP 1. Prints too much, disabled for now
        # print("LOOP 1:")
        # for i in range(1, max_bunch_size):
        #     dyi = absmin[i-1]
        #     if dyi <= bdef:
        #         print(i, dyi, fdmin[i-1], bwmin[i-1], bdmin[i-1])

        # LOOP 2
        r = []
        print("LOOP 2:")
        for i in range(1, max_bunch_size):
            dyi = avbw[i - 1, 1]
            if dyi > 0:
                r.append([
                    i, avbw[i - 1, 0] / dyi, avbd[i - 1] / dyi,
                    avmn[i - 1] / dyi, avf[i - 1] / dyi, avl[i - 1] / dyi
                ])
                print("R1: ", i, avbw[i - 1, 0] / dyi, avbd[i - 1] / dyi)
                print("R2: ", i, avmn[i - 1] / dyi, avf[i - 1] / dyi,
                      avl[i - 1] / dyi)
    r = pd.DataFrame(r)
    r.columns = ["bunch no", "avbw[i - 1, 0] / dyi", "avbd[i - 1] / dyi",
                    "avmn[i - 1] / dyi", "avf[i - 1] / dyi", "avl[i - 1] / dyi"]
    return r

def test():
    bdef = 1.0

    with open("tests/test_input.txt", "r") as f:
        M = int(f.readline())
        Y = np.zeros(M)
        for i in range(Y.shape[0]):
            Y[i] = float(f.readline())

    output = l_stat(Y, 0.0, bdef)
    print("\n\n********Surface statistics (MSI)********\n\n")
    print(
        f'М = {M}\tbdef = {output["bdef"]}\tmind = {output["mind"]}\tmaxbd = {output["maxbd"]}\tmintd = {output["mintd"]}\tmaxtd = {output["maxtd"]}\tavtd = {output["avtd"]}'
    )
    print("\n\n********Surface statistics (MSII)********\n\n")
    h_stat(output["dy"], bdef)




if __name__ == "__main__":
    test()