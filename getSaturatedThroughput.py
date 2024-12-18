import math
import sys
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

## input -----------------------------------------------------------------------------------------
NSS = 2
NSD = 108
IMCS = 7
W = 16
MSMALL = 6
RSMALL = 6
D = 1500
TSYM = 3.6
DELTA = 0.1
## input -----------------------------------------------------------------------------------------

## setting ---------------------------------------------------------------------------------------
BASIC_ACCESS = 1
RTS_CTS_ACCESS = 2
## setting ---------------------------------------------------------------------------------------

## values of IEEE 802.11ac -----------------------------------------------------------------------
TI = 9
TSIFS = 16
TDIFS = 34
## values of IEEE 802.11ac -----------------------------------------------------------------------

## some values decided by other values -----------------------------------------------------------
NVLTF = (0,1,2,3,4,5,6,7,8)
M_R = ((2,1/2),(4,1/2),(4,3/4),(16,1/2),(16,3/4),(64,2/3),(64,3/4),(64,5/6),(256,3/4),(256,5/6)) 
T_CBACK = (68,44,36,32,28,28,28)
T_RTS_CTS = ((52,44),(36,32),(32,28),(28,28),(28,24),(24,24),(24,24))
## some values decided by other values -----------------------------------------------------------

## settings of graph -------------------------------------------------------------------------------
GCOLOR = ('black', 'red', 'blue', 'yellow', 'green')
GSHAPE = ('o', '*', '^', '+', 's')
## settings of graph -------------------------------------------------------------------------------

def getTdata():
    m, r = M_R[IMCS]
    inCeil = (16+8*AP*(38+AS*(14+D))) / (NSS*NSD*r*math.log2(m))
    return (36+NVLTF[NSS]*4) + TSYM * math.ceil(inCeil)

def getTs_c(mode, Tdata):
    Imcs = getValidImcs(IMCS)
    Tcback = T_CBACK[Imcs]
    if mode == BASIC_ACCESS:
        Ts = Tdata + Tcback + TSIFS + TDIFS + 2*DELTA
        Tc = Tdata + TDIFS + DELTA
    elif mode == RTS_CTS_ACCESS:
        Trts, Tcts = T_RTS_CTS[Imcs]
        Ts = Trts + Tcts + Tdata + Tcback + 3*TSIFS + TDIFS + 4*DELTA
        Tc = Trts + TDIFS + DELTA
    else:
        sys.exit(-1)
    return (Ts, Tc)

def getValidImcs(Imcs):
    if not isinstance(Imcs, int) or Imcs < 0 or Imcs > 9:
        sys.exit(-1)
    else:
        if 7 <= Imcs <= 9:
            Imcs = 6
        return Imcs

def func(x):
    u = min(MSMALL, RSMALL)
    return [2 - x[0] * (1 + W * ((((2*x[1])**u-1)/(2*x[1]-1) + (2*x[1])**u * (((2*x[1])**(RSMALL-u+1)-1)/(2*x[1]-1))) / ((x[1]**(RSMALL+1)-1)/(x[1]-1)))), 
            1 - x[1] - (1-x[0])**(n-1)]

def getTau_P():
    return fsolve(func, [0,0])

def getP(tau):
    Pi = (1-tau)**n
    Ps = n * tau * (1-tau)**(n-1)
    Pc = 1 - Pi - Ps
    return (Pi, Ps, Pc)

def getST(Pi_s_c, t):
    Pi, Ps, Pc = Pi_s_c
    Ts, Tc = t
    return (Ps*AS*AP*D) / (Pi*TI+Ps*Ts+Pc*Tc)

while True:
    try:
        input_str = input("input integer for As: ")
        AS = int(input_str)
        break
    except ValueError:
        print("invalid input!!")

while True:
    try:
        input_str = input("input integers for Aps splitting by comma: ")
        APs = [int(x) for x in input_str.split(',')]
        break
    except ValueError:
        print("invalid input!!")

while True:
    try:
        input_str = input("input integer for access mode (BASIC: 1, RTS/CTS: 2): ")
        mode = int(input_str)
        if mode not in (1,2):
            raise ValueError
        break
    except ValueError:
        print("invalid input!!")

for index, AP in enumerate(APs):
    n_values = []
    st_values = []

    for n in range(2, 51):
        Tdata = getTdata()
        Ts_c = getTs_c(mode, Tdata)
        Tau, p = getTau_P()
        Pi_s_c = getP(Tau)
        st = getST(Pi_s_c, Ts_c)

        n_values.append(n)
        st_values.append(st)

    label = f'Ap={AP}'
    plt.plot(n_values, st_values, label=label, color=GCOLOR[index], marker=GSHAPE[index])

plt.title(f'S vs n for As={AS}')
plt.xlabel('n')
plt.ylabel('Saturated Throughput')
plt.legend()

plt.show()