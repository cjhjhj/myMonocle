#!/usr/bin/python


import multiprocessing as mp
from datetime import datetime
import os
import math

def f(x):
    # print("value", x, "PID = ",os.getpid())
    # time.sleep(1)
    return x*x

if __name__ == '__main__':
    p = mp.Pool(2)
    startTime = datetime.now()
    print(p.map(f, range(0,10)))
    endTime = datetime.now()
    print("Total elapsed time", (endTime - startTime).total_seconds())