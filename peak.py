import math, numpy as np
from datetime import datetime

def withinError(theoretical, observed, tol, tolUnit="ppm"):
    # try:
    #     tol = options["tolerance"]
    # except KeyError:
    #     tol = 10    # Default 10 ppm
    # try:
    #     tolUnit = options["tolerance_unit"]
    # except KeyError:
    #     tolUnit = "ppm"
    if tolUnit.lower() == "ppm":
        return abs(theoretical - observed) / theoretical * 1e6 < tol
    if tolUnit.lower() == "da":
        return abs(theoretical - observed) < tol
    return False


def match(scan, targetMz, tol, tolUnit="ppm"):
    # try:
    #     tol = options["tolerance"]
    # except KeyError:
    #     tol = 10    # Default 10 ppm
    # try:
    #     tolUnit = options["tolerance_unit"]
    # except KeyError:
    #     tolUnit = "ppm"
    mzArray = scan["m/z array"]
    nPeaks = len(mzArray)
    i = np.argmin(abs(mzArray - targetMz))

    foundNext = False
    errorNext = 0
    if i < nPeaks:
        foundNext = True
        errorNext = abs(mzArray[i] - targetMz)
    foundPrevious = False
    errorPrevious = 0
    if i > 0:
        foundPrevious = True
        errorPrevious = abs(mzArray[i - 1] - targetMz)

    if not foundNext and not foundPrevious:
        return -1
    if not foundNext or (foundPrevious and errorPrevious < errorNext):
        if withinError(mzArray[i - 1], targetMz, tol, tolUnit):
            return i - 1
    elif not foundPrevious or (foundNext and errorNext < errorPrevious):
        if withinError(mzArray[i], targetMz, tol, tolUnit):
            return i
    return -1


def mostIntenseIndex(scan, targetMz, tol, tolUnit="ppm"):
    # try:
    #     tol = options["tolerance"]
    # except KeyError:
    #     tol = 10    # Default 10 ppm
    # try:
    #     tolUnit = options["tolerance_unit"]
    # except KeyError:
    #     tolUnit = "ppm"
    if tolUnit.lower() == "ppm":
        lL = targetMz - tol * targetMz / 1e6
        uL = targetMz + tol * targetMz / 1e6
    elif tolUnit.lower() == "da":
        lL = targetMz - tol
        uL = targetMz + tol
    mzArray = scan["m/z array"]
    i = nearestIndex(mzArray, lL)
    if mzArray[i] > lL:
        i += 1

    maxIntensity = 0
    bestIndex = -1
    intensityArray = scan["intensity array"]
    for i in range(len(intensityArray)):
        if mzArray[i] > uL:
            break
        else:
            if intensityArray[i] > maxIntensity:
                maxIntensity = intensityArray[i]
                bestIndex = i
    return bestIndex


def nearestIndex(peaks, target):
    # Binary search
    low, mid, high = 0, 0, len(peaks) - 1
    while True:
        if low == high:
            return low
        mid = math.floor((high + low) / 2)
        if peaks[mid] < target:
            low = mid + 1
        else:
            high = mid


def isotopeRange(mass):
    res = {}
    if mass > 2900:
        res["isotopes"] = 14
        res["left"] = -7
        res["compareSize"] = 7
    elif mass > 1200:
        res["isotopes"] = 10
        res["left"] = -5
        res["compareSize"] = 5
    else:
        res["isotopes"] = 7
        res["left"] = -3
        res["compareSize"] = 4
    res["monoisotopicIndex"] = -res["left"]
    return res

"""
    def __init__(self, targetMz, scan, tol, tolUnit = "ppm"):
        self.targetMz = targetMz
        self.scan = scan
        self.tol = tol
        self.tolUnit = tolUnit

    def withinError(self, theoretical, observed, tol, tolUnit):
        if tolUnit.lower() == "ppm":
            return abs(theoretical - observed) / theoretical * 1e6 < tol
        if tolUnit.lower() == "da":
            return abs(theoretical - observed) < tol
        return False

    def match(self):
        mzArray = self.scan["m/z array"]
        nPeaks = len(mzArray)
        i = self.nearestIndex(mzArray, self.targetMz)
        foundNext = False
        errorNext = 0
        if i < nPeaks:
            foundNext = True
            errorNext = abs(mzArray[i] - self.targetMz)
        foundPrevious = False
        errorPrevious = 0
        if i > 0:
            foundPrevious = True
            errorPrevious = abs(mzArray[i] - self.targetMz)

        if not foundNext and not foundPrevious:
            return -1
        if not foundNext or (foundPrevious and errorPrevious < errorNext):
            if self.withinError(mzArray[i - 1], self.targetMz, self.tol, self.tolUnit):
                return i - 1
        elif not foundPrevious or (foundNext and errorNext < errorPrevious):
            if self.withinError(mzArray[i], self.targetMz, self.tol, self.tolUnit):
                return i
        return -1

    def mostIntenseIndex(self):
        if self.tolUnit.lower() == "ppm":
            lL = self.targetMz - self.tol * self.targetMz / 1e6
            uL = self.targetMz + self.tol * self.targetMz / 1e6
        elif self.tolUnit.lower() == "da":
            lL = self.targetMz - self.tol
            uL = self.targetMz + self.tol
        mzArray = self.scan["m/z array"]
        i = self.nearestIndex(mzArray, lL)
        if mzArray[i] > lL:
            i += 1

        maxIntensity = 0
        bestIndex = -1
        intensityArray = self.scan["intensity array"]
        for i in range(len(intensityArray)):
            if mzArray[i] > uL:
                break
            else:
                if intensityArray[i] > maxIntensity:
                    maxIntensity = intensityArray[i]
                    bestIndex = i
        return bestIndex

    def nearestIndex(self, peaks, target):
        # Binary search
        low, mid, high = 0, 0, len(peaks) - 1
        while True:
            if low == high:
                return low
            mid = math.floor((high + low) / 2)
            if peaks[mid] < target:
                low = mid + 1
            else:
                high = mid

class isotopeRange:

    def __init__(self, mass):
        if mass > 2900:
            self.isotopes = 14
            self.left = -7
            self.compareSize = 7
        elif mass > 1200:
            self.isotopes = 10
            self.left = -5
            self.comareSize = 5
        else:
            self.isotopes = 7
            self.left = -3
            self.compareSize = 4
        self.monoisotopicIndex = -self.left

"""


