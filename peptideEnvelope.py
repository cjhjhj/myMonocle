import math, statistics
from scipy.stats import binom
import peak

def getTheoreticalEnvelope(precMz, charge, isotopeRange):
    nCarbons = estimateCarbons(precMz, charge)
    res = [0]
    for i in range(1, isotopeRange["compareSize"]):
        res.append(binom.pmf(i - 1, nCarbons, 0.011))
    return res


def estimateCarbons(mz, z):
    """
    Estimates the number of carbons in a peptide based only on its precursor m/z and charge.
    Original monocle: mz = 111; carbons = 5.1
    Senko et al. 1995: mz = 111.1254; carbons = 4.9384
    DKS Uniprot TREMBL 2019_08: mz = 110.3963; carbons = 4.9243
    """
    protonMass = 1.007276466812
    return math.floor((((mz - protonMass) * z) / 111) * 5.1)


def extract(scans, precMz, charge, isotopeRange):
    diff = 1.00286864    # Averagine difference(?)
    # diff = 1.003355 # C13 - C12 difference
    protonMass = 1.007276466812
    left, compareSize, nIsotopes = isotopeRange["left"], isotopeRange["compareSize"], isotopeRange["isotopes"]
    mzArray = [[] for i in range(nIsotopes)]
    intArray = [[] for i in range(nIsotopes)]
    for scan in scans:
        for i in range(nIsotopes):
            matchMz = precMz + (((i + left) * diff) / charge)
            idx = peak.match(scan, matchMz, 3, "ppm")   # Hard-coded tolerance in extract, 3ppm
            if idx >= 0:
                mz = scan["m/z array"][idx]
                intensity = scan["intensity array"][idx]
                mzArray[i].append(mz)
                intArray[i].append(intensity)
    nMax = 0
    avgIntensity = []
    for intensities in intArray:
        if len(intensities) > nMax:
            nMax = len(intensities)
        if len(intensities) > 0:
            avgIntensity.append(sum(intensities) / len(intensities))
        else:
            avgIntensity.append(0)

    output = {"mz": mzArray, "intensity": intArray, "avgIntensity": avgIntensity, "maxPeakCount": nMax}
    return output


def scaleByPeakCount(x, env, i):
    if env["maxPeakCount"] > 0:
        for j in range(i, i + len(x)):
            x[j - i] *= len(env["mz"][j]) / env["maxPeakCount"]
    return x
