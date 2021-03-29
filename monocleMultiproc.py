#!/usr/bin/python

import collections, re, math
import peak, peptideEnvelope
from pyteomics import mzxml
from multiprocessing import Process, Pool, Manager
from datetime import datetime


def getNearbyScans(reader, scanNumber, nTotScans, options):
    scans = []
    nScans = options["number_nearby_scans"]
    MAX_SCAN_NUM = nTotScans

    # Gets nearby scans around the MS1 scan given by scanNumber
    # 1. Flanking scans before the precursor scan
    n = 0
    i = scanNumber - 1
    while n < nScans:
        if i > 0:
            if reader[str(i)]["msLevel"] == 1:
                scans.append(reader[str(i)])
                n += 1
        else:
            break
        i -= 1
    scans.reverse()
    scans.append(reader[str(scanNumber)])

    # 2. Flanking scans after the precursor scan
    n = 0
    i = scanNumber + 1
    while n < nScans:
        if i <= MAX_SCAN_NUM:
            if reader[str(i)]["msLevel"] == 1:
                scans.append(reader[str(i)])
                n += 1
        else:
            break
        i += 1

    return scans


def dotProduct(a, b):
    num, den1, den2 = 0, 0, 0
    for i in range(min(len(a), len(b))):
        num += a[i] * b[i]
        den1 += a[i] ** 2
        den2 += b[i] ** 2
    if num == 0:
        return 0
    else:
        return num / math.sqrt(den1 * den2)


def monocle(precScan, nearbyScans, options):
    tol = options["tolerance"]
    tolUnit = options["tolerance_unit"]
    # Re-assign the precursor m/z in the precursor scan
    precMz = precScan["precMz"]
    idx = peak.match(precScan, precMz, 50, tolUnit)
    if idx >= 0:
        precMz = precScan["m/z array"][idx]
        precScan["precMz"] = precMz

    # Charge detection
    if precScan["precZ"] is not None:
        chargeLow = precScan["precZ"]
        chargeHigh = precScan["precZ"]
    else:   # Precursor charge is not determined
        # Check presumable charges from 1 to 5
        chargeLow = 2
        chargeHigh = 5

    bestCharge = 0
    bestScore = -1
    bestIndex = 0
    bestPeakMzs = None
    bestPeakIntensities = None

    for charge in range(chargeLow, chargeHigh + 1):
        mass = precMz * charge   # mass = [M + charge * H]
        isotopeRange = peak.isotopeRange(mass)

        # Generate expected relative intensities, and extract envelopes
        expected = peptideEnvelope.getTheoreticalEnvelope(precMz, charge, isotopeRange)
        expected = expected / max(expected)
        envelope = peptideEnvelope.extract(nearbyScans, precMz, charge, isotopeRange)

        # Get the best match using dot product
        for i in range(isotopeRange["isotopes"] - isotopeRange["compareSize"] + 1):
            observed = envelope["avgIntensity"][i: min(i + len(expected), len(envelope["avgIntensity"]))]
            if max(observed) > 0:
                observed = observed / max(observed)
            observed = peptideEnvelope.scaleByPeakCount(observed, envelope, i)
            score = dotProduct(observed, expected)

            if score > bestScore * 1.05:
                bestScore = score
                if score > 0.1:
                    bestIndex = i + 1
                    bestCharge = charge
                    bestPeakMzs = envelope["mz"][bestIndex]
                    bestPeakIntensities = envelope["intensity"][bestIndex]

    # End charge for loop
    if len(bestPeakMzs) > 0:
        num = 0
        den = 0
        for i in range(len(bestPeakMzs)):
            num += bestPeakMzs[i] * bestPeakIntensities[i]
            den += bestPeakIntensities[i]
        precMz = num / den

    return precMz, bestCharge


def decharge(reader, ms2Num, ms1Num, nScans, options):
    # Get nearby scans
    nearbyScans = getNearbyScans(reader, ms1Num, nScans, options)

    # Run Monocole algorithm
    currScan = reader[str(ms2Num)]
    precMz = float(re.search(r"ms2 ([0-9.]+)\@", currScan["filterLine"]).group(1))
    precScan = reader[str(ms1Num)]
    precScan["precMz"] = precMz
    try:
        precScan["precZ"] = currScan["precursorMz"][0]["precursorCharge"]
    except KeyError:
        precScan["precZ"] = None
    precMz, precZ = monocle(precScan, nearbyScans, options)

    # res["ms2"].append(ms2Num)
    # res["precMz"].append(precMz)
    # res["charge"].append(precZ)
    return ms2Num



########
# Main #
########
if __name__ == "__main__":
    inputFile = "FTLD_Batch2_F76.ReAdW.mzXML" # 1~34910
    reader = mzxml.read(inputFile)

    # Options (i.e. parameters)
    options = {"number_nearby_scans": 6, "tolerance": 10, "tolerance_unit": "ppm"}

    # Get the relationship between MS2 and its precursor MS1 scan
    ms2ToMs1 = {}
    nScans = 0
    print("Reading scans")
    n = 0
    for spec in reader:
        n += 1
        nScans += 1
        msLevel = int(spec["msLevel"])
        scanNum = int(spec["num"])
        if msLevel == 1:
            precScanNum = scanNum
        if msLevel == 2:
            ms2ToMs1[scanNum] = precScanNum
        if n > 100:
            break
    print("Finished reading scans")

    pool = Pool(2)
    # m = Manager()
    # res = m.dict()
    # res = {"ms2": [], "precMz": [], "charge": []}
    # pool.starmap_async(decharge, [(reader, k, v, options, res) for k, v in ms2ToMs1.items()])
    res = pool.starmap_async(decharge, [(reader, k, v, n, options) for k, v in ms2ToMs1.items()])
    res.wait()
    pool.close()
    res = res.get()
    # pool.join()
    print()

# # Monocle
# for spec in reader:
#     msLevel = int(spec["msLevel"])
#     if msLevel == 1:
#         precScanNum = int(spec["num"])
#     if msLevel == 2:
#         # Preprocessing
#         precMz = float(re.search(r"ms2 ([0-9.]+)\@", spec["filterLine"]).group(1))
#         precScan = reader[str(precScanNum)]
#         precScan["precMz"] = precMz
#         try:
#             precScan["precZ"] = spec["precursorMz"][0]["precursorCharge"]
#         except KeyError:
#             precScan["precZ"] = None
#
#         # # Get nearby scans
#         # nearbyScans = getNearbyScans(reader, precScanNum, options)
#         #
#         # # Run Monocle algorithm
#         # precMz, precZ = monocle(precScan, nearbyScans, options)
