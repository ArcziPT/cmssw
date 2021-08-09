import sys
import subprocess
import os.path

def execute_job(step, jsonIn, jsonOut, geometry, meanMax=0.2, rmsMax=0.2, treshold = 0):
    print(f"execute {step} step")

    command = "cmsRun Analyzer/DiamondTimingAnalyzer/python/test.py"
    command += " validOOT=\"-1\""
    command += f" calibInput=\"{jsonIn}\""
    command += f" calibOutput=\"{jsonOut}\""
    command += f" geometryFile=\"{geometry}\""
    command += f" loopIndex=\"{step}\""
    command += f" treshold=\"{treshold}\""
    command += f" meanMax=\"{meanMax}\""
    command += f" rmsMax=\"{rmsMax}\""
    
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()

    if process.returncode != 0:
        print("CMSSW error")
        exit(1)


if __name__ == "__main__":
    geometry = f"Geometry.VeryForwardGeometry.geometryRPFromDD_{sys.argv[3]}_cfi"
    meanMax = float(sys.argv[4])
    rmsMax = float(sys.argv[5])

    if sys.argv[1] == "fixed":
        n = int(sys.argv[2])
        for i in range(n):
            jsonIn = f"calib_{i}.json"
            jsonOut = f"calib_{i+1}.json"
            execute_job(i, jsonIn, jsonOut, geometry, meanMax, rmsMax)

    if sys.argv[1] == "diff":
        treshold = float(sys.argv[2])

        i = 0
        while True:
            jsonIn = f"calib_{i}.json"
            jsonOut = f"calib_{i+1}.json"
            execute_job(i, jsonIn, jsonOut, geometry, meanMax, rmsMax, treshold)

            if os.path.isfile("finish"):
                break

            i += 1