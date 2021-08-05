import sys
import subprocess
import os.path

def execute_job(step, jsonIn, jsonOut, geometry, treshold = 0):
    print(f"execute {step} step")

    command = "cmsRun Analyzer/DiamondTimingAnalyzer/python/test.py"
    command += " validOOT=\"-1\""
    command += f" calibInput=\"{jsonIn}\""
    command += f" calibOutput=\"{jsonOut}\""
    command += f" geometryFile=\"{geometry}\""
    command += f" loopIndex=\"{step}\""
    command += f" treshold=\"{treshold}\""
    
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()


if __name__ == "__main__":
    if sys.argv[1] == "fixed":
        geometry = sys.argv[3]

        n = int(sys.argv[2])
        for i in range(n):
            jsonIn = f"calib_{i}.json"
            jsonOut = f"calib_{i+1}.json"
            execute_job(i, jsonIn, jsonOut, geometry)

    if sys.argv[1] == "diff":
        treshold = float(sys.argv[2])
        geometry = sys.argv[3]

        i = 0
        while True:
            jsonIn = f"calib_{i}.json"
            jsonOut = f"calib_{i+1}.json"
            execute_job(i, jsonIn, jsonOut, geometry, treshold)

            if os.path.isfile("finish"):
                break

            i += 1