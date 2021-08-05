import sys
import subprocess
import json

def execute_job(step, jsonIn, jsonOut, geometry):
    print(f"execute {step} step")

    command = "cmsRun Analyzer/DiamondTimingAnalyzer/python/test.py"
    command += " validOOT=\"-1\""
    command += f" calibInput=\"{jsonIn}\""
    command += f" calibOutput=\"{jsonOut}\""
    command += f" geometryFile=\"{geometry}\""
    command += f" loopIndex=\"{step}\""
    
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
        treshold = int(sys.argv[2])
        geometry = sys.argv[3]

        i = 0
        while True:
            jsonIn = f"calib_{i}.json"
            jsonOut = f"calib_{i+1}.json"
            execute_job(i, jsonIn, jsonOut, geometry)

            #check JSON
            #new_res = json.loads(open(jsonOut, "r").read().replace("\n", " "))
            #prev_res = json.loads(open(jsonIn, "r").read().replace("\n", " "))

            n += 1