#!/usr/bin/env python

import os
import subprocess
import sys

def trackJobs(jobs, waittime=5):
    while len(jobs) != 0:
        for jobid in jobs:
            x = subprocess.Popen(['qstat', '-j', jobid], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            std_out, std_err = x.communicate()
            if std_err :
                jobs.remove(jobid)
                break
        os.system("sleep " + str(waittime))
    return

if __name__ == "__main__":
    trackJobs(sys.argv)
