#!/usr/bin/env python

import sys
import htcondor
from htcondor import JobEventType
from os.path import join


def print_and_exit(s):
    print(s)
    exit()


jobID, UUID, clusterID = sys.argv[1].split("_")

jobDir = "/home/mattia/.condor_jobs{}_{}".format(jobID, UUID)
jobLog = join(jobDir, "condor.log")

failed_states = [
    JobEventType.JOB_HELD,
    JobEventType.JOB_ABORTED,
    JobEventType.EXECUTABLE_ERROR,
]

import os
import time

# Log file may not exist yet if condor hasn't flushed it — treat as still running
if not os.path.exists(jobLog):
    print_and_exit("running")

try:
    jel = htcondor.JobEventLog(join(jobLog))
    for event in jel.events(stop_after=1):
        if event.type in failed_states:
            print_and_exit("failed")
        if event.type is JobEventType.JOB_TERMINATED:
            if event["ReturnValue"] == 0:
                print_and_exit("success")
            print_and_exit("failed")
except OSError as e:
    # Log exists but isn't readable yet — still initialising
    print_and_exit("running")

print_and_exit("running")
