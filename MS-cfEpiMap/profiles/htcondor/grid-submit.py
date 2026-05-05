#!/usr/bin/env python3

import sys
import htcondor
from os import makedirs
from os.path import join
from uuid import uuid4

from snakemake.utils import read_job_properties


jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

UUID = uuid4()  # random UUID
jobDir = "/home/mattia/.condor_jobs{}_{}".format(job_properties["jobid"], UUID)
makedirs(jobDir, exist_ok=True)

sub = htcondor.Submit(
    {
        "executable": "/bin/bash",
        "arguments": jobscript,
        "max_retries": "1",
        "log": join(jobDir, "condor.log"),
        "output": join(jobDir, "condor.out"),
        "error": join(jobDir, "condor.err"),
        "getenv": "True",
        "request_cpus": str(job_properties["threads"]),
    }
)

# If Snakemake requests GPUs, tell Condor and restrict to GPU nodes
request_gpus = job_properties["resources"].get("gpu", None)
if request_gpus is not None:
    sub["request_gpus"] = str(request_gpus)
    # Use the full hostnames as Condor reports Machine as FQDN (e.g. monod30.mbb.ki.se)
    sub["requirements"] = '(Machine == "monod30.mbb.ki.se" || Machine == "monod31.mbb.ki.se" || Machine == "monod32.mbb.ki.se" || Machine == "monod33.mbb.ki.se")'

request_memory = job_properties["resources"].get("mem_mb", None)
# If this job is using GPU slots, avoid pinning it to a small Condor slot memory value
# (GPU slots often advertise a low Memory value per GPU slot, which can prevent matching).
if request_memory is not None and request_gpus is None:
    sub["request_memory"] = str(request_memory)

request_disk = job_properties["resources"].get("disk_mb", None)
if request_disk is not None:
    sub["request_disk"] = str(request_disk)

# Support user-requested walltime (seconds) or human-readable time in params
walltime = job_properties["resources"].get("walltime", None)
if walltime is None:
    # Fallback: allow a params.time string like "20:00:00"
    time_param = job_properties.get("params", {}).get("time", None)
    if time_param:
        try:
            # parse HH:MM:SS
            h, m, s = map(int, time_param.split(':'))
            walltime = h * 3600 + m * 60 + s
        except Exception:
            try:
                walltime = int(time_param)
            except Exception:
                walltime = None

if walltime is not None:
    # Add as a custom classad attribute that HTCondor/cluster policies can use
    sub["+MaxWallTime"] = str(int(walltime))

# Add kerberos credentials
# c.f. https://batchdocs.web.cern.ch/local/pythonapi.html
#col = htcondor.Collector()
#credd = htcondor.Credd()
#credd.add_user_cred(htcondor.CredTypes.Kerberos, None)
#sub["MY.SendCredential"] = "True"

schedd = htcondor.Schedd()
clusterID = schedd.submit(sub)

# print jobid for use in Snakemake
print("{}_{}_{}".format(job_properties["jobid"], UUID, clusterID))
