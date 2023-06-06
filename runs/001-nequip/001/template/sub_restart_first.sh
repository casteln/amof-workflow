#!/bin/sh
cd restart{restart_number}
JID_JOB{jid_job}=`sbatch  job.slurm | cut -d " " -f 4`
cd ..