cd restart{restart_number}
JID_JOB{jid_job}=`sbatch  --dependency=afterok:$JID_JOB{jid_job_previous} job.slurm | cut -d " " -f 4`
cd ..