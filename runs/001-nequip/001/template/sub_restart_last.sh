cd restart{restart_number}
sbatch  --dependency=afterok:$JID_JOB{jid_job_previous}  job.slurm
cd ..