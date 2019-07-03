## Reference
- BT Thomas Yeo, et al. "Functional specialization and flexibility in human association cortex." Cerebral Cortex, 2015.

## Instructions

To replicate the 12 cognitive components estimated from the BrainMap database (Yeo et al., Cerebral Cortex, 2015), please follow the steps below.
Note that the instructions are meant for CBIG's collaborators who have access to CBIG's HPC server and are authorized to access the BrainMap database.

1. From `bash` terminal, execute `bash CBIG_AuthorTopicEM_submit_replication_jobs.sh`, which submits 100 pbs jobs for 100 different initializations for estimating the 12 components to CBIG's queue system.
2. Once all the jobs are finished, run `CBIG_AuthorTopicEM_ReplicateAverageOrBestSolutionWrapper` in Matlab to compute the most typical 12-component solution among all initializations.

Check `/mnt/eql/yeo1/CBIG_private_data/replication/stable_projects/meta-analysis/Yeo2015_AuthorTopicEM/results/matlab2018b_results` for reference results.
