#!/bin/bash
#SBATCH --job-name=all_samples
#SBATCH --account=boucher
#SBATCH --qos=boucher
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nbonin@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --time=1:00:00
#SBATCH --output=logs/%j_disp.log
#SBATCH --error=logs/%j_disp.log

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date


##----------------------------------------------------------
# Modules
module load snakemake


##----------------------------------------------------------
# Run

#snakemake --touch --forceall -j 10000 --rerun-incomplete 

snakemake --cluster "sbatch -A {cluster.account} -q {cluster.qos} -c {cluster.cpus-per-task} -N {cluster.Nodes} \
  -t {cluster.runtime} --mem {cluster.mem} -J {cluster.jobname} --mail-type={cluster.mail_type} \
  --mail-user={cluster.mail} --output {cluster.out} --error {cluster.err}" \
  --cluster-config config/cluster.json --jobs 500 --rerun-trigger mtime --latency-wait 20  --use-envmodules --use-conda --conda-frontend mamba
