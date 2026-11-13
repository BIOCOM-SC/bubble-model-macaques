#!/bin/bash
# Slurm job written by runs/launch.py for '{name}' -- do not edit (site settings: hpc_site.toml)
#SBATCH --job-name={job_name}
#SBATCH --account={account}
#SBATCH --qos={qos}
#SBATCH --time={time}
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ntasks_per_node}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --output={log_dir}/slurm-%j.out
#SBATCH --error={log_dir}/slurm-%j.out
{extra_sbatch}
module purge
{modules}
set -eo pipefail          # after `module` (Lmod does not like -u)
export SRUN_CPUS_PER_TASK=${{SLURM_CPUS_PER_TASK}}
export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK}}
{env}
cd {root}
{python} runs/launch.py run {run_args}
