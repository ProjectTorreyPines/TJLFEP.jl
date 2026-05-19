#!/bin/bash -l
#SBATCH -A m4909
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:05:00
#SBATCH -C cpu
#SBATCH -J TJLFEP

module load julia/1.11.7

threads_main=2

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia"
mkdir -p "${JULIA_DEPOT_PATH}/compiled"

cd ~/.julia/dev/TJLFEP/src/DIIIDfiles

# Start node monitor in the background; poll every 5s until the Julia job ends.
MONITOR_LOG="monitor_${SLURM_JOB_ID}.log"
# Resolve nvidia-smi path in the SLURM environment (modules loaded) so it can
# be passed as a literal into pdsh SSH shells which have a stripped PATH.

(
    while true; do
        echo "=== $(date '+%H:%M:%S') ===" >> "$MONITOR_LOG"
        # GPU query via srun so each node runs it within the job's cgroup
        # --overlap allows this step to run concurrently with the worker srun
        srun --overlap --ntasks-per-node=1 bash -c \
            'nvidia-smi --query-gpu=index,memory.used,utilization.gpu,temperature.gpu \
             --format=csv,noheader,nounits 2>/dev/null \
             | awk -F", " -v h=$(hostname -s) \
               '"'"'{printf "  %s GPU%s: mem=%sMiB util=%s%% temp=%s°C\n",h,$1,$2,$3,$4}'"'"'' \
            2>/dev/null | sort >> "$MONITOR_LOG"
        # CPU/thread info via pdsh (SSH session, no GPU access)
        pdsh -w "$SLURM_JOB_NODELIST" "
            julia_procs=\$(pgrep -x julia 2>/dev/null | wc -l)
            total_threads=0
            proc_info=\"\"
            while IFS= read -r pid; do
                t=\$(awk '/^Threads:/{print \$2}' /proc/\$pid/status 2>/dev/null || echo 0)
                c=\$(ps -p \$pid -o %cpu= 2>/dev/null | tr -d ' ')
                proc_info+=\"  PID \$pid threads=\$t cpu%=\$c\n\"
                total_threads=\$(( total_threads + t ))
            done < <(pgrep -x julia 2>/dev/null)
            cpu=\$(mpstat 1 1 2>/dev/null | awk '/^Average:/{printf \"usr=%.1f%% idle=%.1f%%\",\$3,\$12}')
            echo \"\$(hostname -s): julia_procs=\$julia_procs  total_threads=\$total_threads  \$cpu\"
            echo -e \"\$proc_info\"
        " 2>/dev/null | sort >> "$MONITOR_LOG"
        echo "" >> "$MONITOR_LOG"
        sleep 5
    done
) &
MONITOR_PID=$!

stdbuf -o0 -e0 julia --project --threads=$threads_main --sysimage ../../build/noTJLF_TJLFEP_sysimage.so DIIID_juliaValidation.jl
# stdbuf -o0 -e0 julia --project --threads=20 DIIID_juliaValidation.jl

# Stop the monitor once Julia exits
kill $MONITOR_PID 2>/dev/null
echo "Monitor log written to $MONITOR_LOG"

# Move slurm output and monitor log into the run output directory
OUTDIR=$(find . -maxdepth 1 -type d -name "*_${SLURM_JOB_ID}" 2>/dev/null | head -1)
if [ -n "$OUTDIR" ]; then
    [ -f "$MONITOR_LOG" ] && mv "$MONITOR_LOG" "$OUTDIR/"
    [ -f "${SLURM_SUBMIT_DIR}/slurm-${SLURM_JOB_ID}.out" ] && mv "${SLURM_SUBMIT_DIR}/slurm-${SLURM_JOB_ID}.out" "$OUTDIR/"
fi

# Remove the symlink that submit_sweep left in the DIIIDfiles root (no-op if running directly)
JL_SYMLINK="DIIID_juliaValidation.jl"
[ -L "$JL_SYMLINK" ] && rm "$JL_SYMLINK"
