# Complete-stage surrogate hardware comparison

This experiment uses the PengRobinson complete-stage map, the mass-flow coordinates selected by the paired coordinate ablation, and one shared `2 -> 64 -> 64 -> 32 -> 3` ReLU MLP. The Hailo-10H is an inference accelerator (NPU); it uses the same ONNX model compiled to HEF. The exporter emits ONNX, PyTorch weights, and NumPy reference weights. Hailo host-side decoding applies `exp(z * target_log_std + target_log_mean)` after inference.

The model's training-fit metrics are not independent validation. The earlier coordinate-ablation report evaluates a different training experiment and is not the qualification report for this hardware candidate. Use the independent held-out-speed-line qualification below for this exported architecture. The Hailo validator measures quantization deviation on supplied direct-stage rows; that measures quantization error, not speed-line generalization.

Scope: this measures the three stage responses on the supplied valid sample domain. It does not yet interpolate the boundary curves at arbitrary speed or reject out-of-envelope runtime queries; keep those checks in the surrounding operating-envelope logic.

## Independent held-out-speed-line qualification

Training metrics and the Hailo comparison above do not establish speed-line
generalization. This independent leave-one-complete-speed-line-out workflow
fits a fresh copy of the same `2 -> 64 -> 64 -> 32 -> 3` network per fold.
It excludes the held-out line from training and from input and target scaling.
The endpoint folds measure speed extrapolation; interior folds measure speed
interpolation. Per-row predictions and per-speed/aggregate metrics are saved.

Run on the Dell 7710 with the CPU PyTorch environment:

```bash
cd ~/Sync/TurboCompGPU
. .venv/bin/activate
python benchmarks/surrogate/qualify_stage_surrogate.py \
  build/relu-surrogate/complete-stage-boundary-grid/benchmarks/results-complete-stage-characteristic/complete-stage-boundary-grid.csv \
  --model PengRobinson \
  --output-dir build/relu-surrogate/stage-speedline-qualification \
  --cpu-threads 2
```

For Slurm on the 7710:

```bash
cd ~/Sync/TurboCompGPU
PROJECT="$HOME/Sync/TurboCompGPU"
mkdir -p "$PROJECT/build/relu-surrogate/slurm-logs"
sbatch --nodelist=dell7710 --partition=compute --chdir="$PROJECT" \
  --cpus-per-task=2 --mem=4G --time=00:30:00 \
  --output="$PROJECT/build/relu-surrogate/slurm-logs/qualify-%j.out" \
  --error="$PROJECT/build/relu-surrogate/slurm-logs/qualify-%j.err" \
  --wrap='. .venv/bin/activate && python benchmarks/surrogate/qualify_stage_surrogate.py build/relu-surrogate/complete-stage-boundary-grid/benchmarks/results-complete-stage-characteristic/complete-stage-boundary-grid.csv --model PengRobinson --output-dir build/relu-surrogate/stage-speedline-qualification --cpu-threads "$SLURM_CPUS_PER_TASK"'
```

Outputs are `qualification_report.json` and
`held_out_speedline_predictions.csv`. Review interior and edge
metrics separately before using the candidate in the pipeline component. This
only establishes approximation of this numerical model over the sampled
fixed-setup domain; it does not establish experimental accuracy or validity
for other compressor inputs.

## Train and export on the 7710

Use the already-generated C++ complete-stage CSV and submit training/export through Slurm:

```bash
cd ~/Sync/TurboCompGPU
PROJECT="$HOME/Sync/TurboCompGPU"
DATA=build/relu-surrogate/complete-stage-boundary-grid/benchmarks/results-complete-stage-characteristic/complete-stage-boundary-grid.csv
OUT=build/relu-surrogate/hardware/pengrobinson
mkdir -p "$PROJECT/build/relu-surrogate/slurm-logs"

sbatch --nodelist=dell7710 --partition=compute --chdir="$PROJECT" \
  --cpus-per-task=2 --mem=4G --time=00:20:00 \
  --output="$PROJECT/build/relu-surrogate/slurm-logs/export-%j.out" \
  --error="$PROJECT/build/relu-surrogate/slurm-logs/export-%j.err" \
  --wrap='. .venv/bin/activate && python benchmarks/surrogate/train_export_stage_surrogate.py '"$DATA"' --model PengRobinson --output-dir '"$OUT"' --cpu-threads "$SLURM_CPUS_PER_TASK"'
```

## CPU and GPU throughput on the 7710

```bash
cd ~/Sync/TurboCompGPU
PROJECT="$HOME/Sync/TurboCompGPU"
MODEL="$PROJECT/build/relu-surrogate/hardware/pengrobinson"

# CPU single-thread and CPU parallel from the repo's CPU PyTorch environment
sbatch --nodelist=dell7710 --partition=compute --chdir="$PROJECT" \
  --cpus-per-task=8 --mem=3G --time=00:10:00 \
  --output="$PROJECT/build/relu-surrogate/slurm-logs/cpu-bench-%j.out" \
  --error="$PROJECT/build/relu-surrogate/slurm-logs/cpu-bench-%j.err" \
  --wrap='. .venv/bin/activate && python benchmarks/surrogate/benchmark_stage_surrogate.py '"$MODEL"' --cpu-threads "$SLURM_CPUS_PER_TASK" --output '"$MODEL"'/cpu-gpu-benchmark.json'

# GPU measurement on the 7710; use the CUDA-enabled environment already used for GPU benchmarks
srun --nodelist=dell7710 --chdir="$PROJECT" --gres=gpu:quadro_m3000m:1 \
  bash -lc 'source "$HOME/Sync/FrictionLibraryCpp20/.venv-hailo-mlp/bin/activate" && \
    python benchmarks/surrogate/benchmark_stage_surrogate.py \
      build/relu-surrogate/hardware/pengrobinson --cpu-threads 8 \
      --output build/relu-surrogate/hardware/pengrobinson/gpu-benchmark.json'
```

The benchmark reports resident batched function evaluations/s, CPU single-thread, CPU multi-thread, GPU resident, and GPU host-to-device plus inference when CUDA is available. The nominal arithmetic estimate is 6,368 MAC/evaluation (12,736 FLOP/evaluation, excluding bias, ReLU and decoding).

## Stage the files on shared storage

The 7680's Slurm jobs read `/cluster/projects`; copy only the source files and generated inputs needed for compile/runtime from the 7710:

```bash
cd ~/Sync/TurboCompGPU
SHARED=/cluster/projects/TurboCompGPU
mkdir -p "$SHARED/benchmarks/surrogate" "$SHARED/scripts/hailo"
mkdir -p "$SHARED/build/relu-surrogate/hardware"
cp benchmarks/surrogate/*.py "$SHARED/benchmarks/surrogate/"
cp scripts/hailo/compile_stage_surrogate.py scripts/hailo/validate_stage_surrogate.py \
  scripts/hailo/stage_surrogate_gpu.alls "$SHARED/scripts/hailo/"
cp -a build/relu-surrogate/complete-stage-boundary-grid \
  "$SHARED/build/relu-surrogate/"
cp -a build/relu-surrogate/hardware/pengrobinson \
  "$SHARED/build/relu-surrogate/hardware/"
```

## Compile for Hailo-10H

Submit DFC on dell7680 with its RTX 4090 allocation. The container settings and XLA CUDA-library links below match the working DFC 5.4 setup on this cluster:

```bash
SHARED=/cluster/projects/TurboCompGPU
srun --nodelist=dell7680 --chdir=/tmp --gres=gpu:rtx_4090:1 bash -lc '
  set -euo pipefail
  PROJECT=/cluster/projects/TurboCompGPU
  podman --cgroup-manager=cgroupfs run --rm \
    --cgroupns=host --device nvidia.com/gpu=0 \
    --security-opt=label=disable \
    -v "$PROJECT:/workspace" -w /workspace \
    localhost/hailo-dfc:5.4 bash -lc "
      set -e
      source /opt/hailo-dfc/bin/activate
      CUDA_ROOT=/tmp/xla-cuda
      NVIDIA_ROOT=/opt/hailo-dfc/lib/python3.12/site-packages/triton/backends/nvidia
      mkdir -p \$CUDA_ROOT/nvvm/libdevice \$CUDA_ROOT/bin
      ln -sf \$NVIDIA_ROOT/lib/libdevice.10.bc \$CUDA_ROOT/nvvm/libdevice/libdevice.10.bc
      PTXAS=\$(find /opt/hailo-dfc -type f -name ptxas -print -quit)
      if [ -n \"\$PTXAS\" ]; then ln -sf \$PTXAS \$CUDA_ROOT/bin/ptxas; fi
      export XLA_FLAGS=--xla_gpu_cuda_data_dir=\$CUDA_ROOT
      python scripts/hailo/compile_stage_surrogate.py \\
        build/relu-surrogate/hardware/pengrobinson/stage_surrogate.onnx \\
        build/relu-surrogate/hardware/pengrobinson/calibration_inputs.npy \\
        --output-dir build/relu-surrogate/hardware/pengrobinson/hailo \\
        --hw-arch hailo10h \\
        --model-script scripts/hailo/stage_surrogate_gpu.alls
    "
'
```

The compiler writes parsed/optimized HAR files and `stage_surrogate.hef`.

## Hailo accuracy and throughput on the 7680

With the HEF and CSV on shared storage, run both under Slurm:

```bash
PROJECT=/cluster/projects/TurboCompGPU
srun --nodelist=dell7680 --chdir="$PROJECT" --gres=hailo10h:1 \
  bash -lc 'source /cluster/software/hailo10h/pyhailort-5.4.0-py313/bin/activate && \
    python scripts/hailo/validate_stage_surrogate.py \
      build/relu-surrogate/hardware/pengrobinson/hailo/stage_surrogate.hef \
      build/relu-surrogate/hardware/pengrobinson \
      build/relu-surrogate/complete-stage-boundary-grid/benchmarks/results-complete-stage-characteristic/complete-stage-boundary-grid.csv \
      --model PengRobinson --limit 0 \
      --output build/relu-surrogate/hardware/pengrobinson/hailo/validation.json'

srun --nodelist=dell7680 --chdir="$PROJECT" --gres=hailo10h:1 \
  bash -lc 'hailortcli benchmark \
    build/relu-surrogate/hardware/pengrobinson/hailo/stage_surrogate.hef'
```

`hailortcli benchmark` reports the Hailo network FPS, which is one MLP function evaluation per frame. `validate_stage_surrogate.py` also reports serial Python/HailoRT rate and the HEF quantization error relative to the float network. These timing paths differ: resident CPU/GPU rates exclude Python and PCIe call overhead, while serial HailoRT includes it. Report both categories rather than treating them as identical end-to-end latencies.
