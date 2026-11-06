#!/usr/bin/env python3
"""Benchmark one exported stage MLP on CPU (1/N threads) and optional CUDA."""
from __future__ import annotations
import argparse, csv, json, statistics, time
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn

class StageMLP(nn.Module):
    def __init__(self):
        super().__init__(); self.layers=nn.Sequential(nn.Linear(2,64),nn.ReLU(),nn.Linear(64,64),nn.ReLU(),nn.Linear(64,32),nn.ReLU(),nn.Linear(32,3))
    def forward(self,x): return self.layers(x)

def bench(fn, sync, repeats, warmups):
    for _ in range(warmups): fn(); sync()
    times=[]
    for _ in range(repeats):
        sync(); t=time.perf_counter(); fn(); sync(); times.append(time.perf_counter()-t)
    return statistics.median(times)

def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("model_dir",type=Path); ap.add_argument("--evals",type=int,default=32768); ap.add_argument("--cpu-threads",type=int,default=8); ap.add_argument("--repeats",type=int,default=11); ap.add_argument("--warmups",type=int,default=3); ap.add_argument("--output",type=Path)
    a=ap.parse_args(); meta=json.loads((a.model_dir/"stage_surrogate.json").read_text()); base=np.load(a.model_dir/"calibration_inputs.npy")
    if not len(base): raise ValueError("Empty calibration inputs")
    reps=(a.evals+len(base)-1)//len(base); x=np.tile(base,(reps,1))[:a.evals].astype(np.float32)
    rows=[]; model=StageMLP(); model.load_state_dict(torch.load(a.model_dir/"stage_surrogate.pt",map_location="cpu")); model.eval()
    mac=int(meta["mac_per_evaluation"]); flops=int(meta["nominal_flop_per_evaluation"])
    def record(label,seconds):
        rate=a.evals/seconds; rows.append({"backend":label,"evaluations":a.evals,"seconds_median":seconds,"function_evaluations_per_second":rate,"nominal_MFLOP_per_second":rate*flops/1e6,"MAC_per_evaluation":mac})
    for threads,label in ((1,"CPU single-thread"),(a.cpu_threads,f"CPU {a.cpu_threads}-thread")):
        torch.set_num_threads(threads); xt=torch.from_numpy(x)
        with torch.inference_mode():
            sec=bench(lambda:model(xt),lambda:None,a.repeats,a.warmups)
        record(label,sec)
    if torch.cuda.is_available():
        device=torch.device("cuda"); gpu_model=model.to(device); xt=torch.from_numpy(x).to(device)
        torch.cuda.synchronize()
        with torch.inference_mode():
            sec=bench(lambda:gpu_model(xt),torch.cuda.synchronize,a.repeats,a.warmups)
        record(f"GPU resident: {torch.cuda.get_device_name(0)}",sec)
        host=torch.from_numpy(x)
        def h2d():
            inp=host.to(device,non_blocking=False); out=gpu_model(inp); return out.cpu()
        with torch.inference_mode():
            sec=bench(h2d,torch.cuda.synchronize,a.repeats,a.warmups)
        record("GPU host-to-device + inference",sec)
    else:
        print("GPU: unavailable in this Python environment; rerun under Slurm with CUDA-enabled PyTorch.")
    print(f"Model: {meta['model']} | {a.evals} resident evaluations/pass | {mac} MAC/eval, {flops} nominal FLOP/eval")
    for r in rows: print(f"{r['backend']}: {r['function_evaluations_per_second']/1e6:.4f} M eval/s, {r['nominal_MFLOP_per_second']/1e3:.3f} nominal GFLOP/s")
    result={"model":meta["model"],"measurement_scope":"resident batched network forward; CPU/GPU host overhead reported separately","rows":rows}
    if a.output: a.output.write_text(json.dumps(result,indent=2)+"\n"); print(f"Report: {a.output}")
    return 0
if __name__=="__main__": raise SystemExit(main())
