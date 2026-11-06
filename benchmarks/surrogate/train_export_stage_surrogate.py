#!/usr/bin/env python3
"""Train/export the PengRobinson complete-stage surrogate for CPU/GPU/Hailo."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn

TARGETS = ("outlet_total_temperature", "outlet_total_pressure", "specific_work")

class StageMLP(nn.Module):
    def __init__(self):
        super().__init__()
        self.layers = nn.Sequential(
            nn.Linear(2, 64), nn.ReLU(),
            nn.Linear(64, 64), nn.ReLU(),
            nn.Linear(64, 32), nn.ReLU(),
            nn.Linear(32, 3),
        )
    def forward(self, x):
        return self.layers(x)

def val(row, key):
    try: x = float(row[key])
    except (KeyError, ValueError, TypeError) as exc: raise ValueError(f"Invalid {key}: {row.get(key)!r}") from exc
    if not math.isfinite(x): raise ValueError(f"Non-finite {key}")
    return x

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", type=Path)
    ap.add_argument("--model", default="PengRobinson", choices=("PengRobinson", "IdealGas"))
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--epochs", type=int, default=600)
    ap.add_argument("--learning-rate", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=20261106)
    ap.add_argument("--cpu-threads", type=int, default=1)
    args = ap.parse_args()
    if args.epochs < 1 or args.learning_rate <= 0 or args.cpu_threads < 1: ap.error("training settings must be positive")

    with args.csv.open(newline="") as f: rows = [r for r in csv.DictReader(f) if r.get("model") == args.model]
    if not rows: raise ValueError(f"No {args.model} samples in {args.csv}")
    req = {"rotational_speed_rpm", "mass_flow", "converged", "validated", "basic_physics_valid", *TARGETS}
    missing = req.difference(rows[0])
    if missing: raise ValueError(f"Missing columns: {sorted(missing)}")
    for r in rows:
        if any(r[k].strip().lower() not in {"1", "true", "yes"} for k in ("converged", "validated", "basic_physics_valid")):
            raise ValueError("Training CSV contains an unvalidated stage sample")
        if r.get("stage_status", "Complete") != "Complete": raise ValueError("Training CSV contains an incomplete stage")

    rpm = np.array([val(r, "rotational_speed_rpm") for r in rows], dtype=np.float64)
    flow = np.array([val(r, "mass_flow") for r in rows], dtype=np.float64)
    y = np.array([[val(r, t) for t in TARGETS] for r in rows], dtype=np.float64)
    if np.any(y <= 0): raise ValueError("Log-target training requires positive outputs")
    rmin, rmax, mmin, mmax = float(rpm.min()), float(rpm.max()), float(flow.min()), float(flow.max())
    if rmax <= rmin or mmax <= mmin: raise ValueError("Training data must span speed and mass flow")
    x = np.column_stack((2*(rpm-rmin)/(rmax-rmin)-1, 2*(flow-mmin)/(mmax-mmin)-1)).astype(np.float32)

    np.random.seed(args.seed); torch.manual_seed(args.seed); torch.set_num_threads(args.cpu_threads)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log_y = np.log(y).astype(np.float32)
    ymean = log_y.mean(axis=0); ystd = log_y.std(axis=0); ystd[ystd < 1e-8] = 1.0
    yn = (log_y-ymean)/ystd
    xt = torch.as_tensor(x, device=device); yt = torch.as_tensor(yn, device=device)
    model = StageMLP().to(device); opt = torch.optim.Adam(model.parameters(), lr=args.learning_rate); mse = nn.MSELoss()
    for epoch in range(args.epochs):
        opt.zero_grad(set_to_none=True); loss = mse(model(xt), yt); loss.backward(); opt.step()
    model.eval()
    with torch.no_grad():
        pred_log = model(xt).cpu().numpy().astype(np.float64)*ystd+ymean
    pred = np.exp(pred_log)
    rel = (pred-y)/y
    metrics = {}
    for j,t in enumerate(TARGETS):
        a = np.abs(rel[:,j]); metrics[t] = {
            "training_max_relative_error": float(a.max()),
            "training_p99_relative_error": float(np.percentile(a,99)),
            "training_rms_relative_error": float(np.sqrt(np.mean(rel[:,j]**2))),
        }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    state = model.cpu().state_dict()
    torch.save(state, args.output_dir/"stage_surrogate.pt")
    np.savez(args.output_dir/"stage_surrogate_weights.npz", **{
        key.replace("layers.", "layer_").replace(".weight", "_weight").replace(".bias", "_bias"):
            value.detach().cpu().numpy()
        for key, value in state.items()
    })
    metadata = {
        "model": args.model, "architecture": "2->64->64->32->3 Linear/ReLU hidden",
        "input_features": ["normalized rotational speed", "globally normalized mass flow"],
        "speed_min_rpm": rmin, "speed_max_rpm": rmax,
        "mass_flow_min_kg_s": mmin, "mass_flow_max_kg_s": mmax,
        "target_names": list(TARGETS), "target_log_mean": ymean.tolist(), "target_log_std": ystd.tolist(),
        "target_transform": "network predicts standardized ln(target); host decodes exp(z*std+mean)",
        "epochs": args.epochs, "learning_rate": args.learning_rate, "seed": args.seed,
        "training_samples": len(rows), "training_device": str(device),
        "training_fit_metrics_not_independent_validation": metrics,
        "mac_per_evaluation": 6368, "nominal_flop_per_evaluation": 12736,
        "note": "Use leave-one-speed-line-out report for generalization metrics. This artifact is the full-data model for hardware timing.",
    }
    (args.output_dir/"stage_surrogate.json").write_text(json.dumps(metadata,indent=2)+"\n")
    np.save(args.output_dir/"calibration_inputs.npy", x)
    dummy = torch.zeros((1,2),dtype=torch.float32)
    torch.onnx.export(model, dummy, args.output_dir/"stage_surrogate.onnx", export_params=True,
        opset_version=11, do_constant_folding=True, dynamo=False, input_names=["input_coords"], output_names=["scaled_log_targets"])
    print(f"Model: {args.model}; device={device}; samples={len(rows)}; MAC/eval=6368")
    for t,m in metrics.items(): print(f"{t}: train max_rel={m['training_max_relative_error']:.6g} rms_rel={m['training_rms_relative_error']:.6g}")
    print(f"PT: {args.output_dir/'stage_surrogate.pt'}")
    print(f"NumPy reference weights: {args.output_dir/'stage_surrogate_weights.npz'}")
    print(f"ONNX: {args.output_dir/'stage_surrogate.onnx'}")
    print(f"Calibration: {args.output_dir/'calibration_inputs.npy'}")
    print("Training fit metrics are in-sample; use cross-validation for accuracy qualification.")
    return 0

if __name__ == "__main__": raise SystemExit(main())
