#!/usr/bin/env python3
"""Qualify the complete-stage ReLU map by holding out complete speed lines."""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from pathlib import Path
from typing import Any
import numpy as np
import torch
import torch.nn as nn

FEATURES=("rotational_speed_rpm","mass_flow")
TARGETS=("outlet_total_temperature","outlet_total_pressure","specific_work")
FLAGS=("converged","validated","basic_physics_valid")

class StageMLP(nn.Module):
    """Same 2 -> 64 -> 64 -> 32 -> 3 network as the exported candidate."""
    def __init__(self):
        super().__init__()
        self.layers=nn.Sequential(nn.Linear(2,64),nn.ReLU(),nn.Linear(64,64),
            nn.ReLU(),nn.Linear(64,32),nn.ReLU(),nn.Linear(32,3))
    def forward(self,x): return self.layers(x)

def number(row,name,line):
    try: value=float(row[name])
    except (KeyError,TypeError,ValueError) as exc:
        raise ValueError(f"CSV line {line}: invalid {name}") from exc
    if not math.isfinite(value): raise ValueError(f"CSV line {line}: non-finite {name}")
    return value

def read_samples(path,model_name):
    required={"model",*FEATURES,*TARGETS,*FLAGS}
    unique={}; selected=0
    with path.open(newline="",encoding="utf-8") as stream:
        reader=csv.DictReader(stream)
        missing=required.difference(reader.fieldnames or ())
        if missing: raise ValueError(f"{path} missing columns: {sorted(missing)}")
        for line,row in enumerate(reader,start=2):
            if row.get("model","").strip()!=model_name: continue
            selected+=1
            if any(row.get(k,"").strip().lower() not in {"1","true","yes"} for k in FLAGS):
                raise ValueError(f"{path}:{line}: sample failed solver validation")
            if row.get("stage_status","Complete").strip()!="Complete":
                raise ValueError(f"{path}:{line}: sample is not a complete stage")
            speed,flow=(number(row,k,line) for k in FEATURES)
            target=np.array([number(row,k,line) for k in TARGETS],dtype=np.float64)
            if speed<=0 or flow<=0 or np.any(target<=0):
                raise ValueError(f"{path}:{line}: inputs and targets must be positive")
            key=(speed,flow)
            if key in unique:
                if not np.allclose(unique[key],target,rtol=1e-8,atol=1e-8):
                    raise ValueError(f"{path}:{line}: duplicate condition {key} disagrees")
            else: unique[key]=target
    if not unique: raise ValueError(f"No validated {model_name} samples in {path}")
    points=sorted(unique.items())
    x=np.asarray([k for k,_ in points],dtype=np.float64)
    y=np.asarray([v for _,v in points],dtype=np.float64)
    if len(np.unique(x[:,0]))<3: raise ValueError("Need at least three speed lines")
    return x,y,selected

def encode(x,low,high):
    span=high-low
    if np.any(span<=0): raise ValueError("Training fold must span both input features")
    # No clipping: edge folds must expose out-of-range speed coordinates.
    return (2.0*(x-low)/span-1.0).astype(np.float32)

def fit_predict(train_x,train_y,eval_x,seed,epochs,learning_rate,cpu_threads):
    torch.manual_seed(seed); torch.set_num_threads(cpu_threads)
    low,high=train_x.min(axis=0),train_x.max(axis=0)
    tx=torch.from_numpy(encode(train_x,low,high))
    ex=torch.from_numpy(encode(eval_x,low,high))
    log_y=np.log(train_y); mean=log_y.mean(axis=0); std=log_y.std(axis=0)
    std=np.where(std>1e-8,std,1.0)
    ty=torch.from_numpy(((log_y-mean)/std).astype(np.float32))
    model=StageMLP(); optimizer=torch.optim.Adam(model.parameters(),lr=learning_rate)
    loss=nn.MSELoss(); model.train()
    for _ in range(epochs):
        optimizer.zero_grad(set_to_none=True); loss(model(tx),ty).backward(); optimizer.step()
    model.eval()
    with torch.inference_mode(): pred=model(ex).numpy().astype(np.float64)
    return np.exp(pred*std+mean)

def summarize(actual,predicted):
    rel=(predicted-actual)/actual
    return {"sample_count":int(actual.size),
        "max_absolute_relative_error":float(np.max(np.abs(rel))),
        "p99_absolute_relative_error":float(np.percentile(np.abs(rel),99)),
        "mean_absolute_relative_error":float(np.mean(np.abs(rel))),
        "rms_relative_error":float(np.sqrt(np.mean(rel**2))),
        "rmse":float(np.sqrt(np.mean((predicted-actual)**2)))}

def main():
    default=Path("build/relu-surrogate/complete-stage-boundary-grid/"
        "benchmarks/results-complete-stage-characteristic/"
        "complete-stage-boundary-grid.csv")
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("csv",nargs="?",type=Path,default=default)
    p.add_argument("--model",default="PengRobinson",choices=("PengRobinson","IdealGas"))
    p.add_argument("--output-dir",type=Path,
        default=Path("build/relu-surrogate/stage-speedline-qualification"))
    p.add_argument("--epochs",type=int,default=600)
    p.add_argument("--learning-rate",type=float,default=.01)
    p.add_argument("--seed",type=int,default=20261106)
    p.add_argument("--cpu-threads",type=int,default=1)
    a=p.parse_args()
    if a.epochs<1 or a.learning_rate<=0 or a.cpu_threads<1:
        p.error("epochs, learning-rate, and cpu-threads must be positive")
    x,y,selected=read_samples(a.csv,a.model); speeds=np.unique(x[:,0])
    if len(speeds)<3: raise ValueError("Need at least three speed lines")
    records=[]; per_speed={}
    for fold,speed in enumerate(speeds):
        held=x[:,0]==speed; train_x,eval_x=x[~held],x[held]
        train_y,eval_y=y[~held],y[held]
        pred=fit_predict(train_x,train_y,eval_x,a.seed+fold,a.epochs,
                         a.learning_rate,a.cpu_threads)
        kind=("edge_extrapolation" if fold in (0,len(speeds)-1)
              else "interior_interpolation")
        per_speed[str(float(speed))]={t:summarize(eval_y[:,j],pred[:,j])
            for j,t in enumerate(TARGETS)}
        for i in range(len(eval_x)):
            row={"held_out_speed_rpm":float(speed),"fold_type":kind,
                 "mass_flow":float(eval_x[i,1]),"training_rows":int(len(train_x))}
            for j,t in enumerate(TARGETS):
                row[f"{t}_actual"]=float(eval_y[i,j])
                row[f"{t}_predicted"]=float(pred[i,j])
                row[f"{t}_relative_error"]=float((pred[i,j]-eval_y[i,j])/eval_y[i,j])
            records.append(row)
    masks={"all_held_out_lines":np.ones(len(records),dtype=bool),
        "interior_interpolation_lines":np.array([r["fold_type"]=="interior_interpolation" for r in records]),
        "edge_extrapolation_lines":np.array([r["fold_type"]=="edge_extrapolation" for r in records])}
    aggregate={}
    for group,mask in masks.items():
        aggregate[group]={}
        for t in TARGETS:
            actual=np.array([r[f"{t}_actual"] for r in records])[mask]
            pred=np.array([r[f"{t}_predicted"] for r in records])[mask]
            aggregate[group][t]=summarize(actual,pred)
    a.output_dir.mkdir(parents=True,exist_ok=True)
    pred_path=a.output_dir/"held_out_speedline_predictions.csv"
    fields=["held_out_speed_rpm","fold_type","mass_flow","training_rows"]
    for t in TARGETS: fields += [f"{t}_actual",f"{t}_predicted",f"{t}_relative_error"]
    with pred_path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader(); writer.writerows(records)
    report={"qualification":"leave-one-complete-speed-line-out","model":a.model,
        "architecture":"2->64->64->32->3 Linear/ReLU hidden","features":list(FEATURES),
        "targets":list(TARGETS),"dataset":{"path":str(a.csv),
        "sha256":hashlib.sha256(a.csv.read_bytes()).hexdigest(),"selected_csv_rows":selected,
        "unique_conditions":int(len(x)),"speed_lines_rpm":[float(n) for n in speeds],
        "mass_flow_range_kg_s":[float(x[:,1].min()),float(x[:,1].max())]},
        "fit":{"epochs":a.epochs,"learning_rate":a.learning_rate,"base_seed":a.seed,
        "cpu_threads":a.cpu_threads,
        "feature_scaling":"min/max fitted on each fold's training rows only",
        "target_scaling":"log mean/std fitted on each fold's training rows only"},
        "split":{"fold_count":int(len(speeds)),"held_out_rows":len(records),
        "endpoint_folds":"speed extrapolation checks","interior_folds":"speed interpolation checks",
        "duplicates":"deduplicated after target consistency check"},
        "metrics_by_held_out_speed_rpm":per_speed,"aggregate_metrics":aggregate,
        "interpretation":"Approximates the current numerical model over this sampled fixed-setup domain. This is not experimental validation and does not qualify other inlet states, geometries, gas mixtures, or transient conditions.",
        "predictions_csv":str(pred_path)}
    report_path=a.output_dir/"qualification_report.json"
    report_path.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(f"Model: {a.model}; unique conditions: {len(x)}")
    print(f"Folds: {len(speeds)} ({len(speeds)-2} interior, 2 edge)")
    for group,by_target in aggregate.items():
        print(f"{group}:")
        for target,m in by_target.items():
            print(f"  {target}: max|rel|={m['max_absolute_relative_error']:.6g} "
                  f"p99|rel|={m['p99_absolute_relative_error']:.6g} "
                  f"rms_rel={m['rms_relative_error']:.6g}")
    print(f"Predictions: {pred_path}"); print(f"Report: {report_path}")
    return 0

if __name__=="__main__": raise SystemExit(main())
