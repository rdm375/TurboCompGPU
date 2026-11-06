#!/usr/bin/env python3
"""Compare Hailo HEF outputs with the float MLP and direct-stage labels."""
from __future__ import annotations
import argparse,csv,json,math,time
from pathlib import Path
import numpy as np
from hailo_platform import HEF, VDevice, FormatType

TARGETS=("outlet_total_temperature","outlet_total_pressure","specific_work")
def numpy_reference(x, weights_path):
    weights = np.load(weights_path)
    y = np.asarray(x, dtype=np.float32)
    for layer in (0, 2, 4, 6):
        y = y @ weights[f"layer_{layer}_weight"].T + weights[f"layer_{layer}_bias"]
        if layer < 6:
            y = np.maximum(y, 0.0)
    return y.astype(np.float64)
def summary(a,p):
    rel=(p-a)/a; return {"max_relative_error":float(np.abs(rel).max()),"p99_relative_error":float(np.percentile(np.abs(rel),99)),"rms_relative_error":float(np.sqrt(np.mean(rel**2)))}
def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("hef",type=Path); ap.add_argument("model_dir",type=Path); ap.add_argument("samples_csv",type=Path); ap.add_argument("--model",default="PengRobinson"); ap.add_argument("--limit",type=int,default=0); ap.add_argument("--output",type=Path)
    a=ap.parse_args(); meta=json.loads((a.model_dir/"stage_surrogate.json").read_text())
    selected=[]
    with a.samples_csv.open(newline="") as f:
        for r in csv.DictReader(f):
            if r.get("model")==a.model:
                selected.append(r)
    if a.limit>0:selected=selected[:a.limit]
    if not selected:raise ValueError("No matching validation samples")
    rpm_min,rpm_max=meta["speed_min_rpm"],meta["speed_max_rpm"]; mmin,mmax=meta["mass_flow_min_kg_s"],meta["mass_flow_max_kg_s"]
    x=np.array([[2*(float(r["rotational_speed_rpm"])-rpm_min)/(rpm_max-rpm_min)-1,2*(float(r["mass_flow"])-mmin)/(mmax-mmin)-1] for r in selected],dtype=np.float32)
    actual=np.array([[float(r[k]) for k in TARGETS] for r in selected],dtype=np.float64)
    float_scaled = numpy_reference(x, a.model_dir/"stage_surrogate_weights.npz")
    means=np.array(meta["target_log_mean"]); stds=np.array(meta["target_log_std"]); float_values=np.exp(float_scaled*stds+means)
    hef=HEF(str(a.hef)); params=VDevice.create_params()
    t0=time.perf_counter(); outputs=[]
    with VDevice(params) as device:
        # Hailo-10H requires the InferModel API. The older
        # ConfigureParams/VDevice.configure/InferVStreams path is not
        # implemented for this device family (HAILO_NOT_IMPLEMENTED).
        infer_model=device.create_infer_model(str(a.hef))
        infer_model.set_batch_size(1)
        input_model=infer_model.input()
        output_model=infer_model.output()
        # Let HailoRT quantize the normalized float coordinates and
        # dequantize the three network outputs for the reference comparison.
        input_model.set_format_type(FormatType.FLOAT32)
        output_model.set_format_type(FormatType.FLOAT32)
        output_name=output_model.name
        input_shape=tuple(input_model.shape)
        output_shape=tuple(output_model.shape)
        with infer_model.configure() as configured_model:
            for point in x:
                output_buffer=np.empty(output_shape,dtype=np.float32)
                bindings=configured_model.create_bindings(
                    output_buffers={output_name: output_buffer})
                bindings.input().set_buffer(
                    np.ascontiguousarray(point.reshape(input_shape),dtype=np.float32))
                # The no-callback form is synchronous through job.wait().
                job=configured_model.run_async([bindings])
                job.wait(10000)
                outputs.append(output_buffer.reshape(-1)[:3].copy())
    elapsed=time.perf_counter()-t0; quant_scaled=np.asarray(outputs,dtype=np.float64)
    quant_values=np.exp(quant_scaled*stds+means)
    report={"model":a.model,"samples":len(x),"validation_scope":"provided CSV rows; report CV generalization separately","serial_runtime_seconds":elapsed,"serial_function_evaluations_per_second":len(x)/elapsed,"targets":{}}
    for j,k in enumerate(TARGETS):
        report["targets"][k]={"float_vs_direct":summary(actual[:,j],float_values[:,j]),"HEF_vs_direct":summary(actual[:,j],quant_values[:,j]),"HEF_vs_float_quantization_only":summary(float_values[:,j],quant_values[:,j])}
    if a.output:a.output.write_text(json.dumps(report,indent=2)+"\n")
    print(f"Hailo device: {hef.get_network_groups_names() if hasattr(hef,'get_network_groups_names') else 'configured HEF'}")
    print(f"Samples: {len(x)}; serial rate including Python/stream overhead: {len(x)/elapsed:.2f} eval/s")
    for k,m in report["targets"].items():print(f"{k}: HEF-vs-float max quantization rel error={m['HEF_vs_float_quantization_only']['max_relative_error']:.4%}; HEF-vs-direct max rel error={m['HEF_vs_direct']['max_relative_error']:.4%}")
    if a.output:print(f"Report: {a.output}")
    return 0
if __name__=="__main__":raise SystemExit(main())
