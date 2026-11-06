#!/usr/bin/env python3
"""Compile the ONNX stage surrogate to a Hailo-10H HEF using DFC 5.4."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from hailo_sdk_client import ClientRunner

def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("onnx",type=Path); ap.add_argument("calibration",type=Path); ap.add_argument("--output-dir",type=Path,required=True); ap.add_argument("--hw-arch",default="hailo10h"); ap.add_argument("--model-script",type=Path)
    a=ap.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
    calib=np.load(a.calibration).astype(np.float32)
    if calib.ndim!=2 or calib.shape[1]!=2: raise ValueError(f"Calibration expected [N,2], got {calib.shape}")
    runner=ClientRunner(hw_arch=a.hw_arch)
    print("[1/4] Parsing ONNX")
    runner.translate_onnx_model(str(a.onnx),"stage_surrogate",start_node_names=["input_coords"],end_node_names=["scaled_log_targets"])
    runner.save_har(str(a.output_dir/"stage_surrogate_parsed.har"))
    if a.model_script:
        print(f"Loading model script: {a.model_script}"); runner.load_model_script(str(a.model_script))
    print(f"[2/4] Optimizing with {len(calib)} calibration samples")
    runner.optimize(calib)
    runner.save_har(str(a.output_dir/"stage_surrogate_optimized.har"))
    print(f"[3/4] Compiling for {a.hw_arch}")
    hef=runner.compile()
    out=a.output_dir/"stage_surrogate.hef"; out.write_bytes(hef)
    print(f"[4/4] Wrote {out}")
    return 0
if __name__=="__main__": raise SystemExit(main())
