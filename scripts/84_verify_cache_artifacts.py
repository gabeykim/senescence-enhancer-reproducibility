# ---------------------------------------------------------------------------
# 84_verify_cache_artifacts.py
#
# Verify the cache loads: shape, finiteness, all-zero rows, fwd != rc; checkpoint loads.
#
# CONSUMES: ois_cache/region_{fwd,rc}.npy, results/head_ois.pt
# PRODUCES: results/artifact_verification.json
# ---------------------------------------------------------------------------
import numpy as np, torch, json, hashlib, os
from pathlib import Path
C = Path("/workspace/ois_enhancer_run/ois_cache")
R = Path("/workspace/ois_enhancer_run/results")
rep = {}
f = np.load(C/"region_fwd.npy", mmap_mode="r")
r = np.load(C/"region_rc.npy",  mmap_mode="r")
print("fwd", f.shape, f.dtype, "| rc", r.shape, r.dtype)
rep["fwd_shape"]=list(f.shape); rep["rc_shape"]=list(r.shape)
rep["fwd_dtype"]=str(f.dtype); rep["rc_dtype"]=str(r.dtype)
assert f.shape==(111671,1920) and r.shape==(111671,1920), "SHAPE MISMATCH"
# stream in chunks to bound RAM
nz_f=nz_r=0; nonfinite_f=nonfinite_r=0; diff_rows=0; maxabs=0.0
sum_absdiff=0.0; n=f.shape[0]; CH=8192
for i in range(0,n,CH):
    a=np.asarray(f[i:i+CH],np.float32); b=np.asarray(r[i:i+CH],np.float32)
    nonfinite_f+=int((~np.isfinite(a)).sum()); nonfinite_r+=int((~np.isfinite(b)).sum())
    nz_f+=int((np.abs(a).sum(1)==0).sum()); nz_r+=int((np.abs(b).sum(1)==0).sum())
    d=np.abs(a-b); diff_rows+=int((d.sum(1)>0).sum()); sum_absdiff+=float(d.sum())
    maxabs=max(maxabs,float(np.abs(a).max()),float(np.abs(b).max()))
rep.update(nonfinite_fwd=nonfinite_f, nonfinite_rc=nonfinite_r,
           all_zero_rows_fwd=nz_f, all_zero_rows_rc=nz_r,
           rows_where_fwd_differs_from_rc=diff_rows, n_rows=n,
           mean_abs_fwd_minus_rc=sum_absdiff/(n*1920), max_abs_value=maxabs)
print(json.dumps(rep,indent=2))
sd=torch.load(R/"head_ois.pt", map_location="cpu")
print("checkpoint keys:", {k:list(v.shape) for k,v in sd.items()})
rep["checkpoint_keys"]={k:list(v.shape) for k,v in sd.items()}
rep["checkpoint_finite"]=bool(all(torch.isfinite(v).all() for v in sd.values()))
rep["checkpoint_weight_norm"]=float(sum(float((v.float()**2).sum()) for v in sd.values())**0.5)
json.dump(rep, open("/workspace/ois_gates/results/artifact_verification.json","w"), indent=2)
print("VERIFY_DONE finite_ckpt=",rep["checkpoint_finite"])
