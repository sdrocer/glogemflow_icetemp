"""Six-checks validator, Python port (no node available). Checks 2-5 of dataviz/color-formula."""
import numpy as np, itertools, sys

def hex2rgb(h): h=h.lstrip('#'); return np.array([int(h[i:i+2],16)/255 for i in (0,2,4)])
def srgb2lin(c): return np.where(c<=0.04045, c/12.92, ((c+0.055)/1.055)**2.4)
def lin2srgb(c): return np.where(c<=0.0031308, c*12.92, 1.055*np.clip(c,0,None)**(1/2.4)-0.055)

M1=np.array([[0.4122214708,0.5363325363,0.0514459929],
             [0.2119034982,0.6806995451,0.1073969566],
             [0.0883024619,0.2817188376,0.6299787005]])
M2=np.array([[0.2104542553,0.7936177850,-0.0040720468],
             [1.9779984951,-2.4285922050,0.4505937099],
             [0.0259040371,0.7827717662,-0.8086757660]])
def oklab(rgb):
    lms=M1@srgb2lin(rgb); return M2@np.cbrt(np.clip(lms,0,None))
def oklch(rgb):
    L,a,b=oklab(rgb); return L, float(np.hypot(a,b))

# Machado-Oliveira-Fernandes 2009, severity 1.0 (operate on LINEAR rgb)
PROT=np.array([[0.152286,1.052583,-0.204868],[0.114503,0.786281,0.099216],[-0.003882,-0.048116,1.051998]])
DEUT=np.array([[0.367322,0.860646,-0.227968],[0.280085,0.672501,0.047413],[-0.011820,0.042940,0.968881]])
def simulate(rgb,M): return np.clip(lin2srgb(M@srgb2lin(rgb)),0,1)
def dE(a,b): return float(np.linalg.norm(oklab(a)-oklab(b))*100)

def relL(rgb):
    lin=srgb2lin(rgb); return float(0.2126*lin[0]+0.7152*lin[1]+0.0722*lin[2])
def contrast(a,b):
    l1,l2=sorted([relL(a),relL(b)],reverse=True); return (l1+0.05)/(l2+0.05)

def report(hexes, surface='#fcfcfb', mode='light', pairs='all'):
    cols=[hex2rgb(h) for h in hexes]; surf=hex2rgb(surface); ok=True
    lo,hi=(0.43,0.77) if mode=='light' else (0.48,0.67)
    print(f"--- mode={mode} surface={surface} pairs={pairs} ---")
    for h,c in zip(hexes,cols):
        L,C=oklch(c); cr=contrast(c,surf)
        f2='PASS' if lo<=L<=hi else 'FAIL'; f3='PASS' if C>=0.10 else 'FAIL'
        f5='PASS' if cr>=3.0 else 'WARN'
        if 'FAIL' in (f2,f3): ok=False
        print(f"  {h}  L={L:.3f}[{f2}]  C={C:.3f}[{f3}]  contrast={cr:.2f}:1[{f5}]")
    idx=list(itertools.combinations(range(len(cols)),2)) if pairs=='all' else [(i,i+1) for i in range(len(cols)-1)]
    worst_cvd=(1e9,None); worst_nv=(1e9,None)
    for i,j in idx:
        nv=dE(cols[i],cols[j])
        for name,M in (('prot',PROT),('deut',DEUT)):
            d=dE(simulate(cols[i],M),simulate(cols[j],M))
            if d<worst_cvd[0]: worst_cvd=(d,f"{hexes[i]}/{hexes[j]} ({name})")
        if nv<worst_nv[0]: worst_nv=(nv,f"{hexes[i]}/{hexes[j]}")
    c4='PASS' if worst_cvd[0]>=8 else ('WARN(floor)' if worst_cvd[0]>=6 else 'FAIL')
    c4n='PASS' if worst_nv[0]>=15 else 'FAIL'
    if c4=='FAIL' or c4n=='FAIL': ok=False
    print(f"  worst CVD dE  = {worst_cvd[0]:.1f}  [{c4}]  {worst_cvd[1]}   (target>=8, floor>=6)")
    print(f"  worst normal  = {worst_nv[0]:.1f}  [{c4n}]  {worst_nv[1]}   (hard floor>=15)")
    print(f"  => {'OK' if ok else 'HARD FAIL'}\n"); return ok

if __name__=='__main__':
    hexes=sys.argv[1].split(',')
    a=report(hexes,'#fcfcfb','light','all')
    sys.exit(0 if a else 1)
