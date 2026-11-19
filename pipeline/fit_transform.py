import sys, json
import site; sys.path.append(site.getusersitepackages())
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, sys.argv[2])
import numpy as np
import common
common.init(sys.argv[1])
ROOT = sys.argv[1]
st = json.load(open(ROOT + '/vendor/acupuncture-3d/data/structures.json', encoding='utf-8'))['structures']
sk = common.skin()
print('skin', sk.vertices.shape, sk.faces.shape, 'watertight', sk.is_watertight, sk.bounds)
src, dst, labels = [], [], []
for key, meta in st.items():
    nm = key[:-2] if key.endswith('-r') else key
    eng = nm.replace('-', ' ')
    cand = [eng + ' right', eng]
    pid = None
    for c in cand:
        pid = common.by_english(c)
        if pid: break
    if not pid:
        print('missing', key); continue
    m = common.mesh(pid)
    b = np.array(meta['bounds']).reshape(3, 2)
    fc = b.mean(axis=1)
    bc = m.bounds.mean(axis=0)
    src.append(fc); dst.append(bc); labels.append(key)
    # also add bound extremes as correspondences (min/max per axis), mapping under 180deg Z rot
    fb = b; bb = m.bounds.T  # 3x2
    print(f'{key:45s} F_ext={np.ptp(fb,axis=1).round(4)} B_ext={np.ptp(bb,axis=1).round(4)} B_c={bc.round(3)}')
src = np.array(src); dst = np.array(dst)
def umeyama(A, B):
    ma, mb = A.mean(0), B.mean(0)
    A0, B0 = A - ma, B - mb
    U, S, Vt = np.linalg.svd(B0.T @ A0 / len(A))
    D = np.eye(3)
    if np.linalg.det(U) * np.linalg.det(Vt) < 0: D[2, 2] = -1
    R = U @ D @ Vt
    s = np.trace(np.diag(S) @ D) / (A0 ** 2).sum(1).mean()
    t = mb - s * R @ ma
    return s, R, t
s, R, t = umeyama(src, dst)
res = np.linalg.norm((s * (R @ src.T)).T + t - dst, axis=1)
print('scale', s, '\nR\n', R.round(5), '\nt', t, '\nresid mm (BV units/ s*1000 -> furia mm):')
for l, r in sorted(zip(labels, res), key=lambda x: -x[1]): print(f'  {l:45s} {r/s*1000:.1f} mm')
json.dump({'s': s, 'R': R.tolist(), 't': t.tolist()}, open(sys.argv[2] + '/transform_struct.json', 'w'), indent=1)
