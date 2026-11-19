"""Meridian polylines on the BodyVista skin.
usage: python -I paths.py <project_root> <workdir>

Per consecutive acupoint pair: Dijkstra shortest path over the edges of the
(once subdivided) outer skin shell -> resample -> Laplacian smoothing with the
acupoints pinned -> closest-point re-projection onto the skin (5 rounds) ->
lift along the outward normal by LIFT.
"""
import sys, json
import site; sys.path.append(site.getusersitepackages())
sys.stdout.reconfigure(encoding='utf-8')
ROOT, WORK = sys.argv[1], sys.argv[2]
sys.path.insert(0, WORK)
import numpy as np
import trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree
import common
common.init(ROOT)

LIFT = 0.002          # BV units (~1.4 mm real) above the skin
STEP = 0.008          # resampling step, BV units (~5.4 mm real)

data = json.load(open(ROOT + '/data/acupoints.json', encoding='utf-8'))
P = {}
for p in data['points']:
    if p['xyz'] is not None:
        P[(p['code'], p['side'])] = p
MER_ZH = data['points'] and {p['meridian']: p['meridianName'] for p in data['points']}

# Point order per meridian stroke. Originally taken from acupuncture-3d meridians.json,
# with BL39 (between BL38 and BL40) and Yintang GV24+ (between GV24 and GV25) added.
strokes = json.load(open(ROOT + '/data/meridian-strokes.json', encoding='utf-8'))['strokes']
used = {c for m in strokes.values() for s in m for c in s}
all_codes = {p['code'] for p in data['points']}
print('codes not on any stroke:', sorted(all_codes - used))

sk = common.skin_outer()
V, F = sk.vertices, sk.faces
Vs, Fs = trimesh.remesh.subdivide(V, F)
msub = trimesh.Trimesh(Vs, Fs, process=True)
E = msub.edges_unique
L = msub.edges_unique_length
G = coo_matrix((np.r_[L, L], (np.r_[E[:, 0], E[:, 1]], np.r_[E[:, 1], E[:, 0]])),
               shape=(len(msub.vertices),) * 2).tocsr()
tree = cKDTree(msub.vertices)
print('graph nodes', len(msub.vertices), 'edges', len(E))

# outward orientation of face normals: vote with the solver's outward normals
fn = sk.face_normals
votes = []
for p in data['points']:
    if p['normal'] is None:
        continue
    _, _, tri = trimesh.proximity.closest_point(sk, [p['xyz']])
    votes.append(np.dot(fn[tri[0]], p['normal']))
SIGN = 1.0 if np.median(votes) > 0 else -1.0
print('face normal sign', SIGN, 'agree frac', np.mean(np.sign(votes) == SIGN))


def project(X):
    c, d, tri = trimesh.proximity.closest_point(sk, X)
    return c, fn[tri] * SIGN


def geodesic(a, b):
    ia, ib = tree.query(a)[1], tree.query(b)[1]
    chord = np.linalg.norm(b - a)
    if ia == ib:
        return np.array([a, b]), 'chord'
    dist, pred = dijkstra(G, indices=ia, return_predecessors=True, limit=max(4 * chord, 0.1))
    if not np.isfinite(dist[ib]):
        return np.array([a, b]), 'chord_unreachable'
    path = [ib]
    while path[-1] != ia:
        path.append(pred[path[-1]])
    pts = msub.vertices[path[::-1]]
    pts = np.vstack([a, pts[1:-1], b]) if len(pts) > 2 else np.array([a, b])
    return pts, 'geodesic'


def resample(pts, step):
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.r_[0, np.cumsum(seg)]
    if s[-1] < 1e-9:
        return pts[[0, -1]]
    n = max(2, int(np.ceil(s[-1] / step)) + 1)
    t = np.linspace(0, s[-1], n)
    return np.column_stack([np.interp(t, s, pts[:, k]) for k in range(3)])


def surface_segment(a, b):
    pts, how = geodesic(a, b)
    pts = resample(pts, STEP)
    chord = np.linalg.norm(b - a)
    for _ in range(6):
        if len(pts) > 2:
            for _ in range(4):
                pts[1:-1] = 0.5 * pts[1:-1] + 0.25 * (pts[:-2] + pts[2:])
            pts[1:-1], _ = project(pts[1:-1])
        pts = resample(pts, STEP)
        pts[0], pts[-1] = a, b
    length = float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())
    return pts, how, length, chord


out = []
report = []
for mer, sl in strokes.items():
    sides = ['C'] if mer in ('GV', 'CV') else ['L', 'R']
    for side in sides:
        rec = {'meridian': mer, 'name': MER_ZH[mer], 'side': side, 'strokes': []}
        for s in sl:
            codes = [c for c in s if (c, side) in P]
            missing = [c for c in s if (c, side) not in P]
            if missing:
                print('missing points', mer, side, missing)
            allpts, idx = [], []
            for i in range(len(codes) - 1):
                a = np.array(P[(codes[i], side)]['xyz'])
                b = np.array(P[(codes[i + 1], side)]['xyz'])
                pts, how, length, chord = surface_segment(a, b)
                ratio = length / max(chord, 1e-9)
                report.append((mer, side, codes[i], codes[i + 1], how, round(ratio, 2), round(chord / 1.4734 * 1000)))
                if i == 0:
                    idx.append(0)
                    allpts.append(pts)
                else:
                    allpts.append(pts[1:])
                idx.append(sum(len(x) for x in allpts) - 1)
            pts = np.vstack(allpts)
            _, nrm = project(pts)
            lifted = pts + nrm * LIFT
            rec['strokes'].append({'codes': codes, 'acupointIndex': idx,
                                   'points': np.round(lifted, 5).tolist()})
        out.append(rec)

meta = {
    'description': '14 经络在 BodyVista 皮肤表面的折线路径（双侧经分左右，任督为正中 C）',
    'coordinateSystem': data['meta']['coordinateSystem'],
    'method': '相邻穴位之间：在细分一次后的外层皮肤网格边图上求 Dijkstra 最短路，按 0.008 重采样，'
              '固定穴位端点做拉普拉斯平滑并用最近点投影回皮肤（6 轮），最后沿皮肤外法线抬高 0.002（约 1.4 mm）避免 z-fighting。',
    'lift': LIFT, 'step': STEP,
    'strokes': 'strokes 为若干笔画（分支）。codes 为该笔画途经穴位，acupointIndex[i] 是 codes[i] 在 points 中的下标。',
    'orderSource': 'acupuncture-3d meridians.json；补入 BL39（BL38 与 BL40 之间）与印堂 GV24+（GV24 与 GV25 之间）',
    'license': 'CC BY-SA 4.0',
}
json.dump({'meta': meta, 'meridians': out}, open(ROOT + '/data/meridians-paths.json', 'w', encoding='utf-8'),
          ensure_ascii=False)
json.dump(report, open(WORK + '/paths_report.json', 'w'), indent=0)
bad = [r for r in report if r[4] != 'geodesic' or r[5] > 1.8]
print('segments', len(report), 'non-geodesic or ratio>1.8:')
for r in bad:
    print('  ', r)
print('total points', sum(len(s['points']) for m in out for s in m['strokes']))
