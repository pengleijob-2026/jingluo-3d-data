"""Own re-implementation of the anchor geometry (written after reading the
semantics in acupuncture-3d acu3d/anchors.py, rig.py, anatomy.py; no code is
imported from there).

All anchor maths happens in the "Furia frame" of rig_fitted.json
(metres, +X = person's right, +Y = front, +Z = up). The BodyVista skin and
bones are mapped into that frame with the inverse of the fitted similarity
transform, so every number in rig_fitted.json can be used as is.

The left side is computed on a MIRRORED body: BodyVista skin mirrored in X and
BodyVista *left* bones mirrored in X. The anchor logic then always works on
"the right side"; results are mirrored back. Thus each side lands on its own
real skin instead of being a mirror copy of the other.
"""
import json
import math
import numpy as np
import trimesh

import common

MAX_RAY = 0.45
DOUBLE_SKIN = 0.010


def unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])


class Transform:
    def __init__(self, path):
        d = json.load(open(path))
        self.s = float(d['s'])
        self.R = np.array(d['R'])
        self.t = np.array(d['t'])

    def to_bv(self, p):
        p = np.atleast_2d(p)
        return (self.s * (self.R @ p.T)).T + self.t

    def dir_to_bv(self, n):
        n = np.atleast_2d(n)
        return (self.R @ n.T).T

    def to_furia(self, P):
        P = np.atleast_2d(P)
        return (self.R.T @ (P - self.t).T).T / self.s


class Seg:
    def __init__(self, name, p0, p1, cun_length, body, cun_t0=0.0, cun_t1=1.0,
                 cun_lat_m=None, core=False):
        self.name = name
        self.p0 = np.asarray(p0, float)
        self.p1 = np.asarray(p1, float)
        self.cun_length = float(cun_length)
        self.cun_t0, self.cun_t1 = float(cun_t0), float(cun_t1)
        self.cun_lat_m = cun_lat_m
        self.body = body
        self.core = core
        self.e_ax = unit(self.p1 - self.p0)
        ant = np.array([0.0, 1.0, 0.0])
        ant = ant - self.e_ax * ant.dot(self.e_ax)
        if np.linalg.norm(ant) < 1e-6:
            ant = np.array([0.0, 0.0, 1.0]) - self.e_ax * self.e_ax[2]
        self.e_ant = unit(ant)
        self.e_lat = unit(np.cross(self.e_ant, self.e_ax))
        if self.e_lat[0] < 0:
            self.e_lat = -self.e_lat

    @property
    def length(self):
        return float(np.linalg.norm(self.p1 - self.p0))

    @property
    def cun_m(self):
        return self.length * (self.cun_t1 - self.cun_t0) / self.cun_length

    @property
    def cun_lat(self):
        return self.cun_lat_m if self.cun_lat_m else self.cun_m

    def t_from_cun(self, cun, frm='p0'):
        span = self.cun_t1 - self.cun_t0
        if frm == 'p1':
            return self.cun_t1 - span * cun / self.cun_length
        return self.cun_t0 + span * cun / self.cun_length

    def axis_point(self, t):
        t = float(np.clip(t, 0, 1))
        return self.p0 + (self.p1 - self.p0) * t

    def interior(self, t):
        p = self.axis_point(t)
        if self.core:
            p = p.copy()
            p[1] = self.body.core_y(p[2])
            return p
        return p + self.body.section_shift(self.name, t)

    def _limit(self, base, d):
        if abs(self.p0[0]) < 0.02 or abs(d[0]) < 1e-6 or d[0] * base[0] >= 0:
            return MAX_RAY
        return float(np.clip(-base[0] / d[0], 0.01, MAX_RAY))

    def surface_by_azimuth(self, t, az, out=0.0):
        a = math.radians(az)
        d = unit(math.cos(a) * self.e_ant + math.sin(a) * self.e_lat)
        base = self.interior(t)
        p = self.body.cast(base, d, self._limit(base, d))
        return p + d * out * self.cun_lat

    def surface_by_lateral(self, t, lat, face='anterior', out=0.0):
        base = self.interior(t) + self.e_lat * lat * self.cun_lat
        if face == 'posterior':
            d = -self.e_ant
        elif face == 'lateral':
            d = self.e_lat if lat >= 0 else -self.e_lat
            base = self.interior(t)
        else:
            d = self.e_ant
        p = self.body.cast(base, d, self._limit(base, d))
        return p + d * out * self.cun_lat


class Struct:
    def __init__(self, name, centroid, axis, extent, mesh=None):
        self.name = name
        self.centroid = np.asarray(centroid, float)
        self.axis = unit(axis)
        self.extent = tuple(extent)
        self.mesh = mesh

    def at(self, along, ref_dir=None):
        axis, (lo, hi) = self.axis, self.extent
        if ref_dir is not None and axis.dot(ref_dir) < 0:
            axis, lo, hi = -axis, -hi, -lo
        return self.centroid + axis * (lo + (hi - lo) * float(np.clip(along, 0, 1)))


def pca_meta(V):
    c = V.mean(0)
    X = V - c
    w, U = np.linalg.eigh(X.T @ X)
    ax = U[:, -1]
    pr = X @ ax
    return c, ax, (float(pr.min()), float(pr.max()))


class Body:
    LIMBS = ('upper_arm', 'forearm', 'thigh', 'shank')
    SECTION_MAX_R = {'upper_arm': 0.085, 'forearm': 0.085, 'thigh': 0.115, 'shank': 0.080}
    SECTION_MAX_OFF = 0.035
    SECTION_T = np.linspace(0, 1, 41)
    SHORT_RAY = {'hand': (0.045, 0.006), 'foot': (0.065, 0.008)}
    SAGITTAL_CUN = 18.0
    AURICLE = ('helix', 'antihelix', 'concha-of-auricle', 'tragus', 'antitragus',
               'lobule-of-auricle', 'apex-of-auricle', 'anterior-notch-of-auricle')

    def __init__(self, root, workdir, side='R'):
        self.root, self.side = root, side
        self.mirror = side == 'L'
        A = root + '/vendor/acupuncture-3d/data/'
        self.fit = json.load(open(A + 'rig_fitted.json', encoding='utf-8'))
        self.smeta = json.load(open(A + 'structures.json', encoding='utf-8'))['structures']
        self.T = Transform(workdir + '/transform_struct.json')
        self.H = float(self.fit['height'])
        sk = common.skin_outer()
        V = self.T.to_furia(sk.vertices)
        F = sk.faces.copy()
        if self.mirror:
            V[:, 0] *= -1
            F = F[:, ::-1]
        self.skin = trimesh.Trimesh(V, F, process=False)
        self.ray = trimesh.ray.ray_triangle.RayMeshIntersector(self.skin)
        self.flags = []
        rng = np.random.default_rng(7)
        d = rng.normal(size=(7, 3))
        self.parity_dirs = d / np.linalg.norm(d, axis=1, keepdims=True)
        self._structs = {}
        self._frames = {}
        self._bone_cache = {}
        self._edges = {}
        self.cun_trunk_lat = float(self.fit['cun']['trunk_lateral_m'])
        self._core = self._build_core()
        self._ruler_z = [k[0] for k in self.fit['trunk_ruler']['knots']]
        self._ruler_c = [k[1] for k in self.fit['trunk_ruler']['knots']]
        self.vertebra_z = {k: v['z'] for k, v in self.fit['vertebra'].items()}
        self.landmarks = {k: np.array(v) for k, v in self.fit['landmarks'].items()}
        self.sagittal = self._build_sagittal()
        self.levels = self._levels()
        self.segments = self._build_segments()
        self._soft = {}
        for n in self.LIMBS:
            self._soft[n] = self._section_offsets(self.segments[n])
        self._soft['hand'] = self._bone_shift_region('hand')
        self.landmark_cun = {'pubis': 0.0, 'umbilicus': 5.0, 'xiphoid': 13.0,
                             'sternal_notch': 22.0,
                             'nipple': self.cun_of_z(self.levels['nipple']),
                             'iliac_crest': self.cun_of_z(self.levels['iliac_crest'])}

    # ---------------- meshes ----------------
    def bone_mesh(self, name):
        """BodyVista bone mesh for this side, in (possibly mirrored) Furia frame."""
        key = name.lower()
        if key in self._bone_cache:
            return self._bone_cache[key]
        sidew = 'left' if self.mirror else 'right'
        pid = common.by_english(f'{key} {sidew}') or common.by_english(key)
        if pid is None:
            self._bone_cache[key] = None
            return None
        m = common.mesh(pid)
        V = self.T.to_furia(m.vertices)
        F = m.faces.copy()
        if self.mirror:
            V[:, 0] *= -1
            F = F[:, ::-1]
        mm = trimesh.Trimesh(V, F, process=False)
        self._bone_cache[key] = mm
        return mm

    # ---------------- ray helpers ----------------
    def hits(self, origin, d, limit):
        loc, ri, _ = self.ray.intersects_location([origin], [d], multiple_hits=True)
        if len(loc) == 0:
            return np.zeros((0, 3))
        dist = (loc - origin) @ d
        ok = (dist > 1e-7) & (dist <= limit)
        loc, dist = loc[ok], dist[ok]
        return loc[np.argsort(dist)]

    def is_inside(self, p):
        p = np.asarray(p, float)
        origins = np.repeat(p[None], 7, axis=0)
        loc, ri, _ = self.ray.intersects_location(origins, self.parity_dirs, multiple_hits=True)
        cnt = np.bincount(ri, minlength=7) if len(ri) else np.zeros(7, int)
        return int((cnt % 2).sum()) >= 4

    def closest(self, p):
        c, d, tri = trimesh.proximity.closest_point(self.skin, np.atleast_2d(p))
        return c[0], float(d[0]), int(tri[0])

    def cast(self, origin, d, limit=None, probe=0.03):
        d = unit(d)
        lim = MAX_RAY if limit is None else limit
        pts = self.hits(np.asarray(origin, float), d, lim)
        if len(pts):
            near = pts[np.linalg.norm(pts - pts[0], axis=1) <= DOUBLE_SKIN]
            for h in near:
                if not self.is_inside(h + d * 0.002):
                    return h
            return pts[0]
        self.flags.append('ray_miss_closest_point')
        c, _, _ = self.closest(origin + d * probe)
        return c

    def raw_normal(self, p):
        c, _, tri = self.closest(p)
        return self.skin.face_normals[tri], c

    def surface_normal(self, p):
        n, foot = self.raw_normal(p)
        return -n if self.is_inside(foot + n * 0.003) else n

    def walk_on_skin(self, p, direction, distance, step=0.002):
        if abs(distance) < 1e-9:
            return np.asarray(p, float)
        pos = np.asarray(p, float)
        d = np.asarray(direction, float) * np.sign(distance)
        left = abs(distance)
        while left > 1e-6:
            h = min(step, left)
            n, _ = self.raw_normal(pos)
            d = unit(d - n * d.dot(n))
            foot, _, _ = self.closest(pos + d * h)
            d = unit(foot - pos)
            pos = foot
            left -= h
        return pos

    def scalp_side(self, p, forward):
        n, _ = self.raw_normal(p)
        side = np.cross(forward, n)
        side = unit(side - n * side.dot(n))
        return side if side[0] > 0 else -side

    # ---------------- core / rulers ----------------
    def _build_core(self):
        p = self.skin.vertices
        band = p[np.abs(p[:, 0]) < 0.035]
        zs = np.linspace(band[:, 2].min(), band[:, 2].max(), 160)
        ys = []
        for z in zs:
            s = band[np.abs(band[:, 2] - z) < 0.012]
            ys.append(float((s[:, 1].min() + s[:, 1].max()) / 2) if len(s) > 6
                      else (ys[-1] if ys else 0.0))
        ys = np.convolve(np.array(ys), np.ones(5) / 5, mode='same')
        ys[:2], ys[-2:] = ys[2], ys[-3]
        return zs, ys

    def core_y(self, z):
        return float(np.interp(z, *self._core))

    def z_of_cun(self, c):
        return float(np.interp(c, self._ruler_c, self._ruler_z))

    def cun_of_z(self, z):
        return float(np.interp(z, self._ruler_z, self._ruler_c))

    def _first_hit(self, a, b):
        d = b - a
        L = np.linalg.norm(d)
        h = self.hits(a, d / L, L)
        return h[0] if len(h) else None

    def _build_sagittal(self):
        g = self.landmarks['glabella']
        c7 = np.array(self.fit['vertebra']['C7']['spinous'], float)
        center = np.array([0.0, (g[1] + c7[1]) / 2, g[2]])
        a0 = math.atan2(g[2] - center[2], g[1] - center[1])
        a1 = math.atan2(c7[2] - center[2], c7[1] - center[1])
        if a1 <= a0:
            a1 += 2 * math.pi
        c7_skin = self._first_hit(c7, c7 + np.array([0, -MAX_RAY, 0]))
        if c7_skin is None:
            c7_skin = c7
        pts = []
        for a in np.linspace(a0, a1 + 0.2, 480):
            d = np.array([0.0, math.cos(a), math.sin(a)])
            hit = self._first_hit(center, center + d * MAX_RAY)
            if hit is None:
                continue
            if a > a0 + math.pi and hit[2] <= c7_skin[2]:
                break
            pts.append(hit)
        pts = np.array(pts + [c7_skin])
        s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
        knots_s, knots_c = [0.0, s[-1]], [0.0, self.SAGITTAL_CUN]
        apex = self.smeta['apex-of-auricle-r']
        apex_y = apex['centroid'][1]
        over = np.flatnonzero(pts[:, 2] > np.percentile(pts[:, 2], 70))
        k = int(over[np.argmin(np.abs(pts[over, 1] - apex_y))])
        knots_s, knots_c = [0.0, s[k], s[-1]], [0.0, 8.0, self.SAGITTAL_CUN]
        return {'points': pts, 'cun': np.interp(s, knots_s, knots_c)}

    def sagittal_point(self, c):
        sg = self.sagittal
        c = float(np.clip(c, 0, self.SAGITTAL_CUN))
        return np.array([np.interp(c, sg['cun'], sg['points'][:, i]) for i in range(3)])

    def _levels(self):
        lm = self.landmarks
        lv = {k: float(v[2]) for k, v in lm.items()}
        lv['nipple'] = lv['xiphoid'] + (lv['sternal_notch'] - lv['xiphoid']) * 0.42
        lv['iliac_crest'] = float(self.vertebra_z['L4'])
        lv['hairline_ant'] = float(self.sagittal_point(3.0)[2])
        lv['hairline_post'] = float(self.sagittal_point(15.0)[2])
        return lv

    def _build_segments(self):
        f = self.fit['segments']
        segs = {}
        for n in ('upper_arm', 'forearm', 'hand', 'thigh', 'shank', 'foot'):
            s = f[n]
            segs[n] = Seg(n, s['p0'], s['p1'], s['cun_length'], self, cun_t0=s.get('cun_t0', 0.0))
        lv = self.levels
        segs['trunk'] = Seg('trunk', [0, 0, lv['pubis'] - 0.05], [0, 0, lv['sternal_notch'] + 0.03],
                            22.0, self, cun_lat_m=self.cun_trunk_lat, core=True)
        segs['neck'] = Seg('neck', [0, 0, lv['sternal_notch']], [0, 0, lv['chin'] + 0.004],
                           4.0, self, cun_lat_m=self.cun_trunk_lat, core=True)
        segs['head'] = Seg('head', [0, 0, lv['chin']], [0, 0, lv['vertex']],
                           10.0, self, cun_lat_m=0.0125 * self.H, core=True)
        return segs

    def seg(self, n):
        return self.segments[n]

    def ruler_t(self, name, cun, frm='p0'):
        r = self.fit['rulers'][name]
        seg = self.seg(r['seg'])
        t0 = float(np.dot(np.asarray(r['p0']) - seg.p0, seg.e_ax)) / seg.length
        t1 = float(np.dot(np.asarray(r['p1']) - seg.p0, seg.e_ax)) / seg.length
        f = cun / float(r['cun_length'])
        return t1 + (t0 - t1) * f if frm == 'p1' else t0 + (t1 - t0) * f

    # ---------------- skin sections / soft-tissue offsets ----------------
    def skin_loop(self, center, normal):
        try:
            sec = self.skin.section(plane_normal=normal, plane_origin=center)
        except Exception:
            sec = None
        if sec is None:
            return None
        best, bd = None, 1e9
        for loop in sec.discrete:
            P = np.asarray(loop)
            if len(P) < 6:
                continue
            d = np.linalg.norm(P.mean(0) - center)
            if d < bd:
                best, bd = P, d
        return best

    def _section_offsets(self, seg):
        off = np.full((len(self.SECTION_T), 3), np.nan)
        for i, t in enumerate(self.SECTION_T):
            c = seg.axis_point(t)
            loop = self.skin_loop(c, seg.e_ax)
            if loop is None:
                continue
            sc = loop.mean(0)
            radius = np.max(np.linalg.norm(loop - sc, axis=1))
            d = sc - c
            d -= seg.e_ax * d.dot(seg.e_ax)
            if radius < self.SECTION_MAX_R[seg.name] and np.linalg.norm(d) < self.SECTION_MAX_OFF:
                off[i] = d
        ok = ~np.isnan(off[:, 0])
        if not ok.any():
            return np.zeros_like(off)
        idx = np.arange(len(off))
        first, last = int(idx[ok][0]), int(idx[ok][-1])
        for k in range(3):
            off[:, k] = np.interp(idx, idx[ok], off[ok, k])
        if first:
            off[:first] *= np.linspace(0, 1, first + 1)[:-1, None]
        pad = np.pad(off, ((2, 2), (0, 0)), mode='edge')
        out = np.column_stack([np.convolve(pad[:, k], np.ones(5) / 5, mode='valid') for k in range(3)])
        out[last:] = out[last]
        return out

    def section_shift(self, name, t):
        off = self._soft.get(name) if hasattr(self, '_soft') else None
        if off is None or np.ndim(off) != 2:
            return np.zeros(3)
        return np.array([np.interp(t, self.SECTION_T, off[:, k]) for k in range(3)])

    def _bone_shift_region(self, region):
        f = self.frame(region)
        refs = [self.struct(f'{o} metacarpal bone', 'hand') for o in ('Second', 'Third', 'Fourth')]
        pts = np.vstack([r.mesh.vertices for r in refs])
        shifts = []
        for along in (0.3, 0.5, 0.7):
            base = np.mean([r.at(along, f['long']) for r in refs], axis=0)
            loop = self.skin_loop(base, f['long'])
            if loop is None:
                continue
            s = pts[np.abs((pts - base) @ f['long']) < 0.006]
            if len(s) < 5:
                continue
            shifts.append(float(np.max(loop @ f['dorsal'])) - 0.004 - float(np.max(s @ f['dorsal'])))
        if not shifts:
            return np.zeros(3)
        return f['dorsal'] * min(0.0, float(np.median(shifts)))

    def soft_shift(self, region):
        v = self._soft.get(region)
        return np.zeros(3) if v is None or np.ndim(v) != 1 else v

    # ---------------- hand / foot registration ----------------
    HAND_BONES = ['scaphoid bone', 'lunate bone', 'triquetral bone', 'pisiform bone', 'trapezium bone',
                  'trapezoid bone', 'capitate bone', 'hamate bone'] + \
        [f'{o} metacarpal bone' for o in ('first', 'second', 'third', 'fourth', 'fifth')] + \
        [f'{k} phalanx of {o} finger of hand' for o in ('first', 'second', 'third', 'fourth', 'fifth')
         for k in ('proximal', 'middle', 'distal')]
    FOOT_BONES = ['calcaneus', 'talus', 'navicular bone', 'cuboid bone', 'medial cuneiform bone',
                  'intermediate cuneiform bone', 'lateral cuneiform bone'] + \
        [f'{o} metatarsal bone' for o in ('first', 'second', 'third', 'fourth', 'fifth')] + \
        [f'{k} phalanx of {o} finger of foot' for o in ('first', 'second', 'third', 'fourth', 'fifth')
         for k in ('proximal', 'middle', 'distal')]
    DIGITS = ('first', 'second', 'third', 'fourth', 'fifth')

    @staticmethod
    def _kabsch(A, B):
        ca, cb = A.mean(0), B.mean(0)
        U, S, Vt = np.linalg.svd((A - ca).T @ (B - cb))
        D = np.eye(3)
        if np.linalg.det(Vt.T @ U.T) < 0:
            D[2, 2] = -1
        R = Vt.T @ D @ U.T
        return R, cb - R @ ca

    def _icp(self, src, target, iters=40, translate_only=False):
        """Rigid ICP: bone vertices -> nearest skin surface. Returns 4x4."""
        M = np.eye(4)
        X = src.copy()
        for _ in range(iters):
            c, d, _ = trimesh.proximity.closest_point(target, X)
            if translate_only:
                R, t = np.eye(3), (c - X).mean(0)
            else:
                R, t = self._kabsch(X, c)
            X = (R @ X.T).T + t
            step = np.eye(4)
            step[:3, :3] = R
            step[:3, 3] = t
            M = step @ M
            if np.linalg.norm(t) < 1e-5 and abs(np.trace(R) - 3) < 1e-8:
                break
        return M

    def _contain_translation(self, src, target, margin=0.003):
        """Translation keeping the bone vertices >= margin inside the skin:
        minimise sum(max(0, signed_dist + margin)^2), coarse-to-fine grid."""
        fn = target.face_normals

        def cost(t):
            X = src + t
            c, d, tri = trimesh.proximity.closest_point(target, X)
            sd = np.einsum('ij,ij->i', X - c, fn[tri])
            return float(np.sum(np.maximum(0.0, sd + margin) ** 2))
        best, bc = np.zeros(3), cost(np.zeros(3))
        for step, span in ((0.010, 3), (0.004, 2), (0.0015, 2)):
            c0 = best.copy()
            for dx in range(-span, span + 1):
                for dy in range(-span, span + 1):
                    for dz in range(-span, span + 1):
                        t = c0 + step * np.array([dx, dy, dz])
                        v = cost(t)
                        if v < bc:
                            best, bc = t, v
        M = np.eye(4)
        M[:3, 3] = best
        return M

    def _region_skin(self, faces_mask):
        return self.skin.submesh([np.flatnonzero(faces_mask)], append=True, repair=False)

    def _register(self, region):
        """Register the Z-Anatomy hand/foot skeleton into the BodyParts3D skin:
        global rigid ICP of the whole group, then per-digit ICP of the phalanges
        (+ first metacarpal for the thumb)."""
        from scipy.spatial import cKDTree
        names = self.HAND_BONES if region == 'hand' else self.FOOT_BONES
        meshes = {n: self.bone_mesh(n) for n in names}
        meshes = {n: m for n, m in meshes.items() if m is not None}
        C = self.skin.triangles_center
        if region == 'hand':
            W = np.array(self.fit['segments']['forearm']['p1'])
            ax = self.seg('hand').e_ax
            mask = ((C - W) @ ax > -0.015) & (np.linalg.norm(C - W, axis=1) < 0.26) & (C[:, 0] > 0.08)
        else:
            mask = (C[:, 2] < 0.075) & (C[:, 0] > 0.0) & (C[:, 0] < 0.25)
        # keep only the connected skin patch of this hand/foot (the fingers lie
        # next to the thigh, whose skin is in the same half-space)
        idx = np.flatnonzero(mask)
        adj = self.skin.face_adjacency
        keep = np.zeros(len(C), bool)
        keep[idx] = True
        adj = adj[keep[adj[:, 0]] & keep[adj[:, 1]]]
        comps = trimesh.graph.connected_components(adj, nodes=idx, min_len=1)
        ref = np.array(self.fit['segments']['forearm']['p1']) if region == 'hand' else self.seg('foot').axis_point(0.5)
        best = min(comps, key=lambda c: np.min(np.linalg.norm(C[c] - ref, axis=1)))
        mask = np.zeros(len(C), bool)
        mask[best] = True
        target = self._region_skin(mask)
        rng = np.random.default_rng(1)
        allV = np.vstack([m.vertices for m in meshes.values()])
        src = allV[rng.choice(len(allV), min(4000, len(allV)), replace=False)]
        # init: rotate the bone long axis onto the skin long axis about the centroids
        def long_axis(P):
            X = P - P.mean(0)
            w, U = np.linalg.eigh(X.T @ X)
            return U[:, -1]
        sv = target.vertices
        a_b, a_s = long_axis(allV), long_axis(sv)
        if a_b.dot(a_s) < 0:
            a_s = -a_s
        v = np.cross(a_b, a_s)
        sn, cs = np.linalg.norm(v), a_b.dot(a_s)
        if sn > 1e-9:
            K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]) / sn
            ang = np.arctan2(sn, cs)
            R0 = np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K
        else:
            R0 = np.eye(3)
        M0 = np.eye(4)
        M0[:3, :3] = R0
        M0[:3, 3] = sv.mean(0) - R0 @ allV.mean(0)
        if region == 'foot':
            # foot: the shank analysis shows a pure translation (~4 cm) between the
            # Z-Anatomy leg skeleton and the BodyParts3D skin; a free rotation
            # pitches the foot and drags the talus/malleolar region down.
            M0 = np.eye(4)
            Mg = self._contain_translation(src[::2], target)
        else:
            src0 = (R0 @ src.T).T + M0[:3, 3]
            Mg = self._icp(src0, target) @ M0
        T = {n: Mg for n in meshes}
        before = float(np.median(trimesh.proximity.closest_point(target, src)[1]))
        srcg = (Mg[:3, :3] @ src.T).T + Mg[:3, 3]
        after = float(np.median(trimesh.proximity.closest_point(target, srcg)[1]))
        cen = allV.mean(0)
        info = {'global_shift_mm': np.round((Mg[:3, :3] @ cen + Mg[:3, 3] - cen) * 1000, 1).tolist(),
                'global_rot_deg': round(float(np.degrees(np.arccos(np.clip((np.trace(Mg[:3, :3]) - 1) / 2, -1, 1)))), 1),
                'median_bone_skin_mm_before_after': (round(before * 1000, 1), round(after * 1000, 1)),
                'digits': {}}
        kind = 'hand' if region == 'hand' else 'foot'
        if region == 'hand':
            tips = self._skin_fingertips(target, np.array(self.fit['segments']['forearm']['p1']))
            info['skin_tips'] = None if tips is None else np.round(tips, 4).tolist()
            if tips is not None:
                for k, dg in enumerate(self.DIGITS):
                    grp = [n for n in meshes if f'of {dg} finger of hand' in n]
                    if dg == 'first':
                        grp.append('first metacarpal bone')
                    W0 = np.array(self.fit['segments']['forearm']['p1'])
                    base_bone = 'first metacarpal bone' if dg == 'first' else f'proximal phalanx of {dg} finger of hand'
                    tip_bone = f'distal phalanx of {dg} finger of hand'
                    if base_bone not in meshes or tip_bone not in meshes:
                        continue
                    Bv = (Mg[:3, :3] @ meshes[base_bone].vertices.T).T + Mg[:3, 3]
                    c_, ax_, ext_ = pca_meta(Bv)
                    ends = [c_ + ax_ * ext_[0], c_ + ax_ * ext_[1]]
                    pivot = min(ends, key=lambda e: np.linalg.norm(e - W0))
                    Tv = (Mg[:3, :3] @ meshes[tip_bone].vertices.T).T + Mg[:3, 3]
                    btip = Tv[np.argmax(np.linalg.norm(Tv - pivot, axis=1))]
                    tgt = tips[k] - unit(tips[k] - pivot) * 0.004
                    u, w = unit(btip - pivot), unit(tgt - pivot)
                    v = np.cross(u, w)
                    sn, cs = np.linalg.norm(v), float(u.dot(w))
                    R = np.eye(3)
                    if sn > 1e-9:
                        K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]) / sn
                        a = np.arctan2(sn, cs)
                        R = np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K
                    sc = float(np.linalg.norm(tgt - pivot) / np.linalg.norm(btip - pivot))
                    Md = np.eye(4)
                    Md[:3, :3] = sc * R
                    Md[:3, 3] = pivot - sc * R @ pivot
                    info['digits'][dg] = ('aimed', round(float(np.degrees(np.arctan2(sn, cs))), 1),
                                          round(float(np.linalg.norm(tgt - pivot) - np.linalg.norm(btip - pivot)) * 1000, 1))
                    for n in grp:
                        T[n] = Md @ Mg
                T['_global'] = Mg
                return T, info
        for dg in self.DIGITS:
            grp = [n for n in meshes if f'of {dg} finger of {kind}' in n]
            if dg == 'first' and region == 'hand':
                grp.append('first metacarpal bone')
            if not grp:
                continue
            V = np.vstack([(Mg[:3, :3] @ meshes[n].vertices.T).T + Mg[:3, 3] for n in grp])
            dist, _ = cKDTree(V).query(C)
            dmask = mask & (dist < 0.025)
            if dmask.sum() < 20:
                continue
            sub = V[::max(1, len(V) // 1500)]
            tgt = self._region_skin(dmask)
            Md = self._icp(sub, tgt, iters=30)
            vc = V.mean(0)
            ang = float(np.degrees(np.arccos(np.clip((np.trace(Md[:3, :3]) - 1) / 2, -1, 1))))
            if ang > 35:
                Md = self._icp(sub, tgt, iters=30, translate_only=True)
                ang = 0.0
            info['digits'][dg] = (np.round((Md[:3, :3] @ vc + Md[:3, 3] - vc) * 1000, 1).tolist(), round(ang, 1))
            for n in grp:
                T[n] = Md @ Mg
        T['_global'] = Mg
        return T, info

    def _skin_fingertips(self, target, W):
        """5 fingertips of the skin hand = local maxima (r=22 mm) of the distance
        from the wrist, farther than 10 cm; ordered radial (thumb) -> ulnar."""
        from scipy.spatial import cKDTree
        V = target.vertices
        d = np.linalg.norm(V - W, axis=1)
        tree = cKDTree(V)
        tips = [i for i, v in enumerate(V)
                if d[i] > 0.10 and d[i] >= d[tree.query_ball_point(v, 0.022)].max() - 1e-12]
        if len(tips) != 5:
            return None
        T = V[tips]
        radial = self.frame_raw_radial()
        return T[np.argsort(-(T @ radial))]

    def frame_raw_radial(self):
        """Radial direction from the unregistered metacarpals (I minus V)."""
        a = self.bone_mesh('first metacarpal bone').vertices.mean(0)
        b = self.bone_mesh('fifth metacarpal bone').vertices.mean(0)
        return unit(a - b)

    def region_T(self, region):
        key = 'reg_' + region
        if key not in self._edges:
            self._edges[key], self._edges[key + '_info'] = self._register(region)
        return self._edges[key]

    def region_M(self, key, region):
        T = self.region_T(region)
        return T.get(key, T['_global'])

    def region_mesh(self, name, region):
        m = self.bone_mesh(name)
        if m is None or region not in ('hand', 'foot'):
            return m
        mm = m.copy()
        if region == 'foot' and name.lower() in ('tibia', 'fibula'):
            # malleoli belong to the leg bones: use the shank bone-to-skin shift
            mm.apply_translation(self.bone_shift('shank'))
            return mm
        mm.apply_transform(self.region_M(name.lower(), region))
        return mm

    # ---------------- structures ----------------
    def struct(self, name, region=None):
        key = name.strip().lower().replace('(', '').replace(')', '')
        ck = (key, region if region in ('hand', 'foot') else None)
        if ck in self._structs:
            return self._structs[ck]
        slug = '-'.join(w for w in ''.join(ch if ch.isalnum() else ' ' for ch in key).split())
        meta = self.smeta.get(slug + '-r') or self.smeta.get(slug)
        mesh = self.region_mesh(key, ck[1])
        if mesh is not None and len(mesh.vertices) >= 4:
            c, ax, ext = pca_meta(mesh.vertices)
            ref_axis = np.asarray(meta['axis']) if meta is not None else None
            if ref_axis is not None and ck[1] is not None and key not in ('tibia', 'fibula'):
                ref_axis = self.region_M(key, ck[1])[:3, :3] @ ref_axis
            if ref_axis is not None and ax.dot(ref_axis) < 0:
                ax, ext = -ax, (-ext[1], -ext[0])
            ref = Struct(name, c, ax, ext, mesh)
        elif meta is not None:
            # no BodyVista mesh (auricle parts): use the Z-Anatomy metadata.
            # For the mirrored (left) body this assumes a symmetric auricle.
            ref = Struct(name, meta['centroid'], meta['axis'], meta['extent'], None)
            self.flags.append('struct_from_metadata')
        else:
            raise KeyError(name)
        self._structs[ck] = ref
        return ref

    def frame(self, region):
        if region in self._frames:
            return self._frames[region]
        ordn = ('First', 'Second', 'Third', 'Fourth', 'Fifth')
        bone = 'metacarpal' if region == 'hand' else 'metatarsal'
        refs = [self.struct(f'{o} {bone} bone', region) for o in ordn]
        base = np.mean([r.at(0.0) for r in refs], axis=0)
        tip = np.mean([r.at(1.0) for r in refs], axis=0)
        long = unit(tip - base)
        ff = refs[0].centroid - refs[4].centroid
        side = unit(ff - long * ff.dot(long))
        normal = unit(np.cross(side, long))
        if region == 'hand':
            pis = self.struct('Pisiform bone', 'hand').centroid
            if (pis - base).dot(normal) > 0:
                normal = -normal
            f = {'long': long, 'radial': side, 'ulnar': -side, 'dorsal': normal, 'palmar': -normal}
        else:
            if normal[2] < 0:
                normal = -normal
            f = {'long': long, 'medial': side, 'lateral': -side, 'dorsal': normal, 'plantar': -normal}
        f['proximal'], f['distal'] = -long, long
        self._frames[region] = f
        return f

    def direction(self, name, region=None):
        if region and name in self.frame(region):
            return self.frame(region)[name]
        base = {'up': [0, 0, 1], 'down': [0, 0, -1], 'anterior': [0, 1, 0], 'posterior': [0, -1, 0],
                'lateral': [1, 0, 0], 'medial': [-1, 0, 0], 'dorsal': [0, -1, 0], 'palmar': [0, 1, 0]}
        return np.array(base[name], float)

    def cast_side(self, base, d, region=None):
        limit = MAX_RAY
        if abs(base[0]) > 0.02 and abs(d[0]) > 1e-6 and d[0] * base[0] < 0:
            limit = float(np.clip(-base[0] / d[0], 0.01, MAX_RAY))
        if region in self.SHORT_RAY:
            reach, probe = self.SHORT_RAY[region]
            return self.cast(base, d, min(limit, reach), probe)
        return self.cast(base, d, limit)

    # ---------------- bone edges ----------------
    BONE_EDGES = {
        'tibial_crest': ('shank', 'Tibia', 'anterior'),
        'tibia_medial': ('shank', 'Tibia', 'medial'),
        'fibula_anterior': ('shank', 'Fibula', 'anterior'),
        'fibula_posterior': ('shank', 'Fibula', 'posterior'),
        'ulna_posterior': ('forearm', 'Ulna', 'posterior'),
    }
    EDGE_LEVELS = np.linspace(0.02, 0.98, 60)
    BONE_FIT_PROBES = {
        'shank': [('Tibia', 'anterior', 4.0), ('Tibia', 'medial', 4.0),
                  ('Fibula', 'lateral', 7.0), ('Tibia', 'posterior', 25.0)],
        'forearm': [('Ulna', 'posterior', 5.0), ('Radius', 'lateral', 4.0)],
    }
    PALPABLE_MM = 8.0

    def _slice_pts(self, mesh, normal, origin):
        try:
            sec = mesh.section(plane_normal=normal, plane_origin=origin)
        except Exception:
            return np.zeros((0, 3))
        if sec is None:
            return np.zeros((0, 3))
        return np.asarray(sec.vertices)

    def _edge_points(self, seg, mesh, direction, ts=None):
        out = []
        for t in (np.linspace(0.25, 0.8, 10) if ts is None else ts):
            c = seg.interior(t)
            P = self._slice_pts(mesh, seg.e_ax, c)
            if len(P) < 8:
                continue
            v = P - c
            v -= np.outer(v @ seg.e_ax, seg.e_ax)
            out.append(P[int(np.argmax(v @ direction))])
        return np.array(out)

    def bone_shift(self, seg_name):
        key = 'shift_' + seg_name
        if key in self._edges:
            return self._edges[key]
        probes = self.BONE_FIT_PROBES.get(seg_name)
        seg = self.seg(seg_name)
        pts, want = [], []
        for st, side, expect in probes:
            d = {'anterior': seg.e_ant, 'posterior': -seg.e_ant,
                 'medial': -seg.e_lat, 'lateral': seg.e_lat}[side]
            P = self._edge_points(seg, self.struct(st).mesh, d)
            if len(P):
                pts.append(P)
                want.append(expect)

        def err(shift):
            tot = 0.0
            for P, e in zip(pts, want):
                c, dist, _ = trimesh.proximity.closest_point(self.skin, P + shift)
                tot += (float(np.median(dist)) * 1000 - e) ** 2
            return tot
        zero = err(np.zeros(3))
        best, be = np.zeros(3), zero
        if (zero / max(len(pts), 1)) ** 0.5 > 5.0:
            grid = np.arange(-0.045, 0.0451, 0.003)
            for da in grid:
                for dl in grid:
                    sh = seg.e_ant * da + seg.e_lat * dl
                    e = err(sh)
                    if e < be:
                        best, be = sh, e
        self._edges[key] = best
        self._edges[key + '_rms'] = ((be / max(len(pts), 1)) ** 0.5, (zero / max(len(pts), 1)) ** 0.5)
        return best

    def _palpable_border(self, seg, mesh, start_az, direction, t, shift):
        prev = start_az
        for step in range(1, 60):
            az = start_az + direction * step * 3.0
            p = seg.surface_by_azimuth(t, az)
            P = self._slice_pts(mesh, seg.e_ax, seg.interior(t) - shift)
            if len(P) < 4:
                return prev
            if np.min(np.linalg.norm(P + shift - p, axis=1)) * 1000 > self.PALPABLE_MM:
                return prev
            prev = az
        return prev

    def _edge_line(self, ref):
        seg_name, st, side = self.BONE_EDGES[ref]
        seg = self.seg(seg_name)
        mesh = self.struct(st).mesh
        pick = {'anterior': seg.e_ant, 'posterior': -seg.e_ant,
                'medial': -seg.e_lat, 'lateral': seg.e_lat}[side]
        shift = self.bone_shift(seg_name)
        out = np.full(len(self.EDGE_LEVELS), np.nan)
        for i, t in enumerate(self.EDGE_LEVELS):
            c = seg.interior(t)
            P = self._slice_pts(mesh, seg.e_ax, c - shift)
            if len(P) < 8:
                continue
            v = P + shift - c
            v -= np.outer(v @ seg.e_ax, seg.e_ax)
            k = int(np.argmax(v @ pick))
            az = math.degrees(math.atan2(float(v[k] @ seg.e_lat), float(v[k] @ seg.e_ant)))
            if side in ('medial', 'lateral'):
                az = self._palpable_border(seg, mesh, az, -1 if side == 'medial' else 1, t, shift)
            out[i] = az
        ok = ~np.isnan(out)
        out = np.interp(np.arange(len(out)), np.flatnonzero(ok), out[ok])
        pad = np.pad(out, 4, mode='edge')
        return np.convolve(pad, np.ones(9) / 9, mode='valid')

    def bone_edge_azimuth(self, ref, t):
        if ref not in self._edges:
            self._edges[ref] = self._edge_line(ref)
        return float(np.interp(np.clip(t, 0, 1), self.EDGE_LEVELS, self._edges[ref]))

    def _arc_len(self, seg, t, a, b, n=16):
        P = np.array([seg.surface_by_azimuth(t, float(x)) for x in np.linspace(a, b, n)])
        return float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())

    def azimuth_at_surface_cun(self, seg_name, t, az0, cun, tol=0.05):
        seg = self.seg(seg_name)
        target = abs(cun) * seg.cun_m
        step = 1.0 if cun >= 0 else -1.0
        lo, hi = 0.0, 5.0
        for _ in range(40):
            if self._arc_len(seg, t, az0, az0 + step * hi) >= target or hi > 180:
                break
            hi *= 1.6
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            if self._arc_len(seg, t, az0, az0 + step * mid) < target:
                lo = mid
            else:
                hi = mid
            if hi - lo < tol:
                break
        return az0 + step * 0.5 * (lo + hi)

    # ---------------- anchor resolution ----------------
    def _t_for_z(self, seg_name, z):
        seg = self.seg(seg_name)
        dz = seg.p1[2] - seg.p0[2]
        if abs(dz) < 1e-9:
            return 0.5
        raw = (z - seg.p0[2]) / dz
        if raw < -1e-6 or raw > 1 + 1e-6:
            self.flags.append(f'axial_clamped({raw:.2f})')
        return float(np.clip(raw, 0, 1))

    def axial_t(self, a):
        seg = self.seg(a['seg'])
        if 't' in a:
            return float(a['t'])
        if 'ruler' in a:
            return self.ruler_t(a['ruler'], float(a.get('cun', 0)), a.get('from', 'p0'))
        if 'cun' in a:
            return seg.t_from_cun(float(a['cun']), a.get('from', 'p0'))
        if 'z_frac' in a:
            return self._t_for_z(a['seg'], float(a['z_frac']) * self.H)
        if 'z_from' in a:
            c = self.landmark_cun[a['z_from']] + float(a.get('z_cun', 0))
            return self._t_for_z(a['seg'], self.z_of_cun(c))
        if 'vertebra' in a:
            z = self.vertebra_z[a['vertebra']] + float(a.get('dz_cun', 0)) * self.cun_trunk_lat
            return self._t_for_z(a['seg'], z)
        raise ValueError(a)

    def scalp_arc(self, s_cun, lat=0.0, out=0.0):
        head = self.seg('head')
        p = self.sagittal_point(3.0 + s_cun)
        if lat:
            fwd = self.sagittal_point(3.0 + s_cun + 0.05) - self.sagittal_point(3.0 + s_cun - 0.05)
            side = self.scalp_side(p, fwd)
            p = self.walk_on_skin(p, side, lat * head.cun_lat)
        if out:
            p = p + self.surface_normal(p) * out * head.cun_lat
        return p

    def struct_base(self, a):
        region = a.get('region')
        ref = self.struct(a['struct'], region)
        along = float(a.get('along', 0.5))
        long_dir = self.frame(region)['long'] if region else None
        base = ref.at(along, long_dir)
        if a.get('toward'):
            other = self.struct(a['toward'], region)
            base = base + (other.at(along, long_dir) - base) * float(a.get('frac', 0.4))
        for name, mm in (a.get('shift') or {}).items():
            base = base + self.direction(name, region) * (float(mm) / 1000)
        if region:
            base = base + self.soft_shift(region)
        return base

    def resolve(self, a):
        if 'struct' in a:
            base = self.struct_base(a)
            d = self.direction(a.get('dir', 'dorsal'), a.get('region'))
            return self.cast_side(base, d, a.get('region')), base
        if 'arc_cun' in a:
            p = self.scalp_arc(float(a['arc_cun']), float(a.get('lat', 0)), float(a.get('out', 0)))
            return p, self.seg('head').interior(self._t_for_z('head', p[2]))
        seg = self.seg(a['seg'])
        t = self.axial_t(a)
        out = float(a.get('out', 0))
        base = seg.interior(t)
        if 'edge' in a:
            az0 = self.bone_edge_azimuth(a['edge'], t)
            az = self.azimuth_at_surface_cun(a['seg'], t, az0, float(a.get('edge_cun', 0)))
            return seg.surface_by_azimuth(t, az, out), base
        if 'az' in a:
            return seg.surface_by_azimuth(t, float(a['az']), out), base
        if 'lat' in a:
            return seg.surface_by_lateral(t, float(a['lat']), a.get('face', 'anterior'), out), base
        raise ValueError(a)

    def project_midline(self, p):
        """Move a point onto the x=0 line of the skin: ray from the body core,
        inside the sagittal plane, through (0, y, z)."""
        q = np.array([0.0, p[1], p[2]])
        z = p[2]
        o = np.array([0.0, self.core_y(z), z])
        d = q - o
        if np.linalg.norm(d) < 0.01:
            d = np.array([0.0, np.sign(p[1] - o[1]) or 1.0, 0.0])
        d = unit(d)
        h = self.hits(o, d, MAX_RAY)
        if len(h):
            # the hit closest to the original point (skips inner sheets)
            return h[int(np.argmin(np.linalg.norm(h - q, axis=1)))]
        self.flags.append('midline_ray_miss')
        c, _, _ = self.closest(q)
        return c
