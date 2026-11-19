"""Off-screen check renders.  usage: python -I render.py <project_root> <workdir> [prefix]

Skin = outer shell, semi-transparent. Points/lines are drawn only where the skin
faces the camera (normal test), so the far side does not clutter the view.
clip='neg' keeps x<0 (person's right half), clip='pos' keeps x>0 (left half).
White dot inside a point = accuracy 'low'.
"""
import sys, json, os
import site; sys.path.append(site.getusersitepackages())
sys.stdout.reconfigure(encoding='utf-8')
ROOT, WORK = sys.argv[1], sys.argv[2]
PREFIX = sys.argv[3] if len(sys.argv) > 3 else ROOT + '/docs/check/check-'
sys.path.insert(0, WORK)
import numpy as np
import pyvista as pv
import trimesh
import common
common.init(ROOT)

COL = {'LU': '#e6e6e6', 'LI': '#9aa7b8', 'ST': '#f0a500', 'SP': '#c0782a', 'HT': '#ff3b30',
       'SI': '#b0302a', 'BL': '#2f8cff', 'KI': '#1c3f99', 'PC': '#ff5fc8', 'TE': '#9b4fd1',
       'GB': '#2ecc71', 'LR': '#0b7a43', 'GV': '#ffe14d', 'CV': '#ff8a3d'}
pts = json.load(open(ROOT + '/data/acupoints.json', encoding='utf-8'))['points']
paths_file = ROOT + '/data/meridians-paths.json'
paths = json.load(open(paths_file, encoding='utf-8'))['meridians'] if os.path.exists(paths_file) else []

sk = common.skin_outer()
skin = pv.PolyData(np.asarray(sk.vertices), np.hstack([np.full((len(sk.faces), 1), 3), sk.faces]).ravel())


def normals_at(X):
    c, d, tri = trimesh.proximity.closest_point(sk, np.asarray(X))
    return sk.face_normals[tri]


def keep_x(P, clip):
    P = np.asarray(P)
    if clip == 'neg':
        return P[:, 0] < 0.002
    if clip == 'pos':
        return P[:, 0] > -0.002
    return np.ones(len(P), bool)


def scene(cam_pos, focal, scale, fname, title, labels=False, only=None, pt_size=9, clip=None,
          size=(1400, 1800), up=(0, 0, 1)):
    pl = pv.Plotter(off_screen=True, window_size=size)
    pl.set_background('#14161a')
    sm = skin
    if clip == 'neg':
        sm = skin.clip(normal='x', origin=(0.002, 0, 0))
    elif clip == 'pos':
        sm = skin.clip(normal='-x', origin=(-0.002, 0, 0))
    pl.add_mesh(sm, color='#d9b8a6', opacity=0.45, smooth_shading=True, specular=0.15)
    view = np.asarray(focal, float) - np.asarray(cam_pos, float)
    view /= np.linalg.norm(view)
    for m in paths:
        if only and m['meridian'] not in only:
            continue
        for stroke in m['strokes']:
            P = np.array(stroke['points'])
            vis = ((normals_at(P) @ -view) > -0.05) & keep_x(P, clip)
            runs, cur = [], []
            for q, ok in zip(P, vis):
                if ok:
                    cur.append(q)
                elif cur:
                    runs.append(cur)
                    cur = []
            if cur:
                runs.append(cur)
            for r in runs:
                if len(r) >= 2:
                    pl.add_mesh(pv.lines_from_points(np.array(r)), color=COL[m['meridian']], line_width=3.5)
    sel = [p for p in pts if p['xyz'] is not None and (not only or p['meridian'] in only)
           and np.dot(p['normal'], -view) >= -0.05 and keep_x([p['xyz']], clip)[0]]
    by = {}
    for p in sel:
        by.setdefault(p['meridian'], []).append(p)
    for mer, lst in by.items():
        pl.add_points(np.array([p['xyz'] for p in lst]), color=COL[mer], point_size=pt_size,
                      render_points_as_spheres=True)
        low = [p['xyz'] for p in lst if p['accuracy'] == 'low']
        if low:
            pl.add_points(np.array(low), color='white', point_size=pt_size * 0.45, render_points_as_spheres=True)
    if labels and sel:
        pl.add_point_labels(np.array([p['xyz'] for p in sel]), [p['code'] for p in sel], font_size=11,
                            point_size=1, shape_opacity=0.0, text_color='white', always_visible=True,
                            show_points=False)
    pl.camera.position = tuple(cam_pos)
    pl.camera.focal_point = tuple(focal)
    pl.camera.up = up
    pl.camera.parallel_projection = True
    pl.camera.parallel_scale = scale
    pl.add_text(title, font_size=12, color='white')
    entries = [[k, COL[k]] for k in COL if not only or k in only]
    pl.add_legend(entries, bcolor='#14161a', size=(0.12, 0.3), loc='lower right', face='circle')
    pl.screenshot(fname)
    pl.close()
    print('saved', fname)


ARM = ['LU', 'HT', 'PC', 'LI', 'SI', 'TE']
LEG = ['ST', 'SP', 'KI', 'LR', 'GB', 'BL']
scene((0, -6, 0), (0, 0, 0), 1.33, PREFIX + 'front.png', 'Front (camera -Y). 14 meridians')
scene((0, 6, 0), (0, 0, 0), 1.33, PREFIX + 'back.png', 'Back (camera +Y)')
scene((-6, 0, 0), (0, 0, 0), 1.33, PREFIX + 'side-right.png', 'Right side (camera -X), left half clipped', clip='neg')
scene((6, 0, 0), (0, 0, 0), 1.33, PREFIX + 'side-left.png', 'Left side (camera +X), right half clipped', clip='pos')
hc = (0, 0.0, 1.08)
scene((0, -6, 1.08), hc, 0.2, PREFIX + 'head-front.png', 'Head / face, front', labels=True, pt_size=12, size=(1600, 1600))
scene((-6, 0, 1.08), (0, 0, 1.08), 0.2, PREFIX + 'head-side.png', 'Head, right side', labels=True, pt_size=12,
      clip='neg', size=(1600, 1600))
scene((0, 6, 0.98), (0, 0, 0.98), 0.3, PREFIX + 'head-back.png', 'Head / neck / upper back', labels=True,
      pt_size=12, size=(1600, 1600))
scene((-0.35, -6, 0.3), (-0.35, 0, 0.3), 0.55, PREFIX + 'arm-front.png', 'Right arm, palmar side (hand 3 yin)',
      labels=True, pt_size=11, only=ARM, size=(1600, 1600))
scene((-0.35, 6, 0.3), (-0.35, 0, 0.3), 0.55, PREFIX + 'arm-back.png', 'Right arm, dorsal side (hand 3 yang)',
      labels=True, pt_size=11, only=ARM, size=(1600, 1600))
scene((-0.40, -6, -0.12), (-0.40, 0, -0.12), 0.15, PREFIX + 'hand-palmar.png', 'Right hand, palmar', labels=True,
      pt_size=12, only=ARM, size=(1400, 1400))
scene((-0.40, 6, -0.12), (-0.40, 0, -0.12), 0.15, PREFIX + 'hand-dorsal.png', 'Right hand, dorsal', labels=True,
      pt_size=12, only=ARM, size=(1400, 1400))
scene((-6, 0, -0.75), (-0.15, 0, -0.75), 0.6, PREFIX + 'leg-lateral.png', 'Right leg, lateral', labels=True,
      pt_size=11, only=LEG, clip='neg', size=(1600, 1600))
scene((6, 0, -0.75), (-0.15, 0, -0.75), 0.6, PREFIX + 'leg-medial.png', 'Right leg, medial (left half clipped)',
      labels=True, pt_size=11, only=LEG, clip='neg', size=(1600, 1600))
scene((-6, 0, -1.17), (-0.13, 0, -1.17), 0.14, PREFIX + 'foot-lateral.png', 'Right foot, lateral', labels=True,
      pt_size=12, only=LEG, clip='neg', size=(1400, 1400))
scene((6, 0, -1.17), (-0.13, 0, -1.17), 0.14, PREFIX + 'foot-medial.png', 'Right foot, medial', labels=True,
      pt_size=12, only=LEG, clip='neg', size=(1400, 1400))
