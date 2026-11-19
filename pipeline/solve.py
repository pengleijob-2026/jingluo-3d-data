"""Resolve all acupoints on the BodyVista skin.
usage: python -I solve.py <project_root> <workdir>
"""
import sys, json, time, traceback
import site; sys.path.append(site.getusersitepackages())
sys.stdout.reconfigure(encoding='utf-8')
ROOT, WORK = sys.argv[1], sys.argv[2]
sys.path.insert(0, WORK)
import numpy as np
import common
common.init(ROOT)
from body import Body, unit

A = ROOT + '/vendor/acupuncture-3d/data/'
pts = json.load(open(A + 'points.json', encoding='utf-8'))['points']
CH = {'SJ': 'TE', 'REN': 'CV', 'DU': 'GV'}


def gb_code(code):
    ch, n = code.split('-')
    return CH.get(ch, ch) + n


t0 = time.time()
bodies = {'R': Body(ROOT, WORK, 'R'), 'L': Body(ROOT, WORK, 'L')}
print('bodies built', round(time.time() - t0, 1), 's')
for s, b in bodies.items():
    print(s, 'hand soft shift', np.round(b._soft['hand'] * 1000, 1), 'mm')

def mandible_override(b, code):
    """ST5/ST6 from the real mandible mesh (GB text), replacing the ambiguous
    'Mandible along t' anchors whose PCA axis runs left-right."""
    V = b.struct('Mandible').mesh.vertices
    V = V[V[:, 0] > 0.025]
    gonion = V[np.argmax(-V[:, 1] - V[:, 2])]
    if code == 'ST6':      # 下颌角前上方一横指（中指）
        base = gonion + unit([0, 1, 1]) * 0.018 - np.array([0.02, 0, 0])
        b.flags.append('override_mandible')
        return b.cast(base, np.array([1.0, 0, 0])), base, gonion
    # ST5: 下颌角前方，咬肌附着部前缘 - lower border ~1.3 cun anterior to the angle
    ty = gonion[1] + 0.025
    S = V[np.abs(V[:, 1] - ty) < 0.003]
    low = S[np.argmin(S[:, 2])]
    base = low - np.array([0.015, 0, 0]) + np.array([0, 0, 0.004])
    b.flags.append('override_mandible')
    return b.cast(base, unit([1.0, 0, -0.35])), base, gonion


res = {}
for p in pts:
    code = gb_code(p['code'])
    midline = p['side'] == 'midline' or p['channel'] in ('REN', 'DU')
    entry = {'furia_code': p['code'], 'channel': CH.get(p['channel'], p['channel']),
             'furia_accuracy': p['accuracy'], 'anchor': p['anchor'], 'midline': midline,
             'note_pl': p.get('note_pl', ''), 'sides': {}}
    for side in (['R'] if midline else ['R', 'L']):
        b = bodies[side]
        b.flags = []
        try:
            if code in ('ST5', 'ST6'):
                q, base, gon = mandible_override(b, code)
                print('  gonion', side, np.round(gon, 4))
            else:
                q, base = b.resolve(p['anchor'])
            q = np.asarray(q, float)
            if float(p['anchor'].get('out', 0) or 0) != 0:
                q0 = q.copy()
                q = b.closest(q)[0]
                b.flags.append(f'out_projected_mm={np.linalg.norm(q - q0) * 1000:.1f}')
            if midline:
                q0 = q.copy()
                q = b.project_midline(q)
                b.flags.append(f'midline_shift_mm={np.linalg.norm(q - q0) * 1000:.1f}')
            n = b.surface_normal(q)
            _, dist, _ = b.closest(q)
            entry['sides'][side] = {'p': q.tolist(), 'n': n.tolist(), 'base': np.asarray(base).tolist(),
                                    'skin_dist_mm': dist * 1000, 'flags': list(b.flags)}
        except Exception as e:
            traceback.print_exc()
            entry['sides'][side] = {'error': f'{type(e).__name__}: {e}'}
    res[code] = entry
    print(code, {s: (np.round(np.array(v['p']), 4).tolist() if 'p' in v else v['error'])
                 for s, v in entry['sides'].items()},
          {s: v.get('flags') for s, v in entry['sides'].items() if v.get('flags')})

# ---- post-pass: points re-derived from the GB/T 12346-2021 wording where the
#      source anchor contradicts it (found on the check renders)
def setp(code, side, b, q, base, note):
    q = np.asarray(q, float)
    old = res[code]['sides'][side].get('p')
    res[code]['sides'][side] = {'p': q.tolist(), 'n': b.surface_normal(q).tolist(),
                                'base': np.asarray(base).tolist(), 'skin_dist_mm': b.closest(q)[1] * 1000,
                                'flags': list(b.flags) + ['override_gb_text'], 'override_note': note,
                                'moved_mm': None if old is None else float(np.linalg.norm(q - np.array(old)) * 1000)}
    print('override', code, side, np.round(q, 4), 'moved mm', res[code]['sides'][side]['moved_mm'])


def surface_curve(b, a, c, n=60):
    P = np.linspace(a, c, n)
    for _ in range(8):
        P[1:-1] = 0.5 * P[1:-1] + 0.25 * (P[:-2] + P[2:])
        P[1:-1] = np.array([b.closest(x)[0] for x in P[1:-1]])
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    return P, s


def at_frac(P, s, f):
    t = f * s[-1]
    return np.array([np.interp(t, s, P[:, k]) for k in range(3)])


for side in ('R', 'L'):
    b = bodies[side]
    head = b.seg('head')
    pupil_x = b.struct('Anterior segment of eyeball').centroid[0]
    # ST8 头维: 额角发际直上0.5寸，头正中线旁开4.5寸 -> same scalp ruler as GB13/GB15
    b.flags = []
    q = b.scalp_arc(0.5, 4.5)
    setp('ST8', side, b, q, q, '按国标“额角发际直上0.5寸，旁开4.5寸”用与GB13/GB15相同的头皮尺（沿皮肤量4.5寸）重算')
    # GB4/5/6: ST8 -> GB7 arc at 1/4, 1/2, 3/4
    st8 = np.array(res['ST8']['sides'][side]['p'])
    gb7 = np.array(res['GB7']['sides'][side]['p'])
    P, s = surface_curve(b, st8, gb7)
    for code, f in (('GB4', 0.25), ('GB5', 0.5), ('GB6', 0.75)):
        b.flags = []
        q = at_frac(P, s, f)
        setp(code, side, b, q, q, f'按国标“头维(ST8)至曲鬓(GB7)弧形连线”的 {f:g} 等分点重算（连线取皮肤测地线，未模拟鬓发弧度）')
    # GB14 阳白: 眉上1寸，瞳孔直上. 1 cun = (前发际 - 眉) / 3
    b.flags = []
    bl2 = np.array(res['BL2']['sides'][side]['p'])
    z = bl2[2] + (b.levels['hairline_ant'] - bl2[2]) / 3.0
    base = np.array([pupil_x, b.core_y(z), z])
    q = b.cast(base, np.array([0, 1.0, 0]))
    setp('GB14', side, b, q, base, '按国标“眉上1寸，瞳孔直上”重算（原锚点在瞳孔线内侧约1 cm）')
    # ST3 巨髎: 横平鼻翼下缘，瞳孔直下.  鼻翼下缘 ~ subnasale height on the midline profile
    b.flags = []
    gv25 = np.array(res['GV25']['sides']['R']['p'])
    gv26 = np.array(res['GV26']['sides']['R']['p'])
    best = None
    for zz in np.linspace(gv26[2], gv25[2], 40):
        o = np.array([0.0, b.core_y(zz), zz])
        h = b.hits(o, np.array([0, 1.0, 0]), 0.3)
        if len(h) and (best is None or h[-1][1] < best[1]):
            best = h[-1]
    z = best[2]
    base = np.array([pupil_x, b.core_y(z), z])
    q = b.cast(base, np.array([0, 1.0, 0]))
    setp('ST3', side, b, q, base, '按国标“横平鼻翼下缘，瞳孔直下”重算（鼻翼下缘取正中矢状轮廓的鼻下点高度）')
    # ST7 下关: 颧弓下缘中央与下颌切迹之间凹陷中 -> above the mandibular notch
    b.flags = []
    M = b.struct('Mandible').mesh.vertices
    M = M[M[:, 0] > 0.03]
    gon = M[np.argmax(-M[:, 1] - M[:, 2])]
    ram = M[(M[:, 1] < gon[1] + 0.045) & (M[:, 2] > gon[2] + 0.015)]
    ys = np.linspace(ram[:, 1].min(), ram[:, 1].max(), 30)
    top = []
    for y0, y1 in zip(ys[:-1], ys[1:]):
        sl = ram[(ram[:, 1] >= y0) & (ram[:, 1] < y1)]
        top.append(sl[np.argmax(sl[:, 2])] if len(sl) else None)
    top = [t for t in top if t is not None]
    zt_ = np.array([t[2] for t in top])
    i_con = int(np.argmax(zt_[:len(zt_) // 2]))                       # condyle (posterior peak)
    i_cor = len(zt_) // 2 + int(np.argmax(zt_[len(zt_) // 2:]))       # coronoid (anterior peak)
    notch = top[i_con + int(np.argmin(zt_[i_con:i_cor + 1]))]
    Vb0 = np.vstack([b.struct('Temporal bone').mesh.vertices, b.struct('Zygomatic bone').mesh.vertices])
    sl0 = Vb0[(np.abs(Vb0[:, 1] - notch[1]) < 0.004) & (Vb0[:, 2] > notch[2]) & (Vb0[:, 2] < notch[2] + 0.035)]
    arch_low = sl0[sl0[:, 0] > sl0[:, 0].max() - 0.005][:, 2].min()
    base = np.array([notch[0] - 0.015, notch[1], 0.5 * (notch[2] + arch_low)])
    q = b.cast(base, np.array([1.0, 0, 0]))
    print('  mandibular notch', np.round(notch, 4), 'condyle', np.round(top[i_con], 4), 'coronoid', np.round(top[i_cor], 4))
    setp('ST7', side, b, q, base, '按国标“颧弓下缘中央与下颌切迹之间凹陷中”重算：在下颌骨网格上找下颌切迹（原锚点取颧骨体，偏前约2 cm）')
    # GB3 上关: 颧弓上缘中央凹陷中 (directly above ST7 下关)
    b.flags = []
    st7 = np.array(res['ST7']['sides'][side]['p'])
    Vb = np.vstack([b.struct('Temporal bone').mesh.vertices, b.struct('Zygomatic bone').mesh.vertices])
    sl = Vb[(np.abs(Vb[:, 1] - notch[1]) < 0.004) & (Vb[:, 2] > notch[2]) & (Vb[:, 2] < notch[2] + 0.035)]
    arch = sl[sl[:, 0] > sl[:, 0].max() - 0.005]
    zt = arch[:, 2].max()
    base = np.array([arch[:, 0].max() - 0.02, notch[1], zt + 0.004])
    q = b.cast(base, np.array([1.0, 0, 0]))
    print('  GB3 arch z', round(arch[:, 2].min(), 4), '-', round(zt, 4), 'ST7 z', round(st7[2], 4))
    setp('GB3', side, b, q, base, '按国标“颧弓上缘中央凹陷中”重算：取下关(ST7)正上方、颧弓骨网格上缘（原锚点高于眼裂）')
    # GB2 听会: 耳屏间切迹与下颌骨髁突之间 -> same vertical line as TE21 (上切迹) and SI19 (耳屏正中)
    b.flags = []
    te21 = np.array(res['TE21']['sides'][side]['p'])
    si19 = np.array(res['SI19']['sides'][side]['p'])
    guess = si19 + (si19 - te21)
    q = b.closest(guess)[0]
    b.flags.append('override_low')
    setp('GB2', side, b, q, guess, '按国标“耳屏间切迹与下颌骨髁突之间”重算：与耳门(TE21)、听宫(SI19)排成竖线、等距向下（模型无耳廓，耳屏位置取自 Z-Anatomy 元数据）')
    # GB40 丘墟: 外踝的前下方 (source shifts it 'proximal' = behind the malleolus)
    b.flags = []
    fr = b.frame('foot')
    fib = b.struct('Fibula', 'foot').at(1.0, np.array([0, 0, -1.0]))
    base = fib + fr['distal'] * 0.012 + fr['plantar'] * 0.006 - np.array([0.012, 0, 0])
    q = b.cast(base, np.array([1.0, 0, 0]))
    setp('GB40', side, b, q, base, '按国标“外踝的前下方”重算：外踝尖向前 12 mm、向下 6 mm 处向外找皮肤（原锚点向后偏到外踝后方）')
    # KI3 太溪 内踝尖与跟腱之间凹陷中; KI5 水泉 太溪直下1寸; KI6 照海 内踝尖下1寸
    shank = b.seg('shank')
    tib = b.struct('Tibia', 'foot').mesh.vertices
    tip_z = tib[:, 2].min()
    low = tib[tib[:, 2] < tip_z + 0.03]
    prom = low[np.argmin(low[:, 0])]                    # most medial point = 内踝尖

    def interior_at(z):
        t = float(np.clip((z - shank.p0[2]) / (shank.p1[2] - shank.p0[2]), 0, 1))
        o = shank.interior(t)
        return np.array([o[0], o[1], z])

    def cast_toward(target_pt):
        o = interior_at(target_pt[2])
        d = np.asarray(target_pt, float) - o
        d[2] = 0
        return b.cast(o, unit(d)), o
    b.flags = []
    M_skin, _ = cast_toward(prom)
    A_skin = b.cast(interior_at(prom[2]), np.array([0, -1.0, 0]))   # Achilles (posterior) skin
    q3, base3 = cast_toward(0.5 * (M_skin + A_skin))
    setp('KI3', side, b, q3, base3, '按国标“内踝尖与跟腱之间的凹陷中”重算：内踝最突点与跟腱后缘皮肤的中点（原锚点落在内踝骨面上）')
    cun_sh = shank.cun_m
    b.flags = []
    q5, base5 = cast_toward(np.asarray(q3) - np.array([0, 0, cun_sh]))
    setp('KI5', side, b, q5, base5, '按国标“太溪直下1寸”重算（小腿骨度 1 寸 = 膝至外踝 / 16）')
    b.flags = []
    q6, base6 = cast_toward(M_skin - np.array([0, 0, cun_sh]))
    setp('KI6', side, b, q6, base6, '按国标“内踝尖下1寸”重算（内踝尖取胫骨远端最内侧突点）')
    # LI16 巨骨: 锁骨肩峰端与肩胛冈之间凹陷中 (source: 4.2 cun from midline = next to ST12)
    b.flags = []
    clav = b.bone_mesh('clavicle').vertices
    cl_end = clav[np.argmax(clav[:, 0])]
    base = cl_end + np.array([-0.012, -0.018, -0.015])
    q = b.cast(base, np.array([0, 0, 1.0]))
    setp('LI16', side, b, q, base, '按国标“锁骨肩峰端与肩胛冈之间凹陷中”重算：锁骨外端稍内后方向上找肩上皮肤（原锚点距正中线仅4.2寸，落在锁骨上窝附近）')
    # GB21 肩井: 第7颈椎棘突与肩峰最外侧点连线的中点
    b.flags = []
    c7 = np.array(b.fit['vertebra']['C7']['spinous'])
    acr = b.landmarks['acromion']
    mid = 0.5 * (c7 + acr)
    base = mid - np.array([0, 0, 0.03])
    q = b.cast(base, np.array([0, 0, 1.0]))
    setp('GB21', side, b, q, base, '按国标“第7颈椎棘突与肩峰最外侧点连线的中点”重算，向上投到肩上皮肤（原锚点偏外约2.6 cm）')
    # HT1 极泉: 腋窝中央 -> apex of the armpit gap between arm and chest wall
    b.flags = []
    best = None
    yh = b.landmarks['humeral_head'][1]
    for y in np.linspace(yh - 0.025, yh + 0.035, 13):
        for x in np.linspace(0.10, 0.20, 26):
            o = np.array([x, y, 1.17])
            h = b.hits(o, np.array([0, 0, 1.0]), 0.3)
            if not len(h):
                continue
            n_out = b.surface_normal(h[0])
            if n_out[2] > -0.3:          # hit from inside the body (or a side wall)
                continue
            if h[0][2] > b.landmarks['acromion'][2] - 0.03:
                continue
            if best is None or h[0][2] > best[2]:
                best = h[0]
    if best is not None and best[2] < b.landmarks['acromion'][2] - 0.03:
        b.flags.append('override_low')
        setp('HT1', side, b, best, best, '按国标“腋窝中央”重算：取上臂与胸壁之间可见腋窝缝隙的顶点。该皮肤模型上臂与胸壁在腋下互相穿插、腋窝被封住，此点比真实腋窝顶低约 5–10 cm，仅示意')
    else:
        print('HT1 axilla search failed', side, best)

# ---- Yintang (GV24+): midpoint between the two BL2 (inner ends of eyebrows)
b = bodies['R']
b.flags = []
r = np.array(res['BL2']['sides']['R']['p'])
l = np.array(res['BL2']['sides']['L']['p']) * [-1, 1, 1]
mid = (r + l) / 2
q = b.project_midline(mid)
g = b.landmarks['glabella']
gs = b.project_midline(g + [0, 0.02, 0])
res['GV24+'] = {'furia_code': None, 'channel': 'GV', 'furia_accuracy': 'added', 'midline': True,
                'anchor': {'rule': 'midpoint of BL2 L/R projected to sagittal skin line'},
                'note_pl': '', 'sides': {'R': {'p': q.tolist(), 'n': b.surface_normal(q).tolist(),
                                               'base': mid.tolist(), 'skin_dist_mm': 0.0,
                                               'flags': list(b.flags),
                                               'check_glabella_skin_mm': float(np.linalg.norm(q - gs) * 1000)}}}
print('GV24+', q.round(4), 'vs glabella skin', gs.round(4), 'dist mm', np.linalg.norm(q - gs) * 1000)

# edge fit diagnostics
for s, b in bodies.items():
    print(s, {k: (np.round(v * 1000, 1).tolist() if isinstance(v, np.ndarray) and v.shape == (3,) else v)
              for k, v in b._edges.items() if k.startswith('shift_') or k.endswith('_rms')})

json.dump(res, open(WORK + '/solved_furia.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print('done', round(time.time() - t0, 1), 's')
