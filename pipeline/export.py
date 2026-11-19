"""Convert solved points (Furia frame) to BodyVista frame and write acupoints.json.
usage: python -I export.py <project_root> <workdir>
"""
import sys, json
import site; sys.path.append(site.getusersitepackages())
sys.stdout.reconfigure(encoding='utf-8')
ROOT, WORK = sys.argv[1], sys.argv[2]
sys.path.insert(0, WORK)
import numpy as np
import trimesh
import opencc
import common
common.init(ROOT)
from body import Transform

cc = opencc.OpenCC('t2s')
T = Transform(WORK + '/transform_struct.json')
res = json.load(open(WORK + '/solved_furia.json', encoding='utf-8'))
tcm = json.load(open(ROOT + '/data/acupoint-names.json', encoding='utf-8'))['byMeridian']
gb = {}
order = []
for mer, lst in tcm.items():
    for p in lst:
        gb[p['code']] = dict(p, meridian=mer)
        order.append(p['code'])
MER_ZH = {'LU': '手太阴肺经', 'LI': '手阳明大肠经', 'ST': '足阳明胃经', 'SP': '足太阴脾经',
          'HT': '手少阴心经', 'SI': '手太阳小肠经', 'BL': '足太阳膀胱经', 'KI': '足少阴肾经',
          'PC': '手厥阴心包经', 'TE': '手少阳三焦经', 'GB': '足少阳胆经', 'LR': '足厥阴肝经',
          'GV': '督脉', 'CV': '任脉'}

skin = common.skin_outer()
missing_gb = [c for c in order if c not in res]
extra = [c for c in res if c not in gb]
print('GB codes without solution:', missing_gb, ' solved codes not in GB list:', extra)

LIP_INSIDE = {'GV28'}
EAR = {'SI19', 'TE17', 'TE18', 'TE19', 'TE20', 'TE21', 'GB2', 'GB8', 'GB9', 'GB10'}

out = []
for code in order:
    g = gb[code]
    e = res.get(code)
    sides = ['C'] if (e is None or e['midline']) else ['L', 'R']
    for side in sides:
        rec = {'code': code, 'name': cc.convert(g['name']), 'pinyin': g.get('pinyin', ''),
               'meridian': g['meridian'], 'meridianName': MER_ZH[g['meridian']], 'side': side}
        sv = None if e is None else e['sides'].get('R' if side in ('C', 'R') else 'L')
        if sv is None or 'p' not in sv:
            rec.update({'xyz': None, 'normal': None, 'accuracy': 'failed',
                        'accuracyNote': '无锚点数据' if e is None else ('计算出错: ' + sv.get('error', '?'))})
            out.append(rec)
            continue
        p = np.array(sv['p'], float)
        n = np.array(sv['n'], float)
        if side == 'L':                 # computed on the mirrored body
            p = p * [-1, 1, 1]
            n = n * [-1, 1, 1]
        P = T.to_bv(p)[0]
        N = T.dir_to_bv(n)[0]
        N = N / np.linalg.norm(N)
        c, d, tri = trimesh.proximity.closest_point(skin, [P])
        dist_mm = float(d[0]) / T.s * 1000
        notes = []
        level = 'high'
        flags = sv.get('flags', [])
        fa = e['furia_accuracy']
        if 'override_gb_text' in flags:
            fa = 'override'
            level = 'medium'
            notes.append(sv['override_note'])
            if 'override_low' in flags:
                level = 'low'
        if fa == 'structural':
            level = 'medium'; notes.append('按骨性结构定位（Z-Anatomy 骨骼+偏移毫米），非寸法')
        if 'edge' in e['anchor'] and fa != 'override':
            level = 'medium'; notes.append('按骨缘+旁开寸定位，胫腓骨与皮肤源模型错位，已按骨面下软组织厚度把胫腓骨整体平移约 42 mm 校正')
        if 'override_mandible' in flags:
            notes.append('原锚点“下颌骨PCA轴”方向有歧义，改为按国标条文在下颌骨网格上找下颌角重新定位')
        if fa == 'approx':
            level = 'low'; notes.append('原数据即为近似（模型无腋窝凹陷，点放在胸侧壁）')
        if any(f.startswith('axial_clamped') for f in flags):
            if level == 'high':
                level = 'medium'
            notes.append('轴向位置超出体段范围被截断到端点')
        if 'ray_miss_closest_point' in flags:
            dmiss = float(np.linalg.norm(np.array(sv['p']) - np.array(sv['base'])) * 1000)
            if dmiss > 12 or 'struct_from_metadata' in flags:
                level = 'low'
            elif level == 'high':
                level = 'medium'
            notes.append(f'射线未命中皮肤（骨性基点在皮肤外或贴近皮肤），改用最近皮肤点，基点到落点 {dmiss:.0f} mm')
        if 'struct_from_metadata' in flags:
            level = 'low'; notes.append('BodyVista 无耳廓网格，耳廓位置取自 Z-Anatomy 元数据，点落在头侧皮肤上')
        for f in flags:
            if f.startswith('out_projected_mm='):
                notes.append(f'原锚点要求抬离皮肤（out），已投影回皮肤（移动 {float(f.split("=")[1]):.1f} mm）')
        if code in ('CV1', 'GV1'):
            if level == 'high':
                level = 'medium'
            notes.append('会阴区皮肤网格很粗，CV1 可能比“阴囊根部与肛门连线中点”偏前约 1–2 cm，仅示意')
        if code in LIP_INSIDE:
            level = 'low'; notes.append('国标位置在上唇内（唇系带与牙龈交点），模型无口腔内皮肤，放在上唇外表面正中')
        if code == 'GV24+':
            level = 'medium'; notes.append('新增：取两侧攒竹(BL2)连线中点投影到正中线皮肤')
        if e['midline'] and code not in ('GV24+',):
            ms = [float(f.split('=')[1]) for f in flags if f.startswith('midline_shift')]
            if ms and ms[0] > 3:
                notes.append(f'原锚点偏离正中线，已投影到正中线（移动 {ms[0]:.1f} mm）')
        if side in ('L', 'R') and 'L' in e['sides'] and 'p' in e['sides']['L']:
            asym = float(np.linalg.norm(np.array(e['sides']['R']['p']) - np.array(e['sides']['L']['p'])) * 1000)
            rec['asymmetryMm'] = round(asym, 1)
            if asym > 8:
                if level == 'high':
                    level = 'medium'
                notes.append(f'左右镜像差 {asym:.1f} mm（皮肤网格左右不完全对称）')
        if dist_mm > 2:
            level = 'low'; notes.append(f'离皮肤 {dist_mm:.1f} mm')
        rec.update({'xyz': [round(float(x), 5) for x in P], 'normal': [round(float(x), 4) for x in N],
                    'accuracy': level, 'accuracyNote': '；'.join(notes) if notes else '寸法锚点，射线命中皮肤',
                    'skinDistMm': round(dist_mm, 2)})
        out.append(rec)

meta = {
    'description': '362 个国标经穴（GB/T 12346-2021）在 BodyVista full_body.glb 皮肤上的坐标',
    'coordinateSystem': 'BodyVista 世界坐标：Z 朝上，正面朝 -Y，人体左侧在 +X；单位与 GLB 相同（身高约 2.53）',
    'side': 'L/R = 人体左/右侧（左侧 x>0），C = 正中线',
    'normal': '皮肤外法线（单位向量）',
    'accuracy': 'high=寸法锚点且射线命中；medium=骨性结构/骨缘定位、轴向截断或左右差>8mm；low=近似/射线未命中/耳廓缺失/唇内；failed=未能计算',
    'transformFromFuriaRig': {'formula': 'xyz_bv = s * R * xyz_furia + t', 's': T.s,
                              'R': np.round(T.R, 6).tolist(), 't': np.round(T.t, 6).tolist()},
    'sources': ['GB/T 12346-2021 穴名与编码（data/acupoint-names.json，不含定位条文）',
                'acupuncture-3d (FuriaRozkwit) points.json 锚点（MIT）、rig_fitted.json / structures.json / meridians.json（CC BY-SA 4.0）',
                'BodyVista full_body.glb：骨骼等结构来自 Z-Anatomy（CC BY-SA 4.0），皮肤来自 BodyParts3D（CC BY 4.0）'],
    'license': 'CC BY-SA 4.0（继承 Z-Anatomy 网格与 acupuncture-3d 测量数据的许可）',
    'count': len(out),
}
json.dump({'meta': meta, 'points': out}, open(ROOT + '/data/acupoints.json', 'w', encoding='utf-8'),
          ensure_ascii=False, indent=1)
import collections
print('records', len(out), collections.Counter(r['accuracy'] for r in out))
print('unique codes', len(set(r['code'] for r in out)))
