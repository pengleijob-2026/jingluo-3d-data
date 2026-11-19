"""Shared helpers (own code). Loads BodyVista GLB meshes by english name."""
import sys, json, re
import numpy as np
import trimesh

ROOT = None
_scene = None
_names = None

def init(root):
    global ROOT
    ROOT = root

def scene():
    global _scene
    if _scene is None:
        _scene = trimesh.load(ROOT + '/vendor/BodyVista/assets/full_body.glb', force='scene')
    return _scene

def names():
    global _names
    if _names is None:
        cat = json.load(open(ROOT + '/vendor/BodyVista/catalog.json', encoding='utf-8'))
        sysmap = {}
        for s in cat['systems']:
            for p in s['parts']:
                sysmap[p] = s['id']
        _names = {p['id']: (p['english'], p['name'], sysmap.get(p['id'])) for p in cat['parts']}
    return _names

def by_english(eng):
    for k, v in names().items():
        if v[0].lower() == eng.lower():
            return k
    return None

def mesh(pid):
    """World-space trimesh for a part id."""
    sc = scene()
    for node in sc.graph.nodes_geometry:
        T, g = sc.graph[node]
        if g == pid:
            m = sc.geometry[g].copy()
            m.apply_transform(T)
            return m
    raise KeyError(pid)

def mesh_eng(eng):
    pid = by_english(eng)
    if pid is None:
        raise KeyError(eng)
    return mesh(pid)

def skin():
    return mesh(next(k for k, v in names().items() if v[2] == 'skin'))


_outer = None

def skin_outer():
    """BodyVista skin is a double shell (outer + inner ~4.5 mm apart).
    Return the outer shell only (component with the larger bounding box)."""
    global _outer
    if _outer is None:
        sk = skin()
        comps = trimesh.graph.connected_components(sk.face_adjacency, nodes=np.arange(len(sk.faces)))
        comps = sorted(comps, key=len, reverse=True)[:2]
        best = max(comps, key=lambda c: np.prod(np.ptp(sk.triangles_center[c], axis=0)))
        m = sk.submesh([best], append=True, repair=False)
        _outer = m
    return _outer
