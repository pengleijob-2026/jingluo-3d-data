"""Download / verify the upstream files this project is computed from, pinned to exact commits.

usage: python -I pipeline/fetch_upstream.py <project_root>

Files go to <project_root>/vendor/ (ignored by git). Nothing from these repositories
is redistributed here; see THIRD_PARTY_NOTICES.md for their licenses.

BodyVista files are downloaded automatically. The acupuncture-3d repository
(github.com/FuriaRozkwit/acupuncture-3d) was no longer reachable as of 2026-11-19,
so its four files cannot be downloaded: if you have your own copy of commit
7c4fd3feb1b754ef2a90b67aaabf3498114dc932, put the files under
vendor/acupuncture-3d/data/ and this script will check their sha256.
Without them, steps 1-2 (fit_transform.py, solve.py) cannot be rerun; steps 3-5
(export.py, paths.py, render.py) only need BodyVista plus the shipped
pipeline/transform_struct.json, pipeline/solved_furia.json and data/*.json.
"""
import hashlib
import os
import sys
import urllib.request

BODYVISTA = 'https://raw.githubusercontent.com/FreedomIntelligence/BodyVista/82b2ea12f71fd123eb728f465eb413bb8d6a150e/'
FURIA = 'https://raw.githubusercontent.com/FuriaRozkwit/acupuncture-3d/7c4fd3feb1b754ef2a90b67aaabf3498114dc932/'

FILES = [
    # (local path under vendor/, url, sha256)
    ('BodyVista/assets/full_body.glb', BODYVISTA + 'assets/full_body.glb',
     '4d196aab81e263a113ca9ce221feda4e15e1b623d1730505dd35497e1fe20c37'),
    ('BodyVista/catalog.json', BODYVISTA + 'catalog.json',
     '26dee6c277e5ef771e44feed6441b510f90b83976831888b28bc9dfa5107a0e6'),
    # acupuncture-3d: upstream offline, url kept for the record; checked only if present locally
    ('acupuncture-3d/data/points.json', FURIA + 'data/points.json',
     '4adc3d99e02011d64619051c829464457f34c181f6e9812755723fad1e4d3e27'),
    ('acupuncture-3d/data/meridians.json', FURIA + 'data/meridians.json',
     'cf0333eebe2f4dfc0257eb1a16b36325c97c5eada8e5075e5f132bf9e3a9cf68'),
    ('acupuncture-3d/data/rig_fitted.json', FURIA + 'data/rig_fitted.json',
     'c97d298cea42478d6b1bbbe7c8363c3dd50f03b70990b7a20f0f5d8b893a88b7'),
    ('acupuncture-3d/data/structures.json', FURIA + 'data/structures.json',
     'e5d7ee88b2fe8d2804a457631778261d958eea3cb59a5da29166e6c85b2d5acf'),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    vendor = os.path.join(sys.argv[1], 'vendor')
    bad = 0
    for rel, url, want in FILES:
        dst = os.path.join(vendor, rel)
        if os.path.exists(dst) and sha256(dst) == want:
            print('ok      ', rel)
            continue
        if url.startswith(FURIA):
            print('missing ', rel, '(acupuncture-3d is offline; only needed to rerun steps 1-2)')
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        print('download', rel)
        tmp = dst + '.part'
        urllib.request.urlretrieve(url, tmp)
        got = sha256(tmp)
        if got != want:
            os.remove(tmp)
            print('  sha256 mismatch:', got)
            bad += 1
            continue
        os.replace(tmp, dst)
    if bad:
        sys.exit(f'{bad} file(s) failed verification')
    print('done')


if __name__ == '__main__':
    main()
