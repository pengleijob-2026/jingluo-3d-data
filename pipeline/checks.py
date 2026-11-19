import sys, json
import site; sys.path.append(site.getusersitepackages())
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
r = json.load(open(sys.argv[1], encoding='utf-8'))
P = lambda c, s='R': np.array(r[c]['sides'][s]['p'])
for c in ['CV1', 'CV2', 'GV1', 'GV2', 'CV8', 'GV20', 'GV24+', 'HT1']:
    print(c, np.round(P(c), 4))
ct = 0.02069   # trunk lateral cun (m)
cv = 0.02485   # trunk vertical cun
print('BL23 lateral offset from GV4 (cun, want 1.5):', round(abs(P('BL23')[0]) / ct, 2))
print('BL52 lateral (cun, want 3):', round(abs(P('BL52')[0]) / ct, 2))
print('ST25 lateral (cun, want 2):', round(abs(P('ST25')[0]) / ct, 2), ' KI16 (want 0.5):', round(abs(P('KI16')[0]) / ct, 2))
print('CV8 -> CV12 vertical (cun, want 4):', round((P('CV12')[2] - P('CV8')[2]) / cv, 2))
sh = 0.373 / 16
print('ST35 -> ST36 along leg (cun, want 3):', round(np.linalg.norm(P('ST36') - P('ST35')) / sh, 2))
print('SP6 above KI3-level (cun, want ~3):', round((P('SP6')[2] - P('KI3')[2]) / sh, 2))
fa = 0.2446 / 12
print('PC7 -> PC6 (cun, want 2):', round(np.linalg.norm(P('PC6') - P('PC7')) / fa, 2), ' PC7->PC4 (want 5):', round(np.linalg.norm(P('PC4') - P('PC7')) / fa, 2))
print('LI11 -> LI10 (want 2):', round(np.linalg.norm(P('LI10') - P('LI11')) / fa, 2))
