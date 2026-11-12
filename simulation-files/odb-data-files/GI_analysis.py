"""
gi_analysis.py  --  CSA vs Garcia et al. (2024) and gyrification index,
                    outside Abaqus (numpy + matplotlib only).

Reads the ONE results file written by evol_area.py (brain_results.json or
.json.gz).  Legacy <job>_pial_disp.csv (+ <job>_csa.csv) files still work.

Step time is in DAYS (18-day step, t = 0 is E40, t = 18 is P16).  Runs from
the old 1-s step (t = 0 ... 1) are converted to days when they are loaded,
so old and new runs can be drawn in one figure.

Default output: a manuscript-style figure (like the draft's Fig. 8) --
columns Pushing | Pulling, rows CSA vs Garcia and GI vs ferrets, colour =
gamma / gamma_hat, marker = hypothesis -- plus a consistency table.  Runs are
assigned to a column by 'push' / 'pull' in the job name; 'control' runs appear
in both columns.

Usage
-----
    python3 GI_analysis.py Logis_brain_results.json --list        # list running cases
    python3 GI_analysis.py Logis_brain_results.json --gi-only     # exactly the Fig. 8 layout for GI
    python3 GI_analysis.py Logis_brain_results.json --csa-only
    python3 GI_analysis.py Logis_brain_results.json --overlay     # both hypotheses in one panel
    python3 GI_analysis.py Logis_brain_results.json --runs Const_push_Z1_g2p22,Const_pull_Z1_g1
    python3 GI_analysis.py Logis_brain_results.json --save figure8.pdf        # vector for the paper
    python3 GI_analysis.py Logis_brain_results.json --export gi_results.csv   # ONE csv, all runs
    python3 GI_analysis.py Logis_brain_results.json --gi-only --hull outer
    python3 GI_analysis.py Logis_brain_results.json --match Refactor      # one run set
    python3 GI_analysis.py Logis_brain_results.json --drops   # where/why GI falls
    python3 GI_analysis.py Logis_brain_results.json --zcheck  # does the z-slice matter?
    python3 GI_analysis.py Logis_brain_results.json --contours 10 18 --run Logis_pull_Z1_3  # plot deformed pial 


GI definition (Karan's GI script, unchanged; hull='karan'):
    pial nodes on one z-face, off the symmetry planes -> sorted by angle about
    their REFERENCE centroid -> deformed contour (L_oc) -> convex hull with
    interior local-x-minima removed (L_ch) -> GI = L_oc / L_ch.
hull='outer' uses only the outer envelope (hull of contour + origin).
"""

import fnmatch
import gzip
import json
import os
import re
import sys

import numpy as np
import matplotlib

_HEADLESS = sys.platform.startswith('linux') and not os.environ.get('DISPLAY')
if _HEADLESS:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- settings
DAYS_PER_SIM = 1.0      # step time is in days (18-day step)
T_END        = 18.0     # step time at the end of the run (P16)
OLD_STEP_DAYS = 18.0    # days in 1.0 s of the old 1-s step; runs that end at
                        # t ~ 1 are multiplied by this when they are loaded
DAY_OFFSET   = -2.0     # postnatal day at t = 0.  Draft Fig. 1: t = 0 is E40 (P-2),
                        # T1 = 2 ~ P0, T2 = 8 ~ P6, t = 18 is P16.
START_LABEL  = 'E40'    # tick label for t = 0 when it falls before birth
VUMAT_DAY_OFFSET = -2.0 # the offset integ_gray itself uses: t_d = totalTime - 2
T2_SIM       = 8.0      # step time at the end of prescribed growth (P6)
XY_MIN       = 0.01     # [mm] drop pial nodes on the symmetry planes
Z_SLICE      = None     # None -> lowest-z face, 'mid' -> node plane nearest the
                        # middle (what Karan's z_slice = 0.6 did on a multi-layer
                        # mesh), or a z value.  --zcheck compares every plane.
Z_TOL        = 1e-4
HULL         = 'mirror' # 'mirror' (full outline), 'karan' (draft) or 'outer' (+ origin);
                        # the mirrored hull cannot jump at the symmetry walls

# ---------------------------------------------------------------- inputs
DATA_DIR = os.environ.get('DATA_DIR', os.path.dirname(os.path.abspath(__file__)))
RESULTS  = 'brain_results_old_wofg33.json'      # or a list, e.g.
# RESULTS = ['brain_results.json', 'brain_results_old.json']   # both growth laws
RUNS     = None         # None = every run in the file, or e.g.
# RUNS = [('Const_push_Z1_g2p22', 'push',                    {'color': 'C0'}),
#         ('Const_pull_Z1_g1',    r'pull, $\hat\gamma$ = 1', {'color': 'C3'}),
#         ('Const_pull_Z1_g3',    r'pull, $\hat\gamma$ = 3', {'color': 'C1'})]
MATCH    = None         # keep only runs whose name contains this text, e.g. 'Refactor'
                        # or 'Old', or a pattern like '*_push_*Refactor*'; None = all
# Run sets drawn in one figure: (text in the run name, legend label, line style).
# Same gamma -> same colour; the set is told apart by the line style.  A run
# that matches none of these is drawn solid and not named in the legend.
RUN_SETS = [('Const_',   'logistic growth rate', '-'),
            ('Refactor', 'constant growth rate', (0, (5, 2))),
            ('Old',      'constant rate (pre-refactor)', (0, (1, 1.5)))]
MODE     = 'figure'     # 'figure' | 'drops' (contours before/after each GI drop)
                        # | 'overlay' (both hypotheses in one panel)

# Garcia et al. (2024) whole-brain logistic
G_A, G_RG, G_TMID = 8.28, 0.134, 20.5
# cortical surface-area law used in integ_gray (keep in step with the VUMAT)
SA_R, SA_TMID = 0.151, 19.7
# cortical thickness law used in integ_gray_thickness (keep in step with the VUMAT)
TH_H1, TH_HINF, TH_TAU, TH_FEL = 0.421, 0.859, 9.6, 1.0   # mm, mm, day, f_el

# ---------------------------------------------------------------- figure style
# Manuscript-style figure (like the draft's Fig. 8): one column per hypothesis,
# colour = gamma / gamma_hat (same value -> same colour in both columns),
# marker = hypothesis, ferret GI as squares.
CMAP       = 'plasma'
CMAP_RANGE = (0.76, 0.0)           # smallest gamma -> orange, largest -> navy
GAMMA_SHOW = None                  # e.g. (0.05, 0.2, 0.5, 1, 3); None = every run
# Adjacent colours of an 8-step ramp are hard to tell apart (and nearly
# identical for red-green colour-blind readers), so each GI curve also gets its
# gamma value printed at its right-hand end.  Set () to switch off.
END_LABEL_ROWS = ('gi',)
MARKERS    = {'push': 'o', 'pull': '^'}
# control runs: key contains 'control' and ends in _A / _B (Karan's naming).
# Grey + dash pattern keeps them apart from the gamma colour ramp.
CONTROL_STYLE = {'A': ('Control (GM)',      '0.35', (0, (1, 1.5))),
                 'B': ('Control (GZ + GM)', '0.35', (0, (6, 2, 1, 2)))}
# Ferret GI from histology.  READ OFF THE DRAFT FIGURE -- replace with the
# measured values before the figure goes into the paper.
FERRET_GI    = [(10, 1.036), (10, 1.095), (16, 1.255), (16, 1.341), (16, 1.402)]
FERRET_COLOR = '#1f8fff'
XLIM   = (DAY_OFFSET, DAY_OFFSET + DAYS_PER_SIM * T_END + 0.5)   # postnatal days
XTICKS = (0, 4, 8, 12, 16)
# Day the area curves are normalised at (model AND Garcia).  None = start of the
# simulation (E40).  8.0 compares growth only over the window Garcia et al.
# actually measured (their data start at P8); before that it is extrapolation.
CSA_REF_DAY  = None     # None = start of the simulation (E40)
GARCIA_FIRST_DAY = 8.0  # Garcia et al.'s first scan; the fit is drawn dotted before it
SHOW_T2      = True                # dashed line at the end of prescribed growth
FIG_RC = {'font.family': 'serif',
          'font.serif': ['Times New Roman', 'Times', 'Nimbus Roman', 'STIXGeneral',
                         'DejaVu Serif'],
          'mathtext.fontset': 'stix',
          'font.size': 11, 'axes.titlesize': 13, 'axes.labelsize': 12,
          'legend.fontsize': 11, 'xtick.labelsize': 11, 'ytick.labelsize': 11}


def _resolve(path):
    if os.path.isabs(path) or os.path.isfile(path):
        return path
    return os.path.join(DATA_DIR, path)


# ================================================================ geometry
def polyline_length(points):
    pts = np.asarray(points, float)
    if len(pts) < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())


def convex_hull_indices(points, ordered=False):
    """Hull vertex indices, ascending (== np.sort(ConvexHull(P).vertices)),
    or counter-clockwise around the hull if ordered=True."""
    P = np.asarray(points, float)
    if len(P) < 3:
        return np.arange(len(P))
    order = np.lexsort((P[:, 1], P[:, 0]))

    def cross(o, a, b):
        return ((P[a, 0] - P[o, 0]) * (P[b, 1] - P[o, 1])
                - (P[a, 1] - P[o, 1]) * (P[b, 0] - P[o, 0]))

    lower, upper = [], []
    for i in order:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], i) <= 0:
            lower.pop()
        lower.append(int(i))
    for i in order[::-1]:
        while len(upper) >= 2 and cross(upper[-2], upper[-1], i) <= 0:
            upper.pop()
        upper.append(int(i))
    if ordered:
        return np.array(lower[:-1] + upper[:-1], dtype=int)
    return np.array(sorted(set(lower[:-1] + upper[:-1])), dtype=int)


def gyrification_index(ref_xy, def_xy, hull=HULL):
    """(GI, L_oc, L_ch, contour, hull) for one frame."""
    ref_xy = np.asarray(ref_xy, float)
    def_xy = np.asarray(def_xy, float)
    c = ref_xy.mean(axis=0)
    order = np.argsort(np.arctan2(ref_xy[:, 1] - c[1], ref_xy[:, 0] - c[0]))
    contour = def_xy[order]
    if hull == 'mirror':
        # full cross-section: mirror the quarter across x = 0 and y = 0
        q = contour if contour[0, 1] <= contour[-1, 1] else contour[::-1]
        full = np.vstack([q, q[::-1] * [-1, 1], q * [-1, -1], q[::-1] * [1, -1]])
        hv = convex_hull_indices(full, ordered=True)
        H = full[np.r_[hv, hv[:1]]]
        L_oc = polyline_length(np.vstack([full, full[:1]])) / 4.0
        L_ch = polyline_length(H) / 4.0
        return L_oc / L_ch, L_oc, L_ch, q, H
    if hull == 'outer':
        P = np.vstack([contour, [[0.0, 0.0]]])
        hv = convex_hull_indices(P)
        H = contour[hv[hv < len(contour)]]
        H = H[np.argsort(np.arctan2(H[:, 1], H[:, 0]))]
    else:
        H = contour[convex_hull_indices(contour)]
        x = H[:, 0]
        keep = np.ones(len(x), dtype=bool)
        keep[1:-1] &= ~((x[1:-1] < x[:-2]) & (x[1:-1] < x[2:]))
        H = H[keep]
    L_oc, L_ch = polyline_length(contour), polyline_length(H)
    return L_oc / L_ch, L_oc, L_ch, contour, H


# ================================================================ loading
def _read_json(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rb') as fh:
        return json.loads(fh.read().decode('utf-8'))


def _to_days(area, pial):
    """Put the step time of a run in days.  Runs from the old 1-s step end at
    t ~ 1 and are multiplied by OLD_STEP_DAYS; 18-day runs are left as they
    are.  Returns the factor used (1 or OLD_STEP_DAYS)."""
    if area is not None:
        t_end = area['time'][-1]
    elif pial is not None:
        t_end = pial['time'][-1]
    else:
        return 1.0
    k = OLD_STEP_DAYS if t_end < 1.5 else 1.0
    if k != 1.0:
        if area is not None:
            area['time'] = area['time'] * k
        if pial is not None:
            pial['time'] = pial['time'] * k
    return k


def _run_from_json(key, e):
    area = None
    if 'area' in e:
        area = dict((k, np.asarray(v, float)) for k, v in e['area'].items())
    pial = None
    if 'pial' in e:
        p = e['pial']
        t = np.asarray(p['time'], float)
        ref = np.asarray(p['xyz'], float)
        U = np.zeros((len(t), len(ref), 3))
        for c, name in enumerate(('U1', 'U2', 'U3')):
            if name in p:
                U[:, :, c] = np.array(p[name], dtype=float)   # null -> nan
        pial = dict(time=t, ref=ref, U=U)
    tscale = _to_days(area, pial)
    return dict(key=key, label=key, style={}, props=e.get('props', {}),
                area=area, pial=pial, odb_mtime=e.get('odb_mtime'),
                source=e.get('odb', key), tscale=tscale)


_COL = re.compile(r'^\s*U([123])\s*@\s*([-+0-9.eE]+)\s*$')


def _load_pial_csv(path):
    with open(path) as fh:
        header = fh.readline().strip().split(',')
    data = np.atleast_2d(np.genfromtxt(path, delimiter=',', skip_header=1))
    ref = data[:, 1:4]
    cols = {}
    for j, name in enumerate(header):
        m = _COL.match(name)
        if m:
            cols.setdefault(float(m.group(2)), [None] * 3)[int(m.group(1)) - 1] = j
    times = np.array(sorted(cols))
    U = np.zeros((len(times), len(ref), 3))
    for k, t in enumerate(times):
        for c, j in enumerate(cols[t]):
            if j is not None:
                U[k, :, c] = data[:, j]
    return ref, times, U


def _run_from_csv(pial_csv):
    key = re.sub(r'_pial_disp$', '', os.path.splitext(os.path.basename(pial_csv))[0])
    ref, times, U = _load_pial_csv(pial_csv)
    area = None
    csa = re.sub(r'_pial_disp\.csv$', '_csa.csv', pial_csv)
    if os.path.isfile(csa):
        d = np.genfromtxt(csa, delimiter=',', names=True)
        area = dict((n, np.atleast_1d(d[n]).astype(float)) for n in d.dtype.names)
    pial = dict(time=times, ref=ref, U=U)
    tscale = _to_days(area, pial)
    return dict(key=key, label=key, style={}, props={}, area=area,
                pial=pial, odb_mtime=None, source=pial_csv, tscale=tscale)


def load_runs(sources):
    """Runs from any mix of .json / .json.gz results files and pial CSVs."""
    runs = []
    for src in sources:
        p = _resolve(src)
        if not os.path.isfile(p):
            print('not found: %s' % p)
            continue
        if p.endswith('.json') or p.endswith('.json.gz'):
            doc = _read_json(p)
            for key in sorted(doc.get('runs', {})):
                runs.append(_run_from_json(key, doc['runs'][key]))
        elif p.endswith('.csv'):
            runs.append(_run_from_csv(p))
    return runs


def match_runs(runs, pattern):
    """Keep runs whose name contains `pattern` (or matches it, if it has * ? []).
    Several patterns may be given, comma-separated."""
    if not pattern:
        return runs
    pats = [p.strip() for p in pattern.split(',') if p.strip()]

    def ok(key):
        return any(fnmatch.fnmatchcase(key, p) if any(c in p for c in '*?[')
                   else p in key for p in pats)
    out = [r for r in runs if ok(r['key'])]
    if not out:
        print("no run matches '%s'.  Available: %s"
              % (pattern, ', '.join(r['key'] for r in runs)))
    return out


def select_runs(runs, spec):
    """spec: None (all) | list of keys | list of (key, label[, style])."""
    if not spec:
        out = runs
    else:
        by_key = dict((r['key'], r) for r in runs)
        out = []
        for item in spec:
            key = item if isinstance(item, str) else item[0]
            if key not in by_key:
                print("run '%s' not in results.  Available: %s"
                      % (key, ', '.join(sorted(by_key))))
                continue
            r = dict(by_key[key])
            if not isinstance(item, str):
                r['label'] = item[1]
                r['style'] = dict(item[2]) if len(item) > 2 else {}
            out.append(r)
    for i, r in enumerate(out):
        r['style'] = dict(r.get('style', {}))
        r['style'].setdefault('color', 'C%d' % (i % 10))
        r['style'].setdefault('lw', 1.8)
    return out


# ================================================================ analysis
def sim_day(t):
    """Step time (days, after loading) -> postnatal day."""
    return np.asarray(t, float) * DAYS_PER_SIM + DAY_OFFSET


def garcia_csa(day, day_ref=None):
    """(V(day)/V(day_ref))^(2/3); day_ref defaults to the start of the simulation."""
    if day_ref is None:
        day_ref = sim_day(0.0)
    V = lambda d: G_A / (1.0 + np.exp(-G_RG * (np.asarray(d, float) - G_TMID)))
    return (V(day) / V(day_ref)) ** (2.0 / 3.0)


def _csa_ref():
    return sim_day(0.0) if CSA_REF_DAY is None else CSA_REF_DAY


def _z_planes(ref, z_tol=Z_TOL):
    """Distinct reference z-planes of the pial nodes, ascending."""
    zs = np.sort(ref[:, 2])
    planes = [zs[0]]
    for z in zs[1:]:
        if z - planes[-1] > z_tol:
            planes.append(z)
    return np.array(planes)


def _pick_z(ref, z_slice, z_tol=Z_TOL):
    if z_slice is None:
        return ref[:, 2].min()
    planes = _z_planes(ref, z_tol)
    if isinstance(z_slice, str) and z_slice == 'mid':
        return planes[np.argmin(np.abs(planes - 0.5 * (planes[0] + planes[-1])))]
    z = float(z_slice)
    if np.min(np.abs(planes - z)) > z_tol:
        raise ValueError('no pial nodes on z = %g; the planes are %s'
                         % (z, ', '.join('%.4g' % v for v in planes)))
    return z


def gi_series(run, hull=HULL, xy_min=XY_MIN, z_slice=Z_SLICE, z_tol=Z_TOL):
    p = run['pial']
    ref = p['ref']
    z0 = _pick_z(ref, z_slice, z_tol)
    m = ((np.abs(ref[:, 2] - z0) < z_tol)
         & (ref[:, 0] > xy_min) & (ref[:, 1] > xy_min))
    ref_xy = ref[m, :2]
    t_out, gi, loc, lch = [], [], [], []
    for k, t in enumerate(p['time']):
        d = p['U'][k][m, :2]
        if np.isnan(d).any():
            continue
        g, a, b, _, _ = gyrification_index(ref_xy, ref_xy + d, hull)
        t_out.append(t); gi.append(g); loc.append(a); lch.append(b)
    t_out = np.array(t_out)
    return dict(time=t_out, day=sim_day(t_out), GI=np.array(gi),
                L_oc=np.array(loc), L_ch=np.array(lch), n_pial=int(m.sum()))


def _cortex_law(run):
    """Name the cortical growth law from A_gm(end)/A_gm(0) at the end of the run."""
    a = run['area']
    if a is None or abs(a['time'][-1] - T_END) > 0.02 * T_END:
        return float('nan'), 'n/a'
    ratio = a['A_gm'][-1] / a['A_gm'][0]
    k = run.get('tscale', 1.0)                  # props of old 1-s runs are per sim-second
    gr = run['props'].get('GRAY', [])
    T1 = gr[3] * k if len(gr) > 3 else 2.0      # day
    d_end = DAYS_PER_SIM * T_END + VUMAT_DAY_OFFSET
    d_on = DAYS_PER_SIM * T1 + VUMAT_DAY_OFFSET
    SA = lambda d: 1.0 / (1.0 + np.exp(-SA_R * (d - SA_TMID)))
    sa = SA(d_end) / SA(d_on)
    h = lambda d: TH_HINF - (TH_HINF - TH_H1) * np.exp(-(d - 1.0) / TH_TAU)
    thT = max(1.0, h(d_end) / (TH_FEL * TH_H1))
    cands = [('SA law, Fg(3,3)=1', np.sqrt(sa)),
             ('SA + thickness, Fg(3,3)=1', np.sqrt(sa) * thT),
             ('SA law, out-of-plane NOT removed', sa),
             ('SA + thickness, out-of-plane NOT removed', sa * thT)]
    if len(gr) > 2:
        cands.append(('old linear law', (1 + gr[2] / k * (T_END - T1)) * 0.985))
    name, val = min(cands, key=lambda c: abs(np.log(ratio / c[1])))
    if abs(np.log(ratio / val)) > 0.06:
        name = 'unrecognised'
    return ratio, name


def report(runs):
    """Parameters, key numbers and consistency checks for every run."""
    w = max([len(r['key']) for r in runs] + [3])      # never truncate run names
    print('\n%-*s %7s %7s %9s %9s %9s %8s %8s  %s'
          % (w, 'run', 'gamma1', 'gamma', 'CSA(T2)', 'CSA(end)', 'end/T2', 'GI(end)',
             'A_gm x', 'cortex law'))
    t2 = []
    for r in runs:
        wh = r['props'].get('WHIT', [])
        g1 = wh[3] if len(wh) > 4 else float('nan')
        g = wh[4] if len(wh) > 4 else float('nan')
        c_t2 = c_end = float('nan')
        if r['area'] is not None:
            t, A = r['area']['time'], r['area']['A_total'] / r['area']['A_total'][0]
            c_t2, c_end = float(np.interp(T2_SIM, t, A)), A[-1]
            t2.append((r['key'], c_t2))
        gi_end = gi_series(r)['GI'][-1] if r['pial'] is not None else float('nan')
        ratio, law = _cortex_law(r)
        print('%-*s %7.3g %7.3g %9.3f %9.3f %9.3f %8.3f %8.3f  %s'
              % (w, r['key'], g1, g, c_t2, c_end, c_end / c_t2, gi_end, ratio, law))
    print('Garcia (from P%.1f): CSA(T2 = P%.1f) = %.3f   CSA(P%.1f) = %.3f   end/T2 = %.3f'
          '   P8->end = %.3f (the part Garcia et al. measured)'
          % (sim_day(0), sim_day(T2_SIM), garcia_csa(sim_day(T2_SIM)), sim_day(T_END),
             garcia_csa(sim_day(T_END)), garcia_csa(sim_day(T_END), sim_day(T2_SIM)),
             garcia_csa(sim_day(T_END), 8.0)))
    if len(t2) > 1:
        vals = np.array([v for _, v in t2])
        if (vals.max() - vals.min()) / vals.mean() > 0.005:
            print('WARNING: runs differ by %.1f%% at T2 -- they do not share '
                  'Phases 1+2 and the cortex, so push vs pull is not a '
                  'controlled comparison.'
                  % (100 * (vals.max() - vals.min()) / vals.mean()))
    seen = {}
    for r in runs:
        seen.setdefault((hypothesis(r), round(gamma_of(r), 6), run_set(r)[0]),
                        []).append(r['key'])
    dup = [v for v in seen.values() if len(v) > 1]
    if dup:
        print('WARNING: more than one run per hypothesis and gamma -- the figure '
              'overlays them:')
        for v in dup:
            print('    ' + ', '.join(v))
        print("  set MATCH (or --match) to plot one set, e.g. MATCH = 'Refactor'")
    laws = set(_cortex_law(r)[1] for r in runs) - {'n/a'}
    if len(laws) > 1:
        print('WARNING: runs use different cortical growth laws: %s'
              % ', '.join(sorted(laws)))
    print()


# ================================================================ plots
def _t2_line(ax):
    ax.axvline(sim_day(T2_SIM), color='0.6', ls='--', lw=1, zorder=0)


def plot_csa(runs, ax, day_max=None):
    day_max = sim_day(T_END) if day_max is None else day_max
    dd = np.linspace(sim_day(0.0), day_max, 300)
    ax.plot(dd, garcia_csa(dd), 'k-', lw=2.2, zorder=0, label='Garcia et al. (2024)')
    for r in runs:
        if r['area'] is None:
            continue
        a = r['area']
        ax.plot(sim_day(a['time']), a['A_total'] / a['A_total'][0],
                label=r['label'], **r['style'])
    _t2_line(ax)
    ax.set_xlim(sim_day(0.0), day_max)
    ax.set_xlabel('Postnatal day')
    ax.set_ylabel(r'$A(t)/A_0 = (V/V_0)^{2/3}$')
    ax.grid(alpha=0.3)
    ax.legend(frameon=False)


def plot_gi(runs, ax, hull=HULL, day_max=None):
    day_max = sim_day(T_END) if day_max is None else day_max
    for r in runs:
        if r['pial'] is None:
            continue
        s = gi_series(r, hull)
        ax.plot(s['day'], s['GI'], label=r['label'], **r['style'])
    _t2_line(ax)
    ax.set_xlim(sim_day(0.0), day_max)
    ax.set_xlabel('Postnatal day')
    ax.set_ylabel('Gyrification index (GI)')
    ax.grid(alpha=0.3)
    ax.legend(frameon=False)


def plot_summary(runs, hull=HULL):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.8, 4.2))
    plot_csa(runs, a1)
    a1.set_title('(a) Cross-sectional area', loc='left', fontsize=11)
    plot_gi(runs, a2, hull)
    a2.set_title('(b) Gyrification index', loc='left', fontsize=11)
    fig.tight_layout()
    return fig


def plot_contours(run, sim_times=(8.0, 18.0), hull=HULL,
                  xy_min=XY_MIN, z_slice=Z_SLICE, z_tol=Z_TOL):
    p = run['pial']
    ref = p['ref']
    z0 = _pick_z(ref, z_slice, z_tol)
    m = ((np.abs(ref[:, 2] - z0) < z_tol)
         & (ref[:, 0] > xy_min) & (ref[:, 1] > xy_min))
    ref_xy = ref[m, :2]
    fig, axs = plt.subplots(1, len(sim_times),
                            figsize=(4.0 * len(sim_times), 4.0), squeeze=False)
    for ax, ts in zip(axs[0], sim_times):
        k = int(np.argmin(np.abs(p['time'] - ts)))
        g, _, _, cont, H = gyrification_index(ref_xy, ref_xy + p['U'][k][m, :2], hull)
        ax.plot(cont[:, 0], cont[:, 1], '-', color='0.6', lw=2, label='Pial surface')
        ax.plot(H[:, 0], H[:, 1], 'r--', lw=1.8, label='Convex hull')
        ax.set_xlim(min(0, cont[:, 0].min()) - 0.1, cont[:, 0].max() + 0.1)
        ax.set_ylim(min(0, cont[:, 1].min()) - 0.1, cont[:, 1].max() + 0.1)
        ax.set_title('%s\nP%.1f   GI = %.3f'
                     % (run['label'], sim_day(p['time'][k]), g), fontsize=10)
        ax.set_aspect('equal'); ax.axis('off')
    axs[0][-1].legend(frameon=False, loc='lower left')
    fig.tight_layout()
    return fig


def _dist_to_polyline(p, P):
    """Shortest distance from point p to the polyline P."""
    a, b = P[:-1], P[1:]
    ab = b - a
    t = np.clip(np.einsum('ij,ij->i', p - a, ab) / np.maximum(np.einsum('ij,ij->i', ab, ab),
                                                                1e-30), 0, 1)
    return float(np.min(np.linalg.norm(a + t[:, None] * ab - p, axis=1)))


def z_layer_report(runs, hull=HULL):
    """Does the choice of z-slice matter?  For every run: GI at the last frame
    on every z-plane of pial nodes, and the largest in-plane difference in
    position between planes.  In true plane strain every plane is identical."""
    print('\nGI at the last frame on each z-plane (hull = %s)' % hull)
    for r in runs:
        p = r['pial']
        if p is None:
            continue
        ref = p['ref']
        planes = _z_planes(ref)
        gis = []
        for z in planes:
            s = gi_series(r, hull, z_slice=float(z))
            gis.append(s['GI'][-1] if len(s['GI']) else float('nan'))
        # in-plane position difference between planes, matched by (X, Y)
        k = -1
        base = np.abs(ref[:, 2] - planes[0]) < Z_TOL
        xy0 = ref[base, :2]
        u0 = p['U'][k][base, :2]
        dmax = 0.0
        for z in planes[1:]:
            m = np.abs(ref[:, 2] - z) < Z_TOL
            xy, u = ref[m, :2], p['U'][k][m, :2]
            for i in range(len(xy)):
                j = int(np.argmin(np.sum((xy0 - xy[i]) ** 2, axis=1)))
                if np.sum((xy0[j] - xy[i]) ** 2) < 1e-10:
                    dmax = max(dmax, float(np.linalg.norm(u[i] - u0[j])))
        txt = '  '.join('z=%.3g: %.4f' % (z, g) for z, g in zip(planes, gis))
        print('%-36s %d plane(s)   %s   max in-plane diff %.2e mm'
              % (r['key'][:36], len(planes), txt, dmax))
    print('Identical GI on every plane (diff ~0) -> plane strain holds and the '
          'slice does not matter.')


def diagnose_drops(runs, hull=HULL, tol=0.02, max_plots=6, plot=True):
    """Find frames where GI falls by more than `tol` (relative) from one frame
    to the next and say what moved: the pial length L_oc or the hull L_ch.
    Also re-computes the same two frames with the other hull method -- if the
    drop vanishes there, it is the hull algorithm, not the mechanics."""
    other = 'mirror' if hull != 'mirror' else 'karan'
    found = []
    print('\nGI drops > %.1f%% between consecutive frames (hull = %s)' % (100 * tol, hull))
    w = max([len(r['key']) for r in runs] + [3])
    print('%-*s %13s %15s %8s %8s %12s'
          % (w, 'run', 'day', 'GI', 'dL_oc', 'dL_ch', 'GI (%s)' % other))
    for r in runs:
        if r['pial'] is None:
            continue
        s = gi_series(r, hull)
        o = gi_series(r, other)
        g = s['GI']
        for k in range(1, len(g)):
            if g[k] < g[k - 1] * (1.0 - tol):
                dloc = s['L_oc'][k] / s['L_oc'][k - 1] - 1.0
                dlch = s['L_ch'][k] / s['L_ch'][k - 1] - 1.0
                print('%-*s %5.2f->%5.2f %6.3f->%6.3f %+7.1f%% %+7.1f%% %5.3f->%5.3f'
                      % (w, r['key'], s['day'][k - 1], s['day'][k], g[k - 1], g[k],
                         100 * dloc, 100 * dlch, o['GI'][k - 1], o['GI'][k]))
                found.append((g[k - 1] - g[k], r, s['time'][k - 1], s['time'][k]))
    if not found:
        print('  none')
        return None
    print('dL_oc ~ 0 and dL_ch > 0  -> the outline moved out (fold pattern changed,\n'
          '                            or the hull picked different crests)\n'
          'dL_oc < 0                 -> the pial surface itself got shorter\n'
          'drop gone with the other hull -> hull algorithm, not mechanics\n')
    if not plot:
        return found

    # biggest drop per run, largest first
    best = {}
    for d in sorted(found, key=lambda x: -x[0]):
        best.setdefault(d[1]['key'], d)
    picks = sorted(best.values(), key=lambda x: -x[0])[:max_plots]

    fig, axs = plt.subplots(len(picks), 2, figsize=(9.0, 4.2 * len(picks)), squeeze=False)
    for row, (_, r, t0, t1) in zip(axs, picks):
        p = r['pial']
        ref = p['ref']
        z0 = _pick_z(ref, Z_SLICE)
        m = (np.abs(ref[:, 2] - z0) < Z_TOL) & (ref[:, 0] > XY_MIN) & (ref[:, 1] > XY_MIN)
        ref_xy = ref[m, :2]
        k0 = int(np.argmin(np.abs(p['time'] - t0)))
        k1 = int(np.argmin(np.abs(p['time'] - t1)))
        res = [gyrification_index(ref_xy, ref_xy + p['U'][k][m, :2], hull) for k in (k0, k1)]
        # zoom on the place where the hull changed most: the hull vertex (before
        # or after) that lies farthest from the other frame's hull
        H0, H1 = res[0][4], res[1][4]
        d01 = [_dist_to_polyline(v, H0) for v in H1]
        d10 = [_dist_to_polyline(v, H1) for v in H0]
        if max(d01) >= max(d10):
            c, dmax = H1[int(np.argmax(d01))], max(d01)
        else:
            c, dmax = H0[int(np.argmax(d10))], max(d10)
        du = np.array([dmax])
        for ax, zoom in zip(row, (False, True)):
            for (g, _, _, cont, H), col, ls, day in zip(
                    res, ('0.6', 'k'), ('-', '-'), (sim_day(t0), sim_day(t1))):
                ax.plot(cont[:, 0], cont[:, 1], color=col, ls=ls, lw=1.6,
                        label='P%.2f  GI %.3f' % (day, g))
                ax.plot(H[:, 0], H[:, 1], color=col, ls='--', lw=1.0)
            ax.set_aspect('equal')
            if zoom:
                ax.set_xlim(c[0] - 0.8, c[0] + 0.8)
                ax.set_ylim(c[1] - 0.8, c[1] + 0.8)
                ax.set_title('zoom: where the hull changed most (%.3f mm)' % du.max(),
                             fontsize=10)
            else:
                allc = np.vstack([res[0][3], res[1][3]])
                ax.set_xlim(min(0, allc[:, 0].min()) - 0.1, allc[:, 0].max() + 0.1)
                ax.set_ylim(min(0, allc[:, 1].min()) - 0.1, allc[:, 1].max() + 0.1)
                ax.plot(*c, 'o', mfc='none', mec='r', ms=14)
                ax.set_title(r['key'], fontsize=10)
                ax.legend(frameon=False, fontsize=9, loc='lower left')
    fig.tight_layout()
    return fig


# ================================================================ figure (Fig. 8 style)
_TAIL = [(re.compile(r'^(\d+)em(\d+)$'), lambda a, b: int(a) * 10.0 ** -int(b)),
         (re.compile(r'^g?(\d+)p(\d+)$'), lambda a, b: float('%s.%s' % (a, b))),
         (re.compile(r'^g?(\d+(?:\.\d+)?)$'), lambda a: float(a))]


def run_set(run):
    """(set label, line style) from RUN_SETS; ('', '-') if none matches."""
    for text, label, ls in RUN_SETS:
        if text in run['key']:
            return label, ls
    return '', '-'


def hypothesis(run):
    k = run['key'].lower()
    for h in ('push', 'pull', 'control'):
        if h in k:
            return h
    return 'other'


def gamma_of(run):
    """gamma / gamma_hat: WHIT props(5) from the ODB, else parsed from the job
    name (..._5em2 -> 0.05, ..._50em2 -> 0.5, ..._g2p22 -> 2.22, ..._3 -> 3)."""
    wh = run['props'].get('WHIT', [])
    if len(wh) > 4:
        return float(wh[4])
    # '_Old50em2' -> '50em2': drop any letters in front of the number
    tail = re.sub(r'^[A-Za-z]+(?=\d)', '', run['key'].split('_')[-1])
    for rx, f in _TAIL:
        m = rx.match(tail)
        if m:
            return f(*m.groups())
    return float('nan')


def _series(run, row, hull):
    if row == 'csa':
        if run['area'] is None:
            return None
        a = run['area']
        day = sim_day(a['time'])
        return day, a['A_total'] / np.interp(_csa_ref(), day, a['A_total'])
    if run['pial'] is None:
        return None
    s = gi_series(run, hull)
    return s['day'], s['GI']


def _end_labels(ax, ends, ylim, colour, gap=0.045):
    """gamma value at the right end of each curve, nudged apart so labels never
    overlap; a short leader in the curve colour joins a nudged label to its curve."""
    span = ylim[1] - ylim[0]
    ends = sorted(ends)
    ys = [min(max(y, ylim[0]), ylim[1]) for y, _ in ends]
    pos = list(ys)
    for _ in range(50):                       # spread symmetric clusters apart
        moved = False
        for i in range(1, len(pos)):
            if pos[i] - pos[i - 1] < gap * span:
                mid = 0.5 * (pos[i] + pos[i - 1])
                pos[i - 1], pos[i] = mid - 0.5 * gap * span, mid + 0.5 * gap * span
                moved = True
        if not moved:
            break
    lo, hi = ylim[0] + 0.015 * span, ylim[1] - 0.015 * span
    if pos and pos[0] < lo:                   # keep every label inside the axes
        pos = [p + lo - pos[0] for p in pos]
    if pos and pos[-1] > hi:
        pos = [p - (pos[-1] - hi) for p in pos]
    x0, x1 = XLIM
    xl = x1 + 0.02 * (x1 - x0)
    for (y, g), yc, yp in zip(ends, ys, pos):
        ax.plot([x1, xl], [yc, yp], color=colour(g), lw=0.8, clip_on=False, zorder=3)
        ax.text(xl + 0.01 * (x1 - x0), yp, '%g' % g, va='center', ha='left',
                fontsize=8.5, color='0.2', clip_on=False)


def plot_figure(runs, rows=('csa', 'gi'), hull=HULL):
    """Manuscript figure: columns = Pushing | Pulling, rows = CSA and/or GI."""
    from matplotlib.lines import Line2D

    groups = {'push': [], 'pull': [], 'control': []}
    for r in runs:
        h = hypothesis(r)
        if h in groups:
            groups[h].append(r)
    for h in ('push', 'pull'):
        if GAMMA_SHOW:
            keep = [float(g) for g in GAMMA_SHOW]
            groups[h] = [r for r in groups[h]
                         if any(abs(gamma_of(r) - g) < 1e-6 * max(1, g) for g in keep)]
        groups[h].sort(key=gamma_of)
    gam = sorted(set(round(gamma_of(r), 6) for r in groups['push'] + groups['pull']
                     if gamma_of(r) == gamma_of(r)))
    cmap = plt.get_cmap(CMAP)
    # with several run sets in one figure, only the first (RUN_SETS order) gets
    # end-of-curve gamma labels; the others are identified by colour + line style
    order = dict((lab, i) for i, (_, lab, _) in enumerate(RUN_SETS))
    sets_here = sorted(set(run_set(r) for r in groups['push'] + groups['pull']),
                       key=lambda s: order.get(s[0], len(order)))
    first_set = sets_here[0] if sets_here else ('', '-')

    def colour(g):
        if g != g:                                   # gamma unknown
            return '0.45'
        if len(gam) < 2:
            return cmap(0.5 * (CMAP_RANGE[0] + CMAP_RANGE[1]))
        x = gam.index(round(g, 6)) / (len(gam) - 1.0)
        return cmap(CMAP_RANGE[0] + (CMAP_RANGE[1] - CMAP_RANGE[0]) * x)

    def control_style(r):
        for tag, st in CONTROL_STYLE.items():
            if r['key'].rstrip('_').upper().endswith('_' + tag):
                return st
        return r['key'], '0.6', '-'

    # one y-range per row, shared by both columns
    data = {}
    for r in runs:
        for row in rows:
            data[(r['key'], row)] = _series(r, row, hull)
    ymax = {}
    for row in rows:
        vals = [np.nanmax(v[1]) for (k, rw), v in data.items() if rw == row and v is not None]
        if row == 'gi':
            vals += [g for _, g in FERRET_GI]
            top = max(1.6, np.ceil((max(vals) + 0.02) / 0.2) * 0.2) if vals else 1.6
            ymax[row] = (1.0, top)
        else:
            vals += [garcia_csa(XLIM[1], _csa_ref())]
            lo = 1.0
            if _csa_ref() > sim_day(0.0) + 1e-9:   # curves start below 1 before it
                mins = [np.nanmin(v[1]) for (k, rw), v in data.items()
                        if rw == row and v is not None]
                lo = np.floor(min(mins + [garcia_csa(XLIM[0], _csa_ref())]) * 10) / 10
            ymax[row] = (lo, np.ceil(max(vals) * 1.05 / 0.5) * 0.5)

    with plt.rc_context(FIG_RC):
        fig, axs = plt.subplots(len(rows), 2, figsize=(9.2, 3.4 * len(rows)),
                                squeeze=False, sharex=True)
        letters = iter('ABCDEFGH')
        for i, row in enumerate(rows):
            for j, (hyp, title) in enumerate((('push', 'Pushing'), ('pull', 'Pulling'))):
                ax = axs[i][j]
                if row == 'csa':          # dashed on top: stays visible when a run matches it
                    d0 = np.linspace(XLIM[0], GARCIA_FIRST_DAY, 100)
                    d1 = np.linspace(GARCIA_FIRST_DAY, XLIM[1], 200)
                    ax.plot(d0, garcia_csa(d0, _csa_ref()), color='0.35', lw=1.3,
                            ls=(0, (1, 1.5)), zorder=5)
                    ax.plot(d1, garcia_csa(d1, _csa_ref()), color='k', lw=1.6,
                            ls=(0, (5, 3)), zorder=5)
                for r in groups['control']:
                    d = data[(r['key'], row)]
                    if d is not None:
                        _, col, ls = control_style(r)
                        ax.plot(d[0], d[1], color=col, ls=ls, lw=1.5, zorder=2)
                ends = []
                for k, r in enumerate(groups[hyp]):
                    d = data[(r['key'], row)]
                    if d is None:
                        continue
                    step = max(1, len(d[0]) // 10)
                    ax.plot(d[0], d[1], color=colour(gamma_of(r)), lw=1.6,
                            ls=run_set(r)[1], marker=MARKERS[hyp], ms=4.5, mew=0,
                            markevery=((k * step) // max(1, len(groups[hyp])), step),
                            zorder=3)
                    if run_set(r) == first_set:     # label one set only
                        ends.append((d[1][-1], gamma_of(r)))
                if row in END_LABEL_ROWS and ends:
                    _end_labels(ax, ends, ymax[row], colour)
                if not groups[hyp]:
                    ax.text(0.5, 0.5, 'no %s runs' % hyp, transform=ax.transAxes,
                            ha='center', va='center', color='0.5')
                if row == 'gi':
                    ax.plot([d for d, _ in FERRET_GI], [g for _, g in FERRET_GI],
                            ls='none', marker='s', ms=5.5, color=FERRET_COLOR, zorder=4)
                if SHOW_T2:
                    ax.axvline(sim_day(T2_SIM), color='0.75', ls='--', lw=0.9, zorder=0)
                ax.set_xlim(*XLIM)
                ticks = list(XTICKS)
                labels = ['%g' % v for v in ticks]
                if DAY_OFFSET < 0 and START_LABEL:
                    ticks, labels = [DAY_OFFSET] + ticks, [START_LABEL] + labels
                ax.set_xticks(ticks)
                ax.set_xticklabels(labels)
                ax.set_ylim(*ymax[row])
                if row == 'gi':
                    ax.set_yticks(np.arange(1.0, ymax[row][1] + 1e-9, 0.2))
                ax.grid(axis='y', color='0.85', lw=0.8)
                ax.set_axisbelow(True)
                if j == 0:
                    ax.set_ylabel('GI [-]' if row == 'gi'
                                  else (r'Cross-section area $A/A_0$ [-]'
                                        if CSA_REF_DAY is None else
                                        r'Cross-section area $A/A_{\mathrm{P%g}}$ [-]'
                                        % CSA_REF_DAY))
                if i == len(rows) - 1:
                    ax.set_xlabel('Postnatal age [day]')
                if i == 0:
                    ax.set_title(title, fontweight='bold', pad=10)
                ax.text(-0.13, 1.06, next(letters), transform=ax.transAxes,
                        fontsize=15, fontweight='bold', family='sans-serif',
                        va='bottom', ha='left')

        handles = []
        if 'csa' in rows:
            handles.append(Line2D([], [], color='k', lw=1.6, ls=(0, (5, 3)),
                                  label='Garcia et al. (2024)'))
            handles.append(Line2D([], [], color='0.35', lw=1.3, ls=(0, (1, 1.5)),
                                  label='Extrapolated (< P%g)' % GARCIA_FIRST_DAY))
        for r in groups['control']:
            lab, col, ls = control_style(r)
            handles.append(Line2D([], [], color=col, ls=ls, lw=1.5, label=lab))
        for g in gam:
            handles.append(Line2D([], [], color=colour(g), lw=1.8,
                                  label=r'$\gamma,\ \hat{\gamma}$ = %g' % g))
        present = []
        for r in groups['push'] + groups['pull']:
            lab_ls = run_set(r)
            if lab_ls[0] and lab_ls not in present:
                present.append(lab_ls)
        if len(present) > 1:
            for lab, ls in present:
                handles.append(Line2D([], [], color='0.3', ls=ls, lw=1.5, label=lab))
        handles.append(Line2D([], [], color='0.3', lw=1.0, marker='o', ms=4.5,
                              label='Pushing'))
        handles.append(Line2D([], [], color='0.3', lw=1.0, marker='^', ms=4.5,
                              label='Pulling'))
        if 'gi' in rows:
            handles.append(Line2D([], [], ls='none', marker='s', ms=5.5,
                                  color=FERRET_COLOR, label='Ferrets'))
        fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.80, 0.5),
                   frameon=False, handlelength=1.8)
        fig.tight_layout(rect=(0, 0, 0.80, 1))
    return fig


# ================================================================ export
def export(runs, path, hull=HULL):
    """ONE long-format CSV for all runs: run, time, day, CSA, GI, L_oc, L_ch.
    'time' is the step time in days (old 1-s runs already converted)."""
    with open(path, 'w') as fh:
        fh.write('run,time,day,CSA_ratio,GI,L_oc,L_ch\n')
        for r in runs:
            s = gi_series(r, hull) if r['pial'] is not None else None
            if r['area'] is not None:
                a = r['area']
                csa_t, csa = a['time'], a['A_total'] / a['A_total'][0]
            times = s['time'] if s is not None else r['area']['time']
            for i, t in enumerate(times):
                c = float(np.interp(t, csa_t, csa)) if r['area'] is not None else np.nan
                g, lo, lc = ((s['GI'][i], s['L_oc'][i], s['L_ch'][i])
                             if s is not None else (np.nan,) * 3)
                fh.write('%s,%.6g,%.4g,%.6g,%.6g,%.6g,%.6g\n'
                         % (r['key'], t, sim_day(t), c, g, lo, lc))
    print('wrote %s' % os.path.abspath(path))


def show_or_save(outfile=None, dpi=200):
    if outfile is None and not _HEADLESS:
        plt.show()
        return
    outfile = outfile or 'gi_figure.png'
    plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    print('wrote %s' % os.path.abspath(outfile))


# ================================================================ CLI
def main(argv):
    sources, outfile, hull, match = [], None, HULL, MATCH
    mode = {'overlay': 'summary'}.get(MODE, MODE)
    contours, run_key, spec, exp, listing = None, None, RUNS, None, False
    rows = ('csa', 'gi')
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ('--save', '-o'):
            outfile = argv[i + 1]; i += 2
        elif a == '--hull':
            hull = argv[i + 1]; i += 2
        elif a == '--gi-only':
            rows = ('gi',); i += 1
        elif a == '--csa-only':
            rows = ('csa',); i += 1
        elif a == '--overlay':
            mode = 'summary'; i += 1
        elif a == '--drops':
            mode = 'drops'; i += 1
        elif a == '--zcheck':
            mode = 'zcheck'; i += 1
        elif a == '--match':
            match = argv[i + 1]; i += 2
        elif a == '--runs':
            spec = [k.strip() for k in argv[i + 1].split(',') if k.strip()]; i += 2
        elif a == '--run':
            run_key = argv[i + 1]; i += 2
        elif a == '--export':
            exp = argv[i + 1]; i += 2
        elif a == '--list':
            listing = True; i += 1
        elif a == '--contours':
            mode, contours = 'contours', []
            i += 1
            while i < len(argv) and not argv[i].startswith('-'):
                contours.append(float(argv[i])); i += 1
        else:
            sources.append(a); i += 1

    default = list(RESULTS) if isinstance(RESULTS, (list, tuple)) else [RESULTS]
    runs = match_runs(load_runs(sources or default), match)
    if not runs:
        print('No runs loaded.  Put %s in\n  %s\nor pass a file name.'
              % (RESULTS, DATA_DIR))
        return

    if listing:
        w = max([len(r['key']) for r in runs] + [3])
        print('\n%-*s %7s %7s %7s %6s  %-19s %s'
              % (w, 'run', 'frames', 'gamma1', 'gamma', 'step', 'odb written', 'source'))
        for r in runs:
            wh = r['props'].get('WHIT', [])
            nf = len(r['area']['time']) if r['area'] is not None else 0
            print('%-*s %7d %7s %7s %6s  %-19s %s'
                  % (w, r['key'], nf,
                     '%g' % wh[3] if len(wh) > 4 else '-',
                     '%g' % wh[4] if len(wh) > 4 else '-',
                     '1 s' if r.get('tscale', 1.0) != 1.0 else '18 d',
                     r['odb_mtime'] or '-', r['source']))
        return

    runs = select_runs(runs, spec)
    if not runs:
        return
    report(runs)
    if mode == 'zcheck':
        z_layer_report(runs, hull)
        return
    if mode != 'drops':                 # always say whether GI ever falls
        diagnose_drops(runs, hull, plot=False)

    if exp:
        export(runs, exp, hull)
    if mode == 'contours':
        pick = [r for r in runs if r['key'] == run_key] if run_key else runs[:1]
        if not pick:
            print("run '%s' not found" % run_key); return
        plot_contours(pick[0], contours or (8.0, 18.0), hull)
    elif mode == 'summary':
        plot_summary(runs, hull)
    elif mode == 'drops':
        if diagnose_drops(runs, hull) is None:
            return
    else:
        plot_figure(runs, rows, hull)
    show_or_save(outfile)


if __name__ == '__main__':
    main(sys.argv[1:])