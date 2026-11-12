"""
evol_area_all.py  --  ODB extractor, run inside Abaqus.

One pass over the frames of each ODB; every run goes into ONE results file.

    abaqus python extract_odb_all.py                    -> the 10 Logis runs:
                                                         Logis_*_Z1_*.odb minus *_5N_*
                                                         -> Logis_brain_results.json
    abaqus python extract_odb_all.py Logis_*_Z1_*.odb --out Logis_brain_results.json --new
                                                      -> the same 10, file rebuilt
    abaqus python extract_odb_all.py --prune --out Logis_brain_results.json
                                                      -> only REMOVES the *_5N_* runs
                                                         already in the file (no extraction)

ODBs whose name matches EXCLUDE (default '*_5N_*') are skipped, whether they
come from DEFAULT_ODBS or from the command line -- the shell expands
Logis_*_Z1_*.odb to the 5N runs too.  The skipped files are listed.
    --exclude PATTERN       skip these too (repeatable, e.g. --exclude "*_control_*")
    --no-exclude            take every ODB given (also the *_5N_* ones)
    --prune                 drop entries matching EXCLUDE from the results file



Each run is stored under its job name (ODB name without .odb).  Extracting an
ODB whose name is already in the file REPLACES that entry (the script says so,
with the old and new ODB dates); every other entry is kept.  To keep an older
version of a job next to a re-run with the same name, extract with a tag:
    abaqus python extract_odb_all.py Const_*_Z1_*.odb --tag _v2
        -> entries Const_push_Z1_1_v2, ...  (the untagged ones stay)
The output file is RESULTS_JSON unless --out is given; it is printed first.

Run entries:

    props      WHIT / GRAY user-material constants read FROM THE ODB
               (so the gamma values in the file are the ones that ran)
    area       time, A_gm, A_wm, A_total, A_model        [mm^2]
               A_set = sum(EVOL) / depth, depth = z-extent of the brain mesh
               (read from the ODB, so old and new models can differ)
    pial       time, node, xyz (reference), U1[frame][node], U2[...]
               (U3 only if it is non-zero anywhere)
    interface  the same for the gray/white interface nodes (set INTERFACE-FACE),
               used by sulcal_depth.py for the cortical thickness.
               Skipped with a note if the ODB has no such node set.
    odb_mtime  when the ODB was written   } provenance: a stale entry
    extracted  when this entry was made    } cannot pass as a new one

GI_analysis.py reads this file directly (CSA vs Garcia + gyrification index);
sulcal_depth.py reads it for sulcal depth and cortical thickness.
Requires EVOL (areas) and U (node displacements) in the field-output request.

    --no-interface          do not store the interface nodes
    --interface-set NAME    node set of the gray/white interface (default INTERFACE-FACE)
"""

import fnmatch
import glob
import gzip
import json
import math
import os
import sys
import time

# ---------------------------------------------------------------- defaults
# Used when no ODB is named on the command line: every ODB in the current
# directory that matches the pattern.  Currently the 10 OLD runs (draft):
#   QuarterEllipse_push_Z1_5N_Old<gamma>.odb, QuarterEllipse_pull_Z1_5N_Old<gamma>.odb
# Their results go to a separate file so they never mix with the new runs.
DEFAULT_ODBS = 'Logis_*_Z1_*.odb'       # Logis_push_Z1_3.odb, Logis_pull_Z1_5em2.odb, ...
RESULTS_JSON = 'Logis_brain_results.json'
# Skipped wherever the ODB list comes from (fnmatch on the file name).  The
# *_Z1_5N_* jobs are a different model version (push: A_gm x2.7 at P16 instead
# of x7.25) and must not sit in the same results file.
EXCLUDE      = ('*_5N_*',)
# older sets:
#   DEFAULT_ODBS = 'QuarterEllipse_*Refactor*.odb'; RESULTS_JSON = 'brain_results_old.json'
#   DEFAULT_ODBS = 'Const_*_Z1_*.odb';               RESULTS_JSON = 'brain_results.json'
FORMAT       = 'brain_growth_results/1'

# ---------------------------------------------------------------- geometry
MAJ_AXIS = 3.6      # outer semi-major axis of the quarter ellipse [mm]
MIN_AXIS = 3.0      # outer semi-minor axis                        [mm]
DEPTH    = 0.02     # fallback thickness [mm], used only if the mesh is 2-D
                    # (normally measured from the mesh; --depth overrides)

DEFAULT_SETS = (('CORTEX', 'A_gm'), ('SUBCORTEX', 'A_wm'))
PIAL_SET     = 'CORTEX_TOP'             # outer (pial) face of the cortex
# gray/white interface: the first of these node sets that exists is used
INTERFACE_SET = 'INTERFACE-FACE'
INTERFACE_FALLBACK = ('INTERFACE-FACE', 'INTERFACE_FACE', 'INTERFACE', 'INTERFACE-EDGE')


def _r(x, nd=7):
    """Round for storage: 1e-7 mm is far below anything we measure."""
    return round(float(x), nd)


def _iso(t):
    return time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(t))


# ================================================================ ODB helpers
def _find_set(odb, name, kind):
    """Element or node set ('elementSets' / 'nodeSets'), assembly or
    instance level, case-insensitive."""
    key = name.upper()
    ra = odb.rootAssembly
    repo = getattr(ra, kind)
    if key in repo.keys():
        return repo[key]
    for inst in ra.instances.values():
        repo = getattr(inst, kind)
        if key in repo.keys():
            return repo[key]
    avail = list(getattr(ra, kind).keys())
    for inst in ra.instances.values():
        avail += ['%s.%s' % (inst.name, k) for k in getattr(inst, kind).keys()]
    raise KeyError("%s '%s' not in ODB.  Available: %s" % (kind, name, avail))


def _find_first_set(odb, names, kind):
    """(name, set) of the first set in `names` that exists, else (None, None)."""
    for nm in names:
        try:
            return nm, _find_set(odb, nm, kind)
        except KeyError:
            continue
    return None, None


def _flatten(seq):
    seq = list(seq)
    if seq and not hasattr(seq[0], 'label'):
        return [o for sub in seq for o in sub]
    return seq


def _elem_keys(elset):
    return set((e.instanceName, e.label) for e in _flatten(elset.elements))


def _material_props(odb):
    """User-material constants as stored in the ODB (empty if unavailable)."""
    out = {}
    try:
        for name, mat in odb.materials.items():
            try:
                out[name] = [float(c) for c in mat.userMaterial.mechanicalConstants]
            except Exception:
                pass
    except Exception:
        pass
    return out


def _mesh_depth(odb, inst_names):
    """Out-of-plane thickness = z-extent of the nodes of the brain instance(s)
    (the rigid walls are separate instances and are left out)."""
    zmin = zmax = None
    for name in inst_names:
        for n in odb.rootAssembly.instances[name].nodes:
            c = n.coordinates
            if len(c) < 3:
                continue
            z = float(c[2])
            zmin = z if zmin is None else min(zmin, z)
            zmax = z if zmax is None else max(zmax, z)
    return 0.0 if zmin is None else zmax - zmin


# ================================================================ extraction
def extract(odb_path, step_name=None, depth=None, sets=DEFAULT_SETS,
            pial=True, pial_set=PIAL_SET, interface=True,
            interface_set=INTERFACE_SET, verbose=True):
    """Read one ODB and return its run entry (a plain dict)."""
    from odbAccess import openOdb

    odb = openOdb(str(odb_path), readOnly=True)
    try:
        if step_name is None:
            step_name = list(odb.steps.keys())[-1]
        step = odb.steps[step_name]

        labels = [lab for _, lab in sets]
        keys = dict((lab, _elem_keys(_find_set(odb, nm, 'elementSets')))
                    for nm, lab in sets)

        depth_src = 'given'
        if depth is None:
            insts = set(inst for lab in labels for inst, _ in keys[lab])
            depth, depth_src = _mesh_depth(odb, insts), 'mesh'
            if depth < 1e-9:
                depth, depth_src = DEPTH, 'fallback'

        # node sets whose displacements are stored: pial surface and the
        # gray/white interface, both in the same layout
        series = {}
        if pial:
            series['pial'] = {'set_name': pial_set,
                              'set': _find_set(odb, pial_set, 'nodeSets')}
        if interface:
            names = [interface_set] + [n for n in INTERFACE_FALLBACK if n != interface_set]
            nm, ns = _find_first_set(odb, names, 'nodeSets')
            if ns is None:
                print('  note: no interface node set (%s) -- cortical thickness '
                      'will not be available for this run' % ', '.join(names))
            else:
                series['interface'] = {'set_name': nm, 'set': ns}
        for s in series.values():
            s['nodes'] = _flatten(s['set'].nodes)
            s['labels'] = [n.label for n in s['nodes']]
            s.update(time=[], U1=[], U2=[], U3=[], u3_max=0.0)

        area = dict((k, []) for k in ['time'] + labels + ['A_total', 'A_model'])

        for fr in step.frames:
            fo = fr.fieldOutputs

            if 'EVOL' in fo.keys():
                # overlapping output requests can write an element more than
                # once -> de-duplicate on (instance, label) before summing
                evol = {}
                for v in fo['EVOL'].values:
                    evol[(v.instance.name, v.elementLabel)] = v.data
                tot = 0.0
                for lab in labels:
                    a = sum(evol[k] for k in keys[lab] if k in evol) / depth
                    area[lab].append(_r(a, 8))
                    tot += a
                area['A_total'].append(_r(tot, 8))
                area['A_model'].append(_r(sum(evol.values()) / depth, 8))
                area['time'].append(_r(fr.frameValue, 8))

            if series and 'U' in fo.keys():
                for s in series.values():
                    u = {}
                    for v in fo['U'].getSubset(region=s['set']).values:
                        u[v.nodeLabel] = v.data
                    row1, row2, row3 = [], [], []
                    for l in s['labels']:
                        d = u.get(l)
                        if d is None:
                            row1.append(None); row2.append(None); row3.append(None)
                            continue
                        row1.append(_r(d[0])); row2.append(_r(d[1]))
                        z = float(d[2]) if len(d) > 2 else 0.0
                        row3.append(_r(z)); s['u3_max'] = max(s['u3_max'], abs(z))
                    s['U1'].append(row1); s['U2'].append(row2); s['U3'].append(row3)
                    s['time'].append(_r(fr.frameValue, 8))

        if not area['time']:
            raise RuntimeError(
                "no EVOL field output in step '%s'.  Add 'EVOL' to the field "
                "output request and re-run the job." % step_name)

        run = {
            'odb': os.path.basename(str(odb_path)),
            'odb_mtime': _iso(os.path.getmtime(str(odb_path))),
            'extracted': _iso(time.time()),
            'step': step_name,
            'depth': depth,
            'depth_source': depth_src,
            'props': _material_props(odb),
            'area': area,
        }
        for key, s in series.items():
            if not s['time']:
                continue
            run[key] = {
                'set': s['set_name'],
                'time': s['time'],
                'node': s['labels'],
                'xyz': [[_r(c) for c in n.coordinates] for n in s['nodes']],
                'U1': s['U1'],
                'U2': s['U2'],
            }
            if s['u3_max'] > 1e-12:
                run[key]['U3'] = s['U3']

        if verbose:
            _report(odb_path, run, step, labels, depth)
        return run
    finally:
        odb.close()


# ================================================================ output
def _read(path):
    if not os.path.isfile(path):
        return None
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rb') as fh:
        data = fh.read()
    if not isinstance(data, str):
        data = data.decode('utf-8')
    return json.loads(data)


def _write(doc, path):
    """Written to a temporary file first, so an interrupted write never
    corrupts the existing results."""
    doc['updated'] = _iso(time.time())
    data = json.dumps(doc, separators=(',', ':'), sort_keys=True)
    if not isinstance(data, bytes):
        data = data.encode('utf-8')
    tmp = path + '.tmp'
    opener = gzip.open if path.endswith('.gz') else open
    with opener(tmp, 'wb') as fh:
        fh.write(data)
    if os.path.exists(path):
        os.remove(path)
    os.rename(tmp, path)


def save_run(key, run, path=RESULTS_JSON, merge=True):
    """Add/replace one run in the results file."""
    doc = _read(path) if merge else None
    if not doc or 'runs' not in doc:
        doc = {'format': FORMAT, 'runs': {}}
    prev = doc['runs'].get(key)
    doc['runs'][key] = run
    _write(doc, path)
    return doc, prev


# ================================================================ selection
def _excluded(name, patterns):
    """The pattern that `name` (file name or run key) matches, else None."""
    base = os.path.basename(name)
    stem = os.path.splitext(base)[0]
    for pat in patterns:
        if fnmatch.fnmatchcase(base, pat) or fnmatch.fnmatchcase(stem, pat):
            return pat
    return None


def select_odbs(odbs, patterns):
    """(kept, skipped) -- duplicates dropped, order kept."""
    kept, skipped, seen = [], [], set()
    for o in odbs:
        if o in seen:
            continue
        seen.add(o)
        if patterns and _excluded(o, patterns):
            skipped.append(o)
        else:
            kept.append(o)
    return kept, skipped


def _excluded_runs(doc, patterns):
    """Keys of the runs in `doc` whose job name or ODB name matches."""
    hits = []
    for key in sorted((doc or {}).get('runs', {})):
        odb = doc['runs'][key].get('odb', key)
        if _excluded(key, patterns) or _excluded(odb, patterns):
            hits.append(key)
    return hits


def prune(path, patterns):
    """Remove the runs matching `patterns` from the results file."""
    doc = _read(path)
    if not doc or 'runs' not in doc:
        print('prune: %s does not exist or holds no runs' % path)
        return doc
    hits = _excluded_runs(doc, patterns)
    if not hits:
        print('prune: nothing in %s matches %s (%d runs kept)'
              % (os.path.basename(path), ', '.join(patterns), len(doc['runs'])))
        return doc
    for key in hits:
        del doc['runs'][key]
    _write(doc, path)
    print('prune: removed %d run(s) matching %s from %s:'
          % (len(hits), ', '.join(patterns), os.path.basename(path)))
    for key in hits:
        print('    - %s' % key)
    print('  %d run(s) left\n' % len(doc['runs']))
    return doc


def write_legacy_csv(key, run, outdir='.'):
    """<job>_csa.csv and <job>_pial_disp.csv, as before (Karan's layout)."""
    a = run['area']
    names = [k for k in ('A_gm', 'A_wm', 'A_total', 'A_model') if k in a]
    with open(os.path.join(outdir, key + '_csa.csv'), 'w') as fh:
        fh.write(','.join(['time'] + names) + '\n')
        for i, t in enumerate(a['time']):
            fh.write(','.join('%.8g' % x for x in [t] + [a[k][i] for k in names]) + '\n')
    p = run.get('pial')
    if not p:
        return
    comps = [p['U1'], p['U2'], p.get('U3')]
    with open(os.path.join(outdir, key + '_pial_disp.csv'), 'w') as fh:
        hdr = ['Node', 'X', 'Y', 'Z']
        for t in p['time']:
            hdr += ['U%d @ %.6g' % (c, t) for c in (1, 2, 3)]
        fh.write(','.join(hdr) + '\n')
        for i, lab in enumerate(p['node']):
            row = list(p['xyz'][i])
            for k in range(len(p['time'])):
                for c in comps:
                    v = 0.0 if c is None else c[k][i]
                    row.append(float('nan') if v is None else v)
            fh.write(str(lab) + ',' + ','.join('%.8g' % x for x in row) + '\n')


# ================================================================ report
def _report(odb_path, run, step, labels, depth):
    a = run['area']
    a_ref = math.pi * MAJ_AXIS * MIN_AXIS / 4.0
    tot = a['A_total']
    print('%s' % os.path.basename(str(odb_path)))
    print('  step %s, %d frames, t = %.4g .. %.4g,  depth = %.4g mm (%s)'
          % (run['step'], len(a['time']), a['time'][0], a['time'][-1],
             depth, run.get('depth_source', '?')))
    wh = run['props'].get('WHIT', [])
    if len(wh) > 4:
        print('  WHIT props: gamma_1 = %g   gamma/gamma_hat = %g' % (wh[3], wh[4]))
    else:
        mats = ['%s(%d)' % (k, len(v)) for k, v in sorted(run['props'].items())]
        print('  no WHIT material; user materials in ODB: %s'
              % (', '.join(mats) or 'none'))
    print('  t=0   ' + '  '.join('%s=%.4f' % (l, a[l][0]) for l in labels)
          + '  A_total=%.4f  (pi*a*b/4 = %.4f)' % (tot[0], a_ref))
    print('  t=end ' + '  '.join('%s=%.4f' % (l, a[l][-1]) for l in labels)
          + '  A_total=%.4f  (x%.3f)   A_gm x%.3f'
          % (tot[-1], tot[-1] / tot[0], a['A_gm'][-1] / a['A_gm'][0]))
    if abs(tot[0] - a['A_model'][0]) / tot[0] > 1e-6:
        print('  WARNING: sets sum to %.4f but the whole model is %.4f'
              % (tot[0], a['A_model'][0]))
    if abs(tot[0] - a_ref) / a_ref > 0.02:
        print('  WARNING: A_total(0) is %.1f%% off pi*a*b/4 -- check DEPTH'
              % (100.0 * (tot[0] / a_ref - 1.0)))
    for key in ('pial', 'interface'):
        p = run.get(key)
        if p:
            print('  %s: %d nodes x %d frames  (set %s)%s'
                  % (key, len(p['node']), len(p['time']), p['set'],
                     '' if 'U3' not in p else '  (U3 non-zero: stored)'))


def _summarise_file(doc):
    """What the results file holds, grouped by job-name family, so runs from
    different model versions sitting in one file are easy to spot."""
    import re
    groups = {}
    for key, e in doc['runs'].items():
        k = re.sub(r'(push|pull|control)', '*', key)
        m = re.match(r'^(.*?)(\d+(?:em\d+)?)(_[A-Za-z][^_]*)?$', k)
        fam = (m.group(1) + '*' + (m.group(3) or '')) if m else k
        groups.setdefault(fam, []).append(e.get('odb_mtime') or '?')
    print('\nthe file now holds:')
    for fam in sorted(groups):
        dates = sorted(groups[fam])
        print('  %-36s %2d run(s)   ODBs written %s .. %s'
              % (fam, len(dates), dates[0][:16], dates[-1][:16]))
    if len(groups) > 1:
        print('  (several families in one file: GI_analysis plots them together '
              'unless MATCH / --match picks one)')


# ================================================================ CLI
def _explain_no_match():
    """Say what IS here, so a wrong folder and a wrong pattern look different."""
    here = sorted(glob.glob('*.odb'))
    if not here:
        print('There are no .odb files at all in this folder.  The ODBs are '
              'somewhere else\n(or were moved/deleted).  Find them with:\n'
              '    find /groups/mhollan5/commandlab/Chongyi -name "*.odb" 2>/dev/null | head -40\n'
              'then cd there, or pass them with their path:\n'
              '    abaqus python evol_area_all.py /path/to/*.odb --out %s'
              % os.path.abspath(RESULTS_JSON))
        return
    print('This folder has %d .odb file(s), but none match the pattern:' % len(here))
    for f in here[:15]:
        print('    ' + f)
    if len(here) > 15:
        print('    ... and %d more' % (len(here) - 15))
    print('Edit DEFAULT_ODBS at the top of the script, or name them on the '
          'command line.')


def _parse_args(argv):
    if '--' in argv:
        argv = argv[argv.index('--') + 1:]
    odbs = []
    opts = {'depth': None, 'step_name': None, 'pial': True, 'interface': True,
            'interface_set': INTERFACE_SET}
    out, new, csv = RESULTS_JSON, False, False
    sel = {'exclude': list(EXCLUDE), 'prune': False}
    opts_tag = ['']
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--depth':
            opts['depth'] = float(argv[i + 1]); i += 2
        elif a == '--step':
            opts['step_name'] = argv[i + 1]; i += 2
        elif a == '--out':
            out = argv[i + 1]; i += 2
        elif a == '--new':
            new = True; i += 1
        elif a == '--tag':
            opts_tag[0] = argv[i + 1]; i += 2
        elif a == '--csv':
            csv = True; i += 1
        elif a == '--exclude':
            sel['exclude'].append(argv[i + 1]); i += 2
        elif a == '--no-exclude':
            sel['exclude'] = []; i += 1
        elif a == '--prune':
            sel['prune'] = True; i += 1
        elif a == '--no-pial':
            opts['pial'] = False; i += 1
        elif a == '--no-interface':
            opts['interface'] = False; i += 1
        elif a == '--interface-set':
            opts['interface_set'] = argv[i + 1]; i += 2
        elif a in ('-h', '--help'):
            print(__doc__); sys.exit(0)
        elif a.lower().endswith('.odb'):
            # a quoted pattern reaches us unexpanded (e.g. through abaqus cae)
            odbs.extend(sorted(glob.glob(a)) if glob.has_magic(a) else [a]); i += 1
        else:
            i += 1                      # Abaqus launcher flag -- ignore
    return odbs, opts, out, new, csv, opts_tag[0], sel


if __name__ == '__main__':
    odbs, opts, out, new, csv, tag, sel = _parse_args(sys.argv[1:])
    print('results file: %s%s%s\n' % (os.path.abspath(out),
                                        '  (started afresh)' if new else '',
                                        "  (keys tagged '%s')" % tag if tag else ''))
    if sel['prune']:
        if not sel['exclude']:
            print('--prune needs exclude patterns (EXCLUDE or --exclude)')
            sys.exit(1)
        prune(out, sel['exclude'])
        if not odbs:                    # prune only: nothing to extract
            doc = _read(out)
            if doc and doc.get('runs'):
                _summarise_file(doc)
            sys.exit(0)
    if not odbs:
        odbs = sorted(glob.glob(DEFAULT_ODBS))
        print("%d ODB(s) in %s match '%s'"
              % (len(odbs), os.getcwd(), DEFAULT_ODBS))
        if not odbs:
            _explain_no_match()
            sys.exit(1)
    odbs, skipped = select_odbs(odbs, sel['exclude'])
    if skipped:
        print('skipped %d ODB(s) matching %s  (--no-exclude to take them):'
              % (len(skipped), ', '.join(sel['exclude'])))
        for o in skipped:
            print('    - %s' % os.path.basename(o))
    print('extracting %d ODB(s):  %s\n'
          % (len(odbs), ', '.join([os.path.basename(o) for o in odbs])))
    if not odbs:
        print('nothing left to extract.')
        sys.exit(1)
    merge = not new
    done, failed = [], []
    for o in odbs:
        if not os.path.isfile(o):
            print("ERROR: '%s' not found in %s" % (o, os.getcwd()))
            failed.append(o)
            continue
        key = os.path.splitext(os.path.basename(o))[0] + tag
        try:
            run = extract(o, **opts)
        except Exception as exc:        # one bad ODB must not stop the batch
            print('ERROR in %s: %s\n' % (o, exc))
            failed.append(o)
            continue
        doc, prev = save_run(key, run, out, merge=merge)
        done.append(key)
        merge = True                    # later ODBs in this call always add
        if csv:
            write_legacy_csv(key, run)
        if prev is None:
            what = 'added'
        elif prev.get('odb_mtime') == run['odb_mtime']:
            what = 'replaced (same ODB, re-extracted)'
        else:
            what = 'REPLACED: entry was from an ODB written %s, now %s' % (
                prev.get('odb_mtime', '?'), run['odb_mtime'])
        print('  -> %s: %s %s  [%d runs in file]\n'
              % (os.path.basename(out), key, what, len(doc['runs'])))
    print('extracted %d of %d ODB(s) into %s' % (len(done), len(odbs), out))
    if done:
        _summarise_file(doc)
        left = _excluded_runs(doc, sel['exclude'])
        if left:
            print('\n*** the file still holds %d run(s) from an earlier call that'
                  ' match %s:\n    %s\n    remove them with:  abaqus python %s'
                  ' --prune --out %s\n    (or rebuild the file with --new)'
                  % (len(left), ', '.join(sel['exclude']), ', '.join(left),
                     os.path.basename(sys.argv[0]), out))
    for f in failed:
        print('  FAILED: %s' % f)