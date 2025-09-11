#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_for_plot_py38.py  (Python 3.8 compatible)
--------------------------------------------------
- 读取 b1..b8.nml -> 指向的 .out 文件
- 解析 reading_guide.nml
- 返回适合后续画图/编程的数据结构（dict），不写 JSON
- 提供 matplotlib 画图辅助函数（单图、跨 run）

依赖：matplotlib
"""

import re
from dataclasses import dataclass, asdict
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Any, Tuple, Optional

# -------------------- 数据结构 --------------------

@dataclass
class Basics:
    n_obs: int
    total_width: int
    n_bins: int
    n_meas_per_bin: int
    field_width: int
    n_cores: int
    MPI_nblock: int
    MPI_one_block: int
    output_format: str

@dataclass
class ObsDef:
    name: str
    kind: int      # 1 = scalar, 2 = field array
    width: int
    offset_lo: int
    offset_hi: int

# -------------------- 解析 --------------------

def parse_namelist_guide(path: Path) -> Tuple[Basics, List[ObsDef]]:
    txt = path.read_text(errors="ignore")
    basics = {}  # type: Dict[str, Any]
    obs_list = []  # type: List[ObsDef]
    cur_block = None  # type: Optional[str]
    cur_obs = None  # type: Optional[Dict[str, Any]]

    for line in txt.splitlines():
        s = line.strip()
        if not s or s.startswith("!"):
            continue
        if s.startswith("&basics"):
            cur_block = "basics"
            continue
        if s.startswith("&obs"):
            cur_block = "obs"
            cur_obs = {}
            continue
        if s == "/":
            if cur_block == "obs" and isinstance(cur_obs, dict):
                name = str(cur_obs.get("name"))
                kind = int(cur_obs.get("kind", 1))
                width = int(cur_obs.get("width", 1))
                olo = int(cur_obs.get("offset_lo", 0))
                ohi = int(cur_obs.get("offset_hi", 0))
                obs_list.append(ObsDef(name, kind, width, olo, ohi))
            cur_block = None
            cur_obs = None
            continue

        if cur_block == "basics":
            if "=" in s:
                k, v = s.split("=", 1)
                k = k.strip()
                v = v.strip().strip(",").strip('"').strip("'")
                basics[k] = v
        elif cur_block == "obs" and cur_obs is not None:
            if "=" in s:
                k, v = s.split("=", 1)
                k = k.strip()
                v = v.strip().strip(",").strip('"').strip("'")
                try:
                    cur_obs[k] = int(v)
                except ValueError:
                    cur_obs[k] = v

    def _int(k: str) -> int:
        try:
            return int(basics[k])
        except Exception:
            return 0

    b = Basics(
        n_obs=_int("n_obs"),
        total_width=_int("total_width"),
        n_bins=_int("n_bins"),
        n_meas_per_bin=_int("n_meas_per_bin"),
        field_width=_int("field_width"),
        n_cores=_int("n_cores"),
        MPI_nblock=_int("MPI_nblock"),
        MPI_one_block=_int("MPI_one_block"),
        output_format=str(basics.get("output_format", "")).strip(),
    )
    return b, obs_list

def parse_pointer_nml(nml_path: Path) -> Optional[Path]:
    txt = nml_path.read_text(errors="ignore")
    m = re.search(r'filename\s*=\s*"([^"]+)"', txt)
    if not m:
        return None
    fname = m.group(1)
    p = Path(fname)
    candidates = []  # type: List[Path]
    if p.is_absolute():
        candidates.append(p)
    else:
        candidates.append(nml_path.parent / p)
        candidates.append(nml_path.parent / p.name)
        candidates.append(Path.cwd() / p.name)
    for c in candidates:
        if c.exists():
            return c.resolve()
    return None

def parse_out_block_format(path: Path) -> Dict[str, List[Dict[str, Any]]]:
    name_re = re.compile(r'^([A-Za-z0-9_]+)\s+(\d+)\s+([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*$')
    out = defaultdict(list)  # type: Dict[str, List[Dict[str, Any]]]
    for ln in path.read_text(errors="ignore").splitlines():
        m = name_re.match(ln)
        if m:
            nm = m.group(1)
            idx = int(m.group(2))
            val = float(m.group(3))
            err = float(m.group(4))
            out[nm].append({'index': idx, 'value': val, 'stderr': err})
    return out

def assemble_dict(basics: Basics, obs_defs: List[ObsDef], sections: Dict[str, List[Dict[str, Any]]]) -> Dict[str, Any]:
    result = {
        'basics': asdict(basics),
        'observables': {}
    }  # type: Dict[str, Any]
    for od in obs_defs:
        series = sections.get(od.name, [])
        if od.kind == 1:
            val = series[0]['value'] if series else None
            err = series[0]['stderr'] if series else None
            result['observables'][od.name] = {'kind': 1, 'value': val, 'stderr': err}
        else:
            result['observables'][od.name] = {
                'kind': 2,
                'indices': [x['index'] for x in series],
                'values':  [x['value'] for x in series],
                'stderr':  [x['stderr'] for x in series],
            }
    return result

def load_runs(base_dir: Path) -> Dict[str, Dict[str, Any]]:
    base = Path(base_dir)
    guide = base / "reading_guide.nml"
    if not guide.exists():
        raise FileNotFoundError("reading_guide.nml not found in {}".format(base))
    basics, obs_defs = parse_namelist_guide(guide)
    results = {}  # type: Dict[str, Dict[str, Any]]
    for i in range(1, 9):
        nml = base / "b{}.nml".format(i)
        if not nml.exists():
            continue
        target = parse_pointer_nml(nml)
        if not target or not target.exists():
            print("[warn] {}: cannot resolve data file".format(nml.name))
            continue
        sections = parse_out_block_format(target)
        results[nml.stem] = assemble_dict(basics, obs_defs, sections)
    return results

# -------------------- 画图辅助 --------------------

def list_observables(data: Dict[str, Dict[str, Any]], run: str, kind: Optional[int]=None) -> List[str]:
    names = []  # type: List[str]
    obs = data[run]['observables']
    for k, v in obs.items():
        if kind is None or v.get('kind') == kind:
            names.append(k)
    return names

def plot_field(data: Dict[str, Dict[str, Any]], run: str, name: str):
    import matplotlib.pyplot as plt
    obs = data[run]['observables'][name]
    if obs['kind'] != 2:
        raise ValueError("{} is not a field (kind=2)".format(name))
    x = obs['indices']
    y = obs['values']
    yerr = obs['stderr']
    plt.figure()
    plt.errorbar(x, y, yerr=yerr, fmt='o-', capsize=3)
    plt.xlabel('index')
    plt.ylabel(name)
    plt.title('{}: {}'.format(run, name))
    plt.tight_layout()
    plt.show()
def plot_scalar_bar(data: Dict[str, Dict[str, Any]], run: str, names: Optional[List[str]]=None, rotate: int=45):
    import matplotlib.pyplot as plt
    obs = data[run]['observables']
    if names is None:
        names = [k for k, v in obs.items() if v.get('kind') == 1]
    vals = [obs[k]['value'] for k in names]
    err  = [obs[k]['stderr'] for k in names]
    plt.figure()
    xs = list(range(len(names)))
    plt.bar(xs, vals, yerr=err)
    plt.xticks(xs, names, rotation=rotate, ha='right')
    plt.ylabel('value')
    plt.title('{}: scalar observables'.format(run))
    plt.tight_layout()
    plt.show()
    
# -------------------- 跨 run 画图 --------------------

def collect_series_across_runs(data: Dict[str, Dict[str, Any]],
                               name: str,
                               field_index: Optional[int]=None,
                               runs: Optional[List[str]]=None):
    if runs is None:
        runs = sorted([r for r in data.keys() if r.startswith('b')],
                      key=lambda x: (len(x), x))
    runs_used, yvals, yerrs = [], [], []  
    # type: List[str], List[float], List[float]
    for r in runs:
        obs = data.get(r, {}).get('observables', {})
        if name not in obs:
            continue
        item = obs[name]
        k = item.get('kind')
        if k == 1:
            y = item.get('value')
            e = item.get('stderr')
        elif k == 2:
            if field_index is None:
                raise ValueError("{} 是场(kind=2)，需要指定 field_index (1-based).".format(name))
            idxs = item.get('indices', [])
            try:
                pos = idxs.index(field_index)
            except ValueError:
                continue
            y = item['values'][pos]
            e = item['stderr'][pos]
        else:
            continue
        if y is None:
            continue
        runs_used.append(r)
        yvals.append(y)
        yerrs.append(e)
    return runs_used, yvals, yerrs

def plot_across_runs(data: Dict[str, Dict[str, Any]],
                     name: str,
                     field_index: Optional[int]=None,
                     runs: Optional[List[str]]=None,
                     x: Optional[List[float]]=None,
                     xtick_labels: Optional[List[Any]]=None,
                     xlabel: str='x',
                     ylabel: Optional[str]=None,
                     title: Optional[str]=None):
    import matplotlib.pyplot as plt
    runs_used, yvals, yerrs = collect_series_across_runs(data, name, field_index, runs)
    if not runs_used:
        raise ValueError("没有可用数据。")
    if x is not None and len(x) != len(runs_used):
        raise ValueError("x 的长度必须与数据点数量一致。")
    xs = x if x is not None else list(range(len(runs_used)))
    if xtick_labels is None:
        xtick_labels = runs_used
    elif len(xtick_labels) != len(runs_used):
        raise ValueError("xtick_labels 的长度必须与数据点数量一致。")
    plt.figure()
    plt.errorbar(xs, yvals, yerr=yerrs, fmt='o-', capsize=3)
    plt.xticks(xs, [str(t) for t in xtick_labels])
    plt.xlabel(xlabel)
    if ylabel is None:
        ylabel = "{}".format(name) if field_index is None else "{}[{}]".format(name, field_index)
    plt.ylabel(ylabel)
    if title is None:
        title = "{}".format(name) if field_index is None else "{} @ index {}".format(name, field_index)
    plt.title(title)
    
    plt.tight_layout()
    
from typing import List, Optional, Any

def plot_g2_overlay(data: dict,
                    runs: Optional[List[str]] = None,
                    name_substr: str = "G2",
                    with_error: bool = True,
                    xlabel: str = "index",
                    ylabel: Optional[str] = None,
                    title: Optional[str] = None,
                    loglog: bool = False,                 # ← 新增：同时对 x/y 用对数轴
                    drop_non_positive: bool = True):      # ← 新增：丢弃 ≤0 的点，便于 log 轴
    """
    在同一张图上叠加绘制 b1..b4（存在即画）中，名字包含 name_substr（默认 'G2'）
    的场型观测量：x=index，y=value。支持 log-log。
    """
    import matplotlib.pyplot as plt

    if runs is None:
        runs = [r for r in ["b1", "b2", "b3", "b4"] if r in data]
    if not runs:
        raise ValueError("没有找到可用的 b1..b4 运行。")

    fig, ax = plt.subplots()
    used_runs = []
    chosen_name = None

    for r in runs:
        obs_map = data.get(r, {}).get("observables", {})
        if "corf_G2" in obs_map and obs_map["corf_G2"].get("kind") == 2:
            name = "corf_G2"
        else:
            candidates = [k for k, v in obs_map.items()
                          if v.get("kind") == 2 and (name_substr in k)]
            if not candidates:
                continue
            candidates.sort()
            name = candidates[0]

        ob = obs_map[name]
        x = list(ob.get("indices", []))
        x = x[1:20]
        y = list(ob.get("values", []))
        y = y[1:20]
        yerr = list(ob.get("stderr", [])) if ob.get("stderr") is not None else None
        yerr = yerr[1:20] 
        if not x or not y:
            continue

        # log 轴需过滤非正值
        if loglog and drop_non_positive:
            mask = [ (xi > 0) and (yi > 0) for xi, yi in zip(x, y) ]
            x = [xi for xi, m in zip(x, mask) if m]
            y = [yi for yi, m in zip(y, mask) if m]
            if yerr is not None:
                yerr = [ei for ei, m in zip(yerr, mask) if m]
            if not x:
                continue

        if with_error and yerr is not None and len(yerr) == len(y):
            ax.errorbar(x, y, yerr=yerr, fmt='o-', capsize=3, label=r)
        else:
            ax.plot(x, y, 'o-', label=r)

        used_runs.append(r)
        chosen_name = name

    if not used_runs:
        raise ValueError("在所给 runs 中没有找到包含 'G2' 的场型观测量。")

    if ylabel is None:
        ylabel = chosen_name or "G2"

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title is None:
        title = "{} overlay: {}".format(ylabel, ", ".join(used_runs))
    ax.set_title(title)
    ax.legend()

    if loglog:
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.savefig("g2loglog.png", dpi=200, bbox_inches="tight")
    else:
        plt.savefig("g2.png", dpi=200, bbox_inches="tight")
    fig.tight_layout()
    
    return ax
    
# -------------------- 入口（可选） --------------------

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("dir", nargs="?", default=".")
    ap.add_argument("--run", default=None)
    ap.add_argument("--field", default=None)
    ap.add_argument("--scalars", action="store_true")
    args = ap.parse_args()

    data = load_runs(Path(args.dir))
    if not data:
        print("No runs found.")
        return

    run = args.run or sorted(data.keys())[0]
    if args.field is None:
        fields = list_observables(data, run, kind=2)
        if fields:
            args.field = fields[0]

    if args.field:
        plot_field(data, run=run, name=args.field)
    if args.scalars:
        plot_scalar_bar(data, run=run)
    import matplotlib.pyplot as plt
    plt.show()
if __name__ == "__main__":
    main()
