"""
PyMOL 通用共价 warhead 位点筛选脚本 v3.1.3
==========================================

用途
----
本脚本用于重复分析不同蛋白复合物：先在 target 界面寻找可由指定 warhead
化学捕获的亲核残基，再评估 binder 上哪些残基适合突变为 Cys、安装 linker/warhead。
它是结构初筛工具，不预测反应速率、共价收率或体内效果。

加载
----
在 PyMOL 命令行运行：

    run /path/to/warhead_site_rank.py
    warhead_help
    warhead_scan_prompt              # 打开 PyMOL Qt 参数窗口
    warhead_scan_prompt advanced=1   # 参数窗口显示全部高级距离参数

直接运行示例：

    warhead_scan object_name=complex, target_chain=B, binder_chain=A, \
        profile=fluorosulfate, interface_cutoff=5, candidate_cutoff=12, \
        ideal_min=5, ideal_max=10, label_top=3

也支持任意 PyMOL selection；显式 selection 优先于 chain：

    warhead_scan object_name=complex, \
        target_sel="chain B+D and resi 1-200", \
        binder_sel="chain A+C", profile=fluorosulfate

力场/非标准残基名可用内置别名或显式映射：

    warhead_scan ... target_types=CYX+HIE
    warhead_scan ... resn_aliases=CSO:CYS+ABC:ALA

批量模式：

    warhead_batch manifest=/path/jobs.csv, output_dir=/path/results

GUI 交互说明
------------
`warhead_scan_prompt` 会打开 PyMOL 自带的 Qt 参数窗口，不再调用 Python
`input()`。这是因为常见的 Windows/macOS PyMOL 图形界面没有可供脚本读取的
标准输入流。按 Cancel 不会启动分析；如果 Qt 不可用，控制台会打印一条可直接
复制的 `warhead_scan ...` 命令。

v3.1.3 相对 v3.1.2
------------------
- Binder 上已有 Cys 必须先位于结合界面（残基重原子到 target 蛋白 ≤ `interface_cutoff`），
  才会枚举 native+χ1 rotamer 并进入排名。不再对核内/远端 Cys 做构象展开。
- 界面判定用 CB 或 SG 到界面亲核残基的距离做 reach 预筛，避免 native SG 背朝
  target 时把真正的界面 Cys 漏掉。

v3.1.2 相对 v3.1.1
------------------
- 已有 Cys 同时保留 native SG，并枚举三个标准 χ1；与 native 重合的 rotamer 会去重。
- clash、直线路径和复合物暴露度计入同一对象中的有机配体/金属/无机离子（不含溶剂）。
- `best_target` / `target_options` 改为 `B/Tyr134 OH` 这种可读写法。

v3.1.1 相对 v3.1.0
------------------
- N/C 端与肽段断裂端默认 class cap 到 C，避免端基因几何碰巧合适而排进第一轮。
- `tighter_than_preferred` / `longer_than_preferred` 改为用建模 SG（或实际用于
  评分的距离），不再误用突变前锚点距离。
- `class_cap` / `downgrade_reasons` 仅在等级确实被压低时填写。
- 指定的 chain ID 不存在时列出对象里真实可用的链。
- 接触/极性/盐桥计数改用空间网格，判定阈值不变。

v3.1.0 相对 v3.0.1 的轻量增强
----------------------------
1. 残基别名：内置 HIE/HID/CYX/MSE 等映射与 NALA/CGLY 端基规则；化学已改变的
   残基（CSO/SEC/PTR/MLY 等）不静默映射；控制台报告未知残基；支持
   `resn_aliases=`。评分用 canonical 名，CSV 同时保留 `*_resn` 与 `*_resn_std`。
2. 更稳的对象解析：`get_object_list` 返回 None/异常时不再崩溃，并提示 PyMOL
   可能改写对象名。
3. 更可读的 CSV：`A/Ser74`、`A:S74C`、`target_site`/`binder_site`/`mutation`。
   PyMOL 标签保持短格式：`#1 A/Ser74  8.55Å` 与 `B/Tyr134`。

核心原则
--------
1. 5 Å 不是固定阈值。`interface_cutoff`、`candidate_cutoff`、优选距离、接触、
   极性邻近、盐桥、内部接触和 clash 阈值均为独立参数。
2. Target 使用实际反应原子：例如 Tyr-OH、Lys-NZ、His-ND1/NE2、Cys-SG。
3. Binder 不再只用原始 CB 距离排序。脚本构建三个标准 Cys rotamer
   (chi1 = -60, +60, 180 degrees)，评估假想 SG 的距离、方向 orientation、
   clash、暴露度和直线路径阻挡。
4. 原侧链在 SG 暴露度计算前会截短到 Cys 应保留的骨架/CB，因此 Leu、Ile、
   Lys 等大侧链不会因为“自己遮住自己的 CB”而被自动误判。
5. Binder 固有界面特征只计算一次，再用于多个 target pair，适合大界面和批量任务。
6. 输出同时给出 RawScore、UncappedClass、FinalClass、ClassCap 和
   DowngradeReasons；分数高但因结构风险降级时不会显得像排序错误。
7. PyMOL 画面最多标注三个唯一 binder 位点。binder 标签短，例如
   ``#1 A/Ser74  8.55Å``；target 单独标 ``B/Tyr134``；等级用颜色区分
   （A 绿 / B 青 / C 品红 / D 灰）。其他候选只打印并写入 CSV。
8. `label_top=0`、只有 target 命中但没有 binder candidate、完全无命中均可安全运行。

内置 chemistry profile
----------------------
- `fluorosulfate`: TYR:OH, LYS:NZ, HIS:ND1/NE2。
- `broad_nucleophile`: CYS:SG, LYS:NZ, TYR:OH, HIS:ND1/NE2。
- `cys_target`: CYS:SG。
- `custom`: 通过 `reactive_atoms="CYS:SG;LYS:NZ;HIS:ND1|NE2"` 定义。

不同 warhead 对 Cys/Lys/Tyr/His 的反应性不同。不要因为结构上扫描了四类残基，
就假定一种 warhead 能同等捕获四类侧链。

两个主要距离
------------
- `interface_cutoff`: target 反应原子到 binder 任意重原子的最大距离，用于定义
  target 亲核残基是否位于界面。
- `candidate_cutoff`: 模拟 Cys-SG 到 target 反应原子的最大候选距离，用于定义
  linker/warhead 安装位点搜索范围。

`ideal_min`/`ideal_max` 是当前 linker/warhead 家族的优选 SG-to-target 距离，
不是通用化学定律。

排名解释
--------
- A：优先人工复核和第一轮实验。
- B：值得进入第一轮。
- C：备选；存在明确风险或方向/几何信息不足。
- D：通常不建议第一轮使用。

等级先由 raw score 得到，再应用明确 class cap。例如 Pro->Cys、二硫键 Cys、
所有 Cys rotamer 均 clash、极低 SG 暴露度可触发降级。每个降级原因都写入
`downgrade_reasons`。

重要限制
--------
Cys rotamer 和方向性显著优于单纯 CB 距离，但脚本仍未枚举完整
maleimide-linker-warhead 构象。前三名必须继续人工检查：完整 linker 长度与柔性、
warhead 朝向、界面水网络、target 位点质子化、binder 热点、WT/Cys mutant/偶联后
亲和力以及质谱确认。

输出
----
每次扫描写出：
- `<prefix>_targets.csv`
- `<prefix>_pairs.csv`
- `<prefix>_binder_rank.csv`
- `<prefix>_metadata.json`

批量 manifest 常用列：
`job_id,structure,object_name,state,target_chain,binder_chain,target_sel,binder_sel,profile,
interface_cutoff,candidate_cutoff,ideal_min,ideal_max,resn_aliases,output_prefix`。
"""

from dataclasses import dataclass
from collections import defaultdict
import csv
import json
import math
import os
import re
import traceback

from pymol import cmd


HELP_TEXT = (__doc__ or "").strip()
__version__ = "3.1.3"


@dataclass(frozen=True)
class WarheadProfile:
    name: str
    reactive_atoms: dict
    description: str = ""


@dataclass
class ScanConfig:
    interface_cutoff: float = 5.0
    candidate_cutoff: float = 12.0
    ideal_min: float = 5.0
    ideal_max: float = 10.0
    contact_cutoff: float = 4.0
    polar_cutoff: float = 3.5
    salt_bridge_cutoff: float = 4.0
    internal_contact_cutoff: float = 4.5
    clash_cutoff: float = 2.0
    path_radius: float = 1.6
    include_disulfide_targets: bool = False
    compute_rotamers: bool = True

    def validate(self):
        positive = {
            "interface_cutoff": self.interface_cutoff,
            "candidate_cutoff": self.candidate_cutoff,
            "ideal_min": self.ideal_min,
            "ideal_max": self.ideal_max,
            "contact_cutoff": self.contact_cutoff,
            "polar_cutoff": self.polar_cutoff,
            "salt_bridge_cutoff": self.salt_bridge_cutoff,
            "internal_contact_cutoff": self.internal_contact_cutoff,
            "clash_cutoff": self.clash_cutoff,
            "path_radius": self.path_radius,
        }
        for name, value in positive.items():
            if float(value) <= 0:
                raise ValueError("%s must be positive" % name)
        if float(self.ideal_max) < float(self.ideal_min):
            raise ValueError("Require ideal_max >= ideal_min")
        return self

    def to_dict(self):
        return {
            "interface_cutoff": float(self.interface_cutoff),
            "candidate_cutoff": float(self.candidate_cutoff),
            "ideal_min": float(self.ideal_min),
            "ideal_max": float(self.ideal_max),
            "contact_cutoff": float(self.contact_cutoff),
            "polar_cutoff": float(self.polar_cutoff),
            "salt_bridge_cutoff": float(self.salt_bridge_cutoff),
            "internal_contact_cutoff": float(self.internal_contact_cutoff),
            "clash_cutoff": float(self.clash_cutoff),
            "path_radius": float(self.path_radius),
            "include_disulfide_targets": bool(self.include_disulfide_targets),
            "compute_rotamers": bool(self.compute_rotamers),
        }


PROFILES = {
    "fluorosulfate": WarheadProfile(
        "fluorosulfate",
        {"TYR": ("OH",), "LYS": ("NZ",), "HIS": ("ND1", "NE2")},
        "Aryl fluorosulfate/SuFEx-oriented screening.",
    ),
    "broad_nucleophile": WarheadProfile(
        "broad_nucleophile",
        {
            "CYS": ("SG",), "LYS": ("NZ",), "TYR": ("OH",),
            "HIS": ("ND1", "NE2"),
        },
        "Broad structural inventory; chemical compatibility must be reviewed.",
    ),
    "cys_target": WarheadProfile(
        "cys_target", {"CYS": ("SG",)}, "Cys-targeted electrophile screening."
    ),
}


MUTABILITY = {
    "ALA": 3.0, "SER": 2.5, "THR": 1.8, "ASN": 1.5, "GLN": 1.2,
    "VAL": 0.5, "ILE": 0.0, "LEU": 0.0, "MET": 0.3,
    "ASP": -0.8, "GLU": -0.8, "LYS": -1.0, "ARG": -1.2,
    "HIS": -1.0, "PHE": -1.5, "TYR": -1.8, "TRP": -2.0,
    "GLY": -2.2, "PRO": -3.5, "CYS": 3.5,
}

MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0,
    "CYS": 167.0, "GLN": 225.0, "GLU": 223.0, "GLY": 104.0,
    "HIS": 224.0, "ILE": 197.0, "LEU": 201.0, "LYS": 236.0,
    "MET": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0,
    "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

AA_ONE_LETTER = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Only protonation/charge/isosteric variants map onto a standard parent.
# Chemically altered residues (CSO/SEC/PTR/MLY/...) stay unmapped by default.
RESN_ALIASES = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "CYX": "CYS", "CYM": "CYS",
    "LYN": "LYS", "ASH": "ASP", "GLH": "GLU", "AR0": "ARG",
    "MSE": "MET",
    "ILE": "ILE",
}
_TERMINAL_PREFIXES = ("N", "C")

BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}
MUTANT_RETAINED_NAMES = BACKBONE_NAMES | {"CB"}
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
CHARGED = {"ASP", "GLU", "LYS", "ARG"}
POLAR_ATOMS = {
    "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"}, "CYS": {"SG"},
    "ASN": {"OD1", "ND2"}, "GLN": {"OE1", "NE2"},
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"},
    "LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"},
    "HIS": {"ND1", "NE2"}, "TRP": {"NE1"},
}
POSITIVE_ATOMS = {"LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}}
NEGATIVE_ATOMS = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
CLASS_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3}
VDW_RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
CYS_CHI1 = (-60.0, 60.0, 180.0)
CYS_CB_SG_LENGTH = 1.81
CYS_CA_CB_SG_ANGLE = 114.4
NATIVE_ROTAMER_MATCH_A = 0.35


# ---------------------------------------------------------------------------
# Generic parsing and geometry helpers
# ---------------------------------------------------------------------------


def _as_bool(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).strip().lower() not in {"", "0", "false", "no", "off", "none"}


def _as_float(value, name):
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError("%s must be a number, got %r" % (name, value))


def _as_int(value, name):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        raise ValueError("%s must be an integer, got %r" % (name, value))


def _blank(value):
    return value is None or str(value).strip().lower() in {"", "none", "null", "auto"}


def _canonical_resn(resn, extra_aliases=None):
    """Map a residue name onto its standard 3-letter parent for scoring.

    The deposited name must still be used for PyMOL selections and CSV
    ``*_resn`` columns; callers should keep both when writing outputs.
    """
    name = str(resn or "").upper().strip()
    if extra_aliases:
        mapped = extra_aliases.get(name)
        if mapped:
            return str(mapped).upper().strip()
    if name in MUTABILITY:
        return name
    mapped = RESN_ALIASES.get(name)
    if mapped:
        return mapped
    if len(name) == 4 and name[0] in _TERMINAL_PREFIXES:
        stem = name[1:]
        if stem in MUTABILITY:
            return stem
        stem_alias = RESN_ALIASES.get(stem)
        if stem_alias:
            return stem_alias
    return name


def _is_known_resn(resn, extra_aliases=None):
    return _canonical_resn(resn, extra_aliases) in MUTABILITY


def _parse_alias_option(value):
    """Parse ``MSE:MET+HIE:HIS`` style user aliases into a dict."""
    if not value:
        return {}
    if isinstance(value, dict):
        return {
            str(key).upper().strip(): str(val).upper().strip()
            for key, val in value.items()
        }
    out = {}
    for item in re.split(r"[,+;\s]+", str(value).upper()):
        if not item:
            continue
        if ":" not in item:
            raise ValueError(
                "Alias '%s' must be written SOURCE:TARGET, e.g. MSE:MET" % item
            )
        src, dst = item.split(":", 1)
        src, dst = src.strip(), dst.strip()
        if dst not in MUTABILITY:
            raise ValueError(
                "Alias target '%s' is not a standard residue name" % dst
            )
        out[src] = dst
    return out


def _unknown_residue_summary(residue_groups, aliases=None):
    """Return ``{resn: count}`` for residues the scoring tables do not know."""
    unknown = defaultdict(int)
    for key in residue_groups:
        resn = str(key[2] or "").upper().strip()
        if not _is_known_resn(resn, aliases):
            unknown[resn] += 1
    return dict(unknown)


def _residue_title(resn, resi):
    """Readable residue tag such as ``His81`` or ``Ser74``."""
    text = str(resn or "").strip()
    pretty = text[:1].upper() + text[1:].lower() if text else "Unk"
    return "%s%s" % (pretty, resi)


def _distance(a, b):
    return math.sqrt(sum((float(a[i]) - float(b[i])) ** 2 for i in range(3)))


def _d(a, b):
    return _distance(a, b)


def _v_add(a, b):
    return tuple(float(a[i]) + float(b[i]) for i in range(3))


def _v_sub(a, b):
    return tuple(float(a[i]) - float(b[i]) for i in range(3))


def _v_scale(a, s):
    return tuple(float(v) * float(s) for v in a)


def _dot(a, b):
    return sum(float(a[i]) * float(b[i]) for i in range(3))


def _cross(a, b):
    return (
        float(a[1]) * float(b[2]) - float(a[2]) * float(b[1]),
        float(a[2]) * float(b[0]) - float(a[0]) * float(b[2]),
        float(a[0]) * float(b[1]) - float(a[1]) * float(b[0]),
    )


def _norm(a):
    return math.sqrt(_dot(a, a))


def _unit(a):
    n = _norm(a)
    if n < 1e-10:
        raise ValueError("Cannot normalize a zero-length vector")
    return _v_scale(a, 1.0 / n)


def _safe_unit(a, fallback=(1.0, 0.0, 0.0)):
    try:
        return _unit(a)
    except ValueError:
        return tuple(fallback)


def _natural_residue_key(resi):
    text = str(resi or "").strip()
    match = re.match(r"^(-?\d+)(.*)$", text)
    if match:
        number = int(match.group(1))
        insertion = match.group(2).strip()
        # Empty insertion sorts before A/B insertion codes.
        return (0, number, 0 if insertion == "" else 1, insertion.upper(), text)
    return (1, 0, 1, text.upper(), text)


def _residue_sort_key(key):
    chain, resi, resn = key
    return (str(chain), _natural_residue_key(resi), str(resn))


def _atom_element(atom):
    elem = str(getattr(atom, "elem", "") or "").strip().upper()
    if elem:
        return elem
    name = str(getattr(atom, "name", "") or "").strip()
    return name[:1].upper()


def _atom_identity(atom):
    return (
        str(getattr(atom, "model", "")), str(getattr(atom, "segi", "")),
        str(getattr(atom, "chain", "")), str(getattr(atom, "resi", "")),
        str(getattr(atom, "resn", "")), str(getattr(atom, "name", "")),
    )


def _preferred_alt_atoms(atoms):
    chosen = {}
    for atom in atoms:
        key = _atom_identity(atom)
        alt = str(getattr(atom, "alt", "") or "").strip()
        occ = float(getattr(atom, "q", 0.0) or 0.0)
        priority = 3 if alt == "" else (2 if alt == "A" else 1)
        old = chosen.get(key)
        if old is None or (priority, occ) > old[0]:
            chosen[key] = ((priority, occ), atom)
    return [item[1] for item in chosen.values()]


def _atoms(selection, state=1):
    return _preferred_alt_atoms(list(cmd.get_model(selection, int(state)).atom))


def _reskey(atom):
    return (
        str(getattr(atom, "chain", "")), str(getattr(atom, "resi", "")),
        str(getattr(atom, "resn", "")).upper(),
    )


def _group_by_residue(atoms):
    grouped = defaultdict(list)
    for atom in atoms:
        grouped[_reskey(atom)].append(atom)
    return dict(grouped)


def _heavy(atoms):
    return [a for a in atoms if _atom_element(a) != "H"]


def _sidechain(atoms):
    return [a for a in _heavy(atoms) if str(getattr(a, "name", "")) not in BACKBONE_NAMES]


def _find_atom(atoms, name):
    return next((a for a in atoms if str(getattr(a, "name", "")) == name), None)


def _cell_index(coord, inv):
    return (
        int(math.floor(float(coord[0]) * inv)),
        int(math.floor(float(coord[1]) * inv)),
        int(math.floor(float(coord[2]) * inv)),
    )


def _spatial_grid(atoms, cell):
    inv = 1.0 / float(cell)
    grid = defaultdict(list)
    for atom in atoms:
        grid[_cell_index(atom.coord, inv)].append(atom)
    return grid, inv


def _iter_nearby_atoms(atom, grid, inv):
    ix, iy, iz = _cell_index(atom.coord, inv)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                bucket = grid.get((ix + dx, iy + dy, iz + dz))
                if not bucket:
                    continue
                for other in bucket:
                    yield other


def _min_distance(atoms1, atoms2):
    best = float("inf")
    pair = (None, None)
    for atom1 in atoms1:
        for atom2 in atoms2:
            value = _distance(atom1.coord, atom2.coord)
            if value < best:
                best = value
                pair = (atom1, atom2)
    return best, pair


def _is_polar_atom(atom, aliases=None):
    name = str(getattr(atom, "name", ""))
    resn = _canonical_resn(getattr(atom, "resn", ""), aliases)
    return name in {"N", "O", "OXT"} or name in POLAR_ATOMS.get(resn, set())


def _is_positive_atom(atom, aliases=None):
    return str(getattr(atom, "name", "")) in POSITIVE_ATOMS.get(
        _canonical_resn(getattr(atom, "resn", ""), aliases), set()
    )


def _is_negative_atom(atom, aliases=None):
    return str(getattr(atom, "name", "")) in NEGATIVE_ATOMS.get(
        _canonical_resn(getattr(atom, "resn", ""), aliases), set()
    )


def _contact_pairs(atoms1, atoms2, cutoff):
    cutoff = float(cutoff)
    left = _heavy(atoms1)
    right = _heavy(atoms2)
    if not left or not right:
        return 0
    if len(left) * len(right) < 128:
        count = 0
        for atom1 in left:
            for atom2 in right:
                if _distance(atom1.coord, atom2.coord) <= cutoff:
                    count += 1
        return count
    grid, inv = _spatial_grid(right, cutoff)
    count = 0
    for atom1 in left:
        for atom2 in _iter_nearby_atoms(atom1, grid, inv):
            if _distance(atom1.coord, atom2.coord) <= cutoff:
                count += 1
    return count


def _polar_contact_pairs(atoms1, atoms2, cutoff, aliases=None):
    cutoff = float(cutoff)
    left = [atom for atom in _heavy(atoms1) if _is_polar_atom(atom, aliases)]
    right = [atom for atom in _heavy(atoms2) if _is_polar_atom(atom, aliases)]
    if not left or not right:
        return 0
    if len(left) * len(right) < 128:
        count = 0
        for atom1 in left:
            for atom2 in right:
                if _distance(atom1.coord, atom2.coord) <= cutoff:
                    count += 1
        return count
    grid, inv = _spatial_grid(right, cutoff)
    count = 0
    for atom1 in left:
        for atom2 in _iter_nearby_atoms(atom1, grid, inv):
            if _distance(atom1.coord, atom2.coord) <= cutoff:
                count += 1
    return count


def _polar_min_distance(atoms1, atoms2, aliases=None):
    best = None
    for atom1 in _heavy(atoms1):
        if not _is_polar_atom(atom1, aliases):
            continue
        for atom2 in _heavy(atoms2):
            if not _is_polar_atom(atom2, aliases):
                continue
            value = _distance(atom1.coord, atom2.coord)
            if best is None or value < best:
                best = value
    return best


def _salt_bridge_pairs(atoms1, atoms2, cutoff, aliases=None):
    cutoff = float(cutoff)
    left = _heavy(atoms1)
    right = _heavy(atoms2)
    if not left or not right:
        return 0
    grid, inv = _spatial_grid(right, cutoff) if len(left) * len(right) >= 128 else (None, None)
    count = 0
    for atom1 in left:
        neighbors = right if grid is None else _iter_nearby_atoms(atom1, grid, inv)
        for atom2 in neighbors:
            opposite = (
                _is_positive_atom(atom1, aliases) and _is_negative_atom(atom2, aliases)
            ) or (
                _is_negative_atom(atom1, aliases) and _is_positive_atom(atom2, aliases)
            )
            if opposite and _distance(atom1.coord, atom2.coord) <= cutoff:
                count += 1
    return count


def _unique_contact_residues(probe_atoms, residue_groups, cutoff):
    labels = []
    if not probe_atoms:
        return 0, labels
    for key, atoms in residue_groups.items():
        value, _ = _min_distance(_heavy(probe_atoms), _heavy(atoms))
        if value <= float(cutoff):
            labels.append("%s:%s%s" % (key[0], key[2], key[1]))
    return len(labels), labels


def _parse_reactive_atoms_spec(value):
    result = {}
    for item in re.split(r"[;]+", str(value or "")):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError("Reactive atom item must look like TYR:OH: %s" % item)
        resn, atom_text = item.split(":", 1)
        atoms = tuple(x.strip().upper() for x in re.split(r"[|,+/]", atom_text) if x.strip())
        if not atoms:
            raise ValueError("No reactive atom defined for %s" % resn)
        result[resn.strip().upper()] = atoms
    if not result:
        raise ValueError("At least one custom reactive atom definition is required")
    return result


def _parse_target_types(value, aliases=None):
    if _blank(value):
        return None
    if isinstance(value, (list, tuple, set)):
        parts = [str(v).upper() for v in value]
    else:
        parts = [p for p in re.split(r"[,+;/\s]+", str(value).upper()) if p]
    return set(_canonical_resn(part, aliases) for part in parts)


def _resolve_profile(profile="broad_nucleophile", reactive_atoms="", target_types=None, aliases=None):
    if isinstance(profile, WarheadProfile):
        base = profile
    else:
        name = str(profile or "broad_nucleophile").strip().lower()
        if name == "custom":
            raw = _parse_reactive_atoms_spec(reactive_atoms)
            base = WarheadProfile(
                "custom",
                {_canonical_resn(key, aliases): atoms for key, atoms in raw.items()},
                "User-defined profile",
            )
        elif name in PROFILES:
            base = PROFILES[name]
        else:
            raise ValueError("Unknown profile '%s'. Available: %s, custom" % (
                name, ", ".join(sorted(PROFILES))
            ))
    requested = _parse_target_types(target_types, aliases=aliases)
    mapping = dict(base.reactive_atoms)
    if requested is not None:
        unknown = sorted(requested.difference(mapping))
        if unknown:
            raise ValueError("Target type(s) not in selected profile: %s" % ", ".join(unknown))
        mapping = {key: value for key, value in mapping.items() if key in requested}
    return WarheadProfile(base.name, mapping, base.description)


# ---------------------------------------------------------------------------
# Peptide topology and Cys rotamers
# ---------------------------------------------------------------------------


def _terminal_statuses(residue_groups, peptide_min=1.0, peptide_max=1.9):
    by_chain = defaultdict(list)
    for key in residue_groups:
        by_chain[key[0]].append(key)
    statuses = {key: set() for key in residue_groups}
    for chain, keys in by_chain.items():
        ordered = sorted(keys, key=_residue_sort_key)
        if not ordered:
            continue
        statuses[ordered[0]].add("n_terminus")
        statuses[ordered[-1]].add("c_terminus")
        for left, right in zip(ordered, ordered[1:]):
            c_atom = _find_atom(residue_groups[left], "C")
            n_atom = _find_atom(residue_groups[right], "N")
            connected = False
            if c_atom is not None and n_atom is not None:
                value = _distance(c_atom.coord, n_atom.coord)
                connected = float(peptide_min) <= value <= float(peptide_max)
            else:
                # Coordinate topology unavailable: only assume continuity for
                # adjacent natural residue numbers.
                lkey = _natural_residue_key(left[1])
                rkey = _natural_residue_key(right[1])
                connected = lkey[0] == 0 and rkey[0] == 0 and (rkey[1] - lkey[1] in {0, 1})
            if not connected:
                statuses[left].add("fragment_c_terminus")
                statuses[right].add("fragment_n_terminus")
    return {key: tuple(sorted(value)) for key, value in statuses.items()}


def _pseudo_cb_from_backbone(n, ca, c):
    # AlphaFold-style pseudo-CB construction; vectors are scaled by fixed
    # coefficients derived from ideal tetrahedral geometry.
    b = _v_sub(ca, n)
    cc = _v_sub(c, ca)
    a = _cross(b, cc)
    return _v_add(ca, _v_add(_v_scale(a, -0.58273431), _v_add(
        _v_scale(b, 0.56802827), _v_scale(cc, -0.54067466)
    )))


def _dihedral_degrees(a, b, c, d):
    b0 = _v_scale(_v_sub(b, a), -1.0)
    b1 = _v_sub(c, b)
    b2 = _v_sub(d, c)
    b1u = _safe_unit(b1)
    v = _v_sub(b0, _v_scale(b1u, _dot(b0, b1u)))
    w = _v_sub(b2, _v_scale(b1u, _dot(b2, b1u)))
    x = _dot(v, w)
    y = _dot(_cross(b1u, v), w)
    angle = math.degrees(math.atan2(y, x))
    # Normalize +180 to 180 rather than -180 for stable canonical labels.
    if angle <= -179.999999:
        return 180.0
    return angle


def _place_atom(a, b, c, length, angle_deg, dihedral_deg):
    # Place D from A-B-C internal coordinates: |C-D|, angle B-C-D,
    # dihedral A-B-C-D. The phase convention only affects which canonical
    # rotamer receives which sign; all three tetrahedral positions are kept.
    bc = _safe_unit(_v_sub(b, c))  # C -> B
    plane_normal = _cross(_v_sub(b, a), _v_sub(c, b))
    if _norm(plane_normal) < 1e-8:
        trial = (0.0, 0.0, 1.0) if abs(bc[2]) < 0.9 else (0.0, 1.0, 0.0)
        plane_normal = _cross(bc, trial)
    nvec = _unit(plane_normal)
    mvec = _unit(_cross(nvec, bc))
    theta = math.radians(float(angle_deg))
    phi = math.radians(float(dihedral_deg))
    radial = _v_add(_v_scale(mvec, math.cos(phi)), _v_scale(nvec, math.sin(phi)))
    direction = _v_add(_v_scale(bc, math.cos(theta)), _v_scale(radial, math.sin(theta)))
    return _v_add(c, _v_scale(_unit(direction), float(length)))


def _native_cys_chi1(res_atoms):
    n_atom = _find_atom(res_atoms, "N")
    ca_atom = _find_atom(res_atoms, "CA")
    cb_atom = _find_atom(res_atoms, "CB")
    sg_atom = _find_atom(res_atoms, "SG")
    if n_atom is None or ca_atom is None or cb_atom is None or sg_atom is None:
        return None
    return _dihedral_degrees(n_atom.coord, ca_atom.coord, cb_atom.coord, sg_atom.coord)


def _build_cys_rotamers(res_atoms, resn):
    resn = str(resn).upper()
    n_atom = _find_atom(res_atoms, "N")
    ca_atom = _find_atom(res_atoms, "CA")
    cb_atom = _find_atom(res_atoms, "CB")
    rotamers = []
    native_sg = None
    if resn == "CYS":
        sg_atom = _find_atom(res_atoms, "SG")
        if sg_atom is not None and cb_atom is not None:
            native_sg = tuple(sg_atom.coord)
            native_chi1 = _native_cys_chi1(res_atoms)
            rotamers.append({
                "chi1_deg": native_chi1,
                "chi1_label": "native",
                "sg_coord": native_sg,
                "cb_coord": tuple(cb_atom.coord),
                "source": "native_cys",
            })
    if n_atom is None or ca_atom is None:
        return rotamers
    if cb_atom is None:
        c_atom = _find_atom(res_atoms, "C")
        if c_atom is None:
            return rotamers
        cb_coord = _pseudo_cb_from_backbone(n_atom.coord, ca_atom.coord, c_atom.coord)
        source = "pseudo_cb"
    else:
        cb_coord = tuple(cb_atom.coord)
        source = "observed_cb"
    for chi1 in CYS_CHI1:
        sg = _place_atom(
            n_atom.coord, ca_atom.coord, cb_coord,
            CYS_CB_SG_LENGTH, CYS_CA_CB_SG_ANGLE, 180.0 - chi1,
        )
        if native_sg is not None and _distance(sg, native_sg) <= NATIVE_ROTAMER_MATCH_A:
            continue
        rotamers.append({
            "chi1_deg": float(chi1), "chi1_label": "%+g" % chi1,
            "sg_coord": tuple(sg), "cb_coord": tuple(cb_coord), "source": source,
        })
    return rotamers


def _vdw(atom_or_element):
    if isinstance(atom_or_element, str):
        elem = atom_or_element.upper()
    else:
        elem = _atom_element(atom_or_element)
    return VDW_RADII.get(elem, 1.70)


def _fibonacci_points(n=60):
    points = []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for index in range(int(n)):
        y = 1.0 - (index / float(max(n - 1, 1))) * 2.0
        radius = math.sqrt(max(0.0, 1.0 - y * y))
        theta = golden * index
        points.append((math.cos(theta) * radius, y, math.sin(theta) * radius))
    return points


_SPHERE_POINTS = _fibonacci_points(60)


def _truncated_environment(residue_groups, mutated_key, include_mutated_backbone=True):
    atoms = []
    for key, ratoms in residue_groups.items():
        if key != mutated_key:
            atoms.extend(_heavy(ratoms))
        elif include_mutated_backbone:
            atoms.extend([
                atom for atom in _heavy(ratoms)
                if str(getattr(atom, "name", "")) in MUTANT_RETAINED_NAMES
            ])
    return atoms


def _approx_accessible_fraction(center, environment_atoms, probe=1.4):
    shell = _vdw("S") + float(probe)
    # Only atoms whose expanded van der Waals sphere can intersect the SG
    # solvent shell can occlude a sample point. This spatial prefilter keeps
    # the pure-Python calculation practical for large interfaces and batches.
    local_atoms = [
        atom for atom in environment_atoms
        if _distance(center, atom.coord) <= shell + _vdw(atom) + float(probe)
    ]
    accessible = 0
    for direction in _SPHERE_POINTS:
        point = _v_add(center, _v_scale(direction, shell))
        blocked = False
        for atom in local_atoms:
            if _distance(point, atom.coord) < (_vdw(atom) + float(probe)):
                blocked = True
                break
        if not blocked:
            accessible += 1
    return accessible / float(len(_SPHERE_POINTS))


def _min_center_distance(center, atoms):
    if not atoms:
        return float("inf")
    return min(_distance(center, atom.coord) for atom in atoms)


def _point_segment_distance(point, start, end):
    segment = _v_sub(end, start)
    length2 = _dot(segment, segment)
    if length2 < 1e-12:
        return _distance(point, start), 0.0
    t = _dot(_v_sub(point, start), segment) / length2
    t_clamped = max(0.0, min(1.0, t))
    closest = _v_add(start, _v_scale(segment, t_clamped))
    return _distance(point, closest), t_clamped


def _path_obstruction_count(start, end, environment_atoms, radius=1.6):
    count = 0
    for atom in environment_atoms:
        distance, fraction = _point_segment_distance(atom.coord, start, end)
        if 0.12 < fraction < 0.88 and distance < float(radius):
            count += 1
    return count


def _orientation_cosine(cb_coord, sg_coord, target_coord):
    sg_vector = _safe_unit(_v_sub(sg_coord, cb_coord))
    target_vector = _safe_unit(_v_sub(target_coord, cb_coord))
    return max(-1.0, min(1.0, _dot(sg_vector, target_vector)))


def _find_disulfide_residues(*group_sets, aliases=None):
    cys = []
    for groups in group_sets:
        for key, atoms in groups.items():
            if _canonical_resn(key[2], aliases) != "CYS":
                continue
            sg = _find_atom(atoms, "SG")
            if sg is not None:
                cys.append((key, sg))
    bound = set()
    for index, (key1, sg1) in enumerate(cys):
        for key2, sg2 in cys[index + 1:]:
            value = _distance(sg1.coord, sg2.coord)
            if 1.7 <= value <= 2.4:
                bound.add(key1)
                bound.add(key2)
    return bound


# ---------------------------------------------------------------------------
# Binder feature precomputation and pair scoring
# ---------------------------------------------------------------------------


def _anchor_coord(res_atoms, resn):
    if str(resn).upper() == "CYS":
        atom = _find_atom(res_atoms, "SG")
        if atom is not None:
            return tuple(atom.coord), "SG"
    atom = _find_atom(res_atoms, "CB")
    if atom is not None:
        return tuple(atom.coord), "CB"
    n_atom = _find_atom(res_atoms, "N")
    ca_atom = _find_atom(res_atoms, "CA")
    c_atom = _find_atom(res_atoms, "C")
    if n_atom is not None and ca_atom is not None and c_atom is not None:
        return _pseudo_cb_from_backbone(n_atom.coord, ca_atom.coord, c_atom.coord), "pseudo_CB"
    if ca_atom is not None:
        return tuple(ca_atom.coord), "CA"
    return None, ""


def _prefilter_anchor_coords(res_atoms, resn):
    """Coordinates used to decide whether a binder residue is in warhead reach.

    Existing Cys uses both CB and native SG: a side chain pointing away from the
    target can still sit on the interface, and χ1 sampling is only worthwhile
    for those residues.
    """
    coords = []
    if str(resn).upper() == "CYS":
        for name in ("SG", "CB"):
            atom = _find_atom(res_atoms, name)
            if atom is not None:
                coords.append(tuple(atom.coord))
        if coords:
            return coords
    coord, _name = _anchor_coord(res_atoms, resn)
    if coord is not None:
        coords.append(coord)
    return coords


def _residue_to_atoms_min_distance(res_atoms, other_atoms):
    heavy = _heavy(res_atoms)
    if not heavy or not other_atoms:
        return float("inf")
    value, _ = _min_distance(heavy, other_atoms)
    return value


def _geometry_score(distance, ideal_min, ideal_max, candidate_cutoff):
    distance = float(distance)
    if float(ideal_min) <= distance <= float(ideal_max):
        return 3.0
    if distance < float(ideal_min):
        lower = min(2.5, float(ideal_min) - 0.1)
        if distance <= lower:
            return -2.0
        fraction = (distance - lower) / max(float(ideal_min) - lower, 0.1)
        return -2.0 + 5.0 * fraction
    if float(candidate_cutoff) <= float(ideal_max):
        return 0.0
    fraction = (float(candidate_cutoff) - distance) / (
        float(candidate_cutoff) - float(ideal_max)
    )
    return max(0.0, min(3.0, 3.0 * fraction))


def _orientation_score(cosine):
    return 0.75 * max(-1.0, min(1.0, float(cosine)))


def _exposure_score(fraction):
    if fraction is None:
        return 0.0
    fraction = float(fraction)
    if fraction >= 0.45:
        return 1.5
    if fraction >= 0.25:
        return 1.0
    if fraction >= 0.12:
        return 0.3
    if fraction >= 0.05:
        return -0.5
    return -1.5


def _class_from_score(score):
    if float(score) >= 5.0:
        return "A"
    if float(score) >= 3.0:
        return "B"
    if float(score) >= 1.0:
        return "C"
    return "D"


def _worst_class(classes):
    values = [value for value in classes if value in CLASS_ORDER]
    return max(values, key=lambda value: CLASS_ORDER[value]) if values else ""


def _apply_class_caps(uncapped_class, caps):
    cap_class = _worst_class([cap for cap, reason in caps])
    if not cap_class or CLASS_ORDER[cap_class] <= CLASS_ORDER[uncapped_class]:
        return uncapped_class, "", ""
    reasons = [
        reason for cap, reason in caps
        if CLASS_ORDER.get(cap, -1) >= CLASS_ORDER[cap_class]
    ]
    return cap_class, cap_class, ";".join(dict.fromkeys(reasons))


def _compute_binder_feature(
    bkey,
    batoms,
    target_res,
    binder_res,
    all_target_heavy,
    terminal_status,
    disulfide_bound,
    config,
    sasa_free=None,
    sasa_bound=None,
    aliases=None,
    extra_atoms=None,
):
    bch, bresi, bresn = bkey
    bresn_std = _canonical_resn(bresn, aliases)
    extra_atoms = list(extra_atoms or [])
    sidechain = _sidechain(batoms)
    sidechain_contact_pairs = _contact_pairs(sidechain, all_target_heavy, config.contact_cutoff) if sidechain else 0
    polar_contacts = _polar_contact_pairs(sidechain, all_target_heavy, config.polar_cutoff, aliases=aliases) if sidechain else 0
    polar_min_distance = _polar_min_distance(sidechain, all_target_heavy, aliases=aliases) if sidechain else None
    salt_bridges = _salt_bridge_pairs(sidechain, all_target_heavy, config.salt_bridge_cutoff, aliases=aliases) if sidechain else 0
    contacted_count, contacted_labels = _unique_contact_residues(sidechain, target_res, config.contact_cutoff)
    other_binder = {key: value for key, value in binder_res.items() if key != bkey}
    internal_count, _ = _unique_contact_residues(
        sidechain or batoms, other_binder, config.internal_contact_cutoff
    )

    free_asa = (sasa_free or {}).get(bkey)
    bound_asa = (sasa_bound or {}).get(bkey)
    relative_sasa = None
    buried_area = None
    buried_fraction = None
    if free_asa is not None and MAX_ASA.get(bresn_std):
        relative_sasa = max(0.0, float(free_asa) / MAX_ASA[bresn_std])
    if free_asa is not None and bound_asa is not None:
        buried_area = max(0.0, float(free_asa) - float(bound_asa))
        if float(free_asa) > 1e-6:
            buried_fraction = buried_area / float(free_asa)

    contact_penalty = min(sidechain_contact_pairs * 0.08, 1.5) + min(contacted_count * 0.25, 1.0)
    polar_penalty = min(polar_contacts * 1.2, 3.0)
    salt_penalty = min(salt_bridges * 1.8, 3.6)
    residue_penalty = 0.0
    if bresn_std in AROMATIC:
        residue_penalty += 0.7
    if bresn_std in CHARGED:
        residue_penalty += 0.4
    if bresn_std == "PRO":
        residue_penalty += 2.5
    elif bresn_std == "GLY":
        residue_penalty += 1.0
    burial_penalty = 0.0
    if buried_fraction is not None:
        if buried_fraction >= 0.60:
            burial_penalty = 0.8
        elif buried_fraction >= 0.30:
            burial_penalty = 0.4
    elif relative_sasa is None and internal_count >= 8:
        burial_penalty = 0.5
    hotspot_penalty = contact_penalty + polar_penalty + salt_penalty + residue_penalty + burial_penalty

    risk_flags = []
    if bresn_std == "PRO": risk_flags.append("pro_backbone")
    if bresn_std == "GLY": risk_flags.append("gly_backbone")
    if bresn_std in AROMATIC: risk_flags.append("aromatic_hotspot_risk")
    if bresn_std in CHARGED: risk_flags.append("charged_residue")
    if polar_contacts: risk_flags.append("polar_proximity")
    if polar_min_distance is not None and polar_min_distance <= 3.2:
        risk_flags.append("critical_polar_contact")
    if salt_bridges: risk_flags.append("salt_bridge")
    if sidechain_contact_pairs >= 6: risk_flags.append("dense_interface_contacts")
    if buried_fraction is not None and buried_fraction >= 0.60:
        risk_flags.append("strong_interface_burial")
    risk_flags.extend(terminal_status.get(bkey, ()))
    if bresn_std == "CYS":
        if bkey in disulfide_bound:
            risk_flags.append("disulfide_bound_binder_cys")
        else:
            risk_flags.append("existing_free_cys_check")

    rotamers = _build_cys_rotamers(batoms, bresn_std) if config.compute_rotamers else []
    binder_other_atoms = _truncated_environment(binder_res, bkey, include_mutated_backbone=False) + extra_atoms
    binder_exposure_atoms = _truncated_environment(binder_res, bkey, include_mutated_backbone=True)
    target_atoms = list(all_target_heavy) + extra_atoms
    bound_exposure_atoms = binder_exposure_atoms + target_atoms
    evaluated = []
    for rotamer in rotamers:
        sg = rotamer["sg_coord"]
        binder_clearance = _min_center_distance(sg, binder_other_atoms)
        target_clearance = _min_center_distance(sg, target_atoms)
        min_clearance = min(binder_clearance, target_clearance)
        free_fraction = _approx_accessible_fraction(sg, binder_exposure_atoms)
        bound_fraction = _approx_accessible_fraction(sg, bound_exposure_atoms)
        item = dict(rotamer)
        item.update({
            "binder_clearance_A": binder_clearance,
            "target_clearance_A": target_clearance,
            "min_clearance_A": min_clearance,
            "sg_accessible_fraction_free": free_fraction,
            "sg_accessible_fraction_bound": bound_fraction,
            "valid": min_clearance >= config.clash_cutoff,
        })
        evaluated.append(item)

    anchor_coord, anchor_name = _anchor_coord(batoms, bresn_std)
    return {
        "binder_key": bkey,
        "binder_chain": bch,
        "binder_resi": bresi,
        "binder_resn": bresn,
        "binder_resn_std": bresn_std,
        "binder_site": "%s/%s" % (bch, _residue_title(bresn, bresi)),
        "mutation": "%s:%s%sC" % (bch, AA_ONE_LETTER.get(bresn_std, "X"), bresi),
        "anchor_coord": anchor_coord,
        "anchor_atom": anchor_name,
        "mutability_score": MUTABILITY.get(bresn_std, -0.5),
        "sidechain_contact_pairs": sidechain_contact_pairs,
        "contacted_target_residues": contacted_count,
        "contacted_target_labels": ",".join(sorted(contacted_labels)),
        "sidechain_polar_contacts": polar_contacts,
        "sidechain_polar_min_A": polar_min_distance,
        "salt_bridge_contacts": salt_bridges,
        "internal_contact_residues": internal_count,
        "free_residue_sasa_A2": free_asa,
        "bound_residue_sasa_A2": bound_asa,
        "free_relative_sasa": relative_sasa,
        "interface_buried_sasa_A2": buried_area,
        "interface_buried_fraction": buried_fraction,
        "contact_penalty": contact_penalty,
        "polar_penalty": polar_penalty,
        "salt_bridge_penalty": salt_penalty,
        "residue_penalty": residue_penalty,
        "burial_penalty": burial_penalty,
        "hotspot_penalty": hotspot_penalty,
        "risk_flags": risk_flags,
        "terminal_status": ";".join(terminal_status.get(bkey, ())),
        "rotamers": evaluated,
    }


def _precompute_binder_features(
    binder_res, target_res, all_target_heavy, config, disulfide_bound,
    sasa_free=None, sasa_bound=None, candidate_keys=None, aliases=None,
    extra_atoms=None,
):
    terminal_status = _terminal_statuses(binder_res)
    features = {}
    keys = list(candidate_keys) if candidate_keys is not None else list(binder_res)
    for bkey in sorted(keys, key=_residue_sort_key):
        features[bkey] = _compute_binder_feature(
            bkey, binder_res[bkey], target_res, binder_res, all_target_heavy,
            terminal_status, disulfide_bound, config,
            sasa_free=sasa_free, sasa_bound=sasa_bound, aliases=aliases,
            extra_atoms=extra_atoms,
        )
    return features


def _pair_caps(feature, selected_rotamer, rotamers_present, compute_rotamers):
    caps = []
    bresn = feature.get("binder_resn_std") or feature["binder_resn"]
    flags = set(feature["risk_flags"])
    if bresn == "PRO":
        caps.append(("D", "pro_to_cys_backbone_risk"))
    if "disulfide_bound_binder_cys" in flags:
        caps.append(("D", "binder_cys_in_disulfide"))
    if bresn == "GLY":
        caps.append(("C", "gly_to_cys_backbone_risk"))
    if compute_rotamers:
        if not rotamers_present:
            caps.append(("C", "cys_rotamer_geometry_unavailable"))
        elif selected_rotamer is None:
            caps.append(("D", "all_cys_rotamers_clash"))
    else:
        caps.append(("C", "cys_rotamers_not_computed"))
    if feature.get("interface_buried_fraction") is not None and feature["interface_buried_fraction"] >= 0.75:
        caps.append(("C", "strong_interface_burial"))
    if feature["salt_bridge_contacts"] > 0 or feature["sidechain_polar_contacts"] >= 2:
        caps.append(("C", "likely_interface_hotspot"))
    if feature.get("sidechain_polar_min_A") is not None and feature["sidechain_polar_min_A"] <= 3.2:
        caps.append(("C", "critical_polar_contact"))
    if selected_rotamer is not None:
        if selected_rotamer["sg_accessible_fraction_free"] < 0.03:
            caps.append(("D", "modeled_sg_severely_buried"))
        elif selected_rotamer["sg_accessible_fraction_free"] < 0.08:
            caps.append(("C", "modeled_sg_low_exposure"))
    if flags.intersection({"n_terminus", "c_terminus", "fragment_n_terminus", "fragment_c_terminus"}):
        caps.append(("C", "terminal_or_fragment_end"))
    return caps


def _evaluate_pair_for_target_atom(feature, target_key, target_atom, environment_atoms, config, aliases=None):
    fallback_distance = None
    if feature["anchor_coord"] is not None:
        fallback_distance = _distance(feature["anchor_coord"], target_atom.coord)

    rotamers = feature["rotamers"]
    valid_rotamers = [rot for rot in rotamers if rot["valid"]]
    pair_options = []
    excluded_keys = {feature["binder_key"], target_key}
    path_atoms = [atom for atom in environment_atoms if _reskey(atom) not in excluded_keys]
    for rotamer in valid_rotamers:
        sg = rotamer["sg_coord"]
        distance = _distance(sg, target_atom.coord)
        cosine = _orientation_cosine(rotamer["cb_coord"], sg, target_atom.coord)
        obstruction = _path_obstruction_count(sg, target_atom.coord, path_atoms, config.path_radius)
        geometry = _geometry_score(distance, config.ideal_min, config.ideal_max, config.candidate_cutoff)
        orientation = _orientation_score(cosine)
        exposure = _exposure_score(rotamer["sg_accessible_fraction_free"])
        path_penalty = min(obstruction * 0.35, 1.75)
        raw = (
            feature["mutability_score"] + geometry + orientation + exposure
            - feature["hotspot_penalty"] - path_penalty
        )
        option = dict(rotamer)
        option.update({
            "target_atom": target_atom,
            "sg_target_distance_A": distance,
            "orientation_cosine": cosine,
            "geometry_score": geometry,
            "orientation_score": orientation,
            "sg_exposure_score": exposure,
            "path_obstruction_atoms": obstruction,
            "path_penalty": path_penalty,
            "raw_score": raw,
        })
        pair_options.append(option)

    selected = None
    eligible_options = [
        item for item in pair_options
        if item["sg_target_distance_A"] <= config.candidate_cutoff
    ]
    if eligible_options:
        eligible_options.sort(key=lambda item: (
            -item["raw_score"], item["sg_target_distance_A"],
            -item["orientation_cosine"], -item["sg_accessible_fraction_free"],
            str(item["chi1_label"]),
        ))
        selected = eligible_options[0]

    if selected is not None:
        candidate_distance = selected["sg_target_distance_A"]
    else:
        candidate_distance = fallback_distance
    if candidate_distance is None or candidate_distance > config.candidate_cutoff:
        # Keep all-clash sites visible when their unmutated anchor is within
        # reach; otherwise they cannot inform the requested target pair.
        if not (rotamers and not valid_rotamers and fallback_distance is not None and fallback_distance <= config.candidate_cutoff):
            return None

    if selected is None:
        geometry = _geometry_score(fallback_distance, config.ideal_min, config.ideal_max, config.candidate_cutoff)
        raw = feature["mutability_score"] + geometry - feature["hotspot_penalty"]
        option_values = {
            "chi1_deg": None, "chi1_label": "", "sg_coord": None,
            "sg_target_distance_A": fallback_distance, "orientation_cosine": 0.0,
            "geometry_score": geometry, "orientation_score": 0.0,
            "sg_exposure_score": 0.0, "path_obstruction_atoms": None,
            "path_penalty": 0.0, "raw_score": raw,
            "sg_accessible_fraction_free": None,
            "sg_accessible_fraction_bound": None,
            "min_clearance_A": None,
        }
    else:
        option_values = selected

    uncapped = _class_from_score(option_values["raw_score"])
    caps = _pair_caps(feature, selected, bool(rotamers), config.compute_rotamers)
    final_class, class_cap, downgrade_reasons = _apply_class_caps(uncapped, caps)
    risk_flags = list(feature["risk_flags"])
    scored_distance = option_values.get("sg_target_distance_A")
    if scored_distance is None:
        scored_distance = fallback_distance
    if scored_distance is not None and scored_distance < config.ideal_min:
        risk_flags.append("tighter_than_preferred")
    elif scored_distance is not None and scored_distance > config.ideal_max:
        risk_flags.append("longer_than_preferred")
    if option_values.get("path_obstruction_atoms", 0) and option_values["path_obstruction_atoms"] >= 3:
        risk_flags.append("straight_path_obstructed")

    sg_coord = option_values.get("sg_coord")
    binder_resn_std = feature.get("binder_resn_std") or _canonical_resn(feature["binder_resn"], aliases)
    target_resn_std = _canonical_resn(target_key[2], aliases)
    return {
        "binder_chain": feature["binder_chain"],
        "binder_resi": feature["binder_resi"],
        "binder_resn": feature["binder_resn"],
        "binder_resn_std": binder_resn_std,
        "binder_site": feature.get("binder_site") or "%s/%s" % (
            feature["binder_chain"], _residue_title(feature["binder_resn"], feature["binder_resi"])
        ),
        "mutation": feature.get("mutation") or "%s:%s%sC" % (
            feature["binder_chain"], AA_ONE_LETTER.get(binder_resn_std, "X"), feature["binder_resi"]
        ),
        "anchor_atom": feature["anchor_atom"],
        "target_chain": target_key[0],
        "target_resi": target_key[1],
        "target_resn": target_key[2],
        "target_resn_std": target_resn_std,
        "target_site": "%s/%s" % (target_key[0], _residue_title(target_key[2], target_key[1])),
        "target_atom": str(getattr(target_atom, "name", "")),
        "fallback_anchor_distance_A": fallback_distance,
        "sg_target_distance_A": option_values["sg_target_distance_A"],
        "best_sg_x": sg_coord[0] if sg_coord is not None else None,
        "best_sg_y": sg_coord[1] if sg_coord is not None else None,
        "best_sg_z": sg_coord[2] if sg_coord is not None else None,
        "best_chi1_deg": option_values.get("chi1_deg"),
        "best_chi1_label": option_values.get("chi1_label", ""),
        "valid_cys_rotamers": len(valid_rotamers),
        "total_cys_rotamers": len(rotamers),
        "orientation_cosine": option_values.get("orientation_cosine"),
        "path_obstruction_atoms": option_values.get("path_obstruction_atoms"),
        "sg_accessible_fraction_free": option_values.get("sg_accessible_fraction_free"),
        "sg_accessible_fraction_bound": option_values.get("sg_accessible_fraction_bound"),
        "minimum_rotamer_clearance_A": option_values.get("min_clearance_A"),
        "mutability_score": feature["mutability_score"],
        "geometry_score": option_values["geometry_score"],
        "orientation_score": option_values["orientation_score"],
        "sg_exposure_score": option_values["sg_exposure_score"],
        "path_penalty": option_values["path_penalty"],
        "contact_penalty": feature["contact_penalty"],
        "polar_penalty": feature["polar_penalty"],
        "salt_bridge_penalty": feature["salt_bridge_penalty"],
        "residue_penalty": feature["residue_penalty"],
        "burial_penalty": feature["burial_penalty"],
        "hotspot_penalty": feature["hotspot_penalty"],
        "raw_score": option_values["raw_score"],
        "uncapped_class": uncapped,
        "final_class": final_class,
        "class_cap": class_cap,
        "downgrade_reasons": downgrade_reasons,
        "sidechain_contact_pairs": feature["sidechain_contact_pairs"],
        "contacted_target_residues": feature["contacted_target_residues"],
        "contacted_target_labels": feature["contacted_target_labels"],
        "sidechain_polar_contacts": feature["sidechain_polar_contacts"],
        "sidechain_polar_min_A": feature["sidechain_polar_min_A"],
        "salt_bridge_contacts": feature["salt_bridge_contacts"],
        "internal_contact_residues": feature["internal_contact_residues"],
        "free_residue_sasa_A2": feature["free_residue_sasa_A2"],
        "bound_residue_sasa_A2": feature["bound_residue_sasa_A2"],
        "free_relative_sasa": feature["free_relative_sasa"],
        "interface_buried_sasa_A2": feature["interface_buried_sasa_A2"],
        "interface_buried_fraction": feature["interface_buried_fraction"],
        "terminal_status": feature["terminal_status"],
        "risk_flags": ";".join(dict.fromkeys(risk_flags)),
    }


def _pair_sort_key(record):
    return (
        CLASS_ORDER[record["final_class"]],
        -float(record["raw_score"]),
        float(record["sg_target_distance_A"] if record["sg_target_distance_A"] is not None else 999.0),
        str(record["binder_chain"]), _natural_residue_key(record["binder_resi"]),
        str(record["target_chain"]), _natural_residue_key(record["target_resi"]),
        str(record["target_atom"]),
    )


def _pair_target_label(record):
    """Readable unique pair target tag such as ``B/Tyr151 OH``."""
    site = record.get("target_site") or "%s/%s" % (
        record["target_chain"], _residue_title(record["target_resn"], record["target_resi"])
    )
    atom = str(record.get("target_atom") or "").strip()
    return ("%s %s" % (site, atom)).strip()


def _pair_target_title(record):
    return _pair_target_label(record)


def _target_option_sort_key(label):
    text = str(label)
    chain, rest = (text.split("/", 1) + [""])[:2]
    if " " in rest:
        site, atom = rest.rsplit(" ", 1)
    else:
        site, atom = rest, ""
    resi = re.sub(r"^[A-Za-z]+", "", site)
    return (chain, _natural_residue_key(resi or "0"), atom, text)


def _analyze_groups(
    target_res,
    binder_res,
    interface_cutoff=5.0,
    candidate_cutoff=12.0,
    ideal_min=5.0,
    ideal_max=10.0,
    profile="broad_nucleophile",
    target_types=None,
    reactive_atoms="",
    contact_cutoff=4.0,
    polar_cutoff=3.5,
    salt_bridge_cutoff=4.0,
    internal_contact_cutoff=4.5,
    clash_cutoff=2.0,
    path_radius=1.6,
    include_disulfide_targets=False,
    compute_rotamers=True,
    sasa_free=None,
    sasa_bound=None,
    resn_aliases=None,
    extra_atoms=None,
):
    aliases = _parse_alias_option(resn_aliases)
    resolved_profile = _resolve_profile(
        profile, reactive_atoms=reactive_atoms, target_types=target_types, aliases=aliases
    )
    config = ScanConfig(
        interface_cutoff=float(interface_cutoff),
        candidate_cutoff=float(candidate_cutoff),
        ideal_min=float(ideal_min),
        ideal_max=float(ideal_max),
        contact_cutoff=float(contact_cutoff),
        polar_cutoff=float(polar_cutoff),
        salt_bridge_cutoff=float(salt_bridge_cutoff),
        internal_contact_cutoff=float(internal_contact_cutoff),
        clash_cutoff=float(clash_cutoff),
        path_radius=float(path_radius),
        include_disulfide_targets=_as_bool(include_disulfide_targets),
        compute_rotamers=_as_bool(compute_rotamers),
    ).validate()

    all_target_heavy = _heavy([atom for atoms in target_res.values() for atom in atoms])
    all_binder_heavy = _heavy([atom for atoms in binder_res.values() for atom in atoms])
    if not all_target_heavy:
        raise ValueError("Target selection contains no heavy atoms")
    if not all_binder_heavy:
        raise ValueError("Binder selection contains no heavy atoms")

    overlap = {_atom_identity(atom) for atom in all_target_heavy}.intersection(
        {_atom_identity(atom) for atom in all_binder_heavy}
    )
    if overlap:
        raise ValueError("Target and binder selections overlap (%d heavy atom identities)" % len(overlap))

    disulfide_bound = _find_disulfide_residues(target_res, binder_res, aliases=aliases)
    target_records = []
    target_context = []
    for tkey in sorted(target_res, key=_residue_sort_key):
        tatoms = target_res[tkey]
        tch, tresi, tresn = tkey
        tresn_std = _canonical_resn(tresn, aliases)
        wanted = resolved_profile.reactive_atoms.get(tresn_std)
        if not wanted:
            continue
        reactive = [atom for atom in tatoms if str(getattr(atom, "name", "")) in wanted]
        if not reactive:
            continue
        reactive_min, reactive_pair = _min_distance(reactive, all_binder_heavy)
        residue_min, _ = _min_distance(_heavy(tatoms), all_binder_heavy)
        if reactive_min > config.interface_cutoff:
            continue
        is_disulfide = tresn_std == "CYS" and tkey in disulfide_bound
        available = (not is_disulfide) or config.include_disulfide_targets
        record = {
            "target_chain": tch,
            "target_resi": tresi,
            "target_resn": tresn,
            "target_resn_std": tresn_std,
            "target_site": "%s/%s" % (tch, _residue_title(tresn, tresi)),
            "reactive_atoms": ",".join(str(getattr(atom, "name", "")) for atom in reactive),
            "nearest_reactive_atom": str(getattr(reactive_pair[0], "name", "")) if reactive_pair[0] else "",
            "nearest_binder_atom": str(getattr(reactive_pair[1], "name", "")) if reactive_pair[1] else "",
            "reactive_to_binder_min_A": reactive_min,
            "residue_to_binder_min_A": residue_min,
            "disulfide_bound": is_disulfide,
            "available_for_pairing": available,
        }
        target_records.append(record)
        target_context.append((tkey, tatoms, reactive, record))

    # Coarse prefilter: SG can move about one CB-SG bond length from the
    # observed/pseudo CB. Residues farther than candidate_cutoff plus this
    # reach buffer cannot become valid pairs and need no expensive rotamer or
    # exposure calculation. Existing Cys have an extra PPI-interface gate:
    # only Cys already contacting the target protein are χ1-sampled.
    all_reactive_atoms = [
        atom for tkey, tatoms, reactive, record in target_context
        if record["available_for_pairing"] for atom in reactive
    ]
    candidate_keys = []
    skipped_cys_not_interface = 0
    reach_buffer = CYS_CB_SG_LENGTH + 0.5
    if all_reactive_atoms:
        for bkey, batoms in binder_res.items():
            resn_std = _canonical_resn(bkey[2], aliases)
            anchors = _prefilter_anchor_coords(batoms, resn_std)
            if not anchors:
                continue
            nearest = min(
                _distance(coord, atom.coord)
                for coord in anchors
                for atom in all_reactive_atoms
            )
            if nearest > config.candidate_cutoff + reach_buffer:
                continue
            if resn_std == "CYS":
                if _residue_to_atoms_min_distance(batoms, all_target_heavy) > config.interface_cutoff:
                    skipped_cys_not_interface += 1
                    continue
            candidate_keys.append(bkey)

    extra_heavy = _heavy(extra_atoms or [])
    binder_features = _precompute_binder_features(
        binder_res, target_res, all_target_heavy, config, disulfide_bound,
        sasa_free=sasa_free, sasa_bound=sasa_bound, candidate_keys=candidate_keys,
        aliases=aliases, extra_atoms=extra_heavy,
    )
    environment_atoms = all_target_heavy + all_binder_heavy + extra_heavy
    pair_records = []
    for tkey, tatoms, reactive_atoms_list, target_record in target_context:
        if not target_record["available_for_pairing"]:
            continue
        for bkey in sorted(binder_features, key=_residue_sort_key):
            feature = binder_features[bkey]
            best_for_residue = None
            for target_atom in reactive_atoms_list:
                candidate = _evaluate_pair_for_target_atom(
                    feature, tkey, target_atom, environment_atoms, config, aliases=aliases
                )
                if candidate is None:
                    continue
                if best_for_residue is None or _pair_sort_key(candidate) < _pair_sort_key(best_for_residue):
                    best_for_residue = candidate
            if best_for_residue is not None:
                pair_records.append(best_for_residue)

    pair_records.sort(key=_pair_sort_key)
    for index, record in enumerate(pair_records, 1):
        record["pair_rank"] = index

    aggregate = {}
    for record in pair_records:
        bkey = (record["binder_chain"], record["binder_resi"], record["binder_resn"])
        label = _pair_target_label(record)
        if bkey not in aggregate:
            aggregate[bkey] = {
                "binder_chain": record["binder_chain"],
                "binder_resi": record["binder_resi"],
                "binder_resn": record["binder_resn"],
                "binder_resn_std": record.get("binder_resn_std", ""),
                "binder_site": record.get("binder_site", ""),
                "mutation": record.get("mutation", ""),
                "raw_score": record["raw_score"],
                "uncapped_class": record["uncapped_class"],
                "final_class": record["final_class"],
                "class_cap": record["class_cap"],
                "downgrade_reasons": record["downgrade_reasons"],
                "best_target": label,
                "best_sg_target_distance_A": record["sg_target_distance_A"],
                "best_chi1_deg": record["best_chi1_deg"],
                "best_sg_x": record["best_sg_x"],
                "best_sg_y": record["best_sg_y"],
                "best_sg_z": record["best_sg_z"],
                "risk_flags": record["risk_flags"],
                "target_options": [label],
            }
        else:
            aggregate[bkey]["target_options"].append(label)
            # pair_records are already sorted by final class then score; the
            # first pair is the aggregate best and is intentionally retained.

    binder_records = list(aggregate.values())
    for record in binder_records:
        record["target_options"] = ",".join(sorted(
            set(record["target_options"]),
            key=_target_option_sort_key,
        ))
    binder_records.sort(key=lambda record: (
        CLASS_ORDER[record["final_class"]], -float(record["raw_score"]),
        float(record["best_sg_target_distance_A"] if record["best_sg_target_distance_A"] is not None else 999.0),
        str(record["binder_chain"]), _natural_residue_key(record["binder_resi"]),
    ))
    for index, record in enumerate(binder_records, 1):
        record["binder_rank"] = index

    target_records.sort(key=lambda record: (
        float(record["reactive_to_binder_min_A"]), str(record["target_chain"]),
        _natural_residue_key(record["target_resi"]),
    ))
    parameters = config.to_dict()
    parameters.update({
        "profile": resolved_profile.name,
        "reactive_atoms": {
            key: list(value) for key, value in resolved_profile.reactive_atoms.items()
        },
        "resn_aliases": dict(aliases),
        "n_hetero_atoms": len(extra_heavy),
        "skipped_cys_not_interface": int(skipped_cys_not_interface),
    })
    return {
        "target_records": target_records,
        "pair_records": pair_records,
        "binder_records": binder_records,
        "parameters": parameters,
    }


# ---------------------------------------------------------------------------
# PyMOL selection, SASA, output, and visualization
# ---------------------------------------------------------------------------


def _quote_selector(value):
    value = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return '"%s"' % value


def _model_clause(object_name):
    return "model %s" % _quote_selector(object_name)


def _hetero_selection(object_name):
    return (
        "(%s) and (organic or metals or inorganic) and not solvent "
        "and not polymer.protein and not hydro"
        % _model_clause(object_name)
    )


def _load_hetero_atoms(object_name, state=1):
    """Non-solvent ligands/ions in the same object, excluded from protein selections."""
    try:
        atoms = _atoms(_hetero_selection(object_name), int(state))
    except Exception:
        return []
    return _heavy(atoms)


def _parse_chain_ids(value):
    if _blank(value):
        return []
    return [item for item in re.split(r"[,+;/\s]+", str(value).strip()) if item]


def _chain_expression(value):
    chains = _parse_chain_ids(value)
    if not chains:
        raise ValueError("At least one chain ID is required when no explicit selection is supplied")
    return "(" + " or ".join("chain %s" % _quote_selector(chain) for chain in chains) + ")"


def _available_protein_chains(object_name, state=1):
    selection = "(%s) and polymer.protein" % _model_clause(object_name)
    try:
        chains = cmd.get_chains(selection, int(state))
    except TypeError:
        chains = cmd.get_chains(selection)
    except Exception:
        chains = []
    return [str(chain) for chain in (chains or [])]


def _validate_chain_ids(object_name, chain_spec, role, state=1):
    requested = _parse_chain_ids(chain_spec)
    if not requested:
        return requested
    available = _available_protein_chains(object_name, state)
    missing = [chain for chain in requested if chain not in available]
    if missing:
        raise ValueError(
            "%s chain(s) not found: %s. Available protein chains: %s "
            "(IDs are case-sensitive; combine several as H+L)"
            % (role, ", ".join(missing), ", ".join(available) or "none")
        )
    return requested


def _build_scan_selections(object_name, target_chain="", binder_chain="", target_sel="", binder_sel="", state=1):
    base = "(%s) and polymer.protein" % _model_clause(object_name)
    if _blank(target_sel):
        target_chains = _validate_chain_ids(object_name, target_chain, "TARGET", state)
    else:
        target_chains = []
    if _blank(binder_sel):
        binder_chains = _validate_chain_ids(object_name, binder_chain, "BINDER", state)
    else:
        binder_chains = []
    if target_chains and binder_chains:
        overlap = sorted(set(target_chains).intersection(binder_chains))
        if overlap:
            raise ValueError(
                "Target and binder chains must not overlap; both contain: %s"
                % ", ".join(overlap)
            )
    if not _blank(target_sel):
        target = "(%s) and (%s)" % (base, str(target_sel).strip())
    else:
        target = "(%s) and %s" % (base, _chain_expression(target_chain))
    if not _blank(binder_sel):
        binder = "(%s) and (%s)" % (base, str(binder_sel).strip())
    else:
        binder = "(%s) and %s" % (base, _chain_expression(binder_chain))
    return target, binder


def _residue_selection(object_name, chain, resi, resn=None):
    parts = [
        _model_clause(object_name), "chain %s" % _quote_selector(chain),
        "resi %s" % _quote_selector(resi),
    ]
    if resn:
        parts.append("resn %s" % _quote_selector(resn))
    return "(" + " and ".join(parts) + ")"


def _object_list(selection):
    """Normalize ``cmd.get_object_list`` to always return a list.

    Real PyMOL may return None when a selection matches nothing; treating that
    as an empty list keeps "object not found" errors readable.
    """
    try:
        found = cmd.get_object_list(selection)
    except Exception:
        return []
    return list(found) if found else []


def _resolve_object(object_name="auto"):
    value = str(object_name or "auto").strip()
    auto_mode = value.lower() in {"", "auto", "none"}
    if auto_mode:
        objects = _object_list("(polymer.protein)")
    else:
        # Prefer an exact loaded-object hit before querying a selection that
        # can emit a noisy Selector-Error for unknown model names.
        names = list(cmd.get_names("objects") or [])
        if value in names:
            protein_hit = _object_list("(%s) and polymer.protein" % _model_clause(value))
            objects = protein_hit or [value]
        else:
            # Unknown object name: do not query a dead model selection; that
            # only produces a noisy Selector-Error in the PyMOL log.
            objects = []
    if len(objects) == 1:
        return objects[0]
    available = ", ".join(_object_list("(polymer.protein)")) or "none"
    if not objects:
        if auto_mode:
            raise ValueError(
                "No protein object is loaded. Load a structure first, e.g. "
                "load complex.cif, complex"
            )
        raise ValueError(
            "Molecular object '%s' not found or contains no protein. "
            "Available protein objects: %s. "
            "Note PyMOL rewrites object names on load: reserved keywords gain "
            "an underscore ('alt' -> 'alt_') and spaces/invalid characters are "
            "replaced ('my complex' -> 'my_complex')."
            % (value, available)
        )
    raise ValueError(
        "Select exactly one molecular object; matched %d: %s. "
        "Available protein objects: %s"
        % (len(objects), ", ".join(objects), available)
    )


def _sasa_area_maps(atoms):
    residue_out = defaultdict(float)
    for atom in _heavy(atoms):
        residue_out[_reskey(atom)] += max(0.0, float(getattr(atom, "b", 0.0) or 0.0))
    return dict(residue_out)


def _compute_sasa_maps(target_sel, binder_sel, state=1):
    tmp_bound = cmd.get_unused_name("_warhead_bound", 1)
    tmp_free = cmd.get_unused_name("_warhead_free", 1)
    old_dot_solvent = cmd.get("dot_solvent")
    old_dot_density = cmd.get("dot_density")
    try:
        cmd.create(tmp_bound, "(%s) or (%s)" % (target_sel, binder_sel), int(state), 1)
        cmd.create(tmp_free, "(%s)" % binder_sel, int(state), 1)
        cmd.remove("(%s) and hydro" % tmp_bound)
        cmd.remove("(%s) and hydro" % tmp_free)
        cmd.set("dot_solvent", 1)
        cmd.set("dot_density", 3)
        cmd.get_area(tmp_bound, state=1, load_b=1)
        cmd.get_area(tmp_free, state=1, load_b=1)
        return _sasa_area_maps(_atoms(tmp_free, 1)), _sasa_area_maps(_atoms(tmp_bound, 1))
    finally:
        cmd.delete(tmp_bound)
        cmd.delete(tmp_free)
        try:
            cmd.set("dot_solvent", old_dot_solvent)
            cmd.set("dot_density", old_dot_density)
        except Exception:
            pass


def _ensure_parent(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)


def _fmt(value, digits=4):
    if value is None:
        return ""
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, float):
        return ("%%.%df" % digits) % value
    return value


def _write_csv(path, records, fields):
    _ensure_parent(path)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow({field: _fmt(record.get(field)) for field in fields})


TARGET_FIELDS = [
    "target_chain", "target_resi", "target_resn", "target_resn_std", "target_site",
    "reactive_atoms",
    "nearest_reactive_atom", "nearest_binder_atom", "reactive_to_binder_min_A",
    "residue_to_binder_min_A", "disulfide_bound", "available_for_pairing",
]
PAIR_FIELDS = [
    "pair_rank", "raw_score", "uncapped_class", "final_class", "class_cap",
    "downgrade_reasons", "binder_site", "mutation",
    "binder_chain", "binder_resi", "binder_resn", "binder_resn_std",
    "target_site", "target_chain", "target_resi", "target_resn", "target_resn_std",
    "target_atom", "anchor_atom",
    "fallback_anchor_distance_A", "sg_target_distance_A", "best_chi1_deg",
    "best_chi1_label", "best_sg_x", "best_sg_y", "best_sg_z",
    "valid_cys_rotamers", "total_cys_rotamers", "orientation_cosine",
    "path_obstruction_atoms", "sg_accessible_fraction_free",
    "sg_accessible_fraction_bound", "minimum_rotamer_clearance_A",
    "mutability_score", "geometry_score", "orientation_score",
    "sg_exposure_score", "path_penalty", "contact_penalty", "polar_penalty",
    "salt_bridge_penalty", "residue_penalty", "burial_penalty", "hotspot_penalty",
    "sidechain_contact_pairs", "contacted_target_residues", "contacted_target_labels",
    "sidechain_polar_contacts", "sidechain_polar_min_A", "salt_bridge_contacts", "internal_contact_residues",
    "free_residue_sasa_A2", "bound_residue_sasa_A2", "free_relative_sasa",
    "interface_buried_sasa_A2", "interface_buried_fraction", "terminal_status",
    "risk_flags",
]
BINDER_FIELDS = [
    "binder_rank", "raw_score", "uncapped_class", "final_class", "class_cap",
    "downgrade_reasons", "binder_site", "mutation",
    "binder_chain", "binder_resi", "binder_resn", "binder_resn_std",
    "best_target", "best_sg_target_distance_A", "best_chi1_deg", "best_sg_x",
    "best_sg_y", "best_sg_z", "target_options", "risk_flags",
]


def _write_outputs(prefix, result):
    prefix = str(prefix or "warhead_scan")
    if prefix.lower().endswith(".csv"):
        prefix = prefix[:-4]
    paths = (
        prefix + "_targets.csv",
        prefix + "_pairs.csv",
        prefix + "_binder_rank.csv",
        prefix + "_metadata.json",
    )
    _write_csv(paths[0], result.get("target_records", []), TARGET_FIELDS)
    _write_csv(paths[1], result.get("pair_records", []), PAIR_FIELDS)
    _write_csv(paths[2], result.get("binder_records", []), BINDER_FIELDS)
    metadata = {
        "script_version": __version__,
        "parameters": result.get("parameters", {}),
        "counts": {
            "targets": len(result.get("target_records", [])),
            "pairs": len(result.get("pair_records", [])),
            "binder_sites": len(result.get("binder_records", [])),
        },
        "top_candidates": result.get("binder_records", [])[:3],
    }
    _ensure_parent(paths[3])
    with open(paths[3], "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2, default=str)
    return paths


def _mutation_site_label(record):
    site = record.get("binder_site") or "%s/%s" % (
        record["binder_chain"], _residue_title(record["binder_resn"], record["binder_resi"])
    )
    return "%s->Cys" % site


def _target_site_label(record):
    if record.get("target_site"):
        return record["target_site"]
    return "%s/%s" % (
        record["target_chain"], _residue_title(record["target_resn"], record["target_resi"])
    )


def _best_pair_for_binder(result, binder_record):
    key = (
        str(binder_record["binder_chain"]), str(binder_record["binder_resi"]),
        str(binder_record["binder_resn"]),
    )
    fallback = None
    for pair in result.get("pair_records", []):
        pkey = (str(pair["binder_chain"]), str(pair["binder_resi"]), str(pair["binder_resn"]))
        if pkey != key:
            continue
        if fallback is None:
            fallback = pair
        if _pair_target_label(pair) == binder_record.get("best_target"):
            return pair
    return fallback


def _visualize(object_name, result, label_top=3, top_n=None):
    del top_n
    label_top = max(0, min(3, int(label_top)))
    names = ["wh_target_nucleophiles", "wh_top_candidates", "warhead_scan_group"]
    names += ["wh_top_%d" % i for i in range(1, 21)]
    names += ["wh_top_label_%d" % i for i in range(1, 21)]
    names += ["wh_top_target_label_%d" % i for i in range(1, 21)]
    names += ["wh_top_distance_%d" % i for i in range(1, 21)]
    for name in names:
        cmd.delete(name)

    present = []
    zoom_parts = []
    target_parts = [
        _residue_selection(object_name, record["target_chain"], record["target_resi"], record["target_resn"])
        for record in result.get("target_records", [])
    ]
    if target_parts:
        cmd.select("wh_target_nucleophiles", " or ".join(target_parts))
        cmd.show("sticks", "wh_target_nucleophiles")
        cmd.color("yellow", "wh_target_nucleophiles")
        present.append("wh_target_nucleophiles")
        zoom_parts.append("wh_target_nucleophiles")

    top_parts = []
    labelled_targets = set()
    target_label_count = 0
    colors = {"A": "green", "B": "cyan", "C": "magenta", "D": "gray70"}
    for index, binder_record in enumerate(result.get("binder_records", [])[:label_top], 1):
        pair = _best_pair_for_binder(result, binder_record)
        if pair is None:
            continue
        binder_sel = _residue_selection(
            object_name, binder_record["binder_chain"], binder_record["binder_resi"], binder_record["binder_resn"]
        )
        target_sel = _residue_selection(
            object_name, pair["target_chain"], pair["target_resi"], pair["target_resn"]
        )
        reactive_sel = "(%s) and name %s" % (target_sel, pair["target_atom"])
        top_name = "wh_top_%d" % index
        label_name = "wh_top_label_%d" % index
        distance_name = "wh_top_distance_%d" % index
        cmd.select(top_name, binder_sel)
        cmd.show("sticks", top_name)
        cmd.color(colors.get(binder_record.get("final_class", "D"), "gray70"), top_name)

        site = binder_record.get("binder_site") or "%s/%s" % (
            binder_record["binder_chain"],
            _residue_title(binder_record["binder_resn"], binder_record["binder_resi"]),
        )
        distance = pair.get("sg_target_distance_A")
        if distance is None:
            label_text = "#%d %s" % (index, site)
        else:
            label_text = "#%d %s  %.2fÅ" % (index, site, float(distance))
        sg_values = (pair.get("best_sg_x"), pair.get("best_sg_y"), pair.get("best_sg_z"))
        if all(value is not None for value in sg_values):
            cmd.pseudoatom(label_name, pos=list(sg_values), label=label_text)
        else:
            anchor_sel = "(%s) and name %s" % (binder_sel, pair.get("anchor_atom") or "CB")
            cmd.pseudoatom(label_name, selection=anchor_sel, label=label_text)
        try:
            cmd.show("spheres", label_name)
            cmd.show("labels", label_name)
            cmd.set("sphere_scale", 0.25, label_name)
            cmd.set("label_size", 16, label_name)
        except Exception:
            pass

        # The visible text identifies the residue (e.g. B:His81), so ND1
        # and NE2 pairings of the same histidine share one target label. Each
        # distance object still points to the pair-specific reactive atom.
        target_key = (
            str(pair["target_chain"]), str(pair["target_resi"]),
            str(pair["target_resn"]),
        )
        target_label_name = None
        if target_key not in labelled_targets:
            labelled_targets.add(target_key)
            target_label_count += 1
            target_label_name = "wh_top_target_label_%d" % target_label_count
            cmd.pseudoatom(target_label_name, selection=reactive_sel, label=_target_site_label(pair))
            try:
                cmd.show("labels", target_label_name)
                cmd.set("label_size", 16, target_label_name)
            except Exception:
                pass

        cmd.distance(distance_name, label_name, reactive_sel)
        try:
            cmd.hide("labels", distance_name)
        except Exception:
            pass
        top_parts.append(binder_sel)
        present.extend([top_name, label_name, distance_name])
        if target_label_name:
            present.append(target_label_name)

    if top_parts:
        cmd.select("wh_top_candidates", " or ".join(top_parts))
        present.append("wh_top_candidates")
        zoom_parts.append("wh_top_candidates")
    if present:
        cmd.group("warhead_scan_group", " ".join(present))
    if zoom_parts:
        cmd.zoom(" or ".join(zoom_parts), 8)


def _records_for_console(records, top_n=0):
    limit = int(top_n)
    return list(records) if limit <= 0 else list(records)[:limit]


def _print_binder_rank_table(records, top_n=0, label_top=3):
    print("Mark Rank Binder         Mutation   Final Raw   Uncap Cap Best target          SG-dist  Downgrade reason")
    labelled = max(0, min(3, int(label_top)))
    for record in _records_for_console(records, top_n):
        site = record.get("binder_site") or "%s/%s" % (
            record["binder_chain"], _residue_title(record["binder_resn"], record["binder_resi"])
        )
        mutation = record.get("mutation") or ""
        mark = "LABEL" if int(record["binder_rank"]) <= labelled else ""
        distance = record.get("best_sg_target_distance_A")
        distance_text = "" if distance is None else "%.2f" % float(distance)
        print("%-5s %4d %-14s %-10s %5s %5.2f %5s %3s %-20s %7s  %s" % (
            mark, int(record["binder_rank"]), site, mutation, record["final_class"],
            float(record["raw_score"]), record["uncapped_class"],
            record.get("class_cap") or "-", record["best_target"], distance_text,
            record.get("downgrade_reasons", "") or "-",
        ))


# ---------------------------------------------------------------------------
# Public PyMOL commands
# ---------------------------------------------------------------------------


def warhead_scan(
    target_chain="B",
    binder_chain="A",
    target_sel="",
    binder_sel="",
    profile="broad_nucleophile",
    reactive_atoms="",
    target_types="",
    interface_cutoff=5.0,
    candidate_cutoff=12.0,
    ideal_min=5.0,
    ideal_max=10.0,
    contact_cutoff=4.0,
    polar_cutoff=3.5,
    salt_bridge_cutoff=4.0,
    internal_contact_cutoff=4.5,
    clash_cutoff=2.0,
    path_radius=1.6,
    object_name="auto",
    state=1,
    out_prefix="warhead_scan",
    top_n=0,
    label_top=3,
    compute_sasa=1,
    compute_rotamers=1,
    include_disulfide_targets=0,
    visualize=1,
    write_outputs=1,
    resn_aliases="",
):
    object_name = _resolve_object(object_name)
    state = _as_int(state, "state")
    top_n = _as_int(top_n, "top_n")
    label_top = max(0, min(3, _as_int(label_top, "label_top")))
    interface_cutoff = _as_float(interface_cutoff, "interface_cutoff")
    candidate_cutoff = _as_float(candidate_cutoff, "candidate_cutoff")
    ideal_min = _as_float(ideal_min, "ideal_min")
    ideal_max = _as_float(ideal_max, "ideal_max")
    contact_cutoff = _as_float(contact_cutoff, "contact_cutoff")
    polar_cutoff = _as_float(polar_cutoff, "polar_cutoff")
    salt_bridge_cutoff = _as_float(salt_bridge_cutoff, "salt_bridge_cutoff")
    internal_contact_cutoff = _as_float(internal_contact_cutoff, "internal_contact_cutoff")
    clash_cutoff = _as_float(clash_cutoff, "clash_cutoff")
    path_radius = _as_float(path_radius, "path_radius")
    compute_sasa = _as_bool(compute_sasa)
    aliases = _parse_alias_option(resn_aliases)
    visualize = _as_bool(visualize)
    write_outputs = _as_bool(write_outputs)
    resolved_profile = _resolve_profile(
        profile, reactive_atoms=reactive_atoms, target_types=target_types, aliases=aliases
    )
    target_selection, binder_selection = _build_scan_selections(
        object_name, target_chain=target_chain, binder_chain=binder_chain,
        target_sel=target_sel, binder_sel=binder_sel, state=state,
    )
    target_atoms = _atoms(target_selection, state)
    binder_atoms = _atoms(binder_selection, state)
    target_res = _group_by_residue(target_atoms)
    binder_res = _group_by_residue(binder_atoms)
    if not target_res:
        raise ValueError("Target selection contains no protein residues: %s" % target_selection)
    if not binder_res:
        raise ValueError("Binder selection contains no protein residues: %s" % binder_selection)

    extra_atoms = _load_hetero_atoms(object_name, state)

    sasa_free = {}
    sasa_bound = {}
    if compute_sasa:
        try:
            sasa_free, bound_all = _compute_sasa_maps(target_selection, binder_selection, state)
            sasa_bound = {key: value for key, value in bound_all.items() if key in binder_res}
            sasa_free = {key: value for key, value in sasa_free.items() if key in binder_res}
        except Exception as exc:
            print("WARNING: residue SASA failed; continuing with modeled SG accessibility only: %s" % exc)

    result = _analyze_groups(
        target_res, binder_res,
        interface_cutoff=interface_cutoff,
        candidate_cutoff=candidate_cutoff,
        ideal_min=ideal_min,
        ideal_max=ideal_max,
        profile=resolved_profile,
        contact_cutoff=contact_cutoff,
        polar_cutoff=polar_cutoff,
        salt_bridge_cutoff=salt_bridge_cutoff,
        internal_contact_cutoff=internal_contact_cutoff,
        clash_cutoff=clash_cutoff,
        path_radius=path_radius,
        include_disulfide_targets=include_disulfide_targets,
        compute_rotamers=compute_rotamers,
        sasa_free=sasa_free,
        sasa_bound=sasa_bound,
        resn_aliases=aliases,
        extra_atoms=extra_atoms,
    )
    result["parameters"].update({
        "object_name": object_name,
        "target_selection": target_selection,
        "binder_selection": binder_selection,
        "state": state,
        "resn_aliases": dict(aliases),
        "n_hetero_atoms": len(extra_atoms),
    })

    print("\n=== Warhead site scan v%s ===" % __version__)
    print("Object: %s" % object_name)
    print("Target selection: %s" % target_selection)
    print("Binder selection: %s" % binder_selection)
    print("Profile: %s | reactive atoms: %s" % (
        resolved_profile.name,
        ";".join("%s:%s" % (key, "|".join(value)) for key, value in sorted(resolved_profile.reactive_atoms.items())),
    ))
    print("Distances: interface %.2f Å | candidate %.2f Å | preferred %.2f-%.2f Å" % (
        float(interface_cutoff), float(candidate_cutoff), float(ideal_min), float(ideal_max)
    ))
    print("Advanced: contact %.2f | polar %.2f | salt %.2f | internal %.2f | clash %.2f Å" % (
        float(contact_cutoff), float(polar_cutoff), float(salt_bridge_cutoff),
        float(internal_contact_cutoff), float(clash_cutoff)
    ))
    print("Non-protein clash/path atoms: %d (organic/metals/inorganic; solvent excluded)" % len(extra_atoms))
    skipped_cys = int((result.get("parameters") or {}).get("skipped_cys_not_interface") or 0)
    if skipped_cys:
        print(
            "Existing Cys skipped (not on target interface, cutoff %.2f Å): %d"
            % (float(interface_cutoff), skipped_cys)
        )
    if aliases:
        print("Residue aliases (user): " + ", ".join(
            "%s->%s" % (key, value) for key, value in sorted(aliases.items())
        ))
    mapped = sorted({
        str(key[2]).upper()
        for key in list(target_res) + list(binder_res)
        if str(key[2]).upper() not in MUTABILITY and _is_known_resn(str(key[2]).upper(), aliases)
    })
    if mapped:
        print("Non-standard names mapped: " + ", ".join(
            "%s->%s" % (name, _canonical_resn(name, aliases)) for name in mapped
        ))
    unknown_t = _unknown_residue_summary(target_res, aliases)
    unknown_b = _unknown_residue_summary(binder_res, aliases)
    if unknown_t or unknown_b:
        merged = {}
        for src in (unknown_t, unknown_b):
            for name, count in src.items():
                merged[name] = merged.get(name, 0) + count
        listed = ", ".join("%sx%d" % (name, count) for name, count in sorted(merged.items()))
        print("WARNING: unrecognised residue name(s): %s" % listed)
        print("         These score with the unknown-residue default. Map them with")
        print("         resn_aliases=SRC:DST if they are standard residues under another name.")

    print("\n=== Target interface nucleophiles ===")
    if not result["target_records"]:
        print("No requested target reactive atom lies inside interface_cutoff.")
    else:
        for record in result["target_records"]:
            status = "available" if record["available_for_pairing"] else "disulfide-bound/excluded"
            site = record.get("target_site") or "%s/%s" % (
                record["target_chain"], _residue_title(record["target_resn"], record["target_resi"])
            )
            print("  %-12s reactive_atom=%-4s reactive_min=%.2f Å  %s" % (
                site, record["nearest_reactive_atom"],
                float(record["reactive_to_binder_min_A"]), status,
            ))

    print("\n=== Binder X->Cys rank ===")
    if result["binder_records"]:
        _print_binder_rank_table(result["binder_records"], top_n=top_n, label_top=label_top)
    else:
        print("No binder site has a modeled/fallback anchor within candidate_cutoff.")

    paths = []
    if write_outputs:
        paths = list(_write_outputs(out_prefix, result))
        print("\nOutputs:")
        for path in paths:
            print("  " + os.path.abspath(path))
    result["output_paths"] = paths
    if visualize:
        _visualize(object_name, result, label_top=label_top)
    return result


def warhead_help():
    print(HELP_TEXT)
    return HELP_TEXT


def warhead_scan_help():
    return warhead_help()


def warhead_profiles():
    for name in sorted(PROFILES):
        profile = PROFILES[name]
        mapping = ";".join("%s:%s" % (key, "|".join(value)) for key, value in sorted(profile.reactive_atoms.items()))
        print("%-20s %-35s %s" % (name, mapping, profile.description))
    print("custom               user-defined via reactive_atoms=...")



def _prompt_objects():
    """Return loaded protein objects in a deterministic order."""
    objects = _object_list("(polymer.protein)")
    if not objects:
        raise ValueError("No protein molecular object is loaded")
    return sorted(set(str(name) for name in objects))


def _prompt_chains(object_name):
    return sorted(set(_available_protein_chains(object_name)), key=_natural_residue_key)


def _default_chain_pair(chains):
    """Choose convenient editable defaults without assuming they are correct."""
    chains = list(chains)
    if not chains:
        return "B", "A"
    target = "B" if "B" in chains else chains[0]
    if "A" in chains and "A" != target:
        binder = "A"
    else:
        binder = next((chain for chain in chains if chain != target), target)
    return target, binder


def _direct_command_example():
    """Build a copyable command for environments where the Qt dialog is unavailable."""
    try:
        objects = _prompt_objects()
        object_name = objects[0]
        target_chain, binder_chain = _default_chain_pair(_prompt_chains(object_name))
    except Exception:
        object_name, target_chain, binder_chain = "complex", "B", "A"
    return (
        "warhead_scan object_name=%s, target_chain=%s, binder_chain=%s, "
        "profile=fluorosulfate, interface_cutoff=5.0, candidate_cutoff=12.0, "
        "ideal_min=5.0, ideal_max=10.0, label_top=3, out_prefix=warhead_scan"
        % (object_name, target_chain, binder_chain)
    )


def _qt_dialog_exec(dialog):
    """Execute a modal dialog across PyQt/PySide versions bundled by PyMOL."""
    execute = getattr(dialog, "exec", None)
    if execute is None:
        execute = getattr(dialog, "exec_", None)
    if execute is None:
        raise RuntimeError("Qt dialog has neither exec() nor exec_()")
    return execute()


def _show_scan_dialog(advanced=False):
    """Show one modal PyMOL Qt form and return warhead_scan keyword arguments.

    Returning ``None`` means that the user pressed Cancel.  This function
    intentionally never calls Python ``input()`` because GUI PyMOL sessions
    commonly have no usable standard-input stream.
    """
    try:
        from pymol.Qt import QtWidgets
    except Exception as exc:
        raise RuntimeError("PyMOL Qt interface is unavailable: %s" % exc)

    objects = _prompt_objects()
    default_object = objects[0]
    default_chains = _prompt_chains(default_object)
    default_target, default_binder = _default_chain_pair(default_chains)

    dialog = QtWidgets.QDialog()
    dialog.setWindowTitle("Warhead Site Scanner v%s" % __version__)
    dialog.setModal(True)
    try:
        dialog.resize(650, 690 if advanced else 560)
    except Exception:
        pass

    outer = QtWidgets.QVBoxLayout(dialog)
    intro = QtWidgets.QLabel(
        "Select target/binder definitions and chemistry. Distances are editable; "
        "5 Å is only a starting value, not a fixed rule."
    )
    try:
        intro.setWordWrap(True)
    except Exception:
        pass
    outer.addWidget(intro)

    form = QtWidgets.QFormLayout()
    outer.addLayout(form)

    object_box = QtWidgets.QComboBox()
    object_box.addItems(objects)
    object_box.setCurrentText(default_object)
    form.addRow("Protein object", object_box)

    chain_hint = QtWidgets.QLabel(", ".join(default_chains) if default_chains else "(none)")
    form.addRow("Available chains", chain_hint)

    target_chain_edit = QtWidgets.QLineEdit(default_target)
    target_chain_edit.setPlaceholderText("B or B+D")
    form.addRow("TARGET chain(s)", target_chain_edit)

    binder_chain_edit = QtWidgets.QLineEdit(default_binder)
    binder_chain_edit.setPlaceholderText("A or A+C")
    form.addRow("BINDER chain(s)", binder_chain_edit)

    target_sel_edit = QtWidgets.QLineEdit()
    target_sel_edit.setPlaceholderText("optional; overrides TARGET chain(s)")
    form.addRow("TARGET selection", target_sel_edit)

    binder_sel_edit = QtWidgets.QLineEdit()
    binder_sel_edit.setPlaceholderText("optional; overrides BINDER chain(s)")
    form.addRow("BINDER selection", binder_sel_edit)

    profile_box = QtWidgets.QComboBox()
    profile_names = sorted(PROFILES)
    profile_box.addItems(profile_names + ["custom"])
    profile_box.setCurrentText("broad_nucleophile")
    form.addRow("Warhead profile", profile_box)

    reactive_edit = QtWidgets.QLineEdit()
    reactive_edit.setPlaceholderText("e.g. CYS:SG;LYS:NZ;HIS:ND1|NE2")
    reactive_edit.setEnabled(False)
    form.addRow("Custom reactive atoms", reactive_edit)

    def double_box(value, maximum=100.0):
        widget = QtWidgets.QDoubleSpinBox()
        widget.setDecimals(2)
        widget.setRange(0.10, maximum)
        widget.setSingleStep(0.50)
        widget.setValue(float(value))
        widget.setSuffix(" Å")
        return widget

    interface_box = double_box(5.0)
    candidate_box = double_box(12.0)
    ideal_min_box = double_box(5.0)
    ideal_max_box = double_box(10.0)
    form.addRow("Target interface cutoff", interface_box)
    form.addRow("Candidate SG cutoff", candidate_box)
    form.addRow("Preferred SG distance min", ideal_min_box)
    form.addRow("Preferred SG distance max", ideal_max_box)

    advanced_boxes = {}
    if advanced:
        advanced_defaults = [
            ("contact_cutoff", "Side-chain contact cutoff", 4.0),
            ("polar_cutoff", "Polar proximity cutoff", 3.5),
            ("salt_bridge_cutoff", "Salt-bridge cutoff", 4.0),
            ("internal_contact_cutoff", "Internal contact cutoff", 4.5),
            ("clash_cutoff", "Modeled SG severe-clash cutoff", 2.0),
            ("path_radius", "Straight-path obstruction radius", 1.6),
        ]
        for key, label, value in advanced_defaults:
            advanced_boxes[key] = double_box(value)
            form.addRow(label, advanced_boxes[key])

    compute_sasa_box = QtWidgets.QCheckBox("Compute residue SASA")
    compute_sasa_box.setChecked(True)
    form.addRow("SASA", compute_sasa_box)

    compute_rotamers_box = QtWidgets.QCheckBox("Model three Cys rotamers")
    compute_rotamers_box.setChecked(True)
    form.addRow("Cys geometry", compute_rotamers_box)

    include_disulfide_box = QtWidgets.QCheckBox("Allow disulfide-bound target Cys")
    include_disulfide_box.setChecked(False)
    form.addRow("Target disulfides", include_disulfide_box)

    label_top_box = QtWidgets.QSpinBox()
    label_top_box.setRange(0, 3)
    label_top_box.setValue(3)
    form.addRow("Label top binder sites", label_top_box)

    aliases_edit = QtWidgets.QLineEdit()
    aliases_edit.setPlaceholderText("optional; e.g. MSE:MET+HIE:HIS or CSO:CYS")
    form.addRow("Residue aliases", aliases_edit)

    out_prefix_edit = QtWidgets.QLineEdit("warhead_scan")
    out_prefix_edit.setPlaceholderText("path/prefix; four output files will be written")
    form.addRow("Output prefix", out_prefix_edit)

    def update_object_context(object_name):
        chains = _prompt_chains(str(object_name))
        chain_hint.setText(", ".join(chains) if chains else "(none)")
        target_default, binder_default = _default_chain_pair(chains)
        target_chain_edit.setText(target_default)
        binder_chain_edit.setText(binder_default)

    object_box.currentTextChanged.connect(update_object_context)
    profile_box.currentTextChanged.connect(
        lambda name: reactive_edit.setEnabled(str(name).lower() == "custom")
    )

    try:
        standard = QtWidgets.QDialogButtonBox.StandardButton
        buttons = QtWidgets.QDialogButtonBox(standard.Ok | standard.Cancel)
    except AttributeError:
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
    outer.addWidget(buttons)

    def validate_and_accept():
        problems = []
        if not object_box.currentText().strip():
            problems.append("Select a protein object.")
        if not target_sel_edit.text().strip() and not target_chain_edit.text().strip():
            problems.append("Provide TARGET chain(s) or an explicit TARGET selection.")
        if not binder_sel_edit.text().strip() and not binder_chain_edit.text().strip():
            problems.append("Provide BINDER chain(s) or an explicit BINDER selection.")
        if ideal_max_box.value() < ideal_min_box.value():
            problems.append("Preferred maximum distance must be >= preferred minimum distance.")
        if profile_box.currentText().lower() == "custom" and not reactive_edit.text().strip():
            problems.append("Custom profile requires reactive atom definitions.")
        if not out_prefix_edit.text().strip():
            problems.append("Output prefix cannot be blank.")
        if problems:
            QtWidgets.QMessageBox.warning(dialog, "Invalid scan parameters", "\n".join(problems))
            return
        dialog.accept()

    buttons.accepted.connect(validate_and_accept)
    buttons.rejected.connect(dialog.reject)

    if not _qt_dialog_exec(dialog):
        return None

    values = {
        "object_name": str(object_box.currentText()).strip(),
        "target_chain": str(target_chain_edit.text()).strip(),
        "binder_chain": str(binder_chain_edit.text()).strip(),
        "target_sel": str(target_sel_edit.text()).strip(),
        "binder_sel": str(binder_sel_edit.text()).strip(),
        "profile": str(profile_box.currentText()).strip(),
        "reactive_atoms": str(reactive_edit.text()).strip(),
        "interface_cutoff": float(interface_box.value()),
        "candidate_cutoff": float(candidate_box.value()),
        "ideal_min": float(ideal_min_box.value()),
        "ideal_max": float(ideal_max_box.value()),
        "compute_sasa": int(compute_sasa_box.isChecked()),
        "compute_rotamers": int(compute_rotamers_box.isChecked()),
        "include_disulfide_targets": int(include_disulfide_box.isChecked()),
        "label_top": int(label_top_box.value()),
        "resn_aliases": str(aliases_edit.text()).strip(),
        "out_prefix": str(out_prefix_edit.text()).strip(),
    }
    for key, widget in advanced_boxes.items():
        values[key] = float(widget.value())
    return values


def warhead_scan_prompt(advanced=0):
    """Open a PyMOL-native Qt parameter dialog and run ``warhead_scan``.

    GUI PyMOL does not guarantee a usable ``sys.stdin``.  Consequently this
    command deliberately uses a Qt form rather than a sequence of ``input()``
    calls.  Pressing Cancel returns ``None`` without running an analysis.
    """
    advanced = _as_bool(advanced)
    print("\nWarhead site scanner v%s" % __version__)
    try:
        gui_caller = getattr(cmd, "_call_in_gui_thread", None)
        if callable(gui_caller):
            values = gui_caller(lambda: _show_scan_dialog(advanced=advanced))
        else:
            values = _show_scan_dialog(advanced=advanced)
    except Exception as exc:
        print("Could not open the PyMOL parameter dialog: %s" % exc)
        print("Run the scanner directly with a command such as:")
        print("  " + _direct_command_example())
        return None
    if values is None:
        print("Warhead scan cancelled.")
        return None
    return warhead_scan(**values)


BATCH_FIELDS = [
    "job_id", "structure", "object_name", "state", "target_chain", "binder_chain",
    "target_sel", "binder_sel", "profile", "reactive_atoms", "target_types",
    "interface_cutoff", "candidate_cutoff", "ideal_min", "ideal_max",
    "contact_cutoff", "polar_cutoff", "salt_bridge_cutoff",
    "internal_contact_cutoff", "clash_cutoff", "path_radius", "compute_sasa",
    "compute_rotamers", "include_disulfide_targets", "resn_aliases", "output_prefix",
]


def _manifest_value(row, name, default=None, cast=None):
    value = row.get(name, "")
    if value is None or str(value).strip() == "":
        return default
    return cast(value) if cast else value


def warhead_batch(
    manifest,
    output_dir="",
    combined_summary="warhead_batch_summary.csv",
    visualize=0,
    keep_objects=0,
    stop_on_error=0,
):
    manifest = os.path.abspath(str(manifest))
    if not os.path.exists(manifest):
        raise ValueError("Manifest not found: %s" % manifest)
    output_dir = os.path.abspath(str(output_dir or os.path.dirname(manifest) or "."))
    os.makedirs(output_dir, exist_ok=True)
    visualize = _as_bool(visualize)
    keep_objects = _as_bool(keep_objects)
    stop_on_error = _as_bool(stop_on_error)
    summary = []
    with open(manifest, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    for index, row in enumerate(rows, 1):
        job_id = str(row.get("job_id") or "job_%03d" % index)
        structure = str(row.get("structure") or "").strip()
        if structure and not os.path.isabs(structure):
            structure = os.path.abspath(os.path.join(os.path.dirname(manifest), structure))
        object_name = str(row.get("object_name") or "").strip()
        if not object_name:
            try:
                object_name = cmd.get_unused_name("warhead_job_", 1)
            except Exception:
                object_name = "warhead_job_%03d" % index
        record = {
            "job_id": job_id, "structure": structure, "object_name": object_name,
            "status": "error", "error": "", "n_targets": 0, "n_pairs": 0,
            "n_binder_sites": 0, "top1_site": "", "top1_mutation": "",
            "top1_class": "", "top1_raw_score": "", "top1_target": "",
        }
        try:
            if not structure:
                raise ValueError("Missing structure path")
            cmd.load(structure, object_name)
            prefix_value = str(row.get("output_prefix") or job_id)
            if not os.path.isabs(prefix_value):
                prefix_value = os.path.join(output_dir, prefix_value)
            kwargs = {
                "object_name": object_name,
                "state": _manifest_value(row, "state", 1, int),
                "target_chain": _manifest_value(row, "target_chain", "B"),
                "binder_chain": _manifest_value(row, "binder_chain", "A"),
                "target_sel": _manifest_value(row, "target_sel", ""),
                "binder_sel": _manifest_value(row, "binder_sel", ""),
                "profile": _manifest_value(row, "profile", "broad_nucleophile"),
                "reactive_atoms": _manifest_value(row, "reactive_atoms", ""),
                "target_types": _manifest_value(row, "target_types", ""),
                "interface_cutoff": _manifest_value(row, "interface_cutoff", 5.0, float),
                "candidate_cutoff": _manifest_value(row, "candidate_cutoff", 12.0, float),
                "ideal_min": _manifest_value(row, "ideal_min", 5.0, float),
                "ideal_max": _manifest_value(row, "ideal_max", 10.0, float),
                "contact_cutoff": _manifest_value(row, "contact_cutoff", 4.0, float),
                "polar_cutoff": _manifest_value(row, "polar_cutoff", 3.5, float),
                "salt_bridge_cutoff": _manifest_value(row, "salt_bridge_cutoff", 4.0, float),
                "internal_contact_cutoff": _manifest_value(row, "internal_contact_cutoff", 4.5, float),
                "clash_cutoff": _manifest_value(row, "clash_cutoff", 2.0, float),
                "path_radius": _manifest_value(row, "path_radius", 1.6, float),
                "compute_sasa": _manifest_value(row, "compute_sasa", 1, int),
                "compute_rotamers": _manifest_value(row, "compute_rotamers", 1, int),
                "include_disulfide_targets": _manifest_value(row, "include_disulfide_targets", 0, int),
                "resn_aliases": _manifest_value(row, "resn_aliases", ""),
                "visualize": int(visualize),
                "out_prefix": prefix_value,
                "label_top": 3 if visualize else 0,
                "top_n": 0,
            }
            result = warhead_scan(**kwargs)
            record.update({
                "status": "ok",
                "n_targets": len(result.get("target_records", [])),
                "n_pairs": len(result.get("pair_records", [])),
                "n_binder_sites": len(result.get("binder_records", [])),
            })
            if result.get("binder_records"):
                top = result["binder_records"][0]
                record.update({
                    "top1_site": top.get("binder_site") or "%s/%s" % (
                        top["binder_chain"],
                        _residue_title(top["binder_resn"], top["binder_resi"]),
                    ),
                    "top1_mutation": top.get("mutation", ""),
                    "top1_class": top["final_class"],
                    "top1_raw_score": top["raw_score"],
                    "top1_target": top["best_target"],
                })
        except Exception as exc:
            record["error"] = "%s: %s" % (exc.__class__.__name__, exc)
            print("ERROR in batch job %s: %s" % (job_id, record["error"]))
            if stop_on_error:
                raise
        finally:
            if not keep_objects:
                try:
                    cmd.delete(object_name)
                except Exception:
                    pass
        summary.append(record)

    if not os.path.isabs(str(combined_summary)):
        combined_summary = os.path.join(output_dir, str(combined_summary))
    fields = [
        "job_id", "structure", "object_name", "status", "error", "n_targets",
        "n_pairs", "n_binder_sites", "top1_site", "top1_mutation", "top1_class",
        "top1_raw_score", "top1_target",
    ]
    _write_csv(combined_summary, summary, fields)
    print("Batch summary: " + os.path.abspath(combined_summary))
    return summary


def warhead_batch_help():
    text = """warhead_batch manifest=/path/jobs.csv, output_dir=/path/results

Required manifest column: structure
Common columns: job_id, object_name, state, target_chain, binder_chain, target_sel,
binder_sel, profile, reactive_atoms, interface_cutoff, candidate_cutoff,
ideal_min, ideal_max, output_prefix.
Batch visualization is off by default. Failed rows are recorded in the combined summary.
"""
    print(text)
    return text


cmd.extend("warhead_scan", warhead_scan)
cmd.extend("warhead_scan_prompt", warhead_scan_prompt)
cmd.extend("warhead_help", warhead_help)
cmd.extend("warhead_scan_help", warhead_scan_help)
cmd.extend("warhead_profiles", warhead_profiles)
cmd.extend("warhead_batch", warhead_batch)
cmd.extend("warhead_batch_help", warhead_batch_help)

print("Loaded warhead-site ranker v%s." % __version__)
print("Help: warhead_help | Profiles: warhead_profiles | Interactive: warhead_scan_prompt")
print("Batch: warhead_batch_help")
