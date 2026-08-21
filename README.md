# Warhead Site Rank

PyMOL 共价 warhead 位点初筛：在蛋白复合物的 **target 界面亲核残基** 上，给 **binder** 排序适合突变为 Cys、安装 linker/warhead 的位点。

这是结构启发式筛选，**不预测**反应速率、共价收率或体内效果。

---

## 最终版本（请用这个）

| | |
|---|---|
| **脚本** | [`warhead_site_rank.py`](warhead_site_rank.py) |
| **版本** | **v3.1.3** |
| **位置** | 仓库根目录 |

另外两份历史脚本已归档，不要当默认入口：

| 归档路径 | 说明 |
|----------|------|
| [`legacy/v2/pymol_warhead_site_rank_v2_1.py`](legacy/v2/pymol_warhead_site_rank_v2_1.py) | v2 线：CB 锚点 + 截短 SASA |
| [`legacy/v3.0.1/warhead_site_rank_v3_0_1.py`](legacy/v3.0.1/warhead_site_rank_v3_0_1.py) | v3.0.1：profile + Cys rotamer 主干 |

最终版以 v3.0.1 为架构，吸收了 v2 的别名、对象解析和可读输出，并补了界面 Cys 门控、配体 clash、短标签等。

---

## 依赖

- PyMOL 2.x / 3.x（命令 `run` 加载）
- 无需额外 pip 包

## 用法

```pymol
load complex.cif, complex
run /path/to/warhead_site_rank.py
warhead_help
warhead_profiles
warhead_scan_prompt
warhead_scan_prompt advanced=1
```

直接扫描：

```pymol
warhead_scan object_name=complex, target_chain=B, binder_chain=A, \
    profile=fluorosulfate, interface_cutoff=5, candidate_cutoff=12, \
    ideal_min=5, ideal_max=10, label_top=3, out_prefix=warhead_scan
```

多链或任意选择：

```pymol
warhead_scan object_name=complex, binder_chain=H+L, target_chain=B
warhead_scan object_name=complex, \
    target_sel="chain B+D and resi 1-200", \
    binder_sel="chain A+C", profile=fluorosulfate
```

力场残基名：

```pymol
warhead_scan ... target_types=CYX+HIE
warhead_scan ... resn_aliases=CSO:CYS
```

批量：

```pymol
warhead_batch manifest=/path/jobs.csv, output_dir=/path/results
```

## Chemistry profiles

| profile | 反应原子 |
|---------|----------|
| `fluorosulfate` | TYR:OH, LYS:NZ, HIS:ND1/NE2 |
| `broad_nucleophile` | CYS:SG, LYS:NZ, TYR:OH, HIS:ND1/NE2 |
| `cys_target` | CYS:SG |
| `custom` | `reactive_atoms="CYS:SG;LYS:NZ;HIS:ND1\|NE2"` |

## 两个主要距离

- `interface_cutoff`：target 反应原子到 binder 重原子，定义亲核残基是否在界面（默认 5 Å，**不是**化学定律）。
- `candidate_cutoff`：建模 Cys-SG 到 target 反应原子，定义安装位点搜索半径（默认 12 Å）。
- `ideal_min` / `ideal_max`：当前 linker/warhead 家族的优选 SG 距离。

Binder 上**已有 Cys** 还要满足：残基重原子到 target 蛋白 ≤ `interface_cutoff`，才会枚举 native + χ1 rotamer。

## 输出

每次扫描写出：

- `<prefix>_targets.csv`
- `<prefix>_pairs.csv`
- `<prefix>_binder_rank.csv`
- `<prefix>_metadata.json`

PyMOL 最多标 3 个 binder 位点：

- binder：`#1 A/Ser74  8.55Å`
- target：`B/Tyr134`
- 等级只用颜色：A 绿 / B 青 / C 品红 / D 灰

等级：A 优先复核，B 可进第一轮，C 备选，D 一般不作为第一轮。

## 公开命令

`warhead_help` · `warhead_profiles` · `warhead_scan` · `warhead_scan_prompt` · `warhead_batch`

## 限制

脚本不枚举完整 maleimide-linker-warhead 构象。前三名仍需人工检查 rotamer、linker 路径、热点、WT/Cys mutant/偶联后亲和力，以及质谱确认交联位点。
