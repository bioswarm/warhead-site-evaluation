"""
PyMOL 共价 warhead 位点筛选脚本 v2.6
======================================

目的
----
该脚本用于蛋白复合物中的初步工程筛选：

1. 用户指定 target 链和 binder 链（各自可以是多条链）；
2. 在 target 界面寻找指定亲核残基，默认 CYS/LYS/TYR/HIS；
3. 以 target 的实际反应原子为中心，寻找 binder 上可改造成 Cys 的位置；
4. 对每个 ``binder residue -> Cys`` 与 ``target nucleophile`` 组合进行启发式排序；
5. 将完整结果打印到 PyMOL 控制台并写出 CSV；
6. PyMOL 画面只显示和标注综合排名最高的三个唯一 binder 位点，避免界面杂乱。

加载与帮助
----------
在 PyMOL 中运行：

    load complex.cif, complex
    run /path/to/pymol_warhead_site_rank_v2.py
    warhead_help
    warhead_scan_prompt

直接模式示例：

    warhead_scan target_chain=B, binder_chain=A, object_name=complex, \
        interface_cutoff=5.0, candidate_cutoff=12.0, \
        ideal_min=5.0, ideal_max=10.0, \
        target_types=CYS+LYS+TYR+HIS, out_prefix=my_scan

``warhead_help`` 会在 PyMOL 控制台重新打印本说明，因此脚本可独立使用，
不依赖外部 README。

适用范围：任意蛋白-蛋白复合物
--------------------------
脚本不绑定任何特定结构，对任意 PDB/mmCIF 复合物都适用。所有参数（链、cutoff、
target 残基类型、输出前缀）都由用户给出，默认值只是示例结构下的常用起点，
**不是**推荐值——换结构时请按体系重新判断。

v2.5 起明确支持以下常见情况：

1. **多链的 target 或 binder**：用 ``+`` 连接，例如抗体 binder 是重链加轻链::

       warhead_scan target_chain=B, binder_chain=H+L

   链 ID 区分大小写，也支持 mmCIF 的多字符 auth ID（如 ``AAA``）。
   两侧链集合不允许重叠，重叠会直接报错而不是给出无意义的结果。

2. **力场/非标准残基名**：AMBER、CHARMM、GROMACS 导出的模型常把组氨酸写成
   ``HIE`` / ``HID`` / ``HSD``，二硫键半胱氨酸写成 ``CYX``，硒代甲硫氨酸写成
   ``MSE``。v2.5 之前这些名字**静默匹配不到任何 target**（等于白跑一遍）。
   现在这类"同原子、同化学"的别名会自动映射到标准母体，评分表照常适用，
   ``target_types`` 也接受这些写法::

       warhead_scan target_chain=B, binder_chain=A, target_types=CYX+HIE

   同时支持 AMBER 的端基写法（``NALA``、``CGLY`` 等），按规则解析。

3. **化学性质已改变的残基不会被静默当成母体**：``CSO`` / ``CSD``（氧化 Cys）、
   ``SEC``（硒半胱氨酸）、``PTR`` / ``SEP`` / ``TPO``（磷酸化）、``MLY``（修饰
   Lys）等一律**不映射**——它们的反应性与母体不同，悄悄改名会误导判断。
   这些残基会出现在"未识别残基"警告里，可用 ``resn_aliases`` 显式指定::

       warhead_scan ... resn_aliases=CSO:CYS+ABC:ALA

4. **未识别残基有明确提示**：凡是评分表不认识的残基名（配体、核酸、修饰残基），
   都会按未知残基默认值处理（mutability -0.5、无参考 SASA），并在控制台列出
   名称和数量，提醒这部分排名偏保守，而不是让你以为一切正常。

输出中同时保留**原始残基名**（``*_resn``，与结构文件一致，可直接用于 PyMOL
选择）和**规范化母体名**（``*_resn_std``，用于评分），两者都写进 CSV。

两个距离参数必须分开：5 Å 不是固定值
----------------------------------
``interface_cutoff`` 定义：

    target 反应原子到 binder 任意重原子的最大距离。

Target 反应原子：

- Cys：SG
- Lys：NZ
- Tyr：OH
- His：ND1 / NE2

``candidate_cutoff`` 定义：

    binder 未来安装锚点到 target 反应原子的最大搜索半径。

Binder 安装锚点：

- 普通残基：CB
- Gly：CA
- 已有 Cys：SG

因此，target 是否位于直接界面，以及带 linker/warhead 的 binder 位点是否可达，
使用两个独立 cutoff；两者都由用户输入，不共用固定 5 Å。

``ideal_min`` 与 ``ideal_max`` 是用户针对具体 linker/warhead 家族设置的
优选锚点距离区间，只用于连续几何评分，不代表通用化学定律。

``contact_cutoff``（默认 4.0 Å）与 ``polar_cutoff``（默认 3.5 Å）控制侧链接触
与极性接触的判定半径，v2.2 起同样可由用户传入，不再写死在代码里：

    warhead_scan ... contact_cutoff=4.5, polar_cutoff=3.2

排序审查内容
------------
脚本综合检查：

- binder CB/SG/CA 到 target 反应原子的距离；
- binder 原侧链与整条 target 链的接触，而不只看目标亲核残基本身；
- binder 原侧链的极性接触和盐桥风险；
- Pro/Gly 主链构象风险；
- 芳香残基和带电残基的界面热点风险；
- target Cys 是否已形成二硫键；
- binder 安装锚原子的自由态 SASA；
- 锚原子在复合物中是否被界面严重埋藏；
- N/C 端位置；
- 原始残基变成 Cys 的经验性可接受度；
- 锚原子的侧链朝向（``anchor_direction_deg``，始终报告，默认不计入总分）。

SASA 使用临时对象计算，不覆盖原结构 B-factor。若 SASA 计算失败，脚本会
打印警告，并以不含 SASA 的较低置信度模式继续运行。

控制台、标签和输出
------------------
完整 binder 排名默认全部打印；``top_n=0`` 表示打印全部，``top_n=N`` 可限制
控制台只打印前 N 名。无论打印多少，CSV 始终保存完整结果。

只在 PyMOL 中标注总排名前三
--------------------------
PyMOL 画面默认只显示综合排名前三的唯一 binder 位点：

- ``wh_top_1``、``wh_top_2``、``wh_top_3``；
- 每个位点只在安装锚原子处放置一个标签；
- 标签用可读的混合大小写残基名，例如::

      #1 A/Ser74 -> Cys  ->  B/Tyr151 OH  |  8.01 A  |  class A

  即"链/残基名+编号"的写法（``Ser74``、``His81``、``Tyr151``），
  binder 位点、建议突变、配对 target 及其反应原子、锚点距离和等级一并显示，
  不用回头查 CSV 就能看懂画面；
- ``wh_top_distance_1`` 至 ``wh_top_distance_3`` 显示各 binder 的最佳 target 配对；
- 其余 binder 候选不显示、不加标签，只在控制台和 CSV 中报告。

``label_top`` 默认为 3，可设为 0、1、2 或 3；脚本会强制限制为最多 3 个，
保证其余 binder 候选只在控制台和 CSV 中报告。``label_top=0`` 时不放任何标签，
此时只对 target 亲核残基做 zoom（v2.2 修正：旧版本会引用未创建的选择对象而报错）。

控制台表格除排名、分数、等级外，还有两列便于阅读和排查：

- ``Binder`` / ``Mutation``：``A/Ser74`` 与 ``A:S74C`` 两种写法并列；
- ``Capped by``：当某位点的原始分数本该给更高等级、却被风险 flag 压低时，
  这一列写出压制它的 flag。排序是先按等级再按分数，所以高分位点可能排在
  低分位点之后；有了这一列就能看出那是刻意降级，而不是排序出错。

输出文件：

- ``<prefix>_targets.csv``：target 界面亲核残基；
- ``<prefix>_pairs.csv``：每个 binder-target 配对的完整评分与风险；
- ``<prefix>_binder_rank.csv``：按 binder 位点汇总的最佳配对排名。

三个 CSV 都新增了可读列（``target_site``、``binder_site``、``mutation``，
形如 ``Tyr151``、``Ser74``、``S74C``）。原有的 ``*_resi`` / ``*_resn`` /
``best_target`` 等列保持不变，已有的下游解析脚本不受影响。

等级解释
--------
- A：优先进行人工结构复核；
- B：值得进入第一轮；
- C：备选或存在明确风险；
- D：通常不建议作为第一轮位置。

等级是工程启发式，不预测共价反应速率，也不能代替：

1. Cys rotamer 检查；
2. 完整 maleimide-linker-warhead 建模；
3. clash 和 linker 路径检查；
4. WT、未修饰 Cys mutant、修饰后 binder 的结合实验；
5. target 位点突变和竞争对照；
6. 完整质谱确认交联位点。

Target 类型与 warhead 化学
-------------------------
默认扫描 CYS/LYS/TYR/HIS，并不表示同一种 warhead 会同等有效地攻击四类残基。
应根据实际化学类型调整，例如：

    warhead_scan ... target_types=TYR+LYS+HIS

或：

    warhead_scan ... target_types=CYS

对上传示例的独立复核
--------------------
示例设置：

- target：B 链
- binder：A 链
- ``interface_cutoff=5.0 Å``
- ``candidate_cutoff=12.0 Å``
- ``ideal_min=5.0 Å``
- ``ideal_max=10.0 Å``

识别到的 target 界面亲核残基：

- ``B:TYR151``，反应原子最近距离约 2.95 Å；
- ``B:HIS81``，反应原子最近距离约 3.11 Å；
- ``B:TYR134``，反应原子最近距离约 4.75 Å。

加入独立 SASA 复核后，``A:SER74`` 是最清晰的首选候选；``A:ALA80`` 虽几何距离合适，
但锚原子在界面中严重埋藏，因此降为 C；``A:ALA73`` 的 CB 在自由 binder 中也埋藏，因此降为 D。
前三名标签为 ``A:S74C``、``A:L4C``、``A:E3C``，在 ``truncate_anchor_sasa``
开启和关闭两种模式下都一致（两种模式下第 4 名之后的排序不同，见 v2.3 一节）。

v2.2 修正与已知取舍
------------------
v2.2 修掉的问题：

1. ``label_top=0`` 时 zoom 引用未创建的 ``wh_top_candidates`` 而报错；
2. N/C 端判断原先取"文件里原子出现的顺序"的首末残基，重建过 loop 或
   拼接过链的结构会标错端点；现按链内残基编号判定，且每条链各自算端点；
3. 排序的残基号 tie-break 原为字符串比较（``10`` 排在 ``2`` 前），现按数值
   排序，插入码作为次级键；
4. ``contact_cutoff`` / ``polar_cutoff`` 原先写死，现已暴露给用户。

当时列为"有意保留"的三项里，第二项（锚原子 SASA 的自遮挡偏差）已在 v2.3
修掉，见下一节。仍然保留的是：

- 锚原子被判埋藏时，``exposure_score`` 与 ``penalty`` 会各扣一次（合计约
  -4.5），同时 ``anchor_buried`` flag 已把等级硬压到 D。等级不会因此再变，
  但 ``pair_score`` 这个数值在"有 SASA"与"无 SASA"的位点之间不可直接比大小；
- 第三项（缺少方向性项）已在 v2.4 作为**可选**功能补上，默认仍关闭，
  见 v2.4 一节。

v2.3：锚原子 SASA 改为截短侧链后计算
----------------------------------
X->Cys 设计真正要问的是"未来那个 SG 够不够得着"，但结构里还带着 X 的完整
侧链，而大侧链会遮住自己的 CB。v2.2 及之前直接用原构象的 CB SASA，因此
Leu/Ile/Lys/Phe 这类残基即使截短成 Cys 后 SG 朝外，也常被判成埋藏。

v2.3 增加 ``truncate_anchor_sasa``（默认 1，即开启）：先用第一遍分析找出
真正落在 ``candidate_cutoff`` 内的候选，再对**每个候选单独**把它自己的侧链
截到 CB（保留 N/CA/C/O/OXT/CB），重新计算该锚原子的自由态与结合态 SASA，
然后用新值重跑评分。

两个要点：

1. **只截候选残基自己的侧链**，邻近残基的侧链一律保留。单点突变不会让邻居
   消失，它们的遮挡必须继续算作风险；
2. Gly（锚点是 CA，本来没有侧链可截）与已有 Cys（SG 本身就是反应基团）
   不做截短，沿用原值。

额外开销只跟候选数成正比，不随链长增长。``truncate_anchor_sasa=0`` 可退回
v2.2 的"原构象"行为。若截短计算失败，脚本打印警告并保留未截短的值。

在示例结构上开启后的实测效果（真实 Shrake-Rupley 复核）：

- 侧链重原子 >= 4 的残基，锚原子自由态 SASA 平均上升约 27 Å²；
- Ala 变化恰好为 0 Å²（CB 之后本来就没有原子可截），是一个干净的对照；
- ``A:LYS78``、``A:LEU77`` 由 D 升为 B，``A:ILE81`` 由 D 升为 C，
  ``A:LEU4`` 由 B 升为 A；
- 前三名 ``Ser74`` / ``Leu4`` / ``Glu3`` 不变，但 ``Lys78`` 从第 16 名升到
  第 4 名，说明原先的排序确实系统性地冤枉了大侧链残基。

需要注意：截短只解决"自遮挡"这一项偏差，**不代表**该位点在突变后一定可用。
Cys rotamer、linker 路径与碰撞仍需人工检查。

v2.4：可选的方向性项（默认关闭）
------------------------------
在此之前评分只有标量距离：8 Å 且背朝 target 的位点，与 8 Å 正对 target 的
位点完全同分。v2.4 增加 ``anchor_direction_deg``，即 CA->CB 与 CB->反应原子
两个向量的夹角：

- 0°：侧链正朝 target 生长；
- 180°：侧链朝反方向生长，SG 得绕过主链才能够到。

物理依据：Cys 的 CA-CB-SG 键角约 114°，因此 SG 偏离 CA->CB 轴约 66°，chi1
旋转时 SG 在这个半角约 66° 的锥面上扫动。所以：

- 夹角 <= 60°：锥内，某个 chi1 一定能把 SG 指向 target，给 +1.5；
- 夹角 >= 120°：已明显超出锥面加上主链让步的余量，给 -1.5；
- 中间线性过渡（90° 恰好为 0）。

两个阈值都比几何极限保守一点，这是刻意的：夹角是在**突变前**的主链上量的，
真实 Cys rotamer 还有转动自由度，所以这一项的权重（±1.5）明显小于几何项
（±3）和暴露项（±2.5），只用来微调排序，不用来拍板。

``anchor_direction_deg`` 与 ``direction_score`` **始终**写进 pairs CSV，方便
自行核对；但只有 ``use_direction_score=1`` 时才计入总分并产生
``anchor_points_toward`` / ``anchor_points_away`` 标记。**默认为 0**，即默认
行为与 v2.3 完全一致::

    warhead_scan ... use_direction_score=1

Gly（锚点即 CA，没有 CA->CB 向量）等几何未定义的情况返回 ``None``，
``_direction_score`` 对 ``None`` 返回 0.0——既不奖励也不惩罚。

各候选实测夹角（示例结构，B 为 target、A 为 binder）：

- 朝向 target 一侧：``A:SER1`` 19°、``A:PRO75`` 24°、``A:TYR83`` 41°、
  ``A:LYS79`` 42°、``A:LYS78`` 53°，均拿到满额 +1.5；
- 中间带：``A:ALA110`` 62°、``A:LEU4`` 67°、``A:ALA73`` 73°、``A:SER74`` 83°、
  ``A:ALA80`` 85°；
- 背朝 target 一侧：``A:ILE81`` 108°、``A:PHE76`` 118° 被扣分，
  ``A:LEU77`` 142° 是唯一触发 ``anchor_points_away`` 的位点。

在**同一套 SASA 数据**下开关对照（这才是可比的口径），17 个位点里有 3 个
等级变化：``A:LEU4`` B->A、``A:GLU82`` C->B、``A:SER1`` D->C；前两名
``Ser74`` / ``Leu4`` 不变，第三名由 ``Glu3``（96°，小幅扣分）换成
``Glu2``（85°）。

也就是说：方向项对**头部**候选影响有限，主要作用是把"几何距离合适但侧链
背朝 target"的位点往下压。是否采用请自行判断——它修正的是"侧链朝向"这一
真实的几何约束，但用的是突变前主链的近似，**不能**替代真正的 Cys rotamer
建模。

校验范围
--------
当前版本已通过：

- Python 语法编译；
- 48 个回归测试（含端点顺序、每链各自端点、残基号数值排序、
  ``label_top=0`` zoom、contact/polar cutoff 生效、可读列、``Capped by`` 列，
  以及截短 SASA 的 ``dot_solvent`` 保护、只截候选自身侧链、Gly/Cys 跳过，
  方向项的 0°/180°/Gly 边界、阈值平台单调性与"默认不计入总分"，
  以及别名映射、"化学已改变的残基不映射"、多链解析与多链 binder 扫描，
  还有 ``get_object_list`` 返回 None / 选择器抛异常 / 多对象歧义三种情况）；
- 在真实 PyMOL 3.1.7.2 上完成端到端验收，详见下一节；
- 用上传 CIF 进行独立坐标解析和距离复核；
- 用独立 Shrake-Rupley SASA 计算复核锚原子暴露度；
- 用模拟 PyMOL ``cmd`` 接口校验前三名标签、距离对象、完整控制台输出和 CSV 逻辑；
- v2.2 的四项修复：以上述 CIF 重放修改后的分析核心，与修改前的 CSV 逐字段
  比对，等级、排名、risk_flags 全部一致，仅有的差异是 CSV 保留三位小数与
  原始浮点数的显示差别；
- v2.3 的截短 SASA：用真实 Shrake-Rupley 实现的 mock ``get_area``
  （会区分 ``dot_solvent`` 0/1）跑开关对照，确认 Ala 变化为 0、大侧链残基
  平均上升约 27 Å²，并确认 ``dot_solvent`` 被还原、临时对象被清理。
  ``dot_solvent`` 保护还做过一次变异测试：故意去掉该行，对应测试确实失败。
- v2.5 的通用化：默认路径再次与最初的 baseline CSV 逐字段比对，17 个位点
  **0 处差异**；另外把示例结构的 HIS/CYS/MET 批量改名为 HIE/CYX/MSE 后重跑，
  target 数、pair 数与全部排名与改名前**完全一致**，证明别名映射被彻底吸收；
- 多链：把 binder 链 A 拆成 H（resi<=60）与 L（resi>60）后用 ``binder_chain=H+L``
  端到端跑通，H 与 L 的候选都出现在排名里，同号不同链的残基不会混淆；
  链缺失、链重叠、空链参数三种错误输入都给出明确的 ValueError。

v2.6：在真实 PyMOL 中验证（PyMOL 3.1.7.2）
----------------------------------------
此前所有验证都靠模拟 ``cmd`` 接口，docstring 里也一直写着"建议再做一次真实
PyMOL smoke test"。现已在真实 PyMOL 3.1.7.2（``pymol -cq`` 无头模式）跑完，
结论如下。

**发现并修掉一个只有真实 PyMOL 才会暴露的 bug**：``cmd.get_object_list()`` 在
选择匹配不到任何东西时返回 ``None`` 而不是空列表，``len(None)`` 直接抛
``TypeError``，把"对象找不到"变成一个看不懂的崩溃。触发场景很常见——PyMOL
载入时会自己改对象名：

- 保留字加下划线：``load x.cif, alt`` 实际得到 ``alt_``；
- 非法字符被替换：``load x.cif, my complex`` 实际得到 ``my_complex``。

用户按自己写的名字去调用就会踩到。现在 ``_object_list`` 统一把 ``None`` 和
选择器异常都归一化成空列表，并给出可操作的报错（列出实际可用的对象名，
并说明 PyMOL 会改名这件事）。``object_name="auto"`` 在没有载入任何结构时
也改成"请先 load 结构"而不是"找不到名为 auto 的对象"。

**真实 PyMOL 复核的结果**（与模拟环境一致，无回归）：

- 默认跑通：17 个位点、3 个标签、3 个距离对象，前三名仍是
  ``A:S74C`` / ``A:L4C`` / ``A:E3C``；
- 标签文字在视口中即为::

      #1 A/Ser74 -> Cys  ->  B/Tyr151 OH  |  8.01 A  |  class A

- 三个标注位点的选择表达式都能选到真实原子（6 / 8 / 9 个）；
- ``label_top`` 取 3/2/1/0 全部不报错，``label_top=0`` 时不再引用
  ``wh_top_candidates``（v2.2 修复在真实环境确认）；
- ``dot_solvent`` / ``dot_density`` 用后还原（``cmd.get`` 返回的是
  ``'off'`` / ``'2'`` 这样的字符串，回写同一字符串可正常还原）；
- 临时对象无泄漏，原结构 B-factor 未被 SASA 覆盖（A/74 CB 仍为 81.49）；
- **截短 SASA 用真实 ``get_area`` 复核**：侧链重原子 >= 4 的残基平均上升
  27.3 Å²，Ala 恰好 0.00 Å²，等级变化与模拟结果完全相同
  （Leu4 B->A、Lys78 D->B、Leu77 D->B、Ile81 D->C、Ser1 D->C）；
- **别名映射**：用 ``cmd.alter`` 把 HIS/CYS/MET 就地改名为 HIE/CYX/MSE 后，
  target 数、pair 数与全部排名与改名前完全一致；即使 ``target_types`` 用默认
  的标准写法也照样找到这些 target；
- **多链**：用 ``cmd.alter`` 把 A 链拆成 H/L 后 ``binder_chain=H+L`` 跑通，
  排名与"A 作为单链"完全一致，链重叠/链缺失都给出明确 ValueError；
- **插入码**：``resi "52A"`` 的选择表达式在真实解析器下能正确选到 9 个原子
  （引号写法有效），带插入码的结构可完整跑完；
- **altloc**：同一残基 A/B 双构象共 12 个原子，去重后保留 6 个且无重名，
  保留的是占据率较高的 A 构象（q=0.65）；配对不重复计数；
- **多 state**：同一对象 2 个 state 分别扫描，结果一致。

限制
----
排序是初筛，不是自由能、反应速率或共价收率计算。前三名仍需逐个检查 Cys
rotamer、完整 linker 构象、碰撞、结合热点和实验亲和力。
"""

HELP_TEXT = (__doc__ or "").strip()
SCRIPT_HELP = HELP_TEXT
WARHEAD_SCAN_HELP = HELP_TEXT
__version__ = "2.6"
from pymol import cmd
from collections import defaultdict
from math import sqrt
import csv
import os
import re


TARGET_REACTIVE = {
    "CYS": ("SG",),
    "LYS": ("NZ",),
    "TYR": ("OH",),
    "HIS": ("ND1", "NE2"),
}

# Mutation-to-Cys prior. Positive means generally easier to test; negative
# means greater structural/interface risk. These values are deliberately
# heuristic and are only one component of the final rank.
MUTABILITY = {
    "ALA": 3.0, "SER": 2.5, "THR": 1.8, "ASN": 1.5, "GLN": 1.2,
    "VAL": 0.5, "ILE": 0.0, "LEU": 0.0, "MET": 0.3,
    "ASP": -0.8, "GLU": -0.8, "LYS": -1.0, "ARG": -1.2,
    "HIS": -1.0, "PHE": -1.5, "TYR": -1.8, "TRP": -2.0,
    "GLY": -2.2, "PRO": -3.5, "CYS": 3.5,
}

POLAR_ATOMS = {
    "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"}, "CYS": {"SG"},
    "ASN": {"OD1", "ND2"}, "GLN": {"OE1", "NE2"},
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"},
    "LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"},
    "HIS": {"ND1", "NE2"}, "TRP": {"NE1"},
}

POSITIVE_ATOMS = {
    "LYS": {"NZ"},
    "ARG": {"NE", "NH1", "NH2"},
    # Histidine protonation is context-dependent, so it is not treated as a
    # definite positive charge in the automatic salt-bridge count.
}
NEGATIVE_ATOMS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}

AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
CHARGED = {"ASP", "GLU", "LYS", "ARG"}
BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT"}
CLASS_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3}

# Theoretical maximum residue ASA values from the Tien et al. scale. They are
# used only to turn PyMOL's approximate per-residue SASA into a relative value.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0,
    "CYS": 167.0, "GLN": 225.0, "GLU": 223.0, "GLY": 104.0,
    "HIS": 224.0, "ILE": 197.0, "LEU": 201.0,
    "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0,
    "VAL": 174.0,
}

AA_ONE_LETTER = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


# Residue-name aliases seen in real structures. Only protonation states, charge
# states and isosteric substitutions map onto a standard parent: those share the
# parent's atom names and reactive chemistry, so every scoring table applies
# unchanged.
#
# Deliberately NOT mapped, because the chemistry differs from the parent and a
# silent rename would mislead:
#   CSO/CSD/OCS/CME/CSX  oxidised or adducted Cys -- no longer nucleophilic
#   SEC                  selenocysteine -- different pKa and reactivity
#   PTR/SEP/TPO          phosphorylated Tyr/Ser/Thr -- blocked hydroxyl
#   MLY/ALY/M3L          modified Lys -- amine blocked or altered
# These stay unrecognised, are reported by the unknown-residue summary, and can
# be scored explicitly by the user via extra_aliases if that is what they want.
RESN_ALIASES = {
    # Histidine protonation / naming variants (AMBER, CHARMM, GROMACS).
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    # Cysteine charge / disulfide-state variants. CYX is a disulfide-bonded Cys;
    # the geometric disulfide check still runs and will flag it.
    "CYX": "CYS", "CYM": "CYS",
    # Neutral/protonated forms of ionisable residues.
    "LYN": "LYS", "ASH": "ASP", "GLH": "GLU", "AR0": "ARG",
    # Selenomethionine is isosteric with Met and extremely common in PDB entries.
    "MSE": "MET",
    # Alternative Ile CD naming does not change the residue identity.
    "ILE": "ILE",
}

# AMBER/CHARMM terminal residues are the standard name with an N or C prefix
# (NALA, CGLY, ...). Handled by rule rather than by enumerating 40 extra keys.
_TERMINAL_PREFIXES = ("N", "C")


def _unknown_residue_summary(residue_groups, aliases=None):
    """Count residues whose names the scoring tables do not recognise.

    Returns ``{resn: count}``. These are scored with the unknown-residue default
    (mutability -0.5, no MAX_ASA), so a structure full of them silently produces
    a pessimistic ranking. Reporting them lets the user decide whether to supply
    an alias or exclude the residues.
    """
    unknown = defaultdict(int)
    for key in residue_groups:
        resn = str(key[2] or "").upper().strip()
        if not _is_known_resn(resn, aliases):
            unknown[resn] += 1
    return dict(unknown)


def _canonical_resn(resn, extra_aliases=None):
    """Map a residue name onto its standard 3-letter parent for scoring.

    Returns the name unchanged when it is already standard or is not recognised;
    callers use ``_is_known_resn`` to tell those two cases apart. The ORIGINAL
    name must still be used for PyMOL selections and CSV output, because that is
    what the loaded structure actually contains.
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
    # Terminal-variant rule, e.g. NALA/CALA -> ALA. Guarded so that a genuine
    # 4-character ligand code is not mangled into an amino acid.
    if len(name) == 4 and name[0] in _TERMINAL_PREFIXES:
        stem = name[1:]
        if stem in MUTABILITY:
            return stem
        stem_alias = RESN_ALIASES.get(stem)
        if stem_alias:
            return stem_alias
    return name


def _is_known_resn(resn, extra_aliases=None):
    """True when the residue can be scored with the standard tables."""
    return _canonical_resn(resn, extra_aliases) in MUTABILITY


def _parse_alias_option(value):
    """Parse ``MSE:MET+HIE:HIS`` style user aliases into a dict.

    Lets a user score a residue the tables do not know (or override a built-in
    mapping) without editing the script.
    """
    if not value:
        return {}
    if isinstance(value, dict):
        return {str(k).upper().strip(): str(v).upper().strip()
                for k, v in value.items()}
    out = {}
    for item in re.split(r"[,+;\s]+", str(value).upper()):
        if not item:
            continue
        if ":" not in item:
            raise ValueError(
                f"Alias '{item}' must be written SOURCE:TARGET, e.g. MSE:MET"
            )
        src, dst = item.split(":", 1)
        src, dst = src.strip(), dst.strip()
        if dst not in MUTABILITY:
            raise ValueError(
                f"Alias target '{dst}' is not a standard residue name"
            )
        out[src] = dst
    return out


def _residue_title(resn, resi):
    """Readable residue tag such as ``His81`` or ``Ser74``.

    Used for PyMOL labels and console output. CSV columns keep the raw
    upper-case ``resn``/``resi`` fields so existing downstream parsers that
    read the CSVs are unaffected.
    """
    resn = str(resn or "").upper()
    return f"{resn.capitalize()}{resi}"


def _resi_sort_key(resi):
    """Numeric-aware residue-id key so 2 < 9 < 10 < 100.

    A plain string sort puts '10' before '2'. Insertion codes are kept as a
    secondary text key, and anything unparseable sorts last but stays stable.
    """
    match = re.match(r"^\s*(-?\d+)\s*(.*)$", str(resi))
    if match:
        return (0, int(match.group(1)), match.group(2))
    return (1, 0, str(resi))


def _chain_termini(residue_groups):
    """Return (n_terminal_keys, c_terminal_keys) per chain by residue number.

    Residue order must not be taken from the file's atom order: rebuilt loops,
    spliced chains and edited CIFs can list residues out of sequence, which
    would otherwise flag the wrong residues as termini.
    """
    per_chain = defaultdict(list)
    for key in residue_groups:
        per_chain[key[0]].append(key)
    nterm, cterm = set(), set()
    for keys in per_chain.values():
        ordered = sorted(keys, key=lambda k: _resi_sort_key(k[1]))
        nterm.add(ordered[0])
        cterm.add(ordered[-1])
    return nterm, cterm


def _d(a, b):
    return sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _atom_element(atom):
    elem = (getattr(atom, "elem", "") or "").upper().strip()
    if elem:
        return elem
    name = (getattr(atom, "name", "") or "").strip()
    return name[0].upper() if name else ""


def _preferred_alt_atoms(atoms):
    """Keep one conformer per atom identity: blank > A > highest occupancy."""
    chosen = {}
    for atom in atoms:
        key = (
            getattr(atom, "model", ""), getattr(atom, "segi", ""),
            getattr(atom, "chain", ""), str(getattr(atom, "resi", "")),
            getattr(atom, "resn", ""), getattr(atom, "name", ""),
        )
        alt = (getattr(atom, "alt", "") or "").strip()
        occ = float(getattr(atom, "q", 0.0) or 0.0)
        priority = 3 if alt == "" else (2 if alt == "A" else 1)
        old = chosen.get(key)
        if old is None or (priority, occ) > old[0]:
            chosen[key] = ((priority, occ), atom)
    return [item[1] for item in chosen.values()]


def _atoms(selection, state=1):
    return _preferred_alt_atoms(list(cmd.get_model(selection, int(state)).atom))


def _reskey(atom):
    # A scan is restricted to one molecular object, so chain/resi/resn is
    # sufficient and keeps insertion codes intact in the resi string.
    return (getattr(atom, "chain", ""), str(getattr(atom, "resi", "")), getattr(atom, "resn", ""))


def _group_by_residue(atoms):
    out = defaultdict(list)
    for atom in atoms:
        out[_reskey(atom)].append(atom)
    return out


def _heavy(atoms):
    return [a for a in atoms if _atom_element(a) != "H"]


def _sidechain(atoms):
    return [a for a in _heavy(atoms) if getattr(a, "name", "") not in BACKBONE_NAMES]


def _is_polar_atom(atom, aliases=None):
    name = getattr(atom, "name", "")
    if name in {"N", "O", "OXT"}:
        return True
    resn = _canonical_resn(getattr(atom, "resn", ""), aliases)
    return name in POLAR_ATOMS.get(resn, set())


def _is_positive_atom(atom, aliases=None):
    resn = _canonical_resn(getattr(atom, "resn", ""), aliases)
    return getattr(atom, "name", "") in POSITIVE_ATOMS.get(resn, set())


def _is_negative_atom(atom, aliases=None):
    resn = _canonical_resn(getattr(atom, "resn", ""), aliases)
    return getattr(atom, "name", "") in NEGATIVE_ATOMS.get(resn, set())


def _min_distance(atoms1, atoms2):
    best = 999.0
    pair = (None, None)
    for a in atoms1:
        for b in atoms2:
            dd = _d(a.coord, b.coord)
            if dd < best:
                best = dd
                pair = (a, b)
    return best, pair


def _contact_pairs(atoms1, atoms2, cutoff):
    n = 0
    for a in _heavy(atoms1):
        for b in _heavy(atoms2):
            if _d(a.coord, b.coord) <= cutoff:
                n += 1
    return n


def _polar_contact_pairs(atoms1, atoms2, cutoff, aliases=None):
    # Distance-only proxy; hydrogen placement and angular geometry are not
    # inferred here. Only binder side-chain atoms should be passed as atoms1.
    n = 0
    for a in _heavy(atoms1):
        if not _is_polar_atom(a, aliases):
            continue
        for b in _heavy(atoms2):
            if _is_polar_atom(b, aliases) and _d(a.coord, b.coord) <= cutoff:
                n += 1
    return n


def _salt_bridge_pairs(atoms1, atoms2, cutoff=4.0, aliases=None):
    n = 0
    for a in _heavy(atoms1):
        for b in _heavy(atoms2):
            opposite = (
                _is_positive_atom(a, aliases) and _is_negative_atom(b, aliases)
            ) or (
                _is_negative_atom(a, aliases) and _is_positive_atom(b, aliases)
            )
            if opposite and _d(a.coord, b.coord) <= cutoff:
                n += 1
    return n


def _unique_contact_residues(probe_atoms, residue_groups, cutoff):
    count = 0
    labels = []
    if not probe_atoms:
        return count, labels
    for key, atoms in residue_groups.items():
        md, _ = _min_distance(_heavy(probe_atoms), _heavy(atoms))
        if md <= cutoff:
            count += 1
            labels.append(f"{key[0]}:{key[2]}{key[1]}")
    return count, labels


def _reactive_atoms(res_atoms, resn, aliases=None):
    wanted = TARGET_REACTIVE.get(_canonical_resn(resn, aliases), ())
    return [a for a in res_atoms if getattr(a, "name", "") in wanted]


def _anchor_atom(res_atoms, resn, aliases=None):
    """Return the structural origin used for linker-reach screening."""
    heavy = _heavy(res_atoms)
    if _canonical_resn(resn, aliases) == "CYS":
        sg = next((a for a in heavy if getattr(a, "name", "") == "SG"), None)
        if sg is not None:
            return sg
    cb = next((a for a in heavy if getattr(a, "name", "") == "CB"), None)
    if cb is not None:
        return cb
    return next((a for a in heavy if getattr(a, "name", "") == "CA"), None)


def _anchor_direction_angle(res_atoms, resn, anchor, reactive_atom):
    """Angle in degrees between CA->CB and CB->reactive-atom.

    A crude but honest proxy for "would the installed Cys point at the target".
    0 deg means the side chain grows straight toward the reactive atom; 180 deg
    means it grows directly away, so the SG would have to reach around the
    backbone.

    Returns None when the geometry is undefined: Gly (anchor is CA, so there is
    no CA->CB vector), a missing CA, or degenerate coordinates. Callers must
    treat None as "no information", never as a good or bad angle.
    """
    if anchor is None or reactive_atom is None:
        return None
    heavy = _heavy(res_atoms)
    ca = next((a for a in heavy if getattr(a, "name", "") == "CA"), None)
    if ca is None or anchor is ca:
        return None
    # For an existing Cys the anchor is SG; CA->CB remains the right stem.
    stem_end = anchor
    if getattr(anchor, "name", "") == "SG":
        cb = next((a for a in heavy if getattr(a, "name", "") == "CB"), None)
        if cb is None:
            return None
        stem_end = cb

    def _sub(p, q):
        return (p[0] - q[0], p[1] - q[1], p[2] - q[2])

    stem = _sub(stem_end.coord, ca.coord)
    reach = _sub(reactive_atom.coord, stem_end.coord)
    stem_len = sqrt(stem[0] ** 2 + stem[1] ** 2 + stem[2] ** 2)
    reach_len = sqrt(reach[0] ** 2 + reach[1] ** 2 + reach[2] ** 2)
    if stem_len < 1e-6 or reach_len < 1e-6:
        return None
    dot = sum(stem[i] * reach[i] for i in range(3)) / (stem_len * reach_len)
    dot = max(-1.0, min(1.0, dot))
    from math import acos, degrees
    return degrees(acos(dot))


def _direction_score(angle_deg):
    """Map the anchor direction angle onto a -1.5..+1.5 adjustment.

    Deliberately gentler than the geometry and exposure terms: the angle is
    computed on the pre-mutation backbone and a real Cys rotamer has torsional
    freedom, so this nudges the ranking rather than deciding it. Returns 0.0
    when the angle is unavailable, so sites without geometry are neither
    rewarded nor punished.
    """
    if angle_deg is None:
        return 0.0
    if angle_deg <= 60.0:
        return 1.5
    if angle_deg >= 120.0:
        return -1.5
    # Linear between the two plateaus.
    return 1.5 - 3.0 * (angle_deg - 60.0) / 60.0


def _find_disulfide_residues(*residue_group_sets, **kwargs):
    aliases = kwargs.get("aliases")
    cys = []
    for groups in residue_group_sets:
        for key, atoms in groups.items():
            if _canonical_resn(key[2], aliases) != "CYS":
                continue
            sg = next((a for a in atoms if getattr(a, "name", "") == "SG"), None)
            if sg is not None:
                cys.append((key, sg))
    bound = set()
    for i, (key1, sg1) in enumerate(cys):
        for key2, sg2 in cys[i + 1:]:
            dd = _d(sg1.coord, sg2.coord)
            if 1.7 <= dd <= 2.4:
                bound.add(key1)
                bound.add(key2)
    return bound


def _parse_target_types(value, aliases=None):
    """Normalise the requested target types onto canonical residue names.

    Accepts force-field spellings (``HIE``, ``HSD``, ``CYX``) and maps them to
    the standard parent, so a user working from an AMBER-derived model can pass
    the names their file actually contains.
    """
    if value is None:
        return set(TARGET_REACTIVE)
    if isinstance(value, (list, tuple, set)):
        parts = [str(x).upper() for x in value]
    else:
        parts = [x for x in re.split(r"[,+;/\s]+", str(value).upper()) if x]
    result = set()
    unknown = []
    for part in parts:
        canonical = _canonical_resn(part, aliases)
        if canonical in TARGET_REACTIVE:
            result.add(canonical)
        else:
            unknown.append(part)
    if unknown:
        raise ValueError(
            "Unsupported target residue type(s): " + ", ".join(sorted(unknown))
            + ". Supported: " + ", ".join(sorted(TARGET_REACTIVE))
            + " (force-field aliases such as HIE/HSD/CYX are accepted; use "
            "resn_aliases=SRC:DST to map anything else)"
        )
    if not result:
        raise ValueError("At least one target residue type is required")
    return result


def _geometry_score(distance, ideal_min, ideal_max, candidate_cutoff):
    """Continuous, user-parameterized geometry score from -2 to +3."""
    if ideal_min <= distance <= ideal_max:
        return 3.0
    if distance < ideal_min:
        lower = min(2.5, ideal_min - 0.1)
        if distance <= lower:
            return -2.0
        frac = (distance - lower) / max(ideal_min - lower, 0.1)
        return -2.0 + 5.0 * frac
    if candidate_cutoff <= ideal_max:
        return 0.0
    frac = (candidate_cutoff - distance) / (candidate_cutoff - ideal_max)
    return max(0.0, min(3.0, 3.0 * frac))


def _rank_label(score, risk_flags, resn):
    flags = set(risk_flags)
    if "disulfide_bound_binder_cys" in flags or "anchor_buried" in flags or resn == "PRO":
        return "D"
    if score >= 4.5:
        label = "A"
    elif score >= 2.5:
        label = "B"
    elif score >= 0.5:
        label = "C"
    else:
        label = "D"
    # These are not absolute exclusions, but they make a bulky Cys conjugate
    # sufficiently risky that the site should not be presented as first-line.
    if ("anchor_low_exposure" in flags or "anchor_interface_buried" in flags or resn == "GLY") and label in {"A", "B"}:
        return "C"
    return label


def _analyze_groups(
    target_res,
    binder_res,
    interface_cutoff=5.0,
    candidate_cutoff=12.0,
    ideal_min=5.0,
    ideal_max=10.0,
    target_types=None,
    contact_cutoff=4.0,
    polar_cutoff=3.5,
    include_disulfide_targets=False,
    sasa_free=None,
    sasa_bound=None,
    atom_sasa_free=None,
    atom_sasa_bound=None,
    use_direction_score=False,
    resn_aliases=None,
):
    """Pure analysis core used by PyMOL and regression tests."""
    interface_cutoff = float(interface_cutoff)
    candidate_cutoff = float(candidate_cutoff)
    ideal_min = float(ideal_min)
    ideal_max = float(ideal_max)
    contact_cutoff = float(contact_cutoff)
    polar_cutoff = float(polar_cutoff)
    include_disulfide_targets = bool(int(include_disulfide_targets)) if isinstance(include_disulfide_targets, str) else bool(include_disulfide_targets)
    use_direction_score = bool(int(use_direction_score)) if isinstance(use_direction_score, str) else bool(use_direction_score)
    aliases = _parse_alias_option(resn_aliases)
    target_types = _parse_target_types(target_types, aliases)

    if interface_cutoff <= 0 or candidate_cutoff <= 0:
        raise ValueError("Distance cutoffs must be positive")
    if ideal_min <= 0 or ideal_max < ideal_min:
        raise ValueError("Require 0 < ideal_min <= ideal_max")

    all_target_heavy = _heavy([a for atoms in target_res.values() for a in atoms])
    all_binder_heavy = _heavy([a for atoms in binder_res.values() for a in atoms])
    if not all_target_heavy:
        raise ValueError("Target selection contains no heavy atoms")
    if not all_binder_heavy:
        raise ValueError("Binder selection contains no heavy atoms")

    disulfide_bound = _find_disulfide_residues(target_res, binder_res, aliases=aliases)
    sasa_free = sasa_free or {}
    sasa_bound = sasa_bound or {}
    atom_sasa_free = atom_sasa_free or {}
    atom_sasa_bound = atom_sasa_bound or {}

    target_records = []
    target_context = []
    for tkey, tatoms in target_res.items():
        tch, tresi, tresn = tkey
        # Match on the canonical name so HIE/HSD/CYX targets are found, but keep
        # the original name in every record and PyMOL selection.
        tresn_std = _canonical_resn(tresn, aliases)
        if tresn_std not in target_types:
            continue
        r_atoms = _reactive_atoms(tatoms, tresn, aliases)
        if not r_atoms:
            continue
        reactive_min, reactive_pair = _min_distance(r_atoms, all_binder_heavy)
        residue_min, _ = _min_distance(_heavy(tatoms), all_binder_heavy)
        if reactive_min > interface_cutoff:
            continue
        is_disulfide = tkey in disulfide_bound and tresn_std == "CYS"
        available = (not is_disulfide) or include_disulfide_targets
        rec = {
            "target_chain": tch,
            "target_resi": tresi,
            "target_resn": tresn,
            "target_resn_std": tresn_std,
            "target_site": _residue_title(tresn, tresi),
            "reactive_atoms": ",".join(a.name for a in r_atoms),
            "nearest_reactive_atom": reactive_pair[0].name if reactive_pair[0] else "",
            "nearest_binder_atom": reactive_pair[1].name if reactive_pair[1] else "",
            "reactive_to_binder_min_A": reactive_min,
            "residue_to_binder_min_A": residue_min,
            "disulfide_bound": is_disulfide,
            "available_for_pairing": available,
        }
        target_records.append(rec)
        target_context.append((tkey, tatoms, r_atoms, rec))

    # Termini are resolved per chain by residue number, not by file atom order.
    nterm_keys, cterm_keys = _chain_termini(binder_res)

    pair_records = []
    for tkey, tatoms, r_atoms, target_rec in target_context:
        if not target_rec["available_for_pairing"]:
            continue
        tch, tresi, tresn = tkey
        for bkey, batoms in binder_res.items():
            bch, bresi, bresn = bkey
            bresn_std = _canonical_resn(bresn, aliases)
            anchor = _anchor_atom(batoms, bresn, aliases)
            if anchor is None:
                continue
            anchor_distance, anchor_pair = _min_distance([anchor], r_atoms)
            if anchor_distance > candidate_cutoff:
                continue

            sidechain_atoms = _sidechain(batoms)
            min_any_reactive, any_reactive_pair = _min_distance(_heavy(batoms), r_atoms)
            min_residue_distance, _ = _min_distance(_heavy(batoms), _heavy(tatoms))

            sidechain_contact_pairs = _contact_pairs(sidechain_atoms, all_target_heavy, contact_cutoff) if sidechain_atoms else 0
            sidechain_polar_contacts = _polar_contact_pairs(sidechain_atoms, all_target_heavy, polar_cutoff, aliases) if sidechain_atoms else 0
            salt_bridge_contacts = _salt_bridge_pairs(sidechain_atoms, all_target_heavy, 4.0, aliases) if sidechain_atoms else 0
            contacted_target_count, contacted_targets = _unique_contact_residues(sidechain_atoms, target_res, contact_cutoff)

            other_binder = {k: v for k, v in binder_res.items() if k != bkey}
            internal_contact_count, _ = _unique_contact_residues(sidechain_atoms or [anchor], other_binder, 4.5)

            free_asa = sasa_free.get(bkey)
            bound_asa = sasa_bound.get(bkey)
            max_asa = MAX_ASA.get(bresn_std)
            rsa_free = None
            interface_buried_A2 = None
            interface_buried_fraction = None
            if free_asa is not None and max_asa:
                rsa_free = max(0.0, free_asa / max_asa)
            if free_asa is not None and bound_asa is not None:
                interface_buried_A2 = max(0.0, free_asa - bound_asa)
                if free_asa > 1e-6:
                    interface_buried_fraction = interface_buried_A2 / free_asa

            anchor_key = bkey + (anchor.name,)
            anchor_sasa_free = atom_sasa_free.get(anchor_key)
            anchor_sasa_bound = atom_sasa_bound.get(anchor_key)
            anchor_buried_fraction = None
            if anchor_sasa_free is not None and anchor_sasa_bound is not None and anchor_sasa_free > 1e-6:
                anchor_buried_fraction = max(0.0, anchor_sasa_free - anchor_sasa_bound) / anchor_sasa_free

            mut_score = MUTABILITY.get(bresn_std, -0.5)
            geo_score = _geometry_score(anchor_distance, ideal_min, ideal_max, candidate_cutoff)

            # Direction is reported either way so the CSV can be inspected; it
            # only enters the score when the user opts in.
            anchor_angle = _anchor_direction_angle(
                batoms, bresn, anchor, anchor_pair[1]
            )
            direction_score = (
                _direction_score(anchor_angle) if use_direction_score else 0.0
            )

            exposure_score = 0.0
            if anchor_sasa_free is not None:
                # Smooth scoring avoids unstable rank jumps from the small
                # numerical variation inherent to dot-sampled SASA.
                if anchor_sasa_free < 1.0:
                    exposure_score = -2.5
                elif anchor_sasa_free < 5.0:
                    exposure_score = -0.5 + (anchor_sasa_free - 1.0) * 0.25
                elif anchor_sasa_free < 15.0:
                    exposure_score = 0.5 + (anchor_sasa_free - 5.0) * 0.10
                else:
                    exposure_score = 1.5
            elif rsa_free is not None:
                if rsa_free >= 0.50:
                    exposure_score = 1.5
                elif rsa_free >= 0.25:
                    exposure_score = 1.0
                elif rsa_free >= 0.10:
                    exposure_score = 0.0
                else:
                    exposure_score = -1.5

            penalty = 0.0
            penalty += min(sidechain_contact_pairs * 0.08, 1.5)
            penalty += min(contacted_target_count * 0.25, 1.0)
            penalty += min(sidechain_polar_contacts * 1.4, 3.5)
            penalty += min(salt_bridge_contacts * 1.8, 3.6)
            if bresn_std in AROMATIC:
                penalty += 0.7
            if bresn_std in CHARGED:
                penalty += 0.4
            if bresn_std == "PRO":
                penalty += 2.5
            elif bresn_std == "GLY":
                penalty += 1.5
            if anchor_sasa_free is not None:
                if anchor_sasa_free < 1.0:
                    penalty += 2.0
                elif anchor_sasa_free < 5.0:
                    penalty += 1.0
            elif rsa_free is not None and rsa_free < 0.10:
                penalty += 1.0
            elif rsa_free is None and internal_contact_count >= 8:
                penalty += 0.6
            if anchor_buried_fraction is not None:
                if anchor_buried_fraction >= 0.80:
                    penalty += 1.0
                elif anchor_buried_fraction >= 0.50:
                    penalty += 0.5
            elif interface_buried_fraction is not None:
                if interface_buried_fraction >= 0.60:
                    penalty += 0.8
                elif interface_buried_fraction >= 0.30:
                    penalty += 0.4

            flags = []
            if bresn_std == "PRO":
                flags.append("pro_backbone")
            if bresn_std == "GLY":
                flags.append("gly_backbone")
            if bresn_std in AROMATIC:
                flags.append("aromatic_hotspot")
            if bresn_std in CHARGED:
                flags.append("charged_residue")
            if sidechain_polar_contacts:
                flags.append("polar_contact")
            if salt_bridge_contacts:
                flags.append("salt_bridge")
            if sidechain_contact_pairs >= 6:
                flags.append("dense_interface_contacts")
            if anchor_sasa_free is not None:
                if anchor_sasa_free < 1.0:
                    flags.append("anchor_buried")
                elif anchor_sasa_free < 5.0:
                    flags.append("anchor_low_exposure")
                if anchor_buried_fraction is not None and anchor_buried_fraction >= 0.80:
                    flags.append("anchor_interface_buried")
            elif rsa_free is not None and rsa_free < 0.10:
                flags.append("buried_in_binder")
            if interface_buried_fraction is not None and interface_buried_fraction >= 0.50:
                flags.append("strong_interface_burial")
            if bkey in nterm_keys:
                flags.append("n_terminus")
            if bkey in cterm_keys:
                flags.append("c_terminus")
            if bresn_std == "CYS":
                if bkey in disulfide_bound:
                    flags.append("disulfide_bound_binder_cys")
                    penalty += 6.0
                else:
                    flags.append("existing_free_cys_check")
            if anchor_distance < ideal_min:
                flags.append("tighter_than_preferred")
            elif anchor_distance > ideal_max:
                flags.append("longer_than_preferred")

            if anchor_angle is not None and use_direction_score:
                if anchor_angle >= 120.0:
                    flags.append("anchor_points_away")
                elif anchor_angle <= 60.0:
                    flags.append("anchor_points_toward")

            final_score = (
                mut_score + geo_score + exposure_score + direction_score - penalty
            )
            rank_class = _rank_label(final_score, flags, bresn_std)

            pair_records.append({
                "binder_chain": bch,
                "binder_resi": bresi,
                "binder_resn": bresn,
                "binder_site": _residue_title(bresn, bresi),
                "binder_resn_std": bresn_std,
                "mutation": f"{AA_ONE_LETTER.get(bresn_std, bresn)}{bresi}C",
                "anchor_atom": anchor.name,
                "target_chain": tch,
                "target_resi": tresi,
                "target_resn": tresn,
                "target_site": _residue_title(tresn, tresi),
                "target_atom": anchor_pair[1].name if anchor_pair[1] else "",
                "anchor_distance_A": anchor_distance,
                "min_any_binder_to_reactive_A": min_any_reactive,
                "min_binder_to_target_residue_A": min_residue_distance,
                "mutability_score": mut_score,
                "geometry_score": geo_score,
                "exposure_score": exposure_score,
                "anchor_direction_deg": anchor_angle,
                "direction_score": direction_score,
                "hotspot_penalty": penalty,
                "pair_score": final_score,
                "rank_class": rank_class,
                "sidechain_contact_pairs": sidechain_contact_pairs,
                "contacted_target_residues": contacted_target_count,
                "contacted_target_labels": ",".join(contacted_targets),
                "sidechain_polar_contacts": sidechain_polar_contacts,
                "salt_bridge_contacts": salt_bridge_contacts,
                "internal_contact_residues": internal_contact_count,
                "free_residue_sasa_A2": free_asa,
                "bound_residue_sasa_A2": bound_asa,
                "free_relative_sasa": rsa_free,
                "interface_buried_sasa_A2": interface_buried_A2,
                "interface_buried_fraction": interface_buried_fraction,
                "anchor_sasa_free_A2": anchor_sasa_free,
                "anchor_sasa_bound_A2": anchor_sasa_bound,
                "anchor_buried_fraction": anchor_buried_fraction,
                "risk_flags": ";".join(flags),
                "nearest_binder_atom_to_reactive": any_reactive_pair[0].name if any_reactive_pair[0] else "",
            })

    pair_records.sort(key=lambda r: (CLASS_ORDER[r["rank_class"]], -r["pair_score"], r["anchor_distance_A"], r["binder_chain"], _resi_sort_key(r["binder_resi"])))
    for i, rec in enumerate(pair_records, 1):
        rec["pair_rank"] = i

    aggregate = {}
    for rec in pair_records:
        key = (rec["binder_chain"], rec["binder_resi"], rec["binder_resn"])
        agg = aggregate.get(key)
        # Single source of truth: _best_pair_for_binder matches on this exact
        # string, so both sites must build it through the same helper.
        target_label = _pair_target_label(rec)
        if agg is None:
            aggregate[key] = {
                "binder_chain": rec["binder_chain"],
                "binder_resi": rec["binder_resi"],
                "binder_resn": rec["binder_resn"],
                "binder_site": rec["binder_site"],
                "binder_resn_std": rec["binder_resn_std"],
                "mutation": rec["mutation"],
                "best_pair_score": rec["pair_score"],
                "rank_class": rec["rank_class"],
                "best_target": target_label,
                "best_anchor_distance_A": rec["anchor_distance_A"],
                "risk_flags": rec["risk_flags"],
                "target_options": [target_label],
            }
        else:
            agg["target_options"].append(target_label)
            if rec["pair_score"] > agg["best_pair_score"]:
                agg.update({
                    "best_pair_score": rec["pair_score"],
                    "rank_class": rec["rank_class"],
                    "best_target": target_label,
                    "best_anchor_distance_A": rec["anchor_distance_A"],
                    "risk_flags": rec["risk_flags"],
                })

    binder_records = list(aggregate.values())
    for rec in binder_records:
        rec["target_options"] = ",".join(sorted(set(rec["target_options"])))
    binder_records.sort(key=lambda r: (CLASS_ORDER[r["rank_class"]], -r["best_pair_score"], r["best_anchor_distance_A"], r["binder_chain"], _resi_sort_key(r["binder_resi"])))
    for i, rec in enumerate(binder_records, 1):
        rec["binder_rank"] = i

    target_records.sort(key=lambda r: (r["reactive_to_binder_min_A"], r["target_chain"], _resi_sort_key(r["target_resi"])))
    return {
        "target_records": target_records,
        "pair_records": pair_records,
        "binder_records": binder_records,
    }


def _quote_selector(value):
    value = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{value}"'


def _model_clause(object_name):
    return f"model {_quote_selector(object_name)}"


def _chain_clause(chain):
    return f"chain {_quote_selector(chain)}"


def _parse_chain_spec(value):
    """Split a chain option into a list of chain IDs.

    Accepts ``B``, ``H+L``, ``H,L`` or a list. Multi-chain groups matter in
    practice: an antibody binder is usually H plus L, and a target can span
    several chains. Chain IDs are case-sensitive and may be multi-character
    (mmCIF auth IDs like ``AAA``), so only the separators are interpreted.
    """
    if isinstance(value, (list, tuple, set)):
        parts = [str(x).strip() for x in value]
    else:
        parts = [x.strip() for x in re.split(r"[,+;/\s]+", str(value))]
    seen, out = set(), []
    for part in parts:
        if part and part not in seen:
            seen.add(part)
            out.append(part)
    if not out:
        raise ValueError("At least one chain ID is required")
    return out


def _chains_clause(chains):
    """Selection clause matching any of the given chain IDs."""
    chain_list = chains if isinstance(chains, (list, tuple)) else [chains]
    inner = " or ".join(_chain_clause(c) for c in chain_list)
    return f"({inner})"


def _residue_selection(object_name, chain, resi, resn=None):
    parts = [_model_clause(object_name), _chain_clause(chain), f"resi {_quote_selector(resi)}"]
    if resn:
        parts.append(f"resn {resn}")
    return "(" + " and ".join(parts) + ")"



def _pair_target_label(rec):
    """Machine-facing target tag used in CSVs and as the aggregate match key.

    Kept in upper-case ``B:TYR151(OH)`` form for backward compatibility; use
    ``_pair_target_title`` for anything a human reads.
    """
    return (
        f"{rec['target_chain']}:{rec['target_resn']}"
        f"{rec['target_resi']}({rec['target_atom']})"
    )


def _pair_target_title(rec):
    """Readable target tag such as ``B/Tyr151 OH``."""
    return (
        f"{rec['target_chain']}/"
        f"{_residue_title(rec['target_resn'], rec['target_resi'])}"
        f" {rec['target_atom']}"
    )


def _binder_site_title(rec):
    """Readable binder tag such as ``A/Ser74``."""
    return (
        f"{rec['binder_chain']}/"
        f"{_residue_title(rec['binder_resn'], rec['binder_resi'])}"
    )


def _best_pair_for_binder(result, binder_rec):
    """Return the pair record selected as this binder site's aggregate best."""
    key = (
        str(binder_rec["binder_chain"]),
        str(binder_rec["binder_resi"]),
        str(binder_rec["binder_resn"]),
    )
    fallback = None
    for pair in result.get("pair_records", []):
        pair_key = (
            str(pair["binder_chain"]),
            str(pair["binder_resi"]),
            str(pair["binder_resn"]),
        )
        if pair_key != key:
            continue
        if fallback is None:
            fallback = pair
        if _pair_target_label(pair) == binder_rec.get("best_target"):
            return pair
    return fallback


def _object_list(selection):
    """cmd.get_object_list that always returns a list.

    Real PyMOL returns None (not an empty list) when the selection matches
    nothing, which turned a clear "object not found" message into a TypeError.
    """
    try:
        found = cmd.get_object_list(selection)
    except Exception:
        # An unparsable/empty selection is a "no match", not a crash.
        return []
    return list(found) if found else []


def _resolve_object(object_name="auto"):
    object_name = str(object_name or "auto").strip()
    if object_name.lower() in {"auto", "none", ""}:
        objects = _object_list("(polymer.protein)")
    else:
        objects = _object_list(f"({_model_clause(object_name)})")
        if not objects and object_name in (cmd.get_names("objects") or []):
            objects = [object_name]
    if len(objects) != 1:
        available = ", ".join(_object_list("(polymer.protein)")) or "none"
        auto_mode = object_name.lower() in {"auto", "none", ""}
        if not objects:
            if auto_mode:
                raise ValueError(
                    "No protein object is loaded. Load a structure first, e.g. "
                    "load complex.cif, complex"
                )
            raise ValueError(
                f"Molecular object '{object_name}' not found or contains no "
                f"protein. Available protein objects: {available}. "
                "Note PyMOL rewrites object names on load: reserved keywords "
                "gain an underscore ('alt' -> 'alt_') and spaces or other "
                "invalid characters are replaced ('my complex' -> "
                "'my_complex'). Use the name PyMOL actually assigned."
            )
        raise ValueError(
            "Select exactly one molecular object; matched "
            f"{len(objects)}: {', '.join(objects)}. "
            f"Available protein objects: {available}"
        )
    return objects[0]


def _sasa_area_maps(atoms):
    residue_out = defaultdict(float)
    atom_out = {}
    for atom in _heavy(atoms):
        area = max(0.0, float(getattr(atom, "b", 0.0) or 0.0))
        rkey = _reskey(atom)
        residue_out[rkey] += area
        atom_out[rkey + (atom.name,)] = area
    return dict(residue_out), atom_out


def _truncatable_to_cb(resn, anchor_name, aliases=None):
    """True when a residue's side chain can be cut back to CB for SASA.

    Only meaningful when the installation anchor is CB. Gly (anchor CA) has no
    side chain to remove, and an existing Cys already carries the SG that would
    do the reacting, so neither is truncated.
    """
    canonical = _canonical_resn(resn, aliases)
    return anchor_name == "CB" and canonical not in {"GLY", "CYS"}


def _truncated_anchor_sasa(binder_key, anchor_name, target_sel, binder_sel, state=1):
    """Anchor SASA measured after cutting *this* residue back to CB.

    For an X->Cys design the question is whether the future SG is reachable,
    but the deposited structure still carries X's full side chain, which
    shields its own CB. Leu/Ile/Lys/Phe therefore look buried even where a
    truncated Cys would face solvent.

    Only the candidate residue is truncated. Neighbouring side chains are left
    in place on purpose: their shielding survives a single point mutation and
    must keep counting against the site.

    Returns ``(free_area, bound_area)``; either may be None if the anchor atom
    is absent from the rebuilt object.
    """
    chain, resi, resn = binder_key
    residue_clause = (
        f"{_chain_clause(chain)} and resi {_quote_selector(resi)}"
        f" and resn {resn}"
    )
    # Everything past CB; the backbone and CB itself stay.
    strip_clause = f"({residue_clause}) and not name N+CA+C+O+OXT+CB"

    free_area = bound_area = None
    tmp_free = cmd.get_unused_name("_warhead_trunc_free", 1)
    tmp_bound = cmd.get_unused_name("_warhead_trunc_bound", 1)
    # dot_solvent must be 1 here for the same reason as in _compute_sasa_maps:
    # with dot_solvent=0 get_area returns molecular surface, not solvent
    # accessible area, and the truncated values would not be comparable with
    # the untruncated ones they replace.
    old_dot_solvent = cmd.get("dot_solvent")
    old_dot_density = cmd.get("dot_density")
    try:
        cmd.set("dot_solvent", 1)
        cmd.set("dot_density", 3)
        cmd.create(tmp_free, f"({binder_sel})", int(state), 1)
        cmd.create(tmp_bound, f"({target_sel}) or ({binder_sel})", int(state), 1)
        for name in (tmp_free, tmp_bound):
            cmd.remove(f"({name}) and hydro")
            cmd.remove(f"({name}) and {strip_clause}")
            cmd.get_area(name, state=1, load_b=1)
        wanted = (chain, str(resi), resn, anchor_name)
        for name, slot in ((tmp_free, "free"), (tmp_bound, "bound")):
            _, atom_map = _sasa_area_maps(_atoms(name, 1))
            value = atom_map.get(wanted)
            if slot == "free":
                free_area = value
            else:
                bound_area = value
    finally:
        cmd.delete(tmp_free)
        cmd.delete(tmp_bound)
        try:
            cmd.set("dot_solvent", old_dot_solvent)
            cmd.set("dot_density", old_dot_density)
        except Exception:
            pass
    return free_area, bound_area


def _compute_truncated_anchor_maps(candidates, target_sel, binder_sel, state=1,
                                   aliases=None):
    """Replace CB anchor SASA with truncated-side-chain values.

    ``candidates`` is an iterable of ``(binder_key, anchor_atom_name)``. Gly and
    existing Cys anchors are skipped and keep their untruncated values.
    """
    free_out, bound_out = {}, {}
    for binder_key, anchor_name in candidates:
        if not _truncatable_to_cb(binder_key[2], anchor_name, aliases):
            continue
        free_area, bound_area = _truncated_anchor_sasa(
            binder_key, anchor_name, target_sel, binder_sel, state
        )
        atom_key = binder_key + (anchor_name,)
        if free_area is not None:
            free_out[atom_key] = free_area
        if bound_area is not None:
            bound_out[atom_key] = bound_area
    return free_out, bound_out


def _compute_sasa_maps(target_sel, binder_sel, state=1):
    """Return binder-alone and complex-bound residue SASA maps.

    PyMOL get_area is a discrete approximation. Temporary objects are used so
    the original B factors are not overwritten.
    """
    tmp_bound = cmd.get_unused_name("_warhead_bound", 1)
    tmp_free = cmd.get_unused_name("_warhead_free", 1)
    old_dot_solvent = cmd.get("dot_solvent")
    old_dot_density = cmd.get("dot_density")
    try:
        cmd.create(tmp_bound, f"({target_sel}) or ({binder_sel})", int(state), 1)
        cmd.create(tmp_free, f"({binder_sel})", int(state), 1)
        cmd.remove(f"({tmp_bound}) and hydro")
        cmd.remove(f"({tmp_free}) and hydro")
        cmd.set("dot_solvent", 1)
        cmd.set("dot_density", 3)
        cmd.get_area(tmp_bound, state=1, load_b=1)
        cmd.get_area(tmp_free, state=1, load_b=1)
        bound_atoms = _atoms(tmp_bound, 1)
        free_atoms = _atoms(tmp_free, 1)
        free_residue, free_atom = _sasa_area_maps(free_atoms)
        bound_residue, bound_atom = _sasa_area_maps(bound_atoms)
        return free_residue, bound_residue, free_atom, bound_atom
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


def _fmt(value, digits=3):
    if value is None:
        return ""
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return value


def _write_csv(path, records, fields):
    _ensure_parent(path)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for rec in records:
            writer.writerow({field: _fmt(rec.get(field)) for field in fields})


def _write_outputs(prefix, result):
    prefix = str(prefix or "warhead_scan")
    if prefix.lower().endswith(".csv"):
        prefix = prefix[:-4]
    targets_path = prefix + "_targets.csv"
    pairs_path = prefix + "_pairs.csv"
    binder_path = prefix + "_binder_rank.csv"

    _write_csv(targets_path, result["target_records"], [
        "target_chain", "target_site",
        "target_resi", "target_resn", "target_resn_std", "reactive_atoms",
        "nearest_reactive_atom", "nearest_binder_atom",
        "reactive_to_binder_min_A", "residue_to_binder_min_A",
        "disulfide_bound", "available_for_pairing",
    ])
    _write_csv(pairs_path, result["pair_records"], [
        "pair_rank", "rank_class", "pair_score",
        "binder_chain", "binder_site", "mutation",
        "binder_resi", "binder_resn", "binder_resn_std", "anchor_atom",
        "target_chain", "target_site",
        "target_resi", "target_resn", "target_atom",
        "anchor_distance_A", "min_any_binder_to_reactive_A",
        "min_binder_to_target_residue_A", "mutability_score",
        "geometry_score", "exposure_score",
        "anchor_direction_deg", "direction_score", "hotspot_penalty",
        "sidechain_contact_pairs", "contacted_target_residues",
        "contacted_target_labels", "sidechain_polar_contacts",
        "salt_bridge_contacts", "internal_contact_residues",
        "free_residue_sasa_A2", "bound_residue_sasa_A2",
        "free_relative_sasa", "interface_buried_sasa_A2",
        "interface_buried_fraction", "anchor_sasa_free_A2",
        "anchor_sasa_bound_A2", "anchor_buried_fraction", "risk_flags",
    ])
    _write_csv(binder_path, result["binder_records"], [
        "binder_rank", "rank_class", "best_pair_score", "binder_chain",
        "binder_site", "mutation",
        "binder_resi", "binder_resn", "binder_resn_std", "best_target",
        "best_anchor_distance_A", "target_options", "risk_flags",
    ])
    return targets_path, pairs_path, binder_path


def _records_for_console(records, top_n=0):
    """Return all records when top_n <= 0, otherwise the requested prefix."""
    limit = int(top_n)
    return list(records) if limit <= 0 else list(records)[:limit]


def _class_capped_by(rec):
    """Name the risk flag that held this site below its raw-score class.

    The table sorts by class first and score second, so a high-scoring site
    that was capped appears below lower-scoring ones. Without this column that
    looks like a sorting bug rather than an intentional downgrade.
    """
    flags = set(str(rec.get("risk_flags", "")).split(";")) - {""}
    resn = str(rec.get("binder_resn", "")).upper()
    score = float(rec.get("best_pair_score", 0.0))
    raw = "A" if score >= 4.5 else ("B" if score >= 2.5 else ("C" if score >= 0.5 else "D"))
    actual = str(rec.get("rank_class", raw))
    if CLASS_ORDER.get(actual, 0) <= CLASS_ORDER.get(raw, 0):
        return ""
    for flag in ("disulfide_bound_binder_cys", "anchor_buried",
                 "anchor_interface_buried", "anchor_low_exposure"):
        if flag in flags:
            return flag
    if resn in {"PRO", "GLY"}:
        return f"{resn.lower()}_backbone"
    return "risk_flag"


def _print_binder_rank_table(records, top_n=0, label_top=3):
    """Print the complete binder ranking by default; mark the three labels."""
    print("Mark  Rank  Binder      Mutation      Score Class "
          "Capped by                Best target        AnchorDist  Risk flags")
    labelled = max(0, min(3, int(label_top)))
    for rec in _records_for_console(records, top_n):
        binder = _binder_site_title(rec)
        mutation = _mutation_site_label(rec)
        mark = "LABEL" if int(rec["binder_rank"]) <= labelled else ""
        print(
            f"{mark:<5} {int(rec['binder_rank']):>4}  {binder:<11} {mutation:<13} "
            f"{float(rec['best_pair_score']):>5.2f}   {rec['rank_class']:^5} "
            f"{_class_capped_by(rec):<24} {rec['best_target']:<18} "
            f"{float(rec['best_anchor_distance_A']):>7.2f} Å   {rec.get('risk_flags', '')}"
        )


def _mutation_site_label(rec):
    """Compact mutation tag such as ``A:S74C`` (one-letter, for tables).

    Uses the canonical parent name when present, so MSE74 reads as M74C rather
    than falling back to the raw three-letter code.
    """
    resn = str(rec.get("binder_resn_std") or rec["binder_resn"]).upper()
    one = AA_ONE_LETTER.get(resn, str(rec["binder_resn"]))
    return f"{rec['binder_chain']}:{one}{rec['binder_resi']}C"


def _mutation_site_title(rec):
    """Readable mutation tag such as ``A/Ser74 -> Cys``."""
    return f"{_binder_site_title(rec)} -> Cys"


def _visualize(object_name, result, top_n=None, label_top=3):
    """Show exactly the three highest-ranked unique binder sites.

    ``top_n`` is accepted for backward compatibility with v2 calls, but it
    only controls console output in ``warhead_scan``. It never increases the
    number of PyMOL labels beyond three.
    """
    del top_n  # Legacy visualization argument; intentionally ignored.
    label_top = max(0, min(3, int(label_top)))

    cleanup_names = [
        "wh_target_nucleophiles", "wh_candidates_A", "wh_candidates_B",
        "wh_candidates_C", "wh_candidates_D", "wh_top_candidates",
        "warhead_scan_group",
    ]
    cleanup_names.extend(f"wh_top_{i}" for i in range(1, 21))
    cleanup_names.extend(f"wh_top_label_{i}" for i in range(1, 21))
    cleanup_names.extend(f"wh_top_distance_{i}" for i in range(1, 21))
    cleanup_names.extend(f"wh_distance_{i}" for i in range(1, 21))
    for name in cleanup_names:
        cmd.delete(name)

    target_parts = [
        _residue_selection(
            object_name, r["target_chain"], r["target_resi"], r["target_resn"]
        )
        for r in result.get("target_records", [])
    ]
    if target_parts:
        cmd.select("wh_target_nucleophiles", " or ".join(target_parts))
        cmd.show("sticks", "wh_target_nucleophiles")
        cmd.color("yellow", "wh_target_nucleophiles")

    top_records = result.get("binder_records", [])[:label_top]
    top_parts = []
    present = []
    class_colors = {"A": "green", "B": "cyan", "C": "magenta", "D": "gray70"}

    for i, binder_rec in enumerate(top_records, 1):
        pair = _best_pair_for_binder(result, binder_rec)
        if pair is None:
            continue

        binder_sel = _residue_selection(
            object_name,
            binder_rec["binder_chain"],
            binder_rec["binder_resi"],
            binder_rec["binder_resn"],
        )
        target_sel = _residue_selection(
            object_name,
            pair["target_chain"],
            pair["target_resi"],
            pair["target_resn"],
        )
        anchor_sel = f"({binder_sel}) and name {pair['anchor_atom']}"
        reactive_sel = f"({target_sel}) and name {pair['target_atom']}"

        top_name = f"wh_top_{i}"
        label_name = f"wh_top_label_{i}"
        distance_name = f"wh_top_distance_{i}"
        cmd.select(top_name, binder_sel)
        cmd.show("sticks", top_name)
        cmd.color(class_colors.get(binder_rec["rank_class"], "gray70"), top_name)

        distance_value = float(
            pair.get("anchor_distance_A", binder_rec.get("best_anchor_distance_A", 0.0))
        )
        # Readable mixed-case residue names (Ser74 / Tyr151) rather than the
        # upper-case CSV form, so the 3D scene is legible without the tables.
        label_text = (
            f"#{i} {_mutation_site_title(binder_rec)}"
            f"  ->  {_pair_target_title(pair)}"
            f"  |  {distance_value:.2f} A  |  class {binder_rec['rank_class']}"
        )
        cmd.pseudoatom(label_name, selection=anchor_sel, label=label_text)
        try:
            cmd.hide("everything", label_name)
            cmd.show("labels", label_name)
            cmd.set("label_size", 16, label_name)
        except Exception:
            # Old/open-source PyMOL builds differ slightly in display methods;
            # the pseudoatom label itself is still valid.
            pass
        cmd.distance(distance_name, anchor_sel, reactive_sel)

        top_parts.append(binder_sel)
        present.extend([top_name, label_name, distance_name])

    if top_parts:
        cmd.select("wh_top_candidates", " or ".join(top_parts))
        present.append("wh_top_candidates")
    if target_parts:
        present.insert(0, "wh_target_nucleophiles")

    if present:
        cmd.group("warhead_scan_group", " ".join(present))
        # Only zoom on selections that were actually created. With label_top=0
        # no binder selection exists, and naming it here would make PyMOL fail
        # to parse the zoom expression.
        zoom_names = []
        if target_parts:
            zoom_names.append("wh_target_nucleophiles")
        if top_parts:
            zoom_names.append("wh_top_candidates")
        if zoom_names:
            cmd.zoom(" or ".join(zoom_names), 8)

def warhead_scan(
    target_chain="B",
    binder_chain="A",
    interface_cutoff=5.0,
    candidate_cutoff=12.0,
    ideal_min=5.0,
    ideal_max=10.0,
    target_types="CYS+LYS+TYR+HIS",
    object_name="auto",
    state=1,
    out_prefix="warhead_scan",
    top_n=0,
    label_top=3,
    compute_sasa=1,
    include_disulfide_targets=0,
    visualize=1,
    contact_cutoff=4.0,
    polar_cutoff=3.5,
    truncate_anchor_sasa=1,
    use_direction_score=0,
    resn_aliases="",
):
    """Scan target nucleophiles and rank binder X->Cys installation sites."""
    interface_cutoff = float(interface_cutoff)
    candidate_cutoff = float(candidate_cutoff)
    ideal_min = float(ideal_min)
    ideal_max = float(ideal_max)
    contact_cutoff = float(contact_cutoff)
    polar_cutoff = float(polar_cutoff)
    state = int(state)
    top_n = int(top_n)
    label_top = max(0, min(3, int(label_top)))
    compute_sasa = int(compute_sasa)
    include_disulfide_targets = int(include_disulfide_targets)
    visualize = int(visualize)
    truncate_anchor_sasa = int(truncate_anchor_sasa)
    use_direction_score = int(use_direction_score)
    aliases = _parse_alias_option(resn_aliases)

    object_name = _resolve_object(object_name)
    # Either side may span several chains, e.g. an antibody binder is H plus L.
    target_chains = _parse_chain_spec(target_chain)
    binder_chains = _parse_chain_spec(binder_chain)
    overlap = sorted(set(target_chains) & set(binder_chains))
    if overlap:
        raise ValueError(
            "Target and binder chains must not overlap; both contain: "
            + ", ".join(overlap)
        )

    base_sel = f"({_model_clause(object_name)}) and polymer.protein"
    chains = list(cmd.get_chains(base_sel, state))
    missing = [c for c in target_chains + binder_chains if c not in chains]
    if missing:
        raise ValueError(
            "Chain(s) not found: " + ", ".join(missing)
            + f". Available protein chains: {', '.join(chains) or 'none'}"
            + " (chain IDs are case-sensitive; pass several as H+L)"
        )

    target_sel = f"({base_sel}) and {_chains_clause(target_chains)}"
    binder_sel = f"({base_sel}) and {_chains_clause(binder_chains)}"
    target_res = _group_by_residue(_atoms(target_sel, state))
    binder_res = _group_by_residue(_atoms(binder_sel, state))

    sasa_free = {}
    sasa_bound = {}
    atom_sasa_free = {}
    atom_sasa_bound = {}
    if compute_sasa:
        try:
            sasa_free, all_bound, atom_sasa_free, all_atom_bound = _compute_sasa_maps(target_sel, binder_sel, state)
            # The bound maps contain both chains; retain binder keys only.
            sasa_bound = {key: value for key, value in all_bound.items() if key in binder_res}
            atom_sasa_bound = {key: value for key, value in all_atom_bound.items() if key[:3] in binder_res}
            atom_sasa_free = {key: value for key, value in atom_sasa_free.items() if key[:3] in binder_res}
        except Exception as exc:
            print(f"WARNING: SASA calculation failed; continuing without SASA: {exc}")
            sasa_free = {}
            sasa_bound = {}
            atom_sasa_free = {}
            atom_sasa_bound = {}

    def _run(free_atom_map, bound_atom_map):
        return _analyze_groups(
            target_res=target_res,
            binder_res=binder_res,
            interface_cutoff=interface_cutoff,
            candidate_cutoff=candidate_cutoff,
            ideal_min=ideal_min,
            ideal_max=ideal_max,
            target_types=target_types,
            contact_cutoff=contact_cutoff,
            polar_cutoff=polar_cutoff,
            include_disulfide_targets=include_disulfide_targets,
            sasa_free=sasa_free,
            sasa_bound=sasa_bound,
            atom_sasa_free=free_atom_map,
            atom_sasa_bound=bound_atom_map,
            use_direction_score=use_direction_score,
            resn_aliases=aliases,
        )

    result = _run(atom_sasa_free, atom_sasa_bound)

    # Second pass: re-measure the anchor SASA of each candidate with that
    # residue's own side chain cut back to CB, which is the geometry the
    # X->Cys design actually creates. Only candidates from the first pass are
    # re-measured, so the extra get_area cost scales with hits, not chain size.
    truncated_sites = 0
    if compute_sasa and truncate_anchor_sasa and result["pair_records"]:
        candidates = []
        seen = set()
        for rec in result["pair_records"]:
            bkey = (rec["binder_chain"], rec["binder_resi"], rec["binder_resn"])
            if bkey in seen:
                continue
            seen.add(bkey)
            candidates.append((bkey, rec["anchor_atom"]))
        try:
            trunc_free, trunc_bound = _compute_truncated_anchor_maps(
                candidates, target_sel, binder_sel, state, aliases
            )
            if trunc_free or trunc_bound:
                merged_free = dict(atom_sasa_free)
                merged_free.update(trunc_free)
                merged_bound = dict(atom_sasa_bound)
                merged_bound.update(trunc_bound)
                result = _run(merged_free, merged_bound)
                truncated_sites = len(trunc_free)
        except Exception as exc:
            print(
                "WARNING: CB-truncated anchor SASA failed; keeping untruncated "
                f"values: {exc}"
            )

    print("\n=== Scan definition ===")
    print(f"Object: {object_name} | target chain(s): {'+'.join(target_chains)}"
          f" | binder chain(s): {'+'.join(binder_chains)}")
    print(f"Target reactive-atom interface cutoff: {interface_cutoff:.2f} Å")
    print(f"Binder anchor candidate cutoff:       {candidate_cutoff:.2f} Å")
    print(f"Preferred anchor interval:            {ideal_min:.2f}-{ideal_max:.2f} Å")
    print(f"Side-chain contact cutoff:            {contact_cutoff:.2f} Å")
    print(f"Polar contact cutoff:                 {polar_cutoff:.2f} Å")
    print(f"Target residue types:                 {','.join(sorted(_parse_target_types(target_types, aliases)))}")
    if aliases:
        print("Residue aliases (user):               "
              + ", ".join(f"{k}->{v}" for k, v in sorted(aliases.items())))
    applied = sorted({
        str(k[2]).upper() for k in list(target_res) + list(binder_res)
        if str(k[2]).upper() not in MUTABILITY
        and _is_known_resn(str(k[2]).upper(), aliases)
    })
    if applied:
        print("Non-standard names mapped:            "
              + ", ".join(f"{n}->{_canonical_resn(n, aliases)}" for n in applied))
    unknown_t = _unknown_residue_summary(target_res, aliases)
    unknown_b = _unknown_residue_summary(binder_res, aliases)
    if unknown_t or unknown_b:
        merged = defaultdict(int)
        for src in (unknown_t, unknown_b):
            for name, count in src.items():
                merged[name] += count
        listed = ", ".join(f"{n}x{c}" for n, c in sorted(merged.items()))
        print(f"WARNING: unrecognised residue name(s): {listed}")
        print("         These score with the unknown-residue default and have no")
        print("         reference SASA. Map them with resn_aliases=SRC:DST if they")
        print("         are standard residues under another name.")
    if compute_sasa and truncate_anchor_sasa:
        print(f"Anchor SASA mode:                     CB-truncated side chain "
              f"({truncated_sites} site(s) re-measured)")
    elif compute_sasa:
        print("Anchor SASA mode:                     as-deposited side chain")
    direction_mode = "on (CA->CB vs CB->reactive)" if use_direction_score else "off (reported only)"
    print(f"Anchor direction term:                {direction_mode}")
    print(f"PyMOL labels:                        top {max(0, label_top)} unique binder sites")

    print("\n=== Target interface nucleophiles ===")
    if not result["target_records"]:
        print("No requested target nucleophile has its reactive atom inside the selected interface cutoff.")
    else:
        for rec in result["target_records"]:
            status = "available" if rec["available_for_pairing"] else "disulfide-bound/excluded"
            title = _residue_title(rec["target_resn"], rec["target_resi"])
            print(
                f"  {rec['target_chain']}/{title:<10} "
                f"reactive_atom={rec['nearest_reactive_atom']:<4} "
                f"reactive_min={rec['reactive_to_binder_min_A']:.2f} Å  {status}"
            )

    print("\n=== Binder mutation-to-Cys rank (best target per binder site) ===")
    if not result["binder_records"]:
        print("No binder anchor lies inside the candidate cutoff for an available target nucleophile.")
    else:
        _print_binder_rank_table(
            result["binder_records"], top_n=top_n, label_top=label_top
        )

    paths = _write_outputs(out_prefix, result)
    print("\nCSV outputs:")
    for path in paths:
        print("  " + os.path.abspath(path))

    if visualize:
        _visualize(object_name, result, label_top=label_top)

    print("\nInterpretation:")
    print("  A/B = prioritize for manual structural review; C = backup; D = generally avoid first round.")
    print("  The rank is target-pair specific and does not predict chemical reaction rate.")
    print("  Inspect top Cys rotamers, steric clashes and the complete linker/warhead before experiments.")
    print("  Compare WT, unmodified Cys mutant and conjugated binder binding experimentally.")
    print("  Only the top labelled binder sites are shown in PyMOL; all ranked sites are printed/CSV-exported.")
    return result


def warhead_help():
    """Print and return the documentation embedded in this script."""
    print(HELP_TEXT)
    return HELP_TEXT


def warhead_scan_help():
    """Alias for users who expect a scan-specific help command."""
    return warhead_help()


def _prompt_float(label, default):
    value = input(f"{label} [{default}]: ").strip()
    return float(value) if value else float(default)


def _prompt_int(label, default):
    value = input(f"{label} [{default}]: ").strip()
    return int(value) if value else int(default)


def warhead_scan_prompt():
    print("\nCovalent-warhead interface scanner v2.6")
    objects = cmd.get_object_list("(polymer.protein)")
    if not objects:
        raise ValueError("No protein molecular object is loaded")
    print("Protein objects: " + ", ".join(objects))
    default_object = objects[0] if len(objects) == 1 else ""
    object_name = input(f"Object name{f' [{default_object}]' if default_object else ''}: ").strip() or default_object
    object_name = _resolve_object(object_name)

    base_sel = f"({_model_clause(object_name)}) and polymer.protein"
    chains = list(cmd.get_chains(base_sel))
    print("Available protein chains: " + ", ".join(chains))
    print("Chain IDs are case-sensitive; combine several with + (e.g. H+L).")
    target_chain = input("TARGET chain ID(s): ").strip()
    binder_chain = input("BINDER chain ID(s): ").strip()

    interface_cutoff = _prompt_float("Target reactive-atom interface cutoff (Å)", 5.0)
    candidate_cutoff = _prompt_float("Binder anchor maximum search radius (Å)", 12.0)
    ideal_min = _prompt_float("Preferred anchor distance minimum (Å)", 5.0)
    ideal_max = _prompt_float("Preferred anchor distance maximum (Å)", 10.0)
    contact_cutoff = _prompt_float("Side-chain contact cutoff (Å)", 4.0)
    polar_cutoff = _prompt_float("Polar contact cutoff (Å)", 3.5)
    target_types = input("Target residue types [CYS+LYS+TYR+HIS]: ").strip() or "CYS+LYS+TYR+HIS"
    compute_sasa = _prompt_int("Compute approximate SASA (1=yes, 0=no)", 1)
    # Only worth asking when SASA is being computed at all.
    truncate_anchor_sasa = 1
    if compute_sasa:
        truncate_anchor_sasa = _prompt_int(
            "Measure anchor SASA with the side chain cut back to CB (1=yes, 0=no)", 1
        )
    use_direction_score = _prompt_int(
        "Include anchor direction in the score (1=yes, 0=report only)", 0
    )
    resn_aliases = input(
        "Extra residue aliases, e.g. MSE:MET+HIE:HIS [none]: "
    ).strip()
    out_prefix = input("Output prefix [warhead_scan]: ").strip() or "warhead_scan"

    return warhead_scan(
        target_chain=target_chain,
        binder_chain=binder_chain,
        interface_cutoff=interface_cutoff,
        candidate_cutoff=candidate_cutoff,
        ideal_min=ideal_min,
        ideal_max=ideal_max,
        target_types=target_types,
        object_name=object_name,
        out_prefix=out_prefix,
        compute_sasa=compute_sasa,
        contact_cutoff=contact_cutoff,
        polar_cutoff=polar_cutoff,
        truncate_anchor_sasa=truncate_anchor_sasa,
        use_direction_score=use_direction_score,
        resn_aliases=resn_aliases,
    )


cmd.extend("warhead_scan", warhead_scan)
cmd.extend("warhead_scan_prompt", warhead_scan_prompt)
cmd.extend("warhead_help", warhead_help)
cmd.extend("warhead_scan_help", warhead_scan_help)

print("Loaded covalent-warhead scanner v2.6.")
print("Help:        warhead_help  (alias: warhead_scan_help)")
print("Interactive: warhead_scan_prompt")
print("Direct example:")
print("  warhead_scan target_chain=B, binder_chain=A, interface_cutoff=5.0, candidate_cutoff=12.0, ideal_min=5.0, ideal_max=10.0")
