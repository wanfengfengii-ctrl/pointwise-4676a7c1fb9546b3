# Exact observation remapping

Python 3.11 标准库工程，将带仿射密度的三维四面体观测重采样到轴对齐目标盒，并在 SQLite 中保存可追溯成果。

## 输入与业务约定

源版本是一组非空、ID 唯一的 Cell。每个 Cell 有四个非共面的有理坐标顶点，以及全局坐标密度系数 (a,b,c,d)，表示 rho(x,y,z)=a+b*x+c*y+d*z。四个顶点的密度均须非负。不同源单元独立结算，可以表示叠置观测。ID 和名称为非空字符串；数据库版本、计划、成果 ID 为返回的整数。

Target(id, lower, upper) 定义闭盒，三个方向均须有正宽度。一个计划中的目标非空、ID 唯一，盒内部不得相互重叠，共边界合法。坐标和系数接受 int、Fraction 或 Fraction 可解析的有理数字符串，拒绝 bool 和 float。边、点、面接触不产生体积和质量。

源同名登记每次生成新版本并更新 head。save_plan 每次保存一个独立的不可变目标版本，并绑定指定源版本；同名计划不会替换历史计划。更换源 head、修改调用方输入和重开数据库均不改变历史绑定。

几何输出 Moments(volume, first)，first 为全局坐标积分 (积分 x dV, 积分 y dV, 积分 z dV)。质量按 a*volume+b*first[0]+c*first[1]+d*first[2] 精确结算。每个源的完整质量等于其所有目标配对质量与未覆盖质量之和，不对已覆盖部分重新归一化。顶点置换和目标细分不改变相应合计。

## Python API

从 remap 导入 Archive、Cell、Target 和错误类：

```python
from fractions import Fraction
from remap import Archive, Cell, Target

with Archive("observations.sqlite") as a:
    revision = a.put_source("run", [
        Cell("s", ((0,0,0),(1,0,0),(0,1,0),(0,0,1)), (0,1,0,0))
    ])
    plan = a.save_plan("west", revision, [
        Target("w", (0,0,0), (Fraction(1,2),1,1))
    ])
    result = a.execute_plan(plan, "request-001")
    assert result["transfers"][0]["mass"] == Fraction(11,384)
    assert result["source_ledger"][0]["uncovered_mass"] == Fraction(5,384)
    assert a.read_result(result["id"]) == result
```

- put_source(name,cells) -> 源版本 ID；source_head(name)、read_source(revision)、list_sources() 用于查看目录与版本。read_source 返回 id、name、按单元 ID 排序的 Cell 元组。
- save_plan(name,source_revision,targets) -> 计划 ID；read_plan(id) 返回源/目标版本、目标元组及名称；list_plans() 按计划 ID 排序。
- execute_plan(plan_id,request_id,observer=None) -> 完整成果字典。request_id 在库内唯一；同计划重放返回原成果，绑定其他计划时抛 RequestConflict。计算前检查重放，发布事务内再次检查。
- read_result(id) -> 独立的新字典；list_results() 返回按成果 ID 排序的 id、plan_id、request_id 列表。
- remap.geometry.intersection_moments(vertices,lower,upper) 和 tetrahedron_moments(vertices) 是独立精确几何 API，返回不可变 Moments。

成果含 id、plan_id、source_revision、target_revision、receipt，以及以下数组：
- transfers：source、target、volume、first、mass；按 source 再 target 的字符串序排列，所有配对及零值均保留。
- source_ledger：source、total_mass、covered_mass、uncovered_mass。
- targets：target、mass；包括零质量目标。

Python 返回的几何量和质量始终为 Fraction。目录数组按数据库整数 ID 排序；源账与目标汇总按各自公开字符串 ID 排序。读取值的更改不会修改持久状态。

## 错误与发布边界

非法源、重复单元、退化四面体或负顶点密度抛 ValidationError，整次登记不改变目录。非法目标、内部重叠或未知源抛 PlanError，目标版本和计划均不保存。未知读取或执行计划抛 NotFoundError。同请求绑定其他计划抛 RequestConflict。内部负质量或超覆盖抛 BalanceError。

execute_plan 的可选 observer 是同步进度接收函数，每个真实交叠计算后接收 event=overlap、source、target、moments；它不改变结果内容，重放不再次通知。接收函数须只观察事件，不在该连接上重入写入操作；其异常原样传播，成果与回执不发布。几何异常也不发布。所有计算完成后，完整成果与回执在同一 SQLite 事务发布，数据库写入失败整次回滚。源版本与 head、目标版本与计划也分别同事务保存。SQLite 异常原样传播。Archive 为本地连接上下文，每个工作者使用独立实例；本项目不提供并行执行调度。

## CLI 与样例

在项目目录运行。示例数据库必须是新库；samples/plan.json 中 source_revision=1 对应首次登记。

```sh
python3.11 -B -m remap --db observations.sqlite source samples/source.json
python3.11 -B -m remap --db observations.sqlite plan samples/plan.json
python3.11 -B -m remap --db observations.sqlite execute 1 run-001
python3.11 -B -m remap --db observations.sqlite result 1
python3.11 -B -m remap --db observations.sqlite catalog
```

CLI 输出单个 JSON；有理数表示为精确字符串。业务、输入文件或数据库错误输出 error 并退出 2；argparse 参数错误也退出 2。数据库不存在则创建，父目录必须已经存在。CLI 的 head NAME 查询当前源版本。源文件与计划文件格式见 samples/，计划文件须填真实登记返回的源版本。

## 回归与演示

```sh
python3.11 -B -m unittest discover -s tests -v
python3.11 -B demo.py
```

完整回归覆盖独立解析积分、全局矩、顶点置换、穿顶点与沿棱共面裁剪、低维接触、细分合计、仿射多源结算、源和计划拒绝原子性、固定版本、返回值隔离、重开、请求重放和 CLI 子进程。失败回归采用真实进度回调异常与 SQLite 触发器拒绝写入，核对数据库中成果及回执行数，再核对重试。测试仅在自身临时目录创建库；demo 从有限样例真实登记、计划、执行、替换源并重开核对成果，打印新旧结果。运行清单见 .pointwise-runtime.json。

代码职责：sources.py 管理源校验及版本，plans.py 管理目标校验及固定计划，geometry.py 提供几何量，settlement.py 积分与结算，archive.py 管理成果、回执和公开门面，storage.py 管理 SQLite 事务。

## 固定计划逐面通量（新增）

在原有质量接口之外，`remap.Archive` 增加：

- `put_velocity(name, source_revision, fields) -> 场版本 ID`。`fields` 的键必须恰为该源版本的全部单元 ID，每个值为三个仿射分量，分量系数依次表示 `a+b*x+c*y+d*z`。同名登记只追加不可变版本，保存归一化后的副本。
- `execute_flux(plan_id, velocity_revision, request_id) -> 成果字典`。计划与场必须绑定同一源版本；计算完成后在单个事务内发布成果与回执。
- `read_flux(result_id) -> 独立副本`。通量请求与原质量请求使用独立命名空间，同号重放按计划与场版本核对，绑定不同即抛 `RequestConflict`。

成果含 `id`、`plan_id`、`source_revision`、`target_revision`、`velocity_revision`、`receipt`，以及以下数组：

- `faces`：每项 `{id,source,target,vertices,origins,moments,flux,mate}`。`vertices` 为不重复首点、沿交叠体外法向右手绕行的环；`origins.source_faces` 是支撑该片的源面（三个原顶点按坐标字典序排序），`origins.target_faces` 是本目标盒支撑面 `[target_id,axis,side]`，两类可同时非空。`moments` 恰含 `1,x,y,z,xx,xy,xz,yy,yz,zz` 十个键，每个值为 `(∫m·n_x dS, ∫m·n_y dS, ∫m·n_z dS)`；`moments["1"]` 即单位外法向的面积分（有向面积向量），不含根式。`flux` 为 `∫ rho·v·n dS` 的精确值，仿射×仿射按二次精确积分。
- `exchanges`：每对内部公共片一项 `{source,a,b,a_face,b_face,flux}`，`a<b`，`flux` 取 a 侧外向通量；对应 b 侧为相反数。零通量的公共片同样保留。
- `targets`：每个目标一项 `{target,net_outward}`，为该目标全部面通量之和（含内部交换）。
- `source_ledger`：每个源一项 `{source,external_net}`，为该源全部未配对面通量之和。

几何与结算不变量：

1. **面证据重建**：每次执行都从不可变源/目标/场版本精确重建边界，不复用任何历史面或矩。源四面体与目标盒的交叠体若非正体积则不产生面记录；低维接触也不产生。
2. **共面合并**：源面与目标面处于同一半空间时合并为一块物理片，`origins` 同时记录两类来源，但只保留一条记录、只积分一次。
3. **公共面细分**：同一源内相邻目标（盒面贴合且法向相反）的公共区域按实际正面积交叠细分；公共区两侧得到方向相反、矩与通量精确相反的两条记录，`mate` 双向引用。未被任何邻接目标覆盖的余片 `mate=None`。
4. **独立源**：不同源单元即使顶点与场相同也各自生成面记录，不共享、不配对。
5. **结算自检**：发布前核对每片配对对称、同源面通量总和等于 `external_net`、每个交叠体边界闭合（面积和为零）且散度积分等于 `3*体积` 等恒等式；违反抛 `BalanceError`，不发布任何成果或回执。

`verify_flux.py` 是独立解析核账器，不导入 `remap` 的任何几何代码：它用精确半空间枚举重建交叠体、用单纯形参数化与不同推导重算十个矩和通量、核对逐面覆盖与配对、散度积分、目标/源净账、顶点置换、低维接触、无覆盖目标与重开。`python3.11 -B verify_flux.py --check` 运行固定场景与确定性随机场景。

演示：`python3.11 -B demo.py` 现在同时展示部分邻接、双来源面、二次积分与持久重开。
