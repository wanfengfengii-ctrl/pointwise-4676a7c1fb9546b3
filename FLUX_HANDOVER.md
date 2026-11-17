# 逐面通量交付交接

本文件说明新增实现的面身份、方向、支撑来源与结算不变量，便于后续维护和核账。

## 模块分工

| 文件 | 职责 |
| --- | --- |
| `remap/velocities.py` | 仿射速度场校验、不可变版本存储与读取（`velocity_revisions`） |
| `remap/fluxgeom.py` | 有向面积向量、十个有向二阶矩、平面多边形裁剪、共面交叠、凸差分解 |
| `remap/flux.py` | 交叠体边界提取、共面合并、公共面细分与配对、通量积分、发布前自检 |
| `remap/flux_settlement.py` | `execute_flux` 流程：重放、绑定校验、结算、发布 |
| `remap/archive.py` | 公开门面新增 `put_velocity/execute_flux/read_flux`，事务发布 |
| `verify_flux.py` | 不导入 `remap` 几何的独立解析核账器 |
| `verify_all.py` | 原回归 + 独立核账的单入口（运行清单 `verify`） |

新增存储：`velocity_revisions`、`flux_results`、`flux_receipts`，与原 `results/receipts` 完全分离（独立请求命名空间）。

## 面身份

一块面片由「交叠体 + 支撑平面上的区域」唯一确定：

- 交叠体 = (源单元 ID, 目标 ID)。不同源即使坐标相同也是不同交叠体。
- 面片的 `origins.source_faces` 是包含整片的源四面体面（三个原顶点坐标字典序排序），`origins.target_faces` 是包含整片的盒面 `[target_id, axis, side]`。
- 两类来源可同时非空（源面恰好落在盒面上）。共面平面在 `_merged_planes` 中先合并，因此双重来源只产生一条记录。
- 辅助细分线（邻接切分）不作为来源；`origins` 只列真实支撑平面。

## 方向

- 每个交叠体用其外向半空间法向定向：源面法向朝四面体外部，盒面法向 `±e_axis` 背离盒内。
- `vertices` 按该法向右手绕行；`moments["1"] = ∫ n dS` 是面积向量（有向），无根式。
- 公共面两侧：a 侧按 a 的外向法向，b 侧为同区域反向绕行，故 `flux`、`moments` 精确相反，`mate` 双向。

## 支撑来源与配对

- 只把「另一目标的盒面与本目标盒面在同一点集上贴合（法向相反）」视为邻接；共边、共点不配对。
- 公共区域用 2D 凸裁剪求正面积交叠后细分；余片用凸差分解（逐个扣除洞）得到，每个洞减去后保持线性片数。
- 配对判定以“两块目标盒面共享平面且区域正面积重叠”为准，独立源不相互配对。
- 同一交叠体内面片只共边/共点，不重叠；核账器用独立 2D 裁剪逐对验证。

## 结算

`flux = Σ_axis Σ_monomial c[axis][monomial] · moments[monomial][axis]`，其中 `rho·v_axis` 的十个单项式系数由 `_quadratic_terms` 展开（仿射×仿射=二次）。核账器用单纯形参数化的 `∫_Δ u^i v^j du dv = i!j!/(i+j+2)!` 独立重算。

发布前自检（`settle`）：

- 配对对称、同源、通量互为相反数；
- 每源全部面通量之和 = `external_net`（内部交换相消）；
- 每个交叠体：面积向量和为零、`∫x_i n_i dS = 体积`、`∫x_i^2 n_i dS = 2*一阶矩`、六个交叉矩对应分量等于相应一阶矩、`Σ flux = ∫ div(rho·v) dV`。

以上任一失败抛 `BalanceError`，成果与回执都不写入。

## 幂等与不漂移

- 计划、源、目标、场均为不可变版本，执行时按 ID 重新读取，调用方后续修改输入或替换同名 head/场不影响历史。
- 同请求同绑定直接返回原成果；异绑定抛 `RequestConflict`，绝不重算或另发成果。
- 成果与回执同事务发布；发布失败整次回滚，同号可重试；旧质量成果与旧通量成果不受影响。
- 返回值为独立解码副本，读取方修改不影响持久状态。

## 运行

```sh
python3.11 -B verify_all.py          # 原 37 项回归 + 独立核账（含随机场景）
python3.11 -B verify_flux.py --check # 仅独立核账
python3.11 -B demo.py                # 部分邻接/双来源/二次积分/重开
python3.11 -B flux_examples.py flux  # 公开有限样例
```
