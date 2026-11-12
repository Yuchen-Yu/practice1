# 变质作用脱碳、含碳流体迁移和碳固定模拟


**Metamorphic Carbon Reaction–Transport Simulator** 是一个基于 Python 和 Streamlit 开发的一维交互式科学模拟器，用于展示变质作用过程中**脱碳、含碳流体迁移和碳固定**之间的耦合关系。

模型整合了经典矿物反应、实验约束或实验启发的 P–T–XCO₂ 反应边界、简化反应动力学、一维平流–扩散输运以及显式的碳质量守恒。

---

## 1. 模型概述

变质作用能够使碳在岩石、流体和新形成的矿物之间发生重新分配。

本模拟器主要描述三个相互联系的过程：

1. **脱碳作用（Decarbonation）**
2. **含碳变质流体迁移（Carbon transport）**
3. **碳固定 / 碳酸盐化（Carbon fixation / Carbonation）**

整体概念过程为：

```text
含碳酸盐原岩
      ↓
脱碳作用
      ↓
CO₂ 进入变质流体
      ↓
平流 + 扩散迁移
      ↓
重新固定为碳酸盐和/或输出系统
```

一维模型空间被划分为三个区域：

```text
0 km                    9 km                    20 km                  30 km
│------------------------│------------------------│-----------------------│
       脱碳区                   输运区                    碳固定区
 Decarbonation zone        Transport zone          Carbon fixation zone
```

---

## 2. 基本反应

### 2.1 脱碳作用

模型采用经典的方解石–石英–硅灰石反应：

```text
CaCO₃ + SiO₂ ⇌ CaSiO₃ + CO₂
```

即：

```text
Calcite + Quartz ⇌ Wollastonite + CO₂
```

正向反应表示原岩碳酸盐矿物中的碳被释放进入活动性流体碳库。

在模型中：

```text
Original solid carbon → Mobile fluid carbon
```

---

### 2.2 碳迁移

活动性流体碳采用一维平流–扩散–反应方程描述：

```text
∂C/∂t + v ∂C/∂x = D ∂²C/∂x² + R_dec - R_fix
```

其中：

- `C`：流体中的活动性碳浓度
- `v`：流体平流速度
- `D`：有效扩散系数
- `R_dec`：脱碳作用产生的碳
- `R_fix`：碳固定作用消耗的碳

当前模型将流体中的活动性碳统一表示为：

```text
CO₂-equivalent mobile carbon
```

并未进一步区分：

- CO₂(aq)
- HCO₃⁻
- CO₃²⁻
- CH₄
- 其他含碳流体物种

---

### 2.3 碳固定

碳固定采用：

```text
Mg₂SiO₄ + CO₂ ⇌ MgCO₃ + MgSiO₃
```

即：

```text
Forsterite + CO₂ ⇌ Magnesite + Enstatite
```

正向反应表示流体中的 CO₂ 被重新固定进入新的碳酸盐矿物。

在模型中：

```text
Mobile fluid carbon → Newly fixed carbon
```

---

## 3. P–T–XCO₂ 热力学控制

模型不采用简单的：

```text
温度超过某个阈值 → 自动发生反应
```

而是首先根据当前：

```text
Pressure
Temperature
Fluid XCO₂
```

判断正向反应是否在热力学上有利。

当前数值流程为：

```text
P–T–XCO₂ 条件
       ↓
反应平衡边界
       ↓
正向反应在热力学上是否有利？
       ↓
是
       ↓
计算热力学驱动力
       ↓
Arrhenius 型动力学
       ↓
有效反应速率

否
       ↓
正向有效反应速率 = 0
```

P–T 图中的反应边界和数值模型中的反应判据使用的是**同一组函数**。

因此：

```text
图中显示的热力学区域
=
数值模型实际使用的反应区域
```

这可以避免出现：

```text
P–T 点位于热力学不利区
但数值模型仍然继续大量反应
```

这样的逻辑矛盾。

---

## 4. 脱碳反应边界

脱碳反应为：

```text
Calcite + Quartz ⇌ Wollastonite + CO₂
```

Harker & Tuttle (1956) 通过实验约束了这一经典反应的 PCO₂–T 平衡关系。

当前代码采用一个用于教学和可视化的**平滑近似边界**。

参考温度与压力点为：

```python
T_ref = [
    480,
    550,
    600,
    650,
    700,
    750,
    790,
    820
]

P_ref = [
    0.10,
    0.30,
    0.55,
    0.95,
    1.55,
    2.35,
    3.05,
    3.70
]
```

其中：

```text
T = temperature in °C
P = pressure in kbar
```

代码使用三次多项式对这些参考点进行平滑。

需要特别说明：

> 当前代码中的三次多项式并不是 Harker & Tuttle (1956) 原文发表的解析拟合公式，而是为了教学、可视化和数值边界判断，对经典实验反应边界所作的近似表达。

P–T 图中：

- 主要参考温度范围内使用实线
- 超出参考温度范围的部分使用虚线

虚线仅用于展示外推趋势，不代表实验测量结果。

---

## 5. 碳固定反应边界

碳固定采用以下反应的逆反应：

```text
Magnesite + Enstatite ⇌ Forsterite + CO₂
```

Koziol & Newton (1998) 对该反应进行了实验研究。

实验范围约为：

```text
Temperature: 700–1100 °C
Pressure: 6–25 kbar
```

论文给出的实验平衡边界拟合为：

```text
P_eq = -3.215 - 4.784 × 10⁻³ T + 2.542 × 10⁻⁵ T²
```

其中：

```text
P_eq = equilibrium pressure in kbar
T    = temperature in °C
```

该公式在模拟器中直接用于计算碳固定反应边界。

对于正向碳固定反应：

```text
Forsterite + CO₂ → Magnesite + Enstatite
```

在给定温度下，当：

```text
P_fix > P_eq
```

时，模型认为正向碳固定在热力学上有利。

---

## 6. Fluid XCO₂ 的影响

真实变质反应不仅受压力和温度控制，也受到流体组成控制。

例如：

```text
Calcite + Quartz → Wollastonite + CO₂
```

会产生 CO₂。

因此，当流体中 CO₂ 的化学势升高时，进一步向右发生脱碳会变得更加困难。

模拟器提供：

```text
Fluid XCO₂
```

参数，用于表示流体中 CO₂ 的简化摩尔分数。

模型采用如下形式的理想流体化学势修正：

```text
RT ln(XCO₂ / XCO₂,ref)
```

其中：

```text
R          = gas constant
T          = absolute temperature
XCO₂       = selected CO₂ mole fraction
XCO₂,ref   = reference CO₂ mole fraction
```

随后结合局部 Clapeyron slope，对 P–T 反应边界随 XCO₂ 的变化进行近似调整。

经过 XCO₂ 修正后的同一条边界同时用于：

1. P–T 图绘制
2. 数值模型反应判据

因此改变 `Fluid XCO₂` 后：

```text
P–T curve changes
        +
Reaction gate changes
```

两者同步更新。

这一处理属于教学性简化，并不等同于完整的高温高压 H₂O–CO₂ 非理想流体状态方程。

---

## 7. 热力学区域

Panel **a** 中使用三个半透明区域表示不同的热力学状态。

### 红色区域：脱碳有利区

表示：

```text
Decarbonation thermodynamically favourable
```

正向反应：

```text
Calcite + Quartz → Wollastonite + CO₂
```

只有在：

```text
P_decarb < P_eq,decarb
```

时才允许发生。

位于该区域之外时：

```text
R_dec = 0
```

也就是说：

```text
Effective decarbonation rate = 0
```

---

### 绿色区域：碳固定有利区

表示：

```text
Carbonation thermodynamically favourable
```

正向反应：

```text
Forsterite + CO₂ → Magnesite + Enstatite
```

只有在：

```text
P_fix > P_eq,fix
```

时才允许发生。

位于该区域之外时：

```text
R_fix = 0
```

也就是说：

```text
Effective fixation rate = 0
```

---

### 灰色区域

灰色区域表示在当前简化模型下，相应正向反应在热力学上不利。

需要注意：

红色边界和绿色边界分别对应**两个不同的矿物反应体系**。

因此灰色区域不能严格理解为：

```text
一个统一岩石体系中的绝对无反应区
```

它只是将两个独立正向反应的热力学有利范围放在同一个 P–T 图中的教学性表达。

---

## 8. 反应动力学

热力学条件决定的是：

```text
反应是否可以向当前方向进行
```

而不是：

```text
反应进行得有多快
```

当反应进入热力学有利区域后，模型采用 Arrhenius 型关系描述温度对反应速率的影响：

```text
k(T) = k_ref × exp[
    -Ea/R × (1/T - 1/T_ref)
]
```

其中：

```text
k_ref = reference reaction rate
Ea    = activation energy
R     = gas constant
T     = absolute temperature
T_ref = reference temperature
```

当前教学参数为：

```text
Decarbonation Ea = 100 kJ mol⁻¹
Carbonation Ea   = 80 kJ mol⁻¹
Reference T      = 800 °C
```

这些参数主要用于控制模拟中的反应时间尺度，并未针对某一具体天然变质岩体系进行实验标定。

---

## 9. 距离反应边界的影响

模型并不是一进入热力学有利区就立即达到最大反应速率。

以脱碳作用为例，模型使用如下形式的平滑驱动力：

```text
f_dec = 1 - exp[
    -(P_eq - P) / P_scale
]
```

因此：

```text
接近平衡边界
→ 热力学驱动力较弱

深入热力学有利区
→ 热力学驱动力增强

位于热力学不利区
→ 驱动力 = 0
```

最终有效反应速率为：

```text
k_effective = k_Arrhenius × f_thermo
```

这一处理用于简化表达 reaction affinity 的概念。

---

## 10. 碳迁移

流体中的活动性碳通过：

```text
Advection + Diffusion
```

进行一维空间迁移。

### 平流

平流表示含碳流体整体向前移动。

主要控制参数为：

```text
Fluid velocity
```

单位：

```text
km kyr⁻¹
```

---

### 扩散

扩散表示由于碳浓度梯度产生的空间扩散。

主要控制参数为：

```text
Effective diffusivity
```

单位：

```text
km² kyr⁻¹
```

---

## 11. 碳质量守恒

模型显式追踪四个碳库。

### Original solid carbon

表示尚未通过脱碳作用进入流体的原始碳酸盐碳：

```text
C_original
```

---

### Mobile fluid carbon

表示随变质流体迁移的活动性碳：

```text
C_fluid
```

---

### Newly fixed carbon

表示通过碳固定作用进入新形成碳酸盐矿物的碳：

```text
C_fixed
```

---

### Exported carbon

表示通过模型右侧边界离开系统的碳：

```text
C_exported
```

模型满足：

```text
C_initial
=
C_original
+
C_fluid
+
C_fixed
+
C_exported
```

程序在每次模拟后自动检查质量守恒误差。

理想情况下：

```text
Initial carbon ≈ Final accounted carbon
```

数值差值应接近浮点计算误差范围。

---

## 12. 可视化结构

主复合图由三个 panel 组成。

---

### Panel a — P–T reaction boundaries

显示：

- Decarbonation equilibrium boundary
- Carbonation equilibrium boundary
- Decarbonation favourable field
- Carbonation favourable field
- Thermodynamically unfavourable schematic field
- Selected decarbonation-zone P–T state
- Selected fixation-zone P–T state

主要实验或参考温度范围内：

```text
solid line = experimentally constrained / reference interval
```

超出主要约束范围：

```text
dashed line = extrapolated for visualization
```

---

### Panel b — Spatial carbon distribution

显示当前时刻三种碳库的一维空间分布：

```text
Original solid carbon
Mobile fluid carbon
Newly fixed carbon
```

横坐标为：

```text
Distance along metamorphic system (km)
```

三个空间区域分别为：

```text
0–9 km
Decarbonation zone

9–20 km
Transport zone

20–30 km
Carbon fixation zone
```

该 panel 将动态示意图与底层数值浓度直接对应。

---

### Panel c — Dynamic reaction–transport schematic

动态展示：

```text
Decarbonation
      ↓
CO₂ release
      ↓
CO₂-bearing metamorphic fluid
      ↓
Advection + diffusion
      ↓
Carbon fixation and/or export
```

显示的矿物和组分包括：

- Calcite
- Quartz
- Wollastonite
- CO₂
- Forsterite
- Magnesite
- Enstatite

同时显示：

- Fluid-flow arrows
- Reaction progress
- Fixation progress
- Simulation time
- Carbon export

矿物转化程度和 CO₂ 粒子分布均与当前数值模型状态关联，而不是独立播放的预设动画。

---

### Panel d — Carbon reservoirs through time

显示以下碳库随时间的变化：

```text
Original solid carbon
Mobile fluid carbon
Newly fixed carbon
Exported carbon
```

横坐标：

```text
Time (kyr)
```

---

### Panel e — Mass-balance check

比较：

```text
Initial carbon
```

与：

```text
Accounted carbon
```

用于检查整个数值模拟过程中的碳质量守恒。

---

## 13. 交互参数

| 参数 | 含义 |
|---|---|
| Metamorphic temperature | 模型变质温度 |
| Decarbonation-zone pressure | 脱碳区压力 |
| Fixation-zone pressure | 碳固定区压力 |
| Fluid XCO₂ | 简化的流体 CO₂ 摩尔分数 |
| Fluid velocity | 一维流体平流速度 |
| Effective diffusivity | 有效碳扩散系数 |
| Reference decarbonation rate | 脱碳参考动力学速率 |
| Reference fixation rate | 碳固定参考动力学速率 |
| Maximum fixation capacity | 碳固定区最大储碳容量 |
| Simulation duration | 模拟持续时间 |

---

## 14. 运行逻辑

整个模拟可以概括为：

```text
Input:
Temperature
Pressure
Fluid XCO₂
Transport parameters
Kinetic parameters

        ↓

Calculate P–T–XCO₂ reaction boundaries

        ↓

Check thermodynamic favourability

        ↓

If favourable:
    calculate Arrhenius kinetic rate
    calculate distance-from-equilibrium driving force
    calculate effective reaction rate

If unfavourable:
    effective reaction rate = 0

        ↓

Decarbonation releases carbon

        ↓

Advection + diffusion transport mobile carbon

        ↓

Carbonation fixes part of mobile carbon

        ↓

Remaining carbon may leave the model

        ↓

Check carbon mass balance
```

---

## 15. 安装

建议使用：

```text
Python 3.9+
```

克隆或下载仓库后，运行：

```bash
pip install -r requirements.txt
```

也可以直接安装依赖：

```bash
pip install streamlit numpy pandas matplotlib
```

---

## 16. 运行模拟器

假设主程序文件为：

```text
metamorphic_carbon_simulator.py
```

在仓库目录运行：

```bash
streamlit run metamorphic_carbon_simulator.py
```

Streamlit 将启动本地 Web 应用，并通常会自动在默认浏览器中打开。

如果浏览器没有自动打开，可根据终端中显示的 Local URL 手动访问。

---

## 17. 仓库结构

建议 GitHub 仓库采用以下结构：

```text
metamorphic-carbon-simulator/
├── README.md
├── README_CN.md
├── metamorphic_carbon_simulator.py
├── requirements.txt
└── LICENSE
```

其中：

```text
README.md
```

为英文说明；

```text
README_CN.md
```

为中文说明；

```text
metamorphic_carbon_simulator.py
```

为主模拟程序；

```text
requirements.txt
```

记录 Python 依赖；

```text
LICENSE
```

记录代码许可协议。

---

## 18. 模型定位与限制

本模拟器是一个简化的一维变质碳反应–输运模型。

模型在热力学、动力学和空间结构上均进行了有意简化，因此模拟结果应理解为：

```text
conceptual to semi-quantitative simulation
```

而不是对某一具体天然变质系统的精确重建。

当前主要简化包括：

1. 采用一维空间结构；
2. 脱碳区和碳固定区分别使用单一代表性压力；
3. 整个模型采用统一变质温度；
4. 活动性流体碳统一表示为 CO₂-equivalent carbon；
5. 未显式区分 CO₂(aq)、HCO₃⁻、CO₃²⁻、CH₄ 等不同含碳流体组分；
6. 未显式计算氧逸度；
7. 未考虑孔隙度随反应变化；
8. 未考虑渗透率随反应变化；
9. 未考虑矿物反应表面积随时间变化；
10. 未执行完整的多组分 Gibbs free-energy minimization；
11. XCO₂ 对反应边界的影响采用理想化近似；
12. 脱碳 P–T 边界采用经典实验结果的平滑近似，而不是原始实验论文中直接给出的解析多项式；
13. 超出主要实验或参考范围的 P–T 曲线使用虚线，仅用于展示外推趋势；
14. Arrhenius 活化能和参考反应速率主要用于教学模拟，并未针对某一特定天然岩石体系进行标定；
15. 当前模型不考虑岩石变形、裂隙形成以及反应诱导渗透率变化；
16. 当前模型不显式求解多矿物、多流体组分之间的同步化学平衡。

P–T 区域主要用于判断：

```text
当前简化模型中的正向反应在热力学上是否有利
```

而实际反应程度还同时受到：

```text
Temperature-dependent kinetics
Fluid transport
Available reactants
Carbon fixation capacity
Simulation duration
```

等因素控制。

---

## 19. 参考文献

### 脱碳反应与 P–T–fluid equilibrium

Harker, R. I. & Tuttle, O. F. (1956). Experimental data on the PCO₂–T curve for the reaction calcite + quartz = wollastonite + carbon dioxide. *American Journal of Science*, **254**, 239–256.

Greenwood, H. J. (1967). Wollastonite: Stability in H₂O–CO₂ mixtures and occurrence in a contact-metamorphic aureole near Salmo, British Columbia, Canada. *American Mineralogist*, **52**, 1669–1680.

Tanner, S. B., Kerrick, D. M. & Lasaga, A. C. (1985). Experimental kinetic study of the reaction: calcite + quartz = wollastonite + carbon dioxide, from 1 to 3 kilobars and 500° to 850°C. *American Journal of Science*, **285**, 577–620.

DOI:

```text
https://doi.org/10.2475/ajs.285.7.577
```

---

### 碳固定与 Mg-silicate–carbonate equilibrium

Koziol, A. M. & Newton, R. C. (1998). Experimental determination of the reaction: magnesite + enstatite = forsterite + CO₂ in the ranges 6–25 kbar and 700–1100°C. *American Mineralogist*, **83**, 213–219.

DOI:

```text
https://doi.org/10.2138/am-1998-3-403
```

---

### H₂O–CO₂ 流体组成

Aranovich, L. Y. & Newton, R. C. (1999). Experimental determination of CO₂–H₂O activity-composition relations at 600–1000°C and 6–14 kbar by reversed decarbonation and dehydration reactions. *American Mineralogist*, **84**, 1319–1332.

---

### 变质反应动力学

Walther, J. V. & Wood, B. J. (1984). Rate and mechanism in prograde metamorphism. *Contributions to Mineralogy and Petrology*, **88**, 246–259.

Rubie, D. C. & Thompson, A. B. (1985). Kinetics of metamorphic reactions at elevated temperatures and pressures: an appraisal of available experimental data. In A. B. Thompson & D. C. Rubie (eds.), *Metamorphic Reactions: Kinetics, Textures, and Deformation*. Springer.

---

## 20. AI 使用披露

本项目在代码开发和文档整理过程中使用了生成式人工智能工具作为辅助。

OpenAI ChatGPT 被用于协助：

- Python 代码结构设计；
- Streamlit 交互界面设计；
- 部分数值模型代码实现；
- P–T 图和动态过程可视化代码编写；
- 代码调试；
- 模型逻辑检查；
- README 文档结构整理；
- 中英文文本润色；
- Markdown 格式整理。

项目中的以下内容由项目作者进行人工审查和最终确认：

- 科学问题；
- 模型总体结构；
- 脱碳和碳固定反应体系的选择；
- 文献选择；
- 热力学边界的解释；
- 模型参数设置；
- 数值模型逻辑；
- 可视化结构；
- 最终源代码。

生成式 AI 输出不被视为独立科学依据。

涉及以下内容时，应以原始文献为最终依据：

```text
Experimental reaction boundaries
Thermodynamic relationships
Reaction kinetics
Published numerical parameters
Mineral stability fields
```

代码中的教学性近似、参数化方案以及实验范围之外的外推部分已尽可能在 README 和源代码注释中明确标注。

---

## 21. GitHub Issues

如果发现：

- 代码错误；
- 数值稳定性问题；
- 热力学边界问题；
- 可视化问题；
- 参数设置问题；
- 文献引用问题；

或希望提出新的模型功能，可以通过 GitHub Issues 提交。

建议 Issue 中同时提供：

```text
Python version
Operating system
Streamlit version
Parameter settings
Error message
Screenshot
```

以便复现问题。

---

## 22. License

本项目采用 **MIT License**。

详细内容见：

```text
LICENSE
```

MIT License 允许在保留原始版权和许可声明的前提下：

- 使用；
- 复制；
- 修改；
- 合并；
- 发布；
- 分发；
- 再许可。

具体条款以仓库中的 `LICENSE` 文件为准。