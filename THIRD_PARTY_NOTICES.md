# 第三方来源与许可

本仓库的数据由下列上游资源计算得出。本仓库**不转发**其中任何文件（人体模型和 acupuncture-3d 原始数据都不在仓库里）；但计算结果属于它们的衍生作品，所以下列署名必须随本仓库数据一起保留。

## 1. BodyVista（人体模型、结构目录）

- 地址：https://github.com/FreedomIntelligence/BodyVista ，使用的版本是 commit `82b2ea12f71fd123eb728f465eb413bb8d6a150e`
- 用到的文件：`assets/full_body.glb`、`catalog.json`（由 `pipeline/fetch_upstream.py` 下载到本地，不随仓库分发）
- 许可：BodyVista 自有代码是 MIT；**模型不在 MIT 范围内**，按 `catalog.json` 中每个结构的 `source_attribution` 记录处理。2899 个结构标为 `Z-Anatomy CC-BY-SA-4.0 with per-source attribution and NC exceptions; research only`，1 个（皮肤 `bp2_FJ2810`）标为 `BodyParts3D CC BY 4.0`。BodyVista 说明这些标签尚未逐项核实。
- 本仓库的计算只用到 127 个结构（对流水线实际加载的网格逐一记录得出）：皮肤 1 个（BodyParts3D）、骨骼 124 个（颅面骨、锁骨，右侧肱骨/桡骨/尺骨/股骨，双侧胫腓骨、腕骨、掌骨、指骨、跗骨、跖骨、趾骨等），以及左右眼球前段 2 个（定位面部穴位用）。
  其中不含上游列为禁止商用的内耳（Dundee，CC BY-NC-SA 4.0）和肾脏（Lissie Cowley，CC BY-NC 4.0）结构。

## 2. Z-Anatomy（骨骼网格的原始来源）

- 地址：https://github.com/Z-Anatomy/Models-of-human-anatomy
- 许可：CC BY-SA 4.0（上游另有内耳、肾脏两个 NC 例外，本仓库未使用）
- 署名：

> Z-Anatomy, © Gauthier Kervyn, licensed under CC Attribution-Share Alike 4.0 International.

## 3. BodyParts3D（皮肤外壳的原始来源）

- 地址：https://dbarchive.biosciencedbc.jp/en/bodyparts3d/
- 许可：现行 CC BY 4.0（旧版本为 CC BY-SA 2.1 JP）
- 署名：

> BodyParts3D, © The Database Center for Life Science licensed under CC Attribution 4.0 International.

## 4. acupuncture-3d（穴位比例锚点、人体测量数据）

- 原地址：https://github.com/FuriaRozkwit/acupuncture-3d ，作者 Paweł Turzyński，使用的版本是 commit `7c4fd3feb1b754ef2a90b67aaabf3498114dc932`。**该仓库在 2026-11-19 已无法访问。**
- 用到的文件和许可（按该仓库 `LICENSE` 和 `ATTRIBUTION.md`）：
  - `data/points.json`（穴位锚点）、`data/meridians.json`（穴位顺序）：MIT，Copyright (c) 2026 Pawel Turzynski
  - `data/rig_fitted.json`、`data/structures.json`（从 Z-Anatomy / BodyParts3D 网格上量出的数据）：CC BY-SA 4.0
- 本仓库没有使用该项目的任何代码，也没有使用其中受版权保护的书籍摘录（`point_meta.json`）。解析器在 `pipeline/body.py` 中独立重写。
- 本仓库对锚点做的修改：坐标系变换；左右分别计算；骨骼与皮肤错位校正；按 GB/T 12346-2021 重新定位 18 个穴位（ST3、ST5、ST6、ST7、ST8、GB2、GB3、GB4、GB5、GB6、GB14、GB21、GB40、LI16、KI3、KI5、KI6、HT1）；新增印堂 GV24+；把任督脉穴位投到正中线；穴位顺序补入 BL39 和 GV24+。逐穴说明见 `data/acupoints.json` 的 `accuracyNote`。

## 5. GB/T 12346-2021《经穴名称与定位》

- 本仓库只使用穴位名称和编码。**不收录**标准中的定位条文全文；`accuracyNote` 中引用了少量条文短句，用于说明为什么重新定位。

## 许可汇总

| 内容 | 许可 |
|---|---|
| `data/`、`docs/check/`、`pipeline/` 下的 `.json` / `.log` | CC BY-SA 4.0（继承第 1–4 项） |
| `pipeline/*.py` | MIT（本仓库作者自写） |

本文件只记录已知来源，不代表任何上游作者授予额外权限。如需商用，请自行核实各上游来源的许可。
