/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Uniformization.Assemble

/-!
# 一様化(取りまとめ)

★★★★**2026-09-03(第 1456)——292 KB / 325 宣言 / 5,745 行を 12 枚に割った**。

☆割る軸は**ファイル内の見出し**である。論文のセクションでは割れない
——このファイルは [GenEll] §3 の `Lemma 3.5` と `Proposition 3.4` の
2 項目しか持たず、セクションで割っても 146 KB のままだからである(実測)。

★後方参照は 9 件しか無く、すべて末尾の `.src` 台帳への言及であった
——実コードの依存は見出し順に流れているので、順に鎖でつなげる。

☆本ファイルは**既存の import を壊さないための取りまとめ**である。
必要な部分だけを引きたいときは `ABC3.Found.GenEll.Uniformization.<名前>` を直接 import する。

| 枚 | 中身 |
|---|---|
| `Basic` | 一様化の座標・周期性と対称性・`ω`-正規化・同種写像の定義的性質 |
| `VeluAnalytic` | 楕円関数の Liouville・Vélu の公式の解析側・極の打ち消し・正規化定数 |
| `Surjective` | 平行移動が誘導する置換・`℘` は全射・2-捩れ点 |
| `AdditionEntry` | Laurent の入口・加法定理の欠損関数・`z ≡ −w` の側・零点勘定の第一の煉瓦 |
| `AdditionODE` | ODE の一意性で零点勘定を回避する・因数分解・極を埋める道具 |
| `FilledPole` | 極を埋めた `F_w` |
| `AdditionFormula` | `y` 座標の加法公式・mathlib の `Point` への橋・倍加公式 |
| `Phi` | 一様化写像 `Φ : ℂ → E(ℂ)` |
| `GroupIso` | 部分群の原像・群同型 `ℂ/Λ ≅ E(ℂ)`・階数 1 の部分群・指数 |
| `Sublattice` | `Λ′` の基底と行列式・巡回の場合の代表系・代表系の対称性・Laurent |
| `G2G3` | `g₂`・`g₃` の変換・代数側との突き合わせ |
| `Assemble` | 代表系の像は `⟨Q⟩`・Vélu の `B` の和・組み立て・代表系の具体形・出典 |
-/
