---
name: pgc-ramification-naturality-gap
description: pGC Theorem 4.2のΦ構成——裸のΓ_K≅Γ_K′の部分はsorry無しで完成、FilteredGroup.Isoへの持ち上げに要る自然性公理がRamificationFiltrationのInterfaceに無い(2026-09-04)
metadata:
  type: project
---

★★★2026-09-04 続報(完成): 下記で「見込みを確認した」としていた「延長+共役+
選択非依存」の部分を、実際に`Found/PGC/GaloisTransfer.lean`にsorry無しで
完成させた——`galMulEquivOf`(互換な延長ᾱから誘導される裸のΓ_K≅Γ_K′)、
`galMulEquivOf_indep`(2つの異なる延長から誘導される同型が、K′.absGalの
内部自己同型でちょうど繋がる——`FilteredGroup.Iso.setoid`の同値関係と
同じ形)。局所類体論は一切使っていない。★★下記の「別の穴」(FilteredGroup.Iso
への持ち上げに要る自然性公理がRamificationFiltrationに無い)は未解決のまま
残る——本ファイルはその手前(裸の群同型の部分)を完成させたに留まる。

Theorem 4.2(`Skeleton/PGC/Section4.lean`)の `implicitStep`(体の同型 α:K≃K′ を
代数閉包への同型 ᾱ:K̄≃K̄′ へ延長し、共役で誘導される Γ_K≅Γ_K′ が ᾱ の取り方に
よらないことを示す構成、原文の「自然な射」Φ)を実際に埋められるか検討した。

**良いニュース**: 「延長+共役+選択非依存」の部分自体は**局所類体論を経由しない**
——一般の Galois 理論だけで組める見込みを確認した。鍵になる mathlib の道具
(実測、在庫あり):
- `IsAlgClosure.equivOfEquiv (L) (M) (hSR : S≃+*R) : L≃+*M` —— 代数閉包の
  同型を延長する一般定理。`K.closure = AlgebraicClosure K.carrier` は
  `IsAlgClosure K.carrier K.closure` のインスタンスを自動で持つ(`infer_instance`
  で実測)。
- `IsAlgClosure.equivOfEquiv_algebraMap`/`_symm_algebraMap` —— 延長が基礎体上で
  もとの同型と一致することの証明。これで conjugation `g ↦ ᾱ∘g∘ᾱ⁻¹` が
  `K′.absGal` の元(K′ を固定する)であることが出せる。
- 独立性(2つの延長 ᾱ₁,ᾱ₂ の差 `c := ᾱ₂⁻¹∘ᾱ₁ ∈ Γ_K` を取ると、誘導される
  写像は `φ₁ = (conj by φ₂(c)) ∘ φ₂` という関係になり、これはちょうど
  `FilteredGroup.Iso.setoid` の同値関係(後合成による内部自己同型)に一致する)
  も群論だけで出る計算——実際に手で確認した。

**★★見つかった別の穴**: この構成が実際に返すのは Γ_K の**裸の**群同型
(`ContinuousMulEquiv`)であって、`FilteredGroup.Iso`(= 高次分岐群のフィルトレーション
`Gv` を保つ同型)ではない。`Interface.PGC.RamificationFiltration`
(`Interface/PGC/LocalFieldData.lean`)の現在の定義は

```
structure RamificationFiltration (p : ℕ) [Fact p.Prime] where
  Gv : (K : PAdicLocalField p) → ℝ → Subgroup K.absGal
  isClosed : ∀ K v, IsClosed (Gv K v : Set K.absGal)
  isNormal : ∀ K v, (Gv K v).Normal
  antitone : ∀ K, Antitone (Gv K)
```

——**K ごとに完全に独立な `Gv` を許す**抽象データであり、「K の同型から誘導される
共役が `Gv` を保つ」という**自然性の公理が無い**。本物の高次分岐群(付値から
定義される)ならこの自然性は当然成り立つはずだが、現在の Interface はそれを
要求していないので、上の Φ 構成が「裸の同型」から「filtered group の同型」へ
橋渡しできない——`RF` にこの自然性を追加で課すか(逸脱として記録可能、
結論を密輸することにはならない——本物の性質を明示化するだけ)、`RF` を具体的に
構成してから示すかのどちらかが必要。

**教訓**: Cor 3.1/3.3 の「裸の同型 vs filtered group の同型」の不整合
([[pgc-filtered-vs-bare-galois-iso]])と**同根の問題**が Theorem 4.2 側にも
潜んでいた——「フィルトレーションを保つ」という条件を Interface レベルで
自然性公理として持たない限り、どのセクションでも同じ穴に落ちる。次に
`RamificationFiltration` を触るときは、この自然性公理を追加することを
検討する価値がある。

関連: [[pgc-filtered-vs-bare-galois-iso]], [[measure-mathlib-before-skeleton]]
