import ABC3.Found.PGC.TotallyRamifiedAdjoinSupply

/-!
# [pGC] ★★★★★★★★diamond は**橋を掛けるのではなく避けられた**

前波の結論は「次の 1 点は 2 つのスペクトルノルムの一致を橋渡す補題」だった。
★★**これは必要なかった。** 完備な底を `ℚ_[p]` ではなく
★**`K.carrier` 自身**にすれば diamond に触れない。

## ★先に測ったこと（§1、持ち場の提案 3 番目が当たった）

| 測定 | 結果 |
|---|---|
| (1) `CompleteSpace K.carrier` | ★出る |
| (2) `NormedAlgebra K.carrier ↥K.carrier⟮x⟯` | ★★**出る** |
| (3) `Algebra.IsAlgebraic K.carrier ↥K.carrier⟮x⟯` | ★★**出る** |
| (4) `‖(p : K(x))‖ = p⁻¹` | ★`spectralNorm_extends` 2 回で出る（`NormedAlgebra ℚ_[p]` 不要） |

★★前波で見つからなかったのは `NormedAlgebra ℚ_[p] ↥K.carrier⟮x⟯` と
`Algebra.IsAlgebraic ℚ_[p] ↥K.carrier⟮x⟯` であり、
★**底を `K.carrier` にすればどちらも出る**。
`AdjoinPAdicLocalField.lean:44-49` の警告（「両者の一致を橋渡す補題を先に用意する」）は
★**真だが、必須ではなかった**。★スペクトルノルムの一致自体は本波も**測っていない**。

## ★実測したエラーと直し方（linter が diamond を名指しした）

```
⚠ `[Algebra K M]` and `[NormedAlgebra K M]` can be used to infer conflicting versions
of `[SMul K M]`.
💡 Of these, `[Algebra K M]` may be removed.
```

★両方書くと `FiniteDimensional K M` の合成が失敗する。`[Algebra K M]` を消す。

## 出したもの

| 宣言 | 内容 |
|---|---|
| `norm_natCast_p_adjoin` | `‖(p : K(x))‖ = p⁻¹`（`NormedAlgebra ℚ_[p]` なし） |
| ★★★`exists_norm_sub_algebraMap_le_prod_axDecay_of_completeBase` | ★**底が完備なら出口が出る** |
| ★★★★`exists_norm_sub_algebraMap_le_prod_axDecay_of_adjoin` | ★★**`PAdicLocalField` の上での出口** |

## ★★前波の自分の保留を検算した

前波は「これで残るのは `hnK` だけは★まだ真ではありません」と書いた。
★本波で `exists_norm_sub_algebraMap_le_prod_axDecay_of_adjoin` が通ったので、
★★**ノルムの仮説は全部消えた**ことが確かめられた。残るのは

* `ht : IsTotallyRamifiedAdjoin K x`（完全分岐）
* `hnK : [K(x):K] = p^{k+1}`（次数）
* 塔のデータ `g` / `hg` / `τ` / `hτ` / `s` / `hsg`

の 3 群だけである。★これらはすべて「一般の `K` と `x` から塔を作る」の**内容**である。

## ★①③の現在

* ① **不分岐側（`f = p`）** —— ★**依然に生きている。本波も触っていない。**
* ③ **構成側の `adjoin`** —— ★★**配管は全部通った**。
  残るのは上の 3 群（内容）だけである。

## 逸脱の記録

1. 完備な底を `ℚ_[p]` ではなく `K` 自身にした（上の理由）。
   `HarithPAdicSupply` の `ℚ_[p]` 版は残してある（別の宣言）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

namespace SpectralBridge

open ABC3.Skeleton.PGC
open scoped NormedField Valued

section Probe

variable {p : ℕ} [Fact p.Prime]

/-- 測定 (1): `CompleteSpace K.carrier` は出るか。 -/
noncomputable example (K : PAdicLocalField p) : CompleteSpace K.carrier := inferInstance

/-- 測定 (2): `NormedAlgebra K.carrier ↥K.carrier⟮x⟯` は出るか。 -/
noncomputable example (K : PAdicLocalField p) (x : K.closure) :
    NormedAlgebra K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := inferInstance

/-- 測定 (3): `Algebra.IsAlgebraic K.carrier ↥K.carrier⟮x⟯` は出るか。 -/
example (K : PAdicLocalField p) (x : K.closure) :
    Algebra.IsAlgebraic K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := inferInstance

/-- ★★**`‖(p : K(x))‖ = p⁻¹`** —— `NormedAlgebra ℚ_[p] K(x)` なしで出る。

★`spectralNorm_extends` を 2 回使うだけ（`K(x)/K.carrier` と `K.carrier/ℚ_[p]`）。 -/
theorem norm_natCast_p_adjoin (K : PAdicLocalField p) (x : K.closure) :
    ‖((p : ℕ) : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖
      = ((p : ℝ))⁻¹ := by
  have h1 : ((p : ℕ) : IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = algebraMap K.carrier
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ((p : ℕ) : K.carrier) := by
    rw [map_natCast]
  have h3 : ‖algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ((p : ℕ) : K.carrier)‖
        = ‖((p : ℕ) : K.carrier)‖ := spectralNorm_extends _
  have h2 : ((p : ℕ) : K.carrier) = algebraMap ℚ_[p] K.carrier ((p : ℕ) : ℚ_[p]) := by
    rw [map_natCast]
  rw [h1, h3, h2, norm_algebraMap, Padic.norm_p]

end Probe

/-! ## §2 ★★★diamond を**避ける** —— 完備な底を `ℚ_p` ではなく `K` にする -/

section Base

/-- ★★★★★★★**底が完備なノルム体なら出口が出る**。

★★`HarithPAdicSupply.exists_norm_sub_algebraMap_le_prod_axDecay_of_padic` は
`[NormedAlgebra ℚ_[p] M]` を要求し、それが `K.carrier⟮x⟯` で見つからなかった。
★**底を `K` 自身にすれば diamond に触れない**。

★`hiso` は `PureStepSetup.norm_algEquiv_eq`（`PureStepSetup.lean:283`）を
`k := K` で使うだけ。`hnormp` は仮説で受ける（具体層では §1 測定 (4) が与える）。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_completeBase
    {K M : Type*} [NontriviallyNormedField K] [IsUltrametricDist K] [CompleteSpace K]
    [NormedField M] [IsUltrametricDist M] [NormedAlgebra K M]
    [FiniteDimensional K M] [Algebra.IsAlgebraic K M]
    {p : ℕ} [Fact p.Prime] {k : ℕ} {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (τ : M →+ M) (hτ : ∀ z : M, τ z = g z)
    (s : ℕ → (M →+* M)) (hsg : ∀ j, ∀ z : M, s j z = (g ^ p ^ j) z)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (hnormp : ‖((p : ℕ) : M)‖ = ((p : ℝ))⁻¹) (x : M) :
    ∃ y : GainedTowerModel.GaloisTower.twr g p k,
      ‖x - algebraMap (GainedTowerModel.GaloisTower.twr g p k) M y‖
        ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp : p.Prime := Fact.out
  have hq0 : 0 < p ^ (k + 1) := pow_pos hp.pos _
  obtain ⟨e, he, heM⟩ :=
    HarithPAdic.exists_absRamIndex (q := p ^ (k + 1)) hq0 hp.one_lt hπ0 hπ1 hvalK hnormp
  exact HarithAssembly.exists_norm_sub_algebraMap_le_prod_axDecay_of_norm
    (p := p) (e := e) (k := k) (π := π) g hg
    (fun z => PureStepSetup.norm_algEquiv_eq g z) τ hτ s hsg he hnK hvalK hnormp heM x

end Base

/-! ## §3 ★★★★★★★★`PAdicLocalField` の上での出口 —— 残るのは `hnK` だけ -/

section Adjoin

variable {p : ℕ} [Fact p.Prime]

/-- ★★★★★★★★**完全分岐な単純拡大の上での出口**。

★★ノルムの仮説は全部消えた。残るのは
* `ht : IsTotallyRamifiedAdjoin K x`（完全分岐）
* `hnK : [K(x):K] = p^{k+1}`（次数）
* 塔のデータ `g` / `hg` / `τ` / `hτ` / `s` / `hsg`
だけである。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_adjoin
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {k : ℕ}
    (hnK : Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = p ^ (k + 1))
    (g : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (hg : orderOf g = p ^ (k + 1))
    (τ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) →+
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (hτ : ∀ z, τ z = g z)
    (s : ℕ → ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)) →+*
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
    (hsg : ∀ j, ∀ z, s j z = (g ^ p ^ j) z)
    (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ∃ w : GainedTowerModel.GaloisTower.twr g p k,
      ‖y - algebraMap (GainedTowerModel.GaloisTower.twr g p k)
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) w‖
        ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ y - y‖ := by
  obtain ⟨π, hπ0, hπ1, hvalK⟩ := TotRamSupply.valK_of_isTotallyRamifiedAdjoin K x ht
  rw [hnK] at hvalK
  exact exists_norm_sub_algebraMap_le_prod_axDecay_of_completeBase
    (π := π) g hg τ hτ s hsg hπ0 hπ1 hnK hvalK (norm_natCast_p_adjoin K x) y

end Adjoin



/-! ## `.src` と 公理 -/

def exists_norm_sub_algebraMap_le_prod_axDecay_of_adjoin.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms norm_natCast_p_adjoin
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_completeBase
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_adjoin

end SpectralBridge

end ABC3.Found.PGC
