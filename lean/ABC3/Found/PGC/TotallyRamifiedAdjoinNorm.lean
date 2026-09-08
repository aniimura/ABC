import ABC3.Found.PGC.RamificationIndexNorm
import ABC3.Found.PGC.TotallyRamified

/-!
# [pGC] ★★★★残っていた 2 つが両方落ちた —— `hcompat` と `e = n`

## ★先に安い順を測った（持ち場の問い）

| 項目 | 結果 |
|---|---|
| `hcompat` | ★★**`rfl`**（1 行） |
| `IsTotallyRamifiedAdjoin ⇒ e = n` | ★★**3 行**（先例 `TotallyRamified.lean:60` と同じ手） |

両方とも落ちた。★★**持ち場が写した先例がそのまま使えた**。

`hcompat` は `adjoinIntegersAlgebra`（`AdjoinIntegers.lean:319`）が
`adjoinIntegersAlgebraMap` の `toAlgebra` であり、その写像が
`y ↦ ⟨⟨algebraMap K.carrier K.closure y, _⟩, _⟩` だからである
（★本体は本波で初めて読んだ。docstring だけでは分からなかった）。

## 出したもの

| 宣言 | 内容 |
|---|---|
| ★★★`coe_algebraMap_adjoinIntegers` | **`hcompat`**（`rfl`） |
| ★★★`ramificationIndex_eq_finrank_of_isTotallyRamified` | **`e = [K(x):K]`** |
| ★★`exists_zpow_norm_of_dvr` | ★抽象核。DVR → 値群の離散性（ノルム） |
| ★★★★`valK_of_ramificationIdx_eq` | ★**完全分岐（体）⇒ `hvalK`（ノルム）** |

★`exists_zpow_norm_of_dvr` には分岐の語が 1 語も出ない。
★**分数体の仮説は要らなかった** —— 非アルキメデス体では
`‖z‖ ≤ 1` か `‖z⁻¹‖ ≤ 1` のどちらかが必ず成り立つので場分けで済む。

## ★★どこで止まったか（`file:line`）

残るのは ★**§3 を `PAdicLocalField` に代入すること 1 ノードだけ**である。
要る部品はすべて場所が分かっている:

| 部品 | 場所 |
|---|---|
| `hcompat` | ★本ファイル §1 |
| `he : e = n` | ★本ファイル §1（`ramificationIndex` の定義を展開すれば同じ形） |
| `π`（`K(x)` の素元） | `UnramifiedExtension.lean:702` が `isDiscreteValuationRing_adjoinIntegers K x` と `IsDiscreteValuationRing.exists_irreducible` で取っている。★`letI := nontriviallyNormedField_adjoin K x` を先に置くこと（同 `:690-695` の注意書き） |
| `πK`（`K` の素元） | `UnramifiedExtension.lean:236` の `isDiscreteValuationRing_carrierIntegers K` から同様 |
| `hvalM` / `hvalKK` | ★本ファイル §2 `exists_zpow_norm_of_dvr` |
| `hiso : ‖algebraMap K.carrier M a‖ = ‖a‖` | `spectralNorm_extends`（`AdjoinIntegers.lean:102` が同じ手で使っている） |
| `hval : ‖algebraMap πK‖ = ‖π‖^c` | ★`c` を `TotallyRamified.exists_valSub_gen` から取る |

★★**本波はこの代入に降りていない。** 部品は全部揃っているが、
`letI` の順と instance の探索（同 `:690-695` が「`whnf` でタイムアウトする」と警告している）が
未測定である。

## ★前回の①③の現在

* ① **不分岐側（`f = p`）** —— ★**依然に生きている。本波も触っていない。**
  （本ファイルの `IsTotallyRamifiedAdjoin` は `f = 1` の場合である。）
* ③ **構成側の `adjoin`** —— ★★`hcompat` が `rfl` だったことでさらに剥がれた。
  残るのは上の代入 1 ノードだけである。

## 逸脱の記録

1. §3 は `hiso`（包含写像が等長）を仮説で受ける。`spectralNorm_extends` が与えるはずだが
   ★本波未測定である。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

namespace TotRamAdjoin

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-! ## §1 ★★残っていた 2 つのうちの安い方 -/

section Two

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**`hcompat` は `rfl`**。

`RamificationIndexNorm.ramificationIdx_eq_valIndex` が仮説で受けていた
「整数環の間の `algebraMap` は体の間のものの制限」は、
`adjoinIntegersAlgebra`（`AdjoinIntegers.lean:319`）に対して★**定義的に成り立つ**。 -/
theorem coe_algebraMap_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    (y : 𝒪[K.carrier]) :
    ((algebraMap 𝒪[K.carrier] (adjoinIntegers K x) y : adjoinIntegers K x) :
       IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = algebraMap K.carrier
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (y : K.carrier) := rfl

/-- ★★★**`IsTotallyRamifiedAdjoin ⇒ e = [K(x):K]`**。

★先例 `TotallyRamified.lean:60 finrank_eq_one_of_isUnramified_of_isTotallyRamified` と同じ手で
`ramificationIndex_mul_inertiaDegree`（`UnramifiedExtension.lean:444`）に `f = 1` を入れるだけ。 -/
theorem ramificationIndex_eq_finrank_of_isTotallyRamified (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) :
    ramificationIndex K x
      = Module.finrank K.carrier
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  have h := ramificationIndex_mul_inertiaDegree K x
  rw [show inertiaDegree K x = 1 from ht] at h
  omega

end Two

/-! ## §2 ★★DVR から値群の離散性をノルムで出す -/

section Discrete

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★**抽象核** —— `𝒪_M` が DVR で `ϖ` が既約なら、
`M` の値群は `‖ϖ‖^ℤ`。

★分岐の語が 1 語も出ない。`z` が整数のときは
`IsDiscreteValuationRing.associated_pow_irreducible`（`RingTheory/DiscreteValuationRing/Basic.lean:321`）
と `IntegerNorm.isUnit_iff_norm_eq_one`。一般の `z` は `‖z‖ ≤ 1` か `‖z⁻¹‖ ≤ 1` で場分けする
（★分数体の仮説は要らない）。 -/
theorem exists_zpow_norm_of_dvr {ϖ : M} (_hϖ0 : 0 < ‖ϖ‖)
    (hϖmem : ϖ ∈ IntegerNorm.integerSubring M)
    [IsDiscreteValuationRing ↥(IntegerNorm.integerSubring M)]
    (hirr : Irreducible (⟨ϖ, hϖmem⟩ : ↥(IntegerNorm.integerSubring M)))
    {z : M} (hz : z ≠ 0) : ∃ m : ℤ, ‖z‖ = ‖ϖ‖ ^ m := by
  have hint : ∀ y : ↥(IntegerNorm.integerSubring M), y ≠ 0 → ∃ n : ℕ, ‖(y : M)‖ = ‖ϖ‖ ^ n := by
    intro y hy
    obtain ⟨n, u, hu⟩ := IsDiscreteValuationRing.associated_pow_irreducible hy hirr
    refine ⟨n, ?_⟩
    have hcoe := congrArg (fun w : ↥(IntegerNorm.integerSubring M) => ‖(w : M)‖) hu
    simp only [Subring.coe_mul, SubmonoidClass.coe_pow, norm_mul, norm_pow] at hcoe
    have hu1 : ‖(((u : ↥(IntegerNorm.integerSubring M))) : M)‖ = 1 :=
      IntegerNorm.isUnit_iff_norm_eq_one.mp u.isUnit
    rw [hu1, mul_one] at hcoe
    simpa using hcoe
  rcases le_or_gt ‖z‖ 1 with h1 | h1
  · have hzmem : z ∈ IntegerNorm.integerSubring M := h1
    have hy0 : (⟨z, hzmem⟩ : ↥(IntegerNorm.integerSubring M)) ≠ 0 := fun hc => hz (by
      simpa using congrArg (fun w : ↥(IntegerNorm.integerSubring M) => (w : M)) hc)
    obtain ⟨n, hn⟩ := hint _ hy0
    exact ⟨(n : ℤ), by simpa [zpow_natCast] using hn⟩
  · have hinv : ‖z⁻¹‖ ≤ 1 := by
      rw [norm_inv]
      exact inv_le_one_of_one_le₀ (le_of_lt h1)
    have hzi : z⁻¹ ≠ 0 := inv_ne_zero hz
    have hy0 : (⟨z⁻¹, hinv⟩ : ↥(IntegerNorm.integerSubring M)) ≠ 0 := fun hc => hzi (by
      simpa using congrArg (fun w : ↥(IntegerNorm.integerSubring M) => (w : M)) hc)
    obtain ⟨n, hn⟩ := hint _ hy0
    refine ⟨-(n : ℤ), ?_⟩
    have hn' : ‖z⁻¹‖ = ‖ϖ‖ ^ n := by simpa using hn
    rw [norm_inv] at hn'
    rw [zpow_neg, zpow_natCast, ← hn', inv_inv]

end Discrete

/-! ## §3 ★★★★`e = n` から `hvalK` を出す -/

section Assemble

variable {K M : Type*} [NormedField K] [IsUltrametricDist K] [NormedField M]
  [IsUltrametricDist M] [Algebra K M]
  [Algebra ↥(IntegerNorm.integerSubring K) ↥(IntegerNorm.integerSubring M)]

/-- ★★★★**分岐指数が拡大次数に等しければ `hvalK` が出る**。

★これが「完全分岐（体の言葉）⇒ `hvalK`（ノルムの言葉）」の本体である。
`RamificationIndexNorm.ramificationIdx_maximalIdeal_eq` で `e = c` を計算し、
仮説 `he : e = n` と合わせて `c = n`。あとは `hvalKK` を `m` 乗するだけ。 -/
theorem valK_of_ramificationIdx_eq
    (hcompat : ∀ y : ↥(IntegerNorm.integerSubring K),
      ((algebraMap ↥(IntegerNorm.integerSubring K) ↥(IntegerNorm.integerSubring M) y :
        ↥(IntegerNorm.integerSubring M)) : M) = algebraMap K M (y : K))
    {πK : K} {π : M} {c n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hπmem : π ∈ IntegerNorm.integerSubring M)
    (hvalM : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m)
    (hK0 : 0 < ‖πK‖) (hK1 : ‖πK‖ < 1) (hKmem : πK ∈ IntegerNorm.integerSubring K)
    (hvalKK : ∀ y : K, y ≠ 0 → ∃ m : ℤ, ‖y‖ = ‖πK‖ ^ m)
    (hval : ‖algebraMap K M πK‖ = ‖π‖ ^ c)
    (hiso : ∀ a : K, ‖algebraMap K M a‖ = ‖a‖)
    (he : letI := IntegerNorm.isLocalRing_integerSubring (M := K)
          letI := IntegerNorm.isLocalRing_integerSubring (M := M)
          Ideal.ramificationIdx (IsLocalRing.maximalIdeal ↥(IntegerNorm.integerSubring K))
            (IsLocalRing.maximalIdeal ↥(IntegerNorm.integerSubring M)) = n) :
    ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m) := by
  letI := IntegerNorm.isLocalRing_integerSubring (M := K)
  letI := IntegerNorm.isLocalRing_integerSubring (M := M)
  have hcn : c = n := by
    rw [← RamIndexNorm.ramificationIdx_maximalIdeal_eq hcompat hπ0 hπ1 hπmem hvalM
      hK0 hK1 hKmem hvalKK hval]
    exact he
  intro a ha
  obtain ⟨m, hm⟩ := hvalKK a ha
  refine ⟨m, ?_⟩
  rw [hiso a, hm, ← hiso πK, hval, ← zpow_natCast ‖π‖ c, ← zpow_mul, hcn]

end Assemble



/-! ## `.src` と 公理 -/

def valK_of_ramificationIdx_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms coe_algebraMap_adjoinIntegers
#print axioms ramificationIndex_eq_finrank_of_isTotallyRamified
#print axioms exists_zpow_norm_of_dvr
#print axioms valK_of_ramificationIdx_eq

end TotRamAdjoin

end ABC3.Found.PGC
