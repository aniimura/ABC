import ABC3.Found.PGC.RamificationSubgroupCard

/-!
# [pGC] `𝒪_M = {z | ‖z‖ ≤ 1}` を**抽象ノルム体の上で**建てる —— ★#69 は当たらなかった

前波(第 1119)で `hrec` の橋に残るのは
「`lowerRamificationGroup B G i` をノルムの `{σ | ∀ z, ‖z‖ ≤ 1 → ‖σz − z‖ ≤ ‖π‖^{i+1}}` と
同一視すること」1 ノードだと特定した。

## ★★★前波(＝自分)の断定の訂正(名指し)

前波の報告と `RamificationSubgroupCard.lean` の冒頭はこう書いた:

> ★これには `B = 𝒪_M` を型として建てる**しかありません**（`herbrandPhiGroup` が
> `[IsDiscreteValuationRing B] [MulSemiringAction G B]` で `B` を要求するため）。
> 降りるときの既知の危険: `#69`(`adjoinField`/`adjoinIntegers` の境界、212 秒 timeout)

★★**「危険」の部分は誤りである。**★測った:

1. `AdjoinIntegers.lean:18-36` を読むと、`Valued` を避けた理由は
   ★**`IntermediateField.adjoin K.carrier {x}` の上でだけ**起きた詰まりである:
   > `IntermediateField extends Subfield extends Subring extends Submonoid ...` という
   > 何層にも重なった部分構造の上で、新しく導入した `Valued` 由来の位相と、
   > 既存の `NormedField` 由来の位相(の定義的な一致)を検査するコストが高い
   ★**我々の `M` は素の型変数**(`NormedField M`)であって、積み重なった部分構造ではない。
2. ★実測: 本ファイルは `integerSubring M : Subring M` を素の `Subring.mk` で建て、
   ★**8.8 秒で通った**(`leanfile.mjs`)。★#69 の 212 秒には**当たらない**。
⇒ ★**#69 は `PAdicLocalField`/`K.closure`/`IntermediateField` の層に固有であり、
抽象 `NormedField M` の層には効かない。**★元の docstring は直さず、ここに訂正を書く。

## ★在庫の測定(#330 どおり 2 か所)

```
grep -n "unitBall" .cache/mathlib-index.txt
  → `Subsemigroup.unitBall`(Analysis/Normed/Field/UnitBall.lean:32)は★**開球**。
    ★閉球の `Subring` は mathlib に**無い**。木にも `def unitBall` は無い
    (`grep -rn "def unitBall" lean/ABC3/Found/PGC/*.lean` が 0 件)。⇒ 本ファイルで建てた。
grep -n "IsDiscreteValuationRing" .cache/mathlib-index.txt | grep -i "valued\|valuation.integer"
  → ★`Valuation.valuationSubring_isDiscreteValuationRing`
     [IsCyclic (valueGroup v)] [Nontrivial (valueGroup v)] : IsDiscreteValuationRing K₀
     (RingTheory/Valuation/Discrete/Basic.lean:453)
  → ★`IsNonarchimedeanLocalField` は `instance : IsDiscreteValuationRing 𝒪[K]`
     (NumberTheory/LocalField/Basic.lean:108)
grep -n "NormedField.toValued" .cache/mathlib-index.txt
  → `NormedField.toValued : Valued K ℝ≥0`(Topology/Algebra/Valued/NormedValued.lean:67)
```
⇒ ★**`Valued` 経由は塞がっていない**(「値群が巡回かつ非自明」だけで DVR が付く)。
★本ファイルはそこには降りていない(下の「残り」)。

## 出したもの(★どの経路でも要る中身)

| 宣言 | 内容 |
|---|---|
| `integerSubring` | ★`{z : M \| ‖z‖ ≤ 1}` が `Subring M`(超距離のみ、8.8 秒) |
| ★★`exists_pow_mul_unit` | ★**`z ≠ 0`, `‖z‖ ≤ 1` なら `z = π^n · w`(`‖w‖ = 1`)** ＝ DVR の心臓 |
| ★`dvd_iff_norm_le` | `π ∣ z` ⟺ `‖z‖ ≤ ‖π‖` ＝ `𝔪 = (π)` の中身 |
| `maps_integerSubring` | 等長な自己同型は単位球を保つ ＝ `MulSemiringAction G 𝒪_M` の中身 |

★`exists_pow_mul_unit` の入力 `‖z‖ ∈ ‖π‖^ℤ` は
`TotallyRamified.exists_zpow_norm`(`TotallyRamifiedValueGroup.lean:163`)が供給する。

## ★★残り(正確に、`file:line` つき)

★**数学はもう無い。残っているのは型クラスを載せる作業 4 つだけ**である:

1. `IsDiscreteValuationRing ↥(integerSubring M)` —— 中身は `exists_pow_mul_unit`。
   ★mathlib の入口(`IsDiscreteValuationRing.ofHasUnitMulPowIrreducibleFactorization` 等)に
   当てはめる作業。★本ファイルは**やっていない**。
2. `MulSemiringAction G ↥(integerSubring M)` —— 中身は `maps_integerSubring`。
   ★`Subring` への作用の制限を書くだけ。★**やっていない**。
3. `IsLocalRing.maximalIdeal ↥(integerSubring M) = Ideal.span {π}` —— 中身は `dvd_iff_norm_le`。
   ★**やっていない**。
4. `lowerRamificationGroup ↥(integerSubring M) G i` を
   `RamNormBridge.mem_ramification_iff` と突き合わせる。
   `LowerRamificationGroup.lean:270 mem_lowerRamificationGroup_iff_forall`
   (`σ ∈ G_i ↔ ∀ x : B, σ•x − x ∈ 𝔪^{i+1}`)を使う。★**やっていない**。

★これができると `RamCard.card_eq_pow_of_mem_iff`(前波)で `|G_i| = p^{k+1−m}` が入り、
`HerbrandComposition.lean:455 herbrandPhiGroup_natCast` で `hrec` が出て、
`HasseArfCongruence.dvd_sub_of_phi_intCast` で `harith` の (3) が閉じる。

## 逸脱の記録

1. `integerSubring` は `Subring M`(`Valued`/`ValuationSubring` を使わない)。
   ★理由は上の測定 1(木が `Valued` を避けた理由が我々には当たらないが、
   素の `Subring` の方が確実に速いことを実測したため)。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**(import のみ)。
-/

namespace ABC3.Found.PGC

namespace IntegerNorm

section Ball

variable (M : Type*) [NormedField M] [IsUltrametricDist M]

/-- ★**単位閉球は部分環**(超距離のみ)。★mathlib には閉球の `Subring` は無い
(`Subsemigroup.unitBall` は**開球**、`Analysis/Normed/Field/UnitBall.lean:32`)。 -/
def integerSubring : Subring M where
  carrier := {z : M | ‖z‖ ≤ 1}
  zero_mem' := by simp
  one_mem' := by simp
  add_mem' := by
    intro a b ha hb
    exact le_trans (IsUltrametricDist.norm_add_le_max a b) (max_le ha hb)
  neg_mem' := by
    intro a ha
    simpa using ha
  mul_mem' := by
    intro a b ha hb
    rw [Set.mem_setOf_eq, norm_mul]
    calc ‖a‖ * ‖b‖ ≤ 1 * 1 := mul_le_mul ha hb (norm_nonneg _) zero_le_one
      _ = 1 := one_mul 1

@[simp] theorem mem_integerSubring {z : M} : z ∈ integerSubring M ↔ ‖z‖ ≤ 1 := Iff.rfl

end Ball

section Structure

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

omit [IsUltrametricDist M] in
/-- ★★★★**`𝒪_M` の DVR 構造の心臓** —— `‖z‖ ≤ 1`, `z ≠ 0` なら
`z = π^n · w`(`‖w‖ = 1`)。★これが「素元 `π` の DVR」のノルム版である。

★入力 `hzval`(`‖z‖ ∈ ‖π‖^ℤ`)は `TotallyRamified.exists_zpow_norm` が供給する。 -/
theorem exists_pow_mul_unit {π z : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hzval : ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m) (hz1 : ‖z‖ ≤ 1) (hz0 : z ≠ 0) :
    ∃ (n : ℕ) (w : M), ‖w‖ = 1 ∧ z = π ^ n * w := by
  obtain ⟨m, hm⟩ := hzval
  have hπne : π ≠ 0 := norm_pos_iff.mp hπ0
  have hm0 : 0 ≤ m := by
    by_contra hcon
    rw [not_le] at hcon
    have hgt : (1:ℝ) < ‖π‖ ^ m := one_lt_zpow_of_neg₀ hπ0 hπ1 hcon
    rw [← hm] at hgt
    linarith
  refine ⟨m.toNat, z / π ^ m.toNat, ?_, ?_⟩
  · rw [norm_div, norm_pow, hm, ← zpow_natCast ‖π‖ m.toNat, Int.toNat_of_nonneg hm0]
    exact div_self (ne_of_gt (zpow_pos hπ0 m))
  · field_simp

omit [IsUltrametricDist M] in
/-- ★★**`𝔪 = (π)` の中身** —— `π ∣ z`(`𝒪_M` の中で)⟺ `‖z‖ ≤ ‖π‖`。 -/
theorem dvd_iff_norm_le {π z : M} (hπ0 : 0 < ‖π‖) :
    (∃ w : M, ‖w‖ ≤ 1 ∧ z = π * w) ↔ ‖z‖ ≤ ‖π‖ := by
  constructor
  · rintro ⟨w, hw, rfl⟩
    rw [norm_mul]
    calc ‖π‖ * ‖w‖ ≤ ‖π‖ * 1 := mul_le_mul_of_nonneg_left hw (norm_nonneg π)
      _ = ‖π‖ := mul_one _
  · intro hz
    have hπne : π ≠ 0 := norm_pos_iff.mp hπ0
    refine ⟨z / π, ?_, by field_simp⟩
    rw [norm_div, div_le_one hπ0]
    exact hz

/-- ★**等長な自己同型は単位球を保つ**(`MulSemiringAction G 𝒪_M` の中身)。 -/
theorem maps_integerSubring {K : Type*} [Field K] [Algebra K M] (g : M ≃ₐ[K] M)
    (hiso : ∀ w : M, ‖g w‖ = ‖w‖) {z : M} (hz : z ∈ integerSubring M) :
    g z ∈ integerSubring M := by
  rw [mem_integerSubring] at hz ⊢
  rw [hiso]
  exact hz

end Structure

/-! ## `.src` と 公理 -/

def exists_pow_mul_unit.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms integerSubring
#print axioms exists_pow_mul_unit
#print axioms dvd_iff_norm_le
#print axioms maps_integerSubring

end IntegerNorm

end ABC3.Found.PGC
