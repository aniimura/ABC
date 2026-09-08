import ABC3.Found.PGC.TotallyRamifiedCriterion
import ABC3.Found.PGC.UnramifiedExtension

/-!
# [pGC] ★★★★`ramificationIdx = c` をノルムだけで出した

前波で「残る 1 ノードは `ramificationIdx = c`」と測った。本ファイルはそれを埋める。

## ★★まず自分の断定を検査した（持ち場の指摘）

前波の「`IntegerNorm.*` の道具一式は `adjoinIntegers K x` にそのまま移せる」は
★**実際に移して確かめた**。§1 の 4 つの測定はすべて通った:

```lean
-- (a) Ideal の型が defeq
example (K) (x) (I : Ideal ↥(adjoinIntegers K x)) :
    Ideal ↥(IntegerNorm.integerSubring (IntermediateField.adjoin K.carrier {x})) := I
-- (b) 𝒪[K.carrier] の所属がノルムで書ける
example (K) (y : K.carrier) : y ∈ 𝒪[K.carrier] ↔ ‖y‖ ≤ 1 := Iff.rfl
-- (c) 𝒪[K.carrier] の型が defeq
example (K) (y : 𝒪[K.carrier]) : ↥(IntegerNorm.integerSubring K.carrier) := y
-- (d) Subring として等しい
example (K) : (𝒪[K.carrier] : Subring K.carrier) = IntegerNorm.integerSubring K.carrier :=
  Subring.ext (fun _ => Iff.rfl)
```

★★★**(b)(c)(d) は `𝒪[·] = Valued.integer` の層も透明だということである**。
#69 が名指ししていた `Valued` 経由の遅さは、少なくとも
**所属・型・部分環の同一性**には現れない（#341 の続き）。

## 出したもの

| 宣言 | 内容 |
|---|---|
| ★★`ramificationIdx_eq_of_span` | ★抽象核（純可換環論）単項イデアルの分岐指数 |
| `mem_span_pow_iff_norm_le` | `∈ (π^c) ↔ ‖y‖ ≤ ‖π‖^c` |
| ★★★`ramificationIdx_eq_valIndex` | `ramificationIdx (π_K) (π) = c` |
| ★★★★`ramificationIdx_maximalIdeal_eq` | ★**極大イデアルの形**（`UnramifiedExtension.lean:425` と同じ形） |

★抽象核 `ramificationIdx_eq_of_span` には分岐・付値・ノルムの語が 1 語も出ない。
中身は `Ideal.ramificationIdx_spec`
（`NumberTheory/RamificationInertia/Ramification.lean:80`、★同名の
`Ideal.ramificationIdx'` が `RingTheory/…` にあるが**別物**）に
`Ideal.map_span` を差し込むだけである。

定義の字面は `ramificationIdx p P = sSup {n | map f p ≤ P ^ n}` で、
★持ち場の docstring 「何乗まで入るか」は★**正しい**（実際に読んで確かめた）。

## ★★どこで止まったか（`file:line`）

残るのは **2 つ**であり、どちらも★**配管ではなく内容**である。

1. ★`hcompat`（整数環の間の `algebraMap` が体の間のものの制限であること）を
   `AdjoinIntegers.lean:89 adjoinIntegersAlgebraMap` から取ること。
   ★本波はここに降りていない（本定理では仮説で受けている）。
2. ★★**`IsTotallyRamifiedAdjoin K x` から `c = n` を出すこと**。
   `TotallyRamified.lean:54` で `IsTotallyRamifiedAdjoin K x := inertiaDegree K x = 1`、
   `UnramifiedExtension.lean:444 ramificationIndex_mul_inertiaDegree` で
   `e·f = [K(x):K]` なので `e = n` は出る。
   本ファイルの `ramificationIdx_maximalIdeal_eq` は逆に `c` から `e` を計算するので、
   ★**この 2 つを繋げば `IsTotallyRamifiedAdjoin ⇒ c = n ⇒ hvalK`**になる。
   ★ただし `π` を実際に取る（`adjoinIntegers K x` が DVR であること）が先に要る。

## ★前回の①③の現在

* ① **不分岐側（`f = p`）** —— ★**依然に生きている。本波も触っていない。**
* ③ **構成側の `adjoin`** —— ★★**さらに剥がれた**。
  整数環の側（`𝒪[·]` と `adjoinIntegers`）はもう完全に透明である。
  残るのは上の 2 つだけ。

## 逸脱の記録

1. `hcompat` を仮説で受ける（上の 1）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

namespace RamIndexNorm

open ABC3.Skeleton.PGC
open scoped NormedField Valued

section Probe

variable {p : ℕ} [Fact p.Prime]

/-- 測定 (a): `Ideal ↥(adjoinIntegers K x)` と `Ideal ↥(integerSubring _)` は defeq か。 -/
example (K : PAdicLocalField p) (x : K.closure)
    (I : Ideal ↥(adjoinIntegers K x)) :
    Ideal ↥(IntegerNorm.integerSubring
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) := I

/-- 測定 (b): `𝒪[K.carrier]` の所属はノルムで書けるか。 -/
example (K : PAdicLocalField p) (y : K.carrier) :
    y ∈ 𝒪[K.carrier] ↔ ‖y‖ ≤ 1 := Iff.rfl

end Probe

/-! ## §2 ★抽象核（純可換環論）—— `ramificationIdx` の仕様を単項イデアルで -/

section Core

/-- ★★**抽象核** —— `span {a}` と `span {b}` の分岐指数は
`algebraMap a ∈ (b^c)` かつ `∉ (b^{c+1})` で決まる。

★分岐・付値・ノルムの語が 1 語も出ない。
`Ideal.ramificationIdx_spec`（`NumberTheory/RamificationInertia/Ramification.lean:80`）
に `Ideal.map_span` を差し込むだけ。 -/
theorem ramificationIdx_eq_of_span {R S : Type*} [CommRing R] [CommRing S] [Algebra R S]
    {a : R} {b : S} {c : ℕ}
    (hle : algebraMap R S a ∈ Ideal.span ({b ^ c} : Set S))
    (hgt : algebraMap R S a ∉ Ideal.span ({b ^ (c + 1)} : Set S)) :
    Ideal.ramificationIdx (Ideal.span ({a} : Set R)) (Ideal.span ({b} : Set S)) = c := by
  refine Ideal.ramificationIdx_spec ?_ ?_
  · rw [Ideal.map_span, Set.image_singleton, Ideal.span_singleton_pow, Ideal.span_le,
      Set.singleton_subset_iff]
    exact hle
  · intro hcon
    rw [Ideal.map_span, Set.image_singleton, Ideal.span_singleton_pow, Ideal.span_le,
      Set.singleton_subset_iff] at hcon
    exact hgt hcon

end Core

/-! ## §3 ノルム層 —— `∈ (π^c)` をノルムで言う -/

section Norm

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★**`(π^c)` への所属はノルムで書ける**。 -/
theorem mem_span_pow_iff_norm_le {π : M} (hπ0 : 0 < ‖π‖)
    (hπmem : π ∈ IntegerNorm.integerSubring M)
    (y : ↥(IntegerNorm.integerSubring M)) (c : ℕ) :
    y ∈ Ideal.span ({(⟨π, hπmem⟩ : ↥(IntegerNorm.integerSubring M)) ^ c} : Set _)
      ↔ ‖(y : M)‖ ≤ ‖π‖ ^ c := by
  have hp0 : (0 : ℝ) < ‖π ^ c‖ := by rw [norm_pow]; exact pow_pos hπ0 _
  have hpn : ‖π ^ c‖ = ‖π‖ ^ c := norm_pow π c
  rw [Ideal.mem_span_singleton]
  constructor
  · rintro ⟨w, hw⟩
    have hcoe : (y : M) = π ^ c * (w : M) := by
      simpa using congrArg (fun z : ↥(IntegerNorm.integerSubring M) => (z : M)) hw
    rw [← hpn]
    exact (IntegerNorm.dvd_iff_norm_le (π := π ^ c) (z := (y : M)) hp0).mp
      ⟨(w : M), w.2, hcoe⟩
  · intro hy
    rw [← hpn] at hy
    obtain ⟨w, hw, hyw⟩ := (IntegerNorm.dvd_iff_norm_le (π := π ^ c) (z := (y : M)) hp0).mpr hy
    refine ⟨⟨w, hw⟩, ?_⟩
    apply Subtype.ext
    simpa using hyw

end Norm

/-! ## §4 ★★★★`ramificationIdx = c` をノルムだけで -/

section Main

variable {K M : Type*} [NormedField K] [IsUltrametricDist K] [NormedField M]
  [IsUltrametricDist M] [Algebra K M]
  [Algebra ↥(IntegerNorm.integerSubring K) ↥(IntegerNorm.integerSubring M)]

/-- ★★★**分岐指数は値群の指数 `c` に等しい**（単項イデアルの形）。

★`hcompat` は「整数環の間の `algebraMap` は体の間のものの制限」。
`adjoinIntegers` の場合は `AdjoinIntegers.lean` の
`adjoinIntegersAlgebraMap` がこの形で与えるはず（★本波未測定）。 -/
theorem ramificationIdx_eq_valIndex
    (hcompat : ∀ y : ↥(IntegerNorm.integerSubring K),
      ((algebraMap ↥(IntegerNorm.integerSubring K) ↥(IntegerNorm.integerSubring M) y :
        ↥(IntegerNorm.integerSubring M)) : M) = algebraMap K M (y : K))
    {πK : K} {π : M} {c : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hπmem : π ∈ IntegerNorm.integerSubring M)
    (hKmem : πK ∈ IntegerNorm.integerSubring K)
    (hval : ‖algebraMap K M πK‖ = ‖π‖ ^ c) :
    Ideal.ramificationIdx
        (Ideal.span ({⟨πK, hKmem⟩} : Set ↥(IntegerNorm.integerSubring K)))
        (Ideal.span ({⟨π, hπmem⟩} : Set ↥(IntegerNorm.integerSubring M))) = c := by
  refine ramificationIdx_eq_of_span ?_ ?_
  · rw [mem_span_pow_iff_norm_le hπ0 hπmem, hcompat]
    exact le_of_eq hval
  · intro hcon
    rw [mem_span_pow_iff_norm_le hπ0 hπmem, hcompat] at hcon
    rw [hval] at hcon
    have hlt : ‖π‖ ^ (c + 1) < ‖π‖ ^ c :=
      pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 (by omega)
    linarith

/-- ★★★★**極大イデアルの形** —— `UnramifiedExtension.lean:425 ramificationIndex`
と同じ形をノルムだけから出す。 -/
theorem ramificationIdx_maximalIdeal_eq
    (hcompat : ∀ y : ↥(IntegerNorm.integerSubring K),
      ((algebraMap ↥(IntegerNorm.integerSubring K) ↥(IntegerNorm.integerSubring M) y :
        ↥(IntegerNorm.integerSubring M)) : M) = algebraMap K M (y : K))
    {πK : K} {π : M} {c : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hπmem : π ∈ IntegerNorm.integerSubring M)
    (hvalM : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m)
    (hK0 : 0 < ‖πK‖) (hK1 : ‖πK‖ < 1) (hKmem : πK ∈ IntegerNorm.integerSubring K)
    (hvalKK : ∀ y : K, y ≠ 0 → ∃ m : ℤ, ‖y‖ = ‖πK‖ ^ m)
    (hval : ‖algebraMap K M πK‖ = ‖π‖ ^ c) :
    letI := IntegerNorm.isLocalRing_integerSubring (M := K)
    letI := IntegerNorm.isLocalRing_integerSubring (M := M)
    Ideal.ramificationIdx (IsLocalRing.maximalIdeal ↥(IntegerNorm.integerSubring K))
      (IsLocalRing.maximalIdeal ↥(IntegerNorm.integerSubring M)) = c := by
  letI := IntegerNorm.isLocalRing_integerSubring (M := K)
  letI := IntegerNorm.isLocalRing_integerSubring (M := M)
  rw [IntegerNorm.maximalIdeal_eq_span hK0 hK1 hvalKK hKmem,
    IntegerNorm.maximalIdeal_eq_span hπ0 hπ1 hvalM hπmem]
  exact ramificationIdx_eq_valIndex hcompat hπ0 hπ1 hπmem hKmem hval

end Main

/-! ## §5 ★`PAdicLocalField` への代入の測定 -/

section Instantiate

variable {p : ℕ} [Fact p.Prime]

/-- ★★測定 (c): `𝒪[K.carrier]` の**型**は `↥(integerSubring K.carrier)` と defeq か。 -/
example (K : PAdicLocalField p) (y : 𝒪[K.carrier]) :
    ↥(IntegerNorm.integerSubring K.carrier) := y

/-- ★★測定 (d): 環の構造も一致するか（`Subring` として同じか）。 -/
example (K : PAdicLocalField p) :
    (𝒪[K.carrier] : Subring K.carrier) = IntegerNorm.integerSubring K.carrier :=
  Subring.ext (fun _ => Iff.rfl)

end Instantiate




/-! ## `.src` と 公理 -/

def ramificationIdx_maximalIdeal_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms ramificationIdx_eq_of_span
#print axioms mem_span_pow_iff_norm_le
#print axioms ramificationIdx_eq_valIndex
#print axioms ramificationIdx_maximalIdeal_eq

end RamIndexNorm

end ABC3.Found.PGC
