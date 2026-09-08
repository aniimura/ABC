import ABC3.Found.PGC.IntegerRingInstances

/-!
# [pGC] ★★`IsDiscreteValuationRing 𝒪_M` —— 前波が「残り 1 本」と言ったもの

前波(第 1122)で残りは `IsDiscreteValuationRing ↥(integerSubring M)` 1 本と特定した。
★**本ファイルで載った。**★前波の「中身は `exists_pow_mul_unit` で尽きている」は**当たり**であった
(自分の断定を 3 波連続で覆してきたので今回も疑ったが、今回は正しかった)。

## ★入口の測定(3 候補のうちどれを使ったか)

```
sed -n '145,150p' lean/.lake/packages/mathlib/Mathlib/RingTheory/DiscreteValuationRing/Basic.lean
  → def HasUnitMulPowIrreducibleFactorization [CommRing R] : Prop :=
      ∃ p : R, Irreducible p ∧ ∀ {x : R}, x ≠ 0 → ∃ n : ℕ, Associated (p ^ n) x
```

★**これは前波の `IntegerSubringNorm.exists_pow_mul_unit`
(`z ≠ 0`, `‖z‖ ≤ 1` ⇒ `z = π^n · w`, `‖w‖ = 1`)そのものである。**
⇒ 候補 1(`IsDiscreteValuationRing.ofHasUnitMulPowIrreducibleFactorization`、
`RingTheory/DiscreteValuationRing/Basic.lean:293`)を採った。
★候補 2(`of_ufd_of_unique_irreducible`、UFM が要る)と
候補 3(`Valuation.valuationSubring_isDiscreteValuationRing`、`Valued` 経由で
`K₀` と `integerSubring M` の同一視が要る)は**使わなかった**。

★`[IsDomain ↥(integerSubring M)]` は `inferInstance` で降りてくる(体の部分環)。
★`IsLocalRing` は入力に要らない(`ofHasUnitMulPowIrreducibleFactorization` が自分で作る)。

## 出したもの

| 宣言 | 内容 |
|---|---|
| ★`irreducible_pi` | `π` は `𝒪_M` の既約元(`π = ab` で両方が非単元なら `‖π‖ ≤ ‖π‖²`) |
| ★★★`isDiscreteValuationRing_integerSubring` | **`IsDiscreteValuationRing 𝒪_M`** |

★流儀は前波と同じで **`instance` にせず `theorem` ＋ `letI`/`haveI`**
(木も `LowerRamificationGroup.lean:669` で
`attribute [local instance] isDiscreteValuationRing_adjoinIntegers` と貼っている)。

## ★★残り(正確に、`file:line` つき)

`harith` の (3) を閉じるのに残っているのは
`HerbrandComposition.lean:455 herbrandPhiGroup_natCast` の**まだ揃っていない仮説**である。
実測した過不足:

| `herbrandPhiGroup_natCast` の要求 | 現状 |
|---|---|
| `[IsDomain B]` | ★済(`inferInstance`) |
| `[IsDiscreteValuationRing B]` | ★★**本ファイルで済** |
| `[MulSemiringAction G B]` | ★済(第 1122 `integerMulSemiringAction`) |
| `huni : maximalIdeal B = Ideal.span {α}` | ★済(第 1122 `maximalIdeal_eq_span`) |
| `[Algebra A B]`(`A = 𝒪_K`) | ★**未**。`𝒪_K` を建てて `𝒪_K → 𝒪_M` を作る |
| `[SMulCommClass G A B]` | ★**未** |
| `hadj : Algebra.adjoin A {α} = ⊤`(`𝒪_M = 𝒪_K[π]`) | ★**中身は在る**(第 1120 `IntegerSubringNorm.norm_algebraMap_le_one_of_le_one` ＋ `JumpMono.exists_coeff_norm_le`)が、**部分環の statement にしていない** |
| `[Fintype G]` | ★**未**(`G = ⟨g⟩` は有限なので容易) |

★さらに Hasse–Arf 本体
(`HasseArfStrongInduction.lean:447 exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top`)
は `hresA`(`𝒪_M = 𝒪_K + 𝔪`)・`hAinj`・`habel`・`h1 : lowerRamificationGroup B G 1 = ⊤` を要求する。
★`h1` は第 1122 の `mem_lowerRamificationGroup_iff_norm` と
`RamNormBridge.mem_ramification_iff` で出るはずだが★**測っていない**。

⇒ ★**次の 1 ノードは「`𝒪_K` を建てて `Algebra 𝒪_K 𝒪_M` と `𝒪_M = 𝒪_K[π]` を載せる」**である。

## 逸脱の記録

1. `instance` にせず `theorem`(`letI`/`haveI` で貼る)。全体 import への影響を避けるため。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**(import のみ)。
-/

namespace ABC3.Found.PGC

namespace IntegerNorm

section DVR

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★`π` は `𝒪_M` の既約元。 -/
theorem irreducible_pi {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m) (hπmem : π ∈ integerSubring M) :
    Irreducible (⟨π, hπmem⟩ : ↥(integerSubring M)) := by
  constructor
  · intro hu
    have := isUnit_iff_norm_eq_one.mp hu
    simp only at this
    linarith
  · intro a b hab
    by_contra hcon
    rw [not_or] at hcon
    obtain ⟨hna, hnb⟩ := hcon
    have ha1 : ‖(a : M)‖ < 1 := lt_of_le_of_ne a.2 (fun h => hna (isUnit_iff_norm_eq_one.mpr h))
    have hb1 : ‖(b : M)‖ < 1 := lt_of_le_of_ne b.2 (fun h => hnb (isUnit_iff_norm_eq_one.mpr h))
    have hale : ‖(a : M)‖ ≤ ‖π‖ := norm_le_norm_pi_of_lt_one hπ0 hπ1 hval ha1
    have hble : ‖(b : M)‖ ≤ ‖π‖ := norm_le_norm_pi_of_lt_one hπ0 hπ1 hval hb1
    have hmul : ‖π‖ = ‖(a : M)‖ * ‖(b : M)‖ := by
      have := congrArg (fun y : ↥(integerSubring M) => ‖(y : M)‖) hab
      simpa [norm_mul] using this
    nlinarith [norm_nonneg (a : M), norm_nonneg (b : M)]

/-- ★★★★★**`𝒪_M` は離散付値環**。中身は `exists_pow_mul_unit`。

★`instance` にせず `theorem` にして `letI`/`haveI` で貼る(前波の `isLocalRing_integerSubring` と
同じ流儀。★木も `attribute [local instance] isDiscreteValuationRing_adjoinIntegers` で貼っている)。 -/
theorem isDiscreteValuationRing_integerSubring {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m) (hπmem : π ∈ integerSubring M) :
    IsDiscreteValuationRing ↥(integerSubring M) := by
  refine IsDiscreteValuationRing.ofHasUnitMulPowIrreducibleFactorization
    ⟨⟨π, hπmem⟩, irreducible_pi hπ0 hπ1 hval hπmem, ?_⟩
  intro x hx0
  have hxne : (x : M) ≠ 0 := fun h => hx0 (Subtype.ext h)
  obtain ⟨n, w, hw1, hxw⟩ := exists_pow_mul_unit hπ0 hπ1 (hval (x : M) hxne) x.2 hxne
  have hwmem : w ∈ integerSubring M := by rw [mem_integerSubring, hw1]
  have hwunit : IsUnit (⟨w, hwmem⟩ : ↥(integerSubring M)) := isUnit_iff_norm_eq_one.mpr hw1
  obtain ⟨u, hu⟩ := hwunit
  refine ⟨n, u, ?_⟩
  apply Subtype.ext
  have hucoe : ((u : ↥(integerSubring M)) : M) = w := by rw [hu]
  simpa [hucoe] using hxw.symm

end DVR

/-! ## `.src` と 公理 -/

def isDiscreteValuationRing_integerSubring.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms irreducible_pi
#print axioms isDiscreteValuationRing_integerSubring

end IntegerNorm

end ABC3.Found.PGC
