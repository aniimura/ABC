/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.GlobalUnits
import ABC3.Found.Divisor.Hartogs
import ABC3.Found.Divisor.SchemeOrdUnit

/-!
# 代数的 Hartogs のスキーム版(鎖 `normalize` の `global-units` の 1 段目)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> that [since V [L] is a proper normal variety] for A ∈Ob(CV,

## ★★何を閉じるか

環の側の代数的 Hartogs(`Found/Divisor/Hartogs.lean`)を**スキームに載せる**:

```
ord_v(u) ≥ 0 が余次元 1 の点すべてで成り立つ  ⟹  u は大域切断
```

★★これで `global-units` の 3 段(環の Hartogs / 貼り合わせ / 配線)が揃う。

## ★★★配線は 3 本だった

| 本 | 使う在庫 |
|---|---|
| 高さ 1 の素イデアル `p` ↦ `X` の点 `hU.fromSpec p` が余次元 1 | `isLocalization_stalk'` ＋ `IsLocalization.AtPrime.ringKrullDim_eq_height` |
| `ord ≥ 0` ⟹ 茎に入る | `exists_mem_stalk_of_ordPt_nonneg`(実装済み) |
| 茎 `= R_p` から「分母が `p` の外」を取り出す | `IsLocalization.surj` ＋ `functionField_isScalarTower` |

## ★★★★実装上の罠 —— 茎の `Algebra` は明示に入れる

mathlib の `TopCat.Presheaf.algebra_section_stalk` は **`x : U`(部分型の元)**
の形で instance 登録されているので、`hU.fromSpec.base p : X`(素の点)の形では
**instance 探索が届かない**。`letI` で入れてから `IsLocalization.surj` と
`IsScalarTower` を使うこと。★`haveI` では駄目で `letI` が要る
(`isLocalization_stalk'` の主張に出てくる項と**同じ項**でなければ合わない)。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Meta

universe u

/-! ## ★1. 高さ 1 の素イデアルは余次元 1 の点を与える -/

/-- ★★★**アフィン開の高さ 1 の素イデアルは `X` の余次元 1 の点**。

★茎が局所化(`isLocalization_stalk'`)で、局所化の Krull 次元は高さ。 -/
theorem isCodimOnePt_fromSpec {X : Scheme.{u}} {U : X.Opens} (hU : IsAffineOpen U)
    (p : PrimeSpectrum Γ(X, U)) (hp : p.asIdeal.height = 1) :
    IsCodimOnePt X (hU.fromSpec.base p) := by
  have hmem : hU.fromSpec.base p ∈ U := by
    have hr : hU.fromSpec.base p ∈ Set.range hU.fromSpec.base := ⟨p, rfl⟩
    rwa [hU.range_fromSpec] at hr
  have hloc := hU.isLocalization_stalk' p hmem
  have h1 : ringKrullDim (X.presheaf.stalk (hU.fromSpec.base p))
      = (p.asIdeal.height : WithBot ℕ∞) :=
    @IsLocalization.AtPrime.ringKrullDim_eq_height _ _ p.asIdeal p.isPrime
      (X.presheaf.stalk (hU.fromSpec.base p)) _
      (TopCat.Presheaf.algebra_section_stalk X.presheaf ⟨_, hmem⟩) hloc
  rw [IsCodimOnePt, h1, hp]
  rfl

/-! ## ★2. アフィン局所の Hartogs -/

/-- ★★★★**アフィン開ごとの Hartogs** —— `ord ≥ 0` がすべての余次元 1 の点で
成り立つなら、`u` は `Γ(X,U)` から来る。 -/
theorem mem_range_germ_affine {X : Scheme.{u}} [IsIntegral X] [IsLocallyNoetherian X]
    (hnorm : IsNormalScheme X) {U : X.Opens} (hU : IsAffineOpen U) (hne : Nonempty U)
    (hnoeth : IsNoetherianRing Γ(X, U)) (hic : IsIntegrallyClosed Γ(X, U))
    (u : X.functionField)
    (h : ∀ x : PrimeDivisorPt X, 0 ≤ ordPt X hnorm x u) :
    u ∈ Set.range (ConcreteCategory.hom (X.germToFunctionField U)) := by
  rcases eq_or_ne u 0 with rfl | hu
  · exact ⟨0, map_zero _⟩
  haveI := hne; haveI := hnoeth; haveI := hic
  haveI := functionField_isFractionRing_of_isAffineOpen (X := X) U hU
  have hmem : u ∈ Hartogs.Rsub Γ(X, U) X.functionField := by
    refine Hartogs.mem_Rsub_of_forall_heightOnePrime u (fun v => ?_)
    set p : PrimeSpectrum Γ(X, U) := ⟨v.asIdeal, v.isPrime⟩ with hpdef
    haveI : p.asIdeal.IsPrime := p.isPrime
    have hmemU : hU.fromSpec.base p ∈ U := by
      have hr : hU.fromSpec.base p ∈ Set.range hU.fromSpec.base := ⟨p, rfl⟩
      rwa [hU.range_fromSpec] at hr
    have hcodim : IsCodimOnePt X (hU.fromSpec.base p) :=
      isCodimOnePt_fromSpec hU p v.height_eq_one
    obtain ⟨r, hr⟩ := exists_mem_stalk_of_ordPt_nonneg hnorm ⟨_, hcodim⟩ hu (h ⟨_, hcodim⟩)
    -- ★茎の `Algebra` は明示に入れる(素の点では instance 探索が届かない)
    letI algInst : Algebra Γ(X, U) (X.presheaf.stalk (hU.fromSpec.base p)) :=
      TopCat.Presheaf.algebra_section_stalk X.presheaf ⟨_, hmemU⟩
    haveI hloc : IsLocalization.AtPrime (X.presheaf.stalk (hU.fromSpec.base p)) p.asIdeal :=
      hU.isLocalization_stalk' p hmemU
    haveI hst : IsScalarTower Γ(X, U) (X.presheaf.stalk (hU.fromSpec.base p)) X.functionField :=
      functionField_isScalarTower X U ⟨_, hmemU⟩
    obtain ⟨⟨a, b⟩, hab⟩ := IsLocalization.surj (M := p.asIdeal.primeCompl) r
    refine ⟨(b : Γ(X, U)), ?_, b.2⟩
    rw [Hartogs.mem_den]
    refine ⟨a, ?_⟩
    show (algebraMap Γ(X, U) X.functionField) a = (b : Γ(X, U)) • u
    have hpush := congrArg (algebraMap (X.presheaf.stalk (hU.fromSpec.base p)) X.functionField) hab
    rw [map_mul, hr] at hpush
    rw [← IsScalarTower.algebraMap_apply, ← IsScalarTower.algebraMap_apply] at hpush
    rw [← hpush, Algebra.smul_def]
    ring
  obtain ⟨a, ha⟩ := hmem
  exact ⟨a, ha⟩

/-! ## ★3. 本体 —— スキーム版の Hartogs -/

/-- ★★★★★★**代数的 Hartogs のスキーム版** ——
`ord_v(u) ≥ 0` がすべての余次元 1 の点で成り立つなら `u` は**大域切断**。

★アフィン開ごとに環の Hartogs を当て(`mem_range_germ_affine`)、
貼り合わせる(`mem_range_germ_top_of_forall_affine`)。 -/
theorem mem_range_germ_top_of_forall_ordPt_nonneg
    {X : Scheme.{u}} [IsIntegral X] [IsLocallyNoetherian X]
    (hnorm : IsNormalScheme X)
    (hnoeth : ∀ U : X.Opens, IsAffineOpen U → Nonempty U → IsNoetherianRing Γ(X, U))
    (hic : ∀ U : X.Opens, IsAffineOpen U → Nonempty U → IsIntegrallyClosed Γ(X, U))
    (u : X.functionField)
    (h : ∀ x : PrimeDivisorPt X, 0 ≤ ordPt X hnorm x u) :
    u ∈ Set.range (ConcreteCategory.hom (X.germToFunctionField ⊤)) :=
  mem_range_germ_top_of_forall_affine u
    (fun U hU hne => mem_range_germ_affine hnorm hU hne (hnoeth U hU hne) (hic U hU hne) u h)

/-! ## ★4. 出口 —— `𝒪^×(A) = k_L^×` -/

/-- ★★★★★★**[FrdI] Example 6.1** —— proper normal なら
「すべての余次元 1 の点で `ord = 0`」の元は `Γ(X,⊤)` の**単元**である。

★★`Γ(X,⊤)` は `k` 上有限次の体(`globalSections_isField` / `globalSections_finite`)
なので、これが原文の `𝒪^×(A) = 𝒪^▷(A) = k_L^×` である。 -/
theorem exists_unit_of_forall_ordPt_eq_zero
    {X : Scheme.{u}} [IsIntegral X] [IsLocallyNoetherian X] {k : Type u} [Field k]
    (g : X ⟶ Spec (CommRingCat.of k)) (hg : IsProper g)
    (hnorm : IsNormalScheme X)
    (hnoeth : ∀ U : X.Opens, IsAffineOpen U → Nonempty U → IsNoetherianRing Γ(X, U))
    (hic : ∀ U : X.Opens, IsAffineOpen U → Nonempty U → IsIntegrallyClosed Γ(X, U))
    (u : X.functionField) (hu : u ≠ 0)
    (h : ∀ x : PrimeDivisorPt X, ordPt X hnorm x u = 0) :
    ∃ t : Γ(X, ⊤), IsUnit t ∧ (ConcreteCategory.hom (X.germToFunctionField ⊤)) t = u :=
  exists_unit_globalSection g hg hu
    (mem_range_germ_top_of_forall_ordPt_nonneg hnorm hnoeth hic u (fun x => (h x).ge))

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def mem_range_germ_top_of_forall_ordPt_nonneg.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — 代数的 Hartogs のスキーム版",
    sectionId := "frdi-example-6-1" }

def mem_range_germ_top_of_forall_ordPt_nonneg.needs : List ProofObligation :=
  [ .citation "[ABC3]" "代数的 Hartogs(環の側、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.Hartogs.mem_Rsub_of_forall_heightOnePrime") 110,
    .citation "[ABC3]" "貼り合わせ(sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.mem_range_germ_top_of_forall_affine") 110,
    .citation "[ABC3]" "ord ≥ 0 なら茎に入る"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_mem_stalk_of_ordPt_nonneg") 110,
    .citation "[mathlib]" "IsAffineOpen.isLocalization_stalk'(茎は局所化)"
      (.inMathlib "AlgebraicGeometry.IsAffineOpen.isLocalization_stalk'") 110,
    .citation "[mathlib]" "functionField_isScalarTower(Γ(X,U) → 茎 → K(X) の両立)"
      (.inMathlib "AlgebraicGeometry.functionField_isScalarTower") 110,
    .derivation "高さ 1 の素イデアルに対応する点で ord ≥ 0 を使い、分母が素イデアルの外にあることを取り出す" 110 ]

def isCodimOnePt_fromSpec.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — アフィン開の高さ 1 の素は X の余次元 1 の点",
    sectionId := "frdi-example-6-1" }

def isCodimOnePt_fromSpec.needs : List ProofObligation :=
  [ .citation "[mathlib]" "IsLocalization.AtPrime.ringKrullDim_eq_height"
      (.inMathlib "IsLocalization.AtPrime.ringKrullDim_eq_height") 110,
    .citation "[mathlib]" "IsAffineOpen.isLocalization_stalk'"
      (.inMathlib "AlgebraicGeometry.IsAffineOpen.isLocalization_stalk'") 110 ]

def mem_range_germ_affine.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — アフィン開ごとの Hartogs",
    sectionId := "frdi-example-6-1" }

def mem_range_germ_affine.needs : List ProofObligation :=
  [ .citation "[ABC3]" "Hartogs.mem_Rsub_of_forall_heightOnePrime"
      (.inProject "ABC3" "ABC3.Found.Divisor.Hartogs.mem_Rsub_of_forall_heightOnePrime") 110,
    .citation "[mathlib]" "functionField_isFractionRing_of_isAffineOpen"
      (.inMathlib "AlgebraicGeometry.functionField_isFractionRing_of_isAffineOpen") 110,
    .derivation "茎への Algebra は letI で明示に入れる(素の点では instance 探索が届かない)" 110 ]

def exists_unit_of_forall_ordPt_eq_zero.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — O^×(A) = O^▷(A) = k_L^×",
    sectionId := "frdi-example-6-1" }

def exists_unit_of_forall_ordPt_eq_zero.needs : List ProofObligation :=
  [ .citation "[ABC3]" "mem_range_germ_top_of_forall_ordPt_nonneg"
      (.inProject "ABC3" "ABC3.Found.Divisor.mem_range_germ_top_of_forall_ordPt_nonneg") 110,
    .citation "[ABC3]" "exists_unit_globalSection(Γ(X,⊤) は体)"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_unit_globalSection") 110,
    .implicitStep
      "★原文は「[since V [L] is a proper normal variety]」の角括弧ひとつで畳む" 110 ]

end ABC3.Found.Divisor
