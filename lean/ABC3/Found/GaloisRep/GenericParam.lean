/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.PowerSeriesUniversal
import ABC3.Meta.Claim
import Mathlib.RingTheory.Localization.Away.Basic
import Mathlib.RingTheory.Polynomial.Basic

/-!
# 第 1143 ブロック —— **一般の径数 `T` を載せる環**（`Found`、`l = 2` の枝の節点 1）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★これは何か（`l = 2` の枝の節点 1 の土台）

`veluV2_eq_tateDYpair`（`TateODE.lean:89`）は `DX ≠ 0` で割っている。
★`l = 2` では点が **2-捩れ**なので `DX = 2Y + X = 0` になり、この形は使えない。

☆しかし `DY = veluV2` は **`a` を径数と見れば恒等式**である。★そこで

    `U ≔ ℤ[T]` を `T·(1 − T)` で局所化した環

を作る。`U` では `T` も `1 − T` も単元で、しかも `1 + T ≠ 0` である
（`ℤ[T] → U` は単射だから）。☆したがって `PowerSeries U` の中で

    `DX` の定数項 `= T(1+T)·(1−T)⁻³ ≠ 0`

となり、**在庫の `veluV2_eq_tateDYpair` がそのまま効く**。

★あとは `T ↦ a` で特殊化すればよい。`a·(1 − a)` が単元なら
`IsLocalization.lift` が環準同型 `U →+* R` を与える。

☆これは第 1128 と同じ「万有な環を経由する」型である。
-/

namespace ABC3.Found.GaloisRep

open Polynomial

/-- ★★★★★★**一般の径数の環**——`ℤ[T]` を `T(1−T)` で局所化したもの。 -/
noncomputable abbrev GenericParamRing :=
  Localization.Away ((Polynomial.X : Polynomial ℤ) * (1 - Polynomial.X))

/-- ☆`T(1−T) ≠ 0`。 -/
theorem genericParam_ne_zero :
    ((Polynomial.X : Polynomial ℤ) * (1 - Polynomial.X)) ≠ 0 := by
  intro h
  rcases mul_eq_zero.1 h with h1 | h2
  · exact Polynomial.X_ne_zero h1
  · have hc := congrArg (fun p => Polynomial.coeff p 1) h2
    simp [Polynomial.coeff_one] at hc

noncomputable instance : IsDomain GenericParamRing :=
  IsLocalization.isDomain_localization
    (powers_le_nonZeroDivisors_of_noZeroDivisors genericParam_ne_zero)

/-- ☆`ℤ[T] → U` は単射（非零因子での局所化だから）。 -/
theorem genericParam_injective :
    Function.Injective (algebraMap (Polynomial ℤ) GenericParamRing) :=
  IsLocalization.injective _
    (powers_le_nonZeroDivisors_of_noZeroDivisors genericParam_ne_zero)

/-- ★★**径数 `T` そのもの**。 -/
noncomputable def genericT : GenericParamRing :=
  algebraMap (Polynomial ℤ) GenericParamRing Polynomial.X

/-- ★★`T` は単元。 -/
theorem isUnit_genericT : IsUnit genericT := by
  have hf : IsUnit (algebraMap (Polynomial ℤ) GenericParamRing
      ((Polynomial.X : Polynomial ℤ) * (1 - Polynomial.X))) :=
    IsLocalization.Away.algebraMap_isUnit _
  rw [map_mul] at hf
  exact isUnit_of_mul_isUnit_left hf

/-- ★★`1 − T` は単元——★これが「`a` を径数にすると `hu` が自動になる」の中身である。 -/
theorem isUnit_one_sub_genericT : IsUnit (1 - genericT) := by
  have hf : IsUnit (algebraMap (Polynomial ℤ) GenericParamRing
      ((Polynomial.X : Polynomial ℤ) * (1 - Polynomial.X))) :=
    IsLocalization.Away.algebraMap_isUnit _
  rw [map_mul] at hf
  have h2 : IsUnit (algebraMap (Polynomial ℤ) GenericParamRing (1 - Polynomial.X)) :=
    isUnit_of_mul_isUnit_right hf
  rwa [map_sub, map_one] at h2

/-- ★★★★**`1 + T ≠ 0`**——★これが「`a = −1` の場合を径数で迂回できる」理由である。

☆`ℤ[T] → U` は単射で、`1 + T` は `ℤ[T]` の中で `0` でない。 -/
theorem one_add_genericT_ne_zero : (1 : GenericParamRing) + genericT ≠ 0 := by
  intro h
  have h0 : algebraMap (Polynomial ℤ) GenericParamRing (1 + Polynomial.X)
      = algebraMap _ _ 0 := by
    rw [map_add, map_one, map_zero]
    exact h
  have hpoly := genericParam_injective h0
  have hc := congrArg (fun p => Polynomial.coeff p 1) hpoly
  simp [Polynomial.coeff_one] at hc

/-- ★★★★★★★★**特殊化**——`a·(1 − a)` が単元なら `T ↦ a` が環準同型を与える。 -/
noncomputable def genericSpec {R : Type} [CommRing R] (a : R) (ha : IsUnit (a * (1 - a))) :
    GenericParamRing →+* R :=
  IsLocalization.lift
    (M := Submonoid.powers ((Polynomial.X : Polynomial ℤ) * (1 - Polynomial.X)))
    (g := (Polynomial.aeval a).toRingHom.comp (Polynomial.mapRingHom (Int.castRingHom R)))
    (by
      rintro ⟨y, n, hn⟩
      simp only [← hn, map_pow]
      refine IsUnit.pow n ?_
      simpa using ha)

/-- ★★**特殊化は `T` を `a` に送る**。 -/
@[simp] theorem genericSpec_T {R : Type} [CommRing R] (a : R) (ha : IsUnit (a * (1 - a))) :
    genericSpec a ha genericT = a := by
  rw [genericSpec, genericT, IsLocalization.lift_eq]
  simp

/-! ## ★出典の紐付け(`.src`) -/

def GenericParamRing.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(一般の径数を載せる環——ℤ[T] を T(1−T) で局所化)",
    sectionId := "genell-def-3-3" }

def genericSpec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(径数の特殊化 T ↦ a——a(1−a) が単元なら環準同型)",
    sectionId := "genell-def-3-3" }

def genericSpec.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsLocalization.lift(局所化の普遍性)"
      (.inMathlib "IsLocalization.lift") 1,
    .citation "[mathlib]" "IsLocalization.isDomain_localization(非零因子での局所化は整域)"
      (.inMathlib "IsLocalization.isDomain_localization") 1,
    .implicitStep
      ("★★**2026-09-01（第 1143）**——`l = 2` では点が 2-捩れなので `DX = 0` になり、" ++
       "`veluV2_eq_tateDYpair` の「`DX ≠ 0` で割る」段が使えない。" ++
       "☆`a` を径数 `T` に取れば `T`・`1−T` は単元、`1+T ≠ 0` なので " ++
       "`DX` の定数項 `T(1+T)(1−T)⁻³` が `0` でなく、在庫がそのまま効く。" ++
       "★`T ↦ −1` は `(−1)·(1−(−1)) = −2` が単元なら定義される。") 6 ]

end ABC3.Found.GaloisRep
