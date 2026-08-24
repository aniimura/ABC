/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormNormal
import Mathlib.AlgebraicGeometry.Morphisms.Finite
import Mathlib.AlgebraicGeometry.Morphisms.Proper

/-!
# `V[L] → V` は有限、したがって `V[L]` は proper(鎖 `normalize` の `normalization-proper-normal`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> If we write V [L] for the normalization of V in L [so V [L] is also

## ★★正規性は済、残っていたのは proper の側だった

`normObj_isNormalScheme`(`NormNormal.lean`)で **`V[L]` は正規**が閉じている。
本ファイルは **proper** の側、すなわち

```
V[L] → V は有限射   ⟹   V[L] → V は proper   ⟹   V[L] → Spec k は proper
```

を閉じる。★`IsFinite ⟹ IsProper` と proper の合成は mathlib のインスタンス。

## ★★★中身は「整閉包が加群として有限」の 1 本

`normDown` が整射であることは mathlib のインスタンス(`IsIntegralHom fromNormalization`)。
★**有限型であること**が残っていて、それはアフィン局所で

```
Γ(V[L], D ⁻¹ᵁ U) ≅ integralClosure Γ(V,U) L   が Γ(V,U)-加群として有限
```

である。これは `IsIntegralClosure.finite`(mathlib)がそのまま与える。
★★ただしその仮定として **`Γ(V,U)` が整閉かつ Noether**、および
**`L/K` が分離的**が要る —— 原文の「proper normal variety」という仮定と、

原文 (FrdI p.109):
> be thought of as schemes Spec(L), where L ⊆ K is a finite [necessarily separable]

の「necessarily separable」に対応する。★節点としては **`V` 側の仮定を
明示的に引数に取る形**で書く(グラフの note が指示しているとおり)。

## ★★★★実装上の注意 —— 空のアフィン開

`IsFinite` の場は「**すべての**アフィン開 `U` について `f.app U` が有限」であり、
`U = ⊥` も含む。★`Γ(V,⊥)` と `Γ(V[L],⊥)` は自明環なので別扱いで潰す
(`Scheme.Hom.preimage_bot` で逆像が `⊥` になることを使う)。
★`NormNormal.lean` の側は `isNormalScheme_of_exists_affine`(点ごと)なので
この場合分けが要らなかった —— **場の形が違うと同じ道具が使えない**。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI

universe u

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]

/-! ## ★1. 空でないアフィン開の上で切断は有限 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★**`Γ(V[L], D ⁻¹ᵁ U)` は `Γ(V,U)`-加群として有限**(`U` は空でないアフィン開)。

★中身は `IsIntegralClosure.finite` 1 本。繋ぎは `NormNormal.lean` の
`sectionsIso` / `app_comp_sectionsIso` / `integralClosureCongr` をそのまま使う。 -/
theorem finite_normSections (L : FinSub V.functionField Kbar) {U : V.Opens}
    (hUaff : IsAffineOpen U) (hU : (U : Set V).Nonempty)
    (hnoeth : IsNoetherianRing Γ(V, U))
    (hic : IsIntegrallyClosed Γ(V, U))
    (hsep : Algebra.IsSeparable V.functionField L.toIF) :
    ((normDown V L).app U).hom.Finite := by
  haveI : Nonempty U := ⟨⟨hU.some, hU.some_mem⟩⟩
  haveI := L.fin
  haveI := hnoeth; haveI := hic; haveI := hsep
  letI : Algebra Γ(V, U) L.toIF :=
    ((algebraMap V.functionField L.toIF).comp (V.germToFunctionField U).hom).toAlgebra
  haveI : IsScalarTower Γ(V, U) V.functionField L.toIF :=
    IsScalarTower.of_algebraMap_eq' (by
      show ((algebraMap V.functionField L.toIF).comp (V.germToFunctionField U).hom) = _
      rfl)
  haveI := functionField_isFractionRing_of_isAffineOpen (X := V) U hUaff
  -- ★整閉包は加群として有限(mathlib)
  haveI hfin : Module.Finite Γ(V, U) (integralClosure Γ(V, U) L.toIF) :=
    IsIntegralClosure.finite Γ(V, U) V.functionField L.toIF (integralClosure Γ(V, U) L.toIF)
  letI := ((specToV V L).app U).hom.toAlgebra
  -- ★`Γ(Spec L, f ⁻¹ᵁ U) ≃ₐ[Γ(V,U)] L`(`NormNormal.lean` の繋ぎ)
  let e : Γ(Spec (CommRingCat.of L.toIF), (specToV V L) ⁻¹ᵁ U) ≃ₐ[Γ(V, U)] L.toIF :=
    { (sectionsIso V L hU).commRingCatIsoToRingEquiv with
      commutes' := fun r => by
        show (sectionsIso V L hU).hom.hom (((specToV V L).app U).hom r) = _
        exact congrArg (fun t : Γ(V, U) ⟶ CommRingCat.of L.toIF => CommRingCat.Hom.hom t r)
          (app_comp_sectionsIso V L hU) }
  haveI hfinC : Module.Finite Γ(V, U)
      (integralClosure Γ(V, U) Γ(Spec (CommRingCat.of L.toIF), (specToV V L) ⁻¹ᵁ U)) :=
    Module.Finite.equiv (integralClosureCongr e).symm.toLinearEquiv
  -- ★`fromNormalization.app U` は「整閉包への algebraMap」と同型の合成
  have happ := Scheme.Hom.fromNormalization_app (specToV V L) hUaff
  show ((specToV V L).fromNormalization.app U).hom.Finite
  rw [happ]
  show (((specToV V L).normalizationObjIso hUaff).inv.hom.comp
    (algebraMap Γ(V, U) _)).Finite
  exact RingHom.finite_respectsIso.1 _
    (((specToV V L).normalizationObjIso hUaff).commRingCatIsoToRingEquiv).symm
    (RingHom.finite_algebraMap.mpr hfinC)

/-- ★**空のアフィン開の場は自明環**なので有限。 -/
theorem finite_normSections_of_empty (L : FinSub V.functionField Kbar) {U : V.Opens}
    (hempty : (U : Set V) = ∅) :
    ((normDown V L).app U).hom.Finite := by
  have hbot : U = ⊥ := TopologicalSpace.Opens.ext (by simpa using hempty)
  subst hbot
  haveI : Subsingleton Γ(normObj V L, (normDown V L) ⁻¹ᵁ (⊥ : V.Opens)) := by
    rw [Scheme.Hom.preimage_bot]
    infer_instance
  letI := ((normDown V L).app (⊥ : V.Opens)).hom.toAlgebra
  exact ⟨⟨{0}, Subsingleton.elim _ _⟩⟩

/-! ## ★2. `V[L] → V` は有限射 -/

/-- ★★★★★★**`V[L] → V` は有限射**。

★仮定は `V` が正規 Noether(アフィン開の切断が整閉かつ Noether)であることと、
`L/K` が分離的であること —— 原文の「proper normal variety」と
「necessarily separable」に対応する。 -/
theorem isFinite_normDown (L : FinSub V.functionField Kbar)
    (hnoeth : ∀ U : V.Opens, IsAffineOpen U → (U : Set V).Nonempty → IsNoetherianRing Γ(V, U))
    (hic : ∀ U : V.Opens, IsAffineOpen U → (U : Set V).Nonempty → IsIntegrallyClosed Γ(V, U))
    (hsep : Algebra.IsSeparable V.functionField L.toIF) :
    IsFinite (normDown V L) where
  finite_app U hU := by
    rcases Set.eq_empty_or_nonempty (U : Set V) with h | h
    · exact finite_normSections_of_empty V L h
    · exact finite_normSections V L hU h (hnoeth U hU h) (hic U hU h) hsep

/-! ## ★3. `V[L]` は proper -/

/-- ★★**`V[L] → V` は proper**(有限射だから)。 -/
theorem isProper_normDown (L : FinSub V.functionField Kbar)
    (hnoeth : ∀ U : V.Opens, IsAffineOpen U → (U : Set V).Nonempty → IsNoetherianRing Γ(V, U))
    (hic : ∀ U : V.Opens, IsAffineOpen U → (U : Set V).Nonempty → IsIntegrallyClosed Γ(V, U))
    (hsep : Algebra.IsSeparable V.functionField L.toIF) :
    IsProper (normDown V L) := by
  haveI := isFinite_normDown V L hnoeth hic hsep
  infer_instance

/-- ★★★★★★**`V[L]` も proper normal variety である**(原文の角括弧の中身)。

原文 (FrdI p.109):
> a proper normal variety], then let us write DL for the set of prime divisors of V [L]

★正規性は `normObj_isNormalScheme`(`NormNormal.lean`、仮定なし)。
★proper は `V → S` が proper であることと合成する。 -/
theorem isProper_normObj {S : Scheme.{u}} (L : FinSub V.functionField Kbar) (g : V ⟶ S)
    (hg : IsProper g)
    (hnoeth : ∀ U : V.Opens, IsAffineOpen U → (U : Set V).Nonempty → IsNoetherianRing Γ(V, U))
    (hic : ∀ U : V.Opens, IsAffineOpen U → (U : Set V).Nonempty → IsIntegrallyClosed Γ(V, U))
    (hsep : Algebra.IsSeparable V.functionField L.toIF) :
    IsProper (normDown V L ≫ g) ∧ IsNormalScheme (normObj V L) := by
  haveI := isProper_normDown V L hnoeth hic hsep
  haveI := hg
  exact ⟨inferInstance, normObj_isNormalScheme V L⟩

/-! ### ★出典の紐付け(`.src`) -/

/-- ★★★★★★locator —— `Example 6.1` の「`V[L]` は proper normal variety である」。 -/
def isProper_normObj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — V[L] は proper normal variety である",
    sectionId := "frdi-example-6-1" }

def isFinite_normDown.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — V[L] → V は有限射",
    sectionId := "frdi-example-6-1" }

def finite_normSections.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 整閉包は加群として有限",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
