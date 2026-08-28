/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ImmersionGlobalToProj
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★埋め込み性は**部分族**で確かめれば足りる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★★★★★★これは何か —— 段 E3 の**障害を外す**

`§9-913` の `isImmersion_globalToProj` は **`Fin (N+1)` のすべての `i`** について
`haff`（`X_{s_i}` がアフィン）と `hsurj`（チャート写像が全射）を要求していた。

★★★★**そこに穴があった**（2026-08-28 実測）。段 E3 の全射性を出すには
`§9-CommonGluedRatio` の `exists_common_glued_globalRatio` で
**分子の切断 `t_k` を族に足す**必要があるが、

* 分母 `s_i` の非消失軌跡は `IsAmple` の定義からアフィン（`§9-915`）
* ★**分子 `t_k` の非消失軌跡はアフィンとは限らない**

——`X_{t_k} ⊓ X_{s_i}` はアフィンの基本開だからアフィンだが、
`X_{t_k}` 自身はそれらの**和集合**であってアフィンとは限らない。

★★★★★**本ファイルがその穴を塞ぐ**: `haff`・`hsurj` は
**`X` を覆う部分族についてだけ**あればよい。

## ★★★機構 —— 像の入る開へ落とす

`V ≔ ⨆_{i ∈ I₀} D₊(x_i)` と置く。★`§9-913` の `ψ⁻¹(D₊(x_i)) = X_{s_i}` から

    `ψ⁻¹(V) = ⨆_{i ∈ I₀} X_{s_i} = ⊤`

なので `ψ` は `V` を通る（`IsOpenImmersion.lift`）。★★`ψ' : X ⟶ V` について
`{V.ι⁻¹(D₊(x_i))}_{i ∈ I₀}` は `V` の被覆だから、
`IsZariskiLocalAtTarget` が使える。
★★★各チャートでは `morphismRestrict_comp` と `IsImmersion.of_comp` で
`ψ ∣_ D₊(x_i)`（`§9-913` で埋め込み）に戻る。

## ★配管の記録

★★★★`rw [hfac]`（`hfac : ψ' ≫ V.ι = ψ`）は **motive が型付かない**
——`f ∣_ U` の**始域が `f` に依存する**からである。
★**汎化してから `rintro rfl`** すれば通る:

```lean
have hgen : ∀ χ, χ = ψ → IsImmersion (χ ∣_ U) := by rintro χ rfl; …
exact hgen _ hfac
```
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★**埋め込み性は部分族で確かめれば足りる**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`haff`（チャートがアフィン）と `hsurj`（チャート写像が全射）は、
**`X` を覆う部分族 `{i | p i}` についてだけ**あればよい。

★★これで段 E3 の全射性のために**分子の切断を族に足しても**埋め込み性が壊れない
——分子の非消失軌跡はアフィンとは限らないからである。 -/
theorem isImmersion_globalToProj_of_subfamily (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (p : Fin (N + 1) → Prop)
    (hcov₀ : (⨆ i : {i // p i}, nonVanishing M (s i.1)) = ⊤)
    (haff : ∀ i : {i // p i}, IsAffineOpen (nonVanishing M (s i.1)))
    (hsurj : ∀ i : {i // p i}, Function.Surjective (globalAwayHom M hM φ s i.1)) :
    IsImmersion (globalToProj M hM φ s hcov) := by
  set 𝒜 := MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R with h𝒜
  set V : (projSpace N R).Opens :=
    ⨆ i : {i // p i}, Proj.basicOpen 𝒜 (MvPolynomial.X i.1) with hV
  have hpre : globalToProj M hM φ s hcov ⁻¹ᵁ V = ⊤ := by
    rw [hV, Scheme.Hom.preimage_iSup]
    have hb : ∀ i : {i // p i},
        globalToProj M hM φ s hcov ⁻¹ᵁ (Proj.basicOpen 𝒜 (MvPolynomial.X i.1))
        = nonVanishing M (s i.1) := fun i => globalToProj_preimage_basicOpen M hM φ s hcov i.1
    simp only [hb]
    exact hcov₀
  have hrange : Set.range (globalToProj M hM φ s hcov).base ⊆ Set.range V.ι.base := by
    rw [Scheme.Opens.range_ι]
    rintro y ⟨x, rfl⟩
    have hx : x ∈ globalToProj M hM φ s hcov ⁻¹ᵁ V := by rw [hpre]; trivial
    exact hx
  set ψ' := IsOpenImmersion.lift V.ι (globalToProj M hM φ s hcov) hrange with hψ'
  have hfac : ψ' ≫ V.ι = globalToProj M hM φ s hcov := IsOpenImmersion.lift_fac _ _ _
  have himm : IsImmersion ψ' := by
    refine IsZariskiLocalAtTarget.of_iSup_eq_top
      (fun i : {i // p i} => V.ι ⁻¹ᵁ (Proj.basicOpen 𝒜 (MvPolynomial.X i.1))) ?_ ?_
    · rw [← Scheme.Hom.preimage_iSup, ← hV]
      exact V.ι_preimage_self
    · intro i
      have hgen : ∀ (χ : X ⟶ projSpace N R), χ = globalToProj M hM φ s hcov →
          IsImmersion (χ ∣_ (Proj.basicOpen 𝒜 (MvPolynomial.X i.1))) := by
        rintro χ rfl
        haveI : IsImmersion (globalToProj M hM φ s hcov
            ∣_ (Proj.basicOpen 𝒜 (MvPolynomial.X i.1))
            ≫ (Proj.basicOpen 𝒜 (MvPolynomial.X i.1)).ι) := by
          rw [morphismRestrict_ι, globalToProj_preimage_basicOpen M hM φ s hcov i.1,
            ι_globalToProj]
          exact isImmersion_globalChartToProj M hM φ s i.1 (haff i) (hsurj i)
        exact IsImmersion.of_comp _ (Proj.basicOpen 𝒜 (MvPolynomial.X i.1)).ι
      haveI hbig : IsImmersion (ψ' ∣_ (V.ι ⁻¹ᵁ (Proj.basicOpen 𝒜 (MvPolynomial.X i.1)))
          ≫ V.ι ∣_ (Proj.basicOpen 𝒜 (MvPolynomial.X i.1))) := by
        rw [← morphismRestrict_comp]
        exact hgen _ hfac
      exact IsImmersion.of_comp _ (V.ι ∣_ (Proj.basicOpen 𝒜 (MvPolynomial.X i.1)))
  rw [← hfac]
  infer_instance

/-! ## ★出典の紐付け(`.src`) -/

def isImmersion_globalToProj_of_subfamily.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(埋め込み性は部分族で確かめれば足りる)",
    sectionId := "genell-prop-1-4" }

def isImmersion_globalToProj_of_subfamily.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "globalToProj_preimage_basicOpen(ψ⁻¹(D₊(x_i)) = X_{s_i}、§9-913)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalToProj_preimage_basicOpen") 2,
    .citation "[mathlib]" "IsOpenImmersion.lift(像が入る開へ落とす)"
      (.inMathlib "AlgebraicGeometry.IsOpenImmersion.lift") 1,
    .citation "[mathlib]" "morphismRestrict_comp・IsImmersion.of_comp"
      (.inMathlib "AlgebraicGeometry.morphismRestrict_comp") 1,
    .implicitStep
      ("★★★★測定(2026-08-28): §9-913 の形には**穴があった**。" ++
       "段 E3 の全射性を出すには exists_common_glued_globalRatio で" ++
       "**分子の切断 t_k を族に足す**必要があるが、" ++
       "分母 s_i の非消失軌跡は IsAmple からアフィン(§9-915)でも、" ++
       "**分子 t_k の非消失軌跡はアフィンとは限らない**" ++
       "——X_{t_k} ⊓ X_{s_i} はアフィンの基本開だからアフィンだが、" ++
       "X_{t_k} 自身はそれらの和集合である") 5,
    .implicitStep
      ("★★★配管: rw [hfac](hfac : ψ' ≫ V.ι = ψ)は motive が型付かない" ++
       "——f ∣_ U の**始域が f に依存する**からである。" ++
       "汎化してから rintro rfl すれば通る") 2 ]

end ABC3.Found.GenEll
