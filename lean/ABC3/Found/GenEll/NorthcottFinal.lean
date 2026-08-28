/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottGeom
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★`Proposition 1.4, (iv)` —— 幾何の仮定も外れた（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 最後の仮定が消えた

`§9-937` まで、`Proposition 1.4, (iv)` の鎖には**幾何の側の仮定**

    `U`・`hcovU`・`hU`・`eU`——`M` を自明化するアフィン開被覆

が残っていた。★★★★これは `IsLocallyTrivial X M` から**そのまま出る**:

1. `IsLocallyTrivial X M` は `⊤` の上に被覆篩 `S` を与え、`S` の各元の上で `M ≅ 𝒪` である
2. アフィン開は位相の基底である（`Scheme.isBasis_affineOpens`）
3. 篩は前合成で閉じているので、`S` の元に含まれるアフィン開もまた `S` の元である
4. ★したがって「`S` に属するアフィン開」の全体が求める被覆である

★★★★★**有限性は要らない**——`northcott_of_isAmple_geom` の添字型 `ιU` に
有限性の要求がないからである。

## ★★★★★★到達 —— `Proposition 1.4, (iv)` の仮定

```lean
theorem northcott_of_isAmple_final
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)      -- 可逆層
    (hample : IsAmple M) (x₀ : (X : Type))                     -- 豊富、点がある
    (φ : ℤ →+* Γ(X, ⊤)) (d : ℕ)                                -- 構造射、次数の上限
```

——★★★★★★**これだけである**。点の側の仮定は `§9-937` で全部幾何の言葉になり、
幾何の側の仮定は本ファイルで消えた。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open MvPolynomial HomogeneousLocalization NumberField

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★★自明化するアフィン開被覆は存在する -/

/-- ★★★★★★★★★★★★★**可逆層は自明化するアフィン開被覆を持つ**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★これが `§9-927`〜`§9-937` の鎖が受けていた**幾何の側の仮定**そのものであり、
`IsLocallyTrivial X M` から出る。

★機構は 3 つだけ:
* `IsLocallyTrivial` は `⊤` の上に**被覆篩** `S` を与える
* アフィン開は**位相の基底**である（`Scheme.isBasis_affineOpens`）
* 篩は**前合成で閉じている**ので、`S` の元に含まれるアフィン開も `S` の元である -/
theorem exists_affine_trivializing_cover (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) :
    ∃ (ι : Type) (U : ι → X.Opens), (⨆ j, U j) = ⊤ ∧ (∀ j, IsAffineOpen (U j)) ∧
      ∀ j, Nonempty ((restrictPresheafFunctor X (U j)).obj M
        ≅ 𝟙_ (PresheafModulesOn X (U j))) := by
  obtain ⟨S, hS, htriv⟩ := hM ⊤
  refine ⟨{W : X.Opens // IsAffineOpen W ∧ S.arrows (homOfLE (le_top : W ≤ ⊤))},
    fun j => j.1, ?_, fun j => j.2.1, fun j => htriv _ j.2.2⟩
  refine top_le_iff.mp (fun x _ => ?_)
  obtain ⟨V, f, hf, hxV⟩ := hS x trivial
  obtain ⟨W, hWmem, hxW, hWV⟩ :=
    (TopologicalSpace.Opens.isBasis_iff_nbhd.mp (Scheme.isBasis_affineOpens X)) hxV
  have hSW : S.arrows (homOfLE (le_top : W ≤ ⊤)) := by
    have h := S.downward_closed hf (homOfLE hWV)
    have he : (homOfLE (le_top : W ≤ ⊤)) = homOfLE hWV ≫ f := Subsingleton.elim _ _
    rw [he]
    exact h
  exact TopologicalSpace.Opens.mem_iSup.mpr ⟨⟨W, ⟨hWmem, hSW⟩⟩, hxW⟩

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★到達 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**`IsAmple` から Northcott まで一本**
—— `Proposition 1.4, (iv)` の仮定は可逆層・豊富性・点・構造射・次数の上限だけ。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★★★★`§9-937` で点の側の仮定は全部幾何の言葉になり、
本定理で**幾何の側の仮定（自明化するアフィン開被覆）も消えた**。 -/
theorem northcott_of_isAmple_final (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (hample : IsAmple M) (x₀ : (X : Type))
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type)) (d : ℕ) :
    ∃ (L : ℕ) (_ : 0 < L) (N : ℕ)
      (s : Fin (N + 1) → ((presheafTensorPow M L).obj (op ⊤) : Type))
      (hcov : (⨆ k, nonVanishing (presheafTensorPow M L) (s k)) = ⊤),
      (∀ i, IsAffineOpen (nonVanishing (presheafTensorPow M L) (s i))) ∧
      ∀ (D : ArithCartier X)
        (_ : D.divisor = divisorOfSection (presheafTensorPow M L) (s 0))
        (_ : D.green = fun p => greenFS N (p ≫ globalToProj (presheafTensorPow M L)
          (isLocallyTrivial_presheafTensorPow hM L) φ s hcov))
        {P : Type}
        (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
        (_ : ∀ p, Module.finrank ℚ (fld p) ≤ d)
        (xF : ∀ p, haveI := hnf p; specRingOfIntegers (fld p) ⟶ X)
        (x : ∀ p, Fin (N + 1) → NumberField.RingOfIntegers (fld p))
        (_ : ∀ p, x p ≠ 0) (_ : ∀ p, x p 0 ≠ 0)
        (chartOf : ∀ p, ∀ (Q : Ideal (𝓞 (fld p))), Q.IsMaximal → Fin (N + 1))
        (y : ∀ p, haveI := hnf p; ∀ (Q : Ideal (𝓞 (fld p))) (hQ : Q.IsMaximal),
          Spec (CommRingCat.of (Localization.AtPrime Q)) ⟶
            (nonVanishing (presheafTensorPow M L) (s (chartOf p Q hQ))).toScheme)
        (_ : ∀ p, haveI := hnf p; ∀ (Q : Ideal (𝓞 (fld p))) (hQ : Q.IsMaximal),
          Spec.map (CommRingCat.ofHom
              (algebraMap (𝓞 (fld p)) (Localization.AtPrime Q))) ≫ xF p
            = y p Q hQ ≫ (nonVanishing (presheafTensorPow M L) (s (chartOf p Q hQ))).ι)
        (_ : ∀ p, ∀ (Q : Ideal (𝓞 (fld p))) (hQ : Q.IsMaximal),
          Ideal.map (algebraMap (𝓞 (fld p)) (Localization.AtPrime Q))
              (Ideal.span (Set.range (x p)))
            = Ideal.span {algebraMap (𝓞 (fld p)) (Localization.AtPrime Q)
                (x p (chartOf p Q hQ))})
        (_ : ∀ p, haveI := hnf p; ∀ (Q : Ideal (𝓞 (fld p))) (hQ : Q.IsMaximal),
          (Spec.preimage (y p Q hQ ≫ (nonVanishing (presheafTensorPow M L)
              (s (chartOf p Q hQ))).toScheme.toSpecΓ)).hom
            ((nonVanishing (presheafTensorPow M L) (s (chartOf p Q hQ))).topIso.inv.hom
              (globalRatio (presheafTensorPow M L)
                (isLocallyTrivial_presheafTensorPow hM L) (s 0) (s (chartOf p Q hQ))))
            * algebraMap (𝓞 (fld p)) (Localization.AtPrime Q) (x p (chartOf p Q hQ))
            = algebraMap (𝓞 (fld p)) (Localization.AtPrime Q) (x p 0))
        (archChart : ∀ p, haveI := hnf p; InfinitePlace (fld p) → Fin (N + 1))
        (ρ : ∀ p, haveI := hnf p; ∀ v : InfinitePlace (fld p),
          CommRingCat.of (HomogeneousLocalization.Away
            (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ)
            (MvPolynomial.X (archChart p v))) ⟶ CommRingCat.of ℂ)
        (_ : ∀ p, haveI := hnf p; ∀ v : InfinitePlace (fld p),
          archPoint (xF p) v ≫ globalToProj (presheafTensorPow M L)
              (isLocallyTrivial_presheafTensorPow hM L) φ s hcov
            = Spec.map (ρ p v) ≫ chartA N ℤ (archChart p v))
        (_ : ∀ p, haveI := hnf p; ∀ (v : InfinitePlace (fld p)) (k : Fin (N + 1)),
          archRingHom (fld p) v (x p k)
            = (ρ p v).hom (projCoord N ℤ (archChart p v) k)
              * archRingHom (fld p) v (x p (archChart p v)))
        (_ : ∀ p, haveI := hnf p; ∀ v : InfinitePlace (fld p), x p (archChart p v) ≠ 0)
        (v₀ : ∀ p, haveI := hnf p; InfinitePlace (fld p))
        (idx : Fin (N + 1)) (_ : ∀ p, x p idx ≠ 0)
        (_ : ∀ p, haveI := hnf p; ∀ k,
          archRingHom (fld p) (v₀ p) (x p k) = ((x p k : fld p) : ℂ))
        (_ : Function.Injective (fun p => haveI := hnf p;
          archPoint (xF p) (v₀ p) ≫ globalToProj (presheafTensorPow M L)
            (isLocallyTrivial_presheafTensorPow hM L) φ s hcov))
        (C : ℝ),
        {p : P | haveI := hnf p; htArith (fld p) D (xF p) ≤ C}.Finite := by
  obtain ⟨ιU, U, hcovU, hU, htriv⟩ := exists_affine_trivializing_cover M hM
  exact northcott_of_isAmple_geom M hM U hcovU hU (fun j => (htriv j).some) hample x₀ φ d

/-! ## ★出典の紐付け(`.src`) -/

def exists_affine_trivializing_cover.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(可逆層は自明化するアフィン開被覆を持つ)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple_final.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(IsAmple から Northcott まで一本——幾何の仮定も外れた)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple_final.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_of_isAmple_geom(点の側は全部幾何の言葉、§9-937)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_isAmple_geom") 5,
    .citation "[mathlib]" "Scheme.isBasis_affineOpens(アフィン開は位相の基底)"
      (.inMathlib "AlgebraicGeometry.Scheme.isBasis_affineOpens") 2,
    .implicitStep
      ("★★★★★測定(2026-08-29): §9-927〜937 の鎖が受けていた**幾何の側の仮定**" ++
       "(M を自明化するアフィン開被覆 U・hcovU・hU・eU)は " ++
       "IsLocallyTrivial X M からそのまま出る。" ++
       "篩は前合成で閉じているので、S の元に含まれるアフィン開も S の元だからである。" ++
       "★有限性は要らない——northcott_of_isAmple_geom の添字型に有限性の要求がない") 5,
    .implicitStep
      ("★★★★★★これで Proposition 1.4, (iv) の仮定は" ++
       "可逆層(hM)・豊富性(hample)・点(x₀)・構造射(φ)・次数の上限(d)だけになった") 6 ]

end ABC3.Found.GenEll
