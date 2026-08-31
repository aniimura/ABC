/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottLocalAmple
import ABC3.Found.GenEll.InjOfComplexPoints
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★`Proposition 1.4, (iv)` —— 幾何の言葉だけで（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 点の側の仮定が全部翻訳された

`§9-935` までで点の側に残っていたのは

* 大域チャート（`chart`）——`§9-930`〜`§9-935` で**素点ごとのチャート**に置き換わった
* `hinj`（正規化座標が点を分けること）——`§9-936` で**複素点の単射性**に翻訳された

★★★★本ファイルはその 2 つを繋いで、`Proposition 1.4, (iv)` を

    「点の族が相異なる複素点を与え、素点ごとにチャートが取れる ⟹ 高さで有界な点は有限」

という形にする。★★★★★**もはや点の側に算術的な仮定はない**。

## ★★★これで残っている仮定は 1 つだけ

★★★★★**幾何の側の有限アフィン自明化被覆**（`U`・`hcovU`・`hU`・`eU`）である。
`X` が `Spec A` 上分離的・準コンパクトで `M` が可逆なら出るはずのもので、
**Arakelov 理論ではなくスキーム論の一般論**である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open MvPolynomial HomogeneousLocalization NumberField

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★到達 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★★**切断の零点因子についての Northcott**
——点の側の仮定を全部幾何の言葉にした形。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★点の側の仮定は
* **素点ごとのチャート**（有限素点・無限素点それぞれ、★常に取れる）
* **点の族が相異なる複素点を与えること**（`hpt`）

だけである。★★★★`§9-928` の `chart`（イデアル類の障害）も
`§9-935` の `hinj`（正規化座標の単射性）も消えた。 -/
theorem northcott_of_divisorOfSection_geom (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ} (d : ℕ)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (haff : ∀ i, IsAffineOpen (nonVanishing M (s i)))
    {ιU : Type} (U : ιU → X.Opens) (hcovU : (⨆ j, U j) = ⊤)
    (hU : ∀ j, IsAffineOpen (U j))
    (eU : ∀ j, (restrictPresheafFunctor X (U j)).obj M ≅ 𝟙_ (PresheafModulesOn X (U j)))
    (D : ArithCartier X)
    (hDdiv : D.divisor = divisorOfSection M (s 0))
    (hDgreen : D.green = fun p => greenFS N (p ≫ globalToProj M hM φ s hcov))
    {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (xF : ∀ p, haveI := hnf p; specRingOfIntegers (fld p) ⟶ X)
    (x : ∀ p, Fin (N + 1) → NumberField.RingOfIntegers (fld p))
    (hx : ∀ p, x p ≠ 0) (h0 : ∀ p, x p 0 ≠ 0)
    (chartOf : ∀ p, ∀ (Q : Ideal (𝓞 (fld p))), Q.IsMaximal → Fin (N + 1))
    (y : ∀ p, haveI := hnf p; ∀ (Q : Ideal (𝓞 (fld p))) (hQ : Q.IsMaximal),
      Spec (CommRingCat.of (Localization.AtPrime Q)) ⟶
        (nonVanishing M (s (chartOf p Q hQ))).toScheme)
    (hy : ∀ p, haveI := hnf p; ∀ (Q : Ideal (𝓞 (fld p))) (hQ : Q.IsMaximal),
      Spec.map (CommRingCat.ofHom (algebraMap (𝓞 (fld p)) (Localization.AtPrime Q))) ≫ xF p
        = y p Q hQ ≫ (nonVanishing M (s (chartOf p Q hQ))).ι)
    (hspan : ∀ p, ∀ (Q : Ideal (𝓞 (fld p))) (hQ : Q.IsMaximal),
      Ideal.map (algebraMap (𝓞 (fld p)) (Localization.AtPrime Q))
          (Ideal.span (Set.range (x p)))
        = Ideal.span {algebraMap (𝓞 (fld p)) (Localization.AtPrime Q)
            (x p (chartOf p Q hQ))})
    (hw : ∀ p, haveI := hnf p; ∀ (Q : Ideal (𝓞 (fld p))) (hQ : Q.IsMaximal),
      (Spec.preimage (y p Q hQ ≫
          (nonVanishing M (s (chartOf p Q hQ))).toScheme.toSpecΓ)).hom
        ((nonVanishing M (s (chartOf p Q hQ))).topIso.inv.hom
          (globalRatio M hM (s 0) (s (chartOf p Q hQ))))
        * algebraMap (𝓞 (fld p)) (Localization.AtPrime Q) (x p (chartOf p Q hQ))
        = algebraMap (𝓞 (fld p)) (Localization.AtPrime Q) (x p 0))
    (archChart : ∀ p, haveI := hnf p; InfinitePlace (fld p) → Fin (N + 1))
    (ρ : ∀ p, haveI := hnf p; ∀ v : InfinitePlace (fld p),
      CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X (archChart p v)))
      ⟶ CommRingCat.of ℂ)
    (hfac : ∀ p, haveI := hnf p; ∀ v : InfinitePlace (fld p),
      archPoint (xF p) v ≫ globalToProj M hM φ s hcov
        = Spec.map (ρ p v) ≫ chartA N ℤ (archChart p v))
    (hcv : ∀ p, haveI := hnf p; ∀ (v : InfinitePlace (fld p)) (k : Fin (N + 1)),
      archRingHom (fld p) v (x p k)
        = (ρ p v).hom (projCoord N ℤ (archChart p v) k)
          * archRingHom (fld p) v (x p (archChart p v)))
    (hiv : ∀ p, haveI := hnf p; ∀ v : InfinitePlace (fld p), x p (archChart p v) ≠ 0)
    (v₀ : ∀ p, haveI := hnf p; InfinitePlace (fld p))
    (idx : Fin (N + 1)) (hidx : ∀ p, x p idx ≠ 0)
    (hemb : ∀ p, haveI := hnf p; ∀ k,
      archRingHom (fld p) (v₀ p) (x p k) = ((x p k : fld p) : ℂ))
    (hpt : Function.Injective (fun p => haveI := hnf p;
      archPoint (xF p) (v₀ p) ≫ globalToProj M hM φ s hcov))
    (C : ℝ) :
    {p : P | haveI := hnf p; htArith (fld p) D (xF p) ≤ C}.Finite := by
  have hinj := hinj_of_injective_archPoint M hM φ s hcov fld hnf xF x v₀
    (fun p => archChart p (v₀ p)) (fun p => ρ p (v₀ p)) (fun p => hfac p (v₀ p))
    (fun p => hcv p (v₀ p)) (fun p => hiv p (v₀ p)) idx hidx hemb hpt
  exact northcott_of_divisorOfSection_local M hM d φ s hcov haff U hcovU hU eU D
    hDdiv hDgreen fld hnf hdeg xF x hx h0 chartOf y hy hspan hw
    archChart ρ hfac hcv hiv idx hinj C

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★★**`IsAmple` から Northcott まで一本**
——点の側の仮定を全部幾何の言葉にした形。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★★★残るのは**幾何の側の有限アフィン自明化被覆**（`U`・`hcovU`・`hU`・`eU`）だけである
——`X` が `Spec A` 上分離的・準コンパクトで `M` が可逆なら出るはずのもので、
**Arakelov 理論ではなくスキーム論の一般論**である。 -/
theorem northcott_of_isAmple_geom (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {ιU : Type} (U : ιU → X.Opens) (hcovU : (⨆ j, U j) = ⊤)
    (hU : ∀ j, IsAffineOpen (U j))
    (eU : ∀ j, (restrictPresheafFunctor X (U j)).obj M ≅ 𝟙_ (PresheafModulesOn X (U j)))
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
  obtain ⟨L, hL, N, s, hcov, haffs⟩ := exists_fin_cover_of_isAmple M hM hample x₀
  refine ⟨L, hL, N, s, hcov, haffs, ?_⟩
  intro D hDdiv hDgreen P fld hnf hdeg xF x hx h0 chartOf y hy hspan hw
    archChart ρ hfac hcv hiv v₀ idx hidx hemb hpt C
  exact northcott_of_divisorOfSection_geom (presheafTensorPow M L)
    (isLocallyTrivial_presheafTensorPow hM L) d φ s hcov haffs
    U hcovU hU (fun j => tensorPowTriv (eU j) L) D hDdiv hDgreen
    fld hnf hdeg xF x hx h0 chartOf y hy hspan hw archChart ρ hfac hcv hiv
    v₀ idx hidx hemb hpt C

/-! ## ★出典の紐付け(`.src`) -/

def northcott_of_divisorOfSection_geom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(切断の零点因子についての Northcott——点の側は全部幾何の言葉)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple_geom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(IsAmple から Northcott まで一本——点の側は全部幾何の言葉)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple_geom.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_of_divisorOfSection_local(局所チャート版、§9-935)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_divisorOfSection_local") 4,
    .citation "[ABC3]" "hinj_of_injective_archPoint(hinj は複素点の単射性、§9-936)"
      (.inProject "ABC3" "ABC3.Found.GenEll.hinj_of_injective_archPoint") 4,
    .citation "[ABC3]" "exists_fin_cover_of_isAmple(IsAmple から被覆、§9-914・915)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_fin_cover_of_isAmple") 3,
    .implicitStep
      ("★★★★★測定(2026-08-29): Proposition 1.4, (iv) の**点の側の仮定は全部翻訳された**。" ++
       "大域チャート(chart)は素点ごとのチャートに(§9-930〜935)、" ++
       "hinj は複素点の単射性に(§9-936)。" ++
       "★もはや点の側に算術的な仮定はない") 6,
    .implicitStep
      ("★★★残っているのは幾何の側の有限アフィン自明化被覆(U・hcovU・hU・eU)だけである" ++
       "——X が Spec A 上分離的・準コンパクトで M が可逆なら出るはずのもので、" ++
       "**Arakelov 理論ではなくスキーム論の一般論**である") 5 ]

end ABC3.Found.GenEll
