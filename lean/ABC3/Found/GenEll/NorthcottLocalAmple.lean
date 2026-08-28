/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottLocal
import ABC3.Found.GenEll.NorthcottOfIsAmple
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★`IsAmple` から Northcott まで（局所チャート版、`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★★これは何か —— `§9-928` の `chart` を外した形

`§9-927`・`§9-928` は点の側に

    `chart : P → Fin (N+1)`、`xF p : Spec 𝓞_{F_p} ⟶ X_{s_{chart p}}`

——**点が大域的にどれかのチャートを通ること**——を要求していた。
★`§9-928` で測ったとおり、これは**イデアル類が自明でないと成り立たない**。

★★★★本ファイルはそれを `§9-934`（局所チャート版の Northcott）で置き換える:

* 点は `xF p : Spec 𝓞_{F_p} ⟶ X`（★チャートを通らなくてよい）
* 代わりに**素点ごとのチャート**を与える（有限素点・無限素点それぞれ）

## ★★★これで `Proposition 1.4, (iv)` に残った仮定

| 側 | 仮定 |
|---|---|
| 幾何 | 有限アフィン自明化被覆（`U`・`hcovU`・`hU`・`eU`）★スキーム論の一般論 |
| 点 | 素点ごとのチャート ★**常に取れる** |
| 点 | `hinj`（座標が点を分けること）★`ψ` が埋め込みであることの点版 |

★★★★★イデアル類の障害（`§9-928` の測定）は**完全に消えた**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open MvPolynomial HomogeneousLocalization NumberField

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★切断の零点因子は超平面の引き戻しである -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★**`div(s₀)` ＋ Fubini–Study の引き戻し ＝ `ψ^*(超平面の算術因子)`**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-926`（段 E0）と `ArithCartier.ext'` を繋いだだけである。 -/
theorem eq_hyperplaneArith_comap_of_divisorOfSection (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (haff : ∀ i, IsAffineOpen (nonVanishing M (s i)))
    {ι : Type} (U : ι → X.Opens) (hcovU : (⨆ j, U j) = ⊤)
    (hU : ∀ j, IsAffineOpen (U j))
    (eU : ∀ j, (restrictPresheafFunctor X (U j)).obj M ≅ 𝟙_ (PresheafModulesOn X (U j)))
    (D : ArithCartier X)
    (hDdiv : D.divisor = divisorOfSection M (s 0))
    (hDgreen : D.green = fun p => greenFS N (p ≫ globalToProj M hM φ s hcov)) :
    D = (hyperplaneArith N).comap (globalToProj M hM φ s hcov) := by
  refine ArithCartier.ext' ?_ hDgreen
  rw [hDdiv, ArithCartier.comap_divisor]
  exact divisorOfSection_eq_comap_globalToProj M hM φ s hcov haff U hcovU hU eU

/-! ## ★★★★★★★★★★★★★★★★★★★★★★局所チャート版の到達 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★**切断の零点因子についての Northcott**
（局所チャート版）。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★`§9-927` は点が**大域的に**チャートを通ることを要求していたが、
本定理では点は `Spec 𝓞_F ⟶ X` でよく、チャートは**素点ごと**である。 -/
theorem northcott_of_divisorOfSection_local (M : X.PresheafOfModules)
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
    (idx : Fin (N + 1))
    (hinj : Function.Injective (fun (p : P) (k : Fin (N + 1)) =>
      ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)))
    (C : ℝ) :
    {p : P | haveI := hnf p; htArith (fld p) D (xF p) ≤ C}.Finite := by
  rw [eq_hyperplaneArith_comap_of_divisorOfSection M hM φ s hcov haff U hcovU hU eU D
    hDdiv hDgreen]
  exact northcott_globalToProj_of_local M hM d φ s hcov fld hnf hdeg xF x hx h0
    chartOf y hy hspan hw archChart ρ hfac hcv hiv idx hinj C

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★**`IsAmple` から Northcott まで一本**
（局所チャート版）。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★★`§9-928` の `chart` 仮定（点が大域的にチャートを通ること）が
**素点ごとのチャート**に置き換わった形である。 -/
theorem northcott_of_isAmple_local (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
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
        (idx : Fin (N + 1))
        (_ : Function.Injective (fun (p : P) (k : Fin (N + 1)) =>
          ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)))
        (C : ℝ),
        {p : P | haveI := hnf p; htArith (fld p) D (xF p) ≤ C}.Finite := by
  obtain ⟨L, hL, N, s, hcov, haffs⟩ := exists_fin_cover_of_isAmple M hM hample x₀
  refine ⟨L, hL, N, s, hcov, haffs, ?_⟩
  intro D hDdiv hDgreen P fld hnf hdeg xF x hx h0 chartOf y hy hspan hw
    archChart ρ hfac hcv hiv idx hinj C
  exact northcott_of_divisorOfSection_local (presheafTensorPow M L)
    (isLocallyTrivial_presheafTensorPow hM L) d φ s hcov haffs
    U hcovU hU (fun j => tensorPowTriv (eU j) L) D hDdiv hDgreen
    fld hnf hdeg xF x hx h0 chartOf y hy hspan hw archChart ρ hfac hcv hiv idx hinj C

/-! ## ★出典の紐付け(`.src`) -/

def eq_hyperplaneArith_comap_of_divisorOfSection.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(div(s₀) ＋ Fubini–Study の引き戻し ＝ ψ^*(超平面の算術因子))",
    sectionId := "genell-prop-1-4" }

def northcott_of_divisorOfSection_local.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(切断の零点因子についての Northcott——局所チャート版)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple_local.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(IsAmple から Northcott まで一本——局所チャート版)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple_local.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_globalToProj_of_local(局所チャート版の Northcott、§9-934)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_globalToProj_of_local") 4,
    .citation "[ABC3]" "divisorOfSection_eq_comap_globalToProj(段 E0、§9-926)"
      (.inProject "ABC3" "ABC3.Found.GenEll.divisorOfSection_eq_comap_globalToProj") 4,
    .citation "[ABC3]" "exists_fin_cover_of_isAmple(IsAmple から被覆、§9-914・915)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_fin_cover_of_isAmple") 3,
    .implicitStep
      ("★★★★★測定(2026-08-29): §9-928 の chart 仮定(点が大域的にチャートを通ること" ++
       "——イデアル類が自明でないと成り立たない)は**素点ごとのチャート**に置き換わった。" ++
       "有限素点では v_P(x_i) が最小になる i、無限素点では v(x_i) ≠ 0 になる i を" ++
       "取ればよく、どちらも常に取れる。★イデアル類の障害は完全に消えた") 6,
    .implicitStep
      ("★★残るのは (1) 幾何の側の有限アフィン自明化被覆(U・hcovU・hU・eU" ++
       "——X が Spec A 上分離的・準コンパクトで M が可逆なら出るはずのもの)と " ++
       "(2) hinj(座標が点を分けること——ψ が埋め込み(§9-920)であることの点版)である") 4 ]

end ABC3.Found.GenEll
