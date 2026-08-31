/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LocalChartExists
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★`Proposition 1.4, (iv)` —— 束ねた形（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 点ごとの仮定を 2 本に束ねる

`§9-938` の `northcott_of_isAmple_final` は点ごとのデータを
`chartOf`・`y`・`hy`・`hspan`・`hw`・`archChart`・`ρ`・`hfac`・`hcv`・`hiv` と
**10 個の束縛**で受けていた。★★本ファイルはそれを**存在命題 2 本**に束ねる:

    有限素点: `∀ p Q, ∃ i y, (局所化した点が X_{s_i} を通る) ∧ (𝔞_Q = (x_i)) ∧ ((s_0/s_i)(y)·x_i = x_0)`
    無限素点: `∀ p v, ∃ i ρ, (複素点が D₊(x_i) を通る) ∧ (σ_v(x_k) = ρ(x_k/x_i)·σ_v(x_i)) ∧ (x_i ≠ 0)`

★★★★これで `Proposition 1.4, (iv)` は**読める形**になった:

| 仮定 | 意味 |
|---|---|
| `hM`・`hample`・`x₀`・`φ`・`d` | 可逆層・豊富・点・構造射・次数の上限 |
| `x`・`hx`・`h0`・`hidx` | 点の同次座標（`0` でない） |
| 有限素点の整合 | ★`§9-939` でチャートの存在は定理、残るのは座標との整合 |
| 無限素点の整合 | ★同上 |
| `v₀`・`hemb`・`hpt` | 点の族が相異なる複素点を与えること（`§9-936`） |

## ★残っている段（明示）

★★★整合（`𝔞_Q = (x_i)`、`(s_0/s_i)(y)·x_i = x_0`、`σ_v(x_k) = ρ(x_k/x_i)·σ_v(x_i)`）は
**座標 `x` を点から構成する**段と一体である——`x` を点から作れば構成の定義から出る。
★`§9-939` が「チャートは常に取れる」ことを与えているので、
残るのは「そのチャートに合わせて座標を作る」ことだけである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open MvPolynomial HomogeneousLocalization NumberField

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★束ねた到達形 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★**`Proposition 1.4, (iv)`（束ねた形）**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★点ごとのデータ（`§9-938` では 10 個の束縛）を**存在命題 2 本**に束ねた形である。
★★★★`§9-939` により「チャートが取れる」部分は定理なので、
実質的に残っているのは**座標とチャートの整合**だけである。 -/
theorem northcott_of_isAmple_packed (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
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
        -- ★有限素点の整合
        (_ : ∀ p, haveI := hnf p; ∀ (Q : Ideal (𝓞 (fld p))) (_ : Q.IsMaximal),
          ∃ (i : Fin (N + 1)) (y : Spec (CommRingCat.of (Localization.AtPrime Q)) ⟶
            (nonVanishing (presheafTensorPow M L) (s i)).toScheme),
            Spec.map (CommRingCat.ofHom
                (algebraMap (𝓞 (fld p)) (Localization.AtPrime Q))) ≫ xF p
              = y ≫ (nonVanishing (presheafTensorPow M L) (s i)).ι ∧
            Ideal.map (algebraMap (𝓞 (fld p)) (Localization.AtPrime Q))
                (Ideal.span (Set.range (x p)))
              = Ideal.span {algebraMap (𝓞 (fld p)) (Localization.AtPrime Q) (x p i)} ∧
            (Spec.preimage (y ≫ (nonVanishing (presheafTensorPow M L)
                (s i)).toScheme.toSpecΓ)).hom
              ((nonVanishing (presheafTensorPow M L) (s i)).topIso.inv.hom
                (globalRatio (presheafTensorPow M L)
                  (isLocallyTrivial_presheafTensorPow hM L) (s 0) (s i)))
              * algebraMap (𝓞 (fld p)) (Localization.AtPrime Q) (x p i)
              = algebraMap (𝓞 (fld p)) (Localization.AtPrime Q) (x p 0))
        -- ★無限素点の整合
        (_ : ∀ p, haveI := hnf p; ∀ v : InfinitePlace (fld p),
          ∃ (i : Fin (N + 1)) (ρ : CommRingCat.of (HomogeneousLocalization.Away
              (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X i))
            ⟶ CommRingCat.of ℂ),
            archPoint (xF p) v ≫ globalToProj (presheafTensorPow M L)
                (isLocallyTrivial_presheafTensorPow hM L) φ s hcov
              = Spec.map ρ ≫ chartA N ℤ i ∧
            (∀ k, archRingHom (fld p) v (x p k)
              = ρ.hom (projCoord N ℤ i k) * archRingHom (fld p) v (x p i)) ∧
            x p i ≠ 0)
        (v₀ : ∀ p, haveI := hnf p; InfinitePlace (fld p))
        (idx : Fin (N + 1)) (_ : ∀ p, x p idx ≠ 0)
        (_ : ∀ p, haveI := hnf p; ∀ k,
          archRingHom (fld p) (v₀ p) (x p k) = ((x p k : fld p) : ℂ))
        (_ : Function.Injective (fun p => haveI := hnf p;
          archPoint (xF p) (v₀ p) ≫ globalToProj (presheafTensorPow M L)
            (isLocallyTrivial_presheafTensorPow hM L) φ s hcov))
        (C : ℝ),
        {p : P | haveI := hnf p; htArith (fld p) D (xF p) ≤ C}.Finite := by
  obtain ⟨L, hL, N, s, hcov, haffs, key⟩ := northcott_of_isAmple_final M hM hample x₀ φ d
  refine ⟨L, hL, N, s, hcov, haffs, ?_⟩
  intro D hDdiv hDgreen P fld hnf hdeg xF x hx h0 hfin harc v₀ idx hidx hemb hpt C
  choose chartOf y hy hspan hw using hfin
  choose archChart ρ hfac hcv hiv using harc
  exact key D hDdiv hDgreen fld hnf hdeg xF x hx h0 chartOf y hy hspan hw
    archChart ρ hfac hcv hiv v₀ idx hidx hemb hpt C

/-! ## ★出典の紐付け(`.src`) -/

def northcott_of_isAmple_packed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(点ごとの仮定を存在命題 2 本に束ねた形)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple_packed.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_of_isAmple_final(幾何の仮定も外れた形、§9-938)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_isAmple_final") 5,
    .citation "[ABC3]" "exists_localChart_family(チャートは常に取れる、§9-939)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_localChart_family") 3,
    .implicitStep
      ("★★★§9-939 により「チャートが取れる」部分は定理なので、" ++
       "束ねた 2 本の存在命題に実質的に残っているのは**座標とチャートの整合**だけである" ++
       "——𝔞_Q = (x_i)、(s_0/s_i)(y)·x_i = x_0、σ_v(x_k) = ρ(x_k/x_i)·σ_v(x_i)") 5,
    .implicitStep
      ("★★整合は「座標 x を点から構成する」段と一体である" ++
       "——x を点から作れば構成の定義から出る。" ++
       "★これが Proposition 1.4, (iv) に残った最後の段である") 4 ]

end ABC3.Found.GenEll
