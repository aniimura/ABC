/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightLocalArch
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★局所チャートだけで Northcott（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★これは何か —— Northcott も局所化した

`§9-933` で「局所チャートだけで `ht = log H/[F:ℚ]`」が取れた。
★★本ファイルはそれを Northcott まで押し切る:

    `ht` が座標の素朴高さの形で書けて、座標が点を分ける
      ⟹ `{p | ht(p) ≤ C}` は有限

★★★★これで `§9-921` の `northcott_globalToProj` から
**大域のチャート仮定 `chart`・`xF : Spec 𝓞_F ⟶ X_{s_i}` が消えた**
——点は `Spec 𝓞_F ⟶ X` でよく、チャートは**素点ごとに**取ればよい。

## ★★★機構 —— `northcott_of_projModel` に直接渡す

`§9-878` の `northcott_hyperplane` はチャート分解を経由していたが、
`§9-933` は `ht` の値そのものを与えるので、`§9-877`（`northcott_of_projModel`）に
**直接**渡せる。★次数の扱いだけ注意が要る:

* `H(x) = exp([F:ℚ]·ht)` であって `exp(ht)` ではない
* `[F:ℚ] ≤ d` かつ `ht ≥ 0`（`H ≥ 1` だから）なので `[F:ℚ]·ht ≤ d·ht`
* よって `H(x) ≤ exp(d·ht)` となり `northcott_of_projModel` の `hcmp` が埋まる

## ★残っている仮定（明示）

★★★点の側に残るのは **`hinj`（座標が点を分けること）**だけである。
これは `ψ` が埋め込み（`§9-920`）であることの点版であり、
点の族 `P` を座標で添字づけるための道具である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization NumberField
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★高さが素朴高さの形なら Northcott -/

/-- ★★★★★★★★★★★**高さが座標の素朴高さで書けるなら Northcott**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-878`（`northcott_hyperplane`）はチャート分解を経由していたが、
`ht` の値そのものが分かっていれば `§9-877` に**直接**渡せる。
★★次数の扱いが要点である——`H(x) = exp([F:ℚ]·ht)` かつ `[F:ℚ] ≤ d`、`ht ≥ 0` だから
`H(x) ≤ exp(d·ht)` である。 -/
theorem northcott_of_log_mulHeight (d : ℕ) {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    {ι : Type} [Finite ι] (x : ∀ p, ι → NumberField.RingOfIntegers (fld p))
    (ht : P → ℝ)
    (hht : ∀ p, haveI := hnf p; ht p
      = Real.log (Height.mulHeight (fun k => ((x p k : fld p))))
        / (Module.finrank ℚ (fld p) : ℝ))
    (idx : ι)
    (hinj : Function.Injective (fun (p : P) (k : ι) =>
      ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)))
    (C : ℝ) :
    {p : P | ht p ≤ C}.Finite := by
  refine Set.Finite.subset (northcott_of_projModel d (fun p => (d : ℝ) * ht p) fld hnf hdeg
    (fun p k => ((x p k : fld p))) idx 0 (fun p => ?_) hinj ((d : ℝ) * C)) (fun p hp => ?_)
  · haveI := hnf p
    have hpos : (0:ℝ) < Height.mulHeight (fun k => ((x p k : fld p))) :=
      lt_of_lt_of_le zero_lt_one (Height.one_le_mulHeight _)
    have hfr : (0:ℝ) < (Module.finrank ℚ (fld p) : ℝ) := by exact_mod_cast Module.finrank_pos
    have hlognn : 0 ≤ Real.log (Height.mulHeight (fun k => ((x p k : fld p)))) :=
      Real.log_nonneg (Height.one_le_mulHeight _)
    have htnn : 0 ≤ ht p := by rw [hht p]; exact div_nonneg hlognn hfr.le
    have hlog : (Module.finrank ℚ (fld p) : ℝ) * ht p
        = Real.log (Height.mulHeight (fun k => ((x p k : fld p)))) := by
      rw [hht p, mul_div_cancel₀ _ hfr.ne']
    have hle : (Module.finrank ℚ (fld p) : ℝ) * ht p ≤ (d : ℝ) * ht p :=
      mul_le_mul_of_nonneg_right (by exact_mod_cast hdeg p) htnn
    rw [add_zero]
    calc Height.mulHeight (fun k => ((x p k : fld p)))
        = Real.exp (Real.log (Height.mulHeight (fun k => ((x p k : fld p))))) :=
          (Real.exp_log hpos).symm
      _ = Real.exp ((Module.finrank ℚ (fld p) : ℝ) * ht p) := by rw [hlog]
      _ ≤ Real.exp ((d : ℝ) * ht p) := Real.exp_le_exp.2 hle
  · have hp' : ht p ≤ C := hp
    exact mul_le_mul_of_nonneg_left hp' (Nat.cast_nonneg d)

/-! ## ★★★★★★★★★★★★★★★★★★★★★到達 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★**局所チャートだけで Northcott**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★`§9-921` の `northcott_globalToProj` は点が**大域的に**チャートを通ること
（`chart`・`xF : Spec 𝓞_F ⟶ X_{s_i}`）を要求していた。
★★★★本定理では点は `Spec 𝓞_F ⟶ X` でよく、チャートは**素点ごとに**取る
——`§9-928` で測ったイデアル類の障害はここで消える。

★残る点の側の仮定は `hinj`（座標が点を分けること）だけである。 -/
theorem northcott_globalToProj_of_local (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ} (d : ℕ)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
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
    {p : P | haveI := hnf p;
      htArith (fld p) ((hyperplaneArith N).comap (globalToProj M hM φ s hcov)) (xF p)
        ≤ C}.Finite := by
  refine northcott_of_log_mulHeight d fld hnf hdeg x _ (fun p => ?_) idx hinj C
  haveI := hnf p
  exact htArith_globalToProj_eq_log_mulHeight_of_local M hM φ s hcov (fld p) (xF p)
    (x p) (hx p) (h0 p) (chartOf p) (y p) (hy p) (hspan p) (hw p)
    (archChart p) (ρ p) (hfac p) (hcv p) (hiv p)

/-! ## ★出典の紐付け(`.src`) -/

def northcott_of_log_mulHeight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(高さが座標の素朴高さで書けるなら Northcott)",
    sectionId := "genell-prop-1-4" }

def northcott_globalToProj_of_local.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(局所チャートだけで Northcott——大域チャート仮定なし)",
    sectionId := "genell-prop-1-4" }

def northcott_globalToProj_of_local.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htArith_globalToProj_eq_log_mulHeight_of_local(局所形の高さ、§9-933)"
      (.inProject "ABC3"
        "ABC3.Found.GenEll.htArith_globalToProj_eq_log_mulHeight_of_local") 4,
    .citation "[ABC3]" "northcott_of_projModel(hcmp ＋ hinj ⟹ 有限、§9-877)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_projModel") 3,
    .implicitStep
      ("★★★★測定(2026-08-29): これで §9-921 の northcott_globalToProj から" ++
       "**大域のチャート仮定(chart・xF : Spec 𝓞_F ⟶ X_{s_i})が消えた**。" ++
       "点は Spec 𝓞_F ⟶ X でよく、チャートは素点ごとに取ればよい" ++
       "——§9-928 で測ったイデアル類の障害はここで消える") 5,
    .implicitStep
      ("★次数の扱いが要点である——H(x) = exp([F:ℚ]·ht) であって exp(ht) ではない。" ++
       "[F:ℚ] ≤ d かつ ht ≥ 0(H ≥ 1 だから)なので [F:ℚ]·ht ≤ d·ht となり、" ++
       "northcott_of_projModel の hcmp が埋まる") 3,
    .implicitStep
      ("★★点の側に残るのは hinj(座標が点を分けること)だけである" ++
       "——ψ が埋め込み(§9-920)であることの点版であり、" ++
       "点の族 P を座標で添字づけるための道具である") 3 ]

end ABC3.Found.GenEll
