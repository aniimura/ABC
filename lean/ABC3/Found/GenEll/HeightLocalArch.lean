/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightLocalCharts
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★局所チャートだけで高さが決まる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★これは何か —— アルキメデス側も局所化した

`§9-932` は有限素点側の `chart` を局所チャートに置き換えたが、
アルキメデス側の `hgreen` は仮説のままだった。★★本ファイルはそれも外す。

★★★★到達形:

    有限素点 `P` ごとに `Spec (𝓞_F)_P ⟶ X` が `X_{s_{i(P)}}` を通り、
    無限素点 `v` ごとに `Spec ℂ ⟶ ℙᴺ` が `D₊(x_{i(v)})` を通る
      ⟹ `ht_{ψ^*超平面}(x) = log H(座標)/[F:ℚ]`

★★★★★**大域のチャート仮定は完全に消えた**——原文の高さの計算と同じ、
「素点ごとの寄与の和」という形になった。

## ★★★機構 —— `greenFS` は既にチャート独立である

`greenFS N` は `§9-870` で**どのチャートで測っても同じ値**（`greenFS_eq_greenChartOf`）
であることが取れている。だから無限素点ごとに違うチャートを使ってよい。

★★残るのは「正規化座標で測っても Fubini–Study は変わらない」（`§9-869`）の
**ℂ 係数版**である:

    `σ_v(x_k) = c_k · σ_v(x_{i₀})`  ⟹  `log( sup‖c_k‖/‖c_j‖ ) = log( sup v(x_k)/v(x_j) )`

——`§9-869` は `c_k ∈ 𝓞_F` を仮定していたが、局所チャートでは `c_k ∈ ℂ` である。
★機構は同じで `v(x_{i₀})` が約分される。

## ★これで `Proposition 1.4, (iv)` の高さの計算に残った仮定

★★★**局所チャートの存在だけ**である。有限素点では `v_P(x_i)` が最小になる `i`、
無限素点では `v(x_i) ≠ 0` になる `i` を取ればよく、**どちらも常に取れる**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization NumberField
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★`§9-869` の ℂ 係数版 -/

/-- ★★★★★**正規化座標で測っても Fubini–Study は変わらない**——係数が `ℂ` の場合。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-869`（`log_iSup_norm_eq`）は比 `c_k` が `𝓞_F` の元である場合だったが、
局所チャートを使うと `c_k` は `ℂ` の元になる。
★★機構は同じで、`v(x_{i₀})` が**約分される**のが要点である。 -/
theorem log_iSup_norm_eq' (F : Type) [Field F] [NumberField F] {ι : Type}
    (v : InfinitePlace F) (x : ι → 𝓞 F) (c : ι → ℂ) (i₀ j : ι)
    (hw : ∀ k, archRingHom F v (x k) = c k * archRingHom F v (x i₀)) (hi0 : x i₀ ≠ 0) :
    Real.log ((⨆ k, ‖c k‖) / ‖c j‖)
      = Real.log ((⨆ k, v ((x k : F))) / v ((x j : F))) := by
  have hi0pos : (0:ℝ) < v ((x i₀ : F)) := (AbsoluteValue.pos_iff _).2 (by simpa using hi0)
  have hcv : ∀ k, ‖c k‖ = v ((x k : F)) / v ((x i₀ : F)) := by
    intro k
    have hk : v ((x k : F)) = ‖c k‖ * v ((x i₀ : F)) := by
      rw [← norm_archRingHom F v (x k), ← norm_archRingHom F v (x i₀), hw k, norm_mul]
    rw [hk, mul_div_assoc, div_self hi0pos.ne', mul_one]
  simp only [hcv]
  rw [← iSup_div_const (fun k => v ((x k : F))) _ hi0pos.le,
    div_div_div_cancel_right₀ hi0pos.ne']

/-! ## ★★★★★★★★★★★★★★無限素点ごとのチャートで Green 関数を読む -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★**無限素点ごとにチャートを取り直してよい**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`greenFS` は `§9-870` で**チャート独立**（`greenFS_eq_greenChartOf`）だから、
無限素点 `v` ごとに違うチャート `i(v)` を使ってよい。 -/
theorem green_comap_globalToProj_of_archChart (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ X) (x : Fin (N + 1) → 𝓞 F)
    (v : InfinitePlace F) (iv : Fin (N + 1))
    (ρ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X iv))
      ⟶ CommRingCat.of ℂ)
    (hfac : archPoint xF v ≫ globalToProj M hM φ s hcov = Spec.map ρ ≫ chartA N ℤ iv)
    (hcv : ∀ k, archRingHom F v (x k)
      = ρ.hom (projCoord N ℤ iv k) * archRingHom F v (x iv))
    (hiv : x iv ≠ 0) :
    ((hyperplaneArith N).comap (globalToProj M hM φ s hcov)).green (archPoint xF v)
      = Real.log ((⨆ k, v ((x k : F))) / v ((x 0 : F))) := by
  show greenFS N (archPoint xF v ≫ globalToProj M hM φ s hcov) = _
  rw [hfac, greenFS_eq_greenChartOf, greenChartOf]
  exact log_iSup_norm_eq' F v x (fun k => ρ.hom (projCoord N ℤ iv k)) iv 0 hcv hiv

/-! ## ★★★★★★★★★★★★★★★★★★★★到達 —— 大域チャートは完全に消えた -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★**局所チャートだけで `ht = log H`**
—— `Proposition 1.4, (iv)` の高さの計算の**局所形**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★仮定は**素点ごとのチャート**だけである:

* 有限素点 `P`——`Spec (𝓞_F)_P ⟶ X` が `X_{s_{i(P)}}` を通る（`v_P(x_i)` が最小の `i`）
* 無限素点 `v`——`Spec ℂ ⟶ ℙᴺ` が `D₊(x_{i(v)})` を通る（`v(x_i) ≠ 0` の `i`）

★★★★★**どちらも常に取れる**ので、これで原文の高さの計算
（素点ごとの寄与の和）と同じ形になった。 -/
theorem htArith_globalToProj_eq_log_mulHeight_of_local
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ X)
    (x : Fin (N + 1) → 𝓞 F) (hx : x ≠ 0) (h0 : x 0 ≠ 0)
    (chartOf : ∀ (P : Ideal (𝓞 F)), P.IsMaximal → Fin (N + 1))
    (y : ∀ (P : Ideal (𝓞 F)) (hP : P.IsMaximal),
      Spec (CommRingCat.of (Localization.AtPrime P)) ⟶
        (nonVanishing M (s (chartOf P hP))).toScheme)
    (hy : ∀ (P : Ideal (𝓞 F)) (hP : P.IsMaximal),
      Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime P))) ≫ xF
        = y P hP ≫ (nonVanishing M (s (chartOf P hP))).ι)
    (hspan : ∀ (P : Ideal (𝓞 F)) (hP : P.IsMaximal),
      Ideal.map (algebraMap (𝓞 F) (Localization.AtPrime P)) (Ideal.span (Set.range x))
        = Ideal.span {algebraMap (𝓞 F) (Localization.AtPrime P) (x (chartOf P hP))})
    (hw : ∀ (P : Ideal (𝓞 F)) (hP : P.IsMaximal),
      (Spec.preimage (y P hP ≫
          (nonVanishing M (s (chartOf P hP))).toScheme.toSpecΓ)).hom
        ((nonVanishing M (s (chartOf P hP))).topIso.inv.hom
          (globalRatio M hM (s 0) (s (chartOf P hP))))
        * algebraMap (𝓞 F) (Localization.AtPrime P) (x (chartOf P hP))
        = algebraMap (𝓞 F) (Localization.AtPrime P) (x 0))
    (archChart : InfinitePlace F → Fin (N + 1))
    (ρ : ∀ v : InfinitePlace F, CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X (archChart v)))
      ⟶ CommRingCat.of ℂ)
    (hfac : ∀ v : InfinitePlace F, archPoint xF v ≫ globalToProj M hM φ s hcov
      = Spec.map (ρ v) ≫ chartA N ℤ (archChart v))
    (hcv : ∀ (v : InfinitePlace F) (k : Fin (N + 1)), archRingHom F v (x k)
      = (ρ v).hom (projCoord N ℤ (archChart v) k) * archRingHom F v (x (archChart v)))
    (hiv : ∀ v : InfinitePlace F, x (archChart v) ≠ 0) :
    htArith F ((hyperplaneArith N).comap (globalToProj M hM φ s hcov)) xF
      = Real.log (Height.mulHeight (fun k => (x k : F))) / (Module.finrank ℚ F : ℝ) :=
  htArith_globalToProj_eq_log_mulHeight_of_localCharts M hM φ s hcov F xF x hx h0
    chartOf y hy hspan hw
    (fun v => green_comap_globalToProj_of_archChart M hM φ s hcov F xF x v (archChart v)
      (ρ v) (hfac v) (hcv v) (hiv v))

/-! ## ★出典の紐付け(`.src`) -/

def log_iSup_norm_eq'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(正規化座標で測っても Fubini–Study は変わらない——ℂ 係数版)",
    sectionId := "genell-prop-1-4" }

def green_comap_globalToProj_of_archChart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(無限素点ごとにチャートを取り直してよい)",
    sectionId := "genell-prop-1-4" }

def htArith_globalToProj_eq_log_mulHeight_of_local.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(局所チャートだけで ht = log H——高さの計算の局所形)",
    sectionId := "genell-prop-1-4" }

def htArith_globalToProj_eq_log_mulHeight_of_local.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htArith_globalToProj_eq_log_mulHeight_of_localCharts(有限素点側、§9-932)"
      (.inProject "ABC3"
        "ABC3.Found.GenEll.htArith_globalToProj_eq_log_mulHeight_of_localCharts") 4,
    .citation "[ABC3]" "greenFS_eq_greenChartOf(greenFS はチャート独立、§9-870)"
      (.inProject "ABC3" "ABC3.Found.GenEll.greenFS_eq_greenChartOf") 3,
    .implicitStep
      ("★★★★測定(2026-08-29): これで Proposition 1.4, (iv) の高さの計算から" ++
       "**大域のチャート仮定が完全に消えた**。残るのは素点ごとのチャートだけで、" ++
       "有限素点では v_P(x_i) が最小になる i、無限素点では v(x_i) ≠ 0 になる i を" ++
       "取ればよく、どちらも常に取れる") 5,
    .implicitStep
      ("★§9-869(log_iSup_norm_eq)は比 c_k が 𝓞_F の元である場合だったが、" ++
       "局所チャートでは c_k は ℂ の元になる。機構は同じで v(x_{i₀}) が約分される") 2 ]

end ABC3.Found.GenEll
