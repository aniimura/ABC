/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.PullbackChartLocal
import ABC3.Found.GenEll.HeightArithLocalization
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★局所チャートだけで `ht = log H`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★これは何か —— 幾何と算術を局所で継ぐ

* `§9-930`（幾何）——**局所化した点がチャートを通れば引き戻しイデアルの局所化が読める**
* `§9-931`（算術）——**局所条件だけで `ht = log H`**

★★★★本ファイルはこの 2 本を継ぐ:

    各極大イデアル `P` について
      `Spec (𝓞_F)_P ⟶ Spec 𝓞_F ⟶ X` が `X_{s_{i(P)}}` を通り、
      `𝔞·(𝓞_F)_P = (x_{i(P)})`、`(s_0/s_{i(P)})(y_P)·x_{i(P)} = x_0`
    ⟹ `ht_{ψ^*超平面}(x) = log H(座標)/[F:ℚ]`

★★★★★**大域のチャート仮定は消えた**。残るのは局所の因子で、これは
`v_P(x_i)` が最小になる `i` を取れば**常に満たせる**。

## ★★★これまでとの差 —— 何が外れたか

| 段 | 仮定 |
|---|---|
| `§9-921` | `xF : Spec 𝓞_F ⟶ X_{s_i}`（★**大域**のチャート） |
| `§9-929` | ノルムの等式 `N(I)·N(𝔞) = N((x_0))` |
| `§9-931` | 各極大イデアルでの局所条件 |
| **本ファイル** | ★**局所チャート**（常に取れる） |

## ★残っている段（明示）

★★アルキメデス側の `hgreen` は仮説のままである。
大域チャートの場合（`§9-871`）は**同じチャートが全ての素点で使えた**ので自動だったが、
局所チャートの場合は**無限素点ごとにチャートを取り直す**必要がある
——`archPoint xF v : Spec ℂ ⟶ X` が `X_{s_{i(v)}}` を通ること。
★これは体の上の点なので有限素点の側より易しい（座標が 0 でない `i` を選べばよい）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization NumberField
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★★★局所チャートから局所条件を作る -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★**局所チャートから `§9-931` の局所条件が出る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-930` の `map_pullbackIdeal_globalToProj_of_localChart` が
`I·(𝓞_F)_P = ( (s_0/s_i)(y_P) )` を与えるので、
あとは `𝔞·(𝓞_F)_P = (x_i)` と `(s_0/s_i)(y_P)·x_i = x_0` を渡すだけである。 -/
theorem hloc_of_localCharts (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ X)
    (x : Fin (N + 1) → 𝓞 F)
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
        = algebraMap (𝓞 F) (Localization.AtPrime P) (x 0)) :
    ∀ (P : Ideal (𝓞 F)) (_ : P.IsMaximal),
      ∃ (i : Fin (N + 1)) (g : Localization.AtPrime P),
        Ideal.map (algebraMap (𝓞 F) (Localization.AtPrime P)) (Ideal.span (Set.range x))
            = Ideal.span {algebraMap (𝓞 F) (Localization.AtPrime P) (x i)} ∧
          Ideal.map (algebraMap (𝓞 F) (Localization.AtPrime P))
              (pullbackIdeal F ((hyperplaneArith N).comap
                (globalToProj M hM φ s hcov)).divisor xF) = Ideal.span {g} ∧
          g * algebraMap (𝓞 F) (Localization.AtPrime P) (x i)
            = algebraMap (𝓞 F) (Localization.AtPrime P) (x 0) := by
  intro P hP
  refine ⟨chartOf P hP, _, hspan P hP, ?_, hw P hP⟩
  rw [ArithCartier.comap_divisor]
  exact map_pullbackIdeal_globalToProj_of_localChart M hM φ s hcov (chartOf P hP) F
    (Localization.AtPrime P) xF (y P hP) (hy P hP)

/-! ## ★★★★★★★★★★★★★★★★★★到達 —— 大域チャートなしの高さ -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★**局所チャートだけで `ht_{ψ^*超平面} = log H`**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★`§9-921` の `htArith_globalToProj_eq_log_mulHeight` は
点が**大域的に**チャートを通ることを要求していた
（`§9-928` の測定によりイデアル類が自明でないと成り立たない）。
★★★★本定理はそれを**各極大イデアルでの局所チャート**に置き換える
——こちらは常に取れる。 -/
theorem htArith_globalToProj_eq_log_mulHeight_of_localCharts
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
    (hgreen : ∀ v : InfinitePlace F,
      ((hyperplaneArith N).comap (globalToProj M hM φ s hcov)).green (archPoint xF v)
        = Real.log ((⨆ k, v ((x k : F))) / v ((x 0 : F)))) :
    htArith F ((hyperplaneArith N).comap (globalToProj M hM φ s hcov)) xF
      = Real.log (Height.mulHeight (fun k => (x k : F))) / (Module.finrank ℚ F : ℝ) :=
  htArith_eq_log_mulHeight_of_localization F _ xF x hx 0 h0
    (hloc_of_localCharts M hM φ s hcov F xF x chartOf y hy hspan hw) hgreen

/-! ## ★出典の紐付け(`.src`) -/

def hloc_of_localCharts.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(局所チャートから局所条件が出る)",
    sectionId := "genell-prop-1-4" }

def htArith_globalToProj_eq_log_mulHeight_of_localCharts.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(局所チャートだけで ht = log H——大域チャートなし)",
    sectionId := "genell-prop-1-4" }

def htArith_globalToProj_eq_log_mulHeight_of_localCharts.needs :
    List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "map_pullbackIdeal_globalToProj_of_localChart(幾何の側、§9-930)"
      (.inProject "ABC3" "ABC3.Found.GenEll.map_pullbackIdeal_globalToProj_of_localChart") 4,
    .citation "[ABC3]" "htArith_eq_log_mulHeight_of_localization(算術の側、§9-931)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_eq_log_mulHeight_of_localization") 4,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-921 の大域チャート仮定" ++
       "(xF : Spec 𝓞_F ⟶ X_{s_i})は**各極大イデアルでの局所チャート**に置き換わった。" ++
       "局所チャートは v_P(x_i) が最小になる i を取れば**常に取れる**ので、" ++
       "これでイデアル類の障害は消えた") 5,
    .implicitStep
      ("★★アルキメデス側の hgreen は仮説のままである。" ++
       "大域チャートの場合(§9-871)は**同じチャートが全ての素点で使えた**ので自動だったが、" ++
       "局所チャートの場合は無限素点ごとにチャートを取り直す必要がある" ++
       "——archPoint xF v : Spec ℂ ⟶ X が X_{s_{i(v)}} を通ること。" ++
       "★体の上の点なので有限素点の側より易しい(座標が 0 でない i を選べばよい)") 4 ]

end ABC3.Found.GenEll
