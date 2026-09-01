/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Found.GaloisRep.DegInfTateParam
import ABC3.Found.GenEll.JScale

/-!
# 第 954 ブロック —— **★★★★★★★★★★★★★★★★悪い素点の局所データを取り出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——(D3) の (b) の半分

`minDeltaExp_eq_mul_of_veluMu`（第 904）は `C`・`hC`・`hc4ne`・`hc4`、
すなわち**極小モデルとその `c₄` が単元であること**を受ける。

★`SemistableAt p W` の定義は

    `minDeltaExp p W = 0` ∨ `∃ C, IsMinimal (primeSubring p) (C • W) ∧ ∃ h, valAdd p c₄ = 0`

だから、**悪い素点（`jExp p W < 0`）では左の選択肢が落ちる**——
`minDeltaExp p W = max 0 (-jExp p W) = -jExp p W > 0` だからである。

☆したがって半安定性だけから 4 つのデータがそのまま出る。

| 定理 | 内容 |
|---|---|
| `minDeltaExp_pos_of_jExp_neg` | ★`jExp < 0` なら `Δ_min > 0` |
| `exists_minimal_c4_unit_of_jExp_neg` | ★★★★★★★★★★★★★★★★**極小モデルと `c₄` の単元性** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ★**半安定で `jExp < 0` なら `Δ_min > 0`**。 -/
theorem minDeltaExp_pos_of_jExp_neg (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] (hss : SemistableAt p W)
    (hj : jExp p W < 0) : 0 < minDeltaExp p W := by
  rw [minDeltaExp_eq_maxJ_of_semistable p W hss]
  omega

/-- ★★★★★★★★★★★★★★★★**悪い素点では半安定性が
極小モデルと `c₄` の単元性を直に与える**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 954）**——これが
`minDeltaExp_eq_mul_of_veluMu`（第 904）の `C`・`hC`・`hc4ne`・`hc4` の出どころである。
☆`SemistableAt` の左の選択肢（`Δ_min = 0`）は `jExp < 0` と矛盾する。 -/
theorem exists_minimal_c4_unit_of_jExp_neg (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] (hss : SemistableAt p W)
    (hj : jExp p W < 0) :
    ∃ C : WeierstrassCurve.VariableChange L, IsMinimal (primeSubring p) (C • W) ∧
      ∃ h : (C • W).c₄ ≠ 0, valAdd p (Units.mk0 ((C • W).c₄) h) = 0 := by
  rcases hss with h0 | hdata
  · exact absurd h0 (minDeltaExp_pos_of_jExp_neg p W (Or.inl h0) hj).ne'
  · exact hdata

/-! ## ★★★★★★★★★★★★第 973 —— `v_p(Δ) = −jExp`（極小モデルで）

★第 909（`hasMultiplicativeReduction_baseChange`）は `0 < v_p(Δ)` を受ける。
☆極小モデルでは `v_p(c₄) = 0` なので、`j = Δ⁻¹c₄³` より `v_p(j) = −v_p(Δ)`、
すなわち `v_p(Δ) = −jExp`。★悪い素点（`jExp < 0`）ではこれが正である。 -/

/-- ★★★★★★★★★★★★**極小モデルでは `v_p(Δ) = −jExp`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`j = Δ⁻¹c₄³` を単元の等式に直し、`valAdd` の乗法性で割るだけである。 -/
theorem valAdd_Delta_eq_neg_jExp (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0) :
    valAdd p (Units.mk0 W.Δ hΔ) = - jExp p W := by
  have hjeq : W.j = W.Δ⁻¹ * W.c₄ ^ 3 := ABC3.Found.GenEll.j_eq_inv_Delta_mul W
  have hj : W.j ≠ 0 := by
    rw [hjeq]
    exact mul_ne_zero (inv_ne_zero hΔ) (pow_ne_zero 3 hc4ne)
  have hunit : Units.mk0 W.j hj
      = (Units.mk0 W.c₄ hc4ne) ^ 3 * (Units.mk0 W.Δ hΔ)⁻¹ := by
    ext
    show W.j = ((Units.mk0 W.c₄ hc4ne) ^ 3 * (Units.mk0 W.Δ hΔ)⁻¹ : Lˣ)
    push_cast
    rw [hjeq]
    show W.Δ⁻¹ * W.c₄ ^ 3 = W.c₄ ^ 3 * (W.Δ)⁻¹
    ring
  have hJ : jExp p W = valAdd p (Units.mk0 W.j hj) := dif_neg hj
  rw [hJ, hunit, valAdd_mul, valAdd_pow, valAdd_inv, hc4]
  omega

/-- ★★★★★★★★**悪い素点では極小モデルの `v_p(Δ)` は正**——第 909 の仮説である。 -/
theorem valAdd_Delta_pos_of_jExp_neg (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0) (hj : jExp p W < 0) :
    0 < valAdd p (Units.mk0 W.Δ hΔ) := by
  rw [valAdd_Delta_eq_neg_jExp p W hΔ hc4ne hc4]
  omega

def valAdd_Delta_eq_neg_jExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極小モデルでは v_p(Δ) = −jExp。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_Delta_pos_of_jExp_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では極小モデルの v_p(Δ) は正。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★出典の紐付け(`.src`) -/

def minDeltaExp_pos_of_jExp_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半安定で jExp < 0 なら Δ_min > 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_minimal_c4_unit_of_jExp_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では半安定性が極小モデルと c₄ の単元性を与える。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
