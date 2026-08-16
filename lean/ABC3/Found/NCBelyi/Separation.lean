import ABC3.Meta.Claim
import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Data.Complex.Basic
import Mathlib.Analysis.Normed.Field.Basic
import Mathlib.Order.Filter.AtTopBot.Basic

/-!
# [NCBelyi] Lemma 2.3 —— Separation of Collections of Points(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.4。
**260 dpi 目視確認 2026-08-17**。

原文 (NCBelyi p.4):
> Lemma 2.3. (Separation of Collections of Points) Let

## ★★原文の証明は 1 行である

原文 (NCBelyi p.4):
> Indeed, it suffices to take λ such |λ − β| is sufficiently small.

★★**「sufficiently small」を明示した**——それがこのファイルの仕事である。

`δ ≝ min_{α ∈ S ∩ ℂ} ‖α − β‖ > 0` と置き、
**`0 < ε ≤ δ/4` かつ `C·ε ≤ δ/4`** なる実数 `ε` を取って `λ ≝ β + ε` とすればよい。
そのとき

- `‖α − λ‖ ≥ ‖α − β‖ − ε ≥ δ − δ/4 = 3δ/4`
- `C·‖1/(α−λ)‖ = C·ε/(ε·‖α−λ‖) ≤ (δ/4)/(ε·(δ/4))·… `

——要点は **`C·ε ≤ δ/4 ≤ ‖α − λ‖`**、すなわち `C/‖α−λ‖ ≤ 1/ε = ‖1/(β−λ)‖`。

## ★`ℙ¹(ℂ)` の模し方(明示する)

★`ℙ¹(ℂ)` を **`Option ℂ`**(`none` が `∞`)として模した。
本補題が使うのは**集合としての `ℙ¹(ℂ)` と `f` の値**だけなので、これで足りる。

★★**`∞ ∈ S` の場合を落としていない。** `f(∞) = 0` なので
`|f(β)| ≥ C·0` は自明に成り立つが、**落とすと主張が弱くなる**ので含めてある。
-/

namespace ABC3.Found.NCBelyi

/-- `ℙ¹(ℂ)` を `Option ℂ` として模す。★`none` が `∞`。 -/
abbrev P1C := Option ℂ

/-- `f(x) = 1/(x − λ)` の絶対値。★`f(∞) = 0` なので `∞` では `0`。 -/
noncomputable def absInvShift (lam : ℂ) : P1C → ℝ
  | none => 0
  | some x => ‖1 / (x - lam)‖

@[simp] theorem absInvShift_none (lam : ℂ) : absInvShift lam none = 0 := rfl

@[simp] theorem absInvShift_some (lam x : ℂ) :
    absInvShift lam (some x) = ‖1 / (x - lam)‖ := rfl

/-- ★**核となる評価** —— 「sufficiently small」を `ε` で明示した形。

`δ` が `S` の有限点と `β` の距離の下界で、`ε ≤ δ/4` かつ `C·ε ≤ δ/4` なら、
`λ ≝ β + ε` が原文の 3 条件を満たす。 -/
theorem separation_core (S : Finset P1C) (C : ℝ) (hC : 0 < C) (b : ℂ)
    (delta eps : ℝ) (hdelta : 0 < delta) (heps : 0 < eps)
    (hS : ∀ p ∈ S, ∀ x : ℂ, p = some x → delta ≤ ‖x - b‖)
    (h1 : eps ≤ delta / 4) (h2 : C * eps ≤ delta / 4) :
    (b + (eps : ℂ)) ≠ b
      ∧ (some (b + (eps : ℂ)) ∉ S)
      ∧ ∀ p ∈ S, C * absInvShift (b + (eps : ℂ)) p
          ≤ absInvShift (b + (eps : ℂ)) (some b) := by
  set lam : ℂ := b + (eps : ℂ) with hlam
  have hnorm : ‖b - lam‖ = eps := by
    rw [hlam]
    simp only [sub_add_cancel_left, norm_neg, Complex.norm_real, Real.norm_eq_abs]
    exact abs_of_pos heps
  have hne : lam ≠ b := by
    intro hcon
    have : ‖b - lam‖ = 0 := by rw [hcon, sub_self, norm_zero]
    rw [hnorm] at this
    exact absurd this heps.ne'
  -- `S` の有限点は `λ` から `3δ/4` 以上離れている
  have hfar : ∀ p ∈ S, ∀ x : ℂ, p = some x → delta / 4 ≤ ‖x - lam‖ := by
    intro p hp x hpx
    have hd := hS p hp x hpx
    have htri : ‖x - b‖ ≤ ‖x - lam‖ + ‖lam - b‖ := by
      have : x - b = (x - lam) + (lam - b) := by ring
      rw [this]
      exact norm_add_le _ _
    have hlb : ‖lam - b‖ = eps := by rw [← norm_neg]; simpa using hnorm
    rw [hlb] at htri
    linarith
  refine ⟨hne, ?_, ?_⟩
  · intro hmem
    have := hfar _ hmem lam rfl
    rw [sub_self, norm_zero] at this
    linarith
  · intro p hp
    cases p with
    | none => simp only [absInvShift_none, mul_zero, absInvShift_some]; positivity
    | some x =>
      have hx := hfar _ hp x rfl
      have hxlam : x - lam ≠ 0 := by
        intro hcon
        rw [hcon, norm_zero] at hx
        linarith
      have hxpos : (0 : ℝ) < ‖x - lam‖ := norm_pos_iff.2 hxlam
      have hbpos : (0 : ℝ) < ‖b - lam‖ := by rw [hnorm]; exact heps
      simp only [absInvShift_some, norm_div, norm_one, mul_one_div]
      rw [div_le_div_iff₀ hxpos hbpos, hnorm, one_mul]
      -- 残るのは `C·ε ≤ ‖x − λ‖`
      exact le_trans h2 hx

/-- **[NCBelyi] Lemma 2.3**(Separation of Collections of Points)。

原文 (NCBelyi p.4):
> Lemma 2.3. (Separation of Collections of Points) Let

`S ⊆ ℙ¹(ℂ)` を複素点の有限集合、`C > 0`、`β ∈ ℂ\S` とすると、
**`λ ∈ ℂ` が存在して** `f(x) = 1/(x − λ)` が
`f(β) ≠ 0, ∞`(= `λ ≠ β`)、`f(α) ≠ ∞`(= `λ ∉ S`)、
**`|f(β)| ≥ C·|f(α)|`**(すべての `α ∈ S`)を満たす。

★原文の証明は「`|λ − β|` を十分小さく取ればよい」の 1 行。
**本実装はその `ε` を明示的に構成する**(`separation_core`)。 -/
theorem lemma_2_3 (S : Finset P1C) (C : ℝ) (hC : 0 < C) (b : ℂ)
    (hb : (some b) ∉ S) :
    ∃ lam : ℂ, lam ≠ b ∧ (some lam ∉ S)
      ∧ ∀ p ∈ S, C * absInvShift lam p ≤ absInvShift lam (some b) := by
  classical
  -- `δ` を「`S` の点と `b` の距離」の下界として作る(`∞` は `1` に潰す)
  set g : P1C → ℝ := fun p => match p with | none => 1 | some x => ‖x - b‖ with hg
  set D : Finset ℝ := S.image g with hD
  have hgpos : ∀ p ∈ S, 0 < g p := by
    intro p hp
    cases p with
    | none => simp [hg]
    | some x =>
      have hxb : x ≠ b := by
        intro hcon; rw [hcon] at hp; exact hb hp
      simp only [hg]
      exact norm_pos_iff.2 (sub_ne_zero.2 hxb)
  set delta : ℝ := if h : D.Nonempty then D.min' h else 1 with hdel
  have hdelta : 0 < delta := by
    simp only [hdel]
    split_ifs with h
    · obtain ⟨p, hp, hpe⟩ := Finset.mem_image.1 (D.min'_mem h)
      rw [← hpe]
      exact hgpos p hp
    · norm_num
  have hSle : ∀ p ∈ S, ∀ x : ℂ, p = some x → delta ≤ ‖x - b‖ := by
    intro p hp x hpx
    have hmem : g p ∈ D := Finset.mem_image_of_mem g hp
    have hne : D.Nonempty := ⟨g p, hmem⟩
    have : delta ≤ g p := by simp only [hdel, dif_pos hne]; exact D.min'_le _ hmem
    rw [hpx] at this
    simpa [hg] using this
  refine ⟨b + ((min (delta / 4) (delta / (4 * C)) : ℝ) : ℂ), ?_⟩
  have heps : 0 < min (delta / 4) (delta / (4 * C)) :=
    lt_min (by positivity) (by positivity)
  have h1 : min (delta / 4) (delta / (4 * C)) ≤ delta / 4 := min_le_left _ _
  have h2 : C * min (delta / 4) (delta / (4 * C)) ≤ delta / 4 := by
    have := min_le_right (delta / 4) (delta / (4 * C))
    calc C * min (delta / 4) (delta / (4 * C)) ≤ C * (delta / (4 * C)) := by
          exact mul_le_mul_of_nonneg_left this hC.le
      _ = delta / 4 := by field_simp
  exact separation_core S C hC b delta _ hdelta heps hSle h1 h2

/-! ## ★出典の紐付け(`.src`) -/

/-- ★**条なし** —— `Lemma 2.3` は主張全体を実装した。 -/
def lemma_2_3.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4, item := "Lemma 2.3",
    sectionId := "ncbelyi-lemma-2-3" }

end ABC3.Found.NCBelyi
