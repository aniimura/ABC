import ABC3.Found.GenEll.Covolume

/-!
# GenEll 第 350 ブロック —— **★★★★★★不変量は束だけで決まる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★なぜ要るのか——`htFalt` の一意性

第 349 で `archInv = ‖Δ_lat‖·covol⁶` が**相似**で変わらないことを出した。
★しかし `htFalt` を固定するには、もう一段——**基底の取り替えでも変わらない**、
すなわち `archInv` が**周期対ではなく束の関数である**——が要る。

★★`covol` は `ω₁, ω₂` の行列式で定義されているので、これは自明ではない。
★★★同じ束の 2 つの基底は `GL₂(ℤ)` で移り合い、行列式が `±1` なので `|det|` は等しい
——これが `covol_congr` である。

## ★★★★★★★次の障害の**測定**(2026-08-26)

`htFalt` を固定する最後の一歩は

> `C • W = latticeCurve P` と `C' • W = latticeCurve P'` ⟹ `archInv P = archInv P'`

である。★`D := C'C⁻¹` の `u` について `g₂(P') = u⁻⁴g₂(P)`・`g₃(P') = u⁻⁶g₃(P)`、
すなわち `P'` と `scalePair P u` は**同じ `(g₂,g₃)`** を持つ。
★★したがって残るのは **`(g₂,g₃)` が束を決める**ことである。

### ★★★★★★道は測れている

| 段 | 状態 |
|---|---|
| `℘` が同じ ⟹ 束が同じ | ★**本ブロックで完了**(`lattice_eq_of_weierstrassP_eq`) |
| 束が同じ ⟹ `covol`・`archInv` が同じ | ★**本ブロックで完了** |
| `(g₂,g₃)` が同じ ⟹ `℘` が同じ | ★残る 1 つ |

★★★最後の段の道具は **mathlib にそろっている**:
`hasFPowerSeriesAt_weierstrassPExcept` は `℘[L−0]` の 0 での冪級数係数が
**`(i+1)·G(i+2)`** であることを与え(`sumInvPow_zero : sumInvPow 0 = G`)、
`derivWeierstrassP_sq` は `℘'² = 4℘³ − g₂℘ − g₃` を与える。
★★★★`f := ℘[L−0]` と置くと `℘ = f + z⁻²` で、微分方程式は

    z²f'' = 6z²f² + 12f − (g₂/2)z²

という**解析的な恒等式**になる。係数を比べると `a_i = (i+1)G(i+2)` について

    (i−4)(i+3)·a_i = 6·Σ_{j+k=i−2} a_j a_k − (g₂/2)·[i=2]

★★★★★`i = 4` で係数が消えるのは `a_4 ↔ G₆ = g₃/140` が**2 つ目の自由母数**だから
であり、`i ≥ 5` では `(i−4)(i+3) ≠ 0` なので **`a_i` は `g₂, g₃` で決まる**。
★見積もり **4-8 ブロック**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `det2`・`det2_comb` | ★実 2 次元行列式と基底取り替えの規則 |
| `covol_congr` | ★★★★★**共体積は束だけで決まる** |
| `archInv_congr` | ★★★★★★**`archInv` は束だけで決まる** |
| `lattice_eq_discontinuity` | ★★★★★**束は `℘` の不連続点の集合** |
| `lattice_eq_of_weierstrassP_eq` | ★★★★★★**`℘` が同じなら束も同じ** |
| `archInv_eq_of_weierstrassP_eq` | ★★★★★★同じく `archInv` も |
-/

namespace ABC3.Found.GenEll

open Complex Real PeriodPair

/-! ## ★実 2 次元の行列式 -/

/-- ★実 2 次元の行列式 `x.re·y.im − x.im·y.re`。 -/
def det2 (x y : ℂ) : ℝ := x.re * y.im - x.im * y.re

theorem covol_eq_det2 (L : PeriodPair) : covol L = |det2 L.ω₁ L.ω₂| := rfl

/-- ★★**基底を `ℤ` 係数で取り替えると行列式は `ad − bc` 倍**。 -/
theorem det2_comb (a b c d : ℤ) (u v : ℂ) :
    det2 ((a : ℤ) • u + (b : ℤ) • v) ((c : ℤ) • u + (d : ℤ) • v)
      = ((a * d - b * c : ℤ) : ℝ) * det2 u v := by
  simp only [det2, zsmul_eq_mul, Complex.add_re, Complex.add_im, Complex.mul_re, Complex.mul_im,
    Complex.intCast_re, Complex.intCast_im]
  push_cast
  ring

theorem omega_mem (L : PeriodPair) : L.ω₁ ∈ L.lattice ∧ L.ω₂ ∈ L.lattice := by
  constructor <;> (rw [PeriodPair.lattice]; exact Submodule.subset_span (by simp))

theorem det2_ne_zero (L : PeriodPair) : det2 L.ω₁ L.ω₂ ≠ 0 := by
  have := covol_pos L
  rw [covol_eq_det2, abs_pos] at this
  exact this

/-! ## ★★★★★束だけで決まること -/

/-- ★★★★★**共体積は束だけで決まる**——基底の取り替えは `GL₂(ℤ)` で行列式が `±1`。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★両向きの `ℤ` 係数表示を取ると `(ad−bc)(ps−qr) = 1` になり、整数の単元は `±1` である。 -/
theorem covol_congr {L L' : PeriodPair} (h : L.lattice = L'.lattice) : covol L = covol L' := by
  have h1 : L.ω₁ ∈ L'.lattice := h ▸ (omega_mem L).1
  have h2 : L.ω₂ ∈ L'.lattice := h ▸ (omega_mem L).2
  have h1' : L'.ω₁ ∈ L.lattice := h ▸ (omega_mem L').1
  have h2' : L'.ω₂ ∈ L.lattice := h ▸ (omega_mem L').2
  rw [PeriodPair.lattice] at h1 h2 h1' h2'
  obtain ⟨a, b, hab⟩ := Submodule.mem_span_pair.1 h1
  obtain ⟨c, d, hcd⟩ := Submodule.mem_span_pair.1 h2
  obtain ⟨p, q, hpq⟩ := Submodule.mem_span_pair.1 h1'
  obtain ⟨r, s, hrs⟩ := Submodule.mem_span_pair.1 h2'
  have e1 : det2 L.ω₁ L.ω₂ = ((a * d - b * c : ℤ) : ℝ) * det2 L'.ω₁ L'.ω₂ := by
    rw [← hab, ← hcd, det2_comb]
  have e2 : det2 L'.ω₁ L'.ω₂ = ((p * s - q * r : ℤ) : ℝ) * det2 L.ω₁ L.ω₂ := by
    rw [← hpq, ← hrs, det2_comb]
  have hne := det2_ne_zero L
  have h3 : det2 L.ω₁ L.ω₂
      = (((a * d - b * c : ℤ) : ℝ) * ((p * s - q * r : ℤ) : ℝ)) * det2 L.ω₁ L.ω₂ := by
    conv_lhs => rw [e1]
    rw [e2]; ring
  have hprod : ((a * d - b * c : ℤ) : ℝ) * ((p * s - q * r : ℤ) : ℝ) = 1 :=
    mul_right_cancel₀ hne (by rw [one_mul]; exact h3.symm)
  have hz : (a * d - b * c : ℤ) * (p * s - q * r : ℤ) = 1 := by exact_mod_cast hprod
  have habs : |((a * d - b * c : ℤ) : ℝ)| = 1 := by
    rcases Int.eq_one_or_neg_one_of_mul_eq_one hz with h' | h' <;> simp [h']
  rw [covol_eq_det2, covol_eq_det2, e1, abs_mul, habs, one_mul]

/-- ★★★★★★**アルキメデス不変量は束だけで決まる**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★これで `archInv` は周期対ではなく**束の関数**になった。 -/
theorem archInv_congr {L L' : PeriodPair} (h : L.lattice = L'.lattice) :
    archInv L = archInv L' := by
  rw [archInv, archInv, covol_congr h, latticeDisc_congr h]

/-! ## ★★★★★★束は `℘` の不連続点の集合 -/

/-- ★★★★★**束は `℘` の不連続点の集合である**——極は束の点にちょうど乗る。 -/
theorem lattice_eq_discontinuity (L : PeriodPair) :
    (L.lattice : Set ℂ) = {z : ℂ | ¬ ContinuousAt (PeriodPair.weierstrassP L) z} := by
  ext z
  constructor
  · intro hz
    exact L.not_continuousAt_weierstrassP z hz
  · intro hz
    by_contra hnot
    exact hz ((L.analyticOnNhd_weierstrassP z hnot).continuousAt)

/-- ★★★★★★**`℘` が同じなら束も同じ**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★不連続点の集合が束そのものだから。 -/
theorem lattice_eq_of_weierstrassP_eq {L L' : PeriodPair}
    (h : PeriodPair.weierstrassP L = PeriodPair.weierstrassP L') : L.lattice = L'.lattice := by
  apply SetLike.coe_injective
  rw [lattice_eq_discontinuity L, lattice_eq_discontinuity L', h]

/-- ★★★★★★**`℘` が同じならアルキメデス不変量も同じ**。 -/
theorem archInv_eq_of_weierstrassP_eq {L L' : PeriodPair}
    (h : PeriodPair.weierstrassP L = PeriodPair.weierstrassP L') : archInv L = archInv L' :=
  archInv_congr (lattice_eq_of_weierstrassP_eq h)

/-! ## ★出典の紐付け(`.src`) -/

def covol_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def archInv_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def lattice_eq_of_weierstrassP_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
