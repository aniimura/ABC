import ABC3.Found.GenEll.GcdProd
import Mathlib.NumberTheory.ModularForms.EisensteinSeries.QExpansion

/-!
# GenEll 第 340 ブロック —— **★★★★★★格子和の絶対収束**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★mathlib の Eisenstein 級数の収束がそのまま効いた

(i) の残りは解析——`Σ_d d⁻ᵏ` と `Σ_{w 原始} (w₁+w₂τ)⁻ᵏ` の絶対収束であった。

★★★**mathlib は `EisensteinSeries.summable_norm_eisSummand` を持っており、しかも
それは `gammaSet` に制限した和ではなく `Fin 2 → ℤ` **全体**についての絶対収束である**
(2026-08-26 実測):

    3 ≤ k → ∀ z, Summable fun x : Fin 2 → ℤ => ‖eisSummand k x z‖

★★したがって格子和の絶対収束は**そのまま**得られ、部分型への制限も `Summable.subtype`
で出る。★★★当初「絶対収束を自分で積む」と見ていたが、**在庫で足りた**。

## ★★★★索引の突き合わせ

mathlib の `eisSummand k x z = (x₀·z + x₁)^(-k)` に対し、格子側は `(m + nτ)⁻ᵏ`。
★**順序が `(m,n) ↔ (x₁,x₀)` で入れ替わる**ので、
`swapEquiv : ℤ × ℤ ≃ (Fin 2 → ℤ)`(`v ↦ ![v.2, v.1]`)を噛ませる。

## ★★残り((i) の最後)

★`tsum_mul_tsum_of_summable_norm` で積に分け、`Σ_d d⁻ᵏ = ζ(k)` と
`Σ_{w 原始} = 2E_k` を突き合わせる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `summable_natPow_inv` | ★★★`Σ_{d>0} d⁻ᵏ` は可和(`k ≥ 2`) |
| `swapEquiv` | ★索引の入れ替え |
| `summable_lattice_inv` | ★★★★★★**`Σ_{v ∈ ℤ²} (v₁+v₂τ)⁻ᵏ` は可和(`k ≥ 3`)** |
| `summable_primitive_inv`・`summable_nonzero_inv` | ★★★部分型への制限 |
-/

namespace ABC3.Found.GenEll

open EisensteinSeries

/-! ## ★★★`Σ d⁻ᵏ` の可和性 -/

/-- ★★★`Σ_{d>0} d⁻ᵏ` は可和である(`k ≥ 2`)。 -/
theorem summable_natPow_inv (k : ℕ) (hk : 2 ≤ k) :
    Summable (fun d : {d : ℕ // 0 < d} => (((d : ℕ) : ℂ) ^ k)⁻¹) := by
  have h : Summable (fun n : ℕ => 1 / ((n : ℝ) ^ k)) := by
    rw [Real.summable_one_div_nat_pow]
    omega
  refine Summable.of_norm ?_
  have h2 : Summable (fun d : {d : ℕ // 0 < d} => 1 / (((d : ℕ) : ℝ) ^ k)) := h.subtype _
  refine h2.congr (fun d => ?_)
  rw [norm_inv, norm_pow, Complex.norm_natCast, one_div]

/-! ## ★★★★★★格子和の可和性 -/

/-- ★索引の入れ替え `v ↦ ![v.2, v.1]`。 -/
noncomputable def swapEquiv : ℤ × ℤ ≃ (Fin 2 → ℤ) :=
  (Equiv.prodComm ℤ ℤ).trans (piFinTwoEquiv (fun _ => ℤ)).symm

/-- ★★★★★★**格子和は絶対収束する**(`k ≥ 3`)。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★mathlib の `summable_norm_eisSummand` は `Fin 2 → ℤ` 全体についての絶対収束なので、
索引を入れ替えるだけで出る。 -/
theorem summable_lattice_inv (τ : UpperHalfPlane) (K : ℕ) (hK : 3 ≤ K) :
    Summable (fun v : ℤ × ℤ => (((v.1 : ℂ) + (v.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹) := by
  have h := EisensteinSeries.summable_norm_eisSummand (k := (K : ℤ)) (by exact_mod_cast hK) τ
  have h2 : Summable (fun x : Fin 2 → ℤ => EisensteinSeries.eisSummand (K : ℤ) x τ) :=
    Summable.of_norm h
  have h3 := (swapEquiv.summable_iff (f := fun x : Fin 2 → ℤ =>
    EisensteinSeries.eisSummand (K : ℤ) x τ)).2 h2
  refine h3.congr (fun v => ?_)
  show EisensteinSeries.eisSummand (K : ℤ) (swapEquiv v) τ = _
  rw [EisensteinSeries.eisSummand]
  show (((v.2 : ℂ)) * (τ : ℂ) + ((v.1 : ℂ))) ^ (-(K : ℤ)) = _
  rw [zpow_neg, zpow_natCast]
  congr 2
  ring

/-- ★★★原始ベクトルへの制限。 -/
theorem summable_primitive_inv (τ : UpperHalfPlane) (K : ℕ) (hK : 3 ≤ K) :
    Summable (fun w : {w : ℤ × ℤ // Int.gcd w.1 w.2 = 1} =>
      (((w.1.1 : ℂ) + (w.1.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹) :=
  (summable_lattice_inv τ K hK).subtype _

/-- ★★★零でない格子点への制限。 -/
theorem summable_nonzero_inv (τ : UpperHalfPlane) (K : ℕ) (hK : 3 ≤ K) :
    Summable (fun v : {v : ℤ × ℤ // v ≠ (0, 0)} =>
      (((v.1.1 : ℂ) + (v.1.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹) :=
  (summable_lattice_inv τ K hK).subtype _

/-! ## ★出典の紐付け(`.src`) -/

def summable_lattice_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def summable_natPow_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
