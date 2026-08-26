import ABC3.Found.GenEll.LatticeSummable

/-!
# GenEll 第 341 ブロック —— **★★★★★★★格子和が積に分かれた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★到達点

> **`Σ_{v ≠ 0} (v₁ + v₂τ)⁻ᴷ = (Σ_{d>0} d⁻ᴷ) · (Σ_{w 原始} (w₁ + w₂τ)⁻ᴷ)`**(`K ≥ 3`)

★これが第 337-340 の 4 ブロックの帰結である:
gcd 分解(337)→ 全単射(338)→ 積の形(339)→ 絶対収束(340)→ **Fubini(本ブロック)**。

## ★★★★★★残る 1 段——ζ の値との突き合わせ

`g₂ = 60·G₄`、`g₃ = 140·G₆` と `G_k = ζ(k)·(2E_k)` から

    g₂³ − 27g₃² = (120ζ(4))³·E₄³ − 27·(280ζ(6))²·E₆²

★★ここで **`ζ(4) = π⁴/90`・`ζ(6) = π⁶/945`** を入れると

    (120·π⁴/90)³ = (4π⁴/3)³ = 64π¹²/27
    27·(280·π⁶/945)² = 27·(8π⁶/27)² = 64π¹²/27

と**両者が一致する**ので、`g₂³ − 27g₃² = (64π¹²/27)(E₄³ − E₆²) = (2π)¹²·Δ`。
★★★★`Δ = η²⁴ ≠ 0`(mathlib)から **(i) 判別式の非消失が閉じる**。

★係数が一致することが本質的である——一致しなければ `Δ` の非消失だけでは足りない。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `summable_norm_natPow_inv` | ★★`Σ‖d⁻ᴷ‖` は可和 |
| `summable_norm_lattice_inv` | ★★★`Σ‖(v₁+v₂τ)⁻ᴷ‖` は可和 |
| `summable_norm_primitive_inv` | ★★原始への制限 |
| `tsum_nonzero_eq_mul` | ★★★★★★★**格子和 = (Σd⁻ᴷ)·(Σ原始)** |
-/

namespace ABC3.Found.GenEll

open EisensteinSeries

/-! ## ★★ノルムの可和性 -/

/-- ★★`Σ ‖d⁻ᴷ‖` は可和(`K ≥ 2`)。 -/
theorem summable_norm_natPow_inv (k : ℕ) (hk : 2 ≤ k) :
    Summable (fun d : {d : ℕ // 0 < d} => ‖(((d : ℕ) : ℂ) ^ k)⁻¹‖) := by
  have h : Summable (fun n : ℕ => 1 / ((n : ℝ) ^ k)) := by
    rw [Real.summable_one_div_nat_pow]; omega
  have h2 : Summable (fun d : {d : ℕ // 0 < d} => 1 / (((d : ℕ) : ℝ) ^ k)) := h.subtype _
  refine h2.congr (fun d => ?_)
  rw [norm_inv, norm_pow, Complex.norm_natCast, one_div]

/-- ★★★`Σ ‖(v₁+v₂τ)⁻ᴷ‖` は可和(`K ≥ 3`)。 -/
theorem summable_norm_lattice_inv (τ : UpperHalfPlane) (K : ℕ) (hK : 3 ≤ K) :
    Summable (fun v : ℤ × ℤ => ‖(((v.1 : ℂ) + (v.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹‖) := by
  have h := EisensteinSeries.summable_norm_eisSummand (k := (K : ℤ)) (by exact_mod_cast hK) τ
  have h3 := (swapEquiv.summable_iff (f := fun x : Fin 2 → ℤ =>
    ‖EisensteinSeries.eisSummand (K : ℤ) x τ‖)).2 h
  refine h3.congr (fun v => ?_)
  show ‖EisensteinSeries.eisSummand (K : ℤ) (swapEquiv v) τ‖ = _
  congr 1
  rw [EisensteinSeries.eisSummand]
  show (((v.2 : ℂ)) * (τ : ℂ) + ((v.1 : ℂ))) ^ (-(K : ℤ)) = _
  rw [zpow_neg, zpow_natCast]
  congr 2
  ring

/-- ★★原始ベクトルへの制限。 -/
theorem summable_norm_primitive_inv (τ : UpperHalfPlane) (K : ℕ) (hK : 3 ≤ K) :
    Summable (fun w : {w : ℤ × ℤ // Int.gcd w.1 w.2 = 1} =>
      ‖(((w.1.1 : ℂ) + (w.1.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹‖) :=
  (summable_norm_lattice_inv τ K hK).subtype _

/-! ## ★★★★★★★Fubini -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**格子和は `(Σ_{d>0} d⁻ᴷ)·(Σ_{w 原始} (w₁+w₂τ)⁻ᴷ)` に分かれる**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★第 337-340 の帰結:gcd 分解 → 全単射 → 積の形 → 絶対収束 → Fubini。 -/
theorem tsum_nonzero_eq_mul (τ : UpperHalfPlane) (K : ℕ) (hK : 3 ≤ K) :
    ∑' v : {v : ℤ × ℤ // v ≠ (0, 0)}, (((v.1.1 : ℂ) + (v.1.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹
      = (∑' d : {d : ℕ // 0 < d}, (((d : ℕ) : ℂ) ^ K)⁻¹)
        * (∑' w : {w : ℤ × ℤ // Int.gcd w.1 w.2 = 1},
            (((w.1.1 : ℂ) + (w.1.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹) := by
  rw [tsum_mul_tsum_of_summable_norm (summable_norm_natPow_inv K (by omega))
    (summable_norm_primitive_inv τ K hK)]
  rw [← nonzeroProdEquiv.tsum_eq (fun v : {v : ℤ × ℤ // v ≠ (0, 0)} =>
    (((v.1.1 : ℂ) + (v.1.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹)]
  exact tsum_congr (fun q => summand_nonzeroProdEquiv (τ : ℂ) K q)

/-! ## ★出典の紐付け(`.src`) -/

def tsum_nonzero_eq_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def summable_norm_lattice_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
