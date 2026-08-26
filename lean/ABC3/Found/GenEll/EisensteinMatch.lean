import ABC3.Found.GenEll.ZetaTsum

/-!
# GenEll 第 344 ブロック —— **★★★★★★★格子和 = `ζ(K)·2E_K(τ)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★到達点

> **`Σ_{v ≠ 0} (v₁ + v₂τ)⁻ᴷ = ζ(K) · 2·E_K(τ)`**(`K ≥ 3`)

★特に `K = 4, 6` では `ζ(4) = π⁴/90`・`ζ(6) = π⁶/945` を入れた形も取れる。

## ★★★★★2 つの橋

1. **`gammaSet 1 1 0` は原始ベクトルの集合である**(`gammaSet_one_eq`)
   ——`ZMod 1` は自明なので合同条件が消える。
2. **索引の入れ替え**——mathlib の `eisSummand k v z = (v₀z + v₁)^(-k)` に対し
   格子側は `(m + nτ)⁻ᵏ` なので、`swapEquiv`(`w ↦ ![w.2, w.1]`)を噛ませる。
   ★`gcd` は対称なので原始性は保たれる。

★★★`ModularForm.E hk = (1/2) • eisensteinSeriesMF hk 0` の中身は
`.copy` で包まれているが、`simp only [E, copy, eisensteinSeriesMF]` のあと `rfl` で開く。

## ★★★★★★★これで (i) に残るのは 1 つだけ

    g₂ = 60·G₄ = 60·ζ(4)·2E₄ = (4π⁴/3)·E₄       ★本ブロックまでで出る
    g₃ = 140·G₆ = 140·ζ(6)·2E₆ = (8π⁶/27)·E₆    ★同上
    g₂³ − 27g₃² = (64π¹²/27)·(E₄³ − E₆²)         ★係数が一致する(第 341 で確認)
    E₄³ − E₆² = 1728·Δ                           ★**残る 1 つ**(mathlib に無い)
    Δ = η²⁴ ≠ 0                                  ★mathlib

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `gammaSet_one_eq` | ★★★`gammaSet 1 1 0` は原始ベクトルの集合 |
| `swapEquiv_zero`・`swapEquiv_one` | ★索引の成分 |
| `primEquivGammaSet` | ★★★★原始ベクトル ≃ `gammaSet 1 1 0` |
| `tsum_primitive_eq_eisensteinSeries` | ★★★★★和の同定 |
| `E_apply_eq` | ★★★`E_K = (1/2)·eisensteinSeries` |
| `tsum_primitive_eq_two_mul_E` | ★★★★★★**`Σ_{w 原始} = 2E_K`** |
| `tsum_nonzero_eq_zeta_mul_E` | ★★★★★★★**格子和 = `ζ(K)·2E_K`** |
| `tsum_nonzero_four`・`tsum_nonzero_six` | ★★★★`K = 4, 6` の値入り |
-/

namespace ABC3.Found.GenEll

open EisensteinSeries

/-! ## ★★★`gammaSet 1 1 0` は原始ベクトルの集合 -/

/-- ★★★**`gammaSet 1 1 0` は原始ベクトルの集合である**——`ZMod 1` は自明。 -/
theorem gammaSet_one_eq :
    EisensteinSeries.gammaSet 1 1 (0 : Fin 2 → ZMod 1)
      = {v : Fin 2 → ℤ | Int.gcd (v 0) (v 1) = 1} := by
  ext v
  simp only [EisensteinSeries.gammaSet, Set.mem_setOf_eq]
  constructor
  · rintro ⟨-, h⟩; exact h
  · intro h
    refine ⟨?_, h⟩
    funext i
    exact Subsingleton.elim _ _

/-- ★索引の第 0 成分。 -/
theorem swapEquiv_zero (w : ℤ × ℤ) : (swapEquiv w) 0 = w.2 := rfl

/-- ★索引の第 1 成分。 -/
theorem swapEquiv_one (w : ℤ × ℤ) : (swapEquiv w) 1 = w.1 := rfl

/-- ★★★★原始ベクトルと `gammaSet 1 1 0` の全単射。 -/
noncomputable def primEquivGammaSet :
    {w : ℤ × ℤ // Int.gcd w.1 w.2 = 1}
      ≃ (EisensteinSeries.gammaSet 1 1 (0 : Fin 2 → ZMod 1)) :=
  swapEquiv.subtypeEquiv (fun w => by
    rw [gammaSet_one_eq]
    show Int.gcd w.1 w.2 = 1 ↔ Int.gcd ((swapEquiv w) 0) ((swapEquiv w) 1) = 1
    rw [swapEquiv_zero, swapEquiv_one, Int.gcd_comm])

/-! ## ★★★★★★和の同定 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★**原始ベクトル上の和は `eisensteinSeries` である**。 -/
theorem tsum_primitive_eq_eisensteinSeries (τ : UpperHalfPlane) (K : ℕ) :
    ∑' w : {w : ℤ × ℤ // Int.gcd w.1 w.2 = 1}, (((w.1.1 : ℂ) + (w.1.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹
      = eisensteinSeries (0 : Fin 2 → ZMod 1) (K : ℤ) τ := by
  rw [eisensteinSeries, ← primEquivGammaSet.tsum_eq
    (fun x : (EisensteinSeries.gammaSet 1 1 (0 : Fin 2 → ZMod 1)) =>
      EisensteinSeries.eisSummand (K : ℤ) (x : Fin 2 → ℤ) τ)]
  refine tsum_congr (fun w => ?_)
  show _ = EisensteinSeries.eisSummand (K : ℤ) (swapEquiv w.1) τ
  rw [EisensteinSeries.eisSummand, swapEquiv_zero, swapEquiv_one, zpow_neg, zpow_natCast]
  congr 2
  ring

/-- ★★★`E_K = (1/2)·eisensteinSeries`。 -/
theorem E_apply_eq (τ : UpperHalfPlane) (K : ℕ) (hK : 3 ≤ K) :
    ModularForm.E (k := K) hK τ = (1/2 : ℂ) * eisensteinSeries (0 : Fin 2 → ZMod 1) (K : ℤ) τ := by
  simp only [ModularForm.E, ModularForm.copy, ModularForm.eisensteinSeriesMF]
  rfl

/-- ★★★★★★**`Σ_{w 原始} (w₁+w₂τ)⁻ᴷ = 2·E_K(τ)`**。 -/
theorem tsum_primitive_eq_two_mul_E (τ : UpperHalfPlane) (K : ℕ) (hK : 3 ≤ K) :
    ∑' w : {w : ℤ × ℤ // Int.gcd w.1 w.2 = 1}, (((w.1.1 : ℂ) + (w.1.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹
      = 2 * ModularForm.E (k := K) hK τ := by
  rw [tsum_primitive_eq_eisensteinSeries, E_apply_eq τ K hK]
  ring

/-! ## ★★★★★★★格子和 = `ζ(K)·2E_K` -/

/-- ★★★★★★★**`Σ_{v ≠ 0} (v₁+v₂τ)⁻ᴷ = ζ(K)·2·E_K(τ)`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
theorem tsum_nonzero_eq_zeta_mul_E (τ : UpperHalfPlane) (K : ℕ) (hK : 3 ≤ K) :
    ∑' v : {v : ℤ × ℤ // v ≠ (0, 0)}, (((v.1.1 : ℂ) + (v.1.2 : ℂ) * (τ : ℂ)) ^ K)⁻¹
      = (∑' d : {d : ℕ // 0 < d}, (((d : ℕ) : ℂ) ^ K)⁻¹) * (2 * ModularForm.E (k := K) hK τ) := by
  rw [tsum_nonzero_eq_mul τ K hK, tsum_primitive_eq_two_mul_E τ K hK]

/-- ★★★★`K = 4` の値入り。 -/
theorem tsum_nonzero_four (τ : UpperHalfPlane) :
    ∑' v : {v : ℤ × ℤ // v ≠ (0, 0)}, (((v.1.1 : ℂ) + (v.1.2 : ℂ) * (τ : ℂ)) ^ 4)⁻¹
      = ((Real.pi ^ 4 / 90 : ℝ) : ℂ) * (2 * ModularForm.E (k := 4) (by norm_num) τ) := by
  rw [tsum_nonzero_eq_zeta_mul_E τ 4 (by norm_num), tsum_natPow_inv_four]

/-- ★★★★`K = 6` の値入り。 -/
theorem tsum_nonzero_six (τ : UpperHalfPlane) :
    ∑' v : {v : ℤ × ℤ // v ≠ (0, 0)}, (((v.1.1 : ℂ) + (v.1.2 : ℂ) * (τ : ℂ)) ^ 6)⁻¹
      = ((Real.pi ^ 6 / 945 : ℝ) : ℂ) * (2 * ModularForm.E (k := 6) (by norm_num) τ) := by
  rw [tsum_nonzero_eq_zeta_mul_E τ 6 (by norm_num), tsum_natPow_inv_six]

/-! ## ★出典の紐付け(`.src`) -/

def tsum_nonzero_eq_zeta_mul_E.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def gammaSet_one_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
