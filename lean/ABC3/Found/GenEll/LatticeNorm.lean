import ABC3.Found.GenEll.LatticeScale

/-!
# GenEll 第 334 ブロック —— **★★★★★正規化された束と `ℤ²` 上の和**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★モジュラー形式の在庫に載せるための段

第 333 で、判別式の非消失は**正規化された束 `Λ_τ = ℤ + τℤ`** に帰着した。
★その先で使いたいのは mathlib のモジュラー Eisenstein 級数
`ModularForm.E`(= `(1/2)·Σ_{gcd(m,n)=1} (mz+n)⁻ᵏ`)であり、
そちらは **`Fin 2 → ℤ` 上の和**で書かれている。

★★本ブロックは格子側の和を `ℤ × ℤ` 上の和に書き換える:

    G(Λ_τ, k) = Σ_{(m,n) ∈ ℤ²} (m + nτ)⁻ᵏ

★★★mathlib の `latticeEquivProd`(束 `≃ₗ[ℤ] ℤ × ℤ`)と
`latticeEquiv_symm_apply`(`(m,n) ↦ mω₁ + nω₂`)がそのまま効く。
★`(m,n) = (0,0)` の項は `(0^k)⁻¹ = 0` なので、和から除く必要はない。

## ★★★次の段

★`Σ_{(m,n) ∈ ℤ²} (m+nτ)⁻ᵏ = 2ζ(k)·E_k(τ)`——gcd で分解する。
★★`ζ(4) = π⁴/90`・`ζ(6) = π⁶/945` は mathlib にある。
★★★そこまで行けば
`g₂³ − 27g₃² = (2π)¹²·Δ(τ)` と `Δ = η²⁴ ≠ 0` で **(i) が閉じる**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `linearIndependent_one_of_im_ne_zero` | ★`im τ ≠ 0` なら `1, τ` は `ℝ` 上一次独立 |
| `normPair` | ★★正規化された周期対 `⟨1, τ⟩` |
| `G_normPair` | ★★★★★**`G(Λ_τ, k) = Σ_{(m,n)} (m+nτ)⁻ᵏ`** |
-/

namespace ABC3.Found.GenEll

open PeriodPair

/-! ## ★正規化された周期対 -/

/-- ★`im τ ≠ 0` ならば `1` と `τ` は `ℝ` 上一次独立である。 -/
theorem linearIndependent_one_of_im_ne_zero (t : ℂ) (h : t.im ≠ 0) :
    LinearIndependent ℝ ![(1 : ℂ), t] := by
  rw [LinearIndependent.pair_iff]
  intro s u hsu
  have h1 : (s : ℂ) + (u : ℂ) * t = 0 := by
    simpa [Complex.real_smul] using hsu
  have him : u * t.im = 0 := by
    have := congrArg Complex.im h1
    simpa using this
  have hu : u = 0 := (mul_eq_zero.1 him).resolve_right h
  subst hu
  have hs : (s : ℂ) = 0 := by simpa using h1
  exact ⟨by exact_mod_cast hs, rfl⟩

/-- ★★**正規化された周期対** `⟨1, τ⟩`。 -/
noncomputable def normPair (t : ℂ) (h : t.im ≠ 0) : PeriodPair :=
  ⟨1, t, linearIndependent_one_of_im_ne_zero t h⟩

/-! ## ★★★★★`ℤ²` 上の和への書き換え -/

/-- ★★★★★**`G(Λ_τ, k) = Σ_{(m,n) ∈ ℤ²} (m + nτ)⁻ᵏ`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`latticeEquivProd`(束 `≃ₗ[ℤ] ℤ × ℤ`)で添字を取り替えるだけ。
★★`(0,0)` の項は `(0^k)⁻¹ = 0` なので除く必要はない。 -/
theorem G_normPair (t : ℂ) (h : t.im ≠ 0) (k : ℕ) :
    (normPair t h).G k = ∑' x : ℤ × ℤ, (((x.1 : ℂ) + (x.2 : ℂ) * t) ^ k)⁻¹ := by
  rw [PeriodPair.G]
  rw [← (normPair t h).latticeEquivProd.toEquiv.symm.tsum_eq
    (fun l : (normPair t h).lattice => ((l : ℂ) ^ k)⁻¹)]
  congr 1
  funext x
  congr 2
  have hx := (normPair t h).latticeEquiv_symm_apply x
  rw [show (((normPair t h).latticeEquivProd.toEquiv.symm x : (normPair t h).lattice) : ℂ)
      = (((normPair t h).latticeEquivProd.symm x : (normPair t h).lattice) : ℂ) from rfl, hx]
  show (x.1 : ℂ) * 1 + (x.2 : ℂ) * t = (x.1 : ℂ) + (x.2 : ℂ) * t
  ring

/-! ## ★出典の紐付け(`.src`) -/

def G_normPair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
