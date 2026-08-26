import ABC3.Found.GenEll.LatticePoint

/-!
# GenEll 第 337 ブロック —— **★★★★★格子点の gcd 分解**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★一意化の (i) に向けた代数の芯

(i)「束の判別式 `g₂³ − 27g₃² ≠ 0`」は、第 333・334 により

    G(Λ_τ, k) = Σ_{(m,n) ∈ ℤ²} (m + nτ)⁻ᵏ

を mathlib のモジュラー Eisenstein 級数 `E_k`(`= (1/2)Σ_{gcd(m,n)=1} (mz+n)⁻ᵏ`)に
繋ぐことに帰着した。★その繋ぎは古典的な **gcd 分解**

    Σ_{v ≠ 0} f(v) = Σ_{d ≥ 1} Σ_{w 原始} f(d·w) = Σ_{d ≥ 1} d⁻ᵏ · Σ_{w 原始} f(w)
                   = ζ(k) · Σ_{w 原始} f(w) = 2ζ(k)·E_k

である。★★本ブロックはその**代数の芯**——分解の存在と一意性、および被加数の分離——を取る
(和の入れ替え自体は解析の段で、別ブロック)。

## ★★★★★中身

* `exists_unique_gcd_decomp`: 零でない `v ∈ ℤ²` は `v = d·w`(`d > 0`、`w` 原始)と
  **一意に**分解する。★一意性は `gcd(d·w₁, d·w₂) = d·gcd(w₁,w₂) = d` から出る。
* `eisSummand_gcd_factor`: 被加数は `((d·w)₁ + (d·w)₂τ)⁻ᵏ = d⁻ᵏ·(w₁ + w₂τ)⁻ᵏ` と**分離する**。

## ★★残り

★和の入れ替え(`tsum` の分割と絶対収束)と、mathlib の `E_k` の正規化との突き合わせ。
★★`ζ(4) = π⁴/90`・`ζ(6) = π⁶/945` は mathlib にある。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_unique_gcd_decomp` | ★★★★**gcd 分解の存在と一意性** |
| `eisSummand_gcd_factor` | ★★★被加数の分離 |
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★gcd 分解 -/

/-- ★★★★**零でない格子点は「gcd × 原始ベクトル」に一意に分解する**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★一意性は `gcd(d·w₁, d·w₂) = d·gcd(w₁,w₂) = d` から出る。 -/
theorem exists_unique_gcd_decomp (v : ℤ × ℤ) (hv : v ≠ (0, 0)) :
    ∃! p : ℕ × (ℤ × ℤ),
      0 < p.1 ∧ Int.gcd p.2.1 p.2.2 = 1 ∧ v = ((p.1 : ℤ) * p.2.1, (p.1 : ℤ) * p.2.2) := by
  have hg : 0 < Int.gcd v.1 v.2 := by
    rcases Nat.eq_zero_or_pos (Int.gcd v.1 v.2) with h | h
    · exfalso
      rw [Int.gcd_eq_zero_iff] at h
      exact hv (Prod.ext h.1 h.2)
    · exact h
  set d : ℕ := Int.gcd v.1 v.2 with hd
  have hdne : (d : ℤ) ≠ 0 := by exact_mod_cast hg.ne'
  have hd1 : (d : ℤ) ∣ v.1 := Int.gcd_dvd_left v.1 v.2
  have hd2 : (d : ℤ) ∣ v.2 := Int.gcd_dvd_right v.1 v.2
  refine ⟨⟨d, (v.1 / d, v.2 / d)⟩, ⟨hg, Int.gcd_div_gcd_div_gcd hg, ?_⟩, ?_⟩
  · refine Prod.ext ?_ ?_
    · exact (Int.ediv_mul_cancel hd1).symm.trans (mul_comm _ _)
    · exact (Int.ediv_mul_cancel hd2).symm.trans (mul_comm _ _)
  · rintro ⟨e, w⟩ ⟨he, hw, hveq⟩
    have hv1 : v.1 = (e : ℤ) * w.1 := congrArg Prod.fst hveq
    have hv2 : v.2 = (e : ℤ) * w.2 := congrArg Prod.snd hveq
    have hde : d = e := by
      rw [hd, hv1, hv2, Int.gcd_mul_left, hw, mul_one, Int.natAbs_natCast]
    subst hde
    refine Prod.ext rfl (Prod.ext ?_ ?_)
    · show w.1 = v.1 / (d : ℤ)
      rw [hv1, Int.mul_ediv_cancel_left _ hdne]
    · show w.2 = v.2 / (d : ℤ)
      rw [hv2, Int.mul_ediv_cancel_left _ hdne]

/-! ## ★★★被加数の分離 -/

/-- ★★★**被加数は `d⁻ᵏ` と原始部分に分離する**。 -/
theorem eisSummand_gcd_factor (t : ℂ) (d : ℕ) (w : ℤ × ℤ) (k : ℕ) :
    (((((d : ℤ) * w.1 : ℤ) : ℂ) + (((d : ℤ) * w.2 : ℤ) : ℂ) * t) ^ k)⁻¹
      = ((d : ℂ) ^ k)⁻¹ * (((w.1 : ℂ) + (w.2 : ℂ) * t) ^ k)⁻¹ := by
  push_cast
  rw [show (d : ℂ) * (w.1 : ℂ) + (d : ℂ) * (w.2 : ℂ) * t
      = (d : ℂ) * ((w.1 : ℂ) + (w.2 : ℂ) * t) by ring, mul_pow, mul_inv]

/-! ## ★出典の紐付け(`.src`) -/

def exists_unique_gcd_decomp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def eisSummand_gcd_factor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
