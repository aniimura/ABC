import ABC3.Found.GaloisRep.EdsLambda

/-!
# Galois (G1) 第 57 ブロック —— **★★★★★恒等列で普遍環の非零性を出す**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★普遍環 `ℤ[b,c,d]` で `normEDS ≠ 0` を言う道

Ward の定理(第 58 ブロック)の帰納段では

    W(m+n−1)·W(m−n+1) ≠ 0

で割る必要がある。★普遍環 `ℤ[b,c,d]` は整域なので、
**`normEDS X₀ X₁ X₂ j ≠ 0` (`j ≠ 0`)** さえ言えればよい。

## ★★★★★★恒等列へ送るだけでよい

    normEDS 2 3 2 = id

★これは `normEDSRec` で示せる(実際に計算すると)

    偶数段: (m−1)²m(m+2) − (m−2)m(m+1)² = 4m ✅
    奇数段: (m+2)m³ − (m−1)(m+1)³ = 2m+1 ✅

★★したがって `X₀ ↦ 2, X₁ ↦ 3, X₂ ↦ 2` で送れば `normEDS X₀ X₁ X₂ j ↦ j ≠ 0` ✅

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `normEDS_id_nat` / `normEDS_id` | ★★★**`normEDS 2 3 2 = id`** |
| `normEDS_univ_ne_zero` | ★★★★★**普遍環での非零性** |
-/

namespace ABC3.Found.GaloisRep

/-- ★★★**恒等列**——`normEDS 2 3 2 n = n`(自然数版)。 -/
theorem normEDS_id_nat (n : ℕ) : normEDS (2 : ℤ) 3 2 (n : ℤ) = (n : ℤ) := by
  induction n using normEDSRec with
  | zero => simp
  | one => simp
  | two => norm_num
  | three => norm_num
  | four => norm_num
  | even m h1 h2 h3 h4 h5 =>
      have E := normEDS_even (2 : ℤ) 3 2 ((m : ℤ) + 3)
      rw [show (m : ℤ) + 3 - 1 = ((m + 2 : ℕ) : ℤ) by push_cast; ring,
        show (m : ℤ) + 3 = ((m + 3 : ℕ) : ℤ) by push_cast; ring,
        show ((m + 3 : ℕ) : ℤ) + 2 = ((m + 5 : ℕ) : ℤ) by push_cast; ring,
        show ((m + 3 : ℕ) : ℤ) - 2 = ((m + 1 : ℕ) : ℤ) by push_cast; ring,
        show ((m + 3 : ℕ) : ℤ) + 1 = ((m + 4 : ℕ) : ℤ) by push_cast; ring,
        h1, h2, h3, h4, h5] at E
      rw [show ((2 * (m + 3) : ℕ) : ℤ) = 2 * ((m + 3 : ℕ) : ℤ) by push_cast; ring]
      refine mul_right_cancel₀ (by norm_num : (2 : ℤ) ≠ 0) ?_
      rw [E]; push_cast; ring
  | odd m h1 h2 h3 h4 =>
      have O := normEDS_odd (2 : ℤ) 3 2 ((m : ℤ) + 2)
      rw [show (m : ℤ) + 2 + 2 = ((m + 4 : ℕ) : ℤ) by push_cast; ring,
        show (m : ℤ) + 2 = ((m + 2 : ℕ) : ℤ) by push_cast; ring,
        show ((m + 2 : ℕ) : ℤ) - 1 = ((m + 1 : ℕ) : ℤ) by push_cast; ring,
        show ((m + 2 : ℕ) : ℤ) + 1 = ((m + 3 : ℕ) : ℤ) by push_cast; ring,
        h1, h2, h3, h4] at O
      rw [show ((2 * (m + 2) + 1 : ℕ) : ℤ) = 2 * ((m + 2 : ℕ) : ℤ) + 1 by push_cast; ring, O]
      push_cast; ring

/-- ★★★**恒等列**——`normEDS 2 3 2 n = n`。 -/
theorem normEDS_id (n : ℤ) : normEDS (2 : ℤ) 3 2 n = n := by
  rcases le_or_gt 0 n with h | h
  · lift n to ℕ using h; exact normEDS_id_nat n
  · have key := normEDS_id_nat (-n).toNat
    rw [Int.toNat_of_nonneg (by omega)] at key
    have hneg := normEDS_neg (2 : ℤ) 3 2 n
    rw [key] at hneg
    omega

/-- ★★★★★**普遍環では `normEDS X₀ X₁ X₂ n ≠ 0`**(`n ≠ 0`)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★恒等列 `normEDS 2 3 2 = id` へ送れば `n ≠ 0` になるからである。 -/
theorem normEDS_univ_ne_zero (n : ℤ) (hn : n ≠ 0) :
    normEDS (MvPolynomial.X 0 : MvPolynomial (Fin 3) ℤ)
      (MvPolynomial.X 1) (MvPolynomial.X 2) n ≠ 0 := by
  intro h
  let f : MvPolynomial (Fin 3) ℤ →+* ℤ := MvPolynomial.eval₂Hom (RingHom.id ℤ) ![2, 3, 2]
  have hX0 : f (MvPolynomial.X 0) = 2 := by simp [f]
  have hX1 : f (MvPolynomial.X 1) = 3 := by simp [f]
  have hX2 : f (MvPolynomial.X 2) = 2 := by simp [f]
  have hc := congrArg f h
  rw [map_normEDS, hX0, hX1, hX2, normEDS_id, map_zero] at hc
  exact hn hc

/-! ## ★出典の紐付け(`.src`) -/

def normEDS_univ_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(普遍環での分点列の非零性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
