import ABC3.Found.GaloisRep.EdsIdentity

/-!
# Galois (G1) 第 58 ブロック —— **★★★★★★★★Ward の定理**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★mathlib の TODO を閉じる

`Mathlib/NumberTheory/EllipticDivisibilitySequence.lean` の

> TODO: prove that `normEDS` satisfies `IsEllDivSequence`.

のうち、**`IsEllSequence` の側**が本ブロックで閉じる:

    W(m+n)W(m−n)W(r)² = W(m+r)W(m−r)W(n)² − W(n+r)W(n−r)W(m)²

## ★★★★★★機構——積の並べ替え

`r = 1` の場合 `EllAt m n` を **`n` についての 2 歩帰納**で示す。要は

    (S₁D₁)·(S₋₁D₋₁) = (S₁D₋₁)·(S₋₁D₁)     (S_j = W(m+n+j), D_j = W(m−n+j))

という**積の並べ替え**である。★左辺は帰納の仮定 `EllAt (m±1) n` で値が分かり、
右辺の `S₋₁D₁` は `EllAt m (n−1)` で分かるので、`S₁D₋₁` が決まる。

★★そのとき残る差が **Somos-4(第 47)と (Λ)(第 56)** でちょうど消える
——係数は `c·a₀²b₀⁴`、`−c·a₀⁴b₀²`、`−a₀b₀²b₁b₋₁`、`a₀²a₁a₋₁b₀`(§9-396 で実測)。

## ★★★★割り算は普遍環で行う

帰納段は `W(m+n−1)W(m−n+1)` で割る。★**普遍環 `ℤ[b,c,d]` は整域**で、
第 57 ブロックより `W(j) ≠ 0` (`j ≠ 0`) なので割れる。
★★例外(`m = n−1`、`m = 1−n`)は **`normEDS_even` そのもの**であった。

## ★★★★★一般の `r` は `r = 1` から出る

    W(r)²·EllAt(m,n) − W(n)²·EllAt(m,r) + W(m)²·EllAt(n,r)

を取ると、余分な項がすべて消える。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `normEDS_lam` | ★(Λ) の `normEDS` 版 |
| `EllAt` | ★Ward の関係(`r = 1`) |
| `ell_prod` / `ell_prod'` | ★★★★積の並べ替えで出る恒等式 |
| `ell_step` | ★★★★★★帰納段 |
| `ell` | ★★★★★★★`EllAt`——任意の可換環 |
| `isEllSequence_normEDS` | ★★★★★★★★**Ward の定理** |
-/

namespace ABC3.Found.GaloisRep

universe u

variable {R : Type u} [CommRing R] (b c d : R)

/-- ★(Λ) の `normEDS` 版。 -/
theorem normEDS_lam (n : ℤ) :
    c * (b ^ 2 * normEDS b c d n ^ 3
        + normEDS b c d (n + 1) ^ 2 * normEDS b c d (n - 2)
        + normEDS b c d (n + 2) * normEDS b c d (n - 1) ^ 2)
      = (b ^ 4 + d) * normEDS b c d (n - 1) * normEDS b c d n * normEDS b c d (n + 1) := by
  have h := lam (b ^ 4) c d n
  simp only [LamAt] at h
  simp only [normEDS]
  have e2 : Even (n + 2) ↔ Even n := by simp
  have em2 : Even (n - 2) ↔ Even n := by simp
  have e1 : Even (n + 1) ↔ ¬ Even n := Int.even_add_one
  have em1 : Even (n - 1) ↔ ¬ Even n := Int.even_sub_one
  by_cases hn : Even n
  · rw [if_pos hn, if_pos (e2.mpr hn), if_pos (em2.mpr hn),
      if_neg (fun hh => (e1.mp hh) hn), if_neg (fun hh => (em1.mp hh) hn)] at *
    linear_combination b * h
  · rw [if_neg hn, if_neg (fun hh => hn (e2.mp hh)), if_neg (fun hh => hn (em2.mp hh)),
      if_pos (e1.mpr hn), if_pos (em1.mpr hn)] at *
    linear_combination b ^ 2 * h

/-- ★**Ward の関係(`r = 1` の場合)**。 -/
def EllAt (b c d : R) (m n : ℤ) : Prop :=
  normEDS b c d (m + n) * normEDS b c d (m - n)
    = normEDS b c d (m + 1) * normEDS b c d (m - 1) * normEDS b c d n ^ 2
      - normEDS b c d (n + 1) * normEDS b c d (n - 1) * normEDS b c d m ^ 2

theorem ell_zero (m : ℤ) : EllAt b c d m 0 := by
  simp only [EllAt]; norm_num; ring

theorem ell_one (m : ℤ) : EllAt b c d m 1 := by
  simp only [EllAt]; norm_num

/-- ★★★★**積の並べ替えで出る恒等式**——Somos-4 と (Λ) だけで閉じる。 -/
theorem ell_prod (m n : ℤ) :
    (normEDS b c d (m + 2) * normEDS b c d m * normEDS b c d n ^ 2
        - normEDS b c d (n + 1) * normEDS b c d (n - 1) * normEDS b c d (m + 1) ^ 2)
      * (normEDS b c d m * normEDS b c d (m - 2) * normEDS b c d n ^ 2
        - normEDS b c d (n + 1) * normEDS b c d (n - 1) * normEDS b c d (m - 1) ^ 2)
      * c
    = (normEDS b c d (m + 1) * normEDS b c d (m - 1) * normEDS b c d (n + 1) ^ 2
        - normEDS b c d (n + 2) * normEDS b c d n * normEDS b c d m ^ 2)
      * (normEDS b c d (m + 1) * normEDS b c d (m - 1) * normEDS b c d (n - 1) ^ 2
        - normEDS b c d n * normEDS b c d (n - 2) * normEDS b c d m ^ 2)
      * c := by
  have SA := normEDS_somos b c d m
  have SB := normEDS_somos b c d n
  have LA := normEDS_lam b c d m
  have LB := normEDS_lam b c d n
  set a0 := normEDS b c d m
  set a1 := normEDS b c d (m + 1)
  set am1 := normEDS b c d (m - 1)
  set b0 := normEDS b c d n
  set b1 := normEDS b c d (n + 1)
  set bm1 := normEDS b c d (n - 1)
  linear_combination (c * a0^2 * b0^4) * SA - (c * a0^4 * b0^2) * SB
    - (a0 * b0^2 * b1 * bm1) * LA + (a0^2 * a1 * am1 * b0) * LB

theorem ell_prod' [NoZeroDivisors R] (hc : c ≠ 0) (m n : ℤ) :
    (normEDS b c d (m + 2) * normEDS b c d m * normEDS b c d n ^ 2
        - normEDS b c d (n + 1) * normEDS b c d (n - 1) * normEDS b c d (m + 1) ^ 2)
      * (normEDS b c d m * normEDS b c d (m - 2) * normEDS b c d n ^ 2
        - normEDS b c d (n + 1) * normEDS b c d (n - 1) * normEDS b c d (m - 1) ^ 2)
    = (normEDS b c d (m + 1) * normEDS b c d (m - 1) * normEDS b c d (n + 1) ^ 2
        - normEDS b c d (n + 2) * normEDS b c d n * normEDS b c d m ^ 2)
      * (normEDS b c d (m + 1) * normEDS b c d (m - 1) * normEDS b c d (n - 1) ^ 2
        - normEDS b c d n * normEDS b c d (n - 2) * normEDS b c d m ^ 2) :=
  mul_right_cancel₀ hc (ell_prod b c d m n)

/-- ★★例外の場合(`m = n−1`)は `normEDS_even` そのもの。 -/
theorem ell_exc (n : ℤ) : EllAt b c d (n - 1) (n + 1) := by
  simp only [EllAt]
  have E := normEDS_even b c d n
  rw [show n - 1 + (n + 1) = 2 * n by ring, show n - 1 - (n + 1) = -2 by ring,
    show n - 1 + 1 = n by ring, show n - 1 - 1 = n - 2 by ring,
    show n + 1 + 1 = n + 2 by ring, show n + 1 - 1 = n by ring]
  rw [show (-2 : ℤ) = -(2 : ℤ) from rfl, normEDS_neg, normEDS_two]
  linear_combination -E

/-- ★★例外の場合(`m = 1−n`)。 -/
theorem ell_exc' (n : ℤ) : EllAt b c d (1 - n) (n + 1) := by
  simp only [EllAt]
  have E := normEDS_even b c d n
  rw [show 1 - n + (n + 1) = 2 by ring, show 1 - n - (n + 1) = -(2 * n) by ring,
    show 1 - n + 1 = -(n - 2) by ring, show 1 - n - 1 = -n by ring,
    show (1 : ℤ) - n = -(n - 1) by ring, show n + 1 + 1 = n + 2 by ring,
    show n + 1 - 1 = n by ring]
  rw [normEDS_neg, normEDS_neg, normEDS_neg, normEDS_neg, normEDS_two]
  linear_combination -E

/-- ★★★★★★**帰納段**——`EllAt (m±1) n` と `EllAt m (n−1)` から `EllAt m (n+1)`。 -/
theorem ell_step [NoZeroDivisors R] (hc : c ≠ 0)
    (hW : ∀ j : ℤ, j ≠ 0 → normEDS b c d j ≠ 0) (m n : ℤ)
    (h1 : EllAt b c d (m + 1) n) (h2 : EllAt b c d (m - 1) n) (h3 : EllAt b c d m (n - 1)) :
    EllAt b c d m (n + 1) := by
  by_cases he1 : m - n + 1 = 0
  · have hm : m = n - 1 := by omega
    subst hm; exact ell_exc b c d n
  by_cases he2 : m + n - 1 = 0
  · have hm : m = 1 - n := by omega
    subst hm; exact ell_exc' b c d n
  simp only [EllAt] at h1 h2 h3 ⊢
  rw [show m + 1 + n = m + n + 1 by ring, show m + 1 - n = m - n + 1 by ring,
    show m + 1 + 1 = m + 2 by ring, show m + 1 - 1 = m by ring] at h1
  rw [show m - 1 + n = m + n - 1 by ring, show m - 1 - n = m - n - 1 by ring,
    show m - 1 + 1 = m by ring, show m - 1 - 1 = m - 2 by ring] at h2
  rw [show m + (n - 1) = m + n - 1 by ring, show m - (n - 1) = m - n + 1 by ring,
    show n - 1 + 1 = n by ring, show n - 1 - 1 = n - 2 by ring] at h3
  rw [show m + (n + 1) = m + n + 1 by ring, show m - (n + 1) = m - n - 1 by ring,
    show n + 1 + 1 = n + 2 by ring, show n + 1 - 1 = n by ring]
  have hRne : normEDS b c d (m + n - 1) * normEDS b c d (m - n + 1) ≠ 0 :=
    mul_ne_zero (hW _ he2) (hW _ he1)
  refine mul_right_cancel₀ hRne ?_
  have hP := ell_prod' b c d hc m n
  linear_combination (normEDS b c d (m + n - 1) * normEDS b c d (m - n - 1)) * h1
    + (normEDS b c d (m + 2) * normEDS b c d m * normEDS b c d n ^ 2
        - normEDS b c d (n + 1) * normEDS b c d (n - 1) * normEDS b c d (m + 1) ^ 2) * h2
    - (normEDS b c d (m + 1) * normEDS b c d (m - 1) * normEDS b c d (n + 1) ^ 2
        - normEDS b c d (n + 2) * normEDS b c d n * normEDS b c d m ^ 2) * h3
    + hP

theorem ell_nat [NoZeroDivisors R] (hc : c ≠ 0)
    (hW : ∀ j : ℤ, j ≠ 0 → normEDS b c d j ≠ 0) :
    ∀ (n : ℕ) (m : ℤ), EllAt b c d m (n : ℤ) := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    match n with
    | 0 => exact fun m => by simpa using ell_zero b c d m
    | 1 => exact fun m => by simpa using ell_one b c d m
    | (k + 2) =>
        intro m
        have hk1 : ∀ m', EllAt b c d m' ((k + 1 : ℕ) : ℤ) := ih (k + 1) (by omega)
        have hk0 : ∀ m', EllAt b c d m' ((k : ℕ) : ℤ) := ih k (by omega)
        have key := ell_step b c d hc hW m ((k + 1 : ℕ) : ℤ)
          (hk1 (m + 1)) (hk1 (m - 1))
          (by rw [show ((k + 1 : ℕ) : ℤ) - 1 = ((k : ℕ) : ℤ) by push_cast; ring]; exact hk0 m)
        rw [show ((k + 2 : ℕ) : ℤ) = ((k + 1 : ℕ) : ℤ) + 1 by push_cast; ring]
        exact key

theorem ell_int [NoZeroDivisors R] (hc : c ≠ 0)
    (hW : ∀ j : ℤ, j ≠ 0 → normEDS b c d j ≠ 0) (m n : ℤ) : EllAt b c d m n := by
  rcases le_or_gt 0 n with h | h
  · lift n to ℕ using h; exact ell_nat b c d hc hW n m
  · have key := ell_nat b c d hc hW (-n).toNat m
    rw [Int.toNat_of_nonneg (by omega)] at key
    simp only [EllAt] at key ⊢
    rw [show m + -n = m - n by ring, show m - -n = m + n by ring,
      show -n + 1 = -(n - 1) by ring, show -n - 1 = -(n + 1) by ring] at key
    simp only [normEDS_neg] at key
    linear_combination key

/-- ★★★★★★★`EllAt`——任意の可換環(普遍環から送る)。 -/
theorem ell (m n : ℤ) : EllAt b c d m n := by
  have hu := ell_int (MvPolynomial.X 0 : MvPolynomial (Fin 3) ℤ) (MvPolynomial.X 1)
    (MvPolynomial.X 2) (MvPolynomial.X_ne_zero 1) normEDS_univ_ne_zero m n
  simp only [EllAt] at hu ⊢
  let f : MvPolynomial (Fin 3) ℤ →+* R := MvPolynomial.eval₂Hom (Int.castRingHom R) ![b, c, d]
  have hX0 : f (MvPolynomial.X 0) = b := by simp [f]
  have hX1 : f (MvPolynomial.X 1) = c := by simp [f]
  have hX2 : f (MvPolynomial.X 2) = d := by simp [f]
  have h := congrArg f hu
  simp only [map_mul, map_sub, map_pow, map_normEDS, hX0, hX1, hX2] at h
  exact h

/-- ★★★★★★★★**Ward の定理**——`normEDS` は楕円列である。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これは mathlib の
`Mathlib/NumberTheory/EllipticDivisibilitySequence.lean` の TODO
「prove that `normEDS` satisfies `IsEllDivSequence`」の `IsEllSequence` の側である。 -/
theorem isEllSequence_normEDS : IsEllSequence (normEDS b c d) := by
  intro m n r
  have hmn := ell b c d m n
  have hmr := ell b c d m r
  have hnr := ell b c d n r
  simp only [EllAt] at hmn hmr hnr
  linear_combination (normEDS b c d r ^ 2) * hmn - (normEDS b c d n ^ 2) * hmr
    + (normEDS b c d m ^ 2) * hnr

/-! ## ★出典の紐付け(`.src`) -/

def isEllSequence_normEDS.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Ward の定理——normEDS は楕円列)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
