import ABC3.Meta.Claim
import Mathlib.NumberTheory.EllipticDivisibilitySequence
import Mathlib.Tactic.Ring

/-!
# Galois (G5) 第 322 ブロック —— **★★★★★★★楕円列の核心は 2 変数の 3 項恒等式だけ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★到達点

> **`IsEllSequence W` は `r = 1` の場合(2 変数の 3 項恒等式)から `ring` で出る**

★★★これで **mathlib の TODO**(`normEDS` が `IsEllDivSequence` を満たすこと)の
楕円列の側は、**3 変数の恒等式から 2 変数の恒等式へ縮む**。

## ★★★★★★なぜ縮むのか——2 次元の Plücker 関係

`h(x,y) := W(x+y)·W(x−y)` と置くと、3 項恒等式は

    h(m,n) = h(m,1)·h(n,0) − h(n,1)·h(m,0)

すなわち `h(x,y) = det(v_x, v_y)`、`v_x := (W(x+1)W(x−1), W(x)²)` という
**2 次元ベクトルの行列式**である。★★`IsEllSequence` の主張

    h(m,n)·h(r,0) = h(m,r)·h(n,0) − h(n,r)·h(m,0)

は、その行列式表示を代入すれば **2 次元の Plücker 関係**そのものになり、`ring` で閉じる。
★★★同じ理由で **Riemann の 4 項恒等式**(`riemann_four_term`)も出る。

## ★★★★★★残る核心は「差 `3` 以上」だけである

mathlib の 2 本の漸化式は、**3 項恒等式の特別な場合そのもの**である:

| mathlib | 3 項恒等式の場合 | 差 `m − n` |
|---|---|---|
| `normEDS_odd` | `T(m+1, m)` | `1` |
| `normEDS_even` | `T(m+1, m−1)` | `2` |

★さらに `n = 0`・`n = 1` は自明(`normEDS_threeTerm_zero`・`normEDS_threeTerm_one`)、
`m` と `n` の交換は奇関数性から出る(`threeTerm_swap`)。
★★★★★★★**残るのは差が `3` 以上の場合だけ**であり、それが
`Skeleton/GaloisRep/EDSThreeTerm.lean` の `normEDS_isThreeTerm` である。

## ★★これがどこに効くか

(G5) の非退化性は `deg[n] = n²`(第 196・197)に帰着し、それは
`x([n]P) = Φ_n/ΨSq_n`(第 198)に帰着し、それは **EDS 恒等式**に帰着していた
(§9-6xx、`Skeleton/GaloisRep/WeilFunctionField.lean` の測定)。
★本ブロックはその最後の的を **2 変数に縮めた**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `IsThreeTerm` | ★★2 変数の 3 項恒等式(`r = 1` の場合) |
| `isEllSequence_of_isThreeTerm` | ★★★★★★★**`IsEllSequence` は 3 項恒等式から出る** |
| `isThreeTerm_of_isEllSequence` | ★★逆(`W 1 = 1` のとき) |
| `riemann_four_term` | ★★★★**Riemann の 4 項恒等式** |
| `threeTerm_swap` | ★★`m ↔ n` の交換 |
| `normEDS_threeTerm_zero`・`_one` | ★`n = 0, 1` は自明 |
| `normEDS_threeTerm_diff_one` | ★★★★**`normEDS_odd` は差 `1` の場合** |
| `normEDS_threeTerm_diff_two` | ★★★★**`normEDS_even` は差 `2` の場合** |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R]

/-! ## ★★2 変数の 3 項恒等式 -/

/-- ★★**2 変数の 3 項恒等式**——`IsEllSequence` の `r = 1` の場合。

    W(m+n)·W(m−n) = W(m+1)W(m−1)·W(n)² − W(n+1)W(n−1)·W(m)²

★`h(x,y) := W(x+y)W(x−y)` と置けば `h(x,y) = det((W(x+1)W(x−1), W(x)²), (…))` である。 -/
def IsThreeTerm (W : ℤ → R) : Prop :=
  ∀ m n : ℤ, W (m + n) * W (m - n)
    = W (m + 1) * W (m - 1) * W n ^ 2 - W (n + 1) * W (n - 1) * W m ^ 2

/-- ★★★★★★★**`IsEllSequence` は 3 項恒等式から出る**——2 次元の Plücker 関係。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★3 変数の恒等式が 2 変数の恒等式に縮む。`W 1 = 1` すら要らない。 -/
theorem isEllSequence_of_isThreeTerm {W : ℤ → R} (hT : IsThreeTerm W) : IsEllSequence W := by
  intro m n r
  rw [hT m n, hT m r, hT n r]
  ring

/-- ★★逆——`W 1 = 1` のとき `IsEllSequence` は 3 項恒等式を含む。 -/
theorem isThreeTerm_of_isEllSequence {W : ℤ → R} (h1 : W 1 = 1) (hE : IsEllSequence W) :
    IsThreeTerm W := by
  intro m n
  have h := hE m n 1
  rw [h1] at h
  simpa using h

/-- ★★★★**Riemann の 4 項恒等式**——3 項恒等式から `ring` で出る。

★`h(x,y)` が 2 次元ベクトルの行列式であることの言い換えである。 -/
theorem riemann_four_term {W : ℤ → R} (hT : IsThreeTerm W) (p q r s : ℤ) :
    W (p + q) * W (p - q) * (W (r + s) * W (r - s))
      - W (p + r) * W (p - r) * (W (q + s) * W (q - s))
      + W (p + s) * W (p - s) * (W (q + r) * W (q - r)) = 0 := by
  rw [hT p q, hT r s, hT p r, hT q s, hT p s, hT q r]
  ring

/-- ★★`m` と `n` の交換——奇関数性から。 -/
theorem threeTerm_swap {W : ℤ → R} (hodd : ∀ n : ℤ, W (-n) = -W n) {m n : ℤ}
    (h : W (m + n) * W (m - n) = W (m + 1) * W (m - 1) * W n ^ 2
      - W (n + 1) * W (n - 1) * W m ^ 2) :
    W (n + m) * W (n - m) = W (n + 1) * W (n - 1) * W m ^ 2
      - W (m + 1) * W (m - 1) * W n ^ 2 := by
  have h1 : n - m = -(m - n) := by ring
  have h2 : n + m = m + n := by ring
  rw [h1, h2, hodd]
  linear_combination -h

/-! ## ★★★★`normEDS` について既に分かっている場合 -/

variable (b c d : R)

/-- ★`n = 0` は自明。 -/
theorem normEDS_threeTerm_zero (m : ℤ) :
    normEDS b c d (m + 0) * normEDS b c d (m - 0)
      = normEDS b c d (m + 1) * normEDS b c d (m - 1) * normEDS b c d 0 ^ 2
        - normEDS b c d (0 + 1) * normEDS b c d (0 - 1) * normEDS b c d m ^ 2 := by
  have h1 : (0 : ℤ) + 1 = 1 := by ring
  have h2 : (0 : ℤ) - 1 = -1 := by ring
  rw [add_zero, sub_zero, h1, h2, normEDS_zero, normEDS_one,
    show ((-1 : ℤ)) = -(1 : ℤ) from rfl, normEDS_neg, normEDS_one]
  ring

/-- ★`n = 1` は自明。 -/
theorem normEDS_threeTerm_one (m : ℤ) :
    normEDS b c d (m + 1) * normEDS b c d (m - 1)
      = normEDS b c d (m + 1) * normEDS b c d (m - 1) * normEDS b c d 1 ^ 2
        - normEDS b c d (1 + 1) * normEDS b c d (1 - 1) * normEDS b c d m ^ 2 := by
  have h1 : (1 : ℤ) - 1 = 0 := by ring
  rw [h1, normEDS_zero, normEDS_one]
  ring

/-- ★★★★**`normEDS_odd` は 3 項恒等式の差 `1` の場合そのものである**。 -/
theorem normEDS_threeTerm_diff_one (m : ℤ) :
    normEDS b c d (m + 1 + m) * normEDS b c d (m + 1 - m)
      = normEDS b c d (m + 1 + 1) * normEDS b c d (m + 1 - 1) * normEDS b c d m ^ 2
        - normEDS b c d (m + 1) * normEDS b c d (m - 1) * normEDS b c d (m + 1) ^ 2 := by
  have h : m + 1 + m = 2 * m + 1 := by ring
  have h2 : m + 1 - m = 1 := by ring
  have h3 : m + 1 + 1 = m + 2 := by ring
  have h4 : m + 1 - 1 = m := by ring
  rw [h, h2, h3, h4, normEDS_one, normEDS_odd]
  ring

/-- ★★★★**`normEDS_even` は 3 項恒等式の差 `2` の場合そのものである**。 -/
theorem normEDS_threeTerm_diff_two (m : ℤ) :
    normEDS b c d (m + 1 + (m - 1)) * normEDS b c d (m + 1 - (m - 1))
      = normEDS b c d (m + 1 + 1) * normEDS b c d (m + 1 - 1) * normEDS b c d (m - 1) ^ 2
        - normEDS b c d (m - 1 + 1) * normEDS b c d (m - 1 - 1) * normEDS b c d (m + 1) ^ 2 := by
  have h : m + 1 + (m - 1) = 2 * m := by ring
  have h2 : m + 1 - (m - 1) = 2 := by ring
  have h3 : m + 1 + 1 = m + 2 := by ring
  have h4 : m + 1 - 1 = m := by ring
  have h5 : m - 1 + 1 = m := by ring
  have h6 : m - 1 - 1 = m - 2 := by ring
  rw [h, h2, h3, h4, h5, h6, normEDS_two, normEDS_even]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def IsThreeTerm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(楕円列の 3 項恒等式)",
    sectionId := "genell-thm-3-8" }

def isEllSequence_of_isThreeTerm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(楕円列は 3 項恒等式から出る)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
