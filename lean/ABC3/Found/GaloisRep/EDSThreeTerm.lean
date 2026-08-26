import ABC3.Found.GaloisRep.EdsWard

/-!
# Galois (G5) 第 322 ブロック —— **★★★3 項恒等式の抽象化と Riemann の 4 項恒等式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★2026-08-26 の訂正——**本体は第 47-58 ブロックで既に済んでいた**

★当初このファイルは「`normEDS` が楕円列であること」を未解決の葉として立てたが、
**それは第 58 ブロック `EdsWard.lean` の `isEllSequence_normEDS` で閉じていた**
(Somos-4 は第 47、(Λ) は第 56、`W(j) ≠ 0` は第 57)。
★★在庫(`EdsSomos`・`EdsLambda`・`EdsIdentity`・`EdsWard`)を引かずに書いたための重複であり、
立てた Skeleton 節点は撤去した。

## ★★本ファイルに残す中身

`EdsWard` が `normEDS` について具体的に行っていることを、**任意の列 `W : ℤ → R`** に
対する言明として抜き出したものだけを置く:

* `IsThreeTerm W` —— `EllAt` の全称版(`r = 1` の場合)
* `isEllSequence_of_isThreeTerm` —— ★★★**`IsEllSequence` は `r = 1` の場合から `ring` で出る**
* `isThreeTerm_of_isEllSequence` —— 逆(`W 1 = 1` のとき)
* `riemann_four_term` —— ★★★**Riemann の 4 項恒等式**(`EdsWard` には無い)

## ★★★なぜ `r = 1` で足りるのか——2 次元の Plücker 関係

`h(x,y) := W(x+y)·W(x−y)` と置くと、3 項恒等式は

    h(m,n) = h(m,1)·h(n,0) − h(n,1)·h(m,0)

すなわち `h(x,y) = det(v_x, v_y)`、`v_x := (W(x+1)W(x−1), W(x)²)` という
**2 次元ベクトルの行列式**である。★★`IsEllSequence` の主張

    h(m,n)·h(r,0) = h(m,r)·h(n,0) − h(n,r)·h(m,0)

は、その行列式表示を代入すれば **2 次元の Plücker 関係**そのものになる。
★★★同じ理由で 4 つの添字に対する **Riemann の 4 項恒等式**も出る。
★`W 1 = 1` すら要らない。
-/

namespace ABC3.Found.GaloisRep

variable {A : Type} [CommRing A]

/-! ## ★★2 変数の 3 項恒等式 -/

/-- ★★**2 変数の 3 項恒等式**——`IsEllSequence` の `r = 1` の場合。

    W(m+n)·W(m−n) = W(m+1)W(m−1)·W(n)² − W(n+1)W(n−1)·W(m)²

★`normEDS` についてはこれが `EdsWard.lean` の `EllAt`(第 58 で証明済み)である。 -/
def IsThreeTerm (W : ℤ → A) : Prop :=
  ∀ m n : ℤ, W (m + n) * W (m - n)
    = W (m + 1) * W (m - 1) * W n ^ 2 - W (n + 1) * W (n - 1) * W m ^ 2

/-- ★★★**`IsEllSequence` は 3 項恒等式から `ring` で出る**——2 次元の Plücker 関係。

★★3 変数の恒等式が 2 変数の恒等式に縮む。`W 1 = 1` すら要らない。 -/
theorem isEllSequence_of_isThreeTerm {W : ℤ → A} (hT : IsThreeTerm W) : IsEllSequence W := by
  intro m n r
  rw [hT m n, hT m r, hT n r]
  ring

/-- ★★逆——`W 1 = 1` のとき `IsEllSequence` は 3 項恒等式を含む。 -/
theorem isThreeTerm_of_isEllSequence {W : ℤ → A} (h1 : W 1 = 1) (hE : IsEllSequence W) :
    IsThreeTerm W := by
  intro m n
  have h := hE m n 1
  rw [h1] at h
  simpa using h

/-- ★★★**Riemann の 4 項恒等式**——3 項恒等式から `ring` で出る。

★`h(x,y) := W(x+y)W(x−y)` が 2 次元ベクトルの行列式であることの言い換えである。 -/
theorem riemann_four_term {W : ℤ → A} (hT : IsThreeTerm W) (p q r s : ℤ) :
    W (p + q) * W (p - q) * (W (r + s) * W (r - s))
      - W (p + r) * W (p - r) * (W (q + s) * W (q - s))
      + W (p + s) * W (p - s) * (W (q + r) * W (q - r)) = 0 := by
  rw [hT p q, hT r s, hT p r, hT q s, hT p s, hT q r]
  ring

/-- ★★`m` と `n` の交換——奇関数性から。 -/
theorem threeTerm_swap {W : ℤ → A} (hodd : ∀ n : ℤ, W (-n) = -W n) {m n : ℤ}
    (h : W (m + n) * W (m - n) = W (m + 1) * W (m - 1) * W n ^ 2
      - W (n + 1) * W (n - 1) * W m ^ 2) :
    W (n + m) * W (n - m) = W (n + 1) * W (n - 1) * W m ^ 2
      - W (m + 1) * W (m - 1) * W n ^ 2 := by
  rw [show n - m = -(m - n) by ring, show n + m = m + n by ring, hodd]
  linear_combination -h

/-! ## ★★`normEDS` はこれを満たす(第 58 ブロックの言い換え) -/

/-- ★★**`normEDS` は 3 項恒等式を満たす**——第 58 の `ell` そのもの。 -/
theorem normEDS_isThreeTerm (b c d : A) : IsThreeTerm (normEDS b c d) :=
  fun m n => ell b c d m n

/-- ★★`normEDS` についての Riemann の 4 項恒等式。 -/
theorem normEDS_riemann_four_term (b c d : A) (p q r s : ℤ) :
    normEDS b c d (p + q) * normEDS b c d (p - q)
        * (normEDS b c d (r + s) * normEDS b c d (r - s))
      - normEDS b c d (p + r) * normEDS b c d (p - r)
        * (normEDS b c d (q + s) * normEDS b c d (q - s))
      + normEDS b c d (p + s) * normEDS b c d (p - s)
        * (normEDS b c d (q + r) * normEDS b c d (q - r)) = 0 :=
  riemann_four_term (normEDS_isThreeTerm b c d) p q r s

/-! ## ★出典の紐付け(`.src`) -/

def IsThreeTerm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(楕円列の 3 項恒等式)",
    sectionId := "genell-thm-3-8" }

def riemann_four_term.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Riemann の 4 項恒等式)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
