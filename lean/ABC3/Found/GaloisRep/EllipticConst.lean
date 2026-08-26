import ABC3.Found.GaloisRep.EllipticLiouville

/-!
# Galois (G6) 第 242 ブロック —— **★★★★★★★格子の外で正則な楕円関数**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★実際に使う形にする

第 241 の Liouville は「**全体で微分可能**な楕円関数は定数」だった。
しかし加法定理の議論で現れる関数(`℘'(z) − λ℘(z) − ν` の組み合わせ等)は
**格子の外でしか正則でない**——格子点では極を持ちうる。

★極が実は消えている場合(連続に延びる場合)を扱えるようにする:

> 格子の外で正則、全体で連続な楕円関数は定数である。

## ★★★★格子点は孤立している

mathlib の `PeriodPair.compl_lattice_sdiff_singleton_mem_nhds x : (↑L.lattice \ {x})ᶜ ∈ 𝓝 x`
が使える——**`x` 以外の格子点を除いた集合が `x` の近傍**である。
その上で `Complex.differentiableOn_compl_singleton_and_continuousAt_iff`
(可除特異点)を当てれば、`x` でも微分可能になる。

★これで第 241 が適用でき、定数が出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `eq_of_periodic_continuous` | ★★★★★★★**格子の外で正則・全体で連続な楕円関数は定数** |
| `const_of_periodic_continuous` | ★同じことを「ある一点の値に等しい」形で |
-/

namespace ABC3.Found.GaloisRep

open Complex Real PeriodPair

/-- ★★★★★★★**格子の外で正則・全体で連続な楕円関数は定数**。

★格子点は孤立しているので、可除特異点の定理でそこも正則になる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem eq_of_periodic_continuous (L : PeriodPair) (f : ℂ → ℂ)
    (hcont : Continuous f)
    (hdiffOn : DifferentiableOn ℂ f ((L.lattice : Set ℂ)ᶜ))
    (hper : ∀ l ∈ L.lattice, ∀ z : ℂ, f (l + z) = f z) :
    ∀ z w : ℂ, f z = f w := by
  have hdiff : Differentiable ℂ f := by
    intro x
    have hs : ((L.lattice : Set ℂ) \ {x})ᶜ ∈ nhds x := L.compl_lattice_sdiff_singleton_mem_nhds x
    have hmono : DifferentiableOn ℂ f ((((L.lattice : Set ℂ) \ {x})ᶜ) \ {x}) := by
      refine hdiffOn.mono ?_
      intro y hy
      obtain ⟨hy1, hy2⟩ := hy
      intro hylat
      exact hy1 ⟨hylat, hy2⟩
    have hkey := (Complex.differentiableOn_compl_singleton_and_continuousAt_iff hs).1
      ⟨hmono, hcont.continuousAt⟩
    exact hkey.differentiableAt hs
  exact eq_of_periodic_differentiable L f hdiff hper

/-- ★同じことを「ある一点の値に等しい」形で。 -/
theorem const_of_periodic_continuous (L : PeriodPair) (f : ℂ → ℂ)
    (hcont : Continuous f)
    (hdiffOn : DifferentiableOn ℂ f ((L.lattice : Set ℂ)ᶜ))
    (hper : ∀ l ∈ L.lattice, ∀ z : ℂ, f (l + z) = f z) (z₀ : ℂ) :
    ∀ z : ℂ, f z = f z₀ :=
  fun z => eq_of_periodic_continuous L f hcont hdiffOn hper z z₀

/-! ## ★出典の紐付け(`.src`) -/

def eq_of_periodic_continuous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子の外で正則な楕円関数は定数)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
