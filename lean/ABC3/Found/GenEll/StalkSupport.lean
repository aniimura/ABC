import ABC3.Found.GenEll.StalkPullback

/-!
# [GenEll] Definition 1.5, (iv) の足場 —— **茎と台をつなぐ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> as arithmetic divisors ∈ ADiv(F ) that are supported in V(F )non. We shall refer to

## ★★★なぜこれが要るのか

`x_F^* D` を `ADiv(F)` の元として書くには、
**「有限個の素点でしか係数が非零でない」**ことが要る(`Finsupp` だから)。

★★これは幾何的には「`x_F` が `D` の台と有限個の点でしか交わらない」ことであり、
mathlib には **`support_comap`**(`(I.comap f).support = f⁻¹(I.support)`)がある。

★★★**足りないのは「茎」と「台」をつなぐ 1 本**である:

    `I.stalk x = ⊤  ↔  x ∉ I.support`

★これが本ファイルの主定理である。

## ★機構

`x ∈ I.support ↔ x ∈ zeroLocus (I.ideal U)`(`mem_support_iff_of_mem`)であり、
`x ∈ zeroLocus s ↔ ∀ f ∈ s, x ∉ X.basicOpen f`(`mem_zeroLocus_iff`)、
`x ∈ X.basicOpen f ↔ IsUnit (germ f)`(`mem_basicOpen`)である。

★★**茎は局所環**なので、イデアルが `⊤` であることと**単元を含む**ことは同値。
したがって 3 本を繋ぐだけで出る。

## ★★これが何を可能にするか

`comap` の積の保存(`Proposition 1.4, (i)` に残った 1 本)を**迂回する**道が開く:
台が有限なら、各素点での重複度から `Finsupp` を直接作れる。
★重複度は `stalkPullback` から取れ、**`stalkPullback` は積を保つ**
(`StalkPullback.lean`)。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

variable {X : Scheme}

/-! ## ★★★茎が `⊤` であることと台の外にいることは同値 -/

/-- ★★★**イデアル層の茎が `⊤` ⟺ 点が台の外にある**。

★★これが「茎」(本プロジェクトの定義)と「台」(mathlib の定義)をつなぐ 1 本である。

★機構は 3 本の合成:
`mem_support_iff_of_mem` / `mem_zeroLocus_iff` / `mem_basicOpen`。
★★茎が**局所環**であることが効く——イデアルが `⊤` ⟺ 単元を含む。 -/
theorem stalk_eq_top_iff (I : X.IdealSheafData) (x : X) :
    I.stalk x = ⊤ ↔ x ∉ I.support := by
  obtain ⟨_, ⟨U, hU, rfl⟩, hxU, -⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open (Set.mem_univ x) isOpen_univ
  rw [Scheme.IdealSheafData.mem_support_iff_of_mem (U := ⟨U, hU⟩) hxU,
    Scheme.IdealSheafData.stalk_eq_map (U := ⟨U, hU⟩) hxU,
    Scheme.mem_zeroLocus_iff]
  constructor
  · -- ★すべての生成元の芽が非単元なら、像は極大イデアルに入る——`⊤` と矛盾する
    intro htop hall
    have hle : (I.ideal ⟨U, hU⟩).map (X.presheaf.germ U x hxU).hom
        ≤ IsLocalRing.maximalIdeal (X.presheaf.stalk x) := by
      rw [Ideal.map_le_iff_le_comap]
      intro f hf
      refine IsLocalRing.mem_maximalIdeal _ |>.2 fun hcon => ?_
      exact hall f hf ((Scheme.mem_basicOpen X f x hxU).2 hcon)
    rw [htop, top_le_iff] at hle
    exact (IsLocalRing.maximalIdeal.isMaximal (X.presheaf.stalk x)).ne_top hle
  · -- 台の外なら、芽が単元になる生成元が取れる
    intro hns
    rw [Classical.not_forall] at hns
    obtain ⟨f, hf⟩ := hns
    rw [Classical.not_imp, Classical.not_not] at hf
    obtain ⟨hfI, hfb⟩ := hf
    have hunit : IsUnit (X.presheaf.germ U x hxU f) :=
      (Scheme.mem_basicOpen X f x hxU).1 hfb
    exact Ideal.eq_top_of_isUnit_mem _ (Ideal.mem_map_of_mem _ hfI) hunit

/-- ★**台の外では引き戻しも `⊤`**。

★★これが「有限個の素点でしか係数が非零でない」ことの幾何的な内容である。 -/
theorem stalkPullback_eq_top_of_not_mem_support {Y : Scheme} (I : Y.IdealSheafData)
    (f : X ⟶ Y) (x : X) (h : f.base x ∉ I.support) :
    stalkPullback I f x = ⊤ := by
  rw [stalkPullback, (stalk_eq_top_iff I (f.base x)).2 h, Ideal.map_top]

/-- ★★**引き戻しが `⊤` でない点は、台の引き戻しに入る**。

★mathlib の `support_comap`(`(I.comap f).support = f⁻¹(I.support)`)と合わせると、
**その集合は閉集合**である。★★`Spec 𝓞_F` は 1 次元なので、
生成点を含まない閉集合は**有限**——それが `Finsupp` に載る根拠になる。 -/
theorem mem_preimage_support_of_stalkPullback_ne_top {Y : Scheme}
    (I : Y.IdealSheafData) (f : X ⟶ Y) (x : X) (h : stalkPullback I f x ≠ ⊤) :
    f.base x ∈ I.support := by
  by_contra hc
  exact h (stalkPullback_eq_top_of_not_mem_support I f x hc)

/-! ## ★出典の紐付け(`.src`) -/

def stalk_eq_top_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(茎と台の対応のみ)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
