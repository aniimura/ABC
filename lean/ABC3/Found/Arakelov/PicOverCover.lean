import ABC3.Found.Arakelov.PicPresieveIso

/-!
# Arakelov (B1) 第 116 ブロック —— **`over` 位相の被覆は下で判定できる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★生成族の被覆性を空間の側で言う

第 115 ブロックの器具は `Sieve.generate R ∈ (J.over V) W` を要求する。
★これを**空間 `X` の側**(`W.left` の被覆)で判定できれば、
基本開集合が基底であることがそのまま使える。

★★mathlib の `GrothendieckTopology.mem_over_iff` が**まさにそれ**である
(2026-08-24 実測、`CategoryTheory/Sites/Over.lean` 行 242):

    S ∈ (J.over X) Y  ↔  Sieve.overEquiv _ S ∈ J Y.left      (**`rfl`**)

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `mem_over_iff'` | ★`over` の被覆は下で判定できる |
| `symm_mem_over` | ★下の被覆篩を `over` へ持ち上げる |
| `generate_mem_over` | ★★★**生成族の被覆性を下で言う** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (V : X.Opens)

/-- ★**`over` 位相の被覆は下で判定できる**(mathlib の在庫に名を付ける)。 -/
theorem mem_over_iff' (W : Over V) (S : Sieve W) :
    S ∈ ((Opens.grothendieckTopology X).over V) W
      ↔ Sieve.overEquiv W S ∈ (Opens.grothendieckTopology X) W.left :=
  GrothendieckTopology.mem_over_iff _ S

/-- ★**下の被覆篩を `over` へ持ち上げる**。 -/
theorem symm_mem_over (W : Over V) (S : Sieve W.left)
    (hS : S ∈ (Opens.grothendieckTopology X) W.left) :
    (Sieve.overEquiv W).symm S ∈ ((Opens.grothendieckTopology X).over V) W :=
  GrothendieckTopology.overEquiv_symm_mem_over _ W S hS

/-- ★★★**生成族の被覆性を空間の側で言う**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これで「基本開集合が基底である」がそのまま使える。 -/
theorem generate_mem_over (W : Over V) (R : Presieve W)
    (h : Sieve.overEquiv W (Sieve.generate R) ∈ (Opens.grothendieckTopology X) W.left) :
    Sieve.generate R ∈ ((Opens.grothendieckTopology X).over V) W :=
  (mem_over_iff' V W _).2 h

/-! ## ★出典の紐付け(`.src`) -/

def generate_mem_over.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——over 位相の被覆は下で判定できる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
