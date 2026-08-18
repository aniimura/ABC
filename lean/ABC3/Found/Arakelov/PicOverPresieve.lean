import ABC3.Found.Arakelov.PicPresheafId

/-!
# Arakelov (B1) 第 122 ブロック —— **`D(h·g)` の生成族は `over` でも覆う**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性の被覆が `over` 位相で揃った

第 117 ブロックで「`D(h·g)` の生成族は `U` を覆う」が出た。
★それを `Over V` の site へ持ち上げる。

★★mathlib の `Sieve.overEquiv_symm_generate` が**生成族の対応**を与える(実測):

    (overEquiv Y).symm (generate R) = generate (Presieve.functorPullback (Over.forget X) R)

★★★これと第 116(`over` の被覆は下で判定)を合わせて**一発**で通った。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `overBasicPresieve_mem` | ★★★★**`D(h·g)` の生成族は `over` でも覆う** |

## ★★★これで第 115 の被覆の仮定が埋まる

残るのは「その生成族の上で `s` の倍が全単射」——
第 120(制限と局所化の可換図式)+ 第 112(局所化は全単射を保つ)
+ 第 103(生成元の乗法は全単射)である。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (g : (R : Type u))

/-- ★★★★**`D(h·g)` の生成族は `Over V` の site でも覆う**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 117(空間の側の被覆)を `overEquiv_symm_generate` で持ち上げたものである。 -/
theorem overBasicPresieve_mem (V : (Spec R).Opens) (hV : V ≤ PrimeSpectrum.basicOpen g)
    (W : Over V) :
    Sieve.generate (Presieve.functorPullback (Over.forget V)
        (fun (Z : (Spec R).Opens) (_ : Z ⟶ W.left) =>
          ∃ h : (R : Type u), Z = PrimeSpectrum.basicOpen (h * g)))
      ∈ ((Opens.grothendieckTopology (Spec R)).over V) W := by
  rw [← Sieve.overEquiv_symm_generate]
  refine symm_mem_over V W _ ?_
  exact basicOpenMulSieve_mem R g W.left (le_trans (leOfHom W.hom) hV)

/-! ## ★出典の紐付け(`.src`) -/

def overBasicPresieve_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——D(h·g) の生成族は over でも覆う)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
