import ABC3.Found.Arakelov.PicEquivRing
import ABC3.Found.Arakelov.PicInterface
import ABC3.Found.Arakelov.PicPullGroup
import ABC3.Interface.Arakelov.LineBundle

/-!
# Arakelov (B1) 第 146 ブロック —— ★★★★★★★★★★★**`PicardData` の witness**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★★B1 の `PicardData` が閉じた

`Interface` の `PicardData` は 14 欄。★**すべて埋まった。**

| 欄 | witness | 出所 |
|---|---|---|
| `Pic` | `PicSheaf` | 第 62 |
| `group` | `CommGroup (PicSheaf X)` | 第 62 |
| `pullback` | `picPullback` | 第 63 |
| `pullback_mul` | `picPullback_mul` | 第 63(第 18–60 の帰着先) |
| `pullback_id` / `pullback_comp` | 同上 | 第 63 |
| **`equivPicRing`** | **`equivPicRingSheaf`** | ★★**第 145** |
| `sheafOf` 系 7 欄 | `picSheafOf_*` | 第 73 |

## ★★★逸脱の記録(2026-08-18)

`Interface` の `Pic : Scheme.{0} → Type` を **`Type 1`** に緩めた。

★理由: `PicSheaf X = Quotient (InvSheaf.setoid X)` は `InvSheaf X : Type 1`
(層加群の圏の対象を含む)の商なので `Type 1` である。
★★`Type 0` へ落とすには `Small.{0} (PicSheaf X)` が要り、一般のスキームでは
それ自体が別の仕事になる。

★★★影響: `PicardData` を**消費するコードは無い**(2026-08-18 時点で
`Check/` の負の対照は独自の構造を持ち、`PicardData` の instance ではない)。
`equivPicRing` は `MulEquiv` が宇宙をまたげるのでそのまま通る。

★★★★数学的内容は変わらない——`Pic X` は同型類の群であり、
それがどの宇宙に住むかは ABC の議論に影響しない。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory

/-- ★★★★★★★★★★★**B1 の witness**——`PicardData` の 14 欄すべて。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 1 ブロックから 146 ブロックの到達点である。 -/
noncomputable def picardDataWitness : ABC3.Interface.Arakelov.PicardData where
  Pic := fun X => PicSheaf X
  group := fun _ => inferInstance
  pullback := fun f L => picPullback f L
  pullback_mul := fun f L M => picPullback_mul f L M
  pullback_id := fun L => picPullback_id L
  pullback_comp := fun f g L => picPullback_comp f g L
  equivPicRing := fun R => equivPicRingSheaf R
  sheafOf := fun X L => picSheafOf X L
  sheafOf_invertible := fun X L => picSheafOf_invertible X L
  sheafOf_one := fun X => picSheafOf_one X
  sheafOf_mul := fun X L M => picSheafOf_mul X L M
  sheafOf_pullback := fun f L => picSheafOf_pullback f L
  sheafOf_injective := fun X L M h => picSheafOf_injective X L M h
  sheafOf_surjective := fun X F h => picSheafOf_surjective X F h

/-! ## ★出典の紐付け(`.src`) -/

def picardDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——PicardData の witness)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
