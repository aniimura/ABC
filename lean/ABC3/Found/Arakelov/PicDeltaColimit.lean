import ABC3.Found.Arakelov.PicSchemeDelta
import Mathlib.Algebra.Category.ModuleCat.Presheaf.Generator

/-!
# Arakelov (B1) 第 23 ブロック —— **`δ` を同型にするための余極限の部品**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★数学は見えている

§9-41 で訂正した通り、**`f^*` は 𝒪 加群について strong monoidal である**
(平坦性は不要。Stacks 01BJ)。★したがって `δ` は**全対象で同型**である。

★★★Lean で出す道は**生成元による道**である:

| 段 | 主張 | 在庫 |
|---|---|---|
| 1 | 両辺とも各変数で余極限を保つ | ★★本ブロック |
| 2 | 各前層加群は `free (yoneda V)` の余極限 | ★mathlib `isColimitFreeYonedaCoproductsCokernelCofork` |
| 3 | 生成元の上で `δ` は同型 | ★`free(yoneda V) ⊗ free(yoneda W) ≅ free(yoneda (V ⊓ W))` かつ `F` が `⊓` を保つ |

★★★★**段 3 が効く理由**: 前順序では `yoneda V × yoneda W = yoneda (V ⊓ W)` であり、
`F = Opens.map f.base`(逆像)は `⊓` を保つ。
★したがって生成元の上では `δ` は `free(yoneda (FV ⊓ FW))` の恒等射になる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pullbackPre_preservesColimits` | ★★引き戻しは余極限を保つ(左随伴) |
| `tensorLeft_preservesColimits` | ★テンソルは余極限を保つ(mathlib のインスタンス) |
| `freeYoneda_isDetecting` | ★★`free (yoneda -)` は検出族 |

## ★★★実測メモ(2026-08-18)

★`PreservesColimits (tensorRight Q)` は**探索に失敗する**が、
`PreservesColimitsOfSize.{u, u} (tensorRight Q)` は**在る**
(`Presheaf/Monoidal.lean` 末尾)。★★universe の指定が違うだけであった。
——「無いと決める前に測る」の 8 例目である。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-! ## ★★引き戻しは余極限を保つ -/

/-- ★★**引き戻しは余極限を保つ**——左随伴だからである。

★これが生成元による議論の第 1 段である。 -/
noncomputable instance pullbackPre_preservesColimits : PreservesColimits (pullbackPre f) :=
  (PresheafOfModules.pullbackPushforwardAdjunction
    (alphaR (Opens.map f.base) X.presheaf Y.presheaf f.c)).leftAdjoint_preservesColimits

/-! ## ★テンソルは余極限を保つ -/

/-- ★**テンソルは余極限を保つ**(左から)。

★★★mathlib の `Presheaf/Monoidal.lean` 末尾にインスタンスが在る。
★`PreservesColimits`(= `PreservesColimitsOfSize.{v, v}`)では探索に失敗し、
`PreservesColimitsOfSize.{u, u}` なら通る——**universe の指定だけの違い**であった。 -/
theorem tensorLeft_preservesColimits
    (P : PresheafOfModules.{u} (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})) :
    PreservesColimitsOfSize.{u, u} (tensorLeft P) := inferInstance

/-- ★**テンソルは余極限を保つ**(右から)。 -/
theorem tensorRight_preservesColimits
    (Q : PresheafOfModules.{u} (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})) :
    PreservesColimitsOfSize.{u, u} (tensorRight Q) := inferInstance

/-! ## ★★生成元 -/

/-- ★★**`free (yoneda -)` は前層加群の検出族である**。

★★★これがあるので「生成元の上で同型ならば同型」が言える。 -/
theorem freeYoneda_isDetecting :
    ObjectProperty.IsDetecting
      (PresheafOfModules.freeYoneda (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})) :=
  PresheafOfModules.freeYoneda.isDetecting _

/-! ## ★★★逆像は交わりを保つ -/

/-- ★★★**開集合の逆像は交わりを保つ**。

★★これが生成元の上で `δ` が同型になる理由である:
前順序では `yoneda V × yoneda W = yoneda (V ⊓ W)` なので、
`δ` は `free (yoneda (F V ⊓ F W))` の恒等射になる。 -/
theorem opensMap_inf (V W : Y.Opens) :
    (Opens.map f.base).obj (V ⊓ W)
      = (Opens.map f.base).obj V ⊓ (Opens.map f.base).obj W := rfl

/-! ## ★出典の紐付け(`.src`) -/

def pullbackPre_preservesColimits.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——引き戻しが余極限を保つこと)",
    sectionId := "genell-def-1-1-i" }

def opensMap_inf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——逆像が交わりを保つこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
