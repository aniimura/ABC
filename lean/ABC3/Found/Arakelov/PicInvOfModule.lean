import ABC3.Found.Arakelov.PicTildeLocTrivial
import ABC3.Found.Arakelov.PicTildeTensorIso
import ABC3.Found.Arakelov.PicTildeUnit
import ABC3.Found.Arakelov.PicSheafGroup

/-!
# Arakelov (B1) 第 133 ブロック —— ★★★★★★**可逆加群から可逆層へ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★`equivPicRing` の片道が通った

第 132(局所自明性)と第 91(`tilde` はテンソルを保つ)が揃ったので、

    Module.Invertible R M  ⟹  InvSheaf (Spec R)

が**一発で**書ける。

| `InvSheaf` の欄 | 中身 |
|---|---|
| `carrier` | `tilde M` |
| `inv` | `tilde (Mᵛ)` |
| `isInv` | 第 91 + `Module.Invertible.linearEquiv` + 第 82(`tilde R = 𝒪`) |
| `trivial` | ★第 132 |
| `invTrivial` | ★第 132(`Mᵛ` も可逆) |

## ★★これが `PicardData` の最後の欄 `equivPicRing` の**前進側**である

残るのは逆向き(可逆層 ⟹ 可逆加群)と、両者が互いに逆であることである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

/-- ★★★★★★**可逆加群は可逆層を与える**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 132(局所自明性)と第 91(`tilde` はテンソルを保つ)の合流点である。 -/
noncomputable def invSheafOfModule (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u))
    [Module.Invertible (R : Type u) (M : Type u)] : InvSheaf (Spec R) where
  carrier := tilde M
  inv := tilde (ModuleCat.of (R : Type u) (Module.Dual (R : Type u) (M : Type u)))
  isInv := ⟨tildeTensorIso R M _ ≪≫
    (tildeFunctor R).mapIso (LinearEquiv.toModuleIso
      (TensorProduct.comm (R : Type u) (M : Type u)
          (Module.Dual (R : Type u) (M : Type u)) ≪≫ₗ
        Module.Invertible.linearEquiv (R : Type u) (M : Type u))) ≪≫ tildeUnit R⟩
  trivial := isLocallyTrivial_tilde R M
  invTrivial := isLocallyTrivial_tilde R _

/-! ## ★出典の紐付け(`.src`) -/

def invSheafOfModule.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆加群から可逆層へ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
