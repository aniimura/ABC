import ABC3.Found.Arakelov.PicPushMu

/-!
# Arakelov (B1) 第 33 ブロック —— **`δ` を計算するための手がかりが揃った**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★手がかりの一覧(第 31–33 ブロック)

| 手がかり | どこで |
|---|---|
| `adj.homEquiv (δ P Q) = (unit P ⊗ₘ unit Q) ≫ μ` | ★本ブロック(`Equiv.apply_symm_apply`) |
| `unit (free (yoneda V))` の書き下し | ★第 31 ブロック |
| mathlib 側の `homEquiv` は `freeYonedaEquiv` の合成 | ★本ブロック(`rfl`) |
| `μ` は切断の上で純テンソルを純テンソルに送る | ★第 32 ブロック |

★★★**これで `δ` は生成元の切断の上で完全に計算できる。**

## ★★残り

生成元 `V, W` について

    IsIso (δ (free (yoneda V)) (free (yoneda W)))

を示せば、第 30 ブロックの二重持ち上げで**全対象**に伝播する。
★対象の同型は第 27 ブロックで既にある(`pullbackTensorFreeIso`)。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-! ## ★`δ` の adjunct -/

/-- ★★★**`δ` の adjunct は定義式そのものである**。

★`Adjunction.leftAdjointOplaxMonoidal` は
`δ X Y := homEquiv.symm ((unit ⊗ₘ unit) ≫ μ)` と定義しているので、
`homEquiv` を当て直せば取り出せる。 -/
theorem homEquiv_pullbackDelta (P Q : Y.PresheafOfModules) :
    (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).homEquiv _ _
        (pullbackDelta f P Q)
      = ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app P
          ⊗ₘ (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app Q)
        ≫ Functor.LaxMonoidal.μ (PresheafOfModules.pushforward (pullbackPhi f)) _ _ :=
  Equiv.apply_symm_apply _ _

/-! ## ★mathlib 側の余表現 -/

/-- ★★**mathlib 側の余表現の `homEquiv` は `freeYonedaEquiv` の合成である**(`rfl`)。

★★★これが `unit` の書き下し(第 31 ブロック)を**要素の言葉**に翻訳する。 -/
theorem mathlibCorep_homEquiv (V : Y.Opens) (N : X.PresheafOfModules)
    (g : (PresheafOfModules.free (X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
      (yoneda.obj ((Opens.map f.base).obj V)) ⟶ N) :
    (mathlibCorep f V).homEquiv g
      = (PresheafOfModules.freeYonedaEquiv
          (M := (PresheafOfModules.pushforward (pullbackPhi f)).obj N) (X := V)).symm
        (PresheafOfModules.freeYonedaEquiv
          (M := N) (X := (Opens.map f.base).obj V) g) := rfl

/-! ## ★出典の紐付け(`.src`) -/

def homEquiv_pullbackDelta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——δ の adjunct が定義式であること)",
    sectionId := "genell-def-1-1-i" }

def mathlibCorep_homEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——余表現の homEquiv が freeYonedaEquiv の合成であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
