import ABC3.Found.Arakelov.PicUnitVal
import ABC3.Found.Arakelov.PicDualPre

/-!
# Arakelov (B2) 第 173 ブロック —— **評価射のための両側の橋**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★評価射は**両側で綴りが違う**

評価射 `F ⊗ F^∨ → 𝟙_` の双線型性を書くとき、

| 側 | 綴り |
|---|---|
| `φ.app t` の線型性 | `PresheafModulesOn X W` の環(**RingCat 綴り、終対象**) |
| `x : F.obj W` の加群 | `X.presheaf.obj W`(**CommRingCat 綴り**) |

★どちらに寄せても片側が欠ける。

## ★★逃げ道——**両側に型付き恒等関数**([[typed-identity-bridge]] の 2 例目)

    fVal : F.obj W → ((restrictPresheafFunctor X W).obj F).obj t
    rVal : X.presheaf.obj W → (RingCat 綴りの環)

★★これで `fVal (c • x) = rVal c • fVal x` が **`rfl`** で通る。

★★★教訓: **橋は片側では足りない**。値と係数の両方に架ける。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (F : X.PresheafOfModules)


/-- ★切断を制限した前層の終対象での切断と見る。 -/
def fVal (W : (X.Opens)ᵒᵖ) (x : (F.obj W : Type u)) :
    (((restrictPresheafFunctor X W.unop).obj F).obj (op (Over.mk (𝟙 W.unop))) : Type u) := x

/-- ★係数も同じく橋を架ける。 -/
def rVal (W : (X.Opens)ᵒᵖ) (c : ((X.presheaf.obj W) : Type u)) :
    ((((Over.forget W.unop).op ⋙ X.presheaf) ⋙ forget₂ CommRingCat.{u}
      RingCat.{u}).obj (op (Over.mk (𝟙 W.unop))) : Type u) := c

/-- ★★作用も一致する。 -/
theorem fVal_smul (W : (X.Opens)ᵒᵖ) (c : ((X.presheaf.obj W) : Type u))
    (x : (F.obj W : Type u)) :
    fVal F W (c • x) = rVal (X := X) W c • fVal F W x := rfl

/-- ★★加法も一致する。 -/
theorem fVal_add (W : (X.Opens)ᵒᵖ) (x y : (F.obj W : Type u)) :
    fVal F W (x + y) = fVal F W x + fVal F W y := rfl


/-! ## ★出典の紐付け(`.src`) -/

def fVal_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——評価射のための両側の橋)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
