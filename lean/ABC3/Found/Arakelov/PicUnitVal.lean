import ABC3.Found.Arakelov.PicUnitSect

/-!
# Arakelov (B2) 第 172 ブロック —— **`unitMul` の終対象での値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★評価射の構成に要る 1 行

    unitVal ((unitMul U c).app t y) = unitVal y * c

★「`c` 倍」自己射は**本当に掛け算である**ことを言う。
★★評価射 `F ⊗ F^∨ → 𝟙_` の双線型性(第 2 変数)がこれで出る。

## ★★逃げ道——**項で書く**

`rw`/`simp` は `instances` 透明度の型検査で落ちる(3 通り試して全滅):

| 手 | 結果 |
|---|---|
| `show unitVal U (_ • _) = _` から `rw [unitVal_smul]` | パターン不一致 |
| `simp only [unitVal_smul]` | `X.presheaf` の型検査で落ちる |
| **`Eq.trans` の連鎖を項で書く** | ★★**通った** |

★★★[[exact-term-over-rw]] の 6 例目。本 session で繰り返し効いている。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (U : X.Opens)

/-- ★`unitMul` の終対象での値。 -/
theorem unitMul_app_apply (c : (Γ(X, U) : Type u))
    (y : ((𝟙_ (PresheafModulesOn X U)).obj (op (Over.mk (𝟙 U))) : Type u)) :
    unitVal U (((unitMul U c).app (op (Over.mk (𝟙 U)))).hom y) = unitVal U y * c := by
  have hy : y = (y : ((((Over.forget U).op ⋙ X.presheaf)
      ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op (Over.mk (𝟙 U))) : Type u)) • unitOne U :=
    (smul_unitOne U y).symm
  rw [hy, map_smul]
  exact (unitVal_smul U y _).trans
    (((congrArg (fun z => unitVal U y * z) (unitEnd_unitMul U c)).trans
      (congrArg (fun w => w * c) (mul_one (unitVal U y)).symm)).trans
      (congrArg (fun w => w * c) (unitVal_smul U y (unitOne U)).symm))


/-! ## ★出典の紐付け(`.src`) -/

def unitMul_app_apply.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——unitMul の終対象での値)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
