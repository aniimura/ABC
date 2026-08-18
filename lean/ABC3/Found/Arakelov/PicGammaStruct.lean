import ABC3.Found.Arakelov.PicBaseTensor

/-!
# Arakelov (B1) 第 70 ブロック —— **`Γ` の `R` 加群構造の在り処**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★測ったら係数制限は既に中に入っていた

§9-66 で「`Γ(F)` の `R` 加群構造と `F.val(⊤)` の `𝒪(⊤)` 加群構造を
`IsScalarTower` で繋ぐ必要がある」と旗を立てた。

★★**測った結果(2026-08-18)**: `modulesSpecToSheaf` の中に
`ModuleCat.restrictScalars (ΓSpecIso R).inv` が**既に入っている**ので、

    Γ(F) = (modulesSpecToSheaf.obj F).val.obj (op ⊤)      (`rfl`、しかも `ModuleCat R`)

★★★**`IsScalarTower` を自分で繋ぐ必要は無かった。**

## ★★これで前段は `restrictScalars` の `μ` そのもの

    restrictScalars ι A ⊗ restrictScalars ι B  ⟶  restrictScalars ι (A ⊗ B)

は mathlib の `Functor.LaxMonoidal.μ` である(第 18 ブロックで使ったもの)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `gammaEqTopSection` | ★★`Γ(F)` は `modulesSpecToSheaf` の `⊤` 切断(`rfl`) |
| `restrictScalarsMu` | ★前段の射(mathlib の `μ`) |

## ★方法論

★★§9-66 で立てた旗は**杞憂だった**——第 57 ブロックと同じである。
★★★だが旗を立てたから**測った**のであり、測ったから**短く済んだ**。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (R : CommRingCat.{u})

/-- ★★**`Γ(F)` は `modulesSpecToSheaf` の `⊤` での切断である**(`rfl`)。

★★★係数制限は `modulesSpecToSheaf` の中に既に入っている。 -/
theorem gammaEqTopSection (F : (Spec R).Modules) :
    AlgebraicGeometry.moduleSpecΓFunctor.obj F
      = (AlgebraicGeometry.modulesSpecToSheaf.obj F).val.obj (op ⊤) := rfl

/-- ★**係数制限の `μ`**——比較射の前段である(mathlib の在庫)。 -/
noncomputable abbrev restrictScalarsMu (A B : ModuleCat (Γ(Spec R, ⊤) : Type u)) :
    (ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom).obj A
        ⊗ (ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom).obj B
      ⟶ (ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom).obj (A ⊗ B) :=
  Functor.LaxMonoidal.μ _ A B

/-! ## ★出典の紐付け(`.src`) -/

def gammaEqTopSection.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Γ の R 加群構造の在り処)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
