import ABC3.Found.Arakelov.PicInterface

/-!
# Arakelov (B1) 第 74 ブロック —— **切断の比較射 `Γ(F) ⊗_R Γ(G) ⟶ Γ(F ⊗ G)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-66 の旗もまた杞憂だった

§9-66 は「`Γ(F)` の `R` 加群構造と `F.val(⊤)` の `𝒪(⊤)` 加群構造を
`IsScalarTower` で繋ぐのが残りの作業」と書いた。

★★**測った結果(2026-08-18)**: mathlib の `AlgebraicGeometry/Modules/Tilde.lean` に

    instance : Module R Γ(M, U)
    instance : IsScalarTower R Γ(Spec R, U) Γ(M, U)

が**既にある**。★★★しかも `Γ(M, U)` は `Ab` の対象なので、
要素を使った書き方では **`Module 𝒪(⊤) Γ(G,⊤)` の方が見つからない**
——[[ring-instance-two-paths]] と同じ二経路である(実測)。

## ★★★要素を捨てると 2 行で済む

★**`moduleSpecΓFunctor.obj F = (restrictScalars ι).obj (F.val(⊤))` は `rfl`**(実測)。
★★したがって比較射は

| 段 | 射 |
|---|---|
| 前 | `restrictScalars` の **lax monoidal 構造** `μ`(mathlib の在庫) |
| 後 | 第 68 ブロックの `gammaSheafifyUnit` を `restrictScalars` で送ったもの |

の合成である。**要素は 1 つも要らない。**

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `gammaTensorCmp` | ★★★★**`Γ(F) ⊗_R Γ(G) ⟶ Γ(F ⊗ G)`** |
| `gammaEqRestrict` | ★`Γ(F)` は `restrictScalars` の像(`rfl`) |

## ★★★次

これを `tilde ⊣ Γ` の随伴で転置すると

    tilde (M ⊗_R N)  ⟶  tensorModules (tilde M) (tilde N)

が出る——それが `equivPicRing` の乗法の段である。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (R : CommRingCat.{u}) (F G : (Spec R).Modules)

/-- ★**`Γ(F)` は `F.val(⊤)` の係数制限である**(`rfl`)。 -/
theorem gammaEqRestrict :
    AlgebraicGeometry.moduleSpecΓFunctor.obj F
      = (ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom).obj
          (F.val.obj (op ⊤)) := rfl

/-- ★★★★**切断の比較射** `Γ(F) ⊗_R Γ(G) ⟶ Γ(F ⊗ G)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★2 段の合成である:

    Γ(F) ⊗_R Γ(G) --μ--> (F.val ⊗ G.val)(⊤) --層化の unit--> (層化 ..)(⊤) = Γ(F ⊗ G)

★前段は `restrictScalars` の lax monoidal 構造(mathlib)、
後段は第 68 ブロックの `gammaSheafifyUnit` である。 -/
noncomputable def gammaTensorCmp :
    AlgebraicGeometry.moduleSpecΓFunctor.obj F
        ⊗ AlgebraicGeometry.moduleSpecΓFunctor.obj G
      ⟶ AlgebraicGeometry.moduleSpecΓFunctor.obj (tensorModules F G) :=
  Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom)
      (F.val.obj (op ⊤)) (G.val.obj (op ⊤))
    ≫ (ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom).map (gammaSheafifyUnit R F G)

/-! ## ★出典の紐付け(`.src`) -/

def gammaTensorCmp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断の比較射 Γ(F) ⊗ Γ(G) → Γ(F ⊗ G))",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
