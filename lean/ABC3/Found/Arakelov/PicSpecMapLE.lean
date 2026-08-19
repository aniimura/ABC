import ABC3.Found.Arakelov.PicHcompat

/-!
# Arakelov (B2) 第 228 ブロック —— ★★★**`Spec.map` の `appLE` を `ΓSpecIso` で書く**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`hcompat` の右半分

第 227 の左辺に現れる `(Spec.map φ).appLE ⊤ W` を、**環の射 `φ` で書く**:

    (Spec.map φ).appLE ⊤ W e = ΓSpecIso.hom ≫ φ ≫ ΓSpecIso.inv ≫ 制限

★これで第 227 が `hcompat` の形(`appIso`/`ΓSpecIso` の合成)に繋がる。

## ★★`simp` 一発だった

`Scheme.ΓSpecIso_naturality` は **`@[reassoc (attr := simp)]`** が付いているので、
★`simp` が自動で使う。**手で `rw` すると結合の向きで詰まる**が、
`simp` は正規化してから当てるので通る。

★★★摩擦 #3(`simp` の過剰正規化)の**裏返し**である——
`simp` が強すぎて困ることもあれば、`simp` でしか通らないこともある。
★**両方試す**のが実務的である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `specMap_appLE` | ★★★**`Spec.map` の `appLE` を `ΓSpecIso` で書く** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

/-- ★★★**`Spec.map φ` の `appLE` を環の射で書く**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`ΓSpecIso_naturality` は `@[reassoc (attr := simp)]` なので `simp` 一発。 -/
theorem specMap_appLE {R S : CommRingCat.{u}} (φ : R ⟶ S)
    (W : (Spec S).Opens) (e : W ≤ (Spec.map φ) ⁻¹ᵁ ⊤) :
    (Spec.map φ).appLE ⊤ W e
      = (Scheme.ΓSpecIso R).hom ≫ φ ≫ (Scheme.ΓSpecIso S).inv
        ≫ (Spec S).presheaf.map (homOfLE le_top).op := by
  show Scheme.Hom.app (Spec.map φ) ⊤ ≫ _ = _
  simp

/-! ## ★出典の紐付け(`.src`) -/

def specMap_appLE.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Spec.map の appLE を ΓSpecIso で書く)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
