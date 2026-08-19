import ABC3.Found.Arakelov.PicHcompatTop
import ABC3.Found.Arakelov.PicSurjPair

/-!
# Arakelov (B2) 第 232 ブロック —— **`appIso` を `appLE` の言葉に直す**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 221 と第 231 を繋ぐ最後の変換

第 221 の `hcompat` は **`appIso ⊤`** で書いてあり、
第 231 は **`app` / `appLE`** で書いてある。★mathlib の

    appIso_hom' : (f.appIso U).hom = f.appLE (f ''ᵁ U) U _

がこの 2 つを繋ぐ。

## ★★これで座標も言葉も揃った

| 側面 | 第 221 | 第 231 | 変換 |
|---|---|---|---|
| 座標 | `⊤` | ★`⊤`(第 231 で揃えた) | 済 |
| 言葉 | `appIso` | `appLE` | ★**本ブロック** |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `appIso_top_as_app` | ★`appIso ⊤` は `appLE` である |
| `hcompat_as_appLE` | ★★`hcompat` を `appLE` の言葉に直す |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★**`appIso ⊤` は `appLE` である**(mathlib `appIso_hom'`)。 -/
theorem appIso_top_as_app {A : X.Opens} (hA : IsAffineOpen A) :
    (hA.fromSpec.appIso (⊤ : (Spec Γ(X, A)).Opens)).hom
      = hA.fromSpec.appLE (hA.fromSpec ''ᵁ ⊤) ⊤ (Scheme.Hom.preimage_image_eq _ _).ge :=
  Scheme.Hom.appIso_hom' _ _

/-- ★★**`hcompat` を `appLE` の言葉に直す**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで第 221(`appIso`)と第 231(`appLE`)が同じ言葉になる。 -/
theorem hcompat_as_appLE {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B) :
    ((hB.fromSpec.appIso ⊤).hom ≫ ((Scheme.ΓSpecIso (Γ(Y, B))).hom
        ≫ f.appLE B A i ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv)
        ≫ (hA.fromSpec.appIso ⊤).inv)
      = hB.fromSpec.appLE (hB.fromSpec ''ᵁ ⊤) ⊤ (Scheme.Hom.preimage_image_eq _ _).ge
        ≫ ((Scheme.ΓSpecIso (Γ(Y, B))).hom ≫ f.appLE B A i
          ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv)
        ≫ (hA.fromSpec.appIso ⊤).inv := by
  rw [appIso_top_as_app hB]

/-! ## ★出典の紐付け(`.src`) -/

def hcompat_as_appLE.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——appIso を appLE の言葉に直す)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
