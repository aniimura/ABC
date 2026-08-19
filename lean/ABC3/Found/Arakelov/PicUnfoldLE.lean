import ABC3.Found.Arakelov.PicAppIsoAsLE

/-!
# Arakelov (B2) 第 233 ブロック —— **両辺を `app ≫ 制限` に開く**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`hcompat` の両辺を同じ形に開く

第 232 で言葉が揃ったので、あとは `appLE` の定義

    f.appLE U V e = f.app U ≫ X.presheaf.map (homOfLE e).op

を両辺に当てる。★どちらも **`rfl`** である。

## ★★`appLE` は `abbrev` ではなく `def` だが `rfl` で開く

`Scheme.Hom.appLE` は `def` なので `rw [Scheme.Hom.appLE]` は
「equation theorems で書き換えられない」と言われる(§9-260 で実測)。
★★しかし**等式として `rfl` で書けば通る**——`show` も要らない。

★★★これは摩擦 #3'(`rw` が詰まる)の**もう 1 つの逃げ道**である:
**等式を別建ての補題にして `rfl` で証明する**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `lhs_as_appLE` | ★左辺を `app ≫ 制限` に開く |
| `rhs_as_app` | ★右辺を `app ≫ 制限` に開く |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★**左辺を `app ≫ 制限` に開く**(`rfl`)。 -/
theorem lhs_as_appLE {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B)
    (hle : (hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens))
      ≤ f ⁻¹ᵁ (hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens))) :
    f.appLE (hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens))
        (hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens)) hle
      = f.app (hB.fromSpec ''ᵁ ⊤) ≫ X.presheaf.map (homOfLE hle).op :=
  rfl

/-- ★**右辺を `app ≫ 制限` に開く**(`rfl`)。 -/
theorem rhs_as_app {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B) :
    hB.fromSpec.appLE (hB.fromSpec ''ᵁ ⊤) ⊤ (Scheme.Hom.preimage_image_eq _ _).ge
        ≫ ((Scheme.ΓSpecIso (Γ(Y, B))).hom ≫ f.appLE B A i
          ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv)
        ≫ (hA.fromSpec.appIso ⊤).inv
      = (hB.fromSpec.app (hB.fromSpec ''ᵁ ⊤)
          ≫ (Spec Γ(Y, B)).presheaf.map (homOfLE (Scheme.Hom.preimage_image_eq _ _).ge).op)
        ≫ ((Scheme.ΓSpecIso (Γ(Y, B))).hom ≫ f.appLE B A i
          ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv)
        ≫ (hA.fromSpec.appIso ⊤).inv :=
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def lhs_as_appLE.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——両辺を app と制限に開く)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
