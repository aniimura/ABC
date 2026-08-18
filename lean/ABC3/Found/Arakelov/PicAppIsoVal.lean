import ABC3.Found.Arakelov.PicFromSpecApp

/-!
# Arakelov (B2) 第 217 ブロック —— ★★★★**`appIso` を元のレベルで開く**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-242 で「`eqToHom` の輸送が重い」と測った所

射の等式として書くと `simp` が正規化してしまうので、★**元のレベルで書く**。

| 定理 | 内容 |
|---|---|
| `appIso_top_apply` | ★★`(appIso ⊤).hom s` を `app ≫ 制限` に開く |
| `fromSpec_app_image_top` | ★★★`fromSpec ''ᵁ ⊤` の app と `U` の app の関係 |

## ★★★元のレベルなら `rfl` と `congrArg` で済む

★`appIso_top_apply` は `appIso_hom'`(mathlib)を `rw` して **`rfl`**。
★★`fromSpec_app_image_top` は `fromSpec.naturality` に元を当てた **`congrArg` 1 発**。

★★★§9-242 で「2–4 ブロック」と測ったが、
**射の等式を諦めて元の等式にした**ので **1 ブロック**で済んだ。
これは第 181(同型の合成でなく逆写像を直書き)と**同じ判断**である
——[[ring-instance-two-paths]] の「下流の要求に合わせる」。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `appIso_top_apply` | ★★元のレベルの `appIso` |
| `fromSpec_app_image_top` | ★★★app の輸送(自然性) |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} {U : X.Opens} (hU : IsAffineOpen U)

/-- ★元のレベルで同定する。 -/
theorem appIso_top_apply (s : (Γ(X, hU.fromSpec ''ᵁ (⊤ : (Spec Γ(X, U)).Opens)) : Type u)) :
    ((hU.fromSpec.appIso ⊤).hom).hom s
      = ((Spec Γ(X, U)).presheaf.map (eqToHom (Scheme.Hom.preimage_image_eq hU.fromSpec ⊤).symm).op).hom
        ((hU.fromSpec.app (hU.fromSpec ''ᵁ ⊤)).hom s) := by
  rw [Scheme.Hom.appIso_hom']
  rfl




/-- ★`fromSpec ''ᵁ ⊤` での app と `U` での app の関係(自然性そのもの)。 -/
theorem fromSpec_app_image_top
    (s : (Γ(X, hU.fromSpec ''ᵁ (⊤ : (Spec Γ(X, U)).Opens)) : Type u)) :
    (hU.fromSpec.app U).hom
        ((X.presheaf.map (eqToHom (fromSpec_image_top hU).symm).op).hom s)
      = ((Spec Γ(X, U)).presheaf.map
          ((Opens.map hU.fromSpec.base).map
            (eqToHom (fromSpec_image_top hU).symm).op.unop).op).hom
        ((hU.fromSpec.app (hU.fromSpec ''ᵁ ⊤)).hom s) :=
  congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) s)
    (hU.fromSpec.naturality (eqToHom (fromSpec_image_top hU).symm).op)


/-! ## ★出典の紐付け(`.src`) -/

def appIso_top_apply.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——appIso を元のレベルで開く)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
