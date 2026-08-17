import ABC3.Found.Arakelov.PicResScalarsLax

/-!
# Arakelov (B1) 第 19 ブロック —— **押し出しの lax monoidal**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★合成で出る

    pushforward = pushforward₀(同じ環、★mathlib で **strong** monoidal)
                ⋙ restrictScalars(★第 18 ブロックで **lax** monoidal)

★★lax monoidal 関手の合成は lax monoidal なので、**インスタンス探索が自動で見つける**。

## ★★★これで引き戻しへ届く

    (pushforward φ).LaxMonoidal            ★本ブロック
    → Adjunction.leftAdjointOplaxMonoidal  ★mathlib
    → pullback が oplax monoidal
    → f^*(L ⊗ M) ≅ f^*L ⊗ f^*M
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite

variable {C D : Type u} [Category.{u} C] [Category.{u} D] (F : C ⥤ D)
  (R : Dᵒᵖ ⥤ CommRingCat.{u}) (S : Cᵒᵖ ⥤ CommRingCat.{u}) (α : S ⟶ F.op ⋙ R)

/-- ★★**押し出し(`CommRingCat` 版)** = 押し出し(同じ環)⋙ 係数変換。 -/
noncomputable abbrev pushCR :
    PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) ⥤
      PresheafOfModules.{u} (S ⋙ forget₂ CommRingCat.{u} RingCat.{u}) :=
  PresheafOfModules.pushforward₀OfCommRingCat F R ⋙ resScalars α

/-- ★★★★**押し出しは lax monoidal である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが引き戻し `f^*` のモノイダル性への入口である。
★機構は「strong ⋙ lax = lax」——インスタンス探索が自動で見つける。 -/
noncomputable instance pushCRLax : (pushCR F R S α).LaxMonoidal := inferInstance

/-! ## ★出典の紐付け(`.src`) -/

def pushCRLax.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——押し出しの lax monoidal)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
