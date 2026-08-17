import ABC3.Found.Arakelov.PicPushLax

/-!
# Arakelov (B1) 第 20 ブロック —— **引き戻しは oplax monoidal**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★連鎖が繋がった

    (restrictScalars α).LaxMonoidal        ★第 18 ブロック(mathlib に無かった)
    → (pushforward φ).LaxMonoidal          ★第 19 ブロック(strong ⋙ lax)
    → Adjunction.leftAdjointOplaxMonoidal  ★mathlib
    → **pullback が oplax monoidal**       ★本ブロック

★★これで比較射 `δ : f^*(L ⊗ M) ⟶ f^*L ⊗ f^*M` が得られる。
★★★残るのは**それが同型であること**——可逆層なら局所的に恒等なので、
第 15・16 ブロック(局所自明性の保存)が効く。

## ★実装の要点

★`(F.op ⋙ R) ⋙ forget₂` と `F.op ⋙ (R ⋙ forget₂)` は **defeq だが構文が違う**。
★★型注釈つきの別名(`alphaR`)を置くと `pullback` の暗黙引数が決まる。
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite

variable {C D : Type u} [Category.{u} C] [Category.{u} D] (F : C ⥤ D)
  (R : Dᵒᵖ ⥤ CommRingCat.{u}) (S : Cᵒᵖ ⥤ CommRingCat.{u}) (α : S ⟶ F.op ⋙ R)

/-- ★型を合わせた係数射(`RingCat` 側)。 -/
noncomputable abbrev alphaR : (S ⋙ forget₂ CommRingCat.{u} RingCat.{u})
    ⟶ F.op ⋙ (R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) :=
  Functor.whiskerRight α (forget₂ CommRingCat.{u} RingCat.{u})

/-- ★★★★★★**引き戻しは oplax monoidal である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで比較射 `δ : f^*(L ⊗ M) ⟶ f^*L ⊗ f^*M` が得られる。 -/
noncomputable instance pullbackOplax
    [(PresheafOfModules.pushforward.{u} (alphaR F R S α)).IsRightAdjoint] :
    (PresheafOfModules.pullback.{u} (alphaR F R S α)).OplaxMonoidal :=
  haveI : (PresheafOfModules.pushforward.{u} (alphaR F R S α)).LaxMonoidal :=
    pushCRLax F R S α
  (PresheafOfModules.pullbackPushforwardAdjunction (alphaR F R S α)).leftAdjointOplaxMonoidal

/-! ## ★出典の紐付け(`.src`) -/

def pullbackOplax.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——引き戻しが oplax monoidal であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
