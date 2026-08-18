import ABC3.Found.Arakelov.PicBCSquare

/-!
# Arakelov (B1) 第 46 ブロック —— **Beck–Chevalley の mate**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★mate ができた

    (f|)^*_pre ∘ restrict_V  ⟹  restrict_{f⁻¹V} ∘ f^*_pre

★★随伴の下でこれは

    restrict_V  ⟹  restrict_V ∘ f_* ∘ f^*_pre

に対応し、右辺は `restrict_V` を随伴の **unit** に当てたものである。
★★★第 45 ブロック(押し出しの四角形が `rfl` で可換)がこれを可能にした。

## ★★実装の要点(2026-08-18 実測)

★`rw` は透明度で落ちる——`erw` が要る箇所が 2 つある。
★★最後の 1 歩(`map_comp` の往復)は `rw` では当たらないので、
**項として書く**(`map_comp.symm.trans (congrArg _ ...).trans map_comp`)。

## ★★残り

    第 47: mate が生成元の上で同型(第 24 + 第 23 `opensMap_inf`)
    第 48: 第 29 の器具で全体へ持ち上げる
    第 49: 第 41 の仮定を落とす
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (V : Y.Opens)

/-- ★制限した随伴 `(f|)^*_pre ⊣ (f|)_*`。 -/
noncomputable abbrev adjOn :=
  PresheafOfModules.pullbackPushforwardAdjunction
    (alphaR (overPost f V)
      ((Over.forget ((Opens.map f.base).obj V)).op ⋙ X.presheaf)
      ((Over.forget V).op ⋙ Y.presheaf)
      (restrictedC f V))

/-- ★★★★★★**Beck–Chevalley の mate**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが同型であることが、第 41 ブロックの仮定(局所階数 1)を
落とすための鍵である。

★機構は「随伴の下で unit に対応する」——第 45 ブロックの
四角形の可換性(`rfl`)がそれを可能にした。 -/
noncomputable def bcMate :
    restrictPresheafFunctor Y V ⋙ pullbackPreOn f V
      ⟶ pullbackPre f ⋙ restrictPresheafFunctor X ((Opens.map f.base).obj V) where
  app P := ((adjOn f V).homEquiv _ _).symm
    ((restrictPresheafFunctor Y V).map
      ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app P))
  naturality := by
    intro P P' g
    dsimp
    erw [← Adjunction.homEquiv_naturality_left_symm,
      ← Adjunction.homEquiv_naturality_right_symm]
    congr 1
    have hsq : ∀ {A B : X.PresheafOfModules} (h : A ⟶ B),
        (pushOn f V).map ((restrictPresheafFunctor X ((Opens.map f.base).obj V)).map h)
          = (restrictPresheafFunctor Y V).map
            ((PresheafOfModules.pushforward (pullbackPhi f)).map h) := fun _ => rfl
    erw [hsq]
    exact ((restrictPresheafFunctor Y V).map_comp _ _).symm.trans
      ((congrArg (restrictPresheafFunctor Y V).map
        ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.naturality g)).trans
        ((restrictPresheafFunctor Y V).map_comp _ _))

/-! ## ★出典の紐付け(`.src`) -/

def bcMate.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Beck–Chevalley の mate)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
