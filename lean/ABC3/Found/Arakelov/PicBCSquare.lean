import ABC3.Found.Arakelov.PicRestrictedPull

/-!
# Arakelov (B1) 第 45 ブロック —— **Beck–Chevalley の四角形は可換である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★可換性は `rfl` で出た

Beck–Chevalley の mate を作るには、**押し出しの側の四角形が可換**であることが要る:

    restrict_V ∘ f_*  =  (f|)_* ∘ restrict_{f⁻¹V}

★★★どちらも **precomposition** であり、site の四角形

    Over.forget V ⋙ Opens.map f.base = overPost f V ⋙ Over.forget (f⁻¹V)

が `rfl` で可換だから、**両辺は同じ関手である**(2026-08-18 実測)。

## ★★これで mate が作れる

`(f|)^*_pre` は `(f|)_*` の左随伴なので、

    (f|)^*_pre ∘ restrict_V  ⟹  restrict_{f⁻¹V} ∘ f^*_pre

は随伴の下で

    restrict_V  ⟹  (f|)_* ∘ restrict_{f⁻¹V} ∘ f^*_pre  =  restrict_V ∘ f_* ∘ f^*_pre

に対応する。★★★右辺は `restrict_V` を随伴の **unit** に当てたものである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `overSquare_comm` | ★site の四角形が可換(`rfl`) |
| `pushOn` | ★制限した押し出し `(f|)_*` |
| `pushSquare_obj` | ★★**押し出しの四角形が可換**(`rfl`) |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (V : Y.Opens)

/-! ## ★site の四角形 -/

/-- ★**site の四角形は可換である**——どちらも「`W ≤ V` を `f⁻¹W` に送る」。 -/
theorem overSquare_comm :
    Over.forget V ⋙ Opens.map f.base
      = overPost f V ⋙ Over.forget ((Opens.map f.base).obj V) := rfl

/-! ## ★制限した押し出し -/

/-- ★**制限した押し出し `(f|)_*`**。 -/
noncomputable abbrev pushOn :
    PresheafModulesOn X ((Opens.map f.base).obj V) ⥤ PresheafModulesOn Y V :=
  PresheafOfModules.pushforward
    (alphaR (overPost f V)
      ((Over.forget ((Opens.map f.base).obj V)).op ⋙ X.presheaf)
      ((Over.forget V).op ⋙ Y.presheaf)
      (restrictedC f V))

/-! ## ★★押し出しの四角形 -/

/-- ★★★**押し出しの四角形は可換である**(対象の上で `rfl`)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★どちらも precomposition であり、site の四角形が `rfl` で可換だから
**同じものになる**。★これが Beck–Chevalley の mate を作る土台である。 -/
theorem pushSquare_obj (M : X.PresheafOfModules) :
    (restrictPresheafFunctor Y V).obj
        ((PresheafOfModules.pushforward (pullbackPhi f)).obj M)
      = (pushOn f V).obj
        ((restrictPresheafFunctor X ((Opens.map f.base).obj V)).obj M) := rfl

/-! ## ★出典の紐付け(`.src`) -/

def overSquare_comm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——site の四角形が可換であること)",
    sectionId := "genell-def-1-1-i" }

def pushSquare_obj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——押し出しの四角形が可換であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
