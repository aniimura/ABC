import ABC3.Found.Arakelov.PicPresheafTensor
import Mathlib.Algebra.Category.ModuleCat.Presheaf.Sheafification

/-!
# Arakelov (B1) 第 2 ブロック —— **層加群のテンソル積**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★層化で `⊗` を層の側へ渡す

第 1 ブロックで `X.PresheafOfModules` に対称モノイダル構造を渡した。
★★本ブロックは **層化関手**でそれを `X.Modules` へ落とす:

    M ⊗ N := 層化 (M.val ⊗ N.val)

★★★**2026-08-17 の実測で、層化関手はスキーム上で取れることが分かった**
——`PresheafOfModules.sheafification (𝟙 X.ringCatSheaf.obj)` が
必要なインスタンス(局所単射・局所全射・`WEqualsLocallyBijective`・
`HasWeakSheafify`)をすべて自動で見つける。

## ★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `sheafifyFunctor` | 層化関手 `X.PresheafOfModules ⥤ X.Modules` |
| `sheafifyAdjunction` | その随伴(普遍性) |
| `tensorModules` | ★**層加群のテンソル積** |
| `tensorModulesComm` | ★★**可換**(`Pic X` の `CommGroup` のもと) |
| `sheafifyValIso` | ★層をもう一度層化しても変わらない(counit が iso) |

## ★★残る難所

結合律 `sheafify (sheafify (A ⊗ B) ⊗ C) ≅ sheafify (A ⊗ B ⊗ C)` は
「層化がテンソル積と可換」を要し、これは `Localization/Monoidal` の
`W.IsMonoidal`(局所同型がテンソルで安定)に帰着する。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory

/-! ## ★★層化関手 -/

/-- ★★**スキーム上の前層加群の層化**。

★★★mathlib の `PresheafOfModules.sheafification` を恒等係数射に当てる。
`α = 𝟙` なので局所単射・局所全射は自動である。 -/
noncomputable abbrev sheafifyFunctor (X : Scheme.{u}) :
    X.PresheafOfModules ⥤ X.Modules :=
  PresheafOfModules.sheafification (R := X.ringCatSheaf) (𝟙 X.ringCatSheaf.obj)

/-- ★**層化は忘却の左随伴である**(層化の普遍性)。 -/
noncomputable abbrev sheafifyAdjunction (X : Scheme.{u}) :
    sheafifyFunctor X ⊣
      SheafOfModules.forget X.ringCatSheaf ⋙
        PresheafOfModules.restrictScalars (𝟙 X.ringCatSheaf.obj) :=
  PresheafOfModules.sheafificationAdjunction (𝟙 X.ringCatSheaf.obj)

/-- ★★**層をもう一度層化しても変わらない**(counit が同型)。

★★★これが「`M` が既に層なら層化は何もしない」という当たり前を型にしたもので、
単位律 `𝒪 ⊗ M ≅ M` を出すときに効く。 -/
noncomputable def sheafifyValIso {X : Scheme.{u}} (M : X.Modules) :
    (sheafifyFunctor X).obj
        ((SheafOfModules.forget X.ringCatSheaf ⋙
          PresheafOfModules.restrictScalars (𝟙 X.ringCatSheaf.obj)).obj M) ≅ M :=
  (asIso (PresheafOfModules.sheafificationAdjunction
    (R := X.ringCatSheaf) (𝟙 X.ringCatSheaf.obj)).counit).app M

/-! ## ★★★層加群のテンソル積 -/

variable {X : Scheme.{u}}

/-- ★★★**層加群のテンソル積**——前層でテンソルしてから層化する。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★**これが `Pic X` の乗法である。** -/
noncomputable def tensorModules (M N : X.Modules) : X.Modules :=
  (sheafifyFunctor X).obj (M.val ⊗ N.val)

/-- ★テンソル積は射に沿って関手的である。 -/
noncomputable def tensorModulesMap {M₁ M₂ N₁ N₂ : X.Modules}
    (f : M₁ ⟶ M₂) (g : N₁ ⟶ N₂) : tensorModules M₁ N₁ ⟶ tensorModules M₂ N₂ :=
  (sheafifyFunctor X).map ((SheafOfModules.forget X.ringCatSheaf).map f ⊗ₘ
    (SheafOfModules.forget X.ringCatSheaf).map g)

/-- ★同型どうしのテンソル積は同型である(`Pic` が well-defined になる根拠)。 -/
noncomputable def tensorModulesIso {M₁ M₂ N₁ N₂ : X.Modules}
    (e : M₁ ≅ M₂) (f : N₁ ≅ N₂) : tensorModules M₁ N₁ ≅ tensorModules M₂ N₂ :=
  (sheafifyFunctor X).mapIso
    ((SheafOfModules.forget X.ringCatSheaf).mapIso e ⊗ᵢ
      (SheafOfModules.forget X.ringCatSheaf).mapIso f)

/-- ★★★**テンソル積は可換である**。

★★★これが `Pic X` の **`CommGroup`**(可換群)の可換性を出す。
★機構は前層の**対称**モノイダル構造(第 1 ブロック)の braiding を層化で送るだけ。 -/
noncomputable def tensorModulesComm (M N : X.Modules) :
    tensorModules M N ≅ tensorModules N M :=
  (sheafifyFunctor X).mapIso (β_ M.val N.val)

/-! ## ★★★単位律 -/

/-- ★**構造層それ自身を層加群と見たもの**——テンソル積の単位。 -/
noncomputable abbrev unitModules (X : Scheme.{u}) : X.Modules :=
  SheafOfModules.unit X.ringCatSheaf

/-- ★その台の前層は、前層モノイダル構造の単位そのものである。 -/
theorem unitModules_val (X : Scheme.{u}) :
    (unitModules X).val = 𝟙_ X.PresheafOfModules := rfl

/-- ★★★**単位律(左)** `𝒪_X ⊗ M ≅ M`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `Pic X` の**単位元**を与える。
★機構は前層の左単位子 `λ_` を層化で送り、`sheafifyValIso`(層は層化で変わらない)で閉じる。 -/
noncomputable def tensorUnitLeft (M : X.Modules) :
    tensorModules (unitModules X) M ≅ M :=
  (sheafifyFunctor X).mapIso (λ_ M.val) ≪≫ sheafifyValIso M

/-- ★★★**単位律(右)** `M ⊗ 𝒪_X ≅ M`。 -/
noncomputable def tensorUnitRight (M : X.Modules) :
    tensorModules M (unitModules X) ≅ M :=
  (sheafifyFunctor X).mapIso (ρ_ M.val) ≪≫ sheafifyValIso M

/-! ## ★出典の紐付け(`.src`) -/

def sheafifyFunctor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——前層加群の層化関手)",
    sectionId := "genell-def-1-1-i" }

def tensorModules.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層加群のテンソル積)",
    sectionId := "genell-def-1-1-i" }

def tensorUnitLeft.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——テンソル積の単位律)",
    sectionId := "genell-def-1-1-i" }

def tensorModulesComm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——テンソル積の可換性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
