import ABC3.Found.Arakelov.PicWhiskerW
import ABC3.Found.Arakelov.PicSheafTensor
import Mathlib.Algebra.Category.ModuleCat.Sheaf.Localization

/-!
# Arakelov (B1) 第 11 ブロック —— **層化が `f ▷ M` を同型に送る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★`W` から `IsIso` へ

第 10 ブロックで `W (f ▷ M)` を取った。★★mathlib は

    (sheafification α).IsLocalization (J.W.inverseImage (toPresheaf R₀))

を持つ(`ModuleCat/Sheaf/Localization.lean`)ので、
`Localization.inverts` で **`IsIso` に変換できる**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isIso_sheafify_whiskerRight` | 層化は `f ▷ M` を同型に送る |

★★★これで結合子を組む材料が揃う。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

variable (X : Scheme.{u})

/-- ★★★★**層化は `f ▷ M` を同型に送る**(`W` から `IsIso` へ)。 -/
theorem isIso_sheafify_whiskerRight
    {P Q : PresheafOfModules.{u} (X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})}
    (f : P ⟶ Q)
    [PresheafOfModules.IsLocallyInjective (Opens.grothendieckTopology X) f]
    [PresheafOfModules.IsLocallySurjective (Opens.grothendieckTopology X) f]
    (M : PresheafOfModules.{u} (X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}))
    (htriv : ∀ U : X.Opens, ∃ S : Sieve U, S ∈ (Opens.grothendieckTopology X) U ∧
      ∀ ⦃V : X.Opens⦄ (i : V ⟶ U), S i →
        Nonempty ((X.presheaf.obj (op V) : Type u) ≃ₗ[(X.presheaf.obj (op V) : Type u)]
          M.obj (op V))) :
    IsIso ((PresheafOfModules.sheafification (R := X.ringCatSheaf)
      (𝟙 X.ringCatSheaf.obj)).map (f ▷ M)) :=
  Localization.inverts _
    ((Opens.grothendieckTopology X).W.inverseImage
      (PresheafOfModules.toPresheaf X.ringCatSheaf.obj)) _
    (W_whiskerRight_of_basis (Opens.grothendieckTopology X) f M htriv)

/-! ## ★★左からのテンソル(対称性で出る)

★`M ◁ f = (β_ M P).hom ≫ (f ▷ M) ≫ (β_ Q M).hom` なので、
対称モノイダル構造の braiding で右からの場合に帰着する。
★★★ただし `rw` が `sheafification` の `map` に噛まないため
(`@[simps! -isSimp map]` で定義されている)、
**関手の合成保存を明示的に当てる必要がある**。次のブロックで行う。 -/

/-! ## ★出典の紐付け(`.src`) -/

def isIso_sheafify_whiskerRight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層化が f ▷ M を同型に送ること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
