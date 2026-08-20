import ABC3.Found.Arakelov.PicBijGood

/-!
# Arakelov (B2) 第 239 ブロック —— ★★★★★★**比較射の層化は同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★良い開集合は**点ごとに取れる**

`x ∈ U` に対し、次の順で絞る:

| 段 | 取るもの | 根拠 |
|---|---|---|
| 1 | `Y` のアフィン開 `B ∋ f(x)` | ★`Y.isBasis_affineOpens` |
| 2 | `W₁ ∋ x`、`W₁ ≤ U ⊓ f⁻¹ᵁB`、引き戻し側が自明 | ★第 59 + 第 237 |
| 3 | `W₂ ∋ x`、`W₂ ≤ W₁`、逆像因子側が自明 | ★第 202(平坦)+ 第 237 |
| 4 | `V ∋ x`、`V ≤ W₂`、`V` はアフィン | ★`X.isBasis_affineOpens` |

★★`V` は 4 条件をすべて満たす——自明性は第 237 で `W₁, W₂` から降ろす。

★★★★これで第 236(述語版の局所全単射判定)に渡せる。

## ★★★★★★層化は局所同型を同型にする

`W_of_isLocallyBijective` + `inverseImage_W_toPresheaf_eq_inverseImage_isomorphisms`
——第 190 ブロック(積の場合)と**同じ 3 行**である。

| 定理 | 内容 |
|---|---|
| `exists_good` | ★★★★良い開集合は点ごとに取れる |
| `isIso_sheafify_pullIdealHom` | ★★★★★★**比較射の層化は同型** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) [AlgebraicGeometry.Flat f] (D : Y.IdealSheafData)

/-- ★★★★**良い開集合は点ごとに取れる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★平坦性は第 202(`isCartier_comap`)でしか使わない。 -/
theorem exists_good (hD : IsCartier Y D) (U : X.Opens) (x : X) (hx : x ∈ U) :
    ∃ V : X.Opens, x ∈ V ∧ V ≤ U ∧
      Function.Bijective (((pullIdealHom f D).app (op V)).hom) := by
  obtain ⟨B, hB, hfxB, -⟩ := Opens.isBasis_iff_nbhd.1 Y.isBasis_affineOpens
    (U := ⊤) (x := f.base x) trivial
  have hsrc : IsLocallyTrivial X ((pullbackPre f).obj (idealPresheaf D)) :=
    isLocallyTrivial_pullbackPre f (idealPresheaf D) (isLocallyTrivial_idealSheaf D hD)
  have htgt : IsLocallyTrivial X (idealPresheaf (D.comap f)) :=
    isLocallyTrivial_idealSheaf (D.comap f) (isCartier_comap f D hD)
  have hxUB : x ∈ U ⊓ f ⁻¹ᵁ B := ⟨hx, hfxB⟩
  obtain ⟨W₁, hxW₁, hW₁, ⟨e₁⟩⟩ := exists_trivial_nbhd hsrc (U ⊓ f ⁻¹ᵁ B) x hxUB
  obtain ⟨W₂, hxW₂, hW₂, ⟨e₂⟩⟩ := exists_trivial_nbhd htgt W₁ x hxW₁
  obtain ⟨V, hVaff, hxV, hVW₂⟩ :=
    Opens.isBasis_iff_nbhd.1 X.isBasis_affineOpens (U := W₂) (x := x) hxW₂
  refine ⟨V, hxV, le_trans hVW₂ (le_trans hW₂ (le_trans hW₁ inf_le_left)), ?_⟩
  exact bij_of_good f D hVaff hB
    (le_trans hVW₂ (le_trans hW₂ (le_trans hW₁ inf_le_right)))
    (restrict_trivial_of_le (le_trans hVW₂ hW₂) e₁)
    (restrict_trivial_of_le hVW₂ e₂)

/-- ★★★★★★**比較射の層化は同型である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 190 ブロック(積の場合)と同じ形——述語が「アフィン」から
「良い開集合」に変わっただけである。 -/
theorem isIso_sheafify_pullIdealHom (hD : IsCartier Y D) :
    IsIso ((sheafifyFunctor X).map (pullIdealHom f D)) := by
  obtain ⟨hi, hs⟩ := locallyBijective_of_bijective_on_pred
    ((PresheafOfModules.toPresheaf _).map (pullIdealHom f D))
    (fun V => Function.Bijective (((pullIdealHom f D).app (op V)).hom))
    (fun U x hx => by
      obtain ⟨V, hxV, hVU, hbij⟩ := exists_good f D hD U x hx
      exact ⟨V, hbij, hxV, hVU⟩)
    (fun _ hA => hA)
  haveI hi' : Presheaf.IsLocallyInjective (Opens.grothendieckTopology X)
      ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (pullIdealHom f D)) := hi
  haveI hs' : Presheaf.IsLocallySurjective (Opens.grothendieckTopology X)
      ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (pullIdealHom f D)) := hs
  have hW : (MorphismProperty.inverseImage (Opens.grothendieckTopology X).W
      (PresheafOfModules.toPresheaf X.ringCatSheaf.obj)) (pullIdealHom f D) :=
    GrothendieckTopology.W_of_isLocallyBijective _ _
  have heq := PresheafOfModules.inverseImage_W_toPresheaf_eq_inverseImage_isomorphisms
    (R := X.ringCatSheaf) (𝟙 X.ringCatSheaf.obj)
  have h2 : (MorphismProperty.isomorphisms (SheafOfModules X.ringCatSheaf)).inverseImage
      (PresheafOfModules.sheafification (𝟙 X.ringCatSheaf.obj)) (pullIdealHom f D) := heq ▸ hW
  exact h2

/-! ## ★出典の紐付け(`.src`) -/

def isIso_sheafify_pullIdealHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——比較射の層化は同型)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
