import ABC3.Found.Arakelov.ArcAwayLift
import ABC3.Found.Arakelov.ArcLandsInScheme

/-!
# Arakelov (C1) 段 B —— **開部分スキームへの合成は開写像である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★段 A から段 B への道

段 A(`arcTopology_opens_of_affine`)は**アフィン周囲**での chart 埋め込みを与えた。
★★本ファイルはそれを**一般のスキーム**へ上げる。

| 手順 | 内容 |
|---|---|
| 1 | `(· ≫ O.ι)` の像が開(**一般の `Y` で**——`isOpen_landsIn_scheme` の系) |
| 2 | 段 A を `IsOpenMap` に翻訳 |
| 3 | 同型による輸送 |
| 4 | 任意のアフィンスキームへ(`isoSpec` で輸送) |
| 5 | ★★★一般の `Y` へ(`U ⊓ O` を経由) |

## ★★★手順 5 の機構

`arcTopology Y = ⨆`(アフィン chart)なので、開性は chart ごとに確かめればよい。
★★chart `U` での逆像は

    (· ≫ U.ι) ⁻¹' ((· ≫ O.ι) '' V) = (· ≫ homOfLE) '' ((· ≫ homOfLE) ⁻¹' V)

(両辺 `U ⊓ O` を経由)となり、右辺は**アフィンな `U` の上での段 A**で開である。
-/

universe u v

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory TopologicalSpace

/-! ## ★★手順 1 —— 像は開である(一般の `Y`) -/

/-- ★★★**`(· ≫ O.ι)` の像は `Y^arc` の開集合である**(`O` は任意の開集合)。

★`isOpen_range_chart` のアフィン性を落とした版である
——証明はアフィン性を使っていなかった。 -/
theorem isOpen_range_comp_ι {Y : Scheme.{0}} (O : Y.Opens) :
    @IsOpen _ (arcTopology Y)
      (Set.range (fun p : Spec (CommRingCat.of ℂ) ⟶ O.toScheme => p ≫ O.ι)) := by
  letI := arcTopology Y
  have hset : Set.range (fun p : Spec (CommRingCat.of ℂ) ⟶ O.toScheme => p ≫ O.ι)
      = {q : Spec (CommRingCat.of ℂ) ⟶ Y | q.base default ∈ O} := by
    ext q
    rw [mem_range_comp_iff_base_default O.ι q, Set.mem_setOf_eq, Scheme.Opens.range_ι]
    exact Iff.rfl
  rw [hset]
  exact isOpen_landsIn_scheme O

/-! ## ★★手順 2 —— 段 A を `IsOpenMap` に翻訳 -/

/-- ★**誘導位相 + 像が開 ⟹ 開写像**(一般位相の補助)。

★★単射性は要らない——`g '' (g ⁻¹' W) = W ∩ range g` は常に成り立つ。 -/
theorem isOpenMap_of_induced {α β : Type v} (tα : TopologicalSpace α)
    (tβ : TopologicalSpace β) {g : α → β} (hind : tα = TopologicalSpace.induced g tβ)
    (hopen : @IsOpen _ tβ (Set.range g)) : @IsOpenMap _ _ tα tβ g := by
  letI := tβ
  intro V hV
  rw [hind, isOpen_induced_iff] at hV
  obtain ⟨W, hW, rfl⟩ := hV
  rw [Set.image_preimage_eq_inter_range]
  exact hW.inter hopen

/-- ★★★**アフィン周囲では `(· ≫ O.ι)` は開写像である**(段 A の翻訳)。 -/
theorem isOpenMap_comp_ι_of_affine (A : CommRingCat.{0}) (O : (Spec A).Opens) :
    @IsOpenMap _ _ (arcTopology O.toScheme) (arcTopology (Spec A))
      (fun p : Spec (CommRingCat.of ℂ) ⟶ O.toScheme => p ≫ O.ι) :=
  isOpenMap_of_induced _ _ (arcTopology_opens_of_affine A O) (isOpen_range_comp_ι O)

/-! ## ★★手順 3 —— 同型による輸送 -/

/-- ★**同型との合成は開写像**(逆写像が連続だから)。 -/
theorem isOpenMap_comp_iso_hom {X Y : Scheme.{0}} (e : X ≅ Y) :
    @IsOpenMap _ _ (arcTopology X) (arcTopology Y)
      (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ e.hom) := by
  letI := arcTopology X
  letI := arcTopology Y
  intro V hV
  have hset : (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ e.hom) '' V
      = (fun q : Spec (CommRingCat.of ℂ) ⟶ Y => q ≫ e.inv) ⁻¹' V := by
    ext q
    constructor
    · rintro ⟨p, hp, rfl⟩
      show (p ≫ e.hom) ≫ e.inv ∈ V
      rw [Category.assoc, Iso.hom_inv_id, Category.comp_id]
      exact hp
    · intro hq
      refine ⟨q ≫ e.inv, hq, ?_⟩
      show (q ≫ e.inv) ≫ e.hom = q
      rw [Category.assoc, Iso.inv_hom_id, Category.comp_id]
  rw [hset]
  exact (continuous_comp_openImmersion e.inv).isOpen_preimage _ hV

/-- ★**開埋め込みは自身の像への同型を与える**。 -/
noncomputable def opensRangeIso {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f] :
    X ≅ f.opensRange.toScheme :=
  IsOpenImmersion.isoOfRangeEq f f.opensRange.ι (by rw [Scheme.Opens.range_ι]; rfl)

@[simp] theorem opensRangeIso_fac {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f] :
    (opensRangeIso f).hom ≫ f.opensRange.ι = f :=
  IsOpenImmersion.isoOfRangeEq_hom_fac f f.opensRange.ι _

/-! ## ★★手順 4 —— 任意のアフィンスキームへ -/

/-- ★★**`Spec A` への任意の開埋め込みは開写像**。 -/
theorem isOpenMap_comp_of_affine {A : CommRingCat.{0}} {X : Scheme.{0}} (g : X ⟶ Spec A)
    [IsOpenImmersion g] :
    @IsOpenMap _ _ (arcTopology X) (arcTopology (Spec A))
      (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ g) := by
  letI := arcTopology X
  letI := arcTopology (Spec A)
  letI := arcTopology g.opensRange.toScheme
  have hcomp : (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ g)
      = (fun q : Spec (CommRingCat.of ℂ) ⟶ g.opensRange.toScheme => q ≫ g.opensRange.ι)
        ∘ (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ (opensRangeIso g).hom) := by
    funext p
    show p ≫ g = (p ≫ (opensRangeIso g).hom) ≫ g.opensRange.ι
    rw [Category.assoc, opensRangeIso_fac]
  rw [hcomp]
  exact (isOpenMap_comp_ι_of_affine A g.opensRange).comp
    (isOpenMap_comp_iso_hom (opensRangeIso g))

/-- ★★★**アフィンスキームへの任意の開埋め込みは開写像**。

★`Spec A` に限らない——`isoSpec` で輸送する。 -/
theorem isOpenMap_comp_of_isAffine {W X : Scheme.{0}} [IsAffine W] (g : X ⟶ W)
    [IsOpenImmersion g] :
    @IsOpenMap _ _ (arcTopology X) (arcTopology W)
      (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ g) := by
  letI := arcTopology X
  letI := arcTopology W
  letI := arcTopology (Spec (W.presheaf.obj (Opposite.op ⊤)))
  haveI : IsOpenImmersion (g ≫ W.isoSpec.hom) := inferInstance
  have hcomp : (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ g)
      = (fun q : Spec (CommRingCat.of ℂ) ⟶ Spec (W.presheaf.obj (Opposite.op ⊤)) =>
          q ≫ W.isoSpec.inv)
        ∘ (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ (g ≫ W.isoSpec.hom)) := by
    funext p
    show p ≫ g = (p ≫ (g ≫ W.isoSpec.hom)) ≫ W.isoSpec.inv
    rw [Category.assoc, Category.assoc, Iso.hom_inv_id, Category.comp_id]
  rw [hcomp]
  exact (isOpenMap_comp_iso_hom W.isoSpec.symm).comp
    (isOpenMap_comp_of_affine (g ≫ W.isoSpec.hom))

/-! ## ★出典の紐付け(`.src`) -/

def isOpen_range_comp_ι.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——一般のスキームで chart の像が開であること)",
    sectionId := "genell-def-1-1-i" }

def isOpenMap_comp_of_isAffine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィンスキームへの開埋め込みは開写像)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
