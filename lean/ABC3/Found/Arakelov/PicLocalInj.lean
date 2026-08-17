import ABC3.Found.Arakelov.PicLocalSurj

/-!
# Arakelov (B1) 第 7 ブロック —— **局所単射性の基底と輸送**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★結合律の残り 1 条

第 6 ブロックで**局所全射**の側を取った(局所自由性を使わずに)。
残るのは**局所単射**の側である。

★★★一般の `M` では**偽**である(平坦性が要る)。しかし可逆層は
**局所的に構造層**なので、局所的には `M ≅ 𝟙_` となり成り立つ。

## ★本ブロックが取る 2 つ

| 定理 | 内容 |
|---|---|
| `isLocallyInjective_whiskerRight_unit` | ★**基底**: `M = 𝟙_` なら成り立つ |
| `isLocallyInjective_whiskerRight_iso` | ★**輸送**: `M ≅ N` なら移せる |

★★この 2 つを合わせると「局所的に `M ≅ 𝟙_`」から結論が出る
——残るのは**局所性の貼り合わせ**だけになる。

## ★機構

いずれも mathlib の

    isLocallyInjective_of_isLocallyInjective : IsLocallyInjective (φ ≫ ψ) → IsLocallyInjective φ

と、モノイダル圏の**自然性**(右単位子の自然性・交換律)の組み合わせである。
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite

variable {C : Type u} [Category.{u} C] (J : GrothendieckTopology C)
  {R : Cᵒᵖ ⥤ CommRingCat.{u}}
  {P Q : PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})}

/-- ★★★**基底の場合** —— 単位対象とのテンソルは局所単射性を保つ。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★機構は**右単位子の自然性** `(f ▷ 𝟙_) ≫ ρ_ = ρ_ ≫ f` と、
「合成が局所単射なら最初の射も局所単射」。 -/
theorem isLocallyInjective_whiskerRight_unit (f : P ⟶ Q)
    [PresheafOfModules.IsLocallyInjective J f] :
    PresheafOfModules.IsLocallyInjective J
      (f ▷ (𝟙_ (PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})))) := by
  have fac : (f ▷ (𝟙_ _)) ≫ (ρ_ Q).hom = (ρ_ P).hom ≫ f :=
    MonoidalCategory.rightUnitor_naturality f
  have h1 : PresheafOfModules.IsLocallyInjective J ((ρ_ P).hom ≫ f) := by
    dsimp [PresheafOfModules.IsLocallyInjective]
    rw [Functor.map_comp]
    infer_instance
  rw [← fac] at h1
  dsimp [PresheafOfModules.IsLocallyInjective] at h1 ⊢
  rw [Functor.map_comp] at h1
  exact Presheaf.isLocallyInjective_of_isLocallyInjective J _
    ((PresheafOfModules.toPresheaf _).map (ρ_ Q).hom)

/-- ★★★**輸送** —— 同型な対象へ移しても局所単射性は保たれる。

★機構は**交換律** `(f ▷ M) ≫ (Q ◁ e.hom) = (P ◁ e.hom) ≫ (f ▷ N)`。 -/
theorem isLocallyInjective_whiskerRight_iso (f : P ⟶ Q)
    {M N : PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})} (e : M ≅ N)
    [PresheafOfModules.IsLocallyInjective J (f ▷ N)] :
    PresheafOfModules.IsLocallyInjective J (f ▷ M) := by
  have fac : (f ▷ M) ≫ (Q ◁ e.hom) = (P ◁ e.hom) ≫ (f ▷ N) :=
    (whisker_exchange f e.hom).symm
  have h1 : PresheafOfModules.IsLocallyInjective J ((P ◁ e.hom) ≫ (f ▷ N)) := by
    dsimp [PresheafOfModules.IsLocallyInjective]
    rw [Functor.map_comp]
    infer_instance
  rw [← fac] at h1
  dsimp [PresheafOfModules.IsLocallyInjective] at h1 ⊢
  rw [Functor.map_comp] at h1
  exact Presheaf.isLocallyInjective_of_isLocallyInjective J _
    ((PresheafOfModules.toPresheaf _).map (Q ◁ e.hom))

/-! ## ★出典の紐付け(`.src`) -/

def isLocallyInjective_whiskerRight_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——単位対象とのテンソルが局所単射性を保つこと)",
    sectionId := "genell-def-1-1-i" }

def isLocallyInjective_whiskerRight_iso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所単射性の同型による輸送)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
