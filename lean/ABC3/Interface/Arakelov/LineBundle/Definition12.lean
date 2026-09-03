/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.AlgebraicGeometry.IdealSheaf.Functorial
import Mathlib.AlgebraicGeometry.Morphisms.Flat
import Mathlib.RingTheory.ClassGroup.Basic
import Mathlib.RingTheory.PicardGroup
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.AlgebraicGeometry.Modules.Sheaf
import Mathlib.Algebra.Category.ModuleCat.Presheaf.Monoidal
import Mathlib.Algebra.Category.ModuleCat.Presheaf.Sheafification
import ABC3.Interface.Arakelov.LineBundle.Definition11

/-!
# LineBundle —— `[GenEll] Definition 1.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Interface.Arakelov

open ABC3.Meta AlgebraicGeometry CategoryTheory NumberField

/-! ## ★★B3 —— `Spec 𝓞_F` 上での `Pic` -/

/-- **(B3)** `Pic(Spec 𝓞_F) ≅ ClassGroup 𝓞_F`。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★★★**高さの定義が `Spec 𝓞_F` へ引き戻した先で完結するのは、これがあるからである。**
★原文は `ADiv(F)/APrc(F) ⥲ APic(Spec 𝓞_F)` の形で使う——
その有限素点側が本 obligation である。 -/
structure PicSpecData where
  /-- 台となる `Pic`。 -/
  toPicardData : PicardData
  /-- `Pic(Spec 𝓞_F) ≃ ClassGroup 𝓞_F`。 -/
  equivClassGroup : (F : Type) → [Field F] → [NumberField F] →
    toPicardData.Pic (Spec (CommRingCat.of (𝓞 F))) ≃ ClassGroup (𝓞 F)
  /-- ★★★**その同型は (B1) の `equivPicRing` と整合する**。

  ★★★**これが自由な posit を殺す。**mathlib には
  `ClassGroup.equivPic : ClassGroup R ≃* CommRing.Pic R` があるので、
  本条件は `equivClassGroup` を**完全に決めてしまう**——
  すなわち **(B3) は独立の難所ではなく、(B1) の系である**。 -/
  equivClassGroup_compat : ∀ (F : Type) [Field F] [NumberField F]
    (L : toPicardData.Pic (Spec (CommRingCat.of (𝓞 F)))),
    ClassGroup.equivPic (𝓞 F) (equivClassGroup F L)
      = toPicardData.equivPicRing (CommRingCat.of (𝓞 F)) L

def PicSpecData.waiting : WaitingFor :=
  { what := "(B3) Pic(Spec 𝓞_F) ≅ ClassGroup 𝓞_F —— 数体の整数環上の可逆層が類群で書けること"
    trackB := "Found/Arakelov — ★`ClassGroup` は mathlib にある(`RingTheory/ClassGroup/Basic.lean`)。★★無いのは Pic 側(B1)と、その間の橋である。★Dedekind 環では可逆イデアル ≅ 可逆層なので、(B1) が入れば機械的である" }

def PicSpecData.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "Definition 1.2, (i)(層 B——Spec 𝓞_F 上の Pic)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Interface.Arakelov
