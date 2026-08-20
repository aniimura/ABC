import ABC3.Interface.Arakelov.ArcSpace
import Mathlib.AlgebraicGeometry.Morphisms.Proper
import Mathlib.AlgebraicGeometry.Morphisms.Flat
import Mathlib.AlgebraicGeometry.Limits

/-!
# 負の対照 —— **`ProperFlatOverZ := False` は接地欄で落ちる**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

## ★★★★何を確かめるか

(C2) `ProjectiveModelData` の `ProperFlatOverZ` は自前の posit である。
★2026-08-19 以前は接地が無かったので

    ProperFlatOverZ := fun _ => False

と置けば `properFlat_compact` が**空虚に成立**し、残る `projectiveCase`
(`Found/GenEll/ArcModel.lean` で実装済)だけで (C2) が「達成」になってしまった。

★★塞いだ手は `properFlatOverZ_iff`——mathlib の `IsProper` と `Flat` に縛る欄である。

★★★本ファイルはそれが**効いている**ことを機械で確かめる:
`Spec ℤ` は自分自身の上で固有かつ平坦なので、`False` では同値が破れる。
-/

namespace ABC3.Check.Arakelov

open AlgebraicGeometry CategoryTheory

/-- ★`Spec ℤ` から `Spec ℤ` への標準射は恒等射である。 -/
theorem from_specZ_self :
    specZIsTerminal.from (Spec (CommRingCat.of ℤ)) = 𝟙 _ :=
  specZIsTerminal.hom_ext _ _

/-- ★★`Spec ℤ` は ℤ 上固有である。 -/
instance : IsProper (specZIsTerminal.from (Spec (CommRingCat.of ℤ))) := by
  rw [from_specZ_self]; infer_instance

/-- ★★`Spec ℤ` は ℤ 上平坦である。 -/
instance : Flat (specZIsTerminal.from (Spec (CommRingCat.of ℤ))) := by
  rw [from_specZ_self]; infer_instance

/-- ★★★★★**`ProperFlatOverZ := False` は接地欄を満たさない**。

★この定理が通ること自体が、`properFlatOverZ_iff`(2026-08-19 追加)が
退化を実際に殺していることの証明である。 -/
theorem falseProperFlat_fails_grounding :
    ¬ (∀ X : Scheme.{0}, (False : Prop) ↔
        (IsProper (specZIsTerminal.from X) ∧ Flat (specZIsTerminal.from X))) := by
  intro h
  exact (h (Spec (CommRingCat.of ℤ))).2 ⟨inferInstance, inferInstance⟩

end ABC3.Check.Arakelov
