import ABC3.Meta.Claim
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs
import Mathlib.Algebra.Group.Units.Defs
import Mathlib.Tactic.Group

/-!
# [LocProP] §4 —— J-Geometric Sections の posit(`Interface`)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、物理 p.30-31。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

`Lemma 4.1`(`J_X^{(1)}(K)` の幾何点との一致 ⟺ `δ_θ = 0`)・`Lemma 4.2`(共役公式)・
`Definition 4.3`(`J-geometric` の命名)・`Proposition 4.4`(`(∗)^spl` ⟺ `J-geometric`)は
§1-§3 で posit した Malčev completion・weight zero quotient の理論の上に立つ。
★★引き続き mathlib に丸ごと無いので posit する。
-/

namespace ABC3.Interface.LocProP

universe u

/-- ★posit —— §4(J-Geometric Sections)の骨組み。 -/
structure JGeometricSetup where
  /-- `Π_{X_K} → Γ_K` の節 `θ` 全体。 -/
  Theta : Type u
  /-- `J_X^{(1)}(K)` の幾何点。 -/
  Lpt : Type u
  /-- `L ↦ θ_L`。 -/
  thetaOf : Lpt → Theta
  H1 : Type u
  H1Grp : AddCommGroup H1
  /-- `δ_θ ∈ H¹(K, 𝔷_X[0])`。 -/
  delta : Theta → H1
  /-- **[LocProP] Lemma 4.1**。 -/
  lemma_4_1 : ∀ θ, (∃ l, thetaOf l = θ) ↔ delta θ = 0
  /-- `Z_X(K̂̂)`。 -/
  Zhat : Type u
  ZhatGrp : Group Zhat
  GammaK : Type u
  /-- `α` が定める `Γ_K` の `Z_X(K̂̂)` への作用。 -/
  actionAlpha : GammaK → Zhat → Zhat
  /-- `β` が定める作用。 -/
  actionBeta : GammaK → Zhat → Zhat
  /-- `ζ(γ)_Z`。 -/
  zetaZ : GammaK → Zhat
  /-- **[LocProP] Lemma 4.2**。 -/
  lemma_4_2 : ∀ γ φ, actionAlpha γ φ = zetaZ γ * actionBeta γ φ * (zetaZ γ)⁻¹
  /-- `Π_{X_K} → Γ_K` の節 `α'` 全体。 -/
  AlphaPrime : Type u
  /-- **[LocProP] Definition 4.3** —— `α'` が J-geometric であること(posit)。 -/
  isJGeometric : AlphaPrime → Prop
  /-- 条件 `(∗)^spl`(posit)。 -/
  satisfiesSpl : AlphaPrime → Prop
  /-- **[LocProP] Proposition 4.4**。 -/
  prop_4_4 : ∀ a, satisfiesSpl a ↔ isJGeometric a

/-- ★★**非退化性の witness**(具体項)。 -/
@[reducible] def JGeometricSetup.example : JGeometricSetup.{0} where
  Theta := ℤ
  Lpt := ℤ
  thetaOf := fun _ => 0
  H1 := ℤ; H1Grp := inferInstance
  delta := id
  lemma_4_1 := fun θ => ⟨fun ⟨_, hl⟩ => hl.symm, fun h => ⟨0, h.symm⟩⟩
  Zhat := ℤˣ; ZhatGrp := inferInstance
  GammaK := ℤ
  actionAlpha := fun _ φ => φ
  actionBeta := fun _ φ => φ
  zetaZ := fun _ => 1
  lemma_4_2 := fun _ φ => by group
  AlphaPrime := ℤ
  isJGeometric := fun _ => True
  satisfiesSpl := fun _ => True
  prop_4_4 := fun _ => Iff.rfl

@[reducible] def JGeometricSetup.nonvacuous : Nonempty (JGeometricSetup.{0}) := ⟨JGeometricSetup.example⟩

end ABC3.Interface.LocProP
