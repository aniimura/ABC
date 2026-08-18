/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34
import ABC3.Found.FrdI.Prop114

/-!
# [FrdI] `Theorem 3.4, (ii)` の入口 —— 圏論的な述語の移送

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.63。

原文 (FrdI p.63):
> phisms] the fact that Ψ preserves pre-steps follows formally from Proposition 1.14,

## ★原文の括弧書きを開く

原文は `Theorem 3.4, (ii)` の証明で
「**任意の圏同値は明らかに FSM 射と irreducible 射を保つ**」と括弧で済ませ、
そのうえで「pre-step を保つことは `Proposition 1.14, (ii)(iii)` から**形式的に**従う」と書く。

★★**「形式的に」の中身は移送 3 本である**:

| 述語 | 宣言 | どこで要るか |
|---|---|---|
| `IsFSMMorphism` | `isFSMMorphism_map_iff` | `Prop 1.14, (ii)` の右辺 第 1 項 |
| `IsIrreducibleMor` | `isIrreducibleMor_map_iff` | `irredNonPreStep` の第 1 項 |
| `BoundedFSMIFactor` | `boundedFSMIFactor_map_equiv` | `Prop 1.14, (iii)` 経由で `irredNonPreStep` の第 2 項 |

★`Thm34.lean` に**保つ**側(`isFSMMorphism_map` / `isIrreducibleMor_map` / `isFSMI_map`)は既にある。
★★本ファイルが足すのは**反射する**側であり、`Prop 1.14` の同値を**両向きに**運ぶのに要る。

## ★型の罠(記録)

`Ψ.asEquivalence.functor` は `Ψ` と **defeq だが構文的に別物**で、
`(Ψ.asEquivalence.functor ⋙ Ψ.asEquivalence.inverse).obj A` と
`Ψ.inv.obj (Ψ.obj A)` が `rw` の照合で食い違う。
★★そこで `BoundedFSMIFactor` の移送だけは **`Equivalence` を直接受ける形**で書き、
`e.functor` / `e.inverse` / `e.unit` に統一した。
★最後の結合律も `rw` ではなく `congrArg` ＋ `Eq.trans` で組んでいる(同じ理由)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u v2 u2

variable {C : Type u} [Category.{v} C] {D : Type u2} [Category.{v2} D]

/-! ## ★1. FSM 射・irreducible 射の反射 -/

theorem isFSMMorphism_isIso_comp {A A' B : C} (e : A ⟶ A') [IsIso e] {f : A' ⟶ B}
    (h : IsFSMMorphism f) : IsFSMMorphism (e ≫ f) := by
  haveI := h.2
  exact ⟨(isFiberwiseSurjective_of_isIso e).comp h.1, mono_comp e f⟩

theorem isFSMMorphism_comp_isIso {A B E : C} {f : A ⟶ B} (h : IsFSMMorphism f)
    (g : B ⟶ E) [IsIso g] : IsFSMMorphism (f ≫ g) := by
  haveI := h.2
  exact ⟨h.1.comp (isFiberwiseSurjective_of_isIso g), mono_comp f g⟩

variable (Ψ : C ⥤ D) [Ψ.IsEquivalence]

/-- ★★**圏同値は FSM 射を反射する**。

★`⟸` は擬逆 `Ψ⁻¹` に `isFSMMorphism_map` を当て、単位同型で挟んだ形を剥がす
(`isFSMI_map_iff` と同じ定型)。 -/
theorem isFSMMorphism_map_iff {B A : C} (β : B ⟶ A) :
    IsFSMMorphism β ↔ IsFSMMorphism (Ψ.map β) := by
  refine ⟨isFSMMorphism_map Ψ, fun h => ?_⟩
  have h' : IsFSMMorphism (Ψ.inv.map (Ψ.map β)) := isFSMMorphism_map Ψ.inv h
  set eA := Ψ.asEquivalence.unitIso.app A with hEA
  set eB := Ψ.asEquivalence.unitIso.app B with hEB
  have hnat : β ≫ eA.hom = eB.hom ≫ Ψ.inv.map (Ψ.map β) :=
    (Ψ.asEquivalence.unitIso).hom.naturality β
  have hcomp : IsFSMMorphism (β ≫ eA.hom) := by
    rw [hnat]; exact isFSMMorphism_isIso_comp eB.hom h'
  have hfin := isFSMMorphism_comp_isIso hcomp eA.inv
  simpa using hfin

/-- ★★**圏同値は irreducible 射を反射する**。 -/
theorem isIrreducibleMor_map_iff {B A : C} (β : B ⟶ A) :
    IsIrreducibleMor β ↔ IsIrreducibleMor (Ψ.map β) := by
  refine ⟨isIrreducibleMor_map Ψ, fun h => ?_⟩
  have h' : IsIrreducibleMor (Ψ.inv.map (Ψ.map β)) := isIrreducibleMor_map Ψ.inv h
  set eA := Ψ.asEquivalence.unitIso.app A with hEA
  set eB := Ψ.asEquivalence.unitIso.app B with hEB
  have hnat : β ≫ eA.hom = eB.hom ≫ Ψ.inv.map (Ψ.map β) :=
    (Ψ.asEquivalence.unitIso).hom.naturality β
  have hcomp : IsIrreducibleMor (β ≫ eA.hom) := by
    rw [hnat]; exact h'.isIso_comp eB.hom
  have hfin := hcomp.comp_isIso eA.inv
  simpa using hfin

/-! ## ★2. FSMI 鎖と `BoundedFSMIFactor` の移送 -/

/-- ★★**圏同値は FSMI 鎖を保つ**(長さを変えない)。 -/
theorem isFSMIChain_map : ∀ {n : ℕ} {A B : C} {f : A ⟶ B},
    IsFSMIChain n f → IsFSMIChain n (Ψ.map f) := by
  intro n A B f h
  induction h with
  | nil => rw [CategoryTheory.Functor.map_id]; exact IsFSMIChain.nil
  | cons hφ _ ih =>
    rw [CategoryTheory.Functor.map_comp]; exact IsFSMIChain.cons (isFSMI_map Ψ hφ) ih

/-- ★★★**圏同値は `BoundedFSMIFactor` を保つ**。

★★これが `Theorem 3.4, (ii)` で `Proposition 1.14, (iii)` を運ぶ経路である。
★`D` 側の鎖を `e.inverse` で `C` へ戻し、単位同型で挟んで長さを保つ。
★`n = 0` は上界が自明なので分けている。 -/
theorem boundedFSMIFactor_map_equiv (e : C ≌ D) {A B : C} (φ : A ⟶ B)
    (h : BoundedFSMIFactor φ) : BoundedFSMIFactor (e.functor.map φ) := by
  obtain ⟨N, hN⟩ := h
  refine ⟨N, fun n E' ψ' χ' hψ' hchain' hfac' => ?_⟩
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · exact Nat.zero_le N
  subst hfac'
  have hnat : φ ≫ e.unit.app B = e.unit.app A ≫ e.inverse.map (e.functor.map φ) :=
    e.unit.naturality φ
  have hchain₀ : IsFSMIChain n (e.unit.app A ≫ e.inverse.map (e.functor.map φ ≫ ψ')) :=
    (isFSMIChain_map e.inverse hchain').isIso_comp (e.unit.app A) hn
  have hψ₀ : IsFSMI (e.unit.app B ≫ e.inverse.map ψ') :=
    isFSMI_isIso_comp (e.unit.app B) (isFSMI_map e.inverse hψ')
  refine hN n _ (e.unit.app B ≫ e.inverse.map ψ') _ hψ₀ hchain₀ ?_
  rw [e.inverse.map_comp]
  refine Eq.trans (Category.assoc _ _ _).symm ?_
  exact Eq.trans (congrArg (· ≫ e.inverse.map ψ') hnat.symm) (Category.assoc _ _ _)

end ABC3.Found.FrdI
