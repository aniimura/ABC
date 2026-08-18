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


/-! ## ★3. `BoundedFSMIFactor` の同型不変性と反射 -/

/-- ★同型を前に付けても `BoundedFSMIFactor` は変わらない。 -/
theorem boundedFSMIFactor_isIso_comp {A A' B : C} (u : A ⟶ A') [IsIso u] {φ : A' ⟶ B}
    (h : BoundedFSMIFactor φ) : BoundedFSMIFactor (u ≫ φ) := by
  obtain ⟨N, hN⟩ := h
  refine ⟨N, fun n E ψ χ hψ hchain hfac => ?_⟩
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · exact Nat.zero_le N
  refine hN n E ψ (inv u ≫ χ) hψ (hchain.isIso_comp (inv u) hn) ?_
  rw [hfac, Category.assoc, ← Category.assoc (inv u), IsIso.inv_hom_id, Category.id_comp]

/-- ★同型を後ろに付けても `BoundedFSMIFactor` は変わらない。 -/
theorem boundedFSMIFactor_comp_isIso {A B E : C} {φ : A ⟶ B} (h : BoundedFSMIFactor φ)
    (v : B ⟶ E) [IsIso v] : BoundedFSMIFactor (φ ≫ v) := by
  obtain ⟨N, hN⟩ := h
  refine ⟨N, fun n E' ψ χ hψ hchain hfac => ?_⟩
  refine hN n E' (v ≫ ψ) χ (isFSMI_isIso_comp v hψ) hchain ?_
  rw [hfac, Category.assoc]

/-- ★同型を後ろから剥がしても変わらない。 -/
theorem boundedFSMIFactor_of_comp_iso {A B E : C} (φ : A ⟶ B) (w : B ≅ E)
    (h : BoundedFSMIFactor (φ ≫ w.hom)) : BoundedFSMIFactor φ := by
  have h2 := boundedFSMIFactor_comp_isIso h w.inv
  simpa using h2

/-- ★★★**圏同値は `BoundedFSMIFactor` を保ち、かつ反射する**。 -/
theorem boundedFSMIFactor_map_iff (e : C ≌ D) {A B : C} (φ : A ⟶ B) :
    BoundedFSMIFactor φ ↔ BoundedFSMIFactor (e.functor.map φ) := by
  refine ⟨boundedFSMIFactor_map_equiv e φ, fun h => ?_⟩
  have h' : BoundedFSMIFactor (e.inverse.map (e.functor.map φ)) :=
    boundedFSMIFactor_map_equiv e.symm _ h
  have hnat : φ ≫ e.unit.app B = e.unit.app A ≫ e.inverse.map (e.functor.map φ) :=
    e.unit.naturality φ
  have h2 : BoundedFSMIFactor (e.unit.app A ≫ e.inverse.map (e.functor.map φ)) :=
    boundedFSMIFactor_isIso_comp (e.unit.app A) h'
  rw [← hnat] at h2
  exact boundedFSMIFactor_of_comp_iso φ (e.unitIso.app B) h2

/-! ## ★4. ★★★純粋に圏論的な射のクラス `irredBounded`

★★ここが `Theorem 3.4, (ii)` の要である。`Proposition 1.14, (iii)` により
`irredNonPreStep P`(**`Definition 1.3` を見る**)は
`irredBounded`(**圏だけを見る**)と一致する。
★一致が言えれば、pre-step の保存は圏同値の一般論に落ちる。 -/

variable (C) in
/-- ★★**irreducible かつ `BoundedFSMIFactor`** —— 純粋に圏論的。 -/
def irredBounded : MorphismProperty C := fun _ _ f => IsIrreducibleMor f ∧ BoundedFSMIFactor f

theorem irredBounded_isIso_comp {A A' B : C} (u : A ⟶ A') [IsIso u] {f : A' ⟶ B}
    (h : irredBounded C f) : irredBounded C (u ≫ f) :=
  ⟨h.1.isIso_comp u, boundedFSMIFactor_isIso_comp u h.2⟩

theorem irredBounded_comp_isIso {A B E : C} {f : A ⟶ B} (h : irredBounded C f)
    (v : B ⟶ E) [IsIso v] : irredBounded C (f ≫ v) :=
  ⟨h.1.comp_isIso v, boundedFSMIFactor_comp_isIso h.2 v⟩

theorem irredBounded_of_isIso_comp {A A' B : C} (u : A ⟶ A') [IsIso u] {f : A' ⟶ B}
    (h : irredBounded C (u ≫ f)) : irredBounded C f := by
  have h2 := irredBounded_isIso_comp (inv u) h
  simpa using h2

theorem irredBounded_of_comp_isIso {A B E : C} {f : A ⟶ B} (v : B ⟶ E) [IsIso v]
    (h : irredBounded C (f ≫ v)) : irredBounded C f := by
  have h2 := irredBounded_comp_isIso h (inv v)
  simpa using h2

/-- ★★★圏同値は `irredBounded` を**保ち、かつ反射する**。 -/
theorem irredBounded_map_iff (e : C ≌ D) {A B : C} (f : A ⟶ B) :
    irredBounded C f ↔ irredBounded D (e.functor.map f) :=
  and_congr (isIrreducibleMor_map_iff e.functor f) (boundedFSMIFactor_map_iff e f)

/-! ## ★5. 3 段の分解の反射と `mid-adjoint` の移送 -/

/-- ★★**3 段の分解も反射する** —— `Ψ.map φ = γ ≫ β ≫ α` を `C` へ戻す。 -/
theorem exists_factor3_of_map_factor {A B : C} (φ : A ⟶ B) {X Y : D}
    (γ : Ψ.obj A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ Ψ.obj B) (h : Ψ.map φ = γ ≫ β ≫ α) :
    ∃ (X₀ Y₀ : C) (γ₀ : A ⟶ X₀) (β₀ : X₀ ⟶ Y₀) (α₀ : Y₀ ⟶ B)
      (ex : Ψ.obj X₀ ≅ X) (ey : Ψ.obj Y₀ ≅ Y),
      φ = γ₀ ≫ β₀ ≫ α₀ ∧ β = ex.inv ≫ Ψ.map β₀ ≫ ey.hom := by
  obtain ⟨X₀, γ₀, δ₀, ex, hfac1, hγ, hδ⟩ :=
    exists_factor_of_map_factor Ψ φ γ (β ≫ α) h
  have hδ' : Ψ.map δ₀ = (ex.hom ≫ β) ≫ α := by
    rw [Category.assoc]
    exact ((Iso.eq_inv_comp ex).mp hδ).symm
  obtain ⟨Y₀, β₀, α₀, ey, hfac2, hβ, hα⟩ :=
    exists_factor_of_map_factor Ψ δ₀ (ex.hom ≫ β) α hδ'
  refine ⟨X₀, Y₀, γ₀, β₀, α₀, ex, ey, by rw [hfac1, hfac2], ?_⟩
  exact (Iso.eq_inv_comp ex).mpr hβ

/-- ★同じ射に対する 2 つのクラスが一致すれば `mid-adjoint` も一致する。 -/
theorem isMidAdjoint_congr {S S' : MorphismProperty C}
    (h : ∀ {X Y : C} (β : X ⟶ Y), S β ↔ S' β) {A B : C} (φ : A ⟶ B) :
    IsMidAdjoint S φ ↔ IsMidAdjoint S' φ :=
  ⟨fun hm X Y γ β α hfac hs => hm X Y γ β α hfac ((h β).mpr hs),
   fun hm X Y γ β α hfac hs => hm X Y γ β α hfac ((h β).mp hs)⟩

/-- ★★★**`mid-adjoint to irredBounded` は圏同値で移り、かつ反射する**。 -/
theorem isMidAdjoint_irredBounded_map_iff (e : C ≌ D) {A B : C} (φ : A ⟶ B) :
    IsMidAdjoint (irredBounded C) φ ↔ IsMidAdjoint (irredBounded D) (e.functor.map φ) := by
  constructor
  · intro h X' Y' γ' β' α' hfac hS
    obtain ⟨X₀, Y₀, γ₀, β₀, α₀, ex, ey, hfac₀, hβ⟩ :=
      exists_factor3_of_map_factor e.functor φ γ' β' α' hfac
    have hS₀ : irredBounded D (e.functor.map β₀) := by
      refine irredBounded_of_comp_isIso ey.hom (irredBounded_of_isIso_comp ex.inv ?_)
      rw [← Category.assoc, ← hβ] at *
      exact hβ ▸ hS
    haveI : IsIso β₀ := h X₀ Y₀ γ₀ β₀ α₀ hfac₀ ((irredBounded_map_iff e β₀).mpr hS₀)
    rw [hβ]
    infer_instance
  · intro h X Y γ β α hfac hS
    have hfac' : e.functor.map φ
        = e.functor.map γ ≫ e.functor.map β ≫ e.functor.map α := by
      rw [hfac, e.functor.map_comp, e.functor.map_comp]
    haveI := h _ _ (e.functor.map γ) (e.functor.map β) (e.functor.map α) hfac'
      ((irredBounded_map_iff e β).mp hS)
    exact (Functor.FullyFaithful.ofFullyFaithful e.functor).isIso_of_isIso_map β

/-! ## ★6. ★★★★★`Theorem 3.4, (ii)` の主部 —— pre-step の保存と反射 -/

section PreStep

universe w w3 v3 u3 v4 u4

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} (P₁ : PreFrobenioid C₁ Φ₁)
  {D₂ : Type u3} [Category.{v3} D₂] {C₂ : Type u4} [Category.{v4} C₂]
  {Φ₂ : MonoidOn.{v3, u3, w3} D₂} (P₂ : PreFrobenioid C₂ Φ₂)

include P₁ in
/-- ★★★**`irredNonPreStep` は純粋に圏論的な `irredBounded` と一致する**
(`Proposition 1.14, (iii)`)。

★★これが `Theorem 3.4, (ii)` の要である —— 一致が言えれば
`Definition 1.3` を見ずに圏同値で運べる。 -/
theorem irredNonPreStep_iff_irredBounded (F : FrobenioidCore P₁) (G : Frobenioid P₁)
    (hiso : ∀ X : C₁, IsIsotropic P₁ X) (hFSMFF : IsOfFSMFFType D₁)
    (hFrobMono : ∀ {X Y : C₁} (ε : X ⟶ Y), IsFrobeniusType P₁ ε → Mono ε)
    (hFrobFS : ∀ {X Y : C₁} (ε : X ⟶ Y), IsFrobeniusType P₁ ε → IsFiberwiseSurjective ε)
    {A B : C₁} (β : A ⟶ B) : irredNonPreStep P₁ β ↔ irredBounded C₁ β := by
  constructor
  · rintro ⟨hirr, hnps⟩
    exact ⟨hirr, (prop_1_14_iii P₁ F G hiso hFSMFF hFrobMono hFrobFS β hirr).mpr hnps⟩
  · rintro ⟨hirr, hb⟩
    exact ⟨hirr, (prop_1_14_iii P₁ F G hiso hFSMFF hFrobMono hFrobFS β hirr).mp hb⟩

include P₁ P₂ in
/-- ★★★★★**[FrdI] Theorem 3.4, (ii) の主部** —— **圏同値は pre-step を保ち、反射する**。

原文 (FrdI p.63):
> phisms] the fact that Ψ preserves pre-steps follows formally from Proposition 1.14,

★★**逸脱の記録(分類 (B))**: `Proposition 1.14, (iii)` の `⟸` が
`Definition 1.3` から出ない(反例 `𝔽_ℕ ⋉ ∏ℤ/n`、`Gap_1_14_iii`)ため、
原典に無い 2 条 `hFrobMono` / `hFrobFS` を**両側で**仮定に置いている。

★筋は 3 段:
1. `Prop 1.14, (ii)` —— `pre-step ⟺ FSM ＋ mid-adjoint to irredNonPreStep`
2. `Prop 1.14, (iii)` —— `irredNonPreStep = irredBounded`(**純粋に圏論的**)
3. あとは 2 つの述語を圏同値で運ぶだけ(`isFSMMorphism_map_iff` /
   `isMidAdjoint_irredBounded_map_iff`)

★★原文が「follows formally from Proposition 1.14, (ii), (iii)」と 1 行で書く所である。 -/
theorem isPreStep_map_iff (e : C₁ ≌ C₂)
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hFSMFF₁ : IsOfFSMFFType D₁)
    (hFrobMono₁ : ∀ {X Y : C₁} (ε : X ⟶ Y), IsFrobeniusType P₁ ε → Mono ε)
    (hFrobFS₁ : ∀ {X Y : C₁} (ε : X ⟶ Y), IsFrobeniusType P₁ ε → IsFiberwiseSurjective ε)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X) (hFSMFF₂ : IsOfFSMFFType D₂)
    (hFrobMono₂ : ∀ {X Y : C₂} (ε : X ⟶ Y), IsFrobeniusType P₂ ε → Mono ε)
    (hFrobFS₂ : ∀ {X Y : C₂} (ε : X ⟶ Y), IsFrobeniusType P₂ ε → IsFiberwiseSurjective ε)
    {A B : C₁} (φ : A ⟶ B) :
    IsPreStep P₁ φ ↔ IsPreStep P₂ (e.functor.map φ) := by
  rw [prop_1_14_ii P₁ F₁ G₁ hiso₁ hFSMFF₁ φ,
    prop_1_14_ii P₂ F₂ G₂ hiso₂ hFSMFF₂ (e.functor.map φ),
    isMidAdjoint_congr
      (irredNonPreStep_iff_irredBounded P₁ F₁ G₁ hiso₁ hFSMFF₁ hFrobMono₁ hFrobFS₁) φ,
    isMidAdjoint_congr
      (irredNonPreStep_iff_irredBounded P₂ F₂ G₂ hiso₂ hFSMFF₂ hFrobMono₂ hFrobFS₂)
      (e.functor.map φ)]
  exact and_congr (isFSMMorphism_map_iff e.functor φ) (isMidAdjoint_irredBounded_map_iff e φ)

end PreStep

end ABC3.Found.FrdI
