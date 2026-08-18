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


/-! ## ★7. ★★★★★`Theorem 3.4, (ii)` —— 3 主張 -/

include P₁ in
/-- ★`A` から出る co-angular step の存在は、対象の同型で移る(isotropic 型のとき)。 -/
theorem exists_coaStep_of_iso (F₁ : FrobenioidCore P₁) (hiso : ∀ X : C₁, IsIsotropic P₁ X)
    {A A' : C₁} (k : A ≅ A')
    (h : ∃ (B : C₁) (φ : A ⟶ B), IsCoAngular P₁ φ ∧ IsStep P₁ φ) :
    ∃ (B : C₁) (φ : A' ⟶ B), IsCoAngular P₁ φ ∧ IsStep P₁ φ := by
  obtain ⟨B, φ, -, hst⟩ := h
  refine ⟨B, k.inv ≫ φ, isCoAngular_of_isotropic_dom (P := P₁) F₁ (hiso A') _, ?_, ?_⟩
  · exact IsPreStep.comp P₁ (isPreStep_of_isIso P₁ k.inv) hst.1
  · intro hi
    haveI := hi
    refine hst.2 ?_
    have hk : φ = k.hom ≫ (k.inv ≫ φ) := by
      rw [← Category.assoc, Iso.hom_inv_id, Category.id_comp]
    rw [hk]
    infer_instance

include P₁ P₂ in
/-- ★★**pre-step の保存から、出ていく co-angular step の存在が移る**。 -/
theorem exists_coaStep_map (e : C₁ ≌ C₂)
    (F₂ : FrobenioidCore P₂) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f))
    {A : C₁} (h : ∃ (B : C₁) (φ : A ⟶ B), IsCoAngular P₁ φ ∧ IsStep P₁ φ) :
    ∃ (B : C₂) (φ : e.functor.obj A ⟶ B), IsCoAngular P₂ φ ∧ IsStep P₂ φ := by
  obtain ⟨B, φ, -, hst⟩ := h
  refine ⟨e.functor.obj B, e.functor.map φ,
    isCoAngular_of_isotropic_dom (P := P₂) F₂ (hiso₂ _) _, (hps φ).mp hst.1, ?_⟩
  intro hi
  haveI := hi
  exact hst.2 ((Functor.FullyFaithful.ofFullyFaithful e.functor).isIso_of_isIso_map φ)

include P₁ P₂ in
/-- ★★★★**[FrdI] Theorem 3.4, (ii)** —— **圏同値は group-like 対象を保つ**。

原文 (FrdI p.64):
> also Proposition 1.4, (i)], that Ψ preserves group-like objects. This completes the

★`Proposition 1.8, (iii)` が「non-group-like ⟺ 出ていく co-angular step がある」を
与えるので、pre-step の保存からそのまま出る。 -/
theorem isGroupLikeObj_map_iff (e : C₁ ≌ C₂)
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f))
    (hps' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPreStep P₂ f ↔ IsPreStep P₁ (e.inverse.map f))
    (A : C₁) :
    IsGroupLikeObj P₁ A ↔ IsGroupLikeObj P₂ (e.functor.obj A) := by
  letI := coaPreProp_isMultiplicative P₁ G₁.core.coAngularComp
  letI := coaPreProp_isMultiplicative P₂ G₂.core.coAngularComp
  rw [← not_iff_not, prop_1_8_iii_out P₁ (fun X => G₁.coaPreUnderEquiv X) A,
    prop_1_8_iii_out P₂ (fun X => G₂.coaPreUnderEquiv X) (e.functor.obj A)]
  constructor
  · exact exists_coaStep_map P₁ P₂ e F₂ hiso₂ hps
  · intro h
    have h1 := exists_coaStep_map P₂ P₁ e.symm F₁ hiso₁ hps' h
    exact exists_coaStep_of_iso P₁ F₁ hiso₁ (e.unitIso.app A).symm h1

set_option maxHeartbeats 1000000 in
include P₁ P₂ in
/-- ★★★★★**[FrdI] Theorem 3.4, (ii)** —— **3 主張をまとめた形**。

原文 (FrdI p.62):
> (ii) Suppose that C1, C2 are of quasi-isotropic type, and that D1, D2 are of

★★**逸脱の記録(分類 (B))** —— 2 点ある。

1. `Proposition 1.14, (iii)` の `⟸` が `Definition 1.3` から出ない
   (反例 `𝔽_ℕ ⋉ ∏ℤ/n`、`Gap_1_14_iii`)ので、原典に無い 2 条
   `hFrobMono`(Frobenius 型は mono)と `hFrobFS`(同 fiberwise 全射)を
   **両側で**仮定に置いた。
2. 原文は **quasi-isotropic 型**で述べ、`Proposition 1.7, (iv)` と (i) で
   **isotropic 型へ帰着**させる。★我々は**isotropic 型を直接仮定している**
   (帰着の側は未実装)。したがって原文より**狭い**主張である。

★内容は 3 つ: pre-step の保存・co-angular pre-step の保存・group-like 対象の保存。
★isotropic 型では co-angular はすべての射について自動なので、2 番目は 1 番目に帰着する。 -/
theorem thm_3_4_ii (e : C₁ ≌ C₂)
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hFSMFF₁ : IsOfFSMFFType D₁)
    (hFrobMono₁ : ∀ {X Y : C₁} (ε : X ⟶ Y), IsFrobeniusType P₁ ε → Mono ε)
    (hFrobFS₁ : ∀ {X Y : C₁} (ε : X ⟶ Y), IsFrobeniusType P₁ ε → IsFiberwiseSurjective ε)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X) (hFSMFF₂ : IsOfFSMFFType D₂)
    (hFrobMono₂ : ∀ {X Y : C₂} (ε : X ⟶ Y), IsFrobeniusType P₂ ε → Mono ε)
    (hFrobFS₂ : ∀ {X Y : C₂} (ε : X ⟶ Y), IsFrobeniusType P₂ ε → IsFiberwiseSurjective ε) :
    (∀ {A B : C₁} (φ : A ⟶ B), IsPreStep P₁ φ ↔ IsPreStep P₂ (e.functor.map φ)) ∧
    (∀ {A B : C₁} (φ : A ⟶ B),
      (IsCoAngular P₁ φ ∧ IsPreStep P₁ φ) ↔
        (IsCoAngular P₂ (e.functor.map φ) ∧ IsPreStep P₂ (e.functor.map φ))) ∧
    (∀ A : C₁, IsGroupLikeObj P₁ A ↔ IsGroupLikeObj P₂ (e.functor.obj A)) := by
  have hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f) :=
    fun f => isPreStep_map_iff P₁ P₂ e F₁ G₁ hiso₁ hFSMFF₁ hFrobMono₁ hFrobFS₁
      F₂ G₂ hiso₂ hFSMFF₂ hFrobMono₂ hFrobFS₂ f
  have hps' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPreStep P₂ f ↔ IsPreStep P₁ (e.inverse.map f) :=
    fun f => isPreStep_map_iff P₂ P₁ e.symm F₂ G₂ hiso₂ hFSMFF₂ hFrobMono₂ hFrobFS₂
      F₁ G₁ hiso₁ hFSMFF₁ hFrobMono₁ hFrobFS₁ f
  refine ⟨hps, fun φ => ?_, isGroupLikeObj_map_iff P₁ P₂ e F₁ G₁ F₂ G₂ hiso₁ hiso₂ hps hps'⟩
  constructor
  · rintro ⟨-, h⟩
    exact ⟨isCoAngular_of_isotropic_dom (P := P₂) F₂ (hiso₂ _) _, (hps φ).mp h⟩
  · rintro ⟨-, h⟩
    exact ⟨isCoAngular_of_isotropic_dom (P := P₁) F₁ (hiso₁ _) _, (hps φ).mpr h⟩

end PreStep

/-! ## ★8. ★★★★`Theorem 3.4, (iii)` の入口 —— prime-Frobenius 射の保存

原文 (FrdI p.64):
> ces to prove that, for each prime p1 ∈Primes, there exists a prime p2 ∈Primes,

★★原文は「(iii) を示すには、**`Ψ^istr` が p₁-Frobenius 射を p₂-Frobenius 射に写す**
ことだけ示せばよい」と書く。★その入口が本節である。

★★**`Proposition 1.14, (v)`** が
`IsPrimeFrobenius ∧ IsDivIdentity` の**圏論的特徴づけ**を与える。
★その右辺は `IsStep`・`IsIrreducibleMor`・`IsPreStep` だけで書かれており、
**(ii) で得た pre-step の保存がそのまま効く**。 -/

section Step

universe w

variable {C : Type u2} [Category.{v2} C] {D : Type u} [Category.{v} D]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

include P in
theorem isStep_isIso_comp {A A' B : C} (u : A ⟶ A') [IsIso u] {f : A' ⟶ B}
    (h : IsStep P f) : IsStep P (u ≫ f) := by
  refine ⟨IsPreStep.comp P (isPreStep_of_isIso P u) h.1, fun hi => ?_⟩
  haveI := hi
  refine h.2 ?_
  have hf : f = inv u ≫ (u ≫ f) := by
    rw [← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  rw [hf]; infer_instance

include P in
theorem isStep_comp_isIso {A B E : C} {f : A ⟶ B} (h : IsStep P f) (v : B ⟶ E) [IsIso v] :
    IsStep P (f ≫ v) := by
  refine ⟨IsPreStep.comp P h.1 (isPreStep_of_isIso P v), fun hi => ?_⟩
  haveI := hi
  refine h.2 ?_
  have hf : f = (f ≫ v) ≫ inv v := by
    rw [Category.assoc, IsIso.hom_inv_id, Category.comp_id]
  rw [hf]; infer_instance

include P in
theorem not_isPreStep_isIso_comp {A A' B : C} (u : A ⟶ A') [IsIso u] {f : A' ⟶ B}
    (h : ¬ IsPreStep P f) : ¬ IsPreStep P (u ≫ f) := by
  intro hc
  refine h ?_
  have hf : f = inv u ≫ (u ≫ f) := by
    rw [← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  rw [hf]
  exact IsPreStep.comp P (isPreStep_of_isIso P (inv u)) hc

end Step

section PrimeFrob

universe w w4 v5 u5 v6 u6

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} (P₁ : PreFrobenioid C₁ Φ₁)
  {D₂ : Type u5} [Category.{v5} D₂] {C₂ : Type u6} [Category.{v6} C₂]
  {Φ₂ : MonoidOn.{v5, u5, w4} D₂} (P₂ : PreFrobenioid C₂ Φ₂)

include P₁ P₂ in
/-- ★★**step も圏同値で移り、反射する** —— `step = pre-step ＋ 非同型`。 -/
theorem isStep_map_iff (e : C₁ ≌ C₂)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f))
    {A B : C₁} (φ : A ⟶ B) :
    IsStep P₁ φ ↔ IsStep P₂ (e.functor.map φ) := by
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨(hps φ).mp h1, fun hi => ?_⟩
    haveI := hi
    exact h2 ((Functor.FullyFaithful.ofFullyFaithful e.functor).isIso_of_isIso_map φ)
  · rintro ⟨h1, h2⟩
    refine ⟨(hps φ).mpr h1, fun hi => ?_⟩
    haveI := hi
    exact h2 inferInstance

variable (C₁) in
include P₁ in
/-- ★★`Proposition 1.14, (v)` の右辺 —— **純粋に圏論的な条件**
(`IsStep`・`IsIrreducibleMor`・`IsPreStep` だけで書かれている)。

★★これが `IsPrimeFrobenius ∧ IsDivIdentity` と同値であることが `prop_1_14_v`。 -/
def PrimeFrobCond {A : C₁} (φ : A ⟶ A) : Prop :=
  ∀ (B : C₁) (α : A ⟶ B), IsStep P₁ α →
    ∃ (B' : C₁) (ψ : B ⟶ B') (β : B ⟶ B'),
      IsIrreducibleMor ψ ∧ ¬ IsPreStep P₁ ψ ∧ IsStep P₁ β ∧ α ≫ ψ = φ ≫ α ≫ β

include P₁ P₂ in
/-- ★★★★**prime-Frobenius の圏論的条件は圏同値で移る**。

★★`Theorem 3.4, (iii)` は「`Ψ^istr` が p-Frobenius 射を保つ」ことに帰着するので、
`prop_1_14_v` と本補題を組めば、(ii) の pre-step 保存から (iii) の入口が出る。

★手: 余単位で `α₂` を `C₁` へ引き戻し(`ff.preimage`)、
そこで条件を使ってから像を余単位で戻す。★`step` / `irreducible` / `非 pre-step` は
いずれも同型との合成で保たれるので、余単位の付け外しは無害である。 -/
theorem primeFrobCond_map (e : C₁ ≌ C₂)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f))
    {A : C₁} (φ : A ⟶ A) (h : PrimeFrobCond C₁ P₁ φ) :
    PrimeFrobCond C₂ P₂ (e.functor.map φ) := by
  intro B₂ α₂ hα₂
  set ff := Functor.FullyFaithful.ofFullyFaithful e.functor with hff
  set ε := e.counitIso.app B₂ with hε
  set α₀ := ff.preimage (α₂ ≫ ε.inv) with hα₀
  have hmap : e.functor.map α₀ = α₂ ≫ ε.inv := ff.map_preimage _
  have hst₀ : IsStep P₁ α₀ := by
    refine (isStep_map_iff P₁ P₂ e hps α₀).mpr ?_
    rw [hmap]
    exact isStep_comp_isIso P₂ hα₂ ε.inv
  obtain ⟨B'₀, ψ₀, β₀, hirr, hnps, hstβ, heq⟩ := h _ α₀ hst₀
  refine ⟨e.functor.obj B'₀, ε.inv ≫ e.functor.map ψ₀, ε.inv ≫ e.functor.map β₀,
    (isIrreducibleMor_map e.functor hirr).isIso_comp ε.inv,
    not_isPreStep_isIso_comp P₂ ε.inv ((hps ψ₀).not.mp hnps),
    isStep_isIso_comp P₂ ε.inv ((isStep_map_iff P₁ P₂ e hps β₀).mp hstβ), ?_⟩
  have hE := congrArg e.functor.map heq
  rw [e.functor.map_comp, e.functor.map_comp, e.functor.map_comp, hmap] at hE
  refine Eq.trans (Category.assoc _ _ _).symm ?_
  refine Eq.trans hE ?_
  exact congrArg (fun t => e.functor.map φ ≫ t)
    (Category.assoc α₂ ε.inv (e.functor.map β₀))

end PrimeFrob

/-! ### ★共役不変性と両向きの移送 -/

section Conj

universe w7

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w7} D} (P : PreFrobenioid C Φ)

include P in
/-- ★★`PrimeFrobCond` は同型による共役で不変。

★`θ.hom ≫ α` に条件を当て、`θ.inv` を前から掛けて戻すだけ。 -/
theorem primeFrobCond_conj {A A' : C} (θ : A ≅ A') {φ : A ⟶ A}
    (h : PrimeFrobCond C P φ) : PrimeFrobCond C P (θ.inv ≫ φ ≫ θ.hom) := by
  intro B α hα
  obtain ⟨B', ψ, β, hirr, hnps, hstβ, heq⟩ :=
    h B (θ.hom ≫ α) (isStep_isIso_comp P θ.hom hα)
  refine ⟨B', ψ, β, hirr, hnps, hstβ, ?_⟩
  have h1 : θ.inv ≫ ((θ.hom ≫ α) ≫ ψ) = α ≫ ψ := by
    rw [← Category.assoc, ← Category.assoc, Iso.inv_hom_id, Category.id_comp]
  rw [← h1, heq]
  simp only [Category.assoc]

end Conj

section PrimeFrob2

universe w8 v9 u9 v10 u10

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w8} D₁} (P₁ : PreFrobenioid C₁ Φ₁)
  {D₂ : Type u9} [Category.{v9} D₂] {C₂ : Type u10} [Category.{v10} C₂]
  {Φ₂ : MonoidOn.{v9, u9, w8} D₂} (P₂ : PreFrobenioid C₂ Φ₂)

include P₁ P₂ in
/-- ★★★★**prime-Frobenius の圏論的条件は圏同値で移り、かつ反射する**。

★`⟸` は擬逆に同じ補題を当て、単位同型による共役を `primeFrobCond_conj` で剥がす。 -/
theorem primeFrobCond_map_iff (e : C₁ ≌ C₂)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f))
    (hps' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPreStep P₂ f ↔ IsPreStep P₁ (e.inverse.map f))
    {A : C₁} (φ : A ⟶ A) :
    PrimeFrobCond C₁ P₁ φ ↔ PrimeFrobCond C₂ P₂ (e.functor.map φ) := by
  refine ⟨primeFrobCond_map P₁ P₂ e hps φ, fun h => ?_⟩
  have h1 : PrimeFrobCond C₁ P₁ (e.inverse.map (e.functor.map φ)) :=
    primeFrobCond_map P₂ P₁ e.symm hps' _ h
  have hnat : φ ≫ e.unit.app A = e.unit.app A ≫ e.inverse.map (e.functor.map φ) :=
    e.unit.naturality φ
  have h2 := primeFrobCond_conj P₁ (e.unitIso.app A).symm h1
  have heq : (e.unitIso.app A).symm.inv ≫ e.inverse.map (e.functor.map φ)
      ≫ (e.unitIso.app A).symm.hom = φ := by
    show e.unit.app A ≫ e.inverse.map (e.functor.map φ) ≫ (e.unitIso.app A).inv = φ
    refine Eq.trans (Category.assoc _ _ _).symm ?_
    rw [← hnat]
    refine Eq.trans (Category.assoc _ _ _) ?_
    exact Eq.trans (congrArg (fun t => φ ≫ t) ((e.unitIso.app A).hom_inv_id))
      (Category.comp_id φ)
  exact heq ▸ h2


include P₁ P₂ in
/-- ★★★★★**[FrdI] Theorem 3.4, (iii) の入口** ——
**圧同値は prime-Frobenius な Div-identity 自己射を保ち、反射する**。

原文 (FrdI p.64):
> ces to prove that, for each prime p1 ∈Primes, there exists a prime p2 ∈Primes,

★`Proposition 1.14, (v)` を両側で使って圧論的な条件に直し、
`primeFrobCond_map_iff` で運ぶだけ。 -/
theorem primeFrobDivId_map_iff (e : C₁ ≌ C₂)
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hnd₁ : MonoidOn.IsNonDilatingOn Φ₁)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X) (hnd₂ : MonoidOn.IsNonDilatingOn Φ₂)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f))
    (hps' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPreStep P₂ f ↔ IsPreStep P₁ (e.inverse.map f))
    {A : C₁} (hA : ¬ IsGroupLikeObj P₁ A)
    (hA₂ : ¬ IsGroupLikeObj P₂ (e.functor.obj A))
    (φ : A ⟶ A) (hirr : IsIrreducibleMor φ) (hnps : ¬ IsPreStep P₁ φ) :
    (IsPrimeFrobenius P₁ φ ∧ IsDivIdentity P₁ φ) ↔
      (IsPrimeFrobenius P₂ (e.functor.map φ) ∧ IsDivIdentity P₂ (e.functor.map φ)) := by
  rw [prop_1_14_v P₁ F₁ G₁ hiso₁ hnd₁ hA φ hirr hnps,
    prop_1_14_v P₂ F₂ G₂ hiso₂ hnd₂ hA₂ (e.functor.map φ)
      (isIrreducibleMor_map e.functor hirr) ((hps φ).not.mp hnps)]
  exact primeFrobCond_map_iff P₁ P₂ e hps hps' φ

include P₁ P₂ in
/-! ## ★★(F1) の非 group-like の場合の中核

原文 (FrdI p.65):
> exists a base-identity [hence Div-identity] p1-Frobenius endomorphism φ

★原典は base-identity な p₁-Frobenius 自己射を取って、
`Proposition 1.14, (v)` の特徴づけで像も prime-Frobenius になることを出す。 -/

/-- ★★★★**base-identity な prime-Frobenius 自己射は像も prime-Frobenius**。

★`¬ IsPreStep` は**次数が素数だから**自動で出る(1 は素数でない)。 -/
theorem primeFrob_baseId_map (e : C₁ ≌ C₂)
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hnd₁ : MonoidOn.IsNonDilatingOn Φ₁)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X) (hnd₂ : MonoidOn.IsNonDilatingOn Φ₂)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f))
    (hps' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPreStep P₂ f ↔ IsPreStep P₁ (e.inverse.map f))
    {A : C₁} (hA : ¬ IsGroupLikeObj P₁ A)
    (hA₂ : ¬ IsGroupLikeObj P₂ (e.functor.obj A))
    (φ : A ⟶ A) (hirr : IsIrreducibleMor φ)
    (hpf : IsPrimeFrobenius P₁ φ) (hbid : IsBaseIdentity P₁ φ) :
    IsPrimeFrobenius P₂ (e.functor.map φ) := by
  have hnps : ¬ IsPreStep P₁ φ := by
    intro h
    have h1 : P₁.degFr φ = 1 := h.1
    have h2 := hpf.2
    rw [h1] at h2
    exact Nat.not_prime_one (by simpa using h2)
  exact ((primeFrobDivId_map_iff P₁ P₂ e F₁ G₁ hiso₁ hnd₁ F₂ G₂ hiso₂ hnd₂
    hps hps' hA hA₂ φ hirr hnps).mp ⟨hpf, isDivIdentity_of_isBaseIdentity P₁ hbid⟩).1

def primeFrob_baseId_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F1) — base-identity な prime-Frobenius 自己射の像",
    sectionId := "frdi-thm-3-4" }

/-- ★★`hA₂` を `isGroupLikeObj_map_iff` から導いた形。 -/
theorem primeFrobDivId_map_iff' (e : C₁ ≌ C₂)
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hnd₁ : MonoidOn.IsNonDilatingOn Φ₁)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X) (hnd₂ : MonoidOn.IsNonDilatingOn Φ₂)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f))
    (hps' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPreStep P₂ f ↔ IsPreStep P₁ (e.inverse.map f))
    {A : C₁} (hA : ¬ IsGroupLikeObj P₁ A)
    (φ : A ⟶ A) (hirr : IsIrreducibleMor φ) (hnps : ¬ IsPreStep P₁ φ) :
    (IsPrimeFrobenius P₁ φ ∧ IsDivIdentity P₁ φ) ↔
      (IsPrimeFrobenius P₂ (e.functor.map φ) ∧ IsDivIdentity P₂ (e.functor.map φ)) :=
  primeFrobDivId_map_iff P₁ P₂ e F₁ G₁ hiso₁ hnd₁ F₂ G₂ hiso₂ hnd₂ hps hps' hA
    (fun h => hA ((isGroupLikeObj_map_iff P₁ P₂ e F₁ G₁ F₂ G₂ hiso₁ hiso₂ hps hps' A).mpr h))
    φ hirr hnps

end PrimeFrob2



end ABC3.Found.FrdI
