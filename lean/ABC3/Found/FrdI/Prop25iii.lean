import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.Def24

/-!
# [FrdI] Proposition 2.5, (iii) へ向けた 4 重分解

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.49。

原文 (FrdI p.49):
> Finally, we consider assertion (iii). By applying the factorizations of Definition

原文 (FrdI p.49):
> alences of categories of Definition 1.3, (iii), (d)], we conclude that every morphism

## ★この段の内容

`Proposition 2.5, (iii)`(unit-linear Frobenius 函手 `Ψ : 𝒞 ≃ 𝒞(d)`)の
証明が最初に用意する **4 重分解**

  `φ = δ ≫ γ ≫ β ≫ α`
  (`δ` Frobenius 型、`γ` 等長 pre-step、`β` **base-identity な pre-step 自己射**、
   `α` pull-back 射)

をここで取る。★**原文が「`Definition 1.3, (iv), (a)`; `(v), (c)` と
assertion (i) の全単射(および `Definition 1.3, (iii), (d)` の圏同値)を当てる」と
一言で述べる箇所**である。

★★**要点は `β` を「自己射」に、しかも「base-identity」にすること。**
`Definition 1.3` の分解が与えるのは co-angular pre-step `β₀ : X ⟶ Y` までで、
- **metrically trivial 型**が `Y ≅ X` を与えて**自己射**にし、
- **Aut-ample 型**がずれた底の自己同型を打ち消して **base-identity** にする。

★`δ`・`γ`・`α` は**どれも等長射**である(Frobenius 型は LB-invertible、
pull-back も LB-invertible)。★**これが (a)「`Ψ` は等長射の上で恒等」から
`Ψ` が一意に決まる理由**でもある —— 非等長な部分は `β` に集まっている。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

include P in
/-- ★★★**4 重分解** —— metrically trivial ＋ Aut-ample 型のもとで、
どの射も `δ ≫ γ ≫ β ≫ α` と分解する。

| 因子 | 型 |
|---|---|
| `δ` | Frobenius 型 |
| `γ` | 等長 pre-step |
| `β` | **base-identity な pre-step 自己射**(＝ `𝒪^▷(X)` の元) |
| `α` | pull-back 射 | -/
theorem quadFactor (Fc : FrobenioidCore P)
    (hmt : ∀ X : C, IsMetricallyTrivial P X) (haa : IsOfAutAmpleType P)
    {A B : C} (φ : A ⟶ B) :
    ∃ (X Y : C) (δ : A ⟶ X) (γ : X ⟶ Y) (β : Y ⟶ Y) (α : Y ⟶ B),
      φ = δ ≫ γ ≫ β ≫ α ∧ IsFrobeniusType P δ ∧
        IsIsometric P γ ∧ IsPreStep P γ ∧
        IsBaseIdentity P β ∧ IsPreStep P β ∧ IsPullBack P α := by
  -- 手 1: `Definition 1.3, (iv), (a)`
  obtain ⟨X, Z, δ, p, α₀, hfac, hδ, hp, hα₀⟩ := Fc.arbFactor φ
  -- 手 2: `Definition 1.3, (v), (c)` —— pre-step を「等長 ≫ co-angular」に
  obtain ⟨Y, γ, β₀, hpfac, hγi, hγs, hβ₀c, hβ₀s⟩ := Fc.preStepFactor' p hp
  -- 手 3: metrically trivial 型で `β₀` の終域を `Y` に戻す
  obtain ⟨e⟩ := hmt Y Z β₀ hβ₀c hβ₀s
  -- 手 4: Aut-ample 型で底を `𝟙` にする
  have hβ₁s : IsPreStep P (β₀ ≫ e.hom) :=
    IsPreStep.comp P hβ₀s (isPreStep_of_isIso P e.hom)
  haveI hb₁ : IsIso (P.Base (β₀ ≫ e.hom)) := hβ₁s.2
  obtain ⟨u, hui, hub⟩ := haa Y (P.Base (β₀ ≫ e.hom)) hb₁
  haveI := hui
  obtain ⟨w, hw1, hw2⟩ : ∃ w : Y ⟶ Y, (u : Y ⟶ Y) ≫ w = 𝟙 Y ∧ w ≫ (u : Y ⟶ Y) = 𝟙 Y :=
    ⟨inv (u : Y ⟶ Y), IsIso.hom_inv_id _, IsIso.inv_hom_id _⟩
  haveI hwi : IsIso w := ⟨(u : Y ⟶ Y), hw2, hw1⟩
  refine ⟨X, Y, δ, γ, (β₀ ≫ e.hom) ≫ w, ((u : Y ⟶ Y) ≫ e.inv) ≫ α₀,
    ?_, hδ, hγi, hγs, ?_, ?_, ?_⟩
  · -- 分解が元の射に戻る
    have hkey : ((β₀ ≫ e.hom) ≫ w) ≫ (((u : Y ⟶ Y)) ≫ e.inv) ≫ α₀ = β₀ ≫ α₀ := by
      simp only [Category.assoc]
      rw [← Category.assoc w, hw2, Category.id_comp, e.hom_inv_id_assoc]
    rw [hfac, hpfac, hkey]
    simp
  · -- `β` は base-identity
    show P.Base ((β₀ ≫ e.hom) ≫ w) = P.Base (𝟙 Y)
    rw [P.Base_comp, P.Base_id, ← hub, ← P.Base_comp, hw1, P.Base_id]
  · -- `β` は pre-step
    exact IsPreStep.comp P hβ₁s (isPreStep_of_isIso P w)
  · -- `α` は pull-back
    exact IsPullBack.comp P
      (IsPullBack.comp P (isPullBack_of_isIso P (u : Y ⟶ Y)) (isPullBack_of_isIso P e.inv))
      hα₀

include P in
/-- ★**4 重分解の `δ`・`γ`・`α` はどれも等長射**。

★**これが原文の (a)「`Ψ` は対象と等長射の上で恒等」から
`Ψ` が一意に決まる理由**である —— 非等長な部分は `β` にすべて集まる。 -/
theorem quadFactor_isometric (Fc : FrobenioidCore P) {A X Y B : C}
    {δ : A ⟶ X} {γ : X ⟶ Y} {α : Y ⟶ B}
    (hδ : IsFrobeniusType P δ) (hγ : IsIsometric P γ) (hα : IsPullBack P α) :
    IsIsometric P δ ∧ IsIsometric P γ ∧ IsIsometric P α :=
  ⟨hδ.1.2, hγ, (Fc.pullBackLB α hα).1.2⟩

include P in
/-- ★**`β` は `𝒪^▷(Y)` の元**である。 -/
theorem quadFactor_mem_otri {Y : C} {β : Y ⟶ Y}
    (hβb : IsBaseIdentity P β) (hβs : IsPreStep P β) : β ∈ OTri P Y :=
  ⟨hβb, hβs.1⟩

include P in
/-- ★★**分解から `Div` が読める** —— `Div φ = Φ.map (Base (δ ≫ γ)) (Div β)`。

★`δ`・`γ`・`α` が等長で `γ`・`β`・`α` が linear なので、
**`Div` はすべて `β` から来る**。 -/
theorem quadFactor_div (Fc : FrobenioidCore P) {A X Y B : C}
    {δ : A ⟶ X} {γ : X ⟶ Y} {β : Y ⟶ Y} {α : Y ⟶ B}
    (hδ : IsFrobeniusType P δ) (hγi : IsIsometric P γ) (hγs : IsPreStep P γ)
    (hβb : IsBaseIdentity P β) (hβs : IsPreStep P β) (hα : IsPullBack P α) :
    P.Div (δ ≫ γ ≫ β ≫ α) = Φ.map (P.Base (δ ≫ γ)) (P.Div β) := by
  have hαd : P.Div α = 0 := (Fc.pullBackLB α hα).1.2
  have hαdeg : P.degFr α = 1 := (Fc.pullBackLB α hα).2
  have hβbb : P.Base β = 𝟙 _ := by
    have h : P.Base β = P.Base (𝟙 _) := hβb
    rwa [P.Base_id] at h
  have h1 : P.Div (β ≫ α) = P.Div β := by
    rw [P.Div_comp, hαd, hβbb, MonoidOn.map_id, hαdeg]
    simp
  have h2 : P.Div (γ ≫ β ≫ α) = Φ.map (P.Base γ) (P.Div β) := by
    rw [P.Div_comp, h1, hγi]
    simp
  rw [P.Div_comp, h2, show P.Div δ = 0 from hδ.1.2, P.Base_comp]
  simp

end ABC3.Found.FrdI
