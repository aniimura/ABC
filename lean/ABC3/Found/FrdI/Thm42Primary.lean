/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42Div
import ABC3.Found.FrdI.Prop41Cat

/-!
# [FrdI] Theorem 4.2, (i) —— **primary step** の保存

原文 (FrdI p.78):
> B1 → C1 are mapped to primary steps B2 → A2, A2 → C2 with primary composite

## ★★★`Proposition 4.1, (i)` が**圏論的特徴づけ**を与える

`prop_4_1_i`(`Prop41Cat.lean:350`)は、`A` が Div-Frobenius-trivial
(各次数の Div-identity Frobenius 型自己射 `α n` を持つ)で
`Φ(A_𝒟)` が perfect なとき、co-angular step `φ : B ⟶ A` について

```
IsPrimaryElt (Div φ) ↔
  ∀ B' φB φA, IsStep φB → IsStep φA → φ = φB ≫ φA →
    ∃ Y n β' ζ ζ', IsFrobeniusType β' ∧ IsPreStep ζ ∧ IsPreStep ζ' ∧
      ζ = ζ' ≫ φ ∧ β' ≫ ζ = φA ≫ α n
```

★★**右辺に `Div` も `Φ` も現れない** —— `IsStep`・`IsPreStep`・`IsFrobeniusType`
と `α` の族だけである。★`α` の族は `divFrobeniusTrivial_map` が `Ψ` の側にも与える。

★したがって `Proposition 1.14, (v)` のときと**同じ手筋**で右辺が移る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 v3 u3 w3 u4 v4

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u3} [Category.{v3} D₂] {C₂ : Type u4} [Category.{v4} C₂]
  {Φ₂ : MonoidOn.{v3, u3, w3} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★`Proposition 4.1, (i)` の右辺 —— 純粋に圏論的な条件。 -/
def Prop41iRhs (P : PreFrobenioid C₁ Φ₁) {A B : C₁} (φ : B ⟶ A)
    (α : ℕ+ → (A ⟶ A)) : Prop :=
  ∀ (B' : C₁) (φB : B ⟶ B') (φA : B' ⟶ A), IsStep P φB → IsStep P φA → φ = φB ≫ φA →
    ∃ (Y : C₁) (n : ℕ+) (β' : B' ⟶ Y) (ζ : Y ⟶ A) (ζ' : Y ⟶ B),
      IsFrobeniusType P β' ∧ IsPreStep P ζ ∧ IsPreStep P ζ' ∧
      ζ = ζ' ≫ φ ∧ β' ≫ ζ = φA ≫ α n

/-- ★★同じものを第 2 の圏で。 -/
def Prop41iRhs' (P : PreFrobenioid C₂ Φ₂) {A B : C₂} (φ : B ⟶ A)
    (α : ℕ+ → (A ⟶ A)) : Prop :=
  ∀ (B' : C₂) (φB : B ⟶ B') (φA : B' ⟶ A), IsStep P φB → IsStep P φA → φ = φB ≫ φA →
    ∃ (Y : C₂) (n : ℕ+) (β' : B' ⟶ Y) (ζ : Y ⟶ A) (ζ' : Y ⟶ B),
      IsFrobeniusType P β' ∧ IsPreStep P ζ ∧ IsPreStep P ζ' ∧
      ζ = ζ' ≫ φ ∧ β' ≫ ζ = φA ≫ α n

/-- ★★★★★**右辺は `Ψ` で移る**。

★手筋は `prop114vRhs_map` と同じ ——
本質的全射性で `B'` を像の中に取り直し、充満性で `φB`・`φA` を引き戻し、
忠実性で分解式を `𝒞₁` へ移す。 -/
theorem prop41iRhs_map (F₂ : FrobenioidCore P₂) (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f →
      IsFrobeniusType P₂ (Ψ.functor.map f))
    {A B : C₁} (φ : B ⟶ A) (α : ℕ+ → (A ⟶ A)) (h : Prop41iRhs P₁ φ α) :
    Prop41iRhs' P₂ (Ψ.functor.map φ) (fun n => Ψ.functor.map (α n)) := by
  intro B' φB φA hstB hstA hfac
  -- ★段 1: `B'` を像の中に取り直す
  obtain ⟨B'₀, ⟨ε⟩⟩ := Functor.EssSurj.mem_essImage (F := Ψ.functor) B'
  -- ★段 2: 充満性で引き戻す
  obtain ⟨φB₀, hφB₀⟩ := Ψ.functor.map_surjective (φB ≫ ε.inv)
  obtain ⟨φA₀, hφA₀⟩ := Ψ.functor.map_surjective (ε.hom ≫ φA)
  have hpsB : IsPreStep P₂ (φB ≫ ε.inv) :=
    IsPreStep.comp P₂ hstB.1 (isPreStep_of_isIso P₂ ε.inv)
  have hpsA : IsPreStep P₂ (ε.hom ≫ φA) :=
    IsPreStep.comp P₂ (isPreStep_of_isIso P₂ ε.hom) hstA.1
  have hstB₀ : IsStep P₁ φB₀ := by
    refine ⟨hPS' φB₀ (by rw [hφB₀]; exact hpsB), fun hiso => hstB.2 ?_⟩
    haveI := hiso
    haveI : IsIso (Ψ.functor.map φB₀) := Ψ.functor.map_isIso φB₀
    have hb : φB = Ψ.functor.map φB₀ ≫ ε.hom := by
      rw [hφB₀, Category.assoc, ε.inv_hom_id, Category.comp_id]
    rw [hb]; infer_instance
  have hstA₀ : IsStep P₁ φA₀ := by
    refine ⟨hPS' φA₀ (by rw [hφA₀]; exact hpsA), fun hiso => hstA.2 ?_⟩
    haveI := hiso
    haveI : IsIso (Ψ.functor.map φA₀) := Ψ.functor.map_isIso φA₀
    have ha : φA = ε.inv ≫ Ψ.functor.map φA₀ := by
      rw [hφA₀, ← Category.assoc, ε.inv_hom_id, Category.id_comp]
    rw [ha]; infer_instance
  -- ★段 3: 忠実性で分解式を `𝒞₁` へ
  have hfac₀ : φ = φB₀ ≫ φA₀ := by
    refine Ψ.functor.map_injective ?_
    rw [CategoryTheory.Functor.map_comp, hφB₀, hφA₀, hfac, Category.assoc,
      ← Category.assoc ε.inv, ε.inv_hom_id, Category.id_comp]
  -- ★段 4: `𝒞₁` で右辺を使う
  obtain ⟨Y₀, n, β'₀, ζ₀, ζ'₀, hβF, hζ, hζ', hζeq, hsq⟩ := h B'₀ φB₀ φA₀ hstB₀ hstA₀ hfac₀
  refine ⟨Ψ.functor.obj Y₀, n, ε.inv ≫ Ψ.functor.map β'₀, Ψ.functor.map ζ₀,
    Ψ.functor.map ζ'₀, ?_, hPS ζ₀ hζ, hPS ζ'₀ hζ', ?_, ?_⟩
  · exact IsFrobeniusType.comp P₂ F₂ (isFrobeniusType_of_isIso P₂ ε.inv) (hFT β'₀ hβF)
  · rw [← CategoryTheory.Functor.map_comp, ← hζeq]
  · calc (ε.inv ≫ Ψ.functor.map β'₀) ≫ Ψ.functor.map ζ₀
        = ε.inv ≫ Ψ.functor.map (β'₀ ≫ ζ₀) := by
          rw [CategoryTheory.Functor.map_comp, Category.assoc]
      _ = ε.inv ≫ Ψ.functor.map (φA₀ ≫ α n) := by rw [hsq]
      _ = ε.inv ≫ Ψ.functor.map φA₀ ≫ Ψ.functor.map (α n) := by
          rw [CategoryTheory.Functor.map_comp]
      _ = φA ≫ Ψ.functor.map (α n) := by
          rw [hφA₀, ← Category.assoc, ← Category.assoc, ε.inv_hom_id, Category.id_comp]

def prop41iRhs_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 78,
    item := "Theorem 4.2, (i) — Proposition 4.1, (i) の右辺は Ψ で移る",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★★`primary step` の保存 -/

/-- ★★★★★★**[FrdI] Theorem 4.2, (i)** —— `Ψ` は
Div-Frobenius-trivial な対象への co-angular step の **primary 性**を保つ。

★両側で `Proposition 4.1, (i)` を当て、間を `prop41iRhs_map` で繋ぐだけ。 -/
theorem isPrimaryElt_div_map (G₁ : Frobenioid P₁) (G₂ : Frobenioid P₂)
    (F₂ : FrobenioidCore P₂) (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f →
      IsFrobeniusType P₂ (Ψ.functor.map f))
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    {A B : C₁} (α : ℕ+ → (A ⟶ A))
    (hαd : ∀ n, P₁.degFr (α n) = n) (hαD : ∀ n, IsDivIdentity P₁ (α n))
    (hαF : ∀ n, IsFrobeniusType P₁ (α n))
    (hαd₂ : ∀ n, P₂.degFr (Ψ.functor.map (α n)) = n)
    (hαD₂ : ∀ n, IsDivIdentity P₂ (Ψ.functor.map (α n)))
    (hperf₁ : IsPerfectMonoid (Φ₁.val (P₁.toElem.obj A).base))
    (hperf₂ : IsPerfectMonoid (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base))
    (φ : B ⟶ A) (hcφ : IsCoAngular P₁ φ) (hst : IsStep P₁ φ)
    (hcφ₂ : IsCoAngular P₂ (Ψ.functor.map φ)) (hst₂ : IsStep P₂ (Ψ.functor.map φ))
    (h : IsPrimaryElt (P₁.Div φ)) :
    IsPrimaryElt (P₂.Div (Ψ.functor.map φ)) :=
  (prop_4_1_i G₂ hiso₂ (fun n => Ψ.functor.map (α n)) hαd₂ hαD₂
      (fun n => hFT (α n) (hαF n)) hperf₂ (Ψ.functor.map φ) hcφ₂ hst₂).mpr
    (prop41iRhs_map F₂ Ψ hPS hPS' hFT φ α
      ((prop_4_1_i G₁ hiso₁ α hαd hαD hαF hperf₁ φ hcφ hst).mp h))

def isPrimaryElt_div_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 78,
    item := "Theorem 4.2, (i) — Ψ は primary step を保つ",
    sectionId := "frdi-thm-4-2" }

/-- ★★★★★★**[FrdI] Theorem 4.2, (i)** —— `primary pre-step` の保存(語彙どおりの形)。 -/
theorem isPrimaryPreStep_map (G₁ : Frobenioid P₁) (G₂ : Frobenioid P₂)
    (F₂ : FrobenioidCore P₂) (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f →
      IsFrobeniusType P₂ (Ψ.functor.map f))
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    {A B : C₁} (α : ℕ+ → (A ⟶ A))
    (hαd : ∀ n, P₁.degFr (α n) = n) (hαD : ∀ n, IsDivIdentity P₁ (α n))
    (hαF : ∀ n, IsFrobeniusType P₁ (α n))
    (hαd₂ : ∀ n, P₂.degFr (Ψ.functor.map (α n)) = n)
    (hαD₂ : ∀ n, IsDivIdentity P₂ (Ψ.functor.map (α n)))
    (hperf₁ : IsPerfectMonoid (Φ₁.val (P₁.toElem.obj A).base))
    (hperf₂ : IsPerfectMonoid (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base))
    (φ : B ⟶ A) (hcφ : IsCoAngular P₁ φ) (hst : IsStep P₁ φ)
    (hcφ₂ : IsCoAngular P₂ (Ψ.functor.map φ)) (hst₂ : IsStep P₂ (Ψ.functor.map φ))
    (h : IsPrimaryPreStep P₁ φ) :
    IsPrimaryPreStep P₂ (Ψ.functor.map φ) :=
  ⟨hPS φ h.1,
   isPrimaryElt_div_map G₁ G₂ F₂ Ψ hPS hPS' hFT hiso₁ hiso₂ α hαd hαD hαF hαd₂ hαD₂
     hperf₁ hperf₂ φ hcφ hst hcφ₂ hst₂ h.2⟩


/-! ## ★★★★★逆向きの移送(反射) -/

/-- ★★★★★**右辺は `Ψ` で戻る** —— `prop41iRhs_map` の逆向き。

★手筋は同じ ——`𝒞₂` の側で得た `Y`・`β'`・`ζ`・`ζ'` を
本質的全射性と充満性で `𝒞₁` へ引き戻し、忠実性で等式を移す。 -/
theorem prop41iRhs_map_rev (F₂ : FrobenioidCore P₂) (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    (hFT' : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₂ (Ψ.functor.map f) →
      IsFrobeniusType P₁ f)
    {A B : C₁} (φ : B ⟶ A) (α : ℕ+ → (A ⟶ A))
    (h : Prop41iRhs' P₂ (Ψ.functor.map φ) (fun n => Ψ.functor.map (α n))) :
    Prop41iRhs P₁ φ α := by
  intro B' φB φA hstB hstA hfac
  have hstB₂ : IsStep P₂ (Ψ.functor.map φB) := by
    refine ⟨hPS φB hstB.1, fun hiso => hstB.2 ?_⟩
    haveI := hiso
    exact isIso_of_reflects_iso φB Ψ.functor
  have hstA₂ : IsStep P₂ (Ψ.functor.map φA) := by
    refine ⟨hPS φA hstA.1, fun hiso => hstA.2 ?_⟩
    haveI := hiso
    exact isIso_of_reflects_iso φA Ψ.functor
  obtain ⟨Y₂, n, β'₂, ζ₂, ζ'₂, hβF, hζ, hζ', hζeq, hsq⟩ :=
    h (Ψ.functor.obj B') (Ψ.functor.map φB) (Ψ.functor.map φA) hstB₂ hstA₂
      (by rw [hfac, Ψ.functor.map_comp])
  obtain ⟨Y₀, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage (F := Ψ.functor) Y₂
  obtain ⟨β'₀, hβ'₀⟩ := Ψ.functor.map_surjective (β'₂ ≫ e.inv)
  obtain ⟨ζ₀, hζ₀⟩ := Ψ.functor.map_surjective (e.hom ≫ ζ₂)
  obtain ⟨ζ'₀, hζ'₀⟩ := Ψ.functor.map_surjective (e.hom ≫ ζ'₂)
  refine ⟨Y₀, n, β'₀, ζ₀, ζ'₀, ?_, ?_, ?_, ?_, ?_⟩
  · refine hFT' β'₀ ?_
    rw [hβ'₀]
    exact IsFrobeniusType.comp P₂ F₂ hβF (isFrobeniusType_of_isIso P₂ e.inv)
  · refine hPS' ζ₀ ?_
    rw [hζ₀]
    exact IsPreStep.comp P₂ (isPreStep_of_isIso P₂ e.hom) hζ
  · refine hPS' ζ'₀ ?_
    rw [hζ'₀]
    exact IsPreStep.comp P₂ (isPreStep_of_isIso P₂ e.hom) hζ'
  · refine Ψ.functor.map_injective ?_
    rw [Ψ.functor.map_comp, hζ₀, hζ'₀, Category.assoc, ← hζeq]
  · refine Ψ.functor.map_injective ?_
    rw [Ψ.functor.map_comp, Ψ.functor.map_comp, hβ'₀, hζ₀, Category.assoc,
      ← Category.assoc e.inv, e.inv_hom_id, Category.id_comp, hsq]

def prop41iRhs_map_rev.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 78,
    item := "Theorem 4.2, (i) — Proposition 4.1, (i) の右辺は Ψ で戻る",
    sectionId := "frdi-thm-4-2" }


/-! ## ★★★★★★基底 —— Div-Frobenius-trivial な対象へ入る primary の両向き移送 -/

/-- ★★★★★★**基底(step 版)** —— 両側で `Proposition 4.1, (i)` を当てる。 -/
theorem isPrimaryElt_div_map_iff (G₁ : Frobenioid P₁) (G₂ : Frobenioid P₂)
    (F₂ : FrobenioidCore P₂) (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f →
      IsFrobeniusType P₂ (Ψ.functor.map f))
    (hFT' : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₂ (Ψ.functor.map f) →
      IsFrobeniusType P₁ f)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    {A B : C₁} (α : ℕ+ → (A ⟶ A))
    (hαd : ∀ n, P₁.degFr (α n) = n) (hαD : ∀ n, IsDivIdentity P₁ (α n))
    (hαF : ∀ n, IsFrobeniusType P₁ (α n))
    (hαd₂ : ∀ n, P₂.degFr (Ψ.functor.map (α n)) = n)
    (hαD₂ : ∀ n, IsDivIdentity P₂ (Ψ.functor.map (α n)))
    (hperf₁ : IsPerfectMonoid (Φ₁.val (P₁.toElem.obj A).base))
    (hperf₂ : IsPerfectMonoid (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base))
    (φ : B ⟶ A) (hcφ : IsCoAngular P₁ φ) (hst : IsStep P₁ φ)
    (hcφ₂ : IsCoAngular P₂ (Ψ.functor.map φ)) (hst₂ : IsStep P₂ (Ψ.functor.map φ)) :
    IsPrimaryElt (P₁.Div φ) ↔ IsPrimaryElt (P₂.Div (Ψ.functor.map φ)) :=
  ⟨fun h => isPrimaryElt_div_map G₁ G₂ F₂ Ψ hPS hPS' hFT hiso₁ hiso₂ α hαd hαD hαF
      hαd₂ hαD₂ hperf₁ hperf₂ φ hcφ hst hcφ₂ hst₂ h,
   fun h => (prop_4_1_i G₁ hiso₁ α hαd hαD hαF hperf₁ φ hcφ hst).mpr
      (prop41iRhs_map_rev F₂ Ψ hPS hPS' hFT' φ α
        ((prop_4_1_i G₂ hiso₂ (fun n => Ψ.functor.map (α n)) hαd₂ hαD₂
          (fun n => hFT (α n) (hαF n)) hperf₂ (Ψ.functor.map φ) hcφ₂ hst₂).mp h))⟩

/-- ★★★★★★**基底(pre-step 版)** —— 同型の場合は両辺とも偽なので自明。 -/
theorem primaryInIff_base (G₁ : Frobenioid P₁) (G₂ : Frobenioid P₂)
    (F₂ : FrobenioidCore P₂) (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f →
      IsFrobeniusType P₂ (Ψ.functor.map f))
    (hFT' : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₂ (Ψ.functor.map f) →
      IsFrobeniusType P₁ f)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    {A : C₁} (α : ℕ+ → (A ⟶ A))
    (hαd : ∀ n, P₁.degFr (α n) = n) (hαD : ∀ n, IsDivIdentity P₁ (α n))
    (hαF : ∀ n, IsFrobeniusType P₁ (α n))
    (hαd₂ : ∀ n, P₂.degFr (Ψ.functor.map (α n)) = n)
    (hαD₂ : ∀ n, IsDivIdentity P₂ (Ψ.functor.map (α n)))
    (hperf₁ : IsPerfectMonoid (Φ₁.val (P₁.toElem.obj A).base))
    (hperf₂ : IsPerfectMonoid (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base))
    (E : C₁) (ϵ : E ⟶ A) (hsϵ : IsPreStep P₁ ϵ) :
    IsPrimaryElt (P₁.Div ϵ) ↔ IsPrimaryElt (P₂.Div (Ψ.functor.map ϵ)) := by
  classical
  by_cases hst : IsStep P₁ ϵ
  · have hst₂ : IsStep P₂ (Ψ.functor.map ϵ) := by
      refine ⟨hPS ϵ hsϵ, fun hiso => hst.2 ?_⟩
      haveI := hiso
      exact isIso_of_reflects_iso ϵ Ψ.functor
    exact isPrimaryElt_div_map_iff G₁ G₂ F₂ Ψ hPS hPS' hFT hFT' hiso₁ hiso₂ α
      hαd hαD hαF hαd₂ hαD₂ hperf₁ hperf₂ ϵ (prop_1_4_i P₁ _ (fun Y _ => hiso₁ Y)) hst
      (prop_1_4_i P₂ _ (fun Y _ => hiso₂ Y)) hst₂
  · haveI hi : IsIso ϵ := by
      by_contra hc
      exact hst ⟨hsϵ, hc⟩
    haveI : IsIso (Ψ.functor.map ϵ) := Ψ.functor.map_isIso ϵ
    have h1 : P₁.Div ϵ = 0 := isIsometric_of_isIso P₁ ϵ
    have h2 : P₂.Div (Ψ.functor.map ϵ) = 0 := isIsometric_of_isIso P₂ _
    rw [h1, h2]
    exact ⟨fun hp => absurd rfl hp.1, fun hp => absurd rfl hp.1⟩

def primaryInIff_base.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 78,
    item := "Theorem 4.2, (i) — Div-Frobenius-trivial な対象での基底",
    sectionId := "frdi-thm-4-2" }

end ABC3.Found.FrdI
