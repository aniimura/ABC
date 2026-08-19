/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop114
import ABC3.Found.FrdI.Thm34Pre
import Mathlib.Data.PNat.Factors

/-!
# [FrdI] Theorem 4.2, (i) —— `Div-identity` 自己射の保存

原文 (FrdI p.77):
> preserves primary steps, Div-identity endomorphisms, Div-Frobenius-

## ★★★`Proposition 1.14, (v)` が**圏論的特徴づけ**を与える(2026-08-19 に測った)

`prop_1_14_v`(`Prop114.lean:1648`):
```
(IsPrimeFrobenius φ ∧ IsDivIdentity φ) ↔
  ∀ B α, IsStep α → ∃ B' ψ β,
    IsIrreducibleMor ψ ∧ ¬ IsPreStep ψ ∧ IsStep β ∧ α ≫ ψ = φ ≫ α ≫ β
```

★★**右辺に `Div` も `Φ` も現れない** —— `IsStep`・`IsPreStep`・`IsIrreducibleMor`
だけである。★これらは `Theorem 3.4, (ii)(iii)` で `Ψ` が保つ/反射するので、
**右辺がそのまま移る**。

★したがって `Div-identity` の保存は
「pre-step の保存と反射」＋「`Proposition 1.14, (v)` の両側適用」で出る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 v3 u3 w3 u4 v4

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u3} [Category.{v3} D₂] {C₂ : Type u4} [Category.{v4} C₂]
  {Φ₂ : MonoidOn.{v3, u3, w3} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★`Proposition 1.14, (v)` の右辺 —— 純粋に圏論的な条件。 -/
def Prop114vRhs (P : PreFrobenioid C₁ Φ₁) {A : C₁} (φ : A ⟶ A) : Prop :=
  ∀ (B : C₁) (α : A ⟶ B), IsStep P α →
    ∃ (B' : C₁) (ψ : B ⟶ B') (β : B ⟶ B'),
      IsIrreducibleMor ψ ∧ ¬ IsPreStep P ψ ∧ IsStep P β ∧ α ≫ ψ = φ ≫ α ≫ β

/-- ★★同じものを第 2 の圏で。 -/
def Prop114vRhs' (P : PreFrobenioid C₂ Φ₂) {A : C₂} (φ : A ⟶ A) : Prop :=
  ∀ (B : C₂) (α : A ⟶ B), IsStep P α →
    ∃ (B' : C₂) (ψ : B ⟶ B') (β : B ⟶ B'),
      IsIrreducibleMor ψ ∧ ¬ IsPreStep P ψ ∧ IsStep P β ∧ α ≫ ψ = φ ≫ α ≫ β

/-- ★★★★★**右辺は `Ψ` で移る**(像の中の対象について)。

★`IsStep` / `IsPreStep` の保存と反射、`IsIrreducibleMor` の保存(圏論的)だけを使う。 -/
theorem prop114vRhs_map_image (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    {A : C₁} (φ : A ⟶ A) (h : Prop114vRhs P₁ φ)
    (B₀ : C₁) (α : Ψ.functor.obj A ⟶ Ψ.functor.obj B₀) (hα : IsStep P₂ α) :
    ∃ (B' : C₂) (ψ : Ψ.functor.obj B₀ ⟶ B') (β : Ψ.functor.obj B₀ ⟶ B'),
      IsIrreducibleMor ψ ∧ ¬ IsPreStep P₂ ψ ∧ IsStep P₂ β ∧
        α ≫ ψ = Ψ.functor.map φ ≫ α ≫ β := by
  obtain ⟨α₀, rfl⟩ := Ψ.functor.map_surjective α
  have hstep₀ : IsStep P₁ α₀ := by
    refine ⟨hPS' α₀ hα.1, fun hiso => hα.2 ?_⟩
    haveI := hiso
    exact Ψ.functor.map_isIso α₀
  obtain ⟨B'₀, ψ₀, β₀, hirr, hnps, hstepβ, hsq⟩ := h B₀ α₀ hstep₀
  refine ⟨Ψ.functor.obj B'₀, Ψ.functor.map ψ₀, Ψ.functor.map β₀,
    isIrreducibleMor_map Ψ.functor hirr, fun hc => hnps (hPS' ψ₀ hc),
    ⟨hPS β₀ hstepβ.1, fun hiso => hstepβ.2 ?_⟩, ?_⟩
  · haveI := hiso
    exact isIso_of_reflects_iso β₀ Ψ.functor
  · rw [← CategoryTheory.Functor.map_comp, ← CategoryTheory.Functor.map_comp,
      ← CategoryTheory.Functor.map_comp, hsq]

def prop114vRhs_map_image.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (i) — Proposition 1.14, (v) の右辺は Ψ で移る",
    sectionId := "frdi-thm-4-2" }

/-! ## ★2. 像の外へ広げる —— 本質的全射性で `B` を取り直す -/

/-- ★同型を前置しても irreducible。 -/
theorem isIrreducibleMor_iso_comp {X Y Z : C₂} (e : X ⟶ Y) [IsIso e] {f : Y ⟶ Z}
    (h : IsIrreducibleMor f) : IsIrreducibleMor (e ≫ f) := by
  refine ⟨fun hiso => h.1 ?_, ?_⟩
  · haveI := hiso
    have : f = inv e ≫ (e ≫ f) := by rw [← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
    rw [this]; infer_instance
  · intro W β α hfac
    have hfac' : f = (inv e ≫ β) ≫ α := by
      rw [Category.assoc, ← hfac, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
    rcases h.2 W (inv e ≫ β) α hfac' with hβ | hα
    · left
      haveI := hβ
      have : β = e ≫ (inv e ≫ β) := by
        rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
      rw [this]; infer_instance
    · right; exact hα

/-- ★★★★★**右辺は `𝒞₂` の全対象について移る** —— 本質的全射性で `B` を像へ引き戻す。 -/
theorem prop114vRhs_map (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    {A : C₁} (φ : A ⟶ A) (h : Prop114vRhs P₁ φ) :
    Prop114vRhs' P₂ (Ψ.functor.map φ) := by
  intro B α hα
  obtain ⟨B₀, ⟨ε⟩⟩ := Functor.EssSurj.mem_essImage (F := Ψ.functor) B
  have hstep : IsStep P₂ (α ≫ ε.inv) := by
    refine ⟨IsPreStep.comp P₂ hα.1 (isPreStep_of_isIso P₂ ε.inv), fun hiso => hα.2 ?_⟩
    haveI := hiso
    have : α = (α ≫ ε.inv) ≫ ε.hom := by
      rw [Category.assoc, ε.inv_hom_id, Category.comp_id]
    rw [this]; infer_instance
  obtain ⟨B', ψ, β, hirr, hnps, hstepβ, hsq⟩ :=
    prop114vRhs_map_image Ψ hPS hPS' φ h B₀ (α ≫ ε.inv) hstep
  refine ⟨B', ε.inv ≫ ψ, ε.inv ≫ β, isIrreducibleMor_iso_comp ε.inv hirr, ?_, ?_, ?_⟩
  · intro hc
    refine hnps ?_
    have h1 : ψ = ε.hom ≫ (ε.inv ≫ ψ) := by
      rw [← Category.assoc, ε.hom_inv_id, Category.id_comp]
    rw [h1]
    exact IsPreStep.comp P₂ (isPreStep_of_isIso P₂ ε.hom) hc
  · refine ⟨IsPreStep.comp P₂ (isPreStep_of_isIso P₂ ε.inv) hstepβ.1, fun hiso => hstepβ.2 ?_⟩
    haveI := hiso
    have : β = ε.hom ≫ (ε.inv ≫ β) := by
      rw [← Category.assoc, ε.hom_inv_id, Category.id_comp]
    rw [this]; infer_instance
  · calc α ≫ ε.inv ≫ ψ = (α ≫ ε.inv) ≫ ψ := (Category.assoc _ _ _).symm
      _ = Ψ.functor.map φ ≫ (α ≫ ε.inv) ≫ β := hsq
      _ = Ψ.functor.map φ ≫ α ≫ ε.inv ≫ β := by simp only [Category.assoc]

def prop114vRhs_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (i) — 右辺は 𝒞₂ の全対象について移る",
    sectionId := "frdi-thm-4-2" }

/-! ## ★3. ★★★★★★`Theorem 4.2, (i)` の `Div-identity` 保存

原文 (FrdI p.77):
> preserves primary steps, Div-identity endomorphisms, Div-Frobenius-

★★`Proposition 1.14, (v)` を**両側で適用する**だけ ——
`Div` も `Φ` も現れない右辺が `Ψ` で移るから。 -/

/-- ★★★★★★**[FrdI] Theorem 4.2, (i)** —— `Ψ` は
`prime-Frobenius かつ Div-identity` な自己射を保つ。 -/
theorem primeFrobenius_divIdentity_map (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hnd₁ : MonoidOn.IsNonDilatingOn Φ₁)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X) (hnd₂ : MonoidOn.IsNonDilatingOn Φ₂)
    {A : C₁} (hA₁ : ¬ IsGroupLikeObj P₁ A)
    (hA₂ : ¬ IsGroupLikeObj P₂ (Ψ.functor.obj A))
    (φ : A ⟶ A) (hirr : IsIrreducibleMor φ) (hnps : ¬ IsPreStep P₁ φ)
    (h : IsPrimeFrobenius P₁ φ ∧ IsDivIdentity P₁ φ) :
    IsPrimeFrobenius P₂ (Ψ.functor.map φ) ∧ IsDivIdentity P₂ (Ψ.functor.map φ) :=
  (prop_1_14_v P₂ F₂ G₂ hiso₂ hnd₂ hA₂ (Ψ.functor.map φ)
      (isIrreducibleMor_map Ψ.functor hirr) (fun hc => hnps (hPS' φ hc))).mpr
    (prop114vRhs_map Ψ hPS hPS' φ
      ((prop_1_14_v P₁ F₁ G₁ hiso₁ hnd₁ hA₁ φ hirr hnps).mp h))

def primeFrobenius_divIdentity_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (i) — Div-identity 自己射の保存",
    sectionId := "frdi-thm-4-2" }

/-! ## ★4. `Div-Frobenius-trivial` 対象の保存

原文 (FrdI p.77):
> trivial objects, and universally Div-Frobenius-trivial objects.

★★`IsDivFrobeniusTrivial A := ∃ ζ : ℕ≥1 →* End A, 次数が `n` で、各 `ζ n` が
`Div-identity` かつ Frobenius 型`。
★`ζ` は**単系準同型**なので、`Div-identity` を**素数の所だけ**言えば
素因数分解の帰納で全体へ延びる。 -/

/-- ★`Div-identity` は合成で閉じる。 -/
theorem isDivIdentity_comp {A : C₂} {φ ψ : A ⟶ A}
    (hφ : IsDivIdentity P₂ φ) (hψ : IsDivIdentity P₂ ψ) :
    IsDivIdentity P₂ (φ ≫ ψ) := by
  have h1 : Φ₂.map (P₂.Base φ) = Φ₂.map (P₂.Base (𝟙 A)) := hφ
  have h2 : Φ₂.map (P₂.Base ψ) = Φ₂.map (P₂.Base (𝟙 A)) := hψ
  rw [P₂.Base_id] at h1 h2
  show Φ₂.map (P₂.Base (φ ≫ ψ)) = Φ₂.map (P₂.Base (𝟙 A))
  rw [P₂.Base_id, P₂.Base_comp]
  refine AddMonoidHom.ext (fun x => ?_)
  rw [Φ₂.map_comp]
  have e2 : Φ₂.map (P₂.Base ψ) x = Φ₂.map (𝟙 ((P₂.toElem.obj A).base)) x :=
    congrArg (fun t : Φ₂.val _ →+ Φ₂.val _ => t x) h2
  rw [e2, Φ₂.map_id]
  have e1 : Φ₂.map (P₂.Base φ) x = Φ₂.map (𝟙 ((P₂.toElem.obj A).base)) x :=
    congrArg (fun t : Φ₂.val _ →+ Φ₂.val _ => t x) h1
  rw [e1, Φ₂.map_id]

/-- ★`𝟙` は `Div-identity`。 -/
theorem isDivIdentity_id (A : C₂) : IsDivIdentity P₂ (𝟙 A) := rfl

/-- ★★★★`Div-identity` な自己射のなす部分単系(`ℕ≥1` 側)。 -/
def divIdSubmonoid {A : C₂} (ζ : ℕ+ →* End A) : Submonoid ℕ+ where
  carrier := {n : ℕ+ | IsDivIdentity P₂ ((ζ n : End A) : A ⟶ A)}
  one_mem' := by
    show IsDivIdentity P₂ ((ζ 1 : End A) : A ⟶ A)
    rw [map_one]
    exact isDivIdentity_id A
  mul_mem' {m n} hm hn := by
    show IsDivIdentity P₂ ((ζ (m * n) : End A) : A ⟶ A)
    rw [map_mul]
    exact isDivIdentity_comp hn hm

/-- ★★★★★**素数の所だけ言えば全体へ延びる** —— `ℕ≥1` は素数で生成される。 -/
theorem isDivIdentity_of_primes {A : C₂} (ζ : ℕ+ →* End A)
    (hp : ∀ q : ℕ+, (q : ℕ).Prime → IsDivIdentity P₂ ((ζ q : End A) : A ⟶ A))
    (n : ℕ+) : IsDivIdentity P₂ ((ζ n : End A) : A ⟶ A) := by
  have hn : n ∈ divIdSubmonoid ζ := by
    rw [← PNat.prod_factorMultiset n]
    show ((PNat.factorMultiset n : Multiset ℕ+).prod) ∈ divIdSubmonoid ζ
    refine Submonoid.multiset_prod_mem _ _ (fun x hx => ?_)
    exact hp x (PrimeMultiset.coePNat_prime _ x hx)
  exact hn

def isDivIdentity_of_primes.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (i) — Div-identity は素数の所から全体へ延びる",
    sectionId := "frdi-thm-4-2" }

end ABC3.Found.FrdI
