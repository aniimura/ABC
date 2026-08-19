/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34EndBs
import ABC3.Found.FrdI.Prop14

/-!
# [FrdI] Theorem 3.4, (iv) —— `𝒪^▷` の保存を**すべての対象**へ広げる

原文 (FrdI p.63):
> (c) D1, D2 are Frobenius-slim. Then Ψ preserves the submonoids

★★在庫の `thm_3_4_iv_otri_map` は `C₀` が **Frobenius-trivial かつ
Frobenius-normalized** であることを要求していた。原文はそのような制限を置かない。

## ★★★★★測って分かった道 —— `OTriGenCond` を降ろす必要は無い

★以前は「`OTriGenCond` を Frobenius-trivial でない対象へ降ろす」道を探していたが、
`𝔽 = ElemFrob.Standard` は**次数の部分を持つ**ので、
その道は「span の頂点が Frobenius-trivial か」という(成り立たない)問いに帰着していた。

★★★**代わりに span をそのまま使えばよい**:

| 段 | 使うもの |
|---|---|
| 1 | `Definition 1.3, (i), (a)`(`baseSurj`)で Frobenius-trivial な `A₀` を取る |
| 2 | `Definition 1.3, (i), (b)`(`preStepSpan`)で span `X → A₀`、`X → C₀` を取る |
| 3 | isotropic 型では**すべての射が co-angular**(`prop_1_4_i`) |
| 4 | `Definition 1.3, (iii), (c)`(`otriBwd` / `otriFwd`)で `γ` を `X` 経由で `A₀` へ運ぶ |
| 5 | `A₀` では在庫の `thm_3_4_iv_otri_map` が使える |
| 6 | ★**pre-step の四角形で `𝒪^▷` は両方向に伝わる**(底の同型と次数の消去だけ) |

★★段 6 の消去は `oTri_of_otriEndCond` が使っていたものと同じ手筋である。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

section Cancel

variable {Dd : Type u} [Category.{v} Dd] {Cc : Type u2} [Category.{v2} Cc]
  {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)

/-- ★★★**pre-step の四角形で `𝒪^▷` は終域から始域へ降りる**。 -/
theorem oTri_dom_of_preStep_sq {A B : Cc} (φ : A ⟶ B) (hφ : IsPreStep P φ)
    {α : End A} {β : End B} (hβ : β ∈ OTri P B)
    (hsq : φ ≫ (β : B ⟶ B) = (α : A ⟶ A) ≫ φ) : α ∈ OTri P A := by
  haveI : IsIso (P.Base φ) := hφ.2
  have hb := congrArg P.Base hsq
  rw [P.Base_comp, P.Base_comp,
    show P.Base (β : B ⟶ B) = P.Base (𝟙 B) from hβ.1, P.Base_id, Category.comp_id] at hb
  have hd := congrArg P.degFr hsq
  rw [P.degFr_comp, P.degFr_comp, show P.degFr (β : B ⟶ B) = 1 from hβ.2, one_mul] at hd
  refine ⟨?_, ?_⟩
  · show P.Base (α : A ⟶ A) = P.Base (𝟙 A)
    rw [P.Base_id]
    have hb2 : (𝟙 ((P.toElem.obj A).base)) ≫ P.Base φ = P.Base (α : A ⟶ A) ≫ P.Base φ := by
      rw [Category.id_comp]; exact hb
    exact ((cancel_mono (P.Base φ)).mp hb2).symm
  · show P.degFr (α : A ⟶ A) = 1
    have hd2 : P.degFr φ * 1 = P.degFr φ * P.degFr (α : A ⟶ A) := by
      rw [mul_one]; exact hd
    exact (mul_left_cancel hd2).symm

/-- ★★★**pre-step の四角形で `𝒪^▷` は始域から終域へ上がる**。 -/
theorem oTri_cod_of_preStep_sq {A B : Cc} (φ : A ⟶ B) (hφ : IsPreStep P φ)
    {α : End A} {β : End B} (hα : α ∈ OTri P A)
    (hsq : φ ≫ (β : B ⟶ B) = (α : A ⟶ A) ≫ φ) : β ∈ OTri P B := by
  haveI : IsIso (P.Base φ) := hφ.2
  have hb := congrArg P.Base hsq
  rw [P.Base_comp, P.Base_comp,
    show P.Base (α : A ⟶ A) = P.Base (𝟙 A) from hα.1, P.Base_id, Category.id_comp] at hb
  have hd := congrArg P.degFr hsq
  rw [P.degFr_comp, P.degFr_comp, show P.degFr (α : A ⟶ A) = 1 from hα.2, mul_one] at hd
  refine ⟨?_, ?_⟩
  · show P.Base (β : B ⟶ B) = P.Base (𝟙 B)
    rw [P.Base_id]
    have hb2 : P.Base φ ≫ P.Base (β : B ⟶ B) = P.Base φ ≫ 𝟙 ((P.toElem.obj B).base) := by
      rw [Category.comp_id]; exact hb
    exact (cancel_epi (P.Base φ)).mp hb2
  · show P.degFr (β : B ⟶ B) = 1
    have hd2 : P.degFr (β : B ⟶ B) * P.degFr φ = 1 * P.degFr φ := by
      rw [one_mul]; exact hd
    exact mul_right_cancel hd2

def oTri_dom_of_preStep_sq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 66,
    item := "Theorem 3.4, (iv) — pre-step の四角形での 𝒪^▷ の伝播",
    sectionId := "frdi-thm-3-4" }

end Cancel

section Gen

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★★★**[FrdI] Theorem 3.4, (iv)** —— `Ψ` は**すべての対象**で `𝒪^▷` を保つ。

★★在庫版の `IsFrobeniusTrivial C₀` / `IsFrobeniusNormalized C₀` の仮定が落ちる。
★仮定として要るのは原文どおり
「`𝒞₁` は isotropic 型」「`𝒞₁` は Frobenius-normalized 型」「`𝒟₂` は Frobenius-slim」。 -/
theorem thm_3_4_iv_otri_map_general (F₁ : FrobenioidCore P₁) (Fc₂ : FrobenioidCore P₂)
    (hslim₂ : IsFrobeniusSlim D₂)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X)
    (hfn₁ : ∀ A : C₁, IsFrobeniusNormalized P₁ A)
    (Ψ : C₁ ⥤ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    (hGen : ∀ (A : C₁) (α : End A), OTriGenCond P₁ A α →
      OTriGenCond P₂ (Ψ.obj A) (Ψ.map (α : A ⟶ A)))
    {C₀ : C₁} (γ : End C₀) (h : γ ∈ OTri P₁ C₀) :
    Ψ.map (γ : C₀ ⟶ C₀) ∈ OTri P₂ (Ψ.obj C₀) := by
  -- 段 1: Frobenius-trivial な `A₀`
  obtain ⟨A₀, hft, ⟨e⟩⟩ := F₁.baseSurj ((P₁.toElem.obj C₀).base)
  -- 段 2: span
  obtain ⟨X, p, q, hp, hq, -⟩ := F₁.preStepSpan A₀ C₀ e.hom (by infer_instance)
  -- 段 3: isotropic 型ではすべての射が co-angular
  have hpc : IsCoAngular P₁ p := prop_1_4_i P₁ p (fun Y _ => hiso₁ Y)
  have hqc : IsCoAngular P₁ q := prop_1_4_i P₁ q (fun Y _ => hiso₁ Y)
  -- 段 4: `γ` を `X` へ降ろし、`A₀` へ上げる
  obtain ⟨α, ⟨hαm, hαsq⟩, -⟩ := F₁.otriBwd q hqc hq (γ : End C₀) h
  obtain ⟨β, ⟨hβm, hβsq⟩, -⟩ := F₁.otriFwd p hpc hp α hαm
  -- 段 5: `A₀` では在庫が使える
  have hΨβ : Ψ.map (β : A₀ ⟶ A₀) ∈ OTri P₂ (Ψ.obj A₀) :=
    thm_3_4_iv_otri_map Fc₂ hslim₂ Ψ hPS hGen hft (hfn₁ A₀) β hβm
  -- 段 6: pre-step の四角形で伝播させる
  have hΨα : Ψ.map (α : X ⟶ X) ∈ OTri P₂ (Ψ.obj X) :=
    oTri_dom_of_preStep_sq P₂ (Ψ.map p) (hPS p hp) hΨβ
      (by rw [← Ψ.map_comp, ← Ψ.map_comp, hβsq])
  exact oTri_cod_of_preStep_sq P₂ (Ψ.map q) (hPS q hq) hΨα
    (by rw [← Ψ.map_comp, ← Ψ.map_comp, hαsq])

def thm_3_4_iv_otri_map_general.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 63,
    item := "Theorem 3.4, (iv) — Ψ は 𝒪^▷ / 𝒪^× をすべての対象で保つ",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★★原文の仮定だけで `Ψ^un-tr` を作る -/

section UnTrGen

variable (e : C₁ ≌ C₂)
  (hPB : ∀ {X Y : C₁} (f : X ⟶ Y), IsPullBack P₁ f → IsPullBack P₂ (e.functor.map f))
  (hPB' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPullBack P₂ f → IsPullBack P₁ (e.inverse.map f))
  (hBI : ∀ {X Y : C₁} (f : X ⟶ Y),
    IsBaseIsomorphism P₁ f → IsBaseIsomorphism P₂ (e.functor.map f))
  (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (e.functor.map f))

include hPB hPB' hBI hPS in
/-- ★★★★★★**[FrdI] Theorem 3.4, (iv)** —— `Ψ` は `𝒪^▷` を保つ(仮定は原文どおり)。 -/
theorem thm_3_4_iv_otri_map_full (F₁ : FrobenioidCore P₁) (Fc₂ : FrobenioidCore P₂)
    (hslim₂ : IsFrobeniusSlim D₂)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X)
    (hfn₁ : ∀ A : C₁, IsFrobeniusNormalized P₁ A)
    (X : C₁) (γ : End X) (h : γ ∈ OTri P₁ X) :
    e.functor.map (γ : X ⟶ X) ∈ OTri P₂ (e.functor.obj X) :=
  thm_3_4_iv_otri_map_general F₁ Fc₂ hslim₂ hiso₁ hfn₁ e.functor hPS
    (fun A α hα => otriGenCond_map e hPB hPB' A hBI α hα) γ h

include hPB hPB' hBI hPS in
/-- ★★★★★★**[FrdI] Theorem 3.4, (iv)** —— `Ψ^un-tr` を原文の仮定だけで構成する。 -/
noncomputable def psiUnTrFull (F₁ : FrobenioidCore P₁) (Fc₂ : FrobenioidCore P₂)
    (hslim₂ : IsFrobeniusSlim D₂)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X)
    (hfn₁ : ∀ A : C₁, IsFrobeniusNormalized P₁ A)
    (h₁ : IsOfQuasiIsotropicType C₁ P₁) (h₂ : IsOfQuasiIsotropicType C₂ P₂) :
    UnTr P₁ ⥤ UnTr P₂ :=
  psiUnTr e.functor h₁ h₂
    (fun α₁ α₂ hh => toElem_map_congr_of_congr e.functor F₁ hiso₁
      (fun Xx δ hδ =>
        thm_3_4_iv_otri_map_full e hPB hPB' hBI hPS F₁ Fc₂ hslim₂ hiso₁ hfn₁ Xx δ hδ)
      α₁ α₂ hh)

def psiUnTrFull.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 63,
    item := "Theorem 3.4, (iv) — Ψ^un-tr（原文の仮定のみ）",
    sectionId := "frdi-thm-3-4" }

end UnTrGen

end Gen

end ABC3.Found.FrdI
