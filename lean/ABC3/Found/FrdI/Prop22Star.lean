import ABC3.Found.FrdI.Prop21
import ABC3.Found.FrdI.Prop22

/-!
# [FrdI] Proposition 2.2, (ii) —— `𝒟*` 全体への拡張

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.45。

原文 (FrdI p.45):
> (ii) There is a unique contravariant functor

★`Found/FrdI/Prop22.lean` では `(𝒞^istr)^lin` の上で作った(`otriFunctorLin`)。
★**ここではそれを `𝒟*` 全体へ拡張する。**

## ★段取り

**段 A** —— `𝒟*` の射は「pre-step の逆 ∘ linear」で書ける(`dstar_span`)。

**段 B** —— その表示が well-defined であること。★**急所は
`otriLin_base_indep`**(linear 射に沿う `𝒪^▷` の引き戻しは `Base` にしか依らない)。

**段 C** —— 関手性と一意性。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★★段 B の急所 —— `otriLin` は `Base` にしか依らない

★`Definition 1.3, (iii), (c)` の「this bijection depends only on Base(φ)」は
**co-angular pre-step** についての主張である(`F.otriBase`)。
★★**それを linear 射へ広げる**のがここの仕事。

★**筋**: `Proposition 1.7, (iii)` で `ρ = β ≫ α`(pre-step ≫ pull-back)と分解し、
★**`α'` の pull-back 普遍性で `e : X ⟶ X'` を「底を指定して」取る**と
`β ≫ e` と `β'` が**底の等しい pre-step**になる。あとは `F.otriBase` が効く。

★**pull-back 部分 `α'` は共通に取れる**ので、`α'` に沿う引き戻しは
両辺で同じものになり、消約(`α'` は mono とは限らない)が要らない。 -/
theorem otriLin_base_indep {A B : C} (hA : IsIsotropic P A) {ρ ρ' : A ⟶ B}
    (hl : IsLinear P ρ) (hl' : IsLinear P ρ') (hb : P.Base ρ = P.Base ρ') :
    otriLin P F hA hl = otriLin P F hA hl' := by
  obtain ⟨X, β, α, hρ, hβ, hα⟩ := (prop_1_7_iii_linear_factor P F ρ).mp hl
  obtain ⟨X', β', α', hρ', hβ', hα'⟩ := (prop_1_7_iii_linear_factor P F ρ').mp hl'
  haveI hbβ : IsIso (P.Base β) := hβ.2
  haveI hbβ2 : IsIso (P.Base β') := hβ'.2
  have hX' : IsIsotropic P X' := F.isotropicClosed β' hA
  have hα'lin : IsLinear P α' := (F.pullBackLB α' hα').2
  -- 手 1: `α'` の pull-back 普遍性で `e : X ⟶ X'` を底を指定して作る
  have hbase : P.Base α = (inv (P.Base β) ≫ P.Base β') ≫ P.Base α' := by
    rw [Category.assoc, ← P.Base_comp, ← hρ', ← hb, hρ, P.Base_comp,
      ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  obtain ⟨e, he1, he2⟩ := isPullBack_lift P hα' α (inv (P.Base β) ≫ P.Base β') hbase
  have hel : IsLinear P e := by
    have h1 : P.degFr α = P.degFr α' * P.degFr e := by rw [← he1, P.degFr_comp]
    rw [show P.degFr α = 1 from (F.pullBackLB α hα).2,
      show P.degFr α' = 1 from hα'lin, one_mul] at h1
    exact h1.symm
  haveI hbe : IsIso (P.Base e) := by rw [he2]; infer_instance
  -- 手 2: `β ≫ e` は pre-step で、その底は `Base β'` に等しい
  have hβe : IsPreStep P (β ≫ e) := by
    refine ⟨IsLinear.comp P hβ.1 hel, ?_⟩
    show IsIso (P.Base (β ≫ e))
    rw [P.Base_comp]; infer_instance
  have hbβe : P.Base (β ≫ e) = P.Base β' := by
    rw [P.Base_comp, he2, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  have hρ2 : ρ = (β ≫ e) ≫ α' := by rw [hρ, Category.assoc, he1]
  -- 手 3: `α'` に沿う引き戻し `γ₁` は両辺で共通
  refine MonoidHom.ext fun γ => ?_
  set γ₁ : OTri P X' := otriLin P F hX' hα'lin γ with hγ₁
  have hsq₁ : α' ≫ ((γ : End B) : B ⟶ B) = ((γ₁ : End X') : X' ⟶ X') ≫ α' :=
    otriLin_spec P F hX' hα'lin γ
  set δ₁ : OTri P A := otriLin P F hA hβe.1 γ₁ with hδ₁
  have hsq₂ : (β ≫ e) ≫ ((γ₁ : End X') : X' ⟶ X')
      = ((δ₁ : End A) : A ⟶ A) ≫ (β ≫ e) :=
    otriLin_spec P F hA hβe.1 γ₁
  -- 手 4: `δ₁` は `ρ` 側の値でもある
  have hval : otriLin P F hA hl γ = δ₁ := by
    refine (otriLin_uniq P F hA hl γ δ₁ ?_).symm
    rw [hρ2, Category.assoc, hsq₁, ← Category.assoc, hsq₂, Category.assoc]
  -- 手 5: `F.otriBase` で `β'` 側へ移す
  have hcoa : IsCoAngular P (β ≫ e) := isCoAngular_of_isotropic_dom P F hA (β ≫ e)
  have hcoa' : IsCoAngular P β' := isCoAngular_of_isotropic_dom P F hA β'
  have hsq₃ : β' ≫ ((γ₁ : End X') : X' ⟶ X')
      = ((δ₁ : End A) : A ⟶ A) ≫ β' :=
    F.otriBase (β ≫ e) β' hcoa hβe hcoa' hβ' hbβe (δ₁ : End A) δ₁.2
      (γ₁ : End X') γ₁.2 hsq₂
  -- 手 6: `ρ'` 側の四角形を組み立てる
  rw [hval]
  refine otriLin_uniq P F hA hl' γ δ₁ ?_
  rw [hρ', Category.assoc, hsq₁, ← Category.assoc, hsq₃, Category.assoc]

/-! ## ★段 A —— `𝒟*` の射は「pre-step の逆 ∘ linear」で書ける

★`𝒟*` の対象は `𝒞^istr` の対象、射は**底の圏 `𝒟` の射**である。
★★**その射を `𝒞` へ持ち上げるのがここ** —— ただし `𝒞` の 1 本の射にはならず、
**`Z` からの span `(σ, ρ)`(`σ` は pre-step、`ρ` は linear)**になる。

★使う道具は 4 つ:
`plBk_realize`(`Definition 1.3, (i), (c)`)、
`F.isotropicHullExists`(同 (vii)(a))、
`F.preStepSpan`(同 (i)(b))、
`F.pullBackLB`(同 (iv)(b))。 -/
include F in
theorem dstar_span (A B : Istr P)
    (ψ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base) :
    ∃ (Z : C) (_ : IsIsotropic P Z) (σ : Z ⟶ A.obj) (_ : IsPreStep P σ)
      (ρ : Z ⟶ B.obj), IsLinear P ρ ∧ P.Base ρ = P.Base σ ≫ ψ := by
  -- 手 1: `ψ` を底に持つ pull-back を作る
  obtain ⟨X, π, hπ, θ, hbπ⟩ := plBk_realize P F B.obj ψ
  have hπlin : IsLinear P π := (F.pullBackLB π hπ).2
  -- 手 2: `X` の isotropic hull
  obtain ⟨Xi, ξ, hξ⟩ := F.isotropicHullExists X
  haveI hbξ : IsIso (P.Base ξ) := hξ.2.1.2
  -- 手 3: `B.obj` は isotropic なので `π` は `ξ` を経由する
  obtain ⟨π', hπ'e, -⟩ := hξ.2.2.2 B.obj B.property π
  have hπ'lin : IsLinear P π' := by
    have h1 : P.degFr π = P.degFr π' * P.degFr ξ := by rw [hπ'e, P.degFr_comp]
    rw [show P.degFr π = 1 from hπlin, show P.degFr ξ = 1 from hξ.2.1.1, mul_one] at h1
    exact h1.symm
  -- 手 4: `Definition 1.3, (i), (b)` で span を作る
  haveI hiα : IsIso (θ.inv ≫ P.Base ξ) := inferInstance
  obtain ⟨Z₀, σ₀, τ₀, hσ₀, hτ₀, hspan⟩ :=
    F.preStepSpan A.obj Xi (θ.inv ≫ P.Base ξ) hiα
  haveI hbσ₀ : IsIso (P.Base σ₀) := hσ₀.2
  haveI hbτ₀ : IsIso (P.Base τ₀) := hτ₀.2
  -- 手 5: `Z₀` の hull を取って isotropic にする
  obtain ⟨Z, ζ, hζ⟩ := F.isotropicHullExists Z₀
  haveI hbζ : IsIso (P.Base ζ) := hζ.2.1.2
  obtain ⟨σ, hσe, -⟩ := hζ.2.2.2 A.obj A.property σ₀
  obtain ⟨τ, hτe, -⟩ := hζ.2.2.2 Xi hξ.2.2.1 τ₀
  have hζlin : P.degFr ζ = 1 := hζ.2.1.1
  have hσlin : IsLinear P σ := by
    have h1 : P.degFr σ₀ = P.degFr σ * P.degFr ζ := by rw [hσe, P.degFr_comp]
    rw [show P.degFr σ₀ = 1 from hσ₀.1, hζlin, mul_one] at h1
    exact h1.symm
  have hτlin : IsLinear P τ := by
    have h1 : P.degFr τ₀ = P.degFr τ * P.degFr ζ := by rw [hτe, P.degFr_comp]
    rw [show P.degFr τ₀ = 1 from hτ₀.1, hζlin, mul_one] at h1
    exact h1.symm
  have hbσe : P.Base ζ ≫ P.Base σ = P.Base σ₀ := by rw [← P.Base_comp, ← hσe]
  have hbτe : P.Base ζ ≫ P.Base τ = P.Base τ₀ := by rw [← P.Base_comp, ← hτe]
  haveI hbσ : IsIso (P.Base σ) := by
    haveI : IsIso (P.Base ζ ≫ P.Base σ) := by rw [hbσe]; infer_instance
    exact IsIso.of_isIso_comp_left (P.Base ζ) _
  haveI hbτ : IsIso (P.Base τ) := by
    haveI : IsIso (P.Base ζ ≫ P.Base τ) := by rw [hbτe]; infer_instance
    exact IsIso.of_isIso_comp_left (P.Base ζ) _
  -- 手 6: `ρ := τ ≫ π'`
  refine ⟨Z, hζ.2.2.1, σ, ⟨hσlin, hbσ⟩, τ ≫ π', IsLinear.comp P hτlin hπ'lin, ?_⟩
  -- 底の計算
  have h0 : P.Base σ₀ ≫ (θ.inv ≫ P.Base ξ) = P.Base τ₀ := by
    rw [hspan, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  have hb1 : P.Base τ = P.Base σ ≫ θ.inv ≫ P.Base ξ := by
    haveI : Epi (P.Base ζ) := inferInstance
    refine ((cancel_epi (P.Base ζ)).mp ?_).symm
    rw [← Category.assoc, hbσe, h0, hbτe]
  have hb2 : P.Base ξ ≫ P.Base π' = θ.hom ≫ ψ := by
    rw [← P.Base_comp, ← hπ'e]; exact hbπ
  rw [P.Base_comp, hb1, Category.assoc, Category.assoc, hb2]
  simp

/-! ## ★段 B —— span の表す写像は well-defined

★2 つの span が同じ `ψ` を与えるなら、定める写像は一致する。

★**筋**: `Definition 1.3, (i), (b)`(`F.preStepSpan`)で**共通細分** `W` を作り、
hull を取って isotropic にする。★**`W` から見ると 2 つの span は
「底の等しい linear 射」になる**ので `otriLin_base_indep` が効く。
★最後は `σ` 側の `otriLin` が単射(pre-step)であることで締める。 -/
include F in
theorem dstar_map_unique (A B : Istr P)
    (ψ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base)
    {Z₁ : C} (h₁ : IsIsotropic P Z₁) {σ₁ : Z₁ ⟶ A.obj} (hσ₁ : IsPreStep P σ₁)
    {ρ₁ : Z₁ ⟶ B.obj} (hρ₁ : IsLinear P ρ₁) (hb₁ : P.Base ρ₁ = P.Base σ₁ ≫ ψ)
    {Z₂ : C} (h₂ : IsIsotropic P Z₂) {σ₂ : Z₂ ⟶ A.obj} (hσ₂ : IsPreStep P σ₂)
    {ρ₂ : Z₂ ⟶ B.obj} (hρ₂ : IsLinear P ρ₂) (hb₂ : P.Base ρ₂ = P.Base σ₂ ≫ ψ)
    (f₁ f₂ : OTri P B.obj →* OTri P A.obj)
    (hf₁ : (otriLin P F h₁ hσ₁.1).comp f₁ = otriLin P F h₁ hρ₁)
    (hf₂ : (otriLin P F h₂ hσ₂.1).comp f₂ = otriLin P F h₂ hρ₂) :
    f₁ = f₂ := by
  haveI hbσ₁ : IsIso (P.Base σ₁) := hσ₁.2
  haveI hbσ₂ : IsIso (P.Base σ₂) := hσ₂.2
  -- 手 1: 共通細分を作り、hull を取って isotropic にする
  obtain ⟨W₀, μ₁₀, μ₂₀, hμ₁₀, hμ₂₀, hspan⟩ :=
    F.preStepSpan Z₁ Z₂ (P.Base σ₁ ≫ inv (P.Base σ₂)) inferInstance
  haveI hbμ₁₀ : IsIso (P.Base μ₁₀) := hμ₁₀.2
  haveI hbμ₂₀ : IsIso (P.Base μ₂₀) := hμ₂₀.2
  obtain ⟨W, ν, hν⟩ := F.isotropicHullExists W₀
  haveI hbν : IsIso (P.Base ν) := hν.2.1.2
  have hνlin : P.degFr ν = 1 := hν.2.1.1
  have hW : IsIsotropic P W := hν.2.2.1
  obtain ⟨μ₁, hμ₁e, -⟩ := hν.2.2.2 Z₁ h₁ μ₁₀
  obtain ⟨μ₂, hμ₂e, -⟩ := hν.2.2.2 Z₂ h₂ μ₂₀
  -- 手 2: `μᵢ` は pre-step
  have hbμ₁e : P.Base ν ≫ P.Base μ₁ = P.Base μ₁₀ := by rw [← P.Base_comp, ← hμ₁e]
  have hbμ₂e : P.Base ν ≫ P.Base μ₂ = P.Base μ₂₀ := by rw [← P.Base_comp, ← hμ₂e]
  have hμ₁ : IsPreStep P μ₁ := by
    refine ⟨?_, ?_⟩
    · have h1 : P.degFr μ₁₀ = P.degFr μ₁ * P.degFr ν := by rw [hμ₁e, P.degFr_comp]
      rw [show P.degFr μ₁₀ = 1 from hμ₁₀.1, hνlin, mul_one] at h1
      exact h1.symm
    · show IsIso (P.Base μ₁)
      haveI : IsIso (P.Base ν ≫ P.Base μ₁) := by rw [hbμ₁e]; infer_instance
      exact IsIso.of_isIso_comp_left (P.Base ν) _
  have hμ₂ : IsPreStep P μ₂ := by
    refine ⟨?_, ?_⟩
    · have h1 : P.degFr μ₂₀ = P.degFr μ₂ * P.degFr ν := by rw [hμ₂e, P.degFr_comp]
      rw [show P.degFr μ₂₀ = 1 from hμ₂₀.1, hνlin, mul_one] at h1
      exact h1.symm
    · show IsIso (P.Base μ₂)
      haveI : IsIso (P.Base ν ≫ P.Base μ₂) := by rw [hbμ₂e]; infer_instance
      exact IsIso.of_isIso_comp_left (P.Base ν) _
  -- 手 3: 底が一致すること
  have hspan' : P.Base μ₁₀ ≫ P.Base σ₁ = P.Base μ₂₀ ≫ P.Base σ₂ := by
    have h : P.Base μ₁₀ ≫ (P.Base σ₁ ≫ inv (P.Base σ₂)) = P.Base μ₂₀ := by
      rw [hspan, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
    rw [← h, Category.assoc, Category.assoc, IsIso.inv_hom_id, Category.comp_id]
  have hbμσ : P.Base (μ₁ ≫ σ₁) = P.Base (μ₂ ≫ σ₂) := by
    rw [P.Base_comp, P.Base_comp]
    refine (cancel_epi (P.Base ν)).mp ?_
    rw [← Category.assoc, ← Category.assoc, hbμ₁e, hbμ₂e]
    exact hspan'
  have hbμρ : P.Base (μ₁ ≫ ρ₁) = P.Base (μ₂ ≫ ρ₂) := by
    rw [P.Base_comp, P.Base_comp, hb₁, hb₂, ← Category.assoc, ← Category.assoc]
    exact congrArg (· ≫ ψ) (by rw [← P.Base_comp, ← P.Base_comp]; exact hbμσ)
  -- 手 4: `otriLin_base_indep`
  have hμσ₁ : IsPreStep P (μ₁ ≫ σ₁) := IsPreStep.comp P hμ₁ hσ₁
  have hμσ₂ : IsPreStep P (μ₂ ≫ σ₂) := IsPreStep.comp P hμ₂ hσ₂
  have e1 : otriLin P F hW hμσ₁.1 = otriLin P F hW hμσ₂.1 :=
    otriLin_base_indep P F hW hμσ₁.1 hμσ₂.1 hbμσ
  have e2 : otriLin P F hW (IsLinear.comp P hμ₁.1 hρ₁)
      = otriLin P F hW (IsLinear.comp P hμ₂.1 hρ₂) :=
    otriLin_base_indep P F hW _ _ hbμρ
  -- 手 5: 分解して単射性で締める
  have hd₁ : (otriLin P F hW hμσ₁.1).comp f₁
      = otriLin P F hW (IsLinear.comp P hμ₁.1 hρ₁) := by
    rw [otriLin_comp P F hW h₁ hμ₁.1 hσ₁.1, otriLin_comp P F hW h₁ hμ₁.1 hρ₁,
      MonoidHom.comp_assoc, hf₁]
  have hd₂ : (otriLin P F hW hμσ₂.1).comp f₂
      = otriLin P F hW (IsLinear.comp P hμ₂.1 hρ₂) := by
    rw [otriLin_comp P F hW h₂ hμ₂.1 hσ₂.1, otriLin_comp P F hW h₂ hμ₂.1 hρ₂,
      MonoidHom.comp_assoc, hf₂]
  refine MonoidHom.ext fun γ => ?_
  refine otriLin_injective P F hW hμσ₁.1 ?_
  have k1 : (otriLin P F hW hμσ₁.1) (f₁ γ) = (otriLin P F hW hμσ₁.1) (f₂ γ) := by
    have l1 : (otriLin P F hW hμσ₁.1) (f₁ γ)
        = otriLin P F hW (IsLinear.comp P hμ₁.1 hρ₁) γ := DFunLike.congr_fun hd₁ γ
    have l2 : (otriLin P F hW hμσ₁.1) (f₂ γ)
        = otriLin P F hW (IsLinear.comp P hμ₂.1 hρ₂) γ := by
      rw [e1]; exact DFunLike.congr_fun hd₂ γ
    rw [l1, l2, e2]
  exact k1

/-! ## ★段 C —— `𝒟*` 上の関手

★span `(Z, σ, ρ)` が定める写像は `(otriLin σ)⁻¹ ∘ otriLin ρ` である。
★`σ` は pre-step なので `otriLin σ` は**全単射**(`otriLin_bijective_of_preStep`)。 -/

/-- ★pre-step に沿う `𝒪^▷` の引き戻しは**同型**(`Definition 1.3, (iii), (c)`)。 -/
noncomputable def otriPreStepEquiv {A B : C} (hA : IsIsotropic P A) {φ : A ⟶ B}
    (hs : IsPreStep P φ) : OTri P B ≃* OTri P A :=
  MulEquiv.ofBijective (otriLin P F hA hs.1) (otriLin_bijective_of_preStep P F hA hs)

/-- ★span が定める写像。 -/
noncomputable def spanMap {A B Z : C} (hZ : IsIsotropic P Z)
    {σ : Z ⟶ A} (hσ : IsPreStep P σ) {ρ : Z ⟶ B} (hρ : IsLinear P ρ) :
    OTri P B →* OTri P A :=
  ((otriPreStepEquiv P F hZ hσ).symm.toMonoidHom).comp (otriLin P F hZ hρ)

theorem spanMap_spec {A B Z : C} (hZ : IsIsotropic P Z)
    {σ : Z ⟶ A} (hσ : IsPreStep P σ) {ρ : Z ⟶ B} (hρ : IsLinear P ρ) :
    (otriLin P F hZ hσ.1).comp (spanMap P F hZ hσ hρ) = otriLin P F hZ hρ := by
  refine MonoidHom.ext fun γ => ?_
  exact (otriPreStepEquiv P F hZ hσ).apply_symm_apply _

/-- ★★`𝒟*` の射 `ψ` を表す写像であるという条件 —— ★**どの span で見ても同じ**。 -/
def IsDStarMap (A B : Istr P)
    (ψ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base)
    (f : OTri P B.obj →* OTri P A.obj) : Prop :=
  ∀ (Z : C) (hZ : IsIsotropic P Z) (σ : Z ⟶ A.obj) (hσ : IsPreStep P σ)
    (ρ : Z ⟶ B.obj) (hρ : IsLinear P ρ), P.Base ρ = P.Base σ ≫ ψ →
    (otriLin P F hZ hσ.1).comp f = otriLin P F hZ hρ

include F in
/-- ★★★**段 A ＋ 段 B の結論** —— `𝒟*` の各射に対し、条件を満たす写像が
**ただ 1 つ**存在する。 -/
theorem dstar_map_exists_unique (A B : Istr P)
    (ψ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base) :
    ∃! f : OTri P B.obj →* OTri P A.obj, IsDStarMap P F A B ψ f := by
  obtain ⟨Z, hZ, σ, hσ, ρ, hρ, hb⟩ := dstar_span P F A B ψ
  refine ⟨spanMap P F hZ hσ hρ, ?_, ?_⟩
  · intro Z' hZ' σ' hσ' ρ' hρ' hb'
    have heq := dstar_map_unique P F A B ψ hZ hσ hρ hb hZ' hσ' hρ' hb'
      (spanMap P F hZ hσ hρ) (spanMap P F hZ' hσ' hρ')
      (spanMap_spec P F hZ hσ hρ) (spanMap_spec P F hZ' hσ' hρ')
    rw [heq]
    exact spanMap_spec P F hZ' hσ' hρ'
  · intro g hg
    exact dstar_map_unique P F A B ψ hZ hσ hρ hb hZ hσ hρ hb g
      (spanMap P F hZ hσ hρ) (hg Z hZ σ hσ ρ hρ hb) (spanMap_spec P F hZ hσ hρ)

/-- ★★`𝒟*` の射が誘導する `𝒪^▷` の準同型。 -/
noncomputable def dstarMap (A B : Istr P)
    (ψ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base) :
    OTri P B.obj →* OTri P A.obj :=
  (dstar_map_exists_unique P F A B ψ).choose

theorem dstarMap_spec (A B : Istr P)
    (ψ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base) :
    IsDStarMap P F A B ψ (dstarMap P F A B ψ) :=
  (dstar_map_exists_unique P F A B ψ).choose_spec.1

theorem dstarMap_uniq (A B : Istr P)
    (ψ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base)
    (f : OTri P B.obj →* OTri P A.obj) (hf : IsDStarMap P F A B ψ f) :
    f = dstarMap P F A B ψ :=
  (dstar_map_exists_unique P F A B ψ).choose_spec.2 f hf

include F in
/-- ★★**(a)** `𝒞^istr` の linear 射に沿う値は `otriLin` そのもの
(＝`Proposition 1.11, (iv)` の埋め込み) —— 原文の (a)。 -/
theorem dstarMap_eq_otriLin (A B : Istr P) {ρ : A.obj ⟶ B.obj} (hρ : IsLinear P ρ) :
    dstarMap P F A B (P.Base ρ) = otriLin P F A.property hρ := by
  refine (dstarMap_uniq P F A B (P.Base ρ) _ ?_).symm
  intro Z hZ σ hσ ρ' hρ' hb'
  rw [← otriLin_comp P F hZ A.property hσ.1 hρ]
  exact otriLin_base_indep P F hZ _ _ (by rw [P.Base_comp]; exact hb'.symm)

/-- ★**恒等射では恒等準同型**。 -/
theorem dstarMap_id (A : Istr P) :
    dstarMap P F A A (𝟙 ((P.toElem.obj A.obj).base)) = MonoidHom.id (OTri P A.obj) := by
  have h := dstarMap_eq_otriLin P F A A (isLinear_id P A.obj)
  rw [P.Base_id] at h
  rw [h, otriLin_id]

variable (G : Frobenioid P)

include F G in
/-- ★★**合成を(反変に)保つ**。

★**2 つの span を繋ぐのに `Proposition 1.11, (vii)`** が要る ——
`ρ₁` に沿って `σ₂` を引き戻し、共通の `W` を作る。
★hull を取って `W` を isotropic にするのは `otriLin` の定義域条件のため。 -/
theorem dstarMap_comp (A B E : Istr P)
    (ψ₁ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base)
    (ψ₂ : (P.toElem.obj B.obj).base ⟶ (P.toElem.obj E.obj).base) :
    dstarMap P F A E (ψ₁ ≫ ψ₂)
      = (dstarMap P F A B ψ₁).comp (dstarMap P F B E ψ₂) := by
  obtain ⟨Z₁, hZ₁, σ₁, hσ₁, ρ₁, hρ₁, hb₁⟩ := dstar_span P F A B ψ₁
  obtain ⟨Z₂, hZ₂, σ₂, hσ₂, ρ₂, hρ₂, hb₂⟩ := dstar_span P F B E ψ₂
  -- 手 1: `ρ₁` に沿って `σ₂` を引き戻す
  obtain ⟨W₀, γ₀, α₀, -, hγ₀s, hcomm₀⟩ :=
    prop_1_11_vii P F G ρ₁ σ₂ (isCoAngular_of_isotropic_dom P F hZ₂ σ₂) hσ₂
  -- 手 2: hull を取って isotropic にする
  obtain ⟨W, ν, hν⟩ := F.isotropicHullExists W₀
  haveI hbν : IsIso (P.Base ν) := hν.2.1.2
  have hW : IsIsotropic P W := hν.2.2.1
  obtain ⟨γ, hγe, -⟩ := hν.2.2.2 Z₁ hZ₁ γ₀
  obtain ⟨α, hαe, -⟩ := hν.2.2.2 Z₂ hZ₂ α₀
  have hγ : IsPreStep P γ := by
    refine ⟨?_, ?_⟩
    · have h1 : P.degFr γ₀ = P.degFr γ * P.degFr ν := by rw [hγe, P.degFr_comp]
      rw [show P.degFr γ₀ = 1 from hγ₀s.1, show P.degFr ν = 1 from hν.2.1.1,
        mul_one] at h1
      exact h1.symm
    · show IsIso (P.Base γ)
      haveI : IsIso (P.Base ν ≫ P.Base γ) := by
        rw [← P.Base_comp, ← hγe]; exact hγ₀s.2
      exact IsIso.of_isIso_comp_left (P.Base ν) _
  have hcomm : γ ≫ ρ₁ = α ≫ σ₂ := by
    haveI : Epi ν := P.totEpiC _ _ _
    refine (cancel_epi ν).mp ?_
    rw [← Category.assoc, ← hγe, ← Category.assoc, ← hαe]
    exact hcomm₀
  have hαlin : IsLinear P α := by
    have h1 : P.degFr (γ ≫ ρ₁) = P.degFr (α ≫ σ₂) := by rw [hcomm]
    rw [P.degFr_comp, P.degFr_comp, show P.degFr ρ₁ = 1 from hρ₁,
      show P.degFr γ = 1 from hγ.1, show P.degFr σ₂ = 1 from hσ₂.1, one_mul,
      one_mul] at h1
    exact h1.symm
  -- 手 3: 合成 span の底
  have hbγσ : P.Base (α ≫ ρ₂) = P.Base (γ ≫ σ₁) ≫ (ψ₁ ≫ ψ₂) := by
    have e1 : P.Base (γ ≫ ρ₁) = P.Base (α ≫ σ₂) := by rw [hcomm]
    rw [P.Base_comp, P.Base_comp] at e1
    rw [P.Base_comp, hb₂, ← Category.assoc, ← e1, P.Base_comp, hb₁,
      Category.assoc, Category.assoc, Category.assoc]
  -- 手 4: 右辺が合成 span の条件を満たす
  have hγσ : IsPreStep P (γ ≫ σ₁) := IsPreStep.comp P hγ hσ₁
  have s1 := dstarMap_spec P F A B ψ₁ Z₁ hZ₁ σ₁ hσ₁ ρ₁ hρ₁ hb₁
  have s2 := dstarMap_spec P F B E ψ₂ Z₂ hZ₂ σ₂ hσ₂ ρ₂ hρ₂ hb₂
  have hbi := otriLin_base_indep P F hW (IsLinear.comp P hγ.1 hρ₁)
    (IsLinear.comp P hαlin hσ₂.1) (by rw [hcomm])
  have hright : (otriLin P F hW hγσ.1).comp
      ((dstarMap P F A B ψ₁).comp (dstarMap P F B E ψ₂))
      = otriLin P F hW (IsLinear.comp P hαlin hρ₂) := by
    refine MonoidHom.ext fun δ => ?_
    have e1 : otriLin P F hW hγσ.1 (dstarMap P F A B ψ₁ (dstarMap P F B E ψ₂ δ))
        = otriLin P F hW hγ.1
          (otriLin P F hZ₁ hσ₁.1 (dstarMap P F A B ψ₁ (dstarMap P F B E ψ₂ δ))) :=
      DFunLike.congr_fun (otriLin_comp P F hW hZ₁ hγ.1 hσ₁.1) _
    have e2 : otriLin P F hZ₁ hσ₁.1 (dstarMap P F A B ψ₁ (dstarMap P F B E ψ₂ δ))
        = otriLin P F hZ₁ hρ₁ (dstarMap P F B E ψ₂ δ) :=
      DFunLike.congr_fun s1 _
    have e3 : otriLin P F hW hγ.1 (otriLin P F hZ₁ hρ₁ (dstarMap P F B E ψ₂ δ))
        = otriLin P F hW (IsLinear.comp P hγ.1 hρ₁) (dstarMap P F B E ψ₂ δ) :=
      (DFunLike.congr_fun (otriLin_comp P F hW hZ₁ hγ.1 hρ₁) _).symm
    have e4 : otriLin P F hW (IsLinear.comp P hαlin hσ₂.1) (dstarMap P F B E ψ₂ δ)
        = otriLin P F hW hαlin (otriLin P F hZ₂ hσ₂.1 (dstarMap P F B E ψ₂ δ)) :=
      DFunLike.congr_fun (otriLin_comp P F hW hZ₂ hαlin hσ₂.1) _
    have e5 : otriLin P F hZ₂ hσ₂.1 (dstarMap P F B E ψ₂ δ)
        = otriLin P F hZ₂ hρ₂ δ := DFunLike.congr_fun s2 _
    have e6 : otriLin P F hW hαlin (otriLin P F hZ₂ hρ₂ δ)
        = otriLin P F hW (IsLinear.comp P hαlin hρ₂) δ :=
      (DFunLike.congr_fun (otriLin_comp P F hW hZ₂ hαlin hρ₂) _).symm
    show otriLin P F hW hγσ.1 (dstarMap P F A B ψ₁ (dstarMap P F B E ψ₂ δ)) = _
    rw [e1, e2, e3, DFunLike.congr_fun hbi (dstarMap P F B E ψ₂ δ), e4, e5, e6]
  -- 手 5: 段 B で締める
  exact dstar_map_unique P F A E (ψ₁ ≫ ψ₂) hW hγσ (IsLinear.comp P hαlin hρ₂) hbγσ
    hW hγσ (IsLinear.comp P hαlin hρ₂) hbγσ _ _
    (dstarMap_spec P F A E (ψ₁ ≫ ψ₂) W hW (γ ≫ σ₁) hγσ (α ≫ ρ₂)
      (IsLinear.comp P hαlin hρ₂) hbγσ) hright

/-! ### ★★★関手として束ねる —— `𝒟* → Mon` -/

include F G in
/-- ★★★**[FrdI] Proposition 2.2, (ii)** —— 反変関手 `𝒪^▷(−) : 𝒟* → Mon`。

原文 (FrdI p.45):
> (ii) There is a unique contravariant functor

★**対象は `𝒞^istr` の対象、射は底の圏 `𝒟` の射**である。
★★**`(𝒞^istr)^lin` 上の `otriFunctorLin` の拡張になっている**
(`dstarMap_eq_otriLin` が (a) を、`otriLin_bijective_of_preStep` が (b) を与える)。 -/
noncomputable def otriStar : (InducedCategory D (istrBase P))ᵒᵖ ⥤ MonCat.{v2} where
  obj A := MonCat.of (OTri P (A.unop : Istr P).obj)
  map {A B} f := MonCat.ofHom (dstarMap P F (B.unop : Istr P) (A.unop : Istr P)
    ((dStarToD P).map f.unop))
  map_id A := by
    apply MonCat.hom_ext
    exact dstarMap_id P F (A.unop : Istr P)
  map_comp {A B E} f g := by
    apply MonCat.hom_ext
    exact dstarMap_comp P F G (E.unop : Istr P) (B.unop : Istr P) (A.unop : Istr P)
      ((dStarToD P).map g.unop) ((dStarToD P).map f.unop)

include F in
/-- ★★**(b)** pre-step の底での値は**全単射** ——
`Definition 1.3, (iii), (c)` の全単射そのもの(原文の (b))。 -/
theorem dstarMap_bijective_of_preStep (A B : Istr P) {σ : A.obj ⟶ B.obj}
    (hσ : IsPreStep P σ) : Function.Bijective (dstarMap P F A B (P.Base σ)) := by
  rw [dstarMap_eq_otriLin P F A B hσ.1]
  exact otriLin_bijective_of_preStep P F A.property hσ

include F in
/-- ★★★**一意性**(原文の "There is a **unique** contravariant functor")。

★**`𝒞^istr` の linear 射での値(＝原文の (a))と関手性だけで、
`𝒟*` の全射での値が決まってしまう。**

★**理由**: `𝒟*` の射はすべて span `(Z, σ, ρ)` で書ける(`dstar_span`)。
`ψ = (Base σ)⁻¹ ≫ Base ρ` なので、`Base σ` と `Base ρ` での値が決まれば
関手性から `ψ` での値が決まる。 -/
theorem dstarMap_uniq_of_linear
    (m : ∀ X Y : Istr P,
      ((P.toElem.obj X.obj).base ⟶ (P.toElem.obj Y.obj).base) →
      (OTri P Y.obj →* OTri P X.obj))
    (hlin : ∀ (X Y : Istr P) (ρ : X.obj ⟶ Y.obj) (hρ : IsLinear P ρ),
      m X Y (P.Base ρ) = otriLin P F X.property hρ)
    (hcomp : ∀ (X Y Z : Istr P)
      (u : (P.toElem.obj X.obj).base ⟶ (P.toElem.obj Y.obj).base)
      (v : (P.toElem.obj Y.obj).base ⟶ (P.toElem.obj Z.obj).base),
      m X Z (u ≫ v) = (m X Y u).comp (m Y Z v))
    (A B : Istr P)
    (ψ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base) :
    m A B ψ = dstarMap P F A B ψ := by
  refine dstarMap_uniq P F A B ψ _ ?_
  intro Z hZ σ hσ ρ hρ hb
  have h1 : m ⟨Z, hZ⟩ B (P.Base ρ) = (m ⟨Z, hZ⟩ A (P.Base σ)).comp (m A B ψ) := by
    rw [← hcomp ⟨Z, hZ⟩ A B (P.Base σ) ψ, hb]
  rw [hlin ⟨Z, hZ⟩ B ρ hρ, hlin ⟨Z, hZ⟩ A σ hσ.1] at h1
  exact h1.symm

/-! ## ★(iii) の前半 —— `𝒪^×` は部分関手であり `𝒪^▷^±` に等しい

原文 (FrdI p.45):
> functor of the functor of (ii) which is equal to the subfunctor
-/

/-- ★★**(iii)-1** `𝒪^×(A)` は `𝒪^▷(A)` の**単元群**に一致する
(原文の `𝒪^▷(A)^±`)。

★`⟸` は「部分モノイドの単元は元の単元」。
★`⟹` は「`End A` で可逆な `𝒪^▷` の元は、その逆射も `𝒪^▷` に入る」——
底は `Base (inv x) ≫ Base x = 𝟙` から、次数は `degFr_of_isIso` から。 -/
theorem otimes_eq_units {A : C} (x : End A) :
    x ∈ OTimes P A ↔ ∃ u : (OTri P A)ˣ, ((u : OTri P A) : End A) = x := by
  constructor
  · intro hx
    haveI hiso : IsIso ((x : End A) : A ⟶ A) :=
      (CategoryTheory.isUnit_iff_isIso _).mp hx.2
    have hbx : P.Base ((x : End A) : A ⟶ A) = 𝟙 _ := by
      have h : P.Base ((x : End A) : A ⟶ A) = P.Base (𝟙 A) := hx.1.1
      rwa [P.Base_id] at h
    have hinv : (inv ((x : End A) : A ⟶ A)) ∈ OTri P A := by
      refine ⟨?_, degFr_of_isIso P _⟩
      show P.Base (inv ((x : End A) : A ⟶ A)) = P.Base (𝟙 A)
      have h : P.Base (inv ((x : End A) : A ⟶ A)) ≫ P.Base ((x : End A) : A ⟶ A)
          = P.Base (𝟙 A) := by rw [← P.Base_comp, IsIso.inv_hom_id]
      rwa [hbx, Category.comp_id] at h
    exact ⟨⟨⟨x, hx.1⟩, ⟨inv ((x : End A) : A ⟶ A), hinv⟩,
      Subtype.ext (by show inv ((x : End A) : A ⟶ A) ≫ _ = _; simp),
      Subtype.ext (by show ((x : End A) : A ⟶ A) ≫ _ = _; simp)⟩, rfl⟩
  · rintro ⟨u, rfl⟩
    exact ⟨(u : OTri P A).2, IsUnit.map (OTri P A).subtype u.isUnit⟩

include F in
/-- ★★**(iii)-1** `𝒪^×` は (ii) の関手の**部分関手** ——
`dstarMap` は `𝒪^×` を `𝒪^×` へ写す。

★モノイド準同型は単元を単元へ送るので、上の `otimes_eq_units` から直ちに従う。 -/
theorem dstarMap_otimes_mem (A B : Istr P)
    (ψ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base)
    (u : OTri P B.obj) (hu : ((u : End B.obj)) ∈ OTimes P B.obj) :
    ((dstarMap P F A B ψ u : End A.obj)) ∈ OTimes P A.obj := by
  obtain ⟨w, hw⟩ := (otimes_eq_units P ((u : End B.obj))).mp hu
  have hwu : (w : OTri P B.obj) = u := Subtype.ext hw
  refine (otimes_eq_units P _).mpr ⟨Units.map (dstarMap P F A B ψ) w, ?_⟩
  show ((dstarMap P F A B ψ (w : OTri P B.obj) : OTri P A.obj) : End A.obj)
      = ((dstarMap P F A B ψ u : OTri P A.obj) : End A.obj)
  rw [hwu]

/-! ## ★(iii) の後半 —— `Div` は関手的な準同型 `𝒪^▷(A) → Φ(A)`

原文 (FrdI p.45):
> the notation of §0]. Moreover, the operation “Div(−)” determines a functorial
-/

/-- ★★★**(iii)-3** `Div` は `otriLin` と**両立する** ——
`Div (otriLin ρ γ) = ρ*(Div γ)`。

★**四角形 `ρ ≫ γ = δ ≫ ρ` に `Div` を当てるだけ**だが、
`Div ρ` を両辺から消すのに ★**`Φ` が integral(消約的)であること**が要る
(`Definition 1.1, (i)` の pre-divisorial に含まれる)。 -/
theorem otriLin_Div {A B : C} (hA : IsIsotropic P A) {ρ : A ⟶ B} (hl : IsLinear P ρ)
    (γ : OTri P B) :
    P.Div ((otriLin P F hA hl γ : End A) : A ⟶ A)
      = Φ.map (P.Base ρ) (P.Div ((γ : End B) : B ⟶ B)) := by
  haveI := isCancelAdd_of_isIntegralMonoid
    (M := Φ.val (P.toElem.obj A).base) (P.divisorial _).1.1
  have hd := congrArg P.Div (otriLin_spec P F hA hl γ)
  rw [P.Div_comp, P.Div_comp] at hd
  have hbδ : P.Base ((otriLin P F hA hl γ : End A) : A ⟶ A) = 𝟙 _ := by
    have h : P.Base ((otriLin P F hA hl γ : End A) : A ⟶ A) = P.Base (𝟙 A) :=
      (otriLin P F hA hl γ).2.1
    rwa [P.Base_id] at h
  rw [hbδ, MonoidOn.map_id, show P.degFr ((γ : End B) : B ⟶ B) = 1 from γ.2.2,
    show P.degFr ρ = 1 from hl] at hd
  simp only [PNat.one_coe, one_smul] at hd
  refine add_left_cancel (a := P.Div ρ) ?_
  rw [← hd]
  exact add_comm _ _

include F in
/-- ★★★**(iii)-3** `𝒟*` の版 —— `Div` は `dstarMap` と両立する。

★**これが原文の「`Div(−)` determines a functorial homomorphism `𝒪^▷(A) → Φ(A)`」**
の中身である。★span を 1 本取って `otriLin_Div` を 2 回使い、
`Φ.map (Base σ)` の単射性(`Definition 1.1, (ii), (a)`)で締める。 -/
theorem dstarMap_Div (A B : Istr P)
    (ψ : (P.toElem.obj A.obj).base ⟶ (P.toElem.obj B.obj).base) (γ : OTri P B.obj) :
    P.Div ((dstarMap P F A B ψ γ : End A.obj) : A.obj ⟶ A.obj)
      = Φ.map ψ (P.Div ((γ : End B.obj) : B.obj ⟶ B.obj)) := by
  obtain ⟨Z, hZ, σ, hσ, ρ, hρ, hb⟩ := dstar_span P F A B ψ
  have key : otriLin P F hZ hσ.1 (dstarMap P F A B ψ γ) = otriLin P F hZ hρ γ :=
    DFunLike.congr_fun (dstarMap_spec P F A B ψ Z hZ σ hσ ρ hρ hb) γ
  have h1 := otriLin_Div P F hZ hσ.1 (dstarMap P F A B ψ γ)
  have h2 := otriLin_Div P F hZ hρ γ
  refine Φ.map_injective (P.Base σ) ?_
  rw [← h1, key, h2, hb, Φ.map_comp]

/-! ## ★★★`Proposition 2.2` 全体の `.src`

★★**9 主張すべての実装名**(数え落とし防止のため表にする):

| 条 | 主張 | 実装 |
|---|---|---|
| (i) | `𝒟* → 𝒟` は圏同値 | `prop_2_2_i` |
| (ii) | 反変関手 `𝒪^▷(−) : 𝒟* → Mon` の**存在** | `otriStar`(＋`dstarMap` / `dstar_map_exists_unique`) |
| (ii) | その**一意性** | `dstarMap_uniq_of_linear`(＋`dstarMap_uniq`) |
| (ii)(a) | linear 射での値が `Proposition 1.11, (iv)` の埋め込み | `dstarMap_eq_otriLin`(＋`otriLin` / `otriLin_injective`) |
| (ii)(b) | pre-step での値が `Definition 1.3, (iii), (c)` の全単射 | `dstarMap_bijective_of_preStep`(＋`otriLin_bijective_of_preStep` / `otriLin_otriFwd`) |
| (iii) | `𝒪^×` は (ii) の関手の**部分関手** | `dstarMap_otimes_mem`(＋`otriLin_otimes_mem`) |
| (iii) | それが `𝒪^▷(A)^±` に**等しい** | `otimes_eq_units` |
| (iii) | `Div` が**関手的な準同型** `𝒪^▷(A) → Φ(A)`、その核が `𝒪^×` | `otri_Div_mul` / `otri_Div_one` / `otriLin_Div` / `dstarMap_Div`、核は `otri_isUnit_iff_Div_zero` / `otri_div_eq_iff` |
| (iv) | isotropic hull が誘導する `𝒪^▷` の埋め込み | `hullOTriHom` / `hullOTriHom_injective` / `hullOTriHom_mem` |

★**`(𝒞^istr)^lin` への制限**(原文の「the restriction of this functor on D* to (Cistr)lin」)は
`otriFunctorLin`(`Prop22.lean`)であり、`dstarMap_eq_otriLin` がその一致を与える。

★★**追加の仮定は 1 つも足していない。** `otriStar` / `dstarMap_comp` が
`G : Frobenioid P` を取るのは、原文が `Proposition 2.2` を Frobenioid について
述べているためであって、原文より強い仮定ではない。 -/
def prop_2_2.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 45, item := "Proposition 2.2",
    sectionId := "frdi-prop-2-2" }

end ABC3.Found.FrdI
