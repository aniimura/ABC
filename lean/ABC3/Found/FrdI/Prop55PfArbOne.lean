/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfArb

/-!
# [FrdI] Proposition 5.5, (i) を**任意の対象**へ

原文 (FrdI p.104):
> the case of arbitrary A then follows by considering "pairs of pre-steps" as in Theorem

★★`Prop55Pf.lean` の `otriPfEquiv`

  `𝒪^▷(A)^pf ≅ 𝒪^▷(A^pf)`

は **`A` が Frobenius-trivial** のときのものである。原文はそこから
**pre-step の対**で任意の `A` へ渡す。★`Prop55PfArb.lean` で
その道具(`otriPullHom_bijective`、`toPfRoot_isPreStep`)は揃っているので、
本ファイルは**同型そのものを移す**。

## ★移し方

`F.baseSurj` ＋ `F.preStepSpan` で **pre-step の対** `A₀ ⟵ X ⟶ A`
(`A₀` は Frobenius-trivial)を取り、

```
𝒪^▷(A)  ≅  𝒪^▷(X)  ≅  𝒪^▷(A₀)          （𝒞 の側）
𝒪^▷(A^pf) ≅ 𝒪^▷(X^pf) ≅ 𝒪^▷(A₀^pf)     （𝒞^pf の側）
```

を作り、`A₀` での `otriPfEquiv` に挟む。★`Pf` は関手なので `𝒞` の側の同型を
そのまま持ち上げられる(`pfEquivOfEquiv`)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section OTriTransport

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★★**co-angular pre-step に沿って `𝒪^▷` は同型**(`otriPullHom_bijective` の言い換え)。 -/
noncomputable def otriPreStepMulEquiv (F : FrobenioidCore P) {X A : C} (φ : X ⟶ A)
    (hc : IsCoAngular P φ) (hs : IsPreStep P φ) :
    OTri P A ≃* OTri P X :=
  MulEquiv.ofBijective (otriPullHom P F φ hc hs.1) (otriPullHom_bijective P F φ hc hs)

/-- ★★**pre-step の対で `𝒪^▷` を移す**。 -/
noncomputable def otriSpanEquiv (F : FrobenioidCore P) {X A A₀ : C}
    (φ : X ⟶ A₀) (hcφ : IsCoAngular P φ) (hsφ : IsPreStep P φ)
    (ψ : X ⟶ A) (hcψ : IsCoAngular P ψ) (hsψ : IsPreStep P ψ) :
    OTri P A ≃* OTri P A₀ :=
  (otriPreStepMulEquiv P F ψ hcψ hsψ).trans (otriPreStepMulEquiv P F φ hcφ hsφ).symm

end OTriTransport

/-! ## ★2. `Pf` は単系の同型を持ち上げる -/

section PfEquiv

/-- ★★単系の同型は完全化の同型を与える(`gpEquiv` の `Pf` 版)。 -/
noncomputable def pfEquivOfEquiv {M N : Type w} [AddCommMonoid M] [AddCommMonoid N]
    (e : M ≃+ N) : Pf M ≃+ Pf N where
  toFun := Pf.map (e : M →+ N)
  invFun := Pf.map (e.symm : N →+ M)
  left_inv x := by
    induction x using Pf.inductionOn with | _ m a =>
    show Pf.map (e.symm : N →+ M) (Pf.map (e : M →+ N) (Pf.mk m a)) = Pf.mk m a
    rw [Pf.map_mk, Pf.map_mk]
    exact congrArg (fun t => Pf.mk t a) (e.symm_apply_apply m)
  right_inv x := by
    induction x using Pf.inductionOn with | _ n a =>
    show Pf.map (e : M →+ N) (Pf.map (e.symm : N →+ M) (Pf.mk n a)) = Pf.mk n a
    rw [Pf.map_mk, Pf.map_mk]
    exact congrArg (fun t => Pf.mk t a) (e.apply_symm_apply n)
  map_add' x y := map_add _ x y

end PfEquiv

/-! ## ★3. `Proposition 5.5, (i)` を任意の対象へ -/

section ArbOne

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**[FrdI] Proposition 5.5, (i)(任意の対象)** ——
`𝒪^▷(A)^pf ≅ 𝒪^▷(A^pf)`。

原文 (FrdI p.104):
> the case of arbitrary A then follows by considering "pairs of pre-steps" as in Theorem

★`F.baseSurj` で Frobenius-trivial な `A₀`(底は `A` と同型)を取り、
`F.preStepSpan` で pre-step の対 `A₀ ⟵ X ⟶ A` を得て、
`𝒞` の側と `𝒞^pf` の側の両方で `𝒪^▷` を移す。
★`𝒞` / `𝒞^pf` はどちらも isotropic なので pre-step は自動的に co-angular である。 -/
noncomputable def otriPfEquiv_of_span (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfn : ∀ X : C, IsFrobeniusNormalized P X)
    (G : Frobenioid P) (Gpf : Frobenioid (pfRootPre P F))
    {A A₀ X : C} (hA₀ : IsFrobeniusTrivial P A₀)
    (φ : X ⟶ A₀) (hsφ : IsPreStep P φ) (ψ : X ⟶ A) (hsψ : IsPreStep P ψ) :
    @Pf (Additive (OTri P A)) (otriAddCommMonoid (hfn A))
      ≃+ Additive (OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) :=
  letI := otriAddCommMonoid (hfn A)
  letI := otriAddCommMonoid (hfn A₀)
  -- `𝒞` の側は isotropic なので pre-step は co-angular
  haveI hcφ : IsCoAngular P φ := isCoAngular_of_isotropic P G (hiso X) φ hsφ
  haveI hcψ : IsCoAngular P ψ := isCoAngular_of_isotropic P G (hiso X) ψ hsψ
  -- `𝒞^pf` の側も同様
  haveI hsφ' : IsPreStep (pfRootPre P F) ((toPfRoot P F).map φ) := toPfRoot_isPreStep φ hsφ
  haveI hsψ' : IsPreStep (pfRootPre P F) ((toPfRoot P F).map ψ) := toPfRoot_isPreStep ψ hsψ
  haveI hcφ' : IsCoAngular (pfRootPre P F) ((toPfRoot P F).map φ) :=
    isCoAngular_of_isotropic (pfRootPre P F) Gpf
      (pfRoot_isOfIsotropicType (F := F) hfi _) _ hsφ'
  haveI hcψ' : IsCoAngular (pfRootPre P F) ((toPfRoot P F).map ψ) :=
    isCoAngular_of_isotropic (pfRootPre P F) Gpf
      (pfRoot_isOfIsotropicType (F := F) hfi _) _ hsψ'
  -- `A₀` での同型
  (pfEquivOfEquiv (M := Additive (OTri P A)) (N := Additive (OTri P A₀))
      (MulEquiv.toAdditive (otriSpanEquiv P F φ hcφ hsφ ψ hcψ hsψ))).trans
    ((otriPfEquiv (F := F) hiso hA₀ (hfn A₀) (hfn (rtObj P F A₀ 1))
        hA₀.choose hA₀.choose_spec.1 hA₀.choose_spec.2).trans
      (MulEquiv.toAdditive
        (otriSpanEquiv (pfRootPre P F) Gpf.core _ hcφ' hsφ' _ hcψ' hsψ')).symm)

/-- ★★★★★★★**[FrdI] Proposition 5.5, (i)(任意の対象)** ——
`𝒪^▷(A)^pf ≅ 𝒪^▷(A^pf)`。

★`F.baseSurj` ＋ `F.preStepSpan` で pre-step の対を取り、
`otriPfEquiv_of_span` に渡すだけ。 -/
theorem otriPfEquiv_of_arb (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfn : ∀ X : C, IsFrobeniusNormalized P X)
    (G : Frobenioid P) (Gpf : Frobenioid (pfRootPre P F)) (A : C) :
    Nonempty (@Pf (Additive (OTri P A)) (otriAddCommMonoid (hfn A))
      ≃+ Additive (OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F))) := by
  obtain ⟨A₀, hA₀, ⟨e⟩⟩ := F.baseSurj ((P.toElem.obj A).base)
  obtain ⟨X, φ, ψ, hsφ, hsψ, -⟩ := F.preStepSpan A₀ A e.hom (by infer_instance)
  exact ⟨otriPfEquiv_of_span hfi hiso hfn G Gpf hA₀ φ hsφ ψ hsψ⟩

end ArbOne

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Proposition 5.5, (i)` を任意の対象へ渡した段。 -/
def otriPfEquiv_of_arb.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — 𝒪^▷(A)^pf ≅ 𝒪^▷(A^pf)(任意の対象)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
