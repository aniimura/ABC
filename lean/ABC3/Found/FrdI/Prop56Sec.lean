/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Ref
import ABC3.Found.FrdI.Prop25iii

/-!
# [FrdI] Proposition 5.6 —— Frobenius-trivial 対象の base-section

原文 (FrdI p.105):
> that C is of model [hence, in particular, isotropic — cf. Definition 2.7, (iii)] and

★命題の主張は「`(σ, φ)` は `(𝒫, F)` の取り方によらず、`𝒪^×(A)` による**共役を除いて**
一意」である。原文の証明は 2 段に分かれる:

1. **claim(還元)**: 「素数 `p` について `u · φ_p · u⁻¹ = φ'_p` なる `u ∈ 𝒪^×(A)` が
   あれば十分」—— これは `σ` の側が**自動的に**同じ `u` で移ることを言う。
2. **claim の検証**: `𝒞` が **unit-profinite 型**であることを使い、
   `M_p ⊆ 𝒪^×(A)` の逆極限で `u` を作る。

★★**このファイルは 1 を実装する**(2 は pro-`p` 群論が要るので別の葉)。

## ★1 の中身

`v_α := σ'_α · (u σ_α u⁻¹)⁻¹ ∈ 𝒪^×(A)` と置く(底が同じ同型どうしの差なので単元)。
`F`, `F'` の自然性から `σ_α · φ_p = φ_p · σ_α`、`σ'_α · φ'_p = φ'_p · σ'_α`。
★**Frobenius-normalized**(`φ_p · x = x^p · φ_p`)と
★**`𝒪^▷(A)` の可換性**(`Proposition 2.5, (iii)`)で両辺を整理すると

  `v_α · Z = v_α^p · Z`  (`Z := u · φ_p · σ_α · u⁻¹`)

となり、★**全射性**(`𝒞` は totally epimorphic)で `v_α = v_α^p`。
`p = 2` を取れば `v_α = 1`。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★0. `End A` の言葉での道具 -/

/-- ★`𝒪^▷(A)` の可換性を `End A` の言葉で。 -/
theorem otri_comm_end {A : C} (hfn : IsFrobeniusNormalized P A) {x y : End A}
    (hx : x ∈ OTri P A) (hy : y ∈ OTri P A) : x * y = y * x :=
  congrArg Subtype.val (otri_mul_comm P hfn ⟨x, hx⟩ ⟨y, hy⟩)

/-- ★Frobenius-normalized の `End A` 版 —— `φ * α = α ^ degFr φ * φ`。 -/
theorem frobNorm_end {A : C} (hfn : IsFrobeniusNormalized P A) {φ α : End A}
    (hφ : IsBaseIdentity P φ) (hα : α ∈ OTri P A) :
    φ * α = α ^ ((P.degFr ((φ : A ⟶ A)) : ℕ+) : ℕ) * φ :=
  (hfn φ hφ α hα).symm

/-- ★底が恒等な同型は `𝒪^×(A)` の元。 -/
theorem otimes_of_isIso_baseId {A : C} {v : End A} (hv : IsIso ((v : A ⟶ A)))
    (hb : P.Base ((v : A ⟶ A)) = 𝟙 _) : v ∈ OTimes P A := by
  refine ⟨⟨?_, degFr_of_isIso P _⟩, (CategoryTheory.isUnit_iff_isIso v).mpr hv⟩
  show P.Base ((v : A ⟶ A)) = P.Base (𝟙 A)
  rw [hb, P.Base_id]

/-- ★`v ∈ 𝒪^×(A)` で `v^2 = v` なら `v = 1`。 -/
theorem otimes_eq_one_of_sq {A : C} {v : End A} (hv : v ∈ OTimes P A) (h : v ^ 2 = v) :
    v = 1 := by
  have h2 : v * v = v * 1 := by rw [mul_one, ← sq]; exact h
  exact hv.2.mul_left_cancel h2

/-! ## ★1. 主計算 —— `v = v^p` -/

/-- ★★★★**[FrdI] Proposition 5.6 の主計算** ——
`σ'_α = v · (u σ_α u⁻¹)` と `φ' = u φ u⁻¹` が交換関係を満たすなら `v^p = v`。

原文 (FrdI p.106):
> — which [by the total epimorphicity of C] implies that vα = vp
-/
theorem prop_5_6_v_pow (htot : IsTotallyEpimorphic C) {A : C}
    (hfn : IsFrobeniusNormalized P A)
    {u ui v φ s : End A}
    (huT : u ∈ OTri P A) (huiT : ui ∈ OTri P A) (_hui : u * ui = 1) (hiu : ui * u = 1)
    (hvT : v ∈ OTri P A)
    (hφ : IsBaseIdentity P φ)
    (hcomm : s * φ = φ * s)
    (hcomm' : (v * (u * s * ui)) * (u * φ * ui) = (u * φ * ui) * (v * (u * s * ui))) :
    v ^ ((P.degFr ((φ : A ⟶ A)) : ℕ+) : ℕ) = v := by
  set p := ((P.degFr ((φ : A ⟶ A)) : ℕ+) : ℕ) with hp
  have hcuv : Commute u v := otri_comm_end hfn huT hvT
  have hcuiv : ui * v = v * ui := otri_comm_end hfn huiT hvT
  have hL : (v * (u * s * ui)) * (u * φ * ui) = v * (u * (φ * (s * ui))) := by
    calc (v * (u * s * ui)) * (u * φ * ui)
        = v * u * s * (ui * u) * φ * ui := by simp only [mul_assoc]
      _ = v * u * (s * φ) * ui := by rw [hiu]; simp only [mul_assoc, one_mul]
      _ = v * (u * (φ * (s * ui))) := by rw [hcomm]; simp only [mul_assoc]
  have hR : (u * φ * ui) * (v * (u * s * ui)) = u * ((φ * v) * (s * ui)) := by
    calc (u * φ * ui) * (v * (u * s * ui))
        = u * φ * (ui * v * u * (s * ui)) := by simp only [mul_assoc]
      _ = u * ((φ * v) * (s * ui)) := by
          rw [hcuiv, mul_assoc v ui u, hiu, mul_one]; simp only [mul_assoc]
  rw [hL, hR, frobNorm_end hfn hφ hvT, ← hp] at hcomm'
  have hcomm2 : u * v ^ p = v ^ p * u := (hcuv.pow_right p).eq
  have h2 : v * (u * (φ * (s * ui))) = v ^ p * (u * (φ * (s * ui))) := by
    rw [hcomm']
    calc u * (v ^ p * φ * (s * ui)) = (u * v ^ p) * (φ * (s * ui)) := by simp only [mul_assoc]
      _ = (v ^ p * u) * (φ * (s * ui)) := by rw [hcomm2]
      _ = v ^ p * (u * (φ * (s * ui))) := by simp only [mul_assoc]
  haveI : Epi ((u * (φ * (s * ui)) : End A) : A ⟶ A) := htot A A _
  exact ((cancel_epi ((u * (φ * (s * ui)) : End A) : A ⟶ A)).mp h2).symm

/-! ## ★2. `σ` —— base-section が与える `Aut_𝒟(A_𝒟) → Aut_𝒞(A)` -/

/-- ★★**[FrdI] Proposition 5.6 の `σ`**。 -/
noncomputable def BaseSection.sigma (S : BaseSection P) {A : C} (hA : S.objP A)
    (α : (P.toElem.obj A).base ⟶ (P.toElem.obj A).base) : End A := S.lift hA hA α

theorem BaseSection.sigma_base (S : BaseSection P) {A : C} (hA : S.objP A)
    (α : (P.toElem.obj A).base ⟶ (P.toElem.obj A).base) :
    P.Base ((S.sigma hA α : A ⟶ A)) = α := S.lift_base hA hA α

theorem BaseSection.sigma_mul (S : BaseSection P) {A : C} (hA : S.objP A)
    (α β : (P.toElem.obj A).base ⟶ (P.toElem.obj A).base) :
    S.sigma hA β * S.sigma hA α = S.sigma hA (α ≫ β) := by
  show S.lift hA hA α ≫ S.lift hA hA β = S.lift hA hA (α ≫ β)
  rw [S.lift_comp hA hA hA]

theorem BaseSection.sigma_id (S : BaseSection P) {A : C} (hA : S.objP A) :
    S.sigma hA (𝟙 _) = 1 := by
  show S.lift hA hA (𝟙 _) = 𝟙 A
  rw [← P.Base_id A, S.lift_id]

/-- ★`𝒫` の自己射は `End(𝒫 → 𝒞)` の成分と可換(自然性)。 -/
theorem SectionEnd.comm_of_homP {S : BaseSection P} (ε : SectionEnd S) {A : C} (hA : S.objP A)
    {f : End A} (hf : S.homP f) :
    f * ε.app ⟨A, hA⟩ = ε.app ⟨A, hA⟩ * f :=
  (ε.naturality (A := (⟨A, hA⟩ : S.Obj)) (B := (⟨A, hA⟩ : S.Obj)) ⟨f, hf⟩).symm

/-! ## ★3. claim(還元) -/

variable [IsConnected D]

/-- ★★★★★**[FrdI] Proposition 5.6 の「claim」** ——
`F` を共役で移す `u ∈ 𝒪^×(A)` が(次数 2 のところで)1 つあれば、
`σ` も**同じ** `u` の共役で移る。

原文 (FrdI p.106):
> Then I claim that it suffices to show the existence of a
-/
theorem prop_5_6_sigma_conj (htot : IsTotallyEpimorphic C)
    {A : C} (hfn : IsFrobeniusNormalized P A)
    {S S' : BaseSection P} {Fs : ℕ+ →* SectionEnd S} {Fs' : ℕ+ →* SectionEnd S'}
    (hFs : IsFrobeniusSection S Fs)
    (hA : S.objP A) (hA' : S'.objP A)
    {u ui : End A} (huT : u ∈ OTri P A) (huiT : ui ∈ OTri P A)
    (hui : u * ui = 1) (hiu : ui * u = 1)
    (hconj : u * ((Fs 2).app ⟨A, hA⟩) * ui = (Fs' 2).app ⟨A, hA'⟩)
    (α : (P.toElem.obj A).base ⟶ (P.toElem.obj A).base) [IsIso α] :
    u * S.sigma hA α * ui = S'.sigma hA' α := by
  set s := S.sigma hA α with hs
  set s' := S'.sigma hA' α with hs'
  set sinv := S.sigma hA (inv α) with hsinv
  set s'inv := S'.sigma hA' (inv α) with hs'inv
  have hss1 : s * sinv = 1 := by
    rw [hs, hsinv, S.sigma_mul hA, IsIso.inv_hom_id, S.sigma_id]
  have hss2 : sinv * s = 1 := by
    rw [hs, hsinv, S.sigma_mul hA, IsIso.hom_inv_id, S.sigma_id]
  have hs's1 : s' * s'inv = 1 := by
    rw [hs', hs'inv, S'.sigma_mul hA', IsIso.inv_hom_id, S'.sigma_id]
  have hs's2 : s'inv * s' = 1 := by
    rw [hs', hs'inv, S'.sigma_mul hA', IsIso.hom_inv_id, S'.sigma_id]
  set w := u * s * ui with hw
  set winv := u * sinv * ui with hwinv
  have hww1 : w * winv = 1 := by
    calc w * winv = u * s * (ui * u) * sinv * ui := by rw [hw, hwinv]; simp only [mul_assoc]
      _ = u * (s * sinv) * ui := by rw [hiu]; simp only [mul_assoc, one_mul]
      _ = 1 := by rw [hss1, mul_one, hui]
  have hww2 : winv * w = 1 := by
    calc winv * w = u * sinv * (ui * u) * s * ui := by rw [hw, hwinv]; simp only [mul_assoc]
      _ = u * (sinv * s) * ui := by rw [hiu]; simp only [mul_assoc, one_mul]
      _ = 1 := by rw [hss2, mul_one, hui]
  set v := s' * winv with hv
  have hvw : v * w = s' := by rw [hv, mul_assoc, hww2, mul_one]
  have hvunit : IsUnit v := by
    refine ⟨⟨v, w * s'inv, ?_, ?_⟩, rfl⟩
    · calc v * (w * s'inv) = (v * w) * s'inv := (mul_assoc v w s'inv).symm
        _ = 1 := by rw [hvw, hs's1]
    · calc (w * s'inv) * v = w * (s'inv * s') * winv := by rw [hv]; simp only [mul_assoc]
        _ = 1 := by rw [hs's2, mul_one, hww1]
  haveI hviso : IsIso ((v : A ⟶ A)) := (CategoryTheory.isUnit_iff_isIso v).mp hvunit
  have hbu : P.Base ((u : A ⟶ A)) = 𝟙 _ := by
    have := huT.1; rw [show P.Base ((u : A ⟶ A)) = P.Base (𝟙 A) from this, P.Base_id]
  have hbui : P.Base ((ui : A ⟶ A)) = 𝟙 _ := by
    have := huiT.1; rw [show P.Base ((ui : A ⟶ A)) = P.Base (𝟙 A) from this, P.Base_id]
  have hbwinv : P.Base ((winv : A ⟶ A)) = inv α := by
    show P.Base ((ui : A ⟶ A) ≫ (sinv : A ⟶ A) ≫ (u : A ⟶ A)) = inv α
    rw [P.Base_comp, P.Base_comp, hbu, hbui, S.sigma_base hA, Category.comp_id,
      Category.id_comp]
  have hbv : P.Base ((v : A ⟶ A)) = 𝟙 _ := by
    show P.Base ((winv : A ⟶ A) ≫ (s' : A ⟶ A)) = 𝟙 _
    rw [P.Base_comp, hbwinv, S'.sigma_base hA', IsIso.inv_hom_id]
  have hvT : v ∈ OTimes P A := otimes_of_isIso_baseId hviso hbv
  have hnat : s * ((Fs 2).app ⟨A, hA⟩) = ((Fs 2).app ⟨A, hA⟩) * s :=
    (Fs 2).comm_of_homP hA (S.lift_homP hA hA α)
  have hnat' : s' * ((Fs' 2).app ⟨A, hA'⟩) = ((Fs' 2).app ⟨A, hA'⟩) * s' :=
    (Fs' 2).comm_of_homP hA' (S'.lift_homP hA' hA' α)
  have hcomm' : (v * (u * s * ui)) * (u * ((Fs 2).app ⟨A, hA⟩) * ui)
      = (u * ((Fs 2).app ⟨A, hA⟩) * ui) * (v * (u * s * ui)) := by
    rw [← hw, hvw, hconj]; exact hnat'
  have hdeg : P.degFr (((Fs 2).app ⟨A, hA⟩ : A ⟶ A)) = 2 :=
    ((Fs 2).deg_eq ⟨A, hA⟩).symm.trans (hFs.degSection 2)
  have hkey := prop_5_6_v_pow htot hfn huT huiT hui hiu hvT.1
    (hFs.baseIdentity 2 ⟨A, hA⟩) hnat hcomm'
  rw [hdeg] at hkey
  have hv1 : v = 1 := otimes_eq_one_of_sq hvT (by exact hkey)
  rw [← hvw, hv1, one_mul]

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 5.6` の「claim」(還元)の条。
★★命題全体ではない —— `u` の存在(unit-profinite 型を使う pro-`p` の議論)は未実装。 -/
def prop_5_6_sigma_conj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.6 — claim(σ は φ と同じ u の共役で移る)",
    sectionId := "frdi-prop-5-6" }

end ABC3.Found.FrdI
