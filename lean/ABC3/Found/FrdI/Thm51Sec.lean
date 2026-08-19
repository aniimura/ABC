/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm51Hom
import ABC3.Found.FrdI.Def45
import ABC3.Found.FrdI.Prop33v

/-!
# [FrdI] Theorem 5.1, (iv) —— unit-trivial 型は model 型

原文 (FrdI p.99):
> trivial type, it follows immediately [cf., e.g., Proposition 3.3, (iii), (iv)] that given

原文 (FrdI p.99):
> identity endomorphism of Frobenius type of A) is uniquely determined by its projec-

原文 (FrdI p.100):
> Frobenius-section. Moreover, since C is of unit-trivial type, it follows immediately from the structure of an elementary Frobenioid [cf. the description of the kernel in Proposition 4.4, (iii)] that C is of birationally Frobenius-normalized type, hence also of model type, as desired. This completes the proof of assertion (iv).

★3 つの部品からなる:

1. **base-section の構成** —— `𝒞` の同型類の**代表元**のうち Frobenius-trivial なもの。
   ★原文は「`(𝒞^Fr-tr)^pl-bk` の skeletal な部分圏」と言うが、
   `Quotient (isIsomorphicSetoid 𝒞)` の代表元を取れば同じことである。
   ★充満性は `Theorem 5.1, (iii)`(`exists_pullBack_frobTrivial`)、
   忠実性は `Proposition 1.11, (ii)`(unit-trivial なら pull-back 射は底で決まる)。
2. **Frobenius-section の構成** —— unit-trivial 型では base-identity な Frobenius 型
   自己射が次数で一意に決まるので、自然性が**自動**になる。
3. **birationally Frobenius-normalized 性** —— `𝒞 → 𝔽_Φ` の忠実性
   (`Proposition 3.3, (ii)` ＋ unit-trivial)を `𝒞^birat` へ持ち上げる。
   ★`Hom^birat` の 2 元は共通の代表元に揃うので、`𝒞` の忠実性がそのまま効く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. 忠実性 -/

/-- ★★**忠実性** —— isotropic 型 ＋ unit-trivial 型なら `𝒞 → 𝔽_Φ` は忠実。 -/
theorem toElem_faithful_of_isotropic (Fc : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) (hut : IsOfUnitTrivialType P) : P.toElem.Faithful := by
  refine ⟨fun {A B} {f g} h => ?_⟩
  obtain ⟨Cc, γ, β, δ, hδ, h₁, h₂⟩ :=
    prop_3_3_ii_sufficiency P Fc hiso
      (congrArg ElemFrobCat.Hom.deg h) (congrArg ElemFrobCat.Hom.div h)
      (congrArg ElemFrobCat.Hom.base h)
  have hd1 : (δ : Cc ⟶ Cc) = 𝟙 Cc := by
    have hb : δ ∈ (⊥ : Submonoid (End Cc)) := by rw [← hut Cc]; exact hδ
    simpa using hb
  rw [h₁, h₂, hd1]
  simp

/-- ★★★**`𝒞^birat → 𝔽_{Φ^gp}` の忠実性**(isotropic ＋ unit-trivial の下)。

★`Hom^birat` の 2 元は共通の代表元 `(a, φ)`、`(a, ψ)` に揃うので、
底・次数・零因子の一致が `𝒞` の側に降り、`𝒞` の忠実性で `φ = ψ` になる。 -/
theorem birat_faithful_of_unitTrivial (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) (hut : IsOfUnitTrivialType P)
    {A B : C} (f g : HomBirat P G A B)
    (hb : biratBase f = biratBase g) (hd : biratDeg f = biratDeg g)
    (hz : biratDivGp f = biratDivGp g) : f = g := by
  haveI := toElem_faithful_of_isotropic G.core hiso hut
  obtain ⟨W, φ, ψ, hφ, hψ⟩ := HomBirat.exists_rep_pair f g
  subst hφ; subst hψ
  haveI hba : IsIso (P.Base W.unop.hom.hom) := W.unop.hom.property.2.2
  rw [biratBase_mk, biratBase_mk, sliceBaseOf_eq, sliceBaseOf_eq] at hb
  rw [biratDeg_mk, biratDeg_mk] at hd
  rw [biratDivGp_mk, biratDivGp_mk, sliceDivGpOf_eq, sliceDivGpOf_eq, hd] at hz
  have hbase : P.Base φ = P.Base ψ := (cancel_epi (inv (P.Base W.unop.hom.hom))).mp hb
  have hzz := congrArg (gpMap _ (Φ.map (P.Base W.unop.hom.hom))) hz
  rw [gpMap_phi_inv_left, gpMap_phi_inv_left] at hzz
  have hdiv : P.Div φ = P.Div ψ :=
    (P.divisorial _).1.1 (sub_left_injective hzz)
  have hpsi : φ = ψ := P.toElem.map_injective (ElemFrobCat.Hom.ext hbase hdiv hd)
  rw [hpsi]

/-! ## ★2. birationally Frobenius-normalized 性 -/

/-- ★`𝒪^▷` の元のべきの `Div^gp` は次数倍。 -/
theorem biratDivGp_pow_otri (G : Frobenioid P) {X : BiratCat P G} {α : End X}
    (hα : α ∈ OTri (biratPre P G) X) (n : ℕ) :
    biratDivGp (((α ^ n : End X) : X ⟶ X)) = n • biratDivGp ((α : X ⟶ X)) := by
  induction n with
  | zero => simpa using biratDivGp_id X
  | succ k ih =>
    have hmem : (α ^ k : End X) ∈ OTri (biratPre P G) X := Submonoid.pow_mem _ hα k
    have hcomp : ((α ^ (k + 1) : End X) : X ⟶ X)
        = ((α : X ⟶ X) ≫ ((α ^ k : End X) : X ⟶ X)) := by
      rw [pow_succ]; rfl
    rw [hcomp, biratDivGp_comp', gpMap_biratBase_of_baseIdentity hα.1,
      show ((biratDeg ((α ^ k : End X) : X ⟶ X) : ℕ+) : ℕ) = 1 from by
        rw [show biratDeg ((α ^ k : End X) : X ⟶ X) = 1 from hmem.2]; rfl,
      one_smul, ih, succ_nsmul]

/-- ★★★★**unit-trivial 型なら `𝒞^birat` は Frobenius-normalized**。 -/
theorem birat_frobNormalized_of_unitTrivial (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) (hut : IsOfUnitTrivialType P) (X : BiratCat P G) :
    IsFrobeniusNormalized (biratPre P G) X := by
  intro φ hφb α hα
  set d : ℕ := (((biratPre P G).degFr ((φ : End X) : X ⟶ X) : ℕ+) : ℕ) with hddef
  have hβmem : (α ^ d : End X) ∈ OTri (biratPre P G) X := Submonoid.pow_mem _ hα d
  have hbφ : (biratPre P G).Base ((φ : End X) : X ⟶ X) = 𝟙 _ :=
    hφb.trans ((biratPre P G).Base_id X)
  have hbα : (biratPre P G).Base ((α : End X) : X ⟶ X) = 𝟙 _ :=
    hα.1.trans ((biratPre P G).Base_id X)
  have hbβ : (biratPre P G).Base (((α ^ d : End X)) : X ⟶ X) = 𝟙 _ :=
    hβmem.1.trans ((biratPre P G).Base_id X)
  refine birat_faithful_of_unitTrivial G hiso hut _ _ ?_ ?_ ?_
  · show (biratPre P G).Base _ = (biratPre P G).Base _
    rw [(biratPre P G).Base_comp, (biratPre P G).Base_comp, hbφ, hbα, hbβ, Category.id_comp]
  · show (biratPre P G).degFr _ = (biratPre P G).degFr _
    rw [(biratPre P G).degFr_comp, (biratPre P G).degFr_comp,
      show (biratPre P G).degFr ((α : End X) : X ⟶ X) = 1 from hα.2,
      show (biratPre P G).degFr (((α ^ d : End X)) : X ⟶ X) = 1 from hβmem.2,
      one_mul, mul_one]
  · rw [biratDivGp_comp', biratDivGp_comp', gpMap_biratBase_of_baseIdentity hφb,
      gpMap_biratBase_of_baseIdentity hα.1,
      show ((biratDeg (((α ^ d : End X)) : X ⟶ X) : ℕ+) : ℕ) = 1 from by
        rw [show biratDeg (((α ^ d : End X)) : X ⟶ X) = 1 from hβmem.2]; rfl,
      one_smul, biratDivGp_pow_otri G hα d,
      show ((biratDeg ((φ : End X) : X ⟶ X) : ℕ+) : ℕ) = d from rfl]
    abel

theorem isOfBirationallyFrobeniusNormalizedType_of_unitTrivial (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) (hut : IsOfUnitTrivialType P) :
    IsOfBirationallyFrobeniusNormalizedType C P G :=
  fun A => birat_frobNormalized_of_unitTrivial G hiso hut ((toBiratCat P G).obj A)

/-! ## ★3. base-section の構成 -/

/-- ★`𝒫` の対象 —— Frobenius-trivial かつ同型類の代表元。 -/
def FrTrObjP (P : PreFrobenioid C Φ) (A : C) : Prop :=
  IsFrobeniusTrivial P A ∧ (Quotient.mk (isIsomorphicSetoid C) A).out = A

/-- ★★★★**[FrdI] Theorem 5.1, (iv)** —— Frobenius-trivial 対象の骨格が定める **base-section**。

原文 (FrdI p.99):
> [in D] lifts to a pull-back morphism of CFr-tr. Thus, we conclude that the natural
-/
noncomputable def frTrBaseSection (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hut : IsOfUnitTrivialType P) : BaseSection P where
  objP A := FrTrObjP P A
  homP {A B} f := FrTrObjP P A ∧ FrTrObjP P B ∧ IsPullBack P f
  id_mem h := ⟨h, h, isPullBack_of_isIso P _⟩
  comp_mem hf hg := ⟨hf.1, hg.2.1, IsPullBack.comp P hf.2.2 hg.2.2⟩
  isPullBack h := h.2.2
  skeletal := by
    rintro A B f g ⟨hA, hB, -⟩ - h1 h2
    haveI : IsIso f := ⟨g, h1, h2⟩
    have hq : (Quotient.mk (isIsomorphicSetoid C) A) = Quotient.mk (isIsomorphicSetoid C) B :=
      Quotient.sound ⟨asIso f⟩
    rw [← hA.2, ← hB.2, hq]
  frobTrivial h := h.1
  essSurjP := by
    intro X
    obtain ⟨A₀, hA₀, ⟨e⟩⟩ := G.core.baseSurj X
    obtain ⟨ι⟩ := Quotient.mk_out (s := isIsomorphicSetoid C) A₀
    haveI hbι : IsIso (P.Base ι.hom) := base_isIso_of_iso ι
    refine ⟨(Quotient.mk (isIsomorphicSetoid C) A₀).out, ⟨?_, ?_⟩, ⟨asIso (P.Base ι.hom) ≪≫ e⟩⟩
    · exact isFrobeniusTrivial_of_iso P G.core ι.symm hA₀
    · exact congrArg Quotient.out (Quotient.out_eq _)
  fullP := by
    rintro A B hA hB ψ
    obtain ⟨f, hpb, hb⟩ := exists_pullBack_frobTrivial G hiso hA.1 hB.1 ψ
    exact ⟨f, ⟨hA, hB, hpb⟩, hb⟩
  faithfulP := by
    rintro A B f g ⟨-, -, hf⟩ ⟨-, -, hg⟩ hbase
    exact prop_1_11_ii P (hut A) f g hf hg hbase

/-! ## ★4. Frobenius-section の構成 -/

/-- ★★★**unit-trivial 型では、base-identity な Frobenius 型自己射は次数で決まる**。 -/
theorem frobTrivial_endo_unique (G : Frobenioid P) (hut : IsOfUnitTrivialType P)
    {A : C} (φ ψ : A ⟶ A) (hφb : IsBaseIdentity P φ) (hψb : IsBaseIdentity P ψ)
    (hφ : IsFrobeniusType P φ) (hψ : IsFrobeniusType P ψ)
    (hd : P.degFr φ = P.degFr ψ) : φ = ψ := by
  obtain ⟨β, hβiso, hβ⟩ := G.core.frobDegUniq A A A φ ψ hφ hψ hd
  haveI := hβiso
  have hbβ : P.Base β = P.Base (𝟙 A) := by
    have h := congrArg P.Base hβ
    rw [P.Base_comp, hφb, hψb, P.Base_id, Category.id_comp] at h
    rw [h, P.Base_id]
  have hmem : (β : End A) ∈ OTimes P A :=
    ⟨⟨hbβ, degFr_of_isIso P β⟩, (isUnit_iff_isIso (β : End A)).mpr hβiso⟩
  rw [hut A] at hmem
  have hb1 : β = 𝟙 A := hmem
  rw [← hβ, hb1, Category.comp_id]

/-- ★base-section の各対象に付いてくる `ζ`。 -/
noncomputable def BaseSection.zeta (S : BaseSection P) (A : S.Obj) : ℕ+ →* End A.1 :=
  (S.frobTrivial A.2).choose

theorem BaseSection.zeta_degFr (S : BaseSection P) (A : S.Obj) (n : ℕ+) :
    P.degFr ((S.zeta A n : End A.1) : A.1 ⟶ A.1) = n :=
  (S.frobTrivial A.2).choose_spec.1 n

theorem BaseSection.zeta_baseId (S : BaseSection P) (A : S.Obj) (n : ℕ+) :
    IsBaseIdentity P ((S.zeta A n : End A.1) : A.1 ⟶ A.1) :=
  ((S.frobTrivial A.2).choose_spec.2 n).1

theorem BaseSection.zeta_frobType (S : BaseSection P) (A : S.Obj) (n : ℕ+) :
    IsFrobeniusType P ((S.zeta A n : End A.1) : A.1 ⟶ A.1) :=
  ((S.frobTrivial A.2).choose_spec.2 n).2

/-- ★★★**`ζ` は `𝒫` の射に沿って自然** —— unit-trivial 型での一意性から。 -/
theorem BaseSection.zeta_naturality (S : BaseSection P) (G : Frobenioid P)
    (hiso : ∀ Y : C, IsIsotropic P Y) (hut : IsOfUnitTrivialType P)
    {A B : S.Obj} (f : A ⟶ B) (n : ℕ+) :
    f.1 ≫ ((S.zeta B n : End B.1) : B.1 ⟶ B.1)
      = ((S.zeta A n : End A.1) : A.1 ⟶ A.1) ≫ f.1 := by
  have hpb : IsPullBack P f.1 := S.isPullBack f.2
  have hdf : P.degFr f.1 = 1 := (G.core.pullBackLB f.1 hpb).2
  have hzf : P.Div f.1 = 0 := (G.core.pullBackLB f.1 hpb).1.2
  obtain ⟨h, hh1, hh2⟩ := pullBack_lift hpb (f.1 ≫ ((S.zeta B n : End B.1) : B.1 ⟶ B.1))
    (P.Base (𝟙 A.1)) (by
      rw [P.Base_comp, S.zeta_baseId B n, P.Base_id, Category.comp_id, P.Base_id,
        Category.id_comp])
  have hdeg : P.degFr h = n := by
    have h1 := congrArg P.degFr hh1
    rw [P.degFr_comp, P.degFr_comp, hdf, one_mul, mul_one, S.zeta_degFr B n] at h1
    exact h1
  have hdiv : P.Div h = 0 := by
    have h1 := congrArg P.Div hh1
    rw [P.Div_comp, P.Div_comp, hzf, map_zero, smul_zero, add_zero, zero_add,
      show P.Div ((S.zeta B n : End B.1) : B.1 ⟶ B.1) = 0 from (S.zeta_frobType B n).1.2,
      map_zero, hdf] at h1
    rw [← h1]; show _ = ((1 : ℕ+) : ℕ) • P.Div h
    rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  have hft : IsFrobeniusType P h :=
    ⟨⟨prop_1_4_i P _ (fun Y _ => hiso Y), hdiv⟩,
      by show IsIso (P.Base h); rw [hh2, P.Base_id]; infer_instance⟩
  have hkey : h = ((S.zeta A n : End A.1) : A.1 ⟶ A.1) :=
    frobTrivial_endo_unique G hut h _ hh2 (S.zeta_baseId A n) hft (S.zeta_frobType A n)
      (by rw [hdeg, S.zeta_degFr A n])
  rw [← hh1, hkey]

/-- ★★★★**[FrdI] Theorem 5.1, (iv)** —— どの base-section も **Frobenius-section** を持つ。 -/
noncomputable def BaseSection.frobSection (S : BaseSection P) (G : Frobenioid P)
    (hiso : ∀ Y : C, IsIsotropic P Y) (hut : IsOfUnitTrivialType P) :
    ℕ+ →* SectionEnd S where
  toFun n := ⟨fun A => S.zeta A n, fun f => S.zeta_naturality G hiso hut f n⟩
  map_one' := by
    refine SectionEnd.ext ?_
    funext A
    exact map_one (S.zeta A)
  map_mul' m n := by
    refine SectionEnd.ext ?_
    funext A
    exact map_mul (S.zeta A) m n

theorem BaseSection.isFrobeniusSection [IsConnected D] (S : BaseSection P) (G : Frobenioid P)
    (hiso : ∀ Y : C, IsIsotropic P Y) (hut : IsOfUnitTrivialType P) :
    IsFrobeniusSection S (S.frobSection G hiso hut) where
  degSection n := by
    haveI := S.isConnected_obj
    rw [SectionEnd.deg_eq _ (Classical.arbitrary S.Obj)]
    exact S.zeta_degFr _ n
  baseIdentity n A := S.zeta_baseId A n
  frobType n A := S.zeta_frobType A n

/-! ## ★5. ★★★★`Theorem 5.1, (iv)` -/

/-- ★★★★**[FrdI] Theorem 5.1, (iv)** —— unit-trivial 型の Frobenioid は **model 型**。

原文 (FrdI p.97):
> associated Frobenius-section F . Moreover, C is of model type.
-/
theorem thm_5_1_iv [IsConnected D] (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hut : IsOfUnitTrivialType P) : IsOfModelType C P G :=
  ⟨⟨{ sec := frTrBaseSection G hiso hut
      frob := (frTrBaseSection G hiso hut).frobSection G hiso hut
      isFrobSection := (frTrBaseSection G hiso hut).isFrobeniusSection G hiso hut }⟩,
   isOfBirationallyFrobeniusNormalizedType_of_unitTrivial G hiso hut⟩

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Theorem 5.1, (iv)`。 -/
def thm_5_1_iv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 97, item := "Theorem 5.1, (iv)",
    sectionId := "frdi-thm-5-1" }

end ABC3.Found.FrdI
