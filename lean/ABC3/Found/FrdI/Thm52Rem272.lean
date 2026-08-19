/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Ref

/-!
# [FrdI] Remark 2.7.2 の 3 分解 —— 等長な場合は skeleton の仮定が要らない

原文 (FrdI p.52):
> is a base-identity pre-step endomorphism;

`Def27.lean` の `rem_2_7_2_factor` は **`𝒞` が base-trivial 型かつ skeleton**
を仮定する。原文も「WLOG `𝒞` は skeleton としてよい」と書く。

★★**しかし `Theorem 5.2, (iv)` が実際に使うのは `𝒞^Fr-tr` と `𝒞^birat` であり、
そこでは「すべての射が等長」**である。そしてその場合、
**skeleton も base-trivial 型も要らない**:

`φ : A → B`(`A, B ∈ Ob(𝒫)`)に対し
1. `α := 𝒫` の一意な持ち上げ(底 `Base φ`)。これは pull-back 射。
2. その普遍性で `g ≫ α = φ`、`Base g = 𝟙` なる `g : A → A` が一意に取れる。
3. ★**`g` は等長なので Frobenius 型**(isotropic 型なので co-angular は自動)。
4. `Definition 1.3, (ii)` の一意性で `F_n ≫ β = g` なる同型 `β` が取れる。

★一意性(`rem_2_7_2_uniq`)は**もともと skeleton を使っていない**ので、そのまま使える。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} [IsConnected D] {S : BaseSection P}

/-- ★★★★**`Remark 2.7.2` の 3 分解** —— **すべての射が等長**なら
`skeleton` と `base-trivial 型`の仮定なしで出る。

原文 (FrdI p.52):
> where α is P-distinguished; β is a base-identity pre-step endomorphism; γ is F distinguished.
-/
theorem rem_2_7_2_factor_isometric (Fc : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hmet : ∀ {A B : C} (f : A ⟶ B), IsIsometric P f)
    {Fs : ℕ+ →* SectionEnd S} (hFs : IsFrobeniusSection S Fs)
    {A B : C} (hA : S.objP A) (hB : S.objP B) (φ : A ⟶ B) :
    ∃ β : A ⟶ A, IsBaseIdentity P β ∧ IsPreStep P β ∧
      φ = (((Fs (P.degFr φ)).app ⟨A, hA⟩ : End A) : A ⟶ A) ≫ β
            ≫ S.lift hA hB (P.Base φ) := by
  have hαP := S.lift_homP hA hB (P.Base φ)
  have hαpb := S.isPullBack hαP
  have hdegα : P.degFr (S.lift hA hB (P.Base φ)) = 1 := (Fc.pullBackLB _ hαpb).2
  obtain ⟨g, hg1, hg2⟩ := pullBack_lift hαpb φ (P.Base (𝟙 A)) (by
    rw [P.Base_id, Category.id_comp, S.lift_base])
  have hgb : P.Base g = 𝟙 _ := by rw [hg2, P.Base_id]
  have hdeg : P.degFr g = P.degFr φ := by
    have h := congrArg P.degFr hg1
    rw [P.degFr_comp, hdegα, one_mul] at h
    exact h
  have hgF : IsFrobeniusType P g :=
    ⟨⟨prop_1_4_i P g (fun Y _ => hiso Y), hmet g⟩, by
      show IsIso (P.Base g); rw [hgb]; infer_instance⟩
  have hFn : IsFrobeniusType P (((Fs (P.degFr φ)).app ⟨A, hA⟩ : End A) : A ⟶ A) :=
    hFs.frobType _ _
  have hFnb : P.Base (((Fs (P.degFr φ)).app ⟨A, hA⟩ : End A) : A ⟶ A) = 𝟙 _ := by
    have hh := hFs.baseIdentity (P.degFr φ) ⟨A, hA⟩
    show P.Base _ = _
    rw [show P.Base (((Fs (P.degFr φ)).app ⟨A, hA⟩ : End A) : A ⟶ A) = P.Base (𝟙 A) from hh,
      P.Base_id]
  have hdegFn : P.degFr (((Fs (P.degFr φ)).app ⟨A, hA⟩ : End A) : A ⟶ A) = P.degFr φ :=
    (SectionEnd.deg_eq (Fs (P.degFr φ)) ⟨A, hA⟩).symm.trans (hFs.degSection _)
  obtain ⟨β, hβiso, hβ⟩ := Fc.frobDegUniq A A A _ g hFn hgF (hdegFn.trans hdeg.symm)
  haveI := hβiso
  have hbβ : P.Base β = 𝟙 _ := by
    have h := congrArg P.Base hβ
    rw [P.Base_comp, hFnb, Category.id_comp, hgb] at h
    exact h
  refine ⟨β, ?_, ?_, ?_⟩
  · show P.Base β = P.Base (𝟙 A)
    rw [P.Base_id]; exact hbβ
  · exact ⟨degFr_of_isIso P β, by show IsIso (P.Base β); rw [hbβ]; infer_instance⟩
  · rw [← Category.assoc, hβ, hg1]

/-- ★★★**3 分解の `β` は `𝒪^×(A)` の元**(等長な場合)。

★`β` は base-identity な pre-step で、等長なので isotropic 型のもとで**同型**、
したがって `𝒪^×(A)` に入る。★これが `Theorem 5.2, (iv)` の `u` 成分の正体である。 -/
theorem rem_2_7_2_beta_otimes (hiso : ∀ X : C, IsIsotropic P X)
    (hmet : ∀ {A B : C} (f : A ⟶ B), IsIsometric P f)
    {A : C} (β : A ⟶ A) (hb : IsBaseIdentity P β) (hs : IsPreStep P β) :
    (β : End A) ∈ OTimes P A := by
  haveI : IsIso β := hiso A A β (hmet β) hs
  exact ⟨⟨hb, hs.1⟩, (isUnit_iff_isIso (β : End A)).mpr inferInstance⟩

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Remark 2.7.2` の 3 分解(等長版)。 -/
def rem_2_7_2_factor_isometric.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 52,
    item := "Remark 2.7.2 — 3 分解の存在(等長な場合、skeleton 不要)",
    sectionId := "frdi-def-2-7" }

/-! ## ★2. `𝒞^birat` の base-section

★`Remark 2.7.2`(等長版)を `𝒞^birat` に当てるには `𝒞^birat` の base-section が要る。
`𝒫` の像で取れる —— 唯一の非自明な点は
**「`𝒞` の pull-back 射の像は `𝒞^birat` の pull-back 射」**である。 -/

/-- ★`𝒞` の射との合成は代表元の合成。 -/
theorem compBirat_mk_toHomBirat {G : Frobenioid P} (F : FrobenioidCore P) {A B E : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (π : B ⟶ E) :
    compBirat P G F (HomBirat.mk Z φ) (toHomBirat (P := P) (G := G) π)
      = HomBirat.mk Z (φ ≫ π) := by
  rw [toHomBirat, compBirat_mk]
  have hsq2 : biratPullGamma F Z φ (idxBiratOne P G B) ≫ φ
      = biratPullAlpha F Z φ (idxBiratOne P G B) :=
    (biratPull_sq F Z φ (idxBiratOne P G B)).trans (Category.comp_id _)
  rw [← hsq2]
  exact (congrArg (HomBirat.mk (biratPullIdx F Z φ (idxBiratOne P G B)))
      (Category.assoc _ _ _)).trans
    (HomBirat.mk_map (biratPullHom F Z φ (idxBiratOne P G B)) (φ ≫ π))

/-- ★★★**`𝒞` の pull-back 射の像は `𝒞^birat` の pull-back 射**。 -/
theorem birat_isPullBack_of (G : Frobenioid P) {A B : C} (π : A ⟶ B) (hπ : IsPullBack P π) :
    IsPullBack (biratPre P G) ((toBiratCat P G).map π) := by
  intro Z
  constructor
  · intro g g' heq
    obtain ⟨W, φ, φ', hφ, hφ'⟩ := HomBirat.exists_rep_pair g g'
    subst hφ; subst hφ'
    haveI hba : IsIso (P.Base W.unop.hom.hom) := W.unop.hom.property.2.2
    have h1 := congrArg (fun t => t.1) heq
    have hb : biratBase (HomBirat.mk W φ) = biratBase (HomBirat.mk W φ') :=
      congrArg Prod.snd h1
    rw [biratBase_mk, biratBase_mk, sliceBaseOf_eq, sliceBaseOf_eq] at hb
    have hbase : P.Base φ = P.Base φ' := (cancel_epi (inv (P.Base W.unop.hom.hom))).mp hb
    have e1 := compBirat_mk_toHomBirat (G := G) G.core W φ π
    have e2 := compBirat_mk_toHomBirat (G := G) G.core W φ' π
    have hc : HomBirat.mk W (φ ≫ π) = HomBirat.mk W (φ' ≫ π) :=
      e1.symm.trans ((congrArg Prod.fst h1).trans e2)
    obtain ⟨V, u, hu⟩ := HomBirat.eq_iff_same.mp hc
    refine HomBirat.eq_iff_same.mpr ⟨V, u, ?_⟩
    refine pullBack_uniq hπ ?_ ?_
    · exact (Category.assoc _ _ _).trans (hu.trans (Category.assoc _ _ _).symm)
    · exact (P.Base_comp _ _).trans
        ((congrArg (fun t => P.Base u.unop.left.hom ≫ t) hbase).trans (P.Base_comp _ _).symm)
  · rintro ⟨⟨h, k⟩, hcomm⟩
    obtain ⟨W, ψ, hψ⟩ := HomBirat.exists_rep h
    subst hψ
    haveI hba : IsIso (P.Base W.unop.hom.hom) := W.unop.hom.property.2.2
    have hb : inv (P.Base W.unop.hom.hom) ≫ P.Base ψ = k ≫ P.Base π := by
      have h2 : biratBase (HomBirat.mk W ψ) = k ≫ P.Base π :=
        hcomm.trans (congrArg (fun t => k ≫ t) (biratBase_toHomBirat (P := P) (G := G) π))
      rw [biratBase_mk, sliceBaseOf_eq] at h2
      exact h2
    have hb2 : P.Base ψ = (P.Base W.unop.hom.hom ≫ k) ≫ P.Base π := by
      have h3 := congrArg (fun t => P.Base W.unop.hom.hom ≫ t) hb
      simp only [← Category.assoc, IsIso.hom_inv_id, Category.id_comp] at h3
      exact h3.trans (Category.assoc _ _ _).symm
    obtain ⟨φ, hφ1, hφ2⟩ := pullBack_lift hπ ψ (P.Base W.unop.hom.hom ≫ k) hb2
    refine ⟨HomBirat.mk W φ, Subtype.ext (Prod.ext ?_ ?_)⟩
    · exact (compBirat_mk_toHomBirat (G := G) G.core W φ π).trans
        (congrArg (HomBirat.mk W) hφ1)
    · show biratBase (HomBirat.mk W φ) = k
      rw [biratBase_mk, sliceBaseOf_eq, hφ2]
      rw [← Category.assoc, IsIso.inv_hom_id, Category.id_comp]

/-- ★★★**`𝒞^birat` の base-section**(`𝒫` の像)。 -/
noncomputable def biratBaseSection (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (S : BaseSection P) : BaseSection (biratPre P G) where
  objP A := S.objP (biratDown P G A)
  homP {A B} f := ∃ f₀ : biratDown P G A ⟶ biratDown P G B,
    S.homP f₀ ∧ f = (toBiratCat P G).map f₀
  id_mem h := ⟨𝟙 _, S.id_mem h, ((toBiratCat P G).map_id _).symm⟩
  comp_mem := by
    rintro A B E f g ⟨f₀, hf, rfl⟩ ⟨g₀, hg, rfl⟩
    exact ⟨f₀ ≫ g₀, S.comp_mem hf hg, ((toBiratCat P G).map_comp _ _).symm⟩
  isPullBack := by
    rintro A B f ⟨f₀, hf, rfl⟩
    exact birat_isPullBack_of G f₀ (S.isPullBack hf)
  skeletal := by
    rintro A B f g ⟨f₀, hf, rfl⟩ ⟨g₀, hg, rfl⟩ h1 h2
    haveI := toBiratCat_faithful (P := P) (G := G)
    refine S.skeletal hf hg ((toBiratCat P G).map_injective ?_) ((toBiratCat P G).map_injective ?_)
    · rw [(toBiratCat P G).map_comp, (toBiratCat P G).map_id]; exact h1
    · rw [(toBiratCat P G).map_comp, (toBiratCat P G).map_id]; exact h2
  frobTrivial := by
    intro A hA
    obtain ⟨ζ, hdeg, hprop⟩ := S.frobTrivial hA
    have hid : biratBase (𝟙 A) = P.Base (𝟙 (biratDown P G A)) := by
      rw [show (𝟙 A) = toHomBirat (𝟙 (biratDown P G A)) from
        ((toBiratCat P G).map_id _).symm, biratBase_toHomBirat]
    refine ⟨{ toFun := fun n =>
        (((toBiratCat P G).map (((ζ n : End (biratDown P G A))) : _ ⟶ _)) : End A),
              map_one' := ?_, map_mul' := ?_ }, ?_, ?_⟩
    · show (toBiratCat P G).map (((ζ 1 : End (biratDown P G A))) : _ ⟶ _) = 𝟙 _
      rw [map_one]
      exact (toBiratCat P G).map_id _
    · intro m n
      show (toBiratCat P G).map (((ζ (m * n) : End (biratDown P G A))) : _ ⟶ _)
        = (toBiratCat P G).map (((ζ n : End (biratDown P G A))) : _ ⟶ _)
          ≫ (toBiratCat P G).map (((ζ m : End (biratDown P G A))) : _ ⟶ _)
      rw [← (toBiratCat P G).map_comp, map_mul]
      rfl
    · intro n
      exact (biratDeg_toHomBirat (G := G) ((ζ n : End (biratDown P G A)) : _ ⟶ _)).trans (hdeg n)
    · intro n
      refine ⟨?_, ⟨⟨prop_1_4_i (biratPre P G) _ (fun Y _ => birat_isOfIsotropicType hiso Y),
        birat_isIsometric _⟩, ?_⟩⟩
      · exact ((biratBase_toHomBirat (G := G) ((ζ n : End (biratDown P G A)) : _ ⟶ _)).trans
          ((hprop n).1)).trans hid.symm
      · show IsIso (biratBase ((toBiratCat P G).map
            (((ζ n : End (biratDown P G A))) : _ ⟶ _)))
        exact (biratBase_toHomBirat (G := G)
          ((ζ n : End (biratDown P G A)) : _ ⟶ _)).symm ▸ (hprop n).2.2
  essSurjP X := by
    obtain ⟨A, hA, ⟨e⟩⟩ := S.essSurjP X
    exact ⟨A, hA, ⟨e⟩⟩
  fullP := by
    intro A B hA hB ψ
    refine ⟨(toBiratCat P G).map (S.lift hA hB ψ), ⟨_, S.lift_homP hA hB ψ, rfl⟩, ?_⟩
    exact (biratBase_toHomBirat (G := G) (S.lift hA hB ψ)).trans (S.lift_base hA hB ψ)
  faithfulP := by
    rintro A B f g ⟨f₀, hf, rfl⟩ ⟨g₀, hg, rfl⟩ hbase
    have h : P.Base f₀ = P.Base g₀ :=
      ((biratBase_toHomBirat (G := G) f₀).symm.trans hbase).trans
        (biratBase_toHomBirat (G := G) g₀)
    rw [S.faithfulP hf hg h]

/-- ★locator —— `𝒞^birat` の base-section(`Theorem 5.2, (iv)` の準備)。 -/
def biratBaseSection.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 101,
    item := "Theorem 5.2, (iv) — 𝒞^birat の base-section",
    sectionId := "frdi-thm-5-2" }

end ABC3.Found.FrdI
