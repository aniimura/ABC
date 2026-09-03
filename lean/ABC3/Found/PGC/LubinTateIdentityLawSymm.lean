import ABC3.Found.PGC.LubinTateIdentityLaw

/-!
# Lubin-Tate 形式群法則: 単位元則 `F_f(0,Y)=Y`(`Found/PGC/LubinTateIdentityLaw.lean`
の `X_0↔X_1` 対称版)

`F_f(X,0)=X`(`formalGroupLaw_identity`)と全く同じ構造の議論を
`X_0` と `X_1` を入れ替えて繰り返し、`F_f(0,Y)=Y` を確立する——
`emb`(`X_0` 経由の埋め込み)の対 `emb1`(`X_1` 経由)、`restrictR`
(`X_1↦0`)の対 `restrictL`(`X_0↦0`)を新たに用意し、`subst_family_
comp_value`・`hasSubst_const`・`subst_zero_eq_zero`(`Found/PGC/
LubinTateIdentityLaw.lean` の一般補題、`X_0`/`X_1` によらない)を
そのまま再利用する。 -/

namespace ABC3.Found.PGC

variable {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))

/-! ### `X_1` を経由する埋め込み `emb1`(`Found/PGC/LubinTateUniqueness.lean::emb` の対) -/

noncomputable def emb1 (p : PowerSeries A) : MvPowerSeries (Fin 2) A :=
  PowerSeries.subst (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) p

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
theorem hasSubst_X1 : PowerSeries.HasSubst (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) := by
  show IsNilpotent (MvPowerSeries.constantCoeff (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A))
  rw [MvPowerSeries.constantCoeff_X]; exact IsNilpotent.zero

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
theorem coeff_emb1 (p : PowerSeries A) (e : Fin 2 →₀ ℕ) :
    MvPowerSeries.coeff e (emb1 p) =
      if e 0 = 0 then PowerSeries.coeff (e 1) p else 0 := by
  rw [emb1, PowerSeries.coeff_subst hasSubst_X1]
  rw [finsum_eq_sum_of_support_subset _ (s := {e 1}) (fun d hd => by
    simp only [Function.mem_support] at hd
    simp only [Finset.coe_singleton, Set.mem_singleton_iff]
    by_contra hcon
    apply hd
    have : MvPowerSeries.coeff e ((MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) ^ d) = 0 := by
      rw [MvPowerSeries.coeff_X_pow]
      split_ifs with hcase
      · exfalso; apply hcon; rw [hcase, Finsupp.single_eq_same]
      · rfl
    simp [this])]
  by_cases h0 : e 0 = 0
  · rw [if_pos h0]
    have he : e = Finsupp.single (1 : Fin 2) (e 1) := by
      ext i; fin_cases i
      · simp [h0]
      · simp
    rw [Finset.sum_singleton, he, MvPowerSeries.coeff_X_pow]
    simp
  · rw [if_neg h0]
    rw [Finset.sum_singleton]
    have : MvPowerSeries.coeff e ((MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) ^ (e 1)) = 0 := by
      rw [MvPowerSeries.coeff_X_pow]
      split_ifs with hcase
      · exfalso; apply h0; rw [hcase, Finsupp.single_eq_of_ne' (by decide : (1 : Fin 2) ≠ 0)]
      · rfl
    simp [this]

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
theorem emb1_injective : Function.Injective (emb1 (A := A)) := by
  intro p q hpq
  ext n
  have hp := coeff_emb1 (A := A) p (Finsupp.single (1 : Fin 2) n)
  have hq := coeff_emb1 (A := A) q (Finsupp.single (1 : Fin 2) n)
  simp only [Finsupp.single_eq_same] at hp hq
  rw [hpq] at hp
  rw [hp] at hq
  simpa using hq

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
theorem subst_emb1 (p g : PowerSeries A) (hp : PowerSeries.HasSubst p) :
    PowerSeries.subst (emb1 p) g = emb1 (PowerSeries.subst p g) :=
  (PowerSeries.subst_comp_subst_apply hp hasSubst_X1 g).symm

/-! ### `X_0` を `0` に潰す代入 `restrictL`(`restrictR` の対) -/

noncomputable def restrictL : Fin 2 → MvPowerSeries (Fin 2) A :=
  fun i => if i = 1 then (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) else 0

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
theorem hasSubst_restrictL : MvPowerSeries.HasSubst (restrictL (A := A)) := by
  constructor
  · intro i
    show IsNilpotent (MvPowerSeries.constantCoeff (restrictL i))
    unfold restrictL
    split_ifs with h
    · rw [MvPowerSeries.constantCoeff_X]; exact IsNilpotent.zero
    · rw [map_zero]; exact IsNilpotent.zero
  · intro d; exact Set.toFinite _

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
theorem coeff_restrictL_eq_zero_of_ne (Φ : MvPowerSeries (Fin 2) A) (e : Fin 2 →₀ ℕ)
    (he : e 0 ≠ 0) :
    MvPowerSeries.coeff e (MvPowerSeries.subst (restrictL (A := A)) Φ) = 0 := by
  rw [MvPowerSeries.coeff_subst (hasSubst_restrictL (A := A))]
  rw [finsum_eq_sum_of_support_subset _ (s := (∅ : Finset (Fin 2 →₀ ℕ))) (fun d hd => by
    exfalso
    simp only [Function.mem_support] at hd
    apply hd
    by_cases hd0 : d 0 = 0
    · have hprod : d.prod (fun s m => ((restrictL (A := A)) s) ^ m) =
          (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) ^ (d 1) := by
        show (∏ i ∈ d.support, ((restrictL (A := A)) i) ^ (d i)) = _
        rw [Finset.prod_subset (Finset.subset_univ d.support) (fun x _ hx => by
          simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two, hd0]
        simp [restrictL]
      rw [hprod, MvPowerSeries.coeff_X_pow]
      have hne : e ≠ Finsupp.single (1 : Fin 2) (d 1) := by
        intro heq; apply he; rw [heq, Finsupp.single_eq_of_ne' (by decide : (1 : Fin 2) ≠ 0)]
      simp [hne]
    · have hprod0 : d.prod (fun s m => ((restrictL (A := A)) s) ^ m) = 0 := by
        show (∏ i ∈ d.support, ((restrictL (A := A)) i) ^ (d i)) = 0
        apply Finset.prod_eq_zero (i := (0 : Fin 2))
        · simp only [Finsupp.mem_support_iff]; exact hd0
        · simp [restrictL, zero_pow hd0]
      simp [hprod0])]
  simp

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
theorem coeff_restrictL_eq_of_e0_zero (Φ : MvPowerSeries (Fin 2) A) (e : Fin 2 →₀ ℕ)
    (he : e 0 = 0) :
    MvPowerSeries.coeff e (MvPowerSeries.subst (restrictL (A := A)) Φ) =
      MvPowerSeries.coeff e Φ := by
  rw [MvPowerSeries.coeff_subst (hasSubst_restrictL (A := A))]
  rw [finsum_eq_sum_of_support_subset _ (s := ({Finsupp.single (1 : Fin 2) (e 1)} : Finset (Fin 2 →₀ ℕ)))
    (fun d hd => by
      simp only [Function.mem_support] at hd
      simp only [Finset.coe_singleton, Set.mem_singleton_iff]
      by_contra hcon
      apply hd
      by_cases hd0 : d 0 = 0
      · have hprod : d.prod (fun s m => ((restrictL (A := A)) s) ^ m) =
            (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) ^ (d 1) := by
          show (∏ i ∈ d.support, ((restrictL (A := A)) i) ^ (d i)) = _
          rw [Finset.prod_subset (Finset.subset_univ d.support) (fun x _ hx => by
            simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two, hd0]
          simp [restrictL]
        rw [hprod, MvPowerSeries.coeff_X_pow]
        have hd1 : d 1 ≠ e 1 := by
          intro heq1
          apply hcon
          have hd_eq : d = Finsupp.single (1 : Fin 2) (d 1) := by
            ext i; fin_cases i
            · simp [hd0]
            · simp
          rw [hd_eq, heq1]
        have : e ≠ Finsupp.single (1 : Fin 2) (d 1) := by
          intro heq; apply hd1; rw [heq, Finsupp.single_eq_same]
        simp [this]
      · have hprod0 : d.prod (fun s m => ((restrictL (A := A)) s) ^ m) = 0 := by
          show (∏ i ∈ d.support, ((restrictL (A := A)) i) ^ (d i)) = 0
          apply Finset.prod_eq_zero (i := (0 : Fin 2))
          · simp only [Finsupp.mem_support_iff]; exact hd0
          · simp [restrictL, zero_pow hd0]
        simp [hprod0])]
  have he_eq : e = Finsupp.single (1 : Fin 2) (e 1) := by
    ext i; fin_cases i
    · simp [he]
    · simp
  have hprod : (Finsupp.single (1 : Fin 2) (e 1)).prod (fun s m => ((restrictL (A := A)) s) ^ m) =
      (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) ^ (e 1) := by
    show (∏ i ∈ (Finsupp.single (1 : Fin 2) (e 1)).support,
      ((restrictL (A := A)) i) ^ ((Finsupp.single (1 : Fin 2) (e 1)) i)) = _
    rw [Finset.prod_subset (Finset.subset_univ _) (fun x _ hx => by
      simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two]
    simp [restrictL]
  rw [Finset.sum_singleton, hprod, MvPowerSeries.coeff_X_pow]
  rw [show e = Finsupp.single (1 : Fin 2) (e 1) from he_eq] at *
  simp

include hq hπmax hf0 hf1 hfres in
/-- `formalGroupLaw` の次数1・`X_1` の係数も `1`
(`coeff_single01_formalGroupLaw` の対、同じ理由: `φ_0` は次数2)。 -/
theorem coeff_single10_formalGroupLaw :
    MvPowerSeries.coeff (Finsupp.single (1 : Fin 2) 1) (formalGroupLaw hq hπmax f hf0 hf1 hfres) = 1 := by
  have hdeg1 : Finsupp.degree (Finsupp.single (1 : Fin 2) 1) = 1 := by
    rw [finsupp_degree_fin2, Finsupp.single_eq_same,
      Finsupp.single_eq_of_ne' (by decide : (1 : Fin 2) ≠ 0)]
  rw [formalGroupLaw, coeff_LubinTateF_eq_ΦSeq hq hπmax f (hf0' f hf0) hf1 hfres f hf0 hf1 hfres
    (Finsupp.single (1 : Fin 2) 1) 1 (by rw [hdeg1])]
  rw [coeff_ΦSeq_stable hq hπmax f (hf0' f hf0) hf1 hfres f hf0 hf1 hfres
    (Finsupp.single (1 : Fin 2) 1) 0 1 (by omega) (by rw [hdeg1])]
  show MvPowerSeries.coeff (Finsupp.single (1 : Fin 2) 1)
    ((MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) + MvPowerSeries.X 1) = 1
  rw [map_add, MvPowerSeries.coeff_X, MvPowerSeries.coeff_X]
  have h2 : Finsupp.single (1 : Fin 2) 1 ≠ Finsupp.single (0 : Fin 2) 1 := by
    intro heq
    have h01 : (Finsupp.single (1 : Fin 2) 1 : Fin 2 →₀ ℕ) 1 =
        (Finsupp.single (0 : Fin 2) 1 : Fin 2 →₀ ℕ) 1 := by rw [heq]
    rw [Finsupp.single_eq_same, Finsupp.single_eq_of_ne' (by decide : (0 : Fin 2) ≠ 1)] at h01
    exact one_ne_zero h01
  simp [h2]

include hq hπmax hf0 hf1 hfres in
/-- `F_f` を `restrictL` で制限した「1変数の中身」——`ψ_L(Y):=F_f(0,Y)`。 -/
noncomputable def psiL : PowerSeries A :=
  PowerSeries.mk (fun n => MvPowerSeries.coeff (Finsupp.single (1 : Fin 2) n)
    (MvPowerSeries.subst (restrictL (A := A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres)))

include hq hπmax hf0 hf1 hfres in
theorem emb1_psiL : emb1 (psiL hq hπmax f hf0 hf1 hfres) =
    MvPowerSeries.subst (restrictL (A := A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres) := by
  apply MvPowerSeries.ext
  intro e
  rw [coeff_emb1]
  by_cases h0 : e 0 = 0
  · rw [if_pos h0, psiL, PowerSeries.coeff_mk]
    have he : e = Finsupp.single (1 : Fin 2) (e 1) := by
      ext i; fin_cases i
      · simp [h0]
      · simp
    rw [← he]
  · rw [if_neg h0]
    exact (coeff_restrictL_eq_zero_of_ne (formalGroupLaw hq hπmax f hf0 hf1 hfres) e h0).symm

include hq hπmax hf0 hf1 hfres in
theorem constantCoeff_psiL : PowerSeries.constantCoeff (psiL hq hπmax f hf0 hf1 hfres) = 0 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply, psiL, PowerSeries.coeff_mk]
  rw [coeff_restrictL_eq_of_e0_zero _ _ (by simp)]
  rw [show Finsupp.single (1 : Fin 2) 0 = (0 : Fin 2 →₀ ℕ) from Finsupp.single_zero _]
  rw [MvPowerSeries.coeff_zero_eq_constantCoeff_apply]
  exact constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres

include hq hπmax hf0 hf1 hfres in
theorem coeff_one_psiL : PowerSeries.coeff 1 (psiL hq hπmax f hf0 hf1 hfres) = 1 := by
  rw [psiL, PowerSeries.coeff_mk]
  rw [coeff_restrictL_eq_of_e0_zero _ _ (by simp)]
  exact coeff_single10_formalGroupLaw hq hπmax f hf0 hf1 hfres

include hq hπmax hf0 hf1 hfres in
/-- ★★★★★★★★★**`ψ_L:=F_f(0,Y)` は `f` との関数等式を満たす**:
`ψ_L(f(Y))=f(ψ_L(Y))`——`Found/PGC/LubinTateIdentityLaw.lean::
psi_functional_equation` と全く同じ議論を `X_0↔X_1` を入れ替えて行う。 -/
theorem psiL_functional_equation :
    PowerSeries.subst f (psiL hq hπmax f hf0 hf1 hfres) =
      PowerSeries.subst (psiL hq hπmax f hf0 hf1 hfres) f := by
  apply emb1_injective
  have hf0c : PowerSeries.constantCoeff f = 0 := hf0' f hf0
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0c]; exact IsNilpotent.zero
  have hHSψ : PowerSeries.HasSubst (psiL hq hπmax f hf0 hf1 hfres) := by
    show IsNilpotent (PowerSeries.constantCoeff (psiL hq hπmax f hf0 hf1 hfres))
    rw [constantCoeff_psiL]; exact IsNilpotent.zero
  set a : Fin 2 → MvPowerSeries (Fin 2) A :=
    fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f with ha_def
  have hHSa : MvPowerSeries.HasSubst a := hasSubst_g_subst_X f hf0c
  have hembf0 : MvPowerSeries.constantCoeff (emb1 f) = 0 := by
    rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply, coeff_emb1]; simp [hf0]
  have hHSembf : MvPowerSeries.HasSubst (fun _ : Fin 2 => emb1 f) := hasSubst_const hembf0
  have hrestrict_a : ∀ i : Fin 2, MvPowerSeries.subst (restrictL (A := A)) (a i) =
      PowerSeries.subst (restrictL (A := A) i) f := by
    intro i
    rw [ha_def]
    show MvPowerSeries.subst (restrictL (A := A))
      (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) = _
    rw [subst_family_comp_value (hasSubst_restrictL (A := A)) (MvPowerSeries.constantCoeff_X i) f]
    rw [MvPowerSeries.subst_X (hasSubst_restrictL (A := A)) i]
  have has1 : MvPowerSeries.subst (restrictL (A := A)) (a 1) = emb1 f := by
    rw [hrestrict_a]
    show PowerSeries.subst (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) f = emb1 f
    rfl
  have has0 : MvPowerSeries.subst (restrictL (A := A)) (a 0) = 0 := by
    rw [hrestrict_a]
    show PowerSeries.subst (restrictL (A := A) 0) f = 0
    rw [show (restrictL (A := A) 0) = 0 from rfl]
    exact subst_zero_eq_zero hf0
  have hstep1 : MvPowerSeries.subst (restrictL (A := A))
      (MvPowerSeries.subst a (formalGroupLaw hq hπmax f hf0 hf1 hfres)) =
      MvPowerSeries.subst (fun _ : Fin 2 => emb1 f)
        (MvPowerSeries.subst (restrictL (A := A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres)) := by
    rw [MvPowerSeries.subst_comp_subst_apply hHSa (hasSubst_restrictL (A := A))]
    rw [MvPowerSeries.subst_comp_subst_apply (hasSubst_restrictL (A := A)) hHSembf]
    have hpt : ∀ s : Fin 2, MvPowerSeries.subst (restrictL (A := A)) (a s) =
        MvPowerSeries.subst (fun _ : Fin 2 => emb1 f) (restrictL (A := A) s) := by
      intro s
      by_cases hs : s = 1
      · subst hs; rw [has1]; exact (MvPowerSeries.subst_X hHSembf 1).symm
      · have hs0 : s = 0 := by omega
        subst hs0
        rw [has0, show (restrictL (A := A)) 0 = 0 from rfl]
        have hz := MvPowerSeries.subst_sub hHSembf (0 : MvPowerSeries (Fin 2) A) 0
        simpa using hz.symm
    exact congrArg (fun g => MvPowerSeries.subst g (formalGroupLaw hq hπmax f hf0 hf1 hfres))
      (funext hpt)
  have hstep2 : MvPowerSeries.subst (fun _ : Fin 2 => emb1 f) (emb1 (psiL hq hπmax f hf0 hf1 hfres)) =
      emb1 (PowerSeries.subst f (psiL hq hπmax f hf0 hf1 hfres)) := by
    show MvPowerSeries.subst (fun _ : Fin 2 => emb1 f)
      (PowerSeries.subst (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A)
        (psiL hq hπmax f hf0 hf1 hfres)) = _
    rw [subst_family_comp_value hHSembf (MvPowerSeries.constantCoeff_X 1)]
    rw [show MvPowerSeries.subst (fun _ : Fin 2 => emb1 f)
        (MvPowerSeries.X (1 : Fin 2) : MvPowerSeries (Fin 2) A) = emb1 f from
        MvPowerSeries.subst_X hHSembf 1]
    exact subst_emb1 f (psiL hq hπmax f hf0 hf1 hfres) hHSf
  have hstep3 : MvPowerSeries.subst (restrictL (A := A))
      (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hfres) f) =
      emb1 (PowerSeries.subst (psiL hq hπmax f hf0 hf1 hfres) f) := by
    have hFf0 : MvPowerSeries.constantCoeff (formalGroupLaw hq hπmax f hf0 hf1 hfres) = 0 :=
      constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres
    rw [subst_family_comp_value (hasSubst_restrictL (A := A)) hFf0]
    rw [show MvPowerSeries.subst (restrictL (A := A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres) =
        emb1 (psiL hq hπmax f hf0 hf1 hfres) from (emb1_psiL hq hπmax f hf0 hf1 hfres).symm]
    exact subst_emb1 (psiL hq hπmax f hf0 hf1 hfres) f hHSψ
  have hFf := formalGroupLaw_f_isEndomorphism hq hπmax f hf0 hf1 hfres
  calc emb1 (PowerSeries.subst f (psiL hq hπmax f hf0 hf1 hfres))
      = MvPowerSeries.subst (fun _ : Fin 2 => emb1 f) (emb1 (psiL hq hπmax f hf0 hf1 hfres)) :=
        hstep2.symm
    _ = MvPowerSeries.subst (fun _ : Fin 2 => emb1 f)
        (MvPowerSeries.subst (restrictL (A := A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres)) := by
        rw [emb1_psiL]
    _ = MvPowerSeries.subst (restrictL (A := A))
        (MvPowerSeries.subst a (formalGroupLaw hq hπmax f hf0 hf1 hfres)) := hstep1.symm
    _ = MvPowerSeries.subst (restrictL (A := A))
        (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hfres) f) := by rw [hFf]
    _ = emb1 (PowerSeries.subst (psiL hq hπmax f hf0 hf1 hfres) f) := hstep3

include hq hπmax hf0 hf1 hfres in
/-- ★★★★★★★★★★**Lubin-Tate 形式群法則の単位元則(対称版)**: `F_f(0,Y)=Y`。
`Found/PGC/LubinTateIdentityLaw.lean::formalGroupLaw_identity` と同じ議論。 -/
theorem formalGroupLaw_identity_left (hπne0 : π ≠ 0) :
    psiL hq hπmax f hf0 hf1 hfres = PowerSeries.X := by
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0' f hf0]; exact IsNilpotent.zero
  apply powerSeries_uniqueness hπmax hπne0 (hf0' f hf0) hf1
    (constantCoeff_psiL hq hπmax f hf0 hf1 hfres) PowerSeries.constantCoeff_X
  · rw [coeff_one_psiL, PowerSeries.coeff_one_X]
  · exact psiL_functional_equation hq hπmax f hf0 hf1 hfres
  · rw [PowerSeries.X_subst, PowerSeries.subst_X hHSf]

end ABC3.Found.PGC
