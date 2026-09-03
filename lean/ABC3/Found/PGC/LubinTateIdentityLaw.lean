import ABC3.Found.PGC.LubinTateUniqueness

/-!
# Lubin-Tate 形式群法則: 単位元則 `F_f(X,0)=X`(進行中)

`Found/PGC/LubinTateUniqueness.lean::powerSeries_uniqueness`(一意性補題)を
`formalGroupLaw`(`F_f`)に実際に適用し、単位元則 `F_f(X,0)=X` を示す。
-/

namespace ABC3.Found.PGC

variable {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))

/-- `MvPowerSeries` の中で「`X_1` を `0` に潰す」代入。 -/
noncomputable def restrictR : Fin 2 → MvPowerSeries (Fin 2) A :=
  fun i => if i = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) else 0

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
theorem hasSubst_restrictR : MvPowerSeries.HasSubst (restrictR (A := A)) := by
  constructor
  · intro i
    show IsNilpotent (MvPowerSeries.constantCoeff (restrictR i))
    unfold restrictR
    split_ifs with h
    · rw [MvPowerSeries.constantCoeff_X]; exact IsNilpotent.zero
    · rw [map_zero]; exact IsNilpotent.zero
  · intro d; exact Set.toFinite _

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
/-- `Φ` を `restrictR` で制限した結果は `X_1` に依存しない
(`e 1 ≠ 0` の係数は消える)。 -/
theorem coeff_restrictR_eq_zero_of_ne (Φ : MvPowerSeries (Fin 2) A) (e : Fin 2 →₀ ℕ)
    (he : e 1 ≠ 0) :
    MvPowerSeries.coeff e (MvPowerSeries.subst (restrictR (A := A)) Φ) = 0 := by
  rw [MvPowerSeries.coeff_subst (hasSubst_restrictR (A := A))]
  rw [finsum_eq_sum_of_support_subset _ (s := (∅ : Finset (Fin 2 →₀ ℕ))) (fun d hd => by
    exfalso
    simp only [Function.mem_support] at hd
    apply hd
    by_cases hd1 : d 1 = 0
    · have hprod : d.prod (fun s m => ((restrictR (A := A)) s) ^ m) =
          (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) ^ (d 0) := by
        show (∏ i ∈ d.support, ((restrictR (A := A)) i) ^ (d i)) = _
        rw [Finset.prod_subset (Finset.subset_univ d.support) (fun x _ hx => by
          simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two, hd1]
        simp [restrictR]
      rw [hprod, MvPowerSeries.coeff_X_pow]
      have hne : e ≠ Finsupp.single (0 : Fin 2) (d 0) := by
        intro heq; apply he; rw [heq, Finsupp.single_eq_of_ne' (by decide : (0 : Fin 2) ≠ 1)]
      simp [hne]
    · have hprod0 : d.prod (fun s m => ((restrictR (A := A)) s) ^ m) = 0 := by
        show (∏ i ∈ d.support, ((restrictR (A := A)) i) ^ (d i)) = 0
        apply Finset.prod_eq_zero (i := (1 : Fin 2))
        · simp only [Finsupp.mem_support_iff]; exact hd1
        · simp [restrictR, zero_pow hd1]
      simp [hprod0])]
  simp

end ABC3.Found.PGC
