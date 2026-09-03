import ABC3.Found.PGC.LubinTateUniqueness

/-!
# Lubin-Tate 形式群法則: 単位元則 `F_f(X,0)=X`(進行中)

`Found/PGC/LubinTateUniqueness.lean::powerSeries_uniqueness`(一意性補題)を
`formalGroupLaw`(`F_f`)に実際に適用し、単位元則 `F_f(X,0)=X` を示す計画。

## ここまでで確立したもの(`sorry` 無し)

- `restrictR`・`hasSubst_restrictR`・`coeff_restrictR_eq_zero_of_ne`・
  `coeff_restrictR_eq_of_e1_zero`: `X_1↦0` の代入とその係数計算。
- `psi`(`ψ:=F_f(X,0)` の1変数の中身)・`emb_psi`(`emb ψ = F_f.subst
  (restrictR)`)・`constantCoeff_psi`・`coeff_one_psi`(定数項0・次数1の
  係数が1)。
- `constantCoeff_formalGroupLaw`・`coeff_single01_formalGroupLaw`:
  `F_f` 自身の対応する係数(`ΦSeq` まで遡って計算)。
- `subst_zero_eq_zero`: `PowerSeries.subst 0 p = 0`(`p` の定数項が0のとき)。

## 残っている段(未着手、具体的な計画あり)

`ψ` が `f` との関数等式 `ψ.subst(f)=f.subst(ψ)` を満たすことを示す段
(その後 `powerSeries_uniqueness` で `ψ=X` を結論する)。手計算では以下まで
詰めた: `formalGroupLaw_f_isEndomorphism` の両辺に `restrictR`(`Y↦0`)を
合成し、`MvPowerSeries.subst_comp_subst_apply`・`PowerSeries.subst_
comp_subst_apply` で丁寧に追跡すると、最終的に
「`c := fun s=>if s=0 then emb f else 0` を `F_f` に代入した結果」を
「`emb f` を `ψ` に代入した結果の埋め込み」に結び付ける一般補題
(`emb ψ' が X_1 に依存しないとき、任意の族 c による MvPowerSeries.subst
は c 0 による PowerSeries.subst と一致する`)が必要だと分かった。この
一般補題自体は `MvPowerSeries.coeff_subst` の finsum を `d 1=0` へ圧縮
するだけで**数学的には**明らかだが、有限台への相当な reindexing
(`Fin2→₀ℕ` 上の finsum を `ℕ` 上の finsum に単射 `n↦single 0 n` 経由で
戻す)を要る具体的な Lean 化がこのセッションでは詰め切れず、`sorry` を
含む形で残すよりも本ファイルから一旦除いた——次に戻るならここが
出発点(`finsum_comp_injective` 相当の道具を探すか、`Finset.sum` へ
早めに落として直接計算するか)。 -/

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

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
/-- `e 1 = 0` のとき、`restrictR` による制限は `Φ` 自身の同じ係数をそのまま返す。 -/
theorem coeff_restrictR_eq_of_e1_zero (Φ : MvPowerSeries (Fin 2) A) (e : Fin 2 →₀ ℕ)
    (he : e 1 = 0) :
    MvPowerSeries.coeff e (MvPowerSeries.subst (restrictR (A := A)) Φ) =
      MvPowerSeries.coeff e Φ := by
  rw [MvPowerSeries.coeff_subst (hasSubst_restrictR (A := A))]
  rw [finsum_eq_sum_of_support_subset _ (s := ({Finsupp.single (0 : Fin 2) (e 0)} : Finset (Fin 2 →₀ ℕ)))
    (fun d hd => by
      simp only [Function.mem_support] at hd
      simp only [Finset.coe_singleton, Set.mem_singleton_iff]
      by_contra hcon
      apply hd
      by_cases hd1 : d 1 = 0
      · have hprod : d.prod (fun s m => ((restrictR (A := A)) s) ^ m) =
            (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) ^ (d 0) := by
          show (∏ i ∈ d.support, ((restrictR (A := A)) i) ^ (d i)) = _
          rw [Finset.prod_subset (Finset.subset_univ d.support) (fun x _ hx => by
            simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two, hd1]
          simp [restrictR]
        rw [hprod, MvPowerSeries.coeff_X_pow]
        have hd0 : d 0 ≠ e 0 := by
          intro heq0
          apply hcon
          have hd_eq : d = Finsupp.single (0 : Fin 2) (d 0) := by
            ext i; fin_cases i
            · simp
            · simp [hd1]
          rw [hd_eq, heq0]
        have : e ≠ Finsupp.single (0 : Fin 2) (d 0) := by
          intro heq; apply hd0; rw [heq, Finsupp.single_eq_same]
        simp [this]
      · have hprod0 : d.prod (fun s m => ((restrictR (A := A)) s) ^ m) = 0 := by
          show (∏ i ∈ d.support, ((restrictR (A := A)) i) ^ (d i)) = 0
          apply Finset.prod_eq_zero (i := (1 : Fin 2))
          · simp only [Finsupp.mem_support_iff]; exact hd1
          · simp [restrictR, zero_pow hd1]
        simp [hprod0])]
  have he_eq : e = Finsupp.single (0 : Fin 2) (e 0) := by
    ext i; fin_cases i
    · simp
    · simp [he]
  have hprod : (Finsupp.single (0 : Fin 2) (e 0)).prod (fun s m => ((restrictR (A := A)) s) ^ m) =
      (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) ^ (e 0) := by
    show (∏ i ∈ (Finsupp.single (0 : Fin 2) (e 0)).support,
      ((restrictR (A := A)) i) ^ ((Finsupp.single (0 : Fin 2) (e 0)) i)) = _
    rw [Finset.prod_subset (Finset.subset_univ _) (fun x _ hx => by
      simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two]
    simp [restrictR]
  rw [Finset.sum_singleton, hprod, MvPowerSeries.coeff_X_pow]
  rw [show e = Finsupp.single (0 : Fin 2) (e 0) from he_eq] at *
  simp

include hq hπmax hf0 hf1 hfres in
/-- `formalGroupLaw` の定数項は0——`ΦSeq 0`(出発点 `X_0+X_1`)まで遡るだけ。 -/
theorem constantCoeff_formalGroupLaw :
    MvPowerSeries.constantCoeff (formalGroupLaw hq hπmax f hf0 hf1 hfres) = 0 := by
  rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply]
  rw [formalGroupLaw, coeff_LubinTateF_eq_ΦSeq hq hπmax f (hf0' f hf0) hf1 hfres f hf0 hf1 hfres
    (0 : Fin 2 →₀ ℕ) 0 (le_refl 0)]
  show MvPowerSeries.coeff (0 : Fin 2 →₀ ℕ)
    ((MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) + MvPowerSeries.X 1) = 0
  simp

include hq hπmax hf0 hf1 hfres in
/-- `formalGroupLaw` の次数1・`X_0` の係数は `1`——`ΦSeq 1` まで遡ると、
次数1の斉次部分は出発点 `X_0+X_1` からまだ変化していない
(`φ_0` は次数2の斉次式なので、次数1には効かない)。 -/
theorem coeff_single01_formalGroupLaw :
    MvPowerSeries.coeff (Finsupp.single (0 : Fin 2) 1) (formalGroupLaw hq hπmax f hf0 hf1 hfres) = 1 := by
  have hdeg1 : Finsupp.degree (Finsupp.single (0 : Fin 2) 1) = 1 := by
    rw [finsupp_degree_fin2, Finsupp.single_eq_same,
      Finsupp.single_eq_of_ne' (by decide : (0 : Fin 2) ≠ 1)]
    ring
  rw [formalGroupLaw, coeff_LubinTateF_eq_ΦSeq hq hπmax f (hf0' f hf0) hf1 hfres f hf0 hf1 hfres
    (Finsupp.single (0 : Fin 2) 1) 1 (by rw [hdeg1])]
  rw [coeff_ΦSeq_stable hq hπmax f (hf0' f hf0) hf1 hfres f hf0 hf1 hfres
    (Finsupp.single (0 : Fin 2) 1) 0 1 (by omega) (by rw [hdeg1])]
  show MvPowerSeries.coeff (Finsupp.single (0 : Fin 2) 1)
    ((MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) + MvPowerSeries.X 1) = 1
  rw [map_add, MvPowerSeries.coeff_X, MvPowerSeries.coeff_X]
  have h2 : Finsupp.single (0 : Fin 2) 1 ≠ Finsupp.single (1 : Fin 2) 1 := by
    intro heq
    have h01 : (Finsupp.single (0 : Fin 2) 1 : Fin 2 →₀ ℕ) 0 =
        (Finsupp.single (1 : Fin 2) 1 : Fin 2 →₀ ℕ) 0 := by rw [heq]
    rw [Finsupp.single_eq_same, Finsupp.single_eq_of_ne' (by decide : (1 : Fin 2) ≠ 0)] at h01
    exact one_ne_zero h01
  simp [h2]

include hq hπmax hf0 hf1 hfres in
/-- `F_f` を `restrictR` で制限した「1変数の中身」——`ψ(X):=F_f(X,0)`。 -/
noncomputable def psi : PowerSeries A :=
  PowerSeries.mk (fun n => MvPowerSeries.coeff (Finsupp.single (0 : Fin 2) n)
    (MvPowerSeries.subst (restrictR (A := A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres)))

include hq hπmax hf0 hf1 hfres in
theorem emb_psi : emb (psi hq hπmax f hf0 hf1 hfres) =
    MvPowerSeries.subst (restrictR (A := A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres) := by
  apply MvPowerSeries.ext
  intro e
  rw [coeff_emb]
  by_cases h1 : e 1 = 0
  · rw [if_pos h1, psi, PowerSeries.coeff_mk]
    have he : e = Finsupp.single (0 : Fin 2) (e 0) := by
      ext i; fin_cases i
      · simp
      · simp [h1]
    rw [← he]
  · rw [if_neg h1]
    exact (coeff_restrictR_eq_zero_of_ne (formalGroupLaw hq hπmax f hf0 hf1 hfres) e h1).symm

include hq hπmax hf0 hf1 hfres in
theorem constantCoeff_psi : PowerSeries.constantCoeff (psi hq hπmax f hf0 hf1 hfres) = 0 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply, psi, PowerSeries.coeff_mk]
  rw [coeff_restrictR_eq_of_e1_zero _ _ (by simp)]
  rw [show Finsupp.single (0 : Fin 2) 0 = (0 : Fin 2 →₀ ℕ) from Finsupp.single_zero _]
  rw [MvPowerSeries.coeff_zero_eq_constantCoeff_apply]
  exact constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres

include hq hπmax hf0 hf1 hfres in
theorem coeff_one_psi : PowerSeries.coeff 1 (psi hq hπmax f hf0 hf1 hfres) = 1 := by
  rw [psi, PowerSeries.coeff_mk]
  rw [coeff_restrictR_eq_of_e1_zero _ _ (by simp)]
  exact coeff_single01_formalGroupLaw hq hπmax f hf0 hf1 hfres

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
theorem subst_zero_eq_zero {p : PowerSeries A} (hp0 : PowerSeries.coeff 0 p = 0) :
    PowerSeries.subst (0 : MvPowerSeries (Fin 2) A) p = 0 :=
  MvPowerSeries.ext (coeff_subst_zero_eq_zero p hp0)

end ABC3.Found.PGC
