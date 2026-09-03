import ABC3.Found.PGC.LubinTateUniqueness

/-!
# Lubin-Tate 形式群法則: 単位元則 `F_f(X,0)=X`(`sorry` 無しで完成)

`Found/PGC/LubinTateUniqueness.lean::powerSeries_uniqueness`(一意性補題)を
`formalGroupLaw`(`F_f`)に実際に適用し、単位元則 `F_f(X,0)=X`
(`formalGroupLaw_identity`)を確立した。

## 確立したもの(`sorry` 無し)

- `restrictR`・`hasSubst_restrictR`・`coeff_restrictR_eq_zero_of_ne`・
  `coeff_restrictR_eq_of_e1_zero`: `X_1↦0` の代入とその係数計算。
- `psi`(`ψ:=F_f(X,0)` の1変数の中身)・`emb_psi`(`emb ψ = F_f.subst
  (restrictR)`)・`constantCoeff_psi`・`coeff_one_psi`(定数項0・次数1の
  係数が1)。
- `constantCoeff_formalGroupLaw`・`coeff_single01_formalGroupLaw`:
  `F_f` 自身の対応する係数(`ΦSeq` まで遡って計算)。
- `subst_zero_eq_zero`: `PowerSeries.subst 0 p = 0`(`p` の定数項が0のとき)。
- `hasSubst_const`・`subst_family_comp_value`: 1変数の値による代入の後に
  2変数の族による代入を合成する一般補題——`PowerSeries.subst_def`
  (`PowerSeries.subst` が定数族による `MvPowerSeries.subst` そのもので
  あること)を経由して `MvPowerSeries.subst_comp_subst_apply` に帰着する。
  当初 finsum の reindexing が要ると見積もっていたが、`subst_def` の
  存在に気づいてからは不要になった。
- ★★★★★★★★★`psi_functional_equation`: `ψ` が `f` との関数等式
  `ψ.subst(f)=f.subst(ψ)` を満たす——`formalGroupLaw_f_isEndomorphism`
  の両辺に `restrictR` を合成し、`subst_family_comp_value` と `subst_emb`
  を組み合わせて `emb` の中で閉じる。
- ★★★★★★★★★★`formalGroupLaw_identity`: **`F_f(X,0)=X`**(単位元則)。
  `ψ` と `X` が同じ関数等式・同じ次数1の係数を持つことから
  `powerSeries_uniqueness` で結論する。

## まだ無いもの

`F_f(0,Y)=Y`(対称な単位元則、`X_0↦0` 版の `restrictR` で同様に示せる
はず)・結合律・可換律——`formalGroupLaw` が実際に(可換)形式群法則
であることを完全に示すにはこれらが要るが、単位元則の一方はこれで
確立できた。 -/

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

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
/-- 定数族への `HasSubst`: `v` の定数項が0なら `fun _ => v` は `HasSubst`。 -/
theorem hasSubst_const {σ : Type*} [Finite σ] {v : MvPowerSeries (Fin 2) A}
    (hv0 : MvPowerSeries.constantCoeff v = 0) :
    MvPowerSeries.HasSubst (fun _ : σ => v) := by
  constructor
  · intro _; show IsNilpotent (MvPowerSeries.constantCoeff v); rw [hv0]; exact IsNilpotent.zero
  · intro d; exact Set.toFinite _

omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
/-- **1変数の値による代入と、その後に続く2変数の族による代入の合成**:
`c` を(`v` を `p` に代入した結果)に代入するのは、まず `c` を `v` に代入して
から `p` に代入するのと同じ——`PowerSeries.subst_def`(`PowerSeries.subst`が
定数族による `MvPowerSeries.subst` そのものであること)経由で
`MvPowerSeries.subst_comp_subst_apply` に帰着する。 -/
theorem subst_family_comp_value {c : Fin 2 → MvPowerSeries (Fin 2) A}
    (hc : MvPowerSeries.HasSubst c) {v : MvPowerSeries (Fin 2) A}
    (hv0 : MvPowerSeries.constantCoeff v = 0) (p : PowerSeries A) :
    MvPowerSeries.subst c (PowerSeries.subst v p) = PowerSeries.subst (MvPowerSeries.subst c v) p := by
  rw [PowerSeries.subst_def v p, MvPowerSeries.subst_comp_subst_apply (hasSubst_const hv0) hc,
    ← PowerSeries.subst_def]

include hq hπmax hf0 hf1 hfres in
/-- ★★★★★★★★★**`ψ:=F_f(X,0)` は `f` との関数等式を満たす**: `ψ(f(X))=f(ψ(X))`。
鍵は `PowerSeries.subst_def`(`PowerSeries.subst a p = MvPowerSeries.subst
(fun _=>a) p`、`PowerSeries.subst` が定数族による `MvPowerSeries.subst`
そのものであること)——これで `restrictR` を経由する族 `c` を
「`emb f` を送る定数族」経由で書き換えられ、finsum の reindexing を
一切経由せずに閉じた。 -/
theorem psi_functional_equation :
    PowerSeries.subst f (psi hq hπmax f hf0 hf1 hfres) =
      PowerSeries.subst (psi hq hπmax f hf0 hf1 hfres) f := by
  apply emb_injective
  have hf0c : PowerSeries.constantCoeff f = 0 := hf0' f hf0
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0c]; exact IsNilpotent.zero
  have hHSψ : PowerSeries.HasSubst (psi hq hπmax f hf0 hf1 hfres) := by
    show IsNilpotent (PowerSeries.constantCoeff (psi hq hπmax f hf0 hf1 hfres))
    rw [constantCoeff_psi]; exact IsNilpotent.zero
  set a : Fin 2 → MvPowerSeries (Fin 2) A :=
    fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f with ha_def
  have hHSa : MvPowerSeries.HasSubst a := hasSubst_g_subst_X f hf0c
  have hembf0 : MvPowerSeries.constantCoeff (emb f) = 0 := by
    rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply, coeff_emb]; simp [hf0]
  have hHSembf : MvPowerSeries.HasSubst (fun _ : Fin 2 => emb f) := hasSubst_const hembf0
  have hrestrict_a : ∀ i : Fin 2, MvPowerSeries.subst (restrictR (A := A)) (a i) =
      PowerSeries.subst (restrictR (A := A) i) f := by
    intro i
    rw [ha_def]
    show MvPowerSeries.subst (restrictR (A := A))
      (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) = _
    rw [subst_family_comp_value (hasSubst_restrictR (A := A)) (MvPowerSeries.constantCoeff_X i) f]
    rw [MvPowerSeries.subst_X (hasSubst_restrictR (A := A)) i]
  -- 段1: `subst restrictR (subst a F_f) = subst c F_f`(`c s := if s=0 then emb f else 0`)。
  have has0 : MvPowerSeries.subst (restrictR (A := A)) (a 0) = emb f := by
    rw [hrestrict_a]
    show PowerSeries.subst (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) f = emb f
    rfl
  have has1 : MvPowerSeries.subst (restrictR (A := A)) (a 1) = 0 := by
    rw [hrestrict_a]
    show PowerSeries.subst (restrictR (A := A) 1) f = 0
    rw [show (restrictR (A := A) 1) = 0 from rfl]
    exact subst_zero_eq_zero hf0
  have hstep1 : MvPowerSeries.subst (restrictR (A := A))
      (MvPowerSeries.subst a (formalGroupLaw hq hπmax f hf0 hf1 hfres)) =
      MvPowerSeries.subst (fun _ : Fin 2 => emb f)
        (MvPowerSeries.subst (restrictR (A := A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres)) := by
    rw [MvPowerSeries.subst_comp_subst_apply hHSa (hasSubst_restrictR (A := A))]
    rw [MvPowerSeries.subst_comp_subst_apply (hasSubst_restrictR (A := A)) hHSembf]
    have hpt : ∀ s : Fin 2, MvPowerSeries.subst (restrictR (A := A)) (a s) =
        MvPowerSeries.subst (fun _ : Fin 2 => emb f) (restrictR (A := A) s) := by
      intro s
      by_cases hs : s = 0
      · subst hs; rw [has0]; exact (MvPowerSeries.subst_X hHSembf 0).symm
      · have hs1 : s = 1 := by omega
        subst hs1
        rw [has1, show (restrictR (A := A)) 1 = 0 from rfl]
        have hz := MvPowerSeries.subst_sub hHSembf (0 : MvPowerSeries (Fin 2) A) 0
        simpa using hz.symm
    exact congrArg (fun g => MvPowerSeries.subst g (formalGroupLaw hq hπmax f hf0 hf1 hfres))
      (funext hpt)
  -- 段2: `subst (fun_=>emb f) (emb ψ) = emb(subst f ψ)`。
  have hstep2 : MvPowerSeries.subst (fun _ : Fin 2 => emb f) (emb (psi hq hπmax f hf0 hf1 hfres)) =
      emb (PowerSeries.subst f (psi hq hπmax f hf0 hf1 hfres)) := by
    show MvPowerSeries.subst (fun _ : Fin 2 => emb f)
      (PowerSeries.subst (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A)
        (psi hq hπmax f hf0 hf1 hfres)) = _
    rw [subst_family_comp_value hHSembf (MvPowerSeries.constantCoeff_X 0)]
    rw [show MvPowerSeries.subst (fun _ : Fin 2 => emb f)
        (MvPowerSeries.X (0 : Fin 2) : MvPowerSeries (Fin 2) A) = emb f from
        MvPowerSeries.subst_X hHSembf 0]
    exact subst_emb f (psi hq hπmax f hf0 hf1 hfres) hHSf
  -- 段3: RHS 側。
  have hstep3 : MvPowerSeries.subst (restrictR (A := A))
      (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hfres) f) =
      emb (PowerSeries.subst (psi hq hπmax f hf0 hf1 hfres) f) := by
    have hFf0 : MvPowerSeries.constantCoeff (formalGroupLaw hq hπmax f hf0 hf1 hfres) = 0 :=
      constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres
    rw [subst_family_comp_value (hasSubst_restrictR (A := A)) hFf0]
    rw [show MvPowerSeries.subst (restrictR (A := A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres) =
        emb (psi hq hπmax f hf0 hf1 hfres) from (emb_psi hq hπmax f hf0 hf1 hfres).symm]
    exact subst_emb (psi hq hπmax f hf0 hf1 hfres) f hHSψ
  have hFf := formalGroupLaw_f_isEndomorphism hq hπmax f hf0 hf1 hfres
  calc emb (PowerSeries.subst f (psi hq hπmax f hf0 hf1 hfres))
      = MvPowerSeries.subst (fun _ : Fin 2 => emb f) (emb (psi hq hπmax f hf0 hf1 hfres)) :=
        hstep2.symm
    _ = MvPowerSeries.subst (fun _ : Fin 2 => emb f)
        (MvPowerSeries.subst (restrictR (A := A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres)) := by
        rw [emb_psi]
    _ = MvPowerSeries.subst (restrictR (A := A))
        (MvPowerSeries.subst a (formalGroupLaw hq hπmax f hf0 hf1 hfres)) := hstep1.symm
    _ = MvPowerSeries.subst (restrictR (A := A))
        (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hfres) f) := by rw [hFf]
    _ = emb (PowerSeries.subst (psi hq hπmax f hf0 hf1 hfres) f) := hstep3

include hq hπmax hf0 hf1 hfres in
/-- ★★★★★★★★★★**Lubin-Tate 形式群法則の単位元則**: `F_f(X,0)=X`。
`ψ:=F_f(X,0)`(`psi`)と `X`(恒等射)がどちらも `f` との同じ関数等式を
満たし(`psi_functional_equation`・`PowerSeries.X_subst`/`subst_X`)、
同じ次数1の係数(`coeff_one_psi`・`PowerSeries.coeff_one_X`)を持つ
ことから、`powerSeries_uniqueness`(一意性補題)で `ψ=X` が出る——
`emb` を経由して `F_f(X,0)=X` を `MvPowerSeries (Fin 2) A` の中で述べる。 -/
theorem formalGroupLaw_identity (hπne0 : π ≠ 0) :
    psi hq hπmax f hf0 hf1 hfres = PowerSeries.X := by
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0' f hf0]; exact IsNilpotent.zero
  apply powerSeries_uniqueness hπmax hπne0 (hf0' f hf0) hf1 (constantCoeff_psi hq hπmax f hf0 hf1 hfres)
    PowerSeries.constantCoeff_X
  · rw [coeff_one_psi, PowerSeries.coeff_one_X]
  · exact psi_functional_equation hq hπmax f hf0 hf1 hfres
  · rw [PowerSeries.X_subst, PowerSeries.subst_X hHSf]

end ABC3.Found.PGC
