import ABC3.Found.PGC.LubinTateCommutativity

/-!
# Lubin-Tate の一意性補題: 任意個の変数への一般化(結合律への土台)

`Found/PGC/LubinTateCommutativity.lean::mvPowerSeries_uniqueness`(2変数版、
`Fin 2` 固定)の証明を読み直すと、`Fin 2` という具体的な添字型には
どこにも本質的に依存していない——`Finsupp.prod` を `Fin.prod_univ_two` で
2項に展開していた箇所(`coeff_ad_eq_of_degree`・`order_prod_pow_sub_prod_pow_ge`)
だけが `Fin 2` 固有だった。ここではその展開を「有限集合上の望遠鏡和」
(`Finset.induction_on`)に置き換えることで、**任意の有限添字型 `σ`**
(`Fin 3` も含む)に対する一意性補題 `mvPowerSeries_uniqueness'` を得る。

結合律 `F_f(F_f(X,Y),Z)=F_f(X,F_f(Y,Z))` は `σ := Fin 3` として
`G:=F_f(F_f(X,Y),Z)` と `H:=F_f(X,F_f(Y,Z))` にこれを適用することで
示せる見込み——次数1の係数が一致すること(`F_f` の次数1部分が対称
だから)と、関数等式(`F_f` 自身の関数等式を2回合成して得る)を示す
段が残っている。
-/

namespace ABC3.Found.PGC

/-! ### 部品0: 次数 `≤1` の `Finsupp` の分類(任意の添字型) -/

/-- `Finsupp.degree d ≤ 1` ならば `d = 0` か、ある `i` について
`d = Finsupp.single i 1`——`Fin 2` に限らず任意の添字型で成り立つ。
`d.support` の濃度が `≤ 1`(各項が `≥1` を足すので `card ≤ degree`)
であることから出す。 -/
theorem finsupp_degree_le_one_cases {σ : Type*} [DecidableEq σ] (d : σ →₀ ℕ)
    (h : Finsupp.degree d ≤ 1) : d = 0 ∨ ∃ i, d = Finsupp.single i 1 := by
  rw [Finsupp.degree_apply] at h
  by_cases hsupp : d.support = ∅
  · left
    rwa [Finsupp.support_eq_empty] at hsupp
  · right
    obtain ⟨i, hi⟩ := Finset.nonempty_iff_ne_empty.mpr hsupp
    have hcard : d.support.card ≤ ∑ j ∈ d.support, d j := by
      calc d.support.card = ∑ _j ∈ d.support, 1 := by rw [Finset.sum_const, smul_eq_mul, mul_one]
        _ ≤ ∑ j ∈ d.support, d j := by
            apply Finset.sum_le_sum
            intro j hj
            exact Nat.one_le_iff_ne_zero.mpr (Finsupp.mem_support_iff.mp hj)
    have hcard1 : d.support.card = 1 := le_antisymm (le_trans hcard h) (Finset.card_pos.mpr ⟨i, hi⟩)
    obtain ⟨i', hi'⟩ := Finset.card_eq_one.mp hcard1
    have hii' : i = i' := by
      have : i ∈ ({i'} : Finset σ) := hi' ▸ hi
      simpa using this
    subst hii'
    refine ⟨i, ?_⟩
    have hdi : d i = 1 := by
      have hle : d i ≤ ∑ j ∈ d.support, d j := Finset.single_le_sum (fun j _ => Nat.zero_le _) hi
      have hine0 : d i ≠ 0 := Finsupp.mem_support_iff.mp hi
      omega
    ext j
    by_cases hji : j = i
    · rw [hji, hdi, Finsupp.single_eq_same]
    · have hjnotsupp : j ∉ d.support := by rw [hi']; simpa using hji
      have hdj0 : d j = 0 := by
        by_contra hc; exact hjnotsupp (Finsupp.mem_support_iff.mpr hc)
      rw [hdj0, eq_comm]
      exact Finsupp.single_eq_of_ne hji

/-! ### 部品1: 有限集合上の望遠鏡和(`Fin 2` 固有の `order_prod_pow_sub_prod_pow_ge`
の一般化) -/

/-- `order_prod_pow_sub_prod_pow_ge`(`Fin 2` 固定・2項)の一般化: 任意の
有限集合 `s` 上の積の差 `∏_{i∈s}a_i^{k_i} − ∏_{i∈s}a_i'^{k_i}` の次数が
`(∑_{i∈s}k_i)+1` 以上であること。`Finset.induction_on` で1項ずつ剥がし、
`(a_j^{k_j}P − a_j'^{k_j}P') = (a_j^{k_j}-a_j'^{k_j})P + a_j'^{k_j}(P-P')`
という2項の場合と同じ望遠鏡分解を繰り返すだけ——`order_pow_sub_pow_ge'`
(既存、`k` に下限なし)と `MvPowerSeries.le_order_prod`(mathlib既存)の
組み合わせで閉じる。 -/
theorem order_prod_pow_sub_prod_pow_ge_finset {A σ : Type*} [CommRing A] [DecidableEq σ]
    {a a' : σ → MvPowerSeries σ A}
    (ha : ∀ i, 1 ≤ (a i).order) (ha' : ∀ i, 1 ≤ (a' i).order)
    (hdiff : ∀ i, (2 : ℕ∞) ≤ (a i - a' i).order) (k : σ → ℕ) (s : Finset σ) :
    ((∑ i ∈ s, k i : ℕ) : ℕ∞) + 1 ≤
      ((∏ i ∈ s, (a i) ^ (k i)) - (∏ i ∈ s, (a' i) ^ (k i))).order := by
  induction s using Finset.induction_on with
  | empty => simp [MvPowerSeries.order_zero]
  | insert j t hjt ih =>
    rw [Finset.prod_insert hjt, Finset.prod_insert hjt, Finset.sum_insert hjt]
    have htelescope : (a j) ^ (k j) * ∏ i ∈ t, (a i) ^ (k i) - (a' j) ^ (k j) * ∏ i ∈ t, (a' i) ^ (k i) =
        ((a j) ^ (k j) - (a' j) ^ (k j)) * (∏ i ∈ t, (a i) ^ (k i)) +
          (a' j) ^ (k j) * ((∏ i ∈ t, (a i) ^ (k i)) - (∏ i ∈ t, (a' i) ^ (k i))) := by
      ring
    rw [htelescope]
    have hPorder : ((∑ i ∈ t, k i : ℕ) : ℕ∞) ≤ (∏ i ∈ t, (a i) ^ (k i)).order := by
      rw [Nat.cast_sum]
      calc ∑ i ∈ t, (k i : ℕ∞) ≤ ∑ i ∈ t, (a i ^ k i).order := by
            apply Finset.sum_le_sum
            intro i _
            calc (k i : ℕ∞) = k i • (1 : ℕ∞) := by simp
              _ ≤ k i • (a i).order := by gcongr; exact ha i
              _ ≤ (a i ^ k i).order := MvPowerSeries.le_order_pow (k i)
        _ ≤ (∏ i ∈ t, (a i) ^ (k i)).order := MvPowerSeries.le_order_prod _ t
    have hajpow_order : (k j : ℕ∞) ≤ ((a' j) ^ (k j)).order := by
      calc (k j : ℕ∞) = k j • (1 : ℕ∞) := by simp
        _ ≤ k j • (a' j).order := by gcongr; exact ha' j
        _ ≤ ((a' j) ^ (k j)).order := MvPowerSeries.le_order_pow (k j)
    have hjdiff_pow : ((2 : ℕ∞) + ((k j - 1 : ℕ) : ℕ∞)) ≤ ((a j) ^ (k j) - (a' j) ^ (k j)).order :=
      le_trans (add_le_add (hdiff j) le_rfl) (order_pow_sub_pow_ge' (ha j) (ha' j) (k j))
    have hkey : (((k j + ∑ i ∈ t, k i : ℕ) : ℕ∞)) + 1 ≤
        ((2 : ℕ∞) + ((k j - 1 : ℕ) : ℕ∞)) + ((∑ i ∈ t, k i : ℕ) : ℕ∞) := by
      have hnat : (k j + ∑ i ∈ t, k i) + 1 ≤ 2 + (k j - 1) + ∑ i ∈ t, k i := by omega
      exact_mod_cast hnat
    have hterm1 : (((k j + ∑ i ∈ t, k i : ℕ) : ℕ∞)) + 1 ≤
        (((a j) ^ (k j) - (a' j) ^ (k j)) * (∏ i ∈ t, (a i) ^ (k i))).order :=
      calc (((k j + ∑ i ∈ t, k i : ℕ) : ℕ∞)) + 1
          ≤ ((2 : ℕ∞) + ((k j - 1 : ℕ) : ℕ∞)) + ((∑ i ∈ t, k i : ℕ) : ℕ∞) := hkey
        _ ≤ ((a j) ^ (k j) - (a' j) ^ (k j)).order + (∏ i ∈ t, (a i) ^ (k i)).order :=
            add_le_add hjdiff_pow hPorder
        _ ≤ _ := MvPowerSeries.le_order_mul
    have hterm2 : (((k j + ∑ i ∈ t, k i : ℕ) : ℕ∞)) + 1 ≤
        ((a' j) ^ (k j) * ((∏ i ∈ t, (a i) ^ (k i)) - (∏ i ∈ t, (a' i) ^ (k i)))).order := by
      have heq : (((k j + ∑ i ∈ t, k i : ℕ) : ℕ∞)) + 1
          = (k j : ℕ∞) + (((∑ i ∈ t, k i : ℕ) : ℕ∞) + 1) := by push_cast; ring
      rw [heq]
      calc (k j : ℕ∞) + (((∑ i ∈ t, k i : ℕ) : ℕ∞) + 1)
          ≤ ((a' j) ^ (k j)).order + (∏ i ∈ t, (a i) ^ (k i) - ∏ i ∈ t, (a' i) ^ (k i)).order :=
            add_le_add hajpow_order ih
        _ ≤ _ := MvPowerSeries.le_order_mul
    exact le_trans (le_min hterm1 hterm2) MvPowerSeries.min_order_le_add

/-! ### 部品2: 次数 `n+1` の斉次式に対する `a^d` と `a'^d` の一致(任意の添字型) -/

/-- `coeff_ad_eq_of_degree`(`Fin 2` 固定)の一般化: `d.prod` の形のまま
(`Fin.prod_univ_two` 等の展開を経由せず)`order_prod_pow_sub_prod_pow_ge_finset`
を `d.support` に適用するだけで閉じる。 -/
theorem coeff_ad_eq_of_degree_finsupp {A σ : Type*} [CommRing A] [DecidableEq σ] {n : ℕ}
    (a a' : σ → MvPowerSeries σ A)
    (hai_order : ∀ i, 1 ≤ (a i).order) (hai'_order : ∀ i, 1 ≤ (a' i).order)
    (hdiff_order : ∀ i, (2 : ℕ∞) ≤ (a i - a' i).order)
    (e : σ →₀ ℕ) (he : Finsupp.degree e ≤ n + 1)
    (d : σ →₀ ℕ) (hd : Finsupp.degree d = n + 1) :
    MvPowerSeries.coeff e (d.prod fun i m => (a i) ^ m) =
      MvPowerSeries.coeff e (d.prod fun i m => (a' i) ^ m) := by
  have hprod_eq : ∀ b : σ → MvPowerSeries σ A,
      d.prod (fun i m => (b i) ^ m) = ∏ i ∈ d.support, (b i) ^ (d i) := fun b => rfl
  rw [hprod_eq a, hprod_eq a']
  have hge : ((Finsupp.degree d : ℕ) : ℕ∞) + 1 ≤
      ((∏ i ∈ d.support, (a i) ^ (d i)) - (∏ i ∈ d.support, (a' i) ^ (d i))).order := by
    have hkey := order_prod_pow_sub_prod_pow_ge_finset hai_order hai'_order hdiff_order d d.support
    rwa [← Finsupp.degree_apply] at hkey
  rw [hd] at hge
  have hlt : ((Finsupp.degree e : ℕ) : ℕ∞) <
      ((∏ i ∈ d.support, (a i) ^ (d i)) - (∏ i ∈ d.support, (a' i) ^ (d i))).order := by
    calc ((Finsupp.degree e : ℕ) : ℕ∞) ≤ ((n + 1 : ℕ) : ℕ∞) := by exact_mod_cast he
      _ < ((n + 2 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : n + 1 < n + 2)
      _ = ((n + 1 : ℕ) : ℕ∞) + 1 := by push_cast; ring
      _ ≤ _ := hge
  have hz := MvPowerSeries.coeff_of_lt_order hlt
  rw [map_sub] at hz
  exact sub_eq_zero.mp hz

/-! ### 部品3: `g.subst(X_i)` と `π•X_i` のずれ(任意の添字型) -/

/-- `coeff_a_diff_order`(`Fin 2` 固定)の一般化——証明は添字型に依存しない
(`hasSubst_X_i` も既に任意の添字型で使える)。 -/
theorem coeff_a_diff_order_general {A σ : Type*} [CommRing A] (g : PowerSeries A) (π : A)
    (hg0 : PowerSeries.constantCoeff g = 0) (hg1 : PowerSeries.coeff 1 g = π) (i : σ) :
    (2 : ℕ∞) ≤ (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g -
      π • (MvPowerSeries.X i : MvPowerSeries σ A)).order := by
  have hXi : PowerSeries.HasSubst (MvPowerSeries.X i : MvPowerSeries σ A) := hasSubst_X_i i
  have heq : PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g -
      π • (MvPowerSeries.X i : MvPowerSeries σ A)
      = PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) (g - π • PowerSeries.X) := by
    rw [PowerSeries.subst_sub hXi, PowerSeries.subst_smul hXi]
    have hXsubst :
        PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) (PowerSeries.X : PowerSeries A)
        = (MvPowerSeries.X i : MvPowerSeries σ A) := PowerSeries.subst_X hXi
    rw [hXsubst]
  rw [heq]
  have hbound : (2 : ℕ∞) ≤ (g - π • PowerSeries.X : PowerSeries A).order := by
    apply PowerSeries.nat_le_order
    intro k hk
    interval_cases k
    · rw [map_sub]
      have h0 : PowerSeries.coeff 0 g = 0 := by
        rw [PowerSeries.coeff_zero_eq_constantCoeff_apply, hg0]
      rw [h0]; simp
    · rw [map_sub, hg1]; simp [PowerSeries.coeff_one_X]
  calc (2 : ℕ∞) ≤ (g - π • PowerSeries.X : PowerSeries A).order := hbound
    _ ≤ (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) (g - π • PowerSeries.X)).order :=
        PowerSeries.le_order_subst_left (MvPowerSeries.constantCoeff_X i)

/-! ### 部品4: 斉次式のスケーリング(任意の添字型) -/

/-- `homogeneous_subst_const_smul`(`Fin 2` 固定)の一般化——`MvPowerSeries.rescale`
の一般理論(`Function.const σ c`)がそのまま使える。 -/
theorem homogeneous_subst_const_smul_general {A σ : Type*} [CommRing A] {φ : MvPowerSeries σ A}
    {p : ℕ} (hφ : ∀ d : σ →₀ ℕ, Finsupp.degree d ≠ p → MvPowerSeries.coeff d φ = 0) (c : A) :
    MvPowerSeries.subst (fun i => c • (MvPowerSeries.X i : MvPowerSeries σ A)) φ = c ^ p • φ := by
  have hsupp : ∀ d ∈ (φ : MvPowerSeries σ A).support, Finsupp.degree d = p := by
    intro d hd
    by_contra hne
    exact hd (hφ d hne)
  have h1 : MvPowerSeries.rescale (Function.const σ c) φ = c ^ p • φ :=
    MvPowerSeries.rescale_homogeneous_eq_smul hsupp
  rw [← h1, MvPowerSeries.rescale_eq_subst]
  congr 1

/-! ### 部品5: g側の線形化(任意の添字型) -/

/-- `coeff_subst_g_linearize`(`Fin 2` 固定)の一般化。 -/
theorem coeff_subst_g_linearize_general {A σ : Type*} [CommRing A] [Finite σ] [DecidableEq σ]
    {φ : MvPowerSeries σ A} {n : ℕ}
    (hφ : ∀ d : σ →₀ ℕ, Finsupp.degree d ≠ n + 1 → MvPowerSeries.coeff d φ = 0)
    {g : PowerSeries A} (π : A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (e : σ →₀ ℕ) (he : Finsupp.degree e ≤ n + 1) :
    MvPowerSeries.coeff e (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) φ)
      = π ^ (n + 1) * MvPowerSeries.coeff e φ := by
  set a : σ → MvPowerSeries σ A := fun i => PowerSeries.subst (MvPowerSeries.X i) g with ha_def
  set a' : σ → MvPowerSeries σ A :=
    fun i => π • (MvPowerSeries.X i : MvPowerSeries σ A) with ha'_def
  have haHS : MvPowerSeries.HasSubst a := hasSubst_g_subst_X g hg0
  have hHSa' : MvPowerSeries.HasSubst a' := by
    constructor
    · intro i
      show IsNilpotent (MvPowerSeries.constantCoeff (a' i))
      show IsNilpotent (MvPowerSeries.constantCoeff
        (π • (MvPowerSeries.X i : MvPowerSeries σ A)))
      rw [MvPowerSeries.constantCoeff_smul, MvPowerSeries.constantCoeff_X, smul_zero]
      exact IsNilpotent.zero
    · intro d; exact Set.toFinite _
  have hai_order : ∀ i, 1 ≤ (a i).order := by
    intro i
    apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
    exact PowerSeries.constantCoeff_subst_eq_zero (MvPowerSeries.constantCoeff_X i) g hg0
  have hai'_order : ∀ i, 1 ≤ (a' i).order := by
    intro i
    apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
    show MvPowerSeries.constantCoeff (π • (MvPowerSeries.X i : MvPowerSeries σ A)) = 0
    rw [MvPowerSeries.constantCoeff_smul, MvPowerSeries.constantCoeff_X, smul_zero]
  have hdiff_order : ∀ i, (2 : ℕ∞) ≤ (a i - a' i).order := fun i => coeff_a_diff_order_general g π hg0 hg1 i
  have hstep2 :
      MvPowerSeries.coeff e (MvPowerSeries.subst a φ) = MvPowerSeries.coeff e (MvPowerSeries.subst a' φ) := by
    rw [MvPowerSeries.coeff_subst haHS, MvPowerSeries.coeff_subst hHSa']
    refine finsum_congr (fun d => ?_)
    by_cases hd : Finsupp.degree d = n + 1
    · congr 1
      exact coeff_ad_eq_of_degree_finsupp a a' hai_order hai'_order hdiff_order e he d hd
    · rw [hφ d hd]; simp
  have hstep1 : MvPowerSeries.coeff e (MvPowerSeries.subst a' φ) = π ^ (n + 1) * MvPowerSeries.coeff e φ := by
    have hscale := homogeneous_subst_const_smul_general hφ π
    rw [ha'_def]
    show MvPowerSeries.coeff e (MvPowerSeries.subst
      (fun i => π • (MvPowerSeries.X i : MvPowerSeries σ A)) φ) = _
    rw [hscale, MvPowerSeries.coeff_smul]
  rw [ha_def, hstep2, hstep1]

/-! ### 部品6: g側の線形化(order版、任意の添字型) -/

/-- `coeff_subst_g_linearize_order`(`Fin 2` 固定)の一般化。 -/
theorem coeff_subst_g_linearize_order_general {A σ : Type*} [CommRing A] [Finite σ] [DecidableEq σ]
    {φ : MvPowerSeries σ A} {n : ℕ}
    (hφorder : ((n + 1 : ℕ) : ℕ∞) ≤ φ.order)
    {g : PowerSeries A} (π : A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (e : σ →₀ ℕ) (he : Finsupp.degree e ≤ n + 1) :
    MvPowerSeries.coeff e (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g) φ)
      = π ^ (n + 1) * MvPowerSeries.coeff e φ := by
  rcases lt_or_eq_of_le he with hlt | heq
  · have hφe0 : MvPowerSeries.coeff e φ = 0 := by
      apply MvPowerSeries.coeff_of_lt_order
      calc ((Finsupp.degree e : ℕ) : ℕ∞) < ((n + 1 : ℕ) : ℕ∞) := by exact_mod_cast hlt
        _ ≤ φ.order := hφorder
    have haHS := hasSubst_g_subst_X (σ := σ) g hg0
    have hai_order : ∀ i, 1 ≤ (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g).order :=
      fun i => by
        apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
        exact PowerSeries.constantCoeff_subst_eq_zero (MvPowerSeries.constantCoeff_X i) g hg0
    have hsubst_order : ((n + 1 : ℕ) : ℕ∞) ≤
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g) φ).order := by
      have hle := MvPowerSeries.le_order_subst haHS φ
      have hinf : (1 : ℕ∞) ≤ ⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g).order :=
        le_iInf hai_order
      calc ((n + 1 : ℕ) : ℕ∞) = 1 * ((n + 1 : ℕ) : ℕ∞) := by ring
        _ ≤ (⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g).order) *
              ((n + 1 : ℕ) : ℕ∞) := by gcongr
        _ ≤ _ := by gcongr
        _ ≤ _ := hle
    have hlhs0 : MvPowerSeries.coeff e
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g) φ) = 0 := by
      apply MvPowerSeries.coeff_of_lt_order
      calc ((Finsupp.degree e : ℕ) : ℕ∞) < ((n + 1 : ℕ) : ℕ∞) := by exact_mod_cast hlt
        _ ≤ _ := hsubst_order
    rw [hlhs0, hφe0, mul_zero]
  · set φ1 := MvPowerSeries.homogeneousComponent (n + 1) φ with hφ1_def
    set tail := φ - φ1 with htail_def
    have hφ1hom : ∀ d : σ →₀ ℕ, Finsupp.degree d ≠ n + 1 → MvPowerSeries.coeff d φ1 = 0 :=
      fun d hd => MvPowerSeries.IsHomogeneous.coeff_eq_zero
        (MvPowerSeries.isHomogeneous_homogeneousComponent φ (n + 1)) hd
    have hφ1e : MvPowerSeries.coeff e φ1 = MvPowerSeries.coeff e φ := by
      rw [hφ1_def, MvPowerSeries.coeff_homogeneousComponent, if_pos heq]
    have htail_order : ((n + 2 : ℕ) : ℕ∞) ≤ tail.order := by
      apply MvPowerSeries.nat_le_order
      intro d hd
      rw [htail_def, map_sub]
      by_cases hdeq : Finsupp.degree d = n + 1
      · rw [hφ1_def, MvPowerSeries.coeff_homogeneousComponent, if_pos hdeq, sub_self]
      · have hdlt : Finsupp.degree d < n + 1 := by omega
        have hφd0 : MvPowerSeries.coeff d φ = 0 := by
          apply MvPowerSeries.coeff_of_lt_order
          calc ((Finsupp.degree d : ℕ) : ℕ∞) < ((n + 1 : ℕ) : ℕ∞) := by exact_mod_cast hdlt
            _ ≤ φ.order := hφorder
        rw [hφd0, hφ1hom d hdeq, sub_zero]
    have haHS := hasSubst_g_subst_X (σ := σ) g hg0
    have hai_order : ∀ i, 1 ≤ (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g).order :=
      fun i => by
        apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
        exact PowerSeries.constantCoeff_subst_eq_zero (MvPowerSeries.constantCoeff_X i) g hg0
    have htailsubst_order : ((n + 2 : ℕ) : ℕ∞) ≤
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g) tail).order := by
      have hle := MvPowerSeries.le_order_subst haHS tail
      have hinf : (1 : ℕ∞) ≤ ⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g).order :=
        le_iInf hai_order
      calc ((n + 2 : ℕ) : ℕ∞) = 1 * ((n + 2 : ℕ) : ℕ∞) := by ring
        _ ≤ (⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g).order) *
              ((n + 2 : ℕ) : ℕ∞) := by gcongr
        _ ≤ _ := by gcongr
        _ ≤ _ := hle
    have htailcoeff0 : MvPowerSeries.coeff e
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g) tail) = 0 := by
      apply MvPowerSeries.coeff_of_lt_order
      calc ((Finsupp.degree e : ℕ) : ℕ∞) < ((n + 2 : ℕ) : ℕ∞) := by
            rw [heq]; exact_mod_cast (by omega : n + 1 < n + 2)
        _ ≤ _ := htailsubst_order
    have hsplit : φ = φ1 + tail := by rw [htail_def]; ring
    have hφsubst : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g) φ =
        MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g) φ1 +
          MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) g) tail := by
      rw [hsplit, MvPowerSeries.subst_add haHS]
    rw [hφsubst, map_add, htailcoeff0, add_zero]
    rw [coeff_subst_g_linearize_general hφ1hom π hg0 hg1 e (by omega)]
    rw [hφ1e]

/-! ### ★★★★★★★★★★★★本題: 一意性補題の任意の添字型への一般化 -/

/-- ★★★★★★★★★★★★**Lubin-Tate の一意性補題(任意の添字型版)**。
`mvPowerSeries_uniqueness`(`Fin 2` 固定の2変数版)と全く同じ証明が、
`σ := Fin 3` を含む任意の有限添字型に対してそのまま成り立つ——
結合律 `F_f(F_f(X,Y),Z)=F_f(X,F_f(Y,Z))` を示すための土台となる。 -/
theorem mvPowerSeries_uniqueness_general {A σ : Type*} [CommRing A] [Finite σ] [DecidableEq σ]
    [IsLocalRing A] [IsDomain A] {π : A}
    (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    {f : PowerSeries A} (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    {Φ Ψ : MvPowerSeries σ A} (hΦ0 : MvPowerSeries.constantCoeff Φ = 0)
    (hΨ0 : MvPowerSeries.constantCoeff Ψ = 0)
    (hlead : ∀ i : σ, MvPowerSeries.coeff (Finsupp.single i 1) Φ =
      MvPowerSeries.coeff (Finsupp.single i 1) Ψ)
    (heqΦ : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) f) Φ =
      PowerSeries.subst Φ f)
    (heqΨ : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) f) Ψ =
      PowerSeries.subst Ψ f) :
    Φ = Ψ := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0c]; exact IsNilpotent.zero
  have hHSa : MvPowerSeries.HasSubst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) f) :=
    hasSubst_g_subst_X f hf0c
  set δ := Φ - Ψ with hδ_def
  have hδ0 : MvPowerSeries.constantCoeff δ = 0 := by rw [hδ_def, map_sub, hΦ0, hΨ0, sub_zero]
  have hδlead : ∀ i : σ, MvPowerSeries.coeff (Finsupp.single i 1) δ = 0 := by
    intro i; rw [hδ_def, map_sub, hlead i, sub_self]
  have hbase : ((2 : ℕ) : ℕ∞) ≤ δ.order := by
    apply MvPowerSeries.nat_le_order
    intro d hd
    have hddeg : Finsupp.degree d < 2 := by exact_mod_cast hd
    rcases finsupp_degree_le_one_cases d (by omega) with h0 | ⟨i, hi⟩
    · rw [h0]
      rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply] at hδ0
      exact hδ0
    · rw [hi]; exact hδlead i
  have hstep : ∀ n : ℕ, 1 ≤ n → ((n + 1 : ℕ) : ℕ∞) ≤ δ.order →
      ((n + 2 : ℕ) : ℕ∞) ≤ δ.order := by
    intro n hn1 hδorder
    apply MvPowerSeries.nat_le_order
    intro e he
    have hedeg' : Finsupp.degree e < n + 2 := by exact_mod_cast he
    by_cases hcase : Finsupp.degree e < n + 1
    · exact MvPowerSeries.coeff_of_lt_order
        (lt_of_lt_of_le (by exact_mod_cast hcase) hδorder)
    have hedeg : Finsupp.degree e = n + 1 := by omega
    have hΨorder : (1 : ℕ∞) ≤ Ψ.order := MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hΨ0
    have hΨδeq : Ψ + δ = Φ := by rw [hδ_def]; ring
    have hlin := coeff_subst_linearize hΨ0 hδ0 hΨorder hδorder (by omega : 1 ≤ n + 1) f π hf0 hf1 e
      (by rw [hedeg])
    rw [hΨδeq] at hlin
    have hglin := coeff_subst_g_linearize_order_general (φ := δ) hδorder π hf0c hf1 e (by omega)
    have hδsubst : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) f) δ =
        MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) f) Φ -
          MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries σ A) f) Ψ := by
      rw [hδ_def, MvPowerSeries.subst_sub hHSa]
    rw [hδsubst, map_sub, heqΦ, heqΨ] at hglin
    rw [hglin] at hlin
    have hfactor : MvPowerSeries.coeff e δ * (π ^ (n + 1) - π) = 0 := by
      have h2 := hlin
      ring_nf
      ring_nf at h2
      linear_combination h2
    have hπn_mem : π ^ n ∈ IsLocalRing.maximalIdeal A :=
      Ideal.pow_mem_of_mem _ (hπmax ▸ Ideal.mem_span_singleton_self π) n (by omega)
    have hπn_ne_one : π ^ n ≠ 1 := fun heqp => by
      rw [heqp] at hπn_mem
      exact IsLocalRing.maximalIdeal.isMaximal A |>.ne_top
        (Ideal.eq_top_of_isUnit_mem _ hπn_mem isUnit_one)
    have hne : π ^ (n + 1) - π ≠ 0 := by
      intro hcon
      apply hπn_ne_one
      have h3 : π * (π ^ n - 1) = 0 := by ring_nf; ring_nf at hcon; linear_combination hcon
      rcases mul_eq_zero.mp h3 with h4 | h5
      · exact absurd h4 hπne0
      · exact sub_eq_zero.mp h5
    exact (mul_eq_zero.mp hfactor).resolve_right hne
  have hall : ∀ n : ℕ, 1 ≤ n → ((n + 1 : ℕ) : ℕ∞) ≤ δ.order := by
    intro n hn
    induction n, hn using Nat.le_induction with
    | base => exact hbase
    | succ n hn1 ih => exact hstep n hn1 ih
  have horder : δ.order = ⊤ := by
    by_contra hne
    obtain ⟨m, hm⟩ := WithTop.ne_top_iff_exists.mp hne
    have hcontra := hall (m + 1) (by omega)
    rw [← hm] at hcontra
    have : ((m : ℕ) : ℕ∞) < ((m + 1 + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega)
    exact absurd hcontra (not_le.mpr this)
  have hδzero : δ = 0 := MvPowerSeries.order_eq_top_iff.mp horder
  rw [hδ_def, sub_eq_zero] at hδzero
  exact hδzero

end ABC3.Found.PGC
