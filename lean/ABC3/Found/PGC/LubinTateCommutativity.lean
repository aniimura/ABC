import ABC3.Found.PGC.LubinTateIdentityLawSymm

/-!
# Lubin-Tate 形式群法則: 可換律 `F_f(X,Y)=F_f(Y,X)`(進行中)

`Found/PGC/LubinTateUniqueness.lean::powerSeries_uniqueness`(1変数の
一意性補題)の2変数版を確立し、`swap(F_f):=F_f(Y,X)` が `F_f` 自身と
同じ関数等式を満たすことと組み合わせて可換律を示す計画。
-/

namespace ABC3.Found.PGC

variable {A : Type*} [CommRing A]

/-- `coeff_subst_g_linearize` の一般化: `φ` が**次数 `n+1` の斉次式**である
ことを要求せず、**次数 `≥n+1`**(次数 `n+1` 以上の項しか持たない)だけで
足りる——次数 `n+1` を超える部分は代入後の次数を押し上げるので、
求める係数(次数 `e=n+1`)には効かない。`homogeneousComponent` で
次数 `n+1` の部分を取り出し、`coeff_subst_g_linearize` を直接適用し、
残り(`tail`、次数 `≥n+2`)が効かないことを次数の議論で示す。 -/
theorem coeff_subst_g_linearize_order {φ : MvPowerSeries (Fin 2) A} {n : ℕ}
    (hφorder : ((n + 1 : ℕ) : ℕ∞) ≤ φ.order)
    {g : PowerSeries A} (π : A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (e : Fin 2 →₀ ℕ) (he : Finsupp.degree e ≤ n + 1) :
    MvPowerSeries.coeff e (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g) φ)
      = π ^ (n + 1) * MvPowerSeries.coeff e φ := by
  rcases lt_or_eq_of_le he with hlt | heq
  · -- 次数 e < n+1: 両辺とも0。
    have hφe0 : MvPowerSeries.coeff e φ = 0 := by
      apply MvPowerSeries.coeff_of_lt_order
      calc ((Finsupp.degree e : ℕ) : ℕ∞) < ((n + 1 : ℕ) : ℕ∞) := by exact_mod_cast hlt
        _ ≤ φ.order := hφorder
    have haHS := hasSubst_g_subst_X g hg0
    have hai_order : ∀ i, 1 ≤ (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g).order :=
      hai_order_g g hg0
    have hsubst_order : ((n + 1 : ℕ) : ℕ∞) ≤
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g) φ).order := by
      have hle := MvPowerSeries.le_order_subst haHS φ
      have hinf : (1 : ℕ∞) ≤ ⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g).order :=
        le_iInf hai_order
      calc ((n + 1 : ℕ) : ℕ∞) = 1 * ((n + 1 : ℕ) : ℕ∞) := by ring
        _ ≤ (⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g).order) *
              ((n + 1 : ℕ) : ℕ∞) := by gcongr
        _ ≤ _ := by gcongr
        _ ≤ _ := hle
    have hlhs0 : MvPowerSeries.coeff e
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g) φ) = 0 := by
      apply MvPowerSeries.coeff_of_lt_order
      calc ((Finsupp.degree e : ℕ) : ℕ∞) < ((n + 1 : ℕ) : ℕ∞) := by exact_mod_cast hlt
        _ ≤ _ := hsubst_order
    rw [hlhs0, hφe0, mul_zero]
  · -- 次数 e = n+1: `homogeneousComponent` で切り出して `coeff_subst_g_linearize` を適用。
    set φ1 := MvPowerSeries.homogeneousComponent (n + 1) φ with hφ1_def
    set tail := φ - φ1 with htail_def
    have hφ1hom : ∀ d : Fin 2 →₀ ℕ, Finsupp.degree d ≠ n + 1 → MvPowerSeries.coeff d φ1 = 0 :=
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
    have haHS := hasSubst_g_subst_X g hg0
    have hai_order : ∀ i, 1 ≤ (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g).order :=
      hai_order_g g hg0
    have htailsubst_order : ((n + 2 : ℕ) : ℕ∞) ≤
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g) tail).order := by
      have hle := MvPowerSeries.le_order_subst haHS tail
      have hinf : (1 : ℕ∞) ≤ ⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g).order :=
        le_iInf hai_order
      calc ((n + 2 : ℕ) : ℕ∞) = 1 * ((n + 2 : ℕ) : ℕ∞) := by ring
        _ ≤ (⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g).order) *
              ((n + 2 : ℕ) : ℕ∞) := by gcongr
        _ ≤ _ := by gcongr
        _ ≤ _ := hle
    have htailcoeff0 : MvPowerSeries.coeff e
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g) tail) = 0 := by
      apply MvPowerSeries.coeff_of_lt_order
      calc ((Finsupp.degree e : ℕ) : ℕ∞) < ((n + 2 : ℕ) : ℕ∞) := by
            rw [heq]; exact_mod_cast (by omega : n + 1 < n + 2)
        _ ≤ _ := htailsubst_order
    have hsplit : φ = φ1 + tail := by rw [htail_def]; ring
    have hφsubst : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g) φ =
        MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g) φ1 +
          MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g) tail := by
      rw [hsplit, MvPowerSeries.subst_add haHS]
    rw [hφsubst, map_add, htailcoeff0, add_zero]
    rw [coeff_subst_g_linearize hφ1hom π hg0 hg1 e (by omega)]
    rw [hφ1e]

/-- ★★★★★★★★★**Lubin-Tate の一意性補題(2変数版)**。`f≡πX(mod deg2)` に
対し、次数1の係数が一致し(かつ定数項0)、同じ関数等式
`Φ.subst(f,f)=f.subst(Φ)`・`Ψ.subst(f,f)=f.subst(Ψ)` を満たす2つの
2変数冪級数 `Φ,Ψ` は一致する——`Found/PGC/LubinTateUniqueness.lean::
powerSeries_uniqueness`(1変数版)と全く同じ次数ごとの帰納法。`emb` を
経由せず直接 `MvPowerSeries (Fin 2) A` の中で行える——f側の線形化は
既存の `coeff_subst_linearize`、g側(ここでは `f` 自身が両方の役)は
今回一般化した `coeff_subst_g_linearize_order` をそのまま使う。 -/
theorem mvPowerSeries_uniqueness [IsLocalRing A] [IsDomain A] {π : A}
    (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    {f : PowerSeries A} (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    {Φ Ψ : MvPowerSeries (Fin 2) A} (hΦ0 : MvPowerSeries.constantCoeff Φ = 0)
    (hΨ0 : MvPowerSeries.constantCoeff Ψ = 0)
    (hlead : ∀ i : Fin 2, MvPowerSeries.coeff (Finsupp.single i 1) Φ =
      MvPowerSeries.coeff (Finsupp.single i 1) Ψ)
    (heqΦ : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) f) Φ =
      PowerSeries.subst Φ f)
    (heqΨ : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) f) Ψ =
      PowerSeries.subst Ψ f) :
    Φ = Ψ := by
  have hf0c : PowerSeries.constantCoeff f = 0 := hf0' f hf0
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0c]; exact IsNilpotent.zero
  have hHSa : MvPowerSeries.HasSubst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) :=
    hasSubst_g_subst_X f hf0c
  set δ := Φ - Ψ with hδ_def
  have hδ0 : MvPowerSeries.constantCoeff δ = 0 := by rw [hδ_def, map_sub, hΦ0, hΨ0, sub_zero]
  have hδlead : ∀ i : Fin 2, MvPowerSeries.coeff (Finsupp.single i 1) δ = 0 := by
    intro i; rw [hδ_def, map_sub, hlead i, sub_self]
  have hbase : ((2 : ℕ) : ℕ∞) ≤ δ.order := by
    apply MvPowerSeries.nat_le_order
    intro d hd
    have hddeg : Finsupp.degree d < 2 := by exact_mod_cast hd
    have hdcase : Finsupp.degree d = 0 ∨ d = Finsupp.single (0 : Fin 2) 1 ∨
        d = Finsupp.single (1 : Fin 2) 1 := by
      have hd01 : d 0 + d 1 ≤ 1 := by rw [← finsupp_degree_fin2]; omega
      have hde : d = Finsupp.single (0 : Fin 2) (d 0) + Finsupp.single (1 : Fin 2) (d 1) := by
        ext i; fin_cases i <;> simp
      rcases (by omega : d 0 = 0 ∧ d 1 = 0 ∨ d 0 = 1 ∧ d 1 = 0 ∨ d 0 = 0 ∧ d 1 = 1)
        with ⟨h0, h1⟩ | ⟨h0, h1⟩ | ⟨h0, h1⟩
      · left; rw [finsupp_degree_fin2, h0, h1]
      · right; left; rw [hde, h0, h1]; simp
      · right; right; rw [hde, h0, h1]; simp
    rcases hdcase with h | h | h
    · rw [Finsupp.degree_eq_zero_iff] at h
      rw [h]
      rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply] at hδ0
      exact hδ0
    · rw [h]; exact hδlead 0
    · rw [h]; exact hδlead 1
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
    have hglin := coeff_subst_g_linearize_order (φ := δ) hδorder π hf0c hf1 e (by omega)
    have hδsubst : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) δ =
        MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) Φ -
          MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) Ψ := by
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

/-! ### 変数の入れ替え `swap` -/

/-- `X_0` と `X_1` を入れ替える代入 family。 -/
noncomputable def swapFam : Fin 2 → MvPowerSeries (Fin 2) A :=
  fun i => if i = 0 then (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) else MvPowerSeries.X 0

theorem hasSubst_swapFam : MvPowerSeries.HasSubst (swapFam (A := A)) := by
  constructor
  · intro i
    show IsNilpotent (MvPowerSeries.constantCoeff (swapFam i))
    unfold swapFam
    split_ifs <;> (rw [MvPowerSeries.constantCoeff_X]; exact IsNilpotent.zero)
  · intro d; exact Set.toFinite _

/-- `X_0↔X_1` を入れ替える写像。 -/
noncomputable def swap (Φ : MvPowerSeries (Fin 2) A) : MvPowerSeries (Fin 2) A :=
  MvPowerSeries.subst (swapFam (A := A)) Φ

theorem swapFam_zero : (swapFam (A := A)) 0 = (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) := rfl

theorem swapFam_one : (swapFam (A := A)) 1 = (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) := rfl

/-- ★★★★★★★★★**`swap` は関数等式を保つ**: `Φ` が `Φ.subst(a)=f.subst(Φ)`
(`a:=fun i=>f.subst(X_i)`)を満たせば、`swap Φ` も同じ関数等式を満たす。
`swap(Φ).subst(a) = Φ.subst(a') = swap(Φ.subst(a)) = swap(f.subst(Φ))
= f.subst(swap(Φ))`(`a'` は `a` の成分を入れ替えたもの)という
`subst_comp_subst_apply` の連鎖で示す。 -/
theorem swap_preserves_functional_equation {f : PowerSeries A} (hf0 : PowerSeries.coeff 0 f = 0)
    {Φ : MvPowerSeries (Fin 2) A} (hΦ0 : MvPowerSeries.constantCoeff Φ = 0)
    (heqΦ : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) Φ =
      PowerSeries.subst Φ f) :
    MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) (swap Φ) =
      PowerSeries.subst (swap Φ) f := by
  have hf0c : PowerSeries.constantCoeff f = 0 := hf0' f hf0
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0c]; exact IsNilpotent.zero
  have hHSa : MvPowerSeries.HasSubst
      (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) :=
    hasSubst_g_subst_X f hf0c
  set a : Fin 2 → MvPowerSeries (Fin 2) A :=
    fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f with ha_def
  -- 段1: `swap(Φ).subst(a) = Φ.subst(a')`(`a'` は `a` の成分を入れ替えたもの)。
  have hstep1 : MvPowerSeries.subst a (swap Φ) =
      MvPowerSeries.subst (fun i => a (if i = 0 then 1 else 0)) Φ := by
    rw [swap, MvPowerSeries.subst_comp_subst_apply hasSubst_swapFam hHSa]
    congr 1
    funext i
    fin_cases i
    · show MvPowerSeries.subst a (swapFam 0) = a 1
      rw [swapFam_zero, MvPowerSeries.subst_X hHSa 1]
    · show MvPowerSeries.subst a (swapFam 1) = a 0
      rw [swapFam_one, MvPowerSeries.subst_X hHSa 0]
  -- 段2: `Φ.subst(a') = swap(Φ.subst(a))`。
  have hstep2 : MvPowerSeries.subst (fun i => a (if i = 0 then 1 else 0)) Φ =
      swap (MvPowerSeries.subst a Φ) := by
    rw [swap, MvPowerSeries.subst_comp_subst_apply hHSa hasSubst_swapFam]
    congr 1
    funext i
    fin_cases i
    · show a (if (0 : Fin 2) = 0 then 1 else 0) = MvPowerSeries.subst (swapFam (A := A)) (a 0)
      rw [ha_def]
      show PowerSeries.subst (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) f =
        MvPowerSeries.subst (swapFam (A := A)) (PowerSeries.subst (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) f)
      rw [subst_family_comp_value hasSubst_swapFam (MvPowerSeries.constantCoeff_X 0) f,
        MvPowerSeries.subst_X hasSubst_swapFam 0, swapFam_zero]
    · show a (if (1 : Fin 2) = 0 then 1 else 0) = MvPowerSeries.subst (swapFam (A := A)) (a 1)
      rw [ha_def]
      show PowerSeries.subst (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) f =
        MvPowerSeries.subst (swapFam (A := A)) (PowerSeries.subst (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) f)
      rw [subst_family_comp_value hasSubst_swapFam (MvPowerSeries.constantCoeff_X 1) f,
        MvPowerSeries.subst_X hasSubst_swapFam 1, swapFam_one]
  -- 段3: `swap(f.subst(Φ)) = f.subst(swap(Φ))`。
  have hstep3 : swap (PowerSeries.subst Φ f) = PowerSeries.subst (swap Φ) f := by
    rw [swap, swap, subst_family_comp_value hasSubst_swapFam hΦ0 f]
  rw [ha_def] at hstep1 ⊢
  rw [hstep1, hstep2, heqΦ, hstep3]

theorem constantCoeff_swap {Φ : MvPowerSeries (Fin 2) A} (hΦ0 : MvPowerSeries.constantCoeff Φ = 0) :
    MvPowerSeries.constantCoeff (swap Φ) = 0 :=
  MvPowerSeries.constantCoeff_subst_eq_zero hasSubst_swapFam
    (fun i => by unfold swapFam; split_ifs <;> exact MvPowerSeries.constantCoeff_X _) hΦ0

theorem swapFam_order (i : Fin 2) : (1 : ℕ∞) ≤ ((swapFam (A := A)) i).order := by
  unfold swapFam; split_ifs <;> simp [MvPowerSeries.one_le_order_iff_constCoeff_eq_zero,
    MvPowerSeries.constantCoeff_X]

/-- `d.prod(swapFam^·)` の次数は `≥ degree d`——`swapFam` の各成分の次数が
`≥1` であることから。 -/
theorem order_prod_swapFam_ge (d : Fin 2 →₀ ℕ) :
    ((Finsupp.degree d : ℕ) : ℕ∞) ≤ (d.prod (fun s m => ((swapFam (A := A)) s) ^ m)).order := by
  rw [finsupp_degree_fin2]
  have hprod : d.prod (fun s m => ((swapFam (A := A)) s) ^ m) =
      ((swapFam (A := A)) 0) ^ (d 0) * ((swapFam (A := A)) 1) ^ (d 1) := by
    show (∏ j ∈ d.support, ((swapFam (A := A)) j) ^ (d j)) = _
    rw [Finset.prod_subset (Finset.subset_univ d.support) (fun x _ hx => by
      simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two]
  rw [hprod]
  have h0 : ((d 0 : ℕ) : ℕ∞) ≤ (((swapFam (A := A)) 0) ^ (d 0)).order := by
    calc ((d 0 : ℕ) : ℕ∞) = (d 0) • (1 : ℕ∞) := by simp
      _ ≤ (d 0) • ((swapFam (A := A)) 0).order := by gcongr; exact swapFam_order 0
      _ ≤ _ := MvPowerSeries.le_order_pow (d 0)
  have h1 : ((d 1 : ℕ) : ℕ∞) ≤ (((swapFam (A := A)) 1) ^ (d 1)).order := by
    calc ((d 1 : ℕ) : ℕ∞) = (d 1) • (1 : ℕ∞) := by simp
      _ ≤ (d 1) • ((swapFam (A := A)) 1).order := by gcongr; exact swapFam_order 1
      _ ≤ _ := MvPowerSeries.le_order_pow (d 1)
  calc ((d 0 + d 1 : ℕ) : ℕ∞) = ((d 0 : ℕ) : ℕ∞) + ((d 1 : ℕ) : ℕ∞) := by push_cast; ring
    _ ≤ (((swapFam (A := A)) 0) ^ (d 0)).order + (((swapFam (A := A)) 1) ^ (d 1)).order :=
          add_le_add h0 h1
    _ ≤ _ := MvPowerSeries.le_order_mul

/-- `Fin 2 →₀ ℕ` の次数1の元は `single 0 1` か `single 1 1` のどちらか。 -/
theorem degree_one_cases {d : Fin 2 →₀ ℕ} (h : Finsupp.degree d = 1) :
    d = Finsupp.single (0 : Fin 2) 1 ∨ d = Finsupp.single (1 : Fin 2) 1 := by
  rw [finsupp_degree_fin2] at h
  by_cases hd0 : d 0 = 0
  · right
    have hd1 : d 1 = 1 := by omega
    ext j; fin_cases j
    · simp [hd0]
    · simp [hd1]
  · left
    have hd0' : d 0 = 1 := by omega
    have hd1 : d 1 = 0 := by omega
    ext j; fin_cases j
    · simp [hd0']
    · simp [hd1]

theorem single0_ne_single1 : (Finsupp.single (0 : Fin 2) 1 : Fin 2 →₀ ℕ) ≠ Finsupp.single (1 : Fin 2) 1 := by
  intro heq
  have h0 : (Finsupp.single (0 : Fin 2) 1 : Fin 2 →₀ ℕ) 0 = (Finsupp.single (1 : Fin 2) 1 : Fin 2 →₀ ℕ) 0 := by
    rw [heq]
  rw [Finsupp.single_eq_same, Finsupp.single_eq_of_ne' (by decide : (1 : Fin 2) ≠ 0)] at h0
  exact one_ne_zero h0

theorem single_ne_zero (i : Fin 2) : (Finsupp.single i 1 : Fin 2 →₀ ℕ) ≠ 0 := by
  intro heq
  have h0 : (Finsupp.single i 1 : Fin 2 →₀ ℕ) i = (0 : Fin 2 →₀ ℕ) i := by rw [heq]
  rw [Finsupp.single_eq_same, Finsupp.coe_zero, Pi.zero_apply] at h0
  exact one_ne_zero h0

/-- `d` の次数が `≠1` のとき、`d.prod(swapFam^·)` の次数1の係数は0
——次数0(`d=0`)なら定数項1で狙う単項式には効かず、次数`≥2`なら
`order_prod_swapFam_ge` により次数の議論で落ちる。 -/
theorem coeff_single_prod_swapFam_eq_zero (i : Fin 2) {d : Fin 2 →₀ ℕ} (hdeg : Finsupp.degree d ≠ 1) :
    MvPowerSeries.coeff (Finsupp.single i 1) (d.prod (fun s m => ((swapFam (A := A)) s) ^ m)) = 0 := by
  rcases eq_or_ne (Finsupp.degree d) 0 with hd0 | hd0
  · have hde0 : d = 0 := by
      rw [finsupp_degree_fin2] at hd0
      ext j; fin_cases j <;> simp <;> omega
    rw [hde0]
    show MvPowerSeries.coeff (Finsupp.single i 1) (1 : MvPowerSeries (Fin 2) A) = 0
    rw [MvPowerSeries.coeff_one, if_neg (single_ne_zero i)]
  · apply MvPowerSeries.coeff_of_lt_order
    have hge2 : ((Finsupp.degree d : ℕ) : ℕ∞) ≥ 2 := by
      have : Finsupp.degree d ≥ 2 := by omega
      exact_mod_cast this
    have h1 : ((Finsupp.degree (Finsupp.single i 1) : ℕ) : ℕ∞) = 1 := by
      simp
    calc ((Finsupp.degree (Finsupp.single i 1) : ℕ) : ℕ∞) = 1 := h1
      _ < 2 := by norm_num
      _ ≤ ((Finsupp.degree d : ℕ) : ℕ∞) := hge2
      _ ≤ _ := order_prod_swapFam_ge d

/-- `swap Φ` の `X_0` 係数は `Φ` の `X_1` 係数に等しい。 -/
theorem coeff_single0_swap (Φ : MvPowerSeries (Fin 2) A) :
    MvPowerSeries.coeff (Finsupp.single (0 : Fin 2) 1) (swap Φ) =
      MvPowerSeries.coeff (Finsupp.single (1 : Fin 2) 1) Φ := by
  show MvPowerSeries.coeff (Finsupp.single (0 : Fin 2) 1)
    (MvPowerSeries.subst (swapFam (A := A)) Φ) = _
  rw [MvPowerSeries.coeff_subst hasSubst_swapFam]
  rw [finsum_eq_sum_of_support_subset _ (s := ({Finsupp.single (1 : Fin 2) 1} : Finset (Fin 2 →₀ ℕ)))
    (fun d hd => by
      simp only [Function.mem_support] at hd
      simp only [Finset.coe_singleton, Set.mem_singleton_iff]
      by_contra hcon
      apply hd
      by_cases hdeg : Finsupp.degree d = 1
      · rcases degree_one_cases hdeg with hd01 | hd11
        · rw [hd01]
          have hprod01 : (Finsupp.single (0 : Fin 2) 1).prod (fun s m => ((swapFam (A := A)) s) ^ m) =
              (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) := by
            show (∏ j ∈ (Finsupp.single (0 : Fin 2) 1).support,
              ((swapFam (A := A)) j) ^ ((Finsupp.single (0 : Fin 2) 1) j)) = _
            rw [Finset.prod_subset (Finset.subset_univ _) (fun x _ hx => by
              simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two]
            simp [swapFam]
          rw [hprod01, MvPowerSeries.coeff_X, if_neg single0_ne_single1, smul_zero]
        · exact absurd hd11 hcon
      · rw [coeff_single_prod_swapFam_eq_zero 0 hdeg, smul_zero])]
  rw [Finset.sum_singleton]
  have hprod : (Finsupp.single (1 : Fin 2) 1).prod (fun s m => ((swapFam (A := A)) s) ^ m) =
      (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) := by
    show (∏ j ∈ (Finsupp.single (1 : Fin 2) 1).support,
      ((swapFam (A := A)) j) ^ ((Finsupp.single (1 : Fin 2) 1) j)) = _
    rw [Finset.prod_subset (Finset.subset_univ _) (fun x _ hx => by
      simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two]
    simp [swapFam]
  rw [hprod, MvPowerSeries.coeff_X, if_pos rfl, smul_eq_mul, mul_one]

/-- `swap Φ` の `X_1` 係数は `Φ` の `X_0` 係数に等しい。 -/
theorem coeff_single1_swap (Φ : MvPowerSeries (Fin 2) A) :
    MvPowerSeries.coeff (Finsupp.single (1 : Fin 2) 1) (swap Φ) =
      MvPowerSeries.coeff (Finsupp.single (0 : Fin 2) 1) Φ := by
  show MvPowerSeries.coeff (Finsupp.single (1 : Fin 2) 1)
    (MvPowerSeries.subst (swapFam (A := A)) Φ) = _
  rw [MvPowerSeries.coeff_subst hasSubst_swapFam]
  rw [finsum_eq_sum_of_support_subset _ (s := ({Finsupp.single (0 : Fin 2) 1} : Finset (Fin 2 →₀ ℕ)))
    (fun d hd => by
      simp only [Function.mem_support] at hd
      simp only [Finset.coe_singleton, Set.mem_singleton_iff]
      by_contra hcon
      apply hd
      by_cases hdeg : Finsupp.degree d = 1
      · rcases degree_one_cases hdeg with hd01 | hd11
        · exact absurd hd01 hcon
        · rw [hd11]
          have hprod11 : (Finsupp.single (1 : Fin 2) 1).prod (fun s m => ((swapFam (A := A)) s) ^ m) =
              (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) := by
            show (∏ j ∈ (Finsupp.single (1 : Fin 2) 1).support,
              ((swapFam (A := A)) j) ^ ((Finsupp.single (1 : Fin 2) 1) j)) = _
            rw [Finset.prod_subset (Finset.subset_univ _) (fun x _ hx => by
              simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two]
            simp [swapFam]
          rw [hprod11, MvPowerSeries.coeff_X, if_neg single0_ne_single1.symm, smul_zero]
      · rw [coeff_single_prod_swapFam_eq_zero 1 hdeg, smul_zero])]
  rw [Finset.sum_singleton]
  have hprod : (Finsupp.single (0 : Fin 2) 1).prod (fun s m => ((swapFam (A := A)) s) ^ m) =
      (MvPowerSeries.X 1 : MvPowerSeries (Fin 2) A) := by
    show (∏ j ∈ (Finsupp.single (0 : Fin 2) 1).support,
      ((swapFam (A := A)) j) ^ ((Finsupp.single (0 : Fin 2) 1) j)) = _
    rw [Finset.prod_subset (Finset.subset_univ _) (fun x _ hx => by
      simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two]
    simp [swapFam]
  rw [hprod, MvPowerSeries.coeff_X, if_pos rfl, smul_eq_mul, mul_one]

variable [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))

include hq hπmax hf0 hf1 hfres in
/-- ★★★★★★★★★★★**Lubin-Tate 形式群法則の可換律**: `F_f(X,Y)=F_f(Y,X)`。
`swap F_f` は `F_f` と同じ関数等式(`swap_preserves_functional_equation`)・
同じ定数項0(`constantCoeff_swap`)・同じ次数1の係数
(`coeff_single0_swap`・`coeff_single1_swap`、`F_f` の次数1部分が
`X_0+X_1` と対称だから)を満たすので、2変数の一意性補題
(`mvPowerSeries_uniqueness`)で `F_f = swap F_f` が出る。 -/
theorem formalGroupLaw_commutative (hπne0 : π ≠ 0) :
    formalGroupLaw hq hπmax f hf0 hf1 hfres = swap (formalGroupLaw hq hπmax f hf0 hf1 hfres) := by
  have hΦ0 : MvPowerSeries.constantCoeff (formalGroupLaw hq hπmax f hf0 hf1 hfres) = 0 :=
    constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres
  have heqΦ := formalGroupLaw_f_isEndomorphism hq hπmax f hf0 hf1 hfres
  have heqΨ := swap_preserves_functional_equation hf0 hΦ0 heqΦ
  apply mvPowerSeries_uniqueness hπmax hπne0 hf0 hf1 hΦ0 (constantCoeff_swap hΦ0)
  · intro i
    fin_cases i
    · show MvPowerSeries.coeff (Finsupp.single (0 : Fin 2) 1) (formalGroupLaw hq hπmax f hf0 hf1 hfres) =
        MvPowerSeries.coeff (Finsupp.single (0 : Fin 2) 1) (swap (formalGroupLaw hq hπmax f hf0 hf1 hfres))
      rw [coeff_single0_swap, coeff_single01_formalGroupLaw, coeff_single10_formalGroupLaw]
    · show MvPowerSeries.coeff (Finsupp.single (1 : Fin 2) 1) (formalGroupLaw hq hπmax f hf0 hf1 hfres) =
        MvPowerSeries.coeff (Finsupp.single (1 : Fin 2) 1) (swap (formalGroupLaw hq hπmax f hf0 hf1 hfres))
      rw [coeff_single1_swap, coeff_single10_formalGroupLaw, coeff_single01_formalGroupLaw]
  · exact heqΦ
  · exact heqΨ

end ABC3.Found.PGC

