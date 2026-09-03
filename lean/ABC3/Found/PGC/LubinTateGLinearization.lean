import ABC3.Found.PGC.LubinTateHomogeneousScaling
import ABC3.Found.PGC.LubinTateDivisibility
import ABC3.Found.PGC.LubinTateLinearization

/-!
# Lubin-Tate 形式群の存在補題: g 側の線形化(`sorry` 無し)

f 側の線形化(`Found/PGC/LubinTateLinearization.lean`)・g 側の準備
(斉次式のスケーリング、`Found/PGC/LubinTateHomogeneousScaling.lean`)に続き、
g 側の線形化そのもの——`φ` が次数 `n+1` の斉次式のとき、`φ.subst(g,g)` が
`π^{n+1}•φ` と次数 `n+1` の範囲で一致すること——を確立する。これで
次数ごとの帰納法の1ステップに要る2つの線形化(f 側・g 側)が両方揃った。

## 数学的な内容

`g ≡ πX (mod deg 2)` のとき、`a_i := g.subst(X_i)`(`Φ.subst(g,g)` の
i 番目のスロット)は `π•X_i` から次数 `≥2` だけずれている
(`coeff_a_diff_order`)。`φ` が次数 `n+1` の斉次式のとき、`φ.subst(a)` の
`coeff_subst` 展開に効くのは `degree d = n+1` の項だけで、そこでの
`a^d − a'^d`(`a' := π•X`)の次数は `order_prod_pow_sub_prod_pow_ge`
(2変数の telescoping、`(m=2, d0+d1=n+1)` を代入すると `≥ n+2`)によって
押さえられる——`degree e ≤ n+1` の範囲では効かない。`a'` 側は
`homogeneous_subst_const_smul` で厳密に `π^{n+1}•φ` に一致する。

## まだ無いもの

f 側・g 側の線形化を組み合わせて実際に「`Φ_n` の障害が次数 `n` まで消えて
いる」→「`Φ_{n+1}:=Φ_n+φ_{n+1}` の障害が次数 `n+1` まで消えている」という
帰納法の1ステップ全体を組み立てる段、そして `Nat.rec` による無限次までの
構成・極限の検証は、まだ残る。
-/

namespace ABC3.Found.PGC

/-- `Fin 2 →₀ ℕ` の全次数は `d 0 + d 1`。 -/
theorem finsupp_degree_fin2 (d : Fin 2 →₀ ℕ) : Finsupp.degree d = d 0 + d 1 := by
  unfold Finsupp.degree
  simp only [AddMonoidHom.coe_mk, ZeroHom.coe_mk]
  rw [Finset.sum_subset (Finset.subset_univ d.support)]
  · rw [Fin.sum_univ_two]
  · intro x _ hx
    simpa using hx

/-! ### 部品1: `g.subst(X_i)` は `π•X_i` から次数 `≥2` だけずれている -/

theorem coeff_a_diff_order {A : Type*} [CommRing A] (g : PowerSeries A) (π : A)
    (hg0 : PowerSeries.constantCoeff g = 0) (hg1 : PowerSeries.coeff 1 g = π) (i : Fin 2) :
    (2 : ℕ∞) ≤ (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g -
      π • (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)).order := by
  have hXi : PowerSeries.HasSubst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) := hasSubst_X_i i
  have heq : PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g -
      π • (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
      = PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) (g - π • PowerSeries.X) := by
    rw [PowerSeries.subst_sub hXi, PowerSeries.subst_smul hXi]
    have hXsubst :
        PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) (PowerSeries.X : PowerSeries A)
        = (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) := PowerSeries.subst_X hXi
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
    _ ≤ (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) (g - π • PowerSeries.X)).order :=
        PowerSeries.le_order_subst_left (MvPowerSeries.constantCoeff_X i)

/-! ### 部品1.5: `order_pow_sub_pow_ge` の一般化(`k ≥ 1` すべて、`order(a-a')` 任意) -/

/-- `Found/PGC/LubinTateLinearization.lean::order_pow_sub_pow_ge` の一般化版:
`k` に下限を課さず、`Φ` 側の次数が厳密に `1` であることも要求しない——
`order(a),order(a')≥1` だけから `order(a^k-a'^k) ≥ order(a-a')+(k-1)` が
`geom_sum₂_mul` の同じ因数分解でそのまま出る(`k=0` は両辺が自明に成立、
`a-a'` 側の次数を「ちょうど `1`」ではなく変数 `order(a-a')` のまま持ち回る
だけの違い)。 -/
theorem order_pow_sub_pow_ge' {A σ : Type*} [CommRing A] {a a' : MvPowerSeries σ A}
    (ha : 1 ≤ a.order) (ha' : 1 ≤ a'.order) (k : ℕ) :
    ((a - a').order + ((k - 1 : ℕ) : ℕ∞)) ≤ (a ^ k - a' ^ k).order := by
  have hfact : (∑ i ∈ Finset.range k, a ^ i * a' ^ (k - 1 - i)) * (a - a') = a ^ k - a' ^ k :=
    geom_sum₂_mul a a' k
  have hsum : ((k - 1 : ℕ) : ℕ∞) ≤
      (∑ i ∈ Finset.range k, a ^ i * a' ^ (k - 1 - i)).order := by
    apply le_order_sum
    intro i hi
    rw [Finset.mem_range] at hi
    have h1 : (i : ℕ∞) ≤ (a ^ i).order := by
      calc (i : ℕ∞) = i • (1 : ℕ∞) := by simp
        _ ≤ i • a.order := by gcongr
        _ ≤ (a ^ i).order := MvPowerSeries.le_order_pow i
    have h2 : ((k - 1 - i : ℕ) : ℕ∞) ≤ (a' ^ (k - 1 - i)).order := by
      calc ((k - 1 - i : ℕ) : ℕ∞) = (k - 1 - i) • (1 : ℕ∞) := by simp
        _ ≤ (k - 1 - i) • a'.order := by gcongr
        _ ≤ (a' ^ (k - 1 - i)).order := MvPowerSeries.le_order_pow (k - 1 - i)
    have heq : ((k - 1 : ℕ) : ℕ∞) = (i : ℕ∞) + ((k - 1 - i : ℕ) : ℕ∞) := by
      have hnat : k - 1 = i + (k - 1 - i) := by omega
      exact_mod_cast hnat
    rw [heq]
    calc (i : ℕ∞) + ((k - 1 - i : ℕ) : ℕ∞) ≤ (a ^ i).order + (a' ^ (k - 1 - i)).order :=
          add_le_add h1 h2
      _ ≤ (a ^ i * a' ^ (k - 1 - i)).order := MvPowerSeries.le_order_mul
  calc (a - a').order + ((k - 1 : ℕ) : ℕ∞)
      = ((k - 1 : ℕ) : ℕ∞) + (a - a').order := by ring
    _ ≤ (∑ i ∈ Finset.range k, a ^ i * a' ^ (k - 1 - i)).order + (a - a').order :=
        add_le_add hsum le_rfl
    _ ≤ ((∑ i ∈ Finset.range k, a ^ i * a' ^ (k - 1 - i)) * (a - a')).order :=
        MvPowerSeries.le_order_mul
    _ = (a ^ k - a' ^ k).order := by rw [hfact]

/-! ### 部品1.6: 2変数の telescoping(`a_0^{k_0}a_1^{k_1} − a_0'^{k_0}a_1'^{k_1}`) -/

/-- `a_0^{k_0}a_1^{k_1} − a_0'^{k_0}a_1'^{k_1} = (a_0^{k_0}−a_0'^{k_0})·a_1^{k_1}
+ a_0'^{k_0}·(a_1^{k_1}−a_1'^{k_1})` の各項に `order_pow_sub_pow_ge'` を当てて、
`order(a_i−a_i')≥2` から全体の次数が `≥ 2+(k_0+k_1-1)` であることを得る。 -/
theorem order_prod_pow_sub_prod_pow_ge {A : Type*} [CommRing A]
    {a a' : Fin 2 → MvPowerSeries (Fin 2) A}
    (ha : ∀ i, 1 ≤ (a i).order) (ha' : ∀ i, 1 ≤ (a' i).order)
    (hdiff : ∀ i, (2 : ℕ∞) ≤ (a i - a' i).order) (k0 k1 : ℕ) :
    ((2 : ℕ∞) + ((k0 + k1 - 1 : ℕ) : ℕ∞)) ≤
      ((a 0) ^ k0 * (a 1) ^ k1 - (a' 0) ^ k0 * (a' 1) ^ k1).order := by
  have htelescope : (a 0) ^ k0 * (a 1) ^ k1 - (a' 0) ^ k0 * (a' 1) ^ k1 =
      ((a 0) ^ k0 - (a' 0) ^ k0) * (a 1) ^ k1 + (a' 0) ^ k0 * ((a 1) ^ k1 - (a' 1) ^ k1) := by
    ring
  have hpow0 : (a 0 - a' 0).order + ((k0 - 1 : ℕ) : ℕ∞) ≤ ((a 0) ^ k0 - (a' 0) ^ k0).order :=
    order_pow_sub_pow_ge' (ha 0) (ha' 0) k0
  have hpow1 : (a 1 - a' 1).order + ((k1 - 1 : ℕ) : ℕ∞) ≤ ((a 1) ^ k1 - (a' 1) ^ k1).order :=
    order_pow_sub_pow_ge' (ha 1) (ha' 1) k1
  have ha1pow : (k1 : ℕ∞) ≤ ((a 1) ^ k1).order := by
    calc (k1 : ℕ∞) = k1 • (1 : ℕ∞) := by simp
      _ ≤ k1 • (a 1).order := by gcongr; exact ha 1
      _ ≤ ((a 1) ^ k1).order := MvPowerSeries.le_order_pow k1
  have ha0'pow : (k0 : ℕ∞) ≤ ((a' 0) ^ k0).order := by
    calc (k0 : ℕ∞) = k0 • (1 : ℕ∞) := by simp
      _ ≤ k0 • (a' 0).order := by gcongr; exact ha' 0
      _ ≤ ((a' 0) ^ k0).order := MvPowerSeries.le_order_pow k0
  have hterm1 : ((2 : ℕ∞) + ((k0 + k1 - 1 : ℕ) : ℕ∞)) ≤
      (((a 0) ^ k0 - (a' 0) ^ k0) * (a 1) ^ k1).order := by
    have hb : ((2 : ℕ∞) + ((k0 - 1 : ℕ) : ℕ∞)) ≤ ((a 0) ^ k0 - (a' 0) ^ k0).order :=
      le_trans (add_le_add (hdiff 0) le_rfl) hpow0
    have hnat : (k0 + k1 - 1 : ℕ) ≤ (k0 - 1 : ℕ) + k1 := by omega
    calc ((2 : ℕ∞) + ((k0 + k1 - 1 : ℕ) : ℕ∞))
        ≤ (2 : ℕ∞) + (((k0 - 1 : ℕ) + k1 : ℕ) : ℕ∞) := by
          gcongr
      _ = ((2 : ℕ∞) + ((k0 - 1 : ℕ) : ℕ∞)) + (k1 : ℕ∞) := by push_cast; ring
      _ ≤ ((a 0) ^ k0 - (a' 0) ^ k0).order + ((a 1) ^ k1).order := add_le_add hb ha1pow
      _ ≤ (((a 0) ^ k0 - (a' 0) ^ k0) * (a 1) ^ k1).order := MvPowerSeries.le_order_mul
  have hterm2 : ((2 : ℕ∞) + ((k0 + k1 - 1 : ℕ) : ℕ∞)) ≤
      ((a' 0) ^ k0 * ((a 1) ^ k1 - (a' 1) ^ k1)).order := by
    have hb : ((2 : ℕ∞) + ((k1 - 1 : ℕ) : ℕ∞)) ≤ ((a 1) ^ k1 - (a' 1) ^ k1).order :=
      le_trans (add_le_add (hdiff 1) le_rfl) hpow1
    have hnat : (k0 + k1 - 1 : ℕ) ≤ k0 + (k1 - 1 : ℕ) := by omega
    calc ((2 : ℕ∞) + ((k0 + k1 - 1 : ℕ) : ℕ∞))
        ≤ (2 : ℕ∞) + ((k0 + (k1 - 1 : ℕ) : ℕ) : ℕ∞) := by
          gcongr
      _ = (k0 : ℕ∞) + ((2 : ℕ∞) + ((k1 - 1 : ℕ) : ℕ∞)) := by push_cast; ring
      _ ≤ ((a' 0) ^ k0).order + ((a 1) ^ k1 - (a' 1) ^ k1).order := add_le_add ha0'pow hb
      _ ≤ ((a' 0) ^ k0 * ((a 1) ^ k1 - (a' 1) ^ k1)).order := MvPowerSeries.le_order_mul
  rw [htelescope]
  exact le_trans (le_min hterm1 hterm2) MvPowerSeries.min_order_le_add

/-! ### 部品2: 次数 `n+1` の斉次式に対する `a^d` と `a'^d` の一致(`degree e ≤ n+1`) -/

theorem coeff_ad_eq_of_degree {A : Type*} [CommRing A] {n : ℕ}
    (a a' : Fin 2 → MvPowerSeries (Fin 2) A)
    (hai_order : ∀ i, 1 ≤ (a i).order) (hai'_order : ∀ i, 1 ≤ (a' i).order)
    (hdiff_order : ∀ i, (2 : ℕ∞) ≤ (a i - a' i).order)
    (e : Fin 2 →₀ ℕ) (he : Finsupp.degree e ≤ n + 1)
    (d : Fin 2 →₀ ℕ) (hd : Finsupp.degree d = n + 1) :
    MvPowerSeries.coeff e ((a 0) ^ (d 0) * (a 1) ^ (d 1)) =
      MvPowerSeries.coeff e ((a' 0) ^ (d 0) * (a' 1) ^ (d 1)) := by
  have hge : ((2 : ℕ∞) + ((d 0 + d 1 - 1 : ℕ) : ℕ∞)) ≤
      ((a 0) ^ (d 0) * (a 1) ^ (d 1) - (a' 0) ^ (d 0) * (a' 1) ^ (d 1)).order :=
    order_prod_pow_sub_prod_pow_ge hai_order hai'_order hdiff_order (d 0) (d 1)
  have hd01 : d 0 + d 1 - 1 = n := by rw [finsupp_degree_fin2] at hd; omega
  rw [hd01] at hge
  have hlt : ((Finsupp.degree e : ℕ) : ℕ∞) <
      ((a 0) ^ (d 0) * (a 1) ^ (d 1) - (a' 0) ^ (d 0) * (a' 1) ^ (d 1)).order := by
    calc ((Finsupp.degree e : ℕ) : ℕ∞) ≤ ((n + 1 : ℕ) : ℕ∞) := by exact_mod_cast he
      _ < ((2 : ℕ) : ℕ∞) + ((n : ℕ) : ℕ∞) := by
            have hlt2 : ((n + 1 : ℕ) : ℕ∞) < ((n + 2 : ℕ) : ℕ∞) :=
              by exact_mod_cast (by omega : n + 1 < n + 2)
            calc ((n + 1 : ℕ) : ℕ∞) < ((n + 2 : ℕ) : ℕ∞) := hlt2
              _ = ((2 : ℕ) : ℕ∞) + ((n : ℕ) : ℕ∞) := by push_cast; ring
      _ ≤ _ := hge
  have hz := MvPowerSeries.coeff_of_lt_order hlt
  rw [map_sub] at hz
  exact sub_eq_zero.mp hz

/-! ### ★★★g 側の線形化そのもの -/

/-- ★★★**g 側の線形化**: `g ≡ πX (mod deg 2)` のとき、次数 `n+1` の斉次式
`φ` について `φ.subst(g,g)` は `π^{n+1}•φ` と次数 `≤n+1` の範囲で一致する。
`homogeneous_subst_const_smul`(スケーリング)+ `coeff_ad_eq_of_degree`
(`a=g.subst(X)` と `a'=π•X` の一致)を `coeff_subst` の有限和に適用する。 -/
theorem coeff_subst_g_linearize {A : Type*} [CommRing A] {φ : MvPowerSeries (Fin 2) A} {n : ℕ}
    (hφ : ∀ d : Fin 2 →₀ ℕ, Finsupp.degree d ≠ n + 1 → MvPowerSeries.coeff d φ = 0)
    {g : PowerSeries A} (π : A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (e : Fin 2 →₀ ℕ) (he : Finsupp.degree e ≤ n + 1) :
    MvPowerSeries.coeff e (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) φ)
      = π ^ (n + 1) * MvPowerSeries.coeff e φ := by
  set a : Fin 2 → MvPowerSeries (Fin 2) A := fun i => PowerSeries.subst (MvPowerSeries.X i) g
    with ha_def
  set a' : Fin 2 → MvPowerSeries (Fin 2) A :=
    fun i => π • (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) with ha'_def
  have haHS : MvPowerSeries.HasSubst a := by
    constructor
    · intro i
      show IsNilpotent (MvPowerSeries.constantCoeff (a i))
      have heq0 : MvPowerSeries.constantCoeff (a i) = 0 :=
        PowerSeries.constantCoeff_subst_eq_zero (MvPowerSeries.constantCoeff_X i) g hg0
      rw [heq0]; exact IsNilpotent.zero
    · intro d; exact Set.toFinite _
  have hHSa' : MvPowerSeries.HasSubst a' := by
    constructor
    · intro i
      show IsNilpotent (MvPowerSeries.constantCoeff (a' i))
      show IsNilpotent (MvPowerSeries.constantCoeff
        (π • (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)))
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
    show MvPowerSeries.constantCoeff (π • (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)) = 0
    rw [MvPowerSeries.constantCoeff_smul, MvPowerSeries.constantCoeff_X, smul_zero]
  have hdiff_order : ∀ i, (2 : ℕ∞) ≤ (a i - a' i).order := fun i => coeff_a_diff_order g π hg0 hg1 i
  have hstep2 :
      MvPowerSeries.coeff e (MvPowerSeries.subst a φ) = MvPowerSeries.coeff e (MvPowerSeries.subst a' φ) := by
    rw [MvPowerSeries.coeff_subst haHS, MvPowerSeries.coeff_subst hHSa']
    refine finsum_congr (fun d => ?_)
    by_cases hd : Finsupp.degree d = n + 1
    · congr 1
      show MvPowerSeries.coeff e (d.prod fun s m => (a s) ^ m)
          = MvPowerSeries.coeff e (d.prod fun s m => (a' s) ^ m)
      have hdprod : ∀ (b : Fin 2 → MvPowerSeries (Fin 2) A),
          d.prod (fun s m => (b s) ^ m) = (b 0) ^ (d 0) * (b 1) ^ (d 1) := by
        intro b
        show (∏ i ∈ d.support, (b i) ^ (d i)) = (b 0) ^ (d 0) * (b 1) ^ (d 1)
        rw [Finset.prod_subset (Finset.subset_univ d.support) (fun x _ hx => by
          simp only [Finsupp.mem_support_iff, not_not] at hx; simp [hx]), Fin.prod_univ_two]
      rw [hdprod a, hdprod a']
      exact coeff_ad_eq_of_degree a a' hai_order hai'_order hdiff_order e he d hd
    · rw [hφ d hd]; simp
  have hstep1 : MvPowerSeries.coeff e (MvPowerSeries.subst a' φ) = π ^ (n + 1) * MvPowerSeries.coeff e φ := by
    have hscale := homogeneous_subst_const_smul hφ π
    rw [ha'_def]
    show MvPowerSeries.coeff e (MvPowerSeries.subst
      (fun i => π • (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)) φ) = _
    rw [hscale, MvPowerSeries.coeff_smul]
  rw [ha_def, hstep2, hstep1]

end ABC3.Found.PGC
