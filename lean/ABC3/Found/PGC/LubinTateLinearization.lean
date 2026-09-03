import ABC3.Found.PGC.LubinTateStep

/-!
# Lubin-Tate 形式群の存在補題: 線形化(`sorry` 無し)

`Found/PGC/LubinTateStep.lean` で「見当をつけたが未着手」としていた**線形化**
——`f` 側の代入 `f(Φ+φ)` を `f(Φ) + π·φ` で近似したときの誤差が、`φ` の次数が
高いことにより高次(`φ` の次数以上)に押し出されること——を、一般の可換環
`A` について sorry 無しで確立する。

## 数学的な内容

`φ` の次数が(下から)`m ≥ 1` 以上、`Φ` の次数が `1` 以上のとき、`k ≥ 2` について

```
(Φ+φ)^k − Φ^k = φ · (Σ_{i<k} (Φ+φ)^i Φ^{k-1-i})    (`geom_sum₂_mul`)
```

の右辺の和は各項が次数 `≥ k-1`(各因子が次数 `≥1` の `k-1` 個の積)なので、
全体の次数は `≥ m + (k-1) ≥ m+1`(`k≥2` だから)。`f` の `k≥2` の項がこの
`(Φ+φ)^k−Φ^k` を通じてしか現れないので、`f(Φ+φ) − f(Φ) − π·φ` の次数 `≤ m`
の係数はすべて 0 になる——`coeff_subst_linearize` がこれを直接与える。

## 使った mathlib の道具(実測、2026-09-04)

- `MvPowerSeries.order`・`le_order_pow`・`le_order_mul`・`min_order_le_add`
  (`RingTheory/MvPowerSeries/Order.lean`)。
- `geom_sum₂_mul`(`x^n-y^n = (x-y)·Σxⁱyⁿ⁻¹⁻ⁱ`、`Found/PGC/PadicLogSurjective.lean`
  で p 進対数の全射性にも使った古典的因数分解)。
- `PowerSeries.coeff_subst`(finsum 表示)+ `finsum_eq_sum_of_support_subset`
  (有限台への還元、`Found/PGC/LubinTateDivisibility.lean::residue_divides_R`
  と同じ技法)。

## 続報(2026-09-04): g 側の線形化も完成した

`g` 側(`Φ.subst(g,g)` の変化分 `φ.subst(g,g)`)についても、`φ` が**ちょうど**
次数 `n+1` の斉次式であるとき `φ.subst(g,g) ≡ π^{n+1}·φ`(次数 `n+2` 以上を除いて)
となることが、`Found/PGC/LubinTateGLinearization.lean::coeff_subst_g_linearize`
として sorry 無しで示された。これで `Found/PGC/LubinTateStep.lean` の
`exists_step_solution_for_R`(障害方程式の可解性)と合わせて、次数ごとの
帰納法の1ステップに要る**3つの部品**(f側の線形化・g側の線形化・可解性)が
すべて揃った。残るのはこれらを組み合わせて実際に1ステップの帰納法の
補題を組み立て、`Nat.rec` で無限次まで構成する純粋な組み立て作業のみ。
-/

namespace ABC3.Found.PGC

/-! ### 部品1: 有限和の次数は各項の次数の最小値以上 -/

theorem le_order_sum {σ R ι : Type*} [Semiring R] (s : Finset ι) (f : ι → MvPowerSeries σ R)
    {n : ℕ∞} (h : ∀ i ∈ s, n ≤ (f i).order) : n ≤ (∑ i ∈ s, f i).order := by
  induction s using Finset.cons_induction with
  | empty => simp [MvPowerSeries.order_zero]
  | cons a s ha ih =>
    rw [Finset.sum_cons]
    refine le_trans ?_ MvPowerSeries.min_order_le_add
    simp only [le_min_iff]
    exact ⟨h a (Finset.mem_cons_self a s), ih (fun i hi => h i (Finset.mem_cons_of_mem hi))⟩

/-! ### ★★部品2: `(Φ+φ)^k − Φ^k` の次数下界(`k ≥ 2`、`geom_sum₂_mul` 経由) -/

/-- ★★**`(Φ+φ)^k − Φ^k` の次数は `m+1` 以上**(`Φ` の次数 `≥1`・`φ` の次数
`≥m≥1`・`k≥2`)。`geom_sum₂_mul` の因数分解 `(Φ+φ)^k−Φ^k=φ·Σᵢ(Φ+φ)ⁱΦᵏ⁻¹⁻ⁱ`
で、和の各項の次数が `≥k-1`(`k-1` 個の次数 `≥1` の因子の積)であることから
`φ` の次数 `m` と合わせて `m+(k-1)≥m+1` を得る。 -/
theorem order_pow_sub_pow_ge {A σ : Type*} [CommRing A] {Φ φ : MvPowerSeries σ A} {m : ℕ}
    (hΦ : 1 ≤ Φ.order) (hφ : (m : ℕ∞) ≤ φ.order) (hm : 1 ≤ m) (k : ℕ) (hk : 2 ≤ k) :
    ((m : ℕ∞) + 1) ≤ ((Φ + φ) ^ k - Φ ^ k).order := by
  have hfact : (∑ i ∈ Finset.range k, (Φ + φ) ^ i * Φ ^ (k - 1 - i)) * ((Φ + φ) - Φ)
      = (Φ + φ) ^ k - Φ ^ k := geom_sum₂_mul (Φ + φ) Φ k
  have hΦφ : 1 ≤ (Φ + φ).order :=
    le_trans (le_min hΦ (le_trans (by exact_mod_cast hm) hφ)) MvPowerSeries.min_order_le_add
  have hsum : ((k - 1 : ℕ) : ℕ∞) ≤
      (∑ i ∈ Finset.range k, (Φ + φ) ^ i * Φ ^ (k - 1 - i)).order := by
    apply le_order_sum
    intro i hi
    rw [Finset.mem_range] at hi
    have h1 : (i : ℕ∞) ≤ ((Φ + φ) ^ i).order := by
      calc (i : ℕ∞) = i • (1 : ℕ∞) := by simp
        _ ≤ i • (Φ + φ).order := by gcongr
        _ ≤ ((Φ + φ) ^ i).order := MvPowerSeries.le_order_pow i
    have h2 : ((k - 1 - i : ℕ) : ℕ∞) ≤ (Φ ^ (k - 1 - i)).order := by
      calc ((k - 1 - i : ℕ) : ℕ∞) = (k - 1 - i) • (1 : ℕ∞) := by simp
        _ ≤ (k - 1 - i) • Φ.order := by gcongr
        _ ≤ (Φ ^ (k - 1 - i)).order := MvPowerSeries.le_order_pow (k - 1 - i)
    have heq : ((k - 1 : ℕ) : ℕ∞) = (i : ℕ∞) + ((k - 1 - i : ℕ) : ℕ∞) := by
      have hnat : k - 1 = i + (k - 1 - i) := by omega
      exact_mod_cast hnat
    rw [heq]
    calc (i : ℕ∞) + ((k - 1 - i : ℕ) : ℕ∞) ≤ ((Φ + φ) ^ i).order + (Φ ^ (k - 1 - i)).order :=
          add_le_add h1 h2
      _ ≤ ((Φ + φ) ^ i * Φ ^ (k - 1 - i)).order := MvPowerSeries.le_order_mul
  have hdiff : (Φ + φ) - Φ = φ := by ring
  rw [hdiff] at hfact
  have hfinal : ((m : ℕ∞) + 1) ≤ ((k - 1 : ℕ) : ℕ∞) + (m : ℕ∞) := by
    have h1lek1 : (1 : ℕ∞) ≤ ((k - 1 : ℕ) : ℕ∞) := by
      have : 1 ≤ k - 1 := by omega
      exact_mod_cast this
    calc ((m : ℕ∞) + 1) = (m : ℕ∞) + 1 := by ring
      _ ≤ (m : ℕ∞) + ((k - 1 : ℕ) : ℕ∞) := by gcongr
      _ = ((k - 1 : ℕ) : ℕ∞) + (m : ℕ∞) := by ring
  calc ((m : ℕ∞) + 1) ≤ ((k - 1 : ℕ) : ℕ∞) + (m : ℕ∞) := hfinal
    _ ≤ (∑ i ∈ Finset.range k, (Φ + φ) ^ i * Φ ^ (k - 1 - i)).order + φ.order :=
        add_le_add hsum hφ
    _ ≤ ((∑ i ∈ Finset.range k, (Φ + φ) ^ i * Φ ^ (k - 1 - i)) * φ).order :=
        MvPowerSeries.le_order_mul
    _ = ((Φ + φ) ^ k - Φ ^ k).order := by rw [hfact]

/-! ### ★★★線形化そのもの -/

/-- ★★★**線形化**: `f ≡ πX (mod deg 2)` のとき、`f(Φ+φ) − f(Φ)` は次数
`≤ deg φ` の範囲で `π·φ` にちょうど一致する。`f` の `k≥2` 次の項は
`order_pow_sub_pow_ge` によって次数 `deg φ + 1` 以上にしか寄与しないので、
`coeff_subst` の有限和から `k=0,1` の項だけが残る。 -/
theorem coeff_subst_linearize {A : Type*} [CommRing A] {Φ φ : MvPowerSeries (Fin 2) A}
    (hΦ0 : MvPowerSeries.constantCoeff Φ = 0) (hφ0 : MvPowerSeries.constantCoeff φ = 0)
    (hΦ : 1 ≤ Φ.order) {m : ℕ} (hφord : (m : ℕ∞) ≤ φ.order) (hm : 1 ≤ m)
    (f : PowerSeries A) (π : A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (e : Fin 2 →₀ ℕ) (he : ((Finsupp.degree e : ℕ) : ℕ∞) ≤ m) :
    MvPowerSeries.coeff e (PowerSeries.subst (Φ + φ) f) -
        MvPowerSeries.coeff e (PowerSeries.subst Φ f) = π * MvPowerSeries.coeff e φ := by
  have hΦHS : PowerSeries.HasSubst Φ := by
    show IsNilpotent (MvPowerSeries.constantCoeff Φ); rw [hΦ0]; exact IsNilpotent.zero
  have hΦφHS : PowerSeries.HasSubst (Φ + φ) := by
    show IsNilpotent (MvPowerSeries.constantCoeff (Φ + φ))
    rw [map_add, hΦ0, hφ0, add_zero]; exact IsNilpotent.zero
  rw [PowerSeries.coeff_subst hΦφHS, PowerSeries.coeff_subst hΦHS]
  set N := Finsupp.degree e + 2 with hN
  have hΦbound : ∀ d, N ≤ d → MvPowerSeries.coeff e (Φ ^ d) = 0 := by
    intro d hd
    apply MvPowerSeries.coeff_of_lt_order
    calc ((Finsupp.degree e : ℕ) : ℕ∞) < N := by
          rw [hN]; exact_mod_cast (by omega : Finsupp.degree e < Finsupp.degree e + 2)
      _ ≤ (d : ℕ∞) := by exact_mod_cast hd
      _ ≤ d • Φ.order := by
            calc (d : ℕ∞) = d • (1 : ℕ∞) := by simp
              _ ≤ d • Φ.order := by gcongr
      _ ≤ (Φ ^ d).order := MvPowerSeries.le_order_pow d
  have hΦφbound : ∀ d, N ≤ d → MvPowerSeries.coeff e ((Φ + φ) ^ d) = 0 := by
    intro d hd
    have hΦφ1 : 1 ≤ (Φ + φ).order :=
      le_trans (le_min hΦ (le_trans (by exact_mod_cast hm) hφord)) MvPowerSeries.min_order_le_add
    apply MvPowerSeries.coeff_of_lt_order
    calc ((Finsupp.degree e : ℕ) : ℕ∞) < N := by
          rw [hN]; exact_mod_cast (by omega : Finsupp.degree e < Finsupp.degree e + 2)
      _ ≤ (d : ℕ∞) := by exact_mod_cast hd
      _ ≤ d • (Φ + φ).order := by
            calc (d : ℕ∞) = d • (1 : ℕ∞) := by simp
              _ ≤ d • (Φ + φ).order := by gcongr
      _ ≤ ((Φ + φ) ^ d).order := MvPowerSeries.le_order_pow d
  rw [finsum_eq_sum_of_support_subset _ (s := Finset.range N) (fun d hd => by
        simp only [Function.mem_support] at hd
        simp only [Finset.coe_range, Set.mem_Iio]
        by_contra hcon
        exact hd (by simp [hΦφbound d (by omega)])),
      finsum_eq_sum_of_support_subset _ (s := Finset.range N) (fun d hd => by
        simp only [Function.mem_support] at hd
        simp only [Finset.coe_range, Set.mem_Iio]
        by_contra hcon
        exact hd (by simp [hΦbound d (by omega)]))]
  rw [← Finset.sum_sub_distrib]
  have hterms : ∀ d ∈ Finset.range N,
      (PowerSeries.coeff d) f • MvPowerSeries.coeff e ((Φ + φ) ^ d) -
        (PowerSeries.coeff d) f • MvPowerSeries.coeff e (Φ ^ d) =
      if d = 1 then π * MvPowerSeries.coeff e φ else 0 := by
    intro d _
    rw [smul_eq_mul, smul_eq_mul, ← mul_sub]
    rcases d with _ | _ | d
    · simp [hf0]
    · simp [hf1, pow_one]
    · have hdge2 : 2 ≤ d + 2 := by omega
      have hzero : MvPowerSeries.coeff e ((Φ + φ) ^ (d + 2)) -
          MvPowerSeries.coeff e (Φ ^ (d + 2)) = 0 := by
        have hordge : ((m : ℕ∞) + 1) ≤ ((Φ + φ) ^ (d + 2) - Φ ^ (d + 2)).order :=
          order_pow_sub_pow_ge hΦ hφord hm (d + 2) hdge2
        have hlt : ((Finsupp.degree e : ℕ) : ℕ∞) < ((Φ + φ) ^ (d + 2) - Φ ^ (d + 2)).order := by
          have hstep : (m : ℕ∞) < (m : ℕ∞) + 1 := by
            have : ((m : ℕ) : ℕ∞) < ((m + 1 : ℕ) : ℕ∞) := by exact_mod_cast Nat.lt_succ_self m
            simpa using this
          calc ((Finsupp.degree e : ℕ) : ℕ∞) ≤ (m : ℕ∞) := he
            _ < (m : ℕ∞) + 1 := hstep
            _ ≤ _ := hordge
        have := MvPowerSeries.coeff_of_lt_order hlt
        rw [map_sub] at this
        exact this
      simp [hzero]
  rw [Finset.sum_congr rfl hterms, Finset.sum_ite_eq' (Finset.range N) 1
    (fun _ => π * MvPowerSeries.coeff e φ)]
  have h1mem : 1 ∈ Finset.range N := by rw [Finset.mem_range, hN]; omega
  simp [h1mem]

end ABC3.Found.PGC
