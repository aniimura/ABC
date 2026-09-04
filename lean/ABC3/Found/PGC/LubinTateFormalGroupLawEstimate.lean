import ABC3.Found.PGC.LubinTateIdentityLawSymm

/-!
# `F_f(z,w)-z-w` の評価レベルでのノルム評価(`sorry` 無し)

`Found/PGC/LubinTateIdentityLaw.lean::aeval_formalGroupLaw_eq_of_snd_eq_zero`
(`Y=0`のときの厳密な等式`F_f(y₀,0)=y₀`)の**不等式版**——`Y`を`0`に
限定せず、`‖z‖,‖w‖≤1`である任意の評価点`(z,w)`について

  `‖F_f(z,w) - z - w‖ ≤ ‖z‖ * ‖w‖`

を示す。古典的なLubin-Tate理論の教科書的な補題(単位元則
`F_f(X,0)=X`・`F_f(0,Y)=Y`により`F_f`の次数`≤1`部分がちょうど
`X+Y`であることから、残る「次数`≥2`」部分の各単項式が`i≥1`かつ
`j≥1`を満たすことを使う)の形式化。

## 使い道

この評価は、`𝒪_K`加群の作用`a·x`の**単射性**
(`a·x=b·x⟹π^n∣(a-b)`)を、`F_f`の形式逆元・対数といった大掛かりな
構成を一切経由せず、既存の加法公式`lubinTateAction_add`
(`(a+b)·x=F_f(a·x,b·x)`)と組み合わせるだけで直接導くための鍵になる
——`c:=a-b`とおくと`a·x=b·x`という仮定から`F_f(b·x,c·x)=b·x`が出て、
この不等式を`(z,w):=(b·x,c·x)`に適用すれば`‖c·x‖≤‖b·x‖*‖c·x‖`、
`b·x`が捩れ点で`‖b·x‖<1`であることと合わせて`c·x=0`(`w≠0`なら
矛盾)が出る——`Found/PGC/LubinTateActionInjective.lean`で使う。

## 証明の構造

`aeval_formalGroupLaw_eq_of_snd_eq_zero`と同じtruncation-limit手法:
各truncation`trunc' A n F_f`(`n₀≥1∧n₁≥1`)を`(z,w)`で評価した値から
`z+w`を引いたものは、その多項式のsupportから2つの特別な単項式
(`X_0^1`・`X_1^1`、係数がちょうど`1`——`coeff_single0_formalGroupLaw`・
`coeff_single1_formalGroupLaw`)を`Finset.add_sum_erase`で2回取り除いた
残りの有限和に一致し、残りの各単項式は`i≥1∧j≥1`(でなければ係数`0`)
なので`‖coeff‖≤1`・`‖z‖^i≤‖z‖`・`‖w‖^j≤‖w‖`(`pow_le_of_le_one`)から
`‖z‖*‖w‖`で抑えられる(`IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg`)。
この不等式が全ての`n₀≥1∧n₁≥1`で一様に成り立つので、
`MvPowerSeries.continuous_aeval`+`le_of_tendsto`で極限(`n→∞`)へ持ち
上げる。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical MvPowerSeries.WithPiTopology

variable {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))
    [UniformSpace A] [IsUniformAddGroup A] [IsTopologicalSemiring A]

include hq hπmax hπne0 hf0 hf1 hfres in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**`F_f(z,w)-z-w`のノルムは`‖z‖*‖w‖`
以下**——`F_f`の単位元則(両側)から、`z`か`w`のどちらかが`0`のとき
`F_f(z,w)-z-w`は恒等的に`0`になる(残る「次数`≥2`」項はすべて
`i≥1∧j≥1`)ことを使った評価レベルの不等式。 -/
theorem norm_aeval_formalGroupLaw_sub_le
    {S : Type*} [NormedCommRing S] [NormMulClass S] [NormOneClass S] [IsUltrametricDist S]
    [Algebra A S] [CompleteSpace S] [IsLinearTopology S S] [ContinuousSMul A S]
    {y : Fin 2 → S} (hy : MvPowerSeries.HasEval y)
    (hy0 : ‖y 0‖ ≤ 1) (hy1 : ‖y 1‖ ≤ 1) (halg : ∀ c : A, ‖(algebraMap A S) c‖ ≤ 1) :
    ‖MvPowerSeries.aeval hy (formalGroupLaw hq hπmax f hf0 hf1 hfres) - y 0 - y 1‖ ≤ ‖y 0‖ * ‖y 1‖ := by
  set F := formalGroupLaw hq hπmax f hf0 hf1 hfres with hF
  have hab : (Finsupp.single (0 : Fin 2) 1 : Fin 2 →₀ ℕ) ≠ Finsupp.single (1 : Fin 2) 1 := by
    intro heq
    have h0 : (Finsupp.single (0 : Fin 2) 1 : Fin 2 →₀ ℕ) 0 = (Finsupp.single (1 : Fin 2) 1 : Fin 2 →₀ ℕ) 0 := by
      rw [heq]
    simp [Finsupp.single_eq_of_ne' (by decide : (1 : Fin 2) ≠ 0)] at h0
  have hnonneg : (0 : ℝ) ≤ ‖y 0‖ * ‖y 1‖ := mul_nonneg (norm_nonneg _) (norm_nonneg _)
  have htrunc : ∀ n : Fin 2 →₀ ℕ, 1 ≤ n 0 → 1 ≤ n 1 →
      ‖MvPowerSeries.aeval hy ((MvPowerSeries.trunc' A n F : MvPolynomial (Fin 2) A) :
          MvPowerSeries (Fin 2) A) - y 0 - y 1‖ ≤ ‖y 0‖ * ‖y 1‖ := by
    intro n hn0 hn1
    rw [MvPowerSeries.aeval_coe, MvPolynomial.aeval_def, MvPolynomial.eval₂_eq']
    set g : (Fin 2 →₀ ℕ) → S := fun d =>
      (algebraMap A S) (MvPolynomial.coeff d (MvPowerSeries.trunc' A n F)) * ∏ i, y i ^ d i with hg
    have ha : Finsupp.single (0 : Fin 2) 1 ∈ (MvPowerSeries.trunc' A n F).support := by
      rw [MvPolynomial.mem_support_iff, MvPowerSeries.coeff_trunc']
      rw [if_pos]
      · rw [coeff_single0_formalGroupLaw hq hπmax f hf0 hf1 hfres hπne0]
        norm_num
      · intro i; fin_cases i
        · simpa using hn0
        · simp
    have hb : Finsupp.single (1 : Fin 2) 1 ∈ (MvPowerSeries.trunc' A n F).support := by
      rw [MvPolynomial.mem_support_iff, MvPowerSeries.coeff_trunc']
      rw [if_pos]
      · rw [coeff_single1_formalGroupLaw hq hπmax f hf0 hf1 hfres hπne0]
        norm_num
      · intro i; fin_cases i
        · simp
        · simpa using hn1
    have hga : g (Finsupp.single (0 : Fin 2) 1) = y 0 := by
      simp only [hg]
      rw [MvPowerSeries.coeff_trunc', if_pos, coeff_single0_formalGroupLaw hq hπmax f hf0 hf1 hfres hπne0]
      · simp [Fin.prod_univ_two]
      · intro i; fin_cases i
        · simpa using hn0
        · simp
    have hgb : g (Finsupp.single (1 : Fin 2) 1) = y 1 := by
      simp only [hg]
      rw [MvPowerSeries.coeff_trunc', if_pos, coeff_single1_formalGroupLaw hq hπmax f hf0 hf1 hfres hπne0]
      · simp [Fin.prod_univ_two]
      · intro i; fin_cases i
        · simp
        · simpa using hn1
    have hb' : Finsupp.single (1 : Fin 2) 1 ∈
        (MvPowerSeries.trunc' A n F).support.erase (Finsupp.single (0 : Fin 2) 1) :=
      Finset.mem_erase.mpr ⟨fun h => hab h.symm, hb⟩
    have hsplit : ∑ d ∈ (MvPowerSeries.trunc' A n F).support, g d =
        y 0 + (y 1 + ∑ d ∈ ((MvPowerSeries.trunc' A n F).support.erase
          (Finsupp.single (0 : Fin 2) 1)).erase (Finsupp.single (1 : Fin 2) 1), g d) := by
      rw [← Finset.add_sum_erase _ g ha, ← Finset.add_sum_erase _ g hb', hga, hgb]
    have hkey : ∑ d ∈ (MvPowerSeries.trunc' A n F).support, g d - y 0 - y 1 =
        ∑ d ∈ ((MvPowerSeries.trunc' A n F).support.erase
          (Finsupp.single (0 : Fin 2) 1)).erase (Finsupp.single (1 : Fin 2) 1), g d := by
      rw [hsplit]; abel
    show ‖∑ d ∈ (MvPowerSeries.trunc' A n F).support, g d - y 0 - y 1‖ ≤ _
    rw [hkey]
    apply IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hnonneg
    intro d hd
    rw [Finset.mem_erase, Finset.mem_erase] at hd
    obtain ⟨hd1, hd0, hdmem⟩ := hd
    simp only [hg]
    by_cases hd00 : d 0 = 0
    · have hdeq : d = Finsupp.single (1 : Fin 2) (d 1) := by
        ext i; fin_cases i
        · simpa using hd00
        · simp
      have hd1ne1 : d 1 ≠ 1 := by
        intro heq; apply hd1; rw [hdeq, heq]
      have hcoeff0 : MvPolynomial.coeff d (MvPowerSeries.trunc' A n F) = 0 := by
        rw [MvPowerSeries.coeff_trunc']
        by_cases hled : d ≤ n
        · rw [if_pos hled, hdeq, coeff_single1_formalGroupLaw hq hπmax f hf0 hf1 hfres hπne0, if_neg hd1ne1]
        · rw [if_neg hled]
      rw [hcoeff0]
      simp only [map_zero, zero_mul, norm_zero]
      exact hnonneg
    · by_cases hd10 : d 1 = 0
      · have hdeq : d = Finsupp.single (0 : Fin 2) (d 0) := by
          ext i; fin_cases i
          · simp
          · simpa using hd10
        have hd0ne1 : d 0 ≠ 1 := by
          intro heq; apply hd0; rw [hdeq, heq]
        have hcoeff0 : MvPolynomial.coeff d (MvPowerSeries.trunc' A n F) = 0 := by
          rw [MvPowerSeries.coeff_trunc']
          by_cases hled : d ≤ n
          · rw [if_pos hled, hdeq, coeff_single0_formalGroupLaw hq hπmax f hf0 hf1 hfres hπne0, if_neg hd0ne1]
          · rw [if_neg hled]
        rw [hcoeff0]
        simp only [map_zero, zero_mul, norm_zero]
        exact hnonneg
      · have h1 : ‖(algebraMap A S) (MvPolynomial.coeff d (MvPowerSeries.trunc' A n F))‖ ≤ 1 := halg _
        have hd0pos : d 0 ≠ 0 := hd00
        have hd1pos : d 1 ≠ 0 := hd10
        have hprod : ∏ i, y i ^ d i = y 0 ^ d 0 * y 1 ^ d 1 := by rw [Fin.prod_univ_two]
        rw [hprod, norm_mul, norm_mul, norm_pow, norm_pow]
        have hy0le : ‖y 0‖ ^ d 0 ≤ ‖y 0‖ := pow_le_of_le_one (norm_nonneg _) hy0 hd0pos
        have hy1le : ‖y 1‖ ^ d 1 ≤ ‖y 1‖ := pow_le_of_le_one (norm_nonneg _) hy1 hd1pos
        calc ‖(algebraMap A S) (MvPolynomial.coeff d (MvPowerSeries.trunc' A n F))‖ *
              (‖y 0‖ ^ d 0 * ‖y 1‖ ^ d 1)
            ≤ 1 * (‖y 0‖ * ‖y 1‖) := by
              apply mul_le_mul h1 (mul_le_mul hy0le hy1le (by positivity) (norm_nonneg _))
                (by positivity) zero_le_one
          _ = ‖y 0‖ * ‖y 1‖ := one_mul _
  have h1 : Filter.Tendsto (fun n => MvPowerSeries.aeval hy
      ((MvPowerSeries.trunc' A n F : MvPolynomial (Fin 2) A) : MvPowerSeries (Fin 2) A) - y 0 - y 1)
      Filter.atTop (nhds (MvPowerSeries.aeval hy F - y 0 - y 1)) := by
    have hcont : Filter.Tendsto (fun n => MvPowerSeries.aeval hy
        ((MvPowerSeries.trunc' A n F : MvPolynomial (Fin 2) A) : MvPowerSeries (Fin 2) A))
        Filter.atTop (nhds (MvPowerSeries.aeval hy F)) :=
      (MvPowerSeries.continuous_aeval hy).continuousAt.tendsto.comp
        (MvPowerSeries.WithPiTopology.tendsto_trunc'_atTop F)
    exact (hcont.sub_const (y 0)).sub_const (y 1)
  have h2 : Filter.Tendsto (fun n : Fin 2 →₀ ℕ => ‖MvPowerSeries.aeval hy
      ((MvPowerSeries.trunc' A n F : MvPolynomial (Fin 2) A) : MvPowerSeries (Fin 2) A) - y 0 - y 1‖)
      Filter.atTop (nhds ‖MvPowerSeries.aeval hy F - y 0 - y 1‖) :=
    (continuous_norm.continuousAt).tendsto.comp h1
  apply le_of_tendsto h2
  filter_upwards [Filter.eventually_ge_atTop (Finsupp.single (0 : Fin 2) 1 + Finsupp.single (1 : Fin 2) 1)]
    with n hn
  apply htrunc n
  · have := hn 0; simpa using this
  · have := hn 1; simpa using this

end ABC3.Found.PGC
