import ABC3.Found.PGC.LubinTateActionPsiField
import ABC3.Found.PGC.ValuationRingComplete

/-!
# `ψ_n`・`ψ_m`(`n≠m`)は共通根を持たない——Newton polygon(`sorry` 無し)

古典的な Lubin-Tate 理論の Newton polygon 論法を、実在の p進局所体
`K:PAdicLocalField p` について確立する。

## 証明の筋

1. `Found/PGC/LubinTateActionPsi.lean::iteratedLubinTatePsi_coeff_zero_mul`
   (`ψ_n.coeff0・U'_n定数項=π`)の両辺のノルムを取り、`U'_n` が単元で
   ノルムがちょうど `1`(`Valued.integer.isUnit_iff_norm_eq_one`)である
   ことから、`‖ψ_n.coeff0‖=‖π‖`(`n` に依らず一定)を得る
   (`norm_iteratedLubinTatePsi_coeff_zero`)。
2. `K` の代数閉包 `K.closure` の中で `ψ_n`(`K.carrier` へ写したもの)が
   既約(Gauss の補題で橋渡し)・モニックであることから、**任意の**根 `x`
   について `minpoly K.carrier x = ψ_n` となることを示し、mathlib の
   `spectralNorm.spectralNorm_eq_norm_coeff_zero_rpow`
   (既約多項式の根のスペクトルノルムは、その多項式の定数項のノルムの
   `1/次数` 乗)と 1. を組み合わせて
   `spectralNorm K.carrier K.closure x = ‖π‖^(1/(q^n-q^{n-1}))` を得る
   (`spectralNorm_root_iteratedLubinTatePsi`)。
3. `0<‖π‖<1`(`norm_pi_pos_lt_one`)なので、指数 `1/(q^n-q^{n-1})` は
   `n` ごとに異なる値を与える(`torsionDegree_ne`・
   `rpow_ne_rpow_of_base_lt_one`)——**異なる段の根は異なるノルムを持つ**。
4. これで `ψ_n`・`ψ_m`(`n≠m`)が共通根を持たないことが従う
   (`no_common_root_iteratedLubinTatePsi`)——共通根があれば同じ `x` の
   スペクトルノルムが2つの異なる値に等しくなり矛盾する。

これは `D_n` 全体の分離性(異なる段の`ψ_i`たちが互いに素であること)・
正規性(`K(Λ_n)`が`ψ_n`の分解体と一致すること)へ向けた核心的な一歩。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-- ★★★★★★★★★**`ψ_n` の定数項のノルムは `‖π‖`**(`n` に依らず一定)——
`ψ_n.coeff0・U'_n定数項=π`(`iteratedLubinTatePsi_coeff_zero_mul`)の
両辺のノルムを取り、`U'_n` の定数項が単元でありノルムが`1`であること
(`Valued.integer.isUnit_iff_norm_eq_one`)で簡約する。 -/
theorem norm_iteratedLubinTatePsi_coeff_zero {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    ‖(iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0‖ = ‖π‖ := by
  haveI := valuationRing_isDVR K
  have hone := iteratedLubinTatePsi_coeff_zero_mul hq hπmax hπne0 f hf0 hf1 hf n hn
  have hUunit : IsUnit (iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn).constantCoeff :=
    PowerSeries.isUnit_constantCoeff _ (isUnit_iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn)
  have hUnorm : ‖(iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn).constantCoeff‖ = 1 :=
    Valued.integer.isUnit_iff_norm_eq_one.mp hUunit
  have hcongr := congrArg norm hone
  rw [norm_mul, hUnorm, mul_one] at hcongr
  exact hcongr

/-- ★★★★★★★★★★★**`ψ_n` の根のスペクトルノルムは `‖π‖^(1/(q^n-q^{n-1}))`**
(`α` が `ψ_n`(`K.carrier` へ写したもの)の**任意の**根のとき)——`ψ_n` が
既約(`irreducible_iteratedLubinTatePsi` を Gauss の補題で `K.carrier`
上へ橋渡し)・モニックであることから、任意の根 `x` について
`minpoly K.carrier x` が `ψ_n` 自身に一致することを示し
(`minpoly.eq_of_irreducible_of_monic`)、mathlib の
`spectralNorm.spectralNorm_eq_norm_coeff_zero_rpow`
(既約多項式の根のスペクトルノルムは、その多項式の定数項のノルムの
`1/次数` 乗)を `norm_iteratedLubinTatePsi_coeff_zero`(定数項のノルムは
`‖π‖`)と組み合わせて結論する。 -/
theorem spectralNorm_root_iteratedLubinTatePsi {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hx : Polynomial.aeval x (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)) = 0) :
    spectralNorm K.carrier K.closure x =
      ‖π‖ ^ (1 / ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).natDegree : ℝ)) := by
  haveI := valuationRing_isDVR K
  haveI : UniqueFactorizationMonoid (𝒪[K.carrier]) := uniqueFactorizationMonoid_valuationRing K
  have hnormcoeff := norm_iteratedLubinTatePsi_coeff_zero K hq hπmax hπne0 f hf0 hf1 hf n hn
  set g := Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn) with hg_def
  have hmonic := (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
  have hirr : Irreducible g :=
    (Polynomial.IsPrimitive.irreducible_iff_irreducible_map_fraction_map hmonic.isPrimitive).mp
      (irreducible_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)
  have hgmonic : g.Monic := hmonic.map _
  have hgdeg : g.natDegree = (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).natDegree :=
    hmonic.natDegree_map _
  have hminpoly : g = minpoly K.carrier x :=
    minpoly.eq_of_irreducible_of_monic hirr hx hgmonic
  rw [spectralNorm.spectralNorm_eq_norm_coeff_zero_rpow, ← hminpoly]
  have hgcoeff0 : g.coeff 0 = algebraMap (𝒪[K.carrier]) K.carrier
      ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0) := by
    rw [hg_def, Polynomial.coeff_map]
  rw [hgcoeff0]
  have hnormcoe : ‖algebraMap (𝒪[K.carrier]) K.carrier
      ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0)‖ =
      ‖(iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0‖ := by
    show ‖((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0 : K.carrier)‖ = _
    exact (AddSubgroupClass.coe_norm 𝒪[K.carrier] _).symm
  rw [hnormcoe, hnormcoeff, hgdeg]

/-- `ψ_n`(`K.carrier` へ写したもの)は代数閉体 `K.closure` に少なくとも
1つ根を持つ(既約・モニックなので次数は正、`IsAlgClosed`)。 -/
theorem exists_root_iteratedLubinTatePsi {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    ∃ x : K.closure, Polynomial.aeval x (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)) = 0 := by
  haveI := valuationRing_isDVR K
  haveI : UniqueFactorizationMonoid (𝒪[K.carrier]) := uniqueFactorizationMonoid_valuationRing K
  have hmonic := (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
  set g := Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn) with hg_def
  have hgdeg : g.natDegree = (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).natDegree :=
    hmonic.natDegree_map _
  have hgdegpos : 0 < g.natDegree := by
    rw [hgdeg, natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]
    have h2 : 1 < pp ^ ff := hq ▸ Fintype.one_lt_card
    have hlt : (pp ^ ff) ^ (n - 1) < (pp ^ ff) ^ n := by
      apply Nat.pow_lt_pow_right h2; omega
    omega
  exact IsAlgClosed.exists_aeval_eq_zero K.closure g
    (Polynomial.natDegree_pos_iff_degree_pos.mp hgdegpos).ne'

/-- `π` のノルムは `0` より大きく `1` 未満(`π≠0`・`π` は単元でない)。 -/
theorem norm_pi_pos_lt_one {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0) :
    0 < ‖π‖ ∧ ‖π‖ < 1 := by
  haveI := valuationRing_isDVR K
  refine ⟨by rw [norm_pos_iff]; exact hπne0, ?_⟩
  rw [lt_iff_le_and_ne]
  refine ⟨Valued.integer.norm_le_one π, ?_⟩
  intro hcon
  have hunit : IsUnit π := Valued.integer.isUnit_iff_norm_eq_one.mpr hcon
  have hmem : π ∈ IsLocalRing.maximalIdeal (𝒪[K.carrier]) := by
    rw [hπmax]; exact Ideal.subset_span rfl
  rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff] at hmem
  exact hmem hunit

/-- `‖π‖^e` は `0<‖π‖<1` のとき指数 `e` について単射——底が `1` 未満の
冪の狭義単調減少性(`Real.rpow_lt_rpow_left_iff_of_base_lt_one`)から。 -/
theorem rpow_ne_rpow_of_base_lt_one {a e1 e2 : ℝ} (ha0 : 0 < a) (ha1 : a < 1) (he : e1 ≠ e2) :
    a ^ e1 ≠ a ^ e2 := by
  rcases lt_or_gt_of_ne he with h | h
  · exact ((Real.rpow_lt_rpow_left_iff_of_base_lt_one ha0 ha1).mpr h).ne'
  · exact ((Real.rpow_lt_rpow_left_iff_of_base_lt_one ha0 ha1).mpr h).ne

/-- `n≠m`(`n,m≥1`)ならば `q^n-q^{n-1} ≠ q^m-q^{m-1}`——
`q^k-q^{k-1}=q^{k-1}(q-1)` と書き直し、`q^{k-1}` の狭義単調性で従う。 -/
theorem torsionDegree_ne {q n m : ℕ} (hq : 1 < q) (hn : 1 ≤ n) (hm : 1 ≤ m) (hnm : n ≠ m) :
    q ^ n - q ^ (n - 1) ≠ q ^ m - q ^ (m - 1) := by
  have key : ∀ a b : ℕ, 1 ≤ a → 1 ≤ b → a < b → q ^ a - q ^ (a - 1) ≠ q ^ b - q ^ (b - 1) := by
    intro a b ha hb hlt heq
    have h1 : q ^ a - q ^ (a - 1) = q ^ (a - 1) * (q - 1) := by
      rw [Nat.mul_sub_one]
      have hpa : q ^ (a - 1) * q = q ^ a := by rw [← pow_succ]; congr 1; omega
      rw [hpa]
    have h2 : q ^ b - q ^ (b - 1) = q ^ (b - 1) * (q - 1) := by
      rw [Nat.mul_sub_one]
      have hpb : q ^ (b - 1) * q = q ^ b := by rw [← pow_succ]; congr 1; omega
      rw [hpb]
    rw [h1, h2] at heq
    have hq1 : 0 < q - 1 := by omega
    have heq' := Nat.eq_of_mul_eq_mul_right hq1 heq
    have hlt' : q ^ (a - 1) < q ^ (b - 1) := Nat.pow_lt_pow_right hq (by omega)
    omega
  rcases lt_or_gt_of_ne hnm with h | h
  · exact key n m hn hm h
  · exact fun heq => key m n hm hn h heq.symm

/-- ★★★★★★★★★★★★**`ψ_n`・`ψ_m`(`n≠m`)は共通根を持たない**
——古典的な Lubin-Tate 理論の Newton polygon 論法の核心。共通根 `x` が
あれば、その `spectralNorm` は `‖π‖^(1/(q^n-q^{n-1}))` にも
`‖π‖^(1/(q^m-q^{m-1}))` にも等しくなるが、`n≠m` なら指数が異なり
(`torsionDegree_ne`)、`0<‖π‖<1` なのでこの2つの冪の値は異なる
(`rpow_ne_rpow_of_base_lt_one`)——矛盾。これで `ψ_n` たちが(代数閉体上)
互いに素であることの核心部分が確立された。 -/
theorem no_common_root_iteratedLubinTatePsi {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {n m : ℕ} (hn : 1 ≤ n) (hm : 1 ≤ m) (hnm : n ≠ m) (x : K.closure)
    (hxn : Polynomial.aeval x (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)) = 0)
    (hxm : Polynomial.aeval x (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf m hm)) = 0) : False := by
  have hnorm_n := spectralNorm_root_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x hxn
  have hnorm_m := spectralNorm_root_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf m hm x hxm
  rw [natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn] at hnorm_n
  rw [natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf m hm] at hnorm_m
  obtain ⟨hπpos, hπlt1⟩ := norm_pi_pos_lt_one K hπmax hπne0
  have hdegne : (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) ≠ (pp ^ ff) ^ m - (pp ^ ff) ^ (m - 1) :=
    torsionDegree_ne (hq ▸ Fintype.one_lt_card) hn hm hnm
  have hexpne : (1 : ℝ) / (((pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) : ℕ) : ℝ) ≠
      (1 : ℝ) / (((pp ^ ff) ^ m - (pp ^ ff) ^ (m - 1) : ℕ) : ℝ) := by
    intro heq
    apply hdegne
    have hne_n : ((((pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) : ℕ)) : ℝ) ≠ 0 := by
      have hpos : 0 < (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
        have h2 : 1 < pp ^ ff := hq ▸ Fintype.one_lt_card
        have hlt : (pp ^ ff) ^ (n - 1) < (pp ^ ff) ^ n := Nat.pow_lt_pow_right h2 (by omega)
        omega
      exact_mod_cast hpos.ne'
    have hne_m : ((((pp ^ ff) ^ m - (pp ^ ff) ^ (m - 1) : ℕ)) : ℝ) ≠ 0 := by
      have hpos : 0 < (pp ^ ff) ^ m - (pp ^ ff) ^ (m - 1) := by
        have h2 : 1 < pp ^ ff := hq ▸ Fintype.one_lt_card
        have hlt : (pp ^ ff) ^ (m - 1) < (pp ^ ff) ^ m := Nat.pow_lt_pow_right h2 (by omega)
        omega
      exact_mod_cast hpos.ne'
    field_simp at heq
    exact_mod_cast heq.symm
  have hne := rpow_ne_rpow_of_base_lt_one hπpos hπlt1 hexpne
  apply hne
  rw [← hnorm_n, ← hnorm_m]

end ABC3.Found.PGC
