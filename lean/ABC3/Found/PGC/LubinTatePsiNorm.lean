import ABC3.Found.PGC.LubinTateActionPsiField
import ABC3.Found.PGC.ValuationRingComplete

/-!
# `ψ_n` の定数項のノルムは `‖π‖`(`sorry` 無し)

`Found/PGC/LubinTateActionPsi.lean::iteratedLubinTatePsi_coeff_zero_mul`
(`ψ_n.coeff0・U'_n定数項=π`)を、実在の p進局所体 `K:PAdicLocalField p`
の付値(ノルム)へ橋渡しする。`U'_n` が単元であることから、そのノルムは
ちょうど `1`(`Valued.integer.isUnit_iff_norm_eq_one`)——これで
`‖ψ_n.coeff0‖=‖π‖` が、`n` に依らず**一定**であることが分かる。

## 次に繋がる見立て(Newton polygon)

mathlib の `spectralNorm.spectralNorm_eq_norm_coeff_zero_rpow`
(`Analysis/Normed/Unbundled/SpectralNorm.lean`)は、`α` が既約多項式 `g`
の根のとき `spectralNorm K L α = ‖g.coeff 0‖^(1/g.natDegree)` という
公式を与える。`ψ_n` の根 `α` にこれを適用すると
`spectralNorm K L α = ‖π‖^(1/(q^n-q^{n-1}))` となり、指数
`1/(q^n-q^{n-1})` が `n` ごとに異なる(`0<‖π‖<1` なので冪は単射)ため、
**異なる段の根は異なるノルムを持つ**——これは `ψ_i`・`ψ_j`(`i≠j`)が
共通根を持たない(互いに素である)ことの Newton polygon 的な証明の
核心になる。ここではその橋渡しの手前、`ψ_n` 自身の定数項のノルムの
計算までを確立する。
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
——`K` の代数閉包 `K.closure` の中に `ψ_n`(`K.carrier` へ写したもの)の
根 `x` を取り(`IsAlgClosed.exists_aeval_eq_zero`)、それが既約
(`irreducible_iteratedLubinTatePsi` を Gauss の補題で `K.carrier` 上へ
橋渡し)・モニックであることから `minpoly K.carrier x` に一致することを
示し(`minpoly.eq_of_irreducible_of_monic`)、mathlib の
`spectralNorm.spectralNorm_eq_norm_coeff_zero_rpow`
(既約多項式の根のスペクトルノルムは、その多項式の定数項のノルムの
`1/次数` 乗)を `norm_iteratedLubinTatePsi_coeff_zero`(定数項のノルムは
`‖π‖`)と組み合わせて結論する。`ψ_n` は既約なので、この式は`ψ_n`の
**すべての**根について成り立つ(`minpoly` はどの根についても同じ`ψ_n`
になるため)。

`0<‖π‖<1` なので指数 `1/(q^n-q^{n-1})` は `n` ごとに異なる値を与える
——これで**異なる段の捩れ点は異なるノルムを持つ**という Newton polygon
的な事実が確立された。次の一歩は、これから「`ψ_n`・`ψ_m`(`n≠m`)は
共通根を持たない(互いに素)」を導き、`D_n` 全体の分離性・
正規性(`K(Λ_n)`が`ψ_n`の分解体と一致すること)へ繋げること。 -/
theorem exists_root_spectralNorm_iteratedLubinTatePsi {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    ∃ x : K.closure, spectralNorm K.carrier K.closure x =
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
  have hgdegpos : 0 < g.natDegree := by
    rw [hgdeg, natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]
    have h2 : 1 < pp ^ ff := hq ▸ Fintype.one_lt_card
    have hlt : (pp ^ ff) ^ (n - 1) < (pp ^ ff) ^ n := by
      apply Nat.pow_lt_pow_right h2; omega
    omega
  have hdeg : g.degree ≠ 0 := (Polynomial.natDegree_pos_iff_degree_pos.mp hgdegpos).ne'
  obtain ⟨x, hx⟩ := IsAlgClosed.exists_aeval_eq_zero K.closure g hdeg
  have hminpoly : g = minpoly K.carrier x :=
    minpoly.eq_of_irreducible_of_monic hirr hx hgmonic
  refine ⟨x, ?_⟩
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

end ABC3.Found.PGC
