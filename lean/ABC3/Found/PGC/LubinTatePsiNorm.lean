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

end ABC3.Found.PGC
