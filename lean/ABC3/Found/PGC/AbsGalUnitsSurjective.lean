import ABC3.Found.PGC.AbsGalSurjections
import ABC3.Found.PGC.AbsoluteRamification
import ABC3.Found.PGC.LubinTateSeriesExists

/-!
# `Γ_K ↠ 𝒪_K^×` を**仮説なし**にする

`Found/PGC/AbsGalSurjections.lean::exists_surjective_absGal_to_units` は、
Lubin-Tate 相互律の射影極限(節目(5))を単数群の言葉に翻訳したものだが、
以下の仮説を引きずっていた:

* `IsAdicComplete (maximalIdeal 𝒪_K) 𝒪_K`
* 剰余体の `Fintype`・`ExpChar p`・位数 `q = pp^ff`
* 一意化元 `π` と `maximalIdeal = span {π}`・`π ≠ 0`
* Lubin-Tate 級数 `f`(`f ≡ πX mod deg 2`、`f ≡ X^q mod π`)

**これらはすべて本リポジトリで既に構築済みだった**——組み上げるだけで消える:

| 仮説 | 供給元 |
|---|---|
| `IsAdicComplete` | `ValuationRingComplete.lean::isAdicComplete_valuationRing` |
| `Fintype 𝓀`・`ExpChar p` | `ValuationRingDVR.lean` の instance 2 つ |
| `q = p^f` | `AbsoluteRamification.lean::residueCard_eq_pow` |
| `π`(DVR の一意化元) | `ValuationRingDVR.lean::valuationRing_isDVR` + mathlib の `IsDiscreteValuationRing.exists_irreducible`/`irreducible_iff_uniformizer` |
| Lubin-Tate 級数 `f` | `LubinTateSeriesExists.lean::exists_lubinTateSeries` |

結果として **任意の p進局所体 `K` について無条件に `Γ_K ↠ 𝒪_K^×`**。
古典的局所類体論の相互律 `Γ_K^ab ≅ (K^×)^∧` のうち、
`K^× ≅ ℤ × 𝒪_K^×`(`UnitsSplit.lean`)の**完全分岐側の因子**が
`Γ_K` の商として実際に現れる、ということ。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- **★★★★仮説なしの `Γ_K ↠ 𝒪_K^×`**。

`exists_surjective_absGal_to_units` の仮説を、本リポジトリで既に構築済みの
事実(DVR 性・adic 完備性・剰余体の有限性と位数・Lubin-Tate 級数の存在)で
すべて埋めたもの。組み上げだけで、新しい数学は要らなかった。 -/
theorem exists_surjective_absGal_units (K : PAdicLocalField p) :
    ∃ φ : (K.closure ≃ₐ[K.carrier] K.closure) →* (𝒪[K.carrier])ˣ, Function.Surjective φ := by
  haveI := isAdicComplete_valuationRing K
  haveI := valuationRing_isDVR K
  obtain ⟨π, hπirr⟩ := IsDiscreteValuationRing.exists_irreducible (𝒪[K.carrier])
  have hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer π).mp hπirr
  have hπne0 : π ≠ 0 := hπirr.ne_zero
  have hq : Fintype.card 𝓀[K.carrier] = p ^ (absoluteInertiaDegree K) := by
    rw [← Nat.card_eq_fintype_card]
    exact residueCard_eq_pow K
  obtain ⟨f, hf0, hf1, hf⟩ := exists_lubinTateSeries (A := 𝒪[K.carrier]) hq hπmax
  exact exists_surjective_absGal_to_units K hq hπmax hπne0 f hf0 hf1 hf

end ABC3.Found.PGC
