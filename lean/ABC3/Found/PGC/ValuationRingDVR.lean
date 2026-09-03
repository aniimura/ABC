import ABC3.Found.PGC.LocalFieldNorm
import Mathlib.Topology.Algebra.Valued.LocallyCompact

/-!
# 任意の p進局所体の付値環は離散付値環(`sorry` 無し)

Lubin-Tate 形式群の構成(`ResearchPaper/blocked-leaves.json` の
`[pGC] Proposition 1.2 / ...` の `progress2026_09_04f`/`g` で特定した3層の
うちの層 (a))で必要になった、「`𝒪[K.carrier]` の極大イデアルが一意化元1つで
生成される」という事実——**`IsDiscreteValuationRing`**——を、一般の p進局所体
`K`(`PAdicLocalField p`)について確立する。

## ★見積りの訂正

`progress2026_09_04g` では「`RankOne` インスタンスだけでは
`IsCyclic`/`Nontrivial (valueGroup)` が自動的に出ない(値群が離散である保証が
無い)」と記録し、これを解くには局所体論の標準事実(有限次拡大は再び離散付値を
持つ)を独立に証明する必要があると見積もっていた。**実際にはもっと直接の経路が
あった**——mathlib の `Valued.integer.isDiscreteValuationRing_of_compactSpace`
(`Topology/Algebra/Valued/LocallyCompact.lean`)が、`valueGroup` の離散性を
経由せず、**付値環がコンパクトであること**だけから直接 `IsDiscreteValuationRing`
を与える。そして `𝒪[K.carrier]`(= 閉単位球、`Valued.integer.mem_iff` で確認)が
コンパクトであることは、`K.carrier` が(`Found/PGC/LocalFieldNorm.lean` で
既に確立済みの)`ProperSpace`(有限次元 ℚ_p-ベクトル空間だから)であることから
即座に出る——`isCompact_closedBall` を当てるだけでよい。

**教訓**: 「値群の離散性」という抽象的な入口で詰まったら、目的の結論
(`IsDiscreteValuationRing`)への**別の経路**(ここでは「コンパクト性」)を
mathlib で探す価値がある——入口を1つに決め打ちしない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-- ★★**任意の p進局所体 `K` の付値環は離散付値環である**。

`𝒪[K.carrier]` がコンパクト(= 閉単位球、`K.carrier` が有限次元 ℚ_p-ベクトル空間
ゆえ `ProperSpace`)であることから、`Valued.integer.
isDiscreteValuationRing_of_compactSpace` で直接得る。

これにより `𝒪[K.carrier]` は `IsPrincipalIdealRing`(`IsDiscreteValuationRing`
が `extends` する)を持つ——極大イデアルは一意化元1つで生成される。 -/
theorem valuationRing_isDVR {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    IsDiscreteValuationRing (𝒪[K.carrier]) := by
  have hcs : CompactSpace (𝒪[K.carrier]) := by
    have hset : (𝒪[K.carrier] : Set K.carrier) = Metric.closedBall 0 1 := by
      ext x
      simp only [Metric.mem_closedBall, dist_eq_norm, sub_zero]
      exact Valued.integer.mem_iff
    have hcompact : IsCompact (𝒪[K.carrier] : Set K.carrier) := hset ▸ isCompact_closedBall 0 1
    exact isCompact_iff_compactSpace.mp hcompact
  exact Valued.integer.isDiscreteValuationRing_of_compactSpace

/-- **剰余体 `𝓀[K.carrier]` は標数 `p` を持つ**。`(p:K.carrier)` が単数でない
(`‖(p:K.carrier)‖ < 1`、`norm_natCast_p_lt_one`)ことから、その `𝒪[K.carrier]`
への持ち上げが極大イデアルに入り、還元 `residue` で 0 になることによる。 -/
theorem charP_residueField {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    CharP (𝓀[K.carrier]) p := by
  rw [CharP.charP_iff_prime_eq_zero Fact.out]
  set u : 𝒪[K.carrier] :=
    ⟨(p : ℕ), by rw [Valued.integer.mem_iff]; exact le_of_lt (norm_natCast_p_lt_one K)⟩ with hu
  have hnu : (p : 𝓀[K.carrier]) = IsLocalRing.residue _ u := by rw [hu]; rfl
  rw [hnu, IsLocalRing.residue_eq_zero_iff, IsLocalRing.mem_maximalIdeal, mem_nonunits_iff,
    Valued.integer.isUnit_iff_norm_eq_one]
  have hun : ‖u‖ = ‖((p : ℕ) : K.carrier)‖ := rfl
  rw [hun]
  exact ne_of_lt (norm_natCast_p_lt_one K)

/-- 剰余体は有限体であり、標数 `p`(`p` は特にその標数指数)。
`Found/PGC/FrobeniusExpand.lean::mvPowerSeries_pow_card_eq_expand` の
`[ExpChar F p]` 要件をこれで満たせる。 -/
noncomputable instance instFintypeResidueField {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    Fintype (𝓀[K.carrier]) :=
  have := residueField_finite K
  Fintype.ofFinite _

noncomputable instance instExpCharResidueField {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    ExpChar (𝓀[K.carrier]) p :=
  haveI := charP_residueField K
  ExpChar.prime Fact.out

end ABC3.Found.PGC
