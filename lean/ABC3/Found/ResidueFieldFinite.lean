import Mathlib.Topology.Algebra.Valued.LocallyCompact

/-!
# 完備・proper な非アルキメデスノルム体の剰余体は有限

論文にも我々のモデルにも固有でない、**一般の**結果。mathlib へ出せる形で書く。

## 経緯(記録として残す)

当初は mathlib の
`Valued.integer.properSpace_iff_..._finite_residueField` を直接使おうとして、
**ノルムのダイヤモンド**で止まった——その補題の `ProperSpace` は
`Valued.toNormedField`(付値から作ったノルム)に関するもので、
こちらの `ProperSpace` はスペクトルノルムに関するもの。
両者は命題としては等しいが定義的に等しくない。

round-trip 補題(`Valued.toNormedField (NormedField.toValued) ≍ 元のノルム`)を
書きにいったが、**それは不要だった**——`CompactSpace 𝒪[K]` を経由すれば、
問題はノルムでなく**位相**の一致に落ちる。そして `NormedField.toValued` は
`toUniformSpace` に元のノルムのものをそのまま使うので、位相は**定義的に**一致する。

教訓: ダイヤモンドは「等しいことを証明する」より「**触れずに済ませる**」方が安い。
-/

open scoped NormedField Valued

variable {K : Type*} [NontriviallyNormedField K] [IsUltrametricDist K] [ProperSpace K]

omit [ProperSpace K] in
/-- 付値環は閉単位球にほかならない。 -/
theorem integer_eq_closedBall : (𝒪[K] : Set K) = Metric.closedBall 0 1 := by
  ext x
  simp [Valued.integer.mem_iff]

/-- proper なら付値環はコンパクト(閉単位球だから)。 -/
instance compactSpace_integer : CompactSpace 𝒪[K] :=
  isCompact_iff_compactSpace.mp (integer_eq_closedBall (K := K) ▸ isCompact_closedBall 0 1)

/-- **剰余体は有限**。 -/
theorem finite_residueField : Finite 𝓀[K] :=
  letI : (Valued.v : Valuation K NNReal).RankOne :=
    inferInstanceAs (NormedField.valuation (K := K)).RankOne
  (Valued.integer.compactSpace_iff_completeSpace_and_isDiscreteValuationRing_and_finite_residueField
    (K := K) (Γ₀ := NNReal)).mp inferInstance |>.2.2
