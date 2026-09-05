import ABC3.Found.PGC.AdjoinFieldConstruction
import Mathlib.FieldTheory.Galois.Infinite

/-!
# `Interface.PGC.SubgroupCorrespondence` を**構成する**

`Interface/PGC/LocalFieldData.lean` の

```
structure SubgroupCorrespondence where
  field : (K) → (H : Subgroup K.absGal) → IsOpen ↑H → PAdicLocalField p
  field_top : ∀ K h, field K ⊤ h = K
```

は自由なデータだった(`waiting` の `trackB` は「Krull 位相の Galois 対応と、
`ℚ_p` 上の有限次性の接続」と記していた)。本ファイルはそれを構成する。

## 筋

`H ⊆ Γ_K` が開なら `1 ∈ H` なので `H ∈ 𝓝 1`。Krull 位相の近傍基底
(`krullTopology_mem_nhds_one_iff`)から、`K` 上有限次の中間体 `E` で
`E.fixingSubgroup ⊆ H` なるものが取れる。`fixedField` は反変なので

```
fixedField H ≤ fixedField (E.fixingSubgroup) = E
```

(最後の等号は `InfiniteGalois.fixedField_fixingSubgroup`、`K̄/K` が
無限次 Galois——標数 0 なので分離的、代数閉包なので normal)。よって
`fixedField H` は `K` 上有限次、したがって塔で `ℚ_p` 上も有限次。

## ★逸脱の記録: `field_top`

`field_top : field K ⊤ h = K` は `PAdicLocalField p` の**構造としての等式**
なので、`carrier` が**型として** `K.carrier` に一致しなければならない。
`fixedField ⊤ = ⊥` の台は `↥(⊥ : IntermediateField K.carrier K.closure)` で、
`K.carrier` と標準同型ではあるが**型としては別物**。そこで `H = ⊤` の場合
だけ場合分けして `K` そのものを返す。`H ≠ ⊤` では本物の固定体を返すので、
数学的な内容は変わらない(`H = ⊤` の固定体は `K` に他ならない)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NormedField Valued Classical

variable {p : ℕ} [Fact p.Prime]

/-- `K̄/K` は無限次 Galois(標数 0 なので分離的、代数閉包なので normal)。 -/
theorem isGalois_closure (K : PAdicLocalField p) : IsGalois K.carrier K.closure := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI := IsAlgClosure.normal K.carrier K.closure
  exact ⟨⟩

/-- **開部分群の固定体は `K` 上有限次**——Krull 位相の近傍基底から。 -/
theorem finiteDimensional_fixedField_of_isOpen (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    FiniteDimensional K.carrier (IntermediateField.fixedField H) := by
  haveI := isGalois_closure K
  have hmem : (H : Set K.absGal) ∈ nhds 1 := hH.mem_nhds H.one_mem
  obtain ⟨E, hEfin, hEsub⟩ :=
    (krullTopology_mem_nhds_one_iff K.carrier K.closure (H : Set K.absGal)).mp hmem
  haveI := hEfin
  have hle : IntermediateField.fixedField H ≤ E := by
    have h1 : IntermediateField.fixedField H ≤ IntermediateField.fixedField E.fixingSubgroup := by
      intro z hz g
      exact hz ⟨g.1, hEsub g.2⟩
    rw [InfiniteGalois.fixedField_fixingSubgroup E] at h1
    exact h1
  exact Module.Finite.of_injective (IntermediateField.inclusion hle).toLinearMap
    (IntermediateField.inclusion hle).injective

/-- 開部分群の固定体を `PAdicLocalField p` として。 -/
noncomputable def fixedFieldLocalField (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) : PAdicLocalField p where
  carrier := IntermediateField.fixedField H
  isFinite := by
    haveI := finiteDimensional_fixedField_of_isOpen K H hH
    haveI : IsScalarTower ℚ_[p] K.carrier (IntermediateField.fixedField H) := by infer_instance
    exact Module.Finite.trans (R := ℚ_[p]) K.carrier (IntermediateField.fixedField H)

/-- **★★`SubgroupCorrespondence p` の構成**。`H = ⊤` のときだけ `K` 自身を
返す(上の逸脱の記録を参照——固定体は `K` と標準同型だが型としては別物)。 -/
noncomputable def subgroupCorrespondence (p : ℕ) [Fact p.Prime] : SubgroupCorrespondence p where
  field K H hH := if H = ⊤ then K else fixedFieldLocalField K H hH
  field_top K h := by simp


/-! ## 有限部分拡大はすべて単項——`adjoin` の機構が固定体にも届く

標数 0 なので `K̄/K` の有限次中間体はすべて分離的、したがって原始元定理で
**単項**:`E = K(x)`。これで `Found/PGC/UnramifiedExtension.lean` で
`K(x)` について積み上げた理論(`e·f`・不分岐性・Frobenius・一意性)が、
開部分群の固定体 `L_H = fixedField H` にもそのまま適用できる
——原文が「Proposition 1.2 を `(L, H)` に適用する」と書く操作の土台。 -/

/-- **有限次中間体は単項**(原始元定理、標数 0)。 -/
theorem exists_adjoin_eq_of_finiteDimensional (K : PAdicLocalField p)
    (E : IntermediateField K.carrier K.closure) [FiniteDimensional K.carrier E] :
    ∃ x : K.closure, IntermediateField.adjoin K.carrier ({x} : Set K.closure) = E := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI : Algebra.IsSeparable K.carrier E :=
    IntermediateField.isSeparable_tower_bot K.carrier E
  obtain ⟨α, hα⟩ := Field.exists_primitive_element K.carrier E
  refine ⟨(α : K.closure), ?_⟩
  have hle : IntermediateField.adjoin K.carrier ({(α : K.closure)} : Set K.closure) ≤ E := by
    rw [IntermediateField.adjoin_simple_le_iff]
    exact α.2
  refine IntermediateField.eq_of_le_of_finrank_eq hle ?_
  have hint : IsIntegral K.carrier (α : K.closure) :=
    IsAlgebraic.isIntegral (Algebra.IsAlgebraic.isAlgebraic _)
  have hintE : IsIntegral K.carrier α := IsIntegral.of_finite _ _
  have hmp : minpoly K.carrier (α : K.closure) = minpoly K.carrier α :=
    minpoly.algebraMap_eq (A := K.carrier) (B := (E : IntermediateField K.carrier K.closure))
      (B' := K.closure) (IntermediateField.val E).injective α
  rw [IntermediateField.adjoin.finrank hint,
    ← IntermediateField.finrank_top' (F := K.carrier) (E := E), ← hα,
    IntermediateField.adjoin.finrank hintE, hmp]

/-- **開部分群の固定体も単項**——`L_H = K(x)`。 -/
theorem exists_adjoin_eq_fixedField (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    ∃ x : K.closure, IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      = IntermediateField.fixedField H := by
  haveI := finiteDimensional_fixedField_of_isOpen K H hH
  exact exists_adjoin_eq_of_finiteDimensional K (IntermediateField.fixedField H)


/-- **開部分群の指数 = 固定体の次数**——無限次 Galois 対応
(`IntermediateField.finrank_eq_fixingSubgroup_index` と
`InfiniteGalois.fixingSubgroup_fixedField`、開部分群は閉)。
原文の判定条件 `q_L = q^{[Γ_K:H]}` の指数 `[Γ_K:H]` が
`[L_H : K]` であることの型的な裏付け。 -/
theorem finrank_fixedField_eq_index (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    Module.finrank K.carrier (IntermediateField.fixedField H) = H.index := by
  haveI := isGalois_closure K
  have hclosed : IsClosed (H : Set K.absGal) := Subgroup.isClosed_of_isOpen H hH
  have h1 := IntermediateField.finrank_eq_fixingSubgroup_index
    (k := K.carrier) (K := K.closure) (IntermediateField.fixedField H)
  have h2 : (IntermediateField.fixedField H).fixingSubgroup = H :=
    InfiniteGalois.fixingSubgroup_fixedField (⟨H, hclosed⟩ : ClosedSubgroup K.absGal)
  rw [h1, h2]

end ABC3.Found.PGC

/-- **G2 の非空虚 witness**——`SubgroupCorrespondence p` は実際に存在する。

`Interface` は `Found` を import できない(`Interface/PGC/LocalFieldData.lean`
冒頭の規約)ので、実装側から `Interface.PGC` 名前空間へ足す
——`ResidueCardinality.nonvacuous` と同じ配置。

これにより、`SubgroupCorrespondence` を仮説に持つ `Skeleton/` の主張
(`Skeleton/PGC/Section1Cor13.lean`)は**空虚でないことが確定した**。 -/
theorem ABC3.Interface.PGC.SubgroupCorrespondence.nonvacuous (p : ℕ) [Fact p.Prime] :
    Nonempty (ABC3.Interface.PGC.SubgroupCorrespondence p) :=
  ⟨ABC3.Found.PGC.subgroupCorrespondence p⟩
