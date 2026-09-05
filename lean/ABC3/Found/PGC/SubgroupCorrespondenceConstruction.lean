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

end ABC3.Found.PGC
