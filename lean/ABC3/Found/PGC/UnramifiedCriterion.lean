import ABC3.Found.PGC.SubgroupCorrespondenceConstruction
import ABC3.Skeleton.PGC.Section1Cor13

/-!
# 原文の不分岐判定は正しい——`q_L = q^{[Γ_K:H]}` ⟺ `L_H/K` が不分岐

[pGC] p.3 は、開部分群 `H ⊆ Γ_K` に対応する中間体 `L_H` について

> Note, moreover, that L is unramified over K if and only q_L = q^[Γ_K:H].

と述べる。`Skeleton/PGC/Section1Cor13.lean` はこれを
`IsUnramifiedAt RD SC K H hH : RD.card (SC.field K H hH) = RD.card K ^ H.index`
として、**自由なデータ `RD`・`SC`** の上で写した。

本ファイルは、`Found/` で**実際に構成した** `RD = residueCardinality p`・
`SC = subgroupCorrespondence p` を代入したとき、この判定条件が
**本物の不分岐性 `IsUnramifiedAdjoin` と同値**であることを示す。
すなわち原文の判定条件は正しい。

## 使う部品(すべて本セッションで構築済み)

* `residueCardinality p`(`ResidueCardinalityConstruction.lean`)——`q` の実物
* `subgroupCorrespondence p`(`SubgroupCorrespondenceConstruction.lean`)——`L_H` の実物
* `exists_adjoin_eq_fixedField`——`L_H = K(x)`(標数 0 の原始元)
* `card_residueField_adjoinField`——`q_{K(x)} = residueDegree K x`
* `finrank_fixedField_eq_index`——`[L_H : K] = [Γ_K : H]`
* `isUnramifiedAdjoin_iff_residueDegree`——`e = 1 ⟺ f = [K(x):K]`

## ★配管: 中間体の等式に沿った `PAdicLocalField` の移送

`adjoinField K x` と `fixedFieldLocalField K H hH` は、`K(x) = L_H` の
とき同じものだが、`rw` は **`isFinite` の証明項が動かせない**ため
「motive is not type correct」で落ちる。有限性を**明示引数**で取る
`intermediateLocalField` を経由すると `subst` + `rfl`(証明無関係)で通る。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 中間体の等式に沿った移送 -/

/-- 有限次中間体を `PAdicLocalField` として見る(有限性を**明示引数**で取る)。
`adjoinField`・`fixedFieldLocalField` の共通の形。 -/
noncomputable def intermediateLocalField (K : PAdicLocalField p)
    (E : IntermediateField K.carrier K.closure) (hE : FiniteDimensional ℚ_[p] E) :
    PAdicLocalField p where
  carrier := E
  isFinite := hE

/-- 中間体が等しければ `PAdicLocalField` としても等しい(有限性は Prop なので無関係)。 -/
theorem intermediateLocalField_congr (K : PAdicLocalField p)
    {A B : IntermediateField K.carrier K.closure}
    (hA : FiniteDimensional ℚ_[p] A) (hB : FiniteDimensional ℚ_[p] B) (h : A = B) :
    intermediateLocalField K A hA = intermediateLocalField K B hB := by
  subst h; rfl

theorem adjoinField_eq (K : PAdicLocalField p) (x : K.closure) :
    adjoinField K x
      = intermediateLocalField K (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        (adjoinField K x).isFinite := rfl

theorem fixedFieldLocalField_eq (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    fixedFieldLocalField K H hH
      = intermediateLocalField K (IntermediateField.fixedField H)
        (fixedFieldLocalField K H hH).isFinite := rfl

/-- `K(x) = L_H` なら `PAdicLocalField` としても同じもの。 -/
theorem adjoinField_eq_fixedFieldLocalField (K : PAdicLocalField p) (x : K.closure)
    {H : Subgroup K.absGal} (hH : IsOpen (H : Set K.absGal))
    (hx : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      = IntermediateField.fixedField H) :
    adjoinField K x = fixedFieldLocalField K H hH := by
  rw [adjoinField_eq, fixedFieldLocalField_eq]
  exact intermediateLocalField_congr K _ _ hx

/-! ## 判定条件の正しさ -/

/-- **★★★★★原文の不分岐判定は正しい**——構成した `residueCardinality`・
`subgroupCorrespondence` の下で、`q_{L_H} = q^{[Γ_K:H]}` は
`L_H = K(x)` が**本当に不分岐**(`e = 1`)であることと同値。 -/
theorem isUnramifiedAt_iff_isUnramifiedAdjoin (K : PAdicLocalField p) {H : Subgroup K.absGal}
    (hH : IsOpen (H : Set K.absGal)) (hne : H ≠ ⊤) (x : K.closure)
    (hx : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      = IntermediateField.fixedField H) :
    IsUnramifiedAt (residueCardinality p) (subgroupCorrespondence p) K H hH
      ↔ IsUnramifiedAdjoin K x := by
  have hfield : (subgroupCorrespondence p).field K H hH = fixedFieldLocalField K H hH := by
    simp [subgroupCorrespondence, hne]
  have hrank : Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = H.index := by
    rw [hx]; exact finrank_fixedField_eq_index K H hH
  rw [IsUnramifiedAt, hfield, ← adjoinField_eq_fixedFieldLocalField K x hH hx,
    residueCardinality_adjoinField, residueCardinality_card,
    isUnramifiedAdjoin_iff_residueDegree K x, hrank]

/-- `H = Γ_K` は(`L = K` なので)常に不分岐。 -/
theorem isUnramifiedAt_top (K : PAdicLocalField p)
    (h : IsOpen ((⊤ : Subgroup K.absGal) : Set K.absGal)) :
    IsUnramifiedAt (residueCardinality p) (subgroupCorrespondence p) K ⊤ h := by
  rw [IsUnramifiedAt, (subgroupCorrespondence p).field_top K h, Subgroup.index_top, pow_one]

/-! ## 惰性群の下界 -/

/-- 判定条件を満たす `H` の固定体は最大不分岐拡大 `K^ur` に入る。 -/
theorem fixedField_le_unramifiedClosure_of_isUnramifiedAt (K : PAdicLocalField p)
    {H : Subgroup K.absGal} (hH : IsOpen (H : Set K.absGal)) (hne : H ≠ ⊤)
    (hu : IsUnramifiedAt (residueCardinality p) (subgroupCorrespondence p) K H hH) :
    IntermediateField.fixedField H ≤ unramifiedClosure K := by
  obtain ⟨x, hx⟩ := exists_adjoin_eq_fixedField K H hH
  have hux := (isUnramifiedAt_iff_isUnramifiedAdjoin K hH hne x hx).mp hu
  rw [← hx]
  exact adjoin_le_unramifiedClosure K hux

/-- **★★★★`Gal(K̄/K^ur) ≤ I_K`**——`Skeleton` が構成した惰性群
`inertia = sInf {判定条件を満たす開部分群}` は、本物の惰性群
`Gal(K̄/K^ur)` を**含む**。

(逆向きの包含には「閉部分群は自分を含む開部分群すべての共通部分」という
副有限群の事実が要る——次の段。) -/
theorem fixingSubgroup_unramifiedClosure_le_inertia (K : PAdicLocalField p) :
    (unramifiedClosure K).fixingSubgroup
      ≤ inertia (residueCardinality p) (subgroupCorrespondence p) K := by
  haveI := isGalois_closure K
  refine le_sInf ?_
  rintro H ⟨hH, hu⟩
  by_cases hne : H = ⊤
  · rw [hne]; exact le_top
  · have hle := fixedField_le_unramifiedClosure_of_isUnramifiedAt K hH hne hu
    have hclosed : IsClosed (H : Set K.absGal) := Subgroup.isClosed_of_isOpen H hH
    have h2 : (IntermediateField.fixedField H).fixingSubgroup = H :=
      InfiniteGalois.fixingSubgroup_fixedField (⟨H, hclosed⟩ : ClosedSubgroup K.absGal)
    rw [← h2]
    exact IntermediateField.fixingSubgroup_le hle

end ABC3.Found.PGC
