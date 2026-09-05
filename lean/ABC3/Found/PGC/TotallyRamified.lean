import ABC3.Found.PGC.UnramifiedSubextension

/-!
# 完全分岐と、不分岐との交わり

`UnramifiedSubextension.lean` で `K(x) ≤ K^ur ⟺ K(x)/K が不分岐` まで来た。
相互律の全体像(`Γ_K^ab ≅ (K^×)^∧`)に必要な**線型無関係性**
「完全分岐 ∩ K^ur = K」へ向けて、まずその両端を用意する。

* `IsTotallyRamifiedAdjoin K x := inertiaDegree K x = 1`(慣性次数が 1)
* `finrank_eq_one_of_isUnramified_of_isTotallyRamified`——
  不分岐かつ完全分岐なら次数 1(`e·f = [K(x):K]` の直接の帰結)
* `isUnit_adjoinIntegers_iff`——`adjoinIntegers K x` の単数はノルム 1 の元。
  部分拡大に沿った整数環の包含が局所準同型であることを示すのに要る。
* `inertiaDegree_le_of_adjoin_le`——慣性次数の単調性(剰余体は伸びる一方)
* **`isTotallyRamified_of_le`**——完全分岐は部分拡大に遺伝する
* **`finrank_eq_one_of_mem_unramifiedClosure_of_le`**——「完全分岐 ∩ K^ur = K」

## ★配管の記録: 中間体をまたぐ `rfl` は kernel を止める

`K(x) ≤ K(y)` に沿った整数環の包含 `adjoinIntegers K x → adjoinIntegers K y` は
数学的には自明(中間体のノルムは `K.closure` のノルムの制限そのもの)だが、
**Lean では素直に書けない**:

* `IntermediateField.inclusion hle` を使うと、`((inclusion hle z : ↥K⟮y⟯) : K.closure)
  = ((z : ↥K⟮x⟯) : K.closure)` の `rfl` が **kernel deterministic timeout**(実測 60 秒)。
* `⟨(z : K.closure), hle z.2⟩` と素直に書けば、**中間体の元 1 層なら 0.06 秒**で通る
  (`‖⟨(z : K.closure), hle w.2⟩‖ = ‖w‖` は `rfl`)。
* ところが `adjoinIntegers`(中間体の部分環)の元は**2 層**なので、
  同じ形の `def` でも kernel が落ちる(実測 60 秒)。

つまり「1 層は速い・2 層は落ちる」。

**★回避策が効いた**: `def` の中で membership を defeq に頼らせず、
ノルム保存 `norm_mk_of_le` を**先に**補題として用意して `rw` で渡すと、
60 秒の kernel timeout が **0.12 秒**になる。これで包含
`adjoinIntegersIncl` / `adjoinIntegersRingHom` が作れ、剰余体の比較
`𝓀_x ↪ 𝓀_y` から**慣性次数の単調性**が出て、目標の

**`finrank_eq_one_of_mem_unramifiedClosure_of_le`
(= 完全分岐拡大の中の不分岐な部分は自明 = 「完全分岐 ∩ K^ur = K」)**

まで届いた。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-- **完全分岐**——慣性次数が 1(剰余体が伸びない)。 -/
def IsTotallyRamifiedAdjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] : Prop :=
  inertiaDegree K x = 1

/-- **不分岐かつ完全分岐なら次数 1**——`e·f = [K(x):K]` の直接の帰結。
「完全分岐 ∩ K^ur = K」の**片方の端**。 -/
theorem finrank_eq_one_of_isUnramified_of_isTotallyRamified (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (hu : IsUnramifiedAdjoin K x) (ht : IsTotallyRamifiedAdjoin K x) :
    Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = 1 := by
  have h := ramificationIndex_mul_inertiaDegree K x
  rw [show ramificationIndex K x = 1 from hu, show inertiaDegree K x = 1 from ht] at h
  omega

/-- `x ∈ K^ur` かつ `K(x)/K` が完全分岐なら `K(x) = K`。 -/
theorem finrank_eq_one_of_mem_unramifiedClosure_of_isTotallyRamified (K : PAdicLocalField p)
    (x : K.closure) (hx : x ∈ unramifiedClosure K) (ht : IsTotallyRamifiedAdjoin K x) :
    Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = 1 :=
  finrank_eq_one_of_isUnramified_of_isTotallyRamified K x
    ((mem_unramifiedClosure_iff_isUnramified K x).mp hx) ht

/-- **`adjoinIntegers K x` の単数はノルム 1 の元**。
部分拡大に沿った整数環の包含が局所準同型であることを示すのに要る。 -/
theorem isUnit_adjoinIntegers_iff (K : PAdicLocalField p) (x : K.closure)
    (z : adjoinIntegers K x) :
    IsUnit z ↔ ‖(z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ = 1 := by
  constructor
  · rintro ⟨u, rfl⟩
    have h1 : ((u : adjoinIntegers K x) : IntermediateField.adjoin K.carrier
        ({x} : Set K.closure)) * ((↑u⁻¹ : adjoinIntegers K x) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = 1 :=
      congrArg Subtype.val u.mul_inv
    have h2 := congrArg norm h1
    rw [norm_mul, norm_one] at h2
    have hu : ‖((u : adjoinIntegers K x) : IntermediateField.adjoin K.carrier
      ({x} : Set K.closure))‖ ≤ 1 := (u : adjoinIntegers K x).2
    have hv : ‖((↑u⁻¹ : adjoinIntegers K x) : IntermediateField.adjoin K.carrier
      ({x} : Set K.closure))‖ ≤ 1 := (↑u⁻¹ : adjoinIntegers K x).2
    nlinarith [norm_nonneg ((u : adjoinIntegers K x) : IntermediateField.adjoin K.carrier
      ({x} : Set K.closure)), norm_nonneg ((↑u⁻¹ : adjoinIntegers K x) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
  · intro hz
    have hz0 : ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ≠ 0 := by
      intro h0
      rw [h0, norm_zero] at hz
      exact zero_ne_one hz
    refine ⟨⟨z, ⟨((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))⁻¹, ?_⟩, ?_, ?_⟩,
      rfl⟩
    · show ‖((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))⁻¹‖ ≤ 1
      rw [norm_inv, hz, inv_one]
    · apply Subtype.ext
      show ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
        * ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))⁻¹ = 1
      exact mul_inv_cancel₀ hz0
    · apply Subtype.ext
      show ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))⁻¹
        * ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))) = 1
      exact inv_mul_cancel₀ hz0

/-! ## 部分拡大に沿った整数環の包含

★上の docstring に書いた壁の**回避策が効いた**: `def` の中で membership を
defeq に頼らせず、ノルム保存 `norm_mk_of_le` を**先に**補題として用意して
`rw` で渡すと、60 秒の kernel timeout が **0.12 秒**になる。 -/

/-- 中間体をまたぐノルムの保存(1 層なので `rfl` が速い)。 -/
theorem norm_mk_of_le (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (w : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ‖(⟨(w : K.closure), hle w.2⟩ : IntermediateField.adjoin K.carrier
      ({y} : Set K.closure))‖ = ‖w‖ := rfl

/-- `K(x) ≤ K(y)` に沿った整数環の包含(写像)。 -/
noncomputable def adjoinIntegersIncl (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (z : adjoinIntegers K x) : adjoinIntegers K y :=
  ⟨⟨((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure),
      hle (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)).2⟩,
    by
      show ‖(⟨((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure),
        hle (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)).2⟩ :
          IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ ≤ 1
      rw [norm_mk_of_le K hle]
      exact z.2⟩

@[simp] theorem norm_adjoinIntegersIncl (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (z : adjoinIntegers K x) :
    ‖((adjoinIntegersIncl K hle z : adjoinIntegers K y)
      : IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖
      = ‖((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))‖ :=
  norm_mk_of_le K hle _

/-- 環準同型としての包含。 -/
noncomputable def adjoinIntegersRingHom (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    adjoinIntegers K x →+* adjoinIntegers K y where
  toFun := adjoinIntegersIncl K hle
  map_one' := by apply Subtype.ext; apply Subtype.ext; rfl
  map_mul' _ _ := by apply Subtype.ext; apply Subtype.ext; rfl
  map_zero' := by apply Subtype.ext; apply Subtype.ext; rfl
  map_add' _ _ := by apply Subtype.ext; apply Subtype.ext; rfl

@[simp] theorem adjoinIntegersRingHom_apply (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (z : adjoinIntegers K x) :
    adjoinIntegersRingHom K hle z = adjoinIntegersIncl K hle z := rfl

/-- 包含は局所準同型(ノルムを保つから)。 -/
instance isLocalHom_adjoinIntegersRingHom (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    IsLocalHom (adjoinIntegersRingHom K hle) := by
  refine ⟨fun z hz => ?_⟩
  rw [isUnit_adjoinIntegers_iff] at hz ⊢
  rw [adjoinIntegersRingHom_apply, norm_adjoinIntegersIncl] at hz
  exact hz

/-! ## 慣性次数の単調性、そして「完全分岐 ∩ K^ur = K」 -/

/-- **剰余体は伸びる一方**——`K(x) ≤ K(y)` なら `𝓀_x ↪ 𝓀_y`。 -/
theorem card_residueField_le_of_adjoin_le (K : PAdicLocalField p) {x y : K.closure}
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    residueDegree K x ≤ residueDegree K y := by
  haveI := isLocalRing_adjoinIntegers K x
  haveI := isLocalRing_adjoinIntegers K y
  exact Nat.card_le_card_of_injective _
    (IsLocalRing.ResidueField.map (adjoinIntegersRingHom K hle)).injective

/-- **慣性次数は単調**。 -/
theorem inertiaDegree_le_of_adjoin_le (K : PAdicLocalField p) {x y : K.closure}
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    inertiaDegree K x ≤ inertiaDegree K y := by
  have h := card_residueField_le_of_adjoin_le K hle
  rw [residueDegree_eq_residueCard_pow K x, residueDegree_eq_residueCard_pow K y] at h
  have hq2 : 1 < Nat.card 𝓀[K.carrier] := by
    haveI : Fintype 𝓀[K.carrier] := Fintype.ofFinite _
    rw [Nat.card_eq_fintype_card]
    exact Fintype.one_lt_card
  exact (Nat.pow_le_pow_iff_right hq2).mp h

/-- **★★★★★完全分岐は部分拡大に遺伝する**。 -/
theorem isTotallyRamified_of_le (K : PAdicLocalField p) {x y : K.closure}
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (hty : IsTotallyRamifiedAdjoin K y) : IsTotallyRamifiedAdjoin K x := by
  have hle' := inertiaDegree_le_of_adjoin_le K hle
  rw [show inertiaDegree K y = 1 from hty] at hle'
  have hpos : 0 < inertiaDegree K x := by
    have h := ramificationIndex_mul_inertiaDegree K x
    have hn : 0 < Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := Module.finrank_pos
    rcases Nat.eq_zero_or_pos (inertiaDegree K x) with h0 | h0
    · rw [h0, Nat.mul_zero] at h; omega
    · exact h0
  show inertiaDegree K x = 1
  omega

/-- **★★★★★★「完全分岐 ∩ K^ur = K」**——完全分岐拡大の中の不分岐な部分は自明。
相互律の全体像(`Γ_K^ab ≅ (K^×)^∧`)に必要な線型無関係性。 -/
theorem finrank_eq_one_of_mem_unramifiedClosure_of_le (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (hx : x ∈ unramifiedClosure K) (hty : IsTotallyRamifiedAdjoin K y) :
    Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = 1 :=
  finrank_eq_one_of_isUnramified_of_isTotallyRamified K x
    ((mem_unramifiedClosure_iff_isUnramified K x).mp hx)
    (isTotallyRamified_of_le K hle hty)

end ABC3.Found.PGC
