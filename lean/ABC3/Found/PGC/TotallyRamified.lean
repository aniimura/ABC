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

つまり「1 層は速い・2 層は落ちる」。次に触るときは、
`adjoinIntegers` を経由せず `𝒪[(adjoinField K x).carrier]` 側で組む
(`AdjoinFieldConstruction.lean::integers_eq_adjoinIntegers` で移せる)か、
ノルム保存を明示補題として先に用意してから `def` を書くこと。
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

end ABC3.Found.PGC
