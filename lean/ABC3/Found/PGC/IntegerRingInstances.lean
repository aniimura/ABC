import ABC3.Found.PGC.IntegerSubringNorm

/-!
# 骨組み(作業中) —— `𝒪_M` に型クラスを載せる
-/

namespace ABC3.Found.PGC

namespace IntegerNorm

section Units

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★`𝒪_M` の単元 ⟺ ノルムが `1`。 -/
theorem isUnit_iff_norm_eq_one {x : ↥(integerSubring M)} : IsUnit x ↔ ‖(x : M)‖ = 1 := by
  constructor
  · rintro ⟨u, hu⟩
    have h1 : ((u : ↥(integerSubring M)) : M) * ((u⁻¹ : (integerSubring M)ˣ) : M) = 1 := by
      have := u.mul_inv
      exact congrArg (fun y : ↥(integerSubring M) => (y : M)) this
    have hx : ((x : M)) = ((u : ↥(integerSubring M)) : M) := by rw [hu]
    have hnorm : ‖((u : ↥(integerSubring M)) : M)‖ * ‖((u⁻¹ : (integerSubring M)ˣ) : M)‖ = 1 := by
      rw [← norm_mul, h1, norm_one]
    have hle1 : ‖((u : ↥(integerSubring M)) : M)‖ ≤ 1 := (u : ↥(integerSubring M)).2
    have hle2 : ‖((u⁻¹ : (integerSubring M)ˣ) : M)‖ ≤ 1 := (u⁻¹ : (integerSubring M)ˣ).1.2
    rw [hx]
    nlinarith [norm_nonneg ((u : ↥(integerSubring M)) : M),
      norm_nonneg ((u⁻¹ : (integerSubring M)ˣ) : M)]
  · intro hx
    have hx0 : (x : M) ≠ 0 := by
      intro h
      rw [h, norm_zero] at hx
      norm_num at hx
    have hmem : ((x : M))⁻¹ ∈ integerSubring M := by
      rw [mem_integerSubring, norm_inv, hx, inv_one]
    refine isUnit_iff_exists_inv.mpr ⟨⟨((x : M))⁻¹, hmem⟩, ?_⟩
    ext
    simpa using mul_inv_cancel₀ hx0

/-- ★★`𝒪_M` は局所環(超距離だけ)。 -/
theorem isLocalRing_integerSubring : IsLocalRing ↥(integerSubring M) := by
  haveI : Nontrivial ↥(integerSubring M) := inferInstance
  refine IsLocalRing.of_isUnit_or_isUnit_one_sub_self (fun a => ?_)
  have ha1 : ‖(a : M)‖ ≤ 1 := a.2
  rcases eq_or_lt_of_le ha1 with h1 | h1
  · exact Or.inl (isUnit_iff_norm_eq_one.mpr h1)
  · refine Or.inr (isUnit_iff_norm_eq_one.mpr ?_)
    have hco : ((1 - a : ↥(integerSubring M)) : M) = 1 - (a : M) := by push_cast; ring
    rw [hco]
    have hne : ‖(1 : M)‖ ≠ ‖-(a : M)‖ := by
      rw [norm_one, norm_neg]
      exact ne_of_gt h1
    have hmax := IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (x := (1 : M))
      (y := -(a : M)) hne
    rw [← sub_eq_add_neg] at hmax
    rw [hmax, norm_one, norm_neg]
    exact max_eq_left (le_of_lt h1)


end Units

/-! ## §2 作用の制限 -/

section Action

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★★**等長な環作用は `𝒪_M` に制限できる** —— `MulSemiringAction G ↥(integerSubring M)`。

★木の流儀(`lean-idioms.md` #165)に合わせ **`instance` にせず `def`** にして `letI` で貼る
(`hiso` を仮説で受けるので `instance` にはできない)。 -/
@[implicit_reducible] def integerMulSemiringAction {G : Type*} [Monoid G] [MulSemiringAction G M]
    (hiso : ∀ (g : G) (z : M), ‖g • z‖ = ‖z‖) :
    MulSemiringAction G ↥(integerSubring M) where
  smul g x := ⟨g • (x : M), by rw [mem_integerSubring, hiso]; exact x.2⟩
  one_smul x := Subtype.ext (one_smul G (x : M))
  mul_smul g h x := Subtype.ext (mul_smul g h (x : M))
  smul_zero g := Subtype.ext (smul_zero g)
  smul_add g x y := Subtype.ext (smul_add g (x : M) (y : M))
  smul_one g := Subtype.ext (smul_one g)
  smul_mul g x y := Subtype.ext (smul_mul' g (x : M) (y : M))

end Action

end IntegerNorm

end ABC3.Found.PGC
