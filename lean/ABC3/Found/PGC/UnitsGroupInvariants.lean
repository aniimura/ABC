import ABC3.Found.PGC.PrimeToPTorsion
import ABC3.Found.PGC.PrincipalUnitsRank

/-!
# `q` と `[K:ℚ_p]` は**群の同型類だけ**から決まる

[pGC] Proposition 1.2 の論拠は、相互律 `Γ_K^ab ≅ (K^×)^∧` と、そこからの
**計数**である。本ファイルはその計数を、群の言葉だけで完結させる:

* **`residueCard_eq_of_units_mulEquiv`** :
  `𝒪_K^× ≃* 𝒪_{K'}^×` ⟹ `q_K = q_{K'}`
  ——`p` と素な捩れ部分群 `μ^{(p')}` は群論的に定義され、その位数が `q-1`
  (`PrimeToPTorsion.lean`)。
* **`finrank_eq_of_smallPrincipalUnits_mulEquiv`** :
  `1+𝔪_K ≃* 1+𝔪_{K'}`(十分小さい層)⟹ `[K:ℚ_p] = [K':ℚ_p]`
  ——主単数群は `ℤ_p^{[K:ℚ_p]}` そのもので、`p` 乗部分群の指数が
  `p^{[K:ℚ_p]}`。

## 主単数群 = `ℤ_p^{[K:ℚ_p]}`

`PrincipalUnitsLog.lean` の `padicLogUnitsEquiv`(p進対数)で
`1+𝔪^k ≃* Multiplicative(半径 1/4 の球)`、`PrincipalUnitsRank.lean` で
その球が `ℤ_p` 上の階数 `[K:ℚ_p]` の**自由**加群。基底を選べば

```
smallPrincipalUnits K  ≃*  Multiplicative (Fin [K:ℚ_p] → ℤ_[p])
```

(`smallPrincipalUnitsEquivPi`)。階数は群の言葉では `p` 乗部分群の指数
`p^{[K:ℚ_p]}` として読める(`index_powRange_smallPrincipalUnits`)
——`ℤ_p/pℤ_p ≅ ZMod p`(`PadicInt.ker_toZMod`)を成分ごとに使う。

## Proposition 1.2 との距離

`AbsGalUnitsSurjective.lean` で `Γ_K ↠ 𝒪_K^×` は**無条件**になった。
本ファイルで「`𝒪_K^×` の同型類 ⟹ `q`」「主単数群の同型類 ⟹ `[K:ℚ_p]`」も
確定した。残るのは**その商の標準性**——`Γ_K ≅ Γ_{K'}`(位相群として)から
`𝒪_K^× ≃* 𝒪_{K'}^×` を導く一手であり、それがまさに相互律
`Γ_K^ab ≅ (K^×)^∧` そのものである。Proposition 1.2 の未解決部分は、
これで**相互律ちょうど一つ**に絞り込まれた。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found Subgroup
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 主単数群 `≃* ℤ_p^{[K:ℚ_p]}` -/

/-- 台が同じなので、`AddSubgroup` 版と `Submodule` 版は加法群として同型。 -/
def smallBallAddEquiv (K : PAdicLocalField p) : smallBall K ≃+ smallBallSubmodule K where
  toFun y := ⟨y.1, y.2⟩
  invFun y := ⟨y.1, y.2⟩
  left_inv _ := Subtype.ext rfl
  right_inv _ := Subtype.ext rfl
  map_add' _ _ := Subtype.ext rfl

/-- 半径 `1/4` の球の `ℤ_p`-基底(自由・有限・階数 `[K:ℚ_p]`)。 -/
noncomputable def smallBallBasis (K : PAdicLocalField p) :
    Module.Basis (Fin (Module.finrank ℚ_[p] K.carrier)) ℤ_[p] (smallBallSubmodule K) := by
  haveI := module_finite_smallBall K
  haveI := module_free_smallBall K
  exact Module.finBasisOfFinrankEq _ _ (finrank_smallBall K)

/-- **★★★★主単数群は `ℤ_p^{[K:ℚ_p]}` そのもの**——p進対数(`padicLogUnitsEquiv`)で
加法群にし、`ℤ_p`-自由性(`module_free_smallBall`)で基底を取る。 -/
noncomputable def smallPrincipalUnitsEquivPi (K : PAdicLocalField p) :
    smallPrincipalUnits K ≃* Multiplicative (Fin (Module.finrank ℚ_[p] K.carrier) → ℤ_[p]) :=
  (padicLogUnitsEquiv K).trans
    (AddEquiv.toMultiplicative
      ((smallBallAddEquiv K).trans (smallBallBasis K).equivFun.toAddEquiv))

/-! ## `n` 乗部分群の指数——群論的な階数 -/

/-- 群同型は `n` 乗部分群を `n` 乗部分群に写す。 -/
theorem map_powRange_of_mulEquiv {G H : Type*} [CommGroup G] [CommGroup H] (e : G ≃* H) (n : ℕ) :
    Subgroup.map e.toMonoidHom ((powMonoidHom n : G →* G).range)
      = (powMonoidHom n : H →* H).range := by
  ext y
  constructor
  · rintro ⟨x, ⟨z, rfl⟩, rfl⟩
    exact ⟨e z, by simp [powMonoidHom]⟩
  · rintro ⟨w, rfl⟩
    exact ⟨(e.symm w) ^ n, ⟨e.symm w, rfl⟩, by simp [powMonoidHom]⟩

/-- したがって `n` 乗部分群の指数は群の同型不変量。 -/
theorem index_powRange_of_mulEquiv {G H : Type*} [CommGroup G] [CommGroup H] (e : G ≃* H) (n : ℕ) :
    ((powMonoidHom n : G →* G).range).index = ((powMonoidHom n : H →* H).range).index := by
  rw [← map_powRange_of_mulEquiv e n, Subgroup.index_map_of_bijective e.bijective]

/-- 成分ごとの `ℤ_p ↠ ZMod p`。 -/
noncomputable def redToZMod (p d : ℕ) [Fact p.Prime] : (Fin d → ℤ_[p]) →+ (Fin d → ZMod p) :=
  AddMonoidHom.mk' (fun x i => PadicInt.toZMod (x i)) (by intro a b; funext i; simp)

theorem redToZMod_surjective (p d : ℕ) [Fact p.Prime] : Function.Surjective (redToZMod p d) := by
  intro y
  refine ⟨fun i => ((y i).val : ℤ_[p]), ?_⟩
  funext i
  show PadicInt.toZMod (((y i).val : ℕ) : ℤ_[p]) = y i
  rw [map_natCast]
  exact ZMod.natCast_rightInverse _

/-- 核は `p·(ℤ_p^d)`(`PadicInt.ker_toZMod` + `maximalIdeal_eq_span_p` を成分ごとに)。 -/
theorem redToZMod_eq_zero_iff (p d : ℕ) [Fact p.Prime] (x : Fin d → ℤ_[p]) :
    redToZMod p d x = 0 ↔ ∃ y : Fin d → ℤ_[p], x = (p : ℕ) • y := by
  constructor
  · intro h
    have hi : ∀ i, PadicInt.toZMod (x i) = 0 := fun i => congrFun h i
    have hmem : ∀ i, x i ∈ Ideal.span {(p : ℤ_[p])} := by
      intro i
      rw [← PadicInt.maximalIdeal_eq_span_p, ← PadicInt.ker_toZMod]
      exact hi i
    choose c hc using fun i => Ideal.mem_span_singleton'.mp (hmem i)
    exact ⟨c, by funext i; rw [Pi.smul_apply, nsmul_eq_mul, ← hc i]; ring⟩
  · rintro ⟨y, rfl⟩
    funext i
    show PadicInt.toZMod (((p : ℕ) • y) i) = 0
    rw [Pi.smul_apply, nsmul_eq_mul, map_mul, map_natCast]
    simp

/-- **`ℤ_p^d` の `p` 乗部分群の指数は `p^d`**。 -/
theorem index_powRange_pi (p d : ℕ) [Fact p.Prime] :
    ((powMonoidHom p : Multiplicative (Fin d → ℤ_[p]) →* _).range).index = p ^ d := by
  set f : Multiplicative (Fin d → ℤ_[p]) →* Multiplicative (Fin d → ZMod p) :=
    AddMonoidHom.toMultiplicative (redToZMod p d) with hf
  have hfsurj : Function.Surjective f := redToZMod_surjective p d
  have hker : f.ker = (powMonoidHom p : Multiplicative (Fin d → ℤ_[p]) →* _).range := by
    ext x
    constructor
    · intro hx
      have h0 : redToZMod p d (Multiplicative.toAdd x) = 0 := hx
      obtain ⟨y, hy⟩ := (redToZMod_eq_zero_iff p d _).mp h0
      refine ⟨Multiplicative.ofAdd y, Multiplicative.toAdd.injective ?_⟩
      show p • y = Multiplicative.toAdd x
      exact hy.symm
    · rintro ⟨y, rfl⟩
      show redToZMod p d (Multiplicative.toAdd ((powMonoidHom p) y)) = 0
      exact (redToZMod_eq_zero_iff p d _).mpr ⟨Multiplicative.toAdd y, rfl⟩
  rw [← hker]
  have hcard : Nat.card (Multiplicative (Fin d → ℤ_[p]) ⧸ f.ker)
      = Nat.card (Multiplicative (Fin d → ZMod p)) :=
    Nat.card_congr (QuotientGroup.quotientKerEquivOfSurjective f hfsurj).toEquiv
  show Nat.card (Multiplicative (Fin d → ℤ_[p]) ⧸ f.ker) = p ^ d
  rw [hcard]
  show Nat.card (Fin d → ZMod p) = p ^ d
  haveI : NeZero p := ⟨(Fact.out : p.Prime).ne_zero⟩
  simp [Nat.card_eq_fintype_card]

/-- **★★★★`[K:ℚ_p]` は主単数群から群論的に読める**——`p` 乗部分群の指数が
`p^{[K:ℚ_p]}`。 -/
theorem index_powRange_smallPrincipalUnits (K : PAdicLocalField p) :
    ((powMonoidHom p : smallPrincipalUnits K →* _).range).index
      = p ^ Module.finrank ℚ_[p] K.carrier := by
  rw [index_powRange_of_mulEquiv (smallPrincipalUnitsEquivPi K) p]
  exact index_powRange_pi p _

/-! ## 群の同型類だけから決まる二つの量 -/

/-- **★★★★★`[K:ℚ_p]` は主単数群の同型類だけで決まる**。 -/
theorem finrank_eq_of_smallPrincipalUnits_mulEquiv {K K' : PAdicLocalField p}
    (e : smallPrincipalUnits K ≃* smallPrincipalUnits K') :
    Module.finrank ℚ_[p] K.carrier = Module.finrank ℚ_[p] K'.carrier := by
  have h := index_powRange_of_mulEquiv e p
  rw [index_powRange_smallPrincipalUnits, index_powRange_smallPrincipalUnits] at h
  exact Nat.pow_right_injective (Fact.out : p.Prime).two_le h

/-- `μ^{(p')}` は群論的な定義なので群同型で移る。 -/
theorem map_primeToPTorsion {K K' : PAdicLocalField p}
    (e : (𝒪[K.carrier])ˣ ≃* (𝒪[K'.carrier])ˣ) :
    Subgroup.map e.toMonoidHom (primeToPTorsion K) = primeToPTorsion K' := by
  ext y
  constructor
  · rintro ⟨x, ⟨m, hm, hx⟩, rfl⟩
    exact ⟨m, hm, by rw [← map_pow, hx, map_one]⟩
  · rintro ⟨m, hm, hy⟩
    exact ⟨e.symm y, ⟨m, hm, by rw [← map_pow, hy, map_one]⟩, by simp⟩

/-- **★★★★★`q` は `𝒪_K^×` の同型類だけで決まる**。 -/
theorem residueCard_eq_of_units_mulEquiv {K K' : PAdicLocalField p}
    (e : (𝒪[K.carrier])ˣ ≃* (𝒪[K'.carrier])ˣ) :
    Nat.card 𝓀[K.carrier] = Nat.card 𝓀[K'.carrier] := by
  have h1 : primeToPTorsion K ≃* (primeToPTorsion K).map e.toMonoidHom :=
    Subgroup.equivMapOfInjective _ _ e.injective
  have hcard : Nat.card (primeToPTorsion K) = Nat.card (primeToPTorsion K') := by
    rw [Nat.card_congr h1.toEquiv, map_primeToPTorsion e]
  rw [card_primeToPTorsion, card_primeToPTorsion] at hcard
  have h2 : 1 < Nat.card 𝓀[K.carrier] := by
    rw [Nat.card_eq_fintype_card]; exact Fintype.one_lt_card
  have h3 : 1 < Nat.card 𝓀[K'.carrier] := by
    rw [Nat.card_eq_fintype_card]; exact Fintype.one_lt_card
  omega

end ABC3.Found.PGC
