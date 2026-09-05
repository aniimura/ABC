import ABC3.Found.PGC.PadicLogSurjective

/-!
# 主単数群は p 進対数で加法群になる——`1+𝔪^k ≅ (𝔪^k, +)`

`Found/PGC/PadicLogSurjective.lean::padicLog_bijOn`(`padicLog K` は半径 `1/4` の
球からそれ自身への全単射)と `PadicLogMul.lean::padicLog_mul`
(`log((1+x)(1+y)) = log(1+x) + log(1+y)`)を合わせると、

```
{u ∈ K^× | ‖u-1‖ ≤ 1/4}  ≃*  Multiplicative {y ∈ K | ‖y‖ ≤ 1/4}
```

という**群同型**が得られる。左辺は主単数群(の十分小さい層)、右辺は
`(K,+)` の開コンパクト部分群。

## なぜ要るか

`Found/PGC/UnramifiedExtension.lean::nonempty_units_mulEquiv_prod` で
`𝒪_K^× ≅ 𝓀^× × (1+𝔪_K)` を出した。その**第二因子**の構造がこれ:
`1+𝔪_K` は(十分小さい層で)加法群 `𝔪_K` と同型で、`𝔪_K` は
`ℤ_p`-加群として階数 `[K:ℚ_p]` の自由加群。[pGC] Proposition 1.2 が
`Γ_K^ab` の群構造から `[K:ℚ_p]` を読み取るのはこの経路。

第一因子(`𝓀^×`、位数 `q-1` の巡回群)からは `q` が読める
(`Found/Teichmuller.lean`)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- 半径 `1/4` の球——超距離不等式により `(K,+)` の部分群。 -/
def smallBall (K : PAdicLocalField p) : AddSubgroup K.carrier where
  carrier := {y : K.carrier | ‖y‖ ≤ 1 / 4}
  zero_mem' := by simp
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at *
    exact le_trans (IsUltrametricDist.norm_add_le_max a b) (max_le ha hb)
  neg_mem' := by
    intro a ha
    simpa using ha

/-- `‖u-1‖ < 1` なら `‖u‖ = 1`(超距離不等式)。 -/
theorem norm_eq_one_of_norm_sub_one_lt (K : PAdicLocalField p) {u : K.carrier}
    (hu : ‖u - 1‖ < 1) : ‖u‖ = 1 := by
  have h : u = 1 + (u - 1) := by ring
  rw [h, IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (S := K.carrier)]
  · rw [norm_one]; exact max_eq_left (le_of_lt hu)
  · rw [norm_one]; exact ne_of_gt hu

/-- 主単数群(の `‖u-1‖ ≤ 1/4` の層)——`K^×` の部分群。 -/
def smallPrincipalUnits (K : PAdicLocalField p) : Subgroup (K.carrier)ˣ where
  carrier := {u : (K.carrier)ˣ | ‖(u : K.carrier) - 1‖ ≤ 1 / 4}
  one_mem' := by simp
  mul_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at *
    have hexp : ((a * b : (K.carrier)ˣ) : K.carrier) - 1
        = ((a : K.carrier) - 1) + ((b : K.carrier) - 1)
          + ((a : K.carrier) - 1) * ((b : K.carrier) - 1) := by
      push_cast; ring
    rw [hexp]
    refine le_trans (IsUltrametricDist.norm_add_le_max _ _) (max_le ?_ ?_)
    · exact le_trans (IsUltrametricDist.norm_add_le_max _ _) (max_le ha hb)
    · rw [norm_mul]
      calc ‖(a : K.carrier) - 1‖ * ‖(b : K.carrier) - 1‖
          ≤ (1 / 4) * (1 / 4) := mul_le_mul ha hb (norm_nonneg _) (by norm_num)
        _ ≤ 1 / 4 := by norm_num
  inv_mem' := by
    intro a ha
    simp only [Set.mem_setOf_eq] at *
    have hne : ‖(a : K.carrier)‖ = 1 :=
      norm_eq_one_of_norm_sub_one_lt K (lt_of_le_of_lt ha (by norm_num))
    have h1 : ((a : K.carrier)) * ((a⁻¹ : (K.carrier)ˣ) : K.carrier) = 1 := by
      rw [← Units.val_mul, mul_inv_cancel, Units.val_one]
    have hexp : ((a⁻¹ : (K.carrier)ˣ) : K.carrier) - 1
        = -(((a : K.carrier) - 1) * ((a⁻¹ : (K.carrier)ˣ) : K.carrier)) := by
      have h2 : ((a : K.carrier) - 1) * ((a⁻¹ : (K.carrier)ˣ) : K.carrier)
          = 1 - ((a⁻¹ : (K.carrier)ˣ) : K.carrier) := by
        rw [sub_mul, one_mul, h1]
      rw [h2]; ring
    rw [hexp, norm_neg, norm_mul]
    have hinv : ‖((a⁻¹ : (K.carrier)ˣ) : K.carrier)‖ = 1 := by
      have h3 := congrArg norm h1
      rw [norm_mul, hne, one_mul, norm_one] at h3
      exact h3
    rw [hinv, mul_one]
    exact ha

/-- 主単数の対数は球に入る。 -/
theorem mem_smallBall_log (K : PAdicLocalField p) (u : smallPrincipalUnits K) :
    padicLog K (((u : (K.carrier)ˣ) : K.carrier) - 1) ∈ smallBall K := by
  have hu : ‖((u : (K.carrier)ˣ) : K.carrier) - 1‖ ≤ 1 / 4 := u.2
  exact (padicLog_bijOn K).mapsTo hu

/-- **p 進対数が与える群準同型** `1+𝔪^k →* (𝔪^k, +)`。 -/
noncomputable def padicLogUnitsHom (K : PAdicLocalField p) :
    smallPrincipalUnits K →* Multiplicative (smallBall K) where
  toFun u := Multiplicative.ofAdd
    (⟨padicLog K (((u : (K.carrier)ˣ) : K.carrier) - 1), mem_smallBall_log K u⟩ : smallBall K)
  map_one' := by
    apply congrArg Multiplicative.ofAdd
    apply Subtype.ext
    show padicLog K ((((1 : smallPrincipalUnits K) : (K.carrier)ˣ) : K.carrier) - 1) = 0
    have h : (((1 : smallPrincipalUnits K) : (K.carrier)ˣ) : K.carrier) - 1 = 0 := by simp
    rw [h]
    have h2 := padicLog_mul (K := K) (x := 0) (y := 0) (by simp) (by simp)
    simpa using h2
  map_mul' a b := by
    apply congrArg Multiplicative.ofAdd
    apply Subtype.ext
    show padicLog K ((((a * b : smallPrincipalUnits K) : (K.carrier)ˣ) : K.carrier) - 1)
      = padicLog K (((a : (K.carrier)ˣ) : K.carrier) - 1)
        + padicLog K (((b : (K.carrier)ˣ) : K.carrier) - 1)
    have ha : ‖((a : (K.carrier)ˣ) : K.carrier) - 1‖ ≤ 1 / 4 := a.2
    have hb : ‖((b : (K.carrier)ˣ) : K.carrier) - 1‖ ≤ 1 / 4 := b.2
    have hexp : (((a * b : smallPrincipalUnits K) : (K.carrier)ˣ) : K.carrier) - 1
        = (((a : (K.carrier)ˣ) : K.carrier) - 1) + (((b : (K.carrier)ˣ) : K.carrier) - 1)
          + (((a : (K.carrier)ˣ) : K.carrier) - 1) * (((b : (K.carrier)ˣ) : K.carrier) - 1) := by
      push_cast; ring
    rw [hexp]
    exact padicLog_mul (lt_of_le_of_lt ha (by norm_num)) (lt_of_le_of_lt hb (by norm_num))

/-- **全単射**——単射は `padicLog_injOn`、全射は `padicLog_surjOn`
(`1+x` が単数であることは `‖1+x‖ = 1 ≠ 0` から)。 -/
theorem padicLogUnitsHom_bijective (K : PAdicLocalField p) :
    Function.Bijective (padicLogUnitsHom K) := by
  constructor
  · intro a b hab
    have h : padicLog K (((a : (K.carrier)ˣ) : K.carrier) - 1)
        = padicLog K (((b : (K.carrier)ˣ) : K.carrier) - 1) := congrArg Subtype.val hab
    have ha : ((a : (K.carrier)ˣ) : K.carrier) - 1 ∈ {x : K.carrier | ‖x‖ ≤ 1 / 4} := a.2
    have hb : ((b : (K.carrier)ˣ) : K.carrier) - 1 ∈ {x : K.carrier | ‖x‖ ≤ 1 / 4} := b.2
    have hsub := padicLog_injOn K ha hb h
    apply Subtype.ext
    apply Units.ext
    have hval : ((a : (K.carrier)ˣ) : K.carrier) = ((b : (K.carrier)ˣ) : K.carrier) := by
      linear_combination hsub
    exact hval
  · intro w
    obtain ⟨x, hx, hlog⟩ := (padicLog_bijOn K).surjOn w.2
    have hxlt : ‖x‖ < 1 := lt_of_le_of_lt hx (by norm_num)
    have hne : (1 : K.carrier) + x ≠ 0 := by
      have h1 : ‖(1 : K.carrier) + x‖ = 1 := by
        rw [IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (S := K.carrier)]
        · rw [norm_one]; exact max_eq_left (le_of_lt hxlt)
        · rw [norm_one]; exact ne_of_gt hxlt
      intro hcon
      rw [hcon, norm_zero] at h1
      exact zero_ne_one h1
    have hmem : (Units.mk0 ((1 : K.carrier) + x) hne) ∈ smallPrincipalUnits K := by
      show ‖((Units.mk0 ((1 : K.carrier) + x) hne : (K.carrier)ˣ) : K.carrier) - 1‖ ≤ 1 / 4
      simpa using hx
    refine ⟨⟨Units.mk0 ((1 : K.carrier) + x) hne, hmem⟩, ?_⟩
    apply Subtype.ext
    have hs : ((Units.mk0 ((1 : K.carrier) + x) hne : (K.carrier)ˣ) : K.carrier) - 1 = x := by
      simp
    show padicLog K (((Units.mk0 ((1 : K.carrier) + x) hne : (K.carrier)ˣ) : K.carrier) - 1) = _
    rw [hs]
    exact hlog

/-- **★★主単数群は加法群と同型** `1+𝔪^k ≃* Multiplicative (𝔪^k)`。 -/
noncomputable def padicLogUnitsEquiv (K : PAdicLocalField p) :
    smallPrincipalUnits K ≃* Multiplicative (smallBall K) :=
  MulEquiv.ofBijective (padicLogUnitsHom K) (padicLogUnitsHom_bijective K)

end ABC3.Found.PGC
