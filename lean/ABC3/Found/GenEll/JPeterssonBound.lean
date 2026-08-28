/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.PeterssonBound
import ABC3.Found.GenEll.JSurjective
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★`‖j‖·‖Δ‖^{1+ϵ}` は有界（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★★★★★これは何か —— `Prop 3.4` の**第 3 の `≲`**

`§9-980` で `Prop 3.4` の第 1・第 2 の `≲` の合成
（`deg∞ ≲ 12(1+ϵ)·ht^Falt`）は `archSum` の**上界だけ**で取れた。
☆残っていたのは**第 3 の `≲`**（`12(1+ϵ)·ht^Falt ≲ (1+ϵ)·ht∞`）である。

★★★`§9-1000`（第 548）の測定で、それは

    **`h(j) ≤ 12(1+ϵ)·ht^Falt + C`**

に、さらに素点ごとに

    **`log⁺|j(τ)| ≤ (1+ϵ)·log(1/((2π)¹²‖Δ‖_Pet(τ))) + C′`**

まで降ろされた。★★これは指数を払えば

    **`‖j(τ)‖ · (‖Δ‖_Pet(τ))^{1+ϵ}` が有界**

と同じことである。★★★★**本ファイルはそれを証明する。**

## ★機構 —— `j = E₄³/Δ` を使うと綺麗に分かれる

`jFun τ = E₄(τ)³/Δ(τ)`、`peterssonDelta τ = ‖Δ(τ)‖·(Im τ)⁶` なので

    `‖j(τ)‖ · peterssonDelta(τ)^{1+ϵ} = ‖E₄(τ)‖³ · (Im τ)⁶ · peterssonDelta(τ)^ϵ`

★`E₄` は**カスプで有界**（保型形式だから）——`‖E₄‖³` は `O(1)`。
★★`(Im τ)⁶ · peterssonDelta^ϵ → 0`——`Δ` はカスプ形式なので
`peterssonDelta = O(e^{−2πy}y⁶)`、したがって
`y⁶·(C e^{−2πy}y⁶)^ϵ = C^ϵ y^{6+6ϵ} e^{−2πϵy} → 0`。
★★★あとは `§9-355` と**同じ型**——切り詰めた基本領域のコンパクト性と `SL(2,ℤ)` 不変性。

## ★★これで `Prop 3.4` に残るのは

☆`ht∞` の同定（`M̄_ell` の無限遠因子の高さが `h(j)` であること）と、
`h(j)` の有限素点の寄与が `deg∞` であること（半安定のとき `v(j) = −v(Δ_min)`）。
★★どちらも**代数の段**であり、解析はここで尽きた。
-/

namespace ABC3.Found.GenEll

open Complex Real ModularForm MatrixGroups UpperHalfPlane Filter Topology Asymptotics
open ABC3.Found.GaloisRep

/-! ## ★★★★★`Δ` の減衰（`§9-355` の中身を補題に出す） -/

/-- ★★★★★**`peterssonDelta = O(e^{−2πy}·y⁶)`**（`Δ` がカスプ形式であることから）。

★`§9-355`（`tendsto_peterssonDelta_atImInfty`）の証明の中にあったものを、
**係数つきで使えるよう**補題に出した。 -/
theorem isBigO_peterssonDelta :
    peterssonDelta =O[atImInfty] (fun τ : ℍ => Real.exp (-2 * Real.pi * τ.im / 1) * τ.im ^ 6) := by
  have hO := CuspFormClass.exp_decay_atImInfty (CuspForm.discriminant) (h := (1:ℝ))
    one_pos one_mem_strictPeriods_SL
  rw [CuspForm.coe_discriminant] at hO
  exact hO.norm_left.mul (Asymptotics.isBigO_refl (fun τ : ℍ => τ.im ^ 6) atImInfty)

/-! ## ★★★★★★★★★カスプでの極限 -/

/-- ★★★★★★★★★**`(Im τ)⁶ · peterssonDelta^ϵ → 0`**（カスプ）。

★`peterssonDelta = O(e^{−2πy}y⁶)` から
`y⁶·(C e^{−2πy}y⁶)^ϵ = C^ϵ·y^{6+6ϵ}·e^{−2πϵy} → 0`。
★★指数減衰が多項式に勝つ（`tendsto_rpow_mul_exp_neg_mul_atTop_nhds_zero`）。 -/
theorem tendsto_im_pow_mul_pD_rpow (eps : ℝ) (heps : 0 < eps) :
    Tendsto (fun τ : ℍ => τ.im ^ 6 * peterssonDelta τ ^ eps) atImInfty (𝓝 0) := by
  obtain ⟨C, hC⟩ := isBigO_peterssonDelta.bound
  set C' : ℝ := max C 0 with hC'
  have hC'0 : (0:ℝ) ≤ C' := le_max_right _ _
  have hbd : ∀ᶠ τ : ℍ in atImInfty,
      peterssonDelta τ ≤ C' * (Real.exp (-2 * Real.pi * τ.im / 1) * τ.im ^ 6) := by
    filter_upwards [hC] with τ hτ
    have hpos : (0:ℝ) ≤ Real.exp (-2 * Real.pi * τ.im / 1) * τ.im ^ 6 := by positivity
    rw [Real.norm_of_nonneg (peterssonDelta_pos τ).le, Real.norm_of_nonneg hpos] at hτ
    exact le_trans hτ (mul_le_mul_of_nonneg_right (le_max_left _ _) hpos)
  have hgoal : ∀ᶠ τ : ℍ in atImInfty,
      τ.im ^ 6 * peterssonDelta τ ^ eps
        ≤ C' ^ eps * (τ.im ^ (6 + 6 * eps) * Real.exp (-(2 * Real.pi * eps) * τ.im)) := by
    filter_upwards [hbd] with τ hτ
    have hy : (0:ℝ) < τ.im := τ.2
    have hpow : (τ.im : ℝ) ^ (6:ℕ) * τ.im ^ (6 * eps) = τ.im ^ (6 + 6 * eps) := by
      rw [← Real.rpow_natCast τ.im 6, ← Real.rpow_add hy]
      norm_num
    have h1 : peterssonDelta τ ^ eps
        ≤ (C' * (Real.exp (-2 * Real.pi * τ.im / 1) * τ.im ^ 6)) ^ eps :=
      Real.rpow_le_rpow (peterssonDelta_pos τ).le hτ heps.le
    have h2 : (C' * (Real.exp (-2 * Real.pi * τ.im / 1) * τ.im ^ 6)) ^ eps
        = C' ^ eps * (τ.im ^ (6 * eps) * Real.exp (-(2 * Real.pi * eps) * τ.im)) := by
      rw [Real.mul_rpow hC'0 (by positivity), Real.mul_rpow (by positivity) (by positivity)]
      rw [← Real.rpow_natCast τ.im 6, ← Real.rpow_mul hy.le]
      rw [Real.rpow_def_of_pos (Real.exp_pos _), Real.log_exp]
      push_cast
      ring_nf
    calc τ.im ^ 6 * peterssonDelta τ ^ eps
        ≤ τ.im ^ 6 * (C' * (Real.exp (-2 * Real.pi * τ.im / 1) * τ.im ^ 6)) ^ eps :=
          mul_le_mul_of_nonneg_left h1 (by positivity)
      _ = C' ^ eps * ((τ.im ^ (6:ℕ) * τ.im ^ (6 * eps))
            * Real.exp (-(2 * Real.pi * eps) * τ.im)) := by rw [h2]; ring
      _ = C' ^ eps * (τ.im ^ (6 + 6 * eps) * Real.exp (-(2 * Real.pi * eps) * τ.im)) := by
          rw [hpow]
  have hlim : Tendsto
      (fun τ : ℍ => C' ^ eps * (τ.im ^ (6 + 6 * eps) * Real.exp (-(2 * Real.pi * eps) * τ.im)))
      atImInfty (𝓝 0) := by
    have h := tendsto_rpow_mul_exp_neg_mul_atTop_nhds_zero (6 + 6 * eps) (2 * Real.pi * eps)
      (by positivity)
    have h2 := h.const_mul (C' ^ eps)
    rw [mul_zero] at h2
    exact h2.comp tendsto_comap
  refine squeeze_zero' ?_ hgoal hlim
  filter_upwards with τ
  have := peterssonDelta_pos τ
  have hy : (0:ℝ) < τ.im := τ.2
  positivity

/-! ## ★★★★★★`j = E₄³/Δ` による分解 -/

/-- ★★★★★★**`‖j‖·pD^{1+ϵ} = ‖E₄‖³·(Im τ)⁶·pD^ϵ`**。

★`jFun = E₄³/Δ` と `pD = ‖Δ‖·(Im τ)⁶` を代入するだけである。
★★これで「`E₄` の有界性」と「`pD` の減衰」に分かれる。 -/
theorem norm_jFun_mul_peterssonDelta_rpow (eps : ℝ) (τ : ℍ) :
    ‖jFun τ‖ * peterssonDelta τ ^ (1 + eps)
      = ‖ModularForm.E₄ τ‖ ^ 3 * τ.im ^ 6 * peterssonDelta τ ^ eps := by
  have hΔ : ‖ModularForm.discriminant τ‖ ≠ 0 :=
    norm_ne_zero_iff.2 (ModularForm.discriminant_ne_zero τ)
  have hpd : (0:ℝ) < peterssonDelta τ := peterssonDelta_pos τ
  have hsplit : peterssonDelta τ ^ (1 + eps) = peterssonDelta τ * peterssonDelta τ ^ eps := by
    rw [Real.rpow_add hpd, Real.rpow_one]
  have hj : ‖jFun τ‖ = ‖ModularForm.E₄ τ‖ ^ 3 / ‖ModularForm.discriminant τ‖ := by
    rw [jFun, norm_div, norm_pow]
  rw [hsplit, hj, peterssonDelta]
  field_simp

/-! ## ★★★★★★★★★★★★カスプで `0` -/

/-- ★★★★★★★★★★★★**`‖j‖·pD^{1+ϵ} → 0`**（カスプ）。

★`‖E₄‖³` はカスプで有界（保型形式）、`(Im τ)⁶·pD^ϵ → 0`。 -/
theorem tendsto_jFun_peterssonDelta (eps : ℝ) (heps : 0 < eps) :
    Tendsto (fun τ : ℍ => ‖jFun τ‖ * peterssonDelta τ ^ (1 + eps)) atImInfty (𝓝 0) := by
  have hE4 : (fun τ : ℍ => ‖ModularForm.E₄ τ‖ ^ 3) =O[atImInfty] (fun _ : ℍ => (1:ℝ)) := by
    have h : (fun τ : ℍ => ‖ModularForm.E₄ τ‖) =O[atImInfty] (fun _ : ℍ => (1:ℝ)) := by
      have hb := ModularFormClass.bdd_at_infty (ModularForm.E₄)
      rw [UpperHalfPlane.IsBoundedAtImInfty] at hb
      exact hb.norm_left
    simpa using h.pow 3
  have hlim := tendsto_im_pow_mul_pD_rpow eps heps
  have hO : (fun τ : ℍ => ‖ModularForm.E₄ τ‖ ^ 3 * (τ.im ^ 6 * peterssonDelta τ ^ eps))
      =O[atImInfty] (fun τ : ℍ => τ.im ^ 6 * peterssonDelta τ ^ eps) := by
    simpa using hE4.mul (isBigO_refl (fun τ : ℍ => τ.im ^ 6 * peterssonDelta τ ^ eps) atImInfty)
  refine (hO.trans_tendsto hlim).congr (fun τ => ?_)
  rw [norm_jFun_mul_peterssonDelta_rpow eps τ, mul_assoc]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★有界性 -/

/-- ★連続性——`jFun = E₄³/Δ` は `Δ ≠ 0` なので連続、`pD^{1+ϵ}` も正の上で連続。 -/
theorem continuous_normJ_pD (eps : ℝ) :
    Continuous (fun τ : ℍ => ‖jFun τ‖ * peterssonDelta τ ^ (1 + eps)) := by
  refine continuous_jFun.norm.mul ?_
  have h : ContinuousOn (fun x : ℝ => x ^ (1 + eps)) (Set.Ioi 0) :=
    ContinuousOn.rpow_const continuousOn_id (fun x hx => Or.inl (ne_of_gt hx))
  exact ContinuousOn.comp_continuous h continuous_peterssonDelta (fun τ => peterssonDelta_pos τ)

/-- ★`‖j‖·pD^{1+ϵ}` は `SL(2,ℤ)` 不変。 -/
theorem normJ_pD_smul (eps : ℝ) (g : Matrix.SpecialLinearGroup (Fin 2) ℤ) (τ : ℍ) :
    ‖jFun (g • τ)‖ * peterssonDelta (g • τ) ^ (1 + eps)
      = ‖jFun τ‖ * peterssonDelta τ ^ (1 + eps) := by
  rw [jFun_smul, peterssonDelta_smul]

/-- ★★★★★★★★★★★★★★★★★★★★★★**`‖j(τ)‖·(‖Δ‖_Pet(τ))^{1+ϵ}` は有界**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★これが `Proposition 3.4` の**第 3 の `≲`**（`12(1+ϵ)·ht^Falt ≲ (1+ϵ)·ht∞`）の
アルキメデス側の中身である——原文が [Silv2], `Proposition 2.1` に帰す部分。

★段取りは `§9-355`（`exists_bound_peterssonDelta`）と**同じ型**:
カスプで `0`（`tendsto_jFun_peterssonDelta`）＋
切り詰めた基本領域のコンパクト性 ＋ `SL(2,ℤ)` 不変性。 -/
theorem exists_bound_jFun_peterssonDelta (eps : ℝ) (heps : 0 < eps) :
    ∃ M : ℝ, 0 < M ∧ ∀ τ : ℍ, ‖jFun τ‖ * peterssonDelta τ ^ (1 + eps) ≤ M := by
  obtain ⟨y₀, hy₀⟩ := (UpperHalfPlane.atImInfty_mem _).1
    ((tendsto_jFun_peterssonDelta eps heps).eventually_lt_const (zero_lt_one))
  set y₁ := max y₀ 1 with hy₁
  have hcpt := ModularGroup.isCompact_truncatedFundamentalDomain y₁
  have hne : (ModularGroup.truncatedFundamentalDomain y₁).Nonempty := by
    refine ⟨UpperHalfPlane.I, ModularGroup.I_mem_fd, ?_⟩
    simp [hy₁]
  obtain ⟨x₀, hx₀mem, hx₀max⟩ :=
    hcpt.exists_isMaxOn hne (continuous_normJ_pD eps).continuousOn
  refine ⟨max (‖jFun x₀‖ * peterssonDelta x₀ ^ (1 + eps)) 1,
    lt_of_lt_of_le zero_lt_one (le_max_right _ _), fun τ => ?_⟩
  obtain ⟨g, hg⟩ := ModularGroup.exists_smul_mem_fd τ
  have hinv := normJ_pD_smul eps g τ
  rcases le_or_gt ((g • τ).im) y₁ with hle | hgt
  · have hm : ‖jFun (g • τ)‖ * peterssonDelta (g • τ) ^ (1 + eps)
        ≤ ‖jFun x₀‖ * peterssonDelta x₀ ^ (1 + eps) := hx₀max ⟨hg, hle⟩
    rw [hinv] at hm
    exact le_trans hm (le_max_left _ _)
  · have h1 : ‖jFun (g • τ)‖ * peterssonDelta (g • τ) ^ (1 + eps) < 1 :=
      hy₀ (g • τ) (le_of_lt (lt_of_le_of_lt (le_max_left y₀ 1) hgt))
    rw [hinv] at h1
    exact le_trans h1.le (le_max_right _ _)

/-! ## ★出典の紐付け(`.src`) -/

def isBigO_peterssonDelta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(peterssonDelta = O(e^{−2πy}y⁶))",
    sectionId := "genell-prop-3-4" }

def tendsto_jFun_peterssonDelta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(‖j‖·‖Δ‖^{1+ϵ} はカスプで 0)",
    sectionId := "genell-prop-3-4" }

def exists_bound_jFun_peterssonDelta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(第 3 の ≲ のアルキメデス側——‖j‖·‖Δ‖^{1+ϵ} は有界)",
    sectionId := "genell-prop-3-4" }

def exists_bound_jFun_peterssonDelta.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "CuspFormClass.exp_decay_atImInfty(カスプ形式の指数減衰)"
      (.inMathlib "CuspFormClass.exp_decay_atImInfty") 3,
    .citation "[mathlib]" "ModularFormClass.bdd_at_infty(保型形式はカスプで有界)"
      (.inMathlib "ModularFormClass.bdd_at_infty") 2,
    .citation "[mathlib]"
      "tendsto_rpow_mul_exp_neg_mul_atTop_nhds_zero(指数減衰が多項式に勝つ)"
      (.inMathlib "tendsto_rpow_mul_exp_neg_mul_atTop_nhds_zero") 2,
    .citation "[mathlib]"
      "ModularGroup.isCompact_truncatedFundamentalDomain(切り詰めた基本領域のコンパクト性)"
      (.inMathlib "ModularGroup.isCompact_truncatedFundamentalDomain") 3,
    .otherPaper "[Silv2]"
      ("Proposition 2.1——★原文が Prop 3.4 の第 2・第 3 の ≲ の根拠として引く。" ++
       "★★★**そのアルキメデス側の中身が本ファイルである**: " ++
       "‖j(τ)‖·(‖Δ‖_Pet(τ))^{1+ϵ} が有界であること") 9,
    .implicitStep
      ("★★★★★★★★到達点(2026-08-29): §9-980 で第 1・第 2 の ≲ の合成が取れ、" ++
       "第 548 の測定で第 3 の ≲ が『log⁺|j| ≤ (1+ϵ)log(1/((2π)¹²‖Δ‖)) + C』" ++
       "まで降ろされていた。★本ファイルはそれを**指数の形で証明した**。" ++
       "★★機構は j = E₄³/Δ による分解——‖E₄‖³ はカスプで有界(保型形式)、" ++
       "(Im τ)⁶·pD^ϵ → 0(Δ はカスプ形式)。" ++
       "★★★あとは §9-355 と同じ型(基本領域のコンパクト性 ＋ SL(2,ℤ) 不変性)") 8,
    .implicitStep
      ("★★これで Prop 3.4 に残るのは**代数の段だけ**である: " ++
       "ht∞ の同定(M̄_ell の無限遠因子の高さが h(j) であること)と、" ++
       "h(j) の有限素点の寄与が deg∞ であること(半安定のとき v(j) = −v(Δ_min))") 7 ]

end ABC3.Found.GenEll
