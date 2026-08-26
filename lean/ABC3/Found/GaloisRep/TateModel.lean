import ABC3.Found.GaloisRep.TateResidual

/-!
# Galois (G6) 第 309 ブロック —— **★★★★★★★★★★分裂乗法還元の曲線は Tate 曲線である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点

> `W` が分裂乗法還元をもち `Δ ≠ 0` なら、`q ∈ 𝔪`(`q ≠ 0`)と変数変換 `C` が在って
> **`C • (W の整モデル) = E_q`**(`exists_tate_model`)

★★★**G6 の第 3 段がこれで閉じた**。捻りの分類を経由せずに済んだ。

## ★★★★★★★★`a₆` は勝手に合う

第 308 で `a₄` は合わせられた。残る `a₆` は**合わせにいく必要がない**:

    Δ = a₄² − a₆ − 64a₄³ + 72a₄a₆ − 432a₆²             (標準形での判別式)

★`a₆` について 2 次だが、差を取ると

    Δ − Δ' = (a₆ − a₆')·(−1 + 72a₄ − 432(a₆ + a₆'))

と**因数分解する**。★★★★`a₄, a₆, a₆' ∈ 𝔪` なので右の括弧は **`−1 + 𝔪`、すなわち単元**。
したがって **`Δ` と `a₄` が決まれば `a₆` は一意**である(`normal_a6_unique`)。
★★★★★**標数を一切使わない**——これが「捻り」を迂回した道である。

## ★★★★★★`Δ` が合うのは `j` が合うから

`Δ/c₄³` は変数変換で不変(`Delta_c4_ratio`、分母を払った形)。
★第 304 の `q` は `Δ/c₄³` から取ったので、`a₄` を合わせた時点で
`c₄ = 1 − 48a₄` も合い、したがって `Δ` も合う。
★★つまり **`j` の一致が `Δ` の一致に化ける**のは、標準形では `c₄` が `a₄` だけで決まるから。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `normal_c4`・`normal_Delta` | ★★★標準形の `c₄`・`Δ` |
| `Delta_c4_ratio` | ★★★★`Δ/c₄³` の不変性 |
| `normal_a6_unique` | ★★★★★★★★**`a₆` の一意性** |
| `exists_tate_model` | ★★★★★★★★★★**`C • W = E_q`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain

section Normal

variable {R : Type} [CommRing R]

/-- ★★★標準形では `c₄ = 1 − 48a₄`。 -/
theorem normal_c4 (E : WeierstrassCurve R) (h1 : E.a₁ = 1) (h2 : E.a₂ = 0) (h3 : E.a₃ = 0) :
    E.c₄ = 1 - 48 * E.a₄ := by
  simp only [WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄, h1, h2, h3]
  ring

/-- ★★★標準形の判別式。 -/
theorem normal_Delta (E : WeierstrassCurve R) (h1 : E.a₁ = 1) (h2 : E.a₂ = 0) (h3 : E.a₃ = 0) :
    E.Δ = E.a₄ ^ 2 - E.a₆ - 64 * E.a₄ ^ 3 + 72 * E.a₄ * E.a₆ - 432 * E.a₆ ^ 2 := by
  simp only [WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    WeierstrassCurve.b₈, h1, h2, h3]
  ring

/-- ★★★★**`Δ/c₄³` は変数変換で不変**(分母を払った形)。 -/
theorem Delta_c4_ratio (E : WeierstrassCurve R) (C : VariableChange R) :
    (C • E).Δ * E.c₄ ^ 3 = E.Δ * (C • E).c₄ ^ 3 := by
  rw [WeierstrassCurve.variableChange_Δ, WeierstrassCurve.variableChange_c₄]
  ring

end Normal

section TateCoeff

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]

theorem tateCurveAt_a1 (q : R) (hq : q ∈ I) : (tateCurveAt q hq).a₁ = 1 := by
  show evalAdicHom q hq tateCurve.a₁ = 1
  rw [show tateCurve.a₁ = 1 from rfl, map_one]

theorem tateCurveAt_a2 (q : R) (hq : q ∈ I) : (tateCurveAt q hq).a₂ = 0 := by
  show evalAdicHom q hq tateCurve.a₂ = 0
  rw [show tateCurve.a₂ = 0 from rfl, map_zero]

theorem tateCurveAt_a3 (q : R) (hq : q ∈ I) : (tateCurveAt q hq).a₃ = 0 := by
  show evalAdicHom q hq tateCurve.a₃ = 0
  rw [show tateCurve.a₃ = 0 from rfl, map_zero]

end TateCoeff

section Unique

variable {R : Type} [CommRing R] [IsLocalRing R]

/-- ★★★★★★★★**`a₆` は `a₄` と `Δ` から一意に決まる**(`a₆ ∈ 𝔪` の範囲で)。

★差が `(a₆ − a₆')·(−1 + 72a₄ − 432(a₆+a₆'))` と因数分解し、括弧が単元だから。 -/
theorem normal_a6_unique (E E' : WeierstrassCurve R)
    (h1 : E.a₁ = 1) (h2 : E.a₂ = 0) (h3 : E.a₃ = 0)
    (h1' : E'.a₁ = 1) (h2' : E'.a₂ = 0) (h3' : E'.a₃ = 0)
    (h4 : E.a₄ ∈ IsLocalRing.maximalIdeal R)
    (h6 : E.a₆ ∈ IsLocalRing.maximalIdeal R) (h6' : E'.a₆ ∈ IsLocalRing.maximalIdeal R)
    (ha4 : E.a₄ = E'.a₄) (hΔ : E.Δ = E'.Δ) : E.a₆ = E'.a₆ := by
  have hE := normal_Delta E h1 h2 h3
  have hE' := normal_Delta E' h1' h2' h3'
  rw [hE, hE', ← ha4] at hΔ
  have hfac : (E.a₆ - E'.a₆) * (-1 + 72 * E.a₄ - 432 * (E.a₆ + E'.a₆)) = 0 := by
    linear_combination hΔ
  have hunit : IsUnit (-1 + 72 * E.a₄ - 432 * (E.a₆ + E'.a₆)) := by
    have hz : -1 + 72 * E.a₄ - 432 * (E.a₆ + E'.a₆)
        = -(1 + (-72 * E.a₄ + 432 * (E.a₆ + E'.a₆))) := by ring
    rw [hz]
    refine IsUnit.neg (isUnit_one_add_mem ?_)
    exact Ideal.add_mem _ (Ideal.mul_mem_left _ _ h4)
      (Ideal.mul_mem_left _ _ (Ideal.add_mem _ h6 h6'))
  exact sub_eq_zero.1 ((hunit.mul_left_eq_zero).1 hfac)

end Unique

/-! ## ★★★★★★★★★★分裂乗法還元の曲線は Tate 曲線 -/

section Model

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**分裂乗法還元をもつ曲線は、変数変換で Tate 曲線になる**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_tate_model (W : WeierstrassCurve K) [hmul : W.HasMultiplicativeReduction R]
    (hsplit : W.HasSplitMultiplicativeReduction R) (hΔ0 : W.Δ ≠ 0) :
    ∃ (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (C : VariableChange R),
      q ≠ 0 ∧ C • integralModel R W = tateCurveAt q hq := by
  obtain ⟨C₀, n1, n2, n3, n4, n6⟩ := exists_tate_normal_form_of_split W hsplit hΔ0
  have hΔE : (integralModel R W).Δ ∈ IsLocalRing.maximalIdeal R := integralModel_Delta_mem W hΔ0
  have hΔEne : (integralModel R W).Δ ≠ 0 := by
    intro hc
    apply hΔ0
    rw [← WeierstrassCurve.integralModel_Δ_eq R W, hc, map_zero]
  have hΔ₀ : (C₀ • integralModel R W).Δ ∈ IsLocalRing.maximalIdeal R := by
    rw [WeierstrassCurve.variableChange_Δ]
    exact Ideal.mul_mem_left _ _ hΔE
  have hΔ₀ne : (C₀ • integralModel R W).Δ ≠ 0 := by
    rw [WeierstrassCurve.variableChange_Δ]
    exact mul_ne_zero (pow_ne_zero _ (Units.ne_zero _)) hΔEne
  have hc4₀ : IsUnit (C₀ • integralModel R W).c₄ := by
    rw [normal_c4 _ n1 n2 n3]
    have hz : (1 : R) - 48 * (C₀ • integralModel R W).a₄
        = 1 + (-48 * (C₀ • integralModel R W).a₄) := by ring
    rw [hz]
    exact isUnit_one_add_mem (Ideal.mul_mem_left _ _ n4)
  set cinv : R := ((hc4₀.unit⁻¹ : Rˣ) : R) with hcinvdef
  have hcinv : cinv * (C₀ • integralModel R W).c₄ = 1 := by
    have h2 : ((hc4₀.unit⁻¹ : Rˣ) : R) * ((hc4₀.unit : Rˣ) : R) = 1 := by
      rw [← Units.val_mul, inv_mul_cancel]
      rfl
    rw [hc4₀.unit_spec] at h2
    exact h2
  have hcinvne : cinv ≠ 0 := by
    rw [hcinvdef]
    exact Units.ne_zero _
  set t₀ : R := (C₀ • integralModel R W).Δ * cinv ^ 3 with ht₀def
  have ht₀mem : t₀ ∈ IsLocalRing.maximalIdeal R := Ideal.mul_mem_right _ _ hΔ₀
  have ht₀ne : t₀ ≠ 0 := mul_ne_zero hΔ₀ne (pow_ne_zero _ hcinvne)
  have hΔ₀eq : (C₀ • integralModel R W).Δ = t₀ * (C₀ • integralModel R W).c₄ ^ 3 := by
    rw [ht₀def]
    linear_combination (-((C₀ • integralModel R W).Δ)
      * ((cinv * (C₀ • integralModel R W).c₄) ^ 2 + (cinv * (C₀ • integralModel R W).c₄) + 1))
      * hcinv
  obtain ⟨q, hq, hqe⟩ := exists_tateParam (I := IsLocalRing.maximalIdeal R) ht₀mem
  have hΔq : (tateCurveAt q hq).Δ = t₀ * (tateCurveAt q hq).c₄ ^ 3 := by
    rw [tateCurveAt_Delta, tateCurveAt_c4, ← map_pow]
    exact hqe
  have hqne : q ≠ 0 := by
    intro hq0
    obtain ⟨v, hv, hΔeq⟩ := tateCurveAt_Delta_eq_mul_unit q hq
    have hzero : (tateCurveAt q hq).Δ = 0 := by
      rw [hΔeq, hq0, zero_mul]
    rw [hzero] at hΔq
    have hc4u : IsUnit ((tateCurveAt q hq).c₄ ^ 3) := (tateCurveAt_c4_isUnit q hq).pow 3
    rcases mul_eq_zero.1 hΔq.symm with h | h
    · exact ht₀ne h
    · rw [h] at hc4u
      exact not_isUnit_zero hc4u
  obtain ⟨C₁, m1, m2, m3, m4, m6⟩ := exists_residual_match_a4 (C₀ • integralModel R W) n1 n2 n3
    n4 n6 ((tateCurveAt q hq).a₄) (tateCurveAt_a4_mem q hq)
  have hc4eq : (C₁ • (C₀ • integralModel R W)).c₄ = (tateCurveAt q hq).c₄ := by
    rw [normal_c4 _ m1 m2 m3,
      normal_c4 _ (tateCurveAt_a1 q hq) (tateCurveAt_a2 q hq) (tateCurveAt_a3 q hq), m4]
  have hΔ₁ : (C₁ • (C₀ • integralModel R W)).Δ = t₀ * (C₁ • (C₀ • integralModel R W)).c₄ ^ 3 := by
    have hratio := Delta_c4_ratio (C₀ • integralModel R W) C₁
    rw [hΔ₀eq] at hratio
    refine (hc4₀.pow 3).mul_left_cancel ?_
    linear_combination hratio
  have hΔeq : (C₁ • (C₀ • integralModel R W)).Δ = (tateCurveAt q hq).Δ := by
    rw [hΔ₁, hΔq, hc4eq]
  have ha6eq : (C₁ • (C₀ • integralModel R W)).a₆ = (tateCurveAt q hq).a₆ :=
    normal_a6_unique _ _ m1 m2 m3 (tateCurveAt_a1 q hq) (tateCurveAt_a2 q hq)
      (tateCurveAt_a3 q hq) (by rw [m4]; exact tateCurveAt_a4_mem q hq) m6
      (tateCurveAt_a6_mem q hq) m4 hΔeq
  refine ⟨q, hq, C₁ * C₀, hqne, ?_⟩
  rw [mul_smul]
  exact WeierstrassCurve.ext (m1.trans (tateCurveAt_a1 q hq).symm)
    (m2.trans (tateCurveAt_a2 q hq).symm) (m3.trans (tateCurveAt_a3 q hq).symm) m4 ha6eq

end Model

/-! ## ★出典の紐付け(`.src`) -/

def normal_a6_unique.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(標準形の a₆ の一意性)",
    sectionId := "genell-def-3-3" }

def exists_tate_model.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(分裂乗法還元の曲線は Tate 曲線)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
