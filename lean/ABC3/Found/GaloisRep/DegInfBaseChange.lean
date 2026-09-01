/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.PrimeUnder
import ABC3.Found.GaloisRep.FaltingsWitness
import ABC3.Found.GaloisRep.SemistableFin
import Mathlib.NumberTheory.RamificationInertia.Valuation
import ABC3.Meta.Claim
import ABC3.Found.GenEll.LocalHeightRamified

/-!
# `deg∞` と `ht^Falt` の基底変換不変性（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★これは何か

`ResearchPaper/ellmoduli-witness-status.json` の `ellModuliWitnessIsTheBottleneck20260830`
で絞り込んだ通り、`EllModuliData` の `faltingsHeight : EllClass → ℝ`(=`j` の函数)を
作るには「同じ `j` なら同じ `ht^Falt`」が要る。その第一段が**基底変換不変性**である。

## ★★★機構

| 段 | 補題 | 状態 |
|---|---|---|
| (1) 付値は `e` 倍 | `valAdd_algebraMap`（本ファイル） | ★済(mathlib の `valuation_liesOver`) |
| (2) `finsum` は `[L′:L]` 倍 | `finsum_scaling`（`§9-1163`、第 736） | ★済 |
| (3) `deg∞` は不変 | `degInfOf_baseChange_of_minDeltaExp`（本ファイル） | ★済(局所の仮説つき) |
| (4) `ht^Falt` は不変 | `htFaltOf_baseChange_of_minDeltaExp`（本ファイル） | ★済(同上) |
| (5) 局所の仮説 | `minDeltaExp P (E×L′) = e·minDeltaExp p E` | ☆半安定の Tate 母数から |

★(3) の分母の相殺が要点である:

    deg∞(E×L′) = ([L′:L]·N_L) / [L′:ℚ] = ([L′:L]·N_L) / ([L:ℚ]·[L′:L]) = N_L / [L:ℚ] = deg∞(E)
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve

/-! ## ★配管 -/

section Plumbing

variable {R A : Type*} [CommRing R] [CommRing A] [Algebra R A]

/-- ★★底変換した曲線も楕円曲線である（mathlib は `map` の形でしか持っていない）。 -/
instance isElliptic_baseChange (W : WeierstrassCurve R) [W.IsElliptic] :
    (W.baseChange A).IsElliptic :=
  inferInstanceAs ((W.map (algebraMap R A)).IsElliptic)

end Plumbing

/-- ★`e ≥ 0` なら `max 0 (e·x) = e·max 0 x`。 -/
theorem max_zero_mul (e : ℤ) (he : 0 ≤ e) (x : ℤ) : max 0 (e * x) = e * max 0 x := by
  rcases le_total 0 x with h | h
  · rw [max_eq_right h, max_eq_right (mul_nonneg he h)]
  · rw [max_eq_left h, max_eq_left (by nlinarith : e * x ≤ 0), mul_zero]

/-! ## ★★★★★★★★★★半安定なら `v_p(Δ_min) = log⁺|j|_p` -/

section OneField

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★★★★★★★**半安定なら `v_p(Δ_min) = max(0, −v_p(j))`**。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★`≤` は半安定から（`minDeltaExp_le_maxJ`、`§9-1003`）、
`≥` は**無条件**（`maxJ_le_minDeltaExp`、極小モデルの `c₄` が整だから）。
★★★これで `minDeltaExp` が **`j` だけで書けた**——基底変換の議論はここから出る。 -/
theorem minDeltaExp_eq_maxJ_of_semistable (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] (hss : SemistableAt p W) :
    minDeltaExp p W = max 0 (-jExp p W) :=
  le_antisymm (minDeltaExp_le_maxJ p W hss) (maxJ_le_minDeltaExp p W)

/-- ★★`valAdd` の順序は乗法的付値の逆順序。 -/
theorem valAdd_le_iff (p : HeightOneSpectrum (𝓞 L)) (x y : Lˣ) :
    valAdd p x ≤ valAdd p y ↔ (p.valuation L) (y : L) ≤ (p.valuation L) (x : L) := by
  rw [valAdd, valAdd, neg_le_neg_iff, Multiplicative.toAdd_le, ← WithZero.coe_le_coe,
    WithZero.coe_unzero, WithZero.coe_unzero]

/-- ★★★★★★**極小判別式の指数は、どの整モデルの `v_p(Δ)` よりも小さい**。

★`minimal_vAdd_Delta_le`（`DegInfLocal.lean:66`）を `valAdd` の言葉に直したもの。 -/
theorem minDeltaExp_le_of_isIntegral (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    (hΔ : W.Δ ≠ 0) (C : WeierstrassCurve.VariableChange L)
    (hint : WeierstrassCurve.IsIntegral (primeSubring p) (C • W)) :
    minDeltaExp p W
      ≤ valAdd p (Units.mk0 ((C • W).Δ) (variableChange_Delta_ne_zero W hΔ C)) := by
  obtain ⟨C₀, hC₀⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) W
  haveI := hC₀
  have hΔ₀ : (C₀ • W).Δ ≠ 0 := variableChange_Delta_ne_zero W hΔ C₀
  have hsm : (C * C₀⁻¹) • (C₀ • W) = C • W := by rw [mul_smul, inv_smul_smul]
  have hkey := minimal_vAdd_Delta_le (R := primeSubring p) (C₀ • W) hΔ₀ (C * C₀⁻¹)
    (by rw [hsm]; exact hint)
  have hv : (IsDiscreteValuationRing.maximalIdeal (primeSubring p)).valuation L ((C • W).Δ)
      ≤ (IsDiscreteValuationRing.maximalIdeal (primeSubring p)).valuation L ((C₀ • W).Δ) := by
    have h1 := (valuation_le_iff_vAdd_le (R := primeSubring p) _ _).2 hkey
    have h2 : ((C * C₀⁻¹) • (C₀ • W)).Δ = (C • W).Δ := by rw [hsm]
    rw [show (((Units.mk0 (((C * C₀⁻¹) • (C₀ • W)).Δ)
        (variableChange_Delta_ne_zero (C₀ • W) hΔ₀ (C * C₀⁻¹))) : L)) = (C • W).Δ from h2] at h1
    exact h1
  rw [minDeltaExp_eq p W hΔ C₀ hC₀, valAdd_le_iff]
  exact (valuation_isEquiv p _ _).1 hv

/-- ★★★整モデルの係数は付値環に入っている。 -/
theorem mem_primeSubring_of_isIntegral (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] :
    W.a₁ ∈ primeSubring p ∧ W.a₂ ∈ primeSubring p ∧ W.a₃ ∈ primeSubring p ∧
      W.a₄ ∈ primeSubring p ∧ W.a₆ ∈ primeSubring p := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · rw [← WeierstrassCurve.integralModel_a₁_eq (primeSubring p) W]
    exact ((integralModel (primeSubring p) W).a₁).2
  · rw [← WeierstrassCurve.integralModel_a₂_eq (primeSubring p) W]
    exact ((integralModel (primeSubring p) W).a₂).2
  · rw [← WeierstrassCurve.integralModel_a₃_eq (primeSubring p) W]
    exact ((integralModel (primeSubring p) W).a₃).2
  · rw [← WeierstrassCurve.integralModel_a₄_eq (primeSubring p) W]
    exact ((integralModel (primeSubring p) W).a₄).2
  · rw [← WeierstrassCurve.integralModel_a₆_eq (primeSubring p) W]
    exact ((integralModel (primeSubring p) W).a₆).2

def minDeltaExp_eq_maxJ_of_semistable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(半安定なら v_p(Δ_min) = log⁺|j|_p——等号)",
    sectionId := "genell-prop-3-4" }

end OneField

section NumberField

variable (L L' : Type) [Field L] [NumberField L] [Field L'] [NumberField L']
  [Algebra L L'] [Algebra (𝓞 L) L'] [IsScalarTower (𝓞 L) L L']
  [IsScalarTower (𝓞 L) (𝓞 L') L'] [Module.Finite (𝓞 L) (𝓞 L')]

/-! ## ★★★★★★★★付値は `e` 倍になる -/

/-- ★★★★★★★★★★**`v_P(x) = e(P|p)·v_p(x)`**（`x ∈ Lˣ`）。

★mathlib の `IsDedekindDomain.HeightOneSpectrum.valuation_liesOver`
（`Mathlib/NumberTheory/RamificationInertia/Valuation.lean:56`）
`v(x)^e = w(algebraMap x)` を、本プロジェクトの加法付値 `valAdd` の言葉に直したもの。 -/
theorem valAdd_algebraMap (p : HeightOneSpectrum (𝓞 L)) (P : HeightOneSpectrum (𝓞 L'))
    [P.asIdeal.LiesOver p.asIdeal] (x : Lˣ) :
    valAdd P (Units.map (algebraMap L L' : L →* L') x)
      = (p.asIdeal.ramificationIdx P.asIdeal : ℤ) * valAdd p x := by
  have hu : WithZero.unzero (valuationP_ne_zero P (Units.map (algebraMap L L' : L →* L') x))
      = WithZero.unzero (valuationP_ne_zero p x) ^ p.asIdeal.ramificationIdx P.asIdeal := by
    apply WithZero.coe_injective
    rw [WithZero.coe_unzero, WithZero.coe_pow, WithZero.coe_unzero,
      HeightOneSpectrum.valuation_liesOver L' p P (x : L)]
    rfl
  rw [valAdd, valAdd, hu, toAdd_pow, nsmul_eq_mul]
  push_cast
  ring

/-! ## ★★★★★★★★`deg∞` の基底変換不変性 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★
**`deg∞` は基底変換で変わらない**（局所の関係を仮説として受けた形）。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★仮説 `hloc` は「`v_P(Δ_min(E×L′)) = e(P|p)·v_p(Δ_min(E))`」——
**半安定な曲線では極小モデルが基底変換で極小のままである**ことの帰結である
（Tate 母数 `q` の付値がそのまま `e` 倍になる）。 -/
theorem degInfOf_baseChange_of_minDeltaExp (E : WeierstrassCurve L) (hΔ : E.Δ ≠ 0)
    (hloc : ∀ P : HeightOneSpectrum (𝓞 L'),
      (minDeltaExp P (E.baseChange L') : ℝ)
        = ((HeightOneSpectrumUnder (A := 𝓞 L) P).asIdeal.ramificationIdx P.asIdeal : ℝ)
          * (minDeltaExp (HeightOneSpectrumUnder (A := 𝓞 L) P) E : ℝ))
    [IsScalarTower ℚ L L'] :
    degInfOf L' (E.baseChange L') = degInfOf L E := by
  classical
  haveI : FiniteDimensional L L' := FiniteDimensional.right ℚ L L'
  have hfin := minDeltaExp_finite E hΔ
  have hSmem : ∀ p : HeightOneSpectrum (𝓞 L),
      (minDeltaExp p E : ℝ) ≠ 0 → p ∈ hfin.toFinset := by
    intro p hp
    rw [Set.Finite.mem_toFinset]
    intro h
    exact hp (by rw [h]; norm_num)
  have hnum := finsum_scaling L L' (fun p => (minDeltaExp p E : ℝ))
    (fun P => (minDeltaExp P (E.baseChange L') : ℝ)) hloc hfin.toFinset hSmem
  have htower : (Module.finrank ℚ L' : ℝ)
      = (Module.finrank ℚ L : ℝ) * (Module.finrank L L' : ℝ) := by
    rw [← Nat.cast_mul, Module.finrank_mul_finrank ℚ L L']
  have h1 : (Module.finrank ℚ L : ℝ) ≠ 0 := by
    have : 0 < Module.finrank ℚ L := Module.finrank_pos
    positivity
  have h2 : (Module.finrank L L' : ℝ) ≠ 0 := by
    have : 0 < Module.finrank L L' := Module.finrank_pos
    positivity
  rw [degInfOf, degInfOf, hnum, htower]
  field_simp

/-! ## ★★★★★★★★★`ht^Falt` の基底変換不変性 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`ht^Falt` は基底変換で変わらない**（局所の関係を仮説として受けた形）。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★アルキメデス側は `archSum_baseChange`（`§9-1157`、第 725）で既に済んでおり、
有限素点側が本ファイルの `degInfOf_baseChange_of_minDeltaExp` である。
★★これで `faltingsHeight : ℂ → ℝ` を `j` の函数として定義する道の
**片方（体を大きくする向き）が通った**。 -/
theorem htFaltOf_baseChange_of_minDeltaExp (E : WeierstrassCurve L) (hΔ : E.Δ ≠ 0)
    (hloc : ∀ P : HeightOneSpectrum (𝓞 L'),
      (minDeltaExp P (E.baseChange L') : ℝ)
        = ((HeightOneSpectrumUnder (A := 𝓞 L) P).asIdeal.ramificationIdx P.asIdeal : ℝ)
          * (minDeltaExp (HeightOneSpectrumUnder (A := 𝓞 L) P) E : ℝ))
    [IsScalarTower ℚ L L'] :
    htFaltOf L' (E.baseChange L') = htFaltOf L E :=
  htFaltOf_baseChange_of_degInf L L' E
    (degInfOf_baseChange_of_minDeltaExp L L' E hΔ hloc)

/-! ## ★★★★★★★★★★★★整性と極小性の基底変換 -/

/-- ★★★★★★**整なモデルは基底変換しても整**。

★`v_P(f x) = v_p(x)^e` なので `v_p(x) ≤ 1` から `v_P(f x) ≤ 1`。 -/
theorem isIntegral_baseChange_primeSubring (p : HeightOneSpectrum (𝓞 L)) (P : HeightOneSpectrum (𝓞 L'))
    [P.asIdeal.LiesOver p.asIdeal] (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] :
    WeierstrassCurve.IsIntegral (primeSubring P) (W.baseChange L') := by
  have key : ∀ x : L, x ∈ primeSubring p → algebraMap L L' x ∈ primeSubring P := by
    intro x hx
    rw [mem_primeSubring_iff] at hx ⊢
    rw [← HeightOneSpectrum.valuation_liesOver L' p P x]
    exact pow_le_one₀ zero_le' hx
  obtain ⟨h1, h2, h3, h4, h6⟩ := mem_primeSubring_of_isIntegral p W
  exact isIntegral_of_mem (primeSubring P) (W.baseChange L')
    (key _ h1) (key _ h2) (key _ h3) (key _ h4) (key _ h6)

/-- ★★★★★★★★★★★★★★**`v_P(Δ_min(E×L′)) ≤ e(P|p)·v_p(Δ_min(E))`**——★**無条件**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`p` で極小なモデルを基底変換すると `P` で**整**なモデルになる。極小判別式は
どの整モデルの `v(Δ)` よりも小さいので、この向きの不等式は半安定性なしで出る。 -/
theorem minDeltaExp_baseChange_le (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    (E : WeierstrassCurve L) (hΔ : E.Δ ≠ 0) :
    minDeltaExp P (E.baseChange L')
      ≤ (p.asIdeal.ramificationIdx P.asIdeal : ℤ) * minDeltaExp p E := by
  obtain ⟨C₀, hC₀⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) E
  haveI := hC₀
  haveI hint : WeierstrassCurve.IsIntegral (primeSubring p) (C₀ • E) := inferInstance
  have hΔ₀ : (C₀ • E).Δ ≠ 0 := variableChange_Delta_ne_zero E hΔ C₀
  have hΔ' : (E.baseChange L').Δ ≠ 0 := by
    rw [show (E.baseChange L').Δ = algebraMap L L' E.Δ from WeierstrassCurve.map_Δ _ _]
    exact fun h => hΔ ((map_eq_zero_iff _ (RingHom.injective (algebraMap L L'))).1 h)
  have hsm : (C₀.map (algebraMap L L')) • (E.baseChange L') = (C₀ • E).baseChange L' :=
    WeierstrassCurve.map_variableChange E C₀ (algebraMap L L')
  haveI hintP : WeierstrassCurve.IsIntegral (primeSubring P)
      ((C₀.map (algebraMap L L')) • (E.baseChange L')) := by
    rw [hsm]
    exact isIntegral_baseChange_primeSubring L L' p P (C₀ • E)
  have hle := minDeltaExp_le_of_isIntegral P (E.baseChange L') hΔ' _ hintP
  have heq : valAdd P (Units.mk0 (((C₀.map (algebraMap L L')) • (E.baseChange L')).Δ)
        (variableChange_Delta_ne_zero (E.baseChange L') hΔ' (C₀.map (algebraMap L L'))))
      = valAdd P (Units.map (algebraMap L L' : L →* L') (Units.mk0 ((C₀ • E).Δ) hΔ₀)) := by
    refine valAdd_eq_of_valuation_eq P _ _ ?_
    congr 1
    show ((C₀.map (algebraMap L L')) • (E.baseChange L')).Δ = algebraMap L L' ((C₀ • E).Δ)
    rw [hsm]
    exact WeierstrassCurve.map_Δ _ _
  rw [heq, valAdd_algebraMap L L' p P, ← minDeltaExp_eq p E hΔ C₀ hC₀] at hle
  exact hle

/-! ## ★★★★★★★★★★★★半安定なら局所の仮説は自動である -/

/-- ★★★★★★★★**`v_P(j) = e(P|p)·v_p(j)`**。 -/
theorem jExp_baseChange (p : HeightOneSpectrum (𝓞 L)) (P : HeightOneSpectrum (𝓞 L'))
    [P.asIdeal.LiesOver p.asIdeal] (E : WeierstrassCurve L) [E.IsElliptic] :
    jExp P (E.baseChange L') = (p.asIdeal.ramificationIdx P.asIdeal : ℤ) * jExp p E := by
  have hj : (E.baseChange L').j = algebraMap L L' E.j :=
    WeierstrassCurve.map_j (W := E) (f := algebraMap L L')
  by_cases h0 : E.j = 0
  · have hz : (E.baseChange L').j = 0 := by rw [hj, h0, map_zero]
    rw [jExp, dif_pos hz, jExp, dif_pos h0, mul_zero]
  · have hne : (E.baseChange L').j ≠ 0 := by
      rw [hj]
      exact fun h => h0 ((map_eq_zero_iff _ (RingHom.injective (algebraMap L L'))).1 h)
    rw [jExp, dif_neg hne, jExp, dif_neg h0]
    have hunit : Units.mk0 ((E.baseChange L').j) hne
        = Units.map (algebraMap L L' : L →* L') (Units.mk0 E.j h0) := Units.ext hj
    rw [hunit, valAdd_algebraMap L L' p P]

/-- ★★★★★★★★★★★★★★★★★★★★★★
**`v_P(Δ_min(E×L′)) = e(P|p)·v_p(Δ_min(E))`**——★`E` が `p` で半安定なだけでよい。

★★「上でも半安定」を仮定しなくてよいのが要点である:

* `≤` は `minDeltaExp_baseChange_le`（整モデルの比較、★無条件）
* `≥` は `maxJ_le_minDeltaExp`（★無条件）と、下での半安定性 -/
theorem minDeltaExp_baseChange_of_semistableAt (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    (E : WeierstrassCurve L) [E.IsElliptic] (hp : SemistableAt p E) :
    minDeltaExp P (E.baseChange L')
      = (p.asIdeal.ramificationIdx P.asIdeal : ℤ) * minDeltaExp p E := by
  refine le_antisymm (minDeltaExp_baseChange_le L L' p P E E.isUnit_Δ.ne_zero) ?_
  rw [minDeltaExp_eq_maxJ_of_semistable p E hp,
    ← max_zero_mul _ (by positivity : (0:ℤ) ≤ (p.asIdeal.ramificationIdx P.asIdeal : ℤ)),
    mul_neg, ← jExp_baseChange L L' p P E]
  exact maxJ_le_minDeltaExp P (E.baseChange L')

/-- ★★★★★★★★★★★★★★★★★★
**下で半安定なら、上でも `v_P(Δ_min) = max(0, −v_P(j))` が成り立つ**。

★これが `HtFaltJ.lean` の `j`-合同の議論へ渡す形である
——`SemistableAt` そのものではなく、**`maxJ` の等式**だけを持ち上げる。 -/
theorem minDeltaExp_eq_maxJ_baseChange (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    (E : WeierstrassCurve L) [E.IsElliptic] (hp : SemistableAt p E) :
    minDeltaExp P (E.baseChange L') = max 0 (-jExp P (E.baseChange L')) := by
  rw [minDeltaExp_baseChange_of_semistableAt L L' p P E hp,
    minDeltaExp_eq_maxJ_of_semistable p E hp, jExp_baseChange L L' p P E, ← mul_neg,
    max_zero_mul _ (by positivity : (0:ℤ) ≤ (p.asIdeal.ramificationIdx P.asIdeal : ℤ))]

/-- ★★★★★★★★★★★★★★★★★★★★
**半安定なら `v_P(Δ_min(E×L′)) = e(P|p)·v_p(Δ_min(E))`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★これが `degInfOf_baseChange_of_minDeltaExp`（`§9-1164`、第 737）が受けていた
局所の仮説そのものである。**半安定性だけから出る**——Tate 曲線の議論は要らない。

★機構は 2 行:

* 半安定なら `v(Δ_min) = max(0, −v(j))`（`minDeltaExp_eq_maxJ_of_semistable`）
* `v_P(j) = e·v_p(j)`（`jExp_baseChange`）、そして `e ≥ 0` だから `max` は `e` 倍で通る -/
theorem minDeltaExp_baseChange_of_semistable (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    (E : WeierstrassCurve L) [E.IsElliptic]
    (hp : SemistableAt p E) :
    minDeltaExp P (E.baseChange L')
      = (p.asIdeal.ramificationIdx P.asIdeal : ℤ) * minDeltaExp p E :=
  minDeltaExp_baseChange_of_semistableAt L L' p P E hp

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★到達点 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**半安定な楕円曲線の `deg∞` は基底変換で変わらない**——★**無条件**。 -/
theorem degInfOf_baseChange_of_semistable (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    [IsScalarTower ℚ L L'] :
    degInfOf L' (E.baseChange L') = degInfOf L E := by
  refine degInfOf_baseChange_of_minDeltaExp L L' E E.isUnit_Δ.ne_zero (fun P => ?_)
  exact_mod_cast congrArg (fun n : ℤ => (n : ℝ))
    (minDeltaExp_baseChange_of_semistableAt L L' _ P E (hss _))

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**半安定な楕円曲線の `ht^Falt` は基底変換で変わらない**——★**無条件**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★これで `EllModuliData` の `faltingsHeight : EllClass → ℝ` を作る道の
**体を大きくする向きが完全に通った**（残るのは捻り＝同じ `j` の別の曲線の扱いである）。 -/
theorem htFaltOf_baseChange_of_semistable (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    [IsScalarTower ℚ L L'] :
    htFaltOf L' (E.baseChange L') = htFaltOf L E :=
  htFaltOf_baseChange_of_degInf L L' E
    (degInfOf_baseChange_of_semistable L L' E hss)

def htFaltOf_baseChange_of_semistable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(半安定な曲線の ht^Falt は基底変換で不変。★無条件)",
    sectionId := "genell-prop-3-4" }

/-! ## ★★★★★★★★★★★★★★★★上での `Δ_min` の関係を下へ降ろす -/

/-- ★★★★★★★★★★★★★★★★★★★★
**上で `v_P(Δ_min(E')) = l·v_P(Δ_min(E))` なら下でも同じ**——★**無条件**（第 1183）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆両辺とも分岐指数 `e` 倍になるので、`e ≠ 0` で割れば下の関係が出る。
★★★これが `Skeleton/GenEll/LCyclicReading.lean` の**節点 2d-2** である
——安定直線の側では `ζ` が `L_v` に無いので局所の議論を `L_v(ζ_l)` へ上げるが、
そこで得た `Δ_min` の関係は**そのまま降りる**。

☆`Thm38Kummer.lean` の「分岐指数は `l` と素」は**要らなかった**
——必要なのは `e ≠ 0` だけである（第 1182 の見積もりの再訂正）。 -/
theorem minDeltaExp_descend_of_baseChange (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hp : SemistableAt p E) (hp' : SemistableAt p E')
    (he : p.asIdeal.ramificationIdx P.asIdeal ≠ 0) {l : ℕ}
    (h : minDeltaExp P (E'.baseChange L') = l * minDeltaExp P (E.baseChange L')) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  rw [minDeltaExp_baseChange_of_semistableAt L L' p P E' hp',
    minDeltaExp_baseChange_of_semistableAt L L' p P E hp] at h
  have he' : ((p.asIdeal.ramificationIdx P.asIdeal : ℤ)) ≠ 0 := by
    exact_mod_cast he
  refine mul_left_cancel₀ he' ?_
  rw [h]
  ring

/-- ★★★★★★★★★★★★★★**上での `Δ_min` の関係は下へ降りる（仮説なし）**（第 1185）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1183 は `e ≠ 0` を仮説に取っていたが、`P` が `p` の上にあれば
mathlib の `ramificationIdx_ne_zero_of_liesOver` で出る。
★★これで `Skeleton/GenEll/LCyclicReading.lean` の**節点 2d-2 は仮説なし**になった。 -/
theorem minDeltaExp_descend (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hp : SemistableAt p E) (hp' : SemistableAt p E') {l : ℕ}
    (h : minDeltaExp P (E'.baseChange L') = l * minDeltaExp P (E.baseChange L')) :
    minDeltaExp p E' = l * minDeltaExp p E :=
  minDeltaExp_descend_of_baseChange L L' p P E E' hp hp'
    (Ideal.IsDedekindDomain.ramificationIdx_ne_zero_of_liesOver P.asIdeal p.ne_bot) h

/-- ★★★★★★★★**良い還元は底変換で保たれる**（第 1187）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`v_p(Δ_min) = 0` なら上でも `0` である——`minDeltaExp_baseChange_le` で
`≤ e·0 = 0`、`minDeltaExp_nonneg` で `≥ 0`。

★★これが `SemistableAt` の底変換の**第 1 の場合**である
（`Skeleton/GenEll/LCyclicReading.lean` の節点 2d-1、第 1186 の実測）。
☆第 2 の場合（極小モデルと `v(c₄) = 0`）は分岐した `valAdd` のスケーリングが要る。 -/
theorem minDeltaExp_baseChange_eq_zero (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    (E : WeierstrassCurve L) [E.IsElliptic] (h : minDeltaExp p E = 0) :
    minDeltaExp P (E.baseChange L') = 0 := by
  have hle := minDeltaExp_baseChange_le L L' p P E E.isUnit_Δ.ne_zero
  rw [h, mul_zero] at hle
  have hnn := minDeltaExp_nonneg P (E.baseChange L')
  omega

/-- ☆`valAdd` は `LocalHeightRamified.lean` の `ordAt` と同じものである。 -/
theorem valAdd_eq_ordAt (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) :
    valAdd p x = ABC3.Found.GenEll.ordAt p x := rfl

/-- ★★★★★★★★★★★★★★**分岐した拡大での `valAdd` のスケーリング**（第 1188）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`valAdd P (algebraMap x) = e(P|p) · valAdd p x`。
★`vAdd_algebraMap_eq_valAdd`（第 320）は**不分岐を仮定していた**が、
`LocalHeightRamified.lean` の `ordAt_liesOver` は分岐を許す。
☆`valAdd` と `ordAt` は同じ定義なので、そのまま移せる。

★★★これが `Skeleton/GenEll/LCyclicReading.lean` の節点 2d-1 で
第 1186 が名指しした「分岐した `valAdd` のスケーリング」である。 -/
theorem valAdd_algebraMap_liesOver (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal] (x : Lˣ) :
    valAdd P (Units.map (algebraMap L L').toMonoidHom x)
      = (p.asIdeal.ramificationIdx P.asIdeal : ℤ) * valAdd p x :=
  ABC3.Found.GenEll.ordAt_liesOver L' p P x

/-- ★★★★★★★★**付値環は底変換で付値環に入る**（第 1191）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`x ∈ primeSubring p` なら `algebraMap x ∈ primeSubring P`——
`valAdd P (algebraMap x) = e·valAdd p x ≥ 0` だからである（第 1188）。
★これが「整モデルを上げると整モデルになる」の中身である。 -/
theorem algebraMap_mem_primeSubring (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    {x : L} (hx : x ∈ primeSubring p) :
    algebraMap L L' x ∈ primeSubring P := by
  rcases eq_or_ne x 0 with rfl | hx0
  · simpa using (primeSubring P).zero_mem
  · have hx0' : algebraMap L L' x ≠ 0 :=
      (map_ne_zero_iff _ (algebraMap L L').injective).2 hx0
    have hge : 0 ≤ valAdd p (Units.mk0 x hx0) :=
      (valAdd_nonneg_iff p _).2 ((mem_primeSubring_iff p x).1 hx)
    have heq : valAdd P (Units.map (algebraMap L L').toMonoidHom (Units.mk0 x hx0))
        = (p.asIdeal.ramificationIdx P.asIdeal : ℤ) * valAdd p (Units.mk0 x hx0) :=
      valAdd_algebraMap_liesOver L L' p P (Units.mk0 x hx0)
    have hval : valAdd P (Units.mk0 (algebraMap L L' x) hx0')
        = (p.asIdeal.ramificationIdx P.asIdeal : ℤ) * valAdd p (Units.mk0 x hx0) := by
      rw [← heq]
      exact valAdd_eq_of_valuation_eq P _ _ rfl
    have hnn : (0 : ℤ) ≤ (p.asIdeal.ramificationIdx P.asIdeal : ℤ) := by positivity
    have : 0 ≤ valAdd P (Units.mk0 (algebraMap L L' x) hx0') := by
      rw [hval]
      exact mul_nonneg hnn hge
    exact (mem_primeSubring_iff P _).2 ((valAdd_nonneg_iff P _).1 this)

end NumberField

/-! ## ★出典の紐付け(`.src`) -/

def algebraMap_mem_primeSubring.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(付値環は底変換で付値環に入る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_algebraMap_liesOver.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(分岐した拡大での valAdd のスケーリング。★無条件)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_baseChange_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(良い還元は底変換で保たれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_descend.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(上での Δ_min の関係は下へ降りる——仮説なし。★無条件)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_descend_of_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(上での Δ_min の関係は下へ降りる——分岐指数で割るだけ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_algebraMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(付値は e 倍——v_P(x) = e(P|p)·v_p(x)。★無条件)",
    sectionId := "genell-prop-3-4" }

def degInfOf_baseChange_of_minDeltaExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(deg∞ は基底変換で不変——局所の v_P(Δ_min) = e·v_p(Δ_min) を受けた形)",
    sectionId := "genell-prop-3-4" }

def degInfOf_baseChange_of_minDeltaExp.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆残るのは局所の関係 minDeltaExp P (E×L′) = e(P|p)·minDeltaExp p E である" ++
       "——半安定なら Tate 母数 q の付値がそのまま e 倍になる") 6,
    .citation "[ABC3]" "finsum_scaling(finsum を [L′:L] 倍に落とす、§9-1163)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.finsum_scaling") 2 ]

def htFaltOf_baseChange_of_minDeltaExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(ht^Falt は基底変換で不変——局所の関係を受けた形)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GaloisRep
