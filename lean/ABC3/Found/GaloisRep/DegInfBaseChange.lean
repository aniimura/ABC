/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.PrimeUnder
import ABC3.Found.GaloisRep.FaltingsWitness
import ABC3.Found.GaloisRep.SemistableFin
import Mathlib.NumberTheory.RamificationInertia.Valuation
import ABC3.Meta.Claim

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
    (hp : SemistableAt p E) (hP : SemistableAt P (E.baseChange L')) :
    minDeltaExp P (E.baseChange L')
      = (p.asIdeal.ramificationIdx P.asIdeal : ℤ) * minDeltaExp p E := by
  rw [minDeltaExp_eq_maxJ_of_semistable P _ hP, minDeltaExp_eq_maxJ_of_semistable p E hp,
    jExp_baseChange L L' p P E, ← mul_neg,
    max_zero_mul _ (by positivity : (0:ℤ) ≤ (p.asIdeal.ramificationIdx P.asIdeal : ℤ))]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★到達点 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**半安定な楕円曲線の `deg∞` は基底変換で変わらない**——★**無条件**。 -/
theorem degInfOf_baseChange_of_semistable (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    (hss' : ∀ P : HeightOneSpectrum (𝓞 L'), SemistableAt P (E.baseChange L'))
    [IsScalarTower ℚ L L'] :
    degInfOf L' (E.baseChange L') = degInfOf L E := by
  refine degInfOf_baseChange_of_minDeltaExp L L' E E.isUnit_Δ.ne_zero (fun P => ?_)
  exact_mod_cast congrArg (fun n : ℤ => (n : ℝ))
    (minDeltaExp_baseChange_of_semistable L L' _ P E (hss _) (hss' P))

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**半安定な楕円曲線の `ht^Falt` は基底変換で変わらない**——★**無条件**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★これで `EllModuliData` の `faltingsHeight : EllClass → ℝ` を作る道の
**体を大きくする向きが完全に通った**（残るのは捻り＝同じ `j` の別の曲線の扱いである）。 -/
theorem htFaltOf_baseChange_of_semistable (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    (hss' : ∀ P : HeightOneSpectrum (𝓞 L'), SemistableAt P (E.baseChange L'))
    [IsScalarTower ℚ L L'] :
    htFaltOf L' (E.baseChange L') = htFaltOf L E :=
  htFaltOf_baseChange_of_degInf L L' E
    (degInfOf_baseChange_of_semistable L L' E hss hss')

def htFaltOf_baseChange_of_semistable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(半安定な曲線の ht^Falt は基底変換で不変。★無条件)",
    sectionId := "genell-prop-3-4" }

end NumberField

/-! ## ★出典の紐付け(`.src`) -/

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
