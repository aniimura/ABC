/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.PrimeUnder
import ABC3.Found.GaloisRep.FaltingsWitness
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
