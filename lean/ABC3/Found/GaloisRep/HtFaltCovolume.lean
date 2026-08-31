/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.LogProductFormula
import ABC3.Found.GaloisRep.FaltingsWitness
import ABC3.Found.GenEll.ArchInvCovolume
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★`ht^Falt` の共体積表示（`Found`、無条件）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★★★これは何か

    **`12·d·ht^Falt(E) = −12·Σᶠ_p neronExp_p·log N(p) − 12·d·log(2π)`**
    **`                  + 12·Σ_σ log‖u_σ‖ − 6·Σ_σ log(covol P_σ)`**

（`d = [L:ℚ]`、`C_σ • (E ⊗ σ) = latticeCurve P_σ`、`u_σ = (C_σ).u`）。

★`§9-1020`（第 576）が**紙の上で**導いた式の、Lean 上の実体である。
★★`deg∞`（有限素点の `minDeltaExp`）と `archSum`（アルキメデスの `curveArchInv`）が
`Σ_σ log‖σΔ_E‖` を通じて**打ち消し合い**、`neronExp`・`u`・`covol` だけが残る。

## ★機構 —— 3 つの道具

| 道具 | 出どころ |
|---|---|
| 対数版の積公式 `Σ_σ log‖σx‖ = Σᶠ_p v_p(x)·log N(p)` | ★`§9-1021`（第 577） |
| `log(archNorm) = log‖σΔ‖ − 12log‖u‖ + 6log covol` | ★`§9-1022`（第 578） |
| `minDeltaExp = v_p(Δ) − 12·neronExp` | ★`§9-319`（定義） |

★★★どれも本プロジェクトで証明済みであり、本ファイルは**組み立てだけ**である。

## ★★これが何を意味するか —— 同種写像評価の残りが確定した

`Lemma 3.5` に残る唯一の入力 `ht^Falt(E/H) ≤ ht^Falt(E) + 2·log(l)` は、
上の表示のもとで**差を取ると**

    `12·d·(ht^Falt(E′) − ht^Falt(E))`
      `= −12·Σᶠ_p [neronExp_p(E′) − neronExp_p(E)]·log N(p)`
        `+ 12·Σ_σ log‖u′_σ/u_σ‖ − 6·Σ_σ log(covol P′_σ / covol P_σ)`

になる（`log(2π)` の項は打ち消す）。

★**第 3 項は取れている**——解析的な同種写像で `covol` は指数 `l` で `1/l` 倍
（`§9-1017`・`§9-1019`）。
☆残るのは**第 1 項（`neronExp` の差）と第 2 項（`u` の比）**であり、
それがちょうど `[DelSB616] Théorème 2.4` の段 1・段 2
（`ω_{E′} ⊆ ω_E` で余核が `#H` で消える）である。

★★★★**目標は完全に局所化した**: 各素点 `p` で `neronExp` の差が `log(l)` の程度、
かつアルキメデスで `u` の比が制御されること。
-/

namespace ABC3.Found.GaloisRep

open NumberField IsDedekindDomain WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★有限素点側 -/

/-- ★★★★★★**`d·deg∞ = Σᶠ_p v_p(Δ)·log N(p) − 12·Σᶠ_p neronExp_p·log N(p)`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`minDeltaExp p E = v_p(Δ_E) − 12·neronExp_p E`（`§9-319` の定義）を
有限台の和に流すだけである。 -/
theorem finrank_degInfOf_eq (E : WeierstrassCurve L) (hΔ : E.Δ ≠ 0) :
    (Module.finrank ℚ L : ℝ) * degInfOf L E
      = (∑ᶠ p : HeightOneSpectrum (𝓞 L),
          (valAdd p (Units.mk0 E.Δ hΔ) : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        - 12 * (∑ᶠ p : HeightOneSpectrum (𝓞 L),
          (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal)) := by
  have hfinM : (Function.support (fun p : HeightOneSpectrum (𝓞 L) =>
      (minDeltaExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))).Finite := by
    refine (minDeltaExp_finite E hΔ).subset (fun p hp => ?_)
    simp only [Function.mem_support, ne_eq, mul_eq_zero, not_or] at hp
    simp only [Set.mem_setOf_eq]
    exact_mod_cast hp.1
  have hfinN : (Function.support (fun p : HeightOneSpectrum (𝓞 L) =>
      (12:ℝ) * ((neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal)))).Finite := by
    refine (finite_bad_primes' E hΔ).subset (fun p hp => ?_)
    simp only [Function.mem_support, ne_eq, mul_eq_zero, not_or] at hp
    simp only [Set.mem_setOf_eq]
    exact_mod_cast hp.2.1
  have hsplit : (∑ᶠ p : HeightOneSpectrum (𝓞 L),
      (valAdd p (Units.mk0 E.Δ hΔ) : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
      = (∑ᶠ p : HeightOneSpectrum (𝓞 L),
          (minDeltaExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        + ∑ᶠ p : HeightOneSpectrum (𝓞 L),
          (12:ℝ) * ((neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal)) := by
    rw [← finsum_add_distrib hfinM hfinN]
    refine finsum_congr (fun p => ?_)
    rw [minDeltaExp, dif_neg hΔ]
    push_cast
    ring
  rw [finrank_mul_degInfOf, hsplit, mul_finsum]
  ring

/-! ## ★★★★★★アルキメデス側 -/

/-- ★★★★★★**`archSum` の分解**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`§9-1022` の `log_archNorm_eq` を `σ` について足すだけである。 -/
theorem archSum_eq (E : WeierstrassCurve L) [E.IsElliptic]
    (P : (L →+* ℂ) → PeriodPair) (C : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ)) :
    archSum L E
      = 12 * (Module.finrank ℚ L : ℝ) * Real.log (2 * Real.pi)
        + (∑ σ : (L →+* ℂ), Real.log ‖σ E.Δ‖)
        - 12 * (∑ σ : (L →+* ℂ), Real.log ‖((C σ).u : ℂ)‖)
        + 6 * (∑ σ : (L →+* ℂ), Real.log (covol (P σ))) := by
  have hpi : (0:ℝ) < (2 * Real.pi) ^ 12 := by positivity
  have hterm : ∀ σ : (L →+* ℂ), Real.log ((2 * Real.pi) ^ 12 * archNorm E σ)
      = 12 * Real.log (2 * Real.pi)
        + (Real.log ‖σ E.Δ‖ - 12 * Real.log ‖((C σ).u : ℂ)‖ + 6 * Real.log (covol (P σ))) := by
    intro σ
    rw [Real.log_mul (ne_of_gt hpi) (ne_of_gt (archNorm_pos E σ)),
      log_archNorm_eq E σ (hPC σ), Real.log_pow]
    push_cast
    ring
  rw [archSum, Finset.sum_congr rfl (fun σ _ => hterm σ)]
  rw [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ,
    NumberField.Embeddings.card L ℂ, nsmul_eq_mul]
  rw [show (fun σ : (L →+* ℂ) => Real.log ‖σ E.Δ‖ - 12 * Real.log ‖((C σ).u : ℂ)‖
        + 6 * Real.log (covol (P σ)))
      = fun σ => (Real.log ‖σ E.Δ‖ - 12 * Real.log ‖((C σ).u : ℂ)‖)
        + 6 * Real.log (covol (P σ)) from rfl]
  rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum]
  ring

/-! ## ★★★★★★★★★★★★★★★★★★★★共体積表示 -/

/-- ★★★★★★★★★★★★★★★★★★★★**`ht^Falt` の共体積表示**——★**無条件**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

    `12·d·ht^Falt(E) = −12·Σᶠ_p neronExp_p·log N(p) − 12·d·log(2π)`
    `                  + 12·Σ_σ log‖u_σ‖ − 6·Σ_σ log(covol P_σ)`

★`deg∞` と `archSum` が `Σ_σ log‖σΔ_E‖` を通じて**打ち消し合う**のが要点である
（対数版の積公式、`§9-1021`）。

★★★これで同種写像評価 `ht^Falt(E/H) ≤ ht^Falt(E) + 2log(l)` に残るのは
**`neronExp` の差と `u` の比**だけになった——`covol` の比は
`§9-1017`・`§9-1019` で取れている。 -/
theorem twelve_finrank_htFaltOf_eq (E : WeierstrassCurve L) [E.IsElliptic]
    (P : (L →+* ℂ) → PeriodPair) (C : (L →+* ℂ) → VariableChange ℂ)
    (hPC : ∀ σ, C σ • (E.map σ) = latticeCurve (P σ)) :
    12 * (Module.finrank ℚ L : ℝ) * htFaltOf L E
      = -12 * (∑ᶠ p : HeightOneSpectrum (𝓞 L),
            (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        - 12 * (Module.finrank ℚ L : ℝ) * Real.log (2 * Real.pi)
        + 12 * (∑ σ : (L →+* ℂ), Real.log ‖((C σ).u : ℂ)‖)
        - 6 * (∑ σ : (L →+* ℂ), Real.log (covol (P σ))) := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hΔ : E.Δ ≠ 0 := (inferInstance : E.IsElliptic).isUnit.ne_zero
  have hprod := sum_arch_log_eq_finsum_valAdd E.Δ hΔ
  have hdeg := finrank_degInfOf_eq E hΔ
  have harch := archSum_eq E P C hPC
  have hht : 12 * (Module.finrank ℚ L : ℝ) * htFaltOf L E
      = (Module.finrank ℚ L : ℝ) * degInfOf L E - archSum L E := by
    rw [htFaltOf]
    field_simp
  rw [hht, hdeg, harch, ← hprod]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def finrank_degInfOf_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(d·deg∞ を v_p(Δ) と neronExp に分解)",
    sectionId := "genell-prop-3-4" }

def archSum_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(archSum を log‖σΔ‖・log‖u‖・log covol に分解)",
    sectionId := "genell-prop-3-4" }

def twelve_finrank_htFaltOf_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(ht^Falt の共体積表示——deg∞ と archSum が打ち消し合う。★無条件)",
    sectionId := "genell-prop-3-4" }

def twelve_finrank_htFaltOf_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "sum_arch_log_eq_finsum_valAdd(対数版の積公式、§9-1021)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_arch_log_eq_finsum_valAdd") 3,
    .citation "[ABC3]" "log_archNorm_eq(archNorm の共体積分解、§9-1022)"
      (.inProject "ABC3" "ABC3.Found.GenEll.log_archNorm_eq") 3,
    .citation "[ABC3]" "minDeltaExp(v_p(Δ) − 12·neronExp、§9-319)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp") 2,
    .implicitStep
      ("★★★★★★★★★★到達点(2026-08-29、第 579): §9-1020(第 576)が**紙の上で**" ++
       "導いた共体積表示を Lean で証明した。" ++
       "★deg∞(有限素点)と archSum(アルキメデス)が Σ_σ log‖σΔ_E‖ を通じて" ++
       "**打ち消し合う**のが要点であり、それが対数版の積公式(§9-1021)である。" ++
       "★★残るのは neronExp・u・covol の 3 項") 9,
    .implicitStep
      ("☆同種写像評価への帰結: 上の表示で差を取ると " ++
       "12·d·(ht^Falt(E′) − ht^Falt(E)) は " ++
       "(1) neronExp の差 (2) u の比 (3) covol の比 の 3 項になる" ++
       "(log(2π) の項は打ち消す)。" ++
       "★**(3) は取れている**——解析的な同種写像で covol は指数 l で 1/l 倍" ++
       "(§9-1017・§9-1019)。" ++
       "☆残る (1)(2) がちょうど [DelSB616] Théorème 2.4 の段 1・段 2 である" ++
       "(ω_{E′} ⊆ ω_E で余核が #H で消える)") 9 ]

end ABC3.Found.GaloisRep
