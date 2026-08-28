/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Prop17HupCond
import ABC3.Found.GenEll.DifferentTameExact
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★`hlow` が大域で出た —— 局所化は要らなかった（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 道筋が変わった

第 500-502 ブロックでは「`hlow` の残り 1 本は `§9-959`（局所）を大域へ上げる段であり、
`differentIdeal` と局所化の両立（**mathlib に無い**）が要る」と記録していた。

★★★★★**第 503 ブロックの測定でそれが覆った**:
mathlib の `not_dvd_differentIdeal_of_intTrace_not_mem` は

    `{p : Ideal A} (P Q : Ideal B) (hP : P * Q = p.map) (x ∈ Q)
     (hx : intTrace A B x ∉ p) : ¬ P ∣ 𝔡`

であり、★**`P` に `IsMaximal` を要求していない**。
★★したがって `P ≔ 𝔭^e`・`Q ≔ ⊤` と取れば、**大域のまま** `¬ 𝔭^e ∣ 𝔡` が出る。

## ★★★機構 —— 完全分岐なら witness は `1` である

`p·𝓞_K = 𝔭^e`（完全分岐）なら `Q = ⊤` に取れるので、`x ≔ 1` が使える。
★`intTrace(1) = [K:F] = e` なので、条件 `intTrace(1) ∉ p` は
**まさに馴分岐（`p ∤ e`）**である。

★★★★★局所化も、`§9-959` の Eisenstein の議論も**要らなかった**
——mathlib の補題の形が既に大域だったのである。

## ★これで `Proposition 1.7` の different の側は

| 向き | 状態 |
|---|---|
| `hup`（右側の `≲`） | ★★★`§9-964`〜`§9-967` で**原文の `log-cond` の形で定理** |
| `hlow`（左側の `≲`） | ★★★**本ファイルで大域の定理**（完全分岐・馴のとき） |

★★残るのは「一般の分岐（完全分岐でない場合）」と、原文 (ii) の Riemann–Hurwitz である。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★★★★★★★完全分岐・馴なら `𝔭^e ∤ 𝔡`（大域） -/

/-- ★★★★★★★★★★★★★★**完全分岐で馴なら `𝔭^e ∤ 𝔡`**——大域のまま。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

★`p·𝓞_K = 𝔭^e` なら `Q = ⊤` に取れるので witness は `x = 1` でよい。
★★`intTrace(1) = [K:F] = e` だから、条件 `intTrace(1) ∉ p` は**馴分岐そのもの**である。
★★★局所化は要らない——mathlib の `not_dvd_differentIdeal_of_intTrace_not_mem` は
`P` に極大性を要求していないからである（2026-08-29 実測）。 -/
theorem not_pow_dvd_differentIdeal_of_totallyRamified_tame (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    {p : Ideal (𝓞 F)} (Pp : Ideal (𝓞 K)) (e : ℕ)
    (hfac : Pp ^ e = Ideal.map (algebraMap (𝓞 F) (𝓞 K)) p)
    (htame : (Algebra.intTrace (𝓞 F) (𝓞 K)) 1 ∉ p) :
    ¬ Pp ^ e ∣ differentIdeal (𝓞 F) (𝓞 K) := by
  refine not_dvd_differentIdeal_of_intTrace_not_mem (𝓞 F) (Pp ^ e) ⊤ ?_ 1 trivial htame
  rw [Ideal.mul_top]
  exact hfac

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★`hlow` の大域形 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★**`hlow` が大域で出た**
——分岐がすべて完全分岐・馴のとき。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `log-diff(K) − log-diff(F) ≤ (∑_{P∈S} (e_P − 1)·log N(P)) / [K:ℚ]`

★`§9-957`（分岐の言葉での上からの評価）に、本ファイルの `hbound` を渡すだけである。
★★★これで `Proposition 1.7` の different の側は**両方向とも大域の定理**になった。 -/
theorem logDiff_tower_le_sum_of_totallyRamified_tame (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] [DecidableEq (Ideal (𝓞 K))]
    (S : Finset (Ideal (𝓞 K))) (hne : ∀ P ∈ S, P ≠ ⊥)
    (pr : Ideal (𝓞 K) → Ideal (𝓞 F)) (ram : Ideal (𝓞 K) → ℕ)
    (hram1 : ∀ P ∈ S, 1 ≤ ram P)
    (hfac : ∀ P ∈ S, P ^ (ram P) = Ideal.map (algebraMap (𝓞 F) (𝓞 K)) (pr P))
    (htame : ∀ P ∈ S, (Algebra.intTrace (𝓞 F) (𝓞 K)) 1 ∉ pr P)
    (hI : differentIdeal (𝓞 F) (𝓞 K) ≠ 0)
    (hprodne : (∏ P ∈ S, P ^ (ram P - 1)) ≠ 0)
    (hramsupp : ∀ Q : Ideal (𝓞 K), ∀ _ : Q.IsPrime,
      ¬ Algebra.IsUnramifiedAt (𝓞 F) Q → Q ∈ S) :
    logDiffOfField K - logDiffOfField F
      ≤ (∑ P ∈ S, ((ram P - 1 : ℕ) : ℝ) * Real.log (Ideal.absNorm P))
        / (Module.finrank ℚ K : ℝ) := by
  refine logDiff_tower_le_sum_of_ramified F K S hne (fun P => ram P - 1) hI hprodne hramsupp ?_
  intro P hP
  have he : ram P - 1 + 1 = ram P := by have := hram1 P hP; omega
  simp only [he]
  exact not_pow_dvd_differentIdeal_of_totallyRamified_tame F K P (ram P) (hfac P hP) (htame P hP)

/-! ## ★出典の紐付け(`.src`) -/

def not_pow_dvd_differentIdeal_of_totallyRamified_tame.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(完全分岐で馴なら 𝔭^e ∤ 𝔡——大域)",
    sectionId := "genell-prop-1-7" }

def logDiff_tower_le_sum_of_totallyRamified_tame.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(hlow の大域形——分岐がすべて完全分岐・馴のとき)",
    sectionId := "genell-prop-1-7" }

def logDiff_tower_le_sum_of_totallyRamified_tame.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]"
      "not_dvd_differentIdeal_of_intTrace_not_mem(P に極大性を要求しない——2026-08-29 実測)"
      (.inMathlib "IsDedekindDomain.not_dvd_differentIdeal_of_intTrace_not_mem") 4,
    .citation "[ABC3]" "logDiff_tower_le_sum_of_ramified(分岐の言葉での上からの評価、§9-957)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiff_tower_le_sum_of_ramified") 3,
    .implicitStep
      ("★★★★★測定(2026-08-29): 第 500-502 では「hlow の残りは §9-959(局所)を" ++
       "大域へ上げる段で、differentIdeal と局所化の両立(mathlib に無い)が要る」と" ++
       "記録していたが、**それは覆った**。" ++
       "mathlib の not_dvd_differentIdeal_of_intTrace_not_mem は P に極大性を要求せず、" ++
       "P ≔ 𝔭^e・Q ≔ ⊤ と取れば大域のまま ¬ 𝔭^e ∣ 𝔡 が出る" ++
       "——完全分岐なら witness は x = 1 でよく、intTrace(1) = [K:F] = e なので" ++
       "条件はまさに馴分岐である") 6,
    .implicitStep
      ("★★これで Proposition 1.7 の different の側は**両方向とも大域の定理**になった" ++
       "(hup は §9-964〜967 で原文の log-cond の形、hlow は本ファイル)。" ++
       "★残るのは一般の分岐(完全分岐でない場合)と原文 (ii) の Riemann–Hurwitz である") 5 ]

end ABC3.Found.GenEll
