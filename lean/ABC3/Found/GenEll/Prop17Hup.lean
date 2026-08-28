/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DifferentCondBound
import ABC3.Found.GenEll.BDClass
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★`Proposition 1.7` の `hup` は仮定でなくなった（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★これは何か —— `BDle` まで運んだ

`Skeleton/GenEll/Section1.lean` の `prop_1_7` は

    `hup : (1 − 1/e)·cond_E x ≤ (diff_Y x − diff_Z x) + slackUp x`
    `hsu : slackUp x ≤ ∑_{q∈Σ} log q`

を**仮定として受けて** `BDle (diff_Y − diff_Z) ((1 − 1/e)·cond_E)` を出していた。

★★★★本ファイルは、点ごとの体の塔 `F_x ⊆ K_x` と分岐データから
**`hup` を slack ゼロで導く**:

    `(1 − 1/e)·(∑_{P 分岐} log N(P)) / [K_x:ℚ]  ≤  log-diff(K_x) − log-diff(F_x)`

★★したがって `BDle` は**定数 `0`** で成り立つ——`slackUp` も `Σ` も要らない。

## ★★★機構

`§9-964`（`one_sub_inv_mul_cond_le_logDiff`）を点ごとに当てるだけである。
★`BDle α β ≡ ∃ C, ∀ x, β x − α x ≤ C` なので、点ごとの `≤` から `C = 0` で出る。

## ★残っている段（明示）

★★★左辺の `(∑_{P 分岐} log N(P)) / [K_x:ℚ]` が原文の **`log-cond_E(x)`** と一致することが残る
——「被覆 `Y → Z` が `E` の上でだけ分岐する」という**幾何の仮定**が要る。
★これが `Proposition 1.7` の**導手の側**であり、different の側は本日
（`§9-954`〜`§9-959`、`§9-964`）閉じた。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★点ごとの不等式から `BDle` -/

/-- ★**点ごとに `β ≤ α` なら `BDle α β`**（定数 `0`）。

★`BDle α β ≡ ∃ C, ∀ x, β x − α x ≤ C` である。 -/
theorem bdle_of_le {P : Type} (α β : P → ℝ) (h : ∀ x, β x ≤ α x) :
    BDle α β := ⟨0, fun x => by have := h x; linarith⟩

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★`hup` を `BDle` まで -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★**`Proposition 1.7` の `hup`**
—— 分岐データから、slack ゼロで。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    点ごとの体の塔 `F_x ⊆ K_x` と分岐する素点の族 `S x`（`e_P ≥ 2`）について
      `BDle (log-diff(K_x) − log-diff(F_x))
            ((1 − 1/e)·(∑_{P∈S x} log N(P)) / [K_x:ℚ])`

★★`Skeleton` の `prop_1_7` が**仮定として受けていた** `hup`／`hsu` が、
ここでは**定理**になっている（しかも定数は `0`）。
★★★残るのは左辺の和が原文の `log-cond_E(x)` と一致すること
——被覆が `E` の上でだけ分岐するという幾何の仮定である。 -/
theorem hup_pointwise (P : Type) (fldF fldK : P → IntermediateField ℚ ℂ)
    (hnfF : ∀ x, NumberField (fldF x)) (hnfK : ∀ x, NumberField (fldK x))
    (alg : ∀ x, Algebra (fldF x) (fldK x))
    (e : ℕ) (he : 0 < e)
    (S : ∀ x, haveI := hnfK x; Finset (Ideal (𝓞 (fldK x))))
    (hmax : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, Pp.IsMaximal)
    (hne : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, Pp ≠ ⊥)
    (pr : ∀ x, haveI := hnfK x; haveI := hnfF x;
      Ideal (𝓞 (fldK x)) → Ideal (𝓞 (fldF x)))
    (ram : ∀ x, haveI := hnfK x; Ideal (𝓞 (fldK x)) → ℕ)
    (hpmax : ∀ x, haveI := hnfK x; haveI := hnfF x; ∀ Pp ∈ S x, (pr x Pp).IsMaximal)
    (hpne : ∀ x, haveI := hnfK x; haveI := hnfF x; ∀ Pp ∈ S x, pr x Pp ≠ ⊥)
    (hdvd : ∀ x, haveI := hnfK x; haveI := hnfF x; letI := alg x; ∀ Pp ∈ S x,
      Pp ^ (ram x Pp) ∣ (pr x Pp).map (algebraMap (𝓞 (fldF x)) (𝓞 (fldK x))))
    (hram2 : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, 2 ≤ ram x Pp)
    (hI : ∀ x, haveI := hnfK x; haveI := hnfF x; letI := alg x;
      differentIdeal (𝓞 (fldF x)) (𝓞 (fldK x)) ≠ 0) :
    BDle
      (fun x => haveI := hnfK x; haveI := hnfF x;
        logDiffOfField (fldK x) - logDiffOfField (fldF x))
      (fun x => haveI := hnfK x;
        (1 - 1 / (e : ℝ)) * (∑ Pp ∈ S x, Real.log (Ideal.absNorm Pp))
          / (Module.finrank ℚ (fldK x) : ℝ)) := by
  refine bdle_of_le _ _ (fun x => ?_)
  haveI := hnfF x
  haveI := hnfK x
  letI := alg x
  exact one_sub_inv_mul_cond_le_logDiff (fldF x) (fldK x) e he (S x) (hmax x) (hne x)
    (pr x) (ram x) (hpmax x) (hpne x) (hdvd x) (hram2 x) (hI x)

/-! ## ★出典の紐付け(`.src`) -/

def bdle_of_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(点ごとに ≤ なら BDle——定数 0)",
    sectionId := "genell-prop-1-7" }

def hup_pointwise.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(hup を分岐データから——slack ゼロ)",
    sectionId := "genell-prop-1-7" }

def hup_pointwise.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "one_sub_inv_mul_cond_le_logDiff((1 − 1/e)·log-cond ≤ log-diff、§9-964)"
      (.inProject "ABC3" "ABC3.Found.GenEll.one_sub_inv_mul_cond_le_logDiff") 4,
    .implicitStep
      ("★★★★★測定(2026-08-29): Skeleton の prop_1_7 が**仮定として受けていた** hup／hsu は、" ++
       "点ごとの体の塔と分岐データからは**定理である**——しかも slack はゼロで足りる。" ++
       "★Σ の上の食い違いは different の側には現れない") 5,
    .implicitStep
      ("★★★残るのは左辺の和 (∑_{P 分岐} log N(P))/[K_x:ℚ] が" ++
       "原文の log-cond_E(x) と一致することである" ++
       "——「被覆 Y → Z が E の上でだけ分岐する」という幾何の仮定が要る。" ++
       "★これが Proposition 1.7 の**導手の側**であり、" ++
       "different の側は §9-954〜959・§9-964 で閉じた") 5 ]

end ABC3.Found.GenEll
