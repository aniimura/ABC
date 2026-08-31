/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LogCondSupport
import ABC3.Found.GenEll.Prop17Hup
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★`hup` が原文の `log-cond` で書けた（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— `hup` の到達形

`§9-965` は `hup` を **`(∑_{P∈S} log N(P))/[K:ℚ]`** の形で出し、
`§9-966` はその和が **`log-cond`** であることを示した。
★★★★本ファイルは 2 つを繋ぐ:

    `BDle (log-diff(K_x) − log-diff(F_x))  ((1 − 1/e)·log-cond_E(x))`

★★これは原文の `Proposition 1.7, (i)` の**右側の `≲` そのもの**である
（`Skeleton` の `prop_1_7` の第 2 成分）。★定数は `0`。

## ★★★受けている仮定（＝幾何の紐、明示）

点ごとの素点の族 `S x` が**同時に**次の 2 つであること:

| 側 | 条件 |
|---|---|
| 導手 | `S x` が `E` の引き戻しの**台**である（`hmem`・`hsupp`） |
| 分岐 | `S x` の各素点が `K_x/F_x` で **`e_P ≥ 2` で分岐**する（`hdvd`・`hram2`） |

★★★これが原文の「被覆 `Y → Z` が `E` の上でだけ分岐する」の**型での言い方**である。
★他に受けているものは無い——`slackUp` も `Σ` も要らない。

## ★残っている段（明示）

★★`Proposition 1.7` に残るのは **`hlow`（左側の `≲`）の局所→大域**だけである
——`§9-959`（馴分岐で `𝔡 = (λ^{e−1})`）は局所の主張であり、
大域へ上げるには `differentIdeal` と局所化の両立が要る（**mathlib に無い**、
`mathlib-gap.json` の `tame-ramification.measured20260829.newMathlibGap` に記録）。
-/

namespace ABC3.Found.GenEll

open NumberField AlgebraicGeometry

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★`hup` の到達形 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**`Proposition 1.7` の `hup`**
—— 原文の `log-cond` で書いた形。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `BDle (log-diff(K_x) − log-diff(F_x))  ((1 − 1/e)·log-cond_E(x))`

★★受けているのは「点ごとの素点の族 `S x` が**導手の台**であり、かつ
その各素点が **`e_P ≥ 2` で分岐**する」ことだけである
——原文の「被覆が `E` の上でだけ分岐する」の型での言い方。
★★★定数は `0`（`slackUp` も `Σ` も要らない）。 -/
theorem hup_logCond_pointwise {X : Scheme.{0}} (E : X.IdealSheafData)
    (P : Type) (fldF fldK : P → IntermediateField ℚ ℂ)
    (hnfF : ∀ x, NumberField (fldF x)) (hnfK : ∀ x, NumberField (fldK x))
    (alg : ∀ x, Algebra (fldF x) (fldK x))
    (deceq : ∀ x, haveI := hnfK x; DecidableEq (Ideal (𝓞 (fldK x))))
    (e : ℕ) (he : 0 < e)
    (xE : ∀ x, haveI := hnfK x; specRingOfIntegers (fldK x) ⟶ X)
    (S : ∀ x, haveI := hnfK x; Finset (Ideal (𝓞 (fldK x))))
    (hmax : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, Pp.IsMaximal)
    (hne : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, Pp ≠ ⊥)
    -- ★導手の側: `S x` は `E` の引き戻しの台である
    (hI : ∀ x, haveI := hnfK x; pullbackIdeal (fldK x) E (xE x) ≠ 0)
    (hmem : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, Pp ∣ pullbackIdeal (fldK x) E (xE x))
    (hsupp : ∀ x, haveI := hnfK x; ∀ Q : Ideal (𝓞 (fldK x)), Irreducible Q →
      Q ∣ pullbackIdeal (fldK x) E (xE x) → Q ∈ S x)
    -- ★分岐の側: `S x` の素点は `e_P ≥ 2` で分岐する
    (pr : ∀ x, haveI := hnfK x; haveI := hnfF x;
      Ideal (𝓞 (fldK x)) → Ideal (𝓞 (fldF x)))
    (ram : ∀ x, haveI := hnfK x; Ideal (𝓞 (fldK x)) → ℕ)
    (hpmax : ∀ x, haveI := hnfK x; haveI := hnfF x; ∀ Pp ∈ S x, (pr x Pp).IsMaximal)
    (hpne : ∀ x, haveI := hnfK x; haveI := hnfF x; ∀ Pp ∈ S x, pr x Pp ≠ ⊥)
    (hdvd : ∀ x, haveI := hnfK x; haveI := hnfF x; letI := alg x; ∀ Pp ∈ S x,
      Pp ^ (ram x Pp) ∣ (pr x Pp).map (algebraMap (𝓞 (fldF x)) (𝓞 (fldK x))))
    (hram2 : ∀ x, haveI := hnfK x; ∀ Pp ∈ S x, 2 ≤ ram x Pp)
    (hD : ∀ x, haveI := hnfK x; haveI := hnfF x; letI := alg x;
      differentIdeal (𝓞 (fldF x)) (𝓞 (fldK x)) ≠ 0) :
    BDle
      (fun x => haveI := hnfK x; haveI := hnfF x;
        logDiffOfField (fldK x) - logDiffOfField (fldF x))
      (fun x => haveI := hnfK x;
        (1 - 1 / (e : ℝ)) * logCond (fldK x) E (xE x)) := by
  refine bdle_of_le _ _ (fun x => ?_)
  haveI := hnfF x
  haveI := hnfK x
  haveI := deceq x
  letI := alg x
  rw [logCond_eq_sum_log (fldK x) E (xE x) (hI x) (S x) (hmax x) (hne x) (hmem x) (hsupp x),
    ← mul_div_assoc]
  exact one_sub_inv_mul_cond_le_logDiff (fldF x) (fldK x) e he (S x) (hmax x) (hne x)
    (pr x) (ram x) (hpmax x) (hpne x) (hdvd x) (hram2 x) (hD x)

/-! ## ★出典の紐付け(`.src`) -/

def hup_logCond_pointwise.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(hup を原文の log-cond で書いた形)",
    sectionId := "genell-prop-1-7" }

def hup_logCond_pointwise.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "one_sub_inv_mul_cond_le_logDiff((1 − 1/e)·∑ ≤ log-diff、§9-964)"
      (.inProject "ABC3" "ABC3.Found.GenEll.one_sub_inv_mul_cond_le_logDiff") 4,
    .citation "[ABC3]" "logCond_eq_sum_log(導手の等式、§9-966)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logCond_eq_sum_log") 4,
    .implicitStep
      ("★★★★★測定(2026-08-29): Proposition 1.7 の hup(右側の ≲)は" ++
       "**原文の log-cond で書けた**。受けているのは" ++
       "「点ごとの素点の族 S x が導手の台であり、かつ各素点が e_P ≥ 2 で分岐する」" ++
       "——原文の「被覆が E の上でだけ分岐する」の型での言い方——だけであり、" ++
       "★slackUp も Σ も要らない(定数は 0)") 6,
    .implicitStep
      ("★★残るのは hlow(左側の ≲)の**局所→大域**だけである" ++
       "——§9-959(馴分岐で 𝔡 = (λ^{e−1}))は局所の主張であり、" ++
       "大域へ上げるには differentIdeal と局所化の両立が要る" ++
       "(mathlib に無い、mathlib-gap.json に記録済み)") 5 ]

end ABC3.Found.GenEll
