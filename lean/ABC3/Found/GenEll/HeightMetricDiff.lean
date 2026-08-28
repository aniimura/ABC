/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightMetric
import ABC3.Found.GenEll.HyperplaneHeight
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★[GenEll] Proposition 1.4, (iii) —— **差の連続性で足りる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> The BD-class of htL depends only on the isomorphism class of the line bundle LQ on XQ.

## ★★★★★★★★★★★★これは何か —— 測定と、その修正

`Found/GenEll/HeightMetric.lean` の `htArith_sub_abs_le` は

    `hD : Continuous D.green`、`hE : Continuous E.green`

を要求する。★★**しかし本物の因子の Green 関数は連続ではない**——
`g(p) = −log‖s(p)‖` は台 `|D|` の上で `+∞` に発散するからである
（2026-08-28 実測）。実際、`§9-871` の `greenFS`（Fubini–Study）は超平面の上で発散する。

★★★したがって `htArith_sub_abs_le` は**そのままでは超平面因子に当たらない**。
原典が言っているのは「同じ直線束の**2 つの計量**の比が有界」であり、
そこで連続なのは**差** `D.green − E.green` の方である
（特異性が打ち消し合う）。

★★★★本ファイルはその形を取る:

    `hcont : Continuous (fun p => D.green p − E.green p)`
      ⟹ `|ht_D − ht_E| ≤ C`（`F` にも点にも依らない）

## ★★★機構 —— 有限側が消えて、和が線型に割れる

★`D.divisor = E.divisor` なので `htArith_eq_add` の**有限側は打ち消し合う**。
★★残るアルキメデス側は `archADiv` が `g` について線型なので

    `Σ_v mult·D.green(archPoint) − Σ_v mult·E.green(archPoint)
       = Σ_v mult·(D.green − E.green)(archPoint)`

となり、`§9-806` の `archADiv_sum_div_finrank_le/ge` がそのまま効く。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

/-! ## ★`archADiv` は Green 関数について線型である -/

/-- ★**`archADiv` の総和は差について線型である**。 -/
theorem archADiv_sum_sub (F : Type) [Field F] [NumberField F] {X : Scheme.{0}}
    (g g' : GreenFn X) (xF : specRingOfIntegers F ⟶ X) :
    (archADiv F g xF).sum (fun _ r => r) - (archADiv F g' xF).sum (fun _ r => r)
      = (archADiv F (fun p => g p - g' p) xF).sum (fun _ r => r) := by
  rw [archADiv_sum_eq, archADiv_sum_eq, archADiv_sum_eq, ← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl (fun v _ => by ring)

/-- ★★**因子が同じなら、高さの差はアルキメデス側だけで書ける**。

★有限側は `hdiv` で打ち消し合う。 -/
theorem htArith_sub_eq (F : Type) [Field F] [NumberField F] {X : Scheme.{0}}
    (D E : ArithCartier X) (hdiv : D.divisor = E.divisor)
    (xF : specRingOfIntegers F ⟶ X) :
    htArith F D xF - htArith F E xF
      = (archADiv F (fun p => D.green p - E.green p) xF).sum (fun _ r => r)
        / (Module.finrank ℚ F : ℝ) := by
  rw [htArith_eq_add F D xF, htArith_eq_add F E xF, hdiv, ← archADiv_sum_sub]
  ring

/-! ## ★★★★★★★★★★★★差の連続性だけで足りる -/

variable {X : Scheme.{0}} {V : Type} [NormedAddCommGroup V] [NormedSpace ℂ V]
  [FiniteDimensional ℂ V]

/-- ★★★★★★★★★★★★**[GenEll] Proposition 1.4, (iii)——差の連続性版**。

原文 (GenEll p.6):
> The BD-class of htL depends only on the isomorphism class of the line bundle LQ on XQ.

因子が同じで**差**が連続なら、`|ht_D − ht_E|` は一様に有界である
——定数は `F` にも点にも依らない。

★★`HeightMetric.lean` の `htArith_sub_abs_le`（各々の連続性を要求する）の**弱仮定版**であり、
**本物の因子に当てられるのはこちらだけ**である
——`g(p) = −log‖s(p)‖` は台の上で発散するので、個別には連続にならない。 -/
theorem htArith_sub_abs_le_of_diff (M : ArcModel X V) [Nonempty (complexPoints X)]
    (D E : ArithCartier X) (hdiv : D.divisor = E.divisor)
    (hcont : @Continuous _ _ M.topology _ (fun p => D.green p - E.green p)) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
        |htArith F D xF - htArith F E xF| ≤ C := by
  obtain ⟨C, hC, hlo, hhi⟩ := M.exists_bound _ hcont
  refine ⟨C, hC, ?_⟩
  intro F _ _ xF
  rw [htArith_sub_eq F D E hdiv xF, abs_le]
  exact ⟨archADiv_sum_div_finrank_ge F _ xF C hC hlo,
    archADiv_sum_div_finrank_le F _ xF C hhi⟩

/-! ## ★出典の紐付け(`.src`) -/

def archADiv_sum_sub.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(archADiv の総和は差について線型である)",
    sectionId := "genell-prop-1-4" }

def htArith_sub_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(因子が同じなら高さの差はアルキメデス側だけ)",
    sectionId := "genell-prop-1-4" }

def htArith_sub_abs_le_of_diff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(差の連続性版——本物の因子に当てられる形)",
    sectionId := "genell-prop-1-4" }

def htArith_sub_abs_le_of_diff.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htArith_sub_abs_le(各々の連続性を要求する版、§9-806)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_sub_abs_le") 2,
    .citation "[ABC3]" "archADiv_sum_div_finrank_le / _ge(アルキメデス側の一様評価)"
      (.inProject "ABC3" "ABC3.Found.GenEll.archADiv_sum_div_finrank_le") 2,
    .implicitStep
      ("★★★★測定(2026-08-28): htArith_sub_abs_le は各 Green 関数の連続性を要求するが、" ++
       "**本物の因子の Green 関数は連続ではない**" ++
       "——g(p) = −log‖s(p)‖ は台 |D| の上で +∞ に発散する。" ++
       "実際 §9-871 の greenFS(Fubini–Study)は超平面の上で発散する。" ++
       "★したがって htArith_sub_abs_le は**そのままでは超平面因子に当たらない**") 4,
    .implicitStep
      ("★★原典が言っているのは「同じ直線束の**2 つの計量**の比が有界」であり、" ++
       "そこで連続なのは**差** D.green − E.green の方である(特異性が打ち消し合う)。" ++
       "★本補題はその形を取る") 3,
    .implicitStep
      ("★★★段 C2c-3 に残るのは、ℙᴺ の ArcModel(complexPoints (ℙᴺ_ℤ) ↪ ℙ(ℂ^{N+1}))を" ++
       "作る段と、Fubini–Study と与えられた計量の**差**が連続であることの確認である") 4 ]

end ABC3.Found.GenEll
