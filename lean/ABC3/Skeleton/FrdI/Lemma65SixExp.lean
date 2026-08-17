import ABC3.Meta.Claim
import Mathlib.RingTheory.Algebraic.Defs
import Mathlib.Analysis.SpecialFunctions.Complex.Analytic

/-!
# [FrdI] Lemma 6.5, (ii) —— six exponentials theorem(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> Lemma 6.5.

原文 (FrdI p.116):
> linearly independent over Q.

## ★★このファイルの位置づけ —— 「壁」を割って葉を出すために作った

★`Lemma 6.5, (i)` は実装済み(`Found/FrdI/Lemma65.lean` の
`log_primes_linearIndependent`)。(ii) は原文が Lang の定理に送っており、
それは **six exponentials theorem**(6 指数定理)である。

★2026-08-17 まで、ここは「mathlib に超越論の在庫が無い」の一言で止めていた。
CLAUDE.md の**姿勢**に従い、**statement を型で固定して、そこへ至る道を測る**。

## ★★★測定の訂正(2026-08-18)—— 前日の「在庫ゼロ」は探索先が存在しなかった

★2026-08-17 の `Gap/FrdI/Section6.lean` には
「mathlib の `Analysis/Transcendental/` には Liouville 数と `e` の超越性しか無い」
と書いた。★★**mathlib に `Analysis/Transcendental/` というディレクトリは無い。**
実体は **`NumberTheory/Transcendental/`** である。**無いものを探して 0 件を得ていた。**

★★改めて測った在庫(2026-08-18、探索範囲つき):

| 要るもの | mathlib | 判定 |
|---|---|---|
| **Siegel の補題**(整数版) | `Int.Matrix.exists_ne_zero_int_vec_norm_le`(`NumberTheory/SiegelsLemma.lean:181`) | ★★**ある** |
| **数体上の Siegel の補題**(house 版) | `NumberField.house.exists_ne_zero_int_vec_house_le`(`NumberTheory/NumberField/House.lean`) | ★★**ある**(超越性証明の**算術側の心臓部**) |
| 代数的数の house とその評価 | `NumberField.house` / `house_mul_le` / `house_pow_le` / `norm_embedding_le_house` | ★**ある** |
| 0 でない代数的整数の共役下界 | `NumberField.exists_conjugate_one_le_norm` | ★**ある** |
| 最大値原理 | `Complex.norm_le_of_forall_mem_frontier_norm_le`(`Analysis/Complex/AbsMax.lean:400`) | ★**ある** |
| Schwarz の補題(1 点・高位) | `Complex.dist_le_mul_div_pow_of_mapsTo_ball_of_isLittleO`(`Analysis/Complex/Schwarz.lean:144`) | ★**ある** |
| **多零点版の小ささ評価** | `Analysis/Complex/` を目視して該当なし(2026-08-18) | ★**無い** |
| Schneider–Lang / six exponentials | `NumberTheory/Transcendental/` は `Liouville/` と `Lindemann/AnalyticalPart.lean` のみ | ★**無い** |

★★**すなわち「解析的整数論の一分野を入れる作業」ではなかった。**
いちばん重い算術側(数体上の Siegel の補題)は既に在庫にある。
分解の DAG は `ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン
(8 節点・6 層・葉 3 個)。

## ★影響範囲

★`Lemma 6.5, (ii)` は `§6` の**例**(`Example 6.3` の非同型性)を言うために使われ、
Frobenioid の理論そのものには入らない。下流を止めているのは自分 1 件だけである。
-/

namespace ABC3.Skeleton.FrdI

open ABC3.Meta

/-- ★★★**six exponentials theorem**(Lang / Ramachandra)。

`x₁, x₂` が ℚ 上一次独立、`y₁, y₂, y₃` が ℚ 上一次独立な複素数なら、
`exp(xᵢyⱼ)`(6 個)の少なくとも 1 つは超越数である。

★原文 `Lemma 6.5, (ii)` はこれを `{log p₂, log p₄, log p₆}` と
`{1, log p₃/log p₄}` に当てる。 -/
theorem six_exponentials {x : Fin 2 → ℂ} {y : Fin 3 → ℂ}
    (_hx : LinearIndependent ℚ x) (_hy : LinearIndependent ℚ y) :
    ∃ (i : Fin 2) (j : Fin 3), Transcendental ℚ (Complex.exp (x i * y j)) := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def six_exponentials.src : Source :=
  { paper := "FrdI", pdfPage := 116, item := "Lemma 6.5, (ii)",
    sectionId := "frdi-lemma-6-5" }

/-- ★★分解の DAG は `ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
`node tools/frdi-newleaves.mjs` が層と葉を印字する。 -/
def six_exponentials.needs : List ProofObligation :=
  [ .citation "[Lang1]" "Introduction to Transcendental Numbers(原文が送っている先。Baker の本 p.119 も指す)"
      (.absent "0_Source に [Lang1] は無い(papers.json 未登記、2026-08-18)") 116,
    .citation "[mathlib]" "数体上の Siegel の補題(house 版)——補助関数の係数を取る所"
      (.inMathlib "NumberField.house.exists_ne_zero_int_vec_house_le") 116,
    .citation "[mathlib]" "Siegel の補題(整数版)"
      (.inMathlib "Int.Matrix.exists_ne_zero_int_vec_norm_le") 116,
    .citation "[mathlib]" "0 でない代数的整数はノルム 1 以上の共役を持つ(Liouville 不等式の土台)"
      (.inMathlib "NumberField.exists_conjugate_one_le_norm") 116,
    .citation "[mathlib]" "最大値原理(多零点版の小ささ評価の土台)"
      (.inMathlib "Complex.norm_le_of_forall_mem_frontier_norm_le") 116,
    .citation "[mathlib]" "多点で消える正則関数の小ささ(Schwarz の補題の多零点版)"
      (.absent "Analysis/Complex/ の宣言一覧を目視、1 点の高位版(Schwarz.lean:144)はあるが多点版は無い(2026-08-18)") 116,
    .derivation
      "補助関数 F(z) = Σ c_pq exp((p x1 + q x2) z) を Siegel で作り、増大度評価と Liouville 不等式で外挿の帰納を回す(frdi-decomposition.json の sixexp チェーン、葉 3 個)" 116,
    .implicitStep
      "★原文は (ii) の証明を置かず、Lang の定理に送るだけである。適用先の一次独立性(log p の側)は (i) として実装済み(Found/FrdI/Lemma65.lean)" 116 ]

end ABC3.Skeleton.FrdI
