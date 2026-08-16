import ABC3.Meta.Claim
import ABC3.Interface.GenEll.AbcSetup
import ABC3.Found.GenEll.BDClass

/-!
# [GenEll] §2 Bounds on Heights —— Theorem 2.1(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.11。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★この定理の位置づけ —— IUT を経由しない第 2 の頂上

`[IUTchIV]` が [GenEll] を引く回数は `Theorem 2.1` が **10 回**で最多である。
★そして**本定理自身は `[IUTchIV]` を一切使わない**——§1 の枠組みと
noncritical Belyi maps([Mzk1])だけで証明される。

★★**しかしそれは「容易」を意味しない。**
mathlib 実測(2026-08-16): `Belyi` **0 件**、`Arakelov` **0 件**、`LineBundle` **0 件**。
**IUT 非依存であることと、形式化の労力は別の軸である。**

## ★記法の脱落を 2 件実測

★`pdftotext` は (ii) の `ϵ ∈ ℝ_{>0}` を **`∈ R>0`** と出す(`ϵ` が消える)。
★`≲` は**出力に何も残さない**。
★ゆえに**この定理は `.txt` からは書き写せない**。260 dpi 目視から写した。

## ★★向きについて(`Section1.lean` と同じ扱い)

原文 p.5 の定義を逐語で読むと `α ≲ β` は `β ≤ α + C` である。
その読みでは本定理の `ht ≲ (1+ϵ)(log-diff + log-cond)` は
**`(1+ϵ)(…) ≤ ht + C`** となり、abc 予想の向きと**逆**になる。
★**本ファイルは印字どおりに写す**(`BDle`)。
食い違いは `Gap/GenEll/BDDirection.lean` に `GapRecord` と反例つきで記録した。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll ABC3.Found.GenEll

/-- **[GenEll] Theorem 2.1**(Compactly Bounded Subsets and the ABC Conjecture)。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

有限素数集合 `Σ` に対し、次の 2 つは**同値**:

**(i) (Effective Mordell/ABC/Vojta Conjecture)**
`X` を数体上の滑らかな固有幾何的連結曲線、`D ⊆ X` を reduced 因子、`U_X ≝ X\D`、
`d` 正整数、`ϵ ∈ ℝ_{>0}`。`U_X` が双曲的なら、`U_X(ℚ̄)^{≤d}` 上で
`ht_{ω_X(D)} ≲ (1 + ϵ)(log-diff_X + log-cond_D)`。

**(ii) (ABC Conjecture for Σ-Supported Compactly Bounded Subsets)**
`P ≝ ℙ¹_ℚ`、`C` を 3 点 `"0"`, `"1"`, `"∞"` の因子、`U_P ≝ P\C`、
`d` 正整数、`ϵ ∈ ℝ_{>0}`、`K_V ⊆ U_P(ℚ̄)` を support が `Σ` を含む
compactly bounded subset とすると、`K_V ∩ U_P(ℚ̄)^{≤d}` 上で
`ht_{ω_P(C)} ≲ (1 + ϵ)(log-diff_P + log-cond_C)`。

★**(i) ⟹ (ii) は定義から直ちに従う**(原文「immediate from the definitions」)。
実質は **(ii) ⟹ (i)** であり、そこで双曲的曲線の副有限基本群の構造と
`Proposition 1.7` が使われる。 -/
theorem theorem_2_1 (A : AbcSetup) (primes : Finset ℕ) (hprimes : ∀ p ∈ primes, p.Prime) :
    (∀ (X : A.Curve) (dv : (A.data X).Divisor) (d : ℕ) (eps : ℝ),
        A.Reduced X dv → 0 < d → 0 < eps → A.Hyperbolic X dv →
        BDle (fun x : ↥((A.data X).compl dv ∩ (A.data X).degLe d) =>
                (A.data X).ht (A.omegaOf X dv) x.1)
             (fun x : ↥((A.data X).compl dv ∩ (A.data X).degLe d) =>
                (1 + eps) * ((A.data X).logDiff x.1 + (A.data X).logCond dv x.1)))
  ↔
    (∀ (KV : Set (A.data A.projLine).Point) (d : ℕ) (eps : ℝ),
        (A.data A.projLine).CompactlyBounded KV → A.SupportContains KV primes →
        0 < d → 0 < eps →
        BDle (fun x : ↥(KV ∩ (A.data A.projLine).compl A.threePoints
                            ∩ (A.data A.projLine).degLe d) =>
                (A.data A.projLine).ht (A.omegaOf A.projLine A.threePoints) x.1)
             (fun x : ↥(KV ∩ (A.data A.projLine).compl A.threePoints
                            ∩ (A.data A.projLine).degLe d) =>
                (1 + eps) * ((A.data A.projLine).logDiff x.1
                              + (A.data A.projLine).logCond A.threePoints x.1))) := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def theorem_2_1.src : Source :=
  { paper := "GenEll", pdfPage := 11, item := "Theorem 2.1",
    sectionId := "genell-thm-2-1" }

/-- ★原文 p.11–p.13 の証明を通読して数えた。

★**(i) ⟹ (ii) は 1 行**(定義から)。実質はすべて **(ii) ⟹ (i)** の側にある。 -/
def theorem_2_1.needs : List ProofObligation :=
  [ .citation "[GenEll]" "noncritical Belyi maps([Mzk1] = Mochizuki, Noncritical Belyi Maps)"
      (.absent "mathlib に `Belyi` は 0 件(2026-08-16、Mathlib 全体を `Belyi` で検索して 0 件)") 11,
    .otherPaper "[GenEll]" "Proposition 1.7(導手と log-different)" 9,
    .otherPaper "[GenEll]" "Example 1.3, (ii)(compactly bounded subset と support)" 5,
    .otherPaper "[GenEll]" "Remark 1.4.1 / Remark 1.5.1(理論が X_ℚ・(X_ℚ,D_ℚ) だけに依ること)" 8,
    .folklore
      "原文「it follows immediately from the well-known structure of étale fundamental groups of hyperbolic curves over algebraically closed fields of characteristic zero」——任意の正整数 e に対し、E ≝ (D ×_X Y)_red の各点で分岐指数がちょうど e となる連結有限エタール Galois 被覆 U_Y → U_X が存在する。★大きさは未知" 11,
    .implicitStep
      "★statement の語彙(数体上の曲線・標準層 ω_X(D)・双曲性・ℙ¹ の 3 点因子)を Interface/GenEll/AbcSetup.lean に posit した。**我々は作っていない**" 11,
    .implicitStep
      "★★原文 p.5 の ≲ の定義どおりに読むと、本定理の不等式は abc 予想と逆向きになる。**印字どおりに写してある**。Gap/GenEll/BDDirection.lean を参照" 11 ]

end ABC3.Skeleton.GenEll
