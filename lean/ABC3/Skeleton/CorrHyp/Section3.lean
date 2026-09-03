import ABC3.Meta.Claim
import ABC3.Interface.CorrHyp.HyperbolicCurve
import ABC3.Skeleton.CorrHyp.Section1
import ABC3.Skeleton.CorrHyp.Section2

/-!
# [CorrHyp] §3 The Non-arithmetic Case —— 必要 3 件の statement

原典: S. Mochizuki, *Correspondences on Hyperbolic Curves* [CorrHyp]、物理 p.8。
**200dpi 目視確認 2026-09-04**(`ResearchPaper/papers.json`)。

★**`sorry` は「正しい状態」である**——`Skeleton/` は statement 専用トラックだからである。

★`hyperbolic core`(`core` / `coreMap`)は `Interface/CorrHyp/HyperbolicCurve.lean`
に**構成そのものではなく対象と自然な射だけ**を posit してある——`Lemma 5.1` の
内容(非アルキメデス版への一般化を含む構成と普遍性)がまだ `sorry` であることの反映。 -/

namespace ABC3.Skeleton.CorrHyp

open ABC3.Meta ABC3.Interface.CorrHyp

variable (D : HyperbolicCurveData)

/-- [CorrHyp] **Definition 3.1**。

原文 (CorrHyp p.8):
> Suppose that X is not arithmetic. Then we shall refer to Y (respectively,
> Y) as the hyperbolic core of X (respectively, X ).

★`core`/`coreMap` は `Interface` に(非 arithmetic の仮定なしで)posit してある。
本 statement はその対象に「hyperbolic core」という名前を与える段であり、
非 arithmetic 性の仮定を明示するために `h` を引数として持つ(未使用)。 -/
def hyperbolicCore (X : D.Space) (_h : ¬ Arithmetic D (D.Gamma X)) : D.Space := D.core X

/-- [CorrHyp] **Proposition 3.2**。

原文 (CorrHyp p.8):
> We have Γ_Z ⊆ Comm(Γ_X). -/
theorem prop_3_2 (X Z : D.Space) (_hX : ¬ Arithmetic D (D.Gamma X))
    (c : Corr D X Z) (_hC : D.IsConnected c.C) :
    D.Sub (D.Gamma Z) (D.Comm (D.Gamma X)) := sorry

/-- [CorrHyp] **Theorem 3.3**。

原文 (CorrHyp p.8–p.9):
> Suppose that X is not arithmetic. Fix a pair of nonnegative integers
> (g′, r′) such that 2g′ − 2 + r′ > 0. Then there exist (up to isomorphism)
> only finitely many hyperbolic curves Z of type (g′, r′) that are isogenous
> to Z. -/
theorem thm_3_3 (X : D.Space) (_hX : ¬ Arithmetic D (D.Gamma X)) (gr : ℕ × ℕ) :
    FinitelyManyUpTo D (fun Z => D.type Z = gr ∧ IsIsogenous D X Z) := sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def hyperbolicCore.src : Source :=
  { paper := "CorrHyp", pdfPage := 8, item := "Definition 3.1",
    sectionId := "corrhyp-def-3-1" }

def prop_3_2.src : Source :=
  { paper := "CorrHyp", pdfPage := 8, item := "Proposition 3.2",
    sectionId := "corrhyp-prop-3-2" }

/-- ★原文の証明(p.8)は5行、`C` が `Z` 上 Galois という仮定へ帰着してから
`Γ_X ∩ (γ·Γ_X·γ⁻¹)` の有限指数性を直接計算する——外部引用は無い
(`hedge-index.mjs` の合図0件と一致)。 -/
def prop_3_2.needs : List ProofObligation :=
  [ .implicitStep "Cを Z 上 Galois に取り直せる(一般の場合はGalois閉包を取る)、という初手を『we may assume』の1行で済ませている" 8 ]

def thm_3_3.src : Source :=
  { paper := "CorrHyp", pdfPage := 8, item := "Theorem 3.3",
    sectionId := "corrhyp-thm-3-3" }

/-- ★原文の証明(p.9)の核: `Comm(Γ_X)` が有限生成(有限指数部分群 `Γ_X` が
有限生成だから)⟹ 次数を固定した有限被覆は同型を除いて有限個、という
被覆空間論の標準事実。`hedge-index.mjs` の合図0件と一致(この項目には
immediately/well-known 等は現れない——地の文で直接計算されている)。 -/
def thm_3_3.needs : List ProofObligation :=
  [ .folklore "有限生成群の、次数を固定した部分群(=有限被覆)は同型を除いて有限個、という被覆空間論の標準事実" 9,
    .derivation "Γ_Z ⊆ Comm(Γ_X) かつ Comm(Γ_X) が Γ_X を有限指数部分群に持つ(Prop 3.2 と、非arithmeticから Γ_X が Comm(Γ_X) の中で有限指数であること)から Comm(Γ_X) 自身の有限生成性を出す" 9 ]

end ABC3.Skeleton.CorrHyp
