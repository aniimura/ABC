import ABC3.Meta.Claim
import ABC3.Interface.CorrHyp.HyperbolicCurve
import ABC3.Skeleton.CorrHyp.Section1
import ABC3.Skeleton.CorrHyp.Section2
import ABC3.Skeleton.CorrHyp.Section3

/-!
# [CorrHyp] §4 The Main Theorem —— 必要 2 件の statement

原典: S. Mochizuki, *Correspondences on Hyperbolic Curves* [CorrHyp]、物理 p.9–p.10。
**200dpi 目視確認 2026-09-04**(`ResearchPaper/papers.json`)。

★**`sorry` は「正しい状態」である**——`Skeleton/` は statement 専用トラックだからである。 -/

namespace ABC3.Skeleton.CorrHyp

open ABC3.Meta ABC3.Interface.CorrHyp

variable (D : HyperbolicCurveData)

/-- `Ext` の `Corr` への拡張(`Lemma 4.1` の記述に使う)。 -/
def extCorr {A B : D.Space} (c : Corr D A B) : Corr D (D.Ext A) (D.Ext B) :=
  ⟨D.Ext c.C, D.extFEt c.α, D.extFEt c.β⟩

/-- 共通語彙(numbered item ではない)なので、`Lemma 4.1` の `X_K` を指す。 -/
def extCorr.src : Source :=
  { paper := "CorrHyp", pdfPage := 9, item := "Lemma 4.1 直前の地の文(X_K の記法)",
    sectionId := "corrhyp-lemma-4-1" }

/-- [CorrHyp] **Lemma 4.1**。

原文 (CorrHyp p.9):
> Suppose that X_K is isogenous to some hyperbolic curve Z_K over K. Then Z_K
> is the result of base-extending some hyperbolic curve Z over k from k to K,
> and, moreover, any correspondence from X_K to Z_K descends to a
> correspondence from X to Z. -/
theorem lemma_4_1 (X ZK : D.Space) (c : Corr D (D.Ext X) ZK) :
    ∃ (Z : D.Space) (h : ZK = D.Ext Z) (c' : Corr D X Z), h ▸ extCorr D c' = c := sorry

/-- [CorrHyp] **Theorem 4.2**(原文の "first main result"、`Theorem A` の本体)。

原文 (CorrHyp p.10、pdftotextの実測ではこの箇所だけ ′(プライム)が脱落し、
− (U+2212)も素の "-" になる——実際に呼び出したときの抽出結果に合わせて書く):
> Let k be an algebraically closed field of characteristic zero. Let X be a
> hyperbolic curve over k. Let (g, r) be a pair of nonnegative integers satisfying 2g-2+r >
> 0. Then (up to isomorphism) there are only finitely many hyperbolic curves over k of type
> (g, r) that are isogenous to X. -/
theorem thm_4_2 (X : D.Space) (gr' : ℕ × ℕ) :
    FinitelyManyUpTo D (fun Y => D.type Y = gr' ∧ IsIsogenous D X Y) := sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def lemma_4_1.src : Source :=
  { paper := "CorrHyp", pdfPage := 9, item := "Lemma 4.1",
    sectionId := "corrhyp-lemma-4-1" }

def lemma_4_1.needs : List ProofObligation :=
  [ .folklore "有限型の代数的対象は「有限個の方程式」で定義できる、というスキーム論の一般論(降下理論、EGA IV 相当)" 10,
    .implicitStep "étale site の変形に対する剛性(有限エタール射は変形しない)を1行の主張として使っている" 10,
    .implicitStep "接空間の計算(H¹(Z_s, τ_{Z_s}) → H¹(C, τ_C) のpull-back写像がtraceにより単射)で標数0を使う——ここは実際に計算が書かれている(folkloreではない)" 10 ]

def thm_4_2.src : Source :=
  { paper := "CorrHyp", pdfPage := 10, item := "Theorem 4.2",
    sectionId := "corrhyp-thm-4-2" }

/-- ★`hedge-index.mjs --paper CorrHyp --item "Theorem 4.2"` の実測:
`easily×1`(有限生成部分体への還元)・`immediately×1`(arithmeticなら
isogenousな曲線もarithmetic)・`formally×1`(k=Cの場合への還元)。
証明そのものは `Theorem 2.6`(arithmetic の場合)と `Theorem 3.3`
(非arithmetic の場合)への場合分けで**すでに Skeleton にある**——
本項目はその組み合わせだけを畳んでいる。 -/
def thm_4_2.needs : List ProofObligation :=
  [ .implicitStep "任意有限個の曲線を含む有限生成部分体 k′ への還元が「it suffices to show」の1行で済まされている(easily)" 10,
    .implicitStep "k′ の代数閉包を C に埋め込めることから k=C の場合に帰着する段(formally)" 10,
    .derivation "X が arithmetic なら isogenous な曲線もarithmetic、という言明(『it follows easily from the definitions』ではなくimmediatelyに分類)——IsIsogenousとArithmeticの定義から出るはずだが未証明" 10,
    .otherPaper "[CorrHyp]" "Theorem 2.6(X arithmetic の場合)" 7,
    .otherPaper "[CorrHyp]" "Theorem 3.3(X 非arithmetic の場合)" 8 ]

end ABC3.Skeleton.CorrHyp
