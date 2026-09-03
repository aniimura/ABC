import ABC3.Meta.Claim
import ABC3.Interface.CorrHyp.HyperbolicCurve
import ABC3.Skeleton.CorrHyp.Section1
import ABC3.Skeleton.CorrHyp.Section2
import ABC3.Skeleton.CorrHyp.Section3

/-!
# [CorrHyp] §5 Isogenies of General Curves —— 必要 7 件の statement

原典: S. Mochizuki, *Correspondences on Hyperbolic Curves* [CorrHyp]、物理 p.11–p.16。
**200dpi 目視確認 2026-09-04**(`ResearchPaper/papers.json` p.11・p.12 のみ実測、
p.13・p.15 は未目視——`Lemma 5.4`-`5.6`・`Theorem 5.7` は pdftotext 由来)。

★**`sorry` は「正しい状態」である**——`Skeleton/` は statement 専用トラックだからである。 -/

namespace ABC3.Skeleton.CorrHyp

open ABC3.Meta ABC3.Interface.CorrHyp

variable (D : HyperbolicCurveData)

/-- [CorrHyp] **Lemma 5.1**。

原文 (CorrHyp p.11):
> Let k be algebraically closed of characteristic zero. Suppose that k is a
> subfield of C. Let X be a hyperbolic curve over k. Suppose that X_C is not
> arithmetic. Then the morphism X_C → Y_C … descends to some morphism X → Y
> over k. Moreover, X → Y has the universal property that any correspondence
> (C → X, C → Z) over k, fits uniquely into a commutative diagram … Finally,
> the morphism X → Y is independent (up to canonical isomorphism) of the
> embedding of k into C.

★`core`/`coreMap` は `Interface` に(この非 arithmetic の仮定なしで)posit して
ある。本 statement はその普遍性を証明として要求する。「埋め込みに依らない」の
節は、本 Skeleton が `k ↪ ℂ` を型として持たないため対象外とした(逸脱、記録)。 -/
theorem lemma_5_1 (X : D.Space) (_h : ¬ Arithmetic D (D.Gamma (D.Ext X))) :
    ∀ {C Z : D.Space} (γ : D.FEt C X) (δ : D.FEt C Z),
      ∃! ζ : D.FEt Z (D.core X), D.comp δ ζ = D.comp γ (D.coreMap X) := sorry

/-- [CorrHyp] **Definition 5.2**。

原文 (CorrHyp p.12):
> Suppose that we are in the situation of Lemma 5.1. Then we shall refer to
> the stack Y constructed in Lemma 5.1 as the hyperbolic core Y of X. -/
def hyperbolicCoreGeneral (X : D.Space) (_h : ¬ Arithmetic D (D.Gamma (D.Ext X))) :
    D.Space := D.core X

/-- [CorrHyp] **Theorem 5.3**。

原文 (CorrHyp p.12):
> Let k be algebraically closed of characteristic zero. Suppose that k is a
> subfield of C. Fix nonnegative integers g and r such that 2g − 2 + r ≥ 3.
> Then there exists an open dense substack U ⊆ (M_{g,r})_k … with the
> following property: If X is a hyperbolic curve over some extension
> algebraically closed field K of k corresponding to a point ∈ U(K), then the
> hyperbolic core of X is equal to X.

★`2g − 2 + r ≥ 3` は ℕ の減法(切り捨て)を避けて `2g + r ≥ 5`(ℤ上で同値)と書いた。 -/
theorem thm_5_3 (g r : ℕ) (_h : 2 * g + r ≥ 5) :
    ∃ U : Set D.Space, D.IsOpenDenseIn U g r ∧ ∀ X ∈ U, D.Iso (D.core X) X := sorry

/-- [CorrHyp] **Lemma 5.4**。

原文 (CorrHyp p.13):
> The expression in parentheses e_Y is bounded below by an absolute positive
> constant, independent of X, g, and r. -/
theorem lemma_5_4 : ∃ c : ℚ, 0 < c ∧ ∀ X : D.Space, c ≤ (D.stackType (D.core X)).e := sorry

/-- [CorrHyp] **Lemma 5.5**。

原文 (CorrHyp p.13):
> If g and r are fixed, then there is only a finite number of possibilities
> for d, g_Y, r_Y, and Σ_Y. -/
theorem lemma_5_5 (g r : ℕ) :
    Set.Finite {p : ℕ × StackType | ∃ X : D.Space, D.type X = (g, r) ∧
      p = (D.deg (D.coreMap X), D.stackType (D.core X))} := sorry

/-- [CorrHyp] **Lemma 5.6**。

原文 (CorrHyp p.13):
> The locus (inside (M_{g,r})_k) of nonarithmetic curves that are not equal
> to their own hyperbolic cores is constructible (in (M_{g,r})_k). -/
theorem lemma_5_6 (g r : ℕ) :
    D.IsConstructibleIn
      {X : D.Space | D.type X = (g, r) ∧ ¬ Arithmetic D (D.Gamma X) ∧ ¬ D.Iso (D.core X) X}
      g r := sorry

/-- [CorrHyp] **Theorem 5.7**。

原文 (CorrHyp p.15):
> For a "general" (in the same sense as in the statement of Theorem 5.3)
> hyperbolic curve X of type (g, r), the canonical morphism X → Y to the
> hyperbolic core of X may be described as follows: g_Y = 0 and
> (1) If (g, r) = (0, 4), then X → Y has degree 4, r_Y = 1, |Σ_Y| = 3, all the
>     i_σ are 2, and the ramification index at the point at infinity of Y is 1.
> (2) If (g, r) = (1, 1), then X → Y has degree 2, r_Y = 1, |Σ_Y| = 3, all the
>     i_σ are 2, and the ramification index at the point at infinity of Y is 2.
> (3) If (g, r) = (1, 2), then X → Y has degree 2, r_Y = 1, |Σ_Y| = 4, all the
>     i_σ are 2, and the ramification index at the point at infinity of Y is 1.
> (4) If (g, r) = (2, 0), then X → Y has degree 2, r_Y = 0, |Σ_Y| = 6, and all
>     the i_σ are 2.
> Finally, if (g, r) = (0, 3), then X is arithmetic, so the hyperbolic core
> is not defined.

★「ramification index at the point at infinity」は `Σ_Y` の中の特定の1点を
指すため、`Fintype.card` だけでは表せない——本 statement では `|Σ_Y|` と
「全ての `i_σ = 2`」の部分だけを型にし、無限遠点の特定は `Gap` 相当として
`.needs` に残す(逸脱、記録)。 -/
theorem thm_5_7 (g r deg rY sigmaCard : ℕ) (X : D.Space)
    (_hcase : (g, r, deg, rY, sigmaCard) = (0, 4, 4, 1, 3) ∨
              (g, r, deg, rY, sigmaCard) = (1, 1, 2, 1, 3) ∨
              (g, r, deg, rY, sigmaCard) = (1, 2, 2, 1, 4) ∨
              (g, r, deg, rY, sigmaCard) = (2, 0, 2, 0, 6))
    (_hX : D.type X = (g, r)) :
    (D.stackType (D.core X)).g = 0 ∧
    D.deg (D.coreMap X) = deg ∧
    (D.stackType (D.core X)).r = rY ∧
    Fintype.card (D.stackType (D.core X)).Sigma = sigmaCard ∧
    ∀ σ, (D.stackType (D.core X)).i σ = 2 := sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def lemma_5_1.src : Source :=
  { paper := "CorrHyp", pdfPage := 11, item := "Lemma 5.1",
    sectionId := "corrhyp-lemma-5-1" }

/-- ★`hedge-index.mjs --cite` の実測: `clearly×1`(`X'_C → Y_C` が Galois になる
有限 étale Galois 被覆 `X'_C → X_C` が存在すること)・`immediately×1`(普遍性は
`Lemma 4.1` の議論を ℂ から k へ降下させれば「immediately」出る)。 -/
def lemma_5_1.needs : List ProofObligation :=
  [ .folklore "有限次 Galois 閉包の存在(任意の有限エタール被覆はGalois被覆で優越される)" 11,
    .otherPaper "[CorrHyp]" "Lemma 4.1(C から k への降下の議論を流用)" 9,
    .implicitStep "自己同型の剛性(Aut_C(X'_C) = Aut_k(X')、étale coveringの剛性と同型の主張)を1行で済ませている" 11 ]

def hyperbolicCoreGeneral.src : Source :=
  { paper := "CorrHyp", pdfPage := 12, item := "Definition 5.2",
    sectionId := "corrhyp-def-5-2" }

def thm_5_3.src : Source :=
  { paper := "CorrHyp", pdfPage := 12, item := "Theorem 5.3",
    sectionId := "corrhyp-thm-5-3" }

/-- ★`hedge-index.mjs --cite` の実測: `well-known×1`(証明の核となる
Riemann-Hurwitz を使った計算そのものは「well-known」と明言されている)。
本項目が要求するものは `Lemma 5.4`-`5.6` への分解と、Riemann-Hurwitz の式。 -/
def thm_5_3.needs : List ProofObligation :=
  [ .folklore "Riemann-Hurwitz の公式(2g-2+r = d·e_Y)" 12,
    .otherPaper "[CorrHyp]" "Lemma 5.4(e_Y の下界)" 13,
    .otherPaper "[CorrHyp]" "Lemma 5.5(d, g_Y, r_Y, Σ_Y の有限性)" 13,
    .otherPaper "[CorrHyp]" "Lemma 5.6(構成可能性、稠密開集合の構成に使う)" 13 ]

def lemma_5_4.src : Source :=
  { paper := "CorrHyp", pdfPage := 13, item := "Lemma 5.4",
    sectionId := "corrhyp-lemma-5-4" }

def lemma_5_4.needs : List ProofObligation :=
  [ .derivation "gY・rY・|ΣY| の場合分けによる組み合わせ論的な評価(2g-2+r=0,-1,-2の3ケース、各ケースでi_σの制約から下界を出す)" 13 ]

def lemma_5_5.src : Source :=
  { paper := "CorrHyp", pdfPage := 13, item := "Lemma 5.5",
    sectionId := "corrhyp-lemma-5-5" }

/-- ★`hedge-index.mjs --cite` の実測: `immediately×1`(dが有界になれば
可能性が有限個というのが1行で済まされている)。 -/
def lemma_5_5.needs : List ProofObligation :=
  [ .otherPaper "[CorrHyp]" "Lemma 5.4(e_Y の下界、d の有界性の入力)" 13,
    .implicitStep "2g_Y-2+r_Y+(1/2)|Σ_Y| ≤ e_Y という不等式から g_Y, r_Y, |Σ_Y| の有限性を出す段" 13 ]

def lemma_5_6.src : Source :=
  { paper := "CorrHyp", pdfPage := 13, item := "Lemma 5.6",
    sectionId := "corrhyp-lemma-5-6" }

/-- ★`hedge-index.mjs --cite` の実測: `immediately×1`、引用なし
(「引用の無い合図がある」と警告される項目)。証明はモジュライスタックの
理論(有限エタール被覆の変形空間、有限性)を直接使う——本プロジェクトに
モジュライスタックの構成は無い。 -/
def lemma_5_6.needs : List ProofObligation :=
  [ .folklore "有限エタール被覆のモジュライスタック N・N′ の構成(stacksの理論)" 13,
    .implicitStep "N′ → (M_{g,r})_k の像が構成可能である、というChevalleyの定理相当を『hence is constructible』の1行で済ませている" 13 ]

def thm_5_7.src : Source :=
  { paper := "CorrHyp", pdfPage := 15, item := "Theorem 5.7",
    sectionId := "corrhyp-thm-5-7" }

/-- ★`hedge-index.mjs --cite` の実測: `well-known×1`、引用なし。証明(p.15-16)は
4つの例外型それぞれで具体的な被覆(置換で与えられる)を構成する組み合わせ論。
★本 `theorem` は無限遠点での分岐指数を型にしていない(ファイル冒頭の逸脱注記)。 -/
def thm_5_7.needs : List ProofObligation :=
  [ .implicitStep "4つの例外型それぞれで、置換の組(例: (12)(34), (13)(24), (14)(23), id)で被覆を明示的に構成し、慣性群の議論で線形不連結性(=正しい合成)を確認する——場合ごとに異なる組み合わせ論の計算" 15,
    .implicitStep "無限遠点での分岐指数の特定は本 statement の型に入っていない(逸脱として記録済み)" 15 ]

end ABC3.Skeleton.CorrHyp
