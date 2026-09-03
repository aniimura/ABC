import ABC3.Meta.Claim
import ABC3.Interface.CorrHyp.HyperbolicCurve
import ABC3.Skeleton.CorrHyp.Section1

/-!
# [CorrHyp] §2 Review of Results of Margulis and Takeuchi —— 必要 6 件の statement

原典: S. Mochizuki, *Correspondences on Hyperbolic Curves* [CorrHyp]、物理 p.5–p.7。
**200dpi 目視確認 2026-09-04**(`ResearchPaper/papers.json`)。

★**`sorry` は「正しい状態」である**——`Skeleton/` は statement 専用トラックだからである。

★**`Definition 2.2`(Margulis-arithmetic)・`Definition 2.3`(Shimura-arithmetic)は
`Interface/CorrHyp/HyperbolicCurve.lean` の抽象化(ファイル冒頭の docstring)を見よ**。
本節ではその述語を単に名指すだけである。 -/

namespace ABC3.Skeleton.CorrHyp

open ABC3.Meta ABC3.Interface.CorrHyp

variable (D : HyperbolicCurveData)

/-- [CorrHyp] **Definition 2.1**。

原文 (CorrHyp p.5):
> We shall say that X, 𝒳, or Γ has infinitely many correspondences if Γ
> has infinite index in Comm(Γ). -/
def InfinitelyManyCorr (Γ : D.Fuchsian) : Prop := ¬ D.FiniteIndexIn Γ (D.Comm Γ)

/-- [CorrHyp] **Definition 2.2**(Margulis-arithmetic)。

原文 (CorrHyp p.5):
> We shall call X, 𝒳, or Γ Margulis-arithmetic if there exists a connected
> non-commutative almost Q-simple algebraic group G over Q, together with a
> surjection τ : G_R → (PSL₂)_R of algebraic groups over R such that the Lie
> group (Ker τ)(R) is compact, and the subgroups τ(G(Z)) and Γ (of PSL₂(R)⁰)
> are commensurable.

★`Interface` の抽象化(`HyperbolicCurveData.MargulisArithmetic`)をそのまま名指す。
中身(代数群 `G` の構成)は `Proposition 2.4` の証明の中だけで使われる——
本節点自身はそれを要求しない。 -/
def MargulisArithmetic (Γ : D.Fuchsian) : Prop := D.MargulisArithmetic Γ

/-- [CorrHyp] **Definition 2.3**(Shimura-arithmetic)。

原文 (CorrHyp p.5–p.6):
> We shall call X, 𝒳, or Γ Shimura-arithmetic if the following data exist:
> (1) a totally real algebraic number field F; (2) a quaternion algebra A
> over F which is trivial at one of the infinite places of F and nontrivial
> at all the other infinite places; (3) a trivialization of A at the infinite
> place of F at which A is trivial …; (4) an order O_A ⊆ A such that the
> intersection of O_A ⊆ A ⊆ M₂(R) with SL₂(R) ⊆ M₂(R) has image in PSL₂(R)⁰
> commensurable with Γ. -/
def ShimuraArithmetic (Γ : D.Fuchsian) : Prop := D.ShimuraArithmetic Γ

/-- 原文 §2、`Proposition 2.4` の直後の段落: arithmetic = Margulis-arithmetic
∨ Shimura-arithmetic。

原文 (CorrHyp p.6):
> In the future, we shall refer to X, 𝒳, or Γ as arithmetic if it is either
> Margulis-arithmetic or Shimura-arithmetic (since we now know that these two
> notions of arithmeticity are equivalent).

★numbered item ではないが、`Theorem 2.5` 以降の全定理が使う語彙なのでここに置く。 -/
def Arithmetic (Γ : D.Fuchsian) : Prop := D.MargulisArithmetic Γ ∨ D.ShimuraArithmetic Γ

/-- [CorrHyp] **Proposition 2.4**。

原文 (CorrHyp p.6):
> The Riemann surface 𝒳 is Margulis-arithmetic if and only if it is
> Shimura-arithmetic. -/
theorem prop_2_4 (Γ : D.Fuchsian) : D.MargulisArithmetic Γ ↔ D.ShimuraArithmetic Γ := sorry

/-- [CorrHyp] **Theorem 2.5**。

原文 (CorrHyp p.7):
> ([Marg], p. 337, Theorem 27; p. 60, Lemma 3.1.1, (v)) The hyperbolic
> Riemann surface 𝒳 is arithmetic if and only if it has infinitely many
> correspondences (in the sense of Definition 2.1). -/
theorem thm_2_5 (Γ : D.Fuchsian) : Arithmetic D Γ ↔ InfinitelyManyCorr D Γ := sorry

/-- [CorrHyp] **Theorem 2.6**。

原文 (CorrHyp p.7):
> ([Take], Theorem 2.1) There are only finitely many arithmetic X over C of
> a given type (g, r). -/
theorem thm_2_6 (gr : ℕ × ℕ) :
    FinitelyManyUpTo D (fun X => D.type X = gr ∧ Arithmetic D (D.Gamma X)) := sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def InfinitelyManyCorr.src : Source :=
  { paper := "CorrHyp", pdfPage := 5, item := "Definition 2.1",
    sectionId := "corrhyp-def-2-1" }

def MargulisArithmetic.src : Source :=
  { paper := "CorrHyp", pdfPage := 5, item := "Definition 2.2",
    sectionId := "corrhyp-def-2-2" }

def ShimuraArithmetic.src : Source :=
  { paper := "CorrHyp", pdfPage := 5, item := "Definition 2.3",
    sectionId := "corrhyp-def-2-3" }

def prop_2_4.src : Source :=
  { paper := "CorrHyp", pdfPage := 6, item := "Proposition 2.4",
    sectionId := "corrhyp-prop-2-4" }

/-- ★`node tools/hedge-index.mjs --paper CorrHyp --cite` の実測:
`similarly×1`(⇐ 方向の証明の後半、⇒ 方向と同型の議論の再演)・
`easily×1`(「that Shimura-arithmeticity implies Margulis-arithmeticity is
clear」)。証明本体(p.6–p.7)は代数群の almost-simple 因子分解・Galois
コホモロジー・Brauer 群という深い理論を使う——本項目の `sorry` はそこを畳む。 -/
def prop_2_4.needs : List ProofObligation :=
  [ .folklore "⇒ 方向: 「that Shimura-arithmeticity implies Margulis-arithmeticity is clear」(easily に分類されないが同種の省略)" 6,
    .implicitStep
      "⇐ 方向の核: G の almost Q-simple factor 分解(Gal(Q̄/Q) の推移作用)→ Weil restriction → Lie(H_F)_C ≅ sl₂(C) → H¹(F, PGL₂) の非可換 Galois コホモロジー → 四元数環 A の構成、という代数群の理論を丸ごと使う。mathlib に algebraic group の almost-simple 分解・Weil restriction は 2026-09-04 時点で未実測" 6,
    .implicitStep "後半(F が totally real で、量的高々2つの実場所を除いて A が非自明)は「if τ were trivial on this factor」型の議論を similarly で2度畳んでいる(p.7)" 7 ]

def thm_2_5.src : Source :=
  { paper := "CorrHyp", pdfPage := 7, item := "Theorem 2.5",
    sectionId := "corrhyp-thm-2-5" }

def thm_2_5.needs : List ProofObligation :=
  [ .otherPaper "[Marg]" "p. 337, Theorem 27; p. 60, Lemma 3.1.1, (v)" 337 ]

def thm_2_6.src : Source :=
  { paper := "CorrHyp", pdfPage := 7, item := "Theorem 2.6",
    sectionId := "corrhyp-thm-2-6" }

def thm_2_6.needs : List ProofObligation :=
  [ .otherPaper "[Take]" "Theorem 2.1" 1 ]

end ABC3.Skeleton.CorrHyp
