import ABC3.Meta.Claim
import ABC3.Interface.CorrHyp.HyperbolicCurve

/-!
# [CorrHyp] §1 Basic Definitions —— 必要 5 件の statement(`Skeleton`)

原典: S. Mochizuki, *Correspondences on Hyperbolic Curves* [CorrHyp]、物理 p.3–p.4。
**200dpi 目視確認 2026-09-04**(`ResearchPaper/papers.json`)。

★**`sorry` は「正しい状態」である**——`Skeleton/` は statement 専用トラックだからである
(本ファイルには `sorry` は無い——§1 はすべて定義であり、証明を要求しないため)。

`ResearchPaper/corrhyp-goal.md` §1 の 5 件をここで固定する。 -/

namespace ABC3.Skeleton.CorrHyp

open ABC3.Meta ABC3.Interface.CorrHyp

variable (D : HyperbolicCurveData)

/-- [CorrHyp] **Definition 1.1**。

原文 (CorrHyp p.3):
> We shall refer to as a correspondence from X to Y any (ordered) pair
> of finite, étale morphisms α : C → X, β : C → Y, where we assume that C is
> nonempty. -/
structure Corr (X Y : D.Space) where
  C : D.Space
  α : D.FEt C X
  β : D.FEt C Y

/-- [CorrHyp] **Definition 1.2**。

原文 (CorrHyp p.4):
> We shall refer to a correspondence (α : C → X, β : C → Y) from X to Y as
> trivial if there exists a finite étale morphism γ : Y → X such that
> α = γ ◦ β. -/
def Corr.IsTrivial {X Y : D.Space} (c : Corr D X Y) : Prop :=
  ∃ γ : D.FEt Y X, c.α = D.comp c.β γ

/-- [CorrHyp] **Definition 1.3**。

原文 (CorrHyp p.4):
> Given a correspondence (α, β) from X to Y, we shall refer to as the
> transpose correspondence to (α, β) the correspondence (from Y to X) given
> by the pair (β, α). -/
def Corr.transpose {X Y : D.Space} (c : Corr D X Y) : Corr D Y X :=
  ⟨c.C, c.β, c.α⟩

/-- [CorrHyp] **Definition 1.4**。

原文 (CorrHyp p.4):
> Then we shall refer to as the composite of these two correspondences the
> correspondence given by the following pair of morphisms: the first
> morphism C₁ ×_Y C₂ → X is given by composing the projection to C₁ with α₁;
> the second morphism C₁ ×_Y C₂ → Z is given by composing the projection to
> C₂ with β₂. -/
def Corr.comp' {X Y Z : D.Space} (c1 : Corr D X Y) (c2 : Corr D Y Z) : Corr D X Z :=
  ⟨D.pullback c1.β c2.α, D.comp (D.pbFst c1.β c2.α) c1.α, D.comp (D.pbSnd c1.β c2.α) c2.β⟩

/-- [CorrHyp] **Definition 1.5**。

原文 (CorrHyp p.4):
> We shall call two hyperbolic curves X and Y over k isogenous if there
> exists a correspondence from X to Y. -/
def IsIsogenous (X Y : D.Space) : Prop := Nonempty (Corr D X Y)

/-- 「(同型を除いて)有限個」を `Setoid`/`Quotient` を使わずに言い表す言い方。
`P` を満たす任意の `Y` が、有限集合 `S` の中のどれかと同型、という形
(`Theorem 2.6` / `Theorem 3.3` / `Theorem 4.2` で使う共通の語彙、numbered item ではない)。 -/
def FinitelyManyUpTo (P : D.Space → Prop) : Prop :=
  ∃ S : Finset D.Space, ∀ Y, P Y → ∃ Z ∈ S, D.Iso Y Z

/-! ## ★出典の紐付け(`.src`) -/

def Corr.src : Source :=
  { paper := "CorrHyp", pdfPage := 3, item := "Definition 1.1",
    sectionId := "corrhyp-def-1-1" }

def Corr.IsTrivial.src : Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.2",
    sectionId := "corrhyp-def-1-2" }

def Corr.transpose.src : Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.3",
    sectionId := "corrhyp-def-1-3" }

def Corr.comp'.src : Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.4",
    sectionId := "corrhyp-def-1-4" }

def IsIsogenous.src : Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.5",
    sectionId := "corrhyp-def-1-5" }

/-- ★原文 p.4 の一文「it follows immediately that the relation of isogeny is
an equivalence relation」を数えた(`node tools/hedge-index.mjs --paper CorrHyp`
の `immediately×1`)。**未着手**:反射性・対称性・推移性のうち対称性
(`Corr.transpose`)と推移性(`Corr.comp'`)は無条件に構成済み(`sorry` 無し)だが、
反射性(`X` から `X` への自明な correspondence、`C := X`・`α = β := id`)は
`FEt` に恒等射を posit していないので Interface の追加が要る。 -/
def IsIsogenous.needs : List ProofObligation :=
  [ .implicitStep
      "「it follows immediately」の中身は反射性・対称性・推移性の3点。対称性はDefinition 1.3(transpose)、推移性はDefinition 1.4(comp')からsorry無しで出るが、反射性はFEtの恒等射(id : FEt X X)をInterfaceに足す必要があり、まだ足していない" 4 ]

end ABC3.Skeleton.CorrHyp
