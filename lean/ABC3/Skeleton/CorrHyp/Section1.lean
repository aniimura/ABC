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
> of finite, ´etale morphisms α : C →X, β : C →Y , where we assume that C is nonempty. -/
structure Corr (X Y : D.Space) where
  C : D.Space
  α : D.FEt C X
  β : D.FEt C Y

/-- [CorrHyp] **Definition 1.2**。

原文 (CorrHyp p.4):
> We shall refer to a correspondence (α : C →X, β : C →Y ) from X to
> Y as trivial if there exists a finite ´etale morphism γ : Y →X such that α = γ ◦β. -/
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
> the composite of these two correspondences the correspondence given by the following pair
> of morphisms: the first morphism C1 ×Y C2 →X is given by composing the projection to
> C1 with α1; the second morphism C1 ×Y C2 →Z is given by composing the projection to
> C2 with β2. -/
def Corr.comp' {X Y Z : D.Space} (c1 : Corr D X Y) (c2 : Corr D Y Z) : Corr D X Z :=
  ⟨D.pullback c1.β c2.α, D.comp (D.pbFst c1.β c2.α) c1.α, D.comp (D.pbSnd c1.β c2.α) c2.β⟩

/-- [CorrHyp] **Definition 1.5**。

原文 (CorrHyp p.4):
> We shall call two hyperbolic curves X and Y over k isogenous if there
> exists a corresondence from X to Y . -/
def IsIsogenous (X Y : D.Space) : Prop := Nonempty (Corr D X Y)

/-- 原文 p.4 末尾: 「it follows immediately that the relation of isogeny is an
equivalence relation」。反射性は `idFEt`(`Found/CorrHyp/FuchsianGroup.lean` の
`isFiniteIndexIn_refl` が具体化)、対称性は `Corr.transpose`、推移性は
`Corr.comp'` から、sorry 無しで出る。★numbered item ではないが、
`IsIsogenous.needs`(旧版)が「反射性に要る」と記していた欠落はこれで埋まった。 -/
theorem isIsogenous_refl (X : D.Space) : IsIsogenous D X X := ⟨X, D.idFEt X, D.idFEt X⟩

theorem isIsogenous_symm {X Y : D.Space} (h : IsIsogenous D X Y) : IsIsogenous D Y X :=
  h.elim fun c => ⟨Corr.transpose D c⟩

theorem isIsogenous_trans {X Y Z : D.Space} (h1 : IsIsogenous D X Y) (h2 : IsIsogenous D Y Z) :
    IsIsogenous D X Z :=
  h1.elim fun c1 => h2.elim fun c2 => ⟨Corr.comp' D c1 c2⟩

def isIsogenous_refl.src : Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.5 直後(isogenyは同値関係、反射性)",
    sectionId := "corrhyp-def-1-5" }

def isIsogenous_refl.needs : List ProofObligation := []

def isIsogenous_symm.src : Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.5 直後(isogenyは同値関係、対称性)",
    sectionId := "corrhyp-def-1-5" }

def isIsogenous_symm.needs : List ProofObligation := []

def isIsogenous_trans.src : Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.5 直後(isogenyは同値関係、推移性)",
    sectionId := "corrhyp-def-1-5" }

def isIsogenous_trans.needs : List ProofObligation := []

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
の `immediately×1`)。★★**2026-09-04、埋まった**——`idFEt` を `Interface` に
足し、`isIsogenous_refl`/`_symm`/`_trans` を sorry 無しで証明した(このファイル
上部)。 -/
def IsIsogenous.needs : List ProofObligation := []

/-- 共通語彙(numbered item ではない)なので、最初に要る `Definition 1.5` の
箇所を指す。 -/
def FinitelyManyUpTo.src : Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.5(直後、isogenousの語彙を使う共通の言い方)",
    sectionId := "corrhyp-def-1-5" }

end ABC3.Skeleton.CorrHyp
