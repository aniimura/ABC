import ABC3.Meta.Claim
import ABC3.Interface.CorrHyp.HyperbolicCurve
import ABC3.Skeleton.CorrHyp.Section1

/-!
# [CorrHyp] §6 Interpretation of a Theorem of Royden —— 必要 1 件の statement

原典: S. Mochizuki, *Correspondences on Hyperbolic Curves* [CorrHyp]、物理 p.17。
**200dpi 目視確認 2026-09-04**(`ResearchPaper/papers.json`)。

★**`sorry` は「正しい状態」である**——`Skeleton/` は statement 専用トラックだからである。

これで `ResearchPaper/corrhyp-goal.md` の §1–§6、合計 24 件すべてを固定した。 -/

namespace ABC3.Skeleton.CorrHyp

open ABC3.Meta ABC3.Interface.CorrHyp

variable (D : HyperbolicCurveData)

/-- 原文 §6 冒頭が導入する、単一対象上の correspondence
(`Definition 1.1` の 2 対象版とは trivial の意味が違う)。

原文 (CorrHyp p.17):
> Let us refer to as a correspondence on Mg,r an (ordered) pair of finite ´etale morphisms
> α : E →Mg,r, β : E →Mg,r, where E is nonempty. We shall call a correspondence (α, β)
> on Mg,r trivial if α = β. -/
structure CorrOn (M : D.Space) where
  E : D.Space
  α : D.FEt E M
  β : D.FEt E M

def CorrOn.IsTrivial {M : D.Space} (c : CorrOn D M) : Prop := c.α = c.β

def CorrOn.IsTrivial.src : Source :=
  { paper := "CorrHyp", pdfPage := 17, item := "§6 冒頭(correspondence on M_{g,r} の trivial)",
    sectionId := "corrhyp-sec6-corron" }

/-- [CorrHyp] **Theorem 6.1**。

原文 (CorrHyp p.17):
> Suppose that 2g −2 + r ≥3. Then Mg,r is generically a scheme, and
> moreover, does not admit any nontrivial automorphisms or correspondences.

★`2g − 2 + r ≥ 3` は `Theorem 5.3` と同じく `2g + r ≥ 5`(ℤ上同値、ℕ減法回避)
で書いた。「非自明な自己同型を持たない」は `∀ a, a = idAut M`、
「非自明な correspondence を持たない」は `∀ c : CorrOn M, c.IsTrivial` として型にした。 -/
theorem thm_6_1 (g r : ℕ) (_h : 2 * g + r ≥ 5) :
    D.IsGenericallyScheme (D.ModuliStack g r) ∧
    (∀ a : D.Aut (D.ModuliStack g r), a = D.idAut (D.ModuliStack g r)) ∧
    (∀ c : CorrOn D (D.ModuliStack g r), c.IsTrivial) := sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def CorrOn.src : Source :=
  { paper := "CorrHyp", pdfPage := 17, item := "§6 冒頭(correspondence on M_{g,r})",
    sectionId := "corrhyp-sec6-corron" }

def thm_6_1.src : Source :=
  { paper := "CorrHyp", pdfPage := 17, item := "Theorem 6.1",
    sectionId := "corrhyp-thm-6-1" }

/-- ★`hedge-index.mjs --cite` の実測: `well-known×1`、引用なし
(「it is a matter of well-known general nonsense」)。証明の核は
Royden の定理([Gard] §9.2 p.169 Theorem 2、Teichmüller 空間 `𝒯` の
正則自己同型群と `Γ = π₁(M_{g,r})` の同型)であり、本プロジェクトに
Teichmüller 空間の構成は無い。 -/
def thm_6_1.needs : List ProofObligation :=
  [ .folklore "[Gard], §9.2, p. 169, Theorem 2(Roydenの定理——Aut(𝒯) ≅ Γ、外部文献・未登記)" 17,
    .folklore "『well-known general nonsense』(非自明な correspondence の存在 ⟹ Aut(𝒯)∖Γ に、Γ∩(γΓγ⁻¹) が両側で有限指数であるような元 γ が存在する、という一般論。[Marg] p.337 の議論と同型、Theorem 2.5 の証明で使ったのと同じ一般論)" 17 ]

end ABC3.Skeleton.CorrHyp
