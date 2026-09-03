import ABC3.Meta.Claim
import Mathlib.GroupTheory.Commutator.Basic
import Mathlib.Topology.Algebra.Group.Basic

/-!
# [LocProP] §0 —— `p`-derivate(`Δ<i>`)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、物理 p.15。
**260 dpi 目視確認 2026-09-04**(このページに逐語の食い違いは無い)。

原文 (LocProP p.15、`Definition 0.2` 直前の設定):
> Let Δ′ ⊆Δ be the unique normal subgroup of Δ with the following
> property: Δ →Δ/Δ′ is the maximal (topologically) Hausdorffabelian quotient of Δ which
> is annihilated by p. For i ≥0, let Δ<0> def
> = Δ; Δ<i+1> def
> = (Δ<i>)′.

原文 (LocProP p.15):
> Definition 0.2. We shall refer to any one of the Δ<i> as a p-derivate of Δ.

## ★設計 —— 明示的な生成子で `Δ′` を構成する

原文は `Δ′` を**普遍性**(「`Δ/Δ′` が…最大の…商である」)で特徴づける。
★★**本ファイルは同値な明示公式を definition として採用する** ——
「`H` の交換子集合と `p` 乗の集合が生成する正規部分群の位相閉包」。
これは位相群論の標準的事実(`Δ/closure⟨[Δ,Δ], Δ^p⟩` が最大 Hausdorff・可換・
`p` 消去の商であること)であり、**証明はしていない**(下の `.needs` に folklore として記録)。
★★★**この turn で実際に検証したのは「well-defined な降下列であること」の部分**
—— 正規性・閉性・降下性は機械証明済み(`sorry` 無し)。

## mathlib 在庫

- `Subgroup.normalClosure` —— 生成集合の正規閉包(共役不変集合なら通常の閉包と一致)
- `Subgroup.topologicalClosure` / `Subgroup.is_normal_topologicalClosure` /
  `Subgroup.isClosed_topologicalClosure` —— 位相閉包とその正規性・閉性
- `Subgroup.commutator_le_self` —— `⁅H, H⁆ ≤ H`
- `Topology/Algebra/Group/TopologicalAbelianization.lean` の
  `TopologicalAbelianization`(`G ⧸ (commutator G).topologicalClosure`)—— 可換化の部分だけなら既にある。
  ★`p` 消去まで込みの版は無い(2026-09-04 実測、`grep` 0 件)。
-/

namespace ABC3.Found.LocProP

variable {G : Type*} [Group G] [TopologicalSpace G] [IsTopologicalGroup G]

/-- **`derivedStep p H`** —— `H` の交換子集合と `p` 乗の集合が生成する正規部分群の位相閉包。
`H` が `G` で正規なら、この生成集合は既に共役不変なので
(`h[a,b]h⁻¹ = [hah⁻¹, hbh⁻¹]`、`h(g^p)h⁻¹ = (hgh⁻¹)^p`)、
`normalClosure` は通常の `closure` と一致する(明示証明はしていない、下の `.needs`)。 -/
def derivedStep (p : ℕ) (H : Subgroup G) : Subgroup G :=
  Subgroup.topologicalClosure
    (Subgroup.normalClosure ((⁅H, H⁆ : Subgroup G) ∪ {x : G | ∃ g ∈ H, x = g ^ p}))

instance derivedStep_normal (p : ℕ) (H : Subgroup G) : (derivedStep p H).Normal :=
  Subgroup.is_normal_topologicalClosure _

theorem derivedStep_isClosed (p : ℕ) (H : Subgroup G) :
    IsClosed ((derivedStep p H : Subgroup G) : Set G) :=
  Subgroup.isClosed_topologicalClosure _

/-- `H` が正規かつ閉なら `derivedStep p H ≤ H`
(生成子(交換子・`p` 乗)がすべて `H` に入り、`H` は閉なので閉包を取っても超えない)。 -/
theorem derivedStep_le (p : ℕ) (H : Subgroup G) [H.Normal] (hHc : IsClosed (H : Set G)) :
    derivedStep p H ≤ H := by
  apply Subgroup.topologicalClosure_minimal
  · apply Subgroup.normalClosure_le_normal
    rintro x (hx | ⟨g, hg, rfl⟩)
    · exact Subgroup.commutator_le_self H hx
    · exact H.pow_mem hg p
  · exact hHc

/-- **`Δ<i>`**(`pDerivate p i`)。原文どおり `Δ<0> = Δ`、`Δ<i+1> = (Δ<i>)′`。 -/
def pDerivate (p : ℕ) : ℕ → Subgroup G
  | 0 => ⊤
  | i + 1 => derivedStep p (pDerivate p i)

instance pDerivate_normal (p : ℕ) (i : ℕ) : (pDerivate (G := G) p i).Normal := by
  cases i <;> unfold pDerivate <;> infer_instance

theorem pDerivate_isClosed (p : ℕ) (i : ℕ) :
    IsClosed ((pDerivate (G := G) p i : Subgroup G) : Set G) := by
  cases i with
  | zero => simp [pDerivate]
  | succ n => exact derivedStep_isClosed p (pDerivate p n)

/-- ★★**降下列であること**(原文 "descending series … `… ⊆ Δ<i> ⊆ … ⊆ Δ`")。 -/
theorem pDerivate_antitone (p : ℕ) (i : ℕ) :
    pDerivate (G := G) p (i + 1) ≤ pDerivate p i :=
  derivedStep_le p (pDerivate p i) (pDerivate_isClosed p i)

/-- **[LocProP] Definition 0.2** —— `Δ` の "p-derivate" とは、`Δ<i>` のいずれかであること。

★**非退化性**: `i = 0` を取れば `H = ⊤` が常に witness になる
(`Δ` 自身は「0 番目の p-derivate」)。 -/
def IsPDerivate (p : ℕ) (H : Subgroup G) : Prop := ∃ i, H = pDerivate p i

theorem isPDerivate_nonvacuous (p : ℕ) : ∃ H : Subgroup G, IsPDerivate (G := G) p H :=
  ⟨⊤, 0, rfl⟩

def IsPDerivate.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 15, item := "Definition 0.2",
    sectionId := "locprop-def-0-2" }

/-- ★原文はこの定義の直前で `Δ′` を普遍性で特徴づける(証明していない、`folklore`)。 -/
def IsPDerivate.needs : List ABC3.Meta.ProofObligation :=
  [ .folklore
      "Δ/derivedStep(Δ) が「最大の Hausdorff・可換・p 消去の商」であることは古典的事実だが、
       本ファイルは derivedStep を明示公式として definition に採用し、
       この普遍性そのもの(および一意性)は証明していない。機械検証したのは
       正規性・閉性・降下列であることのみ(derivedStep_le / pDerivate_antitone)。" 15 ]

end ABC3.Found.LocProP
