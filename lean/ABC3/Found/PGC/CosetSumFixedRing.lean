import ABC3.Found.PGC.TameQuotientTower
import ABC3.Found.PGC.HerbrandFirstJump
import ABC3.Found.PGC.FirstJumpLedger

/-!
# [pGC] ★★★★10 仮説の配管は**在庫でそろった** —— `u₁ < i_ϖ(σ̄)` が仮説ゼロで出た

## ★★配られた持ち場は「配管」だった。★実際は**新しい補題を 1 本も書かずに**そろった。

前波で「残るのは `coset_sum_truncENat` の 10 仮説を塔で揃える配管だけ」と書いた。
★★**その配管は木の `TameQuotientTower.lean:83-98` が既に一覧表にしていた**
（★本連鎖は本波で初めて開いた。同 `:83` の断定「10 個の仮定は `C := ↥(fixedRing B H)` に
対して**そのまま在庫でそろった**」は★**真だった**）。

| 仮定 | 供給（木の宣言） |
|---|---|
| `hcomp` | `algebraMap_smul_fixedRing`（★`rfl`） |
| `hHtriv` | `smul_fixedRing_eq_self` |
| `hinj` | `fixedRing_injective` |
| `hfixC` | `exists_algebraMap_fixedRing` |
| `hres` | `exists_sub_mem_fixedRing hresA` |
| `hAC` | `exists_algebraMap_fixedRing_eq` |
| `hϖ` | ★**`rfl`** |
| `hadj` | 仮説そのまま |
| `hfix` | `fixedRing_mem_adjoin_uniformizer hresA hAne hϖ` |
| `[MulSemiringAction G C]` / `[IsDiscreteValuationRing C]` | `FixedRingTower.lean:386 fixedRing_isDiscreteValuationRing`（instance） |

⇒ §1 `coset_sum_truncENat_fixedRing` は **`coset_sum_truncENat` への 1 回の代入（項 1 本）**である。

## ★★★成果 —— `FirstJumpLedger.lean:105` が「形式化していない」と書いた 1 本が出た

§2 `lt_ramIndex_quotient_fixedRing`:

★**`u` が第 1 跳びの下（`∀ τ ≠ 1, u < i_{π'}(τ)`）で `σ ∉ H` なら `u < i_ϖ(σ)`。**

★前波の `HerbrandFirstJump.lt_ramIndex_quotient_of_coset_sum` が仮説で受けていた
`hcoset` が★**落ちた**。残っているのは塔の標準的な仮説
（`hπ'` / `huni` / `hresA` / `hadj` / `hAne` / `hϖ`）だけである。

§3 `sub_one_mul_first_jump_le_of_quotient`:

★★`i_ϖ(σ) = i` と `(p−1)·i ≤ p·e_K` から ★**`(p−1)·u₁ ≤ p·e_K`**。
これは `FirstJumpLedger.lean:100-108`（★本連鎖が以前に書いたファイル）が
★「**点 1 は閉じなかった**。足りないのは `u₁ ≤ i(σ̄)` で、これは Herbrand が要る」
と書いた★**その点 1 そのもの**である。
★`ℕ` の糊は木の `sub_one_mul_first_jump_le`（同 `:421`）をそのまま使った。

## ★★それでも `axDecay p k` は出ない（★点 1 が閉じても、である）

`FirstJumpLedger.lean:20-27` の測定 2 が
★**「鎖の台帳が閉じるのは `s_m = 1`（最後の跳びが最大）のときだけ」**
（`chain_ledger_forces_max` / `dvd_of_chain_ledger` / `not_chain_ledger_of_not_dvd`）
を定理にしている。⇒ ★**鎖の道は点 1 を閉じても閉じない。**
★★これは `CyclicLayerDescent` の**鎖**の勘定であって、
`NormalizedTraceDescent.lean:379 FirstJumpRoute` の勘定とは**別**である。
★後者に効くかは★**本波でも測っていない**（前波と同じ。正直に繰り返す）。

## ★穴の現状（本波で増減なし）

①不分岐（`TotallyRamifiedLayer.lean:342`）/ ②′（`ExitLossDepth.not_forall_prod_Icc_le_of_le`）/
③幾何減衰（`JumpGeometricDecay.no_uniform_geometric_of_harith`）/
⑤出口の損失は `axDecay p 1` が限界（`ExitLossDepth.exit_does_not_realize_cEx`）。④は消えたまま。
★本波は**在庫の穴を 2 つ目に埋めた**（前波が Herbrand の 2 本、本波が配管）。

## ★①不分岐側との関係

★本波も独立である。`lt_ramIndex_quotient_fixedRing` は**全分岐を仮定していない**
（`hvalK` も `IsTotallyRamified` も出てこない）。★`huni`（`𝔪_B = (π')`）は
DVR なら常に取れるので分岐の仮定ではない。

## ★在庫の測定（★今回は最初から 2 か所で測った。#330）

```
grep -n "^theorem \|^instance " lean/ABC3/Found/PGC/FixedRingTower.lean
  → :386 instance fixedRing_isDiscreteValuationRing（★前波で「一番高い」と見込んだものは
     ★instance として既に在った）
grep -rn "IsDiscreteValuationRing.*fixedRing" lean/ABC3/Found/PGC/ --include=*.lean
  → FixedRingTower.lean:388 / HasseArfInduction.lean:90（供給表）/ TameQuotientTower.lean:97（供給表）
grep -n "coset_sum_truncENat" .cache/decl-index.txt   → 木にのみ（mathlib 索引には無い）
```
★**前波の見込み「`IsDiscreteValuationRing C` が一番高い」は外れ**（instance で無料だった）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. §1・§2 の `.src` は `HerbrandComposition.lean` と同じ Yoshida08 Lemma 6.10、
   §3 は pGC Cor 3.1（`FirstJumpLedger.lean` と同じ）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★本ファイルは新しい数学を 1 つも書いていない。既存の 2 つの塊
   （`TameQuotientTower` の供給表と前波の `HerbrandFirstJump`）を**繋いだだけ**である。
4. ★`hAne`（`∃ a ∈ 𝔪_A, algebraMap A B a ≠ 0`）は `FixedRingBaseAlgebra.lean:89` が
   「★原典には無い」と自認している逸脱を**そのまま引き継いでいる**。
-/

namespace ABC3.Found.PGC

namespace CosetSumFixedRing

open IsLocalRing IsDiscreteValuationRing

/-! ## §1 配管 —— 10 仮説を `C := fixedRing B H` で埋める（項 1 本） -/

section Plumbing

/-- ★★★**10 仮説の配管** —— `HerbrandComposition.lean:532 coset_sum_truncENat` に
`C := ↥(fixedRing B H)` を代入した形。

★**新しい補題は 1 本も要らなかった**（`TameQuotientTower.lean:83` の断定は真だった）。
★引数は `TameQuotientTower.lean:369 herbrandPhiGroup_comp_fixedRing` と**同じ並び**である。 -/
theorem coset_sum_truncENat_fixedRing
    {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype ↥H]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ) (n : ℝ) (σ : G) :
    ∑ τ : H, truncENat (ramIndex π' (σ * (τ : G))) (n + 1)
      = (Nat.card H : ℝ) * truncENat (ramIndex ϖ σ) (herbrandPhi π' H n + 1) :=
  coset_sum_truncENat (A := A) (fun ρ c => algebraMap_smul_fixedRing ρ c)
    smul_fixedRing_eq_self hπ' fixedRing_injective exists_algebraMap_fixedRing
    (exists_sub_mem_fixedRing hresA) exists_algebraMap_fixedRing_eq rfl hadj
    (fixedRing_mem_adjoin_uniformizer hresA hAne hϖ) n σ

end Plumbing

/-! ## §2 ★★★`hcoset` が落ちた —— `u₁ < i_ϖ(σ)` -/

section FirstJump

/-- ★★★★**`FirstJumpLedger.lean:105` が「形式化していない」と書いた 1 本**。

`u` が第 1 跳びの下（`∀ τ ≠ 1, u < i_{π'}(τ)`）で `σ ∉ H` なら `u < i_ϖ(σ)`。
★前波の `HerbrandFirstJump.lt_ramIndex_quotient_of_coset_sum` の `hcoset` が落ちた。
★中身は「`φ_H` が第 1 跳びの下で恒等」＋「剰余類和」の 2 本だけである。 -/
theorem lt_ramIndex_quotient_fixedRing
    {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype ↥H]
    {π' : B} (hπ' : Irreducible π') (huni : maximalIdeal B = Ideal.span {π'})
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ) {u : ℕ} {σ : G}
    (hfirst : ∀ τ : G, τ ≠ 1 → (u : ℕ∞) < ramIndex π' τ) (hσ : σ ∉ H) :
    (u : ℕ∞) < ramIndex ϖ σ := by
  have hfirstH : ∀ τ : ↥H, (u : ℕ∞) < ramIndex π' ((τ : G)) := by
    intro τ
    by_cases h1 : (τ : G) = 1
    · rw [h1, ramIndex_one]
      exact lt_of_le_of_ne le_top (by simp)
    · exact hfirst _ h1
  have hpos : ∀ τ : ↥H, 0 < ramIndex π' ((τ : G)) := fun τ => pos_ramIndex huni _
  have hphi : herbrandPhi π' H (u : ℝ) = (u : ℝ) :=
    HerbrandFirstJump.herbrandPhi_eq_self_of_first_jump H hpos hfirstH
  have hsum := coset_sum_truncENat_fixedRing (A := A) hπ' hresA hadj hAne hϖ (u : ℝ) σ
  rw [hphi] at hsum
  exact HerbrandFirstJump.lt_of_coset_sum
    (HerbrandFirstJump.forall_lt_ramIndex_of_first_jump hfirst hσ) hsum

end FirstJump

/-! ## §3 ★★`FirstJumpLedger` の「点 1」が閉じた -/

section Payoff

/-- ★★**`WildJumpChain` の点 1 が閉じた** —— `(p−1)·u₁ ≤ p·e_K`。

`FirstJumpLedger.lean:100-108` は「点 1 は閉じなかった。足りないのは `u₁ ≤ i(σ̄)` で、
これは Herbrand が要る」と書いていた。★`ℕ` の糊は同 `:421 sub_one_mul_first_jump_le`。

★★ただし同ファイルの測定 2（`chain_ledger_forces_max`）により、
★**点 1 が閉じても鎖の台帳は閉じない**。 -/
theorem sub_one_mul_first_jump_le_of_quotient
    {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype ↥H]
    {π' : B} (hπ' : Irreducible π') (huni : maximalIdeal B = Ideal.span {π'})
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ) {u : ℕ} {σ : G}
    (hfirst : ∀ τ : G, τ ≠ 1 → (u : ℕ∞) < ramIndex π' τ) (hσ : σ ∉ H)
    {i eK p : ℕ} (hi : ramIndex ϖ σ = (i : ℕ∞)) (hbound : (p - 1) * i ≤ p * eK) :
    (p - 1) * u ≤ p * eK := by
  have hlt : (u : ℕ∞) < (i : ℕ∞) := hi ▸ lt_ramIndex_quotient_fixedRing
    (A := A) hπ' huni hresA hadj hAne hϖ hfirst hσ
  have hui : u ≤ i := le_of_lt (by exact_mod_cast hlt)
  exact sub_one_mul_first_jump_le hui hbound

end Payoff

/-! ## §4 `.src` と 使っている公理の一覧 -/

section Src

def coset_sum_truncENat_fixedRing.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Lemma 6.10", sectionId := "lemma-6-10" }

def lt_ramIndex_quotient_fixedRing.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Lemma 6.10", sectionId := "lemma-6-10" }

def sub_one_mul_first_jump_le_of_quotient.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms coset_sum_truncENat_fixedRing
#print axioms lt_ramIndex_quotient_fixedRing
#print axioms sub_one_mul_first_jump_le_of_quotient

end CosetSumFixedRing

end ABC3.Found.PGC
