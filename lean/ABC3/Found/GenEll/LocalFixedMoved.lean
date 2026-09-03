/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.AlphaLocalInput
import ABC3.Meta.Claim

/-!
# 第 1309 ブロック —— **局所で「固定点」と「動く点」を同時に取る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1308 が受け取る形

第 1308（輸送）は次の 2 つを受け取る:

* `σ_M` が固定する位数 `l` の点
* `σ_M` が動かす `l`-捩れ点

☆本ブロックはそれを**基礎体 `K₀` の 2 つのデータ**から作る:

| 入力 | 出力 |
|---|---|
| `P₀ ∈ E(K₀)` が位数 `l`（`μ_l`、第 1297） | 固定点（第 1299＋`addOrderOf_rhPoint`） |
| `K₀` の `l`-捩れは `l` 個以下（第 1304） | 動く点と `σ_M`（第 1306→1305） |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Classical

variable {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**局所で「固定点」と「動く点」を同時に取る**——★**無条件**（第 1309）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1308（輸送）が受け取る形である。 -/
theorem exists_local_fixed_moved [IsAlgClosed M] [CharZero K₀] [IsGalois K₀ M]
    (W : WeierstrassCurve K₀) [W.IsElliptic]
    [WeierstrassCurve.IsElliptic (W.baseChange M).toAffine]
    (l : ℕ) (hl : l.Prime)
    (P₀ : W.toAffine.Point) (hP₀ : addOrderOf P₀ = l)
    (hcard : ∀ T : Finset (W.toAffine.Point), (∀ p ∈ T, l • p = 0) → T.card ≤ l) :
    ∃ (σM : M ≃ₐ[K₀] M) (Pfix Qmov : (W.baseChange M).toAffine.Point),
      addOrderOf Pfix = l ∧ galPoint W σM Pfix = Pfix ∧
        l • Qmov = 0 ∧ galPoint W σM Qmov ≠ Qmov := by
  haveI hEM : (W.map (algebraMap K₀ M)).IsElliptic :=
    inferInstanceAs ((W.baseChange M).IsElliptic)
  obtain ⟨σM, Qmov, hQ, hmoved⟩ := exists_sigma_moved_of_card_le (M := M) W l hl hcard
  refine ⟨σM, rhPoint (algebraMap K₀ M) W P₀, Qmov, ?_, ?_, hQ, hmoved⟩
  · exact (addOrderOf_rhPoint (algebraMap K₀ M) W P₀).trans hP₀
  · exact galPoint_rhPoint_eq W σM P₀

/-! ## ★出典の紐付け(`.src`) -/

def exists_local_fixed_moved.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所で固定点と動く点を同時に取る。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_local_fixed_moved.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_sigma_moved_of_card_le(第 1307、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_sigma_moved_of_card_le") 1,
    .citation "[ABC3]" "galPoint_rhPoint_eq(第 1299、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galPoint_rhPoint_eq") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1309）**——第 1308（輸送）が受け取る形である。" ++
       "☆入力は「`P₀ ∈ E(K₀)` が位数 `l`」と「`K₀` の `l`-捩れは `l` 個以下」の 2 つだけ" ++
       "——どちらも `K₀` の上の Tate 一意化から出る（第 1297・1304）。") 2 ]

end ABC3.Found.GenEll
