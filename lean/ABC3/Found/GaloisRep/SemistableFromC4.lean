/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.LocalHeightDelta
import ABC3.Found.GaloisRep.SemistableCriterion
import ABC3.Meta.Claim

/-!
# 第 1322 ブロック —— **`c₄` が単元の整モデルは半安定**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——半安定性を出す**一番安い形**

`SemistableAt` の第 2 の選択肢は「極小モデルで `v_p(c₄) = 0`」である。
★整モデルで `v_p(c₄) = 0` なら、そのモデル自身が極小（`isMinimal_of_c4_valAdd_eq_zero`、在庫）
なので、変数変換は `1` でよい。

☆したがって **`c₄` が単元の整モデルがあれば、それだけで半安定**である。

★★★これが `VeluQuotOK` の半安定性を悪い素点で出すときの受け口になる
——Tate モデル（`c₄ = 1 + 240σ₃` は単元）に変数変換で移せばよいからである。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Meta

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★★★★★★★★★★★
**`c₄` が単元の整モデルは半安定**——★**無条件**（第 1322）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`SemistableAt` の第 2 の選択肢に変数変換 `1` を渡すだけである。 -/
theorem semistableAt_of_c4_valAdd_zero (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0) :
    SemistableAt p W := by
  have h1 : ((1 : WeierstrassCurve.VariableChange L) • W) = W := one_smul _ _
  refine Or.inr ⟨1, ?_, ?_⟩
  · simp only [h1]
    exact isMinimal_of_c4_valAdd_eq_zero p W hΔ hc4ne hc4
  · simp only [h1]
    exact ⟨hc4ne, hc4⟩

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_of_c4_valAdd_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(c₄ が単元の整モデルは半安定。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_of_c4_valAdd_zero.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMinimal_of_c4_valAdd_eq_zero(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isMinimal_of_c4_valAdd_eq_zero") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1322）**——半安定性を出す**一番安い形**である。" ++
       "☆`VeluQuotOK` の半安定性を悪い素点で出すときの受け口になる" ++
       "——Tate モデル（`c₄ = 1 + 240σ₃` は単元）に変数変換で移せばよい。") 2 ]

end ABC3.Found.GaloisRep
