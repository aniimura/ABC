/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.C6DenomFree
import ABC3.Found.GaloisRep.TateVelu
import ABC3.Meta.Claim
import ABC3.Found.GaloisRep.VeluTateMuField.Lemma35

/-!
# VeluTateMuField —— `[GenEll] Lemma 3.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Found.GenEll Finset

/-- ★★★★★★**`c₄ ↦ l⁴c₄`・`Δ ↦ l¹²Δ` なら `j` は変わらない**。

☆`j = Δ⁻¹·c₄³` なので `l¹²` が約分される。★これが「`u` 倍の変数変換では `j` が動かない」
という古典的な事実の、必要な形だけを取り出したものである。 -/
theorem j_of_c4_Delta {K : Type} [Field K] (A B : WeierstrassCurve K)
    [A.IsElliptic] [B.IsElliptic] (l : K) (hl : l ≠ 0)
    (h4 : A.c₄ = l ^ 4 * B.c₄) (hΔ : A.Δ = l ^ 12 * B.Δ) :
    A.j = B.j := by
  have hBΔ : B.Δ ≠ 0 := B.isUnit_Δ.ne_zero
  rw [j_eq_inv_Delta_mul, j_eq_inv_Delta_mul, h4, hΔ]
  field_simp

/-- ★★★★★★★★★★★★**`c₄`・`c₆` が `l⁴`・`l⁶` 倍なら `j` は等しい**（第 1132）。

★★これが `p ∣ l` で `j_velu_tate_mu_map`（`hlu` を要求する）の代わりになる形である。 -/
theorem j_eq_of_c4_c6 {K : Type} [Field K] [CharZero K] (A B : WeierstrassCurve K)
    [A.IsElliptic] [B.IsElliptic] (l : K) (hl : l ≠ 0)
    (h4 : A.c₄ = l ^ 4 * B.c₄) (h6 : A.c₆ = l ^ 6 * B.c₆) :
    A.j = B.j :=
  j_of_c4_Delta A B l hl h4 (Delta_of_c4_c6 A B l h4 h6)

def j_eq_of_c4_c6.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(c₄・c₆ が l⁴・l⁶ 倍なら j は等しい。★無条件)",
    sectionId := "genell-lemma-3-2" }

def j_eq_of_c4_c6.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "j_eq_inv_Delta_mul(j = Δ⁻¹·c₄³、在庫)"
      (.inProject "ABC3" "ABC3.Found.GenEll.j_eq_inv_Delta_mul") 1,
    .implicitStep
      ("★★**2026-09-01（第 1132）**——`j = Δ⁻¹·c₄³` なので " ++
       "`c₄ ↦ l⁴c₄`・`Δ ↦ l¹²Δ` では `l¹²` が約分される。" ++
       "☆`p ∣ l` で `j_velu_tate_mu_map`（`hlu` を要求する）の代わりになる形である。") 2 ]

end ABC3.Found.GaloisRep
