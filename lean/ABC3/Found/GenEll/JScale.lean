/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Velu
import ABC3.Meta.Claim

/-!
# `c₄`・`c₆` が尺度倍なら `j` は等しい（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か（2026-08-31、第 834-836）

`Check/GenEll/VeluTateNeedsChange.lean` の測定により、Tate 曲線の `μ_l` による
Vélu の商は `E_{q^l}` と**等しくはない**が、

    `c₄′ = l⁴·c₄(q^l)`,   `c₆′ = −l⁶·c₆(q^l)`

という**尺度倍**の関係にあることが分かった（`l = 5, 7` で `q^21` まで数値確認）。

★★本ファイルはそこから `j′ = j(q^l)` を出す**代数の側**を与える。
☆葉 1（`jExp_velu_bad`）が要求するのは `j` の付値だけなので、これで十分である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve

variable {R : Type*} [CommRing R] [IsDomain R] [CharZero R]

/-- ★★★★★★★★★★★★★★★★
**`c₄`・`c₆` が `n⁴`・`−n⁶` 倍なら `j` は等しい**。

★`Δ` は `n¹²` 倍になり、`j = Δ⁻¹c₄³` では相殺する。 -/
theorem j_eq_of_c4_c6_scale (W W' : WeierstrassCurve R) [W.IsElliptic] [W'.IsElliptic]
    (n : R) (h4 : W.c₄ = n ^ 4 * W'.c₄) (h6 : W.c₆ = -(n ^ 6) * W'.c₆) :
    W.j = W'.j := by
  have h1728 : (1728 : R) ≠ 0 := by
    have : ((1728 : ℕ) : R) ≠ 0 := Nat.cast_ne_zero.2 (by norm_num)
    simpa using this
  have hΔ : W.Δ = n ^ 12 * W'.Δ := by
    refine mul_left_cancel₀ h1728 ?_
    have e1 : (1728 : R) * W.Δ = W.c₄ ^ 3 - W.c₆ ^ 2 := W.c_relation
    have e2 : (1728 : R) * W'.Δ = W'.c₄ ^ 3 - W'.c₆ ^ 2 := W'.c_relation
    rw [e1, h4, h6, show (1728 : R) * (n ^ 12 * W'.Δ) = n ^ 12 * ((1728 : R) * W'.Δ) from by ring,
      e2]
    ring
  have hcoe : ((W.Δ' : Rˣ) : R) = n ^ 12 * ((W'.Δ' : Rˣ) : R) := by
    rw [coe_Δ', coe_Δ', hΔ]
  have hUne : ((W.Δ' : Rˣ) : R) ≠ 0 := (W.Δ').ne_zero
  have hU'ne : ((W'.Δ' : Rˣ) : R) ≠ 0 := (W'.Δ').ne_zero
  have hinvW : ((W.Δ' : Rˣ) : R) * (((W.Δ')⁻¹ : Rˣ) : R) = 1 := by
    rw [← Units.val_mul, mul_inv_cancel, Units.val_one]
  have hinvW' : ((W'.Δ' : Rˣ) : R) * (((W'.Δ')⁻¹ : Rˣ) : R) = 1 := by
    rw [← Units.val_mul, mul_inv_cancel, Units.val_one]
  have hkey : (((W.Δ')⁻¹ : Rˣ) : R) * n ^ 12 = (((W'.Δ')⁻¹ : Rˣ) : R) := by
    refine mul_left_cancel₀ hUne ?_
    calc ((W.Δ' : Rˣ) : R) * ((((W.Δ')⁻¹ : Rˣ) : R) * n ^ 12)
        = (((W.Δ' : Rˣ) : R) * (((W.Δ')⁻¹ : Rˣ) : R)) * n ^ 12 := by ring
      _ = n ^ 12 := by rw [hinvW, one_mul]
      _ = (n ^ 12 * ((W'.Δ' : Rˣ) : R)) * (((W'.Δ')⁻¹ : Rˣ) : R) := by
          rw [mul_assoc, hinvW', mul_one]
      _ = ((W.Δ' : Rˣ) : R) * (((W'.Δ')⁻¹ : Rˣ) : R) := by rw [← hcoe]
  rw [j, j, h4]
  calc (((W.Δ')⁻¹ : Rˣ) : R) * (n ^ 4 * W'.c₄) ^ 3
      = ((((W.Δ')⁻¹ : Rˣ) : R) * n ^ 12) * W'.c₄ ^ 3 := by ring
    _ = (((W'.Δ')⁻¹ : Rˣ) : R) * W'.c₄ ^ 3 := by rw [hkey]

/-- ★★★★★★**符号が `+` の版**——`c₆ = n⁶c₆'` でも `j` は一致する
（`c₆` は 2 乗でしか入らないから）。 -/
theorem j_eq_of_c4_c6_scale_pos (W W' : WeierstrassCurve R) [W.IsElliptic] [W'.IsElliptic]
    (n : R) (h4 : W.c₄ = n ^ 4 * W'.c₄) (h6 : W.c₆ = n ^ 6 * W'.c₆) :
    W.j = W'.j := by
  have h1728 : (1728 : R) ≠ 0 := by
    have : ((1728 : ℕ) : R) ≠ 0 := Nat.cast_ne_zero.2 (by norm_num)
    simpa using this
  have hΔ : W.Δ = n ^ 12 * W'.Δ := by
    refine mul_left_cancel₀ h1728 ?_
    have e1 : (1728 : R) * W.Δ = W.c₄ ^ 3 - W.c₆ ^ 2 := W.c_relation
    have e2 : (1728 : R) * W'.Δ = W'.c₄ ^ 3 - W'.c₆ ^ 2 := W'.c_relation
    rw [e1, h4, h6, show (1728 : R) * (n ^ 12 * W'.Δ) = n ^ 12 * ((1728 : R) * W'.Δ) from by ring,
      e2]
    ring
  have hcoe : ((W.Δ' : Rˣ) : R) = n ^ 12 * ((W'.Δ' : Rˣ) : R) := by
    rw [coe_Δ', coe_Δ', hΔ]
  have hUne : ((W.Δ' : Rˣ) : R) ≠ 0 := (W.Δ').ne_zero
  have hU'ne : ((W'.Δ' : Rˣ) : R) ≠ 0 := (W'.Δ').ne_zero
  have hinvW : ((W.Δ' : Rˣ) : R) * (((W.Δ')⁻¹ : Rˣ) : R) = 1 := by
    rw [← Units.val_mul, mul_inv_cancel, Units.val_one]
  have hinvW' : ((W'.Δ' : Rˣ) : R) * (((W'.Δ')⁻¹ : Rˣ) : R) = 1 := by
    rw [← Units.val_mul, mul_inv_cancel, Units.val_one]
  have hkey : (((W.Δ')⁻¹ : Rˣ) : R) * n ^ 12 = (((W'.Δ')⁻¹ : Rˣ) : R) := by
    refine mul_left_cancel₀ hUne ?_
    calc ((W.Δ' : Rˣ) : R) * ((((W.Δ')⁻¹ : Rˣ) : R) * n ^ 12)
        = (((W.Δ' : Rˣ) : R) * (((W.Δ')⁻¹ : Rˣ) : R)) * n ^ 12 := by ring
      _ = n ^ 12 := by rw [hinvW, one_mul]
      _ = (n ^ 12 * ((W'.Δ' : Rˣ) : R)) * (((W'.Δ')⁻¹ : Rˣ) : R) := by
          rw [mul_assoc, hinvW', mul_one]
      _ = ((W.Δ' : Rˣ) : R) * (((W'.Δ')⁻¹ : Rˣ) : R) := by rw [← hcoe]
  rw [j, j, h4]
  calc (((W.Δ')⁻¹ : Rˣ) : R) * (n ^ 4 * W'.c₄) ^ 3
      = ((((W.Δ')⁻¹ : Rˣ) : R) * n ^ 12) * W'.c₄ ^ 3 := by ring
    _ = (((W'.Δ')⁻¹ : Rˣ) : R) * W'.c₄ ^ 3 := by rw [hkey]

/-! ## ★出典の紐付け(`.src`) -/

def j_eq_of_c4_c6_scale.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(c₄・c₆ が n⁴・−n⁶ 倍なら j は等しい。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GenEll
