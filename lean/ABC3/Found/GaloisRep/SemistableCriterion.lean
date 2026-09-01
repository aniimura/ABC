/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Meta.Claim

/-!
# 第 1269 ブロック —— **`v_p(Δ_min) = max(0, −v_p(j))` なら半安定**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

## ★★★★★★★★★★★★★★★★★★★★これは何か——最後の数学の穴

第 1268 で「残る穴は `v_p(Δ_min) = max(0, −v_p(j)) ⟹ 半安定` の 1 命題」と測った。
★本ブロックはそれを取る。

☆在庫の `valAdd_Delta_eq_neg_jExp`（半安定 ⟹ `v(Δ) = −v(j)`）の**逆**であり、
証明は同じ単元の等式 `j = Δ⁻¹c₄³` を逆向きに読むだけである。

★★★これで `VeluQuotOK` の半安定性が、商の `Δ_min` の関係
（`minDeltaExp_eq_mul_at_bad_prime`、証明済み）から出せるようになる。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Meta

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★★★★★★★
**`v(Δ) = −v(j)` なら `v(c₄) = 0`**——★**無条件**（第 1269）。

☆`valAdd_Delta_eq_neg_jExp`（在庫）の逆向き。
`j = Δ⁻¹c₄³` から `v(j) = 3v(c₄) − v(Δ)` なので、
`v(Δ) = −v(j)` を代入すると `v(c₄) = 0` になる。 -/
theorem valAdd_c4_eq_zero_of_Delta_eq_neg_jExp (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (h : valAdd p (Units.mk0 W.Δ hΔ) = - jExp p W) :
    valAdd p (Units.mk0 W.c₄ hc4ne) = 0 := by
  have hjeq : W.j = W.Δ⁻¹ * W.c₄ ^ 3 := ABC3.Found.GenEll.j_eq_inv_Delta_mul W
  have hj : W.j ≠ 0 := by
    rw [hjeq]
    exact mul_ne_zero (inv_ne_zero hΔ) (pow_ne_zero 3 hc4ne)
  have hunit : Units.mk0 W.j hj
      = (Units.mk0 W.c₄ hc4ne) ^ 3 * (Units.mk0 W.Δ hΔ)⁻¹ := by
    ext
    show W.j = ((Units.mk0 W.c₄ hc4ne) ^ 3 * (Units.mk0 W.Δ hΔ)⁻¹ : Lˣ)
    push_cast
    rw [hjeq]
    show W.Δ⁻¹ * W.c₄ ^ 3 = W.c₄ ^ 3 * (W.Δ)⁻¹
    ring
  have hJ : jExp p W = valAdd p (Units.mk0 W.j hj) := dif_neg hj
  rw [hunit, valAdd_mul, valAdd_pow, valAdd_inv, h] at hJ
  omega

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**`v_p(Δ_min) = max(0, −v_p(j))` なら半安定**——★**無条件**（第 1269）。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

☆`v_p(j) ≥ 0` なら `v_p(Δ_min) = 0`（良還元）。
★`v_p(j) < 0` なら極小モデルで `v(Δ) = −v(j)` になり、
上の補題で `v(c₄) = 0`（乗法還元）。

★★★これが第 1268 で名指しした**最後の数学の穴**であり、
`minDeltaExp_le_maxJ`（在庫、半安定 ⟹ 不等式）の**逆向き**である。 -/
theorem semistableAt_of_minDeltaExp_eq (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] (hΔ : W.Δ ≠ 0)
    (C : WeierstrassCurve.VariableChange L) (hC : IsMinimal (primeSubring p) (C • W))
    (h : minDeltaExp p W = max 0 (- jExp p W)) :
    SemistableAt p W := by
  by_cases hneg : jExp p W < 0
  · right
    refine ⟨C, hC, ?_⟩
    have hΔC : (C • W).Δ ≠ 0 := variableChange_Delta_ne_zero W hΔ C
    have hjC : jExp p (C • W) = jExp p W := jExp_variableChange p W C
    have hjne : (C • W).j ≠ 0 := by
      intro h0
      have : jExp p (C • W) = 0 := dif_pos h0
      omega
    have hc4ne : (C • W).c₄ ≠ 0 := by
      intro h0
      refine hjne ?_
      rw [ABC3.Found.GenEll.j_eq_inv_Delta_mul, h0]
      ring
    refine ⟨hc4ne, ?_⟩
    have hd : valAdd p (Units.mk0 ((C • W).Δ) hΔC) = minDeltaExp p W :=
      (minDeltaExp_eq p W hΔ C hC).symm
    rw [h, max_eq_right (by omega)] at hd
    exact valAdd_c4_eq_zero_of_Delta_eq_neg_jExp p (C • W) hΔC hc4ne (by rw [hd, hjC])
  · left
    rw [h, max_eq_left (by omega)]

/-! ## ★出典の紐付け(`.src`) -/

def valAdd_c4_eq_zero_of_Delta_eq_neg_jExp.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(v(Δ) = −v(j) なら v(c₄) = 0。★無条件)",
    sectionId := "genell-lemma-3-7" }

def semistableAt_of_minDeltaExp_eq.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(v_p(Δ_min) = max(0, −v_p(j)) なら半安定。★無条件)",
    sectionId := "genell-lemma-3-7" }

def semistableAt_of_minDeltaExp_eq.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1269）**——第 1268 で名指しした**最後の数学の穴**である。" ++
       "☆`minDeltaExp_le_maxJ`（在庫、半安定 ⟹ 不等式）の**逆向き**であり、" ++
       "証明は `j = Δ⁻¹c₄³` を逆に読むだけだった。" ++
       "★★★これで `VeluQuotOK` の半安定性が、商の `Δ_min` の関係" ++
       "（`minDeltaExp_eq_mul_at_bad_prime`、証明済み）から出せる。") 3 ]

end ABC3.Found.GaloisRep
