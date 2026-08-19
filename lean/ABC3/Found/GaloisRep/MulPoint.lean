import ABC3.Found.GaloisRep.XYStep

/-!
# Galois (G1) 第 52 ブロック —— **★★★★★★★乗法公式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★ここで乗法公式が閉じる

    ΨSq_k(x) ≠ 0  (1 ≤ k ≤ n)  ⟹  n•P = (Φ_n/ΨSq_n, …)

★§9-393 のとおり、**仮定を帰納に載せた**ので退化の場合分けが要らない。

## ★★帰納の形

| 段 | 使うもの |
|---|---|
| `n = 1` | ★`ΨSq_1 = 1`、`Φ_1 = X`、`preΨ_2 = 1` |
| `n = 2` | ★第 46 ブロック(倍化——第 30・45 ブロック) |
| `n = m+3` | ★★`(m+2)P + P`、第 50(退化なし)+ 第 51(帰納段) |

★★★**2 歩の帰納**である——`x((n+1)P)` は `x((n−1)P)` と `ψ(nP)` から出るので、
`n−1` と `n` の両方が要る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `MulOK` | ★乗法公式の述語(退化なしの版) |
| `mulOK_one` / `mulOK_two` | ★基底 |
| `mulOK_step` | ★★★★★帰納段(点の水準) |
| `mulOK_of_ne` | ★★★★★★★**乗法公式** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★**`n • P` が乗法公式で書けること**(退化なしの版)。 -/
def MulOK {x y : F} (h : W.toAffine.Nonsingular x y) (n : ℕ) : Prop :=
  ∃ (x' y' : F) (h' : W.toAffine.Nonsingular x' y'),
    n • (Point.some x y h) = Point.some x' y' h'
    ∧ x' * (W.ΨSq (n : ℤ)).eval x = (W.Φ (n : ℤ)).eval x
    ∧ (y' - W.toAffine.negY x' y') * (W.ΨSq (n : ℤ)).eval x ^ 2
      = (W.preΨ (2 * (n : ℤ))).eval x * (y - W.toAffine.negY x y)

theorem mulOK_one {x y : F} (h : W.toAffine.Nonsingular x y) : MulOK W h 1 := by
  refine ⟨x, y, h, by rw [one_smul], ?_, ?_⟩
  · norm_num [WeierstrassCurve.ΨSq_one, WeierstrassCurve.Φ_one]
  · norm_num [WeierstrassCurve.ΨSq_one, show (2 : ℤ) * 1 = 2 by norm_num,
      WeierstrassCurve.preΨ_two]

theorem mulOK_two {x y : F} (h : W.toAffine.Nonsingular x y)
    (h2 : (W.ΨSq ((2 : ℕ) : ℤ)).eval x ≠ 0) : MulOK W h 2 := by
  have hy : y ≠ W.toAffine.negY x y := by
    intro hy
    apply h2
    rw [show (((2 : ℕ)) : ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two]
    exact (psi2Sq_eval_eq_zero_iff W h.left).2 hy
  have hQ : (2 : ℕ) • (Point.some x y h)
      = Point.some _ _ (nonsingular_add h h fun hxy => hy hxy.right) := by
    rw [two_nsmul]; exact Point.add_self_of_Y_ne hy
  obtain ⟨-, hform⟩ := mulFormula_two W h
  obtain ⟨ha, hb⟩ := hform _ _ _ hQ
  exact ⟨_, _, _, hQ, ha, hb⟩

/-- ★★★★★**帰納段**——`(m+1)P`, `(m+2)P` から `(m+3)P`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★退化しないこと(第 50)を使って `Point.add_of_X_ne` に持ち込み、
★★体の水準の恒等式(第 51)へ渡す。 -/
theorem mulOK_step {x y : F} (h : W.toAffine.Nonsingular x y) (m : ℕ)
    (hprev : MulOK W h (m + 1)) (hcur : MulOK W h (m + 2))
    (hp : (W.ΨSq ((m : ℤ) + 3)).eval x ≠ 0)
    (hm1 : (W.ΨSq ((m : ℤ) + 1)).eval x ≠ 0)
    (hc : (W.ΨSq ((m : ℤ) + 2)).eval x ≠ 0) :
    MulOK W h (m + 3) := by
  obtain ⟨xm, ym, hm', hPm, hxm, hym⟩ := hcur
  obtain ⟨xmm, ymm, hmm', hPmm, hxmm, hymm⟩ := hprev
  push_cast at hxm hym hxmm hymm
  have hxne : xm ≠ x := by
    refine x_ne_of_PSq_ne W ((m : ℤ) + 2) x xm hxm ?_ ?_
    · rw [show (m : ℤ) + 2 + 1 = (m : ℤ) + 3 by ring]; exact hp
    · rw [show (m : ℤ) + 2 - 1 = (m : ℤ) + 1 by ring]; exact hm1
  have hxne' : x - xm ≠ 0 := sub_ne_zero.mpr (Ne.symm hxne)
  have hsucc : (m + 3 : ℕ) • (Point.some x y h)
      = Point.some xm ym hm' + Point.some x y h := by
    rw [show m + 3 = (m + 2) + 1 from rfl, succ_nsmul, hPm]
  rw [Point.add_of_X_ne hxne] at hsucc
  have key : (m + 2 : ℕ) • (Point.some x y h) + (- (Point.some x y h))
      = (m + 1 : ℕ) • (Point.some x y h) := by
    rw [show (m + 2 : ℕ) = (m + 1) + 1 from rfl, succ_nsmul]
    exact add_neg_cancel_right _ _
  have hpred : Point.some xm ym hm'
        + Point.some x (W.toAffine.negY x y) ((nonsingular_neg ..).mpr h)
      = Point.some xmm ymm hmm' := by
    rw [← Point.neg_some h, ← hPm, ← hPmm]; exact key
  rw [Point.add_of_X_ne hxne] at hpred
  have hxmm_eq : xmm = W.toAffine.addX xm x (W.toAffine.slope xm x ym (W.toAffine.negY x y)) :=
    (Point.some.inj hpred).1.symm
  have hax := W.toAffine.addX_eq_addX_negY_sub (x₁ := xm) (x₂ := x) ym y hxne
  have hnew : (x - xm) ^ 2
      * (xmm - W.toAffine.addX xm x (W.toAffine.slope xm x ym y))
      = (ym - W.toAffine.negY xm ym) * (y - W.toAffine.negY x y) := by
    rw [hax, hxmm_eq]; field_simp; ring
  have hxmm' : xmm * (W.ΨSq ((m : ℤ) + 2 - 1)).eval x = (W.Φ ((m : ℤ) + 2 - 1)).eval x := by
    rw [show (m : ℤ) + 2 - 1 = (m : ℤ) + 1 by ring]; exact hxmm
  have hAmm : (W.ΨSq ((m : ℤ) + 2 - 1)).eval x ≠ 0 := by
    rw [show (m : ℤ) + 2 - 1 = (m : ℤ) + 1 by ring]; exact hm1
  have hxnew : (W.toAffine.addX xm x (W.toAffine.slope xm x ym y))
      * (W.ΨSq ((m : ℤ) + 3)).eval x = (W.Φ ((m : ℤ) + 3)).eval x := by
    rw [show (m : ℤ) + 3 = ((m : ℤ) + 2) + 1 by ring]
    exact xstep W ((m : ℤ) + 2) x xm xmm _ _ _ hxm hxmm' hym
      ((psi2Sq_eval W h.left).symm) hnew hAmm
  have hay := W.toAffine.addY_sub_negY_addY (x₁ := xm) (x₂ := x) ym y hxne
  simp only [] at hay
  have hnewy : ((W.toAffine.addY xm x ym (W.toAffine.slope xm x ym y))
        - W.toAffine.negY (W.toAffine.addX xm x (W.toAffine.slope xm x ym y))
            (W.toAffine.addY xm x ym (W.toAffine.slope xm x ym y))) * (x - xm)
      = (y - W.toAffine.negY x y)
          * (xm - W.toAffine.addX xm x (W.toAffine.slope xm x ym y))
        - (ym - W.toAffine.negY xm ym)
          * (x - W.toAffine.addX xm x (W.toAffine.slope xm x ym y)) := by
    rw [hay]; field_simp
  have hynew := ystep W ((m : ℤ) + 2) x xm _ _ _ _ hxm
    (by rw [show (m : ℤ) + 2 + 1 = (m : ℤ) + 3 by ring]; exact hxnew) hym hnewy hc hxne'
  rw [show (m : ℤ) + 2 + 1 = (m : ℤ) + 3 by ring] at hynew
  exact ⟨_, _, _, hsucc, by push_cast; exact hxnew, by push_cast; exact hynew⟩

/-- ★★★★★★★**乗法公式**——`ΨSq_k(x) ≠ 0` (`1 ≤ k ≤ n`) なら `n • P` が座標で書ける。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが `E[n]` の有限性(G1)の心臓である。 -/
theorem mulOK_of_ne {x y : F} (h : W.toAffine.Nonsingular x y) :
    ∀ n : ℕ, 1 ≤ n → (∀ k : ℕ, 1 ≤ k → k ≤ n → (W.ΨSq (k : ℤ)).eval x ≠ 0) → MulOK W h n := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro hn hne
    match n, hn with
    | 1, _ => exact mulOK_one W h
    | 2, _ => exact mulOK_two W h (hne 2 (by omega) (by omega))
    | (m + 3), _ =>
        have e1 : ((m + 1 : ℕ) : ℤ) = (m : ℤ) + 1 := by push_cast; ring
        have e2 : ((m + 2 : ℕ) : ℤ) = (m : ℤ) + 2 := by push_cast; ring
        have e3 : ((m + 3 : ℕ) : ℤ) = (m : ℤ) + 3 := by push_cast; ring
        refine mulOK_step W h m
          (ih (m + 1) (by omega) (by omega) (fun k hk1 hk2 => hne k hk1 (by omega)))
          (ih (m + 2) (by omega) (by omega) (fun k hk1 hk2 => hne k hk1 (by omega)))
          (e3 ▸ hne (m + 3) (by omega) (by omega))
          (e1 ▸ hne (m + 1) (by omega) (by omega))
          (e2 ▸ hne (m + 2) (by omega) (by omega))

/-! ## ★出典の紐付け(`.src`) -/

def mulOK_of_ne.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式——n 倍点の座標)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
