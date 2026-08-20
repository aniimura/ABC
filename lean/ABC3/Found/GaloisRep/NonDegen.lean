import ABC3.Found.GaloisRep.YSide

/-!
# Galois (G1) 第 50 ブロック —— **★★★★退化しないこと**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★§9-393 の要——仮定を帰納に載せると退化が消える

乗法公式の帰納段 `nP + P` で危ないのは **`x(nP) = x(P)`**(倍化 or `nP = −P`)である。
★ところが `ΨSq_{n±1}(x) ≠ 0` を帰納の仮定に載せておけば、これは**起きない**:

    x(nP) = x  ⟺  preΨ_{n+1}(x)·preΨ_{n−1}(x)·f_n(x) = 0

で、`ΨSq_{n±1} = preΨ_{n±1}²·f_n` だから 3 因子とも `0` でない ✅

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `PSq_eval` | ★`ΨSq` の点での値 |
| `PSq_mul_PSq` | ★★`ΨSq_{n+1}·ΨSq_{n−1} = (preΨ_{n+1}preΨ_{n−1}f_n)²` |
| `x_ne_of_PSq_ne` | ★★★★**退化しない** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

variable {F : Type u} [Field F] (W : WeierstrassCurve F)

/-- ★`ΨSq` の点での値。 -/
theorem PSq_eval (n : ℤ) (x : F) :
    (W.ΨSq n).eval x
      = (W.preΨ n).eval x ^ 2 * (if Even n then W.Ψ₂Sq.eval x else 1) := by
  simp only [WeierstrassCurve.ΨSq]
  by_cases hn : Even n <;> simp [hn]

/-- ★★`ΨSq_{n+1}·ΨSq_{n−1}` は `(preΨ_{n+1}preΨ_{n−1}f_n)²` である。

★`n±1` は同じパリティなので `ΨSq_{n±1} = preΨ_{n±1}²·f_n` となり、
★★積が完全平方になる——これが帰納段で `(x − x_n)²` を消す。 -/
theorem PSq_mul_PSq (m : ℤ) (x : F) :
    (W.ΨSq (m + 1)).eval x * (W.ΨSq (m - 1)).eval x
      = ((W.preΨ (m + 1)).eval x * (W.preΨ (m - 1)).eval x
          * (if Even m then (1 : F) else W.Ψ₂Sq.eval x)) ^ 2 := by
  rw [PSq_eval, PSq_eval]
  by_cases hn : Even m
  · rw [if_neg (fun hh => (Int.even_add_one.mp hh) hn),
      if_neg (fun hh => (Int.even_sub_one.mp hh) hn), if_pos hn]; ring
  · rw [if_pos (Int.even_add_one.mpr hn), if_pos (Int.even_sub_one.mpr hn), if_neg hn]; ring

/-- ★★★★**退化しない**——`ΨSq_{m±1}(x) ≠ 0` なら `x(mP) ≠ x`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが §9-393 の要である。帰納段の場合分けが**丸ごと消える**。 -/
theorem x_ne_of_PSq_ne (m : ℤ) (x xm : F)
    (hxm : xm * (W.ΨSq m).eval x = (W.Φ m).eval x)
    (hp : (W.ΨSq (m + 1)).eval x ≠ 0) (hm : (W.ΨSq (m - 1)).eval x ≠ 0) :
    xm ≠ x := by
  intro hx
  subst hx
  have h := x_mul_PSq_sub_Phi W m xm
  rw [hxm, sub_self] at h
  have hf : (if Even m then (1 : F) else W.Ψ₂Sq.eval xm) ≠ 0 := by
    by_cases hn : Even m
    · rw [if_pos hn]; exact one_ne_zero
    · rw [if_neg hn]
      intro h0
      exact hp (by rw [PSq_eval, if_pos (Int.even_add_one.mpr hn), h0, mul_zero])
  have hpp : (W.preΨ (m + 1)).eval xm ≠ 0 := by
    intro h0; apply hp; rw [PSq_eval, h0]; ring
  have hpm : (W.preΨ (m - 1)).eval xm ≠ 0 := by
    intro h0; apply hm; rw [PSq_eval, h0]; ring
  exact (mul_ne_zero (mul_ne_zero hpp hpm) hf) h.symm

/-! ## ★出典の紐付け(`.src`) -/

def x_ne_of_PSq_ne.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式の帰納で退化が起きないこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
