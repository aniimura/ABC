import ABC3.Found.GaloisRep.LadderThree

/-!
# Galois (G5) 第 203 ブロック —— **★★★★★★★★行列式の公式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★`det ρ = 円分指標` の数学的な中核

`det_galRep_eq_cyclotomic` が使うのは、基底 `P, Q` について

    e_n(aP + cQ, bP + dQ) = e_n(P, Q)^{ad − bc}

である。★引き算を避けるため、**両辺に `e_n(P,Q)^{bc}` を掛けた形**で述べる:

    e_n(aP + cQ, bP + dQ) · e_n(P,Q)^{bc} = e_n(P,Q)^{ad}

★★これなら指数が `ℕ` のままで済む。

### ★★★★★4 つの性質がここで合流する

| 段 | 使うもの |
|---|---|
| `e(aP+cQ, R) = e(aP,R)·e(cQ,R)` | 双線型性(第 1 変数、第 184+191) |
| `e(kP, R) = e(P,R)^k` | 第 195 |
| `e(P, bP+dQ) = e(P,bP)·e(P,dQ)` | 双線型性(第 2 変数、第 195) |
| `e(P,P) = e(Q,Q) = 1` | 交代性(第 190+191) |
| `e(P,Q)·e(Q,P) = 1` | 反対称性(第 195) |

★展開すると `e(P,Q)^{ad} · e(Q,P)^{bc}` になり、反対称性で `e(Q,P)^{bc}` を移す。

## ★★残っているのは非退化性だけ

`det_galRep_eq_cyclotomic` はこの公式に **`e_n(P,Q)` が原始 `n` 乗根であること**
(= 非退化性)を足せば出る。★非退化性は第 197 で `F(E)^{E[n]} = [n]^*F(E)` の 1 つに絞れており、
それは `deg[n] = n²`、すなわち EDS 恒等式(**mathlib 自身の TODO**)に帰着している。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `weilPairingVal_det` | ★★★★★★★★**`e_n(aP+cQ, bP+dQ)·e_n(P,Q)^{bc} = e_n(P,Q)^{ad}`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F] [CharZero F] [IsAlgClosed F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-- ★★★★★★★★**行列式の公式**——`e_n(aP+cQ, bP+dQ)·e_n(P,Q)^{bc} = e_n(P,Q)^{ad}`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★両変数の双線型性・交代性・反対称性がここで合流する。
★★引き算を避けるため `e_n(P,Q)^{bc}` を左辺に掛けた形で述べる。 -/
theorem weilPairingVal_det (n : ℕ) (hn : 1 ≤ n) (P Q : W.Point)
    (hP : n • P = 0) (hQ : n • Q = 0) (a b c d : ℕ) :
    weilPairingVal W n (a • P + c • Q) (b • P + d • Q) * (weilPairingVal W n P Q) ^ (b * c)
      = (weilPairingVal W n P Q) ^ (a * d) := by
  have haP : n • (a • P) = 0 := by rw [smul_comm, hP, smul_zero]
  have hcQ : n • (c • Q) = 0 := by rw [smul_comm, hQ, smul_zero]
  have hbP : n • (b • P) = 0 := by rw [smul_comm, hP, smul_zero]
  have hdQ : n • (d • Q) = 0 := by rw [smul_comm, hQ, smul_zero]
  have hR : n • (b • P + d • Q) = 0 := by rw [nsmul_add, hbP, hdQ, add_zero]
  rw [weilPairingVal_add_left W n hn (a • P) (c • Q) _ haP hcQ hR]
  rw [weilPairingVal_nsmul_left W n hn P _ hP hR a,
    weilPairingVal_nsmul_left W n hn Q _ hQ hR c]
  rw [weilPairingVal_add_right W n hn P (b • P) (d • Q) hP hbP hdQ,
    weilPairingVal_add_right W n hn Q (b • P) (d • Q) hQ hbP hdQ]
  rw [weilPairingVal_nsmul_right W n hn P P hP hP b,
    weilPairingVal_nsmul_right W n hn P Q hP hQ d,
    weilPairingVal_nsmul_right W n hn Q P hQ hP b,
    weilPairingVal_nsmul_right W n hn Q Q hQ hQ d]
  rw [weilPairingVal_alt W n hn P hP, weilPairingVal_alt W n hn Q hQ]
  have hanti := weilPairingVal_antisymm W n hn P Q hP hQ
  rw [one_pow, one_mul, one_pow, mul_one, ← pow_mul, ← pow_mul]
  rw [mul_assoc, ← mul_pow, mul_comm (weilPairingVal W n Q P) (weilPairingVal W n P Q),
    hanti, one_pow, mul_one, Nat.mul_comm]

/-! ## ★出典の紐付け(`.src`) -/

def weilPairingVal_det.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の行列式の公式)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
