import ABC3.Found.GaloisRep.LocalStructure

/-!
# Galois (G5) 第 133 ブロック —— **★★★★★不分岐点の局所構造**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★2 つの場合が揃った

第 132 ブロックで分岐点(`Ψ₂Sq(c) = 0`)の局所構造を出した。★本ブロックは**不分岐**の側:

    s² = Ψ₂Sq(c)  (s ≠ 0)  とすると
      (z − s)(z + s) = z² − s² = Ψ₂Sq(x) − Ψ₂Sq(c) = (x − c)·g(x)

★★極大イデアル `m = (x − c, z − s)` の局所化では `z + s` が単元
(`z ≡ s` かつ `2s ≠ 0` だから `z + s ≡ 2s ≠ 0`)なので

    z − s = (x − c)·g(x)/(z + s) ∈ (x − c)
      ⟹  m = (x − c)

★★★すなわち**不分岐点では `x − c` が一意化元**である。

## ★★★★これで DVR 判定の材料が揃った

| 場合 | 一意化元 | ブロック |
|---|---|---|
| `Ψ₂Sq(c) = 0`(分岐) | `z` | 第 132 |
| `Ψ₂Sq(c) ≠ 0`(不分岐) | `x − c` | **本ブロック** |

★残るのは (i) 極大イデアルが `x` 座標の上にあること、
(ii) `IsDiscreteValuationRing` への接続である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `genZ_sub_factor` | ★★★★★**`(z − s)(z + s) = (x − c)·g(x)`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-- ★★★★★**不分岐点の局所構造**——`(z − s)(z + s) = (x − c)·g(x)`(`s² = Ψ₂Sq(c)`)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★局所化では `z + s` が単元になるので、`z − s ∈ (x − c)`、
すなわち極大イデアルは `x − c` で生成される。 -/
theorem genZ_sub_factor (W : WeierstrassCurve.Affine F) {c s : F}
    (hs : s ^ 2 = W.Ψ₂Sq.eval c) :
    ∃ g : Polynomial F,
      (genZ W - algebraMap F W.CoordinateRing s) * (genZ W + algebraMap F W.CoordinateRing s)
        = (genX W - algebraMap F W.CoordinateRing c)
          * Polynomial.eval₂ (algebraMap F W.CoordinateRing) (genX W) g := by
  obtain ⟨g, hg⟩ := Polynomial.X_sub_C_dvd_sub_C_eval (a := c) (p := W.Ψ₂Sq)
  refine ⟨g, ?_⟩
  have hz := genZ_sq W
  have hexp : (genZ W - algebraMap F W.CoordinateRing s)
      * (genZ W + algebraMap F W.CoordinateRing s)
      = (genZ W) ^ 2 - (algebraMap F W.CoordinateRing s) ^ 2 := by ring
  rw [hexp, hz, ← map_pow, hs]
  have hsub : Polynomial.eval₂ (algebraMap F W.CoordinateRing) (genX W) W.Ψ₂Sq
      - algebraMap F W.CoordinateRing (W.Ψ₂Sq.eval c)
      = Polynomial.eval₂ (algebraMap F W.CoordinateRing) (genX W)
        (W.Ψ₂Sq - Polynomial.C (W.Ψ₂Sq.eval c)) := by
    rw [Polynomial.eval₂_sub, Polynomial.eval₂_C]
  rw [hsub, hg, Polynomial.eval₂_mul, Polynomial.eval₂_sub, Polynomial.eval₂_X,
    Polynomial.eval₂_C]

/-! ## ★出典の紐付け(`.src`) -/

def genZ_sub_factor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——不分岐点での局所構造 (z−s)(z+s) = (x−c)·g(x))",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
