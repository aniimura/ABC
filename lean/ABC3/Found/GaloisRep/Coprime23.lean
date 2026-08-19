import ABC3.Found.GaloisRep.EdsWard

/-!
# Galois (G1) 第 59 ブロック —— **★★★★★`Ψ₂Sq` と `Ψ₃` は共通根を持たない**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★零集合が部分群であることの唯一の例外を潰す

Ward の定理(第 58)から「`ψ` の零集合は位数の倍数全体」が降下で出るが、
★その論法の途中で **`v_K = v_{K+1} = 0`(連続 2 項が消える)** を排除する必要があり、
実測するとそれは `K ≤ 3` に絞られ、さらに

| `K` | 意味 |
|---|---|
| 2 | ★`Ψ₂Sq(x) = 0` かつ `Ψ₃(x) = 0` |
| 3 | ★`Ψ₃(x) = 0` かつ `preΨ₄(x) = 0` |

★★`K = 3` は**乗法公式で潰せる**(位数 3 の点では `preΨ₄(x) = −Ψ₂Sq(x)²`)。
★★★`K = 2` だけが終結式を要する。

## ★★★★終結式の実測

    Res(Ψ₂Sq, Ψ₃) = −Δ²      (b-不変量の関係 4b₈ = b₂b₆ − b₄² のもとで)

★Sylvester 行列の余因子から Bezout 係数を出した(合計 62 項)。
★★点で評価すれば **`R` の中の等式**になるので、`linear_combination` 1 発で済む。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `no_common_root_23` | ★★★★★**共通根があれば `Δ² = 0`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★★★★★**`Ψ₂Sq` と `Ψ₃` に共通根があれば `Δ² = 0`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★係数は Sylvester 行列の余因子(Python で算出、§9-397)。 -/
theorem no_common_root_23 (x : R) (h2 : W.Ψ₂Sq.eval x = 0) (h3 : W.Ψ₃.eval x = 0) :
    W.Δ ^ 2 = 0 := by
  simp only [WeierstrassCurve.Ψ₂Sq, eval_add, eval_mul, eval_pow, eval_C, eval_X] at h2
  simp only [WeierstrassCurve.Ψ₃, eval_add, eval_mul, eval_pow, eval_C, eval_X, eval_ofNat] at h3
  have hb := W.b_relation
  simp only [WeierstrassCurve.Δ]
  linear_combination
    (-(-6*W.b₂^3*W.b₆ + 6*W.b₂^2*W.b₄^2 - 12*W.b₂^2*W.b₈ + 252*W.b₂*W.b₄*W.b₆ - 216*W.b₄^3
        + 288*W.b₄*W.b₈ - 972*W.b₆^2) * x^3
     - (-2*W.b₂^4*W.b₆ + 2*W.b₂^3*W.b₄^2 - W.b₂^3*W.b₈ + 81*W.b₂^2*W.b₄*W.b₆ - 72*W.b₂*W.b₄^3
        - 351*W.b₂*W.b₆^2 + 108*W.b₄^2*W.b₆ + 432*W.b₆*W.b₈) * x^2
     - (W.b₂^4*W.b₈ - 7*W.b₂^3*W.b₄*W.b₆ + 6*W.b₂^2*W.b₄^3 - 50*W.b₂^2*W.b₄*W.b₈
        - 3*W.b₂^2*W.b₆^2 + 288*W.b₂*W.b₄^2*W.b₆ + 168*W.b₂*W.b₆*W.b₈ - 216*W.b₄^4
        + 432*W.b₄^2*W.b₈ - 1134*W.b₄*W.b₆^2 - 192*W.b₈^2) * x
     - (W.b₂^3*W.b₄*W.b₈ - 4*W.b₂^3*W.b₆^2 + 3*W.b₂^2*W.b₄^2*W.b₆ - 7*W.b₂^2*W.b₆*W.b₈
        - 36*W.b₂*W.b₄^2*W.b₈ + 162*W.b₂*W.b₄*W.b₆^2 - 16*W.b₂*W.b₈^2 - 108*W.b₄^3*W.b₆
        + 432*W.b₄*W.b₆*W.b₈ - 729*W.b₆^3)) * h2
    + (-(8*W.b₂^3*W.b₆ - 8*W.b₂^2*W.b₄^2 + 16*W.b₂^2*W.b₈ - 336*W.b₂*W.b₄*W.b₆ + 288*W.b₄^3
        - 384*W.b₄*W.b₈ + 1296*W.b₆^2) * x^2
     - (2*W.b₂^4*W.b₆ - 2*W.b₂^3*W.b₄^2 - 80*W.b₂^2*W.b₄*W.b₆ + 72*W.b₂*W.b₄^3
        + 32*W.b₂*W.b₄*W.b₈ + 360*W.b₂*W.b₆^2 - 144*W.b₄^2*W.b₆ - 576*W.b₆*W.b₈) * x
     - (-W.b₂^4*W.b₈ + 5*W.b₂^3*W.b₄*W.b₆ - 4*W.b₂^2*W.b₄^3 + 48*W.b₂^2*W.b₄*W.b₈
        + W.b₂^2*W.b₆^2 - 204*W.b₂*W.b₄^2*W.b₆ - 176*W.b₂*W.b₆*W.b₈ + 144*W.b₄^4
        - 384*W.b₄^2*W.b₈ + 864*W.b₄*W.b₆^2 + 256*W.b₈^2)) * h3
    + (12*W.b₂^2*W.b₄*W.b₈ + 4*W.b₂^2*W.b₆^2 - 80*W.b₂*W.b₄^2*W.b₆ - 32*W.b₂*W.b₆*W.b₈
        + 64*W.b₄^4 - 112*W.b₄^2*W.b₈ + 324*W.b₄*W.b₆^2 + 64*W.b₈^2) * hb

/-! ## ★出典の紐付け(`.src`) -/

def no_common_root_23.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(2-分点多項式と 3-分点多項式の互いに素性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
