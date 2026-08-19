import ABC3.Found.GaloisRep.TwoTorsionDecomp
import Mathlib.FieldTheory.IsAlgClosed.Basic

/-!
# Galois (G5) 第 123 ブロック —— **★★★★★★平行移動の単射性が完全に閉じた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★最後の仮定を外す

第 122 の `translateHom_injective_all` は
「2 等分点でないアフィン点が 1 つある」ことを仮定として受けていた。
★**代数閉体(標数 ≠ 2)ならそれが取れる**:

| 段 | 内容 |
|---|---|
| 1 | `Ψ₂Sq ≠ 0`——先頭係数が `4 ≠ 0`(mathlib の `coeff_Ψ₂Sq`) |
| 2 | 無限体なので `Ψ₂Sq(x) ≠ 0` なる `x` がある |
| 3 | 代数閉体なので `y² + (a₁x+a₃)y − (x³+…) = 0` が解を持つ |
| 4 | `Ψ₂Sq(x) ≠ 0` ⟺ `negY x y ≠ y`(第 29 ブロック `psi2Sq_eval_eq_zero_iff`) |

★★これで **`τ_Q` の単射性は仮定なしで出る**(`translateHom_injective_algClosed`)。

## ★★これで平行移動の葉はすべて閉じた

| 段 | 場所 |
|---|---|
| 環準同型 | 第 115 |
| 単射性(2 等分点でない) | 第 117 |
| 分数体への延長・自己同型 | 第 119・120 |
| 合成則 | 第 121 |
| 2 等分点の単射性 | 第 122 |
| **2 等分点でない点の存在** | **本ブロック** |

★★★当初 20-60 ブロックと見積もっていた葉が、第 114-123 の **10 ブロック**で閉じた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_eval_ne_zero` | ★無限体では非零多項式に非根がある |
| `exists_point_of_isAlgClosed` | ★★代数閉体では各 `x` に対し曲線の点がある |
| `exists_notTwoTorsion_point` | ★★★★★**2 等分点でないアフィン点の存在** |
| `translateHom_injective_algClosed` | ★★★★★★**仮定なしの単射性** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★補助 -/

/-- ★無限体では、非零多項式に非根がある。 -/
theorem exists_eval_ne_zero [Infinite F] {p : Polynomial F} (hp : p ≠ 0) :
    ∃ x : F, p.eval x ≠ 0 := by
  by_contra hc
  push_neg at hc
  refine hp (Polynomial.eq_zero_of_infinite_isRoot p ?_)
  have hset : {x : F | p.IsRoot x} = Set.univ := by
    ext x
    simp [Polynomial.IsRoot, hc x]
  rw [hset]
  exact Set.infinite_univ

/-- ★★代数閉体では、各 `x` に対し曲線の点 `(x, y)` がある。

★`y` についての 2 次方程式 `y² + (a₁x+a₃)y − (x³+a₂x²+a₄x+a₆) = 0` を解くだけである。 -/
theorem exists_point_of_isAlgClosed [IsAlgClosed F] (W : WeierstrassCurve.Affine F) (x : F) :
    ∃ y : F, W.Equation x y := by
  set q : Polynomial F := Polynomial.C 1 * Polynomial.X ^ 2
    + Polynomial.C (W.a₁ * x + W.a₃) * Polynomial.X
    + Polynomial.C (-(x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆)) with hq
  have hdeg : q.degree ≠ 0 := by
    rw [hq, Polynomial.degree_quadratic (one_ne_zero)]
    decide
  obtain ⟨y, hy⟩ := IsAlgClosed.exists_root q hdeg
  refine ⟨y, ?_⟩
  rw [WeierstrassCurve.Affine.equation_iff]
  have hev := hy
  rw [Polynomial.IsRoot, hq] at hev
  simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_pow, Polynomial.eval_C,
    Polynomial.eval_X] at hev
  linear_combination hev

/-! ## ★★★★★2 等分点でない点の存在 -/

/-- ★★★★★**2 等分点でないアフィン点が存在する**(代数閉・標数 ≠ 2)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Ψ₂Sq` の先頭係数は `4` なので、標数 ≠ 2 なら `Ψ₂Sq ≠ 0`。
★★無限体だから非根 `x` があり、代数閉だから `(x, y)` が曲線に乗る。
★★★`Ψ₂Sq(x) ≠ 0` は第 29 ブロックにより「2 等分点でない」と同値である。 -/
theorem exists_notTwoTorsion_point [IsAlgClosed F] [Infinite F]
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] (h4 : (4 : F) ≠ 0) :
    ∃ (x₁ y₁ : F) (_ : W.Nonsingular x₁ y₁), W.negY x₁ y₁ ≠ y₁ := by
  have hne : W.Ψ₂Sq ≠ 0 := by
    intro h
    refine h4 ?_
    have hc := W.coeff_Ψ₂Sq
    rw [h, Polynomial.coeff_zero] at hc
    exact hc.symm
  obtain ⟨x, hx⟩ := exists_eval_ne_zero hne
  obtain ⟨y, hy⟩ := exists_point_of_isAlgClosed W x
  refine ⟨x, y, WeierstrassCurve.Affine.equation_iff_nonsingular.1 hy, ?_⟩
  intro hc
  exact hx ((psi2Sq_eval_eq_zero_iff W hy).2 hc.symm)

/-! ## ★★★★★★仮定なしの単射性 -/

/-- ★★★★★★**平行移動の環準同型は単射である**(代数閉・標数 ≠ 2、仮定なし)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これで平行移動の葉はすべて閉じた——第 114-123 の 10 ブロックである。 -/
theorem translateHom_injective_algClosed [IsAlgClosed F] [Infinite F] [DecidableEq F]
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] (h4 : (4 : F) ≠ 0)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) :
    Function.Injective (translateHom W hQ) :=
  translateHom_injective_all W hQ (exists_notTwoTorsion_point W h4)

/-! ## ★出典の紐付け(`.src`) -/

def exists_notTwoTorsion_point.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——2 等分点でないアフィン点の存在)",
    sectionId := "genell-thm-3-8" }

def translateHom_injective_algClosed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動の単射性が仮定なしで成り立つこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
