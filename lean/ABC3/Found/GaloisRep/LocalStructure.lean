import ABC3.Found.GaloisRep.Psi2SqSquarefree

/-!
# Galois (G5) 第 132 ブロック —— **★★★★★分岐点の局所構造**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★DVR 判定の中身

§9-443 で **`IsDedekindDomainDvr.isDedekindDomain` が instance である**ことが分かった
——「零でない素イデアルでの局所化がすべて DVR」を示せば Dedekind が出る。

★第 129 で `F[W] = F[x][z]`(`z² = Ψ₂Sq(x)`)、第 131 で `Ψ₂Sq` が squarefree と分かった。
★★極大イデアルは `x` 座標 `c` の上にあり、2 つの場合に分かれる:

| 場合 | 局所構造 |
|---|---|
| `Ψ₂Sq(c) ≠ 0`(不分岐) | `z² − Ψ₂Sq(c)` が相異なる 2 根を持ち、`x − c` が極大イデアルを生成 |
| `Ψ₂Sq(c) = 0`(分岐、= 2 等分点) | `z` が極大イデアルを生成 |

★★★本ブロックは**後者の中身**を出す:

    Ψ₂Sq = (x − c)·e   で   e(c) ≠ 0     (squarefree だから単根)
      ⟹  z² = (x − c)·e(x)  で `e(x)` は局所化で単元
      ⟹  x − c = z²/e(x) ∈ (z)   すなわち  極大イデアル = (z)

★★★★**これが Eisenstein 型の議論**である——`Ψ₂Sq` の squarefree 性が
「単根 ⟹ 分岐指数 2 ⟹ `z` が一意化元」に直結する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_factor_of_squarefree_root` | ★★★squarefree な多項式の根は単根 |
| `genZ_sq_factor` | ★★★★★**`z² = (x − c)·e(x)`、`e(c) ≠ 0`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-- ★★★**squarefree な多項式の根は単根である**——`d = (X−c)·e` で `e(c) ≠ 0`。 -/
theorem exists_factor_of_squarefree_root {d : Polynomial F} (hsq : Squarefree d)
    {c : F} (hc : d.eval c = 0) :
    ∃ e : Polynomial F, d = (Polynomial.X - Polynomial.C c) * e ∧ e.eval c ≠ 0 := by
  obtain ⟨e, he⟩ := (Polynomial.dvd_iff_isRoot.2 hc)
  refine ⟨e, he, ?_⟩
  intro hec
  obtain ⟨f, hf⟩ := (Polynomial.dvd_iff_isRoot.2 hec)
  have hdvd : (Polynomial.X - Polynomial.C c) * (Polynomial.X - Polynomial.C c) ∣ d := by
    rw [he, hf]
    exact ⟨f, by ring⟩
  exact (Polynomial.not_isUnit_X_sub_C c) (hsq _ hdvd)

variable [DecidableEq F]

/-- ★★★★★**分岐点(2 等分点)の局所構造**——`z² = (x − c)·e(x)` で `e(c) ≠ 0`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Ψ₂Sq` が squarefree(第 131)なので `c` は単根であり、`e(x)` は局所化で単元になる。
★★したがって `x − c = z²/e(x) ∈ (z)`、すなわち極大イデアルは `z` で生成される。 -/
theorem genZ_sq_factor (W : WeierstrassCurve.Affine F) [W.IsElliptic] [IsAlgClosed F]
    (h2 : IsUnit (2 : F)) {c : F} (hc : W.Ψ₂Sq.eval c = 0) :
    ∃ e : Polynomial F, e.eval c ≠ 0 ∧
      (genZ W) ^ 2 = (genX W - algebraMap F W.CoordinateRing c)
        * Polynomial.eval₂ (algebraMap F W.CoordinateRing) (genX W) e := by
  obtain ⟨e, hde, hec⟩ := exists_factor_of_squarefree_root (squarefree_Psi2Sq W h2) hc
  refine ⟨e, hec, ?_⟩
  rw [genZ_sq, hde, Polynomial.eval₂_mul, Polynomial.eval₂_sub, Polynomial.eval₂_X,
    Polynomial.eval₂_C]

/-! ## ★出典の紐付け(`.src`) -/

def genZ_sq_factor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——分岐点での局所構造 z² = (x−c)·e(x))",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
