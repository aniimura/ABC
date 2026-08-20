import ABC3.Found.GaloisRep.CoordinateDimension

/-!
# Galois (G5) 第 128 ブロック —— **★★★★★座標環の単元は定数**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★`g_P` の最終段で要る

層 3 の目標は `g_P^n = f_P ∘ [n]` である。★イデアルの言葉で回すと、
まず **`g_P^n` と `μ(f_P)` が同じイデアルを生成する**ところまで出て、
そこから両者が**単元倍で一致する**ことが分かる。

★★その単元が**定数**でなければ、`n` 乗根に吸収できない。
代数閉体なら定数の `n` 乗根は取れるので、これで `g_P` が確定する。

## ★★★★機構——ノルムの次数

mathlib の `CoordinateRing.degree_norm_smul_basis`:

    deg N(p·1 + q·y) = max (2·deg p) (2·deg q + 3)

★単元 `u` のノルムは `F[X]` の単元だから次数 0。
★★右辺の `2·deg q + 3` は `q ≠ 0` なら 3 以上なので、**`q = 0`**。
★★★残った `2·deg p = 0` から `deg p = 0`、すなわち `p` は定数である。

★★★★**幾何的には**「アフィン曲線上いたるところ正則で零点も持たない関数は定数」
——`x`・`y` は無限遠に極を持つので単元になれない、という事実の代数版である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isUnit_coordinateRing` | ★★★★★**単元は定数** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F] {W : WeierstrassCurve.Affine F}

/-- ★★★★★**座標環の単元は定数である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `degree_norm_smul_basis`(ノルムの次数)から出る
——`y` の係数が 0 になり、`x` の係数が定数になる。 -/
theorem isUnit_coordinateRing {u : W.CoordinateRing} (hu : IsUnit u) :
    ∃ c : F, c ≠ 0 ∧ u = algebraMap F W.CoordinateRing c := by
  obtain ⟨p, q, hpq⟩ := CoordinateRing.exists_smul_basis_eq u
  have hnorm : IsUnit (Algebra.norm (Polynomial F) u) := hu.map _
  have hdeg : (Algebra.norm (Polynomial F) u).degree = 0 :=
    Polynomial.degree_eq_zero_of_isUnit hnorm
  rw [← hpq, CoordinateRing.degree_norm_smul_basis] at hdeg
  have hq : q = 0 := by
    by_contra hq0
    have hle : 2 • q.degree + 3 ≤ (0 : WithBot ℕ) := by
      rw [← hdeg]
      exact le_max_right _ _
    rw [Polynomial.degree_eq_natDegree hq0, two_nsmul] at hle
    have hcast : ((q.natDegree + q.natDegree + 3 : ℕ) : WithBot ℕ) ≤ ((0 : ℕ) : WithBot ℕ) := by
      push_cast
      exact hle
    have := Nat.cast_le.1 hcast
    omega
  subst hq
  simp only [Polynomial.degree_zero, two_nsmul] at hdeg
  have hp0 : p ≠ 0 := by
    intro h
    rw [h, Polynomial.degree_zero] at hdeg
    simp at hdeg
  rw [Polynomial.degree_eq_natDegree hp0] at hdeg
  have hcast : ((p.natDegree + p.natDegree : ℕ) : WithBot ℕ) = ((0 : ℕ) : WithBot ℕ) := by
    push_cast
    exact hdeg
  have hd : p.natDegree + p.natDegree = 0 := Nat.cast_injective hcast
  obtain ⟨c, hc⟩ : ∃ c : F, p = Polynomial.C c :=
    ⟨p.coeff 0, Polynomial.eq_C_of_natDegree_eq_zero (by omega)⟩
  refine ⟨c, ?_, ?_⟩
  · intro hc0
    exact hp0 (by rw [hc, hc0, map_zero])
  · rw [← hpq, zero_smul, add_zero, hc, CoordinateRing.smul, mul_one]
    rfl

/-! ## ★出典の紐付け(`.src`) -/

def isUnit_coordinateRing.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——座標環の単元が定数であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
