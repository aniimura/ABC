import ABC3.Found.GenEll.DenominatorBound
import Mathlib.NumberTheory.Height.Projectivization

/-!
# [GenEll] Example 1.3, (i) / Proposition 1.4, (iv) —— 射影空間の Northcott 性(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★mathlib が **TODO** と書いている場所である

`Mathlib/NumberTheory/Height/Northcott.lean` の冒頭:

> We provide instances showing that `K` also satisfies the Northcott property
> * for `logHeight₁`,
> * **(TODO)** for `Projectivization.mulHeight`,
> * **(TODO)** for `Projectivization.logHeight`.
>
> ## TODO
> Add instances for heights on projectivizations.

★`Found/GenEll/DenominatorBound.lean` で `mulHeight₁` の基底を取ったので、
**そこから射影空間の場合が出る**。

## ★機構 —— 座標比に落とす

`P = [x_0 : … : x_n]` に対し `x_j ≠ 0` を 1 つ選び `y_i ≝ x_i / x_j` と置く。
`P = [y_0 : … : y_n]` であり、各座標について

> `H(y_i) = H(x_i / x_j) = H([x_i : x_j]) ≤ H([x_0 : … : x_n]) = H(P) ≤ B`

★2 つ目の等号は `Height.mulHeight₁_div_eq_mulHeight`、
不等号は **`Height.mulHeight_comp_le`**(添字を減らすと高さは減らない)——
**どちらも mathlib にある**。

したがって全座標が有限集合 `{z | H(z) ≤ B}`(こちらは我々が取った)に入り、
`ι` が有限なので組は有限個である。

## ★★これが `Proposition 1.4, (iv)` に対して持つ意味

原文の `ht_L` は `X ⊆ ℙⁿ` に対して `O(1)` の高さの制限であり、
**射影空間の Northcott 性がその基底**である。
★★ただし原文は `X(ℚ̄)^{≤d}`(**次数有界の代数点**)を走るのに対し、
本定理は**固定した数体 `K`** 上の点である。
**基底が取れたことと、全体が取れたことを混同しない。**
-/

namespace ABC3.Found.GenEll

open NumberField

variable {K : Type*} [Field K] [NumberField K] {ι : Type*} [Finite ι]

/-- ★★**射影空間の Northcott 性** —— 高さが `B` 以下の `ℙ(K^ι)` の点は有限個。

★mathlib が TODO と書いている場所である(上の docstring)。 -/
theorem finite_projectivization_mulHeight_le (B : ℝ) :
    {P : Projectivization K (ι → K) | Projectivization.mulHeight P ≤ B}.Finite := by
  classical
  have hS : {z : K | Height.mulHeight₁ z ≤ B}.Finite := finite_mulHeight₁_le B
  have hT : {y : ι → K | ∀ i, Height.mulHeight₁ (y i) ≤ B}.Finite := by
    have heq : {y : ι → K | ∀ i, Height.mulHeight₁ (y i) ≤ B}
        = Set.pi Set.univ (fun _ : ι => {z : K | Height.mulHeight₁ z ≤ B}) := by
      ext y; simp
    rw [heq]
    exact Set.Finite.pi (fun _ => hS)
  have hT' : (Subtype.val ⁻¹' {y : ι → K | ∀ i, Height.mulHeight₁ (y i) ≤ B} :
      Set {y : ι → K // y ≠ 0}).Finite :=
    Set.Finite.preimage (Subtype.val_injective.injOn) hT
  refine Set.Finite.subset
    (hT'.image (fun y : {y : ι → K // y ≠ 0} => Projectivization.mk K y.1 y.2)) ?_
  intro P hP
  simp only [Set.mem_setOf_eq] at hP
  -- 代表元を取る
  set x : ι → K := P.rep with hxdef
  have hx : x ≠ 0 := P.rep_nonzero
  have hPmk : Projectivization.mk K x hx = P := P.mk_rep
  have hPh : Height.mulHeight x ≤ B := by
    rw [← Projectivization.mulHeight_mk hx, hPmk]
    exact hP
  -- 非零な座標を 1 つ選ぶ
  obtain ⟨j, hj⟩ : ∃ j, x j ≠ 0 := by
    by_contra hc
    push_neg at hc
    exact hx (funext hc)
  set y : ι → K := fun i => x i / x j with hy
  have hyj : y j = 1 := by simp [hy, div_self hj]
  have hyne : y ≠ 0 := by
    intro hc
    have : y j = 0 := by rw [hc]; rfl
    rw [hyj] at this
    exact one_ne_zero this
  refine ⟨⟨y, hyne⟩, ?_, ?_⟩
  · -- 各座標の高さが `B` 以下
    intro i
    have hcomp : (x ∘ ![i, j]) = ![x i, x j] := by
      funext k
      fin_cases k <;> rfl
    calc Height.mulHeight₁ (y i) = Height.mulHeight ![x i, x j] := by
          rw [hy]; exact Height.mulHeight₁_div_eq_mulHeight _ _
      _ = Height.mulHeight (x ∘ ![i, j]) := by rw [hcomp]
      _ ≤ Height.mulHeight x := Height.mulHeight_comp_le _ _
      _ ≤ B := hPh
  · -- `[y] = [x] = P`
    rw [← hPmk]
    refine (Projectivization.mk_eq_mk_iff K y x hyne hx).2 ⟨Units.mk0 (x j)⁻¹ (by simpa using hj), ?_⟩
    funext i
    simp [hy, div_eq_inv_mul]

/-- ★**対数版**も同時に出る。 -/
theorem finite_projectivization_logHeight_le (B : ℝ) :
    {P : Projectivization K (ι → K) | Projectivization.logHeight P ≤ B}.Finite := by
  refine Set.Finite.subset (finite_projectivization_mulHeight_le (K := K) (ι := ι)
    (Real.exp B)) ?_
  intro P hP
  simp only [Set.mem_setOf_eq] at hP ⊢
  rw [Projectivization.logHeight_eq_log_mulHeight] at hP
  exact (Real.log_le_iff_le_exp (Projectivization.mulHeight_pos P)).1 hP

/-- ★★★**mathlib の TODO その 1** —— `Northcott (Projectivization.mulHeight)`。 -/
instance northcott_projectivization_mulHeight :
    Northcott (Projectivization.mulHeight (K := K) (ι := ι)) where
  finite_le B := finite_projectivization_mulHeight_le B

/-- ★★★**mathlib の TODO その 2** —— `Northcott (Projectivization.logHeight)`。 -/
instance northcott_projectivization_logHeight :
    Northcott (Projectivization.logHeight (K := K) (ι := ι)) where
  finite_le B := finite_projectivization_logHeight_le B

/-! ## ★出典の紐付け(`.src`) -/

def finite_projectivization_mulHeight_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(固定した数体の射影空間の場合のみ)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
