import ABC3.Found.FrdI.Def31
import Mathlib.GroupTheory.ResiduallyFinite

/-!
# [FrdI] Remark 3.1.2 —— 剰余有限群への `𝔽` からの準同型は次数を経由する

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.58。

原文 (FrdI p.58):
> such that the intersection of the normal subgroups of finite index of G is trivial]:

原文 (FrdI p.58):
> Any homomorphism of monoids F

## ★中身(測定)

**主張は 1 つ**: 剰余有限群 `G` について、任意のモノイド準同型
`𝔽 → G` は自然な全射 `𝔽 ↠ ℕ≥1` を経由する。

## ★証明の骨

★**急所は `𝔽` の中の 1 本の等式**である:

  `(0, d) · (1, 1) = (1, 1)^d · (0, d)`

(左辺・右辺とも `(d, d)`)。これを `f` で送ると
★**`γ^d` が `γ` の共役**になる —— ここで `γ := f(1,1)`。

★**有限群では共役は位数を保つ**ので、`d` を `γ` の位数に取ると
`γ^d = 1` の共役が `γ` になり、`γ = 1` が出る。
★**剰余有限性は「有限商へ落として `γ = 1`」を全ての有限指数正規部分群で
やる**ことに使う。

★`γ = 1` が出れば `f(a, 1) = γ^a = 1` なので、`f` は `(0, n)` の値だけで決まり、
それがちょうど `ℕ≥1` を経由する分解である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

/-! ## ★`𝔽` の中の計算 -/

/-- ★`(1,1)^d = (d,1)` —— `𝔽` の「純 divisor 部分」は `ℕ` そのもの。 -/
theorem elemFrob_gen_pow (d : ℕ) :
    (⟨1, 1⟩ : ElemFrob ℕ) ^ d = ⟨d, 1⟩ := by
  induction d with
  | zero => rfl
  | succ n ih =>
    rw [pow_succ, ih]
    ext
    · show (n : ℕ) + ((1 : ℕ+) : ℕ) • (1 : ℕ) = (n + 1 : ℕ)
      simp
    · show ((1 : ℕ+) * 1 : ℕ+) = 1
      simp

/-- ★★**急所の等式** —— `(0, d) · (1, 1) = (1, 1)^d · (0, d)`(両辺とも `(d, d)`)。

★これを準同型で送ると「`γ^d` は `γ` の共役」が出る。 -/
theorem elemFrob_conj (d : ℕ+) :
    (⟨0, d⟩ : ElemFrob ℕ) * ⟨1, 1⟩ = (⟨1, 1⟩ : ElemFrob ℕ) ^ (d : ℕ) * ⟨0, d⟩ := by
  rw [elemFrob_gen_pow]
  ext
  · show (0 : ℕ) + (d : ℕ) • (1 : ℕ) = (d : ℕ) + ((1 : ℕ+) : ℕ) • (0 : ℕ)
    simp
  · show (d * 1 : ℕ+) = (1 : ℕ+) * d
    rw [mul_one, one_mul]

/-! ## ★有限群の場合 -/

/-- ★★**有限群では `f(1,1) = 1`** —— 共役が位数を保つことから。

★`d` を `γ` の位数に取ると `γ^d = 1`。上の等式で `γ^d` は `γ` の共役なので
`γ` 自身が `1` の共役、すなわち `γ = 1`。 -/
theorem elemFrob_gen_eq_one_of_finite {G : Type*} [Group G] [Finite G]
    (f : ElemFrob ℕ →* G) : f ⟨1, 1⟩ = 1 := by
  set γ : G := f ⟨1, 1⟩ with hγ
  have hpos : 0 < orderOf γ := orderOf_pos γ
  set d : ℕ+ := ⟨orderOf γ, hpos⟩ with hd
  have hkey := congrArg f (elemFrob_conj d)
  rw [map_mul, map_mul, map_pow] at hkey
  -- `hkey : f ⟨0,d⟩ * γ = γ ^ (d : ℕ) * f ⟨0,d⟩`
  have hone : γ ^ ((d : ℕ+) : ℕ) = 1 := pow_orderOf_eq_one γ
  rw [hone, one_mul] at hkey
  -- `hkey : f ⟨0,d⟩ * γ = f ⟨0,d⟩`
  simpa using hkey

/-! ## ★★「`f(1,1) = 1` なら経由する」は**群を使わない**

★★`rem_3_1_2` の段 2–4 は**モノイド構造だけ**で通る。
★括り出しておくと、**別の理由で `f(1,1) = 1` が言えた場合にそのまま使える**
——`Proposition 3.3, (i)` は「可換かつ簡約的」という別の理由でそれを出す。 -/

/-- ★★★**`f(1,1) = 1` なら `𝔽 → N` は `deg` を経由する**(モノイドで十分)。

★段 2: `f(a,1) = f(1,1)^a = 1`。★段 3: `g n := f(0,n)`。
★段 4: `x = (x.div, 1) · (0, x.deg)` に分けて段 2 を当てる。 -/
theorem elemFrob_factors_of_gen_eq_one {N : Type*} [Monoid N]
    (f : ElemFrob ℕ →* N) (hγ : f ⟨1, 1⟩ = 1) :
    ∃ g : ℕ+ →* N, f = g.comp ElemFrob.degHom := by
  have hdiv : ∀ a : ℕ, f ⟨a, 1⟩ = 1 := by
    intro a
    rw [← elemFrob_gen_pow a, map_pow, hγ, one_pow]
  refine ⟨{ toFun := fun n => f ⟨0, n⟩
            map_one' := by
              show f ⟨0, 1⟩ = 1
              exact f.map_one
            map_mul' := fun n m => by
              have hnm : (⟨0, n⟩ : ElemFrob ℕ) * ⟨0, m⟩ = ⟨0, n * m⟩ := by
                ext
                · show (0 : ℕ) + (n : ℕ) • (0 : ℕ) = 0
                  simp
                · rfl
              show f ⟨0, n * m⟩ = f ⟨0, n⟩ * f ⟨0, m⟩
              rw [← hnm, map_mul] }, ?_⟩
  ext x
  show f x = f ⟨0, x.deg⟩
  have hx : x = (⟨x.div, 1⟩ : ElemFrob ℕ) * ⟨0, x.deg⟩ := by
    ext
    · show x.div = x.div + ((1 : ℕ+) : ℕ) • (0 : ℕ)
      simp
    · show x.deg = (1 : ℕ+) * x.deg
      rw [one_mul]
  conv_lhs => rw [hx]
  rw [map_mul, hdiv, one_mul]

/-- ★★★**可換かつ簡約的なモノイドへの準同型は `deg` を経由する**。

原文 (FrdI p.60):
> obtained by considering the Frobenius degree of the induced endomorphism of A --

★★**原文が `Proposition 3.3, (i)` の証明で引く「`ℕ≥1` は可換」の中身**である。
★`elemFrob_conj` を送ると `f(0,d) · γ = γ^d · f(0,d)`。
可換性で右辺を並べ替え、**簡約**して `γ = γ^d`。★`d = 2` に取れば `γ = γ · γ`、
もう一度簡約して `γ = 1`。 -/
theorem elemFrob_factors_of_cancel {N : Type*} [CommMonoid N] [IsCancelMul N]
    (f : ElemFrob ℕ →* N) : ∃ g : ℕ+ →* N, f = g.comp ElemFrob.degHom := by
  refine elemFrob_factors_of_gen_eq_one f ?_
  set γ : N := f ⟨1, 1⟩ with hγdef
  have hkey := congrArg f (elemFrob_conj 2)
  rw [map_mul, map_mul, map_pow] at hkey
  -- `hkey : f ⟨0,2⟩ * γ = γ ^ 2 * f ⟨0,2⟩`
  have h2 : f (⟨0, 2⟩ : ElemFrob ℕ) * γ = f (⟨0, 2⟩ : ElemFrob ℕ) * γ ^ ((2 : ℕ+) : ℕ) := by
    rw [hkey, mul_comm]
  have hpow : γ = γ ^ ((2 : ℕ+) : ℕ) := mul_left_cancel h2
  have h3 : γ * 1 = γ * γ := by
    rw [mul_one]
    calc γ = γ ^ ((2 : ℕ+) : ℕ) := hpow
      _ = γ * γ := by
        show γ ^ (2 : ℕ) = γ * γ
        rw [pow_two]
  exact (mul_left_cancel h3).symm

/-! ## ★剰余有限群の場合 —— 本体 -/

/-- ★★★**[FrdI] Remark 3.1.2** —— 剰余有限群 `G` について、
任意のモノイド準同型 `𝔽 → G` は自然な全射 `𝔽 ↠ ℕ≥1` を経由する。 -/
theorem rem_3_1_2 {G : Type*} [Group G] [Group.ResiduallyFinite G]
    (f : ElemFrob ℕ →* G) :
    ∃ g : ℕ+ →* G, f = g.comp ElemFrob.degHom := by
  -- 段 1: `γ := f(1,1)` は 1 である
  have hγ : f ⟨1, 1⟩ = 1 := by
    refine Group.eq_one_iff_forall_finiteIndexNormalSubroup _ (fun H => ?_)
    haveI : H.toSubgroup.Normal := H.isNormal'
    haveI : H.toSubgroup.FiniteIndex := H.isFiniteIndex'
    haveI : Finite (G ⧸ H.toSubgroup) := H.toSubgroup.finite_quotient_of_finiteIndex
    have h := elemFrob_gen_eq_one_of_finite
      ((QuotientGroup.mk' H.toSubgroup).comp f)
    rw [MonoidHom.comp_apply] at h
    exact (QuotientGroup.eq_one_iff _).mp h
  -- 段 2–4: ★**括り出した一般補題**(群構造を使わない)
  exact elemFrob_factors_of_gen_eq_one f hγ

/-! ## ★系 —— `Aut(E_A → E)` が剰余有限なら `E` は Frobenius-slim

★**これが原文 `Remark 3.1.2` の使いどころ**である(`Proposition 3.2` が引く)。 -/

universe v u

/-- ★★**`Aut(E_A → E)` がすべて剰余有限なら `E` は Frobenius-slim**。 -/
theorem isFrobeniusSlim_of_residuallyFinite (E : Type u) [Category.{v} E]
    (h : ∀ A : E, Group.ResiduallyFinite (Aut (Over.forget A))) : IsFrobeniusSlim E := by
  intro A f
  haveI := h A
  exact rem_3_1_2 f

/-! ## ★★★出典の紐付け(`.src`) -/

/-- ★★★**[FrdI] Remark 3.1.2** —— 1 主張が実装された。

★系として「`Aut(E_A → E)` が剰余有限なら `E` は Frobenius-slim」も出した
(`isFrobeniusSlim_of_residuallyFinite`)。 -/
def rem_3_1_2.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 58, item := "Remark 3.1.2",
    sectionId := "frdi-remark-3-1-2" }

end ABC3.Found.FrdI
