/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateLimit
import ABC3.Meta.Claim

/-!
# `T_l` の第 `n` 射影の核は `l^n · T_l`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か

`galRep` の連続性（`Theorem 3.8` に残る葉 5）の**第 3 段**の核である:

> `σ` が `E[l^n]` に自明に作用するなら `galMat σ ≡ 1 (mod l^n)`

を出すには、`T_l → E[l^n]` の核が `l^n · T_l` であることが要る。

## ★★★証明（`z_n = 0` から `z = l^n · w` を作る）

`w_k ≔ z_{k+n}` と置く。

* `w_k` は `E[l^k]` に入る —— `l^k · z_{k+n} = z_n = 0`（塔の関係を `k` 回）
* `w` は塔の関係を満たす —— `l · w_{k+1} = l · z_{k+1+n} = z_{k+n} = w_k`
* `l^n · w = z` —— `l^n · z_{k+n} = z_k`（塔の関係を `n` 回）

★★逆向き（`l^n · w` の第 `n` 成分は `0`）は `w_n ∈ E[l^n]` から直ちに従う。
-/

namespace ABC3.Found.GaloisRep

universe u

variable {A : Type u} [AddCommGroup A]

/-- ★★塔の関係を `j` 回使った形——`l^j · z_{k+j} = z_k`。 -/
theorem limTors_shift (l : ℕ) (z : limTors A l) (k : ℕ) :
    ∀ j : ℕ, (l ^ j) • ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + j) : A)
      = ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) k : A) := by
  intro j
  induction j with
  | zero => simp
  | succ i ih =>
    have hstep : l • ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + i + 1) : A)
        = ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + i) : A) := z.2 (k + i)
    calc (l ^ (i + 1)) • ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + (i + 1)) : A)
        = (l ^ i) • (l • ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + i + 1) : A)) := by
          rw [smul_smul, ← pow_succ]
          rfl
      _ = (l ^ i) • ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + i) : A) := by rw [hstep]
      _ = ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) k : A) := ih

/-- ★★★★★★★★★★**第 `n` 射影の核は `l^n · T_l`**（片方向）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`z_n = 0` なら `w_k ≔ z_{k+n}` が `z = l^n · w` を与える。 -/
theorem exists_smul_of_proj_zero (l n : ℕ) (z : limTors A l)
    (hz : ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) n : A) = 0) :
    ∃ w : limTors A l, (l ^ n) • w = z := by
  have hmem : ∀ k : ℕ, ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + n) : A)
      ∈ (nsmulHom A (l ^ k)).ker := by
    intro k
    rw [AddMonoidHom.mem_ker]
    show (l ^ k) • ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + n) : A) = 0
    have h := limTors_shift l z n k
    rw [show n + k = k + n from Nat.add_comm n k] at h
    rw [h, hz]
  refine ⟨⟨fun k => ⟨((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + n) : A), hmem k⟩, ?_⟩, ?_⟩
  · intro k
    show l • ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + 1 + n) : A)
      = ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + n) : A)
    have := z.2 (k + n)
    rw [show k + 1 + n = k + n + 1 from by omega]
    exact this
  · refine Subtype.ext (funext fun k => Subtype.ext ?_)
    show (l ^ n) • ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) (k + n) : A)
      = ((z : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) k : A)
    exact limTors_shift l z k n

/-- ★★★★逆向き——`l^n · w` の第 `n` 成分は `0`。 -/
theorem proj_smul_eq_zero (l n : ℕ) (w : limTors A l) :
    (((l ^ n) • w : limTors A l) : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) n
      = (0 : (nsmulHom A (l ^ n)).ker) := by
  refine Subtype.ext ?_
  show (l ^ n) • ((w : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) n : A) = 0
  exact ((w : ∀ m : ℕ, (nsmulHom A (l ^ m)).ker) n).2

/-! ## ★出典の紐付け(`.src`) -/

def exists_smul_of_proj_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(T_l の第 n 射影の核は l^n·T_l。★無条件)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
