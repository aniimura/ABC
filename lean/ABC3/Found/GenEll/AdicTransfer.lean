/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.AdicCompletion.Basic
import Mathlib.RingTheory.Ideal.Operations
import ABC3.Meta.Claim

/-!
# 第 1359 ブロック —— **`I`-進完備性は線型同型で移る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1357 の道の第 3 段

「完備 DVR `R` の有限整閉包 `R′` は完備」を出すには

1. `R′` は `R` 上**有限自由**（`IsIntegralClosure.module_free`、在庫）
2. したがって `R′ ≃ₗ[R] (Fin n → R)`
3. **`IsAdicComplete` は線型同型で移る**（★本ブロック）
4. `Fin n → R` の完備性（次のブロック）

☆本ブロックは 3 を作る。★`I^n • ⊤` は全射線型写像で `I^n • ⊤` に移る
（`Submodule.map_smul''` ＋ `Submodule.map_top`）ので、
`SModEq` がそのまま移る。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

variable {R : Type*} [CommRing R] (I : Ideal R)
variable {M N : Type*} [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]

/-- ★★★★★★**全射線型写像は `I^n • ⊤` を `I^n • ⊤` に移す**（第 1359）。 -/
theorem map_pow_smul_top (e : M ≃ₗ[R] N) (n : ℕ) :
    Submodule.map (e : M →ₗ[R] N) ((I ^ n • ⊤ : Submodule R M)) = (I ^ n • ⊤ : Submodule R N) := by
  rw [Submodule.map_smul'', Submodule.map_top]
  congr 1
  exact LinearMap.range_eq_top.2 e.surjective

/-- ★★★★★★★★**`SModEq` は線型同型で移る**（第 1359）。 -/
theorem sModEq_map (e : M ≃ₗ[R] N) (n : ℕ) {x y : M} :
    x ≡ y [SMOD (I ^ n • ⊤ : Submodule R M)] ↔ e x ≡ e y [SMOD (I ^ n • ⊤ : Submodule R N)] := by
  rw [SModEq.sub_mem, SModEq.sub_mem, ← map_pow_smul_top I e n, ← map_sub]
  constructor
  · intro h
    exact Submodule.mem_map_of_mem h
  · rintro ⟨z, hz, hze⟩
    have : z = x - y := e.injective hze
    rwa [this] at hz

/-- ★★★★★★★★**`IsHausdorff` は線型同型で移る**（第 1359）。 -/
theorem isHausdorff_of_linearEquiv (e : M ≃ₗ[R] N) [IsHausdorff I M] : IsHausdorff I N where
  haus' := by
    intro y hy
    have hx : e.symm y = 0 := by
      refine IsHausdorff.haus' (I := I) (M := M) (e.symm y) ?_
      intro n
      rw [sModEq_map I e n]
      simpa using hy n
    have := congrArg e hx
    simpa using this

/-- ★★★★★★★★**`IsPrecomplete` は線型同型で移る**（第 1359）。 -/
theorem isPrecomplete_of_linearEquiv (e : M ≃ₗ[R] N) [IsPrecomplete I M] :
    IsPrecomplete I N where
  prec' := by
    intro f hf
    have hcauchy : ∀ {m n : ℕ}, m ≤ n →
        e.symm (f m) ≡ e.symm (f n) [SMOD (I ^ m • ⊤ : Submodule R M)] := by
      intro m n hmn
      rw [sModEq_map I e m]
      simpa using hf hmn
    obtain ⟨L, hL⟩ := IsPrecomplete.prec' (I := I) (M := M) (fun n => e.symm (f n)) hcauchy
    refine ⟨e L, fun n => ?_⟩
    have := (sModEq_map I e n).1 (hL n)
    simpa using this

/-- ★★★★★★★★★★★★**`IsAdicComplete` は線型同型で移る**——★**無条件**（第 1359）。 -/
theorem isAdicComplete_of_linearEquiv (e : M ≃ₗ[R] N) [IsAdicComplete I M] :
    IsAdicComplete I N where
  toIsHausdorff := isHausdorff_of_linearEquiv I e
  toIsPrecomplete := isPrecomplete_of_linearEquiv I e

/-! ## ★出典の紐付け(`.src`) -/

def isAdicComplete_of_linearEquiv.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(I-進完備性は線型同型で移る。★無条件)",
    sectionId := "genell-thm-3-8" }

def isAdicComplete_of_linearEquiv.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1359）**——第 1357 の道の**第 3 段**である。" ++
       "☆`R′` は `R` 上有限自由（`IsIntegralClosure.module_free`、在庫）なので、" ++
       "`Fin n → R` の完備性さえ出れば `R′` の完備性が出る。") 19 ]

end ABC3.Found.GenEll
