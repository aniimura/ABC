/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.Polynomial.Cyclotomic.Roots
import Mathlib.RingTheory.AdjoinRoot
import ABC3.Meta.Claim

/-!
# 第 1376 ブロック —— **原始 `l` 乗根を足す拡大は次数 `≤ l−1`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`l ∤ e` の出どころ

第 1372（`exists_h2_h1_of_bad_prime_ram`）は `ζ_l ∈ L_v` を要求する。
☆`L_p` に原始 `l` 乗根が無いときは拡大を取るしかないが、
そのとき分岐指数 `e` が `l` で割れると `¬ l ∣ v(q)` が壊れる。

★★★**円分多項式 `Φ_l` の既約因子 `f` を取れば `[L_p(root f) : L_p] = deg f ≤ l−1`**
であり、第 1369（`e ≤ [L′:L]`）と合わせて `1 ≤ e ≤ l−1 < l`、
すなわち **`l ∤ e`** が出る。

☆標数 `0` なので `Φ_l` の根は原始 `l` 乗根である（`isRoot_cyclotomic_iff`）。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta Polynomial

section CycExt

variable (F : Type) [Field F] [CharZero F] {l : ℕ}

omit [CharZero F] in
/-- ★★★★★★**`Φ_l` は既約因子をもつ**（第 1376）。 -/
theorem exists_cyclotomic_irreducible_factor (hl : l.Prime) :
    ∃ f : F[X], Irreducible f ∧ f ∣ cyclotomic l F := by
  refine WfDvdMonoid.exists_irreducible_factor ?_ (cyclotomic_ne_zero l F)
  intro hu
  have hl2 : 2 ≤ l := hl.two_le
  have hdeg : (cyclotomic l F).natDegree = l - 1 := by
    rw [natDegree_cyclotomic, Nat.totient_prime hl]
  have h0 : (cyclotomic l F).natDegree = 0 := natDegree_eq_zero_of_isUnit hu
  rw [hdeg] at h0
  omega

/-- ★`Φ_l` の既約因子（第 1376）。 -/
noncomputable def cycIrrFactor (hl : l.Prime) : F[X] :=
  (exists_cyclotomic_irreducible_factor F hl).choose

theorem cycIrrFactor_irreducible (hl : l.Prime) : Irreducible (cycIrrFactor F hl) :=
  (exists_cyclotomic_irreducible_factor F hl).choose_spec.1

theorem cycIrrFactor_dvd (hl : l.Prime) : cycIrrFactor F hl ∣ cyclotomic l F :=
  (exists_cyclotomic_irreducible_factor F hl).choose_spec.2

instance factCycIrrFactor (hl : l.Prime) : Fact (Irreducible (cycIrrFactor F hl)) :=
  ⟨cycIrrFactor_irreducible F hl⟩

/-- ★★★★★★★★**原始 `l` 乗根を足した局所体**（第 1376）。

☆`abbrev`（reducible）にしてあるので `Field`・`Algebra` などの
インスタンスは `AdjoinRoot` からそのまま降りてくる。 -/
@[reducible] noncomputable def cycExt (hl : l.Prime) : Type :=
  AdjoinRoot (cycIrrFactor F hl)

/-- ★`cycExt` に足した根（第 1376）。 -/
noncomputable def cycRoot (hl : l.Prime) : cycExt F hl :=
  AdjoinRoot.root (cycIrrFactor F hl)

instance instCharZeroCycExt (hl : l.Prime) : CharZero (cycExt F hl) :=
  charZero_of_injective_algebraMap (algebraMap F (cycExt F hl)).injective

/-- ★★★★★★★★★★★★**足した根は原始 `l` 乗根**（第 1376）。 -/
theorem isPrimitiveRoot_cycRoot (hl : l.Prime) :
    IsPrimitiveRoot (cycRoot F hl) l := by
  haveI : NeZero l := ⟨hl.ne_zero⟩
  haveI : NeZero ((l : ℕ) : cycExt F hl) := NeZero.charZero
  obtain ⟨g, hg⟩ := cycIrrFactor_dvd F hl
  have hself : (Polynomial.aeval (cycRoot F hl)) (cycIrrFactor F hl) = 0 := by
    rw [cycRoot, AdjoinRoot.aeval_eq, AdjoinRoot.mk_self]
  have h1 : (Polynomial.aeval (cycRoot F hl)) (cyclotomic l F) = 0 := by
    rw [hg, map_mul, hself, zero_mul]
  have h2 : (Polynomial.aeval (cycRoot F hl)) (cyclotomic l F)
      = Polynomial.eval (cycRoot F hl) (cyclotomic l (cycExt F hl)) := by
    rw [Polynomial.aeval_def, Polynomial.eval₂_eq_eval_map, Polynomial.map_cyclotomic]
  rw [h2] at h1
  exact isRoot_cyclotomic_iff.1 h1

/-- ★★★★★★★★★★★★★★★★
**円分拡大の次数は `≤ l−1`**——★**無条件**（第 1376）。

★★★これが `l ∤ e` の出どころである
——第 1369 が `1 ≤ e ≤ [L′:L]` を与えるので `0 < e < l`。 -/
theorem finrank_cycExt_le (hl : l.Prime) :
    Module.finrank F (cycExt F hl) ≤ l - 1 := by
  have hne : cycIrrFactor F hl ≠ 0 := (cycIrrFactor_irreducible F hl).ne_zero
  have hfr : Module.finrank F (cycExt F hl) = (cycIrrFactor F hl).natDegree := by
    have h := (AdjoinRoot.powerBasis (f := cycIrrFactor F hl) hne).finrank
    rw [AdjoinRoot.powerBasis_dim] at h
    exact h
  rw [hfr]
  have hdvd : (cycIrrFactor F hl).natDegree ≤ (cyclotomic l F).natDegree :=
    Polynomial.natDegree_le_of_dvd (cycIrrFactor_dvd F hl) (cyclotomic_ne_zero l F)
  rw [natDegree_cyclotomic, Nat.totient_prime hl] at hdvd
  exact hdvd

/-- ★★★★★★**円分拡大は有限次元**（第 1376）。 -/
instance instFiniteDimensionalCycExt (hl : l.Prime) :
    FiniteDimensional F (cycExt F hl) :=
  (AdjoinRoot.powerBasis (f := cycIrrFactor F hl)
    (cycIrrFactor_irreducible F hl).ne_zero).finite

/-- ★★★★★★★★**`l ∤ e`**——分岐指数は `l` で割れない（第 1376）。 -/
theorem not_dvd_of_le_finrank (hl : l.Prime) {e : ℕ} (he1 : 1 ≤ e)
    (he : e ≤ Module.finrank F (cycExt F hl)) : ¬ (l ∣ e) := by
  have hle := le_trans he (finrank_cycExt_le F hl)
  have hl2 : 2 ≤ l := hl.two_le
  intro hdvd
  have := Nat.le_of_dvd (by omega) hdvd
  omega

end CycExt

/-! ## ★出典の紐付け(`.src`) -/

def isPrimitiveRoot_cycRoot.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(円分拡大に足した根は原始 l 乗根。★無条件)",
    sectionId := "genell-thm-3-8" }

def finrank_cycExt_le.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(円分拡大の次数は ≤ l−1。★無条件)",
    sectionId := "genell-thm-3-8" }

def not_dvd_of_le_finrank.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(分岐指数は l で割れない。★無条件)",
    sectionId := "genell-thm-3-8" }

def finrank_cycExt_le.needs : List ProofObligation :=
  [ .citation "[ABC3]" "le_finrank_of_pow_eq_map_maximalIdeal(第 1369、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.le_finrank_of_pow_eq_map_maximalIdeal") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1376）**——`l ∤ e` の出どころである。" ++
       "☆`Φ_l` の既約因子 `f` を取れば `[F(root f) : F] = deg f ≤ l−1`。" ++
       "★第 1369 の `1 ≤ e ≤ [L′:L]` と合わせて `0 < e < l`。") 19 ]

end ABC3.Found.GenEll
