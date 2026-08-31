import ABC3.Check.GenEll.Prop17AxiomGap
import ABC3.Skeleton.GenEll.Section2

/-!
# ★★★[GenEll] `Theorem 2.1` も現在の `Interface` の下では**偽**である

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.11。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★病は §1 だけではなかった

`RemarkAxiomGap.lean` と `Prop17AxiomGap.lean` で §1 の 4 件が偽であることを示した。
★**`Theorem 2.1` にも同じ病がある**。

## ★★理由 —— 今度は `Reduced` が空になれる

`AbcSetup` の `Reduced : (X : Curve) → (data X).Divisor → Prop` は
**公理を持たない述語のフィールド**である。

★`Reduced := fun _ _ => False` と置けば、**(i) は空虚に真**になる。
一方 (ii) は `Reduced` を要求しないので、`logDiff` を点に線形に依らせれば**偽**にできる。
★★ゆえに同値(`↔`)は成り立たない。

★★★これは「不透明な述語は制約にならない」の**もう 1 つの形**である
——`Prop17AxiomGap.lean` では `hyp := True` で仮定を**無力化**したが、
ここでは `Reduced := False` で片側を**空虚化**する。
★★★★**posit した述語は、強くも弱くも取れてしまう。**

## ★★★★★では何をすべきか

`Remark 1.5.1` と同じ治療である——**構成に載せ替える**。
`Theorem 2.1` の場合、載せ替え先は `ArcModel` / `ArithCartier` / `htArith` の側であり、
`Reduced` は `ADivRed` の不動点(`Conductor.lean` の `adivRed_idem`)として**構成できる**。
-/

namespace ABC3.Check.GenEll

open ABC3.Interface.GenEll ABC3.Skeleton.GenEll ABC3.Found.GenEll

/-! ## ★`Reduced` が空な setup -/

/-- ★**`Reduced := False` の `AbcSetup`** —— (i) を空虚にし、(ii) を破る。

★`data` は `linDataDiff`(`logDiff x = x`、`ht L x = L·x`)。
`omegaOf := 0` なので (ii) の左辺は `0`、右辺は `(1+ε)·x` で発散する。 -/
noncomputable def badAbc : AbcSetup where
  Curve := Unit
  data := fun _ => linDataDiff
  Reduced := fun _ _ => False
  omegaOf := fun _ _ => (0 : ℝ)
  Hyperbolic := fun _ _ => True
  projLine := ()
  threePoints := ()
  SupportContains := fun _ _ => True

/-! ## ★★★反例 -/

/-- ★★★★**`Theorem 2.1` の statement は任意の `AbcSetup` については偽**。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

`badAbc` では `Reduced` が `False` なので (i) は空虚に真だが、
(ii) は `KV = ⊤`、`d = 1`、`ε = 1` で `2·x ≤ C` を要求するので偽である。

★★**この `sorry` も「まだ証明していない」ではなく「証明できない」である**。 -/
theorem theorem_2_1_false :
    ¬ (∀ (A : AbcSetup) (primes : Finset ℕ), (∀ p ∈ primes, p.Prime) →
        ((∀ (X : A.Curve) (dv : (A.data X).Divisor) (d : ℕ) (eps : ℝ),
            A.Reduced X dv → 0 < d → 0 < eps → A.Hyperbolic X dv →
            BDle (fun x : ↥((A.data X).compl dv ∩ (A.data X).degLe d) =>
                    (A.data X).ht (A.omegaOf X dv) x.1)
                 (fun x : ↥((A.data X).compl dv ∩ (A.data X).degLe d) =>
                    (1 + eps) * ((A.data X).logDiff x.1 + (A.data X).logCond dv x.1)))
        ↔
          (∀ (KV : Set (A.data A.projLine).Point) (d : ℕ) (eps : ℝ),
            (A.data A.projLine).CompactlyBounded KV → A.SupportContains KV primes →
            0 < d → 0 < eps →
            BDle (fun x : ↥(KV ∩ (A.data A.projLine).compl A.threePoints
                                ∩ (A.data A.projLine).degLe d) =>
                    (A.data A.projLine).ht (A.omegaOf A.projLine A.threePoints) x.1)
                 (fun x : ↥(KV ∩ (A.data A.projLine).compl A.threePoints
                                ∩ (A.data A.projLine).degLe d) =>
                    (1 + eps) * ((A.data A.projLine).logDiff x.1
                                  + (A.data A.projLine).logCond A.threePoints x.1))))) := by
  intro h
  have hiff := h badAbc ∅ (by simp)
  have hi : ∀ (X : badAbc.Curve) (dv : (badAbc.data X).Divisor) (d : ℕ) (eps : ℝ),
      badAbc.Reduced X dv → 0 < d → 0 < eps → badAbc.Hyperbolic X dv →
      BDle (fun x : ↥((badAbc.data X).compl dv ∩ (badAbc.data X).degLe d) =>
              (badAbc.data X).ht (badAbc.omegaOf X dv) x.1)
           (fun x : ↥((badAbc.data X).compl dv ∩ (badAbc.data X).degLe d) =>
              (1 + eps) * ((badAbc.data X).logDiff x.1
                            + (badAbc.data X).logCond dv x.1)) :=
    fun _ _ _ _ hred _ _ _ => hred.elim
  obtain ⟨C, hC⟩ := hiff.1 hi Set.univ 1 1 trivial trivial Nat.one_pos one_pos
  obtain ⟨n, hn⟩ := exists_nat_gt C
  have hmem : (n : ℕ) ∈ (Set.univ : Set ℕ) ∩ (badAbc.data badAbc.projLine).compl
      badAbc.threePoints ∩ (badAbc.data badAbc.projLine).degLe 1 := by
    exact ⟨⟨Set.mem_univ n, Set.mem_univ n⟩, Set.mem_univ n⟩
  have hb := hC ⟨n, hmem⟩
  simp only [badAbc, linDataDiff, zero_mul, sub_zero] at hb
  linarith

end ABC3.Check.GenEll
