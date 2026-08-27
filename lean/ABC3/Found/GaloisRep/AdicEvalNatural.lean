/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicEval

/-!
# 冪級数の **adic 値は自然**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★`galois-equivariant-tate-uniformization` の段 2c の核

台帳の 7 段のうち段 2c は「`tatePtPair` の自然性（**級数と `σ` の交換**）」であり、
そこが同変性の本体である。★本ファイルはその**核**を取る:

> **`σ (evalAdic f q) = evalAdic f (σ q)`**（`σ I ⊆ I` なら）

## ★★★★★★またしても「一意 ⟹ 自然」

`evalAdic f q hq` は `Classical.choose` で定義されているが、
**一意性補題 `evalAdic_unique` がある**——
「部分和と `I^n` を法として合同なものは `evalAdic f q hq` に一致する」。

★機構は 3 段:

1. `partialEval` は**多項式**（`∑ coeff·q^n`、係数は `ℤ`）なので `σ` と可換。
2. `σ I ⊆ I` から `σ (I^n) ⊆ I^n`、したがって `SModEq` が移る。
3. `evalAdic_unique` を当てる。

★★★**`normRep`（段 2a）・`tateAOf`（段 2b）と同じ形**である
——ABC3 の Tate 構成は `.choose` だらけだが、
**どれにも一意性補題が付いているので自然性は機械的に出る**。

## ★残っている段（明示）

★`tateXtail` / `tateXterm` / `tateXpairE` の自然性は本補題と環演算の組み合わせで出る。
★★その先の `tatePtPair`（`Point.some` の合同）と、
段 3（有限拡大）・段 4（帰納極限）が残る。
-/

namespace ABC3.Found.GaloisRep

/-! ## ★部分和は環準同型と可換 -/

/-- ★**部分和は多項式なので `σ` と可換**（係数は `ℤ`）。 -/
theorem partialEval_map {R : Type} [CommRing R] (σ : R →+* R)
    (f : PowerSeries ℤ) (q : R) (n : ℕ) :
    σ (partialEval f q n) = partialEval f (σ q) n := by
  unfold partialEval
  rw [map_sum]
  refine Finset.sum_congr rfl (fun k _ => ?_)
  rw [map_mul, map_pow, map_intCast]

/-! ## ★★`SModEq` は `σ I ⊆ I` で移る -/

/-- ★★**`σ I ⊆ I` なら `I^n` を法とする合同は `σ` で保たれる**。

★`Ideal.map σ (I^n) = (Ideal.map σ I)^n ≤ I^n` による。 -/
theorem smodEq_map {R : Type} [CommRing R] {I : Ideal R} (σ : R →+* R)
    (hσI : ∀ x ∈ I, σ x ∈ I) {n : ℕ} {a b : R}
    (h : a ≡ b [SMOD (I ^ n • ⊤ : Submodule R R)]) :
    σ a ≡ σ b [SMOD (I ^ n • ⊤ : Submodule R R)] := by
  rw [SModEq.sub_mem] at h ⊢
  rw [← map_sub]
  have hI : (I ^ n • ⊤ : Submodule R R) = (I ^ n : Ideal R) := by simp
  rw [hI] at h ⊢
  have hmapI : Ideal.map σ I ≤ I := Ideal.map_le_iff_le_comap.2 (fun x hx => hσI x hx)
  have hpow : Ideal.map σ (I ^ n) ≤ I ^ n := by
    rw [Ideal.map_pow]
    exact Ideal.pow_right_mono hmapI n
  exact hpow (Ideal.mem_map_of_mem σ h)

/-! ## ★★★★★★★★adic 値の自然性 -/

/-- ★★★★★★★★**冪級数の adic 値は環準同型で自然**（`σ I ⊆ I` なら）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

> **`evalAdic f (σ q) = σ (evalAdic f q)`**

★★機構は**一意性 1 本**（`evalAdic_unique`）である
——`σ (evalAdic f q)` が部分和 `partialEval f (σ q) n` と `I^n` を法として
合同であることを見れば済む。

★★★`σ q = q` の場合は **`σ (evalAdic f q) = evalAdic f q`**、
すなわち**級数の値は `σ` で固定される**——これが Tate 一意化の同変性の核である。 -/
theorem evalAdic_map {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I)
    (f : PowerSeries ℤ) (q : R) (hq : q ∈ I) (hq' : σ q ∈ I) :
    evalAdic f (σ q) hq' = σ (evalAdic f q hq) := by
  refine evalAdic_unique f (σ q) hq' _ (fun n => ?_)
  rw [← partialEval_map σ f q n]
  exact smodEq_map σ hσI (evalAdic_spec f q hq n)

/-- ★★★★★★★★★**母数を固定する `σ` は級数の値も固定する**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★これが `G_K` が Tate 曲線の係数に自明に作用することの中身である。 -/
theorem evalAdic_fixed {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    (σ : R →+* R) (hσI : ∀ x ∈ I, σ x ∈ I)
    (f : PowerSeries ℤ) (q : R) (hq : q ∈ I) (hσq : σ q = q) :
    σ (evalAdic f q hq) = evalAdic f q hq := by
  have hq' : σ q ∈ I := by rw [hσq]; exact hq
  have h := evalAdic_map σ hσI f q hq hq'
  rw [← h]
  congr 1

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`（Tate 一意化）の**段 2c の核**であって、
`tatePtPair` 全体の自然性ではない。 -/

def evalAdic_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——冪級数の adic 値の自然性)",
    sectionId := "genell-def-3-3" }

def evalAdic_fixed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——母数を固定する σ は級数の値も固定する)",
    sectionId := "genell-def-3-3" }

def evalAdic_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "evalAdic_unique(adic 極限の一意性)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.evalAdic_unique") 15,
    .implicitStep
      ("★★★★またしても「**一意 ⟹ 自然**」である。" ++
       "★ABC3 の Tate 構成は .choose だらけだが、" ++
       "normRep(段 2a)・tateAOf(段 2b)・evalAdic(段 2c)の**どれにも" ++
       "一意性補題が付いている**ので自然性は機械的に出る") 15,
    .implicitStep
      ("★★残る段: tateXtail / tateXterm / tateXpairE の自然性" ++
       "(本補題と環演算の組み合わせ)と tatePtPair(Point.some の合同)、" ++
       "その先の段 3(有限拡大)・段 4(帰納極限)") 15 ]

end ABC3.Found.GaloisRep
