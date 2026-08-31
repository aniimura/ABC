/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayGenerated
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★段 C2a-2 —— `A⁰_{x_i} ≅ R[x]/(x_i − 1)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★これは何か —— 標準チャートは `𝔸ᴺ` である

`§9-850`（`AwayGenerated.lean`）は「`A⁰_{x_i}` からの環準同型は定数と `x_j/x_i` での値で決まる」
＝**生成**の半分を取った。★段 C2c-1 の逆向きには**自由性**（関係が無いこと）が要る。

★★本ファイルはその両方を一度に取る:

> **`ker (awayEval) = (x_i − 1)`**、すなわち `A⁰_{x_i} ≅ R[x_0,…,x_N]/(x_i − 1)`。

★★★右辺は `x_i` を `1` と置いた多項式環、つまり `R[x_j/x_i : j ≠ i]` である
——**`D₊(x_i) ≅ 𝔸ᴺ`** の代数側の言い方である。

## ★★★★機構 —— 両向きの写像を作って合成を見る

| 向き | 写像 | 中身 |
|---|---|---|
| `→` | `awayEval` | `x_j ↦ x_j/x_i`（`eval₂Hom`） |
| `←` | `awayQuotHom` | `a/x_i^k ↦ [a]`（分母は `[x_i] = 1` で消える） |

★`awayEval` が**全射**であることは `awayEval_mk`（`a` が `k` 次斉次なら `awayEval a = a/x_i^k`）から出る
——`§9-850` の単項式分解をそのまま使う。
★★`awayQuotHom ∘ awayEval = [·]` は `MvPolynomial.ringHom_ext`（`C` と `X` で見る）で出る。
★★★核が `(x_i − 1)` に一致するのはその 2 本からである。

## ★測定の記録

★`Proj` の標準チャートが `𝔸ᴺ` であることは mathlib に**無い**（2026-08-28 実測）。
★★等式はすべて `HomogeneousLocalization.val_injective` で `Localization.mk` の計算に落ちる。
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★定数を `A⁰_{x_i}` へ送る環準同型 -/

/-- ★**定数を `A⁰_{x_i}` へ送る環準同型**（`§9-838` の `awayConst` を束ねた形）。 -/
noncomputable def awayConstHom (N : ℕ) (R : Type) [CommRing R] [Nontrivial R]
    (i : Fin (N + 1)) :
    R →+* HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i) where
  toFun := awayConst N R i
  map_one' := by
    refine HomogeneousLocalization.val_injective _ ?_
    rw [awayConst, Away.val_mk, HomogeneousLocalization.val_one]
    show Localization.mk (MvPolynomial.C (1 : R)) _ = 1
    rw [map_one]
    rw [show (⟨(MvPolynomial.X i : MvPolynomial (Fin (N+1)) R) ^ 0, ⟨0, rfl⟩⟩ :
      Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N+1)) R)) = 1 from by ext; simp]
    exact Localization.mk_one
  map_zero' := by
    refine HomogeneousLocalization.val_injective _ ?_
    rw [awayConst, Away.val_mk, HomogeneousLocalization.val_zero]
    show Localization.mk (MvPolynomial.C (0 : R)) _ = 0
    rw [map_zero, Localization.mk_zero]
  map_mul' := fun a b => by
    refine HomogeneousLocalization.val_injective _ ?_
    rw [awayConst, Away.val_mk, HomogeneousLocalization.val_mul, awayConst, awayConst,
      Away.val_mk, Away.val_mk, Localization.mk_mul, map_mul]
    congr 1
    ext
    simp
  map_add' := fun a b => by
    refine HomogeneousLocalization.val_injective _ ?_
    rw [awayConst, Away.val_mk, HomogeneousLocalization.val_add, awayConst, awayConst,
      Away.val_mk, Away.val_mk, map_add, Localization.add_mk_self]

/-! ## ★★★多項式を `A⁰_{x_i}` へ評価する -/

/-- ★★★**多項式を `A⁰_{x_i}` へ評価する**（`x_j ↦ x_j/x_i`）。 -/
noncomputable def awayEval (N : ℕ) (R : Type) [CommRing R] [Nontrivial R] (i : Fin (N + 1)) :
    MvPolynomial (Fin (N + 1)) R →+* HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i) :=
  MvPolynomial.eval₂Hom (awayConstHom N R i) (fun j => projCoord N R i j)

open scoped Classical in
/-- ★★★★★**斉次式の評価はちょうど `a/x_i^k` である** —— これが全射性を与える。

★`§9-850` の単項式分解（`away_mk_monomial` ＋ `awayMkAddHom`）をそのまま使う。 -/
theorem awayEval_mk (N : ℕ) (R : Type) [CommRing R] [Nontrivial R] (i : Fin (N + 1))
    (k : ℕ) (a : MvPolynomial (Fin (N+1)) R)
    (ha : a ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1)) :
    awayEval N R i a
      = Away.mk (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R)
          (MvPolynomial.isHomogeneous_X R i) k a ha := by
  have hdeg : ∀ m ∈ a.support, (∑ j ∈ m.support, m j) = k := by
    intro m hm
    rw [MvPolynomial.mem_homogeneousSubmodule] at ha
    have h := ha (MvPolynomial.mem_support_iff.1 hm)
    rw [Finsupp.weight_apply] at h
    simpa [Finsupp.sum] using h
  have hmono : ∀ m ∈ a.support,
      (MvPolynomial.monomial m (MvPolynomial.coeff m a) : MvPolynomial (Fin (N+1)) R)
        ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1) := by
    intro m hm
    rw [MvPolynomial.mem_homogeneousSubmodule]
    have h := MvPolynomial.isHomogeneous_monomial (d := m) (MvPolynomial.coeff m a)
      (n := k) (hdeg m hm)
    simpa using h
  have h1 : (⟨a, ha⟩ : MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1))
      = ∑ m ∈ a.support.attach,
        ⟨MvPolynomial.monomial m.1 (MvPolynomial.coeff m.1 a), hmono m.1 m.2⟩ := by
    apply Subtype.ext
    rw [Submodule.coe_sum]
    show a = ∑ m ∈ a.support.attach, MvPolynomial.monomial m.1 (MvPolynomial.coeff m.1 a)
    rw [Finset.sum_attach a.support (fun m => MvPolynomial.monomial m (MvPolynomial.coeff m a))]
    exact MvPolynomial.as_sum a
  have hR : Away.mk (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R)
        (MvPolynomial.isHomogeneous_X R i) k a ha
      = ∑ m ∈ a.support.attach, (awayConst N R i (MvPolynomial.coeff m.1 a)
          * ∏ j ∈ (m.1).support, (projCoord N R i j) ^ (m.1 j)) := by
    show awayMkAddHom R i k ⟨a, ha⟩ = _
    rw [h1, map_sum]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    exact away_mk_monomial N R i m.1 (MvPolynomial.coeff m.1 a) k (hdeg m.1 m.2) (hmono m.1 m.2)
  rw [hR]
  conv_lhs => rw [MvPolynomial.as_sum a,
    ← Finset.sum_attach a.support (fun m => MvPolynomial.monomial m (MvPolynomial.coeff m a))]
  rw [map_sum]
  refine Finset.sum_congr rfl (fun m _ => ?_)
  rw [awayEval, MvPolynomial.eval₂Hom_monomial]
  rfl

/-- ★★**`awayEval` は全射である**。 -/
theorem awayEval_surjective (N : ℕ) (R : Type) [CommRing R] [Nontrivial R] (i : Fin (N + 1)) :
    Function.Surjective (awayEval N R i) := by
  intro z
  obtain ⟨k, a, ha, rfl⟩ := Away.mk_surjective
    (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.isHomogeneous_X R i) z
  exact ⟨a, awayEval_mk N R i k a ha⟩

/-! ## ★★★逆向き —— `a/x_i^k ↦ [a]` -/

theorem quot_X_eq_one (N : ℕ) (R : Type) [CommRing R] (i : Fin (N + 1)) :
    (Ideal.Quotient.mk (Ideal.span {(MvPolynomial.X i - 1 :
      MvPolynomial (Fin (N+1)) R)})) (MvPolynomial.X i) = 1 := by
  have h : (MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) R)
      ∈ Ideal.span {(MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) R)} :=
    Ideal.subset_span rfl
  have h2 := (Ideal.Quotient.eq_zero_iff_mem).2 h
  rw [map_sub, map_one, sub_eq_zero] at h2
  exact h2

/-- ★★**`a/x_i^k ↦ [a]`**（分母は `[x_i] = 1` で消える）。 -/
noncomputable def awayQuot (N : ℕ) (R : Type) [CommRing R] [Nontrivial R] (i : Fin (N + 1))
    (z : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i)) :
    MvPolynomial (Fin (N+1)) R ⧸ Ideal.span {(MvPolynomial.X i - 1 :
      MvPolynomial (Fin (N+1)) R)} :=
  Quotient.liftOn' z
    (fun c => Ideal.Quotient.mk _ (c.num : MvPolynomial (Fin (N+1)) R))
    (fun c c' h => by
      obtain ⟨j, hj⟩ := exists_mul_eq_of_embedding_eq i c c' h
      have h2 := congrArg (Ideal.Quotient.mk (Ideal.span
        {(MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) R)})) hj
      simpa [map_mul, map_pow, quot_X_eq_one] using h2)

theorem quot_den (N : ℕ) (R : Type) [CommRing R] [Nontrivial R] (i : Fin (N + 1))
    (c : NumDenSameDeg (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R))) :
    (Ideal.Quotient.mk (Ideal.span {(MvPolynomial.X i - 1 :
      MvPolynomial (Fin (N+1)) R)})) (c.den : MvPolynomial (Fin (N+1)) R) = 1 := by
  obtain ⟨k, hk, -⟩ := exists_pow_of_numDenSameDeg i c
  rw [hk, map_pow, quot_X_eq_one, one_pow]

/-- ★★★**逆向きの環準同型** `A⁰_{x_i} →+* R[x]/(x_i − 1)`。 -/
noncomputable def awayQuotHom (N : ℕ) (R : Type) [CommRing R] [Nontrivial R]
    (i : Fin (N + 1)) :
    HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i)
      →+* (MvPolynomial (Fin (N+1)) R ⧸ Ideal.span {(MvPolynomial.X i - 1 :
        MvPolynomial (Fin (N+1)) R)}) where
  toFun := awayQuot N R i
  map_one' := by
    show (Ideal.Quotient.mk _) ((NumDenSameDeg.num (1 : NumDenSameDeg
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R))) :
        MvPolynomial (Fin (N+1)) R)) = 1
    rw [NumDenSameDeg.num_one]
    exact map_one _
  map_zero' := by
    show (Ideal.Quotient.mk _) ((NumDenSameDeg.num (0 : NumDenSameDeg
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R))) :
        MvPolynomial (Fin (N+1)) R)) = 0
    rw [NumDenSameDeg.num_zero, ZeroMemClass.coe_zero]
    exact map_zero _
  map_mul' := by
    rintro ⟨c⟩ ⟨c'⟩
    show (Ideal.Quotient.mk _) (((c * c').num : MvPolynomial (Fin (N+1)) R)) = _
    rw [NumDenSameDeg.num_mul, map_mul]
    rfl
  map_add' := by
    rintro ⟨c⟩ ⟨c'⟩
    show (Ideal.Quotient.mk _) (((c + c').num : MvPolynomial (Fin (N+1)) R)) = _
    rw [NumDenSameDeg.num_add, map_add, map_mul, map_mul,
      quot_den N R i c, quot_den N R i c', one_mul, one_mul]
    show _ = (Ideal.Quotient.mk _) (c.num : MvPolynomial (Fin (N+1)) R)
      + (Ideal.Quotient.mk _) (c'.num : MvPolynomial (Fin (N+1)) R)
    ring

/-! ## ★★★★★★★★★段 C2a-2 の本体 -/

/-- ★★★★★**合成は商写像である** —— `MvPolynomial.ringHom_ext` で `C` と `X` を見るだけ。 -/
theorem awayQuotHom_comp_awayEval (N : ℕ) (R : Type) [CommRing R] [Nontrivial R]
    (i : Fin (N + 1)) :
    (awayQuotHom N R i).comp (awayEval N R i)
      = Ideal.Quotient.mk (Ideal.span {(MvPolynomial.X i - 1 :
          MvPolynomial (Fin (N+1)) R)}) := by
  refine MvPolynomial.ringHom_ext (fun r => ?_) (fun j => ?_)
  · show awayQuot N R i (awayEval N R i (MvPolynomial.C r)) = _
    rw [awayEval, MvPolynomial.eval₂Hom_C]
    rfl
  · show awayQuot N R i (awayEval N R i (MvPolynomial.X j)) = _
    rw [awayEval, MvPolynomial.eval₂Hom_X']
    rfl

/-- ★★★★★★★★★**段 C2a-2** —— `ker (awayEval) = (x_i − 1)`。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`awayEval` は全射（`awayEval_surjective`）なので、これは

    `A⁰_{x_i} ≅ R[x_0,…,x_N]/(x_i − 1) ≅ R[x_j/x_i : j ≠ i]`

すなわち **`D₊(x_i) ≅ 𝔸ᴺ`** の代数側の言い方である。

★★これが段 C2c-1 の逆向き（`ker ≤ span {x_0/x_i}`）に要る**自由性**である
——`§9-850` は生成の半分しか取っていなかった。 -/
theorem ker_awayEval (N : ℕ) (R : Type) [CommRing R] [Nontrivial R] (i : Fin (N + 1)) :
    RingHom.ker (awayEval N R i)
      = Ideal.span {(MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) R)} := by
  apply le_antisymm
  · intro p hp
    have h : (awayQuotHom N R i) (awayEval N R i p) = 0 := by
      rw [RingHom.mem_ker] at hp
      rw [hp, map_zero]
    rw [← RingHom.comp_apply, awayQuotHom_comp_awayEval] at h
    exact (Ideal.Quotient.eq_zero_iff_mem).1 h
  · rw [Ideal.span_le]
    rintro z rfl
    rw [SetLike.mem_coe, RingHom.mem_ker, map_sub, map_one]
    rw [show (awayEval N R i) (MvPolynomial.X i) = projCoord N R i i from by
      rw [awayEval, MvPolynomial.eval₂Hom_X']]
    rw [projCoord_self, sub_self]

/-! ## ★出典の紐付け(`.src`) -/

def awayEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(多項式を A⁰_{x_i} へ評価する)",
    sectionId := "genell-prop-1-4" }

def awayEval_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(awayEval は全射である)",
    sectionId := "genell-prop-1-4" }

def awayQuotHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(逆向きの環準同型 A⁰_{x_i} → R[x]/(x_i − 1))",
    sectionId := "genell-prop-1-4" }

def ker_awayEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(段 C2a-2——A⁰_{x_i} ≅ R[x]/(x_i − 1))",
    sectionId := "genell-prop-1-4" }

def ker_awayEval.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "away_mk_monomial / awayMkAddHom(単項式分解、§9-850)"
      (.inProject "ABC3" "ABC3.Found.GenEll.away_mk_monomial") 2,
    .citation "[ABC3]" "exists_mul_eq_of_embedding_eq(§9-814)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_mul_eq_of_embedding_eq") 2,
    .citation "[mathlib]" "MvPolynomial.ringHom_ext(C と X で見る)"
      (.inMathlib "MvPolynomial.ringHom_ext") 2,
    .implicitStep
      ("★Proj の標準チャートが 𝔸ᴺ であることは mathlib に**無い**(2026-08-28 実測)。" ++
       "★★等式はすべて HomogeneousLocalization.val_injective で " ++
       "Localization.mk の計算に落ちる") 3,
    .implicitStep
      ("★★★これが段 C2c-1 の逆向き(ker ≤ span {x_0/x_i})に要る**自由性**である" ++
       "——§9-850 は生成の半分しか取っていなかった") 5 ]

end ABC3.Found.GenEll
