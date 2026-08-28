/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayEvalGenMk
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★`ker (awayEvalOf f) = (f − 1)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★これは何か

`§9-859`（段 C2a-2）は `f = x_i` について `ker (awayEval) = (x_i − 1)`、すなわち

    `A⁰_{x_i} ≅ R[x]/(x_i − 1)`

を取った。★本ファイルはそれを**任意の 1 次斉次形式 `f`** に広げる:

    `A⁰_f ≅ R[x]/(f − 1)`

★★これが段 C2c-1 の (c) に要る——(a) の四角の**右下の辺**が
`awayEvalOf (g f)` であり、その核を知る必要があるからである。

## ★★★★★測定の記録 —— 一般形の方が証明が短い

★`§9-859` は `awayQuot` の**well-defined 性**に
`exists_pow_of_numDenSameDeg`（分母は `x_i^{deg}` であり指数は次数に等しい、`§9-813`）
を使っていた。これは「斉次な非零元の次数は一意」（`IsHomogeneous.inj_right`）に依存し、
`x_i^k ≠ 0` が要るので `[Nontrivial R]` が要った。

★★しかし**次数は要らなかった**（2026-08-28 実測）。必要なのは

> `R[x]/(f − 1)` の中で `[f] = 1`、したがって `powers f` の元はすべて `[·] = 1`

だけである。分母が `f` の何乗かを**知る必要はない**。
★★★結果として一般形は `[Nontrivial R]` も要らず、`f = 0` でも成り立つ
（両辺とも零環になる）。
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {σ : Type}

/-! ## ★商の中で `powers f` はすべて `1` になる -/

/-- ★**`R[x]/(f − 1)` の中で `[f] = 1`**。 -/
theorem quot_f_eq_one (R : Type) [CommRing R] (f : MvPolynomial σ R) :
    (Ideal.Quotient.mk (Ideal.span {(f - 1 : MvPolynomial σ R)})) f = 1 := by
  have h : (f - 1 : MvPolynomial σ R) ∈ Ideal.span {(f - 1 : MvPolynomial σ R)} :=
    Ideal.subset_span rfl
  have h2 := (Ideal.Quotient.eq_zero_iff_mem).2 h
  rw [map_sub, map_one, sub_eq_zero] at h2
  exact h2

/-- ★★**したがって `powers f` の元はすべて `[·] = 1`** —— 指数を知る必要はない。 -/
theorem quot_mem_powers (R : Type) [CommRing R] (f u : MvPolynomial σ R)
    (hu : u ∈ Submonoid.powers f) :
    (Ideal.Quotient.mk (Ideal.span {(f - 1 : MvPolynomial σ R)})) u = 1 := by
  obtain ⟨k, hk⟩ := hu
  rw [← hk, map_pow, quot_f_eq_one, one_pow]

/-! ## ★★★逆向き —— `a/f^k ↦ [a]` -/

/-- ★★**`a/f^k ↦ [a]`**（分母は `[f] = 1` で消える）。 -/
noncomputable def awayQuotOf (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (z : HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule σ R) f) :
    MvPolynomial σ R ⧸ Ideal.span {(f - 1 : MvPolynomial σ R)} :=
  Quotient.liftOn' z
    (fun c => Ideal.Quotient.mk _ (c.num : MvPolynomial σ R))
    (fun c c' hrel => by
      have h : NumDenSameDeg.embedding _ _ c = NumDenSameDeg.embedding _ _ c' := hrel
      rw [NumDenSameDeg.embedding, NumDenSameDeg.embedding, Localization.mk_eq_mk_iff,
        Localization.r_iff_exists] at h
      obtain ⟨u, hu⟩ := h
      have h2 := congrArg
        (Ideal.Quotient.mk (Ideal.span {(f - 1 : MvPolynomial σ R)})) hu
      simp only [map_mul] at h2
      rwa [quot_mem_powers R f _ u.2, quot_mem_powers R f _ c.den_mem,
        quot_mem_powers R f _ c'.den_mem, one_mul, one_mul, one_mul, one_mul] at h2)

theorem quot_den_of (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (c : NumDenSameDeg (MvPolynomial.homogeneousSubmodule σ R) (Submonoid.powers f)) :
    (Ideal.Quotient.mk (Ideal.span {(f - 1 : MvPolynomial σ R)}))
      (c.den : MvPolynomial σ R) = 1 :=
  quot_mem_powers R f _ c.den_mem

/-- ★★★**逆向きの環準同型** `A⁰_f →+* R[x]/(f − 1)`。 -/
noncomputable def awayQuotHomOf (R : Type) [CommRing R] (f : MvPolynomial σ R) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule σ R) f
      →+* (MvPolynomial σ R ⧸ Ideal.span {(f - 1 : MvPolynomial σ R)}) where
  toFun := awayQuotOf R f
  map_one' := by
    show (Ideal.Quotient.mk _) ((NumDenSameDeg.num (1 : NumDenSameDeg
      (MvPolynomial.homogeneousSubmodule σ R) (Submonoid.powers f)) :
        MvPolynomial σ R)) = 1
    rw [NumDenSameDeg.num_one]
    exact map_one _
  map_zero' := by
    show (Ideal.Quotient.mk _) ((NumDenSameDeg.num (0 : NumDenSameDeg
      (MvPolynomial.homogeneousSubmodule σ R) (Submonoid.powers f)) :
        MvPolynomial σ R)) = 0
    rw [NumDenSameDeg.num_zero, ZeroMemClass.coe_zero]
    exact map_zero _
  map_mul' := by
    rintro ⟨c⟩ ⟨c'⟩
    show (Ideal.Quotient.mk _) (((c * c').num : MvPolynomial σ R)) = _
    rw [NumDenSameDeg.num_mul, map_mul]
    rfl
  map_add' := by
    rintro ⟨c⟩ ⟨c'⟩
    show (Ideal.Quotient.mk _) (((c + c').num : MvPolynomial σ R)) = _
    rw [NumDenSameDeg.num_add, map_add, map_mul, map_mul,
      quot_den_of R f c, quot_den_of R f c', one_mul, one_mul]
    show _ = (Ideal.Quotient.mk _) (c.num : MvPolynomial σ R)
      + (Ideal.Quotient.mk _) (c'.num : MvPolynomial σ R)
    ring

/-! ## ★★★★★★★★★本体 -/

/-- ★★★★★**合成は商写像である** —— `MvPolynomial.ringHom_ext` で `C` と `X` を見るだけ。 -/
theorem awayQuotHomOf_comp_awayEvalOf (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) :
    (awayQuotHomOf R f).comp (awayEvalOf R f hf)
      = Ideal.Quotient.mk (Ideal.span {(f - 1 : MvPolynomial σ R)}) := by
  refine MvPolynomial.ringHom_ext (fun r => ?_) (fun j => ?_)
  · show awayQuotOf R f (awayEvalOf R f hf (MvPolynomial.C r)) = _
    rw [awayEvalOf_C]
    rfl
  · show awayQuotOf R f (awayEvalOf R f hf (MvPolynomial.X j)) = _
    rw [awayEvalOf_X]
    rfl

/-- ★**`awayEvalOf f (f) = 1`** —— `f/f = 1`。 -/
theorem awayEvalOf_self (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) :
    awayEvalOf R f hf f = 1 := by
  rw [awayEvalOf_mk R f hf 1 f (by simpa using hf)]
  refine HomogeneousLocalization.val_injective _ ?_
  rw [Away.val_mk]
  simp

/-- ★★★★★★★★★**`ker (awayEvalOf f) = (f − 1)`**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`awayEvalOf` は全射（`§9-861c`）なので、これは

    `A⁰_f ≅ R[x]/(f − 1)`

すなわち「1 次形式 `f` の非零部分は `f = 1` と正規化したアフィン空間である」ことである。

★★`§9-859` の `ker_awayEval`（`f = x_i`）の一般化であり、
**`[Nontrivial R]` も次数の一意性も要らない**——分母が `f` の何乗かを知る必要がないから。 -/
theorem ker_awayEvalOf (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) :
    RingHom.ker (awayEvalOf R f hf) = Ideal.span {(f - 1 : MvPolynomial σ R)} := by
  apply le_antisymm
  · intro p hp
    have h : (awayQuotHomOf R f) (awayEvalOf R f hf p) = 0 := by
      rw [RingHom.mem_ker] at hp
      rw [hp, map_zero]
    rw [← RingHom.comp_apply, awayQuotHomOf_comp_awayEvalOf] at h
    exact (Ideal.Quotient.eq_zero_iff_mem).1 h
  · rw [Ideal.span_le]
    rintro z rfl
    rw [SetLike.mem_coe, RingHom.mem_ker, map_sub, map_one, awayEvalOf_self, sub_self]

/-! ## ★出典の紐付け(`.src`) -/

def awayQuotHomOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(逆向きの環準同型 A⁰_f →+* R[x]/(f − 1))",
    sectionId := "genell-prop-1-4" }

def awayEvalOf_self.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(awayEvalOf f (f) = 1)",
    sectionId := "genell-prop-1-4" }

def ker_awayEvalOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ker (awayEvalOf f) = (f − 1))",
    sectionId := "genell-prop-1-4" }

def ker_awayEvalOf.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ker_awayEval(f = x_i の場合、段 C2a-2、§9-859)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ker_awayEval") 2,
    .citation "[ABC3]" "awayEvalOf_mk(斉次式の評価は a/f^k、§9-861c)"
      (.inProject "ABC3" "ABC3.Found.GenEll.awayEvalOf_mk") 2,
    .implicitStep
      ("★★★★★測定: §9-859 は awayQuot の well-defined 性に " ++
       "exists_pow_of_numDenSameDeg(分母の指数は次数に等しい、§9-813)を使っていたが、" ++
       "**次数は要らなかった**(2026-08-28 実測)。必要なのは " ++
       "「R[x]/(f − 1) の中で powers f の元はすべて [·] = 1」だけで、" ++
       "分母が f の何乗かを**知る必要はない**。" ++
       "★結果として一般形は [Nontrivial R] も要らず、f = 0 でも成り立つ(両辺とも零環)") 3,
    .implicitStep
      ("★★これで段 C2c-1 の (c) の材料が揃った: " ++
       "(a) の四角 ＋ (b) ker σ = (x₀) ＋ 本補題 ker (awayEvalOf) = (f − 1)。" ++
       "★組み立ては Ideal.map_comap_of_surjective / Ideal.comap_map_of_surjective で行う") 4 ]

end ABC3.Found.GenEll
