import ABC3.Found.GaloisRep.MulByNCoordX
import ABC3.Found.GaloisRep.PhiDegree

/-!
# Galois (G5) 第 324 ブロック —— **★★★★★★★`x([n]P)` の超越性と `μ` の単射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★到達点

> **`x([n]P)` は `F` 上超越的**であり、したがって `[n]` の引き戻し `μ` は**単射**である

★★★これで `NondegStep` の `eq_zero_of_weilPairingVal_trivial` が要求する一式
(`hinj`・`hμF`・`hμP`・`hμx`・`hμy`)が**すべて揃った**。残るのは `hfix` だけである。

## ★★★★★★超越性は「次数を数える」だけで出た

原文の道(`n` 等分点での 1 点評価)ではなく、**代数的整数論の 1 行**で出た:

> `z·B(x) = A(x)`、`A` はモニックで `deg B < deg A`、`x` は超越的 ⟹ `z` は超越的

★★`z` が代数的だと仮定すると、`x` は `F[z]` 上**整**である
——`A − z·B` がモニック(次数は `deg A`、`deg B < deg A` だから)。
★★★`F[z]` は `F` 上整だから `isIntegral_trans` で `x` が `F` 上整になり、
`x` の超越性(第 116)に矛盾する。

★★★★★**`A = Φ_n`(モニック、次数 `n²`——第 198)、`B = ΨSq_n`(次数 `≤ n²−1`)**
がちょうどこの形である。★当初の見積もり「第 117 と同じ型、5-15 ブロック」は
**1 ブロック**で済んだ——`Φ_n` の次数を第 198 で既に測ってあったからである。

## ★★これが (G5) のどこに効くか

    (G5) 非退化性
      → eq_zero_of_weilPairingVal_trivial(第 197)の仮定
         ・hμP・hμx・hμy(第 118・125)      ✅
         ・hinj(μ の単射性)                 ✅ ★本ブロック
         ・hμF(μ は F 上恒等)               ✅ 第 118 の `pointHom_algebraMap`
         ・hfix : F(E)^{E[n]} = [n]^*F(E)    ❌ **残る**
      → hfix は deg[n] = n²、すなわち体の拡大次数の帳簿だけ

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `transcendental_of_monic_rel` | ★★★★★**分子の次数が大きい有理式の値は超越的** |
| `transcendental_coordX` | ★第 116 の `Transcendental` 版 |
| `transcendental_mulByN_coordX` | ★★★★★★**`x([n]P)` は超越的** |
| `exists_mulByNHom_full` | ★★★★★★★**`NondegStep` が要求する一式** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

/-! ## ★★★★★分子の次数が大きい有理式の値は超越的 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★**分子の次数が分母より大きい有理式の値は超越的である**。

`z·B(x) = A(x)`、`A` はモニックで `deg B < deg A`、`x` は超越的 ⟹ `z` は超越的。

★`z` が代数的だと仮定すると `x` は `F[z]` 上整(モニック多項式 `A − z·B`)。
★★`F[z]` は `F` 上整だから `x` が `F` 上整になり、`x` の超越性に矛盾する。 -/
theorem transcendental_of_monic_rel {F K : Type} [Field F] [Field K] [Algebra F K]
    (x z : K) (hx : Transcendental F x) (A B : Polynomial F)
    (hA : A.Monic) (hdeg : B.natDegree < A.natDegree)
    (h : z * Polynomial.aeval x B = Polynomial.aeval x A) :
    Transcendental F z := by
  intro hz
  refine hx ?_
  haveI : Algebra.IsIntegral F (Algebra.adjoin F ({z} : Set K)) :=
    Algebra.IsIntegral.adjoin (fun y hy => by
      rw [Set.mem_singleton_iff] at hy
      subst hy
      exact hz.isIntegral)
  have hzmem : z ∈ Algebra.adjoin F ({z} : Set K) := Algebra.subset_adjoin rfl
  set S := Algebra.adjoin F ({z} : Set K) with hS
  set f : F →+* S := algebraMap F S with hf
  have hAm : (A.map f).Monic := hA.map f
  have hdm : (A.map f).natDegree = A.natDegree := hA.natDegree_map f
  have h1 : (Polynomial.C (⟨z, hzmem⟩ : S) * B.map f).natDegree ≤ B.natDegree := by
    refine Polynomial.natDegree_mul_le.trans ?_
    rw [Polynomial.natDegree_C, zero_add]
    exact Polynomial.natDegree_map_le
  have hlt : (Polynomial.C (⟨z, hzmem⟩ : S) * B.map f).natDegree < (A.map f).natDegree := by
    omega
  have hmonic : (A.map f - Polynomial.C (⟨z, hzmem⟩ : S) * B.map f).Monic := by
    rw [sub_eq_add_neg]
    refine Polynomial.Monic.add_of_left hAm ?_
    rw [Polynomial.degree_neg]
    exact Polynomial.degree_lt_degree hlt
  have hint : IsIntegral S x := by
    refine ⟨A.map f - Polynomial.C (⟨z, hzmem⟩ : S) * B.map f, hmonic, ?_⟩
    rw [Polynomial.eval₂_sub, Polynomial.eval₂_mul, Polynomial.eval₂_C, Polynomial.eval₂_map,
      Polynomial.eval₂_map]
    have hcomp : ((algebraMap S K).comp f) = algebraMap F K := by
      ext c
      show (algebraMap S K) (f c) = algebraMap F K c
      rw [hf, ← IsScalarTower.algebraMap_apply]
    rw [hcomp]
    show Polynomial.aeval x A - (z : K) * Polynomial.aeval x B = 0
    rw [h]
    ring
  exact (isIntegral_trans x hint).isAlgebraic

/-! ## ★★★★★★`x([n]P)` の超越性 -/

variable {F : Type} [Field F]

/-- ★第 116 の超越性を `Transcendental` の形で。 -/
theorem transcendental_coordX (W : WeierstrassCurve.Affine F) : Transcendental F (coordX W) := by
  intro halg
  obtain ⟨p, hp, hp0⟩ := halg
  exact coordX_transcendental W hp hp0

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**`x([n]P)` は `F` 上超越的である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Φ_n` はモニックで次数ちょうど `n²`(第 198)、`ΨSq_n` は次数 `≤ n²−1` なので、
`transcendental_of_monic_rel` がそのまま当たる。 -/
theorem transcendental_mulByN_coordX (W : WeierstrassCurve.Affine F) (n : ℕ) (hn : 1 ≤ n)
    {x' : W.FunctionField}
    (hx' : x' * Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.ΨSq (n : ℤ))
        = Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.Φ (n : ℤ))) :
    Transcendental F x' := by
  refine transcendental_of_monic_rel (coordX W) x' (transcendental_coordX W)
    (W.Φ (n : ℤ)) (W.ΨSq (n : ℤ)) (monic_Φ W (n : ℤ)) ?_ hx'
  have h1 : (W.Φ (n : ℤ)).natDegree = n ^ 2 := by
    rw [natDegree_Φ_eq]; simp
  have h2 : (W.ΨSq (n : ℤ)).natDegree ≤ n ^ 2 - 1 := by
    have h := W.natDegree_ΨSq_le (n : ℤ)
    simpa using h
  have h3 : 1 ≤ n ^ 2 := Nat.one_le_pow 2 n (by omega)
  omega

/-! ## ★★★★★★★`NondegStep` が要求する一式 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**`[n]` の引き戻しは単射で `F` 上恒等**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`NondegStep.eq_zero_of_weilPairingVal_trivial` の仮定 `hinj`・`hμF`・`hμP`・`hμx`・`hμy`
がこれで全部そろう。残るのは `hfix` だけである。 -/
theorem exists_mulByNHom_full (W : WeierstrassCurve.Affine F) [W.IsElliptic] [DecidableEq F]
    [CharZero F] (n : ℕ) (hn : 1 ≤ n) :
    ∃ (x' y' : W.FunctionField) (h' : (W.map (algebraMap F W.FunctionField)).Nonsingular x' y')
      (μ : W.CoordinateRing →+* W.FunctionField),
      n • genericPoint W = Point.some x' y' h'
        ∧ μ (genX W) = x' ∧ μ (genY W) = y'
        ∧ Function.Injective μ
        ∧ (∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
        ∧ x' * Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.ΨSq (n : ℤ))
            = Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.Φ (n : ℤ)) := by
  classical
  obtain ⟨x', y', h', heq, hx, _⟩ :=
    mulOK_of_ne (W.map (algebraMap F W.FunctionField)) (nonsingular_coord W) n hn
      (fun k hk1 _ => psiSq_coordX_ne_zero W k hk1)
  rw [WeierstrassCurve.map_ΨSq, Polynomial.eval_map, WeierstrassCurve.map_Φ,
    Polynomial.eval_map] at hx
  have htr : Transcendental F x' := transcendental_mulByN_coordX W n hn hx
  refine ⟨x', y', h', pointHom W h'.1, heq, pointHom_genX W h'.1, pointHom_genY W h'.1, ?_,
    pointHom_algebraMap W h'.1, hx⟩
  refine pointHom_injective_of_transcendental W h'.1 (fun p hp hp0 => htr ⟨p, hp, ?_⟩)
  exact hp0

/-! ## ★出典の紐付け(`.src`) -/

def transcendental_mulByN_coordX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——x([n]P) の超越性)",
    sectionId := "genell-thm-3-8" }

def exists_mulByNHom_full.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——乗法 [n] の引き戻しの単射性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
