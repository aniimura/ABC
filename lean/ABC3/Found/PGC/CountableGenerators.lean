import ABC3.Found.PGC.CoherentFunctional
import ABC3.Found.PGC.LocalFieldNorm
import Mathlib.Analysis.Normed.Field.Krasner

/-!
# 可算個の生成元 —— Krasner の補題で `HasCofinalGaloisTower` を埋める

★**本ファイルで [pGC] Proposition 2.1 が閉じる**(`prop_2_1`、`sorry` なし・`axiom` なし)。

## 何を示すか

`K̄` の**各元**は、ある**可算集合 `D ⊆ K̄`** の元 1 つで `K` 上生成される
(`hasCountableGenerators`)。前段 `Found/PGC/CoherentFunctional.lean` は
これから可算共終 Galois 塔を作り、Maschke 上げで**両立汎関数**を作るので、
`HasCoherentFunctional K` が**無条件に**従い、`SmoothModelCarrier K` を経て
Proposition 2.1(`RecoverableAsAddModule (fun K => K.closure)`)が出る。

## ★★本体の見立てより安かった点(2026-09-08)

本体は残る穴を「**ℚ 上代数的**な元で生成される」(= Krasner + ℚ の ℚ_p での稠密性 +
`Algebraic.countable`)と見立てていた。★**その形は要らない。**
必要なのは「可算集合で生成される」だけで、係数は
**`K` 自身の可算稠密部分集合**から取ればよい(`TopologicalSpace.exists_countable_dense`。
★`K` は局所体なので可分であり、mathlib のインスタンスが**そのまま効く**)。
これで ℚ ⊆ ℚ_p ⊆ K の稠密性の持ち上げ(基底の取り替えが要る)が**丸ごと消える**。

★在庫の測り方(2026-09-08):`NormedField K.closure` は**本木に在った**
(`ABC3.Found.PGC.closureNormedField`、`Found/PGC/LocalFieldNorm.lean:150`、scoped instance)。
`IsKrasner K.carrier K.closure` は `Mathlib/Analysis/Normed/Field/Krasner.lean` の
`IsKrasner.of_completeSpace` が**インスタンスとして直接見つかる**
(`example : IsKrasner K.carrier K.closure := inferInstance` が通る)。
★この 2 つが在ったので、Krasner の道が一気に開いた。

## 内容

* §1 抽象核 —— **単項多項式の値が小さければ根の 1 つが近い**
  (`exists_root_nnnorm_sub_lt`。★体も分岐も Galois も出ない。ノルム付き体だけ)。
* §2 抽象核 —— 係数を `e : ℕ → F` から取る**可算族の単項多項式**と、その値の評価。
* §3 具体層 —— `K` の可算稠密列で近似し、Krasner で `x ∈ K⟮β⟯` を出す。

## 逸脱の記録

無し。原典 [pGC] は Proposition 2.1 の証明でこの補題を明示していないが
(「a p-adic field has only countably many finite extensions」は暗黙)、
本ファイルはその暗黙の一歩を明示にしただけで、主張は変えていない。
-/

namespace ABC3.Found.PGC

namespace ApxRoot

open Polynomial

/-! ## 1. 抽象核 —— 単項多項式の値が小さければ根の 1 つが近い -/

theorem exists_root_nnnorm_sub_lt {L : Type*} [NormedField L] {g : L[X]} (hg : g.Monic)
    (hsplit : g.Splits) {x : L} {δ : NNReal}
    (h : ‖g.eval x‖₊ < δ ^ g.natDegree) :
    ∃ β ∈ g.roots, ‖x - β‖₊ < δ := by
  by_contra hcon
  have hcon' : ∀ β ∈ g.roots, δ ≤ ‖x - β‖₊ := by
    intro β hβ
    by_contra hlt
    exact hcon ⟨β, hβ, not_le.1 hlt⟩
  have hprod : g.eval x = (g.roots.map (fun β => x - β)).prod := by
    conv_lhs => rw [hsplit.eq_prod_roots_of_monic hg]
    rw [Polynomial.eval_multiset_prod, Multiset.map_map]
    simp
  have hcard : Multiset.card g.roots = g.natDegree := Polynomial.splits_iff_card_roots.1 hsplit
  have h2 : δ ^ Multiset.card (g.roots.map (fun β => ‖x - β‖₊))
      ≤ (g.roots.map (fun β => ‖x - β‖₊)).prod := by
    refine Multiset.pow_card_le_prod ?_
    intro y hy
    obtain ⟨β, hβ, rfl⟩ := Multiset.mem_map.1 hy
    exact hcon' β hβ
  have hge : δ ^ g.natDegree ≤ ‖g.eval x‖₊ := by
    rw [hprod, ← hcard]
    have h1 : ‖(g.roots.map (fun β => x - β)).prod‖₊
        = ((g.roots.map (fun β => x - β)).map (fun z => ‖z‖₊)).prod :=
      map_multiset_prod (nnnormHom (α := L)) _
    rw [h1, Multiset.map_map]
    simpa using h2
  exact absurd h (not_lt.2 hge)

/-! ## 2. 抽象核 —— 可算族の単項多項式 -/

noncomputable def apxCoeff {F : Type*} [Field F] (e : ℕ → F) (i : Σ d : ℕ, Fin d → ℕ) (j : ℕ) : F :=
  if h : j < i.1 then e (i.2 ⟨j, h⟩) else 0

noncomputable def apxPoly {F : Type*} [Field F] (e : ℕ → F) (i : Σ d : ℕ, Fin d → ℕ) : F[X] :=
  X ^ i.1 + ∑ j ∈ Finset.range i.1, C (apxCoeff e i j) * X ^ j

variable {F : Type*} [Field F] (e : ℕ → F) (i : Σ d : ℕ, Fin d → ℕ)

theorem degree_apxPoly_aux :
    (∑ j ∈ Finset.range i.1, C (apxCoeff e i j) * X ^ j).degree < (i.1 : WithBot ℕ) := by
  refine lt_of_le_of_lt (Polynomial.degree_sum_le _ _) ?_
  refine (Finset.sup_lt_iff (by exact_mod_cast WithBot.bot_lt_coe i.1)).2 ?_
  intro j hj
  refine lt_of_le_of_lt (Polynomial.degree_C_mul_X_pow_le _ _) ?_
  exact_mod_cast Finset.mem_range.1 hj

theorem apxPoly_monic : (apxPoly e i).Monic :=
  Polynomial.monic_X_pow_add (degree_apxPoly_aux e i)

theorem apxPoly_natDegree : (apxPoly e i).natDegree = i.1 := by
  have h1 : (X ^ i.1 : F[X]).degree = (i.1 : WithBot ℕ) := Polynomial.degree_X_pow i.1
  have h2 := degree_apxPoly_aux e i
  rw [apxPoly, Polynomial.natDegree_add_eq_left_of_degree_lt (by rw [h1]; exact h2),
    Polynomial.natDegree_X_pow]

/-- ★近似多項式の値の評価。 -/
theorem nnnorm_aeval_apxPoly_le {L : Type*} [NormedField L] [Algebra F L]
    (f : F[X]) (hf : f.Monic) (x : L) (hx : Polynomial.aeval x f = 0)
    (ε : NNReal) (n : ℕ → ℕ)
    (hn : ∀ j ∈ Finset.range f.natDegree,
      ‖(algebraMap F L) (e (n j)) - (algebraMap F L) (f.coeff j)‖₊ ≤ ε) :
    ‖Polynomial.aeval x (apxPoly e ⟨f.natDegree, fun j => n (j : ℕ)⟩)‖₊
      ≤ ε * ∑ j ∈ Finset.range f.natDegree, ‖x‖₊ ^ j := by
  set i : (Σ d : ℕ, Fin d → ℕ) := ⟨f.natDegree, fun j => n (j : ℕ)⟩ with hi
  have hcoeff : ∀ j ∈ Finset.range f.natDegree, apxCoeff e i j = e (n j) := by
    intro j hj
    simp only [apxCoeff, hi]
    rw [dif_pos (Finset.mem_range.1 hj)]
  have h1 : Polynomial.aeval x (apxPoly e i)
      = x ^ f.natDegree
        + ∑ j ∈ Finset.range f.natDegree, (algebraMap F L) (apxCoeff e i j) * x ^ j := by
    simp [apxPoly, hi, map_sum]
  have h2 : (0 : L)
      = x ^ f.natDegree
        + ∑ j ∈ Finset.range f.natDegree, (algebraMap F L) (f.coeff j) * x ^ j := by
    rw [← hx]
    conv_lhs => rw [hf.as_sum]
    simp [map_sum]
  have hxd : (x : L) ^ f.natDegree
      = - ∑ j ∈ Finset.range f.natDegree, (algebraMap F L) (f.coeff j) * x ^ j :=
    eq_neg_of_add_eq_zero_left h2.symm
  have hB : ∑ j ∈ Finset.range f.natDegree, (algebraMap F L) (apxCoeff e i j) * x ^ j
      = ∑ j ∈ Finset.range f.natDegree, (algebraMap F L) (e (n j)) * x ^ j :=
    Finset.sum_congr rfl (fun j hj => by rw [hcoeff j hj])
  have hsub : ∑ j ∈ Finset.range f.natDegree,
        ((algebraMap F L) (e (n j)) - (algebraMap F L) (f.coeff j)) * x ^ j
      = (∑ j ∈ Finset.range f.natDegree, (algebraMap F L) (e (n j)) * x ^ j)
        - ∑ j ∈ Finset.range f.natDegree, (algebraMap F L) (f.coeff j) * x ^ j := by
    rw [← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl (fun j _ => sub_mul _ _ _)
  have h3 : Polynomial.aeval x (apxPoly e i)
      = ∑ j ∈ Finset.range f.natDegree,
        ((algebraMap F L) (e (n j)) - (algebraMap F L) (f.coeff j)) * x ^ j := by
    rw [hsub, h1, hB, hxd]
    abel
  rw [h3]
  refine le_trans (nnnorm_sum_le _ _) ?_
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum (fun j hj => ?_)
  rw [nnnorm_mul, nnnorm_pow]
  exact mul_le_mul' (hn j hj) le_rfl

end ApxRoot

/-! ## 3. 具体層 -/

section Concrete

open ABC3.Skeleton.PGC ApxRoot Polynomial

variable {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)

/-- ★★★**可算個の生成元**:`K̄` の各元は、可算集合 `D` の元 1 つで `K` 上生成される。 -/
theorem hasCountableGenerators (K : PAdicLocalField p) : HasCountableGenerators K := by
  classical
  obtain ⟨s, hscount, hsdense⟩ := TopologicalSpace.exists_countable_dense K.carrier
  obtain ⟨e, he⟩ := hscount.exists_eq_range hsdense.nonempty
  have hdense : DenseRange e := by
    rw [DenseRange, ← he]
    exact hsdense
  have hclose : ∀ (y : K.carrier) (ε : NNReal), 0 < ε → ∃ m : ℕ, ‖e m - y‖₊ < ε := by
    intro y ε hε
    obtain ⟨m, hm⟩ := Metric.denseRange_iff.1 hdense y ε (by exact_mod_cast hε)
    refine ⟨m, ?_⟩
    rw [← NNReal.coe_lt_coe]
    calc (‖e m - y‖₊ : ℝ) = dist (e m) y := by rw [coe_nnnorm, dist_eq_norm]
      _ = dist y (e m) := dist_comm _ _
      _ < ε := hm
  set φ := algebraMap K.carrier K.closure with hφ
  refine ⟨⋃ i : (Σ d : ℕ, Fin d → ℕ),
      {α : K.closure | α ∈ ((apxPoly e i).map φ).roots}, Set.countable_iUnion (fun i => ?_), ?_⟩
  · exact (Set.Finite.ofFinset (((apxPoly e i).map φ).roots.toFinset)
      (by intro y; simp)).countable
  · intro x
    have hint : IsIntegral K.carrier x := Algebra.IsIntegral.isIntegral x
    have hfm : (minpoly K.carrier x).Monic := minpoly.monic hint
    have hfsep : (minpoly K.carrier x).Separable := Algebra.IsSeparable.isSeparable K.carrier x
    have hfsplit : ((minpoly K.carrier x).map φ).Splits := IsAlgClosed.splits _
    have hfmap0 : (minpoly K.carrier x).map φ ≠ 0 := (hfm.map φ).ne_zero
    have hd : 0 < (minpoly K.carrier x).natDegree := minpoly.natDegree_pos hint
    -- 共役根からの距離の下界 δ
    obtain ⟨δ, hδpos, hδ⟩ : ∃ δ : NNReal, 0 < δ ∧
        ∀ x' ∈ ((minpoly K.carrier x).map φ).roots, x' ≠ x → δ ≤ ‖x - x'‖₊ := by
      rcases (((minpoly K.carrier x).map φ).roots.toFinset.erase x).eq_empty_or_nonempty with
        hTe | hTn
      · refine ⟨1, one_pos, ?_⟩
        intro x' hx' hne
        exfalso
        have hmem : x' ∈ ((minpoly K.carrier x).map φ).roots.toFinset.erase x :=
          Finset.mem_erase.2 ⟨hne, Multiset.mem_toFinset.2 hx'⟩
        rw [hTe] at hmem
        exact absurd hmem (Finset.notMem_empty x')
      · obtain ⟨y, hy, hmin⟩ :=
          (((minpoly K.carrier x).map φ).roots.toFinset.erase x).exists_min_image
            (fun x' => ‖x - x'‖₊) hTn
        refine ⟨‖x - y‖₊, nnnorm_pos.2 (sub_ne_zero.2 (Ne.symm (Finset.mem_erase.1 hy).1)), ?_⟩
        intro x' hx' hne
        exact hmin x' (Finset.mem_erase.2 ⟨hne, Multiset.mem_toFinset.2 hx'⟩)
    -- 誤差の設定
    set S : NNReal := ∑ j ∈ Finset.range (minpoly K.carrier x).natDegree, ‖x‖₊ ^ j with hS
    have hδ2 : δ / 2 < δ := NNReal.half_lt_self (ne_of_gt hδpos)
    have hδ2pos : 0 < δ / 2 := by
      rw [pos_iff_ne_zero]
      exact div_ne_zero (ne_of_gt hδpos) two_ne_zero
    set ε : NNReal := (δ / 2) ^ (minpoly K.carrier x).natDegree / (S + 1) with hε
    have hS1 : S + 1 ≠ 0 := by
      intro hcon
      simpa using congrArg (fun z : NNReal => z) hcon
    have hεpos : 0 < ε := by
      rw [hε, pos_iff_ne_zero]
      exact div_ne_zero (ne_of_gt (pow_pos hδ2pos _)) hS1
    choose n hn using fun j : ℕ => hclose ((minpoly K.carrier x).coeff j) ε hεpos
    set i : (Σ d : ℕ, Fin d → ℕ) :=
      ⟨(minpoly K.carrier x).natDegree, fun j => n (j : ℕ)⟩ with hi
    have hbound : ‖Polynomial.aeval x (apxPoly e i)‖₊ ≤ ε * S := by
      refine nnnorm_aeval_apxPoly_le e (minpoly K.carrier x) hfm x (minpoly.aeval _ _) ε n ?_
      intro j _
      rw [← map_sub, nnnorm_algebraMap, nnnorm_one, mul_one]
      exact le_of_lt (hn j)
    have hmid : ε * S ≤ (δ / 2) ^ (minpoly K.carrier x).natDegree := by
      calc ε * S ≤ ε * (S + 1) := by gcongr; exact le_self_add
        _ = (δ / 2) ^ (minpoly K.carrier x).natDegree := div_mul_cancel₀ _ hS1
    have hstrict : (δ / 2) ^ (minpoly K.carrier x).natDegree
        < δ ^ (minpoly K.carrier x).natDegree :=
      pow_lt_pow_left₀ hδ2 zero_le hd.ne'
    have hlt : ‖Polynomial.aeval x (apxPoly e i)‖₊
        < δ ^ (minpoly K.carrier x).natDegree :=
      lt_of_le_of_lt (le_trans hbound hmid) hstrict
    -- 根を取る
    have hgm : ((apxPoly e i).map φ).Monic := (apxPoly_monic e i).map φ
    have hgdeg : ((apxPoly e i).map φ).natDegree = (minpoly K.carrier x).natDegree := by
      rw [(apxPoly_monic e i).natDegree_map, apxPoly_natDegree, hi]
    have hgeval : ((apxPoly e i).map φ).eval x = Polynomial.aeval x (apxPoly e i) := by
      rw [Polynomial.eval_map, ← Polynomial.aeval_def]
    obtain ⟨β, hβroot, hβclose⟩ :=
      exists_root_nnnorm_sub_lt hgm (IsAlgClosed.splits _) (x := x) (δ := δ)
        (by rw [hgeval, hgdeg]; exact hlt)
    refine ⟨β, Set.mem_iUnion.2 ⟨i, hβroot⟩, ?_⟩
    have hkr : x ∈ IntermediateField.adjoin K.carrier {β} := by
      have := IsKrasner.krasner (K := K.carrier) (L := K.closure) (x := x) (y := β)
        hfsep hfsplit (Algebra.IsIntegral.isIntegral β) ?_
      · simpa using this
      · intro x' hconj hne
        have hx'root : x' ∈ ((minpoly K.carrier x).map φ).roots := by
          rw [Polynomial.mem_roots hfmap0]
          show ((minpoly K.carrier x).map φ).eval x' = 0
          rw [Polynomial.eval_map, ← Polynomial.aeval_def]
          exact hconj.aeval_eq_zero
        have hge := hδ x' hx'root (Ne.symm hne)
        calc ‖x - β‖ = ((‖x - β‖₊ : NNReal) : ℝ) := rfl
          _ < ((δ : NNReal) : ℝ) := by exact_mod_cast hβclose
          _ ≤ ((‖x - x'‖₊ : NNReal) : ℝ) := by exact_mod_cast hge
          _ = ‖x - x'‖ := rfl
    exact hkr

/-- ★★★**可算共終 Galois 塔の存在**(Krasner)。 -/
theorem hasCofinalGaloisTower (K : PAdicLocalField p) : HasCofinalGaloisTower K :=
  hasCofinalGaloisTower_of_countableGenerators K (hasCountableGenerators K)

/-- ★★★**両立汎関数の存在**。 -/
theorem hasCoherentFunctional (K : PAdicLocalField p) : HasCoherentFunctional K :=
  hasCoherentFunctional_of_cofinalTower K (hasCofinalGaloisTower K)

/-- ★★★**[pGC] Proposition 2.1**。 -/
theorem prop_2_1 : RecoverableAsAddModule (p := p) (fun K => K.closure) :=
  prop_2_1_of_hasCoherentFunctional (fun K => hasCoherentFunctional K)

def prop_2_1.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 4, item := "Proposition 2.1", sectionId := "prop-2-1" }

end Concrete

end ABC3.Found.PGC
