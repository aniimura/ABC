/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.ModularExample

/-!
# [CorrHyp] §2 の非空虚性——`SL(2,ℤ)` は Hecke 対角共役で「無限個の correspondence を持つ」

`Definition 2.1`(`InfinitelyManyCorr`)は「`Γ` が `Comm(Γ)` の中で無限指数」
という主張——`SL(2,ℤ)` は古典的に arithmetic で、これが成り立つことが
知られている(commensurator は `PGL(2,ℚ)` で稠密)。ここでは
`g := diag(2, 1/2) ∈ SL(2,ℝ)`(Hecke 型の対角共役)を使い、
`g ∈ Comm(Γ_SL2Z)` かつ `{g^n : n ∈ ℕ}` が `Γ_SL2Z` の中で相異なる
無限個の剰余類を生む(`g^n ∉ Γ_SL2Z` for `n ≥ 1`、成分計算)ことを示し、
`Definition 2.1` を `corrHypInstance` において非空虚に実現する。

数学的な筋: `g⁻¹ M g = [[a,b/4],[4c,d]]`・`g M g⁻¹ = [[a,4b],[c/4,d]]`
(`M=[[a,b],[c,d]]`)を成分計算で確認し、`Γ(4)`(mathlib に `FiniteIndex`
インスタンス済み)が**両方向**の有限指数下界になることを示す——これで
`g ∈ Comm(Γ_SL2Z)` が言え、`g^n`(`n∈ℕ`)の (1,1) 成分が `2⁻ⁿ ∈ (0,1)`
で整数になりえないことから `g^n ∉ Γ_SL2Z`(`n≥1`)、ゆえに
`{g^n Γ_SL2Z}` が可算無限個の相異なる剰余類を与える。 -/

namespace ABC3.Found.CorrHyp

open Matrix
open scoped Pointwise

/-- Hecke 型の対角行列 `diag(2, 1/2) ∈ SL(2,ℝ)`。 -/
noncomputable def g2 : Matrix.SpecialLinearGroup (Fin 2) ℝ :=
  ⟨!![2, 0; 0, 1/2], by simp [Matrix.det_fin_two]⟩

/-- `Γ(4)`(mod 4 で単位行列に合同)の `SL(2,ℝ)` への像。`Γ_SL2Z` の中で
両方向の有限指数下界として使う(`Γ(2)` の例と同じ、`ker(φ₂)=⊥` を
`Subgroup.relIndex_map_map` で消し込むだけ)。 -/
noncomputable def Γ_Gamma4 : Subgroup (Matrix.SpecialLinearGroup (Fin 2) ℝ) :=
  (CongruenceSubgroup.Gamma 4).map φ₂

theorem Γ_Gamma4_finiteIndex_in_SL2Z : (Γ_Gamma4.subgroupOf Γ_SL2Z).FiniteIndex := by
  rw [Subgroup.finiteIndex_iff]
  show Γ_Gamma4.relIndex Γ_SL2Z ≠ 0
  have hker : φ₂.ker = ⊥ := (MonoidHom.ker_eq_bot_iff φ₂).mpr φ₂_injective
  have hrange : Γ_SL2Z = Subgroup.map φ₂ ⊤ := MonoidHom.range_eq_map φ₂
  rw [hrange]
  unfold Γ_Gamma4
  rw [Subgroup.relIndex_map_map, hker, sup_bot_eq, sup_bot_eq, Subgroup.relIndex_top_right]
  exact (CongruenceSubgroup.instFiniteIndexGamma 4).index_ne_zero

/-- `g⁻¹ φ₂(M) g = [[a, b/4],[4c, d]]`(`M=[[a,b],[c,d]]`)——
`g := g2` による内向きの Hecke 対角共役。成分計算で直接確認する。 -/
theorem conj_inv_formula (M : Matrix.SpecialLinearGroup (Fin 2) ℤ) :
    (g2⁻¹ * φ₂ M * g2 : Matrix (Fin 2) (Fin 2) ℝ) =
      !![(M 0 0 : ℝ), (M 0 1 : ℝ)/4; 4*(M 1 0 : ℝ), (M 1 1 : ℝ)] := by
  simp only [g2, φ₂]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Matrix.vecMul, Fin.sum_univ_two,
      Matrix.SpecialLinearGroup.map, Matrix.SpecialLinearGroup.coe_inv,
      Matrix.adjugate_fin_two, Matrix.vecHead, Matrix.vecTail] <;> ring

/-- `g φ₂(M) g⁻¹ = [[a, 4b],[c/4, d]]`——外向きの Hecke 対角共役。 -/
theorem conj_formula (M : Matrix.SpecialLinearGroup (Fin 2) ℤ) :
    (g2 * φ₂ M * g2⁻¹ : Matrix (Fin 2) (Fin 2) ℝ) =
      !![(M 0 0 : ℝ), 4*(M 0 1 : ℝ); (M 1 0 : ℝ)/4, (M 1 1 : ℝ)] := by
  simp only [g2, φ₂]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Matrix.vecMul, Fin.sum_univ_two,
      Matrix.SpecialLinearGroup.map, Matrix.SpecialLinearGroup.coe_inv,
      Matrix.adjugate_fin_two, Matrix.vecHead, Matrix.vecTail] <;> ring

/-- `Γ(4) ≤ g2 • Γ_SL2Z`——`M ∈ Γ(4)` なら `4|b(M)` なので
`g2⁻¹ φ₂(M) g2`(`conj_inv_formula`)は整数行列になる。 -/
theorem Γ_Gamma4_le_smul : Γ_Gamma4 ≤ ConjAct.toConjAct g2 • Γ_SL2Z := by
  rintro x ⟨M, hM, rfl⟩
  set A : Matrix (Fin 2) (Fin 2) ℤ := (M : Matrix (Fin 2) (Fin 2) ℤ) with hA
  have hAdet : A.det = 1 := M.2
  rw [Subgroup.mem_pointwise_smul_iff_inv_smul_mem, ConjAct.smul_def]
  simp only [ConjAct.ofConjAct_inv, ConjAct.ofConjAct_toConjAct, inv_inv]
  obtain ⟨ha, hb, hc, hd⟩ := CongruenceSubgroup.Gamma_mem.mp hM
  rw [← hA] at ha hb hc hd
  have h4 : (4:ℤ) ∣ A 0 1 := by rwa [ZMod.intCast_zmod_eq_zero_iff_dvd] at hb
  obtain ⟨k, hk⟩ := h4
  have hAdet' : A 0 0 * A 1 1 - A 0 1 * A 1 0 = 1 := by rwa [Matrix.det_fin_two] at hAdet
  rw [hk] at hAdet'
  set N : Matrix (Fin 2) (Fin 2) ℤ := !![A 0 0, k; 4 * A 1 0, A 1 1] with hN
  have hNdet : N.det = 1 := by
    rw [hN, Matrix.det_fin_two_of]
    linarith [hAdet', mul_assoc (4:ℤ) k (A 1 0), mul_comm (4:ℤ) k]
  refine ⟨⟨N, hNdet⟩, ?_⟩
  apply Matrix.SpecialLinearGroup.ext
  intro i
  have hgoal := congrFun (conj_inv_formula M) i
  intro j
  fin_cases i <;> fin_cases j <;>
    simp_all [φ₂, Matrix.SpecialLinearGroup.map, Matrix.SpecialLinearGroup.coe_mk, N]

/-- `Γ(4) ≤ g2⁻¹ • Γ_SL2Z`——`M ∈ Γ(4)` なら `4|c(M)` なので
`g2 φ₂(M) g2⁻¹`(`conj_formula`)は整数行列になる。 -/
theorem Γ_Gamma4_le_inv_smul : Γ_Gamma4 ≤ ConjAct.toConjAct g2⁻¹ • Γ_SL2Z := by
  rintro x ⟨M, hM, rfl⟩
  set A : Matrix (Fin 2) (Fin 2) ℤ := (M : Matrix (Fin 2) (Fin 2) ℤ) with hA
  have hAdet : A.det = 1 := M.2
  rw [Subgroup.mem_pointwise_smul_iff_inv_smul_mem, ConjAct.smul_def]
  simp only [ConjAct.ofConjAct_inv, ConjAct.ofConjAct_toConjAct, inv_inv]
  obtain ⟨ha, hb, hc, hd⟩ := CongruenceSubgroup.Gamma_mem.mp hM
  rw [← hA] at ha hb hc hd
  have h4 : (4:ℤ) ∣ A 1 0 := by rwa [ZMod.intCast_zmod_eq_zero_iff_dvd] at hc
  obtain ⟨k, hk⟩ := h4
  have hAdet' : A 0 0 * A 1 1 - A 0 1 * A 1 0 = 1 := by rwa [Matrix.det_fin_two] at hAdet
  rw [hk] at hAdet'
  set N : Matrix (Fin 2) (Fin 2) ℤ := !![A 0 0, 4 * A 0 1; k, A 1 1] with hN
  have hNdet : N.det = 1 := by
    rw [hN, Matrix.det_fin_two_of]
    linarith [hAdet', mul_assoc (A 0 1) 4 k, mul_comm (A 0 1) 4]
  refine ⟨⟨N, hNdet⟩, ?_⟩
  apply Matrix.SpecialLinearGroup.ext
  intro i
  have hgoal := congrFun (conj_formula M) i
  intro j
  fin_cases i <;> fin_cases j <;>
    simp_all [φ₂, Matrix.SpecialLinearGroup.map, Matrix.SpecialLinearGroup.coe_mk, N]

/-- `relIndex` の単調性: `H ≤ K` かつ `H` が `L` の中で有限指数なら `K` もそう。 -/
theorem relIndex_mono {G : Type} [Group G] (H K L : Subgroup G) (h1 : H ≤ K)
    (h3 : H.relIndex L ≠ 0) : K.relIndex L ≠ 0 := by
  have hsub : H.subgroupOf L ≤ K.subgroupOf L := Subgroup.subgroupOf_mono L h1
  have hH : (H.subgroupOf L).FiniteIndex := Subgroup.finiteIndex_iff.mpr h3
  exact (Subgroup.finiteIndex_of_le hsub).index_ne_zero

/-- **`g2` はモジュラー群の commensurator に属する**——`Γ(4)` を両方向の
有限指数下界として使い、`Commensurable` の2条件をそれぞれ示す。

★**sorry 無し**。標準3公理のみ。 -/
theorem g2_mem_commensurator : g2 ∈ Subgroup.Commensurable.commensurator Γ_SL2Z := by
  rw [Subgroup.Commensurable.commensurator_mem_iff]
  refine ⟨?_, ?_⟩
  · exact relIndex_mono Γ_Gamma4 _ Γ_SL2Z Γ_Gamma4_le_smul
      (Subgroup.finiteIndex_iff.mp Γ_Gamma4_finiteIndex_in_SL2Z)
  · have key : (ConjAct.toConjAct g2⁻¹ • Γ_SL2Z).relIndex Γ_SL2Z ≠ 0 :=
      relIndex_mono Γ_Gamma4 _ Γ_SL2Z Γ_Gamma4_le_inv_smul
        (Subgroup.finiteIndex_iff.mp Γ_Gamma4_finiteIndex_in_SL2Z)
    have heq := Subgroup.relIndex_pointwise_smul (ConjAct.toConjAct g2)
      (ConjAct.toConjAct g2⁻¹ • Γ_SL2Z) Γ_SL2Z
    rw [smul_smul, ← map_mul, mul_inv_cancel, map_one, one_smul] at heq
    rw [heq]
    exact key

theorem g2_pow_eq (k : ℕ) : (g2^k : Matrix (Fin 2) (Fin 2) ℝ) = !![(2:ℝ)^k, 0; 0, (2⁻¹:ℝ)^k] := by
  induction k with
  | zero => ext i j; fin_cases i <;> fin_cases j <;> simp
  | succ n ih =>
    rw [pow_succ, ih]
    have hg2 : (g2 : Matrix (Fin 2) (Fin 2) ℝ) = !![2, 0; 0, 2⁻¹] := by simp [g2]
    rw [hg2]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Matrix.mul_apply, Fin.sum_univ_two, pow_succ]

/-- `g2^k ∉ Γ_SL2Z`(`k ≥ 1`)——`(1,1)` 成分 `2⁻ᵏ` が `(0,1)` の中にあり
整数になりえない。 -/
theorem g2_pow_not_mem (k : ℕ) (hk : 1 ≤ k) : (g2^k : Matrix.SpecialLinearGroup (Fin 2) ℝ) ∉ Γ_SL2Z := by
  rintro ⟨N, hN⟩
  have hN' : Matrix.SpecialLinearGroup.map (Int.castRingHom ℝ) N = g2^k := hN
  have hentry := congrArg (fun x : Matrix.SpecialLinearGroup (Fin 2) ℝ => (x : Matrix (Fin 2) (Fin 2) ℝ) 1 1) hN'
  simp [Matrix.SpecialLinearGroup.map] at hentry
  rw [g2_pow_eq] at hentry
  simp only [Matrix.of_apply, Matrix.cons_val_one, Matrix.cons_val_zero] at hentry
  have hlt1 : (2⁻¹:ℝ)^k < 1 := by
    apply pow_lt_one₀ (by norm_num) (by norm_num) (by omega)
  have hgt0 : (0:ℝ) < (2⁻¹:ℝ)^k := by positivity
  rw [← hentry] at hlt1 hgt0
  have hlt1' : (N 1 1 : ℤ) < 1 := by exact_mod_cast hlt1
  have hgt0' : (0:ℤ) < (N 1 1 : ℤ) := by exact_mod_cast hgt0
  omega

theorem g2_pow_mem_commensurator (n : ℕ) : g2^n ∈ Subgroup.Commensurable.commensurator Γ_SL2Z :=
  Subgroup.pow_mem _ g2_mem_commensurator n

noncomputable def g2pow (n : ℕ) : Subgroup.Commensurable.commensurator Γ_SL2Z :=
  ⟨g2^n, g2_pow_mem_commensurator n⟩

theorem pow_inv_mul_pow_eq_pow_sub {G : Type} [Group G] (g : G) (n m : ℕ) (h : n ≤ m) :
    (g^n)⁻¹ * g^m = g^(m-n) := by
  rw [eq_comm, ← Nat.sub_add_cancel h, pow_add]
  group
  congr 1
  omega

/-- `n ↦ [g2^n]` は `Γ_SL2Z` を法とした剰余類の中で単射——`g2^(m-n) ∉ Γ_SL2Z`
(`n<m` のとき、`g2_pow_not_mem`)を使う。 -/
theorem g2pow_injective_mod :
    Function.Injective (fun n : ℕ => (QuotientGroup.mk (g2pow n) :
      (Subgroup.Commensurable.commensurator Γ_SL2Z) ⧸
        (Γ_SL2Z.subgroupOf (Subgroup.Commensurable.commensurator Γ_SL2Z)))) := by
  have key : ∀ n m : ℕ, n ≤ m →
      (QuotientGroup.mk (g2pow n) : (Subgroup.Commensurable.commensurator Γ_SL2Z) ⧸
        (Γ_SL2Z.subgroupOf (Subgroup.Commensurable.commensurator Γ_SL2Z))) =
        QuotientGroup.mk (g2pow m) → n = m := by
    intro n m h heq
    rw [QuotientGroup.eq] at heq
    simp only [Subgroup.mem_subgroupOf] at heq
    have hval : ((g2pow n)⁻¹ * g2pow m : Subgroup.Commensurable.commensurator Γ_SL2Z).1 =
        (g2^n)⁻¹ * g2^m := rfl
    rw [hval, pow_inv_mul_pow_eq_pow_sub g2 n m h] at heq
    by_contra hne
    exact g2_pow_not_mem (m - n) (by omega) heq
  intro n m heq
  rcases le_total n m with h | h
  · exact key n m h heq
  · exact (key m n h heq.symm).symm

/-- **`Γ_SL2Z` は自身の commensurator の中で無限指数**——`{g2^n}` が
可算無限個の相異なる剰余類を与える(`g2pow_injective_mod`)ので、
`Nat.card` の意味で指数は `0`(有限でない)。

★**sorry 無し**。標準3公理のみ。 -/
theorem Γ_SL2Z_not_finiteIndex_in_commensurator :
    ¬ (Γ_SL2Z.subgroupOf (Subgroup.Commensurable.commensurator Γ_SL2Z)).FiniteIndex := by
  rw [Subgroup.finiteIndex_iff, not_not, Subgroup.index_eq_zero_iff_infinite]
  exact Infinite.of_injective _ g2pow_injective_mod

open ABC3.Skeleton.CorrHyp in
/-- **[CorrHyp] `Definition 2.1`(`InfinitelyManyCorr`)の `corrHypInstance`
における非空虚性**——モジュラー群 `SL(2,ℤ)` は(古典的に arithmetic なので)
自身の commensurator の中で無限指数を持つ、という事実を Hecke 型の対角共役
`g2 = diag(2,1/2)` で具体的に構成した。

★**sorry 無し**。標準3公理のみ。 -/
theorem infinitelyManyCorr_witness : InfinitelyManyCorr corrHypInstance Γ_SL2Z := by
  intro ⟨_, hfi⟩
  exact Γ_SL2Z_not_finiteIndex_in_commensurator hfi

def infinitelyManyCorr_witness.src : ABC3.Meta.Source :=
  { paper := "CorrHyp", pdfPage := 5, item := "Definition 2.1", sectionId := "corrhyp-def-2-1" }

end ABC3.Found.CorrHyp
