/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.GalRep
import ABC3.Found.GaloisRep.TranslateEquiv
import ABC3.Interface.GaloisRep.Torsion
import Mathlib.FieldTheory.KrullTopology
import ABC3.Found.GaloisRep.TateLevel
import ABC3.Found.GaloisRep.PadicLinear
import ABC3.Meta.Claim

/-!
# `galRep` の連続性へ —— 座標を固定する `σ` は点を動かさない（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★これは何か

`Found/GenEll/GalRepClosed.lean`（`§9-1191`、第 765）で、
`Theorem 3.8` に残る位相の側は **`galRep` の連続性**だけになった。
その連続性の証明は

1. `σ` が `E[l^n]` の座標をすべて固定するなら `σ` は `E[l^n]` に自明に作用する ←**本ファイル**
2. `L(E[l^n])` は `L` の有限拡大であり、それを固定する部分群は Krull 位相で開
3. したがって `galMat σ ≡ 1 (mod l^n)` は開集合の上で成り立ち、`galRep` は連続

という 3 段である。本ファイルは**第 1 段**を取る。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve WeierstrassCurve.Affine
open scoped Pointwise

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-- ★★★★★★**座標を固定する `σ` は点を動かさない**——★**無条件**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`galPoint W σ = Point.map σ` なので、`some x y` の像は `some (σ x) (σ y)` である。 -/
theorem galPoint_eq_of_fixed (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L)
    (P : (W.baseChange L).toAffine.Point)
    (h : ∀ (x y : L) (hns : (W.baseChange L).toAffine.Nonsingular x y),
      P = Point.some x y hns → σ x = x ∧ σ y = y) :
    galPoint W σ P = P := by
  cases P with
  | zero => exact Point.map_zero _
  | some x y hns =>
    obtain ⟨hx, hy⟩ := h x y hns rfl
    rw [galPoint, Point.map_some]
    exact point_some_congr hx hy _ _

/-- ★★★★★**すべての点で座標を固定するなら `galPoint` は恒等写像**。 -/
theorem galPoint_eq_id_of_fixed (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L)
    (S : Set ((W.baseChange L).toAffine.Point))
    (h : ∀ P ∈ S, ∀ (x y : L) (hns : (W.baseChange L).toAffine.Nonsingular x y),
      P = Point.some x y hns → σ x = x ∧ σ y = y) :
    ∀ P ∈ S, galPoint W σ P = P :=
  fun P hP => galPoint_eq_of_fixed W σ P (h P hP)

/-! ## ★★★★★★★★捩れ点の座標が生成する有限拡大 -/

/-- ★点の座標の集合（`0` では空）。 -/
def ptCoordSet (W : WeierstrassCurve K) : (W.baseChange L).toAffine.Point → Set L
  | .zero => ∅
  | .some x y _ => {x, y}

theorem finite_ptCoordSet (W : WeierstrassCurve K) (P : (W.baseChange L).toAffine.Point) :
    (ptCoordSet W P).Finite := by
  cases P with
  | zero => exact Set.finite_empty
  | some x y _ => exact (Set.finite_singleton y).insert x

/-- ★★`σ` が座標集合を固定するなら点を動かさない。 -/
theorem galPoint_eq_of_fixed_coordSet (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L)
    (P : (W.baseChange L).toAffine.Point) (h : ∀ z ∈ ptCoordSet W P, σ z = z) :
    galPoint W σ P = P := by
  refine galPoint_eq_of_fixed W σ P ?_
  intro x y hns hP
  subst hP
  refine ⟨h x ?_, h y ?_⟩
  · exact Or.inl rfl
  · exact Or.inr rfl

/-- ★★★**`n`-捩れ点の座標全体**。 -/
def torsionCoordSet (W : WeierstrassCurve K) (n : ℕ) : Set L :=
  ⋃ P ∈ (torsionPoints (W.baseChange L) n : Set ((W.baseChange L).toAffine.Point)),
    ptCoordSet W P

theorem finite_torsionCoordSet (W : WeierstrassCurve K) (n : ℕ)
    (hfin : (torsionPoints (W.baseChange L) n : Set ((W.baseChange L).toAffine.Point)).Finite) :
    (torsionCoordSet W (L := L) n).Finite :=
  hfin.biUnion (fun P _ => finite_ptCoordSet (L := L) W P)

/-- ★★★★★★★★★★**捩れ点の座標が生成する有限拡大**——その固定部分群は
`n`-捩れに自明に作用する。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが `galRep` の連続性の**第 2 段**である
（Krull 位相の開部分群は有限次中間体の `fixingSubgroup` だから）。 -/
theorem exists_finiteDimensional_fixing_torsion [Algebra.IsAlgebraic K L]
    (W : WeierstrassCurve K) (n : ℕ)
    (hfin : (torsionPoints (W.baseChange L) n : Set ((W.baseChange L).toAffine.Point)).Finite) :
    ∃ F : IntermediateField K L, FiniteDimensional K F ∧
      ∀ σ ∈ F.fixingSubgroup, ∀ P ∈ torsionPoints (W.baseChange L) n,
        galPoint W σ P = P := by
  haveI : Finite (torsionCoordSet W (L := L) n) := (finite_torsionCoordSet W n hfin).to_subtype
  refine ⟨IntermediateField.adjoin K (torsionCoordSet W (L := L) n), ?_, ?_⟩
  · exact IntermediateField.finiteDimensional_adjoin
      (fun x _ => Algebra.IsIntegral.isIntegral x)
  · intro σ hσ P hP
    refine galPoint_eq_of_fixed_coordSet W σ P (fun z hz => ?_)
    have hmem : z ∈ IntermediateField.adjoin K (torsionCoordSet W (L := L) n) := by
      refine IntermediateField.subset_adjoin _ _ ?_
      exact Set.mem_biUnion hP hz
    exact (IntermediateField.mem_fixingSubgroup_iff _ _).1 hσ z hmem

/-! ## ★★★★★★★★★★★★`E[l^n]` を固定するなら `galMat ≡ 1 (mod l^n)` -/

/-- ★★★★★★★★★★★★★★
**`σ` が `E[l^n]` を固定するなら `galMat σ ≡ 1 (mod l^n)`**——★**無条件**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが `galRep` の連続性の**第 3 段**である。機構:

1. `σ` が `E[l^n]` を固定 ⟹ `galTate σ y − y` の第 `n` 射影は `0`
2. `ker(proj_n) = l^n·T_l`（`§9-1194`、第 768）⟹ `galTate σ y − y = l^n·w`
3. `e` で写して `(galMat σ − 1)·(e y) = l^n·(e w)`
4. `e` は全射なので `v = Pi.single j 1` を取れば列ごとに `l^n` で割れる -/
theorem galMat_sub_one_dvd (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime] (n : ℕ)
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ : L ≃ₐ[K] L)
    (hfix : ∀ P ∈ torsionPoints (W.baseChange L) (l ^ n), galPoint W σ P = P) (i j : Fin 2) :
    ((l : ℤ_[l]) ^ n) ∣ (galMat W l e σ i j - (1 : Matrix (Fin 2) (Fin 2) ℤ_[l]) i j) := by
  -- (1)(2) 任意の `y` について `e (galTate σ y) − e y ∈ l^n · (Fin 2 → ℤ_[l])`
  have hkey : ∀ y : tateModule (W.baseChange L) l,
      ∃ u : Fin 2 → ℤ_[l], e (galTate W l σ y) - e y = (l ^ n) • u := by
    intro y
    have hproj : (((galTate W l σ y - y :
        tateModule (W.baseChange L) l) :
        ∀ m : ℕ, torsionPoints (W.baseChange L) (l ^ m)) n :
        (W.baseChange L).toAffine.Point) = 0 := by
      show galPoint W σ (((y : ∀ m : ℕ, torsionPoints (W.baseChange L) (l ^ m)) n :
          (W.baseChange L).toAffine.Point))
        - (((y : ∀ m : ℕ, torsionPoints (W.baseChange L) (l ^ m)) n :
          (W.baseChange L).toAffine.Point)) = 0
      rw [hfix _ ((y : ∀ m : ℕ, torsionPoints (W.baseChange L) (l ^ m)) n).2, sub_self]
    obtain ⟨w, hw⟩ := exists_smul_of_proj_zero l n (galTate W l σ y - y) hproj
    refine ⟨e (w : tateModule (W.baseChange L) l), ?_⟩
    calc e (galTate W l σ y) - e y
        = e (galTate W l σ y - y) := (map_sub e _ _).symm
      _ = e ((l ^ n) • (w : tateModule (W.baseChange L) l)) := congrArg e hw.symm
      _ = (l ^ n) • e (w : tateModule (W.baseChange L) l) :=
          map_nsmul e (l ^ n) (w : tateModule (W.baseChange L) l)
  -- (3)(4) `e` の全射性で基底ベクトルを取る
  obtain ⟨y, hy⟩ := e.surjective (Pi.single j (1 : ℤ_[l]))
  obtain ⟨u, hu⟩ := hkey y
  have hmat : e (galTate W l σ y) = Matrix.mulVec (galMat W l e σ) (e y) :=
    matOf_apply e (galTate W l σ) y
  rw [hmat, hy] at hu
  have hcomp := congrFun hu i
  have hmv : Matrix.mulVec (galMat W l e σ) (Pi.single j (1 : ℤ_[l])) i
      = galMat W l e σ i j := by
    simp [Matrix.mulVec, dotProduct, Pi.single_apply]
  have hone : (Pi.single j (1 : ℤ_[l]) : Fin 2 → ℤ_[l]) i
      = (1 : Matrix (Fin 2) (Fin 2) ℤ_[l]) i j := by
    simp [Pi.single_apply, Matrix.one_apply]
  refine ⟨u i, ?_⟩
  have := hcomp
  simp only [Pi.sub_apply, Pi.smul_apply, hmv, hone] at this
  rw [this]
  push_cast
  ring

/-! ## ★★★★★★★★★★★★局所定数性 -/

/-- ★★★★★★★★★★★★
**`galMat` は開部分群の剰余類の上で `mod l^n` 一定である**——★**無条件**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`galMat (σ₀·τ) = galMat σ₀ · galMat τ` と `galMat τ ≡ 1 (mod l^n)`（`§9-1195`、第 769）から

    galMat (σ₀·τ) − galMat σ₀ = galMat σ₀ · (galMat τ − 1) ≡ 0  (mod l^n)

★★これが `galRep` の連続性の**位相の側の核**である
——`galMat` は局所定数（`mod l^n` の意味で）である。 -/
theorem galMat_sub_dvd_of_fix (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime] (n : ℕ)
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ₀ τ : L ≃ₐ[K] L)
    (hfix : ∀ P ∈ torsionPoints (W.baseChange L) (l ^ n), galPoint W τ P = P) (i j : Fin 2) :
    ((l : ℤ_[l]) ^ n) ∣ (galMat W l e (σ₀ * τ) i j - galMat W l e σ₀ i j) := by
  classical
  -- `galMat τ − 1` の各成分は `l^n` で割れる
  have hdvd := galMat_sub_one_dvd W l n e τ hfix
  choose N hN using hdvd
  have hmul : galMat W l e (σ₀ * τ) = galMat W l e σ₀ * galMat W l e τ :=
    galMat_mul W l e σ₀ τ
  have hexp : galMat W l e (σ₀ * τ) i j - galMat W l e σ₀ i j
      = ∑ k : Fin 2, galMat W l e σ₀ i k
          * (galMat W l e τ k j - (1 : Matrix (Fin 2) (Fin 2) ℤ_[l]) k j) := by
    rw [hmul]
    simp only [Matrix.mul_apply, mul_sub, Finset.sum_sub_distrib]
    congr 1
    rw [← Matrix.mul_apply]
    simp
  refine ⟨∑ k : Fin 2, galMat W l e σ₀ i k * N k j, ?_⟩
  rw [hexp, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun k _ => ?_)
  rw [hN k j]
  ring

/-! ## ★★★★★★★★★★★★★★`galMat` の各成分は連続 -/

/-- ★★★★★★★★★★★★★★★★
**`galMat` の各成分は連続**——★**無条件**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★機構: `ε > 0` に対し `l^{-n} < ε` なる `n` を取り、
`E[l^n]` の座標が生成する有限次中間体 `F` の `fixingSubgroup` の左剰余類 `σ₀·V`
（Krull 位相で開）を近傍に取る。`§9-1196`（第 770）により
`σ ∈ σ₀·V` なら `l^n ∣ galMat σ − galMat σ₀`、すなわち距離は `l^{-n} < ε`。 -/
theorem galMat_entry_continuous [Algebra.IsAlgebraic K L]
    (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (i j : Fin 2)
    (hfin : ∀ n : ℕ,
      (torsionPoints (W.baseChange L) (l ^ n) : Set ((W.baseChange L).toAffine.Point)).Finite) :
    Continuous (fun σ : L ≃ₐ[K] L => galMat W l e σ i j) := by
  rw [continuous_iff_continuousAt]
  intro σ₀
  rw [ContinuousAt, Metric.tendsto_nhds]
  intro ε hε
  -- `l^{-n} < ε` なる `n`
  have hl2 : (2 : ℝ) ≤ (l : ℝ) := by exact_mod_cast (Fact.out (p := l.Prime)).two_le
  have hl1 : (1 : ℝ) < (l : ℝ) := by linarith
  have hlinv : (l : ℝ)⁻¹ < 1 := inv_lt_one_of_one_lt₀ hl1
  have hlinv0 : (0 : ℝ) ≤ (l : ℝ)⁻¹ := by positivity
  obtain ⟨n, hn⟩ := exists_pow_lt_of_lt_one hε hlinv
  -- 有限次中間体を取る
  obtain ⟨F, hFfin, hFfix⟩ :=
    exists_finiteDimensional_fixing_torsion (K := K) (L := L) W (l ^ n) (hfin n)
  haveI := hFfin
  have hopen : IsOpen ((σ₀ : L ≃ₐ[K] L) • (F.fixingSubgroup : Set (L ≃ₐ[K] L))) :=
    (IntermediateField.fixingSubgroup_isOpen F).leftCoset σ₀
  have hmem : σ₀ ∈ (σ₀ : L ≃ₐ[K] L) • (F.fixingSubgroup : Set (L ≃ₐ[K] L)) := by
    refine ⟨1, ?_, by simp⟩
    exact Subgroup.one_mem _
  filter_upwards [hopen.mem_nhds hmem] with σ hσ
  obtain ⟨τ, hτ, hστ⟩ := hσ
  have hdvd := galMat_sub_dvd_of_fix W l n e σ₀ τ (hFfix τ hτ) i j
  have hσeq : σ = σ₀ * τ := hστ.symm
  rw [dist_eq_norm, hσeq]
  calc ‖galMat W l e (σ₀ * τ) i j - galMat W l e σ₀ i j‖
      ≤ ((l : ℝ)) ^ (-n : ℤ) :=
        (PadicInt.norm_le_pow_iff_mem_span_pow _ n).2 (Ideal.mem_span_singleton.2 hdvd)
    _ = ((l : ℝ)⁻¹) ^ n := by rw [zpow_neg, zpow_natCast, inv_pow]
    _ < ε := hn

/-! ## ★出典の紐付け(`.src`) -/

def galMat_entry_continuous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galMat の各成分は連続。★無条件)",
    sectionId := "genell-thm-3-8" }

def galMat_sub_dvd_of_fix.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galMat は開部分群の剰余類の上で mod l^n 一定。★無条件)",
    sectionId := "genell-thm-3-8" }

def galMat_sub_one_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(E[l^n] を固定するなら galMat ≡ 1 (mod l^n)。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_finiteDimensional_fixing_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(捩れ点の座標が生成する有限拡大——その固定部分群は捩れに自明に作用)",
    sectionId := "genell-thm-3-8" }

def galPoint_eq_of_fixed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(座標を固定する σ は点を動かさない。★無条件)",
    sectionId := "genell-thm-3-8" }

def galPoint_eq_of_fixed.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆連続性の残り 2 段: (2) L(E[l^n]) は L の有限拡大でありそれを固定する部分群は " ++
       "Krull 位相で開、(3) したがって galMat σ ≡ 1 (mod l^n) は開集合の上で成り立つ") 6 ]

end ABC3.Found.GenEll
