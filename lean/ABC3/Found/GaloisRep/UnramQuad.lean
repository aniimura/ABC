/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.BadPrimeData
import Mathlib.FieldTheory.KummerExtension
import Mathlib.RingTheory.AdicCompletion.Basic
import Mathlib.RingTheory.AdjoinRoot
import Mathlib.RingTheory.AdicCompletion.LocalRing
import Mathlib.RingTheory.Polynomial.IrreducibleRing
import Mathlib.RingTheory.DiscreteValuationRing.TFAE

/-!
# 第 1007 ブロック —— **★★★★★★★★不分岐 2 次拡大の葉**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か

第 1005／1006 に残った 2 本の仮説のうち、`hsplit`（完備化で**分裂**乗法還元）は
**`p ∣ 2` で非分裂の場合**だけが未処理である。

☆古典的な処理は「不分岐 2 次拡大に上げれば分裂になる」であり、
第 992／993 が捻り `d`（完備化の整数環の**単元**）を与えているので、
上げる先は **`Lv(√d)`** でよい。

★★その第一歩が本ファイルである——`d` が `R` で平方でなければ
**`K` でも平方でない**（したがって `X² − d` は既約で、`K(√d)` は 2 次拡大）。

☆証明は付値環の二分法だけである: `a² = d` なら `a` か `a⁻¹` が `R` の元で、
どちらの場合も `d` の平方根が `R` に取れてしまう。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve

/-! ## ★★★★★★★★`R` で平方でなければ `K` でも平方でない -/

/-- ★★★★★★★★**離散付値環で平方でない元は、その商体でも平方でない**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1007）**——不分岐 2 次拡大の最初の葉である。
☆`a² = d` とすると、付値環の二分法より `a ∈ R` か `a⁻¹ ∈ R`。

* `a = b ∈ R` なら `b² = d` でそのまま矛盾。
* `a⁻¹ = c ∈ R` なら `c²d = 1` なので `(cd)² = (c²d)d = d` でやはり矛盾。

★これで `X² − d` が `K[X]` で既約（根を持たない 2 次式）であることが言え、
`K(√d)` が 2 次拡大であることの土台になる。 -/
theorem not_isSquare_in_fractionField {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (d : R) (hns : ∀ b : R, b * b ≠ d) (a : K) :
    a * a ≠ algebraMap R K d := by
  intro hcon
  have hd0 : d ≠ 0 := fun h => hns 0 (by simp [h])
  have hdne : algebraMap R K d ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R K)).2 hd0
  have hane : a ≠ 0 := by
    intro h; rw [h, zero_mul] at hcon; exact hdne hcon.symm
  rcases ValuationRing.isInteger_or_isInteger R a with ⟨b, hb⟩ | ⟨c, hc⟩
  · exact hns b (IsFractionRing.injective R K (by rw [map_mul, hb, hcon]))
  · have hcdR : c * c * d = 1 := by
      refine IsFractionRing.injective R K ?_
      rw [map_mul, map_mul, hc, ← hcon, map_one]
      field_simp
    refine hns (c * d) ?_
    calc c * d * (c * d) = (c * c * d) * d := by ring
      _ = d := by rw [hcdR, one_mul]

/-! ## ★★★★★★★★`X² − d` は既約 -/

open Polynomial in
/-- ★★★★★★★★**`d` が `R` で平方でなければ `X² − d` は `K[X]` で既約**。

★★★★**2026-09-01（第 1008）**——mathlib の
`X_pow_sub_C_irreducible_of_prime`（`p = 2`）に第 1007 を当てるだけである。
☆これで `AdjoinRoot (X² − d)` が**体**になり、`K(√d)` が建つ。 -/
theorem irreducible_X_sq_sub_C_fractionField {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (d : R) (hns : ∀ b : R, b * b ≠ d) :
    Irreducible (X ^ 2 - C (algebraMap R K d)) :=
  X_pow_sub_C_irreducible_of_prime Nat.prime_two
    (fun b hb => not_isSquare_in_fractionField d hns b (by rw [← hb]; ring))

/-! ## ★★★★★★★★★★★★第 1010-1012 —— `R[√d]` は `I`-進完備

★不分岐 2 次拡大の整数環 `R′ = R[X]/(f)` は `R` 上**階数 2 の自由加群**である。
☆したがって `R` が `I`-進完備なら `R′` も `I`-進完備である。

★mathlib には**積の `IsAdicComplete` が無い**（2026-09-01 実測）ので、
`IsHausdorff`／`IsPrecomplete` を成分ごとに作るところから積む。 -/

section AdicProd

variable {R : Type} [CommRing R]

/-- ★★★★**`I • ⊤` は積と可換する**。 -/
theorem smul_top_prod (I : Ideal R) {M N : Type}
    [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N] :
    (I • (⊤ : Submodule R (M × N)))
      = (I • (⊤ : Submodule R M)).prod (I • (⊤ : Submodule R N)) := by
  refine le_antisymm (Submodule.smul_le.2 (fun r hr x _ => ?_)) ?_
  · exact ⟨Submodule.smul_mem_smul hr trivial, Submodule.smul_mem_smul hr trivial⟩
  · rintro ⟨x, y⟩ ⟨hx, hy⟩
    have h1 : ((x, 0) : M × N) ∈ I • (⊤ : Submodule R (M × N)) := by
      have he := Submodule.map_smul'' I (⊤ : Submodule R M) (LinearMap.inl R M N)
      have hmem : ((LinearMap.inl R M N) x)
          ∈ (I • (⊤ : Submodule R M)).map (LinearMap.inl R M N) :=
        Submodule.mem_map_of_mem hx
      rw [he] at hmem
      exact Submodule.smul_mono (fun ⦃_⦄ a => a) le_top hmem
    have h2 : ((0, y) : M × N) ∈ I • (⊤ : Submodule R (M × N)) := by
      have he := Submodule.map_smul'' I (⊤ : Submodule R N) (LinearMap.inr R M N)
      have hmem : ((LinearMap.inr R M N) y)
          ∈ (I • (⊤ : Submodule R N)).map (LinearMap.inr R M N) :=
        Submodule.mem_map_of_mem hy
      rw [he] at hmem
      exact Submodule.smul_mono (fun ⦃_⦄ a => a) le_top hmem
    have hxy : ((x, y) : M × N) = (x, 0) + (0, y) := by simp
    rw [hxy]
    exact Submodule.add_mem _ h1 h2

/-- ★★★★**`IsHausdorff` は積で保たれる**。 -/
theorem isHausdorff_prod (I : Ideal R) {M N : Type}
    [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
    [IsHausdorff I M] [IsHausdorff I N] : IsHausdorff I (M × N) where
  haus' := by
    rintro ⟨x, y⟩ h
    have hmem : ∀ n : ℕ, x ∈ I ^ n • (⊤ : Submodule R M)
        ∧ y ∈ I ^ n • (⊤ : Submodule R N) := by
      intro n
      have h1 := (SModEq.sub_mem).1 (h n)
      rw [sub_zero, smul_top_prod] at h1
      exact h1
    have hx := IsHausdorff.haus' (I := I) x
      (fun n => (SModEq.sub_mem).2 (by rw [sub_zero]; exact (hmem n).1))
    have hy := IsHausdorff.haus' (I := I) y
      (fun n => (SModEq.sub_mem).2 (by rw [sub_zero]; exact (hmem n).2))
    rw [hx, hy]
    rfl

/-- ★★★★**`IsPrecomplete` は積で保たれる**。 -/
theorem isPrecomplete_prod (I : Ideal R) {M N : Type}
    [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
    [IsPrecomplete I M] [IsPrecomplete I N] : IsPrecomplete I (M × N) where
  prec' := by
    intro f hf
    have hsplit : ∀ {m n : ℕ}, m ≤ n →
        ((f m).1 - (f n).1 ∈ I ^ m • (⊤ : Submodule R M)
          ∧ (f m).2 - (f n).2 ∈ I ^ m • (⊤ : Submodule R N)) := by
      intro m n hmn
      have h1 := (SModEq.sub_mem).1 (hf hmn)
      rw [smul_top_prod] at h1
      exact h1
    obtain ⟨L1, hL1⟩ := IsPrecomplete.prec' (I := I) (fun n => (f n).1)
      (fun {m n} hmn => (SModEq.sub_mem).2 (hsplit hmn).1)
    obtain ⟨L2, hL2⟩ := IsPrecomplete.prec' (I := I) (fun n => (f n).2)
      (fun {m n} hmn => (SModEq.sub_mem).2 (hsplit hmn).2)
    refine ⟨(L1, L2), fun n => (SModEq.sub_mem).2 ?_⟩
    rw [smul_top_prod]
    exact ⟨(SModEq.sub_mem).1 (hL1 n), (SModEq.sub_mem).1 (hL2 n)⟩

/-- ★★★★★★★★**`IsAdicComplete` は積で保たれる**（mathlib に無い、2026-09-01 実測）。 -/
theorem isAdicComplete_prod (I : Ideal R) {M N : Type}
    [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
    [IsAdicComplete I M] [IsAdicComplete I N] : IsAdicComplete I (M × N) :=
  { toIsHausdorff := isHausdorff_prod I, toIsPrecomplete := isPrecomplete_prod I }

/-- ★★★★**`IsHausdorff` は線型同型で移る**。 -/
theorem isHausdorff_of_linearEquiv (I : Ideal R) {M N : Type}
    [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
    (e : M ≃ₗ[R] N) [IsHausdorff I M] : IsHausdorff I N where
  haus' := by
    intro y hy
    have hx : ∀ n : ℕ, (e.symm y) ≡ 0 [SMOD I ^ n • (⊤ : Submodule R M)] := by
      intro n
      refine (SModEq.sub_mem).2 ?_
      rw [sub_zero]
      have h1 := (SModEq.sub_mem).1 (hy n)
      rw [sub_zero] at h1
      have h2 : (I ^ n • (⊤ : Submodule R M)).map (e : M →ₗ[R] N)
          = I ^ n • (⊤ : Submodule R N) := by
        rw [Submodule.map_smul'', Submodule.map_top, LinearEquiv.range]
      have h3 : y ∈ (I ^ n • (⊤ : Submodule R M)).map (e : M →ₗ[R] N) := by
        rw [h2]; exact h1
      obtain ⟨x, hxm, hxe⟩ := h3
      have hxx : e.symm y = x := by rw [← hxe]; simp
      rw [hxx]; exact hxm
    have hz := IsHausdorff.haus' (I := I) (e.symm y) hx
    have hy0 : y = e (e.symm y) := by simp
    rw [hy0, hz, map_zero]

/-- ★★★★**`IsPrecomplete` は線型同型で移る**。 -/
theorem isPrecomplete_of_linearEquiv (I : Ideal R) {M N : Type}
    [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
    (e : M ≃ₗ[R] N) [IsPrecomplete I M] : IsPrecomplete I N where
  prec' := by
    intro f hf
    have hmapN : ∀ n : ℕ, (I ^ n • (⊤ : Submodule R M)).map (e : M →ₗ[R] N)
        = I ^ n • (⊤ : Submodule R N) := by
      intro n; rw [Submodule.map_smul'', Submodule.map_top, LinearEquiv.range]
    have hmapM : ∀ n : ℕ, (I ^ n • (⊤ : Submodule R N)).map (e.symm : N →ₗ[R] M)
        = I ^ n • (⊤ : Submodule R M) := by
      intro n; rw [Submodule.map_smul'', Submodule.map_top, LinearEquiv.range]
    obtain ⟨L, hL⟩ := IsPrecomplete.prec' (I := I) (fun n => e.symm (f n)) (by
      intro m n hmn
      refine (SModEq.sub_mem).2 ?_
      have h1 := (SModEq.sub_mem).1 (hf hmn)
      have hsub : e.symm (f m) - e.symm (f n) = e.symm (f m - f n) := by simp
      rw [hsub, ← hmapM m]
      exact Submodule.mem_map_of_mem h1)
    refine ⟨e L, fun n => (SModEq.sub_mem).2 ?_⟩
    have h1 := (SModEq.sub_mem).1 (hL n)
    have h2 : f n - e L = e (e.symm (f n) - L) := by simp
    rw [h2, ← hmapN n]
    exact Submodule.mem_map_of_mem h1

/-- ★★★★★★★★**`IsAdicComplete` は線型同型で移る**。 -/
theorem isAdicComplete_of_linearEquiv (I : Ideal R) {M N : Type}
    [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
    (e : M ≃ₗ[R] N) [IsAdicComplete I M] : IsAdicComplete I N :=
  { toIsHausdorff := isHausdorff_of_linearEquiv I e,
    toIsPrecomplete := isPrecomplete_of_linearEquiv I e }

open Polynomial in
/-- ★★★★★★**2 次のモニック多項式の剰余環は `R × R` と線型同型**。 -/
noncomputable def adjoinRootEquivProd {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2) :
    AdjoinRoot f ≃ₗ[R] (R × R) := by
  have hb := (AdjoinRoot.powerBasis' hf).basis
  rw [AdjoinRoot.powerBasis'_dim hf, hdeg] at hb
  exact (hb.repr.trans (Finsupp.linearEquivFunOnFinite R R (Fin 2))).trans
    (LinearEquiv.finTwoArrow R R)

open Polynomial in
/-- ★★★★★★★★★★★★**`R[X]/(f)` は `I`-進完備**（`f` は 2 次のモニック）。

★★★★**2026-09-01（第 1012）**——不分岐 2 次拡大の整数環が完備であることの中身。
☆`R × R` に落として第 1010 を当てるだけである。 -/
theorem isAdicComplete_adjoinRoot (I : Ideal R) [IsAdicComplete I R]
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2) :
    IsAdicComplete I (AdjoinRoot f) :=
  haveI : IsAdicComplete I (R × R) := isAdicComplete_prod I
  isAdicComplete_of_linearEquiv I (adjoinRootEquivProd hf hdeg).symm

end AdicProd

/-! ## ★★★★★★★★★★★★第 1013-1015 —— `R[X]/(f)` は局所環

★第 1012 は `IsAdicComplete I (AdjoinRoot f)` を **`R`-加群として**与える。
☆mathlib の `isLocalRing_of_isAdicComplete_maximal` が欲しいのは
**`AdjoinRoot f` 自身のイデアル**に対する完備性なので、その橋が要る（第 1013）。

★極大性は `R′/𝔪R′ ≅ k[X]/(f̄)` であり、`f̄` が既約なら体である（第 1014）。
☆両方を合わせれば `R′` は局所環である（第 1015）。 -/

section LocalAdjoinRoot

variable {R : Type} [CommRing R]

/-- ★★★★**`I^n • ⊤`（`R`-加群）と `(I.map φ)^n • ⊤`（`R′`-加群）は同じ集合**。 -/
theorem mem_smul_top_iff_map {R' : Type} [CommRing R'] [Algebra R R']
    (I : Ideal R) (n : ℕ) (x : R') :
    x ∈ (I.map (algebraMap R R')) ^ n • (⊤ : Submodule R' R')
      ↔ x ∈ I ^ n • (⊤ : Submodule R R') := by
  have h1 : ((I.map (algebraMap R R')) ^ n) • (⊤ : Submodule R' R')
      = (I ^ n).map (algebraMap R R') := by
    rw [← Ideal.map_pow, Ideal.smul_top_eq_map (S := R') ((I ^ n).map (algebraMap R R'))]
    ext y
    simp [Ideal.map_id]
  rw [h1, Ideal.smul_top_eq_map (S := R') (I ^ n)]
  rfl

theorem isHausdorff_map_algebraMap {R' : Type} [CommRing R'] [Algebra R R']
    (I : Ideal R) [IsHausdorff I R'] : IsHausdorff (I.map (algebraMap R R')) R' where
  haus' := by
    intro x hx
    refine IsHausdorff.haus' (I := I) (M := R') x (fun n => (SModEq.sub_mem).2 ?_)
    exact (mem_smul_top_iff_map I n _).1 ((SModEq.sub_mem).1 (hx n))

theorem isPrecomplete_map_algebraMap {R' : Type} [CommRing R'] [Algebra R R']
    (I : Ideal R) [IsPrecomplete I R'] : IsPrecomplete (I.map (algebraMap R R')) R' where
  prec' := by
    intro f hf
    obtain ⟨L, hL⟩ := IsPrecomplete.prec' (I := I) (M := R') f (by
      intro m n hmn
      exact (SModEq.sub_mem).2 ((mem_smul_top_iff_map I m _).1 ((SModEq.sub_mem).1 (hf hmn))))
    exact ⟨L, fun n => (SModEq.sub_mem).2
      ((mem_smul_top_iff_map I n _).2 ((SModEq.sub_mem).1 (hL n)))⟩

/-- ★★★★★★★★**`R`-加群としての `I`-進完備性は
`R′` のイデアル `I.map φ` に対する完備性に移る**（第 1013）。 -/
theorem isAdicComplete_map_algebraMap {R' : Type} [CommRing R'] [Algebra R R']
    (I : Ideal R) [IsAdicComplete I R'] : IsAdicComplete (I.map (algebraMap R R')) R' :=
  { toIsHausdorff := isHausdorff_map_algebraMap I,
    toIsPrecomplete := isPrecomplete_map_algebraMap I }

open Polynomial in
/-- ★★★★★★★★**`f̄` が既約なら `I · R[X]/(f)` は極大**（第 1014）。

☆`AdjoinRoot f / I ≅ k[X]/(f̄)` であり、`k` が体なら `k[X]` は PID、
既約元が生成するイデアルは極大である。 -/
theorem isMaximal_map_adjoinRoot (I : Ideal R) [I.IsMaximal]
    {f : R[X]} (hirr : Irreducible (f.map (Ideal.Quotient.mk I))) :
    (I.map (AdjoinRoot.of f)).IsMaximal := by
  letI : Field (R ⧸ I) := Ideal.Quotient.field I
  haveI hmax : (Ideal.span {f.map (Ideal.Quotient.mk I)}).IsMaximal :=
    PrincipalIdealRing.isMaximal_of_irreducible hirr
  refine Ideal.Quotient.maximal_of_isField _ ?_
  refine MulEquiv.isField ?_ (AdjoinRoot.quotAdjoinRootEquivQuotPolynomialQuot I f).toMulEquiv
  exact (Ideal.Quotient.maximal_ideal_iff_isField_quotient _).1 hmax

open Polynomial in
/-- ★★★★★★★★★★★★**`R[X]/(f)` は局所環**（第 1015）。

★`f` は 2 次のモニックで、剰余体での還元 `f̄` が既約であるとする。
☆第 1012（完備）＋第 1013（橋）＋第 1014（極大）を
mathlib の `isLocalRing_of_isAdicComplete_maximal` に流す。 -/
theorem isLocalRing_adjoinRoot [IsLocalRing R]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R)))) :
    IsLocalRing (AdjoinRoot f) := by
  haveI : IsAdicComplete (IsLocalRing.maximalIdeal R) (AdjoinRoot f) :=
    isAdicComplete_adjoinRoot (IsLocalRing.maximalIdeal R) hf hdeg
  haveI hcomp : IsAdicComplete
      ((IsLocalRing.maximalIdeal R).map (algebraMap R (AdjoinRoot f))) (AdjoinRoot f) :=
    isAdicComplete_map_algebraMap _
  haveI hmax : ((IsLocalRing.maximalIdeal R).map (algebraMap R (AdjoinRoot f))).IsMaximal := by
    rw [AdjoinRoot.algebraMap_eq]
    exact isMaximal_map_adjoinRoot (IsLocalRing.maximalIdeal R) hirr
  exact isLocalRing_of_isAdicComplete_maximal
    ((IsLocalRing.maximalIdeal R).map (algebraMap R (AdjoinRoot f)))

end LocalAdjoinRoot

/-! ## ★★★★★★★★★★★★★★★★第 1016-1018 —— `R[X]/(f)` は完備離散付値環

★第 1015 で局所環になった。☆残るのは**整域**と**離散付値環**である。

| 段 | 道具 |
|---|---|
| `f̄` 既約 → `f` 既約 | mathlib `Monic.irreducible_of_irreducible_map_of_isPrime_nilradical` |
| `f` 既約 → `AdjoinRoot f` 整域 | `AdjoinRoot.isDomain_of_prime`（UFD で既約＝素） |
| 極大が `𝔪R′` で主 → 離散付値環 | mathlib `IsDiscreteValuationRing.TFAE` の (5) ⟹ (1) |

☆`𝔪R′` が主であることは `𝔪 = (π)` の像が `(φ π)` だからである。 -/

section DvrAdjoinRoot

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]

open Polynomial in
/-- ★★★★★★★★**剰余体で既約なモニック多項式は `R` でも既約**（第 1016）。 -/
theorem irreducible_of_irreducible_residue {f : R[X]} (hf : f.Monic)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R)))) :
    Irreducible f := by
  letI : Field (R ⧸ IsLocalRing.maximalIdeal R) :=
    Ideal.Quotient.field (IsLocalRing.maximalIdeal R)
  haveI hnil : (nilradical R).IsPrime := by
    have h0 : nilradical R = ⊥ := nilradical_eq_zero R
    rw [h0]; exact Ideal.isPrime_bot
  exact Polynomial.Monic.irreducible_of_irreducible_map_of_isPrime_nilradical
    (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R)) f hf hirr

open Polynomial in
/-- ★★★★★★★★**`R[X]/(f)` は整域**（第 1017）。 -/
theorem isDomain_adjoinRoot {f : R[X]} (hf : f.Monic)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R)))) :
    IsDomain (AdjoinRoot f) :=
  AdjoinRoot.isDomain_of_prime
    (UniqueFactorizationMonoid.irreducible_iff_prime.1
      (irreducible_of_irreducible_residue hf hirr))

open Polynomial in
/-- ★★★★★★★★★★★★★★★★**`R[X]/(f)` は離散付値環**（第 1018）。

★`f` は 2 次のモニックで、剰余体での還元 `f̄` が既約であるとする。
☆極大イデアルは `𝔪R′`（第 1014／1015）であり、`𝔪 = (π)` の像は `(φ π)` なので**主**。
★`IsDiscreteValuationRing.TFAE` の (5) ⟹ (1) がそのまま当たる。 -/
theorem isDiscreteValuationRing_adjoinRoot
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] :
    IsDiscreteValuationRing (AdjoinRoot f) := by
  haveI hloc : IsLocalRing (AdjoinRoot f) := isLocalRing_adjoinRoot hf hdeg hirr
  haveI : Module.Finite R (AdjoinRoot f) := Monic.finite_adjoinRoot hf
  haveI : IsNoetherianRing (AdjoinRoot f) := AdjoinRoot.instIsNoetherianRing
  have hmaxI : ((IsLocalRing.maximalIdeal R).map (algebraMap R (AdjoinRoot f))).IsMaximal := by
    rw [AdjoinRoot.algebraMap_eq]
    exact isMaximal_map_adjoinRoot _ hirr
  have hmeq : IsLocalRing.maximalIdeal (AdjoinRoot f)
      = (IsLocalRing.maximalIdeal R).map (algebraMap R (AdjoinRoot f)) :=
    (IsLocalRing.eq_maximalIdeal hmaxI).symm
  obtain ⟨π, hπ⟩ := (IsPrincipalIdealRing.principal (IsLocalRing.maximalIdeal R))
  have hinj : Function.Injective (algebraMap R (AdjoinRoot f)) := by
    rw [AdjoinRoot.algebraMap_eq]
    refine AdjoinRoot.of.injective_of_degree_ne_zero ?_
    rw [Polynomial.degree_eq_natDegree hf.ne_zero, hdeg]
    simp
  have hπ0 : π ≠ 0 := by
    intro h
    rw [h] at hπ
    have hb : IsLocalRing.maximalIdeal R = ⊥ := by rw [hπ]; simp
    exact (IsDiscreteValuationRing.not_isField R)
      (IsLocalRing.isField_iff_maximalIdeal_eq.2 hb)
  have hprin : Submodule.IsPrincipal (IsLocalRing.maximalIdeal (AdjoinRoot f)) := by
    rw [hmeq, hπ, Ideal.map_span]
    simp only [Set.image_singleton]
    exact ⟨⟨algebraMap R (AdjoinRoot f) π, rfl⟩⟩
  have hnf : ¬ IsField (AdjoinRoot f) := by
    intro hfield
    have h1 : IsLocalRing.maximalIdeal (AdjoinRoot f) = ⊥ :=
      IsLocalRing.isField_iff_maximalIdeal_eq.1 hfield
    have h2 : algebraMap R (AdjoinRoot f) π ∈ IsLocalRing.maximalIdeal (AdjoinRoot f) := by
      rw [hmeq]
      exact Ideal.mem_map_of_mem _ (by rw [hπ]; exact Ideal.subset_span rfl)
    rw [h1, Ideal.mem_bot] at h2
    exact hπ0 (hinj (by rw [h2, map_zero]))
  exact ((IsDiscreteValuationRing.TFAE (AdjoinRoot f) hnf).out 4 0).1 hprin

end DvrAdjoinRoot

/-! ## ★出典の紐付け(`.src`) -/

def not_isSquare_in_fractionField.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(離散付値環で平方でない元は商体でも平方でない。★無条件)",
    sectionId := "genell-lemma-3-5" }

def irreducible_X_sq_sub_C_fractionField.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(d が平方でなければ X² − d は既約。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isAdicComplete_prod.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(IsAdicComplete は積で保たれる——mathlib に無い。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isAdicComplete_of_linearEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(IsAdicComplete は線型同型で移る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isAdicComplete_adjoinRoot.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2 次のモニック多項式の剰余環は I-進完備。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
