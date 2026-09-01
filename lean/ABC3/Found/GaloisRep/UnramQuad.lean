/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Found.GenEll.QuadTwist
import ABC3.Found.GaloisRep.CompletionValuationBridge
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

/-! ## ★★★★★★★★★★★★第 1021 —— 商体の側を建てる

★環の側（完備・局所・整域・離散付値環）は第 1012-1018 で建った。
☆商体は mathlib の `FractionRing (AdjoinRoot f)` をそのまま使えばよく、
`IsFractionRing` は自動でつく。

★要るのは **`K → Frac (AdjoinRoot f)`** の埋め込みである。
☆`R → AdjoinRoot f → Frac (AdjoinRoot f)` が単射なので
`IsFractionRing.lift` で `K` から延びる。 -/

section QuadField

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
variable {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

open Polynomial in
/-- ★★★★★★★★★★★★**`K` から `Frac (R[X]/(f))` への埋め込み**（第 1021）。

★`f` は 2 次のモニックなので `R → R[X]/(f)` は単射、
`R[X]/(f) → Frac` も単射（整域）。☆その合成を `K` に延ばす。 -/
noncomputable def quadFieldHom {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    [IsDomain (AdjoinRoot f)] : K →+* FractionRing (AdjoinRoot f) :=
  IsFractionRing.lift (A := R) (K := K)
    (g := (algebraMap (AdjoinRoot f) (FractionRing (AdjoinRoot f))).comp
      (algebraMap R (AdjoinRoot f)))
    (by
      rw [RingHom.coe_comp]
      refine Function.Injective.comp
        (IsFractionRing.injective (AdjoinRoot f) (FractionRing (AdjoinRoot f))) ?_
      rw [AdjoinRoot.algebraMap_eq]
      refine AdjoinRoot.of.injective_of_degree_ne_zero ?_
      rw [Polynomial.degree_eq_natDegree hf.ne_zero, hdeg]
      simp)

open Polynomial in
/-- ★★★★**埋め込みは `R` の上で `R → R[X]/(f) → Frac` と一致する**。 -/
theorem quadFieldHom_algebraMap {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    [IsDomain (AdjoinRoot f)] (x : R) :
    quadFieldHom (K := K) hf hdeg (algebraMap R K x)
      = algebraMap (AdjoinRoot f) (FractionRing (AdjoinRoot f)) (algebraMap R (AdjoinRoot f) x) :=
  IsFractionRing.lift_algebraMap _ x

open Polynomial in
/-- ★★★★**体からの環準同型なので単射**。 -/
theorem quadFieldHom_injective {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    [IsDomain (AdjoinRoot f)] :
    Function.Injective (quadFieldHom (K := K) hf hdeg) :=
  (quadFieldHom (K := K) hf hdeg).injective

end QuadField

/-! ## ★★★★★★★★★★★★★★★★第 1022-1024 —— 付値の橋（`e = 1`）

★不分岐 2 次拡大の要点は「`R` の素元 `π` が `R′` でも素元」＝ **`e = 1`** である。
☆第 1018 の証明で `maximalIdeal R′ = 𝔪R′` を示したので、
`𝔪 = (π)` なら `maximalIdeal R′ = (φ π)` である（第 1022）。

★★したがって `K` 上の 2 つの付値——`R` の付値と、`R′` の付値を `K` に引き戻したもの——は
**同値であり、両方全射**なので等しい（第 964 の `eq_of_isEquiv_of_surjective`）。

☆同値性の片側は `IsDiscreteValuationRing.exists_lift_of_le_one`、
もう片側は「`R′` の元で `K` に入るものは `R` の元」＝ `R` が整閉であること（第 1023）。 -/

section ValuationBridge

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
variable {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

open Polynomial IsDedekindDomain in
/-- ★★★★**`R[X]/(f)` の極大イデアルは `𝔪` の像**（第 1022）。 -/
theorem maximalIdeal_adjoinRoot_eq_map [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsLocalRing (AdjoinRoot f)] :
    IsLocalRing.maximalIdeal (AdjoinRoot f)
      = (IsLocalRing.maximalIdeal R).map (algebraMap R (AdjoinRoot f)) := by
  have hmaxI : ((IsLocalRing.maximalIdeal R).map (algebraMap R (AdjoinRoot f))).IsMaximal := by
    rw [AdjoinRoot.algebraMap_eq]
    exact isMaximal_map_adjoinRoot _ hirr
  exact (IsLocalRing.eq_maximalIdeal hmaxI).symm

open Polynomial IsDedekindDomain in
/-- ★★★★★★★★**`π` は `R[X]/(f)` でも素元**——これが `e = 1` である（第 1022）。 -/
theorem maximalIdeal_adjoinRoot_eq_span [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsLocalRing (AdjoinRoot f)] {π : R}
    (hπ : IsLocalRing.maximalIdeal R = Ideal.span {π}) :
    IsLocalRing.maximalIdeal (AdjoinRoot f)
      = Ideal.span {algebraMap R (AdjoinRoot f) π} := by
  rw [maximalIdeal_adjoinRoot_eq_map hf hdeg hirr, hπ, Ideal.map_span]
  simp

open Polynomial IsDedekindDomain in
/-- ★★★★**`R → R[X]/(f)` は単射**（`f` は 2 次のモニック）。 -/
theorem algebraMap_adjoinRoot_injective {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2) :
    Function.Injective (algebraMap R (AdjoinRoot f)) := by
  rw [AdjoinRoot.algebraMap_eq]
  refine AdjoinRoot.of.injective_of_degree_ne_zero ?_
  rw [Polynomial.degree_eq_natDegree hf.ne_zero, hdeg]
  simp

open Polynomial IsDedekindDomain in
/-- ★★★★★★★★★★★★**`R′` の付値で `≤ 1` なら `K` の元は `R` の元**（第 1023）。

☆`R′` は `R` 上有限なので整、`R` は整閉なので降りる。 -/
theorem exists_lift_of_quadFieldHom_le_one
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)] (y : K)
    (hy : (HeightOneSpectrum.valuation (FractionRing (AdjoinRoot f))
      (IsDiscreteValuationRing.maximalIdeal (AdjoinRoot f))) (quadFieldHom hf hdeg y) ≤ 1) :
    ∃ a : R, algebraMap R K a = y := by
  obtain ⟨b, hb⟩ := IsDiscreteValuationRing.exists_lift_of_le_one hy
  haveI : Module.Finite R (AdjoinRoot f) := Monic.finite_adjoinRoot hf
  have hbint : IsIntegral R b := IsIntegral.of_finite R b
  letI : Algebra K (FractionRing (AdjoinRoot f)) := (quadFieldHom hf hdeg).toAlgebra
  haveI : IsScalarTower R K (FractionRing (AdjoinRoot f)) :=
    IsScalarTower.of_algebraMap_eq (fun x => (quadFieldHom_algebraMap hf hdeg x).symm)
  have hint : IsIntegral R (algebraMap K (FractionRing (AdjoinRoot f)) y) := by
    rw [show algebraMap K (FractionRing (AdjoinRoot f)) y = quadFieldHom hf hdeg y from rfl, ← hb]
    exact hbint.map (IsScalarTower.toAlgHom R (AdjoinRoot f) (FractionRing (AdjoinRoot f)))
  exact IsIntegrallyClosed.isIntegral_iff.1 (IsIntegral.tower_bot_of_field hint)

open Polynomial IsDedekindDomain in
/-- ★★★★★★★★**`π` の像の付値は `exp(-1)`**——`e = 1` の付値版（第 1023）。 -/
theorem valuation_quadFieldHom_uniformizer [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)]
    {π : R} (hπ : IsLocalRing.maximalIdeal R = Ideal.span {π}) (hπ0 : π ≠ 0) :
    (HeightOneSpectrum.valuation (FractionRing (AdjoinRoot f))
        (IsDiscreteValuationRing.maximalIdeal (AdjoinRoot f)))
      (quadFieldHom hf hdeg (algebraMap R K π)) = WithZero.exp (-1 : ℤ) := by
  rw [quadFieldHom_algebraMap hf hdeg, HeightOneSpectrum.valuation_of_algebraMap]
  refine HeightOneSpectrum.intValuation_singleton _ ?_ ?_
  · intro hc
    exact hπ0 (algebraMap_adjoinRoot_injective hf hdeg (by rw [hc, map_zero]))
  · show IsLocalRing.maximalIdeal (AdjoinRoot f) = _
    exact maximalIdeal_adjoinRoot_eq_span hf hdeg hirr hπ

open Polynomial IsDedekindDomain in
/-- ★★★★★★★★★★★★★★★★**付値の橋**——`K` の付値は拡大で変わらない（第 1024）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1024）**——これが不分岐（`e = 1`）の中身である。
☆第 964 の `eq_of_isEquiv_of_surjective` に、
同値性（第 1023）と全射性（`π` の像が素元）を流す。 -/
theorem valuation_quadFieldHom [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)] (y : K) :
    (HeightOneSpectrum.valuation (FractionRing (AdjoinRoot f))
        (IsDiscreteValuationRing.maximalIdeal (AdjoinRoot f))) (quadFieldHom hf hdeg y)
      = (HeightOneSpectrum.valuation K (IsDiscreteValuationRing.maximalIdeal R)) y := by
  obtain ⟨π, hπ⟩ := (IsPrincipalIdealRing.principal (IsLocalRing.maximalIdeal R))
  have hπ0 : π ≠ 0 := by
    intro h
    rw [h] at hπ
    have hb : IsLocalRing.maximalIdeal R = ⊥ := by rw [hπ]; simp
    exact (IsDiscreteValuationRing.not_isField R)
      (IsLocalRing.isField_iff_maximalIdeal_eq.2 hb)
  have hkey : (HeightOneSpectrum.valuation K (IsDiscreteValuationRing.maximalIdeal R))
      = (HeightOneSpectrum.valuation (FractionRing (AdjoinRoot f))
          (IsDiscreteValuationRing.maximalIdeal (AdjoinRoot f))).comap
        (quadFieldHom hf hdeg) := by
    refine eq_of_isEquiv_of_surjective ?_ ?_ ?_
    · refine Valuation.isEquiv_of_val_le_one (fun x => ?_)
      constructor
      · intro hx
        obtain ⟨a, ha⟩ := IsDiscreteValuationRing.exists_lift_of_le_one hx
        show (HeightOneSpectrum.valuation _ _) (quadFieldHom hf hdeg x) ≤ 1
        rw [← ha, quadFieldHom_algebraMap]
        exact HeightOneSpectrum.valuation_le_one _ _
      · intro hx
        obtain ⟨a, ha⟩ := exists_lift_of_quadFieldHom_le_one hf hdeg x hx
        rw [← ha]
        exact HeightOneSpectrum.valuation_le_one _ _
    · exact HeightOneSpectrum.valuation_surjective _ _
    · intro t
      rcases eq_or_ne t 0 with rfl | ht
      · exact ⟨0, by simp⟩
      · obtain ⟨n, rfl⟩ : ∃ n : ℤ, t = WithZero.exp n :=
          ⟨WithZero.log t, (WithZero.exp_log ht).symm⟩
        refine ⟨(algebraMap R K π) ^ (-n), ?_⟩
        show (HeightOneSpectrum.valuation _ _)
          (quadFieldHom hf hdeg ((algebraMap R K π) ^ (-n))) = _
        rw [map_zpow₀, map_zpow₀, valuation_quadFieldHom_uniformizer hf hdeg hirr hπ hπ0,
          exp_neg_one_zpow, neg_neg]
  rw [hkey]
  rfl

end ValuationBridge

/-! ## ★★★★★★★★★★★★第 1025 —— 剰余体では `f` は分裂する

★`R′ = R[X]/(f)` の剰余体には **`f` の根がある**（`root f` の像）。
☆2 次式は根を 1 つ持てば分裂するので、`f` は `ResidueField R′` で分裂する。

★★これが「非分裂だった 2 次式が拡大で分裂する」ことの中身である——
曲線の分裂性を決める 2 次式を（`c₄` で割って）モニックにし、
その持ち上げを `f` に取れば、この定理がそのまま `hsplit` を与える。 -/

section ResidueSplits

open Polynomial in
/-- ★★★★★★**根を 1 つ持つ 2 次式は分裂する**。 -/
theorem splits_of_isRoot_natDegree_two {k : Type} [Field k] {p : k[X]}
    (hdeg : p.natDegree = 2) {r : k} (hr : p.IsRoot r) : p.Splits := by
  have hmon : (X - C r).Monic := monic_X_sub_C r
  have hfac : (X - C r) * (p /ₘ (X - C r)) = p := mul_divByMonic_eq_iff_isRoot.2 hr
  have hp0 : p ≠ 0 := fun h => by simp [h] at hdeg
  have hq0 : p /ₘ (X - C r) ≠ 0 := by
    intro h
    rw [h, mul_zero] at hfac
    exact hp0 hfac.symm
  rw [← hfac, splits_mul hmon.ne_zero hq0]
  refine ⟨Splits.of_natDegree_le_one (by simp), Splits.of_natDegree_le_one ?_⟩
  rw [natDegree_divByMonic p hmon, natDegree_X_sub_C, hdeg]

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]

open Polynomial in
/-- ★★★★★★★★**`R[X]/(f)` の剰余体には `f` の根がある**——`root f` の像である。 -/
theorem exists_isRoot_residueField (f : R[X]) [IsLocalRing (AdjoinRoot f)] :
    ∃ r : IsLocalRing.ResidueField (AdjoinRoot f),
      (f.map ((IsLocalRing.residue (AdjoinRoot f)).comp
        (algebraMap R (AdjoinRoot f)))).IsRoot r := by
  refine ⟨IsLocalRing.residue (AdjoinRoot f) (AdjoinRoot.root f), ?_⟩
  show Polynomial.eval _ (f.map _) = 0
  rw [Polynomial.eval_map,
    ← Polynomial.hom_eval₂ f (algebraMap R (AdjoinRoot f)) (IsLocalRing.residue _)
      (AdjoinRoot.root f), AdjoinRoot.algebraMap_eq, AdjoinRoot.eval₂_root, map_zero]

open Polynomial in
/-- ★★★★★★★★★★★★**`f` は `R[X]/(f)` の剰余体で分裂する**（第 1025）。

★これが不分岐 2 次拡大の最後の段である——
非分裂だった 2 次式は、拡大の剰余体では根を持つので分裂する。 -/
theorem splits_map_residueField {f : R[X]} (hdeg : f.natDegree = 2)
    [IsLocalRing (AdjoinRoot f)]
    (hmap : (f.map ((IsLocalRing.residue (AdjoinRoot f)).comp
      (algebraMap R (AdjoinRoot f)))).natDegree = 2) :
    (f.map ((IsLocalRing.residue (AdjoinRoot f)).comp
      (algebraMap R (AdjoinRoot f)))).Splits := by
  obtain ⟨r, hr⟩ := exists_isRoot_residueField f
  exact splits_of_isRoot_natDegree_two hmap hr

end ResidueSplits

/-! ## ★★★★★★★★★★★★第 1029 —— 拡大の上での局所データ

★第 1028 が受ける `E ⊗ Lv′` の側の仮説のうち、
**楕円性は mathlib＋在庫（`isElliptic_baseChange`）で自動**である。
☆極大性と乗法還元は第 973／976 が `hp` を受ける形で書いてあるので、
拡大の付値の橋（第 1024）を渡せばそのまま出る。

★★残るのは **`Splits`（分裂性）ただ 1 つ**になる。 -/

section ExtLocalData

variable {L : Type} [Field L] [NumberField L]
variable {Lv : Type} [Field Lv] [Algebra L Lv]
variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
variable [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]

open Polynomial IsDedekindDomain in
/-- ★★★★★★★★★★★★**拡大の上でも付値の橋は成り立つ**（第 1029）。

☆`hp`（`Lv` の側）と第 1024（`Lv → Lv′` の側）を合成するだけである。 -/
theorem valuation_algebraMap_ext (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom hf hdeg (algebraMap L Lv x))
    (x : L) :
    (HeightOneSpectrum.valuation (FractionRing (AdjoinRoot f))
        (IsDiscreteValuationRing.maximalIdeal (AdjoinRoot f)))
      (algebraMap L (FractionRing (AdjoinRoot f)) x)
      = (HeightOneSpectrum.valuation L p) x := by
  rw [halg, valuation_quadFieldHom hf hdeg hirr, hp]

open Polynomial IsDedekindDomain WeierstrassCurve in
/-- ★★★★★★★★**極小性は拡大にも移る**（第 1029）。

☆第 973 が `hp` を受ける形なので、第 1029 の橋を渡すだけ。 -/
theorem isMinimal_baseChange_ext (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom hf hdeg (algebraMap L Lv x))
    (E : WeierstrassCurve L) [E.IsElliptic]
    (C : WeierstrassCurve.VariableChange L)
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0) :
    WeierstrassCurve.IsMinimal (AdjoinRoot f)
      ((C • E).baseChange (FractionRing (AdjoinRoot f))) :=
  isMinimal_baseChange_at_bad_prime p
    (valuation_algebraMap_ext p hp hf hdeg hirr halg) E C hC hc4ne hc4

open Polynomial IsDedekindDomain WeierstrassCurve in
/-- ★★★★★★★★**乗法還元も拡大に移る**（第 1029）。

☆第 976 が `hp` を受ける形なので、第 1029 の橋を渡すだけ。 -/
theorem hasMultiplicativeReduction_ext (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom hf hdeg (algebraMap L Lv x))
    (E : WeierstrassCurve L) [E.IsElliptic]
    (C : WeierstrassCurve.VariableChange L)
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (hj : jExp p E < 0)
    [WeierstrassCurve.IsMinimal (AdjoinRoot f)
      ((C • E).baseChange (FractionRing (AdjoinRoot f)))] :
    WeierstrassCurve.HasMultiplicativeReduction (AdjoinRoot f)
      ((C • E).baseChange (FractionRing (AdjoinRoot f))) :=
  hasMultiplicativeReduction_at_bad_prime p
    (valuation_algebraMap_ext p hp hf hdeg hirr halg) E C hC hc4ne hc4 hj

end ExtLocalData

/-! ## ★★★★★★★★★★★★第 1030 —— 整モデルは拡大で係数ごとに移る

★`HasSplitMultiplicativeReduction` の `Splits` 条件は
**整モデルの係数**で書かれている。
☆したがって `integralModel R′ (W ⊗ Lv′)` が
`(integralModel R (W ⊗ Lv)).map φ` であることが要る。

★mathlib の `integralModel` は `.choose` で定義されているが、
第 925 の `integralModel_eq_of_map_eq`（底変換が一致すれば一意）で決まる。 -/

section IntegralModelExt

open Polynomial IsDedekindDomain WeierstrassCurve in
/-- ★★★★★★★★★★★★**整モデルは拡大で係数ごとに移る**（第 1030）。 -/
theorem integralModel_ext {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv]
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    [IsDomain (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom hf hdeg (algebraMap L Lv x))
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral R (W.baseChange Lv)]
    [WeierstrassCurve.IsIntegral (AdjoinRoot f)
      (W.baseChange (FractionRing (AdjoinRoot f)))] :
    WeierstrassCurve.integralModel (AdjoinRoot f)
        (W.baseChange (FractionRing (AdjoinRoot f)))
      = (WeierstrassCurve.integralModel R (W.baseChange Lv)).map
        (algebraMap R (AdjoinRoot f)) := by
  refine ABC3.Found.GenEll.integralModel_eq_of_map_eq (IsFractionRing.injective _ _) _ _ ?_
  rw [WeierstrassCurve.map_map]
  have hcomp : (algebraMap (AdjoinRoot f) (FractionRing (AdjoinRoot f))).comp
      (algebraMap R (AdjoinRoot f))
      = (quadFieldHom (K := Lv) hf hdeg).comp (algebraMap R Lv) := by
    ext x
    simp only [RingHom.coe_comp, Function.comp_apply]
    exact (quadFieldHom_algebraMap hf hdeg x).symm
  rw [hcomp, ← WeierstrassCurve.map_map]
  have h1 : (WeierstrassCurve.integralModel R (W.baseChange Lv)).map (algebraMap R Lv)
      = W.baseChange Lv := WeierstrassCurve.baseChange_integralModel_eq R (W.baseChange Lv)
  rw [h1]
  show (W.map (algebraMap L Lv)).map (quadFieldHom hf hdeg) = W.map (algebraMap L _)
  rw [WeierstrassCurve.map_map]
  congr 1
  ext x
  simp only [RingHom.coe_comp, Function.comp_apply]
  exact (halg x).symm

end IntegralModelExt

/-! ## ★★★★★★★★★★★★第 1031 —— 2 次式を monic 化する

★`HasSplitMultiplicativeReduction` の `Splits` 条件の 2 次式は

    `C c₄ · X² + C (a₁ c₄) · X − C (54 b₆ − 3 b₂ b₄ + a₂ c₄)`

であり、乗法還元では `c₄` が**単元**なので `c₄` で括り出して monic にできる。
☆その monic な因子を `f` に取れば、第 1025（剰余体で分裂）がそのまま効く。 -/

section SplitQuad

open Polynomial WeierstrassCurve in
/-- ★★★★★★★★**分裂性の 2 次式の monic 化**（第 1031）。 -/
noncomputable def splitQuadPoly {R : Type} [CommRing R] (I : WeierstrassCurve R)
    (hA : IsUnit I.c₄) : Polynomial R :=
  X ^ 2 + C I.a₁ * X - C ((54 * I.b₆ - 3 * I.b₂ * I.b₄ + I.a₂ * I.c₄)
    * ((hA.unit⁻¹ : Rˣ) : R))

open Polynomial WeierstrassCurve in
/-- ★★★★**monic である**。 -/
theorem monic_splitQuadPoly {R : Type} [CommRing R] [Nontrivial R] (I : WeierstrassCurve R)
    (hA : IsUnit I.c₄) : (splitQuadPoly I hA).Monic := by
  rw [splitQuadPoly]
  monicity!

open Polynomial WeierstrassCurve in
/-- ★★★★**次数は 2**。 -/
theorem natDegree_splitQuadPoly {R : Type} [CommRing R] [Nontrivial R] (I : WeierstrassCurve R)
    (hA : IsUnit I.c₄) : (splitQuadPoly I hA).natDegree = 2 := by
  rw [splitQuadPoly]
  compute_degree!

open Polynomial WeierstrassCurve in
/-- ★★★★★★★★**もとの 2 次式は `c₄ ×` monic である**（第 1031）。 -/
theorem quad_eq_c4_mul_splitQuadPoly {R : Type} [CommRing R] (I : WeierstrassCurve R)
    (hA : IsUnit I.c₄) :
    C I.c₄ * X ^ 2 + C (I.a₁ * I.c₄) * X
        - C (54 * I.b₆ - 3 * I.b₂ * I.b₄ + I.a₂ * I.c₄)
      = C I.c₄ * splitQuadPoly I hA := by
  have hu : I.c₄ * ((hA.unit⁻¹ : Rˣ) : R) = 1 := by
    have h1 : ((hA.unit : Rˣ) : R) * ((hA.unit⁻¹ : Rˣ) : R) = 1 := by
      rw [← Units.val_mul, mul_inv_cancel, Units.val_one]
    rwa [hA.unit_spec] at h1
  have key : I.c₄ * ((54 * I.b₆ - 3 * I.b₂ * I.b₄ + I.a₂ * I.c₄)
      * ((hA.unit⁻¹ : Rˣ) : R)) = 54 * I.b₆ - 3 * I.b₂ * I.b₄ + I.a₂ * I.c₄ := by
    rw [← mul_assoc, mul_comm I.c₄, mul_assoc, hu, mul_one]
  rw [splitQuadPoly, mul_sub, mul_add, ← C_mul, key]
  simp only [C_mul]
  ring

open Polynomial in
/-- ★★★★★★★★★★★★**単元 × monic は剰余体で分裂する**（第 1031）。

☆monic の側は第 1025、単元の側は次数 0 なので自明に分裂する。 -/
theorem splits_of_unit_mul_monic {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R]
    {f : Polynomial R} (hf : f.Monic) (hdeg : f.natDegree = 2)
    [IsLocalRing (AdjoinRoot f)]
    (A : R) (hA : IsUnit A) (q : Polynomial R) (hq : q = C A * f) :
    (q.map ((IsLocalRing.residue (AdjoinRoot f)).comp
      (algebraMap R (AdjoinRoot f)))).Splits := by
  set ψ : R →+* IsLocalRing.ResidueField (AdjoinRoot f) :=
    (IsLocalRing.residue (AdjoinRoot f)).comp (algebraMap R (AdjoinRoot f)) with hψ
  have hAu : IsUnit (ψ A) := hA.map ψ
  have hfm : (f.map ψ).Monic := hf.map ψ
  have hdeg' : (f.map ψ).natDegree = 2 := by rw [hf.natDegree_map, hdeg]
  rw [hq, Polynomial.map_mul, Polynomial.map_C]
  refine (splits_mul ?_ hfm.ne_zero).2 ⟨Splits.of_natDegree_le_one (by simp), ?_⟩
  · simpa using hAu.ne_zero
  · exact splits_map_residueField hdeg hdeg'

end SplitQuad

/-! ## ★★★★★★★★★★★★★★★★第 1032 —— 貼り合わせ

★第 1030（整モデルの同定）＋第 1031（monic 化）＋第 1025（剰余体で分裂）を繋ぐ。
☆これで **`Splits`（分裂性）が出る**——非分岐 2 次拡大の最後の一手である。 -/

section SplitsGlue

open Polynomial IsDedekindDomain NumberField WeierstrassCurve in
/-- ★★★★★★★★★★★★★★★★**拡大の整モデルの 2 次式は分裂する**（第 1032）。 -/
theorem splits_integralModel_quad_ext {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv]
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral R (W.baseChange Lv)]
    (hA : IsUnit (WeierstrassCurve.integralModel R (W.baseChange Lv)).c₄)
    {f : Polynomial R}
    (hfdef : f = splitQuadPoly (WeierstrassCurve.integralModel R (W.baseChange Lv)) hA)
    [IsDomain (AdjoinRoot f)] [IsLocalRing (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom (hfdef ▸ monic_splitQuadPoly _ hA)
        (hfdef ▸ natDegree_splitQuadPoly _ hA) (algebraMap L Lv x))
    [WeierstrassCurve.IsIntegral (AdjoinRoot f)
      (W.baseChange (FractionRing (AdjoinRoot f)))] :
    (Polynomial.map (algebraMap (AdjoinRoot f) (IsLocalRing.ResidueField (AdjoinRoot f)))
      (C (WeierstrassCurve.integralModel (AdjoinRoot f)
            (W.baseChange (FractionRing (AdjoinRoot f)))).c₄ * X ^ 2
        + C ((WeierstrassCurve.integralModel (AdjoinRoot f)
              (W.baseChange (FractionRing (AdjoinRoot f)))).a₁
            * (WeierstrassCurve.integralModel (AdjoinRoot f)
              (W.baseChange (FractionRing (AdjoinRoot f)))).c₄) * X
        - C (54 * (WeierstrassCurve.integralModel (AdjoinRoot f)
                (W.baseChange (FractionRing (AdjoinRoot f)))).b₆
            - 3 * (WeierstrassCurve.integralModel (AdjoinRoot f)
                  (W.baseChange (FractionRing (AdjoinRoot f)))).b₂
                * (WeierstrassCurve.integralModel (AdjoinRoot f)
                  (W.baseChange (FractionRing (AdjoinRoot f)))).b₄
            + (WeierstrassCurve.integralModel (AdjoinRoot f)
                  (W.baseChange (FractionRing (AdjoinRoot f)))).a₂
                * (WeierstrassCurve.integralModel (AdjoinRoot f)
                  (W.baseChange (FractionRing (AdjoinRoot f)))).c₄))).Splits := by
  subst hfdef
  set I := WeierstrassCurve.integralModel R (W.baseChange Lv) with hI
  have hI' : WeierstrassCurve.integralModel (AdjoinRoot (splitQuadPoly I hA))
      (W.baseChange (FractionRing (AdjoinRoot (splitQuadPoly I hA))))
      = I.map (algebraMap R (AdjoinRoot (splitQuadPoly I hA))) :=
    integralModel_ext (monic_splitQuadPoly I hA) (natDegree_splitQuadPoly I hA) halg W
  rw [hI']
  simp only [WeierstrassCurve.map_c₄, WeierstrassCurve.map_a₁, WeierstrassCurve.map_b₂,
    WeierstrassCurve.map_b₄, WeierstrassCurve.map_b₆, WeierstrassCurve.map_a₂]
  have hpoly : (C ((algebraMap R (AdjoinRoot (splitQuadPoly I hA))) I.c₄) * X ^ 2
      + C ((algebraMap R (AdjoinRoot (splitQuadPoly I hA))) I.a₁
          * (algebraMap R (AdjoinRoot (splitQuadPoly I hA))) I.c₄) * X
      - C (54 * (algebraMap R (AdjoinRoot (splitQuadPoly I hA))) I.b₆
          - 3 * (algebraMap R (AdjoinRoot (splitQuadPoly I hA))) I.b₂
              * (algebraMap R (AdjoinRoot (splitQuadPoly I hA))) I.b₄
          + (algebraMap R (AdjoinRoot (splitQuadPoly I hA))) I.a₂
              * (algebraMap R (AdjoinRoot (splitQuadPoly I hA))) I.c₄))
      = (C I.c₄ * X ^ 2 + C (I.a₁ * I.c₄) * X
        - C (54 * I.b₆ - 3 * I.b₂ * I.b₄ + I.a₂ * I.c₄)).map
        (algebraMap R (AdjoinRoot (splitQuadPoly I hA))) := by
    simp [Polynomial.map_add, Polynomial.map_sub, Polynomial.map_mul, Polynomial.map_C,
      Polynomial.map_X, Polynomial.map_pow, map_ofNat]
  rw [hpoly, Polynomial.map_map]
  exact splits_of_unit_mul_monic (monic_splitQuadPoly I hA) (natDegree_splitQuadPoly I hA)
    I.c₄ hA _ (quad_eq_c4_mul_splitQuadPoly I hA)

open Polynomial IsDedekindDomain NumberField WeierstrassCurve in
/-- ★★★★★★★★★★★★★★★★★★★★**拡大の上では分裂乗法還元をもつ**（第 1033）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1033）**——不分岐 2 次拡大の**到達点**である。
☆`f` を「整モデルの 2 次式を monic 化したもの」に取れば、
`R[X]/(f)` の剰余体には `f` の根があるので 2 次式が分裂する。
★これで非分裂の枝が閉じる。 -/
theorem hasSplitMultiplicativeReduction_ext {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv]
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral R (W.baseChange Lv)]
    (hA : IsUnit (WeierstrassCurve.integralModel R (W.baseChange Lv)).c₄)
    {f : Polynomial R}
    (hfdef : f = splitQuadPoly (WeierstrassCurve.integralModel R (W.baseChange Lv)) hA)
    [IsDomain (AdjoinRoot f)] [IsLocalRing (AdjoinRoot f)]
    [IsDiscreteValuationRing (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom (hfdef ▸ monic_splitQuadPoly _ hA)
        (hfdef ▸ natDegree_splitQuadPoly _ hA) (algebraMap L Lv x))
    [WeierstrassCurve.HasMultiplicativeReduction (AdjoinRoot f)
      (W.baseChange (FractionRing (AdjoinRoot f)))] :
    WeierstrassCurve.HasSplitMultiplicativeReduction (AdjoinRoot f)
      (W.baseChange (FractionRing (AdjoinRoot f))) :=
  ⟨splits_integralModel_quad_ext W hA hfdef halg⟩

end SplitsGlue

/-! ## ★★★★★★★★★★★★第 1035 —— 分裂しない monic 2 次式は既約

★場合分けの「非分裂」側で要るのは
「2 次式が分裂しない ⟹ その monic 因子は既約」である。
☆体上の monic 2 次式が可約なら 1 次×1 次で、どちらも分裂するから全体も分裂する。 -/

open Polynomial in
/-- ★★★★★★★★**分裂しない monic 2 次式は既約**（第 1035）。 -/
theorem irreducible_of_not_splits_natDegree_two {k : Type} [Field k] {q : Polynomial k}
    (hm : q.Monic) (hdeg : q.natDegree = 2) (hns : ¬ q.Splits) : Irreducible q := by
  have hq0 : q ≠ 0 := hm.ne_zero
  refine ⟨Polynomial.not_isUnit_of_natDegree_pos q (by omega), ?_⟩
  intro a b hab
  by_contra hcon
  rw [not_or] at hcon
  obtain ⟨ha, hb⟩ := hcon
  have ha0 : a ≠ 0 := by rintro rfl; simp at hab; exact hq0 hab
  have hb0 : b ≠ 0 := by rintro rfl; simp at hab; exact hq0 hab
  have hda : 0 < a.natDegree := by
    rcases Nat.eq_zero_or_pos a.natDegree with h | h
    · exact absurd (Polynomial.isUnit_iff_degree_eq_zero.2
        (by rw [Polynomial.degree_eq_natDegree ha0, h]; rfl)) ha
    · exact h
  have hdb : 0 < b.natDegree := by
    rcases Nat.eq_zero_or_pos b.natDegree with h | h
    · exact absurd (Polynomial.isUnit_iff_degree_eq_zero.2
        (by rw [Polynomial.degree_eq_natDegree hb0, h]; rfl)) hb
    · exact h
  have hsum : a.natDegree + b.natDegree = 2 := by
    rw [← hdeg, hab, Polynomial.natDegree_mul ha0 hb0]
  exact hns (by
    rw [hab, splits_mul ha0 hb0]
    exact ⟨Splits.of_natDegree_le_one (by omega), Splits.of_natDegree_le_one (by omega)⟩)

open Polynomial in
/-- ★★★★**単元倍しても分裂性は変わらない**（片側）。 -/
theorem splits_unit_mul {k : Type} [Field k] {u : k} (hu : u ≠ 0) {g : Polynomial k}
    (hg : g ≠ 0) (h : g.Splits) : (C u * g).Splits := by
  refine (splits_mul (by simpa using hu) hg).2 ⟨Splits.of_natDegree_le_one (by simp), h⟩

open Polynomial WeierstrassCurve in
/-- ★★★★★★★★★★★★**monic 化した 2 次式が分裂すれば分裂乗法還元**（第 1036）。

☆場合分けの対偶を取るための形——`¬ HasSplit` から `¬ Splits (f̄)` が直ちに出る。 -/
theorem hasSplit_of_splits_splitQuadPoly {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K)
    [WeierstrassCurve.HasMultiplicativeReduction R W]
    (hA : IsUnit (WeierstrassCurve.integralModel R W).c₄)
    (hsp : ((splitQuadPoly (WeierstrassCurve.integralModel R W) hA).map
      (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))).Splits) :
    WeierstrassCurve.HasSplitMultiplicativeReduction R W := by
  refine ⟨?_⟩
  rw [quad_eq_c4_mul_splitQuadPoly _ hA, Polynomial.map_mul, Polynomial.map_C]
  exact splits_unit_mul (hA.map (algebraMap R (IsLocalRing.ResidueField R))).ne_zero
    (Polynomial.Monic.ne_zero (Polynomial.Monic.map _ (monic_splitQuadPoly _ hA))) hsp

open Polynomial in
/-- ★★★★★★★★**剰余体で分裂しない monic 2 次式は既約**（第 1036）。

☆配管の注意: `Ideal.Quotient.mk (maximalIdeal R)` のままだと
`Field (R ⧸ 𝔪)` のインスタンス探索が爆発する。
★`show` で `IsLocalRing.residue R` の形に言い直してから当てる（0.04 秒になる）。 -/
theorem irreducible_map_residue_of_not_splits {R : Type} [CommRing R] [IsLocalRing R]
    {f : Polynomial R} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hns : ¬ (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))).Splits) :
    Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))) := by
  show Irreducible (f.map (IsLocalRing.residue R))
  refine irreducible_of_not_splits_natDegree_two (hf.map (IsLocalRing.residue R)) ?_ hns
  rw [hf.natDegree_map, hdeg]

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
