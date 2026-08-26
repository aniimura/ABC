/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.HeightOneDVR
import ABC3.Meta.Claim

/-!
# 代数的 Hartogs —— 正規 Noether 整域は高さ 1 の局所化の交わり

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> that [since V [L] is a proper normal variety] for A ∈Ob(CV,

## ★★これは `Example 6.1` の壁だった

原文が "(we observe that)" で畳んだ 1 行 `𝒪^×(A) = 𝒪^▷(A) = k_L^×` を出すのに
2 本要り、そのうち 1 本(`proper` なら `Γ(X,𝒪)` が `k` 上有限な体)は
**mathlib にあった**(`AlgebraicGeometry.isField_of_universallyClosed` /
`finite_appTop_of_universallyClosed`)。★残る 1 本が本ファイルである。

## ★★★証明の骨(4 段)

`x ∈ K` が**すべての高さ 1 の素で局所化に入る**なら `x ∈ R` であること。

| 段 | 中身 |
|---|---|
| 1 | **分母イデアル** `den R y = (R : y)` を `Submodule.colon` で置く |
| 2 | `x ∉ R` なら、`den R x` を含む**極大な分母イデアル** `p = den R y` が取れる(Noether) |
| 3 | ★そのような `p` は**素**である(`(R : ay) ⊇ (R : y)` と極大性) |
| 4 | ★★**二分**: `p·y ⊆ p` なら `y` は整だから `y ∈ R`(正規性)——矛盾。<br>そうでなければ `c ∈ p`, `d ∉ p` があって `d = c·y`、このとき **`p·R_p = c·R_p` は主**。<br>主なら `IsDiscreteValuationRing.TFAE` で `R_p` の非零素は 1 つだけ、すなわち `p` は高さ 1 |

★最後に仮定を `p` に当てると `den R x ⊄ p` だが `den R x ⊆ p` なので矛盾する。

## ★在庫の測定(2026-08-25)

| 要るもの | mathlib | 判定 |
|---|---|---|
| `MaximalSpectrum.iInf_localization_eq_bot` | ある | ★**極大**イデアル版のみ(高さ 1 版ではない) |
| `IsKrullDomain` | 無い | `Mathlib/RingTheory/` に `KrullDimension` しかない |
| `isIntegral_of_smul_mem_submodule` | ある | ★段 4 の「整である」の中身 |
| `IsIntegrallyClosed.isIntegral_iff` | ある | ★正規性 |
| `IsDiscreteValuationRing.TFAE` | ある | ★段 4 の「主 ⟹ 非零素は 1 つ」 |
-/

namespace ABC3.Found.Divisor.Hartogs

open Submodule IsLocalRing

variable {R : Type*} [CommRing R] [IsDomain R] {K : Type*} [Field K] [Algebra R K]
  [IsFractionRing R K]

/-! ## ★1. 分母イデアル -/

/-- ★`K` の中の `R`(像としての部分加群)。 -/
abbrev Rsub (R K : Type*) [CommRing R] [Field K] [Algebra R K] : Submodule R K :=
  LinearMap.range (Algebra.linearMap R K)

/-- ★★`y` の**分母イデアル** `(R : y) = {r ∈ R | r·y ∈ R}`。 -/
def den (R : Type*) [CommRing R] {K : Type*} [Field K] [Algebra R K] (y : K) : Ideal R :=
  (Rsub R K).colon {y}

omit [IsDomain R] [IsFractionRing R K] in
theorem mem_den {y : K} {r : R} : r ∈ den R y ↔ r • y ∈ Rsub R K := by
  rw [den, Submodule.mem_colon]
  exact ⟨fun h => h y rfl, fun h s hs => by rw [Set.mem_singleton_iff] at hs; subst hs; exact h⟩

omit [IsDomain R] [IsFractionRing R K] in
theorem den_ne_top_iff {y : K} : den R y ≠ ⊤ ↔ y ∉ Rsub R K := by
  constructor
  · intro h hy
    exact h (Ideal.eq_top_of_isUnit_mem _ (mem_den.mpr (by simpa using hy)) isUnit_one)
  · intro hy h
    exact hy (by simpa using (mem_den (R := R) (y := y) (r := 1)).mp (h ▸ Submodule.mem_top))

/-- ★分母イデアルは 0 でない(分数体の元は分母を持つ)。 -/
theorem den_ne_bot (y : K) : den R y ≠ ⊥ := by
  obtain ⟨a, b, hb, hab⟩ := IsFractionRing.div_surjective (A := R) (K := K) y
  have hb0 : (algebraMap R K) b ≠ 0 :=
    IsFractionRing.to_map_ne_zero_of_mem_nonZeroDivisors hb
  have hmem : b ∈ den R y := by
    refine mem_den.mpr ⟨a, ?_⟩
    show (algebraMap R K) a = b • y
    rw [← hab, Algebra.smul_def]
    field_simp
  intro h
  rw [h, Ideal.mem_bot] at hmem
  exact (nonZeroDivisors.ne_zero hb) hmem

omit [IsDomain R] [IsFractionRing R K] in
theorem den_le_den_smul (a : R) (y : K) : den R y ≤ den R (a • y) := by
  intro r hr
  refine mem_den.mpr ?_
  have h : r • (a • y) = a • (r • y) := by rw [smul_comm]
  rw [h]
  exact Submodule.smul_mem _ a (mem_den.mp hr)

/-! ## ★2. 極大な分母イデアル -/

variable [IsNoetherianRing R]

omit [IsDomain R] [IsFractionRing R K] in
/-- ★★分母イデアルの中で極大なものが取れる(Noether)。 -/
theorem exists_maximal_den (x : K) (hx : x ∉ Rsub R K) :
    ∃ y : K, y ∉ Rsub R K ∧ den R x ≤ den R y ∧
      ∀ z : K, z ∉ Rsub R K → den R x ≤ den R z → den R y ≤ den R z → den R z = den R y := by
  set S : Set (Ideal R) := {J | ∃ y : K, y ∉ Rsub R K ∧ J = den R y ∧ den R x ≤ J} with hS
  have hne : S.Nonempty := ⟨den R x, x, hx, rfl, le_refl _⟩
  obtain ⟨p, hpS, hmax⟩ :=
    (IsWellFounded.wf (r := (· > · : Ideal R → Ideal R → Prop))).has_min S hne
  obtain ⟨y, hy, hpy, hxp⟩ := hpS
  refine ⟨y, hy, hpy ▸ hxp, ?_⟩
  intro z hz hxz hyz
  by_contra hne'
  exact hmax (den R z) ⟨z, hz, rfl, hxz⟩ (lt_of_le_of_ne (hpy ▸ hyz) fun h => hne' (hpy ▸ h.symm))

omit [IsDomain R] [IsNoetherianRing R] [IsFractionRing R K] in
/-- ★★★**極大な分母イデアルは素イデアル**。

★`a ∉ p` なら `ay ∉ R` で `(R : ay) ⊇ p`、極大性から `(R : ay) = p`。 -/
theorem isPrime_of_maximal_den {x y : K} (hy : y ∉ Rsub R K) (hxy : den R x ≤ den R y)
    (hmax : ∀ z : K, z ∉ Rsub R K → den R x ≤ den R z → den R y ≤ den R z → den R z = den R y) :
    (den R y).IsPrime := by
  refine ⟨den_ne_top_iff.mpr hy, ?_⟩
  intro a b hab
  by_contra hcon
  rw [not_or] at hcon
  obtain ⟨ha, hb⟩ := hcon
  have hay : a • y ∉ Rsub R K := fun h => ha (mem_den.mpr h)
  have h1 : den R y ≤ den R (a • y) := den_le_den_smul a y
  have h2 : den R (a • y) = den R y := hmax _ hay (hxy.trans h1) h1
  refine hb ?_
  rw [← h2]
  refine mem_den.mpr ?_
  have h3 : b • (a • y) = (a * b) • y := by rw [smul_smul, mul_comm]
  rw [h3]
  exact mem_den.mp hab

/-! ## ★3. 二分の第 1 の場合 —— 整だから `R` に入る -/

variable [IsIntegrallyClosed R]

omit [IsDomain R] in
/-- ★★★★**分母イデアルを自分に写すなら整** —— 正規性で `y ∈ R`。 -/
theorem mem_Rsub_of_smul_mem {y : K} {p : Ideal R} (hp0 : p ≠ ⊥)
    (h : ∀ n ∈ Submodule.map (Algebra.linearMap R K) p,
      y • n ∈ Submodule.map (Algebra.linearMap R K) p) :
    y ∈ Rsub R K := by
  have hinj : Function.Injective (algebraMap R K) := IsFractionRing.injective R K
  have hN0 : Submodule.map (Algebra.linearMap R K) p ≠ ⊥ := by
    intro hc
    apply hp0
    rw [eq_bot_iff]
    intro a ha
    have hmem : (Algebra.linearMap R K) a ∈ Submodule.map (Algebra.linearMap R K) p := ⟨a, ha, rfl⟩
    rw [hc, Submodule.mem_bot] at hmem
    exact hinj (by simpa using hmem)
  have hNfg : (Submodule.map (Algebra.linearMap R K) p).FG := (IsNoetherian.noetherian p).map _
  obtain ⟨a, ha⟩ :=
    IsIntegrallyClosed.isIntegral_iff.mp (isIntegral_of_smul_mem_submodule _ hN0 hNfg y h)
  exact ⟨a, ha⟩

/-! ## ★4. 二分の第 2 の場合 —— 局所化の極大イデアルが主になる -/

omit [IsNoetherianRing R] [IsIntegrallyClosed R] in
/-- ★★★★★**`c ∈ p`, `d ∉ p`, `d = c·y` なら `p·R_p = c·R_p`**。

★`r ∈ p` に対して `r·y = e ∈ R` を取ると `r·d = c·e` であり、
`d` は `R_p` で単元なので `r = c·e·d⁻¹ ∈ c·R_p`。 -/
theorem maximalIdeal_eq_span {p : Ideal R} [hp : p.IsPrime] {y : K} (hpy : p = den R y)
    {c : R} (hc : c ∈ p) {d : R} (hd : d ∉ p)
    (hcd : (algebraMap R K) d = (algebraMap R K) c * y) :
    IsLocalRing.maximalIdeal (Localization.AtPrime p)
      = Ideal.span {(algebraMap R (Localization.AtPrime p)) c} := by
  set S := Localization.AtPrime p
  have hinj : Function.Injective (algebraMap R K) := IsFractionRing.injective R K
  have hdu : IsUnit ((algebraMap R S) d) :=
    (IsLocalization.AtPrime.isUnit_to_map_iff S p d).mpr hd
  refine le_antisymm ?_ ?_
  · rw [← IsLocalization.AtPrime.map_eq_maximalIdeal p S, Ideal.map_le_iff_le_comap]
    intro r hr
    obtain ⟨e, he⟩ := mem_den.mp (hpy ▸ hr)
    have heK : (algebraMap R K) e = (algebraMap R K) r * y := by
      have h : (Algebra.linearMap R K) e = r • y := he
      simpa [Algebra.smul_def] using h
    have hRR : r * d = c * e := by
      refine hinj ?_
      rw [map_mul, map_mul, hcd, heK]; ring
    have hSS : (algebraMap R S) r * (algebraMap R S) d
        = (algebraMap R S) c * (algebraMap R S) e := by
      rw [← map_mul, ← map_mul, hRR]
    obtain ⟨u, hu⟩ := hdu
    refine Ideal.mem_span_singleton'.mpr ⟨(algebraMap R S) e * (↑u⁻¹ : S), ?_⟩
    have hru : (algebraMap R S) r = (algebraMap R S) r * (algebraMap R S) d * (↑u⁻¹ : S) := by
      rw [← hu, mul_assoc, Units.mul_inv, mul_one]
    rw [hru, hSS]; ring
  · rw [Ideal.span_le, Set.singleton_subset_iff,
      ← IsLocalization.AtPrime.map_eq_maximalIdeal p S]
    exact Ideal.mem_map_of_mem _ hc

omit [IsIntegrallyClosed R] in
/-- ★★★★★**局所化の極大イデアルが主なら、`p` は極小非零素**(高さ 1)。

★`IsDiscreteValuationRing.TFAE` の「主 ⟺ 整閉かつ非零素が唯一」を使い、
素イデアルの対応(`orderIsoOfPrime`)で `R` の側へ降ろす。 -/
theorem minimal_of_maximalIdeal_isPrincipal {p : Ideal R} [hp : p.IsPrime] (hp0 : p ≠ ⊥)
    (hprin : Submodule.IsPrincipal (IsLocalRing.maximalIdeal (Localization.AtPrime p))) :
    ∀ q : Ideal R, q.IsPrime → q ≠ ⊥ → q ≤ p → q = p := by
  set S := Localization.AtPrime p with hSdef
  haveI : IsNoetherianRing S := IsLocalization.isNoetherianRing p.primeCompl S inferInstance
  haveI : IsLocalRing S := IsLocalization.AtPrime.isLocalRing S p
  have hnf : ¬ IsField S := IsLocalization.AtPrime.not_isField R hp0 S
  obtain ⟨-, huniq⟩ := ((IsDiscreteValuationRing.TFAE S hnf).out 4 3).mp hprin
  intro q hq hq0 hqp
  haveI := hq
  set Q := (IsLocalization.AtPrime.orderIsoOfPrime S p).symm ⟨q, ⟨hq, hqp⟩⟩ with hQ
  have hfwd : (Q.1).comap (algebraMap R S) = q :=
    congrArg Subtype.val
      ((IsLocalization.AtPrime.orderIsoOfPrime S p).apply_symm_apply ⟨q, ⟨hq, hqp⟩⟩)
  have hQ0 : Q.1 ≠ ⊥ := by
    intro hc
    apply hq0
    rw [← hfwd, hc, Ideal.comap_bot_of_injective _
      (IsLocalization.injective S p.primeCompl_le_nonZeroDivisors)]
  have hmax0 : (IsLocalRing.maximalIdeal S) ≠ ⊥ :=
    IsLocalRing.isField_iff_maximalIdeal_eq.not.mp hnf
  have hQm := huniq.unique ⟨hQ0, Q.2⟩ ⟨hmax0, inferInstance⟩
  have hund : Ideal.comap (algebraMap R S) (IsLocalRing.maximalIdeal S) = p :=
    IsLocalization.AtPrime.under_maximalIdeal S p
  rw [← hfwd, hQm, hund]

/-! ## ★5. 本体 -/

/-- ★★★★★★**代数的 Hartogs** ——
`x ∈ K` が**すべての高さ 1 の素で局所化に入る**なら `x ∈ R`。

★「`x ∈ R_p`」は「`den R x` が `p` に含まれない」で表す。 -/
theorem mem_Rsub_of_forall_heightOne (x : K)
    (h : ∀ p : Ideal R, p.IsPrime → p ≠ ⊥ →
      (∀ q : Ideal R, q.IsPrime → q ≠ ⊥ → q ≤ p → q = p) →
      ∃ b ∈ den R x, b ∉ p) :
    x ∈ Rsub R K := by
  by_contra hx
  obtain ⟨y, hy, hxy, hmax⟩ := exists_maximal_den x hx
  haveI hp : (den R y).IsPrime := isPrime_of_maximal_den hy hxy hmax
  have hp0 : den R y ≠ ⊥ := den_ne_bot y
  by_cases hcase : ∀ n ∈ Submodule.map (Algebra.linearMap R K) (den R y),
      y • n ∈ Submodule.map (Algebra.linearMap R K) (den R y)
  · exact hy (mem_Rsub_of_smul_mem hp0 hcase)
  · rw [not_forall] at hcase
    obtain ⟨n, hcase⟩ := hcase
    rw [Classical.not_imp] at hcase
    obtain ⟨hn, hn'⟩ := hcase
    obtain ⟨c, hc, hcn⟩ := hn
    obtain ⟨d, hd⟩ := mem_den.mp hc
    have hdK : (algebraMap R K) d = (algebraMap R K) c * y := by
      have h1 : (Algebra.linearMap R K) d = c • y := hd
      simpa [Algebra.smul_def] using h1
    have hdp : d ∉ den R y := by
      intro hdmem
      refine hn' ?_
      have h2 : y • n = (algebraMap R K) d := by
        rw [← hcn]
        show y * (Algebra.linearMap R K) c = _
        rw [hdK]
        show y * (algebraMap R K) c = (algebraMap R K) c * y
        ring
      rw [h2]
      exact ⟨d, hdmem, rfl⟩
    have hprin : Submodule.IsPrincipal
        (IsLocalRing.maximalIdeal (Localization.AtPrime (den R y))) :=
      ⟨⟨_, maximalIdeal_eq_span rfl hc hdp hdK⟩⟩
    obtain ⟨b, hb, hbp⟩ :=
      h (den R y) hp hp0 (minimal_of_maximalIdeal_isPrincipal hp0 hprin)
    exact hbp (hxy hb)

/-- ★★**`HeightOnePrime` の言葉での言い換え**。 -/
theorem mem_Rsub_of_forall_heightOnePrime (x : K)
    (h : ∀ v : ABC3.Found.Divisor.HeightOnePrime R, ∃ b ∈ den R x, b ∉ v.asIdeal) :
    x ∈ Rsub R K :=
  mem_Rsub_of_forall_heightOne x fun p hp hp0 hmin => h ⟨p, hp, hp0, hmin⟩

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Example 6.1` の代数的 Hartogs。 -/
def mem_Rsub_of_forall_heightOne.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — 代数的 Hartogs(正規 Noether 整域は高さ 1 の局所化の交わり)",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor.Hartogs
