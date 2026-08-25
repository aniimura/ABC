/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithFunctor

/-!
# 素点の基本等式（相対版）—— `Σ_{w | v} localDeg(w) = [M : L]`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

## ★★何のために要るか

`Theorem 6.4, (i)` の 5 圏のうち **`𝒞^rlf`** を算術で出すには
`Φ^rlf : MonoidOn 𝒟` の条件 (a)（節点 `rlf-flat`）が要る。

★★在庫の `isCharacteristicallyInjective_scMap_of_nsmul_retraction`
（`Def24ScTransport.lean`）により、これは

    h ∘ (Φ.map α) = k · id   なる h の構成

に帰着している。算術では **「上にある素点についての和」**

    h(y)(v) := Σ_{V | v} y(V)

を取れば `h(f(x))(v) = (Σ_{V | v} localDeg V) · x(v)` なので、
★**要るのは相対版の基本等式 `Σ_{V | v} localDeg(V) = [M : L]` ただ 1 本**である。

## ★★★中身（無限素点の側）

★mathlib の相対版（`InfinitePlace/Ramification.lean` の `card_isUnramified` ほか）は
**Galois 専用**で、一般の `L ⊆ M` 用の等式は無い（2026-08-25 に確認）。
そこで埋め込みの二重数え上げで作る:

    S := {ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v}

| 数え方 | 結果 | 使う在庫 |
|---|---|---|
| `mk` のファイバー | `Σ_{w \| v} mult w` | `card_filter_mk_eq`（mathlib） |
| `ψ ↦ ψ ∘ algebraMap` のファイバー | `mult v · [M:L]` | `AlgHom.card`（mathlib） |

★★見積りは **150〜300 行**としていたが、mathlib に
`card_filter_mk_eq` / `comap_mk` / `AlgHom.card` が揃っていたので
**約 60 行**で済んだ。
-/

namespace ABC3.Found.Divisor

open NumberField InfinitePlace Finset

/-! ## ★1. 与えられた埋め込みの延長は `[M:L]` 個 -/

open scoped Classical in
/-- ★★`φ : L →+* ℂ` の `M` への延長はちょうど `[M : L]` 個。

★`ℂ` を `φ` で `L`-代数と見れば `M →ₐ[L] ℂ` そのものなので、
mathlib の `AlgHom.card` がそのまま効く。 -/
theorem card_filter_comp_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (φ : L →+* ℂ) :
    #{ψ : M →+* ℂ | ψ.comp (algebraMap L M) = φ} = Module.finrank L M := by
  classical
  letI : Algebra L ℂ := φ.toAlgebra
  have e : {ψ : M →+* ℂ // ψ.comp (algebraMap L M) = φ} ≃ (M →ₐ[L] ℂ) := by
    refine ⟨fun p => ⟨p.1, ?_⟩, fun f => ⟨f.toRingHom, ?_⟩, ?_, ?_⟩
    · intro x
      show p.1 (algebraMap L M x) = algebraMap L ℂ x
      rw [← RingHom.comp_apply, p.2]
      rfl
    · ext x
      show f (algebraMap L M x) = φ x
      exact f.commutes x
    · intro p; rfl
    · intro f; rfl
  rw [← Fintype.card_subtype, Fintype.card_congr e]
  exact AlgHom.card L M ℂ

/-! ## ★2. 二重数え上げ -/

open scoped Classical in
/-- ★★★★`v` の上にある埋め込みの個数は `mult v · [M:L]`。

★`ψ ↦ ψ ∘ algebraMap L M` のファイバーで数える。 -/
theorem card_filter_comap_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (v : InfinitePlace L) :
    #{ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v}
      = v.mult * Module.finrank L M := by
  classical
  rw [Finset.card_eq_sum_card_fiberwise
    (f := fun ψ : M →+* ℂ => ψ.comp (algebraMap L M))
    (t := {φ : L →+* ℂ | mk φ = v})
    (by
      intro ψ hψ
      rw [Finset.mem_coe, mem_filter_univ] at hψ
      rw [Finset.mem_coe, mem_filter_univ, ← comap_mk]
      exact hψ)]
  have hinner : ∀ φ ∈ ({φ : L →+* ℂ | mk φ = v} : Finset (L →+* ℂ)),
      #({a ∈ {ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v} |
          a.comp (algebraMap L M) = φ}) = (fun _ : L →+* ℂ => Module.finrank L M) φ := by
    intro φ hφ
    rw [mem_filter_univ] at hφ
    have hset : ({a ∈ {ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v} |
        a.comp (algebraMap L M) = φ} : Finset (M →+* ℂ))
        = ({ψ : M →+* ℂ | ψ.comp (algebraMap L M) = φ} : Finset (M →+* ℂ)) := by
      ext ψ
      simp only [mem_filter, mem_univ, true_and]
      constructor
      · rintro ⟨-, h2⟩; exact h2
      · intro h2
        refine ⟨?_, h2⟩
        rw [comap_mk, h2]; exact hφ
    show #({a ∈ {ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v} |
        a.comp (algebraMap L M) = φ} : Finset (M →+* ℂ)) = Module.finrank L M
    rw [hset]
    exact card_filter_comp_eq L M φ
  rw [Finset.sum_congr rfl hinner, Finset.sum_const, smul_eq_mul, card_filter_mk_eq]

open scoped Classical in
/-- ★★★★★★**無限素点の相対基本等式** —— `Σ_{w | v} mult w = mult v · [M:L]`。

★同じ集合を `mk` のファイバーで数え直すだけである。 -/
theorem sum_mult_comap_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (v : InfinitePlace L) :
    ∑ w ∈ {w : InfinitePlace M | w.comap (algebraMap L M) = v}, w.mult
      = v.mult * Module.finrank L M := by
  classical
  rw [← card_filter_comap_eq L M v,
    Finset.card_eq_sum_card_fiberwise
      (f := fun ψ : M →+* ℂ => mk ψ)
      (t := {w : InfinitePlace M | w.comap (algebraMap L M) = v})
      (by
        intro ψ hψ
        rw [Finset.mem_coe, mem_filter_univ] at hψ
        rw [Finset.mem_coe, mem_filter_univ]
        exact hψ)]
  refine (Finset.sum_congr rfl ?_).symm
  intro w hw
  rw [mem_filter_univ] at hw
  have hset : ({a ∈ {ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v} | mk a = w}
      : Finset (M →+* ℂ)) = ({ψ : M →+* ℂ | mk ψ = w} : Finset (M →+* ℂ)) := by
    ext ψ
    simp only [mem_filter, mem_univ, true_and]
    constructor
    · rintro ⟨-, h2⟩; exact h2
    · intro h2
      exact ⟨by rw [h2]; exact hw, h2⟩
  rw [hset, card_filter_mk_eq]

/-! ## ★3. `localDeg` の形 -/

open scoped Classical in
/-- ★★★★★★**無限素点の相対基本等式**（`localDeg` の形）——
`Σ_{W | v} localDeg(W) = [M : L]`。

★`localDeg W = mult W / mult v` で `mult v ∣ mult W`（`mult_resInf_dvd`）なので、
上の等式を `mult v` で割ればよい。 -/
theorem sum_localDeg_inf_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (v : InfinitePlace L) :
    ∑ W ∈ {W : InfinitePlace M | resInf (L := L) W = v},
        localDeg (L := L) (Sum.inr W) = Module.finrank L M := by
  classical
  have hkey := sum_mult_comap_eq L M v
  have hmul : (∑ W ∈ {W : InfinitePlace M | resInf (L := L) W = v},
      localDeg (L := L) (Sum.inr W)) * v.mult
      = ∑ W ∈ {W : InfinitePlace M | resInf (L := L) W = v}, W.mult := by
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl (fun W hW => ?_)
    rw [mem_filter_univ] at hW
    have h1 : localDeg (L := L) (Sum.inr W) * (resInf (L := L) W).mult = W.mult :=
      localDeg_mul_mult W
    rw [← hW]
    exact h1
  have hset : ({W : InfinitePlace M | resInf (L := L) W = v} : Finset (InfinitePlace M))
      = ({W : InfinitePlace M | W.comap (algebraMap L M) = v} : Finset (InfinitePlace M)) := rfl
  rw [hset] at hmul
  rw [hkey] at hmul
  have hv : 0 < v.mult := by
    rcases mult_eq_one_or_two v with h | h <;> omega
  rw [hset]
  exact Nat.eq_of_mul_eq_mul_right hv (hmul.trans (by ring))

/-! ## ★4. 有限素点の側

★こちらは mathlib の `Ideal.sum_ramification_inertia`（`Σ e·f = [M:L]`）が
そのまま使える。★`FinitePlace M` は**有限型ではない**（素点は無限個ある）ので、
ファイバーを `finite_fiber_fin` の `Set.Finite.toFinset` で取る。 -/

open IsDedekindDomain in
open scoped Classical in
/-- ★★★★★**有限素点の相対基本等式** —— `Σ_{W | v} localDeg(W) = [M : L]`。

★`W ↦ W.maximalIdeal.asIdeal` が `primesOverFinset` への全単射を与える。 -/
theorem sum_localDeg_fin_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (v : FinitePlace L) :
    ∑ W ∈ (finite_fiber_fin (L := L) (M := M) v).toFinset,
        localDeg (L := L) (Sum.inl W) = Module.finrank L M := by
  classical
  have hne : (v.maximalIdeal.asIdeal) ≠ ⊥ := v.maximalIdeal.ne_bot
  have hmem : ∀ W : FinitePlace M, W ∈ (finite_fiber_fin (L := L) (M := M) v).toFinset
      ↔ resHOS (L := L) W.maximalIdeal = v.maximalIdeal := by
    intro W
    rw [Set.Finite.mem_toFinset]
    show resFin (L := L) W = v ↔ _
    constructor
    · intro h; rw [← h, resFin, FinitePlace.maximalIdeal_mk]
    · intro h; rw [resFin, h, FinitePlace.mk_maximalIdeal]
  rw [← Ideal.sum_ramification_inertia (p := v.maximalIdeal.asIdeal) (𝓞 M) L M hne]
  refine Finset.sum_bij'
    (i := fun W _ => W.maximalIdeal.asIdeal)
    (j := fun P hP =>
      haveI := ((IsDedekindDomain.mem_primesOverFinset_iff hne (𝓞 M)).mp hP).1
      haveI := ((IsDedekindDomain.mem_primesOverFinset_iff hne (𝓞 M)).mp hP).2
      FinitePlace.mk ⟨P, inferInstance, Ideal.ne_bot_of_liesOver_of_ne_bot hne P⟩)
    ?_ ?_ ?_ ?_ ?_
  · intro W hW
    rw [hmem] at hW
    refine (IsDedekindDomain.mem_primesOverFinset_iff hne (𝓞 M)).mpr ⟨inferInstance, ⟨?_⟩⟩
    exact (congrArg HeightOneSpectrum.asIdeal hW).symm
  · intro P hP
    haveI := ((IsDedekindDomain.mem_primesOverFinset_iff hne (𝓞 M)).mp hP).1
    haveI := ((IsDedekindDomain.mem_primesOverFinset_iff hne (𝓞 M)).mp hP).2
    rw [hmem, FinitePlace.maximalIdeal_mk]
    refine HeightOneSpectrum.ext ?_
    show P.under (𝓞 L) = v.maximalIdeal.asIdeal
    exact (Ideal.LiesOver.over (P := P) (p := v.maximalIdeal.asIdeal)).symm
  · intro W hW
    show FinitePlace.mk _ = W
    rw [show (⟨W.maximalIdeal.asIdeal, _, _⟩ : HeightOneSpectrum (𝓞 M)) = W.maximalIdeal from rfl,
      FinitePlace.mk_maximalIdeal]
  · intro P hP
    show (FinitePlace.mk _).maximalIdeal.asIdeal = P
    rw [FinitePlace.maximalIdeal_mk]
  · intro W hW
    rw [hmem] at hW
    show ramIdx (L := L) W * inertDeg (L := L) W = _
    rw [ramIdx, inertDeg, hW]

/-! ## ★5. 統合 —— `ArithPlace` の相対基本等式 -/

open scoped Classical in
/-- ★★★★★★★**相対基本等式** —— `Σ_{V | v} localDeg(V) = [M : L]`。

★★これが `Φ^rlf` の条件 (a)（節点 `rlf-flat`）を算術で埋めるための鍵である ——
「上にある素点についての和」`h(y)(v) := Σ_{V | v} y(V)` が
`h ∘ (Φ.map α) = [M:L] · id` を与える。 -/
theorem sum_localDeg_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (v : ArithPlace L) :
    ∑ V ∈ (finite_fiber (L := L) (M := M) v).toFinset,
        localDeg (L := L) V = Module.finrank L M := by
  classical
  cases v with
  | inl v₀ =>
      have hset : (finite_fiber (L := L) (M := M) (Sum.inl v₀)).toFinset
          = (finite_fiber_fin (L := L) (M := M) v₀).toFinset.image Sum.inl := by
        ext V
        rw [Set.Finite.mem_toFinset, Finset.mem_image]
        constructor
        · intro hV
          cases V with
          | inl W =>
              refine ⟨W, ?_, rfl⟩
              rw [Set.Finite.mem_toFinset]
              simpa using hV
          | inr W => simp at hV
        · rintro ⟨W, hW, rfl⟩
          rw [Set.Finite.mem_toFinset] at hW
          show resPlace (L := L) (Sum.inl W) = Sum.inl v₀
          rw [resPlace_inl, hW]
      rw [hset, Finset.sum_image (fun _ _ _ _ h => Sum.inl_injective h)]
      exact sum_localDeg_fin_eq L M v₀
  | inr v₀ =>
      have hset : (finite_fiber (L := L) (M := M) (Sum.inr v₀)).toFinset
          = ({W : InfinitePlace M | resInf (L := L) W = v₀} : Finset _).image Sum.inr := by
        ext V
        rw [Set.Finite.mem_toFinset, Finset.mem_image]
        constructor
        · intro hV
          cases V with
          | inl W => simp at hV
          | inr W =>
              refine ⟨W, ?_, rfl⟩
              rw [mem_filter_univ]
              simpa using hV
        · rintro ⟨W, hW, rfl⟩
          rw [mem_filter_univ] at hW
          show resPlace (L := L) (Sum.inr W) = Sum.inr v₀
          rw [resPlace_inr, hW]
      rw [hset, Finset.sum_image (fun _ _ _ _ h => Sum.inr_injective h)]
      exact sum_localDeg_inf_eq L M v₀

/-! ## ★6. トレース写像 —— `Finsupp.mapDomain resPlace`

★★★`Φ^rlf` の条件 (a) を出すのに要る「`k` 倍の引き戻し」は、
**素点の制限に沿った `Finsupp.mapDomain`**（＝上にある素点についての和）である。 -/

open scoped Classical in
/-- ★`Finsupp.mapDomain` の値は**像の点でなくても**ファイバー上の和で書ける。

★mathlib の `mapDomain_apply_eq_sum` は `f a` の形の点だけなので、一般の点で置き直す。 -/
theorem mapDomain_apply_eq_sum' {α β N : Type} [DecidableEq β] [AddCommMonoid N] (f : α → β)
    (x : α →₀ N) (v : β) : (x.mapDomain f) v = ∑ i ∈ x.support with f i = v, x i := by
  classical
  simp [Finsupp.mapDomain, Finsupp.sum, Finsupp.single_apply, Finset.sum_ite]

open scoped Classical in
/-- ★★★★★★**トレース写像は `arithExtend` の `[M:L]` 倍の引き戻し**。

    mapDomain resPlace ∘ arithExtend = [M : L] · id

★これが節点 `rlf-flat`（`S = ℝ≥0`）を算術で埋める鍵である
（`isCharacteristicallyInjective_scMap_of_nsmul_retraction` に流す）。 -/
theorem mapDomain_arithExtend (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (x : ArithPlace L →₀ ℝ) :
    Finsupp.mapDomain (resPlace (L := L) (M := M)) (arithExtend (L := L) (M := M) x)
      = (Module.finrank L M : ℕ) • x := by
  classical
  ext v
  rw [mapDomain_apply_eq_sum']
  have hsub : ((arithExtend (L := L) (M := M) x).support.filter
      (fun V => resPlace (L := L) V = v))
      ⊆ (finite_fiber (L := L) (M := M) v).toFinset := by
    intro V hV
    rw [Finset.mem_filter] at hV
    rw [Set.Finite.mem_toFinset]
    exact hV.2
  have hzero : ∀ V ∈ (finite_fiber (L := L) (M := M) v).toFinset,
      V ∉ ((arithExtend (L := L) (M := M) x).support.filter
        (fun V => resPlace (L := L) V = v)) → (arithExtend (L := L) (M := M) x) V = 0 := by
    intro V hV hVn
    rw [Set.Finite.mem_toFinset] at hV
    by_contra hne
    exact hVn (Finset.mem_filter.mpr ⟨Finsupp.mem_support_iff.mpr hne, hV⟩)
  rw [Finset.sum_subset hsub hzero]
  have hterm : ∀ V ∈ (finite_fiber (L := L) (M := M) v).toFinset,
      (arithExtend (L := L) (M := M) x) V = (localDeg (L := L) V : ℝ) * x v := by
    intro V hV
    rw [Set.Finite.mem_toFinset] at hV
    show (localDeg (L := L) V : ℝ) * x (resPlace (L := L) V) = _
    rw [hV]
  rw [Finset.sum_congr rfl hterm, ← Finset.sum_mul]
  have hsum : ∑ V ∈ (finite_fiber (L := L) (M := M) v).toFinset, ((localDeg (L := L) V : ℝ))
      = ((Module.finrank L M : ℕ) : ℝ) := by
    rw [← Nat.cast_sum, sum_localDeg_eq L M v]
  rw [hsum]
  show _ = (Module.finrank L M : ℕ) • x v
  rw [nsmul_eq_mul]

/-! ## ★7. 相対ノルムの冪表示

★★トレース写像が `arithDivGroup`（有限素点の係数が `log N𝔭` の整数倍）を保つには

    N(𝔓) = N(𝔭)^{f(𝔓/𝔭)}

が要る。★mathlib の `Ideal.absNorm_eq_pow_inertiaDeg` は**底が ℤ**なので、
`ℤ ⊆ 𝓞 L ⊆ 𝓞 M` の塔で `f` の乗法性（`Ideal.inertiaDeg_algebra_tower`）を挟む。 -/

open NumberField Ideal in
/-- ★`𝓞 L` の非零素イデアルの下には有理素数がある。 -/
theorem exists_rational_prime_under (L : Type) [Field L] [NumberField L]
    (P : Ideal (𝓞 L)) [hp : P.IsPrime] (hne : P ≠ ⊥) :
    ∃ p : ℕ, p.Prime ∧ P.LiesOver (Ideal.span {(p : ℤ)}) := by
  classical
  haveI : (P.under ℤ).IsPrime := Ideal.IsPrime.under _ _
  have hqne : (P.under ℤ) ≠ ⊥ := Ideal.IsIntegral.comap_ne_bot ℤ hne
  obtain ⟨n, hspan⟩ := (IsPrincipalIdealRing.principal (P.under ℤ)).1
  have hsp : P.under ℤ = Ideal.span {n} := hspan
  have hn0 : n ≠ 0 := by
    rintro rfl
    refine hqne ?_
    rw [hsp, Ideal.span_singleton_eq_bot]
  have hnprime : Prime n := by
    rw [← Ideal.span_singleton_prime hn0, ← hsp]
    infer_instance
  refine ⟨n.natAbs, Int.prime_iff_natAbs_prime.mp hnprime, ⟨?_⟩⟩
  rw [hsp, Int.span_natAbs]

open NumberField Ideal in
open scoped Classical in
/-- ★★★★**相対ノルムの冪表示** —— `N(𝔓) = N(𝔭)^{f(𝔓/𝔭)}`。 -/
theorem absNorm_eq_pow_inertiaDeg_rel (L M : Type) [Field L] [Field M]
    [NumberField L] [NumberField M] [Algebra L M]
    (P : Ideal (𝓞 M)) [P.IsPrime] (hPne : P ≠ ⊥) :
    Ideal.absNorm P
      = (Ideal.absNorm (P.under (𝓞 L))) ^ ((P.under (𝓞 L)).inertiaDeg P) := by
  classical
  haveI hZ : IsScalarTower ℤ (𝓞 L) (𝓞 M) := inferInstance
  obtain ⟨p, hp, hlo⟩ := exists_rational_prime_under M P hPne
  haveI := hlo
  have hpz : Prime ((p : ℤ)) := Nat.prime_iff_prime_int.mp hp
  set q : Ideal (𝓞 L) := P.under (𝓞 L) with hqdef
  haveI : q.IsPrime := Ideal.IsPrime.under _ _
  haveI : P.LiesOver q := ⟨rfl⟩
  haveI hql : q.LiesOver (Ideal.span {(p : ℤ)}) :=
    ⟨by rw [hqdef, Ideal.under_under]; exact Ideal.LiesOver.over⟩
  haveI : P.IsMaximal := Ideal.IsPrime.isMaximal inferInstance hPne
  haveI : (Ideal.span {(p : ℤ)} : Ideal ℤ).IsMaximal :=
    PrincipalIdealRing.isMaximal_of_irreducible hpz.irreducible
  have hspan : (Ideal.span {(p : ℤ)} : Ideal ℤ) ≠ ⊥ := by simp [hp.ne_zero]
  have hqne : q ≠ ⊥ := Ideal.ne_bot_of_liesOver_of_ne_bot hspan q
  haveI : q.IsMaximal := Ideal.IsPrime.isMaximal inferInstance hqne
  have h1 : Ideal.absNorm P = p ^ ((Ideal.span {(p : ℤ)}).inertiaDeg P) := by
    have h := Ideal.absNorm_eq_pow_inertiaDeg (R := 𝓞 M) P hpz
    simpa using h
  have h2 : Ideal.absNorm q = p ^ ((Ideal.span {(p : ℤ)}).inertiaDeg q) := by
    have h := Ideal.absNorm_eq_pow_inertiaDeg (R := 𝓞 L) q hpz
    simpa using h
  have htow := Ideal.inertiaDeg_algebra_tower (Ideal.span {(p : ℤ)}) q P
  rw [h1, h2, ← pow_mul, htow]

/-! ## ★8. トレース写像は `arithDivGroup` を保つ -/

open NumberField Finset Finsupp Ideal in
open scoped Classical in
/-- ★トレース写像の有限素点での値は、上にある有限素点についての和。 -/
theorem mapDomain_apply_inl (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (y : ArithPlace M →₀ ℝ) (w : FinitePlace L) :
    (Finsupp.mapDomain (resPlace (L := L) (M := M)) y) (Sum.inl w)
      = ∑ W ∈ (finite_fiber_fin (L := L) (M := M) w).toFinset, y (Sum.inl W) := by
  classical
  rw [mapDomain_apply_eq_sum']
  have hsub : (y.support.filter (fun V => resPlace (L := L) V = Sum.inl w))
      ⊆ (finite_fiber_fin (L := L) (M := M) w).toFinset.image Sum.inl := by
    intro V hV
    rw [Finset.mem_filter] at hV
    obtain ⟨hVs, hVr⟩ := hV
    cases V with
    | inl W =>
        refine Finset.mem_image.mpr ⟨W, ?_, rfl⟩
        rw [Set.Finite.mem_toFinset]
        simpa using hVr
    | inr W => simp at hVr
  have hzero : ∀ V ∈ (finite_fiber_fin (L := L) (M := M) w).toFinset.image Sum.inl,
      V ∉ (y.support.filter (fun V => resPlace (L := L) V = Sum.inl w)) → y V = 0 := by
    intro V hV hVn
    obtain ⟨W, hW, rfl⟩ := Finset.mem_image.mp hV
    rw [Set.Finite.mem_toFinset] at hW
    by_contra hne
    refine hVn (Finset.mem_filter.mpr ⟨Finsupp.mem_support_iff.mpr hne, ?_⟩)
    show resPlace (L := L) (Sum.inl W) = Sum.inl w
    rw [resPlace_inl]
    exact congrArg Sum.inl hW
  rw [Finset.sum_subset hsub hzero,
    Finset.sum_image (fun _ _ _ _ h => Sum.inl_injective h)]

open NumberField Finset Finsupp Ideal IsDedekindDomain in
open scoped Classical in
/-- ★上にある素点の係数は `f` 倍だけずれる（相対ノルムの冪表示から）。 -/
theorem trace_term_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (y : ArithPlace M →₀ ℝ) (n : FinitePlace M → ℤ)
    (hn : ∀ W : FinitePlace M, y (Sum.inl W)
      = (n W : ℝ) * Real.log (Ideal.absNorm (FinitePlace.maximalIdeal W).asIdeal : ℝ))
    (w : FinitePlace L) (W : FinitePlace M) (hW : resFin (L := L) W = w) :
    y (Sum.inl W) = (n W : ℝ)
      * ((w.maximalIdeal.asIdeal.inertiaDeg W.maximalIdeal.asIdeal : ℤ) : ℝ)
      * Real.log (Ideal.absNorm (FinitePlace.maximalIdeal w).asIdeal : ℝ) := by
  classical
  have hres : resHOS (L := L) W.maximalIdeal = w.maximalIdeal := by
    rw [← hW, resFin, FinitePlace.maximalIdeal_mk]
  have hunder : W.maximalIdeal.asIdeal.under (𝓞 L) = w.maximalIdeal.asIdeal :=
    congrArg IsDedekindDomain.HeightOneSpectrum.asIdeal hres
  have hPne : W.maximalIdeal.asIdeal ≠ ⊥ := W.maximalIdeal.ne_bot
  haveI : W.maximalIdeal.asIdeal.IsPrime := W.maximalIdeal.isPrime
  have hnorm := absNorm_eq_pow_inertiaDeg_rel L M W.maximalIdeal.asIdeal hPne
  rw [hunder] at hnorm
  rw [hn W, hnorm]
  push_cast
  rw [Real.log_pow]
  ring

open NumberField Finset Finsupp Ideal in
open scoped Classical in
/-- ★★★★★**トレース写像は `arithDivGroup` を保つ**。

★有限素点の係数が `log N𝔭` の整数倍であることは、上にある素点では
`log N𝔓 = f · log N𝔭`（`absNorm_eq_pow_inertiaDeg_rel`）なので保たれる。 -/
theorem mapDomain_mem_arithDivGroup (L M : Type) [Field L] [Field M]
    [NumberField L] [NumberField M] [Algebra L M]
    (y : ArithPlace M →₀ ℝ) (hy : y ∈ arithDivGroup M) :
    Finsupp.mapDomain (resPlace (L := L) (M := M)) y ∈ arithDivGroup L := by
  classical
  choose n hn using hy
  intro w
  rw [mapDomain_apply_inl]
  refine ⟨∑ W ∈ (finite_fiber_fin (L := L) (M := M) w).toFinset,
    n W * ((w.maximalIdeal.asIdeal.inertiaDeg W.maximalIdeal.asIdeal : ℤ)), ?_⟩
  push_cast
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl (fun W hW => ?_)
  rw [Set.Finite.mem_toFinset] at hW
  have h := trace_term_eq L M y n hn w W hW
  push_cast at h
  exact h

open NumberField Finset Finsupp in
open scoped Classical in
/-- ★トレース写像は非負性を保つ（非負元の有限和だから）。 -/
theorem mapDomain_nonneg (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] {y : ArithPlace M →₀ ℝ} (hy : ∀ s, 0 ≤ y s) :
    ∀ v, 0 ≤ (Finsupp.mapDomain (resPlace (L := L) (M := M)) y) v := by
  classical
  intro v
  rw [mapDomain_apply_eq_sum']
  exact Finset.sum_nonneg (fun V _ => hy V)

/-! ### ★出典の紐付け -/

def mapDomain_mem_arithDivGroup.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — トレース写像は arithDivGroup を保つ",
    sectionId := "frdi-thm-6-4" }

def absNorm_eq_pow_inertiaDeg_rel.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 相対ノルムの冪表示 N(𝔓) = N(𝔭)^f",
    sectionId := "frdi-thm-6-4" }

def mapDomain_arithExtend.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — トレース写像は arithExtend の [M:L] 倍の引き戻し",
    sectionId := "frdi-thm-6-4" }

def sum_localDeg_fin_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 有限素点の相対基本等式(Σ_{W|v} e·f = [M:L])",
    sectionId := "frdi-thm-6-4" }

def sum_localDeg_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 相対基本等式(Σ_{V|v} localDeg V = [M:L])",
    sectionId := "frdi-thm-6-4" }

def sum_mult_comap_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 無限素点の相対基本等式(Σ_{w|v} mult w = mult v · [M:L])",
    sectionId := "frdi-thm-6-4" }

def sum_localDeg_inf_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 無限素点の相対基本等式(localDeg の形)",
    sectionId := "frdi-thm-6-4" }

def sum_localDeg_inf_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "AlgHom.card(#(M →ₐ[L] ℂ) = [M:L])"
      (.inMathlib "AlgHom.card") 114,
    .citation "[mathlib]" "InfinitePlace.card_filter_mk_eq(#{φ | mk φ = w} = mult w)"
      (.inMathlib "NumberField.InfinitePlace.card_filter_mk_eq") 114,
    .citation "[mathlib]" "InfinitePlace.comap_mk((mk ψ).comap f = mk (ψ ∘ f))"
      (.inMathlib "NumberField.InfinitePlace.comap_mk") 114,
    .implicitStep
      ("★mathlib の相対版(card_isUnramified ほか)は Galois 専用なので、" ++
       "埋め込みの二重数え上げで作った") 114 ]

end ABC3.Found.Divisor
