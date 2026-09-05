import ABC3.Found.Falt1.AlmostDerivation

/-!
# [Falt1] Chapter I §3 —— almost 直和因子からの降下(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §3、
`Theorem 3.2`(物理 p.11 = 印字 p.264)。

> *... Then `A∞` is almost étale over `R∞`, and **`R∞` is almost a direct
> summand in it (the image of the trace form contains `m·R∞`)**. It follows
> that `B∞` is almost isomorphic to `S∞ ⊗_{R∞} A∞` as this is almost étale
> over `S∞`. So enlarging `ε` a little we see that `p^ε·e_{S∞/R∞}` becomes
> integral after `⊗_R A∞`. **As `R∞` is almost a direct summand in `A∞`,
> this already holds before tensoring.** The assertion now follows as `ε`
> can be arbitrarily small.*

## ロードマップの訂正(2026-09-05)

`falt1-goal.md` はこれまで「§3 は `Theorem 2.2`/`2.3` の壁を継承する」と
記録していたが、原典を逐語で読み直すとそうではない:

* `Theorem 3.1` は `dim(R)` の帰納法で、底(`dim = 1`)で使うのは
  **`Theorem 1.2`**(*"the sequence `Rₙ` kills ramification (complete and
  use part 1, Theorem 1.2)"*)である。
* `Theorem 3.2` は `Definition 2.1`・trace 形式の像・「almost 直和因子」
  の議論だけを使う。

本ファイルはそのうち**再利用可能な代数的核**——「almost 直和因子からの
降下」——を切り出す。

## 内容

`τ : A → R` を `R`-線形写像とし `τ(1) = c` とする(trace 形式の像が
`c` を含む場合に `τ := tr(a₀ · −)` と取ればよい)。このとき

    retr : M ⊗_R A → M,  m ⊗ a ↦ τ(a)·m

は `retr(x ⊗ 1) = c·x` を満たし、任意の `f : N → M` について
`retr ∘ (f ⊗ id) = f ∘ retr` が成り立つ。したがって `x ⊗ 1` が
`f ⊗ id` の像に入れば `c·x` は `f` の像に入る——これが原文の
*"this already holds before tensoring"* の内容である。
-/

namespace ABC3.Found.Falt1

universe u

/-- **almost 直和因子からの降下**(Faltings `Theorem 3.2` の
*"As `R∞` is almost a direct summand in `A∞`, this already holds before
tensoring"*)。`τ(1) = c` なら、`x ⊗ 1` が `f ⊗ id` の像に入るとき
`c·x` は `f` の像に入る。 -/
theorem almost_descent_range {R A M N : Type u} [CommRing R] [CommRing A] [Algebra R A]
    [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
    (τ : A →ₗ[R] R) (c : R) (hτ : τ 1 = c)
    (f : N →ₗ[R] M) (x : M)
    (hx : x ⊗ₜ[R] (1 : A)
      ∈ LinearMap.range (TensorProduct.map f (LinearMap.id : A →ₗ[R] A))) :
    c • x ∈ LinearMap.range f := by
  set retrM : TensorProduct R M A →ₗ[R] M :=
    (TensorProduct.rid R M).toLinearMap ∘ₗ TensorProduct.map LinearMap.id τ with hrM
  set retrN : TensorProduct R N A →ₗ[R] N :=
    (TensorProduct.rid R N).toLinearMap ∘ₗ TensorProduct.map LinearMap.id τ with hrN
  obtain ⟨z, hz⟩ := hx
  refine ⟨retrN z, ?_⟩
  have hcomp : ∀ u : TensorProduct R N A,
      retrM (TensorProduct.map f (LinearMap.id : A →ₗ[R] A) u) = f (retrN u) := by
    intro u
    induction u using TensorProduct.induction_on with
    | zero => simp
    | tmul n a => simp [hrM, hrN, map_smul]
    | add u v hu hv => rw [map_add, map_add, hu, hv, map_add, map_add]
  have h1 := congrArg retrM hz
  rw [hcomp] at h1
  have h2 : retrM (x ⊗ₜ[R] (1:A)) = c • x := by simp [hrM, hτ]
  rw [h1, h2]

/-- trace 形式の像が `c` を含む場合の言い換え——原文の
*"the image of the trace form contains `m·R∞`"* にそのまま対応する。 -/
theorem almost_descent_range_of_trace {R A M N : Type u} [CommRing R] [CommRing A]
    [Algebra R A] [Module.Finite R A] [Module.Free R A]
    [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
    (c : R) (a₀ : A) (hc : Algebra.trace R A a₀ = c)
    (f : N →ₗ[R] M) (x : M)
    (hx : x ⊗ₜ[R] (1 : A)
      ∈ LinearMap.range (TensorProduct.map f (LinearMap.id : A →ₗ[R] A))) :
    c • x ∈ LinearMap.range f :=
  almost_descent_range ((Algebra.trace R A) ∘ₗ (LinearMap.mulLeft R a₀)) c
    (by simp only [LinearMap.coe_comp, Function.comp_apply, LinearMap.mulLeft_apply,
      mul_one]; exact hc) f x hx

/-! ## 「almost 直和因子である」を almost étale 性から導く

原文は *"`A∞` is almost étale over `R∞`, and **`R∞` is almost a direct
summand in it (the image of the trace form contains `m·R∞`)**"* と、
almost 直和因子性を almost étale 性の帰結として使う。

`AlmostEtale.trace_ideal_pow_mem_traceIdeal`(既証)は
`span{p^{n·rank}} ≤ span(range tr)` を与えるが、`Algebra.trace` は
`A`-線型なのでその像は既に `A`-部分加群であり、生成するイデアルと一致する。
したがって**実際に `p^{n·rank} = tr(z)` となる `z` が取れる**——
これが `almost_descent_range_of_trace` の入力にそのまま入る。 -/

attribute [local instance] FractionRing.liftAlgebra in
/-- **trace 形式の像が `p^{n·rank}` を含む**(原文の
*"the image of the trace form contains `m·R∞`"*)。 -/
theorem exists_trace_eq_pow {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Finite A B] [Module.Free A B] (p : A)
    (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    [IsDedekindDomain A] [IsDedekindDomain B] [Module.IsTorsionFree A B]
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (e : TensorProduct A B B)
    (he : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p e
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p))) :
    ∃ z : B, Algebra.trace A B z
      = p ^ (n * Module.finrank (FractionRing A) (FractionRing B)) := by
  have hsub := trace_ideal_pow_mem_traceIdeal p hAE hf0inj n e he
  have hmem : p ^ (n * Module.finrank (FractionRing A) (FractionRing B))
      ∈ (Ideal.span (Set.range (Algebra.trace A B)) : Ideal A) :=
    hsub (Ideal.mem_span_singleton_self _)
  have hspan : (Ideal.span (Set.range (Algebra.trace A B)) : Ideal A)
      = LinearMap.range (Algebra.trace A B) :=
    Submodule.span_eq (LinearMap.range (Algebra.trace A B))
  rw [hspan] at hmem
  exact hmem

attribute [local instance] FractionRing.liftAlgebra in
/-- **`Theorem 3.2` の降下段**——*"As `R∞` is almost a direct summand in
`A∞`, this already holds before tensoring."*

★「almost 直和因子である」を**仮定せず、almost étale 性から導いて**いる
(`exists_trace_eq_pow`)のが要点。 -/
theorem thm_3_2_descent {A B M N : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Finite A B] [Module.Free A B] (p : A)
    (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    [IsDedekindDomain A] [IsDedekindDomain B] [Module.IsTorsionFree A B]
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (e : TensorProduct A B B)
    (he : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p e
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p)))
    [AddCommGroup M] [Module A M] [AddCommGroup N] [Module A N]
    (f : N →ₗ[A] M) (x : M)
    (hx : x ⊗ₜ[A] (1 : B)
      ∈ LinearMap.range (TensorProduct.map f (LinearMap.id : B →ₗ[A] B))) :
    (p ^ (n * Module.finrank (FractionRing A) (FractionRing B))) • x
      ∈ LinearMap.range f := by
  obtain ⟨z, hz⟩ := exists_trace_eq_pow p hAE hf0inj n e he
  exact almost_descent_range_of_trace _ z hz f x hx

/-- 非空虚性——`R = ℤ[X]`、`A = Fin 2 → R`(非自明な有限自由拡大)、
`τ =` 第 0 成分への射影(`τ 1 = 1`)。仮定が空虚に真になっていないことの対照。 -/
example : (1 : Polynomial ℤ) • (1 : Polynomial ℤ) ∈
    LinearMap.range (LinearMap.id : Polynomial ℤ →ₗ[Polynomial ℤ] Polynomial ℤ) :=
  almost_descent_range
    (LinearMap.proj 0 : (Fin 2 → Polynomial ℤ) →ₗ[Polynomial ℤ] Polynomial ℤ)
    1 rfl LinearMap.id 1 ⟨(1 : Polynomial ℤ) ⊗ₜ[Polynomial ℤ] 1, by simp⟩

end ABC3.Found.Falt1
