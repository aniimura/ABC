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

/-! ## witness の押し出し——almost étale 性を別の環へ移す

`Theorem 3.2` の *"`B∞` is almost isomorphic to `S∞ ⊗_{R∞} A∞` as this is
almost étale over `S∞`"* のように、almost étale 性を別の環へ移す段で使う。
`Definition 2.1` 条件(iii)の内容(annihilation と augmentation)は
**全射 `A`-代数写像で押し出せる**。 -/

/-- `lmul'` は `A`-代数写像と可換。 -/
theorem lmul'_map {A B B' : Type u} [CommRing A] [CommRing B] [CommRing B']
    [Algebra A B] [Algebra A B'] (f : B →ₐ[A] B') (w : TensorProduct A B B) :
    Algebra.TensorProduct.lmul' A (Algebra.TensorProduct.map f f w)
      = f (Algebra.TensorProduct.lmul' A w) := by
  induction w using TensorProduct.induction_on with
  | zero => simp
  | tmul u v =>
    rw [Algebra.TensorProduct.map_tmul, Algebra.TensorProduct.lmul'_apply_tmul,
      Algebra.TensorProduct.lmul'_apply_tmul, map_mul]
  | add x y hx hy => rw [map_add, map_add, map_add, hx, hy, map_add]

/-- **witness は全射 `A`-代数写像で押し出せる**。`f : B ↠ B'` が全射なら
`w' := (f ⊗ f)(w)` は**同じ水準で** annihilation と augmentation を満たす。 -/
theorem witness_pushforward {A B B' : Type u} [CommRing A] [CommRing B] [CommRing B']
    [Algebra A B] [Algebra A B'] (f : B →ₐ[A] B') (hf : Function.Surjective f)
    (c : A) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = c • (1 : B)) :
    (∀ q' : B', (1 ⊗ₜ[A] q' - q' ⊗ₜ[A] 1) * (Algebra.TensorProduct.map f f w) = 0)
    ∧ Algebra.TensorProduct.lmul' A (Algebra.TensorProduct.map f f w) = c • (1 : B') := by
  refine ⟨fun q' => ?_, ?_⟩
  · obtain ⟨q, rfl⟩ := hf q'
    have h := congrArg (Algebra.TensorProduct.map f f) (hw_ann q)
    rw [map_mul, map_zero, map_sub, Algebra.TensorProduct.map_tmul,
      Algebra.TensorProduct.map_tmul, map_one] at h
    exact h
  · rw [lmul'_map f w, hw_aug, map_smul, map_one]

/-! ## witness の底変換——almost étale 性は `⊗_R` で保たれる

`Theorem 3.2` の *"It follows that `B∞` is almost isomorphic to
`S∞ ⊗_{R∞} A∞` **as this is almost étale over `S∞`**"* の括弧内が
これである:`A∞/R∞` が almost étale なら、その底変換
`S∞ ⊗_{R∞} A∞` は `S∞` 上 almost étale。

`witness_pushforward`(全射で押し出す)と対になる操作で、こちらは
**任意の** `A`-代数 `A'` へ移せる。水準 `c` は
`algebraMap A A' c` に移るだけで悪化しない。 -/

open scoped TensorProduct in
/-- 底変換 `B ⊗_A B → (A' ⊗_A B) ⊗_{A'} (A' ⊗_A B)`。
`b₁ ⊗ b₂ ↦ (1 ⊗ b₁) ⊗ (1 ⊗ b₂)`。 -/
noncomputable def bcWitnessMap (A A' B : Type u) [CommRing A] [CommRing A'] [CommRing B]
    [Algebra A A'] [Algebra A B] :
    (TensorProduct A B B) →ₐ[A]
      (TensorProduct A' (TensorProduct A A' B) (TensorProduct A A' B)) :=
  Algebra.TensorProduct.lift
    (((Algebra.TensorProduct.includeLeft :
        (TensorProduct A A' B) →ₐ[A']
          TensorProduct A' (TensorProduct A A' B) (TensorProduct A A' B)).restrictScalars A).comp
      (Algebra.TensorProduct.includeRight : B →ₐ[A] TensorProduct A A' B))
    (((Algebra.TensorProduct.includeRight :
        (TensorProduct A A' B) →ₐ[A']
          TensorProduct A' (TensorProduct A A' B) (TensorProduct A A' B)).restrictScalars A).comp
      (Algebra.TensorProduct.includeRight : B →ₐ[A] TensorProduct A A' B))
    (fun _ _ => Commute.all _ _)

open scoped TensorProduct in
@[simp] theorem bcWitnessMap_tmul (A A' B : Type u) [CommRing A] [CommRing A'] [CommRing B]
    [Algebra A A'] [Algebra A B] (b₁ b₂ : B) :
    bcWitnessMap A A' B (b₁ ⊗ₜ[A] b₂)
      = ((1 : A') ⊗ₜ[A] b₁) ⊗ₜ[A'] ((1 : A') ⊗ₜ[A] b₂) := by
  simp [bcWitnessMap, Algebra.TensorProduct.includeLeft_apply,
    Algebra.TensorProduct.includeRight_apply, Algebra.TensorProduct.tmul_mul_tmul]

open scoped TensorProduct in
/-- augmentation は底変換と可換——`lmul'` を取ってから `A'` を足しても、
足してから `lmul'` を取っても同じ。 -/
theorem bcWitnessMap_lmul' (A A' B : Type u) [CommRing A] [CommRing A'] [CommRing B]
    [Algebra A A'] [Algebra A B] (w : TensorProduct A B B) :
    Algebra.TensorProduct.lmul' A' (bcWitnessMap A A' B w)
      = Algebra.TensorProduct.includeRight (Algebra.TensorProduct.lmul' A w) := by
  induction w using TensorProduct.induction_on with
  | zero => simp
  | tmul b₁ b₂ =>
    rw [bcWitnessMap_tmul, Algebra.TensorProduct.lmul'_apply_tmul,
      Algebra.TensorProduct.lmul'_apply_tmul, Algebra.TensorProduct.tmul_mul_tmul,
      Algebra.TensorProduct.includeRight_apply, one_mul]
  | add x y hx hy => rw [map_add, map_add, hx, hy, map_add, map_add]

open scoped TensorProduct in
/-- **witness は任意の底変換で押し出せる**——`Definition 2.1` 条件(iii)の
内容(annihilation と augmentation)は `A → A'` で移り、水準は
`c ↦ algebraMap A A' c` になるだけ。

原文 *"`B∞` is almost isomorphic to `S∞ ⊗_{R∞} A∞` **as this is almost
étale over `S∞`**"* の括弧内。 -/
theorem witness_baseChange (A A' B : Type u) [CommRing A] [CommRing A'] [CommRing B]
    [Algebra A A'] [Algebra A B]
    (c : A) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, ((1 : B) ⊗ₜ[A] q - q ⊗ₜ[A] (1 : B)) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = c • (1 : B)) :
    (∀ q' : TensorProduct A A' B,
        ((1 : TensorProduct A A' B) ⊗ₜ[A'] q'
          - q' ⊗ₜ[A'] (1 : TensorProduct A A' B)) * bcWitnessMap A A' B w = 0)
    ∧ Algebra.TensorProduct.lmul' A' (bcWitnessMap A A' B w)
        = (algebraMap A A' c) • (1 : TensorProduct A A' B) := by
  constructor
  · intro q'
    induction q' using TensorProduct.induction_on with
    | zero => simp
    | tmul a' b =>
      have hq : (a' ⊗ₜ[A] b : TensorProduct A A' B) = a' • ((1 : A') ⊗ₜ[A] b) := by
        rw [TensorProduct.smul_tmul', smul_eq_mul, mul_one]
      have h0 : ((1 : TensorProduct A A' B) ⊗ₜ[A'] ((1 : A') ⊗ₜ[A] b)
            - ((1 : A') ⊗ₜ[A] b) ⊗ₜ[A'] (1 : TensorProduct A A' B))
          * bcWitnessMap A A' B w = 0 := by
        have e1 : bcWitnessMap A A' B ((1 : B) ⊗ₜ[A] b)
            = (1 : TensorProduct A A' B) ⊗ₜ[A'] ((1 : A') ⊗ₜ[A] b) := by
          rw [bcWitnessMap_tmul]; rfl
        have e2 : bcWitnessMap A A' B (b ⊗ₜ[A] (1 : B))
            = ((1 : A') ⊗ₜ[A] b) ⊗ₜ[A'] (1 : TensorProduct A A' B) := by
          rw [bcWitnessMap_tmul]; rfl
        rw [← e1, ← e2, ← map_sub, ← map_mul, hw_ann b, map_zero]
      have hE : ((1 : TensorProduct A A' B) ⊗ₜ[A'] (a' ⊗ₜ[A] b)
            - (a' ⊗ₜ[A] b) ⊗ₜ[A'] (1 : TensorProduct A A' B))
          = a' • ((1 : TensorProduct A A' B) ⊗ₜ[A'] ((1 : A') ⊗ₜ[A] b)
            - ((1 : A') ⊗ₜ[A] b) ⊗ₜ[A'] (1 : TensorProduct A A' B)) := by
        rw [smul_sub,
          ← TensorProduct.tmul_smul a' (1 : TensorProduct A A' B) ((1 : A') ⊗ₜ[A] b),
          TensorProduct.smul_tmul' a' ((1 : A') ⊗ₜ[A] b) (1 : TensorProduct A A' B),
          ← hq]
      rw [hE, smul_mul_assoc, h0, smul_zero]
    | add x y hx hy =>
      have hsplit : ((1 : TensorProduct A A' B) ⊗ₜ[A'] (x + y)
            - (x + y) ⊗ₜ[A'] (1 : TensorProduct A A' B))
          = ((1 : TensorProduct A A' B) ⊗ₜ[A'] x - x ⊗ₜ[A'] (1 : TensorProduct A A' B))
            + ((1 : TensorProduct A A' B) ⊗ₜ[A'] y
              - y ⊗ₜ[A'] (1 : TensorProduct A A' B)) := by
        rw [TensorProduct.tmul_add, TensorProduct.add_tmul]; ring
      rw [hsplit, add_mul, hx, hy, add_zero]
  · rw [bcWitnessMap_lmul', hw_aug, map_smul, map_one, algebraMap_smul]

open scoped TensorProduct in
/-- **分裂エタール代数 `ι → A` の分離冪等元**——`Definition 2.1` の
witness 条件(annihilation・augmentation)を `c = 1` で満たす具体例。

`w = Σᵢ eᵢ ⊗ eᵢ`(`eᵢ = Pi.single i 1`)。`q·eᵢ = q(i)·eᵢ` なので
`eᵢ ⊗ (q·eᵢ) = q(i)·(eᵢ ⊗ eᵢ) = (q·eᵢ) ⊗ eᵢ`。

★これは Falt1 全体で「almost étale の仮定が空虚でない」ことを示すための
**具体的な witness の在庫**であり、非空虚性の対照に繰り返し使える。 -/
theorem pi_witness {A : Type u} [CommRing A] {ι : Type u} [Fintype ι] [DecidableEq ι] :
    (∀ q : ι → A, ((1 : ι → A) ⊗ₜ[A] q - q ⊗ₜ[A] (1 : ι → A))
        * (∑ i : ι, (Pi.single i (1 : A)) ⊗ₜ[A] (Pi.single i (1 : A))) = 0)
    ∧ Algebra.TensorProduct.lmul' A
        (∑ i : ι, (Pi.single i (1 : A)) ⊗ₜ[A] (Pi.single i (1 : A)))
      = (1 : A) • (1 : ι → A) := by
  have hmul : ∀ (q : ι → A) (i : ι), q * Pi.single i 1 = q i • Pi.single i (1 : A) := by
    intro q i; funext j; by_cases h : j = i <;> simp [h]
  have hsum : ∑ i : ι, Pi.single i (1 : A) = 1 := by
    funext j; simp [Finset.sum_apply]
  constructor
  · intro q
    rw [Finset.mul_sum]
    refine Finset.sum_eq_zero (fun i _ => ?_)
    rw [sub_mul, Algebra.TensorProduct.tmul_mul_tmul, Algebra.TensorProduct.tmul_mul_tmul,
      one_mul, hmul q i, TensorProduct.tmul_smul, ← TensorProduct.smul_tmul', sub_self]
  · rw [map_sum, one_smul]
    simp only [Algebra.TensorProduct.lmul'_apply_tmul]
    calc ∑ i : ι, Pi.single i (1 : A) * Pi.single i (1 : A)
        = ∑ i : ι, Pi.single i (1 : A) :=
          Finset.sum_congr rfl (fun i _ => by funext j; by_cases h : j = i <;> simp [h])
      _ = 1 := hsum

open scoped TensorProduct in
/-- 非空虚性——`A = ℤ[X]`、`B = Fin 2 → A`(分裂エタール、`pi_witness`)を
`A' = A[Y]` へ底変換する。`witness_baseChange` の仮定が空虚に真に
なっていないことの対照。 -/
example :
    Algebra.TensorProduct.lmul' (Polynomial (Polynomial ℤ))
      (bcWitnessMap (Polynomial ℤ) (Polynomial (Polynomial ℤ)) (Fin 2 → Polynomial ℤ)
        (∑ i : Fin 2, (Pi.single i (1 : Polynomial ℤ)) ⊗ₜ[Polynomial ℤ]
          (Pi.single i (1 : Polynomial ℤ))))
      = algebraMap (Polynomial ℤ) (Polynomial (Polynomial ℤ)) 1 • 1 :=
  (witness_baseChange (Polynomial ℤ) (Polynomial (Polynomial ℤ)) (Fin 2 → Polynomial ℤ)
    1 _ pi_witness.1 pi_witness.2).2

/-- 非空虚性——`R = ℤ[X]`、`A = Fin 2 → R`(非自明な有限自由拡大)、
`τ =` 第 0 成分への射影(`τ 1 = 1`)。仮定が空虚に真になっていないことの対照。 -/
example : (1 : Polynomial ℤ) • (1 : Polynomial ℤ) ∈
    LinearMap.range (LinearMap.id : Polynomial ℤ →ₗ[Polynomial ℤ] Polynomial ℤ) :=
  almost_descent_range
    (LinearMap.proj 0 : (Fin 2 → Polynomial ℤ) →ₗ[Polynomial ℤ] Polynomial ℤ)
    1 rfl LinearMap.id 1 ⟨(1 : Polynomial ℤ) ⊗ₜ[Polynomial ℤ] 1, by simp⟩

end ABC3.Found.Falt1
