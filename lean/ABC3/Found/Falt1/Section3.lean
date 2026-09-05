import ABC3.Found.Falt1.AlmostDerivation
import ABC3.Found.Falt1.Section1

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

/-! ## ★単元の `N` 乗根を添加したときの微分(2026-09-05)

原文 §3(a):

> *We let `Rₙ` denote the normalization of `R ⊗_V Vₙ` in the extension
> generated by the `p^n`-th roots of the `uᵢ`. **It is clear that**
> `Rₙ ≅ (R ⊗_V Vₙ)[X₁,…,X_d]/(Xᵢ^{p^n} − uᵢ)`.*

その微分の計算をしておく。`x^N = u`(単元)から
`0 = D(u) = N·x^{N−1}·Dx` で、`x^{N−1}` は単元だから **`N·Dx = 0`**。
`Ω` が `Dx` たちで生成される(`kaehler_span_of_adjoin`)ので、
`N` は `Ω[B⁄A]` 全体を零化する。

§4(c) の *"modulo `p` any element of `R∞` is a `p`-power"* や
固有空間分解の計算はこの形の Kummer 拡大の微分を使う。 -/

/-- **単元の `N` 乗根を添加すると微分は `N` で消える**(1 元版)。 -/
theorem kaehler_annihilated_of_root_of_unit {A B : Type u} [CommRing A] [CommRing B]
    [Algebra A B] (x : B) (N : ℕ) (hN : 0 < N) (u : A)
    (hxu : x ^ N = algebraMap A B u) (huniv : IsUnit (algebraMap A B u))
    (hcyc : Submodule.span B ({KaehlerDifferential.D A B x} : Set (Ω[B⁄A])) = ⊤) :
    ∀ ω : Ω[B⁄A], (N : ℕ) • ω = 0 := by
  have hsplit : x * x ^ (N - 1) = x ^ N := by
    rw [← pow_succ']
    congr 1
    omega
  have hxunit : IsUnit x := by
    refine isUnit_of_mul_isUnit_left (y := x ^ (N - 1)) ?_
    rw [hsplit, hxu]
    exact huniv
  obtain ⟨v, hv⟩ := hxunit.pow (N - 1)
  have hDx : (N : ℕ) • KaehlerDifferential.D A B x = 0 := by
    have h0 : (N : ℕ) • x ^ (N - 1) • KaehlerDifferential.D A B x = 0 := by
      rw [← (KaehlerDifferential.D A B).leibniz_pow x N, hxu]
      exact (KaehlerDifferential.D A B).map_algebraMap u
    have hinv : ((v⁻¹ : Bˣ) : B) * x ^ (N - 1) = 1 := by rw [← hv]; exact v.inv_mul
    have hkey := congrArg (fun z => ((v⁻¹ : Bˣ) : B) • z) h0
    simp only [smul_zero] at hkey
    rw [smul_comm, ← mul_smul, hinv, one_smul] at hkey
    exact hkey
  intro ω
  have hmem : ω ∈ Submodule.span B ({KaehlerDifferential.D A B x} : Set (Ω[B⁄A])) := by
    rw [hcyc]; exact Submodule.mem_top
  obtain ⟨b, hb⟩ := Submodule.mem_span_singleton.mp hmem
  rw [← hb, smul_comm, hDx, smul_zero]

/-- **単元の `N` 乗根を `d` 個添加した版**——原文の
`Rₙ = (R ⊗_V Vₙ)[X₁,…,X_d]/(Xᵢ^{p^n} − uᵢ)` の形。
`hspan` は `kaehler_span_of_adjoin`(`Section1.lean`)が与える。 -/
theorem kaehler_annihilated_of_roots_of_units {A B ι : Type u} [CommRing A] [CommRing B]
    [Algebra A B] (x : ι → B) (N : ℕ) (hN : 0 < N) (u : ι → A)
    (hxu : ∀ i, (x i) ^ N = algebraMap A B (u i))
    (huniv : ∀ i, IsUnit (algebraMap A B (u i)))
    (hspan : Submodule.span B (KaehlerDifferential.D A B '' (Set.range x)) = ⊤) :
    ∀ ω : Ω[B⁄A], (N : ℕ) • ω = 0 := by
  have hDx : ∀ i, (N : ℕ) • KaehlerDifferential.D A B (x i) = 0 := by
    intro i
    have hsplit : (x i) * (x i) ^ (N - 1) = (x i) ^ N := by
      rw [← pow_succ']
      congr 1
      omega
    have hxunit : IsUnit (x i) := by
      refine isUnit_of_mul_isUnit_left (y := (x i) ^ (N - 1)) ?_
      rw [hsplit, hxu i]
      exact huniv i
    obtain ⟨v, hv⟩ := hxunit.pow (N - 1)
    have h0 : (N : ℕ) • (x i) ^ (N - 1) • KaehlerDifferential.D A B (x i) = 0 := by
      rw [← (KaehlerDifferential.D A B).leibniz_pow (x i) N, hxu i]
      exact (KaehlerDifferential.D A B).map_algebraMap (u i)
    have hinv : ((v⁻¹ : Bˣ) : B) * (x i) ^ (N - 1) = 1 := by rw [← hv]; exact v.inv_mul
    have hkey := congrArg (fun z => ((v⁻¹ : Bˣ) : B) • z) h0
    simp only [smul_zero] at hkey
    rw [smul_comm, ← mul_smul, hinv, one_smul] at hkey
    exact hkey
  intro ω
  have hmem : ω ∈ Submodule.span B (KaehlerDifferential.D A B '' (Set.range x)) := by
    rw [hspan]; exact Submodule.mem_top
  refine Submodule.span_induction ?_ ?_ ?_ ?_ hmem
  · rintro _ ⟨_, ⟨i, rfl⟩, rfl⟩
    exact hDx i
  · simp
  · intro y z _ _ hy hz; rw [smul_add, hy, hz, add_zero]
  · intro c y _ hy; rw [smul_comm, hy, smul_zero]

/-! ### ★原文 §3(a) の `Rₙ` を実際に構成する(2026-09-05)

原文が *"It is clear that `Rₙ ≅ (R ⊗_V Vₙ)[X₁,…,X_d]/(Xᵢ^{p^n} − uᵢ)`"*
と書く環を、Lean で**そのまま定義**して微分を計算する。 -/

/-- **原文 §3(a) の `Rₙ`**——`S[X₁,…,X_d]/(Xᵢ^N − uᵢ)`。 -/
noncomputable abbrev kummerAlgebra (S : Type u) [CommRing S] (ι : Type u) (N : ℕ) (u : ι → S) :
    Type u :=
  MvPolynomial ι S ⧸ Ideal.span (Set.range fun i =>
    (MvPolynomial.X i : MvPolynomial ι S) ^ N - MvPolynomial.C (u i))

/-- 生成元 `xᵢ = Xᵢ` の像。 -/
noncomputable def kummerGen (S : Type u) [CommRing S] (ι : Type u) (N : ℕ) (u : ι → S)
    (i : ι) : kummerAlgebra S ι N u :=
  Ideal.Quotient.mk _ (MvPolynomial.X i)

/-- `xᵢ^N = uᵢ`。 -/
theorem kummerGen_pow (S : Type u) [CommRing S] (ι : Type u) (N : ℕ) (u : ι → S) (i : ι) :
    (kummerGen S ι N u i) ^ N = algebraMap S (kummerAlgebra S ι N u) (u i) := by
  rw [kummerGen, ← map_pow]
  have hrw : (MvPolynomial.X i : MvPolynomial ι S) ^ N
      = MvPolynomial.C (u i) + ((MvPolynomial.X i) ^ N - MvPolynomial.C (u i)) := by ring
  have hmem : ((MvPolynomial.X i : MvPolynomial ι S) ^ N - MvPolynomial.C (u i))
      ∈ Ideal.span (Set.range fun j => (MvPolynomial.X j : MvPolynomial ι S) ^ N
        - MvPolynomial.C (u j)) := Ideal.subset_span ⟨i, rfl⟩
  rw [hrw, map_add, Ideal.Quotient.eq_zero_iff_mem.mpr hmem, add_zero]
  rfl

/-- `xᵢ` たちが `S`-代数として生成する。 -/
theorem kummerGen_adjoin (S : Type u) [CommRing S] (ι : Type u) (N : ℕ) (u : ι → S) :
    Algebra.adjoin S (Set.range (kummerGen S ι N u)) = ⊤ := by
  have hrange : Set.range (kummerGen S ι N u)
      = (Ideal.Quotient.mkₐ S (Ideal.span (Set.range fun i =>
          (MvPolynomial.X i : MvPolynomial ι S) ^ N - MvPolynomial.C (u i)))) ''
        (Set.range MvPolynomial.X) := by
    ext z
    constructor
    · rintro ⟨i, rfl⟩; exact ⟨MvPolynomial.X i, ⟨i, rfl⟩, rfl⟩
    · rintro ⟨_, ⟨i, rfl⟩, rfl⟩; exact ⟨i, rfl⟩
  rw [hrange, Algebra.adjoin_image, MvPolynomial.adjoin_range_X, Algebra.map_top,
    AlgHom.range_eq_top]
  exact Ideal.Quotient.mkₐ_surjective S _

/-- **★`Rₙ` の微分は `N` で消える**——原文 §3(a) の `Rₙ` そのものに対する
主張。`uᵢ` が単元であることだけを使う。 -/
theorem kummerAlgebra_kaehler_annihilated (S : Type u) [CommRing S] (ι : Type u) (N : ℕ)
    (hN : 0 < N) (u : ι → S) (hu : ∀ i, IsUnit (u i)) :
    ∀ ω : Ω[(kummerAlgebra S ι N u)⁄S], (N : ℕ) • ω = 0 :=
  kaehler_annihilated_of_roots_of_units (kummerGen S ι N u) N hN u
    (kummerGen_pow S ι N u)
    (fun i => (hu i).map (algebraMap S (kummerAlgebra S ι N u)))
    (kaehler_span_of_adjoin _ (kummerGen_adjoin S ι N u))

/-- **塔の遷移写像 `Rₙ → R_{nm}`**——`xᵢ ↦ yᵢ^M`。原文の
`Rₙ ⊂ R_{n+1}`(`M = p` の場合)そのもの。 -/
noncomputable def kummerTransition (S : Type u) [CommRing S] (ι : Type u) (N M : ℕ)
    (u : ι → S) :
    kummerAlgebra S ι N u →ₐ[S] kummerAlgebra S ι (N * M) u := by
  refine Ideal.Quotient.liftₐ _
    (MvPolynomial.aeval (fun i => (kummerGen S ι (N * M) u i) ^ M)) ?_
  intro a ha
  refine Submodule.span_induction ?_ ?_ ?_ ?_ ha
  · rintro _ ⟨i, rfl⟩
    simp only [map_sub, map_pow, MvPolynomial.aeval_X, MvPolynomial.aeval_C]
    rw [← pow_mul, mul_comm M N, kummerGen_pow, sub_self]
  · simp
  · intro x y _ _ hx hy; rw [map_add, hx, hy, add_zero]
  · intro c x _ hx; rw [smul_eq_mul, map_mul, hx, mul_zero]

@[simp] theorem kummerTransition_gen (S : Type u) [CommRing S] (ι : Type u) (N M : ℕ)
    (u : ι → S) (i : ι) :
    kummerTransition S ι N M u (kummerGen S ι N u i) = (kummerGen S ι (N * M) u i) ^ M := by
  rw [kummerTransition, kummerGen, Ideal.Quotient.liftₐ_apply, Ideal.Quotient.lift_mk]
  simp

/-- 非空虚性——`S = ℤ[T]`、`ι = Fin 1`、`N = 5`、`u ≡ 1`
(`S[X]/(X⁵−1)`、1 の 5 乗根を添加した環)。その微分は `5` で消える。 -/
example : ∀ ω : Ω[(kummerAlgebra (Polynomial ℤ) (Fin 1) 5 (fun _ => 1))⁄(Polynomial ℤ)],
    (5 : ℕ) • ω = 0 :=
  kummerAlgebra_kaehler_annihilated (Polynomial ℤ) (Fin 1) 5 (by norm_num) (fun _ => 1)
    (fun _ => isUnit_one)

/-! ## ★almost étale は塔で合成できる(2026-09-05)

Faltings は `A → B → C` の塔を常用する(§3 の `Rₙ ⊂ Sₙ`・`R∞ ⊂ A∞ ⊂ B∞`
等)。その基礎になるのが「almost étale の合成」である。

古典的なエタールの合成は分離冪等元の積 `e_{C/A} = e_{C/B}·ι(e_{B/A})`
で作る。almost 版でも同じ式が効く——鍵は

    `C ⊗_A C → C ⊗_B C` の核(= `B` の元による相対対角イデアル
     `relDiag`)が `ι(w_B)` で消える

ことで、これは `(1⊗b − b⊗1)·w_B = 0` を `ι` で送るだけである。
水準は `c_B·c_C` になる(塔があるので回収できる)。 -/

/-- **相対対角イデアル**——`C ⊗_A C → C ⊗_B C` の核を生成する元たち。 -/
def relDiag (A B C : Type u) [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra B C] [Algebra A C] [IsScalarTower A B C] :
    Ideal (TensorProduct A C C) :=
  Ideal.span {x | ∃ b : B, x = (1 : C) ⊗ₜ[A] (algebraMap B C b)
    - (algebraMap B C b) ⊗ₜ[A] (1 : C)}

/-- **相対対角イデアルは `ι(w_B)` で消える**——`(1⊗b − b⊗1)·w_B = 0` を
`ι` で送るだけ。 -/
theorem relDiag_mul_eq_zero (A B C : Type u) [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra B C] [Algebra A C] [IsScalarTower A B C]
    (wB : TensorProduct A B B)
    (hB_ann : ∀ b : B, ((1 : B) ⊗ₜ[A] b - b ⊗ₜ[A] (1 : B)) * wB = 0)
    (x : TensorProduct A C C) (hx : x ∈ relDiag A B C) :
    x * (Algebra.TensorProduct.map (IsScalarTower.toAlgHom A B C)
      (IsScalarTower.toAlgHom A B C) wB) = 0 := by
  refine Submodule.span_induction ?_ ?_ ?_ ?_ hx
  · rintro _ ⟨b, rfl⟩
    have h : ((1 : C) ⊗ₜ[A] (algebraMap B C b) - (algebraMap B C b) ⊗ₜ[A] (1 : C))
        = Algebra.TensorProduct.map (IsScalarTower.toAlgHom A B C)
            (IsScalarTower.toAlgHom A B C) ((1 : B) ⊗ₜ[A] b - b ⊗ₜ[A] (1 : B)) := by
      rw [map_sub, Algebra.TensorProduct.map_tmul, Algebra.TensorProduct.map_tmul]
      simp
    rw [h, ← map_mul, hB_ann b, map_zero]
  · rw [zero_mul]
  · intro y z _ _ hy hz; rw [add_mul, hy, hz, add_zero]
  · intro a y _ hy; rw [smul_mul_assoc, hy, smul_zero]

/-- **★almost étale は塔で合成できる**(witness 版)——
`w := w_C · ι(w_B)` が `C/A` の witness(水準 `c_B·c_C`)。

`hC_ann` は「`w_C` が `C/B` の witness を `C ⊗_A C` へ持ち上げたもの」
という条件で、`C ⊗_A C → C ⊗_B C` の核の記述を使わずに済む形にしてある。 -/
theorem witness_comp (A B C : Type u) [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra B C] [Algebra A C] [IsScalarTower A B C]
    (cB cC : A) (wB : TensorProduct A B B) (wC : TensorProduct A C C)
    (hB_ann : ∀ b : B, ((1 : B) ⊗ₜ[A] b - b ⊗ₜ[A] (1 : B)) * wB = 0)
    (hB_aug : Algebra.TensorProduct.lmul' A wB = cB • (1 : B))
    (hC_ann : ∀ c : C, ((1 : C) ⊗ₜ[A] c - c ⊗ₜ[A] (1 : C)) * wC ∈ relDiag A B C)
    (hC_aug : Algebra.TensorProduct.lmul' A wC = cC • (1 : C)) :
    (∀ c : C, ((1 : C) ⊗ₜ[A] c - c ⊗ₜ[A] (1 : C)) *
        (wC * Algebra.TensorProduct.map (IsScalarTower.toAlgHom A B C)
          (IsScalarTower.toAlgHom A B C) wB) = 0)
    ∧ Algebra.TensorProduct.lmul' A
        (wC * Algebra.TensorProduct.map (IsScalarTower.toAlgHom A B C)
          (IsScalarTower.toAlgHom A B C) wB) = (cB * cC) • (1 : C) := by
  constructor
  · intro c
    rw [← mul_assoc]
    exact relDiag_mul_eq_zero A B C wB hB_ann _ (hC_ann c)
  · rw [map_mul, hC_aug, lmul'_map, hB_aug, map_smul, map_one,
      Algebra.smul_def, Algebra.smul_def, Algebra.smul_def, map_mul]
    ring

/-- **witness の族も塔で合成できる**——水準は `ϖ(k+1)² ∣ ϖ k` で回収する。 -/
theorem hasAlmostWitnesses_comp {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra B C] [Algebra A C] [IsScalarTower A B C]
    {q : ℕ} (hq : 2 ≤ q) (T : PDivTower A q)
    (hB : T.HasAlmostWitnesses B)
    (hC : ∀ k : ℕ, ∃ wC : TensorProduct A C C,
      (∀ c : C, ((1 : C) ⊗ₜ[A] c - c ⊗ₜ[A] (1 : C)) * wC ∈ relDiag A B C) ∧
      Algebra.TensorProduct.lmul' A wC = T.ϖ k • (1 : C)) :
    T.HasAlmostWitnesses C := by
  intro k
  obtain ⟨wB, hBann, hBaug⟩ := hB (k + 1)
  obtain ⟨wC, hCann, hCaug⟩ := hC (k + 1)
  have hdvd : (T.ϖ (k + 1) * T.ϖ (k + 1)) ∣ T.ϖ k := by
    refine ⟨(T.ϖ (k + 1)) ^ (q - 2), ?_⟩
    rw [← T.ϖ_succ k, ← pow_two, ← pow_add]
    congr 1
    omega
  obtain ⟨e, he⟩ := hdvd
  obtain ⟨h1, h2⟩ := witness_comp A B C (T.ϖ (k + 1)) (T.ϖ (k + 1)) wB wC
    hBann hBaug hCann hCaug
  refine ⟨e • (wC * Algebra.TensorProduct.map (IsScalarTower.toAlgHom A B C)
    (IsScalarTower.toAlgHom A B C) wB), ?_, ?_⟩
  · intro c
    rw [mul_smul_comm, h1 c, smul_zero]
  · rw [map_smul, h2, smul_smul, he, mul_comm]

/-- `B = A` のときは相対対角イデアルは `⊥`——合成則がもとの witness を
そのまま返すことの対照(非空虚性)。 -/
theorem relDiag_self (A C : Type u) [CommRing A] [CommRing C] [Algebra A C] :
    relDiag A A C = ⊥ := by
  rw [relDiag, Ideal.span_eq_bot]
  rintro _ ⟨a, rfl⟩
  have h : ((1 : C) ⊗ₜ[A] (algebraMap A C a)) = ((algebraMap A C a) ⊗ₜ[A] (1 : C)) := by
    rw [Algebra.algebraMap_eq_smul_one, TensorProduct.tmul_smul, TensorProduct.smul_tmul']
  rw [h, sub_self]

/-! ## 塔の言葉へ——「`ε` は任意に小さくできる」

`PDivTower.HasAlmostWitnesses` は「**各 `k` について** 水準 `ϖ k` の
witness がある」という述語で、これが原文の *"the assertion now follows
as `ε` can be arbitrarily small"* に対応する(`ε = 1/q^k`)。

上の 2 つの操作(全射での押し出し・底変換)はどちらも**水準を変えない**
ので、そのまま塔の言葉へ持ち上がる。 -/

namespace PDivTower

/-- 塔は環準同型で押し出せる。`ϖ' k := f (ϖ k)`。 -/
def mapRingHom {A A' : Type u} [CommRing A] [CommRing A'] {q : ℕ}
    (T : PDivTower A q) (f : A →+* A') : PDivTower A' q where
  ϖ k := f (T.ϖ k)
  ϖ_succ k := by rw [← map_pow, T.ϖ_succ]

@[simp] theorem mapRingHom_ϖ {A A' : Type u} [CommRing A] [CommRing A'] {q : ℕ}
    (T : PDivTower A q) (f : A →+* A') (k : ℕ) :
    (T.mapRingHom f).ϖ k = f (T.ϖ k) := rfl

end PDivTower

/-- **witness の族は底変換で保たれる**——`A'` 上の塔を
`ϖ' k = algebraMap A A' (ϖ k)` と取れば、`A' ⊗_A B` は `A'` 上
同じ水準の witness の族を持つ。

原文 *"`B∞` is almost isomorphic to `S∞ ⊗_{R∞} A∞` **as this is almost
étale over `S∞`**"* の括弧内を、塔(= 任意に小さい `ε`)の言葉で
述べたもの。 -/
theorem hasAlmostWitnesses_baseChange {A A' B : Type u} [CommRing A] [CommRing A']
    [CommRing B] [Algebra A A'] [Algebra A B] {q : ℕ} (T : PDivTower A q)
    (hwit : T.HasAlmostWitnesses B) :
    (T.mapRingHom (algebraMap A A')).HasAlmostWitnesses (TensorProduct A A' B) := by
  intro k
  obtain ⟨w, hann, haug⟩ := hwit k
  obtain ⟨h1, h2⟩ := witness_baseChange A A' B (T.ϖ k) w hann haug
  exact ⟨bcWitnessMap A A' B w, h1, h2⟩

/-- **witness の族は全射 `A`-代数写像で押し出せる**(`witness_pushforward`
の塔版)。 -/
theorem hasAlmostWitnesses_of_surjective {A B B' : Type u} [CommRing A] [CommRing B]
    [CommRing B'] [Algebra A B] [Algebra A B'] {q : ℕ} (T : PDivTower A q)
    (f : B →ₐ[A] B') (hf : Function.Surjective f)
    (hwit : T.HasAlmostWitnesses B) : T.HasAlmostWitnesses B' := by
  intro k
  obtain ⟨w, hann, haug⟩ := hwit k
  obtain ⟨h1, h2⟩ := witness_pushforward f hf (T.ϖ k) w hann haug
  exact ⟨Algebra.TensorProduct.map f f w, h1, h2⟩

/-- **分裂エタール代数は任意の塔について witness の族を持つ**——
`pi_witness`(水準 `1`)を `ϖ k` 倍するだけ。

★`PDivTower.HasAlmostWitnesses` が空虚に真な述語ではないことの対照。
これで `kaehler_isAlmostZero`・`hochschild_H2_isAlmostCoboundary` など
この述語を仮定する結果すべてに具体的な入力が付いた。 -/
theorem hasAlmostWitnesses_pi {A : Type u} [CommRing A] {q : ℕ} (T : PDivTower A q)
    {ι : Type u} [Fintype ι] [DecidableEq ι] : T.HasAlmostWitnesses (ι → A) := by
  intro k
  obtain ⟨hann, haug⟩ := (pi_witness (A := A) (ι := ι))
  refine ⟨T.ϖ k • (∑ i : ι, (Pi.single i (1 : A)) ⊗ₜ[A] (Pi.single i (1 : A))), ?_, ?_⟩
  · intro b
    rw [mul_smul_comm, hann b, smul_zero]
  · rw [map_smul, haug, one_smul]

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
