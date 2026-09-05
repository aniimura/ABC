import ABC3.Found.Falt1.AlmostDerivation
import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Algebra.Subalgebra.Directed

/-!
# [Falt1] `Theorem 2.3` への第 1 段——almost 冪等元の持ち上げ(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §2、
`Theorem 2.3`(物理 p.7-8 = 印字 p.260-261)。

> *2.3. Theorem. Suppose `I ⊂ A` is a nilpotent ideal, `B` an almost étale
> covering of `Ā = A/I`, `B = Ā + mB`. Then there exists an almost étale
> covering `B̃` of `A` such that `B̃ ⊗_A Ā ≅ B` modulo `p`-torsion.*

原文の証明は 5 段からなる:

1. `I² = 0` としてよい。
2. **`B` は自由加群 `Ā^r` の直和因子「`p^ε` を除いて」**——`ē² = p^ε ē` を
   満たす `r×r` 行列 `ē` があって `B = Ā^r/(p^ε − ē)(Ā^r)`(`p`-捩れを除く)。
3. **"Tripling `ε`"**——`ē` は `A` 上の行列 `e` に `e² = p^ε e` を満たす形で
   持ち上がる。★**本ファイルはこの段を閉じる。**
4. `B_ε := A^r/(p^ε − e)(A^r)` は `p^ε` を除いて `A`-射影的、
   乗法 `m` は `m_ε : B_ε ⊗_A B_ε → B_ε` に持ち上がる。
5. 結合律の障害は `H³(B/Ā)` の類。`ε` を細かくした `C_σ := A + p^σ B_ε` の
   増大列の合併が求める `B̃`。

## 本ファイルの内容(第 3 段)

鍵は **可換環 `A[X]` の 1 つの恒等式**である:

    (3sX² − 2X³)² − s³·(3sX² − 2X³) = (4X² − 4sX − 3s²)·(X² − sX)²

つまり `e := 3s·u² − 2u³`(古典的な冪等元持ち上げ `3u² − 2u³` に `s` を
掛けたもの)は、「冪等からのずれ」`n := u² − s·u` の**平方**を法として
`e² = s³·e` を満たす。`I² = 0` なら `n² = 0` なので、これは厳密な等式になる
——`s = p^ε` と取れば、これがまさに Faltings の *"Tripling `ε`"* である
(指数が `ε` から `3ε` に増える)。

`u` は行列環(非可換)の元なので、恒等式は `Polynomial.eval₂RingHom'`
(係数と可換な元での評価)を通して移す。
-/

namespace ABC3.Found.Falt1

universe u

open Polynomial in
/-- `e := 3sX² − 2X³` について
`e² − s³e = (4X²−4sX−3s²)·(X²−sX)²`。可換環 `A[X]` での恒等式。 -/
theorem almost_idem_poly {A : Type u} [CommRing A] (s : A) :
    (3 * C s * X^2 - 2 * X^3)^2 - C s^3 * (3 * C s * X^2 - 2 * X^3)
      = (4 * X^2 - 4 * C s * X - 3 * C s^2) * (X^2 - C s * X)^2 := by
  ring

open Polynomial in
/-- **almost 冪等元の持ち上げ(核心の計算)**。`u` の「冪等からのずれ」
`n := u² − s·u` が `n² = 0` を満たすなら、`e := 3s·u² − 2u³` は
`e² = s³·e` を**厳密に**満たす。`S` は非可換でよい(行列環)。 -/
theorem almost_idem_lift {A S : Type u} [CommRing A] [Ring S] [Algebra A S]
    (s : A) (u : S) (hn : (u^2 - algebraMap A S s * u)^2 = 0) :
    (3 * algebraMap A S s * u^2 - 2 * u^3)^2
      = algebraMap A S (s^3) * (3 * algebraMap A S s * u^2 - 2 * u^3) := by
  set φ : Polynomial A →+* S :=
    Polynomial.eval₂RingHom' (algebraMap A S) u
      (fun r => Algebra.commute_algebraMap_left r u) with hφ
  have h := congrArg φ (almost_idem_poly (A := A) s)
  simp only [hφ, Polynomial.eval₂RingHom'_apply, map_sub, map_mul, map_pow, map_ofNat,
    Polynomial.eval₂_C, Polynomial.eval₂_X] at h
  rw [hn, mul_zero, sub_eq_zero] at h
  rw [h, map_pow]

/-- **almost 冪等元は平方零イデアルに沿って持ち上がる**——Faltings の
*"Tripling `ε`"*(物理 p.7 = 印字 p.260)。

`ē² = s·ē` を満たす `ē` は、`e² = s³·e` を満たす `e` に持ち上がり、
`e` の像はちょうど `s²·ē` になる。`s = p^ε` と取れば「`ε` を 3 倍すれば
`ē` は持ち上がる」という原文の主張そのもの。 -/
theorem exists_almost_idem_lift {A S S' : Type u} [CommRing A] [Ring S] [Ring S']
    [Algebra A S] [Algebra A S'] (π : S →ₐ[A] S') (hπ : Function.Surjective π)
    (hker : ∀ x y : S, π x = 0 → π y = 0 → x * y = 0)
    (s : A) (ebar : S') (hebar : ebar^2 = algebraMap A S' s * ebar) :
    ∃ e : S, e^2 = algebraMap A S (s^3) * e
      ∧ π e = algebraMap A S' (s^2) * ebar := by
  obtain ⟨u, hu⟩ := hπ ebar
  have hn0 : π (u^2 - algebraMap A S s * u) = 0 := by
    rw [map_sub, map_pow, map_mul, AlgHom.commutes, hu, hebar, sub_self]
  have hn : (u^2 - algebraMap A S s * u)^2 = 0 := by
    rw [sq]; exact hker _ _ hn0 hn0
  refine ⟨3 * algebraMap A S s * u^2 - 2 * u^3, almost_idem_lift s u hn, ?_⟩
  have h3 : ebar^3 = (algebraMap A S' s)^2 * ebar := by
    rw [pow_succ, hebar, mul_assoc, ← sq, hebar, ← mul_assoc, ← sq]
  simp only [map_sub, map_mul, map_pow, map_ofNat, AlgHom.commutes, hu]
  rw [hebar, h3]
  simp only [sq, mul_assoc]
  rw [← sub_mul, show (3:S') - 2 = 1 by norm_num, one_mul]

/-- 非空虚性——`ℤ/4 ↠ ℤ/2`(核 `(2)` は平方零)。 -/
example : ∃ e : ZMod 4, e^2 = algebraMap ℤ (ZMod 4) ((1:ℤ)^3) * e
    ∧ (ZMod.castHom (by norm_num) (ZMod 2)).toIntAlgHom e
      = algebraMap ℤ (ZMod 2) ((1:ℤ)^2) * (1 : ZMod 2) := by
  refine exists_almost_idem_lift ((ZMod.castHom (by norm_num) (ZMod 2)).toIntAlgHom)
    ?_ ?_ 1 1 (by decide)
  · intro x
    exact ⟨(ZMod.val x : ZMod 4), by revert x; decide⟩
  · decide

/-! ## 第 2 段——almost 射影性はそのまま almost 冪等行列を与える

原文の *"First `B` is a direct summand in a free module `Ā^r` up to some
`p^ε`. That is, there exists an `r×r` matrix `ē` such that `ē² = p^ε ē` and
`B = Ā^r/(p^ε − ē)(Ā^r)` modulo `p`-torsion."*

前半(`ē` の存在)は `AlmostProjective.almost_projective_factor` から
**1 行で**出る:`g∘f = c·id` なら `E := f∘g` について
`E² = f g f g = f (c·id) g = c·E`。 -/

/-- **`Theorem 2.3` の第 2 段**——`g∘f = c·id` なら `E := f∘g` は
`E² = c·E` を満たす。 -/
theorem almost_idem_of_factor {A B : Type u} [CommRing A] [AddCommGroup B] [Module A B]
    {ι : Type u} [Fintype ι] (c : A) (f : B →ₗ[A] (ι → A)) (g : (ι → A) →ₗ[A] B)
    (hgf : ∀ b : B, g (f b) = c • b) :
    (f ∘ₗ g) ∘ₗ (f ∘ₗ g) = c • (f ∘ₗ g) := by
  ext x
  simp only [LinearMap.coe_comp, Function.comp_apply, LinearMap.smul_apply]
  rw [hgf, map_smul]

/-- 第 2 段を行列の言葉へ。 -/
theorem almost_idem_matrix_of_factor {A B : Type u} [CommRing A] [AddCommGroup B] [Module A B]
    {ι : Type u} [Fintype ι] [DecidableEq ι] (c : A)
    (f : B →ₗ[A] (ι → A)) (g : (ι → A) →ₗ[A] B)
    (hgf : ∀ b : B, g (f b) = c • b) :
    (LinearMap.toMatrixAlgEquiv' (f ∘ₗ g))^2
      = algebraMap A (Matrix ι ι A) c * LinearMap.toMatrixAlgEquiv' (f ∘ₗ g) := by
  have h2 : (f ∘ₗ g) * (f ∘ₗ g) = c • (f ∘ₗ g) := almost_idem_of_factor c f g hgf
  have h3 := congrArg (LinearMap.toMatrixAlgEquiv' (R := A) (n := ι)) h2
  rw [map_mul, map_smul] at h3
  rw [sq, h3, Algebra.smul_def]

/-- 平方零イデアル `I` に沿った行列環の還元は、積を零化する。 -/
theorem matrix_sq_zero {ι A : Type u} [Fintype ι] [DecidableEq ι] [CommRing A]
    (I : Ideal A) (hI : ∀ u ∈ I, ∀ v ∈ I, u * v = 0) (x y : Matrix ι ι A)
    (hx : (Ideal.Quotient.mkₐ A I).mapMatrix x = 0)
    (hy : (Ideal.Quotient.mkₐ A I).mapMatrix y = 0) : x * y = 0 := by
  ext i j
  rw [Matrix.mul_apply, Matrix.zero_apply]
  refine Finset.sum_eq_zero (fun k _ => ?_)
  refine hI _ ?_ _ ?_
  · rw [← Ideal.Quotient.eq_zero_iff_mem]
    exact congrFun (congrFun (congrArg (fun m => (m : Matrix ι ι (A ⧸ I))) hx) i) k
  · rw [← Ideal.Quotient.eq_zero_iff_mem]
    exact congrFun (congrFun (congrArg (fun m => (m : Matrix ι ι (A ⧸ I))) hy) k) j

/-- **`Theorem 2.3` の第 3 段(行列版)**——`Ā = A/I`(`I² = 0`)上の
almost 冪等行列 `ē`(`ē² = s·ē`)は `A` 上の行列 `e`(`e² = s³·e`)に
持ち上がり、その像はちょうど `s²·ē`。

第 2 段(`almost_idem_matrix_of_factor`)が `Ā` 上の `ē` を供給するので、
この 2 つを繋げば原文の *"Tripling `ε` we may assume that `ē` lifts to an
`r×r` matrix `e` with `e² = p^ε e`"* がそのまま出る。 -/
theorem exists_matrix_almost_idem_lift {ι A : Type u} [Fintype ι] [DecidableEq ι] [CommRing A]
    (I : Ideal A) (hI : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (s : A) (ebar : Matrix ι ι (A ⧸ I))
    (hebar : ebar^2 = algebraMap A (Matrix ι ι (A ⧸ I)) s * ebar) :
    ∃ e : Matrix ι ι A,
      e^2 = algebraMap A (Matrix ι ι A) (s^3) * e ∧
      (Ideal.Quotient.mkₐ A I).mapMatrix e
        = algebraMap A (Matrix ι ι (A ⧸ I)) (s^2) * ebar := by
  refine exists_almost_idem_lift (Ideal.Quotient.mkₐ A I).mapMatrix ?_
    (matrix_sq_zero I hI) s ebar hebar
  intro m
  refine ⟨Matrix.of (fun i j => (Ideal.Quotient.mk_surjective (m i j)).choose), ?_⟩
  ext i j
  exact (Ideal.Quotient.mk_surjective (m i j)).choose_spec

/-! ## 第 4 段(前半)——`B_ε` の almost 射影性

原文の *"Let `B_ε = A^r/(p^ε − e)(A^r)`. Then `B_ε` is `A`-projective up to
`p^ε`, that is, `p^ε` annihilates all higher Ext-groups `Ext^i(B_ε, M)`."*

`E² = s·E` から `E·(s·1 − E) = s·E − E² = 0` なので `E` は商
`N := (ι→A)/im(s·1−E)` を経由する。その `α : N → (ι→A)` と商写像
`β : (ι→A) → N` について `β∘α = s·id_N`——**`s·id_N` が自由加群を経由して
分解する**。これがそのまま「`s` を除いて射影的」の内容であり、
`AlmostProjective.almost_lift_of_surjective` にそのまま渡せる
(Ext の言葉に直すなら `ext_smul_eq_zero_of_almost_split`——
`tools/lean-idioms.md` #54 の「almost split ⟹ almost vanishing」)。 -/

/-- `E² = s·E` から作る商加群 `N := (ι→A)/im(s·1 − E)` について、
`s·id_N` が自由加群 `ι→A` を経由して分解する。 -/
theorem almostIdemQuotient_factor {A : Type u} [CommRing A] {ι : Type u} [Fintype ι]
    (s : A) (E : Module.End A (ι → A)) (hE : E * E = s • E) :
    ∃ (α : ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E)) →ₗ[A] (ι → A))
      (β : (ι → A) →ₗ[A] ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E))),
      ∀ x, β (α x) = s • x := by
  set P := LinearMap.range (s • (1 : Module.End A (ι → A)) - E) with hP
  have hEker : P ≤ LinearMap.ker E := by
    rintro z ⟨y, rfl⟩
    have hz : E * (s • (1 : Module.End A (ι → A)) - E) = 0 := by
      rw [mul_sub, mul_smul_comm, mul_one, hE, sub_self]
    show E _ = 0
    exact congrArg (fun (f : Module.End A (ι → A)) => f y) hz
  refine ⟨P.liftQ E hEker, P.mkQ, ?_⟩
  intro x
  induction x using Submodule.Quotient.induction_on with
  | H y =>
    show P.mkQ ((P.liftQ E hEker) (Submodule.Quotient.mk y)) = s • Submodule.Quotient.mk y
    rw [P.liftQ_apply, ← Submodule.mkQ_apply, ← map_smul, ← sub_eq_zero, ← map_sub,
      Submodule.mkQ_apply, Submodule.Quotient.mk_eq_zero]
    refine ⟨-y, ?_⟩
    simp only [LinearMap.sub_apply, LinearMap.smul_apply, Module.End.one_apply, map_neg]
    abel

/-- **`Theorem 2.3` の第 4 段(前半)**——`B_ε := (ι→A)/im(s·1−E)` は
`s` を除いて射影的:任意の全射に沿って `s·φ` が持ち上がる。 -/
theorem almostIdemQuotient_almost_lift {A C D : Type u} [CommRing A]
    [AddCommGroup C] [Module A C] [AddCommGroup D] [Module A D]
    {ι : Type u} [Fintype ι]
    (s : A) (E : Module.End A (ι → A)) (hE : E * E = s • E)
    (π : C →ₗ[A] D) (hπ : Function.Surjective π)
    (φ : ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E)) →ₗ[A] D) :
    ∃ ψ : ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E)) →ₗ[A] C,
      ∀ x, π (ψ x) = s • φ x := by
  obtain ⟨α, β, hαβ⟩ := almostIdemQuotient_factor s E hE
  exact almost_lift_of_surjective s ⟨ι, inferInstance, α, β, hαβ⟩ π hπ φ

/-! ## 第 4 段(後半への足場)——テンソル積も almost 射影的

原文の *"Let `m : B ⊗_Ā B → B` denote the multiplication map. `p^{3ε}·m`
lifts to `m_ε : B_ε ⊗_A B_ε → B_ε`."*

`B_ε` の almost 射影性(上)から `B_ε ⊗_A B_ε` の almost 射影性が出る
(水準は積 `s·s`)。ここで中間加群は `(ι→A) ⊗_A (ι→A)`——有限自由では
あるが `μ → A` の形には**書いていない**ので、`almost_lift_of_surjective`
を「中間加群は射影的でありさえすればよい」形へ一般化しておく。 -/

/-- `almost_lift_of_surjective` の一般化——中間加群は有限自由でなくとも
**射影的**でありさえすればよい。 -/
theorem almost_lift_of_projective_factor {A B C D P : Type u} [CommRing A]
    [AddCommGroup B] [Module A B] [AddCommGroup C] [Module A C]
    [AddCommGroup D] [Module A D] [AddCommGroup P] [Module A P] [Module.Projective A P]
    (c : A) (f : B →ₗ[A] P) (g : P →ₗ[A] B) (hgf : ∀ b : B, g (f b) = c • b)
    (π : C →ₗ[A] D) (hπ : Function.Surjective π) (φ : B →ₗ[A] D) :
    ∃ ψ : B →ₗ[A] C, ∀ b : B, π (ψ b) = c • φ b := by
  obtain ⟨χ, hχ⟩ := Module.projective_lifting_property π (φ ∘ₗ g) hπ
  refine ⟨χ ∘ₗ f, fun b => ?_⟩
  have h := LinearMap.congr_fun hχ (f b)
  simp only [LinearMap.coe_comp, Function.comp_apply] at h ⊢
  rw [h, hgf b, map_smul]

/-- **almost 射影加群のテンソル積は almost 射影的**——水準は積になる。 -/
theorem almost_projective_tensor {A M N P Q : Type u} [CommRing A]
    [AddCommGroup M] [Module A M] [AddCommGroup N] [Module A N]
    [AddCommGroup P] [Module A P] [AddCommGroup Q] [Module A Q]
    (c d : A)
    (f : M →ₗ[A] P) (g : P →ₗ[A] M) (hgf : ∀ x : M, g (f x) = c • x)
    (f' : N →ₗ[A] Q) (g' : Q →ₗ[A] N) (hgf' : ∀ x : N, g' (f' x) = d • x)
    (z : TensorProduct A M N) :
    TensorProduct.map g g' (TensorProduct.map f f' z) = (c * d) • z := by
  induction z using TensorProduct.induction_on with
  | zero => simp
  | tmul x y =>
    rw [TensorProduct.map_tmul, TensorProduct.map_tmul, hgf, hgf',
      TensorProduct.tmul_smul, TensorProduct.smul_tmul', TensorProduct.smul_tmul',
      smul_smul, mul_comm d c]
  | add u v ihu ihv => rw [map_add, map_add, ihu, ihv, smul_add]

/-- **`B_ε ⊗_A B_ε` も almost 射影的**(水準 `s²`)——乗法 `m` は
`s²` を除いて `B_ε ⊗_A B_ε → B_ε` に持ち上がる。 -/
theorem almostIdemQuotient_tensor_almost_lift {A C D : Type u} [CommRing A]
    [AddCommGroup C] [Module A C] [AddCommGroup D] [Module A D]
    {ι : Type u} [Fintype ι]
    (s : A) (E : Module.End A (ι → A)) (hE : E * E = s • E)
    (π : C →ₗ[A] D) (hπ : Function.Surjective π)
    (φ : TensorProduct A
        ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E))
        ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E)) →ₗ[A] D) :
    ∃ ψ : TensorProduct A
        ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E))
        ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E)) →ₗ[A] C,
      ∀ x, π (ψ x) = (s * s) • φ x := by
  obtain ⟨α, β, hαβ⟩ := almostIdemQuotient_factor s E hE
  exact almost_lift_of_projective_factor (P := TensorProduct A (ι → A) (ι → A)) (s * s)
    (TensorProduct.map α α) (TensorProduct.map β β)
    (fun z => almost_projective_tensor s s α β hαβ α β hαβ z) π hπ φ

/-! ## 第 5 段——結合律の障害は "five-term identity" でコサイクルになる

原文の *"The obstruction to `m_ε` being associative is given by
`b₀[m_ε(b₁, m_ε(b₂,b₃)) − m_ε(m_ε(b₁,b₂),b₃)]b₄ ∈ Kern(B_ε → B)`
and defines a class in `H³(B/Ā)` with values in this module.
**(That we get a cycle amounts to a five-term identity well known in the
study of associativity.)**"*

その "five-term identity" は、**任意の**双線形な「乗法」`m` について
恒等的に成り立つ:結合子 `F(a,b,c) := m a (m b c) − m (m a b) c` は

    m w (F x y z) − F (m w x) y z + F w (m x y) z − F w x (m y z) + m (F w x y) z = 0

を満たす(展開すると 10 項がちょうど打ち消し合う)。
`H²` のときの `HochschildLowDegree.obstruction_identity` に対応する段で、
これと「核の上では `m` が加群作用に一致する」ことを合わせれば
Hochschild 3-コサイクル条件になり、`hochschild_H3_almost_coboundary` が
効く。 -/

/-- **結合律の "five-term identity"**。任意の双線形 `m` について恒等的に
成り立つ——「障害がコサイクルである」ことの代数的な核。 -/
theorem five_term_identity {A N : Type u} [CommRing A] [AddCommGroup N] [Module A N]
    (m : N →ₗ[A] N →ₗ[A] N) (w x y z : N) :
    m w (m x (m y z) - m (m x y) z)
      - (m (m w x) (m y z) - m (m (m w x) y) z)
      + (m w (m (m x y) z) - m (m w (m x y)) z)
      - (m w (m x (m y z)) - m (m w x) (m y z))
      + m (m w (m x y) - m (m w x) y) z = 0 := by
  simp only [map_sub, LinearMap.sub_apply]
  abel

/-- 結合子を三重線形写像として束ねたもの
(`hochschild_H3_almost_coboundary` にそのまま渡せる形)。 -/
noncomputable def associatorTri {A N : Type u} [CommRing A] [AddCommGroup N] [Module A N]
    (m : N →ₗ[A] N →ₗ[A] N) : N →ₗ[A] N →ₗ[A] (N →ₗ[A] N) :=
  LinearMap.mk₂ A (fun a b => (m a) ∘ₗ (m b) - m (m a b))
    (by intro a a' b; ext c; simp; abel)
    (by intro r a b; ext c; simp [smul_sub])
    (by intro a b b'; ext c; simp; abel)
    (by intro r a b; ext c; simp [smul_sub])

@[simp] theorem associatorTri_apply {A N : Type u} [CommRing A] [AddCommGroup N] [Module A N]
    (m : N →ₗ[A] N →ₗ[A] N) (a b c : N) :
    associatorTri m a b c = m a (m b c) - m (m a b) c := rfl

/-! ## 第 5 段の後半——**結合律への補正は倍化なしで厳密に閉じる**

原文の *"It follows that we can modify `p^{4ε}m_ε` such that the
associative law holds."*

`t·(結合子) = δh` が `hochschild_H3_almost_coboundary` から得られたとき、
補正 `ν := t·μ − h` の結合子を展開すると

    F' = t²·F + t·δh + [h(x,h(y,z)) − h(h(x,y),z)]
       = t²·F − t²·F + [二次の項]

となり、**二次の項は `h` の値が「核」`J` に入り、かつ `h` が `J` の引数で
消えることから `0`** になる。すなわち `ν` は**厳密に**結合的である
——`Theorem 2.2` の第 2 段で必要だった「`ε` を倍にする」
(`doubling_multiplicative`)がここでは**要らない**。

「`h` が `J` の引数で消える」のは、`h(x,y) = Σᵢ xᵢ·F(yᵢ,x,y)` の形をして
いて、**結合子 `F` 自体が `J` の引数で消える**からである
(`associator_eq_zero_of_act`・`associator_eq_zero_of_act'`)——
`J` が `π : N → B` を通して `B`-加群であること(`μ x (ι j) = ι (π x • j)`)
と `B` の可換性・結合性だけから出る。 -/

/-- **結合律への補正**。`t·(結合子) = δh` で `h` が二次で消えるなら、
`ν := t·μ − h` は**厳密に**結合的。 -/
theorem associativity_correction {A N : Type u} [CommRing A] [AddCommGroup N] [Module A N]
    (μ : N →ₗ[A] N →ₗ[A] N) (t : A) (h : N →ₗ[A] N →ₗ[A] N)
    (hv2 : ∀ x y z : N, h x (h y z) = 0)
    (hv1 : ∀ x y z : N, h (h x y) z = 0)
    (hcob : ∀ x y z : N, t • (μ x (μ y z) - μ (μ x y) z)
      = μ x (h y z) - h (μ x y) z + h x (μ y z) - μ (h x y) z)
    (x y z : N) :
    (t • μ - h) x ((t • μ - h) y z) = (t • μ - h) ((t • μ - h) x y) z := by
  have hL : (t • μ - h) x ((t • μ - h) y z)
      = t • (t • μ x (μ y z)) - t • μ x (h y z) - t • h x (μ y z) := by
    show (t • μ x - h x) (t • μ y z - h y z) = _
    simp only [LinearMap.sub_apply, LinearMap.smul_apply, map_sub, map_smul]
    rw [hv2]
    module
  have hR : (t • μ - h) ((t • μ - h) x y) z
      = t • (t • μ (μ x y) z) - t • μ (h x y) z - t • h (μ x y) z := by
    show (t • μ (t • μ x y - h x y) - h (t • μ x y - h x y)) z = _
    simp only [map_sub, map_smul, LinearMap.sub_apply, LinearMap.smul_apply]
    rw [hv1]
    module
  rw [hL, hR]
  have hc := hcob x y z
  simp only [smul_sub] at hc ⊢
  linear_combination (norm := module) t • hc

/-- **核の引数では結合子が消える**(第 3 引数)。`J` が `π` を通して
`B`-加群であり `μ x (ι j) = ι (π x • j)` なら、結合子は `0`。 -/
theorem associator_eq_zero_of_act {A N B J : Type u} [CommRing A] [AddCommGroup N] [Module A N]
    [CommRing B] [Algebra A B] [AddCommGroup J] [Module A J] [Module B J] [IsScalarTower A B J]
    (μ : N →ₗ[A] N →ₗ[A] N) (ι : J →ₗ[A] N) (π : N →ₗ[A] B)
    (hact : ∀ (x : N) (j : J), μ x (ι j) = ι (π x • j))
    (hmulπ : ∀ x y : N, π (μ x y) = π x * π y)
    (x y : N) (j : J) :
    μ x (μ y (ι j)) - μ (μ x y) (ι j) = 0 := by
  rw [hact, hact, hact, hmulπ, mul_smul, sub_self]

/-- 第 2 引数が `J` のときも同様(`μ` の可換性を使う)。 -/
theorem associator_eq_zero_of_act' {A N B J : Type u} [CommRing A] [AddCommGroup N] [Module A N]
    [CommRing B] [Algebra A B] [AddCommGroup J] [Module A J] [Module B J] [IsScalarTower A B J]
    (μ : N →ₗ[A] N →ₗ[A] N) (ι : J →ₗ[A] N) (π : N →ₗ[A] B)
    (hcomm : ∀ x y : N, μ x y = μ y x)
    (hact : ∀ (x : N) (j : J), μ x (ι j) = ι (π x • j))
    (x z : N) (j : J) :
    μ x (μ (ι j) z) - μ (μ x (ι j)) z = 0 := by
  rw [hcomm (ι j) z, hact z j, hact x (π z • j), hact x j,
    hcomm (ι (π x • j)) z, hact z (π x • j), smul_comm, sub_self]

/-! ## 第 4 段(後半)——`B_ε ⊗_A Ā` と `B` の差は `p^ε` で消える

原文の *"`B_ε` contains `B ⊗_Ā A/(p-torsion)`, the quotient being
annihilated by `p^ε`."*

"Tripling ε" の帰結 `π e = s²·ē` から

    π(s³·1 − e) = s³·1 − s²·ē = s²·(s·1 − ē)

なので、`B_ε ⊗_A Ā = coker(s²·(s·1 − ē))` と `B = coker(s·1 − ē)` を
比べる自然な全射ができ、**その核は `s²` で零化される**——`coker(d·f)`
から `coker(f)` への全射の核が常に `d` で消えるという初等的な事実。 -/

/-- `coker(d·f) ↠ coker(f)` の**核は `d` で零化される**。 -/
theorem quot_comparison_kernel {R M : Type u} [CommRing R] [AddCommGroup M] [Module R M]
    (f : M →ₗ[R] M) (d : R) :
    ∃ q : (M ⧸ LinearMap.range (d • f)) →ₗ[R] (M ⧸ LinearMap.range f),
      Function.Surjective q ∧ ∀ x, q x = 0 → d • x = 0 := by
  have hle : LinearMap.range (d • f) ≤ LinearMap.range f := by
    rintro _ ⟨y, rfl⟩
    exact ⟨d • y, by simp [LinearMap.smul_apply, map_smul]⟩
  refine ⟨Submodule.mapQ _ _ LinearMap.id hle, ?_, ?_⟩
  · intro y
    induction y using Submodule.Quotient.induction_on with
    | H m => exact ⟨Submodule.Quotient.mk m, rfl⟩
  · intro x hx
    induction x using Submodule.Quotient.induction_on with
    | H m =>
      have hm : m ∈ LinearMap.range f := by
        have hq : Submodule.Quotient.mk (p := LinearMap.range f) m = 0 := hx
        exact (Submodule.Quotient.mk_eq_zero _).mp hq
      obtain ⟨y, rfl⟩ := hm
      show Submodule.Quotient.mk (d • f y) = 0
      rw [Submodule.Quotient.mk_eq_zero]
      exact ⟨y, by simp [LinearMap.smul_apply]⟩

/-- **"Tripling ε" の帰結を商の比較に繋ぐ**。`π e = s²·ē` なので
`π(s³·1 − e) = s²·(s·1 − ē)`——右辺が `quot_comparison_kernel` の
`d·f`(`d = s²`, `f = s·1 − ē`)の形になる。 -/
theorem reduction_of_idem_lift {A S S' : Type u} [CommRing A] [Ring S] [Ring S']
    [Algebra A S] [Algebra A S'] (π : S →ₐ[A] S') (s : A) (e : S) (ebar : S')
    (he : π e = algebraMap A S' (s^2) * ebar) :
    π (algebraMap A S (s^3) - e) = algebraMap A S' (s^2) * (algebraMap A S' s - ebar) := by
  rw [map_sub, AlgHom.commutes, he, mul_sub, ← map_mul]
  congr 2
  ring

/-- **`Theorem 2.3` の第 3・4 段(まとめ)**——`I² = 0` の上で
almost 冪等行列 `ē`(`ē² = s·ē`)は `e`(`e² = s³·e`)に持ち上がり、
しかも `s³·1 − e` の還元がちょうど `s²·(s·1 − ē)` になる。

前者が原文の *"Tripling `ε` we may assume that `ē` lifts to an `r×r`
matrix `e` with `e² = p^ε e`"*、後者が `quot_comparison_kernel` と
合わせて *"`B_ε` contains `B ⊗_Ā A/(p-torsion)`, the quotient being
annihilated by `p^ε`"* に当たる。 -/
theorem thm_2_3_lift_and_compare {A : Type u} [CommRing A] {ι : Type u} [Fintype ι]
    [DecidableEq ι]
    (I : Ideal A) (hI : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (s : A) (ebar : Matrix ι ι (A ⧸ I))
    (hebar : ebar^2 = algebraMap A (Matrix ι ι (A ⧸ I)) s * ebar) :
    ∃ e : Matrix ι ι A,
      e^2 = algebraMap A (Matrix ι ι A) (s^3) * e ∧
      (Ideal.Quotient.mkₐ A I).mapMatrix (algebraMap A (Matrix ι ι A) (s^3) - e)
        = algebraMap A (Matrix ι ι (A ⧸ I)) (s^2)
          * (algebraMap A (Matrix ι ι (A ⧸ I)) s - ebar) := by
  obtain ⟨e, he2, hered⟩ := exists_matrix_almost_idem_lift I hI s ebar hebar
  exact ⟨e, he2, reduction_of_idem_lift (Ideal.Quotient.mkₐ A I).mapMatrix s e ebar hered⟩

/-- 非空虚性——`A = ℤ/4`、`I = (2)`(**`I ≠ 0` かつ `I² = 0`**)、
`s = 1`、`ē = 1`。仮定(特に平方零イデアル)が空虚に真になっていない
ことの対照。 -/
example :
    ∃ e : Matrix (Fin 1) (Fin 1) (ZMod 4),
      e^2 = algebraMap (ZMod 4) (Matrix (Fin 1) (Fin 1) (ZMod 4)) ((1 : ZMod 4)^3) * e ∧
      (Ideal.Quotient.mkₐ (ZMod 4) (Ideal.span {(2 : ZMod 4)})).mapMatrix
          (algebraMap (ZMod 4) (Matrix (Fin 1) (Fin 1) (ZMod 4)) ((1 : ZMod 4)^3) - e)
        = algebraMap (ZMod 4) (Matrix (Fin 1) (Fin 1) (ZMod 4 ⧸ Ideal.span {(2 : ZMod 4)}))
            ((1 : ZMod 4)^2)
          * (algebraMap (ZMod 4) (Matrix (Fin 1) (Fin 1) (ZMod 4 ⧸ Ideal.span {(2 : ZMod 4)}))
              (1 : ZMod 4) - 1) := by
  refine thm_2_3_lift_and_compare (Ideal.span {(2 : ZMod 4)}) ?_ 1 1 (by simp)
  intro u hu v hv
  obtain ⟨a, rfl⟩ := Ideal.mem_span_singleton.mp hu
  obtain ⟨b, rfl⟩ := Ideal.mem_span_singleton.mp hv
  have h4 : (4 : ZMod 4) = 0 := rfl
  linear_combination (a * b) * h4

/-! ## 乗法の持ち上げと結合子の所在

原文の

> *Let `m : B ⊗_Ā B → B` denote the multiplication map. `p^{3ε}·m` lifts
> to `m_ε : B_ε ⊗_A B_ε → B_ε`. The obstruction to `m_ε` being associative
> is given by
> `b₀[m_ε(b₁, m_ε(b₂,b₃)) − m_ε(m_ε(b₁,b₂),b₃)]b₄ ∈ Kern(B_ε → B)`*

の 2 段。前段は `B_ε ⊗_A B_ε` の almost 射影性を全射 `ρ : B_ε ↠ B` に
沿って使うだけ、後段は `B` の結合律だけから出る。 -/

/-- **`ρ : B_ε ↠ B` の構成**——`red(g x) = c·f(red x)` なら、還元 `red` は
`coker g → coker f` を誘導し、全射になる。

`Theorem 2.3` では `g := s³·1 − e`(`A` 上)、`f := s·1 − ē`(`Ā` 上)、
`c := s²` と取る——`reduction_of_idem_lift` がちょうど `hrel` を与える。
これで乗法の持ち上げ `thm_2_3_mul_lift` に渡す全射 `ρ` が手に入る。 -/
theorem exists_quot_map {A M N : Type u} [CommRing A] [AddCommGroup M] [Module A M]
    [AddCommGroup N] [Module A N]
    (g : M →ₗ[A] M) (f : N →ₗ[A] N) (red : M →ₗ[A] N) (c : A)
    (hrel : ∀ x, red (g x) = c • f (red x)) (hsurj : Function.Surjective red) :
    ∃ ρ : (M ⧸ LinearMap.range g) →ₗ[A] (N ⧸ LinearMap.range f),
      Function.Surjective ρ
      ∧ ∀ x, ρ (Submodule.Quotient.mk x) = Submodule.Quotient.mk (red x) := by
  have hle : LinearMap.range g ≤ LinearMap.ker ((LinearMap.range f).mkQ ∘ₗ red) := by
    rintro _ ⟨x, rfl⟩
    show (LinearMap.range f).mkQ (red (g x)) = 0
    rw [hrel x, Submodule.mkQ_apply, Submodule.Quotient.mk_eq_zero]
    exact Submodule.smul_mem _ c ⟨red x, rfl⟩
  refine ⟨Submodule.liftQ _ ((LinearMap.range f).mkQ ∘ₗ red) hle, ?_, fun x => rfl⟩
  intro y
  induction y using Submodule.Quotient.induction_on with
  | H n =>
    obtain ⟨x, rfl⟩ := hsurj n
    exact ⟨Submodule.Quotient.mk x, rfl⟩

/-- **`Theorem 2.3` の乗法の持ち上げ**——原文の
*"`p^{3ε}·m` lifts to `m_ε : B_ε ⊗_A B_ε → B_ε`"*。 -/
theorem thm_2_3_mul_lift {A B : Type u} [CommRing A] [AddCommGroup B] [Module A B]
    {ι : Type u} [Fintype ι]
    (s : A) (E : Module.End A (ι → A)) (hE : E * E = s • E)
    (ρ : ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E)) →ₗ[A] B)
    (hρ : Function.Surjective ρ)
    (m : B →ₗ[A] B →ₗ[A] B) :
    ∃ μ : TensorProduct A
        ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E))
        ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E))
        →ₗ[A] ((ι → A) ⧸ LinearMap.range (s • (1 : Module.End A (ι → A)) - E)),
      ∀ z, ρ (μ z) = (s * s) •
        ((TensorProduct.lift m) (TensorProduct.map ρ ρ z)) := by
  obtain ⟨μ, hμ⟩ := almostIdemQuotient_tensor_almost_lift s E hE ρ hρ
    ((TensorProduct.lift m) ∘ₗ (TensorProduct.map ρ ρ))
  exact ⟨μ, hμ⟩

/-- **`Theorem 2.3`——結合子は `ρ` の核に入る**。原文の
*"The obstruction to `m_ε` being associative is given by
`b₀[...]b₄ ∈ Kern(B_ε → B)`"*。`B` の結合律だけから出る。 -/
theorem associator_mem_ker {A N B : Type u} [CommRing A] [AddCommGroup N] [Module A N]
    [AddCommGroup B] [Module A B]
    (ρ : N →ₗ[A] B) (m : B →ₗ[A] B →ₗ[A] B) (c : A)
    (μ : N →ₗ[A] N →ₗ[A] N)
    (hμ : ∀ x y : N, ρ (μ x y) = c • m (ρ x) (ρ y))
    (hassoc : ∀ a b d : B, m a (m b d) = m (m a b) d)
    (x y z : N) :
    ρ (μ x (μ y z) - μ (μ x y) z) = 0 := by
  simp only [map_sub, hμ, map_smul, LinearMap.smul_apply, smul_smul]
  rw [hassoc, sub_self]

/-- **交換子は導分**(Faltings `Theorem 2.3`: *"It also follows that the
product is commutative, as **the commutator with a fixed element of `B_ε`
is a derivation**"*)。結合律だけから出る:

    c_a(xy) = a(xy) − (xy)a = x·c_a(y) + c_a(x)·y

これと「導分が almost 消える」(`Ω[B⁄A]` の almost 消滅、
`AlmostDifferentials.kaehler_almost_zero`)を合わせると
`μ` の almost 可換性が出る。 -/
theorem commutator_leibniz {A N : Type u} [CommRing A] [AddCommGroup N] [Module A N]
    (μ : N →ₗ[A] N →ₗ[A] N)
    (hassoc : ∀ x y z : N, μ x (μ y z) = μ (μ x y) z)
    (a x y : N) :
    μ a (μ x y) - μ (μ x y) a
      = μ x (μ a y - μ y a) + μ (μ a x - μ x a) y := by
  rw [map_sub, map_sub, LinearMap.sub_apply]
  rw [hassoc x a y, hassoc x y a, hassoc a x y]
  abel

/-! ## ★補正した積から**本物の環**を作る(2026-09-05)

`associativity_correction` は「補正した積 `ν = t·μ − h` が**厳密に**
結合的」を与えた。可換性と単位元が揃えば、これで加群 `N` に
**本物の可換環構造**が入り、`A`-代数になる——原文の `B_ε` の構成
そのものである(`B_ε = A^r/(p^ε−e)(A^r)` は加群としては商加群で、
積は持ち上げた `m_ε` を補正したもの)。

Lean で「双線形写像から環を作る」道具は mathlib に無いので置く。 -/

/-- **双線形写像から可換環**——結合的・可換・単位元つきの双線形写像
`μ` は `N` に可換環構造を与える(積は `x * y = μ x y`)。 -/
@[reducible] noncomputable def commRingOfBilin {A N : Type u} [CommRing A] [AddCommGroup N]
    [Module A N] (μ : N →ₗ[A] N →ₗ[A] N) (one : N)
    (hassoc : ∀ x y z : N, μ x (μ y z) = μ (μ x y) z)
    (hcomm : ∀ x y : N, μ x y = μ y x)
    (hone : ∀ x : N, μ one x = x) : CommRing N :=
  { (inferInstance : AddCommGroup N) with
    mul := fun x y => μ x y
    mul_assoc := fun x y z => (hassoc x y z).symm
    one := one
    one_mul := hone
    mul_one := fun x => by show μ x one = x; rw [hcomm]; exact hone x
    left_distrib := fun x y z => by show μ x (y + z) = μ x y + μ x z; rw [map_add]
    right_distrib := fun x y z => by
      show μ (x + y) z = μ x z + μ y z
      rw [map_add]; rfl
    zero_mul := fun x => by show μ 0 x = 0; rw [map_zero]; rfl
    mul_zero := fun x => by show μ x 0 = 0; rw [map_zero]
    mul_comm := hcomm }

/-- **その環は `A`-代数**——`μ` が `A`-双線形なので
`Algebra.ofModule` がそのまま通る。 -/
@[reducible] noncomputable def algebraOfBilin {A N : Type u} [CommRing A] [AddCommGroup N]
    [Module A N] (μ : N →ₗ[A] N →ₗ[A] N) (one : N)
    (hassoc : ∀ x y z : N, μ x (μ y z) = μ (μ x y) z)
    (hcomm : ∀ x y : N, μ x y = μ y x)
    (hone : ∀ x : N, μ one x = x) :
    letI := commRingOfBilin μ one hassoc hcomm hone
    Algebra A N :=
  letI := commRingOfBilin μ one hassoc hcomm hone
  Algebra.ofModule
    (fun r x y => by show μ (r • x) y = r • μ x y; rw [map_smul]; rfl)
    (fun r x y => by show μ x (r • y) = r • μ x y; rw [map_smul])

/-- **★補正した積から環を作る**——`Theorem 2.3` 第 4 段の帰結を
そのまま環にする。`hcob` は `hochschild_H3_almost_coboundary_act` が
与えるコバウンダリ表示、`hv1`/`hv2` は `associator_eq_zero_of_act(')`
から出る「`h` が二次で消える」性質。 -/
@[reducible] noncomputable def commRingOfCorrectedMul {A N : Type u} [CommRing A]
    [AddCommGroup N] [Module A N]
    (μ : N →ₗ[A] N →ₗ[A] N) (t : A) (h : N →ₗ[A] N →ₗ[A] N) (one : N)
    (hv2 : ∀ x y z : N, h x (h y z) = 0)
    (hv1 : ∀ x y z : N, h (h x y) z = 0)
    (hcob : ∀ x y z : N, t • (μ x (μ y z) - μ (μ x y) z)
      = μ x (h y z) - h (μ x y) z + h x (μ y z) - μ (h x y) z)
    (hcomm : ∀ x y : N, (t • μ - h) x y = (t • μ - h) y x)
    (hone : ∀ x : N, (t • μ - h) one x = x) : CommRing N :=
  commRingOfBilin (t • μ - h) one
    (associativity_correction μ t h hv2 hv1 hcob) hcomm hone

/-- 非空虚性——`N := A`、`μ := (· * ·)`、`one := 1` で
`commRingOfBilin` はもとの可換環をそのまま与える。仮定が空虚に
真になっていないことの対照。 -/
example : (1 : Polynomial ℤ) = 1 := by
  letI := commRingOfBilin (LinearMap.mul (Polynomial ℤ) (Polynomial ℤ)) (1 : Polynomial ℤ)
    (fun x y z => (mul_assoc x y z).symm) mul_comm one_mul
  rfl

/-! ## 最終段——`C_σ = A + p^{5σ}·B_σ` の増大列とその合併

原文の

> *Consider the subalgebras `C_σ = A + p^{5σ}B_σ ⊆ B̃[1/p]`. If `σ < ε/2`
> then `C_ε ⊆ C_σ`. ... If we choose a sequence `εₙ` with `ε_{n+1} < εₙ/2`,
> the corresponding `C_{εₙ}`'s form an increasing sequence. If `B̃` denotes
> their union, `B̃ ⊗_A B̃` contains `m·e_{B̃/A}`, hence `B̃` is the required
> lifting of `B`.*

の骨格。`A + c·B` が部分代数になるのは
`(a + cb)(a' + cb') = aa' + c(ab' + a'b + cbb')` から。増大性は
`c ∈ c'·B` から出る(`σ < ε/2` のときの `p^{5ε} ∈ p^{5σ}·B` に対応)。 -/

/-- **`C_σ := A + c·B`**——Faltings `Theorem 2.3` 最終段の部分代数。 -/
def scaledSubalgebra {A B : Type u} [CommRing A] [CommRing B] [Algebra A B] (c : B) :
    Subalgebra A B where
  carrier := {x | ∃ (a : A) (b : B), x = algebraMap A B a + c * b}
  mul_mem' := by
    rintro _ _ ⟨a, b, rfl⟩ ⟨a', b', rfl⟩
    exact ⟨a * a', algebraMap A B a * b' + algebraMap A B a' * b + c * (b * b'),
      by rw [map_mul]; ring⟩
  add_mem' := by
    rintro _ _ ⟨a, b, rfl⟩ ⟨a', b', rfl⟩
    exact ⟨a + a', b + b', by rw [map_add]; ring⟩
  algebraMap_mem' := fun a => ⟨a, 0, by rw [mul_zero, add_zero]⟩

@[simp] theorem mem_scaledSubalgebra {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (c x : B) :
    x ∈ scaledSubalgebra (A := A) c ↔ ∃ (a : A) (b : B), x = algebraMap A B a + c * b :=
  Iff.rfl

/-- **増大性**——`c ∈ c'·B` なら `A + c·B ⊆ A + c'·B`
(原文の *"If `σ < ε/2` then `C_ε ⊆ C_σ`"*)。 -/
theorem scaledSubalgebra_mono {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (c c' : B) (h : ∃ d : B, c = c' * d) :
    scaledSubalgebra (A := A) c ≤ scaledSubalgebra (A := A) c' := by
  obtain ⟨d, rfl⟩ := h
  rintro _ ⟨a, b, rfl⟩
  exact ⟨a, d * b, by ring⟩

/-- **増大列の合併**——原文の *"the corresponding `C_{εₙ}`'s form an
increasing sequence. If `B̃` denotes their union..."*。 -/
theorem mem_iSup_scaledSubalgebra {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (c : ℕ → B) (hmono : ∀ n, ∃ d : B, c n = c (n+1) * d) (x : B) :
    x ∈ (⨆ n, scaledSubalgebra (A := A) (c n)) ↔ ∃ n, x ∈ scaledSubalgebra (A := A) (c n) := by
  have hmono' : Monotone (fun n => scaledSubalgebra (A := A) (c n)) :=
    monotone_nat_of_le_succ (fun n => scaledSubalgebra_mono _ _ (hmono n))
  have hcoe := Subalgebra.coe_iSup_of_directed hmono'.directed_le
  constructor
  · intro hx
    have hx' : x ∈ (↑(⨆ n, scaledSubalgebra (A := A) (c n)) : Set B) := hx
    rw [hcoe] at hx'
    simpa using hx'
  · rintro ⟨n, hn⟩
    exact le_iSup (fun n => scaledSubalgebra (A := A) (c n)) n hn

/-- **像がテンソルへ伝わる**——`C` の像が `c·B` を含めば、`C ⊗ C` の像は
`c²·(B ⊗ B)` を含む。原文の *"The image of `C_σ` in `B` contains
`p^{6σ}·B`. It follows that `p^{7σ}·e_{B/Ā}` lies in `C_σ ⊗_A C_σ`"* の
代数的核(条件(iii)の witness が `C_σ ⊗ C_σ` に入ることを言う段)。 -/
theorem tensor_range_of_range {A B C : Type u} [CommRing A]
    [AddCommGroup B] [Module A B] [AddCommGroup C] [Module A C]
    (f : C →ₗ[A] B) (c : A) (h : ∀ b : B, c • b ∈ LinearMap.range f)
    (z : TensorProduct A B B) :
    (c * c) • z ∈ LinearMap.range (TensorProduct.map f f) := by
  induction z using TensorProduct.induction_on with
  | zero => simp
  | tmul b b' =>
    obtain ⟨x, hx⟩ := h b
    obtain ⟨y, hy⟩ := h b'
    refine ⟨x ⊗ₜ[A] y, ?_⟩
    rw [TensorProduct.map_tmul, hx, hy, TensorProduct.smul_tmul', TensorProduct.tmul_smul,
      ← smul_smul]
    exact (TensorProduct.smul_tmul' c (c • b) b').symm
  | add u v hu hv =>
    rw [smul_add]
    exact Submodule.add_mem _ hu hv

/-! ## エタール持ち上げの一意性——*"all `B_ε[1/p]` are isomorphic"*

原文の

> *Now each `B_ε[1/p]` is a lifting of the étale `A[1/p]`-algebra `B[1/p]`
> and is a projective `A[1/p]`-module. **It follows that all `B_ε[1/p]` are
> isomorphic to a fixed algebra `B̃[1/p]`** (say)*

は古典的な事実(エタール代数は冪零イデアルに沿って一意に持ち上がる)。
mathlib に `Algebra.FormallySmooth.exists_lift`(存在)と
`Algebra.FormallyUnramified.ext`(一意性)があるので、両方を使うだけで
出る——**`Theorem 2.2` はこれの almost 版**にあたる。 -/

/-- **エタール持ち上げの存在と一意性**(古典版)。 -/
theorem etale_lift_existsUnique {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C] [Algebra.FormallyEtale A B]
    (I : Ideal C) (hI : IsNilpotent I) (φ : B →ₐ[A] (C ⧸ I)) :
    ∃! ψ : B →ₐ[A] C, (Ideal.Quotient.mkₐ A I).comp ψ = φ := by
  obtain ⟨ψ, hψ⟩ := Algebra.FormallySmooth.exists_lift I hI φ
  refine ⟨ψ, hψ, fun ψ' hψ' => ?_⟩
  refine Algebra.FormallyUnramified.ext (I := I) hI (fun x => ?_)
  have h1 := congrArg (fun f : B →ₐ[A] (C ⧸ I) => f x) hψ'
  have h2 := congrArg (fun f : B →ₐ[A] (C ⧸ I) => f x) hψ
  simpa using h1.trans h2.symm

/-- **エタール持ち上げは同型を除いて一意**——原文の
*"It follows that all `B_ε[1/p]` are isomorphic to a fixed algebra
`B̃[1/p]`"*。 -/
theorem etale_lift_iso {A B₁ B₂ : Type u} [CommRing A] [CommRing B₁] [CommRing B₂]
    [Algebra A B₁] [Algebra A B₂] [Algebra.FormallyEtale A B₁] [Algebra.FormallyEtale A B₂]
    (I₁ : Ideal B₁) (I₂ : Ideal B₂) (hI₁ : IsNilpotent I₁) (hI₂ : IsNilpotent I₂)
    (e : (B₁ ⧸ I₁) ≃ₐ[A] (B₂ ⧸ I₂)) :
    ∃ (f : B₁ →ₐ[A] B₂) (g : B₂ →ₐ[A] B₁),
      g.comp f = AlgHom.id A B₁ ∧ f.comp g = AlgHom.id A B₂ := by
  obtain ⟨f, hf⟩ := Algebra.FormallySmooth.exists_lift (R := A) (A := B₁) I₂ hI₂
    (e.toAlgHom.comp (Ideal.Quotient.mkₐ A I₁))
  obtain ⟨g, hg⟩ := Algebra.FormallySmooth.exists_lift (R := A) (A := B₂) I₁ hI₁
    (e.symm.toAlgHom.comp (Ideal.Quotient.mkₐ A I₂))
  have hf' : ∀ b, Ideal.Quotient.mk I₂ (f b) = e (Ideal.Quotient.mk I₁ b) :=
    fun b => congrArg (fun h : B₁ →ₐ[A] (B₂ ⧸ I₂) => h b) hf
  have hg' : ∀ c, Ideal.Quotient.mk I₁ (g c) = e.symm (Ideal.Quotient.mk I₂ c) :=
    fun c => congrArg (fun h : B₂ →ₐ[A] (B₁ ⧸ I₁) => h c) hg
  refine ⟨f, g, ?_, ?_⟩
  · refine Algebra.FormallyUnramified.ext (I := I₁) hI₁ (fun b => ?_)
    show Ideal.Quotient.mk I₁ (g (f b)) = Ideal.Quotient.mk I₁ b
    rw [hg' (f b), hf' b, AlgEquiv.symm_apply_apply]
  · refine Algebra.FormallyUnramified.ext (I := I₂) hI₂ (fun c => ?_)
    show Ideal.Quotient.mk I₂ (f (g c)) = Ideal.Quotient.mk I₂ c
    rw [hf' (g c), hg' c, AlgEquiv.apply_symm_apply]

end ABC3.Found.Falt1
