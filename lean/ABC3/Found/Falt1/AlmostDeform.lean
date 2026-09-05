import ABC3.Found.Falt1.AlmostDerivation
import Mathlib.Data.ZMod.Basic

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

end ABC3.Found.Falt1
