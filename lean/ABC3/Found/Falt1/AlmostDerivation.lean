import ABC3.Found.Falt1.AlmostLifting

/-!
# [Falt1] 捩れを持つ加群への持ち上げ——`Theorem 2.4(i)` の核側へ(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §2、
`Theorem 2.4(i)` の証明(物理 p.8 = 印字 p.261)。

> *For (i) we have to show that for any `B`-module `M` (**which is allowed
> to have `p`-torsion contrary to our general assumptions**), the restriction
> map on derivations `Der(B, M) → Der(A, M)` is an almost isomorphism.
> This follows in the usual way using Hochschild cohomology `H*(B/A, M)`.*

**原文が明示的に「`M` は `p`-捩れを持ってよい」と断っている**のがこの
ファイルの出発点である。`AlmostLifting.lean` の `exists_multiplicative_lift`
は `htors`(`C` の `c`-捩れ無し)を使っていたので、そのままでは
`C := B ⊕ M`(`M` が捩れを持つ)に適用できない。

## 捩れ無し仮定の外し方

`htors` が使われていたのは 1 箇所だけ——`obstruction_isCocycle` で
`c·(δOb) = 0` から `c` を割って `δOb = 0` を出すところである
(`obstruction_identity` が与えるのは `c` 倍の等式だけ)。

**割らずに `c•Ob` を取れば、それ自体が honest な 2-コサイクルになる**
(`obstruction_isCocycle_smul`)。水準が `c` 1 つ分だけ悪くなるが、
almost な主張では `ε` を小さく取り直せば済む。

これにより `exists_multiplicative_lift_nt`——`htors` を仮定しない
乗法的持ち上げ——が得られる。
-/

namespace ABC3.Found.Falt1

universe u

/-- `exists_multiplicative_lift` の中核を、コバウンダリ表示 `(h, hh)` を
入力として切り出したもの。**捩れ無し仮定 `htors` を使わない**——
`htors` はコサイクル性を出すためだけに使われていたので、コバウンダリ
表示を外から渡せば不要になる。 -/
theorem exists_multiplicative_lift_of_coboundary {A B C : Type u}
    [CommRing A] [CommRing B] [CommRing C] [Algebra A B] [Algebra A C]
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (φ : B →ₐ[A] (C ⧸ I))
    (c : A) (ψ : B →ₗ[A] C) (hψ : ∀ b, Ideal.Quotient.mk I (ψ b) = c • φ b)
    (t : A) :
    letI : Module (C ⧸ I) I := sqZeroModule I hsq
    letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
    ∀ (h : B →ₗ[A] I),
      (∀ b₁ b₂ : B, t • (obsBil I φ c ψ hψ) b₁ b₂
        = b₁ • h b₂ - h (b₁*b₂) + b₂ • h b₁) →
      ∃ ψ₂ : B →ₗ[A] C,
        (∀ b, Ideal.Quotient.mk I (ψ₂ b) = ((c*t)*(c*t)) • φ b) ∧
        (∀ x y, ((c*t)*(c*t)) • ψ₂ (x*y) = ψ₂ x * ψ₂ y) := by
  letI : Module (C ⧸ I) I := sqZeroModule I hsq
  letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
  intro h hh
  set Ob := obsBil I φ c ψ hψ with hObdef
  have hObC : ∀ b₁ b₂ : B, ((Ob b₁ b₂ : I) : C) = c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂ :=
    fun b₁ b₂ => rfl
  have hact := sqZero_act_eq I hsq φ c ψ hψ
  set ĥ : B →ₗ[A] C := t • ((I.subtype.restrictScalars A) ∘ₗ h) with hĥ
  have hIsq2 : ∀ b₁ b₂ : B, ĥ b₁ * ĥ b₂ = 0 := by
    intro b₁ b₂
    show (t • ((h b₁ : I) : C)) * (t • ((h b₂ : I) : C)) = 0
    have hz := hsq _ (h b₁).2 _ (h b₂).2
    rw [Algebra.smul_def, Algebra.smul_def]
    linear_combination (algebraMap A C t * algebraMap A C t) * hz
  have hcob : ∀ b₁ b₂ : B,
      (c*t) • ((c*t) • ((t • ψ) (b₁ * b₂)) - (t • ψ) b₁ * (t • ψ) b₂)
      = (t • ψ) b₁ * ĥ b₂ - (c*t) • ĥ (b₁ * b₂) + (t • ψ) b₂ * ĥ b₁ := by
    intro b₁ b₂
    have hh12 := congrArg (fun x : I => (x : C)) (hh b₁ b₂)
    have ha1 := hact b₁ (h b₂)
    have ha2 := hact b₂ (h b₁)
    have hO := hObC b₁ b₂
    simp only [hĥ, LinearMap.smul_apply, LinearMap.coe_comp, Function.comp_apply,
      LinearMap.coe_restrictScalars, Submodule.coe_subtype,
      Submodule.coe_add, Submodule.coe_sub, Submodule.coe_smul_of_tower,
      Algebra.smul_def, map_mul] at hh12 ha1 ha2 hO ⊢
    linear_combination (-(algebraMap A C c * (algebraMap A C t)^3)) * hO
      + (algebraMap A C c * (algebraMap A C t)^2) * hh12
      - ((algebraMap A C t)^2) * ha1 - ((algebraMap A C t)^2) * ha2
  refine ⟨(c*t) • (t • ψ) + ĥ, fun b => ?_, fun x y => ?_⟩
  · show Ideal.Quotient.mk I ((c*t) • (t • ψ b) + t • ((h b : I) : C)) = ((c*t)*(c*t)) • φ b
    have hzero : Ideal.Quotient.mk I ((h b : I) : C) = 0 :=
      (Ideal.Quotient.eq_zero_iff_mem).mpr (h b).2
    have hmk : ∀ x : A, Ideal.Quotient.mk I (algebraMap A C x) = algebraMap A (C ⧸ I) x :=
      fun _ => rfl
    simp only [Algebra.smul_def, map_add, map_mul, hmk, hzero, mul_zero, add_zero, hψ b]
    ring
  · exact doubling_multiplicative (c*t) (t • ψ) ĥ hIsq2 hcob x y

/-- 平方零イデアル `I` 上では、`B` 作用と `A` 作用が可換
(`A → B` を経由しているから)。 -/
theorem sqZero_smul_comm {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (φ : B →ₐ[A] (C ⧸ I)) (c : A) :
    letI : Module (C ⧸ I) I := sqZeroModule I hsq
    letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
    ∀ (v : B) (z : I), v • (c • z) = c • (v • z) := by
  letI : Module (C ⧸ I) I := sqZeroModule I hsq
  letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
  haveI : IsScalarTower A B I := sqZeroTower I hsq φ
  intro v z
  rw [← algebraMap_smul B c z, ← algebraMap_smul B c (v • z), smul_smul, smul_smul, mul_comm]

/-- **障害の `c` 倍は、捩れ無し仮定なしで honest な 2-コサイクル**。

`obstruction_isCocycle` は `obstruction_identity` の与える `c·(δOb) = 0`
から `htors` で `c` を割って `δOb = 0` を出していた。割らずに `c•Ob` を
取れば、同じ等式がそのままコサイクル条件になる。 -/
theorem obstruction_isCocycle_smul {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (φ : B →ₐ[A] (C ⧸ I)) (c : A) (ψ : B →ₗ[A] C)
    (hψ : ∀ b : B, Ideal.Quotient.mk I (ψ b) = c • φ b) :
    letI : Module (C ⧸ I) I := sqZeroModule I hsq
    letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
    ∀ (Ob : B →ₗ[A] B →ₗ[A] I),
      (∀ b₁ b₂ : B, ((Ob b₁ b₂ : I) : C) = c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂) →
      ∀ v b₁ b₂ : B, v • (c • Ob) b₁ b₂ - (c • Ob) (v * b₁) b₂
        + (c • Ob) v (b₁ * b₂) - b₂ • (c • Ob) v b₁ = 0 := by
  letI : Module (C ⧸ I) I := sqZeroModule I hsq
  letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
  haveI : IsScalarTower A B I := sqZeroTower I hsq φ
  intro Ob hOb v b₁ b₂
  have hcomm := sqZero_smul_comm I hsq φ c
  have hkey : v • (c • Ob) b₁ b₂ - (c • Ob) (v * b₁) b₂
      + (c • Ob) v (b₁ * b₂) - b₂ • (c • Ob) v b₁
      = c • (v • Ob b₁ b₂ - Ob (v * b₁) b₂ + Ob v (b₁ * b₂) - b₂ • Ob v b₁) := by
    simp only [LinearMap.smul_apply, hcomm, smul_sub, smul_add]
  rw [hkey]
  apply Subtype.ext
  have hact := sqZero_act_eq I hsq φ c ψ hψ
  have hid := obstruction_identity c ψ v b₁ b₂
  have e1 := (hact v (Ob b₁ b₂)).symm
  have e2 := (hact b₂ (Ob v b₁)).symm
  have h1 := hOb b₁ b₂
  have h2 := hOb (v * b₁) b₂
  have h3 := hOb v (b₁ * b₂)
  have h4 := hOb v b₁
  simp only [Submodule.coe_sub, Submodule.coe_add, Submodule.coe_smul_of_tower,
    Submodule.coe_zero, Algebra.smul_def] at hid e1 e2 h1 h2 h3 h4 ⊢
  rw [← mul_assoc v b₁ b₂] at hid h3
  linear_combination e1 - e2 + (ψ v) * h1 - (algebraMap A C c) * h2
    + (algebraMap A C c) * h3 - (ψ b₂) * h4 + hid

/-- **`exists_multiplicative_lift` の捩れ無し版**。`C` に `c`-捩れが
あってもよい——`Theorem 2.4(i)` が要求する「`M` は `p`-捩れを持ってよい」
という状況にそのまま使える。水準は `c·t·c` になる(元の版は `c·t`)。 -/
theorem exists_multiplicative_lift_nt {A B C : Type u} [CommRing A] [CommRing B] [CommRing C]
    [Algebra A B] [Algebra A C]
    (I : Ideal C) (hsq : ∀ u ∈ I, ∀ v ∈ I, u * v = 0)
    (φ : B →ₐ[A] (C ⧸ I))
    (c : A) (ψ : B →ₗ[A] C) (hψ : ∀ b, Ideal.Quotient.mk I (ψ b) = c • φ b)
    (t : A) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = t • (1 : B)) :
    ∃ ψ₂ : B →ₗ[A] C,
      (∀ b, Ideal.Quotient.mk I (ψ₂ b) = ((c*(t*c))*(c*(t*c))) • φ b) ∧
      (∀ x y, ((c*(t*c))*(c*(t*c))) • ψ₂ (x*y) = ψ₂ x * ψ₂ y) := by
  letI : Module (C ⧸ I) I := sqZeroModule I hsq
  letI : Module B I := Module.compHom _ (φ.toRingHom : B →+* (C ⧸ I))
  haveI : IsScalarTower A B I := sqZeroTower I hsq φ
  set Ob := obsBil I φ c ψ hψ with hObdef
  have hObC : ∀ b₁ b₂ : B, ((Ob b₁ b₂ : I) : C) = c • ψ (b₁ * b₂) - ψ b₁ * ψ b₂ :=
    fun b₁ b₂ => rfl
  have hcoc := obstruction_isCocycle_smul I hsq φ c ψ hψ Ob hObC
  obtain ⟨h, hh⟩ := hochschild_H2_almost_coboundary t w hw_ann hw_aug (c • Ob) hcoc
  refine exists_multiplicative_lift_of_coboundary I hsq φ c ψ hψ (t*c) h (fun b₁ b₂ => ?_)
  have hb := hh b₁ b₂
  rwa [LinearMap.smul_apply, LinearMap.smul_apply, smul_smul] at hb

/-! ## `Theorem 2.4(i)`——導分への還元(Faltings 自身の論法)

> *For (i) we have to show that for any `B`-module `M` ..., the restriction
> map on derivations `Der(B, M) → Der(A, M)` is an almost isomorphism.
> This follows in the usual way using Hochschild cohomology `H*(B/A, M)`.*

これを**完全に初等的に**実行する。鍵は次の1つの構成である:

almost dual basis `(fᵢ, yᵢ)`(`Σᵢ fᵢ(b)·yᵢ = c·b`)と `d ∈ Der_R(A,M)` から

    u(b) := Σᵢ yᵢ · d(fᵢ b)

と置くと、`fᵢ` の `A`-線型性と `d` の Leibniz 則だけから

    u(a·b) = a·u(b) + (c·b)·d a        (`c`-ねじれ半線型性)

が出る。このねじれのおかげで**コバウンダリ `δu` が `A`-双線型**になり、
`δ²=0` から Hochschild 2-コサイクルになる(捩れ無し仮定は一切要らない
——原文が「`M` は `p`-捩れを持ってよい」と断っている状況をそのまま扱える)。
そこで `hochschild_H2_almost_coboundary` により `t·δu = δh` と書ければ

    D := t·u - h

は `δD = 0`、すなわち **`Der_R(B,M)` の元**であり、`D 1 = 0` から
`D|_A = (t·c)·d` が従う。

最後に `M := B ⊗_A Ω[A⁄R]`、`d := (a ↦ 1 ⊗ d_A a)` と取れば、
`D` の誘導する `g : Ω[B⁄R] → B ⊗_A Ω[A⁄R]` は `mapBaseChange` の
**almost な左逆**になり、核が `t·c` で零化されることが出る。
余核側は既証(`AlmostDifferentials.thm_2_4_i_cokernel`)。
-/

section Derivation

variable {R A B M : Type u} [CommRing R] [CommRing A] [CommRing B]
  [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
  [AddCommGroup M] [Module R M] [Module A M] [Module B M]
  [IsScalarTower A B M] [IsScalarTower R A M]

/-- `A` の作用は `B` の作用と可換(`A → B` を経由しているから)。 -/
theorem smul_A_comm (a : A) (m : M) (b : B) : b • (a • m) = (a • b) • m := by
  rw [← algebraMap_smul B a m, ← algebraMap_smul B a b, smul_smul, smul_eq_mul, mul_comm]

/-- `u b := Σᵢ yᵢ • d(fᵢ b)`——almost dual basis と `d` から作る 1-コチェイン。 -/
noncomputable def twistCochain {ι : Type u} [Fintype ι]
    (fi : ι → (B →ₗ[A] A)) (y : ι → B) (d : Derivation R A M) : B → M :=
  fun b => ∑ i, y i • d (fi i b)

theorem twistCochain_add {ι : Type u} [Fintype ι]
    (fi : ι → (B →ₗ[A] A)) (y : ι → B) (d : Derivation R A M) (b b' : B) :
    twistCochain fi y d (b + b') = twistCochain fi y d b + twistCochain fi y d b' := by
  simp only [twistCochain, map_add]
  rw [← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl (fun i _ => by rw [smul_add])

/-- **`u` の `c`-ねじれ `A`-半線型性**: `u(a·b) = a·u(b) + (c·b)·(d a)`。
`fᵢ` の `A`-線型性と `d` の Leibniz 則だけから出る。 -/
theorem twistCochain_smul {ι : Type u} [Fintype ι]
    (fi : ι → (B →ₗ[A] A)) (y : ι → B) (d : Derivation R A M)
    (c : A) (hfy : ∀ b : B, ∑ i, (fi i b) • y i = c • b) (a : A) (b : B) :
    twistCochain fi y d (a • b) = a • twistCochain fi y d b + (c • b) • d a := by
  have hstep : ∀ i : ι, y i • d (fi i (a • b))
      = a • (y i • d (fi i b)) + ((fi i b) • y i) • d a := by
    intro i
    rw [map_smul, smul_eq_mul, d.leibniz, smul_add,
      smul_A_comm a (d (fi i b)) (y i), smul_A_comm (fi i b) (d a) (y i), smul_assoc]
  simp only [twistCochain]
  rw [Finset.sum_congr rfl (fun i _ => hstep i), Finset.sum_add_distrib,
    ← Finset.smul_sum, ← Finset.sum_smul, hfy]

theorem twistCob_add1 {ι : Type u} [Fintype ι]
    (fi : ι → (B →ₗ[A] A)) (y : ι → B) (d : Derivation R A M) (b₁ b₁' b₂ : B) :
    (b₁ + b₁') • twistCochain fi y d b₂ - twistCochain fi y d ((b₁ + b₁') * b₂)
      + b₂ • twistCochain fi y d (b₁ + b₁')
    = (b₁ • twistCochain fi y d b₂ - twistCochain fi y d (b₁ * b₂)
        + b₂ • twistCochain fi y d b₁)
      + (b₁' • twistCochain fi y d b₂ - twistCochain fi y d (b₁' * b₂)
        + b₂ • twistCochain fi y d b₁') := by
  rw [add_smul, add_mul, twistCochain_add, twistCochain_add, smul_add]
  abel

theorem twistCob_symm {ι : Type u} [Fintype ι]
    (fi : ι → (B →ₗ[A] A)) (y : ι → B) (d : Derivation R A M) (b₁ b₂ : B) :
    b₁ • twistCochain fi y d b₂ - twistCochain fi y d (b₁ * b₂)
      + b₂ • twistCochain fi y d b₁
    = b₂ • twistCochain fi y d b₁ - twistCochain fi y d (b₂ * b₁)
      + b₁ • twistCochain fi y d b₂ := by
  rw [mul_comm b₂ b₁]
  abel

theorem twistCob_smul1 {ι : Type u} [Fintype ι]
    (fi : ι → (B →ₗ[A] A)) (y : ι → B) (d : Derivation R A M)
    (c : A) (hfy : ∀ b : B, ∑ i, (fi i b) • y i = c • b) (a : A) (b₁ b₂ : B) :
    (a • b₁) • twistCochain fi y d b₂ - twistCochain fi y d ((a • b₁) * b₂)
      + b₂ • twistCochain fi y d (a • b₁)
    = a • (b₁ • twistCochain fi y d b₂ - twistCochain fi y d (b₁ * b₂)
        + b₂ • twistCochain fi y d b₁) := by
  rw [smul_mul_assoc, twistCochain_smul fi y d c hfy a (b₁ * b₂),
    twistCochain_smul fi y d c hfy a b₁, smul_add]
  have h1 : (a • b₁) • twistCochain fi y d b₂ = a • (b₁ • twistCochain fi y d b₂) :=
    smul_assoc a b₁ _
  have h2 : b₂ • ((c • b₁) • d a) = (c • (b₁ * b₂)) • d a := by
    rw [smul_smul, mul_smul_comm, mul_comm b₂ b₁]
  have h3 : b₂ • (a • twistCochain fi y d b₁) = a • (b₂ • twistCochain fi y d b₁) := by
    rw [smul_A_comm, smul_assoc]
  rw [h1, h3, h2]
  simp only [smul_add, smul_sub]
  abel

/-- コバウンダリ `δu(b₁,b₂) = b₁·u b₂ - u(b₁b₂) + b₂·u b₁` は `A`-**双線型**
——`u` 自身は `A`-線型でないが、`c`-ねじれ項がちょうど打ち消し合う。 -/
noncomputable def twistCoboundary {ι : Type u} [Fintype ι]
    (fi : ι → (B →ₗ[A] A)) (y : ι → B) (d : Derivation R A M)
    (c : A) (hfy : ∀ b : B, ∑ i, (fi i b) • y i = c • b) : B →ₗ[A] B →ₗ[A] M :=
  LinearMap.mk₂ A (fun b₁ b₂ => b₁ • twistCochain fi y d b₂
      - twistCochain fi y d (b₁ * b₂) + b₂ • twistCochain fi y d b₁)
    (fun b₁ b₁' b₂ => twistCob_add1 fi y d b₁ b₁' b₂)
    (fun a b₁ b₂ => twistCob_smul1 fi y d c hfy a b₁ b₂)
    (fun b₁ b₂ b₂' => by
      rw [twistCob_symm fi y d b₁ (b₂ + b₂'), twistCob_add1 fi y d b₂ b₂' b₁,
        twistCob_symm fi y d b₂ b₁, twistCob_symm fi y d b₂' b₁])
    (fun a b₁ b₂ => by
      rw [twistCob_symm fi y d b₁ (a • b₂), twistCob_smul1 fi y d c hfy a b₂ b₁,
        twistCob_symm fi y d b₂ b₁])

/-- `δ² = 0`——`δu` は Hochschild 2-コサイクル(捩れ無し仮定は要らない)。 -/
theorem twistCoboundary_isCocycle {ι : Type u} [Fintype ι]
    (fi : ι → (B →ₗ[A] A)) (y : ι → B) (d : Derivation R A M)
    (c : A) (hfy : ∀ b : B, ∑ i, (fi i b) • y i = c • b) (v b₁ b₂ : B) :
    v • twistCoboundary fi y d c hfy b₁ b₂ - twistCoboundary fi y d c hfy (v * b₁) b₂
      + twistCoboundary fi y d c hfy v (b₁ * b₂)
      - b₂ • twistCoboundary fi y d c hfy v b₁ = 0 := by
  set U := twistCochain fi y d with hU
  show v • (b₁ • U b₂ - U (b₁ * b₂) + b₂ • U b₁)
      - ((v * b₁) • U b₂ - U ((v * b₁) * b₂) + b₂ • U (v * b₁))
      + (v • U (b₁ * b₂) - U (v * (b₁ * b₂)) + (b₁ * b₂) • U v)
      - b₂ • (v • U b₁ - U (v * b₁) + b₁ • U v) = 0
  rw [← mul_assoc v b₁ b₂]
  simp only [smul_sub, smul_add, smul_smul]
  rw [mul_comm b₂ v, mul_comm b₂ b₁]
  abel

/-- **Faltings `Theorem 2.4(i)` の核心——導分の almost 延長**。
`d ∈ Der_R(A,M)` に対し `(t·c)·d` は `Der_R(B,M)` に延びる。
`M` が `p`-捩れを持っていてもよい。 -/
theorem exists_derivation_extension {ι : Type u} [Fintype ι]
    (fi : ι → (B →ₗ[A] A)) (y : ι → B)
    (c : A) (hfy : ∀ b : B, ∑ i, (fi i b) • y i = c • b)
    (t : A) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = t • (1 : B))
    (d : Derivation R A M) :
    ∃ D : Derivation R B M, ∀ a : A, D (algebraMap A B a) = (t * c) • d a := by
  obtain ⟨h, hh⟩ := hochschild_H2_almost_coboundary t w hw_ann hw_aug
    (twistCoboundary fi y d c hfy) (twistCoboundary_isCocycle fi y d c hfy)
  have hhU : ∀ b₁ b₂ : B,
      t • (b₁ • twistCochain fi y d b₂ - twistCochain fi y d (b₁ * b₂)
        + b₂ • twistCochain fi y d b₁)
      = b₁ • h b₂ - h (b₁ * b₂) + b₂ • h b₁ := hh
  have hRsmul : ∀ (r : R) (b : B), twistCochain fi y d (r • b) = r • twistCochain fi y d b := by
    intro r b
    rw [← algebraMap_smul A r b, twistCochain_smul fi y d c hfy,
      Derivation.map_algebraMap, smul_zero, add_zero, algebraMap_smul]
  have hone : t • twistCochain fi y d 1 - h 1 = 0 := by
    have hk := hhU 1 1
    simp only [one_mul, one_smul, smul_sub, smul_add] at hk
    linear_combination (norm := abel) hk
  refine ⟨{
    toLinearMap := {
      toFun := fun b => t • twistCochain fi y d b - h b
      map_add' := fun b b' => by
        rw [twistCochain_add, map_add, smul_add]; abel
      map_smul' := fun r b => by
        simp only [RingHom.id_apply]
        rw [hRsmul, h.map_smul_of_tower, smul_comm t r, smul_sub]
    }
    map_one_eq_zero' := hone
    leibniz' := ?_
  }, ?_⟩
  case refine_1 =>
    intro b₁ b₂
    show t • twistCochain fi y d (b₁ * b₂) - h (b₁ * b₂)
      = b₁ • (t • twistCochain fi y d b₂ - h b₂) + b₂ • (t • twistCochain fi y d b₁ - h b₁)
    have hk := hhU b₁ b₂
    simp only [smul_sub, smul_add] at hk ⊢
    rw [smul_comm b₁ t, smul_comm b₂ t]
    linear_combination (norm := abel) -hk
  case refine_2 =>
    intro a
    show t • twistCochain fi y d (algebraMap A B a) - h (algebraMap A B a) = (t * c) • d a
    have hz : a • (t • twistCochain fi y d 1) - a • h 1 = 0 := by
      rw [← smul_sub, hone, smul_zero]
    rw [Algebra.algebraMap_eq_smul_one, twistCochain_smul fi y d c hfy,
      h.map_smul, smul_add, smul_comm t a,
      ← Algebra.algebraMap_eq_smul_one, algebraMap_smul]
    simp only [mul_smul]
    linear_combination (norm := abel) hz

end Derivation

/-- 導分の almost 延長から作った `g : Ω[B⁄R] → B ⊗_A Ω[A⁄R]` は
`mapBaseChange` の **almost な左逆**(`g ∘ mapBaseChange = s·id`)。 -/
theorem retraction_of_derivation {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    (s : A)
    (D : Derivation R B (TensorProduct A B Ω[A⁄R]))
    (hD : ∀ a : A, D (algebraMap A B a) = s • ((1 : B) ⊗ₜ[A] (KaehlerDifferential.D R A) a))
    (ξ : TensorProduct A B Ω[A⁄R]) :
    D.liftKaehlerDifferential (KaehlerDifferential.mapBaseChange R A B ξ) = s • ξ := by
  induction ξ using TensorProduct.induction_on with
  | zero => simp
  | tmul b x =>
    rw [KaehlerDifferential.mapBaseChange_tmul, map_smul]
    have hx : ∀ z : Ω[A⁄R],
        D.liftKaehlerDifferential ((KaehlerDifferential.map R R A B) z)
          = s • ((1 : B) ⊗ₜ[A] z) := by
      intro z
      have hspan : z ∈ Submodule.span A (Set.range (KaehlerDifferential.D R A)) := by
        rw [KaehlerDifferential.span_range_derivation]; trivial
      induction hspan using Submodule.span_induction with
      | mem u hu =>
        obtain ⟨a, rfl⟩ := hu
        rw [KaehlerDifferential.map_D, Derivation.liftKaehlerDifferential_comp_D, hD]
      | zero => simp
      | add u v _ _ ihu ihv => rw [map_add, map_add, ihu, ihv, TensorProduct.tmul_add, smul_add]
      | smul a u _ ihu =>
        rw [map_smul, LinearMap.map_smul_of_tower, ihu, TensorProduct.tmul_smul,
          smul_comm, TensorProduct.smul_tmul']
    rw [hx, smul_comm, TensorProduct.smul_tmul', smul_eq_mul, mul_one]
  | add u v ihu ihv => rw [map_add, map_add, ihu, ihv, smul_add]

/-- **`Theorem 2.4(i)` の核側**——`Ω[A⁄R] ⊗_A B → Ω[B⁄R]` の核は `t·c` で
零化される。 -/
theorem thm_2_4_i_kernel {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    {ι : Type u} [Fintype ι] (fi : ι → (B →ₗ[A] A)) (y : ι → B)
    (c : A) (hfy : ∀ b : B, ∑ i, (fi i b) • y i = c • b)
    (t : A) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = t • (1 : B))
    (ξ : TensorProduct A B Ω[A⁄R])
    (hξ : KaehlerDifferential.mapBaseChange R A B ξ = 0) :
    (t * c) • ξ = 0 := by
  obtain ⟨D, hD⟩ := exists_derivation_extension (R := R) (M := TensorProduct A B Ω[A⁄R])
    fi y c hfy t w hw_ann hw_aug
    ((TensorProduct.mk A B Ω[A⁄R] 1).compDer (KaehlerDifferential.D R A))
  have h := retraction_of_derivation (t * c) D (fun a => hD a) ξ
  rw [hξ, map_zero] at h
  exact h.symm

/-- `almost_dual_basis` の出力を `(fi, y)` の形に直す橋渡し。 -/
theorem dualBasis_family {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Finite A B] [Module.Free A B]
    (c : A) (S : Finset (B × B))
    (hS : ∀ b : B, c • b = ∑ i ∈ S, (Algebra.trace A B (b * i.1)) • i.2) :
    ∀ b : B, ∑ i : {x // x ∈ S},
      ((Algebra.trace A B ∘ₗ LinearMap.mulRight A ((i : B × B).1)) b) • ((i : B × B).2)
      = c • b := by
  intro b
  rw [hS b]
  exact Finset.sum_coe_sort S (fun i => (Algebra.trace A B (b * i.1)) • i.2)

/-- **`Theorem 2.4(i)`(Faltings, Ch.I §2)——`Ω_A ⊗_A B → Ω_B` は almost 同型**。

> *(i) The map `Ω_A ⊗_A B → Ω_B` is an almost isomorphism, that is, its
> kernel and cokernel are annihilated by `m`.*(物理 p.8 = 印字 p.261)

核は `p^{2n}`、余核は `p^n` で零化される。`n` は `Definition 2.1` 条件(iii)
の witness の水準なので、いくらでも小さく取れる——これが「`m` が零化する」
の内容である。 -/
theorem thm_2_4_i {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    [Module.Finite A B] [Module.Free A B]
    (p : A) (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p))) :
    (∀ ξ : TensorProduct A B Ω[A⁄R],
        KaehlerDifferential.mapBaseChange R A B ξ = 0 → ((p^n) * (p^n)) • ξ = 0) ∧
    (∀ x : Ω[B⁄R], (algebraMap A B (p^n)) • x
        ∈ LinearMap.range (KaehlerDifferential.mapBaseChange R A B)) := by
  have hann := almost_swap_annihilate p hAE hf0inj (p ^ n) w hw
  have haug := almost_swap_augment p hAE hf0inj (p ^ n) w hw
  obtain ⟨S, hS⟩ := almost_dual_basis p hAE hf0inj (p ^ n) w hw
  refine ⟨fun ξ hξ => ?_, fun x => thm_2_4_i_cokernel (p ^ n) w hann haug x⟩
  exact thm_2_4_i_kernel (R := R)
    (fun i : {x // x ∈ S} => Algebra.trace A B ∘ₗ LinearMap.mulRight A ((i : B × B).1))
    (fun i : {x // x ∈ S} => ((i : B × B).2)) (p ^ n)
    (dualBasis_family (p ^ n) S hS) (p ^ n) w hann haug ξ hξ

/-! ## 非空虚性の対照

`R = ℤ`、`A = ℤ[X]`(**`Ω[A⁄R] ≠ 0`** なので主張が退化しない)、
`B = Fin 2 → A`(階数 2 の非自明な有限エタール拡大)、`p = 5`(真の非単元)。 -/

theorem thm_2_4_i_nonvacuous_hf0inj :
    Function.Injective (algebraMap (Fin 2 → Polynomial ℤ)
      (Localization.Away (algebraMap (Polynomial ℤ) (Fin 2 → Polynomial ℤ)
        (5 : Polynomial ℤ)))) := by
  refine IsLocalization.injective
    (M := Submonoid.powers
      (algebraMap (Polynomial ℤ) (Fin 2 → Polynomial ℤ) (5 : Polynomial ℤ))) _ ?_
  rw [Submonoid.powers_le, mem_nonZeroDivisors_iff]
  have key : ∀ z : Fin 2 → Polynomial ℤ,
      (algebraMap (Polynomial ℤ) (Fin 2 → Polynomial ℤ)) 5 * z = 0 → z = 0 := by
    intro z hz
    funext j
    have hj := congrFun hz j
    simp only [Pi.mul_apply, Pi.zero_apply] at hj
    rw [show (algebraMap (Polynomial ℤ) (Fin 2 → Polynomial ℤ)) 5 j = 5 by simp] at hj
    rcases mul_eq_zero.mp hj with h | h
    · exact absurd h (by norm_num)
    · exact h
  exact ⟨key, fun z hz => key z (by rw [mul_comm]; exact hz)⟩

/-- **`thm_2_4_i` の非空虚性**——仮定を全て満たす具体例が実在する。 -/
example :
    (∀ ξ : TensorProduct (Polynomial ℤ) (Fin 2 → Polynomial ℤ) Ω[Polynomial ℤ⁄ℤ],
        KaehlerDifferential.mapBaseChange ℤ (Polynomial ℤ) (Fin 2 → Polynomial ℤ) ξ = 0 →
        (((5:Polynomial ℤ)^1) * ((5:Polynomial ℤ)^1)) • ξ = 0) ∧
    (∀ x : Ω[(Fin 2 → Polynomial ℤ)⁄ℤ],
        (algebraMap (Polynomial ℤ) (Fin 2 → Polynomial ℤ) ((5:Polynomial ℤ)^1)) • x
        ∈ LinearMap.range
          (KaehlerDifferential.mapBaseChange ℤ (Polynomial ℤ) (Fin 2 → Polynomial ℤ))) := by
  have hAE : IsAlmostEtaleCovering (A := Polynomial ℤ) (B := Fin 2 → Polynomial ℤ)
      (5 : Polynomial ℤ) := isAlmostEtaleCovering_of_etale_general _
  refine thm_2_4_i (R := ℤ) (5 : Polynomial ℤ) hAE thm_2_4_i_nonvacuous_hf0inj 1
    ((5 : Polynomial ℤ)^1 • Algebra.FormallyUnramified.elem
      (Polynomial ℤ) (Fin 2 → Polynomial ℤ)) ?_
  letI := awayAlgebra (5 : Polynomial ℤ) (A := Polynomial ℤ) (B := Fin 2 → Polynomial ℤ)
  haveI := hAE.2.2.1
  haveI := (hAE.2.1 : Module.Finite _ _)
  rw [map_smul, diagonalCompare_elem_eq]

/-! ## 原文どおりの「`m` が零化する」形へ

`thm_2_4_i` は witness の水準 `p^n` で述べているが、原文は
*"its kernel and cokernel are annihilated by **`m`**"* と書く。
`p`-可除な塔 `PDivTower` の上では条件(iii)の witness が**各水準 `ϖ k` で
取れる**ので、その形にそのまま強められる。

核側で `ϖ k` そのもの(`ϖ k` の平方ではなく)が効くのは、`q ≥ 2` のとき
`ϖ k = (ϖ(k+1))^q` が `(ϖ(k+1))²` で割れるから——塔があることの
ご利益がここでも効く。 -/

open KaehlerDifferential in
/-- `thm_2_4_i` の一般スカラー版(witness の水準を `c` に一般化)。 -/
theorem thm_2_4_i_gen {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    [Module.Finite A B] [Module.Free A B]
    (p : A) (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (c : A) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = c • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p))) :
    (∀ ξ : TensorProduct A B Ω[A⁄R],
        KaehlerDifferential.mapBaseChange R A B ξ = 0 → (c * c) • ξ = 0) ∧
    (∀ x : Ω[B⁄R], (algebraMap A B c) • x
        ∈ LinearMap.range (KaehlerDifferential.mapBaseChange R A B)) := by
  have hann := almost_swap_annihilate p hAE hf0inj c w hw
  have haug := almost_swap_augment p hAE hf0inj c w hw
  obtain ⟨S, hS⟩ := almost_dual_basis p hAE hf0inj c w hw
  refine ⟨fun ξ hξ => ?_, fun x => thm_2_4_i_cokernel c w hann haug x⟩
  exact thm_2_4_i_kernel (R := R)
    (fun i : {x // x ∈ S} => Algebra.trace A B ∘ₗ LinearMap.mulRight A ((i : B × B).1))
    (fun i : {x // x ∈ S} => ((i : B × B).2)) c
    (dualBasis_family c S hS) c w hann haug ξ hξ

open KaehlerDifferential in
/-- **`Theorem 2.4(i)`、塔の上での形**——核も余核も**各水準 `ϖ k` で**
零化される。 -/
theorem thm_2_4_i_tower {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    [Module.Finite A B] [Module.Free A B]
    {q : ℕ} (hq : 2 ≤ q) (T : PDivTower A q)
    (hAET : IsAlmostEtaleCoveringTower (A := A) (B := B) T)
    (hf0inj : letI := awayAlgebra (T.ϖ 0) (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B (T.ϖ 0))))) :
    (∀ (k : ℕ) (ξ : TensorProduct A B Ω[A⁄R]),
        KaehlerDifferential.mapBaseChange R A B ξ = 0 → T.ϖ k • ξ = 0) ∧
    (∀ (k : ℕ) (x : Ω[B⁄R]), (algebraMap A B (T.ϖ k)) • x
        ∈ LinearMap.range (KaehlerDifferential.mapBaseChange R A B)) := by
  have hAE := isAlmostEtaleCovering_of_tower T hAET
  obtain ⟨_, _, _, _, hwit⟩ := hAET
  constructor
  · intro k ξ hξ
    obtain ⟨wk, hwk⟩ := hwit (k+1)
    have hker := (thm_2_4_i_gen (R := R) (T.ϖ 0) hAE hf0inj (T.ϖ (k+1)) wk hwk).1 ξ hξ
    have hdvd : T.ϖ k = (T.ϖ (k+1))^(q - 2) * (T.ϖ (k+1) * T.ϖ (k+1)) := by
      rw [← T.ϖ_succ k, ← pow_two, ← pow_add]
      congr 1
      omega
    rw [hdvd, ← smul_smul, hker, smul_zero]
  · intro k x
    obtain ⟨wk, hwk⟩ := hwit k
    exact (thm_2_4_i_gen (R := R) (T.ϖ 0) hAE hf0inj (T.ϖ k) wk hwk).2 x

open KaehlerDifferential in
/-- **`Theorem 2.4(i)`——「核は `m` で零化される」**(原文の言葉そのもの)。 -/
theorem thm_2_4_i_m_annihilates {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    [Module.Finite A B] [Module.Free A B]
    {q : ℕ} (hq : 2 ≤ q) (T : PDivTower A q)
    (hAET : IsAlmostEtaleCoveringTower (A := A) (B := B) T)
    (hf0inj : letI := awayAlgebra (T.ϖ 0) (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B (T.ϖ 0))))) :
    T.IsAlmostZero (LinearMap.ker (KaehlerDifferential.mapBaseChange R A B))
    ∧ (∀ (a : A), a ∈ T.m → ∀ ξ : TensorProduct A B Ω[A⁄R],
        KaehlerDifferential.mapBaseChange R A B ξ = 0 → a • ξ = 0) := by
  obtain ⟨hker, hcoker⟩ := thm_2_4_i_tower (R := R) hq T hAET hf0inj
  have hz : T.IsAlmostZero (LinearMap.ker (KaehlerDifferential.mapBaseChange R A B)) := by
    intro k x
    apply Subtype.ext
    show T.ϖ k • (x : TensorProduct A B Ω[A⁄R]) = 0
    exact hker k x x.2
  refine ⟨hz, fun a ha ξ hξ => ?_⟩
  have h := (T.isAlmostZero_iff_m_smul_eq_bot
    (LinearMap.ker (KaehlerDifferential.mapBaseChange R A B))).mp hz ⟨ξ, hξ⟩ a ha
  exact congrArg Subtype.val h

/-! ## 「almost 同型」の言葉

Faltings は §3・§4 で *"the map ... is an **almost isomorphism**"* を
繰り返し使う(`Theorem 2.4(i)`・`4.1(i)`・`4.2(ii)(iii)(vi)`・`4.5(ii)` 等)。
その形を定義として切り出し、**合成則**まで用意しておく——水準は積になる。 -/

/-- **水準 `c` での almost 同型**——核が `c` で零化され、余核も `c` で
零化される(= `c·N ⊆ im f`)。 -/
def AlmostIsoAt {A M N : Type u} [CommRing A] [AddCommGroup M] [Module A M]
    [AddCommGroup N] [Module A N] (c : A) (f : M →ₗ[A] N) : Prop :=
  (∀ x : M, f x = 0 → c • x = 0) ∧ (∀ y : N, c • y ∈ LinearMap.range f)

theorem almostIsoAt_id {A M : Type u} [CommRing A] [AddCommGroup M] [Module A M] :
    AlmostIsoAt (1 : A) (LinearMap.id : M →ₗ[A] M) :=
  ⟨fun x hx => by rw [one_smul]; exact hx, fun y => ⟨y, by simp⟩⟩

/-- **almost 同型の合成**——水準は積になる。 -/
theorem almostIsoAt_comp {A M N P : Type u} [CommRing A] [AddCommGroup M] [Module A M]
    [AddCommGroup N] [Module A N] [AddCommGroup P] [Module A P]
    (c d : A) (f : M →ₗ[A] N) (g : N →ₗ[A] P)
    (hf : AlmostIsoAt c f) (hg : AlmostIsoAt d g) :
    AlmostIsoAt (c * d) (g ∘ₗ f) := by
  refine ⟨fun x hx => ?_, fun z => ?_⟩
  · have h1 : d • f x = 0 := hg.1 (f x) hx
    have h2 : f (d • x) = 0 := by rw [map_smul]; exact h1
    have h3 : c • (d • x) = 0 := hf.1 _ h2
    rw [smul_smul] at h3
    exact h3
  · obtain ⟨y, hy⟩ := hg.2 z
    obtain ⟨w, hw⟩ := hf.2 y
    refine ⟨w, ?_⟩
    show g (f w) = (c * d) • z
    rw [hw, map_smul, hy, smul_smul, mul_comm]

/-- 水準は割り切れる方へ緩められる。 -/
theorem almostIsoAt_of_dvd {A M N : Type u} [CommRing A] [AddCommGroup M] [Module A M]
    [AddCommGroup N] [Module A N] {c c' : A} (h : c ∣ c') (f : M →ₗ[A] N)
    (hf : AlmostIsoAt c f) : AlmostIsoAt c' f := by
  obtain ⟨e, rfl⟩ := h
  refine ⟨fun x hx => ?_, fun y => ?_⟩
  · rw [mul_comm, ← smul_smul, hf.1 x hx, smul_zero]
  · rw [mul_comm, ← smul_smul]
    exact Submodule.smul_mem _ e (hf.2 y)

open KaehlerDifferential in
/-- **`Theorem 2.4(i)` を「almost 同型」の言葉で**——原文の
*"The map `Ω_A ⊗_A B → Ω_B` is an **almost isomorphism**"* そのもの。 -/
theorem thm_2_4_i_almostIso {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    [Module.Finite A B] [Module.Free A B]
    (p : A) (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p))) :
    AlmostIsoAt ((p^n) * (p^n))
      ((KaehlerDifferential.mapBaseChange R A B).restrictScalars A) := by
  obtain ⟨hker, hcoker⟩ := thm_2_4_i (R := R) p hAE hf0inj n w hw
  refine ⟨fun x hx => hker x hx, fun y => ?_⟩
  have h1 := hcoker y
  rw [algebraMap_smul] at h1
  obtain ⟨z, hz⟩ := h1
  refine ⟨(p^n) • z, ?_⟩
  show KaehlerDifferential.mapBaseChange R A B ((p^n) • z) = ((p^n) * (p^n)) • y
  rw [(KaehlerDifferential.mapBaseChange R A B).map_smul_of_tower, hz, smul_smul]

/-! ### 塔の上での almost 同型——**合成で閉じる**

各水準 `ϖ k` で almost 同型であることが Faltings の
「核と余核が `m` で零化される」の内容である。塔があるおかげで
`ϖ(k+1)² ∣ ϖ k`(`q ≥ 2`)が使え、**合成しても水準が回収できる**
——これが almost mathematics が「圏」として機能する理由である。 -/

/-- **塔の上での almost 同型**——各水準 `ϖ k` で almost 同型。 -/
def AlmostIso {A M N : Type u} [CommRing A] [AddCommGroup M] [Module A M]
    [AddCommGroup N] [Module A N] {q : ℕ} (T : PDivTower A q) (f : M →ₗ[A] N) : Prop :=
  ∀ k : ℕ, AlmostIsoAt (T.ϖ k) f

/-- **almost 同型は合成で閉じる**。 -/
theorem almostIso_comp {A M N P : Type u} [CommRing A] [AddCommGroup M] [Module A M]
    [AddCommGroup N] [Module A N] [AddCommGroup P] [Module A P]
    {q : ℕ} (hq : 2 ≤ q) (T : PDivTower A q) (f : M →ₗ[A] N) (g : N →ₗ[A] P)
    (hf : AlmostIso T f) (hg : AlmostIso T g) : AlmostIso T (g ∘ₗ f) := by
  intro k
  have hdvd : (T.ϖ (k+1) * T.ϖ (k+1)) ∣ T.ϖ k := by
    refine ⟨(T.ϖ (k+1))^(q-2), ?_⟩
    rw [← T.ϖ_succ k, ← pow_two, ← pow_add]
    congr 1
    omega
  exact almostIsoAt_of_dvd hdvd _ (almostIsoAt_comp _ _ f g (hf (k+1)) (hg (k+1)))

/-- almost 同型で源に `m`-捩れが無ければ**単射**——`Theorem 4.1(ii)` の
証明の形(*"its kernel is annihilated by `m`. But ... has no `m`-torsion"*)。 -/
theorem injective_of_almostIso_of_noTorsion {A M N : Type u} [CommRing A]
    [AddCommGroup M] [Module A M] [AddCommGroup N] [Module A N]
    {q : ℕ} (T : PDivTower A q) (f : M →ₗ[A] N) (hf : AlmostIso T f)
    (hnotors : ∀ (k : ℕ) (x : M), T.ϖ k • x = 0 → x = 0) :
    Function.Injective f := by
  rw [injective_iff_map_eq_zero]
  intro x hx
  exact hnotors 0 x ((hf 0).1 x hx)

end ABC3.Found.Falt1
