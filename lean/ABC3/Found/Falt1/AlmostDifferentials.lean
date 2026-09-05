import ABC3.Found.Falt1.AlmostEtale
import Mathlib.RingTheory.Kaehler.JacobiZariski

/-!
# [Falt1] Theorem 2.4(i) —— `Ω_A⊗_AB → Ω_B` の almost 同型性、余核側(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §2、物理 p.8
(印字 p.261)。

内容 (Falt1 p.8、260dpi 目視): Suppose B is an almost étale covering of A.
(i) The map Ω_A⊗_A B → Ω_B is an almost isomorphism, that is, its kernel
and cokernel are annihilated by m.

## この節の到達点

**余核側は完全に閉じた**。Jacobi–Zariski 完全列
`B⊗_AΩ[A⁄R] --mapBaseChange--> Ω[B⁄R] --map--> Ω[B⁄A] --> 0`
(mathlib の `KaehlerDifferential.range_mapBaseChange`)により、余核は
`Ω[B⁄A]` そのものなので、**`p^n` が `Ω[B⁄A]` を零化すること**を示せばよい。
これは `Definition 2.1` 条件(iii)の witness `w` から直接出る:

`Ω[B⁄A] = I/I²`(`I := ker(μ : B⊗_AB → B)`、mathlib の
`KaehlerDifferential.ideal` はまさにこれ)であり、`x∈I` について
`v := w - algebraMap A (B⊗_AB) (p^n)` と置くと `μ(v)=0` すなわち `v∈I`
(`almost_swap_augment`)、かつ `x·w = 0`(`almost_swap_annihilate`)なので

    p^n·x = (w - v)·x = w·x - v·x = -(v·x) ∈ I·I = I².

古典的な「分離的なら `e` が冪等元で `x = x·e + x·(1-e) = x·(1-e) ∈ I²`」
という1行の議論の、almost 版そのものである。

## 残っている部分(正直な記録)

**核側**は `Algebra.H1Cotangent.exact_δ_mapBaseChange`
(mathlib、`ker(mapBaseChange) = range(δ)`)により
「`p^n` が `Algebra.H1Cotangent A B` を零化する」に帰着するが、これは
**almost formal smoothness**(余接複体 `L_{B/A}` の almost 消滅)であって、
`Ω` の消滅より一段深い。mathlib は `Algebra.FormallySmooth` から
`H1Cotangent` の消滅を出す道具は持つが、その almost 版は無い。
Faltings 自身も「証明はいつも通り Hochschild cohomology を使う」の一文で
畳んでおり、`Theorem 2.2`(nilpotent ideal に沿った lifting)と同じ
議論を要する見込みである。
-/

namespace ABC3.Found.Falt1

universe u

/-- **鍵の計算**: `p^n·(1⊗b - b⊗1) ∈ I²`(`I = ker μ`)。
`w` を条件(iii)の witness とし、`v := w - p^n` と置くと `v∈I` かつ
`(1⊗b-b⊗1)·w = 0` なので `p^n·x = -(v·x) ∈ I²`。 -/
theorem kaehler_key {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (p : A) (n : ℕ) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = p ^ n • (1 : B)) (b : B) :
    (algebraMap A (TensorProduct A B B) (p ^ n)) * ((1:B) ⊗ₜ[A] b - b ⊗ₜ[A] (1:B))
      ∈ (KaehlerDifferential.ideal A B) ^ 2 := by
  set I := KaehlerDifferential.ideal A B with hI
  set x : TensorProduct A B B := (1:B) ⊗ₜ[A] b - b ⊗ₜ[A] (1:B) with hx
  set v : TensorProduct A B B := w - algebraMap A (TensorProduct A B B) (p ^ n) with hv
  have hxI : x ∈ I := by
    rw [hI, hx]
    show Algebra.TensorProduct.lmul' A ((1:B) ⊗ₜ[A] b - b ⊗ₜ[A] (1:B)) = 0
    rw [map_sub, Algebra.TensorProduct.lmul'_apply_tmul, Algebra.TensorProduct.lmul'_apply_tmul]
    ring
  have hvI : v ∈ I := by
    rw [hI, hv]
    show Algebra.TensorProduct.lmul' A (w - algebraMap A (TensorProduct A B B) (p ^ n)) = 0
    rw [map_sub, hw_aug, AlgHom.commutes, Algebra.smul_def, mul_one, sub_self]
  have hkey : (algebraMap A (TensorProduct A B B) (p ^ n)) * x = -(v * x) := by
    have h1 : w * x = 0 := by rw [mul_comm]; exact hw_ann b
    rw [hv, sub_mul, h1]
    ring
  rw [hkey, sq]
  exact neg_mem (Ideal.mul_mem_mul hvI hxI)

/-- **`p^n` は `Ω[B⁄A]` を零化する**(almost unramified 性)。
`Ω[B⁄A]` は `D b` たちで `B`-生成される(`span_range_derivation`)ので、
生成元について `kaehler_key` + `Ideal.toCotangent_eq_zero` で示せばよい。 -/
theorem kaehler_almost_zero {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    (p : A) (n : ℕ) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = p ^ n • (1 : B))
    (x : Ω[B⁄A]) : (p ^ n) • x = 0 := by
  have hbridge : ∀ (z : Ω[B⁄A]), (algebraMap A B (p ^ n)) • z
      = (algebraMap A (TensorProduct A B B) (p ^ n)) • z := by
    intro z; rw [algebraMap_smul, algebraMap_smul]
  have hgen : ∀ b : B, (algebraMap A B (p ^ n)) • (KaehlerDifferential.D A B b) = 0 := by
    intro b
    have hmem : (1:B) ⊗ₜ[A] b - b ⊗ₜ[A] (1:B) ∈ KaehlerDifferential.ideal A B := by
      show Algebra.TensorProduct.lmul' A ((1:B) ⊗ₜ[A] b - b ⊗ₜ[A] (1:B)) = 0
      rw [map_sub, Algebra.TensorProduct.lmul'_apply_tmul,
        Algebra.TensorProduct.lmul'_apply_tmul]
      ring
    have hzero : (algebraMap A (TensorProduct A B B) (p ^ n)) •
        ((KaehlerDifferential.ideal A B).toCotangent
          (⟨(1:B) ⊗ₜ[A] b - b ⊗ₜ[A] (1:B), hmem⟩ : KaehlerDifferential.ideal A B)) = 0 := by
      rw [← map_smul, Ideal.toCotangent_eq_zero]
      show (algebraMap A (TensorProduct A B B) (p ^ n)) • ((1:B) ⊗ₜ[A] b - b ⊗ₜ[A] (1:B))
        ∈ (KaehlerDifferential.ideal A B) ^ 2
      rw [Algebra.smul_def]
      exact kaehler_key p n w hw_ann hw_aug b
    rw [hbridge, KaehlerDifferential.D_apply]
    exact hzero
  have hsub : Submodule.span B (Set.range (KaehlerDifferential.D A B))
      ≤ LinearMap.ker (LinearMap.lsmul B (Ω[B⁄A]) (algebraMap A B (p ^ n))) := by
    rw [Submodule.span_le]
    rintro _ ⟨b, rfl⟩
    exact hgen b
  rw [KaehlerDifferential.span_range_derivation] at hsub
  have hx := hsub (Submodule.mem_top (x := x))
  rw [← algebraMap_smul B (p ^ n) x]
  exact hx

/-- `Definition 2.1` の仮定から直接述べた形(`almost_swap_annihilate`・
`almost_swap_augment` を噛ませただけ)。 -/
theorem kaehler_almost_zero_of_isAlmostEtale {A B : Type u} [CommRing A] [CommRing B] [Algebra A B]
    [Module.Free A B] (p : A)
    (hAE : IsAlmostEtaleCovering (A := A) (B := B) p)
    (hf0inj : letI := awayAlgebra p (A := A) (B := B)
      Function.Injective (algebraMap B (Localization.Away (algebraMap A B p))))
    (n : ℕ) (w : TensorProduct A B B)
    (hw : letI := awayAlgebra p (A := A) (B := B)
      haveI := hAE.2.2.1
      haveI := (hAE.2.1 : Module.Finite _ _)
      diagonalCompare p w
        = p ^ n • Algebra.FormallyUnramified.elem (Localization.Away p)
            (Localization.Away (algebraMap A B p)))
    (x : Ω[B⁄A]) : (p ^ n) • x = 0 :=
  kaehler_almost_zero p n w
    (fun q => almost_swap_annihilate p hAE hf0inj n w hw q)
    (almost_swap_augment p hAE hf0inj n w hw) x

/-- **`Theorem 2.4(i)` の余核側**: `Ω[A⁄R]⊗_AB → Ω[B⁄R]` の余核は `p^n` で
零化される。すなわち任意の `x : Ω[B⁄R]` について `p^n·x` は像に入る。
Jacobi–Zariski 完全列(`KaehlerDifferential.range_mapBaseChange`——
像がちょうど `Ω[B⁄R] → Ω[B⁄A]` の核)と `kaehler_almost_zero` から。 -/
theorem thm_2_4_i_cokernel {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    (p : A) (n : ℕ) (w : TensorProduct A B B)
    (hw_ann : ∀ q : B, (1 ⊗ₜ[A] q - q ⊗ₜ[A] 1) * w = 0)
    (hw_aug : Algebra.TensorProduct.lmul' A w = p ^ n • (1 : B))
    (x : Ω[B⁄R]) :
    (algebraMap A B (p ^ n)) • x ∈ LinearMap.range (KaehlerDifferential.mapBaseChange R A B) := by
  rw [KaehlerDifferential.range_mapBaseChange, LinearMap.mem_ker, map_smul]
  have h := kaehler_almost_zero p n w hw_ann hw_aug (KaehlerDifferential.map R A B B x)
  rw [← algebraMap_smul B (p ^ n) (KaehlerDifferential.map R A B B x)] at h
  exact h

/-- **`Theorem 2.4(i)` の核側は `H1Cotangent` の almost 消滅にちょうど
帰着する**。`Algebra.H1Cotangent.exact_δ_mapBaseChange`(mathlib、
`ker(mapBaseChange) = range(δ)`)により、`p^n` が `H1Cotangent A B` を
零化しさえすれば核側が出る——逆に言えば、**残っているのはその 1 本だけ**
である(= `B` が `A` 上 almost formally smooth であること)。 -/
theorem thm_2_4_i_kernel_of_h1Cotangent {R A B : Type u} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    (p : A) (n : ℕ)
    (hH1 : ∀ y : Algebra.H1Cotangent A B, (algebraMap A B (p ^ n)) • y = 0)
    (y : TensorProduct A B (Ω[A⁄R]))
    (hy : KaehlerDifferential.mapBaseChange R A B y = 0) :
    (algebraMap A B (p ^ n)) • y = 0 := by
  obtain ⟨z, rfl⟩ := (Algebra.H1Cotangent.exact_δ_mapBaseChange R A B) y |>.mp hy
  rw [← map_smul, hH1 z, map_zero]

end ABC3.Found.Falt1
