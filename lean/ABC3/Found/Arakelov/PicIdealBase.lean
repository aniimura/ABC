import Mathlib.RingTheory.PicardGroup
import Mathlib.RingTheory.Flat.Basic
import Mathlib.LinearAlgebra.TensorProduct.Tower
import ABC3.Meta.Claim

/-!
# Arakelov (B2) 第 183 ブロック —— **平坦底変換で可逆イデアルは可逆**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`isCartierDivisor_affine` の ⟸ 向きに要る唯一の欠落

§9-203 で測った通り、`Spec R` の上で

    Module.Invertible Γ(⊤) (D.ideal ⊤)  ⟹  ∀ A アフィン開, Module.Invertible Γ(A) (D.ideal A)

を出すのに要る部品は 5 個で、**4 個は mathlib に在り**、欠落は本ブロックの 1 個だけだった:

    S ⊗[R] I  ≃ₗ[S]  I.map (algebraMap R S)     (S が R 上平坦なら)

## ★★機構

| 段 | 内容 |
|---|---|
| 1 | `bcMap : S ⊗[R] I → S`(`s ⊗ i ↦ s · φ(i)`) |
| 2 | ★**単射**——`S` が平坦だから `lTensor` が単射性を保つ |
| 3 | ★**像は `I.map φ`**——生成元で両包含 |
| 4 | `Module.Invertible S (S ⊗[R] I)`(mathlib の底変換 instance)を移す |

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `bcMap` / `bcMap_tmul` | ★底変換の写像 |
| `bcMap_injective` | ★★平坦なら単射 |
| `bcMap_range` | ★★像はイデアルの拡大 |
| `idealBaseChangeEquiv` | ★★★`S ⊗[R] I ≃ₗ[S] I.map φ` |
| `invertible_map_of_flat` | ★★★★**平坦底変換で可逆イデアルは可逆** |
-/

universe u

namespace ABC3.Found.Arakelov

open TensorProduct

variable (R S : Type u) [CommRing R] [CommRing S] [Algebra R S] (I : Ideal R)

/-- ★**底変換の写像** `S ⊗[R] I → S`。 -/
noncomputable def bcMap : (S ⊗[R] I) →ₗ[S] S :=
  (TensorProduct.AlgebraTensorModule.rid R S S).toLinearMap ∘ₗ
    LinearMap.baseChange S I.subtype

variable {R S I}

@[simp] theorem bcMap_tmul (s : S) (i : I) :
    bcMap R S I (s ⊗ₜ[R] i) = s * algebraMap R S (i : R) := by
  simp [bcMap, Algebra.smul_def, mul_comm]

/-- ★★**平坦なら底変換の写像は単射である**。 -/
theorem bcMap_injective [Module.Flat R S] : Function.Injective (bcMap R S I) := by
  have h : Function.Injective (LinearMap.lTensor S (I.subtype : I →ₗ[R] R)) :=
    Module.Flat.lTensor_preserves_injective_linearMap _ Subtype.val_injective
  intro a b hab
  refine h ?_
  have key : ∀ t : S ⊗[R] I, LinearMap.lTensor S (I.subtype : I →ₗ[R] R) t
      = (TensorProduct.AlgebraTensorModule.rid R S S).symm (bcMap R S I t) := by
    intro t
    simp [bcMap, LinearMap.baseChange_eq_ltensor]
  rw [key a, key b, hab]

/-- ★★**底変換の写像の像はイデアルの拡大である**。 -/
theorem bcMap_range :
    LinearMap.range (bcMap R S I) = (I.map (algebraMap R S) : Ideal S) := by
  apply le_antisymm
  · rintro _ ⟨t, rfl⟩
    induction t using TensorProduct.induction_on with
    | zero => simp
    | tmul s i =>
        rw [bcMap_tmul]
        exact Ideal.mul_mem_left _ _ (Ideal.mem_map_of_mem _ i.2)
    | add a b ha hb =>
        rw [map_add]
        exact Submodule.add_mem _ ha hb
  · rw [Ideal.map_le_iff_le_comap]
    intro i hi
    exact ⟨(1 : S) ⊗ₜ[R] (⟨i, hi⟩ : I), by rw [bcMap_tmul]; simp⟩

/-- ★★★**平坦なら底変換はイデアルの拡大と一致する**。 -/
noncomputable def idealBaseChangeEquiv [Module.Flat R S] :
    (S ⊗[R] I) ≃ₗ[S] (I.map (algebraMap R S) : Ideal S) :=
  (LinearEquiv.ofInjective (bcMap R S I) bcMap_injective) ≪≫ₗ
    LinearEquiv.ofEq _ _ bcMap_range

/-- ★★★★**平坦底変換で可逆イデアルは可逆イデアルに移る**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `isCartierDivisor_affine` の ⟸ 向きの核である
——アフィン開の制限は平坦(`Scheme.Hom.flat_appLE`)だから。 -/
theorem invertible_map_of_flat [Module.Flat R S] [Module.Invertible R I] :
    Module.Invertible S (I.map (algebraMap R S) : Ideal S) :=
  Module.Invertible.congr idealBaseChangeEquiv

/-! ## ★出典の紐付け(`.src`) -/

def invertible_map_of_flat.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——平坦底変換で可逆イデアルは可逆)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
