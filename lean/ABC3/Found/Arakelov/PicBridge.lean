import ABC3.Found.Arakelov.PicConstSection

/-!
# Arakelov (B1) 第 88 ブロック —— **橋渡しの写像**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★型の運搬(§9-90 で旗を立てたもの)

`tensorSectionMap` の定義域は **`𝒪(U)` 上のテンソル**、
第 87 の恒等式は **`R` 上の元**(`m ⊗ₜ n`)についてである。

★★両者を繋ぐ **`R` 線型写像**

    M ⊗_R N  ⟶  (tilde M)(U) ⊗_{𝒪(U)} (tilde N)(U)
    m ⊗ₜ n   ↦   const m ⊗ₜ const n

を作るのが本ブロックである。

## ★★機構

★`R` の作用は `IsScalarTower R 𝒪(U) (tilde M)(U)`(mathlib の instance)を通す。
★★双線型性は `algebraMap_smul` で `𝒪(U)` の作用に直してから
`smul_tmul'` / `tmul_smul` を当てる——第 80 ブロックと同じ形。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `toOpenTensor` | ★★★★**橋渡しの写像** |
| `tensorSectionMap_comp_toOpenTensor` | ★★★★**合成は定数切断そのもの** |

## ★★★これで第 86 の器具が使える

    tensorSectionMap ∘ toOpenTensor = toOpenₗ (M ⊗_R N)

★両辺とも `M ⊗_R N` からの写像で、`U = D(f)` なら**どちらも局所化**である。
★★したがって `tensorLocalization_ext` で `tensorSectionMap` が同型と決まる。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace TensorProduct StructureSheaf

variable (R : Type u) [CommRing R] (M N : Type u)
  [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
  (U : Opens (PrimeSpectrum.Top R))

/-- ★★★★**橋渡しの写像** `M ⊗_R N ⟶ (tilde M)(U) ⊗_{𝒪(U)} (tilde N)(U)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`m ⊗ₜ n` を「定数切断のテンソル」へ送る。 -/
noncomputable def toOpenTensor :
    M ⊗[R] N →ₗ[R]
      TensorProduct ((structureSheafInType R R).1.obj (op U) : Type u)
        ((structureSheafInType R M).1.obj (op U) : Type u)
        ((structureSheafInType R N).1.obj (op U) : Type u) :=
  TensorProduct.lift
    { toFun := fun m =>
        { toFun := fun n => toOpenₗ R M U m ⊗ₜ toOpenₗ R N U n
          map_add' := fun n n' => by rw [map_add, TensorProduct.tmul_add]
          map_smul' := fun r n => by
            rw [map_smul, RingHom.id_apply,
              ← IsScalarTower.algebraMap_smul
                ((structureSheafInType R R).1.obj (op U) : Type u) r (toOpenₗ R N U n),
              TensorProduct.tmul_smul, IsScalarTower.algebraMap_smul] }
      map_add' := fun m m' => by
        ext n; simp only [LinearMap.coe_mk, AddHom.coe_mk, LinearMap.add_apply]
        rw [map_add, TensorProduct.add_tmul]
      map_smul' := fun r m => by
        ext n; simp only [LinearMap.coe_mk, AddHom.coe_mk, RingHom.id_apply,
          LinearMap.smul_apply]
        rw [map_smul, ← IsScalarTower.algebraMap_smul
            ((structureSheafInType R R).1.obj (op U) : Type u) r (toOpenₗ R M U m),
          TensorProduct.smul_tmul', IsScalarTower.algebraMap_smul] }

@[simp] theorem toOpenTensor_tmul (m : M) (n : N) :
    toOpenTensor R M N U (m ⊗ₜ[R] n) = toOpenₗ R M U m ⊗ₜ toOpenₗ R N U n := rfl

/-- ★★★★**合成は定数切断そのものである**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで第 86 ブロックの `tensorLocalization_ext` が使えるようになる。 -/
theorem tensorSectionMap_comp_toOpenTensor :
    (tensorSectionMap R M N U).restrictScalars R ∘ₗ toOpenTensor R M N U
      = toOpenₗ R (M ⊗[R] N) U := by
  refine TensorProduct.ext' fun m n => ?_
  exact tensorSection_toOpenL R M N U m n

/-! ## ★出典の紐付け(`.src`) -/

def toOpenTensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——橋渡しの写像)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
