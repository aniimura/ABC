import ABC3.Found.Arakelov.PicIdealFlat

/-!
# Arakelov (B2) 第 185 ブロック —— **可逆イデアルは積で閉じる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 183 と**同じ形**で出た

    I ⊗[R] J  --(1 ⊗ J↪R)-->  I ⊗[R] R  ≅  I  ⊆  R

★★この合成の**像は `I * J`**、**単射性は `I` の平坦性**から出る
(可逆 ⟹ 射影的 ⟹ 平坦、mathlib の instance)。
したがって `I ⊗ J ≅ I * J` であり、`Module.Invertible R (I ⊗ J)` を移せばよい。

★★★mathlib の `Submodule.tensorEquivMul` は `(Submodule R A)ˣ`(**単元**)を要求するので
イデアルには直に当たらない——`Module.Invertible R ↥I` は
「`I` が `Submodule R R` の単元」を**意味しない**からである(可逆イデアルの逆は
`R` の中には無い)。★そこで平坦性の筋を採った。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `mulTensorMap` | ★掛け算の写像 `I ⊗ J → R` |
| `mulTensorMap_injective` | ★★`I` が平坦だから単射 |
| `mulTensorMap_range` | ★★像は `I * J` |
| `invertible_mul` | ★★★★**可逆イデアルの積は可逆** |
| `invertible_top` | ★`⊤` は可逆 |
-/

universe u

namespace ABC3.Found.Arakelov

open TensorProduct

variable {R : Type u} [CommRing R] (I J : Ideal R)

/-- ★**掛け算の写像** `I ⊗[R] J → R`。 -/
noncomputable def mulTensorMap : (I ⊗[R] J) →ₗ[R] R :=
  I.subtype ∘ₗ (TensorProduct.rid R I).toLinearMap ∘ₗ LinearMap.lTensor I J.subtype

@[simp] theorem mulTensorMap_tmul (i : I) (j : J) :
    mulTensorMap I J (i ⊗ₜ[R] j) = (i : R) * (j : R) := by
  simp [mulTensorMap, mul_comm]

variable {I J}

/-- ★★**`I` が可逆(ゆえに平坦)なら掛け算の写像は単射である**。 -/
theorem mulTensorMap_injective [Module.Invertible R I] :
    Function.Injective (mulTensorMap I J) := by
  have h1 : Function.Injective
      (LinearMap.lTensor (I : Type u) (J.subtype : (J : Type u) →ₗ[R] R)) :=
    Module.Flat.lTensor_preserves_injective_linearMap _ Subtype.val_injective
  intro a b hab
  exact h1 ((TensorProduct.rid R (I : Type u)).injective (Subtype.val_injective hab))

/-- ★★**掛け算の写像の像は `I * J` である**。 -/
theorem mulTensorMap_range :
    LinearMap.range (mulTensorMap I J) = I * J := by
  apply le_antisymm
  · rintro _ ⟨t, rfl⟩
    induction t using TensorProduct.induction_on with
    | zero => simp
    | tmul i j =>
        rw [mulTensorMap_tmul]
        exact Ideal.mul_mem_mul i.2 j.2
    | add a b ha hb =>
        rw [map_add]
        exact Submodule.add_mem _ ha hb
  · refine Submodule.mul_le.2 (fun i hi j hj => ?_)
    exact ⟨(⟨i, hi⟩ : I) ⊗ₜ[R] (⟨j, hj⟩ : J), by rw [mulTensorMap_tmul]⟩

/-- ★★★★**可逆イデアルの積は可逆である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `CartierPicData.isCartierDivisor_mul` の中身である。 -/
theorem invertible_mul [Module.Invertible R I] [Module.Invertible R J] :
    Module.Invertible R (I * J : Ideal R) :=
  Module.Invertible.congr
    ((LinearEquiv.ofInjective _ mulTensorMap_injective) ≪≫ₗ
      LinearEquiv.ofEq _ _ mulTensorMap_range)

/-- ★**`⊤` は可逆である**。 -/
theorem invertible_top : Module.Invertible R (⊤ : Ideal R) :=
  Module.Invertible.congr (Submodule.topEquiv (R := R) (M := R)).symm

/-! ## ★出典の紐付け(`.src`) -/

def invertible_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆イデアルは積で閉じる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
