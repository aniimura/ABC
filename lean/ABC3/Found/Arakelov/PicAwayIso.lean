import ABC3.Found.Arakelov.PicBridge

/-!
# Arakelov (B1) 第 89 ブロック —— **基本開集合で比較射は同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★局所化の一意性で押し切る

    tensorSectionMap ∘ toOpenTensor = toOpenₗ (M ⊗_R N)      (第 88 ブロック)

★`U = D(f)` のとき**両辺の入力側はどちらも `M ⊗_R N` の `powers f` での局所化**である:

| 写像 | 局所化であること |
|---|---|
| `toOpenₗ R (M ⊗ N) (D f)` | ★mathlib の instance |
| `toOpenTensor R M N (D f)` | ★★本ブロックで示す(第 79 の連鎖) |

★★★したがって `tensorSectionMap` は**2 つの局所化の間の標準同型**であり、
`linearMap_ext` で一意に決まる。

## ★★機構

`toOpenTensor` は

    M ⊗_R N --map(toOpenₗ, toOpenₗ)--> M' ⊗_R N' --moduleTensorEquiv.symm--> M' ⊗_{𝒪(D f)} N'

と分解する(**純テンソルで `rfl`**、実測)。★前半は mathlib の instance、
後半は線型同値なので `of_linearEquiv` で局所化性が伝わる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isLocalizedModule_toOpenTensor` | ★★`toOpenTensor` は局所化 |
| `tensorSectionMap_bijective` | ★★★★★★**基本開集合で比較射は全単射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace TensorProduct StructureSheaf

variable (R : Type u) [CommRing R] (M N : Type u)
  [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N] (f : R)

/-- ★`D(f)` での構造層の切断。 -/
noncomputable abbrev sheafO :=
  ((structureSheafInType R R).1.obj (op (PrimeSpectrum.basicOpen f)) : Type u)

/-- ★`D(f)` での `tilde M` の切断。 -/
noncomputable abbrev sheafM :=
  ((structureSheafInType R M).1.obj (op (PrimeSpectrum.basicOpen f)) : Type u)

/-- ★★係数を `R` から `𝒪(D f)` へ移すテンソルの同値。 -/
noncomputable abbrev awayTensorEquiv :
    sheafM R M f ⊗[sheafO R f] sheafM R N f ≃ₗ[sheafO R f]
      sheafM R M f ⊗[R] sheafM R N f :=
  IsLocalization.moduleTensorEquiv (Submonoid.powers f) (sheafO R f)
    (sheafM R M f) (sheafM R N f)

/-- ★★**`toOpenTensor` は `powers f` での局所化である**。 -/
theorem isLocalizedModule_toOpenTensor : IsLocalizedModule (Submonoid.powers f)
    (toOpenTensor R M N (PrimeSpectrum.basicOpen f)) := by
  have h : toOpenTensor R M N (PrimeSpectrum.basicOpen f)
      = LinearMap.restrictScalars R (awayTensorEquiv R M N f).symm.toLinearMap
        ∘ₗ TensorProduct.map (toOpenₗ R M (PrimeSpectrum.basicOpen f))
            (toOpenₗ R N (PrimeSpectrum.basicOpen f)) :=
    TensorProduct.ext' fun _ _ => rfl
  rw [h]
  exact IsLocalizedModule.of_linearEquiv _ _
    (LinearEquiv.restrictScalars R (awayTensorEquiv R M N f).symm)

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**基本開集合の上で比較射は全単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが「`tilde` はテンソルを保つ」ことの中核である
——**`M`・`N` の可逆性は要らない**。 -/
theorem tensorSectionMap_bijective :
    Function.Bijective (tensorSectionMap R M N (PrimeSpectrum.basicOpen f)) := by
  haveI := isLocalizedModule_toOpenTensor R M N f
  have heq : LinearMap.restrictScalars R (tensorSectionMap R M N (PrimeSpectrum.basicOpen f))
      = (IsLocalizedModule.linearEquiv (Submonoid.powers f)
          (toOpenTensor R M N (PrimeSpectrum.basicOpen f))
          (toOpenₗ R (M ⊗[R] N) (PrimeSpectrum.basicOpen f))).toLinearMap := by
    refine IsLocalizedModule.linearMap_ext (Submonoid.powers f)
      (toOpenTensor R M N (PrimeSpectrum.basicOpen f))
      (toOpenₗ R (M ⊗[R] N) (PrimeSpectrum.basicOpen f)) ?_
    refine LinearMap.ext fun z => ?_
    exact (LinearMap.congr_fun (tensorSectionMap_comp_toOpenTensor R M N
      (PrimeSpectrum.basicOpen f)) z).trans
      (IsLocalizedModule.linearEquiv_apply _ _ _ z).symm
  have : Function.Bijective
      (LinearMap.restrictScalars R (tensorSectionMap R M N (PrimeSpectrum.basicOpen f))) := by
    rw [heq]; exact (IsLocalizedModule.linearEquiv _ _ _).bijective
  exact this

/-! ## ★出典の紐付け(`.src`) -/

def tensorSectionMap_bijective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——基本開集合で比較射は全単射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
