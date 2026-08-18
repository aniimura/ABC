import ABC3.Found.Arakelov.PicAwaySection

/-!
# Arakelov (B1) 第 86 ブロック —— **テンソル局所化の一意性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★比較射の `D(f)` 成分を同定する器具

第 83 の前層射の `D(f)` 成分

    (tilde M)(D f) ⊗ (tilde N)(D f)  ⟶  (tilde (M ⊗_R N))(D f)

が同型であることを示したい。★**直接計算せず、局所化の一意性で押す**:

| 段 | 内容 |
|---|---|
| 1 | 両辺とも `M ⊗_R N` の `powers f` での局所化である |
| 2 | 局所化からの写像は**像で決まる**(`linearMap_ext`) |
| 3 | 2 つの局所化の間には**標準の線型同値**がある(`linearEquiv`) |

★★★2 と 3 から、「純テンソルで一致する」ことだけ見れば同型が出る。

## ★★機構 —— mathlib の在庫

| 在庫 | 内容 |
|---|---|
| `instance IsLocalizedModule S (TensorProduct.map fM fN)` | ★**テンソルは局所化を保つ** |
| `IsLocalizedModule.linearMap_ext` | ★像で決まる |
| `IsLocalizedModule.linearEquiv` | ★2 つの局所化の間の同値 |

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `isLocalizedModule_tensorMap` | ★テンソルは局所化を保つ(名付け) |
| `tensorLocalization_ext` | ★★**像で決まる** |
| `tensorLocalizationEquiv` | ★★★★**2 つの局所化の間の線型同値** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace TensorProduct

variable {R : Type u} [CommRing R] (S : Submonoid R)
  {M N M' N' : Type u} [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
  [AddCommGroup M'] [Module R M'] [AddCommGroup N'] [Module R N']
  (fM : M →ₗ[R] M') (fN : N →ₗ[R] N')
  [IsLocalizedModule S fM] [IsLocalizedModule S fN]

/-- ★**テンソルは局所化を保つ**——mathlib の在庫に名を付ける。 -/
theorem isLocalizedModule_tensorMap :
    IsLocalizedModule S (TensorProduct.map fM fN) := inferInstance

/-- ★★**局所化からの写像は像で決まる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで「純テンソルで一致する」ことだけ見れば写像が決まる。 -/
theorem tensorLocalization_ext {P : Type u} [AddCommGroup P] [Module R P]
    (g : M ⊗[R] N →ₗ[R] P) [IsLocalizedModule S g]
    (h h' : (M' ⊗[R] N') →ₗ[R] P)
    (e : h ∘ₗ TensorProduct.map fM fN = h' ∘ₗ TensorProduct.map fM fN) : h = h' :=
  IsLocalizedModule.linearMap_ext S (TensorProduct.map fM fN) g e

/-- ★★★★**2 つの局所化の間の線型同値**。

★★★これが比較射の `D(f)` 成分が同型であることの中身になる。 -/
noncomputable def tensorLocalizationEquiv {P : Type u} [AddCommGroup P] [Module R P]
    (g : M ⊗[R] N →ₗ[R] P) [IsLocalizedModule S g] :
    (M' ⊗[R] N') ≃ₗ[R] P :=
  IsLocalizedModule.linearEquiv S (TensorProduct.map fM fN) g

/-! ## ★出典の紐付け(`.src`) -/

def tensorLocalizationEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——テンソル局所化の一意性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
