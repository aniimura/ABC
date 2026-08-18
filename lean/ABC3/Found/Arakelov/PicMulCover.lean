import ABC3.Found.Arakelov.PicIdealSq
import ABC3.Found.Arakelov.PicAffCover

/-!
# Arakelov (B2) 第 161 ブロック —— **`D(h·g)` 形の被覆**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★第 117 の一般スキーム版

局所自明性の組み立てでは、被覆の各員が `D(h·g)` の形である必要がある
——第 158(全単射性)がその形で述べてあるからである。

★`D(h) ⊆ U ≤ D(g)` なら `D(h·g) = D(h)` なので、
基本開集合の基底(`IsAffineOpen.exists_basicOpen_le`)がそのまま使える。

| 定理 | 内容 |
|---|---|
| `affineBasicOpenMulSieve_mem` | ★`D(h·g)` 形で被覆篩を生成 |
| `overAffineBasicMulPresieve_mem` | ★★`Over V` の site でも覆う |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★`D(h·g)` の形の基本開集合で覆える。 -/
theorem affineBasicOpenMulSieve_mem (A : X.affineOpens) (g : (Γ(X, A.1) : Type u))
    (U : X.Opens) (hU : U ≤ X.basicOpen g) :
    Sieve.generate (fun (Z : X.Opens) (_ : Z ⟶ U) =>
        ∃ h : (Γ(X, A.1) : Type u), Z = X.basicOpen (h * g))
      ∈ (Opens.grothendieckTopology X) U := by
  intro x hx
  obtain ⟨h, hle, hxh⟩ := A.2.exists_basicOpen_le (V := U) ⟨x, hx⟩
    (le_trans hU (X.basicOpen_le g) hx)
  have hgh : X.basicOpen (h * g) = X.basicOpen h := by
    rw [Scheme.basicOpen_mul]
    exact inf_eq_left.2 (le_trans hle hU)
  refine ⟨X.basicOpen h, homOfLE hle, ?_, hxh⟩
  exact ⟨X.basicOpen h, 𝟙 _, homOfLE hle, ⟨h, hgh.symm⟩, rfl⟩

/-- ★★`Over V` の site でも覆う。 -/
theorem overAffineBasicMulPresieve_mem (A : X.affineOpens) (g : (Γ(X, A.1) : Type u))
    (V : X.Opens) (hV : V ≤ X.basicOpen g) (W : Over V) :
    Sieve.generate (Presieve.functorPullback (Over.forget V)
        (fun (Z : X.Opens) (_ : Z ⟶ W.left) =>
          ∃ h : (Γ(X, A.1) : Type u), Z = X.basicOpen (h * g)))
      ∈ ((Opens.grothendieckTopology X).over V) W := by
  rw [← Sieve.overEquiv_symm_generate]
  refine symm_mem_over V W _ ?_
  exact affineBasicOpenMulSieve_mem A g W.left (le_trans (leOfHom W.hom) hV)


/-! ## ★出典の紐付け(`.src`) -/

def overAffineBasicMulPresieve_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——D(h·g) 形の被覆)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
