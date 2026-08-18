import ABC3.Found.Arakelov.PicIdealFree
import ABC3.Found.Arakelov.PicOverCover

/-!
# Arakelov (B2) 第 155 ブロック —— **アフィン開の基本開集合による被覆**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★第 117・122 を**一般のスキーム**へ

第 117(`basicOpenMulSieve_mem`)と第 122(`overBasicPresieve_mem`)は
`Spec R` 専用だった。★本ブロックはそれを

    A : X.affineOpens、U ≤ A、基本開集合 X.basicOpen h(h : Γ(X,A))

の形で述べ直す。★★第 116(`symm_mem_over`)は**もともと一般**だったので、
`Over V` への持ち上げはそのまま使える。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `affineBasicOpenSieve_mem` | ★基本開集合は被覆篩を生成する |
| `overAffineBasicPresieve_mem` | ★★`Over V` の site でも覆う |

★★★要は `IsAffineOpen.exists_basicOpen_le`(mathlib)一発である。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★アフィン開の基本開集合は被覆篩を生成する。 -/
theorem affineBasicOpenSieve_mem (A : X.affineOpens) (U : X.Opens) (hU : U ≤ A.1) :
    Sieve.generate (fun (Z : X.Opens) (_ : Z ⟶ U) =>
        ∃ h : (Γ(X, A.1) : Type u), Z = X.basicOpen h)
      ∈ (Opens.grothendieckTopology X) U := by
  intro x hx
  obtain ⟨h, hle, hxh⟩ := A.2.exists_basicOpen_le (V := U) ⟨x, hx⟩ (hU hx)
  refine ⟨X.basicOpen h, homOfLE hle, ?_, hxh⟩
  exact ⟨X.basicOpen h, 𝟙 _, homOfLE hle, ⟨h, rfl⟩, rfl⟩

/-- ★★`Over V` の site でも覆う。 -/
theorem overAffineBasicPresieve_mem (A : X.affineOpens) (V : X.Opens) (hV : V ≤ A.1)
    (W : Over V) :
    Sieve.generate (Presieve.functorPullback (Over.forget V)
        (fun (Z : X.Opens) (_ : Z ⟶ W.left) =>
          ∃ h : (Γ(X, A.1) : Type u), Z = X.basicOpen h))
      ∈ ((Opens.grothendieckTopology X).over V) W := by
  rw [← Sieve.overEquiv_symm_generate]
  refine symm_mem_over V W _ ?_
  exact affineBasicOpenSieve_mem A W.left (le_trans (leOfHom W.hom) hV)


/-! ## ★出典の紐付け(`.src`) -/

def overAffineBasicPresieve_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——アフィン開の基本開集合による被覆)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
