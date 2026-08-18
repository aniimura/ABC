import ABC3.Found.Arakelov.PicOverCover

/-!
# Arakelov (B1) 第 117 ブロック —— **`D(h·g)` の生成族は覆う**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性の被覆

`U ≤ D(g)` のとき、**`D(h·g)` の形の基本開集合**が `U` を覆う。

★★機構: `U` の点 `x` に対し基本開集合の基底から `D(h) ⊆ U` を取る。
★★★`D(h) ⊆ U ≤ D(g)` なので `D(h·g) = D(h) ⊓ D(g) = D(h)` である。

★これで第 99 ブロック(自由性は `D(g)` から `D(t·g)` へ運べる)が
**そのまま当たる形**の被覆が得られた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `basicOpenMulSieve_mem` | ★★★★**`D(h·g)` の生成族は `U` を覆う** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (g : (R : Type u))

/-- ★★★★**`D(h·g)` の形の基本開集合が `U` を覆う**(`U ≤ D(g)` のとき)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`D(h) ⊆ U ≤ D(g)` なら `D(h·g) = D(h)` だから、
基本開集合の基底がそのまま `D(h·g)` の族になる。 -/
theorem basicOpenMulSieve_mem (U : (Spec R).Opens) (hU : U ≤ PrimeSpectrum.basicOpen g) :
    Sieve.generate (fun (Z : (Spec R).Opens) (_ : Z ⟶ U) =>
        ∃ h : (R : Type u), Z = PrimeSpectrum.basicOpen (h * g))
      ∈ (Opens.grothendieckTopology (Spec R)) U := by
  intro x hx
  obtain ⟨W, ⟨h, rfl⟩, hxW, hWU⟩ :=
    Opens.isBasis_iff_nbhd.1 PrimeSpectrum.isBasis_basic_opens hx
  have hgW : PrimeSpectrum.basicOpen (h * g) = PrimeSpectrum.basicOpen h := by
    rw [PrimeSpectrum.basicOpen_mul]
    exact inf_eq_left.2 (le_trans hWU hU)
  refine ⟨PrimeSpectrum.basicOpen h, homOfLE hWU, ?_, hxW⟩
  exact ⟨PrimeSpectrum.basicOpen h, 𝟙 _, homOfLE hWU, ⟨h, hgW.symm⟩, rfl⟩

/-! ## ★出典の紐付け(`.src`) -/

def basicOpenMulSieve_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——D(h·g) の生成族は覆う)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
