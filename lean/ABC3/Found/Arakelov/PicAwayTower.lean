import ABC3.Found.Arakelov.PicAwayTrans

/-!
# Arakelov (B1) 第 95 ブロック —— **`R → R_g → R_{t·g}` は足場をなす**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★自由性を小さい開集合へ運ぶ土台

第 94 ブロックの環準同型 `R_g → R_{t·g}` を**代数構造**にし、
`R → R_g → R_{t·g}` が `IsScalarTower` であることを示す。

★★これが `Module.free_of_isLocalizedModule`(自由性は局所化で保たれる)を
当てるための前提である。

## ★★機構

★`IsScalarTower.of_algebraMap_eq` + `IsLocalization.lift_eq`。
★★`awayToAwayMul` は `IsLocalization.Away.lift` なので、
`algebraMap` との合成が元の `algebraMap` に戻ることは**普遍性の等式そのもの**である。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `awayMulAlgebra` | ★`R_{t·g}` は `R_g` 代数 |
| `awayMulTower` | ★★**`R → R_g → R_{t·g}` は足場** |
-/

universe u

namespace ABC3.Found.Arakelov

variable (R : Type u) [CommRing R] (g t : R)

/-- ★**`R_{t·g}` は `R_g` 代数である**。 -/
noncomputable instance awayMulAlgebra : Algebra (Localization (Submonoid.powers g))
    (Localization (Submonoid.powers (t * g))) := (awayToAwayMul R g t).toAlgebra

/-- ★★**`R → R_g → R_{t·g}` は足場をなす**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで自由性を `D(g)` から `D(t·g)` へ運べる。 -/
instance awayMulTower : IsScalarTower R (Localization (Submonoid.powers g))
    (Localization (Submonoid.powers (t * g))) := by
  refine IsScalarTower.of_algebraMap_eq fun r => ?_
  show algebraMap R (Localization (Submonoid.powers (t * g))) r
    = awayToAwayMul R g t (algebraMap R (Localization (Submonoid.powers g)) r)
  rw [awayToAwayMul, IsLocalization.Away.lift]
  exact (IsLocalization.lift_eq _ r).symm

/-! ## ★出典の紐付け(`.src`) -/

def awayMulTower.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所化の足場)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
