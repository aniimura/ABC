/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal

/-!
# 第 953 ブロック —— **★★★★★★★★★★★★`v_p(l) = 0` なら `l` は完備化の整数環で単元**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——(D3) の (c)

`tateParam_quot_velu_of_torsion`（第 948）は `hlu : IsUnit ((l : R))` を受ける。
原文の仮定は「l is prime to … the residue characteristics」であるから、
これは **`v_p(l) = 0`** から出なければならない。

★道は 2 段である:

1. `vAdd = 0` なら `x` も `x⁻¹` も環の元（`dvr_mem_of_nonneg`、第 895）——だから単元
2. `vAdd (l) = valAdd p (l)`（`vAdd_algebraMap_eq_valAdd`、第 893）

| 定理 | 内容 |
|---|---|
| `isUnit_of_vAdd_eq_zero` | ★★★★★★付値が `0` の元は環の単元 |
| `isUnit_natCast_of_valAdd_eq_zero` | ★★★★★★★★★★★★**`v_p(n) = 0` なら `n` は単元** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField
open scoped Classical

section UnitVal

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★★★★★**付値が `0` の元は環の単元**。

☆`x` も `x⁻¹` も環に入るので、その積は `1` である。 -/
theorem isUnit_of_vAdd_eq_zero (x : Kˣ) (h : vAdd (tateDvrVal R K) x = 0)
    (y : R) (hy : algebraMap R K y = (x : K)) : IsUnit y := by
  obtain ⟨z, hz⟩ := dvr_mem_of_nonneg (R := R) x⁻¹ (by rw [vAdd_inv, h]; simp)
  have hyz : y * z = 1 := by
    refine IsFractionRing.injective R K ?_
    rw [map_mul, hy, hz, map_one, Units.val_inv_eq_inv_val]
    field_simp
  exact ⟨⟨y, z, hyz, by rw [mul_comm]; exact hyz⟩, rfl⟩

end UnitVal

section NatCast

variable {L : Type} [Field L] [NumberField L]
  {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★★★★★★★★★★★★**`v_p(n) = 0` なら `n` は完備化の整数環で単元**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 953）**——これが `tateParam_quot_velu_of_torsion`（第 948）の
`hlu : IsUnit ((l : R))` の出どころである。
☆原文の「`l` は残余標数と互いに素」が `v_p(l) = 0` に当たる。 -/
theorem isUnit_natCast_of_valAdd_eq_zero (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x)
    (n : ℕ) (hn : (n : L) ≠ 0) (hnv : (n : Lv) ≠ 0)
    (hval : valAdd p (Units.mk0 ((n : L)) hn) = 0) :
    IsUnit ((n : R)) := by
  have hmapL : algebraMap L Lv ((n : L)) = (n : Lv) := by
    rw [map_natCast]
  have hne' : algebraMap L Lv ((n : L)) ≠ 0 := by rw [hmapL]; exact hnv
  have hbridge := vAdd_algebraMap_eq_valAdd (R := R) p hp ((n : L)) hn hne'
  rw [hval] at hbridge
  refine isUnit_of_vAdd_eq_zero (R := R) (Units.mk0 (algebraMap L Lv ((n : L))) hne')
    hbridge ((n : R)) ?_
  show algebraMap R Lv ((n : R)) = algebraMap L Lv ((n : L))
  rw [map_natCast, hmapL]

end NatCast

/-! ## ★出典の紐付け(`.src`) -/

def isUnit_of_vAdd_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(付値が 0 の元は環の単元。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isUnit_natCast_of_valAdd_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(n) = 0 なら n は完備化の整数環で単元。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
