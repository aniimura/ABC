import ABC3.Found.GaloisRep.Support
import Mathlib.RingTheory.DedekindDomain.AdicValuation

/-!
# Galois (G5) 第 152 ブロック —— **★★★★★★還元写像**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★剰余体を経由しない還元写像

D2 の最後の段(点の還元)の土台を作る。★**剰余体 `κ(v)` を構成せずに済む**:

    w t ≤ 1   ⟹   ∃! d ∈ F,  w(t − d) < 1

★★一意性は易しい——`d ≠ d'` なら `d − d'` は `F` の単元だから付値は 1 である。
★★★存在は mathlib の `valuationSubringAtPrime` が効く:

| mathlib(2026-08-20 実測) | 内容 |
|---|---|
| `HeightOneSpectrum.valuationSubringAtPrime K v` | 付値環 |
| `valuationSubringAtPrime_eq_valuationSubring` | それが付値環と一致 |
| `IsLocalization (v.asIdeal.primeCompl) (valuationSubringAtPrime K v)` | **局所化と一致**(instance) |

★したがって `w t ≤ 1` なら `t = a/b`(`a ∈ F[W]`、`b ∉ v.asIdeal`)と書け、
`d := ev(a)/ev(b)` と置けばよい(`ev` は第 117 の評価写像)。

## ★★★★評価写像の核

    ker(ev_{(c,y₀)}) = XYIdeal(c, y₀)

★包含は生成元を送るだけ。逆は `XYIdeal` が極大(第 138)で `ker ≠ ⊤` だから。

## ★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `evalPt_ker` | ★★★★**評価写像の核は `XYIdeal`** |
| `evalPt_ne_zero` | ★`b ∉ XYIdeal` なら `ev(b) ≠ 0` |
| `exists_const_of_valuation_le_one` | ★★★★★★**還元の存在** |
| `const_unique` | ★★★★還元先の定数は一意 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)

/-! ## ★★★★評価写像の核 -/

/-- ★★★★**評価写像の核は `XYIdeal` である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★包含は生成元を送るだけ。逆は `XYIdeal` が極大(第 138)で `ker ≠ ⊤` だから。 -/
theorem evalPt_ker {c y₀ : F} (h : W.Equation c y₀) :
    RingHom.ker (evalPt W h) = CoordinateRing.XYIdeal W c (Polynomial.C y₀) := by
  have hle : CoordinateRing.XYIdeal W c (Polynomial.C y₀) ≤ RingHom.ker (evalPt W h) := by
    rw [CoordinateRing.XYIdeal, Ideal.span_le]
    rintro z (rfl | rfl)
    · rw [SetLike.mem_coe, RingHom.mem_ker, XClass_eq, map_sub, evalPt_genX,
        evalPt_algebraMap, sub_self]
    · rw [SetLike.mem_coe, RingHom.mem_ker, YClass_eq, map_sub, evalPt_genY,
        evalPt_algebraMap, sub_self]
  have hne : RingHom.ker (evalPt W h) ≠ ⊤ := by
    intro htop
    have h1 : (1 : W.CoordinateRing) ∈ RingHom.ker (evalPt W h) := by rw [htop]; trivial
    rw [RingHom.mem_ker, map_one] at h1
    exact one_ne_zero h1
  exact ((xyIdeal_isMaximal W h).eq_of_le hne hle).symm

/-- ★`b ∉ XYIdeal` なら `ev(b) ≠ 0`。 -/
theorem evalPt_ne_zero {c y₀ : F} (h : W.Equation c y₀) {b : W.CoordinateRing}
    (hb : b ∉ CoordinateRing.XYIdeal W c (Polynomial.C y₀)) : evalPt W h b ≠ 0 := by
  intro h0
  refine hb ?_
  rw [← evalPt_ker W h]
  exact RingHom.mem_ker.2 h0

/-! ## ★★★★★★還元写像 -/

variable [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★★★★★**還元の存在**——`w t ≤ 1` なら `w(t − d) < 1` なる定数 `d ∈ F` がある。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `valuationSubringAtPrime` が**局所化と一致する**ことを使う。
したがって `t = a/b`(`b ∉ v.asIdeal`)と書け、`d := ev(a)/ev(b)` が求めるもの。
★★**剰余体 `κ(v)` を構成する必要は無い。** -/
theorem exists_const_of_valuation_le_one (v : HeightOneSpectrum W.CoordinateRing)
    {c y₀ : F} (h : W.Equation c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))
    {t : W.FunctionField} (ht : v.valuation W.FunctionField t ≤ 1) :
    ∃ d : F, v.valuation W.FunctionField (t - algebraMap F W.FunctionField d) < 1 := by
  set w := v.valuation W.FunctionField with hw
  set S := HeightOneSpectrum.valuationSubringAtPrime W.FunctionField v with hS
  have htS : t ∈ S := by
    rw [hS, HeightOneSpectrum.valuationSubringAtPrime_eq_valuationSubring,
      Valuation.mem_valuationSubring_iff]
    exact ht
  obtain ⟨⟨a, bb⟩, hab⟩ := IsLocalization.mk'_surjective (v.asIdeal.primeCompl)
    (S := S) (⟨t, htS⟩ : S)
  have hab' : IsLocalization.mk' (S : Type _) a bb = (⟨t, htS⟩ : S) := hab
  have hspec := IsLocalization.mk'_spec (S := (S : Type _)) a bb
  rw [hab'] at hspec
  obtain ⟨b, hbmem⟩ := bb
  have hcast : t * algebraMap W.CoordinateRing W.FunctionField b
      = algebraMap W.CoordinateRing W.FunctionField a := by
    have hq := congrArg (fun z : S => (z : W.FunctionField)) hspec
    simpa [IsScalarTower.algebraMap_apply W.CoordinateRing S W.FunctionField] using hq
  have hbnot : b ∉ CoordinateRing.XYIdeal W c (Polynomial.C y₀) := by
    rw [← hv]; exact hbmem
  have hevb : evalPt W h b ≠ 0 := evalPt_ne_zero W h hbnot
  refine ⟨evalPt W h a / evalPt W h b, ?_⟩
  set d := evalPt W h a / evalPt W h b with hd
  have hmem : a - algebraMap F W.CoordinateRing d * b
      ∈ CoordinateRing.XYIdeal W c (Polynomial.C y₀) := by
    rw [← evalPt_ker W h, RingHom.mem_ker, map_sub, map_mul, evalPt_algebraMap, hd]
    field_simp
    ring
  have hker : a - algebraMap F W.CoordinateRing d * b ∈ v.asIdeal := hv.symm ▸ hmem
  have hbval : w (algebraMap W.CoordinateRing W.FunctionField b) = 1 := by
    refine le_antisymm (HeightOneSpectrum.valuation_le_one v _) ?_
    by_contra hlt
    rw [not_le] at hlt
    exact hbnot (hv ▸ (HeightOneSpectrum.valuation_lt_one_iff_mem v _).1 hlt)
  have hprod : (t - algebraMap F W.FunctionField d)
      * algebraMap W.CoordinateRing W.FunctionField b
      = algebraMap W.CoordinateRing W.FunctionField
        (a - algebraMap F W.CoordinateRing d * b) := by
    rw [map_sub, map_mul, sub_mul, hcast,
      ← IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
  have hfin := congrArg w hprod
  rw [Valuation.map_mul, hbval, mul_one] at hfin
  rw [hfin]
  exact (HeightOneSpectrum.valuation_lt_one_iff_mem v _).2 hker

/-- ★★★★**還元先の定数は一意である**——`d ≠ d'` なら `d − d'` は単元で付値が 1。 -/
theorem const_unique (v : HeightOneSpectrum W.CoordinateRing) {t : W.FunctionField} {d d' : F}
    (hd : v.valuation W.FunctionField (t - algebraMap F W.FunctionField d) < 1)
    (hd' : v.valuation W.FunctionField (t - algebraMap F W.FunctionField d') < 1) : d = d' := by
  set w := v.valuation W.FunctionField with hw
  by_contra hne
  have hsub : algebraMap F W.FunctionField (d' - d)
      = (t - algebraMap F W.FunctionField d) - (t - algebraMap F W.FunctionField d') := by
    rw [map_sub]; ring
  have hlt : w (algebraMap F W.FunctionField (d' - d)) < 1 := by
    rw [hsub, sub_eq_add_neg]
    refine lt_of_le_of_lt (Valuation.map_add w _ _) (max_lt hd ?_)
    rwa [Valuation.map_neg]
  have hunit : w (algebraMap F W.FunctionField (d' - d)) = 1 := by
    rw [IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
    exact valuation_algebraMap_isUnit v ((isUnit_iff_ne_zero.2
      (sub_ne_zero.2 (Ne.symm hne))).map (algebraMap F W.CoordinateRing))
  rw [hunit] at hlt
  exact absurd hlt (lt_irrefl 1)

/-! ## ★出典の紐付け(`.src`) -/

def exists_const_of_valuation_le_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——付値環の元が定数へ還元されること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
