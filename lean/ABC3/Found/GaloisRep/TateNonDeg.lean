import ABC3.Found.GaloisRep.TateInvert

/-!
# Galois (G6) 第 101 ブロック —— **★★★`E_q` は非退化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★`Δ(E_q) = q·(単元)` の帰結

第 98 ブロックで `Δ(E_{q}) = q · (単元)` を得た。★整域では

    q ≠ 0  ⟹  Δ(E_q) ≠ 0

すなわち **`E_q` は退化していない**。★★離散付値環なら

    v(Δ(E_q)) = v(q) > 0

——★★★これが原文 `Definition 3.3` の「局所高さ `v_K(q_E) ∈ ℤ_{>0}`」である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateCurveAt_Delta_ne_zero` | ★★`q ≠ 0` なら `Δ(E_q) ≠ 0` |
| `addVal_tateDelta` | ★★★**`v(Δ(E_q)) = v(q)`** |
| `addVal_tateDelta_pos` | ★★★**局所高さは正** |
-/

namespace ABC3.Found.GaloisRep

open PowerSeries

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★`q ≠ 0` なら `Δ(E_q) ≠ 0`——`E_q` は非退化である。 -/
theorem tateCurveAt_Delta_ne_zero [IsAdicComplete I R] [Nontrivial R] [NoZeroDivisors R]
    {q : R} (hq : q ∈ I) (hq0 : q ≠ 0) : (tateCurveAt q hq).Δ ≠ 0 := by
  obtain ⟨v, hv, hΔ⟩ := tateCurveAt_Delta_eq_mul_unit q hq
  rw [hΔ]
  exact mul_ne_zero hq0 hv.ne_zero

/-! ## ★★★離散付値環での局所高さ -/

variable {A : Type} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A]

/-- ★★★**`v(Δ(E_q)) = v(q)`**——局所高さは判別式の付値である。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem addVal_tateDelta [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {q : A} (hq : q ∈ IsLocalRing.maximalIdeal A) (hq0 : q ≠ 0) :
    IsDiscreteValuationRing.addVal A (tateCurveAt q hq).Δ
      = IsDiscreteValuationRing.addVal A q := by
  obtain ⟨v, hv, hΔ⟩ := tateCurveAt_Delta_eq_mul_unit q hq
  rw [hΔ, IsDiscreteValuationRing.addVal_mul,
    IsDiscreteValuationRing.addVal_eq_zero_iff.2 hv, add_zero]

/-- ★★★**局所高さは正である**(原文が `∈ ℤ_{>0}` と書いているもの)。 -/
theorem addVal_tateDelta_pos [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {q : A} (hq : q ∈ IsLocalRing.maximalIdeal A) (hq0 : q ≠ 0) :
    0 < IsDiscreteValuationRing.addVal A (tateCurveAt q hq).Δ := by
  rw [addVal_tateDelta hq hq0]
  have hnu : ¬ IsUnit q := (IsLocalRing.mem_maximalIdeal q).1 hq
  have hne : IsDiscreteValuationRing.addVal A q ≠ 0 := fun h =>
    hnu (IsDiscreteValuationRing.addVal_eq_zero_iff.1 h)
  exact pos_iff_ne_zero.2 hne

/-! ## ★出典の紐付け(`.src`) -/

def addVal_tateDelta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(局所高さ v_K(q_E) が正であること)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
