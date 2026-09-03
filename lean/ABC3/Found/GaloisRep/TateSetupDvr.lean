/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateDoubling
import ABC3.Found.GaloisRep.TateDvrSetup

/-!
# Galois (G6) 第 899 ブロック —— **★★★★★★★★★★★★★★★★★★★★完備 DVR なら
Tate 一意化は無条件に存在する**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★これは何か

`tatePhiAddEquivAll`（`TateDoubling.lean`）は `TateSetup` に加えて 6 つの仮説
（`hloc`・`hI`・`hquad`・`hvalring`・`hvR`・`hvI`）を受けている。

★本ブロックは**その 6 つをすべて DVR の事実から作る**。
したがって、完備な DVR `R` とその分数体 `K`、`0 ≠ q ∈ 𝔪` さえあれば

    **`Additive (Kˣ/q^ℤ) ≃+ E_q(K)`**

が**無条件に**得られる。

| 仮説 | どこから |
|---|---|
| `hloc` | `IsLocalRing`（単元でないなら `𝔪` の元） |
| `hI` | `𝔪` は極大だから素 |
| `hquad` | ★DVR は**整閉**（`IsIntegrallyClosed`） |
| `hvalring` | ★DVR は**付値環**（`dvr_mem_of_nonneg` / `dvr_mem_max_of_pos`） |
| `hvR`・`hvI` | `tateDvrVal_nonneg_of_mem` / `tateDvrVal_pos_of_mem_max` |

☆これで `Lemma 3.5` の残る 3 つのうち**「Tate 一意化 `Φ` の存在」が消える**。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K] [IsFractionRing R K]

/-! ## ★★★★★★DVR から `TateSetup` を作る -/

/-- ★★★★★★**完備な DVR と `0 ≠ q ∈ 𝔪` から `TateSetup` を作る**。 -/
noncomputable def mkTateSetup (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0) :
    TateSetup R (IsLocalRing.maximalIdeal R) K where
  v := tateDvrVal R K
  Q := Units.mk0 (algebraMap R K q)
    ((map_ne_zero_iff _ (IsFractionRing.injective R K)).2 hq0)
  hQ := tateDvrVal_pos_of_mem_max _ ⟨q, hq, rfl⟩
  q := q
  hq := hq
  hQq := rfl
  hinj := IsFractionRing.injective R K
  hmem := fun x hx => dvr_mem_max_of_pos x hx
  hmem0 := fun x hx => dvr_mem_of_nonneg x hx

@[simp] theorem mkTateSetup_q (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0) :
    (mkTateSetup (K := K) q hq hq0).q = q := rfl

@[simp] theorem mkTateSetup_v (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0) :
    (mkTateSetup (K := K) q hq hq0).v = tateDvrVal R K := rfl

/-! ## ★★★★★★ 6 つの仮説を DVR の事実から作る -/

/-- ★`R` は局所環なので、元は単元か `𝔪` の元である。 -/
theorem dvrTate_hloc : ∀ x : R, IsUnit x ∨ x ∈ IsLocalRing.maximalIdeal R := by
  intro x
  by_cases h : IsUnit x
  · exact Or.inl h
  · exact Or.inr ((IsLocalRing.mem_maximalIdeal x).2 h)

/-- ★`𝔪` は素イデアルである。 -/
theorem dvrTate_hI : (IsLocalRing.maximalIdeal R).IsPrime :=
  (IsLocalRing.maximalIdeal.isMaximal R).isPrime

/-- ★★**DVR は整閉**——モニックな 2 次式の根は環の元である。 -/
theorem dvrTate_hquad : ∀ (t : K) (b c : R),
    t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t := by
  intro t b c h
  have hint : IsIntegral R t := by
    refine ⟨Polynomial.X ^ 2 + Polynomial.C b * Polynomial.X + Polynomial.C c, ?_, ?_⟩
    · monicity!
    · simp [h]
  exact IsIntegrallyClosed.isIntegral_iff.1 hint

/-- ★★**DVR は付値環**——`t` か `t⁻¹` のどちらかが環に入る。 -/
theorem dvrTate_hvalring : ∀ t : K, t ≠ 0 →
    (∃ r : R, algebraMap R K r = t)
      ∨ (∃ r ∈ IsLocalRing.maximalIdeal R, algebraMap R K r = t⁻¹) := by
  intro t ht
  rcases le_or_gt 0 (vAdd (tateDvrVal R K) (Units.mk0 t ht)) with h | h
  · exact Or.inl (dvr_mem_of_nonneg (Units.mk0 t ht) h)
  · refine Or.inr ?_
    have hinv : 0 < vAdd (tateDvrVal R K) (Units.mk0 t ht)⁻¹ := by
      rw [vAdd_inv]
      omega
    obtain ⟨y, hy, hval⟩ := dvr_mem_max_of_pos _ hinv
    exact ⟨y, hy, by rw [hval]; rfl⟩

theorem dvrTate_hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) →
    0 ≤ vAdd (tateDvrVal R K) t := fun t h => tateDvrVal_nonneg_of_mem t h

theorem dvrTate_hvI : ∀ t : Kˣ,
    (∃ r ∈ IsLocalRing.maximalIdeal R, algebraMap R K r = (t : K)) →
    0 < vAdd (tateDvrVal R K) t := fun t h => tateDvrVal_pos_of_mem_max t h

/-! ## ★★★★★★★★★★★★★★★★★★★★結論——Tate 一意化は無条件に存在する -/

/-- ★★★★★★★★★★★★★★★★★★★★**完備な DVR なら Tate 一意化は無条件に存在する**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★引数は `q` と `Δ ≠ 0` だけである——**6 つの仮説はすべて内製した**。 -/
noncomputable def dvrTatePhiAddEquiv [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    Additive (Kˣ ⧸ Subgroup.zpowers (mkTateSetup (K := K) q hq hq0).Q)
      ≃+ ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
        (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Point :=
  tatePhiAddEquivAll (mkTateSetup q hq hq0) dvrTate_hloc dvrTate_hI dvrTate_hquad
    dvrTate_hvalring dvrTate_hvR dvrTate_hvI hq0 hΔ

/-! ## ★出典の紐付け(`.src`) -/

def mkTateSetup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化—完備 DVR から TateSetup を作る。★無条件)",
    sectionId := "genell-def-3-3" }

def dvrTatePhiAddEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化—完備 DVR なら無条件に存在する。★無条件)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
