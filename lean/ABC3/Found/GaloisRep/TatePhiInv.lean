import ABC3.Found.GaloisRep.TatePhiSurj

/-!
# Galois (G6) 第 291 ブロック —— **★★★★★★★★★`Φ(c⁻¹) = −Φ(c)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点

> **`Φ(c⁻¹) = −Φ(c)`**(`tatePhi_inv`)

★★★群構造の**半分**である。残りは `Φ(c₁c₂) = Φ(c₁) + Φ(c₂)` の側。

## ★★★★★★★逆元の代表元は 2 通り

類 `c` の正規化代表元を `u`(対 `(a,w)`)とすると、`c⁻¹` の代表元は

| `v(u)` | `c⁻¹` の対 | 使う反転則 |
|---|---|---|
| `> 0`(環帯) | `(w, a)`——**相方に入れ替えるだけ** | `tateXpair_symm`(第 210) |
| `= 0`(単元・原点近傍) | `(1/a, q·a)` | `tateXK_ringInverse`(第 281) |

★★★★どちらの対も**基本領域に入っている**ことを確かめるのが仕事の半分である:
`v(w) = v(Q) − v(u) ∈ (0, v(Q))`、`v(1/u) = −0 = 0`。

★★`tateAOf_eq_of_pair`——「正規化した対とその類が分かれば `tateAOf`・`tateWOf` が
決まる」——を作っておくと、両方の場合が 4 行で済む。

## ★★★点の水準の反転則を類の水準へ

第 281・第 282 で作った 2 つの反転則(`tatePtPair_swap`・`tatePtPair_ringInverse`)が、
ここでちょうど 2 つの場合に対応する。★★★**`Kˣ/q^ℤ` の元の代表元が 2 通りあることが、
反転則が 2 つある理由だった**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `mk_Q_eq_one`・`vAdd_eq_zero_of_isUnit` | ★★★下ごしらえ |
| `tateAOf_eq_of_pair` | ★★★★★★正規化した対からの同定 |
| `tateAOf_inv_annulus`・`tateAOf_inv_unit` | ★★★★★★★逆元の代表元 |
| `tatePhi_inv` | ★★★★★★★★★**`Φ(c⁻¹) = −Φ(c)`** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-! ## ★★★下ごしらえ -/

theorem mk_Q_eq_one (Q : Kˣ) : (QuotientGroup.mk Q : Kˣ ⧸ Subgroup.zpowers Q) = 1 :=
  (QuotientGroup.eq_one_iff _).2 (Subgroup.mem_zpowers Q)

/-- ★★★★`R` の単元の付値は 0。 -/
theorem vAdd_eq_zero_of_isUnit (S : TateSetup R I K)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (a : R) (ha : IsUnit a) (hne : algebraMap R K a ≠ 0) :
    vAdd S.v (Units.mk0 (algebraMap R K a) hne) = 0 := by
  set u : Kˣ := Units.mk0 (algebraMap R K a) hne with hu
  have h1 : 0 ≤ vAdd S.v u := hvR u ⟨a, rfl⟩
  have hinv : algebraMap R K (Ring.inverse a) = ((u⁻¹ : Kˣ) : K) := by
    have h2 := congrArg (algebraMap R K) (Ring.inverse_mul_cancel a ha)
    rw [map_mul, map_one] at h2
    rw [Units.val_inv_eq_inv_val, hu]
    show algebraMap R K (Ring.inverse a) = (algebraMap R K a)⁻¹
    exact eq_inv_of_mul_eq_one_left h2
  have h2 : 0 ≤ vAdd S.v u⁻¹ := hvR u⁻¹ ⟨Ring.inverse a, hinv⟩
  rw [vAdd_inv] at h2
  omega

/-! ## ★★★★★★正規化した対からの同定 -/

/-- ★★★★★★正規化した対から `tateAOf`・`tateWOf` を同定する。 -/
theorem tateAOf_eq_of_pair (S : TateSetup R I K) (d : Kˣ ⧸ Subgroup.zpowers S.Q)
    (a w : R) (u : Kˣ) (hu : algebraMap R K a = (u : K)) (haw : a * w = S.q)
    (h0 : 0 ≤ vAdd S.v u) (h1 : vAdd S.v u < vAdd S.v S.Q)
    (hcls : QuotientGroup.mk u = d) : tateAOf S d = a ∧ tateWOf S d = w :=
  pair_eq_of_same_class' S.hinj S.v S.Q S.hQ S.q S.hQq (tateAOf_spec S d).1 hu
    (tateAOf_mul_tateWOf S d) haw (normRep_nonneg _ _ _ _) (normRep_lt _ _ _ _) h0 h1
    (by rw [normRep_mk, hcls])

/-! ## ★★★★★★★逆元の代表元 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**環帯の類の逆元は相方**。 -/
theorem tateAOf_inv_annulus (S : TateSetup R I K)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0) (c : Kˣ ⧸ Subgroup.zpowers S.Q) (ha : tateAOf S c ∈ I) :
    tateAOf S c⁻¹ = tateWOf S c ∧ tateWOf S c⁻¹ = tateAOf S c := by
  have haw : tateAOf S c * tateWOf S c = S.q := tateAOf_mul_tateWOf S c
  have ha0 : tateAOf S c ≠ 0 := by
    intro h
    rw [h, zero_mul] at haw
    exact hq0 haw.symm
  have hw0 : tateWOf S c ≠ 0 := by
    intro h
    rw [h, mul_zero] at haw
    exact hq0 haw.symm
  have hAK : algebraMap R K (tateAOf S c) ≠ 0 := fun h => ha0 (S.hinj (by rw [h, map_zero]))
  have hWK : algebraMap R K (tateWOf S c) ≠ 0 := fun h => hw0 (S.hinj (by rw [h, map_zero]))
  set ua : Kˣ := Units.mk0 (algebraMap R K (tateAOf S c)) hAK with hua
  set uw : Kˣ := Units.mk0 (algebraMap R K (tateWOf S c)) hWK with huw
  have huaw : ua * uw = S.Q := by
    refine Units.ext ?_
    show algebraMap R K (tateAOf S c) * algebraMap R K (tateWOf S c) = (S.Q : K)
    rw [← map_mul, haw, S.hQq]
  have hmka : (QuotientGroup.mk ua : Kˣ ⧸ Subgroup.zpowers S.Q) = c := by
    have h : ua = normRep S.v S.Q S.hQ c := Units.ext (tateAOf_spec S c).1
    rw [h, normRep_mk]
  have hmkw : (QuotientGroup.mk uw : Kˣ ⧸ Subgroup.zpowers S.Q) = c⁻¹ := by
    have h : uw = ua⁻¹ * S.Q := by
      rw [← huaw]; group
    rw [h, QuotientGroup.mk_mul, QuotientGroup.mk_inv, hmka, mk_Q_eq_one]
    exact mul_one (c⁻¹)
  refine tateAOf_eq_of_pair S c⁻¹ (tateWOf S c) (tateAOf S c) uw rfl (by rw [← haw]; ring)
    (hvR uw ⟨tateWOf S c, rfl⟩) ?_ hmkw
  have hpos : 0 < vAdd S.v ua := hvI ua ⟨tateAOf S c, ha, rfl⟩
  have hsum : vAdd S.v S.Q = vAdd S.v ua + vAdd S.v uw := by rw [← huaw, vAdd_mul]
  omega

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**単元の類の逆元は `1/u`**。 -/
theorem tateAOf_inv_unit (S : TateSetup R I K)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hq0 : S.q ≠ 0) (c : Kˣ ⧸ Subgroup.zpowers S.Q) (ha : IsUnit (tateAOf S c)) :
    tateAOf S c⁻¹ = Ring.inverse (tateAOf S c)
      ∧ tateWOf S c⁻¹ = S.q * tateAOf S c := by
  have haw : tateAOf S c * tateWOf S c = S.q := tateAOf_mul_tateWOf S c
  have ha0 : tateAOf S c ≠ 0 := by
    intro h
    rw [h, zero_mul] at haw
    exact hq0 haw.symm
  have hAK : algebraMap R K (tateAOf S c) ≠ 0 := fun h => ha0 (S.hinj (by rw [h, map_zero]))
  set ua : Kˣ := Units.mk0 (algebraMap R K (tateAOf S c)) hAK with hua
  have hmka : (QuotientGroup.mk ua : Kˣ ⧸ Subgroup.zpowers S.Q) = c := by
    have h : ua = normRep S.v S.Q S.hQ c := Units.ext (tateAOf_spec S c).1
    rw [h, normRep_mk]
  have hinvK : algebraMap R K (Ring.inverse (tateAOf S c)) = ((ua⁻¹ : Kˣ) : K) := by
    have h2 := congrArg (algebraMap R K) (Ring.inverse_mul_cancel (tateAOf S c) ha)
    rw [map_mul, map_one] at h2
    rw [Units.val_inv_eq_inv_val, hua]
    show algebraMap R K (Ring.inverse (tateAOf S c)) = (algebraMap R K (tateAOf S c))⁻¹
    exact eq_inv_of_mul_eq_one_left h2
  have hv0 : vAdd S.v ua = 0 := vAdd_eq_zero_of_isUnit S hvR (tateAOf S c) ha hAK
  refine tateAOf_eq_of_pair S c⁻¹ (Ring.inverse (tateAOf S c)) (S.q * tateAOf S c) ua⁻¹
    hinvK ?_ ?_ ?_ ?_
  · calc Ring.inverse (tateAOf S c) * (S.q * tateAOf S c)
        = S.q * (Ring.inverse (tateAOf S c) * tateAOf S c) := by ring
      _ = S.q := by rw [Ring.inverse_mul_cancel _ ha, mul_one]
  · rw [vAdd_inv, hv0]
    omega
  · rw [vAdd_inv, hv0]
    simpa using S.hQ
  · rw [QuotientGroup.mk_inv, hmka]

/-! ## ★★★★★★★★★`Φ(c⁻¹) = −Φ(c)` -/

section Inv

variable [DecidableEq K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★**`Φ(c⁻¹) = −Φ(c)`**——類の水準の反転則。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_inv (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c : Kˣ ⧸ Subgroup.zpowers S.Q) : tatePhi S hΔ c⁻¹ = -tatePhi S hΔ c := by
  by_cases hc : c = 1
  · subst hc
    rw [inv_one, tatePhi_one, neg_zero]
  · have hc' : c⁻¹ ≠ 1 := fun h => hc (by rw [← inv_inv c, h, inv_one])
    have e1 : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.a₁ = 1 := by
      show algebraMap R K ((tateCurveAt S.q S.hq).a₁) = 1
      rw [show (tateCurveAt S.q S.hq).a₁ = 1 by
        simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_one]
    have e3 : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.a₃ = 0 := by
      show algebraMap R K ((tateCurveAt S.q S.hq).a₃) = 0
      rw [show (tateCurveAt S.q S.hq).a₃ = 0 by
        simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
    rw [tatePhi_eq S hΔ hc, tatePhi_eq S hΔ hc', tatePtPair, tatePtPair, Point.neg_some]
    simp only [Point.some.injEq]
    have hawc : tateAOf S c * tateWOf S c = S.q := tateAOf_mul_tateWOf S c
    rcases hloc (tateAOf S c) with hu | hm
    · obtain ⟨h1, h2⟩ := tateAOf_inv_unit S hvR hq0 c hu
      have hw : wOf S.q (tateAOf S c) = tateWOf S c := wOf_eq_of_mul hu hawc
      have hne : algebraMap R K (1 - tateAOf S c) ≠ 0 := tateAOf_ne_one S hc
      constructor
      · rw [h1, h2, ← hw]
        exact tateXK_ringInverse (tateAOf S c) S.q S.hq hu hne
      · rw [WeierstrassCurve.Affine.negY, e1, e3, h1, h2, ← hw,
          tateYK_ringInverse (tateAOf S c) S.q S.hq hu hne]
        ring
    · obtain ⟨h1, h2⟩ := tateAOf_inv_annulus S hvR hvI hq0 c hm
      have hua : IsUnit (1 - tateAOf S c) := isUnit_one_sub hm
      have huw : IsUnit (1 - tateWOf S c) := isUnit_one_sub (tateWOf_mem S c)
      constructor
      · rw [h1, h2, tateXK_eq _ _ _ _ huw, tateXK_eq _ _ _ _ hua, tateXpair_symm]
      · rw [WeierstrassCurve.Affine.negY, e1, e3, h1, h2,
          tateYK_eq _ _ _ _ huw, tateXK_eq _ _ _ _ hua, tateYK_eq _ _ _ _ hua,
          tateYpair_swap]
        simp only [map_sub, map_neg]
        ring

end Inv

/-! ## ★出典の紐付け(`.src`) -/

def tateAOf_eq_of_pair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——正規化した対からの同定)",
    sectionId := "genell-def-3-3" }

def tatePhi_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Phi(c^{-1}) = -Phi(c))",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
