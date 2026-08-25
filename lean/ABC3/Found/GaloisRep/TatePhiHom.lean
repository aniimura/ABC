import ABC3.Found.GaloisRep.TatePhiInv

/-!
# Galois (G6) 第 292 ブロック —— **★★★★★★★★★★一般の位置での準同型性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点

> `a(c₁)·a(c₂)·a(c₃) = q` で 3 つの類が一般の位置なら
> **`Φ(c₁) + Φ(c₂) + Φ(c₃) = 0`**(`tatePhi_add_add_eq_zero`)

★★★葉 (c)(共線性)と群法則(第 279)が、類の水準に上がった。

## ★★★★★★★★相異性は `±` の分岐だけ

3 点の `x` 座標が相異なることが群法則の要である。**`Φ` が単射**(第 285)と
**`Φ(c⁻¹) = −Φ(c)`**(第 291)があれば、曲線の式の差から

    X(c) = X(d)  ⟹  (Y(c) − Y(d))·(Y(c) + Y(d) + X(c)) = 0
                 ⟹  Y(c) = Y(d)(第 1 の枝)または `P(d) = −P(c)`(第 2 の枝)
                 ⟹  **`c = d` または `c·d = 1`**(`tatePhi_X_eq_imp`)

★★★★★領域ごとの単射性を**もう一度なぞる必要はない**——`Φ` の単射性という
1 つの事実に集約されている。第 285 で 5 通りを貼り合わせておいた甲斐があった。

## ★★★★相方は掛け算で出る

`a₁·w₁ = q = a₁·(a₂a₃)` と `a₁ ≠ 0` から **`w₁ = a₂·a₃`**。
★したがって「相方が `𝔪` に入る」(群法則の仮定)は `w₁ ∈ 𝔪`(定義)そのものである。
★★付値の議論は要らなかった。

## ★★残っているもの

`a(c₁)a(c₂)a(c₃) = q` は仮定として受けている。`c₁c₂c₃ = 1` からこれを出すには
正規化代表元の付値の和が `v(q)` であることを見る(和は `0`・`v(q)`・`2v(q)` の 3 通り)。
また退化した場合(`c₁ = c₂` など)は補助母数で別に扱う。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tatePhi_X_eq_imp` | ★★★★★★★★**`X` の一致は `c = d` か `c·d = 1`** |
| `tatePhi_add_add_eq_zero` | ★★★★★★★★★★**一般の位置での準同型性** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

/-! ## ★★★★★★★★`X` の一致は `c = d` か `c·d = 1` -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**`X` が一致するのは `c = d` か `c·d = 1` のときだけ**。

★`Φ` の単射性(第 285)と反転則(第 291)だけで出る。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_X_eq_imp (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {c d : Kˣ ⧸ Subgroup.zpowers S.Q} (hc : c ≠ 1) (hd : d ≠ 1)
    (hX : tateXK (tateAOf S c) (tateWOf S c) S.q S.hq
      = (tateXK (tateAOf S d) (tateWOf S d) S.q S.hq : K)) :
    c = d ∨ c * d = 1 := by
  have e1 : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.a₁ = 1 := by
    show algebraMap R K ((tateCurveAt S.q S.hq).a₁) = 1
    rw [show (tateCurveAt S.q S.hq).a₁ = 1 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_one]
  have e3 : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.a₃ = 0 := by
    show algebraMap R K ((tateCurveAt S.q S.hq).a₃) = 0
    rw [show (tateCurveAt S.q S.hq).a₃ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  have hec := tate_equationK (K := K) (tateAOf S c) (tateWOf S c) S.q S.hq
    (tateAOf_mul_tateWOf S c) (isUnit_one_sub (tateWOf_mem S c)) (tateAOf_ne_one S hc)
  have hed := tate_equationK (K := K) (tateAOf S d) (tateWOf S d) S.q S.hq
    (tateAOf_mul_tateWOf S d) (isUnit_one_sub (tateWOf_mem S d)) (tateAOf_ne_one S hd)
  rw [← hX] at hed
  have hfac : (tateYK (K := K) (tateAOf S c) (tateWOf S c) S.q S.hq
        - tateYK (K := K) (tateAOf S d) (tateWOf S d) S.q S.hq)
      * (tateYK (K := K) (tateAOf S c) (tateWOf S c) S.q S.hq
        + tateYK (K := K) (tateAOf S d) (tateWOf S d) S.q S.hq
        + tateXK (tateAOf S c) (tateWOf S c) S.q S.hq) = 0 := by
    linear_combination hec - hed
  rcases mul_eq_zero.1 hfac with h | h
  · refine Or.inl (tatePhi_injective S hloc hΔ ?_)
    rw [tatePhi_eq S hΔ hc, tatePhi_eq S hΔ hd, tatePtPair, tatePtPair]
    simp only [Point.some.injEq]
    exact ⟨hX, sub_eq_zero.1 h⟩
  · refine Or.inr ?_
    have hdc : tatePhi S hΔ d = tatePhi S hΔ c⁻¹ := by
      rw [tatePhi_inv S hloc hvR hvI hq0 hΔ c, tatePhi_eq S hΔ hd, tatePhi_eq S hΔ hc,
        tatePtPair, tatePtPair, Point.neg_some]
      simp only [Point.some.injEq]
      refine ⟨hX.symm, ?_⟩
      rw [WeierstrassCurve.Affine.negY, e1, e3]
      linear_combination h
    have hdd := tatePhi_injective S hloc hΔ hdc
    rw [hdd, mul_inv_cancel]

/-! ## ★★★★★★★★★★一般の位置での準同型性 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**一般の位置での準同型性**——3 つの類の和は 0。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_add_add_eq_zero (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c₁ c₂ c₃ : Kˣ ⧸ Subgroup.zpowers S.Q)
    (haaa : tateAOf S c₁ * tateAOf S c₂ * tateAOf S c₃ = S.q)
    (hn1 : c₁ ≠ 1) (hn2 : c₂ ≠ 1) (hn3 : c₃ ≠ 1)
    (h12 : c₁ ≠ c₂) (h13 : c₁ ≠ c₃) (h23 : c₂ ≠ c₃)
    (hp12 : c₁ * c₂ ≠ 1) (hp13 : c₁ * c₃ ≠ 1) (hp23 : c₂ * c₃ ≠ 1) :
    tatePhi S hΔ c₁ + tatePhi S hΔ c₂ + tatePhi S hΔ c₃ = 0 := by
  have haw1 : tateAOf S c₁ * tateWOf S c₁ = S.q := tateAOf_mul_tateWOf S c₁
  have haw2 : tateAOf S c₂ * tateWOf S c₂ = S.q := tateAOf_mul_tateWOf S c₂
  have haw3 : tateAOf S c₃ * tateWOf S c₃ = S.q := tateAOf_mul_tateWOf S c₃
  have ha10 : tateAOf S c₁ ≠ 0 := by
    intro h; rw [h, zero_mul] at haw1; exact hq0 haw1.symm
  have ha20 : tateAOf S c₂ ≠ 0 := by
    intro h; rw [h, zero_mul] at haw2; exact hq0 haw2.symm
  have ha30 : tateAOf S c₃ ≠ 0 := by
    intro h; rw [h, zero_mul] at haw3; exact hq0 haw3.symm
  have hw1 : tateWOf S c₁ = tateAOf S c₂ * tateAOf S c₃ := by
    refine mul_left_cancel₀ ha10 ?_
    rw [haw1, ← haaa]; ring
  have hw2 : tateWOf S c₂ = tateAOf S c₁ * tateAOf S c₃ := by
    refine mul_left_cancel₀ ha20 ?_
    rw [haw2, ← haaa]; ring
  have hw3 : tateWOf S c₃ = tateAOf S c₁ * tateAOf S c₂ := by
    refine mul_left_cancel₀ ha30 ?_
    rw [haw3, ← haaa]; ring
  have hvw : tateAOf S c₂ * tateAOf S c₃ ∈ I := hw1 ▸ tateWOf_mem S c₁
  have huw : tateAOf S c₁ * tateAOf S c₃ ∈ I := hw2 ▸ tateWOf_mem S c₂
  have huv : tateAOf S c₁ * tateAOf S c₂ ∈ I := hw3 ▸ tateWOf_mem S c₃
  have hx12 : tateXK (tateAOf S c₁) (tateAOf S c₂ * tateAOf S c₃) S.q S.hq
      ≠ (tateXK (tateAOf S c₂) (tateAOf S c₁ * tateAOf S c₃) S.q S.hq : K) := by
    rw [← hw1, ← hw2]
    intro h
    rcases tatePhi_X_eq_imp S hloc hvR hvI hq0 hΔ hn1 hn2 h with h' | h'
    · exact h12 h'
    · exact hp12 h'
  have hx13 : tateXK (tateAOf S c₁) (tateAOf S c₂ * tateAOf S c₃) S.q S.hq
      ≠ (tateXK (tateAOf S c₃) (tateAOf S c₁ * tateAOf S c₂) S.q S.hq : K) := by
    rw [← hw1, ← hw3]
    intro h
    rcases tatePhi_X_eq_imp S hloc hvR hvI hq0 hΔ hn1 hn3 h with h' | h'
    · exact h13 h'
    · exact hp13 h'
  have hx23 : tateXK (tateAOf S c₂) (tateAOf S c₁ * tateAOf S c₃) S.q S.hq
      ≠ (tateXK (tateAOf S c₃) (tateAOf S c₁ * tateAOf S c₂) S.q S.hq : K) := by
    rw [← hw2, ← hw3]
    intro h
    rcases tatePhi_X_eq_imp S hloc hvR hvI hq0 hΔ hn2 hn3 h with h' | h'
    · exact h23 h'
    · exact hp23 h'
  have hkey := tate_points_add_eq_zero_K (K := K) (tateAOf S c₁) (tateAOf S c₂) (tateAOf S c₃)
    S.q S.hq haaa hvw huw huv
    (tateAOf_ne_one S hn1) (tateAOf_ne_one S hn2) (tateAOf_ne_one S hn3) hx12 hx13 hx23
    (nonsingular_tateK _ _ S.q S.hq (by rw [← haaa]; ring) (isUnit_one_sub hvw)
      (tateAOf_ne_one S hn1) hΔ)
    (nonsingular_tateK _ _ S.q S.hq (by rw [← haaa]; ring) (isUnit_one_sub huw)
      (tateAOf_ne_one S hn2) hΔ)
    (nonsingular_tateK _ _ S.q S.hq (by rw [← haaa]; ring) (isUnit_one_sub huv)
      (tateAOf_ne_one S hn3) hΔ)
  rw [tatePhi_eq S hΔ hn1, tatePhi_eq S hΔ hn2, tatePhi_eq S hΔ hn3, tatePtPair,
    tatePtPair, tatePtPair]
  simp only [hw1, hw2, hw3]
  exact hkey

/-! ## ★出典の紐付け(`.src`) -/

def tatePhi_X_eq_imp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——X の一致は c = d か cd = 1)",
    sectionId := "genell-def-3-3" }

def tatePhi_add_add_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——一般の位置での準同型性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
