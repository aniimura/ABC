import ABC3.Found.GaloisRep.CollGroupTwo

/-!
# Galois (G6) 第 296 ブロック —— **★★★★★★★★★★準同型性が一般の位置で閉じた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点

> `c₁c₂c₃ = 1`、3 類が互いに異なり単位元でない —— このとき
> **`Φ(c₁) + Φ(c₂) + Φ(c₃) = 0`**(`tatePhi_add_add_eq_zero_all`)

★★★正規化代表元の積の場合分け(`n = 0, 1, 2`)が**すべて済んだ**。
残るのは退化した場合(`cᵢ = cⱼ` または `cᵢ = 1`)だけである。

| `n` | ブロック |
|---|---|
| `0` | **本ブロック** |
| `1` | 第 293 |
| `2` | 第 294(逆元で `n = 1` に還元) |

## ★★★★★★★★`n = 0` は `P(u) + P(v) = P(uv)` そのもの

`a₁a₂a₃ = 1` なので `a₁a₂ = a₃⁻¹`、すなわち **`(a₁a₂, q a₃)` は類 `c₃⁻¹` の正規化した対**
(第 291 の `tateAOf_inv_unit` がそう言っている)。したがって第 295 の

    P(a₁, ·) + P(a₂, ·) = P(a₁a₂, q a₃)

は **`Φ(c₁) + Φ(c₂) = Φ(c₃⁻¹) = −Φ(c₃)`**、すなわち求める式そのものである。
★★★★★「第 3 の点を入れ替えた向きで書く」という技術的な選択が、
**そのまま `c₃⁻¹` の正規化代表元と一致していた**——偶然ではなく、
基本領域の取り方が両方を決めているからである。

## ★★★★相異性の 3 本

| 組 | 破れると |
|---|---|
| `X(c₁) ≠ X(c₂)` | `c₁ = c₂` または `c₃ = 1` |
| `X(c₁) ≠ X(c₃⁻¹)` | `c₂ = 1` または `c₁ = c₃` |
| `X(c₂) ≠ X(c₃⁻¹)` | `c₁ = 1` または `c₂ = c₃` |

★どれも仮定で除いてある。★★`cᵢcⱼ = 1 ⟺ c_k = 1` の補題(`third_eq_one_of_mul` ほか)
を切り出しておくと、3 本が 1 行ずつで書ける。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tate_points_add_eq_K_two` | ★★★★★★★★★★`P(u) + P(v) = P(uv)`(単元 2 つ) |
| `third_eq_one_of_mul` ほか | ★★`cᵢcⱼ = 1 ⟺ c_k = 1` |
| `tatePhi_add_add_eq_zero_zero` | ★★★★★★★★★★**準同型性(`n = 0`)** |
| `tatePhi_add_add_eq_zero_all` | ★★★★★★★★★★**準同型性(一般の位置、全場合)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

/-! ## ★★`cᵢcⱼ = 1 ⟺ c_k = 1` -/

theorem third_eq_one_of_mul {G : Type} [CommGroup G] {c₁ c₂ c₃ : G}
    (hprod : c₁ * c₂ * c₃ = 1) (h : c₁ * c₂ = 1) : c₃ = 1 := by
  have h2 := hprod
  rw [h] at h2
  exact (one_mul c₃).symm.trans h2

theorem second_eq_one_of_mul {G : Type} [CommGroup G] {c₁ c₂ c₃ : G}
    (hprod : c₁ * c₂ * c₃ = 1) (h : c₁ * c₃ = 1) : c₂ = 1 := by
  have h2 : c₂ * (c₁ * c₃) = 1 := by rw [← hprod, mul_comm c₁ c₂, mul_assoc]
  rw [h] at h2
  exact (mul_one c₂).symm.trans h2

theorem first_eq_one_of_mul {G : Type} [CommGroup G] {c₁ c₂ c₃ : G}
    (hprod : c₁ * c₂ * c₃ = 1) (h : c₂ * c₃ = 1) : c₁ = 1 := by
  have h2 : c₁ * (c₂ * c₃) = 1 := by rw [← hprod, mul_assoc]
  rw [h] at h2
  exact (mul_one c₁).symm.trans h2

/-! ## ★★★★★★★★★★`P(u) + P(v) = P(uv)`(単元 2 つ) -/

section Two

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K]
  [DecidableEq K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**`P(u) + P(v) = P(uv)`(単元 2 つの場合)**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_points_add_eq_K_two (u v w q : R) (hq : q ∈ I) (huvw : u * v * w = q)
    (hvw : v * w ∈ I) (huw : u * w ∈ I) (hw : w ∈ I)
    (hu : algebraMap R K (1 - u) ≠ 0) (hv : algebraMap R K (1 - v) ≠ 0)
    (huv : algebraMap R K (1 - u * v) ≠ 0)
    (h12 : tateXK u (v * w) q hq ≠ (tateXK v (u * w) q hq : K))
    (h13 : tateXK u (v * w) q hq ≠ (tateXK (u * v) w q hq : K))
    (h23 : tateXK v (u * w) q hq ≠ (tateXK (u * v) w q hq : K))
    (n₁ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (tateXK u (v * w) q hq) (tateYK (K := K) u (v * w) q hq))
    (n₂ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (tateXK v (u * w) q hq) (tateYK (K := K) v (u * w) q hq))
    (n₃ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (tateXK (u * v) w q hq) (tateYK (K := K) (u * v) w q hq)) :
    Point.some _ _ n₁ + Point.some _ _ n₂ = Point.some _ _ n₃ := by
  have e1 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₁ = 1 := by
    show algebraMap R K ((tateCurveAt q hq).a₁) = 1
    rw [show (tateCurveAt q hq).a₁ = 1 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_one]
  have e3 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₃ = 0 := by
    show algebraMap R K ((tateCurveAt q hq).a₃) = 0
    rw [show (tateCurveAt q hq).a₃ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  have hneg : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.negY
      (tateXK (u * v) w q hq) (tateYK (K := K) (u * v) w q hq)
      = -tateXK (u * v) w q hq - tateYK (K := K) (u * v) w q hq := by
    rw [WeierstrassCurve.Affine.negY, e1, e3]
    ring
  have hd := collDet_K_eq_zero_two (K := K) u v w q hq huvw hvw huw hw hu hv huv
  rw [← hneg] at hd
  have hkey : Point.some _ _ n₁ + Point.some _ _ n₂
      + (-Point.some (tateXK (u * v) w q hq) (tateYK (K := K) (u * v) w q hq) n₃) = 0 := by
    rw [Point.neg_some]
    exact add_add_eq_zero_of_collDet n₁ n₂ _ hd h12 h13 h23
  have h2 : Point.some _ _ n₁ + Point.some _ _ n₂
      - Point.some (tateXK (u * v) w q hq) (tateYK (K := K) (u * v) w q hq) n₃ = 0 := by
    rw [sub_eq_add_neg]
    exact hkey
  exact sub_eq_zero.1 h2

end Two

/-! ## ★★★★★★★★★★準同型性(`n = 0`) -/

section Zero

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**準同型性(`n = 0` の場合)**——`P(u) + P(v) = P(uv)` そのもの。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_add_add_eq_zero_zero (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c₁ c₂ c₃ : Kˣ ⧸ Subgroup.zpowers S.Q) (hprod : c₁ * c₂ * c₃ = 1)
    (h0 : normRep S.v S.Q S.hQ c₁ * normRep S.v S.Q S.hQ c₂ * normRep S.v S.Q S.hQ c₃ = 1)
    (hn1 : c₁ ≠ 1) (hn2 : c₂ ≠ 1) (hn3 : c₃ ≠ 1)
    (h12 : c₁ ≠ c₂) (h13 : c₁ ≠ c₃) (h23 : c₂ ≠ c₃) :
    tatePhi S hΔ c₁ + tatePhi S hΔ c₂ + tatePhi S hΔ c₃ = 0 := by
  have haaa : tateAOf S c₁ * tateAOf S c₂ * tateAOf S c₃ = 1 := by
    refine S.hinj ?_
    rw [map_mul, map_mul, (tateAOf_spec S c₁).1, (tateAOf_spec S c₂).1, (tateAOf_spec S c₃).1,
      ← Units.val_mul, ← Units.val_mul, h0, map_one]
    rfl
  have hvsum : vAdd S.v (normRep S.v S.Q S.hQ c₁) + vAdd S.v (normRep S.v S.Q S.hQ c₂)
      + vAdd S.v (normRep S.v S.Q S.hQ c₃) = 0 := by
    rw [← vAdd_mul, ← vAdd_mul, h0, vAdd_one]
  have f1 := normRep_nonneg S.v S.Q S.hQ c₁
  have f2 := normRep_nonneg S.v S.Q S.hQ c₂
  have f3 := normRep_nonneg S.v S.Q S.hQ c₃
  have hz3 : vAdd S.v (normRep S.v S.Q S.hQ c₃) = 0 := by omega
  have hu3 : IsUnit (tateAOf S c₃) := by
    rcases hloc (tateAOf S c₃) with h | h
    · exact h
    · exfalso
      have := hvI (normRep S.v S.Q S.hQ c₃) ⟨tateAOf S c₃, h, (tateAOf_spec S c₃).1⟩
      omega
  obtain ⟨hi3a, hi3w⟩ := tateAOf_inv_unit S hvR hq0 c₃ hu3
  have hinv : tateAOf S c₁ * tateAOf S c₂ = Ring.inverse (tateAOf S c₃) := by
    have h := Ring.inverse_mul_cancel (tateAOf S c₃) hu3
    calc tateAOf S c₁ * tateAOf S c₂
        = tateAOf S c₁ * tateAOf S c₂ * (Ring.inverse (tateAOf S c₃) * tateAOf S c₃) := by
          rw [h, mul_one]
      _ = Ring.inverse (tateAOf S c₃) * (tateAOf S c₁ * tateAOf S c₂ * tateAOf S c₃) := by ring
      _ = Ring.inverse (tateAOf S c₃) := by rw [haaa, mul_one]
  have hA : tateAOf S c₃⁻¹ = tateAOf S c₁ * tateAOf S c₂ := by rw [hi3a, hinv]
  have huvw : tateAOf S c₁ * tateAOf S c₂ * (S.q * tateAOf S c₃) = S.q := by
    linear_combination S.q * haaa
  have hqa3 : S.q * tateAOf S c₃ ∈ I := Ideal.mul_mem_right _ _ S.hq
  have hvw : tateAOf S c₂ * (S.q * tateAOf S c₃) ∈ I := Ideal.mul_mem_left _ _ hqa3
  have huw : tateAOf S c₁ * (S.q * tateAOf S c₃) ∈ I := Ideal.mul_mem_left _ _ hqa3
  have hn3' : c₃⁻¹ ≠ 1 := fun h => hn3 (by rw [← inv_inv c₃, h, inv_one])
  have ha1 : tateAOf S c₁ ≠ 0 := by
    intro h
    have h2 := tateAOf_mul_tateWOf S c₁
    rw [h, zero_mul] at h2
    exact hq0 h2.symm
  have ha2 : tateAOf S c₂ ≠ 0 := by
    intro h
    have h2 := tateAOf_mul_tateWOf S c₂
    rw [h, zero_mul] at h2
    exact hq0 h2.symm
  have hw1 : tateWOf S c₁ = tateAOf S c₂ * (S.q * tateAOf S c₃) := by
    refine mul_left_cancel₀ ha1 ?_
    linear_combination tateAOf_mul_tateWOf S c₁ - huvw
  have hw2 : tateWOf S c₂ = tateAOf S c₁ * (S.q * tateAOf S c₃) := by
    refine mul_left_cancel₀ ha2 ?_
    linear_combination tateAOf_mul_tateWOf S c₂ - huvw
  have hne3' : algebraMap R K (1 - (tateAOf S c₁ * tateAOf S c₂)) ≠ 0 := by
    rw [← hA]; exact tateAOf_ne_one S hn3'
  have hx12 : tateXK (tateAOf S c₁) (tateAOf S c₂ * (S.q * tateAOf S c₃)) S.q S.hq
      ≠ (tateXK (tateAOf S c₂) (tateAOf S c₁ * (S.q * tateAOf S c₃)) S.q S.hq : K) := by
    rw [← hw1, ← hw2]
    intro h
    rcases tatePhi_X_eq_imp S hloc hvR hvI hq0 hΔ hn1 hn2 h with h' | h'
    · exact h12 h'
    · exact hn3 (third_eq_one_of_mul hprod h')
  have hx13 : tateXK (tateAOf S c₁) (tateAOf S c₂ * (S.q * tateAOf S c₃)) S.q S.hq
      ≠ (tateXK (tateAOf S c₁ * tateAOf S c₂) (S.q * tateAOf S c₃) S.q S.hq : K) := by
    rw [← hw1, ← hA, ← hi3w]
    intro h
    rcases tatePhi_X_eq_imp S hloc hvR hvI hq0 hΔ hn1 hn3' h with h' | h'
    · exact hn2 (second_eq_one_of_mul hprod (by rw [h']; exact inv_mul_cancel c₃))
    · exact h13 ((eq_inv_of_mul_eq_one_left h').trans (inv_inv c₃))
  have hx23 : tateXK (tateAOf S c₂) (tateAOf S c₁ * (S.q * tateAOf S c₃)) S.q S.hq
      ≠ (tateXK (tateAOf S c₁ * tateAOf S c₂) (S.q * tateAOf S c₃) S.q S.hq : K) := by
    rw [← hw2, ← hA, ← hi3w]
    intro h
    rcases tatePhi_X_eq_imp S hloc hvR hvI hq0 hΔ hn2 hn3' h with h' | h'
    · exact hn1 (first_eq_one_of_mul hprod (by rw [h']; exact inv_mul_cancel c₃))
    · exact h23 ((eq_inv_of_mul_eq_one_left h').trans (inv_inv c₃))
  have hkey := tate_points_add_eq_K_two (K := K) (tateAOf S c₁) (tateAOf S c₂)
    (S.q * tateAOf S c₃) S.q S.hq huvw hvw huw hqa3
    (tateAOf_ne_one S hn1) (tateAOf_ne_one S hn2) hne3' hx12 hx13 hx23
    (nonsingular_tateK _ _ _ _ (by linear_combination huvw) (isUnit_one_sub hvw)
      (tateAOf_ne_one S hn1) hΔ)
    (nonsingular_tateK _ _ _ _ (by linear_combination huvw) (isUnit_one_sub huw)
      (tateAOf_ne_one S hn2) hΔ)
    (nonsingular_tateK _ _ _ _ huvw (isUnit_one_sub hqa3) hne3' hΔ)
  have hphi3 : tatePhi S hΔ c₃ = -tatePhi S hΔ c₃⁻¹ := by
    have h := tatePhi_inv S hloc hvR hvI hq0 hΔ c₃⁻¹
    exact (congrArg (tatePhi S hΔ) (inv_inv c₃)).symm.trans h
  rw [tatePhi_eq S hΔ hn1, tatePhi_eq S hΔ hn2, hphi3, tatePhi_eq S hΔ hn3',
    tatePtPair, tatePtPair, tatePtPair]
  simp only [hw1, hw2, hA, hi3w]
  rw [hkey]
  exact add_neg_cancel _

/-! ## ★★★★★★★★★★準同型性(一般の位置、全場合) -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**準同型性(一般の位置、全場合)**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_add_add_eq_zero_all (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c₁ c₂ c₃ : Kˣ ⧸ Subgroup.zpowers S.Q) (hprod : c₁ * c₂ * c₃ = 1)
    (hn1 : c₁ ≠ 1) (hn2 : c₂ ≠ 1) (hn3 : c₃ ≠ 1)
    (h12 : c₁ ≠ c₂) (h13 : c₁ ≠ c₃) (h23 : c₂ ≠ c₃) :
    tatePhi S hΔ c₁ + tatePhi S hΔ c₂ + tatePhi S hΔ c₃ = 0 := by
  obtain ⟨n, hge, hle, hn⟩ := normRep_prod_zpow S c₁ c₂ c₃ hprod
  interval_cases n
  · exact tatePhi_add_add_eq_zero_zero S hloc hvR hvI hq0 hΔ c₁ c₂ c₃ hprod
      (by simpa using hn) hn1 hn2 hn3 h12 h13 h23
  · exact tatePhi_add_add_eq_zero' S hloc hvR hvI hq0 hΔ c₁ c₂ c₃ hprod
      (by simpa using hn) hn1 hn2 hn3 h12 h13 h23
  · exact tatePhi_add_add_eq_zero_two S hloc hvR hvI hq0 hΔ c₁ c₂ c₃ hprod hn
      hn1 hn2 hn3 h12 h13 h23

end Zero

/-! ## ★出典の紐付け(`.src`) -/

def tate_points_add_eq_K_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——P(u) + P(v) = P(uv))",
    sectionId := "genell-def-3-3" }

def tatePhi_add_add_eq_zero_all.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——準同型性(一般の位置、全場合))",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
