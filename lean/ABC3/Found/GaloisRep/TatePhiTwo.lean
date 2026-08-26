import ABC3.Found.GaloisRep.TatePhiProd

/-!
# Galois (G6) 第 294 ブロック —— **★★★★★★★★★★`n = 2` は逆元で `n = 1` に還る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点

> 正規化代表元の積が `Q²` でも **`Φ(c₁) + Φ(c₂) + Φ(c₃) = 0`**
> (`tatePhi_add_add_eq_zero_two`)

★★★これで準同型性の 3 通り(`n = 0, 1, 2`)のうち **2 通りが済んだ**。

## ★★★★★★★`n = 2` なら 3 つとも環帯

`Σ v(uᵢ) = 2v(Q)` で各 `v(uᵢ) < v(Q)`。★もし `v(u₁) = 0` なら残り 2 つの和が `2v(Q)` で
各々 `< v(Q)` だから和は `< 2v(Q)`——矛盾。したがって **3 つとも `v(uᵢ) > 0`**。

★★★これが効く:逆元の代表元は `v(uᵢ) > 0` のとき `q/uᵢ`(相方)で、その付値は
`v(Q) − v(uᵢ)`。★★★★したがって

    Σ v(normRep cᵢ⁻¹) = 3v(Q) − 2v(Q) = v(Q)

すなわち逆元の 3 つ組は **`n = 1` の場合**である。

## ★★★★★★逆元の代表元の付値

`tateAOf_inv_annulus`(第 291)が `a(c⁻¹) = w(c)` を与えるので、`u·u_w = Q` から
`v(normRep c⁻¹) = v(Q) − v(normRep c)`(`vAdd_normRep_inv`)。
★第 291 で代表元そのものを同定しておいたので、付値の計算は 3 行で済む。

## ★★実例のダイヤモンド、再び

`c₁⁻¹ = c₂⁻¹` から `c₁ = c₂` を出すのに `rw [inv_inv]` は当たらない。
`exact inv_inv c₂` と項の水準で書く(第 293 と同じ)。
★点の群 `Point` は `AddCommGroup` なので `linear_combination` は使えない——`abel` を使う。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateAOf_mem_of_pos` | ★★★★★★`v(u) > 0` なら母数は `I` の元 |
| `vAdd_normRep_inv` | ★★★★★★★逆元の代表元の付値 |
| `tatePhi_add_add_eq_zero_two` | ★★★★★★★★★★**準同型性(`n = 2` の場合)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-! ## ★★★★★★逆元の代表元の付値 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★`v(u) > 0` なら母数は `I` の元。 -/
theorem tateAOf_mem_of_pos (S : TateSetup R I K) (c : Kˣ ⧸ Subgroup.zpowers S.Q)
    (hpos : 0 < vAdd S.v (normRep S.v S.Q S.hQ c)) : tateAOf S c ∈ I := by
  obtain ⟨r, hrI, hr⟩ := S.hmem (normRep S.v S.Q S.hQ c) hpos
  have h : r = tateAOf S c := S.hinj (by rw [hr, (tateAOf_spec S c).1])
  rwa [h] at hrI

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**逆元の代表元の付値**(環帯の場合)。 -/
theorem vAdd_normRep_inv (S : TateSetup R I K)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0) (c : Kˣ ⧸ Subgroup.zpowers S.Q)
    (hpos : 0 < vAdd S.v (normRep S.v S.Q S.hQ c)) :
    vAdd S.v (normRep S.v S.Q S.hQ c⁻¹)
      = vAdd S.v S.Q - vAdd S.v (normRep S.v S.Q S.hQ c) := by
  have ha : tateAOf S c ∈ I := tateAOf_mem_of_pos S c hpos
  obtain ⟨h1, _⟩ := tateAOf_inv_annulus S hvR hvI hq0 c ha
  have haw : tateAOf S c * tateWOf S c = S.q := tateAOf_mul_tateWOf S c
  have ha0 : tateAOf S c ≠ 0 := by
    intro h; rw [h, zero_mul] at haw; exact hq0 haw.symm
  have hw0 : tateWOf S c ≠ 0 := by
    intro h; rw [h, mul_zero] at haw; exact hq0 haw.symm
  have hAK : algebraMap R K (tateAOf S c) ≠ 0 := fun h => ha0 (S.hinj (by rw [h, map_zero]))
  have hWK : algebraMap R K (tateWOf S c) ≠ 0 := fun h => hw0 (S.hinj (by rw [h, map_zero]))
  set ua : Kˣ := Units.mk0 (algebraMap R K (tateAOf S c)) hAK with hua
  set uw : Kˣ := Units.mk0 (algebraMap R K (tateWOf S c)) hWK with huw
  have hna : normRep S.v S.Q S.hQ c = ua := (Units.ext (tateAOf_spec S c).1).symm
  have hnw : normRep S.v S.Q S.hQ c⁻¹ = uw := by
    refine (Units.ext ?_).symm
    show algebraMap R K (tateWOf S c) = (normRep S.v S.Q S.hQ c⁻¹ : K)
    rw [← h1, (tateAOf_spec S c⁻¹).1]
  have huaw : ua * uw = S.Q := by
    refine Units.ext ?_
    show algebraMap R K (tateAOf S c) * algebraMap R K (tateWOf S c) = (S.Q : K)
    rw [← map_mul, haw, S.hQq]
  have hsum : vAdd S.v S.Q = vAdd S.v ua + vAdd S.v uw := by rw [← huaw, vAdd_mul]
  rw [hna, hnw]
  omega

/-! ## ★★★★★★★★★★準同型性(`n = 2` の場合) -/

section Two

variable [DecidableEq K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**準同型性(`n = 2` の場合)**——逆元を取って `n = 1` に還元。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_add_add_eq_zero_two (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c₁ c₂ c₃ : Kˣ ⧸ Subgroup.zpowers S.Q) (hprod : c₁ * c₂ * c₃ = 1)
    (h2 : normRep S.v S.Q S.hQ c₁ * normRep S.v S.Q S.hQ c₂ * normRep S.v S.Q S.hQ c₃
      = S.Q ^ (2 : ℤ))
    (hn1 : c₁ ≠ 1) (hn2 : c₂ ≠ 1) (hn3 : c₃ ≠ 1)
    (h12 : c₁ ≠ c₂) (h13 : c₁ ≠ c₃) (h23 : c₂ ≠ c₃) :
    tatePhi S hΔ c₁ + tatePhi S hΔ c₂ + tatePhi S hΔ c₃ = 0 := by
  have hvsum : vAdd S.v (normRep S.v S.Q S.hQ c₁) + vAdd S.v (normRep S.v S.Q S.hQ c₂)
      + vAdd S.v (normRep S.v S.Q S.hQ c₃) = 2 * vAdd S.v S.Q := by
    rw [← vAdd_mul, ← vAdd_mul, h2, vAdd_zpow]
  have g1 := normRep_lt S.v S.Q S.hQ c₁
  have g2 := normRep_lt S.v S.Q S.hQ c₂
  have g3 := normRep_lt S.v S.Q S.hQ c₃
  have f1 := normRep_nonneg S.v S.Q S.hQ c₁
  have f2 := normRep_nonneg S.v S.Q S.hQ c₂
  have f3 := normRep_nonneg S.v S.Q S.hQ c₃
  have p1 : 0 < vAdd S.v (normRep S.v S.Q S.hQ c₁) := by omega
  have p2 : 0 < vAdd S.v (normRep S.v S.Q S.hQ c₂) := by omega
  have p3 : 0 < vAdd S.v (normRep S.v S.Q S.hQ c₃) := by omega
  have i1 := vAdd_normRep_inv S hvR hvI hq0 c₁ p1
  have i2 := vAdd_normRep_inv S hvR hvI hq0 c₂ p2
  have i3 := vAdd_normRep_inv S hvR hvI hq0 c₃ p3
  have hprodinv : c₁⁻¹ * c₂⁻¹ * c₃⁻¹ = 1 := by
    rw [← mul_inv, ← mul_inv, hprod, inv_one]
  obtain ⟨n, _, _, hn⟩ := normRep_prod_zpow S c₁⁻¹ c₂⁻¹ c₃⁻¹ hprodinv
  have hvn : vAdd S.v (normRep S.v S.Q S.hQ c₁⁻¹) + vAdd S.v (normRep S.v S.Q S.hQ c₂⁻¹)
      + vAdd S.v (normRep S.v S.Q S.hQ c₃⁻¹) = n * vAdd S.v S.Q := by
    rw [← vAdd_mul, ← vAdd_mul, hn, vAdd_zpow]
  have hQ := S.hQ
  have hn1' : n = 1 := by
    rw [i1, i2, i3] at hvn
    have hz : (n - 1) * vAdd S.v S.Q = 0 := by linarith
    rcases mul_eq_zero.1 hz with h | h
    · omega
    · omega
  rw [hn1'] at hn
  have hkey := tatePhi_add_add_eq_zero' S hloc hvR hvI hq0 hΔ c₁⁻¹ c₂⁻¹ c₃⁻¹ hprodinv
    (by simpa using hn)
    (fun h => hn1 (by rw [← inv_inv c₁, h, inv_one]))
    (fun h => hn2 (by rw [← inv_inv c₂, h, inv_one]))
    (fun h => hn3 (by rw [← inv_inv c₃, h, inv_one]))
    (fun h => h12 (by rw [← inv_inv c₁, h]; exact inv_inv c₂))
    (fun h => h13 (by rw [← inv_inv c₁, h]; exact inv_inv c₃))
    (fun h => h23 (by rw [← inv_inv c₂, h]; exact inv_inv c₃))
  rw [tatePhi_inv S hloc hvR hvI hq0 hΔ c₁, tatePhi_inv S hloc hvR hvI hq0 hΔ c₂,
    tatePhi_inv S hloc hvR hvI hq0 hΔ c₃] at hkey
  have h3 : -(tatePhi S hΔ c₁ + tatePhi S hΔ c₂ + tatePhi S hΔ c₃) = 0 := by
    rw [← hkey]; abel
  exact neg_eq_zero.1 h3

end Two

/-! ## ★出典の紐付け(`.src`) -/

def vAdd_normRep_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——逆元の代表元の付値)",
    sectionId := "genell-def-3-3" }

def tatePhi_add_add_eq_zero_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——準同型性(n = 2 の場合))",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
