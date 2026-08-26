import ABC3.Found.GaloisRep.TatePhiZero

/-!
# Galois (G6) 第 297 ブロック —— **★★★★★★★★★★残るのは倍化だけ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点

> `c₁c₂c₃ = 1` で **3 類が互いに異なれば** `Φ(c₁) + Φ(c₂) + Φ(c₃) = 0`
> (`tatePhi_add_add_eq_zero_of_distinct`)

★★★単位元を含む場合も込みである。**残るのは `cᵢ = cⱼ`(倍化)だけ**になった。

## ★★★★★★単位元を含む場合は反転則だけ

`c₃ = 1` なら `c₁c₂ = 1` すなわち `c₂ = c₁⁻¹` なので

    Φ(c₁) + Φ(c₁⁻¹) + Φ(1) = Φ(c₁) − Φ(c₁) + 0 = 0

★共線性も群法則も要らない——**第 291 の反転則と第 284 の `Φ(1) = 0` だけ**である。

## ★★★2 位の点も自動で片づく

`c₁ = c₂` かつ `c₁² = 1` なら `c₃ = c₁⁻² = 1` なので、上の場合に落ちる。
★★★★したがって**2 位の類は倍化の例外にならない**——`Φ(c) = Φ(c⁻¹) = −Φ(c)` が
そのまま `2Φ(c) = 0` を与えている。

## ★★残っているもの——真の倍化

`c₁ = c₂ = c`(`c² ≠ 1`)の場合。補助母数 `e` を取って

    Φ(c) + Φ(e) = Φ(ce)、  Φ(ce) + Φ(c) = Φ(c²e)、  Φ(c²) + Φ(e) = Φ(c²e)

の 3 本(いずれも一般の位置)から `Φ(c²) = 2Φ(c)` を出す。
★`e` は有限個の類を避ければよく、`Kˣ/q^ℤ` は無限(`1 + 𝔪` を含む)なので取れる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tatePhi_add_add_eq_zero_of_one` | ★★★★★★単位元を含む 3 つ組 |
| `tatePhi_add_add_eq_zero_of_distinct` | ★★★★★★★★★★**3 類が互いに異なるとき** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

/-! ## ★★★★★★単位元を含む 3 つ組 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**単位元を含む 3 つ組**——反転則だけで済む。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_add_add_eq_zero_of_one (S : TateSetup R I K)
    (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c₁ c₂ c₃ : Kˣ ⧸ Subgroup.zpowers S.Q) (hprod : c₁ * c₂ * c₃ = 1)
    (hone : c₁ = 1 ∨ c₂ = 1 ∨ c₃ = 1) :
    tatePhi S hΔ c₁ + tatePhi S hΔ c₂ + tatePhi S hΔ c₃ = 0 := by
  rcases hone with h | h | h
  · subst h
    have h2 : c₂ * c₃ = 1 := (congrArg (· * c₃) (one_mul c₂)).symm.trans hprod
    have h3 : c₃ = c₂⁻¹ := eq_inv_of_mul_eq_one_right h2
    rw [tatePhi_one, h3, tatePhi_inv S hloc hvR hvI hq0 hΔ c₂]
    abel
  · subst h
    have h2 : c₁ * c₃ = 1 := (congrArg (· * c₃) (mul_one c₁)).symm.trans hprod
    have h3 : c₃ = c₁⁻¹ := eq_inv_of_mul_eq_one_right h2
    rw [tatePhi_one, h3, tatePhi_inv S hloc hvR hvI hq0 hΔ c₁]
    abel
  · subst h
    have h2 : c₁ * c₂ = 1 := (mul_one (c₁ * c₂)).symm.trans hprod
    have h3 : c₂ = c₁⁻¹ := eq_inv_of_mul_eq_one_right h2
    rw [tatePhi_one, h3, tatePhi_inv S hloc hvR hvI hq0 hΔ c₁]
    abel

/-! ## ★★★★★★★★★★3 類が互いに異なるとき -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**準同型性(3 類が互いに異なるとき)**——単位元の場合も込み。

★残るのは `cᵢ = cⱼ`(倍化)だけである。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_add_add_eq_zero_of_distinct (S : TateSetup R I K)
    (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c₁ c₂ c₃ : Kˣ ⧸ Subgroup.zpowers S.Q) (hprod : c₁ * c₂ * c₃ = 1)
    (h12 : c₁ ≠ c₂) (h13 : c₁ ≠ c₃) (h23 : c₂ ≠ c₃) :
    tatePhi S hΔ c₁ + tatePhi S hΔ c₂ + tatePhi S hΔ c₃ = 0 := by
  by_cases hn1 : c₁ = 1
  · exact tatePhi_add_add_eq_zero_of_one S hloc hvR hvI hq0 hΔ c₁ c₂ c₃ hprod (Or.inl hn1)
  by_cases hn2 : c₂ = 1
  · exact tatePhi_add_add_eq_zero_of_one S hloc hvR hvI hq0 hΔ c₁ c₂ c₃ hprod
      (Or.inr (Or.inl hn2))
  by_cases hn3 : c₃ = 1
  · exact tatePhi_add_add_eq_zero_of_one S hloc hvR hvI hq0 hΔ c₁ c₂ c₃ hprod
      (Or.inr (Or.inr hn3))
  exact tatePhi_add_add_eq_zero_all S hloc hvR hvI hq0 hΔ c₁ c₂ c₃ hprod hn1 hn2 hn3 h12 h13 h23

/-! ## ★出典の紐付け(`.src`) -/

def tatePhi_add_add_eq_zero_of_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——単位元を含む 3 つ組)",
    sectionId := "genell-def-3-3" }

def tatePhi_add_add_eq_zero_of_distinct.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——3 類が互いに異なるときの準同型性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
