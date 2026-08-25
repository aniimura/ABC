import ABC3.Found.GaloisRep.TatePhiOne

/-!
# Galois (G6) 第 298 ブロック —— **★★★★★★★★★★★倍化さえあれば `≃+` が建つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★★到達点

> **倍化 `Φ(c²) = 2Φ(c)` を仮定すれば `Kˣ/q^ℤ ≃+ E_q(K)`**
> (`tatePhiAddEquiv`・`tate_uniformization_of_doubling`)

★★★葉 (b)(c)(d)(e) と準同型性が、**ただ一つの残りの仮定**に集約された。

## ★★★★★★★★倍化だけが残る理由

`Φ(ab) = Φ(a) + Φ(b)` を 3 つ組 `(a, b, (ab)⁻¹)` の共線性から出すには、
3 類が互いに異なることが要る(第 297)。破れるのは

| 破れ方 | 意味 | 還元先 |
|---|---|---|
| `a = b` | 倍化 | ——(残り) |
| `a = (ab)⁻¹` | `b = a⁻²` | `Φ(a²) = 2Φ(a)` |
| `b = (ab)⁻¹` | `a = b⁻²` | `Φ(b²) = 2Φ(b)` |

★★★★下 2 つは**倍化そのものに落ちる**:`b = a⁻²` なら `ab = a⁻¹` なので
求める式は `−Φa = Φa − Φ(a²)`、すなわち倍化である。
★★★★★したがって **3 通りの退化がすべて 1 つの主張に集まる**。

## ★★★★★★倍化は補助母数で還元できる

補助類 `e` を取って 3 本の一般の位置の関係

    Φ(c) + Φ(e) = Φ(ce)、  Φ(c) + Φ(ce) = Φ(c²e)、  Φ(c²) + Φ(e) = Φ(c²e)

を並べると `Φ(c²) + Φ(e) = Φ(c) + Φ(c) + Φ(e)`、すなわち `Φ(c²) = 2Φ(c)`
(`tatePhi_doubling`)。★要るのは `e` が 9 つの相異条件を満たすことだけである。

★★`Kˣ/q^ℤ` は無限(`1 + 𝔪` を含む)なので、有限個の条件を避ける `e` は取れる
——それが最後に残った一点である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tatePhi_mul` | ★★★★★★★★★★`Φ(ab) = Φa + Φb`(一般の位置) |
| `tatePhi_doubling` | ★★★★★★★★★★倍化の補助母数による還元 |
| `tatePhi_mul_of_doubling` | ★★★★★★★★★★**倍化さえあれば準同型性は全部出る** |
| `tatePhiAddEquiv` | ★★★★★★★★★★★**`Kˣ/q^ℤ ≃+ E_q(K)`** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

/-! ## ★★★★★★★★★★`Φ(ab) = Φa + Φb`(一般の位置) -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**`Φ(ab) = Φ(a) + Φ(b)`(一般の位置)**。 -/
theorem tatePhi_mul (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (a b : Kˣ ⧸ Subgroup.zpowers S.Q) (h1 : a ≠ b) (h2 : a ≠ (a * b)⁻¹)
    (h3 : b ≠ (a * b)⁻¹) :
    tatePhi S hΔ (a * b) = tatePhi S hΔ a + tatePhi S hΔ b := by
  have hprod : a * b * (a * b)⁻¹ = 1 := mul_inv_cancel _
  have hk := tatePhi_add_add_eq_zero_of_distinct S hloc hvR hvI hq0 hΔ a b (a * b)⁻¹
    hprod h1 h2 h3
  rw [tatePhi_inv S hloc hvR hvI hq0 hΔ (a * b)] at hk
  have h4 : tatePhi S hΔ a + tatePhi S hΔ b - tatePhi S hΔ (a * b) = 0 := by
    rw [← hk]; abel
  exact (sub_eq_zero.1 h4).symm

/-! ## ★★★★★★★★★★倍化の補助母数による還元 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**倍化——補助母数 `e` による還元**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_doubling (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c e : Kˣ ⧸ Subgroup.zpowers S.Q)
    (t11 : c ≠ e) (t12 : c ≠ (c * e)⁻¹) (t13 : e ≠ (c * e)⁻¹)
    (t21 : c ≠ c * e) (t22 : c ≠ (c * (c * e))⁻¹) (t23 : c * e ≠ (c * (c * e))⁻¹)
    (t31 : c * c ≠ e) (t32 : c * c ≠ (c * c * e)⁻¹) (t33 : e ≠ (c * c * e)⁻¹) :
    tatePhi S hΔ (c * c) = tatePhi S hΔ c + tatePhi S hΔ c := by
  have e1 := tatePhi_mul S hloc hvR hvI hq0 hΔ c e t11 t12 t13
  have e2 := tatePhi_mul S hloc hvR hvI hq0 hΔ c (c * e) t21 t22 t23
  have e3 := tatePhi_mul S hloc hvR hvI hq0 hΔ (c * c) e t31 t32 t33
  have hassoc : c * (c * e) = c * c * e := (mul_assoc c c e).symm
  rw [hassoc] at e2
  have h5 : tatePhi S hΔ (c * c) + tatePhi S hΔ e
      = tatePhi S hΔ c + tatePhi S hΔ c + tatePhi S hΔ e := by
    rw [← e3, e2, e1]
    abel
  exact add_right_cancel h5

/-! ## ★★★★★★★★★★倍化さえあれば準同型性は全部出る -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**倍化さえあれば準同型性は全部出る**。

★3 通りの退化(`a = b`、`a = (ab)⁻¹`、`b = (ab)⁻¹`)がすべて倍化に落ちる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_mul_of_doubling (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (hdbl : ∀ x : Kˣ ⧸ Subgroup.zpowers S.Q, tatePhi S hΔ (x * x)
      = tatePhi S hΔ x + tatePhi S hΔ x)
    (a b : Kˣ ⧸ Subgroup.zpowers S.Q) :
    tatePhi S hΔ (a * b) = tatePhi S hΔ a + tatePhi S hΔ b := by
  by_cases hab : a = b
  · subst hab
    exact hdbl a
  by_cases h2 : a = (a * b)⁻¹
  · have h3 : a * (a * b) = 1 := by
      calc a * (a * b) = (a * b)⁻¹ * (a * b) := by rw [← h2]
        _ = 1 := inv_mul_cancel _
    have h4 : a * a * b = 1 := (mul_assoc a a b) ▸ h3
    have h5 : b = (a * a)⁻¹ := eq_inv_of_mul_eq_one_right h4
    subst h5
    have h6 : a * (a * a)⁻¹ = a⁻¹ :=
      eq_inv_of_mul_eq_one_left ((mul_right_comm a ((a * a)⁻¹) a).trans (mul_inv_cancel _))
    rw [h6, tatePhi_inv S hloc hvR hvI hq0 hΔ a,
      tatePhi_inv S hloc hvR hvI hq0 hΔ (a * a), hdbl a]
    abel
  by_cases h3 : b = (a * b)⁻¹
  · have h4 : b * (a * b) = 1 := by
      calc b * (a * b) = (a * b)⁻¹ * (a * b) := by rw [← h3]
        _ = 1 := inv_mul_cancel _
    have h5 : b * b * a = 1 := ((mul_right_comm b b a).trans (mul_assoc b a b)).trans h4
    have h6 : a = (b * b)⁻¹ := eq_inv_of_mul_eq_one_right h5
    subst h6
    have h7 : (b * b)⁻¹ * b = b⁻¹ :=
      eq_inv_of_mul_eq_one_left ((mul_assoc ((b * b)⁻¹) b b).trans (inv_mul_cancel _))
    rw [h7, tatePhi_inv S hloc hvR hvI hq0 hΔ b,
      tatePhi_inv S hloc hvR hvI hq0 hΔ (b * b), hdbl b]
    abel
  exact tatePhi_mul S hloc hvR hvI hq0 hΔ a b hab h2 h3

/-! ## ★★★★★★★★★★★`Kˣ/q^ℤ ≃+ E_q(K)` -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★**Tate 一意化** `Kˣ/q^ℤ ≃+ E_q(K)`(倍化を仮定)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def tatePhiAddEquiv (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hI : I.IsPrime)
    (hquad : ∀ (t : K) (b c : R),
      t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t)
    (hvalring : ∀ t : K, t ≠ 0 →
      (∃ r : R, algebraMap R K r = t) ∨ (∃ r ∈ I, algebraMap R K r = t⁻¹))
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (hdbl : ∀ x : Kˣ ⧸ Subgroup.zpowers S.Q, tatePhi S hΔ (x * x)
      = tatePhi S hΔ x + tatePhi S hΔ x) :
    Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point :=
  AddEquiv.mk'
    (Equiv.ofBijective (fun c : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) =>
      tatePhi S hΔ (Additive.toMul c))
      ⟨tatePhi_injective S hloc hΔ,
        tatePhi_surjective S hloc hI hquad hvalring hvR hvI hq0 hΔ⟩)
    (fun _ _ => tatePhi_mul_of_doubling S hloc hvR hvI hq0 hΔ hdbl _ _)

/-- ★★★★★★★★★★★**Tate 一意化**(界面の向き、倍化を仮定)。 -/
theorem tate_uniformization_of_doubling (S : TateSetup R I K)
    (hloc : ∀ x : R, IsUnit x ∨ x ∈ I) (hI : I.IsPrime)
    (hquad : ∀ (t : K) (b c : R),
      t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t)
    (hvalring : ∀ t : K, t ≠ 0 →
      (∃ r : R, algebraMap R K r = t) ∨ (∃ r ∈ I, algebraMap R K r = t⁻¹))
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (hdbl : ∀ x : Kˣ ⧸ Subgroup.zpowers S.Q, tatePhi S hΔ (x * x)
      = tatePhi S hΔ x + tatePhi S hΔ x) :
    Nonempty (((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point
      ≃+ Additive (Kˣ ⧸ Subgroup.zpowers S.Q)) :=
  ⟨(tatePhiAddEquiv S hloc hI hquad hvalring hvR hvI hq0 hΔ hdbl).symm⟩

/-! ## ★出典の紐付け(`.src`) -/

def tatePhi_mul_of_doubling.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——倍化から準同型性)",
    sectionId := "genell-def-3-3" }

def tatePhiAddEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——K^x/q^Z ≃+ E_q(K))",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
