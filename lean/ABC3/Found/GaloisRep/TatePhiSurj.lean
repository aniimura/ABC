import ABC3.Found.GaloisRep.TateSurjAll

/-!
# Galois (G6) 第 290 ブロック —— **★★★★★★★★★★`Φ` は全単射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点——集合としての一意化

> `Φ : Kˣ/q^ℤ → E_q(K)` は**全単射**(`tatePhi_bijective`)、したがって
> `Kˣ/q^ℤ ≃ E_q(K)`(`tatePhiEquiv`)

★★★葉 (d)(第 285)と葉 (e)(第 289)が類の水準で合流した。
**残るのは群構造だけ**である。

## ★★★★★★★★正規化した対はある類から来る

点の水準の全射性(第 289)は対 `(a,w)` を返す。それが**類から来る**ことを言うのに
基本領域の条件 `0 ≤ v(u) < v(Q)` が要る:

| 条件 | 根拠 |
|---|---|
| `0 ≤ v(u)` | `a ∈ R` |
| `v(u) < v(Q)` | `v(Q) − v(u) = v(w) > 0`(`w ∈ 𝔪`) |

★★★★**相方 `w` が `𝔪` にいることが、そのまま `u` が基本領域にいることを意味する**。
第 279 で見た「正規化した 3 つ組では相方は必ず `𝔪`」の裏返しである。

## ★★★付値は「`R` から来るか」の言い換えとしてだけ使う

`hvR`・`hvI` は「`R` の元は `v ≥ 0`」「`I` の元は `v > 0`」——`TateSetup` の
`hmem`・`hmem0` の逆向きである。★付値の**順序**は使うが、**ultrametric 不等式は使わない**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_class_of_pair` | ★★★★★★★★正規化した対はある類から来る |
| `tate_equation_of_nonsingular` | ★★★★★★曲線の式を明示した形で |
| `tatePhi_surjective` | ★★★★★★★★★★**`Φ` は全射** |
| `tatePhi_bijective`・`tatePhiEquiv` | ★★★★★★★★★★**`Kˣ/q^ℤ ≃ E_q(K)`** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-! ## ★★★★★★★★正規化した対はある類から来る -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**正規化した対はある類から来る**。

★`v(Q) − v(u) = v(w) > 0` が基本領域の上限を与える。 -/
theorem exists_class_of_pair (S : TateSetup R I K)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0) (a w : R) (haw : a * w = S.q) (hwI : w ∈ I) :
    ∃ c : Kˣ ⧸ Subgroup.zpowers S.Q, tateAOf S c = a ∧ tateWOf S c = w := by
  have ha0 : a ≠ 0 := by
    intro h
    rw [h, zero_mul] at haw
    exact hq0 haw.symm
  have hw0 : w ≠ 0 := by
    intro h
    rw [h, mul_zero] at haw
    exact hq0 haw.symm
  have hAK : algebraMap R K a ≠ 0 := fun h => ha0 (S.hinj (by rw [h, map_zero]))
  have hWK : algebraMap R K w ≠ 0 := fun h => hw0 (S.hinj (by rw [h, map_zero]))
  set u : Kˣ := Units.mk0 (algebraMap R K a) hAK with hu
  set uw : Kˣ := Units.mk0 (algebraMap R K w) hWK with huw
  have huuw : u * uw = S.Q := by
    refine Units.ext ?_
    show algebraMap R K a * algebraMap R K w = (S.Q : K)
    rw [← map_mul, haw, S.hQq]
  have h0 : 0 ≤ vAdd S.v u := hvR u ⟨a, rfl⟩
  have h1 : vAdd S.v u < vAdd S.v S.Q := by
    have hpos : 0 < vAdd S.v uw := hvI uw ⟨w, hwI, rfl⟩
    have hsum : vAdd S.v S.Q = vAdd S.v u + vAdd S.v uw := by rw [← huuw, vAdd_mul]
    omega
  refine ⟨QuotientGroup.mk u, ?_, ?_⟩ <;>
  · obtain ⟨h1', h2'⟩ := pair_eq_of_same_class' S.hinj S.v S.Q S.hQ S.q S.hQq
      (a := tateAOf S (QuotientGroup.mk u)) (w := tateWOf S (QuotientGroup.mk u))
      (a' := a) (w' := w) (u := normRep S.v S.Q S.hQ (QuotientGroup.mk u)) (u' := u)
      (tateAOf_spec S _).1 rfl (tateAOf_mul_tateWOf S _) haw
      (normRep_nonneg S.v S.Q S.hQ _) (normRep_lt S.v S.Q S.hQ _) h0 h1
      (by rw [normRep_mk])
    first
    | exact h1'
    | exact h2'

/-! ## ★★★★★★曲線の式を明示した形で -/

set_option maxHeartbeats 1600000 in
theorem tate_equation_of_nonsingular (q : R) (hq : q ∈ I) (x y : K)
    (n : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular x y) :
    y ^ 2 + x * y = x ^ 3 + algebraMap R K ((tateCurveAt q hq).a₄) * x
      + algebraMap R K ((tateCurveAt q hq).a₆) := by
  have h := n.1
  rw [WeierstrassCurve.Affine.equation_iff] at h
  have e1 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₁ = 1 := by
    show algebraMap R K ((tateCurveAt q hq).a₁) = 1
    rw [show (tateCurveAt q hq).a₁ = 1 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_one]
  have e2 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₂ = 0 := by
    show algebraMap R K ((tateCurveAt q hq).a₂) = 0
    rw [show (tateCurveAt q hq).a₂ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  have e3 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₃ = 0 := by
    show algebraMap R K ((tateCurveAt q hq).a₃) = 0
    rw [show (tateCurveAt q hq).a₃ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  have e4 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₄
      = algebraMap R K ((tateCurveAt q hq).a₄) := rfl
  have e6 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₆
      = algebraMap R K ((tateCurveAt q hq).a₆) := rfl
  rw [e1, e2, e3, e4, e6] at h
  linear_combination h

/-! ## ★★★★★★★★★★`Φ` は全単射 -/

section Bij

variable [DecidableEq K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**`Φ` は全射**——葉 (e) が類の水準で閉じた。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_surjective (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hI : I.IsPrime)
    (hquad : ∀ (t : K) (b c : R),
      t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t)
    (hvalring : ∀ t : K, t ≠ 0 →
      (∃ r : R, algebraMap R K r = t) ∨ (∃ r ∈ I, algebraMap R K r = t⁻¹))
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    Function.Surjective (tatePhi S hΔ) := by
  intro P
  cases P with
  | zero => exact ⟨1, tatePhi_one S hΔ⟩
  | some x y n =>
    have he := tate_equation_of_nonsingular S.q S.hq x y n
    obtain ⟨a, w, haw, hwI, hne, hX, hY⟩ :=
      tate_point_surjective S.hinj hloc hI hquad hvalring S.q S.hq x y he
    obtain ⟨c, hac, hwc⟩ := exists_class_of_pair S hvR hvI hq0 a w haw hwI
    refine ⟨c, ?_⟩
    rw [tatePhi, dif_neg (by rw [hac]; exact hne), tatePtPair]
    simp only [Point.some.injEq]
    rw [hac, hwc]
    exact ⟨hX, hY⟩

/-- ★★★★★★★★★★**`Φ` は全単射**。 -/
theorem tatePhi_bijective (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hI : I.IsPrime)
    (hquad : ∀ (t : K) (b c : R),
      t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t)
    (hvalring : ∀ t : K, t ≠ 0 →
      (∃ r : R, algebraMap R K r = t) ∨ (∃ r ∈ I, algebraMap R K r = t⁻¹))
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    Function.Bijective (tatePhi S hΔ) :=
  ⟨tatePhi_injective S hloc hΔ,
    tatePhi_surjective S hloc hI hquad hvalring hvR hvI hq0 hΔ⟩

/-- ★★★★★★★★★★**集合としての Tate 一意化** `Kˣ/q^ℤ ≃ E_q(K)`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def tatePhiEquiv (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hI : I.IsPrime)
    (hquad : ∀ (t : K) (b c : R),
      t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t)
    (hvalring : ∀ t : K, t ≠ 0 →
      (∃ r : R, algebraMap R K r = t) ∨ (∃ r ∈ I, algebraMap R K r = t⁻¹))
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    (Kˣ ⧸ Subgroup.zpowers S.Q) ≃ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point :=
  Equiv.ofBijective _ (tatePhi_bijective S hloc hI hquad hvalring hvR hvI hq0 hΔ)

end Bij

/-! ## ★出典の紐付け(`.src`) -/

def tatePhi_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Phi は全射)",
    sectionId := "genell-def-3-3" }

def tatePhiEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——集合としての一意化)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
