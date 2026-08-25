import ABC3.Found.GaloisRep.TateAux

/-!
# Galois (G6) 第 300 ブロック —— **★★★★★★★★★★★Tate 一意化(仮定なしの倍化)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★★到達点

> `Φ(c²) = 2Φ(c)` が**仮定なしで**成り立つ(`tatePhi_doubling_all`)。
> したがって `Φ` は**仮定なしの加法同型** `Kˣ/q^ℤ ≃+ E_q(K)`(`tatePhiAddEquivAll`)。

★★★第 298 で唯一残っていた仮定 `hdbl` が、第 299 の補助母数で消えた。

## ★★★★★★★★「避けるべき類」は 9 本の条件から機械的に読める

倍化には 3 本の一般の位置の関係(第 298)が要る。その 9 つの条件
`t11, …, t33` はいずれも **`e` か `e²` が特定の類に等しい**という形にほどける:

| 破れる条件 | ほどくと | 避け先 |
|---|---|---|
| `c = (ce)⁻¹` | `e = (c²)⁻¹` | `A` |
| `e = (ce)⁻¹` | `e² = c⁻¹` | `B` |
| `c = ce` | `e = 1` | `A` |
| `c = (c(ce))⁻¹` | `e = (c³)⁻¹` | `A` |
| `ce = (c(ce))⁻¹` | `e² = (c³)⁻¹` | `B` |
| `c² = (c²e)⁻¹` | `e = (c⁴)⁻¹` | `A` |
| `e = (c²e)⁻¹` | `e² = (c²)⁻¹` | `B` |

★★★★**`e` について 6 本、`e²` について 3 本**——第 299 の `exists_aux` の
`A.card ≤ 6`・`B.card ≤ 3` は、この表を数えて決めた数である。

## ★★★★★可換群の側は抽象のまま解く

`Kˣ ⧸ Subgroup.zpowers Q` の上では `simp`/`rw` が具体化の壁に当たる
(`tools\lean-idioms.md` の項)。そこで**抽象の可換群 `G` で補題を作り**、
商群には項の水準で当てる。★`mul_comm c e` は商群でもそのまま型が付く。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `aux_e_eq_one` 他 6 本 | ★★★★★退けるべき類の同定(抽象可換群) |
| `tatePhi_doubling_all` | ★★★★★★★★★★**仮定なしの倍化** |
| `tatePhi_mul_all` | ★★★★★★★★★★**仮定なしの準同型性** |
| `tatePhiAddEquivAll` | ★★★★★★★★★★★**`Kˣ/q^ℤ ≃+ E_q(K)`** |
| `tate_uniformization_all` | ★★★★★★★★★★★**Tate 一意化(界面の向き)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

/-! ## ★★★★★退けるべき類の同定(抽象可換群) -/

section GroupAux

variable {G : Type} [CommGroup G]

/-- ★★★★★`c = ce` なら `e = 1`。 -/
theorem aux_e_eq_one (c e : G) (h : c = c * e) : e = 1 := by
  have h1 : c * 1 = c * e := (mul_one c).trans h
  exact (mul_left_cancel h1).symm

/-- ★★★★★`c = (ce)⁻¹` なら `e = (c²)⁻¹`。 -/
theorem aux_e_eq_c2 (c e : G) (h : c = (c * e)⁻¹) : e = (c * c)⁻¹ := by
  have h1 : c * (c * e) = 1 := mul_eq_one_iff_eq_inv.2 h
  refine eq_inv_of_mul_eq_one_right ?_
  simpa [mul_comm, mul_left_comm, mul_assoc] using h1

/-- ★★★★★`e = (ce)⁻¹` なら `e² = c⁻¹`。 -/
theorem aux_e2_eq_c (c e : G) (h : e = (c * e)⁻¹) : e * e = c⁻¹ := by
  have h1 : e * (c * e) = 1 := mul_eq_one_iff_eq_inv.2 h
  refine eq_inv_of_mul_eq_one_right ?_
  simpa [mul_comm, mul_left_comm, mul_assoc] using h1

/-- ★★★★★`c = (c(ce))⁻¹` なら `e = (c³)⁻¹`。 -/
theorem aux_e_eq_c3 (c e : G) (h : c = (c * (c * e))⁻¹) : e = (c * c * c)⁻¹ := by
  have h1 : c * (c * (c * e)) = 1 := mul_eq_one_iff_eq_inv.2 h
  refine eq_inv_of_mul_eq_one_right ?_
  simpa [mul_comm, mul_left_comm, mul_assoc] using h1

/-- ★★★★★`ce = (c(ce))⁻¹` なら `e² = (c³)⁻¹`。 -/
theorem aux_e2_eq_c3 (c e : G) (h : c * e = (c * (c * e))⁻¹) : e * e = (c * c * c)⁻¹ := by
  have h1 : (c * e) * (c * (c * e)) = 1 := mul_eq_one_iff_eq_inv.2 h
  refine eq_inv_of_mul_eq_one_right ?_
  simpa [mul_comm, mul_left_comm, mul_assoc] using h1

/-- ★★★★★`c² = (c²e)⁻¹` なら `e = (c⁴)⁻¹`。 -/
theorem aux_e_eq_c4 (c e : G) (h : c * c = (c * c * e)⁻¹) : e = (c * c * c * c)⁻¹ := by
  have h1 : (c * c) * (c * c * e) = 1 := mul_eq_one_iff_eq_inv.2 h
  refine eq_inv_of_mul_eq_one_right ?_
  simpa [mul_comm, mul_left_comm, mul_assoc] using h1

/-- ★★★★★`e = (c²e)⁻¹` なら `e² = (c²)⁻¹`。 -/
theorem aux_e2_eq_c2 (c e : G) (h : e = (c * c * e)⁻¹) : e * e = (c * c)⁻¹ := by
  have h1 : e * (c * c * e) = 1 := mul_eq_one_iff_eq_inv.2 h
  refine eq_inv_of_mul_eq_one_right ?_
  simpa [mul_comm, mul_left_comm, mul_assoc] using h1

end GroupAux

/-! ## ★★濃度の上からの評価 -/

section Card

variable {G : Type} [DecidableEq G]

/-- ★★6 つ並べた `Finset` の濃度は高々 6。 -/
theorem card_le_six (a b c d e f : G) : ({a, b, c, d, e, f} : Finset G).card ≤ 6 := by
  have h1 : ({a, b, c, d, e, f} : Finset G).card ≤ ({b, c, d, e, f} : Finset G).card + 1 :=
    Finset.card_insert_le _ _
  have h2 : ({b, c, d, e, f} : Finset G).card ≤ ({c, d, e, f} : Finset G).card + 1 :=
    Finset.card_insert_le _ _
  have h3 : ({c, d, e, f} : Finset G).card ≤ ({d, e, f} : Finset G).card + 1 :=
    Finset.card_insert_le _ _
  have h4 : ({d, e, f} : Finset G).card ≤ ({e, f} : Finset G).card + 1 :=
    Finset.card_insert_le _ _
  have h5 : ({e, f} : Finset G).card ≤ ({f} : Finset G).card + 1 := Finset.card_insert_le _ _
  have h6 : ({f} : Finset G).card = 1 := Finset.card_singleton f
  omega

/-- ★★3 つ並べた `Finset` の濃度は高々 3。 -/
theorem card_le_three (a b c : G) : ({a, b, c} : Finset G).card ≤ 3 := by
  have h1 : ({a, b, c} : Finset G).card ≤ ({b, c} : Finset G).card + 1 := Finset.card_insert_le _ _
  have h2 : ({b, c} : Finset G).card ≤ ({c} : Finset G).card + 1 := Finset.card_insert_le _ _
  have h3 : ({c} : Finset G).card = 1 := Finset.card_singleton c
  omega

end Card

/-! ## ★★★★★★★★★★仮定なしの倍化 -/

section Doubling

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**仮定なしの倍化** `Φ(c²) = 2Φ(c)`。

★第 299 の補助母数 `e = [1 + q^{n+1}]` を、6 つの類と 3 つの平方から外して取る。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_doubling_all (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    tatePhi S hΔ (c * c) = tatePhi S hΔ c + tatePhi S hΔ c := by
  classical
  obtain ⟨n, hA, hB⟩ := exists_aux S hvR hq0
    ({c, (c * c)⁻¹, 1, (c * c * c)⁻¹, c * c, (c * c * c * c)⁻¹} :
      Finset (Kˣ ⧸ Subgroup.zpowers S.Q))
    ({c⁻¹, (c * c * c)⁻¹, (c * c)⁻¹} : Finset (Kˣ ⧸ Subgroup.zpowers S.Q))
    (card_le_six _ _ _ _ _ _) (card_le_three _ _ _)
  set e := auxCls S n with he
  simp only [Finset.mem_insert, Finset.mem_singleton, not_or] at hA hB
  obtain ⟨n1, n2, n3, n4, n5, n6⟩ := hA
  obtain ⟨m1, m2, m3⟩ := hB
  exact tatePhi_doubling S hloc hvR hvI hq0 hΔ c e
    (fun h => n1 h.symm)
    (fun h => n2 (aux_e_eq_c2 c e h))
    (fun h => m1 (aux_e2_eq_c c e h))
    (fun h => n3 (aux_e_eq_one c e h))
    (fun h => n4 (aux_e_eq_c3 c e h))
    (fun h => m2 (aux_e2_eq_c3 c e h))
    (fun h => n5 h.symm)
    (fun h => n6 (aux_e_eq_c4 c e h))
    (fun h => m3 (aux_e2_eq_c2 c e h))

/-! ## ★★★★★★★★★★仮定なしの準同型性 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**仮定なしの準同型性** `Φ(ab) = Φa + Φb`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePhi_mul_all (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (a b : Kˣ ⧸ Subgroup.zpowers S.Q) :
    tatePhi S hΔ (a * b) = tatePhi S hΔ a + tatePhi S hΔ b :=
  tatePhi_mul_of_doubling S hloc hvR hvI hq0 hΔ
    (fun x => tatePhi_doubling_all S hloc hvR hvI hq0 hΔ x) a b

/-! ## ★★★★★★★★★★★仮定なしの一意化 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★**`Kˣ/q^ℤ ≃+ E_q(K)`**——倍化の仮定を落とした形。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def tatePhiAddEquivAll (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hI : I.IsPrime)
    (hquad : ∀ (t : K) (b c : R),
      t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t)
    (hvalring : ∀ t : K, t ≠ 0 →
      (∃ r : R, algebraMap R K r = t) ∨ (∃ r ∈ I, algebraMap R K r = t⁻¹))
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point :=
  tatePhiAddEquiv S hloc hI hquad hvalring hvR hvI hq0 hΔ
    (fun x => tatePhi_doubling_all S hloc hvR hvI hq0 hΔ x)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★**Tate 一意化**(界面の向き、仮定なし)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_uniformization_all (S : TateSetup R I K) (hloc : ∀ x : R, IsUnit x ∨ x ∈ I)
    (hI : I.IsPrime)
    (hquad : ∀ (t : K) (b c : R),
      t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t)
    (hvalring : ∀ t : K, t ≠ 0 →
      (∃ r : R, algebraMap R K r = t) ∨ (∃ r ∈ I, algebraMap R K r = t⁻¹))
    (hvR : ∀ t : Kˣ, (∃ r : R, algebraMap R K r = (t : K)) → 0 ≤ vAdd S.v t)
    (hvI : ∀ t : Kˣ, (∃ r ∈ I, algebraMap R K r = (t : K)) → 0 < vAdd S.v t)
    (hq0 : S.q ≠ 0)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    Nonempty (((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point
      ≃+ Additive (Kˣ ⧸ Subgroup.zpowers S.Q)) :=
  ⟨(tatePhiAddEquivAll S hloc hI hquad hvalring hvR hvI hq0 hΔ).symm⟩

end Doubling

/-! ## ★出典の紐付け(`.src`) -/

def tatePhi_doubling_all.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——仮定なしの倍化)",
    sectionId := "genell-def-3-3" }

def tatePhi_mul_all.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——仮定なしの準同型性)",
    sectionId := "genell-def-3-3" }

def tate_uniformization_all.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——仮定なしの加法同型)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
