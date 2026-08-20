import ABC3.Found.GaloisRep.QTorsion

/-!
# Galois (G6) 第 106 ブロック —— **★★★★`q^ℤ` の基本領域**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★どの代表元で級数を組むか

第 105 ブロックで `X(u,q)` の片側 2 本の尾を作った。★そこには条件が付いていた:

    qᵐu  も  qᵐu⁻¹  も  I^m に入ること(m ≥ 1)

★★これは `u` を勝手に取ると成り立たない——`v(u)` が大きすぎても小さすぎても崩れる。
★★★**`0 ≤ v(u) < v(q)` に正規化すれば必ず成り立つ**:

    v(qᵐu)   = m·v(q) + v(u)   ≥ v(q) > 0
    v(qᵐu⁻¹) = m·v(q) − v(u)   ≥ v(q) − (v(q)−1) = 1 > 0

## ★★★★そしてこれは `Kˣ/q^ℤ` の**基本領域**である

各類にちょうど一つ、正規化された代表元がある(`exists_unique_normalized_rep`)。
★★これは葉 (d)(核が `q^ℤ`)の土台でもある——類の代表を一意に取れるからである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `vAdd` | ★`v(u)` を `ℤ` の側で見る |
| `exists_normalized_rep` | ★★★`u = q^j·u₀`、`0 ≤ v(u₀) < v(q)` |
| `normalized_rep_unique` | ★★★正規化された代表元は一意 |
| `exists_unique_normalized_rep` | ★★★★**基本領域** |
| `vAdd_pos_of_normalized` | ★★**級数の全項の付値が正**(第 105 の条件) |
-/

namespace ABC3.Found.GaloisRep

open QuotientGroup

variable {K : Type} [Field K]

/-! ## ★付値を `ℤ` の側で見る -/

/-- ★`v(u)` を `ℤ` の側で見たもの。 -/
noncomputable def vAdd (v : Kˣ →* Multiplicative ℤ) (u : Kˣ) : ℤ := Multiplicative.toAdd (v u)

theorem vAdd_mul (v : Kˣ →* Multiplicative ℤ) (u w : Kˣ) :
    vAdd v (u * w) = vAdd v u + vAdd v w := by
  simp [vAdd]

theorem vAdd_zpow (v : Kˣ →* Multiplicative ℤ) (u : Kˣ) (n : ℤ) :
    vAdd v (u ^ n) = n * vAdd v u := by
  simp [vAdd, map_zpow, toAdd_zpow]

theorem vAdd_inv (v : Kˣ →* Multiplicative ℤ) (u : Kˣ) : vAdd v u⁻¹ = - vAdd v u := by
  simp [vAdd]

/-! ## ★★★正規化 -/

/-- ★★★**正規化された代表元の存在**——`u = q^j·u₀` で `0 ≤ v(u₀) < v(q)`。

★`j = v(u) / v(q)`(切り捨て除算)、`v(u₀) = v(u) % v(q)` である。 -/
theorem exists_normalized_rep (v : Kˣ →* Multiplicative ℤ) (q : Kˣ)
    (hq : 0 < vAdd v q) (u : Kˣ) :
    ∃ (j : ℤ) (u₀ : Kˣ), u = q ^ j * u₀ ∧ 0 ≤ vAdd v u₀ ∧ vAdd v u₀ < vAdd v q := by
  set h := vAdd v q with hh
  set a := vAdd v u with ha
  have hmod : a - a / h * h = a % h := by
    rw [Int.emod_def]
    ring
  have h0 : 0 ≤ a % h := Int.emod_nonneg a (ne_of_gt hq)
  have h2 : a % h < h := Int.emod_lt_of_pos a hq
  refine ⟨a / h, u * (q ^ (a / h))⁻¹, ?_, ?_, ?_⟩
  · rw [mul_comm (q ^ (a / h)), mul_assoc, inv_mul_cancel, mul_one]
  · rw [vAdd_mul, vAdd_inv, vAdd_zpow, ← ha, ← hh]
    omega
  · rw [vAdd_mul, vAdd_inv, vAdd_zpow, ← ha, ← hh]
    omega

/-- ★★★**正規化された代表元は一意である**。 -/
theorem normalized_rep_unique (v : Kˣ →* Multiplicative ℤ) (q : Kˣ)
    (hq : 0 < vAdd v q) {u w : Kˣ} (k : ℤ) (hk : w = q ^ k * u)
    (hu0 : 0 ≤ vAdd v u) (hu1 : vAdd v u < vAdd v q)
    (hw0 : 0 ≤ vAdd v w) (hw1 : vAdd v w < vAdd v q) : w = u := by
  have hv : vAdd v w = k * vAdd v q + vAdd v u := by
    rw [hk, vAdd_mul, vAdd_zpow]
  have hk0 : k = 0 := by
    rcases lt_trichotomy k 0 with h | h | h
    · nlinarith
    · exact h
    · nlinarith
  rw [hk, hk0, zpow_zero, one_mul]

/-! ## ★★★★基本領域 -/

/-- ★★★★**基本領域**——`Kˣ/q^ℤ` の各類にちょうど一つ
`0 ≤ v(u) < v(q)` なる代表元がある。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_unique_normalized_rep (v : Kˣ →* Multiplicative ℤ) (q : Kˣ)
    (hq : 0 < vAdd v q) (c : Kˣ ⧸ Subgroup.zpowers q) :
    ∃! u : Kˣ, QuotientGroup.mk u = c ∧ 0 ≤ vAdd v u ∧ vAdd v u < vAdd v q := by
  obtain ⟨w, hw⟩ := QuotientGroup.mk_surjective (s := Subgroup.zpowers q) c
  obtain ⟨j, u₀, hju, h0, h1⟩ := exists_normalized_rep v q hq w
  have hmku : QuotientGroup.mk (s := Subgroup.zpowers q) u₀ = c := by
    rw [← hw]
    refine (QuotientGroup.eq (s := Subgroup.zpowers q)).2 ⟨j, ?_⟩
    rw [hju, mul_comm (q ^ j) u₀, ← mul_assoc, inv_mul_cancel, one_mul]
  refine ⟨u₀, ⟨hmku, h0, h1⟩, ?_⟩
  rintro y ⟨hy, hy0, hy1⟩
  have hmk : QuotientGroup.mk (s := Subgroup.zpowers q) y
      = QuotientGroup.mk (s := Subgroup.zpowers q) u₀ := by rw [hy, hmku]
  obtain ⟨k, hkk⟩ := (QuotientGroup.eq (s := Subgroup.zpowers q)).1 hmk
  have hkk' : q ^ k = y⁻¹ * u₀ := hkk
  have hyk : u₀ = q ^ k * y := by
    rw [hkk', mul_assoc, mul_comm u₀ y, ← mul_assoc, inv_mul_cancel, one_mul]
  exact (normalized_rep_unique v q hq k hyk hy0 hy1 h0 h1).symm

/-! ## ★★級数の項の付値 -/

/-- ★★**正規化された代表元では、`m ≥ 1` のとき `qᵐu` も `qᵐu⁻¹` も付値が正**。

★★★これが第 105 ブロックの `tateXtail` / `tateYtail` が使える条件である。 -/
theorem vAdd_pos_of_normalized (v : Kˣ →* Multiplicative ℤ) (q u : Kˣ)
    (hq : 0 < vAdd v q) (hu0 : 0 ≤ vAdd v u) (hu1 : vAdd v u < vAdd v q)
    (m : ℕ) (hm : 1 ≤ m) :
    0 < vAdd v (q ^ (m : ℤ) * u) ∧ 0 < vAdd v (q ^ (m : ℤ) * u⁻¹) := by
  have hm1 : (1 : ℤ) ≤ (m : ℤ) := by exact_mod_cast hm
  refine ⟨?_, ?_⟩
  · rw [vAdd_mul, vAdd_zpow]
    nlinarith
  · rw [vAdd_mul, vAdd_zpow, vAdd_inv]
    nlinarith

/-! ## ★出典の紐付け(`.src`) -/

def exists_unique_normalized_rep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(q^ℤ の基本領域——級数を組む代表元の取り方)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
