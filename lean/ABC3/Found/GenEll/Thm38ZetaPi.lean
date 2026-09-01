/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.GroupTheory.OrderOfElement
import Mathlib.GroupTheory.QuotientGroup.Basic
import ABC3.Found.GaloisRep.TatePhi
import Mathlib.RingTheory.RootsOfUnity.PrimitiveRoots
import ABC3.Meta.Claim

/-!
# 第 1161 ブロック —— **`ζ` と `π` は `Lˣ/qℤ` の中で独立**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19–p.20。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`AlphaBridge` の節点 1 の核

`Skeleton/GenEll/AlphaBridge.lean` の節点 1 は

> Tate 一意化の下で `E[l] ⊆ Lˣ/qℤ` は `ζ` と `π`（`πˡ = q`）で生成される

であった。★**その「独立性」の側を無条件に取る**——これが
`α = (1 1 / 0 1)` という**具体的な行列**を出す唯一の道である。

## ★★★★★★★★機構

`Lˣ/⟨Q⟩` の中で `[ζ]` も `[π]` も `l` で殺される（`ζˡ = 1`、`πˡ = Q ≡ 1`）。
☆独立性は次のように出る:

    `ζᵃ·πᵇ = Qⁿ`  を `l` 乗すると  `Qᵇ = Q^{ln}`
    `Q` は無限位数なので  `b = ln`、すなわち  `l ∣ b`
    すると `πᵇ = (πˡ)ⁿ = Qⁿ` なので  `ζᵃ = 1`、`ζ` の位数は `l` なので  `l ∣ a`

★`Q` が無限位数であること（`hQinf`）は Tate 母数が `0 < v(q)` を満たすことから来る。
☆`ζ` の位数がちょうど `l` であることは原始 `l` 乗根であることそのものである。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★★★`ζ` の位数はちょうど `l` -/

/-- ☆`ζˡ = 1` かつ `0 < n < l` で `ζⁿ ≠ 1` なら `ζ` の位数はちょうど `l`。 -/
theorem orderOf_eq_of_primitive {G : Type*} [Group G] {l : ℕ} (hl : 0 < l) {ζ : G}
    (hζl : ζ ^ l = 1) (hζprim : ∀ n : ℕ, 0 < n → n < l → ζ ^ n ≠ 1) :
    orderOf ζ = l := by
  have hdvd : orderOf ζ ∣ l := orderOf_dvd_of_pow_eq_one hζl
  have hne : orderOf ζ ≠ 0 := by
    intro h0
    rw [h0] at hdvd
    exact absurd (Nat.eq_zero_of_zero_dvd hdvd) hl.ne'
  have hpos : 0 < orderOf ζ := Nat.pos_of_ne_zero hne
  have hle : orderOf ζ ≤ l := Nat.le_of_dvd hl hdvd
  rcases lt_or_eq_of_le hle with hlt | heq
  · exact absurd (pow_orderOf_eq_one ζ) (hζprim _ hpos hlt)
  · exact heq

/-! ## ★★★★★★★★★★★★★★★★独立性 -/

/-- ★★★★★★★★★★★★★★★★
**`ζ` と `π` は `Lˣ/⟨Q⟩` の中で独立**——★**無条件**（第 1161）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`ζᵃ·πᵇ ∈ ⟨Q⟩` なら `l ∣ a` かつ `l ∣ b` である。
★★これが `E[l]` の基底を `([ζ], [π])` に取れることの中身であり、
`sigma_acts_as_alpha`（第 993）の指数の動きを行列 `α = (1 1 / 0 1)` として
読むための土台である。 -/
theorem zeta_pi_indep {G : Type*} [CommGroup G] {l : ℕ} (hl : 0 < l) {Q ζ π : G}
    (hζl : ζ ^ l = 1) (hζprim : ∀ n : ℕ, 0 < n → n < l → ζ ^ n ≠ 1)
    (hπl : π ^ l = Q) (hQinf : ∀ j : ℤ, Q ^ j = 1 → j = 0)
    (a b n : ℤ) (h : ζ ^ a * π ^ b = Q ^ n) :
    ((l : ℤ) ∣ a) ∧ ((l : ℤ) ∣ b) := by
  have hζz : ζ ^ (l : ℤ) = 1 := by rw [zpow_natCast]; exact hζl
  have hπz : π ^ (l : ℤ) = Q := by rw [zpow_natCast]; exact hπl
  -- ☆`l` 乗すると `Qᵇ = Q^{ln}`
  have hpow : (ζ ^ a * π ^ b) ^ (l : ℤ) = (Q ^ n) ^ (l : ℤ) := by rw [h]
  rw [mul_zpow, ← zpow_mul, ← zpow_mul, ← zpow_mul,
    mul_comm a (l : ℤ), zpow_mul, hζz, one_zpow, one_mul,
    mul_comm b (l : ℤ), zpow_mul, hπz] at hpow
  -- hpow : Q ^ b = Q ^ (n * l)
  have hb0 : b - n * (l : ℤ) = 0 := by
    refine hQinf _ ?_
    rw [zpow_sub, hpow]
    simp
  have hbn : b = n * (l : ℤ) := sub_eq_zero.mp hb0
  have hb' : b = (l : ℤ) * n := hbn.trans (mul_comm n (l : ℤ))
  -- ☆`πᵇ = Qⁿ` なので `ζᵃ = 1`
  have hπb : π ^ b = Q ^ n := by
    rw [hb', zpow_mul, hπz]
  have hζa : ζ ^ a = 1 := by
    have := h
    rw [hπb] at this
    exact mul_right_cancel (b := Q ^ n) (by rw [this, one_mul])
  have hord : orderOf ζ = l := orderOf_eq_of_primitive hl hζl hζprim
  refine ⟨?_, ⟨n, hb'⟩⟩
  have := orderOf_dvd_iff_zpow_eq_one.2 hζa
  rwa [hord] at this

/-! ## ★★★★★★★★★★★★核はちょうど `lℤ × lℤ` -/

/-- ☆逆向き——`l ∣ a` かつ `l ∣ b` なら `ζᵃπᵇ ∈ ⟨Q⟩`。 -/
theorem zeta_pi_mem_of_dvd {G : Type*} [CommGroup G] {l : ℕ} {Q ζ π : G}
    (hζl : ζ ^ l = 1) (hπl : π ^ l = Q)
    {a b : ℤ} (ha : (l : ℤ) ∣ a) (hb : (l : ℤ) ∣ b) :
    ∃ n : ℤ, ζ ^ a * π ^ b = Q ^ n := by
  obtain ⟨a', rfl⟩ := ha
  obtain ⟨b', rfl⟩ := hb
  refine ⟨b', ?_⟩
  have hζz : ζ ^ (l : ℤ) = 1 := by rw [zpow_natCast]; exact hζl
  have hπz : π ^ (l : ℤ) = Q := by rw [zpow_natCast]; exact hπl
  rw [zpow_mul, hζz, one_zpow, one_mul, zpow_mul, hπz]

/-- ★★★★★★★★★★★★★★★★
**座標の核はちょうど `lℤ × lℤ`**——★**無条件**（第 1162）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`ζᵃ·πᵇ ∈ ⟨Q⟩ ⟺ l ∣ a かつ l ∣ b`。
★★これが「`([ζ], [π])` は `E[l]` の `ℤ/l`-基底である」ことの正確な形であり、
`AlphaBridge` の節点 2 が消費する界面である。 -/
theorem zeta_pi_mem_zpowers_iff {G : Type*} [CommGroup G] {l : ℕ} (hl : 0 < l) {Q ζ π : G}
    (hζl : ζ ^ l = 1) (hζprim : ∀ n : ℕ, 0 < n → n < l → ζ ^ n ≠ 1)
    (hπl : π ^ l = Q) (hQinf : ∀ j : ℤ, Q ^ j = 1 → j = 0)
    (a b : ℤ) :
    (∃ n : ℤ, ζ ^ a * π ^ b = Q ^ n) ↔ ((l : ℤ) ∣ a ∧ (l : ℤ) ∣ b) := by
  constructor
  · rintro ⟨n, hn⟩
    exact zeta_pi_indep hl hζl hζprim hπl hQinf a b n hn
  · rintro ⟨ha, hb⟩
    exact zeta_pi_mem_of_dvd hζl hπl ha hb

/-! ## ★★★★★★★★★★★★同じ類 ⟺ 座標が `l` を法として一致 -/

/-- ★★★★★★★★★★★★★★★★
**`ζᵃπᵇ` の類が一致するのは座標が `l` を法として一致するとき**——★**無条件**（第 1163）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`Lˣ/⟨Q⟩` の中で `[ζ^{a₁}π^{b₁}] = [ζ^{a₂}π^{b₂}] ⟺ a₁ ≡ a₂, b₁ ≡ b₂ (mod l)`。
★★これで座標は `(ZMod l)²` の元として**一意に**決まる。
☆`AlphaBridge` の節点 2 はこの形を直接使う——`σ` の作用が
`(a,b) ↦ (a+b, b)`（`sigma_acts_as_alpha`、第 993）であることと合わせると、
`σ` の行列がちょうど `α = (1 1 / 0 1)` になる。 -/
theorem zeta_pi_coord_eq_iff {G : Type*} [CommGroup G] {l : ℕ} (hl : 0 < l) {Q ζ π : G}
    (hζl : ζ ^ l = 1) (hζprim : ∀ n : ℕ, 0 < n → n < l → ζ ^ n ≠ 1)
    (hπl : π ^ l = Q) (hQinf : ∀ j : ℤ, Q ^ j = 1 → j = 0)
    (a₁ b₁ a₂ b₂ : ℤ) :
    (∃ n : ℤ, ζ ^ a₁ * π ^ b₁ = ζ ^ a₂ * π ^ b₂ * Q ^ n)
      ↔ ((l : ℤ) ∣ (a₁ - a₂) ∧ (l : ℤ) ∣ (b₁ - b₂)) := by
  have key := zeta_pi_mem_zpowers_iff hl hζl hζprim hπl hQinf (a₁ - a₂) (b₁ - b₂)
  have h1 : ζ ^ (a₁ - a₂) * π ^ (b₁ - b₂)
      = (ζ ^ a₁ * π ^ b₁) * (ζ ^ a₂ * π ^ b₂)⁻¹ := by
    rw [zpow_sub, zpow_sub, mul_inv]
    simp only [mul_assoc, mul_comm, mul_left_comm]
  rw [← key]
  constructor
  · rintro ⟨n, hn⟩
    refine ⟨n, ?_⟩
    rw [h1, hn, mul_comm (ζ ^ a₂ * π ^ b₂) (Q ^ n), mul_assoc, mul_inv_cancel, mul_one]
  · rintro ⟨n, hn⟩
    refine ⟨n, ?_⟩
    rw [h1] at hn
    rw [mul_inv_eq_iff_eq_mul.mp hn, mul_comm]

/-! ## ★★★★★★★★★★★★★★★★`σ` の座標での作用は `α = (1 1 / 0 1)` -/

/-- ★★★★★★★★**`σ(ζ) = ζ`・`σ(π) = ζπ` なら `σ(ζᵃπᵇ) = ζ^{a+b}πᵇ`**——★**無条件**。

☆`sigma_acts_as_alpha`（第 993）の**群の水準の形**である
——あちらは体の単数群に特化していた。 -/
theorem sigma_zeta_pi {G : Type*} [CommGroup G] {ζ π : G} (σ : G →* G)
    (hζ : σ ζ = ζ) (hπ : σ π = ζ * π) (a b : ℤ) :
    σ (ζ ^ a * π ^ b) = ζ ^ (a + b) * π ^ b := by
  rw [map_mul, map_zpow, map_zpow, hζ, hπ, mul_zpow, zpow_add, mul_assoc]

/-- ★★★★★★★★★★★★★★★★★★★★
**`σ` の行列はちょうど `α = (1 1 / 0 1)`**——★**無条件**（第 1164）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`σ(ζᵃπᵇ)` を座標 `(a', b')` で書いたなら、必ず
`a' ≡ a + b`、`b' ≡ b`（`mod l`）である。
★★★これが「mod `l` 像が `α` を含む」ことの**座標の側の全体**である
——残るのは Tate 一意化でこの `(ζ, π)` を `E[l]` の基底として実現する段だけである。 -/
theorem sigma_coord_alpha {G : Type*} [CommGroup G] {l : ℕ} (hl : 0 < l) {Q ζ π : G}
    (hζl : ζ ^ l = 1) (hζprim : ∀ n : ℕ, 0 < n → n < l → ζ ^ n ≠ 1)
    (hπl : π ^ l = Q) (hQinf : ∀ j : ℤ, Q ^ j = 1 → j = 0)
    (σ : G →* G) (hζ : σ ζ = ζ) (hπ : σ π = ζ * π)
    (a b a' b' : ℤ)
    (h : ∃ n : ℤ, σ (ζ ^ a * π ^ b) = ζ ^ a' * π ^ b' * Q ^ n) :
    ((l : ℤ) ∣ (a + b - a')) ∧ ((l : ℤ) ∣ (b - b')) := by
  rw [sigma_zeta_pi σ hζ hπ a b] at h
  exact (zeta_pi_coord_eq_iff hl hζl hζprim hπl hQinf (a + b) b a' b').mp h

/-! ## ★★★★★★★★★★★★全射の側——`l`-捩れは `ζ` と `π` で尽きる -/

/-- ★★★★★★★★★★★★★★★★
**`xˡ ∈ ⟨Q⟩` なら `x = ζᵃπᵇ`**——★**無条件**（第 1165）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`xˡ = Qᵐ` なら `(x·π^{−m})ˡ = 1` なので `x·π^{−m} ∈ μ_l = ⟨ζ⟩`。
★★★これが `E[l]` が `([ζ], [π])` で**生成される**ことの中身であり、
第 1161-1164 の独立性・一意性と合わせて
**`([ζ], [π])` は `E[l]` の `ℤ/l`-基底である**が完全になる。

☆仮説 `hμ` は「`l` 乗して 1 になる元は `ζ` の冪」——体の中では
`μ_l` が巡回群であることそのものである。 -/
theorem zeta_pi_span {G : Type*} [CommGroup G] {l : ℕ} {Q ζ π : G}
    (hπl : π ^ l = Q)
    (hμ : ∀ y : G, y ^ l = 1 → ∃ a : ℤ, y = ζ ^ a)
    (x : G) (m : ℤ) (hx : x ^ l = Q ^ m) :
    ∃ a : ℤ, x = ζ ^ a * π ^ m := by
  have hπz : π ^ (l : ℤ) = Q := by rw [zpow_natCast]; exact hπl
  have hpi : (π ^ m) ^ l = Q ^ m := by
    rw [← zpow_natCast (π ^ m) l, ← zpow_mul, mul_comm m (l : ℤ), zpow_mul, hπz]
  have hy : (x * (π ^ m)⁻¹) ^ l = 1 := by
    rw [mul_pow, inv_pow, hx, hpi]
    simp
  obtain ⟨a, ha⟩ := hμ _ hy
  refine ⟨a, ?_⟩
  rw [← ha, mul_assoc]
  simp

/-- ★★★★★★★★★★★★★★★★★★★★
**`([ζ], [π])` は `l`-捩れの基底**——★**無条件**（第 1165）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆生成（`zeta_pi_span`）と一意性（`zeta_pi_coord_eq_iff`）を 1 本にまとめた形。
★★これが `AlphaBridge` の節点 1 の**結論そのもの**である。 -/
theorem zeta_pi_basis {G : Type*} [CommGroup G] {l : ℕ} (hl : 0 < l) {Q ζ π : G}
    (hζl : ζ ^ l = 1) (hζprim : ∀ n : ℕ, 0 < n → n < l → ζ ^ n ≠ 1)
    (hπl : π ^ l = Q) (hQinf : ∀ j : ℤ, Q ^ j = 1 → j = 0)
    (hμ : ∀ y : G, y ^ l = 1 → ∃ a : ℤ, y = ζ ^ a)
    (x : G) (m : ℤ) (hx : x ^ l = Q ^ m) :
    (∃ a b : ℤ, x = ζ ^ a * π ^ b)
      ∧ (∀ a₁ b₁ a₂ b₂ : ℤ, ζ ^ a₁ * π ^ b₁ = ζ ^ a₂ * π ^ b₂ →
          ((l : ℤ) ∣ (a₁ - a₂) ∧ (l : ℤ) ∣ (b₁ - b₂))) := by
  refine ⟨?_, ?_⟩
  · obtain ⟨a, ha⟩ := zeta_pi_span hπl hμ x m hx
    exact ⟨a, m, ha⟩
  · intro a₁ b₁ a₂ b₂ heq
    refine (zeta_pi_coord_eq_iff hl hζl hζprim hπl hQinf a₁ b₁ a₂ b₂).mp ⟨0, ?_⟩
    rw [heq, zpow_zero, mul_one]

/-! ## ★★★★★★★★★★★★商群の `l`-捩れ -/

/-- ☆`G ⧸ ⟨Q⟩` で `[x]ˡ = 1` は `xˡ ∈ ⟨Q⟩` と同じ。 -/
theorem quotient_pow_eq_one_iff {G : Type*} [CommGroup G] (Q : G) (l : ℕ) (x : G) :
    (QuotientGroup.mk x : G ⧸ Subgroup.zpowers Q) ^ l = 1 ↔ ∃ m : ℤ, x ^ l = Q ^ m := by
  rw [← QuotientGroup.mk_pow, QuotientGroup.eq_one_iff, Subgroup.mem_zpowers_iff]
  constructor
  · rintro ⟨m, hm⟩
    exact ⟨m, hm.symm⟩
  · rintro ⟨m, hm⟩
    exact ⟨m, hm.symm⟩

/-- ★★★★★★★★★★★★★★★★
**`G ⧸ ⟨Q⟩` の `l`-捩れはすべて `[ζᵃπᵇ]`**——★**無条件**（第 1170）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`zeta_pi_span`（第 1165）を商群の言葉に直したもの。
★★★Tate 一意化 `Φ : G ⧸ ⟨Q⟩ ≃+ E.Point` は加法同型なので、
`E[l]` の元は `Φ [ζᵃπᵇ]` の形に書ける——これが `AlphaBridge` の節点 2 の
**群論の側の全体**である。 -/
theorem exists_zeta_pi_of_torsion {G : Type*} [CommGroup G] {l : ℕ} {Q ζ π : G}
    (hπl : π ^ l = Q) (hμ : ∀ y : G, y ^ l = 1 → ∃ a : ℤ, y = ζ ^ a)
    (x : G) (hx : (QuotientGroup.mk x : G ⧸ Subgroup.zpowers Q) ^ l = 1) :
    ∃ a b : ℤ, (QuotientGroup.mk x : G ⧸ Subgroup.zpowers Q)
      = QuotientGroup.mk (ζ ^ a * π ^ b) := by
  obtain ⟨m, hm⟩ := (quotient_pow_eq_one_iff Q l x).mp hx
  obtain ⟨a, ha⟩ := zeta_pi_span hπl hμ x m hm
  exact ⟨a, m, by rw [ha]⟩

/-! ## ★★★★★★★★★★★★★★★★加法同型で運ぶ（Tate 一意化の側） -/

/-- ★★★★★★★★★★★★★★★★★★★★
**加法同型の下でも `l`-捩れは `[ζᵃπᵇ]` の像**——★**無条件**（第 1171）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`Φ : Additive (G ⧸ ⟨Q⟩) ≃+ P` を Tate 一意化 `Kˣ/qℤ ≃ E(K)` と読む。
★`l • Φ(c) = 0` なら `cˡ = 1` なので、第 1170 より `c = [ζᵃπᵇ]`。

★★★これが `AlphaBridge` の節点 2 の**完成形**である
——`E[l]` の元はすべて `Φ [ζᵃπᵇ]` の形に書ける。 -/
theorem torsion_eq_phi_zeta_pi {G P : Type*} [CommGroup G] [AddCommGroup P]
    {Q ζ π : G} {l : ℕ} (hπl : π ^ l = Q)
    (hμ : ∀ y : G, y ^ l = 1 → ∃ a : ℤ, y = ζ ^ a)
    (Φ : Additive (G ⧸ Subgroup.zpowers Q) ≃+ P)
    (x : G) (hx : l • Φ (Additive.ofMul (QuotientGroup.mk x)) = 0) :
    ∃ a b : ℤ, Φ (Additive.ofMul (QuotientGroup.mk x))
      = Φ (Additive.ofMul (QuotientGroup.mk (ζ ^ a * π ^ b))) := by
  have h1 : Φ (l • Additive.ofMul (QuotientGroup.mk x : G ⧸ Subgroup.zpowers Q)) = 0 := by
    rw [map_nsmul]
    exact hx
  have h2 : (l • Additive.ofMul (QuotientGroup.mk x : G ⧸ Subgroup.zpowers Q)) = 0 :=
    Φ.injective (by rw [h1, map_zero])
  have h3 : (QuotientGroup.mk x : G ⧸ Subgroup.zpowers Q) ^ l = 1 := h2
  obtain ⟨a, b, hab⟩ := exists_zeta_pi_of_torsion hπl hμ x h3
  exact ⟨a, b, by rw [hab]⟩

/-! ## ★★★★★★★★★★★★Tate 母数は無限位数（`hQinf` の中身） -/

/-- ★★★★★★★★★★★★★★**Tate 母数は無限位数**——★**無条件**（第 1172）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

☆`TateSetup` は `0 < v(q)` を持っている。`v(qʲ) = j·v(q)` なので
`qʲ = 1` なら `j·v(q) = 0`、したがって `j = 0`。

★★これが第 1161-1171 のすべての補題が受けている `hQinf` の**中身**である
——`TateSetup` からそのまま出るので、仮説ではなく**定理**になる。 -/
theorem tateSetup_Q_zpow_eq_one {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K]
    [Algebra R K] (S : ABC3.Found.GaloisRep.TateSetup R I K) (j : ℤ) (h : S.Q ^ j = 1) :
    j = 0 := by
  have hv : ABC3.Found.GaloisRep.vAdd S.v (S.Q ^ j)
      = j * ABC3.Found.GaloisRep.vAdd S.v S.Q :=
    ABC3.Found.GaloisRep.vAdd_zpow S.v S.Q j
  rw [h] at hv
  have h1 : ABC3.Found.GaloisRep.vAdd S.v (1 : Kˣ) = 0 := by
    simp [ABC3.Found.GaloisRep.vAdd]
  rw [h1] at hv
  rcases mul_eq_zero.mp hv.symm with hj | hq
  · exact hj
  · exact absurd hq S.hQ.ne'

/-! ## ★★★★★★★★★★★★`μ_l = ⟨ζ⟩`（`hμ` の中身） -/

/-- ★★★★★★★★★★★★★★**原始 `l` 乗根なら `μ_l = ⟨ζ⟩`**——★**無条件**（第 1173）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆整域の中では `l` 乗して 1 になる元はちょうど `l` 個で、原始根の冪で尽きる。
★★これが第 1165・1170・1171 が受けている `hμ` の**中身**である
——原始 `l` 乗根があれば仮説ではなく**定理**になる。 -/
theorem exists_zpow_of_pow_eq_one {K : Type*} [CommRing K] [IsDomain K] {l : ℕ} [NeZero l]
    {ζ : Kˣ} (hζ : IsPrimitiveRoot (ζ : K) l) (y : Kˣ) (hy : y ^ l = 1) :
    ∃ a : ℤ, y = ζ ^ a := by
  have hyK : (y : K) ^ l = 1 := by
    have := congrArg (fun u : Kˣ => (u : K)) hy
    simpa using this
  obtain ⟨i, _, hi⟩ := hζ.eq_pow_of_pow_eq_one hyK
  refine ⟨(i : ℤ), ?_⟩
  apply Units.ext
  rw [zpow_natCast]
  simpa using hi.symm

/-! ## ★★★★★★★★★★★★★★★★★★★★Tate 設定での組み立て -/

/-- ☆原始 `l` 乗根は単数群でも `l` 乗して 1。 -/
theorem units_pow_eq_one_of_isPrimitiveRoot {K : Type*} [CommRing K] {l : ℕ} {ζ : Kˣ}
    (hζ : IsPrimitiveRoot (ζ : K) l) : ζ ^ l = 1 := by
  apply Units.ext
  simpa using hζ.pow_eq_one

/-- ☆原始 `l` 乗根は `0 < n < l` で `ζⁿ ≠ 1`（単数群の側）。 -/
theorem units_pow_ne_one_of_isPrimitiveRoot {K : Type*} [CommRing K] {l : ℕ} {ζ : Kˣ}
    (hζ : IsPrimitiveRoot (ζ : K) l) (n : ℕ) (hn : 0 < n) (hnl : n < l) : ζ ^ n ≠ 1 := by
  intro hcon
  refine hζ.pow_ne_one_of_pos_of_lt hn.ne' hnl ?_
  have := congrArg (fun u : Kˣ => (u : K)) hcon
  simpa using this

/-- ★★★★★★★★★★★★★★★★★★★★
**Tate 設定での `α`**——★**無条件**（第 1174）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`TateSetup S` と原始 `l` 乗根 `ζ`、`πˡ = q` なる `π`、そして
`σ(ζ) = ζ`・`σ(π) = ζπ` なる `σ` があれば、`σ` は座標を `α = (1 1 / 0 1)` で動かす。

★★★**受けている仮説はこれだけ**である——`hQinf` は `TateSetup` から（第 1172）、
`hζ ^ l = 1` と原始性は `IsPrimitiveRoot` から（本ファイル）出るので、
`zeta_pi_*` の系列が要求していた仮説は**すべて消えた**。 -/
theorem tate_sigma_coord_alpha {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K]
    [Algebra R K] (S : ABC3.Found.GaloisRep.TateSetup R I K) {l : ℕ} (hl : 0 < l)
    {ζ π : Kˣ} (hζ : IsPrimitiveRoot (ζ : K) l) (hπl : π ^ l = S.Q)
    (σ : Kˣ →* Kˣ) (hσζ : σ ζ = ζ) (hσπ : σ π = ζ * π)
    (a b a' b' : ℤ)
    (h : ∃ n : ℤ, σ (ζ ^ a * π ^ b) = ζ ^ a' * π ^ b' * S.Q ^ n) :
    ((l : ℤ) ∣ (a + b - a')) ∧ ((l : ℤ) ∣ (b - b')) :=
  sigma_coord_alpha hl (units_pow_eq_one_of_isPrimitiveRoot hζ)
    (units_pow_ne_one_of_isPrimitiveRoot hζ) hπl
    (tateSetup_Q_zpow_eq_one S) σ hσζ hσπ a b a' b' h

/-- ★★★★★★★★★★★★★★★★★★★★
**Tate 設定での `l`-捩れの記述**——★**無条件**（第 1174）。

☆`Φ` を Tate 一意化とすると `E[l]` の元はすべて `Φ [ζᵃπᵇ]` の形である。
★★受けている仮説は「`ζ` が原始 `l` 乗根」「`πˡ = q`」だけである。 -/
theorem tate_torsion_eq_phi_zeta_pi {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K]
    [Algebra R K] (S : ABC3.Found.GaloisRep.TateSetup R I K) {l : ℕ} [NeZero l]
    {ζ π : Kˣ} (hζ : IsPrimitiveRoot (ζ : K) l) (hπl : π ^ l = S.Q)
    {P : Type*} [AddCommGroup P]
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (x : Kˣ) (hx : l • Φ (Additive.ofMul (QuotientGroup.mk x)) = 0) :
    ∃ a b : ℤ, Φ (Additive.ofMul (QuotientGroup.mk x))
      = Φ (Additive.ofMul (QuotientGroup.mk (ζ ^ a * π ^ b))) :=
  torsion_eq_phi_zeta_pi hπl (fun y hy => exists_zpow_of_pow_eq_one hζ y hy) Φ x hx

/-! ## ★出典の紐付け(`.src`) -/

def orderOf_eq_of_primitive.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(原始 l 乗根の位数はちょうど l。★無条件)",
    sectionId := "genell-thm-3-8" }

def zeta_pi_indep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ζ と π は Lˣ/⟨Q⟩ の中で独立。★無条件)",
    sectionId := "genell-thm-3-8" }

def tate_sigma_coord_alpha.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 設定での α——受けているのは ζ が原始 l 乗根で πˡ = q なことだけ。★無条件)",
    sectionId := "genell-thm-3-8" }

def tate_torsion_eq_phi_zeta_pi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 設定での l-捻れの記述——E[l] はすべて Φ [ζᵃπᵇ]。★無条件)",
    sectionId := "genell-thm-3-8" }

def tate_sigma_coord_alpha.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1174）**——`zeta_pi_*` の系列が要求していた仮説" ++
       "（`hζl`・`hζprim`・`hQinf`・`hμ`）は**すべて消えた**——" ++
       "`TateSetup` と `IsPrimitiveRoot` から出るからである。" ++
       "☆残るのはこの `ζ`・`π`・`σ` を `L_v(ζ_l, q^{1/l})` で実際に取る段と、" ++
       "`galRep` の行列に読み替える段である。") 1 ]

def exists_zpow_of_pow_eq_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(原始 l 乗根なら μ_l = ⟨ζ⟩——hμ の中身。★無条件)",
    sectionId := "genell-thm-3-8" }

def tateSetup_Q_zpow_eq_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 母数は無限位数——hQinf の中身。★無条件)",
    sectionId := "genell-def-3-3" }

def torsion_eq_phi_zeta_pi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(加法同型の下でも l-捻れは [ζᵃπᵇ] の像。★無条件)",
    sectionId := "genell-thm-3-8" }

def torsion_eq_phi_zeta_pi.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1171）**——`AlphaBridge` の節点 2 の**完成形**である。" ++
       "☆`Φ` を Tate 一意化 `Kˣ/qℤ ≃ E(K)` と読めば、" ++
       "`E[l]` の元はすべて `Φ [ζᵃπᵇ]` の形に書ける。" ++
       "★これで `AlphaBridge` の 3 節点はすべて揃った——" ++
       "残るのは具体の `TateSetup`・`galRep` に当てはめる配管だけである。") 1 ]

def quotient_pow_eq_one_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G ⧸ ⟨Q⟩ で [x]ˡ = 1 は xˡ ∈ ⟨Q⟩ と同じ。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_zeta_pi_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G ⧸ ⟨Q⟩ の l-捻れはすべて [ζᵃπᵇ]。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_zeta_pi_of_torsion.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1170）**——`AlphaBridge` の節点 2 の" ++
       "**群論の側の全体**である。☆Tate 一意化 `Φ : G ⧸ ⟨Q⟩ ≃+ E.Point` は" ++
       "加法同型なので、`E[l]` の元は `Φ [ζᵃπᵇ]` の形に書ける。" ++
       "★残るのは `Φ` を介して `E[l]` の言葉に直す配管だけである。") 2 ]

def zeta_pi_span.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(xˡ ∈ ⟨Q⟩ なら x = ζᵃπᵇ——全射の側。★無条件)",
    sectionId := "genell-thm-3-8" }

def zeta_pi_basis.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(([ζ], [π]) は l-捻れの基底——生成と一意性。★無条件)",
    sectionId := "genell-thm-3-8" }

def zeta_pi_basis.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1165）**——`AlphaBridge` の**節点 1 の結論そのもの**である。" ++
       "☆仮説 `hμ`（`l` 乗して 1 になる元は `ζ` の冪）は体の中で `μ_l` が巡回群である" ++
       "ことそのもの。★残るのは Tate 一意化で `E[l]` を `Lˣ/qℤ` の `l`-捻れと同定する段だけである。") 3 ]

def sigma_zeta_pi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ(ζ)=ζ・σ(π)=ζπ なら σ(ζᵃπᵇ) = ζ^{a+b}πᵇ——群の水準。★無条件)",
    sectionId := "genell-thm-3-8" }

def sigma_coord_alpha.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ の行列はちょうど α = (1 1 / 0 1)——座標の側の全体。★無条件)",
    sectionId := "genell-thm-3-8" }

def sigma_coord_alpha.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1164）**——これで `AlphaBridge` の**座標の側が全部揃った**。" ++
       "☆残るのは Tate 一意化で `(ζ, π)` を `E[l]` の基底として実現する段" ++
       "（節点 2 の残り）と、分解群経由で大域の像へ移す段（節点 3）である。") 3 ]

def zeta_pi_coord_eq_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ζᵃπᵇ の類が一致するのは座標が mod l で一致するとき。★無条件)",
    sectionId := "genell-thm-3-8" }

def zeta_pi_mem_zpowers_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(座標の核はちょうど lℤ × lℤ。★無条件)",
    sectionId := "genell-thm-3-8" }

def zeta_pi_indep.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-01（第 1161）**——`Skeleton/GenEll/AlphaBridge.lean` の" ++
       "節点 1 の**核**である。☆`E[l]` の基底を `([ζ], [π])` に取れることが、" ++
       "`sigma_acts_as_alpha`（第 993）の指数の動きを行列 `α = (1 1 / 0 1)` として" ++
       "読むための土台である。" ++
       "★残るのは Tate 一意化の下で `E[l]` が実際にこの 2 元で生成されること" ++
       "（全射の側）と、`Q` が無限位数であること（`0 < v(q)` から）の紐付けである。") 5 ]

end ABC3.Found.GenEll
