/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.GroupTheory.OrderOfElement
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

/-! ## ★出典の紐付け(`.src`) -/

def orderOf_eq_of_primitive.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(原始 l 乗根の位数はちょうど l。★無条件)",
    sectionId := "genell-thm-3-8" }

def zeta_pi_indep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ζ と π は Lˣ/⟨Q⟩ の中で独立。★無条件)",
    sectionId := "genell-thm-3-8" }

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
