import ABC3.Found.PGC.HerbrandIntegralNorm

/-!
# [pGC] `hrec` —— Herbrand 関数の 1 段の漸化式

## ★★★まず測ったこと（1）添字の対応

持ち場は「`herbrandPhiGroup_natCast` の `φ_G(n)` と
`dvd_sub_of_phi_intCast` の `φ m` の添字が同じものか未測定」と言っていた。
★**同じものではない。** 正しい対応は

    φ m := herbrandPhiGroup G π ((u m : ℕ) : ℝ)        ★すなわち n = u m

である（`n` は分岐の添字、`m` は跡びの段）。この対応の下でのみ
`dvd_sub_of_phi_intCast` が要求する分母 `p^{m+1}` が出る:

    φ_G(u(m+1)) − φ_G(u m) = (u(m+1) − u m)·p^{k−m}/p^{k+1} = (u(m+1) − u m)/p^{m+1}.

★分母の `p^{m+1}` は「`|G_i| = p^{k−m}` と `|G| = p^{k+1}` の比」から出ている。
`k` が消えるのが鍵で、これが合わなければ Hasse–Arf の合同は出ない。

## ★★★まず測ったこと（2）★**配られた形の一部は偽だった**（訂正）

`RamCard.card_eq_pow_of_mem_iff`(`RamificationSubgroupCard.lean:174-178`)の

```
    (hmem : ∀ s : ℕ, s ≤ k + 1 → (g ^ p ^ s ∈ Gr ↔ i ≤ u s))
```

をノルムで供給しようとすると、★**`s = k + 1` で破綻する**。
`orderOf g = p^{k+1}` なので `g^{p^{k+1}} = 1` であり、

    ‖(g^{p^{k+1}})π − π‖ = ‖0‖ = 0 ≠ ‖π‖^{u(k+1)+1}   （右辺は常に正）

⇒ ★★**跡びの列 `u` が定義されるのは `s ≤ k` までであり、
`s = k+1` を含む `hjump` は満たせない。**
（`RamificationSubgroupCard.lean` 自体は真である——仮説が強いだけ。
他ファイルの docstring は書き換えないのでここに訂正を書く。）

本ファイルの `card_eq_pow_of_mem_Ioc` は `hmem` を ★**`∀ s ≤ k`** に弱め、
`s = k+1` の場合を「`g^{p^{k+1}} = 1 ∈ Gr`（常に真）」と
`u` の `k` での建て直し（`u' s := if s ≤ k then u s else u k`）で**内部で**処理する。
その分 `m + 1 ≤ k` を要求するが、`harith` (3) は `∀ m < k` なのでちょうど合う。

## 出したもの

| 宣言 | 内容 |
|---|---|
| `sum_Icc_eq_add_sum_Ioc` | ★抽象核（純 `Finset`）`Σ_{Icc 1 b} = Σ_{Icc 1 a} + Σ_{Ioc a b}` |
| `sum_Ioc_const` | ★抽象核 区間で定数なら和は `(b−a)·C` |
| ★★`herbrand_step_abstract` | ★抽象核 `φ(b) = φ(a) + (b−a)·C/D` |
| ★★`card_eq_pow_of_mem_Ioc` | ★抽象核（純群論）`u m < i ≤ u(m+1) ⇒ |G_i| = p^{k−m}` |
| `mem_lowerRamificationGroup_iff_jump` | `g^{p^s} ∈ G_i ↔ i ≤ u s`（ノルムで供給）|
| `card_eq_orderOf_of_forall_mem_zpowers` | `Nat.card G = orderOf g` |
| `step_arith` | `p^{k−m}/p^{k+1} = 1/p^{m+1}` |
| ★★★★★★`herbrandPhiGroup_step` | **`hrec` 本体** |

## ★在庫の測定（コマンドを残す）

```
grep -n "Ioc_consecutive" .cache/mathlib-index.txt
  → Finset.prod_Ioc_consecutive  Algebra/BigOperators/Intervals.lean:61 だけが出る。
  ★だが加法版 `Finset.sum_Ioc_consecutive` は **実際に存在する**（`to_additive` 生成名、
  索引の「無い」は嘘の 1 つ目の型）。本ファイルはそれを使って通した。
grep -n "Subgroup.card_top\|Nat.card_zpowers" .cache/mathlib-index.txt
  → Algebra/Group/Subgroup/Finite.lean:111 / Data/ZMod/QuotientGroup.lean:161（両方在る）
grep -n "theorem card_eq_pow_of_mem_iff" -A 8 lean/ABC3/Found/PGC/RamificationSubgroupCard.lean
grep -n "theorem herbrandPhiGroup_natCast" -A 8 lean/ABC3/Found/PGC/HerbrandComposition.lean
```

★「無いと思ったが在った」の道具: `Subgroup.eq_top_iff'` は在るが
★**`Subgroup.eq_top_iff'.mpr` と書くと定数名としてパースされる**:
`Unknown constant `Subgroup.eq_top_iff'.mpr``。`(Subgroup.eq_top_iff' _).mpr` と括弧が要る。

## ★★残り（正確に、`file:line` つき）

`harith` (3) `p^{m+1} ∣ u(m+1) − u m` を閉じるには
`HasseArfCongruence.lean:112 dvd_sub_of_phi_intCast` の 2 入力が要る。

* `hrec` —— ★**本ファイルで閉じた**（`herbrandPhiGroup_step`）。
* `hint` —— `HerbrandIntegralNorm.lean:exists_herbrandPhiGroup_natCast_of_norm` が供給するが、
  その仕様に `hne : G_{u m} ≠ G_{u m + 1}` がある。これは
  `|G_{u m}| = p^{k+1−m}`、`|G_{u m + 1}| = p^{k−m}` から出るはずだが、
  ★**前者の形（`i = u m` ちょうどのときの `|G_i|`）は本ファイルにまだ無い**。
  必要なのは `card_eq_pow_of_mem_Ioc` の「下端を開けない版」
  （`∀ s < m, u s < i` と `i ≤ u m` から `p^{k+1−m}`）である。
* `dvd_sub_of_phi_intCast` は `hint`/`hrec` を ★**`∀ m`** で要求するが、
  本ファイルの `hrec` は `m + 1 ≤ k` しか出ない。
  ★`u` を `fun s => u (min s k)` に切り詰めれば `m ≥ k` では差が 0 になるので
  両方とも `∀ m` に伸びる（★未測定。本波ではここに降りていない）。

## 逸脱の記録

1. `G := M ≃ₐ[K] M` に固定（`RamNormBridge.mem_ramification_iff` がこの形）。
2. `u : ℕ → ℕ`（`RamCard` 側は `ℕ → ℤ`）。ノルム側の跡びは `ℕ` で出るため。
3. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

namespace IntegerNorm

/-! ## §1 抽象核（純 `Finset`／ℝ） -/

section Core

/-- ★**抽象核** —— `Icc 1 b` の和は `Icc 1 a` の和と `Ioc a b` の和に割れる。 -/
theorem sum_Icc_eq_add_sum_Ioc {A : Type*} [AddCommMonoid A] (c : ℕ → A) {a b : ℕ} (hab : a ≤ b) :
    ∑ i ∈ Finset.Icc 1 b, c i = (∑ i ∈ Finset.Icc 1 a, c i) + ∑ i ∈ Finset.Ioc a b, c i := by
  have h1 : ∀ n : ℕ, Finset.Icc 1 n = Finset.Ioc 0 n := by
    intro n; ext x; simp [Nat.lt_iff_add_one_le]
  rw [h1, h1]
  exact (Finset.sum_Ioc_consecutive c (Nat.zero_le a) hab).symm

/-- ★**抽象核** —— `Ioc a b` の上で定数なら和は `(b − a)·C`。 -/
theorem sum_Ioc_const {c : ℕ → ℝ} {a b : ℕ} (hab : a ≤ b) {C : ℝ}
    (hc : ∀ i, a < i → i ≤ b → c i = C) :
    ∑ i ∈ Finset.Ioc a b, c i = ((b : ℝ) - (a : ℝ)) * C := by
  have hcong : ∀ i ∈ Finset.Ioc a b, c i = C := by
    intro i hi
    rw [Finset.mem_Ioc] at hi
    exact hc i hi.1 hi.2
  rw [Finset.sum_congr rfl hcong, Finset.sum_const, Nat.card_Ioc, nsmul_eq_mul]
  congr 1
  have : ((b - a : ℕ) : ℝ) = (b : ℝ) - (a : ℝ) := by
    rw [Nat.cast_sub hab]
  exact this

/-- ★★**抽象核** —— 1 段の Herbrand 差分。

`φ(n) = (Σ_{i∈Icc 1 n} c i)/D` に対し、`c` が `Ioc a b` の上で定数 `C` なら

    φ(b) = φ(a) + (b − a)·C/D. -/
theorem herbrand_step_abstract {c : ℕ → ℝ} {a b : ℕ} (hab : a ≤ b) {C D : ℝ}
    (hc : ∀ i, a < i → i ≤ b → c i = C) :
    (∑ i ∈ Finset.Icc 1 b, c i) / D
      = (∑ i ∈ Finset.Icc 1 a, c i) / D + ((b : ℝ) - (a : ℝ)) * C / D := by
  rw [sum_Icc_eq_add_sum_Ioc c hab, sum_Ioc_const hab hc, add_div]

end Core

/-! ## §2 抽象核（純群論）—— `u m < i ≤ u (m+1)` の区間で `|G_i|` は定数 -/

section Group

variable {G : Type*} [Group G] [Finite G]

/-- ★★**抽象核（純群論）** —— `u m < i ≤ u (m+1)` なら `Nat.card Gr = p^{k−m}`。

`RamCard.card_eq_pow_of_mem_iff`(`RamificationSubgroupCard.lean:174`)の包み。
★分岐・付値・ノルムの語が 1 語も出ない。

★★**本波の検算で見つけた穴（訂正）**: 元の `hmem` は `∀ s ≤ k + 1` を要求するが、
★**`s = k + 1` では `g^{p^{k+1}} = 1` なので `‖(g^{p^s})π − π‖ = ‖π‖^{u s + 1}` は満たせない**
（左辺が 0、右辺が正）。すなわち「跡びの列 `u` を `s ≤ k+1` で与える」のは偽であり、
跡びが定義されるのは `s ≤ k` までである。
⇒ 本補題は `hmem` を `∀ s ≤ k` に弱め、`s = k+1` の場合を
`g^{p^{k+1}} = 1 ∈ Gr`（常に真）と `u` の `k` での建て直しで**内部で**処理する。 -/
theorem card_eq_pow_of_mem_Ioc {p k : ℕ} (hp : p.Prime) {g : G}
    (hg : orderOf g = p ^ (k + 1)) (htop : ∀ σ : G, σ ∈ Subgroup.zpowers g)
    {Gr : Subgroup G} {u : ℕ → ℕ} (hmono : ∀ s t : ℕ, s ≤ t → u s ≤ u t)
    {i m : ℕ} (hmk : m + 1 ≤ k)
    (hmem : ∀ s : ℕ, s ≤ k → (g ^ p ^ s ∈ Gr ↔ i ≤ u s))
    (hlo : u m < i) (hhi : i ≤ u (m + 1)) :
    Nat.card Gr = p ^ (k - m) := by
  classical
  have hik : i ≤ u k := le_trans hhi (hmono (m + 1) k hmk)
  set u' : ℕ → ℤ := fun s => if s ≤ k then ((u s : ℕ) : ℤ) else ((u k : ℕ) : ℤ) with hu'
  have hu'le : ∀ s : ℕ, s ≤ k → u' s = ((u s : ℕ) : ℤ) := by
    intro s hs; simp [hu', hs]
  have hone : g ^ p ^ (k + 1) = 1 := by rw [← hg]; exact pow_orderOf_eq_one g
  have hkey := RamCard.card_eq_pow_of_mem_iff (G := G) hp hg htop
    (Gr := Gr) (u := u') (i := (i : ℤ)) (m := m + 1) (by omega)
    (fun s hs => by
      rcases Nat.lt_or_ge s (k + 1) with hsk | hsk
      · have hsk' : s ≤ k := by omega
        rw [hu'le s hsk', hmem s hsk']
        exact_mod_cast Iff.rfl
      · have hsk' : s = k + 1 := by omega
        subst hsk'
        have hnot : ¬ (k + 1 ≤ k) := by omega
        simp only [hu', hnot, if_false]
        constructor
        · intro _; exact_mod_cast hik
        · intro _; rw [hone]; exact Subgroup.one_mem Gr)
    (by rw [hu'le (m + 1) hmk]; exact_mod_cast hhi)
    (fun s _ hs => by
      by_contra hcon
      have hsm : s ≤ m := by omega
      have hsk : s ≤ k := by omega
      rw [hu'le s hsk] at hs
      have h1 : i ≤ u s := by exact_mod_cast hs
      have h2 : u s ≤ u m := hmono s m hsm
      omega)
  rw [hkey]
  congr 1
  omega

end Group

/-! ## §3 具体層 —— `hmem` をノルムで供給する -/

section Norm

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]
  [FiniteDimensional K M]

/-- ★**`hmem` の供給** —— `g^{p^s} ∈ G_i ↔ i ≤ u s`。

`IntegerRingInstances.mem_lowerRamificationGroup_iff_norm`（環→ノルム）と
`RamNormBridge.mem_ramification_iff`（ノルム→跡び）を繋ぐだけ。 -/
theorem mem_lowerRamificationGroup_iff_jump
    {π : M} {n : ℕ} {u : ℕ → ℕ} {p : ℕ} {g : M ≃ₐ[K] M}
    (hiso : ∀ (σ : M ≃ₐ[K] M) (z : M), ‖σ • z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ j : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * j))
    (hval : ∀ z : M, z ≠ 0 → ∃ j : ℤ, ‖z‖ = ‖π‖ ^ j) (hπmem : π ∈ integerSubring M)
    (i s : ℕ) (hjump : ‖(g ^ p ^ s) π - π‖ = ‖π‖ ^ (u s + 1)) :
    letI := isLocalRing_integerSubring (M := M)
    letI := integerMulSemiringAction hiso
    (g ^ p ^ s ∈ lowerRamificationGroup ↥(integerSubring M) (M ≃ₐ[K] M) i ↔ i ≤ u s) := by
  letI := isLocalRing_integerSubring (M := M)
  letI := integerMulSemiringAction hiso
  rw [mem_lowerRamificationGroup_iff_norm hiso hπ0 hπ1 hval hπmem]
  exact RamNormBridge.mem_ramification_iff (i := i) (g ^ p ^ s)
    (fun w => hiso (g ^ p ^ s) w) hπ0 hπ1 hn hvalK hjump

/-- ★巡回群の位数。 -/
theorem card_eq_orderOf_of_forall_mem_zpowers {G : Type*} [Group G] {g : G}
    (htop : ∀ σ : G, σ ∈ Subgroup.zpowers g) : Nat.card G = orderOf g := by
  have h : (Subgroup.zpowers g) = ⊤ := (Subgroup.eq_top_iff' _).mpr htop
  rw [← Nat.card_zpowers g, h, Subgroup.card_top]

end Norm

/-! ## §4 ★★★`hrec` 本体 -/

section Rec

/-- ★算術の段。`p^{k−m}/p^{k+1} = 1/p^{m+1}`。 -/
theorem step_arith {p k m : ℕ} (hp : 0 < p) (hmk : m ≤ k) (X : ℝ) :
    X * ((p ^ (k - m) : ℕ) : ℝ) / ((p ^ (k + 1) : ℕ) : ℝ)
      = X / (p : ℝ) ^ (m + 1) := by
  have hp0 : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have hcast1 : ((p ^ (k - m) : ℕ) : ℝ) = (p : ℝ) ^ (k - m) := by push_cast; ring
  have hcast2 : ((p ^ (k + 1) : ℕ) : ℝ) = (p : ℝ) ^ (k + 1) := by push_cast; ring
  rw [hcast1, hcast2, div_eq_div_iff (by positivity) (by positivity), mul_assoc,
    ← pow_add]
  congr 2
  omega

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]
  [FiniteDimensional K M]

/-- ★★★★★★**`hrec`** —— Herbrand 関数の 1 段の漸化式。

    φ_G(u(m+1)) = φ_G(u m) + (u(m+1) − u m)/p^{m+1}.

★★**添字の対応**（本波で最初に測ったこと）:
`herbrandPhiGroup_natCast` の `n` は**分岐の添字**、
`dvd_sub_of_phi_intCast` の `φ m` の `m` は**跡びの段**であり、
★**同じものではない**。正しい対応は

    φ m := herbrandPhiGroup G π (u m)      ★即ち n = u m

である。この対応の下でのみ分母が `p^{m+1}` になる（下の計算）。

計算: `u m < i ≤ u(m+1)` の区間で `|G_i| = p^{k−m}` なので

    φ_G(u(m+1)) − φ_G(u m) = (u(m+1) − u m)·p^{k−m}/p^{k+1}
                            = (u(m+1) − u m)/p^{m+1}. -/
theorem herbrandPhiGroup_step
    {π : M} {n p k m : ℕ} {u : ℕ → ℕ} {g : M ≃ₐ[K] M} (hp : p.Prime)
    (hiso : ∀ (σ : M ≃ₐ[K] M) (z : M), ‖σ • z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ j : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * j))
    (hval : ∀ z : M, z ≠ 0 → ∃ j : ℤ, ‖z‖ = ‖π‖ ^ j) (hπmem : π ∈ integerSubring M)
    (hjump : ∀ s : ℕ, s ≤ k → ‖(g ^ p ^ s) π - π‖ = ‖π‖ ^ (u s + 1))
    (hmono : ∀ s t : ℕ, s ≤ t → u s ≤ u t)
    (hg : orderOf g = p ^ (k + 1)) (htop : ∀ σ : M ≃ₐ[K] M, σ ∈ Subgroup.zpowers g)
    (hmk : m + 1 ≤ k) :
    letI := integerMulSemiringAction hiso
    letI := isDiscreteValuationRing_integerSubring hπ0 hπ1 hval hπmem
    herbrandPhiGroup (M ≃ₐ[K] M) (⟨π, hπmem⟩ : ↥(integerSubring M)) ((u (m + 1) : ℕ) : ℝ)
      = herbrandPhiGroup (M ≃ₐ[K] M) (⟨π, hπmem⟩ : ↥(integerSubring M)) ((u m : ℕ) : ℝ)
        + (((u (m + 1) : ℕ) : ℝ) - ((u m : ℕ) : ℝ)) / (p : ℝ) ^ (m + 1) := by
  letI := isLocalRing_integerSubring (M := M)
  letI := integerMulSemiringAction hiso
  letI := isDiscreteValuationRing_integerSubring hπ0 hπ1 hval hπmem
  letI := baseAlgebra K M
  letI := smulCommClass_base (K := K) hiso (fun σ a => σ.commutes a)
  have huni := maximalIdeal_eq_span hπ0 hπ1 hval hπmem
  have hadj := adjoin_pi_eq_top (K := K) hπ0 hπ1 hn hvalK hπmem
  have hcard : Nat.card (M ≃ₐ[K] M) = p ^ (k + 1) := by
    rw [card_eq_orderOf_of_forall_mem_zpowers htop, hg]
  have hab : u m ≤ u (m + 1) := hmono m (m + 1) (by omega)
  have hc : ∀ i : ℕ, u m < i → i ≤ u (m + 1) →
      (Nat.card ↥(lowerRamificationGroup ↥(integerSubring M) (M ≃ₐ[K] M) i) : ℝ)
        = ((p ^ (k - m) : ℕ) : ℝ) := by
    intro i h1 h2
    have := card_eq_pow_of_mem_Ioc (G := M ≃ₐ[K] M) hp hg htop hmono hmk
      (fun s hs => mem_lowerRamificationGroup_iff_jump hiso hπ0 hπ1 hn hvalK hval hπmem i s
        (hjump s hs))
      h1 h2
    exact_mod_cast congrArg (fun x : ℕ => (x : ℝ)) this
  rw [herbrandPhiGroup_natCast (A := ↥(baseIntegerSubring K M)) huni hadj (u (m + 1)),
    herbrandPhiGroup_natCast (A := ↥(baseIntegerSubring K M)) huni hadj (u m),
    herbrand_step_abstract (c := fun i =>
      (Nat.card ↥(lowerRamificationGroup ↥(integerSubring M) (M ≃ₐ[K] M) i) : ℝ)) hab hc,
    hcard]
  congr 1
  exact step_arith hp.pos (by omega) _

end Rec




/-! ## `.src` と 公理 -/

def herbrandPhiGroup_step.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Lemma 6.10", sectionId := "lemma-6-10" }

def card_eq_pow_of_mem_Ioc.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms sum_Icc_eq_add_sum_Ioc
#print axioms sum_Ioc_const
#print axioms herbrand_step_abstract
#print axioms card_eq_pow_of_mem_Ioc
#print axioms mem_lowerRamificationGroup_iff_jump
#print axioms card_eq_orderOf_of_forall_mem_zpowers
#print axioms step_arith
#print axioms herbrandPhiGroup_step

end IntegerNorm

end ABC3.Found.PGC
