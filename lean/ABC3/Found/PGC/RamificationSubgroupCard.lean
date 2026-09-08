import ABC3.Found.PGC.RamificationGroupNormBridge

/-!
# [pGC] (b) `|G_i|` を跳びの列 `u` で書き下す —— ★**純群論で閉じた**

前波(第 1118)で `hrec` の橋は (a)(b)(c) に割れ、(a) は閉じ (c) は在庫だった。
残った (b)「`|G_i| = p^{k+1−m_i}`」を、★**環も型も建てずに**閉じる。

## ★★見積り(着手前)と結果

★「`𝒪_M` を型として建てる」道(#69 の危険区間、`PAdicLocalField` の同型不変性 D13/D30–D32)を
**取らなかった**。理由は着手前の見積り:

> (b) の中身は「`G_i` が `⟨g^{p^{m}}⟩` のどれかであること」で、
> それは `Gr` が**部分群である**ことと `⟨g^j⟩ = ⟨g^{p^{v_p(j)}}⟩` だけから出る。
> ★環も付値もノルムも要らない。

★**当たった。**本ファイルは `Subgroup G` と `orderOf` だけで閉じている
(★分岐・付値・ノルム・環の語彙が 1 語も出ない)。

## 出したもの

| 宣言 | 内容 |
|---|---|
| `orderOf_pow_prime_pow` | `orderOf g = p^{k+1}` ⇒ `orderOf (g^{p^m}) = p^{k+1−m}` |
| `zpowers_pow_eq_of_not_dvd` | ★`p ∤ j` なら `⟨h^j⟩ = ⟨h⟩`(`h` の位数が `p` 冪)。modular inverse |
| `zpowers_pow_eq_zpowers_pow_padicValNat` | ★`⟨g^j⟩ = ⟨g^{p^{v_p(j)}}⟩` |
| ★★`eq_zpowers_of_mem_iff` | `g^{p^s} ∈ Gr ⟺ i ≤ u s` ⇒ **`Gr = ⟨g^{p^m}⟩`**(`m` は最小) |
| ★★★`card_eq_pow_of_mem_iff` | **`Nat.card Gr = p^{k+1−m}`** |

★入力の `hmem : ∀ s ≤ k+1, (g^{p^s} ∈ Gr ↔ i ≤ u s)` は、
前波 `RamNormBridge.mem_ramification_iff`
(`(∀ z, ‖z‖ ≤ 1 → ‖h z − z‖ ≤ ‖π‖^{i+1}) ↔ i ≤ t`)が★**ノルムの言葉で供給する**
(`h := g^{p^s}`、`t := u s`)。

## ★★残り(正確に、`file:line` つき)

`hrec` の橋に残っているのは★**1 つだけ**:

★**`herbrandPhiGroup` の `Nat.card (lowerRamificationGroup B G i)` と、
本ファイルの `Nat.card Gr` を同じ `Gr` にすること。**
すなわち `lowerRamificationGroup B G i`(`LowerRamificationGroup.lean:265`、
`B` は DVR、`𝔪_B^{i+1}` の inertia)を、ノルムの言葉の
`{σ | ∀ z, ‖z‖ ≤ 1 → ‖σ z − z‖ ≤ ‖π‖^{i+1}}` と**同一視**すること。

★これには `B = 𝒪_M` を型として建てるしかない(`herbrandPhiGroup` の
`[IsDiscreteValuationRing B] [MulSemiringAction G B]` が `B` を要求するため)。
★**持ち場が挙げた 3 つはどれも「無い」ではなかった**(本体の追測どおり)が、
本ファイルはそこに**降りていない**。降りるときの既知の危険:

* `lean-idioms.md #69`(`adjoinField`/`adjoinIntegers` の境界、212 秒 timeout)
* `PAdicLocalField` の同型不変性(D13/D30–D32)
* `#323`(`NormedField` のまま `Algebra.adjoin` の所属を書くと `isDefEq` が焼き切れる)

⇒ ★★**残りは「数学」ではなく「`𝒪_M` の型を建てて 4 つのインスタンスを載せ、
`lowerRamificationGroup` を `mem_ramification_iff` と突き合わせる」1 ノードである。**

## 逸脱の記録

1. `Gr` は**任意の部分群**として受ける(`lowerRamificationGroup` を参照しない)。
   ★これが「環を建てずに済む」ようにした逸脱である。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**(import のみ)。
-/

namespace ABC3.Found.PGC

namespace RamCard

/-! ## §1 巡回 `p` 群の部分群 -/

section Cyclic

variable {G : Type*} [Group G]

/-- `orderOf g = p^{k+1}` なら `orderOf (g^{p^m}) = p^{k+1−m}`(`m ≤ k+1`)。 -/
theorem orderOf_pow_prime_pow {p k m : ℕ} (hp : p.Prime) {g : G}
    (hg : orderOf g = p ^ (k + 1)) (hm : m ≤ k + 1) :
    orderOf (g ^ p ^ m) = p ^ (k + 1 - m) := by
  rw [orderOf_pow' _ (pow_ne_zero m hp.pos.ne'), hg,
    Nat.gcd_eq_right (pow_dvd_pow p hm), Nat.pow_div hm hp.pos]

/-- ★`p ∤ j` なら `⟨h^j⟩ = ⟨h⟩`(`h` の位数が `p` 冪)。★modular inverse の道。 -/
theorem zpowers_pow_eq_of_not_dvd {p r : ℕ} (hp : p.Prime) {h : G}
    (hh : orderOf h = p ^ r) {j : ℕ} (hj : ¬ (p ∣ j)) :
    Subgroup.zpowers (h ^ j) = Subgroup.zpowers h := by
  haveI : Fact p.Prime := ⟨hp⟩
  refine le_antisymm (Subgroup.zpowers_le.mpr (pow_mem (Subgroup.mem_zpowers h) j)) ?_
  refine Subgroup.zpowers_le.mpr ?_
  have hcop : Nat.Coprime j (p ^ r) :=
    (Nat.Coprime.pow_right r ((hp.coprime_iff_not_dvd.mpr hj).symm))
  set c : ℕ := ((j : ZMod (p ^ r))⁻¹).val with hcdef
  have hmod : j * c ≡ 1 [MOD p ^ r] := by
    rw [← ZMod.natCast_eq_natCast_iff]
    push_cast
    rw [hcdef, ZMod.natCast_val, ZMod.cast_id]
    exact ZMod.coe_mul_inv_eq_one j hcop
  have hh1 : h ^ (j * c) = h ^ 1 := (pow_eq_pow_iff_modEq).mpr (by rw [hh]; exact hmod)
  have hgen : (h ^ j) ^ c = h := by rw [← pow_mul, hh1, pow_one]
  exact ⟨(c : ℤ), by simpa using hgen⟩

/-- ★`⟨g^j⟩ = ⟨g^{p^{v_p(j)}}⟩`(`j ≠ 0`、`v_p(j) ≤ k`)。 -/
theorem zpowers_pow_eq_zpowers_pow_padicValNat {p k : ℕ} (hp : p.Prime) {g : G}
    (hg : orderOf g = p ^ (k + 1)) {j : ℕ} (hj0 : j ≠ 0)
    (hjk : padicValNat p j ≤ k + 1) :
    Subgroup.zpowers (g ^ j) = Subgroup.zpowers (g ^ p ^ (padicValNat p j)) := by
  haveI : Fact p.Prime := ⟨hp⟩
  set s : ℕ := padicValNat p j with hs
  set j' : ℕ := j / p ^ s with hj'
  have hjeq : p ^ s * j' = j := Nat.mul_div_cancel' pow_padicValNat_dvd
  have hnd : ¬ (p ∣ j') := by
    intro ⟨d, hd⟩
    refine pow_succ_padicValNat_not_dvd (p := p) hj0 ⟨d, ?_⟩
    calc j = p ^ s * j' := hjeq.symm
      _ = p ^ s * (p * d) := by rw [hd]
      _ = p ^ (s + 1) * d := by rw [pow_succ]; ring
  have hord : orderOf (g ^ p ^ s) = p ^ (k + 1 - s) := orderOf_pow_prime_pow hp hg hjk
  have hrw : g ^ j = (g ^ p ^ s) ^ j' := by rw [← pow_mul, hjeq]
  rw [hrw]
  exact zpowers_pow_eq_of_not_dvd hp hord hnd

end Cyclic

/-! ## §2 `|G_i|` を跳びの列で書き下す -/

section Card

variable {G : Type*} [Group G] [Finite G]

/-- ★★★★★**`G_i` は跳びの列が決める** ——
`G = ⟨g⟩`(位数 `p^{k+1}`)で `g^{p^s} ∈ Gr ⟺ i ≤ u s` なら `Gr = ⟨g^{p^m}⟩`
(`m` は `i ≤ u m` なる最小)。★分岐・付値・ノルムの語彙が 1 語も出ない。 -/
theorem eq_zpowers_of_mem_iff {p k : ℕ} (hp : p.Prime) {g : G}
    (hg : orderOf g = p ^ (k + 1)) (htop : ∀ σ : G, σ ∈ Subgroup.zpowers g)
    {Gr : Subgroup G} {u : ℕ → ℤ} {i : ℤ} {m : ℕ} (hmk : m ≤ k + 1)
    (hmem : ∀ s : ℕ, s ≤ k + 1 → (g ^ p ^ s ∈ Gr ↔ i ≤ u s))
    (hm : i ≤ u m) (hmin : ∀ s : ℕ, s ≤ k + 1 → i ≤ u s → m ≤ s) :
    Gr = Subgroup.zpowers (g ^ p ^ m) := by
  haveI : Fact p.Prime := ⟨hp⟩
  refine le_antisymm ?_ (Subgroup.zpowers_le.mpr ((hmem m hmk).mpr hm))
  intro σ hσ
  obtain ⟨j, hj⟩ := mem_powers_iff_mem_zpowers.mpr (htop σ)
  have hjσ : g ^ j = σ := hj
  by_cases hdvd : p ^ (k + 1) ∣ j
  · have hone : σ = 1 := by
      rw [← hjσ]
      exact orderOf_dvd_iff_pow_eq_one.mp (by rw [hg]; exact hdvd)
    rw [hone]
    exact Subgroup.one_mem _
  · have hj0 : j ≠ 0 := by
      intro h
      exact hdvd (by rw [h]; exact dvd_zero _)
    have hsk : padicValNat p j ≤ k := by
      by_contra hcon
      rw [not_le] at hcon
      exact hdvd (dvd_trans (pow_dvd_pow p (by omega : k + 1 ≤ padicValNat p j))
        pow_padicValNat_dvd)
    have hzp : Subgroup.zpowers (g ^ j) = Subgroup.zpowers (g ^ p ^ (padicValNat p j)) :=
      zpowers_pow_eq_zpowers_pow_padicValNat hp hg hj0 (by omega)
    have hmemGr : g ^ p ^ (padicValNat p j) ∈ Gr := by
      have hle : Subgroup.zpowers (g ^ j) ≤ Gr :=
        Subgroup.zpowers_le.mpr (by rw [hjσ]; exact hσ)
      rw [hzp] at hle
      exact hle (Subgroup.mem_zpowers _)
    have hms : m ≤ padicValNat p j :=
      hmin _ (by omega) ((hmem _ (by omega)).mp hmemGr)
    have hdvd2 : p ^ m ∣ j :=
      dvd_trans (pow_dvd_pow p hms) pow_padicValNat_dvd
    obtain ⟨c, hc⟩ := hdvd2
    refine ⟨(c : ℤ), ?_⟩
    rw [← hjσ, hc, pow_mul]
    simp

/-- ★★★★★★**`|G_i| = p^{k+1−m}`**。 -/
theorem card_eq_pow_of_mem_iff {p k : ℕ} (hp : p.Prime) {g : G}
    (hg : orderOf g = p ^ (k + 1)) (htop : ∀ σ : G, σ ∈ Subgroup.zpowers g)
    {Gr : Subgroup G} {u : ℕ → ℤ} {i : ℤ} {m : ℕ} (hmk : m ≤ k + 1)
    (hmem : ∀ s : ℕ, s ≤ k + 1 → (g ^ p ^ s ∈ Gr ↔ i ≤ u s))
    (hm : i ≤ u m) (hmin : ∀ s : ℕ, s ≤ k + 1 → i ≤ u s → m ≤ s) :
    Nat.card Gr = p ^ (k + 1 - m) := by
  rw [eq_zpowers_of_mem_iff hp hg htop hmk hmem hm hmin, Nat.card_zpowers,
    orderOf_pow_prime_pow hp hg hmk]

end Card

/-! ## `.src` と 公理 -/

def eq_zpowers_of_mem_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def card_eq_pow_of_mem_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms orderOf_pow_prime_pow
#print axioms zpowers_pow_eq_of_not_dvd
#print axioms zpowers_pow_eq_zpowers_pow_padicValNat
#print axioms eq_zpowers_of_mem_iff
#print axioms card_eq_pow_of_mem_iff

end RamCard

end ABC3.Found.PGC
