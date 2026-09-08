import ABC3.Found.PGC.HerbrandRecurrence

/-!
# [pGC] ★★★★★★★`harith` (3) が閉じた —— `p^{m+1} ∣ u(m+1) − u m`

前波で `hrec`（`HerbrandRecurrence.lean`）が閉じ、残りは `hint` の入力
`hne : G_{u m} ≠ G_{u m + 1}` だけだった。本ファイルはそれを埋め、
`harith` (3) を ★**ノルムだけの仮説から**出す。

## ★★本波で覚したもの（訂正。他ファイルの docstring は書き換えない）

1. ★**「`dvd_sub_of_phi_intCast` の `hint`/`hrec` を `∀ m` で供給する」は不要。**
   `HasseArfCongruence.lean:112` は `∀ m` で受けるが、実際に要るのは **2 点だけ**である。
   本波の設定では跡びの列 `u` が `s ≤ k` までしか定義されないので
   `∀ m` 版は**そのままでは使えない**。⇒ 点列版 `dvd_sub_of_phi_eq` を本ファイルに立てた。
   （木の版を否定しているのではない。仮説が強いだけであり、両方真である。）
2. ★前波の自分の見積もり「残りは `|G_{u m}| = p^{k+1−m}` の形だけ」は
   半分しか当たっていなかった。`m = k` では `G_{u k + 1}` が**自明群**であり、
   `card_eq_pow_of_mem_Ioc`（上端に `u (m+1)` を使う）では届かない。
   ⇒ `card_eq_pow_of_minimal` で ★**`s = k+1` の埋めを `u k` ではなく `i` にする**と
   `m = k+1`（自明群、`p^0 = 1`）まで一本で扱える。

## 出したもの

| 宣言 | 内容 |
|---|---|
| ★★★`card_eq_pow_of_minimal` | ★抽象核（純群論）`|G_i| = p^{k+1−m}`、`m = k+1` も含む |
| ★★`subgroup_ne_of_jump` | ★抽象核（純群論）`hne : G_{u m} ≠ G_{u m + 1}` |
| ★★`dvd_sub_of_phi_eq` | ★抽象核（純実数・整数）1 点での合同 |
| ★★★★★★★`dvd_sub_jump_of_norm` | **`harith` (3)** |

## ★残り（正確に、`file:line` つき）

`JumpFromValueGroup.lean:272-276` の `harith` 4 条件のうち

* (2) `∀ m<k, u m < u (m+1)` —— `JumpStrictMono.lean:308 jump_lt_succ`（(1) から従う）
* (3) `p^{m+1} ∣ u(m+1) − u m` —— ★**本ファイル**
* (1) `1 ≤ u 0` と (4) `(p−1)·u k ≤ p^{k+1}·e` —— ★**残る**。
  `k = 0`（次数 `p`）では `JumpStrictMono.lean:333 one_le_jump_of_zero` と
  `WildBreakUpperBound.lean sub_one_mul_le_of_totallyRamified` で定理化済み。
  ★**`p^{k+1}` 次への一般化が未測定**である（本波では降りていない）。

★本ファイルの `hjump` は `∀ s ≤ k` であることに注意。
`s = k+1` では `g^{p^{k+1}} = 1` なので成り立たない（`HerbrandRecurrence.lean` の訂正を見よ）。

## 逸脱の記録

1. `G := M ≃ₐ[K] M` に固定。
2. 跡びの列は `u : ℕ → ℕ`。結論だけ `ℤ` に上げる（差を取るため）。
3. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

namespace IntegerNorm

/-! ## §1 抽象核（純群論）—— `i = u m` ちょうどのときの `|G_i|` -/

section Group

variable {G : Type*} [Group G] [Finite G]

/-- ★★★**抽象核（純群論）** —— `m` が「`i ≤ u m` なる最小の段」なら
`Nat.card Gr = p^{k+1−m}`。

`RamCard.card_eq_pow_of_mem_iff`(`RamificationSubgroupCard.lean:174`)の包みで、
★**`m = k+1`（すなわち `Gr` が自明群）まで含めて一本にした**もの。

★★`s = k+1` の建て直しは `u' s := if s ≤ k then u s else i` ——
`HerbrandRecurrence.lean` では `u k` で埋めたが、ここでは **`i` で埋める**。
こうすると `m = k+1` も扱え、`|G_{u k + 1}| = p^0 = 1` が出る。
（★`g^{p^{k+1}} = 1 ∈ Gr` は常に真なので、右辺も常に真にしておく必要がある。） -/
theorem card_eq_pow_of_minimal {p k : ℕ} (hp : p.Prime) {g : G}
    (hg : orderOf g = p ^ (k + 1)) (htop : ∀ σ : G, σ ∈ Subgroup.zpowers g)
    {Gr : Subgroup G} {u : ℕ → ℕ} {i m : ℕ} (hmk : m ≤ k + 1)
    (hmem : ∀ s : ℕ, s ≤ k → (g ^ p ^ s ∈ Gr ↔ i ≤ u s))
    (hhi : m ≤ k → i ≤ u m)
    (hlo : ∀ s : ℕ, s < m → s ≤ k → u s < i) :
    Nat.card Gr = p ^ (k + 1 - m) := by
  classical
  set u' : ℕ → ℤ := fun s => if s ≤ k then ((u s : ℕ) : ℤ) else ((i : ℕ) : ℤ) with hu'
  have hu'le : ∀ s : ℕ, s ≤ k → u' s = ((u s : ℕ) : ℤ) := by
    intro s hs; simp [hu', hs]
  have hone : g ^ p ^ (k + 1) = 1 := by rw [← hg]; exact pow_orderOf_eq_one g
  exact RamCard.card_eq_pow_of_mem_iff (G := G) hp hg htop
    (Gr := Gr) (u := u') (i := (i : ℤ)) (m := m) hmk
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
        · intro _; exact le_rfl
        · intro _; rw [hone]; exact Subgroup.one_mem Gr)
    (by
      rcases Nat.lt_or_ge k m with hm | hm
      · have hmk' : m = k + 1 := by omega
        subst hmk'
        have hnot : ¬ (k + 1 ≤ k) := by omega
        simp only [hu', hnot, if_false]
        exact le_rfl
      · rw [hu'le m hm]; exact_mod_cast hhi hm)
    (fun s _ hs => by
      by_contra hcon
      have hsm : s < m := by omega
      have hsk : s ≤ k := by omega
      rw [hu'le s hsk] at hs
      have h1 : i ≤ u s := by exact_mod_cast hs
      have h2 : u s < i := hlo s hsm hsk
      omega)

/-- ★★**抽象核（純群論）** —— `G_{u m} ≠ G_{u m + 1}`（`hne`）。

位数が `p^{k+1−m}` と `p^{k−m}` で違うから。★`p ≥ 2` が効いている。
★★**`m = k` でも成り立つ**（`|G_{u k}| = p`、`|G_{u k + 1}| = 1`）。 -/
theorem subgroup_ne_of_jump {p k : ℕ} (hp : p.Prime) {g : G}
    (hg : orderOf g = p ^ (k + 1)) (htop : ∀ σ : G, σ ∈ Subgroup.zpowers g)
    {Gr Gr' : Subgroup G} {u : ℕ → ℕ} (hsmono : ∀ s t : ℕ, s < t → u s < u t)
    {m : ℕ} (hmk : m ≤ k)
    (hmem : ∀ s : ℕ, s ≤ k → (g ^ p ^ s ∈ Gr ↔ u m ≤ u s))
    (hmem' : ∀ s : ℕ, s ≤ k → (g ^ p ^ s ∈ Gr' ↔ u m + 1 ≤ u s)) :
    Gr ≠ Gr' := by
  have h1 : Nat.card Gr = p ^ (k + 1 - m) :=
    card_eq_pow_of_minimal hp hg htop (by omega) hmem (fun _ => le_rfl)
      (fun s hs _ => hsmono s m hs)
  have h2 : Nat.card Gr' = p ^ (k + 1 - (m + 1)) :=
    card_eq_pow_of_minimal hp hg htop (by omega) hmem'
      (fun h => by have := hsmono m (m + 1) (by omega); omega)
      (fun s hs _ => by
        rcases Nat.eq_or_lt_of_le (Nat.lt_succ_iff.mp hs) with h | h
        · subst h; omega
        · have := hsmono s m h; omega)
  intro heq
  rw [heq, h2] at h1
  have hinj := Nat.pow_right_injective hp.two_le h1
  omega

end Group

/-! ## §2 抽象核（純実数・整数）—— 1 点での合同 -/

section Arith

/-- ★★**抽象核** —— `x, y` が整数値で `y = x + (b−a)/p^{m+1}` なら `p^{m+1} ∣ b−a`。

★★木の `HasseArfCongruence.lean:112 dvd_sub_of_phi_intCast` は
`hint`・`hrec` を **`∀ m`** で要求するが、実際に要るのは **2 点だけ**である。
本波では `∀ m` 版の入力が `m > k` で作れない（跡びの列は `s ≤ k` まで）ので、
★**点列版をここに立てる**。中身は実数の除算と `Int` への戻しだけである。 -/
theorem dvd_sub_of_phi_eq {p m : ℕ} (hp : p ≠ 0) {a b : ℕ} {ja jb : ℤ} {x y : ℝ}
    (hx : x = (ja : ℝ)) (hy : y = (jb : ℝ))
    (hrec : y = x + ((b : ℝ) - (a : ℝ)) / ((p : ℝ) ^ (m + 1))) :
    ((p : ℤ) ^ (m + 1)) ∣ ((b : ℤ) - (a : ℤ)) := by
  have hpR : ((p : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hp
  have hpow : ((p : ℝ)) ^ (m + 1) ≠ 0 := pow_ne_zero _ hpR
  rw [hx, hy] at hrec
  have h : ((b : ℝ) - (a : ℝ)) = ((jb : ℝ) - (ja : ℝ)) * ((p : ℝ) ^ (m + 1)) := by
    field_simp at hrec
    linarith
  have h2 : ((((b : ℤ) - (a : ℤ)) : ℤ) : ℝ)
      = ((((jb - ja) * (p : ℤ) ^ (m + 1)) : ℤ) : ℝ) := by
    push_cast
    linarith
  have h3 : ((b : ℤ) - (a : ℤ)) = (jb - ja) * (p : ℤ) ^ (m + 1) := by exact_mod_cast h2
  exact ⟨jb - ja, by rw [h3]; ring⟩

end Arith

/-! ## §3 具体層 —— `harith` (3) をノルムだけの仮説から -/

section Norm

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]
  [FiniteDimensional K M]

/-- ★★★★★★★**`harith` (3)** —— Hasse–Arf の合同 `p^{m+1} ∣ u(m+1) − u m`。

仮説は `K`・`M`・`π`・`π_K`・`g` とノルムの不等式だけ。
★環・イデアル・分岐群の語を 1 つも仮説に持たない。 -/
theorem dvd_sub_jump_of_norm
    {π : M} {πK : K} {n p k m : ℕ} {u : ℕ → ℕ} {g : M ≃ₐ[K] M}
    (hn1 : 1 < n) (hp : p.Prime)
    (hiso : ∀ (σ : M ≃ₐ[K] M) (z : M), ‖σ • z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ j : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * j))
    (hval : ∀ z : M, z ≠ 0 → ∃ j : ℤ, ‖z‖ = ‖π‖ ^ j) (hπmem : π ∈ integerSubring M)
    (hjump : ∀ s : ℕ, s ≤ k → ‖(g ^ p ^ s) π - π‖ = ‖π‖ ^ (u s + 1))
    (hsmono : ∀ s t : ℕ, s < t → u s < u t) (hu0 : 1 ≤ u 0)
    (hg : orderOf g = p ^ (k + 1)) (htop : ∀ σ : M ≃ₐ[K] M, σ ∈ Subgroup.zpowers g)
    (hK0 : 0 < ‖algebraMap K M πK‖) (hK1 : ‖algebraMap K M πK‖ < 1)
    (hvalKK : ∀ a : K, a ≠ 0 → ∃ j : ℤ,
      ‖algebraMap K M a‖ = ‖algebraMap K M πK‖ ^ j)
    (hKmem : πK ∈ baseIntegerSubring K M) (hpM : ‖((p : ℕ) : M)‖ < 1)
    (hmk : m + 1 ≤ k) :
    ((p : ℤ) ^ (m + 1)) ∣ ((u (m + 1) : ℕ) : ℤ) - ((u m : ℕ) : ℤ) := by
  letI := isLocalRing_integerSubring (M := M)
  letI := integerMulSemiringAction hiso
  letI := isDiscreteValuationRing_integerSubring hπ0 hπ1 hval hπmem
  have hmono : ∀ s t : ℕ, s ≤ t → u s ≤ u t := by
    intro s t hst
    rcases Nat.eq_or_lt_of_le hst with h | h
    · exact le_of_eq (by rw [h])
    · exact le_of_lt (hsmono s t h)
  have hgen : ∀ h : M ≃ₐ[K] M, ∃ j : ℤ, h = g ^ j := by
    intro h
    obtain ⟨j, hj⟩ := htop h
    exact ⟨j, hj.symm⟩
  have hbr : ‖g π - π‖ = ‖π‖ ^ (u 0 + 1) := by
    have := hjump 0 (by omega)
    simpa using this
  have hmemG : ∀ i s : ℕ, s ≤ k →
      (g ^ p ^ s ∈ lowerRamificationGroup ↥(integerSubring M) (M ≃ₐ[K] M) i ↔ i ≤ u s) :=
    fun i s hs => mem_lowerRamificationGroup_iff_jump hiso hπ0 hπ1 hn hvalK hval hπmem i s
      (hjump s hs)
  have hne : ∀ r : ℕ, r ≤ k →
      lowerRamificationGroup ↥(integerSubring M) (M ≃ₐ[K] M) (u r) ≠
        lowerRamificationGroup ↥(integerSubring M) (M ≃ₐ[K] M) (u r + 1) :=
    fun r hr => subgroup_ne_of_jump hp hg htop hsmono hr
      (fun s hs => hmemG (u r) s hs) (fun s hs => hmemG (u r + 1) s hs)
  obtain ⟨ja, hja⟩ := exists_herbrandPhiGroup_natCast_of_norm hn1 hiso hπ0 hπ1 hn hvalK hval
    hπmem hbr hu0 hgen hK0 hK1 hvalKK hKmem hp hpM (hne m (by omega))
  obtain ⟨jb, hjb⟩ := exists_herbrandPhiGroup_natCast_of_norm hn1 hiso hπ0 hπ1 hn hvalK hval
    hπmem hbr hu0 hgen hK0 hK1 hvalKK hKmem hp hpM (hne (m + 1) (by omega))
  have hrec := herbrandPhiGroup_step hp hiso hπ0 hπ1 hn hvalK hval hπmem hjump hmono hg htop hmk
  exact dvd_sub_of_phi_eq (p := p) (m := m) hp.ne_zero (ja := (ja : ℤ)) (jb := (jb : ℤ))
    (by rw [hja]; push_cast; ring) (by rw [hjb]; push_cast; ring) hrec

end Norm


/-! ## `.src` と 公理 -/

def dvd_sub_jump_of_norm.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def card_eq_pow_of_minimal.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 14, item := "Definition 6.1", sectionId := "def-6-1" }

#print axioms card_eq_pow_of_minimal
#print axioms subgroup_ne_of_jump
#print axioms dvd_sub_of_phi_eq
#print axioms dvd_sub_jump_of_norm

end IntegerNorm

end ABC3.Found.PGC
