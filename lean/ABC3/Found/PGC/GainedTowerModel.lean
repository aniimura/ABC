import ABC3.Found.PGC.GainedTowerStep

/-!
# [pGC] `k ≥ 1` の塔の模型 —— `hu`(Hasse–Arf)と `E` の族を供給する

`Found/PGC/GainedTowerStep.lean` の最終出口
`exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_jumps` は、
`k ≥ 1` でも空虚でない形で `hstep` を供給し切っている。残っていたのは仮説 4 本:

| 仮説 | 内容 | 本ファイル |
|---|---|---|
| `hu` | `p^m ≤ u m`(`m ≤ k`) | ★§2 で**落とした**(Hasse–Arf の合同から) |
| `hjump` | `‖τ^{p^j}π − π‖ ≤ ‖π‖^{u j + 1}` | ★§3 で**底の 1 本に落とした** |
| `hfixj` | `τ^{p^j}` が `E j` を固定 | ★§5 で**落とした**(固定体を取る) |
| `hdegj` | `[M : E j] = p^{k+1−j}` | ★§5 で `htopj` から出した(Artin) |

★残るのは `htopj`(`M = (E j)(π)`)と `hvalj`(値群)で、これは
**全分岐という付値の情報**そのものであり、群論では出ない(§6 に理由を書く)。

## ★★配られた字面の検算 —— `(p,k,e,u) = (3,1,2,(2,8))` は本物か

`GainedTowerStep.Numeric.tz` の docstring が言う実データを、
**原典を見ずに古典的定理(Serre, *Corps Locaux* IV §4)から独立に再計算した**:

`K = ℚ₃(ζ₂₇)` の `Gal(K/ℚ₃) ≅ (ℤ/27)^×`(位数 18、巡回)の**下付き**分岐群は

> `G_i = Gal(K/ℚ_p(ζ_{p^m}))`  (`p^{m−1} ≤ i ≤ p^m − 1`)

なので `p = 3`, `n = 3` では

| `i` | `0` | `1,2` | `3,…,8` | `≥ 9` |
|---|---|---|---|---|
| `G_i` | 位数 18 | ★位数 9 = `Gal(K/ℚ₃(ζ₃))` | ★位数 3 = `Gal(K/ℚ₃(ζ₉))` | `1` |

下付き番号は部分群と両立する(`H_i = H ∩ G_i`)ので、`E 0 = ℚ₃(ζ₃)`、
`H = Gal(K/E 0)`(位数 `9 = p^{k+1}`、`k = 1`)の生成元 `τ` について

* `τ ∈ H_2 \ H_3` ⟹ ★`u 0 = 2`
* `τ³ ∈ H_8 \ H_9` ⟹ ★`u 1 = 8`

`e(K/ℚ₃) = 18 = p^{k+1}·e` ⟹ ★`e = 2`(= `e(ℚ₃(ζ₃)/ℚ₃)`)。
★★**3 つとも `tz` の docstring と一致した**(§4 で機械検算)。

★★★さらに `u 1 = 8` は**上下 2 本の独立な不等式に挟まれて一意に決まる**(§4):

* 下から `hstep` の縮小 `p(u 0 + 1) ≤ u 1 + 1` ⟹ `u 1 ≥ 8`
* 上から `hbnd` `(p−1)·u 1 ≤ p^{k+1}·e = 18` ⟹ `u 1 ≤ 9`
* Hasse–Arf の合同 `p ∣ u 1 − u 0` ⟹ `u 1 ≡ 2 (mod 3)` ⟹ ★`u 1 = 8`

★これは 2 経路の一致(`(u₀,u₁) = (2,8)`)を**3 経路目**で裏づける。
-/

open Finset

namespace ABC3.Found.PGC

namespace GainedTowerModel

/-! ## §1 抽象核(ℤ だけ) —— 「合同つきで真に増える列」は指数的に増える

★分岐・付値・Galois・体の語彙が 1 語も出ない。 -/

section IntCore

/-- 真に増えて差が `d` で割れるなら、差は少なくとも `d` ある。
★`0 < d` は要らない(`d ≤ 0` なら `a + d ≤ a < b`)。 -/
theorem add_le_of_lt_of_dvd {a b d : ℤ} (hlt : a < b) (hdvd : d ∣ b - a) :
    a + d ≤ b := by
  have h : d ≤ b - a := Int.le_of_dvd (by omega) hdvd
  omega

/-- ★★**抽象核** —— `u 0 ≥ 1` と `u m + p^{m+1} ≤ u (m+1)` から `p^m ≤ u m`。

★これが `GainedJumpSeq.exists_jump_seq` の `hu` の中身であり、
入力は「1 から始まって段差が `p^{m+1}` 以上」だけである。 -/
theorem pow_le_of_gaps {p : ℤ} (hp : 0 ≤ p) {u : ℕ → ℤ} (h0 : 1 ≤ u 0)
    (hgap : ∀ m, u m + p ^ (m + 1) ≤ u (m + 1)) : ∀ m, p ^ m ≤ u m := by
  intro m
  induction m with
  | zero => simpa using h0
  | succ n ih =>
      have hpos : (0 : ℤ) ≤ p ^ n := pow_nonneg hp n
      have := hgap n
      linarith

/-- `k` までしか段差を仮定しない版(実際に要るのはこれ)。 -/
theorem pow_le_of_gaps_le {p : ℤ} (hp : 0 ≤ p) {u : ℕ → ℤ} {k : ℕ} (h0 : 1 ≤ u 0)
    (hgap : ∀ m, m < k → u m + p ^ (m + 1) ≤ u (m + 1)) : ∀ m, m ≤ k → p ^ m ≤ u m := by
  intro m
  induction m with
  | zero => intro _; simpa using h0
  | succ n ih =>
      intro hn
      have hnk : n < k := by omega
      have hpos : (0 : ℤ) ≤ p ^ n := pow_nonneg hp n
      have := hgap n hnk
      have := ih (by omega)
      linarith

/-- ★★★**抽象核** —— Hasse–Arf の合同から `hu` が出る。

入力は純粋に ℤ の列の性質:
`1 ≤ u 0`、`u m < u (m+1)`、`p^{m+1} ∣ u (m+1) − u m`(`m < k`)。

★★分岐の言葉では「下付きの跳び `u m` は真に増え、
`m+1` 番目の跳びは `p^{m+1}` を法として前と合同(Hasse–Arf ⇔ 上付きが整数)」。
★この 2 つだけで `p^m ≤ u m` が出る。原典が Hasse–Arf を引くのはここである。 -/
theorem pow_le_of_strictMono_of_dvd {p : ℤ} (hp : 0 ≤ p) {u : ℕ → ℤ} {k : ℕ} (h0 : 1 ≤ u 0)
    (hlt : ∀ m, m < k → u m < u (m + 1))
    (hdvd : ∀ m, m < k → p ^ (m + 1) ∣ u (m + 1) - u m) :
    ∀ m, m ≤ k → p ^ m ≤ u m :=
  pow_le_of_gaps_le hp h0 fun m hm => add_le_of_lt_of_dvd (hlt m hm) (hdvd m hm)

end IntCore

/-! ## §2 点 4 —— `hu`(Hasse–Arf)を落とした出口 -/

section DropHu

variable {M : Type*} [NormedField M] [IsUltrametricDist M]
variable {F : Type*} [Field F] [Algebra F M]

/-- ★★★★★★点 4 —— `GainedTowerStep` の最終出口から仮説 `hu` を落とした形。

`hu`(`p^m ≤ u m`)の代わりに**跳びの列そのものの性質**を取る:

* `hu0 : 1 ≤ u 0`(底の跳びは正 ＝ 全分岐が野性的)
* `hult : u m < u (m+1)`(跳びは真に増える)
* `hudvd : p^{m+1} ∣ u (m+1) − u m`(★Hasse–Arf の合同)

★他の仮説は `..._of_tower_jumps` と同一である。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_hasseArf
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M} {u : ℕ → ℤ}
    {E : ℕ → Type*} [∀ j, Field (E j)] [∀ j, Algebra (E j) M]
    (τ : M →+ M) (σ : M →ₐ[F] M) (s : ℕ → (M →+* M)) (he : 0 < e)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z)
    (hs : ∀ j, j < k → ∀ z : M, s j z = (⇑τ)^[p ^ j] z)
    (hu0 : 1 ≤ u 0)
    (hult : ∀ m, m < k → u m < u (m + 1))
    (hudvd : ∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m)
    (hbnd : ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ))
    (hjump : ∀ j, j < k → ‖s j π - π‖ ≤ ‖π‖ ^ (u j + 1))
    (hdegj : ∀ j, j < k → (minpoly (E j) π).natDegree = p ^ (k + 1 - j))
    (htopj : ∀ j, j < k → Algebra.adjoin (E j) ({π} : Set M) = ⊤)
    (hvalj : ∀ j, j < k → ∀ d : E j, d ≠ 0 →
        ∃ m : ℤ, ‖algebraMap (E j) M d‖ = ‖π‖ ^ (((p ^ (k + 1 - j) : ℕ) : ℤ) * m))
    (hfixj : ∀ j, j < k → ∀ d : E j, s j (algebraMap (E j) M d) = algebraMap (E j) M d)
    (hti : (i : ℤ) = u k)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ :=
  GainedTowerStep.exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_jumps
    τ σ s he hσ hs
    (pow_le_of_strictMono_of_dvd (by positivity) hu0 hult hudvd)
    hbnd hjump hdegj htopj hvalj hfixj hti hnormp heM hbreak hval hdeg htop x

end DropHu

/-! ## §3 点 1 の数値層 —— `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` を**独立に再計算**して字面を検算する

★ここは体もノルムも使わない。`(p,k,e,u) = (3,1,2,(2,8))` という**配られた字面**が
古典的定理(Serre IV §4 の `ℚ_p(ζ_{p^n})` の分岐群)と合うかだけを機械で確かめる。 -/

namespace Zeta27

/-- `Gal(ℚ₃(ζ₂₇)/ℚ₃)` の下付き分岐群 `G_i` の位数(Serre IV §4:
`G_i = Gal(K/ℚ_p(ζ_{p^m}))` for `p^{m−1} ≤ i ≤ p^m − 1`)。

`i = 0` : 18、`i = 1,2` : 9、`i = 3,…,8` : 3、`i ≥ 9` : 1。 -/
def gOrd (i : ℕ) : ℕ :=
  if i = 0 then 18 else if i ≤ 2 then 9 else if i ≤ 8 then 3 else 1

/-- `H = Gal(K/ℚ₃(ζ₃))`(位数 `9 = p^{k+1}`)の下付き分岐群の位数。
★下付き番号は部分群と両立する(`H_i = H ∩ G_i`)ので `min (gOrd i) 9` になる。 -/
def hOrd (i : ℕ) : ℕ := min (gOrd i) 9

/-- `H_i = H ∩ G_i` を `min` で書けること(位数が全部 `9` の約数で入れ子だから)。 -/
theorem hOrd_eq (i : ℕ) : hOrd i = if i ≤ 8 then (if i ≤ 2 then 9 else 3) else 1 := by
  unfold hOrd gOrd
  split_ifs <;> omega

/-- ★`τ`(`H` の生成元、位数 9)の跳び: `τ ∈ H_2 \ H_3` なので `u 0 = 2`。 -/
theorem hOrd_two_ne_hOrd_three : hOrd 2 = 9 ∧ hOrd 3 = 3 := by decide

/-- ★`τ³`(位数 3)の跳び: `τ³ ∈ H_8 \ H_9` なので `u 1 = 8`。 -/
theorem hOrd_eight_ne_hOrd_nine : hOrd 8 = 3 ∧ hOrd 9 = 1 := by decide

/-- `e(K/ℚ₃) = [K:ℚ₃] = φ(27) = 18`(全分岐)。 -/
theorem totient_27 : Nat.totient 27 = 18 := by decide

/-- ★★`heM` の指数 `p^{k+1}·e` が実際の `e(K/ℚ₃)` と合う: `3^{1+1}·2 = 18 = φ(27)`。
★これで `e = 2`(= `e(ℚ₃(ζ₃)/ℚ₃)`)が**強制される**。 -/
theorem eM_eq : 3 ^ (1 + 1) * 2 = Nat.totient 27 := by decide

/-- Hilbert の差積公式 `d = Σ_i (|G_i| − 1)` を表から計算すると `45`。 -/
theorem different_from_table : ∑ i ∈ Finset.range 9, (gOrd i - 1) = 45 := by decide

/-- ★同じ `45` が古典的な閉じた形 `p^{n−1}(pn − n − 1)`(`p = 3`, `n = 3`)からも出る。
★★**これが分岐群の表 `gOrd` の独立な検算である**(表を 1 段でも間違えると合わない)。 -/
theorem different_classical : 3 ^ (3 - 1) * (3 * 3 - 3 - 1) = 45 := by decide

/-- ★★塔の差積公式 `d(K/ℚ₃) = d(K/E₀) + e(K/E₀)·d(E₀/ℚ₃)` も閉じる:
`d(K/E₀) = Σ_i (|H_i| − 1) = 36`、`e(K/E₀) = 9`、`d(E₀/ℚ₃) = e − 1 = 1`(順分岐)。
★`36 + 9·1 = 45`。★**3 本目の独立な検算**。 -/
theorem different_tower :
    (∑ i ∈ Finset.range 9, (hOrd i - 1)) + 9 * (2 - 1) = 45 := by decide

/-- 実データの跳びの列 `u = (2, 8, …)`。`k = 1` なので使うのは `u 0`, `u 1` だけ。 -/
def uz : ℕ → ℤ
  | 0 => 2
  | 1 => 8
  | (n + 2) => 8 * 3 ^ (n + 1)

theorem uz_zero : uz 0 = 2 := rfl
theorem uz_one : uz 1 = 8 := rfl

/-- ★§2 の仮説 `hu0`。 -/
theorem uz_hu0 : (1 : ℤ) ≤ uz 0 := by norm_num [uz]

/-- ★§2 の仮説 `hult`(`k = 1`)。 -/
theorem uz_hult : ∀ m, m < 1 → uz m < uz (m + 1) := by
  intro m hm
  interval_cases m
  · norm_num [uz]

/-- ★§2 の仮説 `hudvd`(Hasse–Arf の合同、`k = 1`): `3 ∣ 8 − 2`。 -/
theorem uz_hudvd : ∀ m, m < 1 → (3 : ℤ) ^ (m + 1) ∣ uz (m + 1) - uz m := by
  intro m hm
  interval_cases m
  · norm_num [uz]

/-- ★`hbnd`(頂点の古典的上界): `(3−1)·8 = 16 ≤ 3^{1+1}·2 = 18`。 -/
theorem uz_hbnd : ((3 : ℤ) - 1) * uz 1 ≤ (3 : ℤ) ^ (1 + 1) * (2 : ℤ) := by norm_num [uz]

/-- ★§2 が実際に `hu` を出すこと(`3^m ≤ uz m`, `m ≤ 1`)を機械で確かめる。 -/
theorem uz_hu : ∀ m, m ≤ 1 → (3 : ℤ) ^ m ≤ uz m :=
  pow_le_of_strictMono_of_dvd (by norm_num) uz_hu0 uz_hult uz_hudvd

/-- ★★★**3 経路目** —— `u 1 = 8` は 3 本の不等式・合同で**一意に決まる**。

* 下から: `hstep` の縮小 `p·(u₀+1) ≤ u₁+1`(`9 ≤ u₁+1`)
* 上から: `hbnd` `(p−1)·u₁ ≤ p^{k+1}·e = 18`(`u₁ ≤ 9`)
* Hasse–Arf: `p ∣ u₁ − u₀`(`u₁ ≡ 2 mod 3`)

★★これは既存の 2 経路(実測の 2 系列)とは独立で、`(u₀,u₁) = (2,8)` を裏づける。 -/
theorem uz_one_forced {u1 : ℤ} (hlow : (3 : ℤ) * (2 + 1) ≤ u1 + 1)
    (hhigh : ((3 : ℤ) - 1) * u1 ≤ (3 : ℤ) ^ (1 + 1) * 2) (hcong : (3 : ℤ) ∣ u1 - 2) :
    u1 = 8 := by
  have h9 : ((3 : ℤ) ^ (1 + 1)) = 9 := by norm_num
  rw [h9] at hhigh
  omega

end Zeta27

/-! ## §4 抽象核(純群論) —— 位数 `p^{k+1}` の元の `p^j` 乗の位数 -/

section GroupCore

/-- ★★**抽象核** —— `orderOf g = p^{k+1}` なら `orderOf (g^{p^j}) = p^{k+1−j}`。
★分岐・付値・Galois・体の語彙が 1 語も出ない。 -/
theorem orderOf_pow_pow {G : Type*} [Group G] {p k j : ℕ} (hp : 0 < p) (g : G)
    (hg : orderOf g = p ^ (k + 1)) (hj : j ≤ k + 1) :
    orderOf (g ^ p ^ j) = p ^ (k + 1 - j) := by
  have hne : p ^ j ≠ 0 := (pow_pos hp j).ne'
  have hdvd : p ^ j ∣ orderOf g := by
    rw [hg]; exact pow_dvd_pow p hj
  rw [orderOf_pow_of_dvd hne hdvd, hg, Nat.pow_div hj hp]

/-- ★`zpowers (g^{p^j})` の位数。 -/
theorem card_zpowers_pow_pow {G : Type*} [Group G] {p k j : ℕ} (hp : 0 < p) (g : G)
    (hg : orderOf g = p ^ (k + 1)) (hj : j ≤ k + 1) :
    Nat.card (Subgroup.zpowers (g ^ p ^ j)) = p ^ (k + 1 - j) := by
  rw [Nat.card_zpowers, orderOf_pow_pow hp g hg hj]

end GroupCore

/-! ## §5 点 2 —— 巡回 `p^{k+1}` 次拡大から `E : ℕ → Type*` を作る

★★`IntermediateField` を**使う**が、#59(2 層をまたぐ `rfl`)には当たらない。
理由: 層の間の関係を `rfl` で潰さず、**すべて底 `K` から測る**(`twr j` は
`M ≃ₐ[K] M` の部分群の固定体で、層どうしは比較しない)から。 -/

namespace GaloisTower

variable {K M : Type*} [Field K] [Field M] [Algebra K M]

/-- 塔の層 `E j` —— `⟨g^{p^j}⟩` の固定体。★`[M : E j] = p^{k+1−j}`。 -/
def twr (g : M ≃ₐ[K] M) (p j : ℕ) : IntermediateField K M :=
  IntermediateField.fixedField (Subgroup.zpowers (g ^ p ^ j))

/-- `M ≃ₐ[K] M` の冪は反復合成。 -/
theorem algEquiv_pow_apply (g : M ≃ₐ[K] M) : ∀ (n : ℕ) (z : M), (g ^ n) z = (⇑g)^[n] z := by
  intro n
  induction n with
  | zero => intro z; simp
  | succ m ih =>
      intro z
      rw [pow_succ', Function.iterate_succ_apply', AlgEquiv.mul_apply, ih]

/-- ★★点 2 の `hfixj` —— `g^{p^j}` は `E j` を**点ごとに固定する**(定義から)。 -/
theorem fix_twr (g : M ≃ₐ[K] M) (p j : ℕ) (d : twr g p j) :
    (g ^ p ^ j) (algebraMap (twr g p j) M d) = algebraMap (twr g p j) M d := by
  have hd : (d : M) ∈ IntermediateField.fixedField (Subgroup.zpowers (g ^ p ^ j)) := d.2
  rw [IntermediateField.mem_fixedField_iff] at hd
  exact hd _ (Subgroup.mem_zpowers _)

/-- ★★点 2 の `[M : E j] = p^{k+1−j}`(Artin) -/
theorem finrank_twr [FiniteDimensional K M] {p k j : ℕ} (hp : 0 < p) (g : M ≃ₐ[K] M)
    (hg : orderOf g = p ^ (k + 1)) (hj : j ≤ k + 1) :
    Module.finrank (twr g p j) M = p ^ (k + 1 - j) := by
  rw [twr, IntermediateField.finrank_fixedField_eq_card,
    card_zpowers_pow_pow hp g hg hj]

/-- ★★点 2 の `htopj` —— `π` が底 `K` 上で原始元なら、**どの層の上でも**原始元。
★`IntermediateField.adjoin_eq_top_of_adjoin_eq_top` の直接の帰結で、
層どうしを比較しない(#59 を踏まない)。 -/
theorem topIF_twr [FiniteDimensional K M] (g : M ≃ₐ[K] M) (p j : ℕ) {π : M}
    (htop : Algebra.adjoin K ({π} : Set M) = ⊤) :
    IntermediateField.adjoin (twr g p j) ({π} : Set M) = ⊤ :=
  IntermediateField.adjoin_eq_top_of_adjoin_eq_top (F := K) (E := twr g p j)
    (IntermediateField.adjoin_eq_top_of_algebra K ({π} : Set M) htop)

theorem top_twr [FiniteDimensional K M] (g : M ≃ₐ[K] M) (p j : ℕ) {π : M}
    (htop : Algebra.adjoin K ({π} : Set M) = ⊤) :
    Algebra.adjoin (twr g p j) ({π} : Set M) = ⊤ :=
  (IntermediateField.adjoin_eq_top_iff (F := twr g p j) (E := M)).1 (topIF_twr g p j htop)

/-- ★★★点 2 の `hdegj` —— `htopj` と Artin の次数から出る。 -/
theorem deg_twr [FiniteDimensional K M] {p k j : ℕ} (hp : 0 < p) (g : M ≃ₐ[K] M)
    (hg : orderOf g = p ^ (k + 1)) (hj : j ≤ k + 1) {π : M}
    (htop : Algebra.adjoin K ({π} : Set M) = ⊤) :
    (minpoly (twr g p j) π).natDegree = p ^ (k + 1 - j) := by
  have h := (Field.primitive_element_iff_minpoly_natDegree_eq
    (F := twr g p j) (E := M) π).1 (topIF_twr g p j htop)
  rw [h, finrank_twr hp g hg hj]

end GaloisTower

/-! ## §6 点 2 の出口 —— 巡回 `p^{k+1}` 次拡大だけから塔を作った形

★★★仮説から消えたもの: `E`(体の族)、`hdegj`、`htopj`、`hfixj`、
さらに ★**上の体 `F` 自身**(`= twr g p k`、`[M:F] = p`)と `hdeg` / `htop`。
★残るのは `hvalj`(値群 ＝ 全分岐)だけである。 -/

section CyclicExit

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

open GaloisTower

/-- ★★★★★★点 2 —— **入力は巡回群 1 個**(`orderOf g = p^{k+1}`)と `π` の原始性、
そして値群の条件 `hvalj` だけ。体の族 `E` は `⟨g^{p^j}⟩` の固定体として**構成される**。

* `F := twr g p k`(`[M : F] = p`)は**結論の側にも現れる**。
* `hdeg` / `htop` / `hdegj` / `htopj` / `hfixj` / `hu` は**すべて落ちた**。
* ★`hvalj`(層 `j` の値群が `p^{k+1−j}·ℤ`)だけが残る。これは全分岐そのもので、
  群論からは出ない(§7 参照)。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic
    {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M]
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M} {u : ℕ → ℤ}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (τ : M →+ M) (hτ : ∀ z : M, τ z = g z)
    (s : ℕ → (M →+* M)) (hsg : ∀ j, ∀ z : M, s j z = (g ^ p ^ j) z)
    (σ : M →ₐ[twr g p k] M) (hσg : ∀ z : M, σ z = (g ^ p ^ k) z)
    (he : 0 < e)
    (hu0 : 1 ≤ u 0)
    (hult : ∀ m, m < k → u m < u (m + 1))
    (hudvd : ∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m)
    (hbnd : ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ))
    (hjump : ∀ j, j < k → ‖s j π - π‖ ≤ ‖π‖ ^ (u j + 1))
    (hvalj : ∀ j, j ≤ k → ∀ d : twr g p j, d ≠ 0 →
        ∃ m : ℤ, ‖algebraMap (twr g p j) M d‖ = ‖π‖ ^ (((p ^ (k + 1 - j) : ℕ) : ℤ) * m))
    (hti : (i : ℤ) = u k)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (htop : Algebra.adjoin K ({π} : Set M) = ⊤) (x : M) :
    ∃ y : twr g p k, ‖x - algebraMap (twr g p k) M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have hτf : (⇑τ : M → M) = ⇑g := funext hτ
  have hs : ∀ j, j < k → ∀ z : M, s j z = (⇑τ)^[p ^ j] z := by
    intro j _ z
    rw [hsg j z, hτf, algEquiv_pow_apply]
  have hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z := by
    intro z
    rw [hσg z, hτf, algEquiv_pow_apply]
  have hvalF : ∀ c : twr g p k, c ≠ 0 →
      ∃ m : ℤ, ‖algebraMap (twr g p k) M c‖ = ‖π‖ ^ ((p : ℤ) * m) := by
    intro c hc
    obtain ⟨m, hm⟩ := hvalj k le_rfl c hc
    refine ⟨m, ?_⟩
    rw [hm]
    congr 2
    simp
  have hdegF : (minpoly (twr g p k) π).natDegree = p := by
    have := deg_twr (p := p) (k := k) (j := k) hp'.pos g hg (by omega) htop
    simpa using this
  refine exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_hasseArf
    (E := fun j => twr g p j) τ σ s he hσ hs hu0 hult hudvd hbnd hjump
    (fun j hj => deg_twr hp'.pos g hg (by omega) htop)
    (fun j _ => top_twr g p j htop)
    (fun j hj => hvalj j (le_of_lt hj))
    (fun j _ d => ?_)
    hti hnormp heM hbreak hvalF hdegF (top_twr g p k htop) x
  rw [hsg j (algebraMap (twr g p j) M d)]
  exact fix_twr g p j d

end CyclicExit

/-! ## §7 `k = 1` と ★`k = 2` の数値模型を**もう 2 つ**独立に作る

★`Zeta27` の 1 例だけでは「たまたま合った」を排除できないので、
`ℚ₃(ζ₈₁)` の塔から `(p,k,e,u)` をもう 2 組取って同じ検算を通す。
★★**`k = 2` が取れる**(`ℚ₃(ζ₈₁)/ℚ₃(ζ₃)`、`[M:K] = 27 = p^{k+1}`)。 -/

namespace Zeta81

/-- `Gal(ℚ₃(ζ₈₁)/ℚ₃)` の下付き分岐群の位数(`n = 4`)。
`i=0` : 54、`i=1,2` : 27、`i=3..8` : 9、`i=9..26` : 3、`i≥27` : 1。 -/
def gOrd (i : ℕ) : ℕ :=
  if i = 0 then 54 else if i ≤ 2 then 27 else if i ≤ 8 then 9 else if i ≤ 26 then 3 else 1

/-- `e(ℚ₃(ζ₈₁)/ℚ₃) = φ(81) = 54`。 -/
theorem totient_81 : Nat.totient 81 = 54 := by decide

/-- Hilbert の差積 `Σ_i(|G_i|−1)` は `189`。 -/
theorem different_from_table : ∑ i ∈ Finset.range 27, (gOrd i - 1) = 189 := by decide

/-- ★同じ `189` が閉じた形 `p^{n−1}(pn−n−1)`(`p=3`, `n=4`)からも出る。 -/
theorem different_classical : 3 ^ (4 - 1) * (3 * 4 - 4 - 1) = 189 := by decide

/-! ### `k = 1` の 2 例目 —— `K = ℚ₃(ζ₉)`(`e = 6`)、`u = (8, 26)` -/

/-- `H = Gal(M/ℚ₃(ζ₉))`(位数 `9 = p^{k+1}`, `k = 1`)の下付き分岐群の位数。 -/
def hOrd9 (i : ℕ) : ℕ := min (gOrd i) 9

/-- ★`τ ∈ H_8 \ H_9`、`τ³ ∈ H_26 \ H_27` ⟹ `u = (8, 26)`。 -/
theorem hOrd9_breaks : hOrd9 8 = 9 ∧ hOrd9 9 = 3 ∧ hOrd9 26 = 3 ∧ hOrd9 27 = 1 := by decide

/-- ★塔の差積: `d(M/K) = 108`、`e(M/K)·d(K/ℚ₃) = 9·9 = 81`、和は `189`。 -/
theorem different_tower9 :
    (∑ i ∈ Finset.range 27, (hOrd9 i - 1)) + 9 * (3 ^ (2 - 1) * (3 * 2 - 2 - 1)) = 189 := by
  decide

/-- `heM` の指数: `p^{k+1}·e = 3^{1+1}·6 = 54 = φ(81)`。 -/
theorem eM_eq9 : 3 ^ (1 + 1) * 6 = Nat.totient 81 := by decide

/-- `hbnd`: `(3−1)·26 = 52 ≤ 3^{1+1}·6 = 54`。 -/
theorem hbnd9 : ((3 : ℤ) - 1) * 26 ≤ (3 : ℤ) ^ (1 + 1) * 6 := by norm_num

/-- ★Hasse–Arf の合同 `3 ∣ 26 − 8`。 -/
theorem hudvd9 : (3 : ℤ) ^ (0 + 1) ∣ (26 : ℤ) - 8 := by norm_num

/-- ★★**2 例目でも頂点が一意に決まる** —— 下 `3(8+1) ≤ u₁+1`、上 `2u₁ ≤ 54`、
合同 `u₁ ≡ 8 (mod 3)` ⟹ `u₁ = 26`。 -/
theorem u_one_forced9 {u1 : ℤ} (hlow : (3 : ℤ) * (8 + 1) ≤ u1 + 1)
    (hhigh : ((3 : ℤ) - 1) * u1 ≤ (3 : ℤ) ^ (1 + 1) * 6) (hcong : (3 : ℤ) ∣ u1 - 8) :
    u1 = 26 := by
  have h9 : ((3 : ℤ) ^ (1 + 1)) = 9 := by norm_num
  rw [h9] at hhigh
  omega

/-! ### ★★`k = 2` —— `K = ℚ₃(ζ₃)`(`e = 2`)、`u = (2, 8, 26)` -/

/-- `H = Gal(M/ℚ₃(ζ₃))`(位数 `27 = p^{k+1}`, ★`k = 2`)の下付き分岐群の位数。 -/
def hOrd27 (i : ℕ) : ℕ := min (gOrd i) 27

/-- ★`u = (2, 8, 26)`(3 段の跳び)。 -/
theorem hOrd27_breaks :
    hOrd27 2 = 27 ∧ hOrd27 3 = 9 ∧ hOrd27 8 = 9 ∧ hOrd27 9 = 3 ∧
      hOrd27 26 = 3 ∧ hOrd27 27 = 1 := by decide

/-- ★塔の差積: `d(M/K) = 162`、`e(M/K)·d(K/ℚ₃) = 27·1 = 27`、和は `189`。 -/
theorem different_tower27 :
    (∑ i ∈ Finset.range 27, (hOrd27 i - 1)) + 27 * (2 - 1) = 189 := by decide

/-- `heM` の指数: `p^{k+1}·e = 3^{2+1}·2 = 54 = φ(81)`。★`k = 2`。 -/
theorem eM_eq27 : 3 ^ (2 + 1) * 2 = Nat.totient 81 := by decide

/-- `k = 2` の跳びの列。 -/
def uz : ℕ → ℤ
  | 0 => 2
  | 1 => 8
  | 2 => 26
  | (n + 3) => 26 * 3 ^ (n + 1)

/-- ★`k = 2` で §2 の仮説がすべて成り立つ。 -/
theorem uz_hu0 : (1 : ℤ) ≤ uz 0 := by norm_num [uz]

theorem uz_hult : ∀ m, m < 2 → uz m < uz (m + 1) := by
  intro m hm; interval_cases m <;> norm_num [uz]

theorem uz_hudvd : ∀ m, m < 2 → (3 : ℤ) ^ (m + 1) ∣ uz (m + 1) - uz m := by
  intro m hm; interval_cases m <;> norm_num [uz]

/-- `hbnd`(`k = 2`): `(3−1)·26 = 52 ≤ 3^{2+1}·2 = 54`。 -/
theorem uz_hbnd : ((3 : ℤ) - 1) * uz 2 ≤ (3 : ℤ) ^ (2 + 1) * (2 : ℤ) := by norm_num [uz]

/-- ★★§2 が `k = 2` でも `hu` を出す: `3^m ≤ uz m`(`m ≤ 2`)。 -/
theorem uz_hu : ∀ m, m ≤ 2 → (3 : ℤ) ^ m ≤ uz m :=
  pow_le_of_strictMono_of_dvd (by norm_num) uz_hu0 uz_hult uz_hudvd

end Zeta81

/-! ## §8 点 3 —— 跳びの**伝播**(`p` 冪の得)

★★配られた見立て「`hjump` は `CyclicJumpNorm` の桁展開から出る」は**外れ**である
(§9 に理由)。桁展開は同じ `σ` について `‖σπ−π‖ ≤ ‖π‖^{c+1}` を
`‖σz−z‖ ≤ ‖π‖^c‖z‖` に**言い換える**だけで、`σ` から `σ^p` へは渡らない。

★★渡すのに要るのは**二項係数の可除性 `p ∣ C(p,r)`**(`1 ≤ r ≤ p−1`)である:

```
σ^p − 1 = (1 + D)^p − 1 = Σ_{r=1}^{p} C(p,r) D^r      (D z = σ z − z)
‖D^r π‖ ≤ c^{r−1}·b                                   (縮小率 c を r−1 回)
‖C(p,r)‖ ≤ ‖p‖  (1 ≤ r ≤ p−1)                        ★ここだけが「得」の源
```

⇒ `‖σ^p π − π‖ ≤ max(‖p‖·b, c^{p−1}·b)`。 -/

section BinomialGain

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- 環準同型の加法部分を `ℤ`-線形自己準同型として見る。 -/
def endOf (σ : M →+* M) : Module.End ℤ M := σ.toAddMonoidHom.toIntLinearMap

omit [IsUltrametricDist M] in
@[simp] theorem endOf_apply (σ : M →+* M) (z : M) : endOf σ z = σ z := rfl

omit [IsUltrametricDist M] in
theorem endOf_pow_apply (σ : M →+* M) : ∀ (n : ℕ) (z : M), (endOf σ ^ n) z = (⇑σ)^[n] z := by
  intro n
  induction n with
  | zero => intro z; simp
  | succ m ih =>
      intro z
      rw [pow_succ', Function.iterate_succ_apply', Module.End.mul_apply, ih]
      simp

/-- ★★★**抽象核** —— `σ` の `p` 回反復は `p` 冪の分だけ得をする。

★分岐・付値・Galois の語彙が 1 語も出ない。入力は
`σ : M →+* M`、縮小率 `c ≤ 1`(`‖σz−z‖ ≤ c‖z‖`)、素元の跳び `b`(`‖σπ−π‖ ≤ b`)だけ。 -/
theorem norm_iterate_prime_sub_le {p : ℕ} (hp : p.Prime) {σ : M →+* M} {π : M} {b c : ℝ}
    (hb : 0 ≤ b) (hc0 : 0 ≤ c) (hc1 : c ≤ 1)
    (hstep : ∀ z : M, ‖σ z - z‖ ≤ c * ‖z‖)
    (hw : ‖σ π - π‖ ≤ b) :
    ‖(⇑σ)^[p] π - π‖ ≤ max (‖(p : M)‖ * b) (c ^ (p - 1) * b) := by
  set D : Module.End ℤ M := endOf σ - 1 with hDdef
  have hD : ∀ z : M, D z = σ z - z := by
    intro z; simp [hDdef]
  -- `‖D^{r+1} π‖ ≤ c^r · b`
  have hDpow : ∀ r : ℕ, ‖(D ^ (r + 1)) π‖ ≤ c ^ r * b := by
    intro r
    induction r with
    | zero => simpa [hD] using hw
    | succ n ih =>
        have hstepD : (D ^ (n + 2)) π = D ((D ^ (n + 1)) π) := by
          rw [pow_succ', Module.End.mul_apply]
        rw [hstepD, hD]
        calc ‖σ ((D ^ (n + 1)) π) - (D ^ (n + 1)) π‖ ≤ c * ‖(D ^ (n + 1)) π‖ := hstep _
          _ ≤ c * (c ^ n * b) := by gcongr
          _ = c ^ (n + 1) * b := by ring
  -- 二項展開
  have hsplit : endOf σ = 1 + D := by rw [hDdef]; abel
  have hexp : (⇑σ)^[p] π - π = ∑ m ∈ Finset.range p, (p.choose m : ℕ) • ((D ^ (p - m)) π) := by
    have h1 : (⇑σ)^[p] π = ((endOf σ) ^ p) π := (endOf_pow_apply σ p π).symm
    rw [h1, hsplit, (Commute.one_left D).add_pow p, LinearMap.sum_apply,
      Finset.sum_range_succ]
    simp only [one_pow, one_mul, Module.End.mul_apply, Module.End.natCast_apply, map_nsmul]
    simp
  rw [hexp]
  refine GainedTowerStep.norm_sum_le_of_forall_le ?_ ?_
  · exact le_max_of_le_right (by positivity)
  intro m hm
  have hmp : m < p := Finset.mem_range.mp hm
  have hr : p - m = (p - m - 1) + 1 := by omega
  have hnorm : ‖(p.choose m : ℕ) • ((D ^ (p - m)) π)‖
      = ‖((p.choose m : ℕ) : M)‖ * ‖(D ^ (p - m)) π‖ := by
    rw [nsmul_eq_mul, norm_mul]
  rw [hnorm, hr]
  rcases Nat.eq_zero_or_pos m with hm0 | hm0
  · -- `m = 0`: 二項係数は 1、`D^p π` の分だけ得をする
    subst hm0
    have hc : ‖((p.choose 0 : ℕ) : M)‖ = 1 := by simp
    rw [hc, one_mul]
    refine le_max_of_le_right ?_
    have := hDpow (p - 0 - 1)
    calc ‖(D ^ (p - 0 - 1 + 1)) π‖ ≤ c ^ (p - 0 - 1) * b := this
      _ = c ^ (p - 1) * b := by norm_num
  · -- `1 ≤ m < p`: `p ∣ C(p,m)`
    refine le_max_of_le_left ?_
    obtain ⟨t, ht⟩ := hp.dvd_choose_self (by omega) hmp
    have hcp : ‖((p.choose m : ℕ) : M)‖ ≤ ‖(p : M)‖ := by
      rw [ht]
      push_cast
      rw [norm_mul]
      calc ‖(p : M)‖ * ‖(t : M)‖ ≤ ‖(p : M)‖ * 1 :=
            mul_le_mul_of_nonneg_left (IsUltrametricDist.norm_natCast_le_one M t) (norm_nonneg _)
        _ = ‖(p : M)‖ := mul_one _
    have hDb : ‖(D ^ (p - m - 1 + 1)) π‖ ≤ b := by
      refine le_trans (hDpow (p - m - 1)) ?_
      calc c ^ (p - m - 1) * b ≤ 1 * b := by
            refine mul_le_mul_of_nonneg_right (pow_le_one₀ hc0 hc1) hb
        _ = b := one_mul b
    exact mul_le_mul hcp hDb (norm_nonneg _) (norm_nonneg _)

/-- ★★★★同じことを塔の指数(`‖π‖` の `ℤ` 冪)で書いた形。

`‖σπ−π‖ ≤ ‖π‖^{t+1}` と縮小率 `‖π‖^t`、`‖p‖ = ‖π‖^E` から

  `‖σ^p π − π‖ ≤ ‖π‖^{min(E + t + 1, p·t + 1)}`

★すなわち跳びは `u ↦ min(E + u, p·u)` 以上に伸びる。 -/
theorem norm_iterate_prime_sub_le_zpow {p : ℕ} (hp : p.Prime) {σ : M →+* M} {π : M} {t E : ℤ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (ht0 : 0 ≤ t)
    (hpn : ‖(p : M)‖ = ‖π‖ ^ E)
    (hstep : ∀ z : M, ‖σ z - z‖ ≤ ‖π‖ ^ t * ‖z‖)
    (hw : ‖σ π - π‖ ≤ ‖π‖ ^ (t + 1)) :
    ‖(⇑σ)^[p] π - π‖ ≤ ‖π‖ ^ (min (E + t + 1) ((p : ℤ) * t + 1)) := by
  have hne : ‖π‖ ≠ 0 := ne_of_gt hπ0
  have hc1 : ‖π‖ ^ t ≤ 1 := by
    have := zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) ht0
    simpa using this
  have hmain := norm_iterate_prime_sub_le (M := M) hp (π := π) (b := ‖π‖ ^ (t + 1))
    (c := ‖π‖ ^ t) (by positivity) (by positivity) hc1 hstep hw
  refine le_trans hmain (max_le ?_ ?_)
  · have hA : ‖(p : M)‖ * ‖π‖ ^ (t + 1) = ‖π‖ ^ (E + t + 1) := by
      rw [hpn, ← zpow_add₀ hne]; congr 1; ring
    rw [hA]
    exact zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) (min_le_left _ _)
  · have hcast : ((p - 1 : ℕ) : ℤ) = (p : ℤ) - 1 := by
      have := hp.two_le; omega
    have hB : (‖π‖ ^ t) ^ (p - 1) * ‖π‖ ^ (t + 1) = ‖π‖ ^ ((p : ℤ) * t + 1) := by
      rw [← zpow_natCast (‖π‖ ^ t) (p - 1), ← zpow_mul, ← zpow_add₀ hne, hcast]
      congr 1; ring
    rw [hB]
    exact zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) (min_le_right _ _)

omit [IsUltrametricDist M] in
/-- 層 `j` の `p` 回反復は層 `j+1` である(`τ^{p^j}` を `p` 回で `τ^{p^{j+1}}`)。 -/
theorem iterate_layer_succ {τ : M → M} {s : ℕ → (M →+* M)} {p j : ℕ}
    (hj : ∀ z, s j z = τ^[p ^ j] z) (hj1 : ∀ z, s (j + 1) z = τ^[p ^ (j + 1)] z) (z : M) :
    (⇑(s j))^[p] z = s (j + 1) z := by
  have hfun : (⇑(s j) : M → M) = τ^[p ^ j] := funext hj
  rw [hfun, ← Function.iterate_mul, hj1, pow_succ]

/-- ★★★★★点 3 の 1 段 —— 層 `j` の跳び `u j` と縮小率から**層 `j+1` の跳び**が出る。

★これが「跳びの伝播」の核である。★★配られた見立て(桁展開だけで出る)は外れで、
使うのは `p ∣ C(p,r)`(`norm_iterate_prime_sub_le`)である。 -/
theorem jump_succ_of_jump_of_step {p : ℕ} (hp : p.Prime) {π : M} {t E : ℤ}
    {τ : M → M} {s : ℕ → (M →+* M)} {j : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (ht0 : 0 ≤ t)
    (hpn : ‖(p : M)‖ = ‖π‖ ^ E)
    (hj : ∀ z, s j z = τ^[p ^ j] z) (hj1 : ∀ z, s (j + 1) z = τ^[p ^ (j + 1)] z)
    (hstep : ∀ z : M, ‖s j z - z‖ ≤ ‖π‖ ^ t * ‖z‖)
    (hw : ‖s j π - π‖ ≤ ‖π‖ ^ (t + 1)) :
    ‖s (j + 1) π - π‖ ≤ ‖π‖ ^ (min (E + t + 1) ((p : ℤ) * t + 1)) := by
  rw [← iterate_layer_succ hj hj1 π]
  exact norm_iterate_prime_sub_le_zpow hp hπ0 hπ1 ht0 hpn hstep hw

end BinomialGain

/-! ## §9 ★★点 3 の鎖は点 4 を**含む** —— 底の跳び 1 本から `hu` が出る -/

section ChainToHu

/-- ★★★**抽象核**(ℤ だけ) —— 伝播の鎖 `min(E + u j, p·u j) ≤ u (j+1)` と
予算 `(p−1)·u j ≤ E` から `p^m ≤ u m` が出る。

★これは §8 の `jump_succ_of_jump_of_step` が出す不等式そのものであり、
★★**底の跳び `u 0 ≥ 1` だけから `hu`(Hasse–Arf)が導ける**ことを意味する。
すなわち点 3 が閉じれば点 4 は**要らなくなる**。 -/
theorem pow_le_of_chain {p E : ℤ} (hp : 2 ≤ p) {u : ℕ → ℤ} {k : ℕ} (h0 : 1 ≤ u 0)
    (hchain : ∀ j, j < k → min (E + u j) (p * u j) ≤ u (j + 1))
    (hbudget : ∀ j, j < k → (p - 1) * u j ≤ E) : ∀ m, m ≤ k → p ^ m ≤ u m := by
  intro m
  induction m with
  | zero => intro _; simpa using h0
  | succ n ih =>
      intro hn
      have hnk : n < k := by omega
      have hun : p ^ n ≤ u n := ih (by omega)
      have hpos : (1 : ℤ) ≤ p ^ n := one_le_pow₀ (by linarith)
      have hmin : min (E + u n) (p * u n) = p * u n := by
        have hb := hbudget n hnk
        exact min_eq_right (by linarith)
      have := hchain n hnk
      rw [hmin] at this
      calc (p : ℤ) ^ (n + 1) = p * p ^ n := by ring
        _ ≤ p * u n := by nlinarith
        _ ≤ u (n + 1) := this

end ChainToHu

/-! ## §10 鎖と実データの突き合わせ —— 鎖は**正しいが緩い** -/

namespace ChainCheck

/-- `ℚ₃(ζ₂₇)`: `E = e_M = 18`、`u 0 = 2` ⇒ 鎖の下界は `min(18+2, 3·2) = 6`。
★実データは `u 1 = 8` なので鎖は**成り立つが等号ではない**(`6 ≤ 8`)。 -/
theorem zeta27_chain : min ((18 : ℤ) + 2) (3 * 2) = 6 ∧ (6 : ℤ) ≤ 8 := by norm_num

/-- `ℚ₃(ζ₈₁)`(`k = 1`, `K = ℚ₃(ζ₉)`): `E = 54`、`u 0 = 8` ⇒ 下界 `24 ≤ 26`。 -/
theorem zeta81_chain : min ((54 : ℤ) + 8) (3 * 8) = 24 ∧ (24 : ℤ) ≤ 26 := by norm_num

/-- ★予算 `(p−1)·u j ≤ E` は実データで成り立つ(`4 ≤ 18`, `16 ≤ 54`)。
⇒ §9 が使えて `hu` が出る。 -/
theorem budgets : ((3 : ℤ) - 1) * 2 ≤ 18 ∧ ((3 : ℤ) - 1) * 8 ≤ 54 := by norm_num

end ChainCheck

/-! ## §11 ★★閉じていないもの と その理由(正確に)

### (A) 点 1 —— 体・ノルム込みの模型は**未構成**

★★**配られた見立て「既存の数値層(`Zeta27.*` が 4 ファイルに散在)を集めれば作れる」は外れ**である。

測定(`grep -n "Zeta27" lean/ABC3/Found/PGC/*.lean`、4 ファイル 40 行):
`DeepDescentPairDirect` / `DeepDescentRepair` / `EquivariantProjectionDescent` /
`JumpDefectTradeoff` / `WildDescentMultiStep` の `Zeta27*` は
★**すべて `ℝ` の数値主張**(`axDecay` / `rpow` / 予算の不等式)であり、
`NormedField` のインスタンスも `Padic` も `IsCyclotomicExtension` も 1 つも現れない。
⇒ 集めても `ℝ` の算術が増えるだけで、`M`(体)も `‖·‖` も出てこない。

必要なのは `ℚ₃` の 18 次全分岐拡大の**構成**と、そこでの
`(minpoly _ π).natDegree`、`Algebra.adjoin _ {π} = ⊤`、★`‖σπ − π‖` の**等式**である。
在庫の測定(`grep -in "ramificationGroup\|IsTotallyRamified\|	LocalField" .cache/mathlib-index.txt`):

* 高次分岐群(`G_i`)は mathlib に**無い**(`RingTheory/Valuation/RamificationGroup.lean` は
  分解群と惰性群だけ)。
* `IsTotallyRamified` は**無い**(0 件)。
* `differentIdeal` と `Ideal.ramificationIdx` は**在る**が、`ℚ_p(ζ_{p^n})` の
  分岐の跳びを与える補題は見当たらない。

★★代わりに本ファイルは **ℤ/ℕ 側の模型を 3 組**閉じた(§3・§7):
`(p,k,e,u) = (3,1,2,(2,8))` / `(3,1,6,(8,26))` / ★`(3,2,2,(2,8,26))`。
いずれも**差積の 3 通りの計算**(表・閉じた形・塔の公式)が一致する。

### (B) `hvalj`(値群)は落ちない —— これは**剰余体の情報**である

`hvalj` は `Γ_{E_j} ⊆ p^{k+1−j}·ℤ`、すなわち `e(M/E_j) = [M : E_j]`(全分岐 ＝ `f = 1`)そのもの。

★`hdegj` / `htopj` / `hfixj` からは**出ない**: 同じ次数・同じ原始元・同じ固定体を持つ
**不分岐**拡大では `Γ_{E_j} = Γ_M = ℤ` で `hvalj` は偽になる。
⇒ 群論と超距離だけでは決まらない(反例が同じ仮説を満たす)。

落とすには `d_j`(値群の指数)について 2 本要る:

1. `d_j ≤ [M : E_j]` —— 超距離の 1 次独立(`1,π,…,π^{d_j−1}`)で出る。★1 層で済む。
2. `d_j / d_{j+1} ≤ [E_{j+1} : E_j]` —— ★**2 層をまたぐ**(#59 の危険帯)。

`d_0 = p^{k+1}`(底が全分岐)と 2 を `j` 回つなぐと `d_j ≥ p^{k+1−j}` が出て、
1 と合わせて `d_j = p^{k+1−j}` すなわち `hvalj` になる。★次の波はここを取るとよい。

### (C) 点 3 の残り —— 1 段は閉じたが `∀ j` の配線が (B) に依存する

`jump_succ_of_jump_of_step`(§8)は**層 `j` の縮小率 `hstep` を入力に取る**。
その `hstep` は各層で `exists_expansion_of_adjoin`(桁展開)から作るが、
桁展開の入力が `hvalj` である。⇒ (B) が閉じるまで `∀ j` には回らない。

★★なお §9 の `pow_le_of_chain` が示す通り、点 3 が閉じれば
**点 4(`hu`)は自動的に落ちる**(底の跳び `u 0 ≥ 1` だけになる)。

### (D) 配られた見立ての当否(測定つき)

| 見立て | 判定 |
|---|---|
| 1. 既存の `Zeta27.*` を集めれば体・ノルム込みの模型が作れる | ★**外れ**(全部 `ℝ` の数値主張) |
| 2. #314 の型の族の構え方をそのまま使えば #59/#69 に当たらない | ★**当たり**(§5・§6 は一度も当たっていない) |
| 3. `hjump` は `CyclicJumpNorm` の桁展開から出る | ★**外れ**(桁展開は同じ `σ` を言い換えるだけ。`σ → σ^p` には `p ∣ C(p,r)` が要る) |
| 4. `hu` は Hasse–Arf から | ★**当たり**(§2 で落ちた) |
-/

/-! ## §12 `.src`(原典の対応箇所) -/

def pow_le_of_strictMono_of_dvd.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_hasseArf.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_iterate_prime_sub_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def pow_le_of_chain.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §13 使っている公理の一覧 -/

#print axioms add_le_of_lt_of_dvd
#print axioms pow_le_of_gaps
#print axioms pow_le_of_gaps_le
#print axioms pow_le_of_strictMono_of_dvd
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_hasseArf
#print axioms Zeta27.hOrd_eq
#print axioms Zeta27.different_from_table
#print axioms Zeta27.different_classical
#print axioms Zeta27.different_tower
#print axioms Zeta27.uz_hu
#print axioms Zeta27.uz_one_forced
#print axioms orderOf_pow_pow
#print axioms card_zpowers_pow_pow
#print axioms GaloisTower.algEquiv_pow_apply
#print axioms GaloisTower.fix_twr
#print axioms GaloisTower.finrank_twr
#print axioms GaloisTower.top_twr
#print axioms GaloisTower.deg_twr
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic
#print axioms Zeta81.different_from_table
#print axioms Zeta81.different_tower9
#print axioms Zeta81.different_tower27
#print axioms Zeta81.u_one_forced9
#print axioms Zeta81.uz_hu
#print axioms endOf_pow_apply
#print axioms norm_iterate_prime_sub_le
#print axioms norm_iterate_prime_sub_le_zpow
#print axioms iterate_layer_succ
#print axioms jump_succ_of_jump_of_step
#print axioms pow_le_of_chain
#print axioms ChainCheck.zeta27_chain

end GainedTowerModel

end ABC3.Found.PGC
