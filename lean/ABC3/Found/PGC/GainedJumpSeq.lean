import ABC3.Found.PGC.GainedBridgeSupply
import ABC3.Found.PGC.RamificationJumpBound

/-!
# [pGC] ★★★★★★`GainedBridgeSupply` の `ℤ` 側の仮説を落とす —— 跳びの列は作れる

`Found/PGC/GainedBridgeSupply.lean` の
`exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps` は仮説を 13 本残していた。
持ち場はそれを 2 種類に分けていた:

* (A) `ℤ` 側 —— `hstepZ`(`p·t_{j+1} ≤ t_{j+2}`)と `hlayerZ`(`(p−1)t_{j+1} ≤ p^{j+1}e`、∀ j)
* (B) 入力の言い換え —— `hval` / `hbreak` / `hti` / `hnormp` / `heM` / `hdeg` / `htop`

★★本ファイルは (A) を**丸ごと**落とす。加えて `ht1` と `hbreak` と、`hjump` の
`j > k` の分も落ちる。

## ★★★核心 —— 結論に `t` は出てこない

`exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps` の結論は

```
∃ y ∈ F,  ‖x − y‖ ≤ (∏_{j∈[1,k+1]} axDecay p j) · ‖τ x − x‖
```

で、★**跳びの列 `t` を含まない**。したがって `t` は `∀`(仮説の側)にあり、
「都合のよい `t` を 1 本作って代入する」ことが許される。
★これに気づくと (A) は数学ではなく**構成の問題**になる。

作る列は「上から `p` で割る鎖」である。頂点 `t_{k+1} := u_k`(ここは `i` の定義なので動かせない)、
そこから下へ

```
t_{j} := min( u_{j−1} , ⌊ t_{j+1} / p ⌋ )          (j ≤ k)
t_{j+1} := p^{j−k} · u_k                            (j > k)
```

とすると、`min` の第 2 引数から `p·t_{j+1} ≤ t_{j+2}`(`hstepZ`)が**定義から**出て、
`min` の第 1 引数から `t_{j+1} ≤ u_j`(`hjump` に代入できる)が出る。
さらに `p^{k−j}·t_{j+1} ≤ u_k` なので `hlayerZ` は**頂点の 1 本**
`(p−1)u_k ≤ p^{k+1}e` から全層ぶん出る。
★`1 ≤ t_{j+1}` が保たれるための条件がちょうど `p^m ≤ u_m`(`m ≤ k`)である。

## ★★★落ちた仮説と、代わりに入ったもの

| 落ちた仮説 | どこから出るか |
|---|---|
| `ht1 : 1 ≤ t_{j+1}` | `hu`(`p^m ≤ u_m`) |
| ★★`hstepZ : p·t_{j+1} ≤ t_{j+2}` | ★構成(`min` の第 2 引数)。仮定ではなくなった |
| ★★`hlayerZ : (p−1)t_{j+1} ≤ p^{j+1}e`(∀ j) | ★頂点 1 本 `hbnd : (p−1)u_k ≤ p^{k+1}e` |
| `hjump` の `j > k` の分 | ★`norm_iterate_pow_sub_self_le_lt`(帰納を書き直した) |
| ★`hbnd` 自身(§9) | ★`heM` ＋ 最小多項式の判別式(`RamificationJumpBound`、木に在った) |
| ★`hbreak`(1 つの `σ` の跳び) | ★`hconj`(全共役の跳び)＋ `σπ ≠ π` |

★★★入るのは `hu : ∀ m ≤ k, p^m ≤ u_m` **1 本だけ**である
(§9 ではさらに `M = F(π)` が `p` 次で分離・分解する、という構造の 4 本)。

## ★★`hu` は `hstepZ` より真に弱い(★これが本ファイルの主張の要点)

`hstepZ` を実際の跳び `t = u` に当てると `p·u_j ≤ u_{j+1}` で、`u_0 ≥ 1` と併せて
`p^m ≤ u_m` が出る。逆は出ない。★つまり
**`hu` ⟸ `ht1` ＋ `hstepZ`、`hbnd` ⟸ `hlayerZ`(`j = k`)** であり、
本ファイルの出口は `GainedBridgeSupply` の出口より**真に強い**。

★数値で見ると: `p = 3`, `k = 1`, 実際の下付き跳び `u = (u₀,u₁) = (2,5)` のとき
`hstepZ` を `t = u` に当てると `3·2 = 6 ≤ 5` で**偽**だが、
本ファイルの構成は `t = (min(2, ⌊5/3⌋), 5) = (1,5)` を作り `3·1 = 3 ≤ 5` で通る
(`Numeric.seq_two_five` で機械検算した)。
★`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` の `u = (2,8)` では構成は `t = (2,8)` をそのまま返す
(`Numeric.seq_zeta27`)—— 古典的な場合は列を壊さない。

## ★★★測定 —— 木に在ったもの / 無かったもの(コマンドつき)

| 部品 | 在庫 | 測定 |
|---|---|---|
| 頂点の跳びの上界 `(p−1)i ≤ v_M(p)` | ★★**木に在った** | `grep -n "sub_one_mul_le_of_norm_natCast_eq_pow" lean/ABC3/Found/PGC/RamificationJumpBound.lean`(334 行)＋`norm_natCast_le_pow_of_splits`(316 行) |
| `σ^{p^j} ∈ G_1 ⟹ σ^{p^{j+k}} ∈ G_{1+k}` | 木に在るが**使えない** | `RamificationJumpDivisibility.lean:344`。`u_m ≥ 1+m` しか出ず `p^m ≤ u_m` に届かない |
| Hasse–Arf(`p^{j+1} ∣ u_{j+1}−u_j`) | ★木は**別の設定**で持つ | `HasseArf.lean` / `HasseArfStrongInduction.lean` は `herbrandPhi` が自然数値であることを `[CommRing A] [CommRing B]` の枠で述べる。本ファイルのノルム言語(`F M : Field/NormedField`)に橋が無い |
| `Int.le_ediv_of_mul_le` / `Int.ediv_mul_le` | mathlib に在った | `#check @Int.le_ediv_of_mul_le`(`0 < c → a * c ≤ b → a ≤ b / c`) |
| 上付き/下付き分岐群・Herbrand | mathlib に**無い** | 木が自前で持つ |

## ★★★閉じていないもの(★正確に)

1. ★★**`hu : ∀ m ≤ k, p^m ≤ u_m` は仮説のまま**である。これは
   「下付き跳びが `u_{j+1} ≥ u_j + p^{j+1}` で伸びる」= ★**Hasse–Arf の中身**で、
   木の `HasseArf*.lean` は同じ内容を `herbrandPhi` の自然数性として持つが、
   ★**そちらは `[CommRing A] [CommRing B]` の設定**で、本ファイルの
   `[Field F] [NormedField M]` に繋ぐ橋が無い(測定した。上表)。
   ★この 1 本が pGC の Ax–Sen–Tate に残る `ℤ` 側の穴である。
2. `hnormp` / `heM` / `hval` / `hdeg` / `htop` は (B) の入力で、本ファイルは触っていない。
   ★これらは「`M/F_0` が `p^{k+1}` 次の全分岐塔で `π` が `M` の素元」という**設定**であり、
   `Interface/PGC/LocalFieldData.lean` には対応する `structure` が無い
   (`grep -n "^structure" lean/ABC3/Interface/PGC/LocalFieldData.lean` は
   `ResidueCardinality` / `SubgroupCorrespondence` / `RamificationFiltration` の 3 本のみ)。
   ★持ち場の見立て「`Interface/PGC` の `PAdicLocalField` の構造から出るはず」は**外れ**で、
   `PAdicLocalField` は `Skeleton/PGC/Setup.lean:40` にあり、素元も分岐指数も持っていない。
3. `hjump` は `j < k` に減ったが消えてはいない。これは下付き分岐群の定義そのものである。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. `GainedBridgeSupply` / `GainedDescentBridge` / `RamificationJumpBound` は**読むだけ**で
   1 行も書き換えていない。弱めた版はすべて本ファイルに別名で置いた。
2. ★`norm_iterate_pow_sub_self_le_lt` は `GainedBridge.norm_iterate_pow_sub_self_le` の
   帰納を書き直したもの(`θ' j = if j < k then θ j else 1` の代入では `hstep` の
   `j ≥ k` が消せなかったため)。中身は同じで、仮説を `∀ j < k` にしただけである。
3. 原典(Serre *Corps Locaux* IV)は付値言語。本ファイルはノルム言語
   (`v_M(z) = m` は `‖z‖ = ‖π‖^m`)。`GainedDescentBridge` と同じ約束である。
-/
open Finset

namespace ABC3.Found.PGC

namespace GainedJumpSeq

/-! ## §1 抽象核 -/

/-- 深さ `d` の鎖 `chain p v d = min_{m ≤ d} ⌊v m / p^{d−m}⌋`。 -/
def chain (p : ℤ) (v : ℕ → ℤ) : ℕ → ℤ
  | 0 => v 0
  | (d + 1) => min (v (d + 1)) (chain p v d / p)

theorem chain_zero (p : ℤ) (v : ℕ → ℤ) : chain p v 0 = v 0 := rfl

theorem chain_succ (p : ℤ) (v : ℕ → ℤ) (d : ℕ) :
    chain p v (d + 1) = min (v (d + 1)) (chain p v d / p) := rfl

theorem chain_le_self (p : ℤ) (v : ℕ → ℤ) : ∀ d, chain p v d ≤ v d
  | 0 => le_refl _
  | (_ + 1) => min_le_left _ _

theorem mul_chain_succ_le {p : ℤ} (hp : 0 < p) (v : ℕ → ℤ) (d : ℕ) :
    p * chain p v (d + 1) ≤ chain p v d := by
  have h1 : chain p v (d + 1) * p ≤ chain p v d / p * p :=
    mul_le_mul_of_nonneg_right (min_le_right _ _) hp.le
  have h2 : chain p v d / p * p ≤ chain p v d := Int.ediv_mul_le _ (ne_of_gt hp)
  linarith

theorem pow_mul_chain_le {p : ℤ} (hp : 0 < p) (v : ℕ → ℤ) :
    ∀ d, p ^ d * chain p v d ≤ v 0
  | 0 => by simp [chain_zero]
  | (d + 1) => by
      have ih := pow_mul_chain_le hp v d
      have h := mul_chain_succ_le hp v d
      have hpd : (0 : ℤ) ≤ p ^ d := le_of_lt (pow_pos hp d)
      calc p ^ (d + 1) * chain p v (d + 1) = p ^ d * (p * chain p v (d + 1)) := by ring
        _ ≤ p ^ d * chain p v d := mul_le_mul_of_nonneg_left h hpd
        _ ≤ v 0 := ih

theorem le_chain {p : ℤ} (hp : 0 < p) {v : ℕ → ℤ} :
    ∀ (d : ℕ) (c : ℤ), (∀ m, m ≤ d → p ^ (d - m) * c ≤ v m) → c ≤ chain p v d
  | 0, c, h => by
      have := h 0 (le_refl 0)
      simpa [chain_zero] using this
  | (d + 1), c, h => by
      have ih : p * c ≤ chain p v d := by
        refine le_chain hp d (p * c) ?_
        intro m hm
        have hh := h m (Nat.le_succ_of_le hm)
        have he : d + 1 - m = (d - m) + 1 := by omega
        rw [he] at hh
        calc p ^ (d - m) * (p * c) = p ^ (d - m + 1) * c := by ring
          _ ≤ v m := hh
      rw [chain_succ]
      refine le_min ?_ (Int.le_ediv_of_mul_le hp ?_)
      · have := h (d + 1) (le_refl _)
        simpa using this
      · linarith


/-! ## §2 跳びの列の構成 -/

/-- 実際の跳び `u` から作る列。`j ≤ k` では上から `p` で割りながら `u j` と `min` を取り、
`j > k` では `p^{j−k}·u k` に延ばす。 -/
def seq (p : ℤ) (u : ℕ → ℤ) (k j : ℕ) : ℤ :=
  if j ≤ k then chain p (fun d => u (k - d)) (k - j) else p ^ (j - k) * u k

theorem seq_top (p : ℤ) (u : ℕ → ℤ) (k : ℕ) : seq p u k k = u k := by
  simp [seq, chain_zero]

theorem seq_le_of_le (p : ℤ) (u : ℕ → ℤ) {k j : ℕ} (hj : j ≤ k) : seq p u k j ≤ u j := by
  have h := chain_le_self p (fun d => u (k - d)) (k - j)
  simp only [seq, if_pos hj]
  simpa [Nat.sub_sub_self hj] using h

theorem pow_mul_seq_le (p : ℤ) (u : ℕ → ℤ) {k j : ℕ} (hp : 0 < p) (hj : j ≤ k) :
    p ^ (k - j) * seq p u k j ≤ u k := by
  have h := pow_mul_chain_le hp (fun d => u (k - d)) (k - j)
  simp only [seq, if_pos hj]
  simpa using h

theorem one_le_seq {p : ℤ} (hp : 2 ≤ p) {u : ℕ → ℤ} {k : ℕ}
    (hu : ∀ m, m ≤ k → p ^ m ≤ u m) (j : ℕ) : 1 ≤ seq p u k j := by
  have hp0 : (0 : ℤ) < p := by linarith
  have hp1 : (1 : ℤ) ≤ p := by linarith
  by_cases hj : j ≤ k
  · simp only [seq, if_pos hj]
    refine le_chain hp0 (k - j) 1 ?_
    intro m hm
    have hmk : k - m ≤ k := Nat.sub_le _ _
    have h1 : p ^ (k - m) ≤ u (k - m) := hu _ hmk
    have h2 : p ^ (k - j - m) ≤ p ^ (k - m) := pow_le_pow_right₀ hp1 (by omega)
    simpa using le_trans h2 h1
  · simp only [seq, if_neg hj]
    have h1 : (1 : ℤ) ≤ p ^ (j - k) := one_le_pow₀ hp1
    have h2 : (1 : ℤ) ≤ u k := le_trans (one_le_pow₀ hp1) (hu k (le_refl k))
    nlinarith

theorem mul_seq_le_seq_succ {p : ℤ} (hp : 2 ≤ p) (u : ℕ → ℤ) (k j : ℕ) :
    p * seq p u k j ≤ seq p u k (j + 1) := by
  have hp0 : (0 : ℤ) < p := by linarith
  rcases lt_trichotomy j k with h | h | h
  · have hj : j ≤ k := le_of_lt h
    have hj1 : j + 1 ≤ k := h
    simp only [seq, if_pos hj, if_pos hj1]
    have he : k - j = (k - (j + 1)) + 1 := by omega
    rw [he]
    exact mul_chain_succ_le hp0 _ _
  · subst h
    have h1 : seq p u j j = u j := seq_top p u j
    have h2 : seq p u j (j + 1) = p ^ (j + 1 - j) * u j := by
      simp only [seq, if_neg (by omega : ¬ j + 1 ≤ j)]
    rw [h1, h2]
    have h3 : j + 1 - j = 1 := by omega
    rw [h3]
    ring_nf
    exact le_refl _
  · have hj : ¬ j ≤ k := by omega
    have hj1 : ¬ j + 1 ≤ k := by omega
    simp only [seq, if_neg hj, if_neg hj1]
    have he : j + 1 - k = (j - k) + 1 := by omega
    rw [he]
    ring_nf
    exact le_refl _

theorem sub_one_mul_seq_le {p : ℤ} (hp : 2 ≤ p) {u : ℕ → ℤ} {e : ℤ} {k : ℕ}
    (hbnd : (p - 1) * u k ≤ p ^ (k + 1) * e) (j : ℕ) :
    (p - 1) * seq p u k j ≤ p ^ (j + 1) * e := by
  have hp0 : (0 : ℤ) < p := by linarith
  by_cases hj : j ≤ k
  · have h := pow_mul_seq_le p u hp0 hj
    have hpos : (0 : ℤ) < p ^ (k - j) := pow_pos hp0 _
    have h1 : p ^ (k - j) * ((p - 1) * seq p u k j) ≤ (p - 1) * u k := by
      have h' := mul_le_mul_of_nonneg_left h (by linarith : (0:ℤ) ≤ p - 1)
      nlinarith [h']
    have h2 : p ^ (k - j) * ((p - 1) * seq p u k j) ≤ p ^ (k - j) * (p ^ (j + 1) * e) := by
      refine le_trans h1 (le_trans hbnd (le_of_eq ?_))
      rw [← mul_assoc, ← pow_add]
      congr 2
      omega
    exact le_of_mul_le_mul_left h2 hpos
  · simp only [seq, if_neg hj]
    have hpos : (0 : ℤ) ≤ p ^ (j - k) := le_of_lt (pow_pos hp0 _)
    have h1 : p ^ (j - k) * ((p - 1) * u k) ≤ p ^ (j - k) * (p ^ (k + 1) * e) :=
      mul_le_mul_of_nonneg_left hbnd hpos
    have h2 : p ^ (j - k) * (p ^ (k + 1) * e) = p ^ (j + 1) * e := by
      rw [← mul_assoc, ← pow_add]
      congr 2
      omega
    calc (p - 1) * (p ^ (j - k) * u k) = p ^ (j - k) * ((p - 1) * u k) := by ring
      _ ≤ p ^ (j - k) * (p ^ (k + 1) * e) := h1
      _ = p ^ (j + 1) * e := h2

/-! ## §3 抽象核の出口 -/

/-- ★★★**抽象核**: 「上から `p` で割る鎖」で、跳びの列 `u` から
`GainedTowerDescent` が要求する `ℤ` 側の条件をすべて満たす列 `t` を作る。

★入力は `p^m ≤ u m`(`m ≤ k`)と**頂点だけ**の上界 `(p−1)·u k ≤ p^{k+1}·e` の 2 本。
★出力の `t` は `t_{k+1} = u_k`(頂点は動かさない)、`t_{j+1} ≤ u_j`、
`1 ≤ t_{j+1}`、`p·t_{j+1} ≤ t_{j+2}`、`(p−1)·t_{j+1} ≤ p^{j+1}·e`。

★分岐・付値・Galois の語彙が 1 語も出ない。 -/
theorem exists_jump_seq {p : ℤ} (hp : 2 ≤ p) {u : ℕ → ℤ} {e : ℤ} (k : ℕ)
    (hu : ∀ m, m ≤ k → p ^ m ≤ u m)
    (hbnd : (p - 1) * u k ≤ p ^ (k + 1) * e) :
    ∃ t : ℕ → ℤ,
      t (k + 1) = u k ∧
      (∀ j, 1 ≤ t (j + 1)) ∧
      (∀ j, j ≤ k → t (j + 1) ≤ u j) ∧
      (∀ j, p * t (j + 1) ≤ t (j + 2)) ∧
      (∀ j, (p - 1) * t (j + 1) ≤ p ^ (j + 1) * e) := by
  refine ⟨fun n => seq p u k (n - 1), ?_, ?_, ?_, ?_, ?_⟩
  · simpa using seq_top p u k
  · intro j; simpa using one_le_seq hp hu j
  · intro j hj; simpa using seq_le_of_le p u hj
  · intro j
    have h := mul_seq_le_seq_succ hp u k j
    have e1 : j + 1 - 1 = j := by omega
    have e2 : j + 2 - 1 = j + 1 := by omega
    simp only [e1, e2]
    exact h
  · intro j; simpa using sub_one_mul_seq_le hp hbnd j

/-! ## §4 抽象核(ノルム) —— `hstep` / `hlayer` を `∀ j < k` に落とす -/

section NormCore

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

/-- ★★`GainedBridge.norm_iterate_pow_sub_self_le` の**両方**の仮説を `∀ j < k` に弱めた形。

★`GainedBridgeSupply.norm_iterate_pow_sub_self_le_of_lt` は `hlayer` だけを弱めており、
`hstep` は `∀ j`(`j ≥ k` も)必要だった。その `j ≥ k` の分は
具体層で `τ^{p^j}`(`j > k`)の縮小率を要求するので、`τ` の位数が `p^{k+1}` であることを
仮定しないと供給できない。★本補題は帰納を書き直してその要求を消す
(帰納が実際に使うのは `j < k` の層だけである)。 -/
theorem norm_iterate_pow_sub_self_le_lt {p : ℕ} (hp : p.Prime) {τ : L →+ L} {θ : ℕ → ℝ} :
    ∀ (k : ℕ), (∀ j, j < k → 0 ≤ θ j) → (∀ j, j < k → θ j ≤ 1) →
      (∀ j, j < k → ∀ z : L, ‖(⇑τ)^[p ^ j] z - z‖ ≤ θ j * ‖z‖) →
      (∀ j, j < k → ‖(p : L)‖ ≤ θ j ^ (p - 1)) →
      ∀ x : L, ‖(⇑τ)^[p ^ k] x - x‖ ≤ (∏ j ∈ Finset.range k, θ j ^ (p - 1)) * ‖τ x - x‖ := by
  intro k
  induction k with
  | zero => intro _ _ _ _ x; simp
  | succ k ih =>
    intro hθ0 hθ1 hstep hlayer x
    have hg : ⇑((GainedBridge.toEnd τ) ^ (p ^ k)) = (⇑τ)^[p ^ k] :=
      AddMonoid.End.coe_pow L (GainedBridge.toEnd τ) (p ^ k)
    have hiter : (⇑τ)^[p ^ (k + 1)] x = (⇑((GainedBridge.toEnd τ) ^ (p ^ k)))^[p] x := by
      rw [hg, ← Function.iterate_mul, ← pow_succ]
    have hD : ∀ z : L, ‖((GainedBridge.toEnd τ) ^ (p ^ k)) z - z‖ ≤ θ k * ‖z‖ := by
      intro z; rw [hg]; exact hstep k (Nat.lt_succ_self k) z
    have h1 := GainedBridge.norm_iterate_prime_sub_self_le hp
      (τ := (GainedBridge.toEnd τ) ^ (p ^ k))
      (hθ0 k (Nat.lt_succ_self k)) (hθ1 k (Nat.lt_succ_self k)) hD x
    rw [max_eq_left (hlayer k (Nat.lt_succ_self k))] at h1
    rw [hiter]
    have h2 : ‖((GainedBridge.toEnd τ) ^ (p ^ k)) x - x‖
        ≤ (∏ j ∈ Finset.range k, θ j ^ (p - 1)) * ‖τ x - x‖ := by
      rw [hg]
      exact ih (fun j hj => hθ0 j (Nat.lt_succ_of_lt hj))
        (fun j hj => hθ1 j (Nat.lt_succ_of_lt hj))
        (fun j hj => hstep j (Nat.lt_succ_of_lt hj))
        (fun j hj => hlayer j (Nat.lt_succ_of_lt hj)) x
    have hpow0 : (0 : ℝ) ≤ θ k ^ (p - 1) := pow_nonneg (hθ0 k (Nat.lt_succ_self k)) _
    calc ‖(⇑((GainedBridge.toEnd τ) ^ (p ^ k)))^[p] x - x‖
        ≤ θ k ^ (p - 1) * ‖((GainedBridge.toEnd τ) ^ (p ^ k)) x - x‖ := h1
      _ ≤ θ k ^ (p - 1) * ((∏ j ∈ Finset.range k, θ j ^ (p - 1)) * ‖τ x - x‖) :=
          mul_le_mul_of_nonneg_left h2 hpow0
      _ = (∏ j ∈ Finset.range (k + 1), θ j ^ (p - 1)) * ‖τ x - x‖ := by
          rw [Finset.prod_range_succ]; ring

end NormCore

/-! ## §5 体の層への橋 —— `hstep` を `∀ j < k` にした写し -/

section Bridge

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- `GainedBridgeSupply.norm_sub_digit_zero_le_gain_of_lt` の `hstep` も `∀ j < k` にした形。 -/
theorem norm_sub_digit_zero_le_gain_lt {p i k : ℕ} (hp : p.Prime) {π : M} {a : ℕ → F}
    (τ : M →+ M) (σ : M →ₐ[F] M) {θ : ℕ → ℝ}
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z)
    (hθ0 : ∀ j, j < k → 0 ≤ θ j) (hθ1 : ∀ j, j < k → θ j ≤ 1)
    (hstep : ∀ j, j < k → ∀ z : M, ‖(⇑τ)^[p ^ j] z - z‖ ≤ θ j * ‖z‖)
    (hlayer : ∀ j, j < k → ‖(p : M)‖ ≤ θ j ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i)
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m)) :
    ‖digitSum p π a - algebraMap F M (a 0)‖
      ≤ ‖π‖ ^ (-(i : ℤ)) * (∏ j ∈ Finset.range k, θ j ^ (p - 1))
          * ‖τ (digitSum p π a) - digitSum p π a‖ := by
  have heq := norm_sub_digit_zero_eq_zpow_mul (a := a) σ hp.pos hπ0 hπ1 hi hbreak hchar hval
  have hchain := norm_iterate_pow_sub_self_le_lt hp k hθ0 hθ1 hstep hlayer (digitSum p π a)
  rw [heq, hσ (digitSum p π a)]
  have hz0 : (0 : ℝ) ≤ ‖π‖ ^ (-(i : ℤ)) := le_of_lt (zpow_pos hπ0 _)
  calc ‖π‖ ^ (-(i : ℤ)) * ‖(⇑τ)^[p ^ k] (digitSum p π a) - digitSum p π a‖
      ≤ ‖π‖ ^ (-(i : ℤ)) * ((∏ j ∈ Finset.range k, θ j ^ (p - 1))
          * ‖τ (digitSum p π a) - digitSum p π a‖) := mul_le_mul_of_nonneg_left hchain hz0
    _ = _ := by ring

/-- `GainedBridgeSupply.norm_sub_digit_zero_le_gainedLoss_of_lt` の `hstep` も `∀ j < k` にした形。 -/
theorem norm_sub_digit_zero_le_gainedLoss_lt {p i k : ℕ} (hp : p.Prime) {π : M} {a : ℕ → F}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ j, j < k → ∀ z : M, ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hlayer : ∀ j, j < k → ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m)) :
    ‖digitSum p π a - algebraMap F M (a 0)‖
      ≤ ‖π‖ ^ (-(GainedDescent.gainedLoss (p : ℤ) t (k + 1)))
          * ‖τ (digitSum p π a) - digitSum p π a‖ := by
  have h := norm_sub_digit_zero_le_gain_lt (a := a) (θ := fun j => ‖π‖ ^ (t (j + 1)))
    hp τ σ hσ (fun j _ => le_of_lt (zpow_pos hπ0 _))
    (fun j _ => zpow_le_one₀ hπ0 (le_of_lt hπ1) (ht0 j)) hstep hlayer hπ0 hπ1 hi hbreak hchar hval
  refine le_trans h (mul_le_mul_of_nonneg_right ?_ (norm_nonneg _))
  rw [GainedBridge.prod_theta_eq hπ0 hp.one_lt.le t k, ← zpow_add₀ (ne_of_gt hπ0), hti]
  refine zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) ?_
  have := GainedBridge.topDefect_le_gainedLoss (p : ℤ) t k
  omega

/-- `GainedBridgeSupply.exists_norm_sub_algebraMap_le_gainedLoss_of_lt` の `hstep` も
`∀ j < k` にした形。 -/
theorem exists_norm_sub_algebraMap_le_gainedLoss_lt {p i k : ℕ} (hp : p.Prime) {π : M}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ j, j < k → ∀ z : M, ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hlayer : ∀ j, j < k → ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ ‖π‖ ^ (-(GainedDescent.gainedLoss (p : ℤ) t (k + 1))) * ‖τ x - x‖ := by
  have hint : IsIntegral F π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    exact hp.ne_zero hdeg.symm
  obtain ⟨a, rfl⟩ := exists_digitSum_of_adjoin_eq_top hint hdeg htop x
  exact ⟨a 0, norm_sub_digit_zero_le_gainedLoss_lt hp τ σ t hσ ht0 hstep hlayer hπ0 hπ1 hi hti
    hbreak hchar hval⟩

end Bridge

/-! ## §6 出口(素元の正規化) —— `hstep` を `∀ j < k` にした写し -/

section Exit

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- `GainedBridgeSupply.exists_norm_sub_algebraMap_le_prod_axDecay_of_uniformizer` の
`hstep` を `∀ j < k` にした形。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_lt
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ) (he : 0 < e)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ j, j < k → ∀ z : M, ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hstepZ : ∀ j, (p : ℤ) * t (j + 1) ≤ t (j + 2))
    (hlayerZ : ∀ j, ((p : ℤ) - 1) * t (j + 1) ≤ (p : ℤ) ^ (j + 1) * (e : ℤ))
    (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp'.one_lt
  have hpR0 : (0 : ℝ) < (p : ℝ) := lt_trans zero_lt_one hpR
  have hEne : p ^ (k + 1) * e ≠ 0 := Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne'
  have hπpow : ‖π‖ ^ (p ^ (k + 1) * e) = ((p : ℝ))⁻¹ := by rw [← heM]; exact hnormp
  have hplt : ‖(p : M)‖ < 1 := by rw [hnormp, inv_lt_one_iff₀]; exact Or.inr hpR
  have hπ0 : 0 < ‖π‖ := by
    rcases lt_or_eq_of_le (norm_nonneg π) with h | h
    · exact h
    · exfalso
      rw [← h, zero_pow hEne] at hπpow
      exact absurd hπpow.symm (ne_of_gt (inv_pos.mpr hpR0))
  have hπ1 : ‖π‖ < 1 := by
    refine GainedBridgeSupply.lt_one_of_pow_lt_one (E := p ^ (k + 1) * e) ?_
    rw [hπpow, inv_lt_one_iff₀]; exact Or.inr hpR
  have hπE : ‖π‖ = (p : ℝ) ^ (-(1 : ℝ) / ((p ^ (k + 1) * e : ℕ) : ℝ)) :=
    GainedBridgeSupply.eq_rpow_neg_one_div_of_pow_eq_inv hπ0 (le_of_lt hpR0) hEne hπpow
  have hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1 :=
    fun j hj0 hjp => norm_natCast_eq_one_of_lt_prime hp' hplt hj0 hjp
  have hlayer : ∀ j, j < k → ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1) :=
    fun j hj => GainedBridgeSupply.norm_natCast_le_zpow_pow_of_layerZ hp'.one_lt.le hπ0
      (le_of_lt hπ1) heM hlayerZ j (le_of_lt hj)
  obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_gainedLoss_lt hp' τ σ t hσ ht0
    hstep hlayer hπ0 hπ1 hi hti hbreak hchar hval hdeg htop x
  refine ⟨y, le_trans hy (mul_le_mul_of_nonneg_right ?_ (norm_nonneg _))⟩
  have hΛ0 : 0 ≤ GainedDescent.gainedLoss (p : ℤ) t (k + 1) :=
    GainedBridge.gainedLoss_nonneg (by simpa using ht0 0) hstepZ k
  have hJ : ((GainedDescent.gainedLoss (p : ℤ) t (k + 1)).toNat : ℤ)
      = GainedDescent.gainedLoss (p : ℤ) t (k + 1) := Int.toNat_of_nonneg hΛ0
  have hmain := GainedDescent.rpow_le_prod_axDecay (p := p) (e := e) (k := k + 1)
    (J := (GainedDescent.gainedLoss (p : ℤ) t (k + 1)).toNat) (t := t) he ht0 hstepZ hlayerZ
    (le_of_eq hJ)
  have hcast : ((GainedDescent.gainedLoss (p : ℤ) t (k + 1) : ℤ) : ℝ)
      = (((GainedDescent.gainedLoss (p : ℤ) t (k + 1)).toNat : ℕ) : ℝ) := by
    exact_mod_cast hJ.symm
  rw [GainedBridge.zpow_eq_rpow_div hπE, hcast]
  exact hmain

end Exit

/-! ## §7 出口 —— `hstepZ` / `hlayerZ`(∀ j) / `ht1` が消える -/

section Exit2

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- ★★★★★**`GainedBridgeSupply` の出口から `ℤ` 側の仮説 3 本を落とした形**。

落ちたのは `ht1`(`1 ≤ t_{j+1}`)、★`hstepZ`(`p·t_{j+1} ≤ t_{j+2}`、∀ j)、
★`hlayerZ`(`(p−1)t_{j+1} ≤ p^{j+1}e`、∀ j)の 3 本。代わりに入るのは

* `hu : ∀ m ≤ k, p^m ≤ u m` —— ★跳びの列の下からの伸び、
* `hbnd : (p−1)·u k ≤ p^{k+1}·e` —— ★**頂点 1 点だけ**の上界。

★★どちらも元の仮説より**真に弱い**:
`hu` は `ht1`＋`hstepZ`(を `t = u` に当てたもの)の帰結であり、
`hbnd` は `hlayerZ` を `j = k` に制限したものである。
★列 `t` は `exists_jump_seq` が `u` から作る(結論に `t` は出てこないので自由に取れる)。

★さらに `hjump` が `∀ j`(`j > k` も)から `∀ j < k` になった。
`j > k` の分は `τ` の位数が `p^{k+1}` であることを仮定しないと供給できなかったので、
これも仮説の削減である(`norm_iterate_pow_sub_self_le_lt` を書き直したのが効いている)。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M} {u : ℕ → ℤ}
    (τ : M →ₐ[F] M) (he : 0 < e)
    (hu : ∀ m, m ≤ k → (p : ℤ) ^ m ≤ u m)
    (hjump : ∀ j, j < k → ‖(τ ^ p ^ j) π - π‖ ≤ ‖π‖ ^ (u j + 1))
    (hbnd : ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ))
    (hti : (i : ℤ) = u k)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbreak : ‖(τ ^ p ^ k) π - π‖ = ‖π‖ ^ (i + 1))
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have hpZ : (2 : ℤ) ≤ (p : ℤ) := by exact_mod_cast hp'.two_le
  obtain ⟨t, htk, ht1, htu, hstepZ, hlayerZ⟩ := exists_jump_seq hpZ (e := (e : ℤ)) k hu hbnd
  have ht0 : ∀ j, 0 ≤ t (j + 1) := fun j => le_trans zero_le_one (ht1 j)
  have hti' : (i : ℤ) = t (k + 1) := by rw [htk]; exact hti
  have hi : 0 < i := by
    have h1 : (1 : ℤ) ≤ u k := le_trans (one_le_pow₀ (by linarith)) (hu k (le_refl k))
    omega
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp'.one_lt
  have hpR0 : (0 : ℝ) < (p : ℝ) := lt_trans zero_lt_one hpR
  have hEne : p ^ (k + 1) * e ≠ 0 := Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne'
  have hπpow : ‖π‖ ^ (p ^ (k + 1) * e) = ((p : ℝ))⁻¹ := by rw [← heM]; exact hnormp
  have hplt : ‖(p : M)‖ < 1 := by rw [hnormp, inv_lt_one_iff₀]; exact Or.inr hpR
  have hπ0 : 0 < ‖π‖ := by
    rcases lt_or_eq_of_le (norm_nonneg π) with h | h
    · exact h
    · exfalso
      rw [← h, zero_pow hEne] at hπpow
      exact absurd hπpow.symm (ne_of_gt (inv_pos.mpr hpR0))
  have hπ1 : ‖π‖ < 1 := by
    refine GainedBridgeSupply.lt_one_of_pow_lt_one (E := p ^ (k + 1) * e) ?_
    rw [hπpow, inv_lt_one_iff₀]; exact Or.inr hpR
  have hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1 :=
    fun j hj0 hjp => norm_natCast_eq_one_of_lt_prime hp' hplt hj0 hjp
  set τ' : M →+ M := AddMonoidHom.mk' (⇑τ) (fun a b => map_add τ a b) with hτ'
  have hcoe : ⇑τ' = ⇑τ := rfl
  have hpowcoe : ∀ j : ℕ, (⇑τ')^[p ^ j] = ⇑(τ ^ p ^ j) := by
    intro j; rw [hcoe, AlgHom.coe_pow]
  have hstep : ∀ j, j < k → ∀ z : M, ‖(⇑τ')^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖ := by
    intro j hj z
    rw [hpowcoe j]
    refine GainedBridgeSupply.norm_algHom_sub_le_mul_of_break hp'.pos (τ ^ p ^ j) hπ0 hπ1
      (zpow_lt_one₀ hπ0 hπ1 (lt_of_lt_of_le zero_lt_one (ht1 j))) hchar hval ?_ hdeg htop z
    refine le_trans (hjump j hj) ?_
    rw [← zpow_add_one₀ (ne_of_gt hπ0)]
    refine zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) ?_
    have := htu j (le_of_lt hj)
    omega
  have hσ : ∀ z : M, (τ ^ p ^ k) z = (⇑τ')^[p ^ k] z := by intro z; rw [hpowcoe k]
  obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_prod_axDecay_lt
    (F := F) (M := M) (p := p) (e := e) (i := i) (k := k) (π := π) τ' (τ ^ p ^ k) t he hσ ht0
    hstep hstepZ hlayerZ hi hti' hnormp heM hbreak hval hdeg htop x
  exact ⟨y, hy⟩

end Exit2

/-! ## §8 `hbnd` は追加の仮定ではない —— 最小多項式の判別式から出る -/

section Bound

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- ★★★**頂点の跳びの上界 `(p−1)·i ≤ v_M(p)` は木に在った**。

`RamificationJumpBound.norm_natCast_le_pow_of_splits`(`f'(π) = ∏_{a≠π}(π−a)` と
`‖n·π^{n−1}‖ ≤ ‖f'(π)‖`)＋ `sub_one_mul_le_of_norm_natCast_eq_pow` を繋いだだけである。

★測定: `grep -n "sub_one_mul_le_of_norm_natCast_eq_pow" lean/ABC3/Found/PGC/RamificationJumpBound.lean`
(334 行)。★これで `exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps` の `hbnd` は
**`heM` と `M = F(π)` が `p` 次分離・分解すること**から供給される。 -/
theorem sub_one_mul_jump_le_of_splits {p i E : ℕ} (hp : 0 < p) {π : M}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hint : IsIntegral F π) (hdeg : (minpoly F π).natDegree = p)
    (hsep : (minpoly F π).Separable)
    (hsp : ((minpoly F π).map (algebraMap F M)).Splits)
    (hconj : ∀ a : M, ((minpoly F π).map (algebraMap F M)).IsRoot a → a ≠ π →
      ‖π - a‖ = ‖π‖ ^ (i + 1))
    (heM : ‖(p : M)‖ = ‖π‖ ^ E) :
    (p - 1) * i ≤ E :=
  sub_one_mul_le_of_norm_natCast_eq_pow hπ0 hπ1 heM
    (norm_natCast_le_pow_of_splits hp hπ0 hπ1 hval (minpoly.monic hint) hdeg hsep
      (minpoly.aeval F π) hsp hconj)

omit [IsUltrametricDist M] in
/-- `σ π` は最小多項式の根なので、`hconj`(全共役)から `hbreak`(1 つの `σ`)が出る。 -/
theorem norm_algHom_sub_eq_of_conj {i : ℕ} {π : M} (σ : M →ₐ[F] M)
    (hconj : ∀ a : M, ((minpoly F π).map (algebraMap F M)).IsRoot a → a ≠ π →
      ‖π - a‖ = ‖π‖ ^ (i + 1))
    (hne : σ π ≠ π) : ‖σ π - π‖ = ‖π‖ ^ (i + 1) := by
  have hroot : ((minpoly F π).map (algebraMap F M)).IsRoot (σ π) := by
    simp only [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def,
      Polynomial.aeval_algHom_apply, minpoly.aeval, map_zero]
  rw [← norm_sub_rev]
  exact hconj (σ π) hroot hne

end Bound

/-! ## §9 出口 —— `hbnd` と `hbreak` も落とした形 -/

section Exit3

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- ★★★★★★**`GainedBridgeSupply` の出口から `ℤ` 側の仮説を全部落とした形**。

`hjumps` 版からさらに `hbnd`(頂点の跳びの上界)と `hbreak`(1 つの `σ` の跳び)が消え、
代わりに `M = F(π)` が **`p` 次で分離・分解する**という構造の仮説
(`hsep` / `hsp` / `hconj` / `hne`)が入る。

★残るのは
`he`(`0 < e`)、★`hu`(`p^m ≤ u m`、跳びの列の伸び＝**Hasse–Arf の中身**)、
`hjump`(下付き分岐群の定義)、`hti`(`i` の定義)、
`hnormp` / `heM`(正規化)、`hval`(値群)、`hdeg` / `htop` / `hsep` / `hsp` / `hconj` / `hne`
(`M = F(π)` が `p` 次 Galois)である。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_conj
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M} {u : ℕ → ℤ}
    (τ : M →ₐ[F] M) (he : 0 < e)
    (hu : ∀ m, m ≤ k → (p : ℤ) ^ m ≤ u m)
    (hjump : ∀ j, j < k → ‖(τ ^ p ^ j) π - π‖ ≤ ‖π‖ ^ (u j + 1))
    (hti : (i : ℤ) = u k)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤)
    (hsep : (minpoly F π).Separable)
    (hsp : ((minpoly F π).map (algebraMap F M)).Splits)
    (hconj : ∀ a : M, ((minpoly F π).map (algebraMap F M)).IsRoot a → a ≠ π →
      ‖π - a‖ = ‖π‖ ^ (i + 1))
    (hne : (τ ^ p ^ k) π ≠ π) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have hp1 : 1 ≤ p := hp'.one_lt.le
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp'.one_lt
  have hpR0 : (0 : ℝ) < (p : ℝ) := lt_trans zero_lt_one hpR
  have hEne : p ^ (k + 1) * e ≠ 0 := Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne'
  have hπpow : ‖π‖ ^ (p ^ (k + 1) * e) = ((p : ℝ))⁻¹ := by rw [← heM]; exact hnormp
  have hπ0 : 0 < ‖π‖ := by
    rcases lt_or_eq_of_le (norm_nonneg π) with h | h
    · exact h
    · exfalso
      rw [← h, zero_pow hEne] at hπpow
      exact absurd hπpow.symm (ne_of_gt (inv_pos.mpr hpR0))
  have hπ1 : ‖π‖ < 1 := by
    refine GainedBridgeSupply.lt_one_of_pow_lt_one (E := p ^ (k + 1) * e) ?_
    rw [hπpow, inv_lt_one_iff₀]; exact Or.inr hpR
  have hint : IsIntegral F π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    exact hp'.ne_zero hdeg.symm
  have hbreak : ‖(τ ^ p ^ k) π - π‖ = ‖π‖ ^ (i + 1) :=
    norm_algHom_sub_eq_of_conj (τ ^ p ^ k) hconj hne
  have hnat : (p - 1) * i ≤ p ^ (k + 1) * e :=
    sub_one_mul_jump_le_of_splits hp'.pos hπ0 hπ1 hval hint hdeg hsep hsp hconj heM
  have hcast : ((p - 1 : ℕ) : ℤ) = (p : ℤ) - 1 := by
    push_cast [Nat.cast_sub hp1]; ring
  have hbnd : ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ) := by
    have h3 : (((p - 1) * i : ℕ) : ℤ) ≤ ((p ^ (k + 1) * e : ℕ) : ℤ) := by exact_mod_cast hnat
    rw [Nat.cast_mul, hcast, Nat.cast_mul, Nat.cast_pow] at h3
    rw [← hti]
    exact h3
  exact exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps τ he hu hjump hbnd hti hnormp heM
    hbreak hval hdeg htop x

end Exit3

/-! ## §10 数値の検算(★構成が古典的な列を壊さないこと / 壊れる列を直すこと) -/

namespace Numeric

/-- `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`: `p = 3`, `k = 1`, 下付き跳び `u = (u₀,u₁) = (2,8)`。
★構成は `t = (2,8)` をそのまま返す(古典的な場合は列を壊さない)。 -/
theorem seq_zeta27 :
    seq 3 (fun m => if m = 0 then (2 : ℤ) else 8) 1 0 = 2 ∧
    seq 3 (fun m => if m = 0 then (2 : ℤ) else 8) 1 1 = 8 := by
  constructor <;> norm_num [seq, chain]

/-- ★`u = (2,5)`, `p = 3`: `hstepZ` を実際の跳び `t = u` に当てると `3·2 = 6 ≤ 5` で**偽**。
★構成は `t = (min(2, ⌊5/3⌋), 5) = (1,5)` を作り、`3·1 = 3 ≤ 5` で通る。 -/
theorem seq_two_five :
    ¬ (3 * (2 : ℤ) ≤ 5) ∧
    seq 3 (fun m => if m = 0 then (2 : ℤ) else 5) 1 0 = 1 ∧
    3 * seq 3 (fun m => if m = 0 then (2 : ℤ) else 5) 1 0
      ≤ seq 3 (fun m => if m = 0 then (2 : ℤ) else 5) 1 1 := by
  refine ⟨by norm_num, by norm_num [seq, chain], by norm_num [seq, chain]⟩

end Numeric

/-! ## §11 `.src`(原典の対応箇所) -/

def exists_jump_seq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_iterate_pow_sub_self_le_lt.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def sub_one_mul_jump_le_of_splits.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_conj.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §12 使っている公理の一覧 -/

#print axioms chain
#print axioms le_chain
#print axioms pow_mul_chain_le
#print axioms exists_jump_seq
#print axioms norm_iterate_pow_sub_self_le_lt
#print axioms norm_sub_digit_zero_le_gainedLoss_lt
#print axioms exists_norm_sub_algebraMap_le_gainedLoss_lt
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_lt
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps
#print axioms sub_one_mul_jump_le_of_splits
#print axioms norm_algHom_sub_eq_of_conj
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_conj
#print axioms Numeric.seq_zeta27
#print axioms Numeric.seq_two_five
end GainedJumpSeq
end ABC3.Found.PGC
