import ABC3.Found.PGC.AxEpsilonDecay
import ABC3.Found.PGC.TotallyRamifiedLayer
import ABC3.Found.PGC.WildDepthFieldDescent

/-!
# [pGC] Sylow 降下 → `AxWildDescent` —— ★段数は合う。合わないのは**定数**である

## ★配られた字面の検算（(a) は偽、(b) は既に木に在った）

★**(a)「`exists_pgroup_descent_cyclic` が出すのは 1 段（`|Q/P| = p`）だが、
出口は `Nat.card = p^{k+1}` を要求するので段数が合わない」は偽である。**

`AxWildDescent`（`Found/PGC/AxTowerDecay.lean:473`）の結論は
`∃ x', wildDepth K x' < wildDepth K x ∧ …` であって、★**1 段しか要求しない**。
塔の再帰は `axLemma_of_wildDescent`（`AxTowerDecay.lean:491`）の側にある。
⇒ 使うべき出口は `p^{k+1}` の
`TowerDataFromCyclic.exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_adjoin` ではなく、
★`k = 0` の `TotallyRamifiedLayer.lean:305`
`exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p` である。
群の側の「1 段」は本ファイル `exists_padicValNat_index_lt` で
「`v_p(Q.index)` が `v_p(H.index)` より真に小さい」として書いた（`omega` 1 行）。

★**(b)「`exists_pgroup_descent` を `AxWildDescent` に繋ぐ」はすでに木にある。**
`WildDepthFieldDescent.lean:141 axWildDescent_prime : AxWildDescent K (fun _ => (p:ℝ))`
が**無条件**に成り立っている。★したがって本波の持ち場が新しく買うものは
**定数だけ**である。★この点は `AxEpsilonDecay.lean:177` の断定と一致した
（★断定をそのまま信じず、定義と上流の形を読んで確かめた）。

## ★★★本波の新しい測定（否定的な結果）

次数 `p` の層 1 枚の出口（`TotallyRamifiedLayer.lean:305`）が与える損失は
`axDecay p 1 · ‖g x − x‖` であり、★**深さ `k` に依らない**。
一方 `axLemma_of_axDecay`（`AxEpsilonDecay.lean:532`）が要求するのは
`∀ k, c k ≤ axDecay p k` である。本ファイルはこの隔たりを測った。

| 宣言 | 内容 |
|---|---|
| `axDecay_succ` | ★`axDecay p (k+1) = (axDecay p 1) ^ ((1/p)^k)`（閉じた形） |
| `axDecay_succ_pow` | `(axDecay p (k+1)) ^ (p^k) = axDecay p 1` |
| `axDecay_one_le_p` | `axDecay p 1 ≤ p`（★定数の改良は真である） |
| `axDecay_succ_lt_axDecay_one` | ★`k ≥ 1` で `axDecay p (k+1) < axDecay p 1` |
| `exists_axDecay_lt` | ★★`1 < C` なら `axDecay p (k+1) < C` となる `k` が在る |
| `not_forall_le_axDecay` | ⇒ ★**どんな一様定数 `C > 1` も `axLemma_of_axDecay` の `hdecay` を満たさない** |
| `not_forall_prod_const_le` | ⇒ `axLemma_of_wildDescent` の積の仮説も満たさない |
| `not_forall_prod_Icc_const_le` | ⇒ ★木で一番弱い `axLemma_of_wildDescent_Icc`（`AxEpsilonDecay.lean:440`）でも満たさない |
| `axDecay_two_one` | 非空虚性の数値検算: `axDecay 2 1 = 2`（積は `2^n`） |

⇒ ★★★**一様定数の 1 段降下は、定数をいくら良くしても `AxLemma` に届かない**
（`const_descent_no_go`。★積の形・`Icc` の形・`axDecay` の形の 3 通りとも塞がっている）。
`axWildDescent_prime` の `p` を `axDecay p 1 = p^{1/(p−1)}` に絞るだけでは**足りない**。

★足りない分はちょうど指数の `(1/p)^k` であり（`axDecay_succ`）、
それは「1 段の損失」ではなく「その段で使う `g` について `‖g x − x‖` が
`ε` より幾何的に小さいこと」から来る
（`AxEpsilonDecay.lean:294 norm_iterate_pow_sub_self_le`、
`AxEpsilonDecay.lean:427 norm_smul_pow_prime_sub_le`）。
★★`exists_pgroup_descent_cyclic`（`TotallyRamifiedLayer.lean:460`）は
`Q/P` が位数 `p` の巡回であることしか言わず、
★その生成元が絶対 Galois 群の元の `p^k` 乗であることは**言わない**。
⇒ ★★★**ここが止まる場所である。**

## ★①不分岐側との関係（持ち場の指示）

★**独立に生きている。** `TotallyRamifiedLayer.lean:342 not_valK_of_norm_eq` は
「次数 `p` の層が不分岐なら出口の `hvalK` は偽」をすでに定理にしている
（★本波で読んで確かめた。今も真である）。
★本ファイルの否定は**全分岐の側でも**成り立つ（定数の話であって分岐の話ではない）。
⇒ ★**不分岐側を閉じても定数の穴は残る。2 つは別の穴である。**
（`TotallyRamifiedLayer.lean:37-41` の「`p ∤ deg minpoly` と『層が不分岐』を
同一視しないこと」という警告も、本ファイルは踏んでいない。）

## ★訂正（他ファイルの docstring は書き換えない。名指しでここに書く）

`Found/PGC/TowerDataFromCyclic.lean` 冒頭の
「★★★残るのは『`p`-群の中に巡回な商の列を取る』ノードである」は★**足りない**
（★本連鎖が前波に自分で書いた断定である）。
巡回な商の列を取っても、上のとおり**定数が深さに依らない**ので `AxLemma` には届かない。
★正しくは「巡回な商の列 ＋ **その生成元が `σ^{p^k}` の形であること**」である。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. 本ファイルの中身（定数列の収束の勘定）は原典 (Ax 1970 / pGC Cor 3.1) に
   対応する文が無い。`.src` は `AxEpsilonDecay.lean` と同じ項目を指す。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. §1 の抽象核は `ℝ` と `Finset` しか使わず、分岐・付値・Galois の語彙が 1 語も出ない。
   §2 は純群論。分岐の語彙が出るのは docstring だけである。
-/

namespace ABC3.Found.PGC

namespace PGroupToAxWild

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-! ## §1 抽象核 —— 定数列の積は有界でない（`ℝ` と `Finset` だけ） -/

section Core

/-- ★抽象核: `1 < C` なら、定数 `C` の有限積は `B` で一様に抑えられない。

★これが `axLemma_of_wildDescent`（`AxTowerDecay.lean:491`）の仮説
`∀ s : Finset ℕ, ∏ k ∈ s, c k ≤ C` が定数列で満たせないことの中身である。 -/
theorem not_forall_prod_const_le {C B : ℝ} (hC : 1 < C) :
    ¬ (∀ s : Finset ℕ, ∏ _k ∈ s, C ≤ B) := by
  intro h
  obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt B hC
  have := h (Finset.range n)
  rw [Finset.prod_const, Finset.card_range] at this
  linarith

end Core

/-! ## §2 群の段 —— Sylow 降下は `v_p(index)` を真に下げる（★段数は合う） -/

section Group

variable {G : Type*} [Group G] [Finite G] {p : ℕ} [Fact p.Prime]

/-- ★★**段数の噛み合わせ（本波の測定その 1）** —— Sylow 降下 1 回で
`v_p(index)` は**真に**下がる。

`wildDepth K x` は `Stab(x)` の指数の `v_p` として測られる
（`WildDepthFieldDescent.lean` の配管）ので、これがそのまま
`AxWildDescent`（`AxTowerDecay.lean:473`）が要求する
`wildDepth K x' < wildDepth K x` の**1 段**に対応する。
⇒ ★**「1 段しか出ないから段数が合わない」は誤りである**（`AxWildDescent` は 1 段しか要らない）。 -/
theorem exists_padicValNat_index_lt (H : Subgroup G) (k : ℕ)
    (hk : padicValNat p H.index = k + 1) :
    ∃ P Q : Subgroup G, P ≤ H ∧ P ≤ Q ∧ (P.subgroupOf Q).Normal ∧
      Nat.card (↥Q ⧸ P.subgroupOf Q) = p ∧
      padicValNat p Q.index < padicValNat p H.index := by
  obtain ⟨P, Q, hPH, hPQ, _hPv, hQv, hnorm, hcard⟩ :=
    TotallyRamifiedLayer.exists_pgroup_descent_cyclic (p := p) H k hk
  exact ⟨P, Q, hPH, hPQ, hnorm, hcard, by rw [hQv, hk]; omega⟩

end Group

/-! ## §3 定数 —— `axDecay` の閉じた形と大小 -/

section Const

variable {p : ℕ} [Fact p.Prime]

theorem one_lt_axDecay_one : 1 < axDecay p 1 := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  rw [axDecay_one]
  refine Real.one_lt_rpow_iff_of_pos (by linarith) |>.mpr ?_
  refine Or.inl ⟨h1, ?_⟩
  apply div_pos one_pos; linarith

omit [Fact p.Prime] in
/-- ★★**`axDecay` の閉じた形** —— 深さ `k+1` の予算は深さ `1` の予算の `(1/p)^k` 乗。

`axDecay p (k+1) = (axDecay p 1) ^ ((1/p)^k)`（右辺は `Real.rpow`）。
★これが「次数 `p` の層 1 枚（`TotallyRamifiedLayer.lean:305`）が与える `axDecay p 1`」と
「`axLemma_of_axDecay`（`AxEpsilonDecay.lean:532`）が要求する `axDecay p (k+1)`」の
**隔たりを閉じた形で表した**式である。 -/
theorem axDecay_succ (k : ℕ) :
    axDecay p (k + 1) = (axDecay p 1) ^ ((1 / (p:ℝ)) ^ k) := by
  have hp0 : (0:ℝ) ≤ (p:ℝ) := Nat.cast_nonneg p
  simp only [axDecay, Nat.add_sub_cancel, Nat.sub_self, pow_zero, mul_one]
  rw [← Real.rpow_mul hp0]

theorem axDecay_one_le_p : axDecay p 1 ≤ (p:ℝ) := by
  have hp2 : (2:ℝ) ≤ (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).two_le
  have h : (p:ℝ) ^ (1 / ((p:ℝ) - 1)) ≤ (p:ℝ) ^ (1:ℝ) := by
    refine Real.rpow_le_rpow_of_exponent_le (by linarith) ?_
    rw [div_le_one (by linarith)]
    linarith
  rw [axDecay_one]
  simpa using h

/-- ★★深さ `2` 以上では、次数 `p` の層 1 枚が与える定数は**大きすぎる**。 -/
theorem axDecay_succ_lt_axDecay_one {k : ℕ} (hk : 0 < k) :
    axDecay p (k + 1) < axDecay p 1 := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hlt : ((1:ℝ) / (p:ℝ)) ^ k < 1 := by
    refine pow_lt_one₀ (by positivity) ?_ (by omega)
    rw [div_lt_one (by linarith)]; linarith
  have h : axDecay p 1 ^ (((1:ℝ) / (p:ℝ)) ^ k) < axDecay p 1 ^ (1:ℝ) :=
    Real.rpow_lt_rpow_of_exponent_lt one_lt_axDecay_one hlt
  rw [axDecay_succ]
  simpa using h

theorem not_forall_prod_axDecay_one_le (B : ℝ) :
    ¬ (∀ s : Finset ℕ, ∏ _k ∈ s, axDecay p 1 ≤ B) :=
  not_forall_prod_const_le one_lt_axDecay_one

end Const

/-! ## §4 `AxWildDescent` は定数について単調 -/

section Mono

variable {p : ℕ} [Fact p.Prime]

/-- `AxWildDescent` は定数列について単調（★定数を**緩める**方向にしか動かない）。

⇒ `axWildDescent_prime`（定数 `p`）から `axDecay p` 版を得るには
`∀ k, (p:ℝ) ≤ axDecay p k` が要るが、それは `not_forall_le_axDecay` で**偽**である。 -/
theorem axWildDescent_mono {K : PAdicLocalField p} {c d : ℕ → ℝ}
    (hcd : ∀ k, c k ≤ d k) (h : AxWildDescent K c) : AxWildDescent K d := by
  intro ε hε x hx hdvd
  obtain ⟨x', hlt, h1, h2⟩ := h ε hε x hx hdvd
  exact ⟨x', hlt, h1.trans (mul_le_mul_of_nonneg_right (hcd _) hε),
    fun σ => (h2 σ).trans (mul_le_mul_of_nonneg_right (hcd _) hε)⟩

end Mono

/-! ## §5 `Finset.Icc` 版（木で一番弱い仮説でも塞がる） -/

section Icc

variable {p : ℕ} [Fact p.Prime]

omit [Fact p.Prime] in
theorem prod_Icc_const {C : ℝ} (n : ℕ) : ∏ _k ∈ Finset.Icc 1 n, C = C ^ n := by
  rw [Finset.prod_const, Nat.card_Icc, Nat.add_sub_cancel]

omit [Fact p.Prime] in
theorem not_forall_prod_Icc_const_le {C B : ℝ} (hC : 1 < C) :
    ¬ (∀ n : ℕ, ∏ _k ∈ Finset.Icc 1 n, C ≤ B) := by
  intro h
  obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt B hC
  have := h n
  rw [prod_Icc_const] at this
  linarith

theorem not_forall_prod_Icc_axDecay_one_le (B : ℝ) :
    ¬ (∀ n : ℕ, ∏ _k ∈ Finset.Icc 1 n, axDecay p 1 ≤ B) :=
  not_forall_prod_Icc_const_le one_lt_axDecay_one

/-- ★非空虚性の数値検算（`p = 2`）: `axDecay 2 1 = 2`。
⇒ `∏_{k ∈ Icc 1 n} axDecay 2 1 = 2^n`（`prod_Icc_const`）で、確かに非有界である。 -/
theorem axDecay_two_one : axDecay 2 1 = 2 := by
  haveI : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  rw [axDecay_one]
  norm_num

end Icc

/-! ## §6 ★★`axDecay` は 1 に収束する —— 下から一様定数で押さえられない -/

section Limit

variable {p : ℕ} [Fact p.Prime]

theorem axDecay_succ_pow (k : ℕ) : (axDecay p (k + 1)) ^ (p ^ k) = axDecay p 1 := by
  have hA0 : (0:ℝ) ≤ axDecay p 1 := le_trans zero_le_one (one_le_axDecay p 1)
  have hp0 : (0:ℝ) < (p : ℝ) := by
    exact_mod_cast (Fact.out : p.Prime).pos
  have hmul : ((1:ℝ) / (p:ℝ)) ^ k * ((p ^ k : ℕ) : ℝ) = 1 := by
    push_cast
    rw [← mul_pow, one_div, inv_mul_cancel₀ (ne_of_gt hp0), one_pow]
  rw [axDecay_succ, ← Real.rpow_natCast _ (p ^ k), ← Real.rpow_mul hA0, hmul,
    Real.rpow_one]

/-- ★★**`axDecay p` は 1 に収束する** —— `1` より大きいどんな定数も、
十分深いところで `axDecay` に追い越される。

★証明は `log` を使わない: `axDecay_succ_pow` で `p^k` 乗して
`axDecay p 1 < C ^ (p^n)` に帰着する（`Real.log` / `Real.logb` は
本ファイルの import では `Unknown constant` になる。#68）。 -/
theorem exists_axDecay_lt {C : ℝ} (hC : 1 < C) : ∃ k : ℕ, axDecay p (k + 1) < C := by
  have hp1 : 1 < p := (Fact.out : p.Prime).one_lt
  obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt (axDecay p 1) hC
  have hnp : n ≤ p ^ n := le_of_lt (Nat.lt_pow_self hp1)
  have hstep : C ^ n ≤ C ^ (p ^ n) := pow_le_pow_right₀ (le_of_lt hC) hnp
  refine ⟨n, ?_⟩
  have hpow : (axDecay p (n + 1)) ^ (p ^ n) < C ^ (p ^ n) := by
    rw [axDecay_succ_pow]
    linarith
  exact lt_of_pow_lt_pow_left₀ (p ^ n) (le_of_lt (lt_trans zero_lt_one hC)) hpow

end Limit

/-! ## §7 ★★★一様定数の降下の no-go（3 通りとも塞がる） -/

section NoGo

variable {p : ℕ} [Fact p.Prime]

theorem not_forall_le_axDecay {C : ℝ} (hC : 1 < C) :
    ¬ (∀ k : ℕ, C ≤ axDecay p k) := by
  intro h
  obtain ⟨k, hk⟩ := exists_axDecay_lt (p := p) hC
  exact absurd (h (k + 1)) (not_le.mpr hk)

/-- ★★★**本ファイルの主結果** —— 一様定数 `C > 1` の 1 段降下は、
木にある 3 通りの出口の**どれにも**乗らない。

1. `axLemma_of_wildDescent`（`AxTowerDecay.lean:491`）の `∀ s : Finset ℕ, ∏ ≤ C`
2. `axLemma_of_wildDescent_Icc`（`AxEpsilonDecay.lean:440`）の `∀ n, ∏_{Icc 1 n} ≤ C`
3. `axLemma_of_axDecay`（`AxEpsilonDecay.lean:532`）の `∀ k, c k ≤ axDecay p k`

⇒ ★`WildDepthFieldDescent.lean:141 axWildDescent_prime`（定数 `p`）を
次数 `p` の層の出口（`TotallyRamifiedLayer.lean:305`、定数 `axDecay p 1`）で
置き換えても**足りない**。★必要なのは**深さに依る**定数であり、
その源は「その段の `g` について `‖g x − x‖` が `ε` より幾何的に小さいこと」
（`AxEpsilonDecay.lean:294` / `:427`）である。 -/
theorem const_descent_no_go {C : ℝ} (hC : 1 < C) :
    (∀ B : ℝ, ¬ (∀ s : Finset ℕ, ∏ _k ∈ s, C ≤ B))
      ∧ (∀ B : ℝ, ¬ (∀ n : ℕ, ∏ _k ∈ Finset.Icc 1 n, C ≤ B))
      ∧ ¬ (∀ k : ℕ, C ≤ axDecay p k) :=
  ⟨fun _ => not_forall_prod_const_le hC, fun _ => not_forall_prod_Icc_const_le hC,
    not_forall_le_axDecay hC⟩

end NoGo

/-! ## §8 `.src`（原典の対応箇所）と 使っている公理の一覧 -/

def not_forall_prod_const_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_padicValNat_index_lt.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def axDecay_succ.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_axDecay_lt.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def const_descent_no_go.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms not_forall_prod_const_le
#print axioms exists_padicValNat_index_lt
#print axioms one_lt_axDecay_one
#print axioms axDecay_succ
#print axioms axDecay_one_le_p
#print axioms axDecay_succ_lt_axDecay_one
#print axioms not_forall_prod_axDecay_one_le
#print axioms axWildDescent_mono
#print axioms prod_Icc_const
#print axioms not_forall_prod_Icc_const_le
#print axioms not_forall_prod_Icc_axDecay_one_le
#print axioms axDecay_two_one
#print axioms axDecay_succ_pow
#print axioms exists_axDecay_lt
#print axioms not_forall_le_axDecay
#print axioms const_descent_no_go

end PGroupToAxWild

end ABC3.Found.PGC
