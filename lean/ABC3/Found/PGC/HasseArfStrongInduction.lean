import ABC3.Found.PGC.FixedRingBaseAlgebra

/-!
# Hasse-Arf 段 2 —— `hind`(原文「by inductive hypothesis」)の除去

★★★**結論を先に書く: 帰納は要らなかった。**

Y15b `Found/PGC/HasseArfInduction.lean` と Y20 `Found/PGC/FixedRingBaseAlgebra.lean` が
残した最後の仮定

    hind : herbrandPhi π' H (n : ℝ) = (k : ℝ)      -- 原文「We have φ_H(n) ∈ Z≥0
                                                   --       by inductive hypothesis」

は、**強帰納法を回さずに** `§3` の `exists_natCast_herbrandPhi_of_jump` で
Proposition 6.9 から**直接**出る。したがって

* 本ファイルに `Nat.strong_induction_on` は 1 つも無い。
* 段 2 は `§6` の `exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top`
  で **`G` 可換 + `G = G_1` + `G_n ≠ G_{n+1}` だけから** `φ_G(n) ∈ ℤ≥0` を出す形になった。

★ファイル名は持ち場の名前(`HasseArfStrongInduction`)をそのまま使っている。
実際に入っているのは「強帰納法」ではなく「帰納法が要らないことの証明」である。

## ★★★原典の証明の穴(逸脱ではなく、原典側の不備の記録)

原文 (Yoshida08 p.16, 段 2):

> For j > 1, if Gn ̸= Gn+1 we can find H with G/H ∼= Z/pmiZ, and GnH/H ̸= Gn+1H/H.
> We have φH(n) ∈Z≥0 by inductive hypothesis, and
> (G/H)φH(n) ̸= (G/H)φH(n+1) = (G/H)φH(n)+1 by Proposition 6.9.

ここで「by inductive hypothesis」が引く帰納法の仮定は、Theorem 6.11 を `H`(巡回成分が
`j−1` 個)に当てたもの、すなわち

    H_n ≠ H_{n+1} ⟹ φ_H(n) ∈ ℤ≥0

である。ところが `H_i = H ∩ G_i` なので、**`H_n ≠ H_{n+1}` は成り立たないことがある**。

★具体例(反例ではなく「原文の論法が届かない例」):
`G ≅ (ℤ/p)²`、`G = G_1`、跳びが 2 つ `n₁ < n₂` で
`G_i = G (i ≤ n₁)`, `G_i = C (n₁ < i ≤ n₂)`, `G_i = 1 (i > n₂)`(`|C| = p`)。
`n = n₂` で原文の `H` は「`|H| = p` かつ `G_{n₂}H/H ≠ 1`」すなわち `C ⊄ H`、
つまり `H ≠ C`、したがって `H ∩ C = 1` を満たす。すると
`H_{n₂} = H ∩ C = 1 = H_{n₂+1}` となり、帰納法の仮定の前提が**成り立たない**。
しかも `φ_H(n₂) = n₁ + (n₂ − n₁)/p` で、これが整数であることは
`φ_G(n₂) = n₁ + (n₂ − n₁)/p ∈ ℤ`(いま示したい結論そのもの)と同値である。
★つまり原文の「by inductive hypothesis」は、この場合には**結論と同値のものを引いている**。

★★**本ファイルが与える正しい理由**(`§3`):
`x ∈ G_n ∖ H` と `G_{n+1} ≤ H` があれば、Proposition 6.9 から

* `x ∈ (G/H)_{φ_H(n)}`、`x ∉ (G/H)_{φ_H(n+1)}`

が出る。分岐群の定義 `(G/H)_r = {σ | i_ϖ(σ) ≥ r+1}` より、**整数** `t := i_ϖ(x) − 1` が

    φ_H(n) ≤ t < φ_H(n+1)

を満たす。Lemma 6.10(i) で `φ_H(n) = (Σ_{i=1}^{n} |H_i|)/|H|`、
`φ_H(n+1) − φ_H(n) = |H_{n+1}|/|H|` だから、`N := |H_{n+1}|` とおくと

    0 ≤ |H|·t − Σ_{i=1}^{n} |H_i| < N.

ところが `H_{n+1} ≤ H_i`(`i ≤ n`)と Lagrange より `N ∣ |H_i|`、また `N ∣ |H|` なので
左辺は `N` の倍数である。長さ `N` の半開区間に入る `N` の倍数は `0` だけなので
`φ_H(n) = t ∈ ℤ`。★**帰納法は 1 度も使っていない。**
★働いているのは「中間環 `𝒪_{K″} = 𝒪_{K′}^H` の分岐指数 `i_ϖ(x)` が整数である」という
**整数性だけ**である。

## 節の構成

| 節 | 何 | 語彙 |
|---|---|---|
| §1 | 抽象核 3 本 | ★**純算術 2 本 + 純群論 1 本。分岐・付値・Galois が 1 語も無い** |
| §2 | `φ_H` の整数値判定(Lemma 6.10(i) を通す) | Herbrand 関数 |
| §3 | ★★**`hind` の除去** —— Proposition 6.9 から整数の証人 `t` を作る | 分岐群 |
| §4 | §3 を固定環 `C = B^H` / 底の設定に代入 | 固定環の塔(Y14/Y16/Y17/Y20) |
| §5 | Y20 の段 2 から `hind` が消えた形 | |
| §6 | ★★★**段 2 全体** —— `G` 可換 + `G = G_1` + `G_n ≠ G_{n+1}` ⟹ `φ_G(n) ∈ ℤ≥0` | |

## 退化の自己検査

* ★★**`G` 可換を落とすと偽**(Hasse-Arf はアーベルでのみ成り立つ)。
  Y13/Y15/Y15b の流儀に合わせ `habel : ∀ x y : G, x * y = y * x` の**命題**で渡す
  (`CommGroup` 構造で渡すと具体層のインスタンスに当たらない)。
  `habel` は §6 で(a) 跳びを分離する巡回商 `H` の存在(Y15b `exists_cyclic_quotient_not_mem`)
  と (b) `H.Normal` の 2 箇所に効いている。
* ★★**`hne : G_n ≠ G_{n+1}` を落とすと空虚**。★原典のテキストは `pdftotext` で
  `Gn = Gn+1` にしか見えない(斜線が落ちる)。ここは `≠` で書いてある。
  §3–§5 では `hne` を「証人 `x ∈ G_n ∖ H`」の形で使う(`hxn` + `hxH`)。
* ★**`h1 : G_1 = ⊤`(原文「First assume G = G_1」)を落とすと段 2 は成り立たない。**
  これは Y15b の中で `G` が `p` 群であること・`φ_H(1) = 1` に効いている。
  ★**§3(`hind` の除去)だけは `h1` を使わない** —— 純粋に整除性だけで動く。
* ★**`hle : G_{n+1} ≤ H` を落とすと §3 が壊れる**。`x ∉ (G/H)_{φ_H(n+1)}` が言えなくなる。
* ★**`hxH : x ∉ H` を落とすと §3 が壊れる**(同上)。
* ★`ℕ∞` の切り詰め引き算・除算は 1 箇所も書いていない(`lean-idioms.md` #102)。
  `ℕ` の引き算は §1 の `eq_of_le_of_lt_add_of_dvd` の中の `b - a` だけで、
  そこは `a ≤ b` を仮定している。

## ★帰納の測度について(持ち場の問いへの答え)

持ち場は「`Nat.card H` に関する強帰納法」を見当としていたが、
★**測度は不要だった**(帰納自体が不要)。したがって
「`H ⊊ G` ⟹ `|H| < |G|`」も「`H` について `H = H_1`」も**必要にならなかった**。
参考までに後者は正しい: `H ≤ G = G_1` なら `H_1 = H ∩ G_1 = H`。

## ★まだ閉じていないもの(Theorem 6.11 全体)

段 1(`HasseArf.lean`)・段 2(本ファイル §6)は閉じた。段 3
(`exists_natCast_herbrandPhiGroup_of_tame_quotient`)は
`hk : φ_{G_1}(n) ∈ ℤ≥0` を仮定として持っている。これを段 1+2 から供給するには
★**「部分群 `H` を群として見る橋」**(木にまだ無い)が要る:

* `MulSemiringAction ↥H B`(制限作用。mathlib にインスタンスが無い)
* `herbrandPhi α H n = herbrandPhiGroup ↥H α n`
* `lowerRamificationGroup B ↥H i = (lowerRamificationGroup B G i).subgroupOf H`
* `FaithfulSMul ↥H B` / `SMulCommClass ↥H A B` の移送

★これは新しいノードである(本ファイルの持ち場ではない)。
段 3 の `hcomp` / `htopC` / `honeC`(順分岐商 `K′^{G_1}/K` の塔)も別に要る。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing Pointwise

def exists_natCast_herbrandPhi_of_jump.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

def exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top.src :
    ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

/-! ## §1 抽象核

★★**分岐・付値・Galois・Herbrand の語彙が 1 つも出てこない。**
最初の 2 本は純算術、3 本目は純群論(Lagrange)である。 -/

section AbstractCore

/-- ★★★**抽象核(本ファイルの心臓)** —— 長さ `N` の半開区間 `[a, a+N)` に入る
`N` の倍数は左端だけである。

原文が「by inductive hypothesis」で畳んだところは、実質この 1 行だった。 -/
theorem eq_of_le_of_lt_add_of_dvd {N a b : ℕ} (hab : a ≤ b) (hlt : b < a + N)
    (ha : N ∣ a) (hb : N ∣ b) : a = b := by
  have hd : N ∣ b - a := (Nat.dvd_sub_iff_left hab ha).mpr hb
  have hz : b - a = 0 := Nat.eq_zero_of_dvd_of_lt hd (by omega)
  omega

/-- ★★**抽象核** —— 有理数 `S/m` と自然数 `t` について
`S/m ≤ t < (S+N)/m` かつ `N ∣ S`・`N ∣ m` なら `S/m = t`。

★`N ∣ m` は `N ∣ m·t` に、`N ∣ S` はそのまま §1 の 1 本目に渡る。
★★`N = 0` は `hlt : b < a + 0` が偽になるので自動的に排除される
(具体層では `N = |H_{n+1}| ≥ 1` なので問題にならない)。 -/
theorem div_eq_natCast_of_dvd_of_le_of_lt {S m N t : ℕ} (hm : 0 < m) (hNS : N ∣ S) (hNm : N ∣ m)
    (h1 : (S : ℝ) / (m : ℝ) ≤ (t : ℝ))
    (h2 : (t : ℝ) < ((S : ℝ) + (N : ℝ)) / (m : ℝ)) :
    (S : ℝ) / (m : ℝ) = (t : ℝ) := by
  have hmR : (0 : ℝ) < (m : ℝ) := by exact_mod_cast hm
  have hle : S ≤ t * m := by exact_mod_cast (div_le_iff₀ hmR).1 h1
  have hlt : t * m < S + N := by exact_mod_cast (lt_div_iff₀ hmR).1 h2
  have hE := eq_of_le_of_lt_add_of_dvd hle hlt hNS (hNm.mul_left t)
  rw [hE]
  push_cast
  field_simp

/-- ★**抽象核(純群論)** —— 部分群の族 `F` と添字集合 `s` について、
`F j` がすべての `F i`(`i ∈ s`)に含まれるなら `|F j|` は `Σ_{i∈s} |F i|` を割る。

★Lagrange(`Subgroup.card_dvd_of_le`)を `Finset.dvd_sum` で束ねただけ。 -/
theorem natCard_dvd_sum_natCard {H : Type*} [Group H] {ι : Type*} (F : ι → Subgroup H)
    (s : Finset ι) (j : ι) (h : ∀ i ∈ s, F j ≤ F i) :
    Nat.card (F j) ∣ ∑ i ∈ s, Nat.card (F i) :=
  Finset.dvd_sum fun i hi => Subgroup.card_dvd_of_le (h i hi)

end AbstractCore

/-! ## §2 中間層 —— `φ_H(n)` が整数になる十分条件

Lemma 6.10(i)(`herbrandPhi_natCast`)で `φ_H` を「`|H_i|` の和 ÷ `|H|`」に開き、
§1 に渡すだけである。 -/

/-- `φ_H(n+1)` を `φ_H(n)` の分子に `|H_{n+1}|` を足した形で書く(Lemma 6.10(i) の 1 段)。 -/
theorem herbrandPhi_succ_natCast {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (H : Subgroup G) [Fintype H] (n : ℕ) :
    herbrandPhi α H ((n : ℝ) + 1)
      = ((∑ i ∈ Finset.Icc 1 n,
            (Nat.card ((lowerRamificationGroup B G i).subgroupOf H) : ℝ))
          + (Nat.card ((lowerRamificationGroup B G (n + 1)).subgroupOf H) : ℝ))
        / (Nat.card H : ℝ) := by
  have hc : ((n : ℝ) + 1) = ((n + 1 : ℕ) : ℝ) := by push_cast; ring
  rw [hc, herbrandPhi_natCast (A := A) huni hadj H (n + 1),
    Finset.sum_Icc_succ_top (by omega)]

/-- ★★★**`φ_H(n) ∈ ℤ` の判定** —— 区間 `[φ_H(n), φ_H(n+1))` に自然数 `t` が入るなら
`φ_H(n) = t` である。

★★理由は `§1` の 2 本:`|H|·φ_H(n) = Σ_{i=1}^{n} |H_i|` と
`|H|·(φ_H(n+1) − φ_H(n)) = |H_{n+1}|` があり、
`|H_{n+1}|` は `|H_i|`(`i ≤ n`、フィルトレーションが減少列なので Lagrange)と
`|H|` の両方を割るから。

★★**`φ_H` の増分が 1 未満であること(`|H_{n+1}| ≤ |H|`)は使っていない。**
効いているのは整除性だけである。 -/
theorem exists_natCast_herbrandPhi_of_mem_Ico {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (H : Subgroup G) [Fintype H] {n t : ℕ}
    (h1 : herbrandPhi α H (n : ℝ) ≤ (t : ℝ))
    (h2 : (t : ℝ) < herbrandPhi α H ((n : ℝ) + 1)) :
    herbrandPhi α H (n : ℝ) = (t : ℝ) := by
  classical
  set F : ℕ → Subgroup H := fun i => (lowerRamificationGroup B G i).subgroupOf H with hF
  set S : ℕ := ∑ i ∈ Finset.Icc 1 n, Nat.card (F i) with hS
  set N : ℕ := Nat.card (F (n + 1)) with hN
  have hSR : (∑ i ∈ Finset.Icc 1 n, (Nat.card (F i) : ℝ)) = (S : ℝ) := by
    rw [hS, Nat.cast_sum]
  have hphi : herbrandPhi α H (n : ℝ) = (S : ℝ) / (Nat.card H : ℝ) := by
    rw [herbrandPhi_natCast (A := A) huni hadj H n, hSR]
  have hphi1 : herbrandPhi α H ((n : ℝ) + 1) = ((S : ℝ) + (N : ℝ)) / (Nat.card H : ℝ) := by
    rw [herbrandPhi_succ_natCast (A := A) huni hadj H n, hSR]
  have hmono : ∀ i ∈ Finset.Icc 1 n, F (n + 1) ≤ F i := by
    intro i hi
    exact Subgroup.comap_mono (lowerRamificationGroup_antitone B G
      (by simp only [Finset.mem_Icc] at hi; omega))
  have hNS : N ∣ S := natCard_dvd_sum_natCard F (Finset.Icc 1 n) (n + 1) hmono
  have hNm : N ∣ Nat.card H := Subgroup.card_subgroup_dvd_card (F (n + 1))
  rw [hphi]
  refine div_eq_natCast_of_dvd_of_le_of_lt Nat.card_pos hNS hNm ?_ ?_
  · rw [← hphi]; exact h1
  · rw [← hphi1]; exact h2

/-! ## §3 具体層 —— ★★`hind` の除去

Proposition 6.9(`herbrand_coe_mul_coe_eq`、Y11)を仮定 `h69` として受け取り、
跳びの証人 `x ∈ G_n ∖ H`(`G_{n+1} ≤ H`)から**整数** `t = i_ϖ(x) − 1` を作る。 -/

/-- ★★★**原文「We have φH(n) ∈Z≥0 by inductive hypothesis」の正しい理由の前半**
—— 区間 `[φ_H(n), φ_H(n+1))` に自然数が入る。

★その自然数は `t := i_ϖ(x) − 1`、すなわち**中間環 `C` における `x` の分岐指数**である。
`i_ϖ` は `ℕ∞` に値を取るので `t` は自動的に整数であり、ここに帰納法は要らない。

* `x ∈ (G/H)_{φ_H(n)}` : `x ∈ G_n` と `h69` から(`x = x · 1`)。
* `x ∉ (G/H)_{φ_H(n+1)}` : `G_{n+1} ≤ H` なので `h69` の左辺は `H` そのものになり、
  `x ∉ H` から従う。
* `i_ϖ(x) ≠ ⊤` : ⊤ なら 2 つ目に反する。

★`φ_H(n) ≥ 0`(`φ_H(0) = 0` と狭義単調性)から `i_ϖ(x) ≥ 1` なので `t : ℕ` が取れる。 -/
theorem exists_nat_mem_Ico_herbrandPhi_of_jump {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {H : Subgroup G} [Fintype H]
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    {π' : B} (huni : maximalIdeal B = Ideal.span {π'})
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤) {ϖ : C}
    (h69 : ∀ r : ℝ, ((ramificationGroupReal (G := G) π' r : Subgroup G) : Set G) * (H : Set G)
      = ((ramificationGroupReal (G := G) ϖ (herbrandPhi π' H r) : Subgroup G) : Set G))
    {n : ℕ} {x : G} (hxn : x ∈ lowerRamificationGroup B G n) (hxH : x ∉ H)
    (hle : lowerRamificationGroup B G (n + 1) ≤ H) :
    ∃ t : ℕ, herbrandPhi π' H (n : ℝ) ≤ (t : ℝ) ∧ (t : ℝ) < herbrandPhi π' H ((n : ℝ) + 1) := by
  have hcast : ((n + 1 : ℕ) : ℝ) = (n : ℝ) + 1 := by push_cast; ring
  -- `φ_H(n) ≥ 0`
  have hnonneg : (0 : ℝ) ≤ herbrandPhi π' H (n : ℝ) := by
    have h0 : herbrandPhi π' H 0 = 0 := herbrandPhi_zero huni H
    have hm := (strictMono_herbrandPhi π' H).monotone
      (show (0 : ℝ) ≤ (n : ℝ) from Nat.cast_nonneg n)
    linarith [h0 ▸ hm]
  -- 原文「(G/H)_{φ_H(n)} ∋ x[bar]」
  have hmem : x ∈ ramificationGroupReal (G := G) ϖ (herbrandPhi π' H (n : ℝ)) := by
    have hx : x ∈ ((ramificationGroupReal (G := G) π' ((n : ℕ) : ℝ) : Subgroup G) : Set G)
        * (H : Set G) := by
      refine ⟨x, ?_, 1, H.one_mem, mul_one x⟩
      rw [SetLike.mem_coe, ramificationGroupReal_natCast (A := A) huni hadj n]
      exact hxn
    rw [h69 ((n : ℕ) : ℝ)] at hx
    exact hx
  -- 原文「(G/H)_{φ_H(n+1)} = G_{n+1}H/H = 1」⇒ `x[bar] ∉ (G/H)_{φ_H(n+1)}`
  have hnmem : x ∉ ramificationGroupReal (G := G) ϖ (herbrandPhi π' H ((n : ℝ) + 1)) := by
    have heq := h69 (((n + 1 : ℕ)) : ℝ)
    rw [ramificationGroupReal_natCast (A := A) huni hadj (n + 1), coe_mul_coe_eq_coe hle,
      hcast] at heq
    intro hc
    have hxmem : x ∈ (H : Set G) := by rw [heq]; exact hc
    exact hxH hxmem
  -- `i_ϖ(x)` は有限
  have htop : ramIndex ϖ x ≠ ⊤ := by
    intro hc
    exact hnmem (by rw [mem_ramificationGroupReal, hc]; exact RealLeENat.top)
  set m : ℕ := (ramIndex ϖ x).toNat with hmdef
  have hmeq : ramIndex ϖ x = (m : ℕ∞) := (ENat.coe_toNat htop).symm
  have hlow : herbrandPhi π' H (n : ℝ) + 1 ≤ (m : ℝ) := by
    have hh := hmem
    rw [mem_ramificationGroupReal, hmeq, realLeENat_coe] at hh
    exact hh
  have hhigh : (m : ℝ) < herbrandPhi π' H ((n : ℝ) + 1) + 1 := by
    by_contra hc
    exact hnmem (by rw [mem_ramificationGroupReal, hmeq, realLeENat_coe]; linarith)
  have hm1 : 1 ≤ m := by
    have h1R : (1 : ℝ) ≤ (m : ℝ) := by linarith
    exact_mod_cast h1R
  have hcast1 : ((m - 1 : ℕ) : ℝ) = (m : ℝ) - 1 := by
    rw [Nat.cast_sub hm1]; norm_num
  exact ⟨m - 1, by rw [hcast1]; linarith, by rw [hcast1]; linarith⟩

/-- ★★★★★**原文「We have φH(n) ∈Z≥0 by inductive hypothesis」——
帰納法を使わずに証明した形**。

原文 (Yoshida08 p.16):
> Theorem 6.11 (Hasse-Arf). If G is abelian, n ∈ Z[bb]_≥0 and G_n = G_n+1, then φ_G(n) ∈ Z[bb]_≥0.

★★逐語の `G_n = G_n+1` は `pdftotext` が `≠` の斜線を落とした形である。

本補題が消すのは段 2 の 1 文
「We have φH(n) ∈Z≥0 by inductive hypothesis」である。

★★ファイル冒頭に書いたとおり、原文が引く帰納法の仮定は `H_n ≠ H_{n+1}` を要求するが、
それは成り立たないことがある。ここでは §3 の整数の証人と §2 の整除性だけで結論する。

仮定:
* `h69` : Proposition 6.9(Y11 `herbrand_coe_mul_coe_eq`)。
* `hxn` / `hxH` / `hle` : 原文「GnH/H ̸= Gn+1H/H」が与える跳びの証人
  (Y15b `exists_cyclic_quotient_of_lowerRamificationGroup_ne` が供給する)。

★★**`G` の可換性も `G = G_1` も `p` 群であることも使っていない。** -/
theorem exists_natCast_herbrandPhi_of_jump {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {H : Subgroup G} [Fintype H]
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    {π' : B} (huni : maximalIdeal B = Ideal.span {π'})
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤) {ϖ : C}
    (h69 : ∀ r : ℝ, ((ramificationGroupReal (G := G) π' r : Subgroup G) : Set G) * (H : Set G)
      = ((ramificationGroupReal (G := G) ϖ (herbrandPhi π' H r) : Subgroup G) : Set G))
    {n : ℕ} {x : G} (hxn : x ∈ lowerRamificationGroup B G n) (hxH : x ∉ H)
    (hle : lowerRamificationGroup B G (n + 1) ≤ H) :
    ∃ k : ℕ, herbrandPhi π' H (n : ℝ) = (k : ℝ) := by
  obtain ⟨t, h1, h2⟩ :=
    exists_nat_mem_Ico_herbrandPhi_of_jump (A := A) (C := C) huni hadj h69 hxn hxH hle
  exact ⟨t, exists_natCast_herbrandPhi_of_mem_Ico (A := A) huni hadj H h1 h2⟩

/-! ## §4 固定環 `C = B^H` に代入する

Proposition 6.9 の塔の適合条件は、Y14/Y16/Y17 が `C = fixedRing B H` に対して供給する
(Y15b の `..._fixedRing` と同じ供給の仕方)。 -/

/-- ★★§3 を `C = 𝒪_{K′}^H` に代入した形。★塔の適合条件 7 本が消える。 -/
theorem exists_natCast_herbrandPhi_of_jump_fixedRing {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B]
    {H : Subgroup G} [H.Normal] [Fintype H]
    {π' : B} (hπ' : Irreducible π') {ϖ : ↥(fixedRing B H)}
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) →
      y ∈ Algebra.adjoin A ({algebraMap (↥(fixedRing B H)) B ϖ} : Set B))
    {n : ℕ} {x : G} (hxn : x ∈ lowerRamificationGroup B G n) (hxH : x ∉ H)
    (hle : lowerRamificationGroup B G (n + 1) ≤ H) :
    ∃ k : ℕ, herbrandPhi π' H (n : ℝ) = (k : ℝ) :=
  exists_natCast_herbrandPhi_of_jump (A := A) (C := ↥(fixedRing B H))
    ((irreducible_iff_uniformizer π').mp hπ') hadj
    (fun r => herbrand_coe_mul_coe_eq (A := A) algebraMap_smul_fixedRing smul_fixedRing_eq_self
      hπ' fixedRing_injective exists_algebraMap_fixedRing (exists_sub_mem_fixedRing hresA)
      (fun a => exists_algebraMap_fixedRing_eq a) rfl hadj hfix r)
    hxn hxH hle

/-- ★★★§3 を**底の設定だけ**で述べた形。★Y17 `fixedRing_mem_adjoin_uniformizer` で
`hfix` も消える。残るのは原典 §6.1 の設定(`hresA`:完全分岐、`hadj`:Lemma 5.11)だけ。 -/
theorem exists_natCast_herbrandPhi_of_jump_base {A B : Type*} [CommRing A] [IsLocalRing A]
    [CommRing B] [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B]
    {H : Subgroup G} [H.Normal] [Fintype H]
    {π' : B} (hπ' : Irreducible π') {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    {n : ℕ} {x : G} (hxn : x ∈ lowerRamificationGroup B G n) (hxH : x ∉ H)
    (hle : lowerRamificationGroup B G (n + 1) ≤ H) :
    ∃ k : ℕ, herbrandPhi π' H (n : ℝ) = (k : ℝ) :=
  exists_natCast_herbrandPhi_of_jump_fixedRing (A := A) hπ' hresA hadj
    (fixedRing_mem_adjoin_uniformizer hresA hAne hϖ) hxn hxH hle

/-! ## §5 段 2 から `hind` が消えた形

★Y20 `exists_natCast_herbrandPhiGroup_of_cyclic_quotient_base_of_injective` は
**書き換えない**。ここに新しい名前で `hind` の無い形を作る。 -/

/-- ★★★★**Yoshida 2008 Theorem 6.11 の段 2 —— `hind` が消えた形**。

原文 (Yoshida08 p.16):
> Theorem 6.11 (Hasse-Arf). If G is abelian, n ∈ Z[bb]_≥0 and G_n = G_n+1, then φ_G(n) ∈ Z[bb]_≥0.

★★逐語の `G_n = G_n+1` は `pdftotext` が `≠` の斜線を落とした形である。ここは
「跳びの証人 `x ∈ G_n ∖ H`」の形で `≠` 側を使っている。

Y20 の `..._base_of_injective` から `hind`(と `k`)が消えている。それ以外の仮定は同一:

* `hπ'` : `π′` は `𝒪_{K′}` の素元、`hϖ` : `ϖ` は `𝒪_{K″} = B^H` の素元。
* `hresA` : `K′/K` は完全分岐。`hadj` : `𝒪_{K′} = 𝒪_K[π′]`(Lemma 5.11)。
* `hAinj` : `𝒪_K → 𝒪_{K′}` が単射。
* `h1` : 原文「First assume G = G_1」。
* `hgen` / `hxn` / `hxH` / `hle` : 原文「we can find H with G/H ≅ Z/p^{m_i}Z,
  and G_nH/H ≠ G_{n+1}H/H」が与えるデータ。 -/
theorem exists_natCast_herbrandPhiGroup_of_cyclic_quotient_of_jump
    {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype H] [Fintype (G ⧸ H)]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π') {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h1 : lowerRamificationGroup B G 1 = ⊤)
    {σ : G} (hgen : ∀ g : G, ∃ t : ℤ, (σ ^ t)⁻¹ * g ∈ H)
    {n : ℕ} {x : G} (hxn : x ∈ lowerRamificationGroup B G n) (hxH : x ∉ H)
    (hle : lowerRamificationGroup B G (n + 1) ≤ H) :
    ∃ j : ℕ, herbrandPhiGroup G π' (n : ℝ) = (j : ℝ) := by
  obtain ⟨k, hk⟩ := exists_natCast_herbrandPhi_of_jump_base (A := A) hπ' hϖ hresA hadj
    (exists_mem_maximalIdeal_map_ne_zero hAinj) hxn hxH hle
  exact exists_natCast_herbrandPhiGroup_of_cyclic_quotient_base_of_injective p hp hπ' hϖ hresA
    hadj hAinj h1 hgen hk hxn hxH hle

/-! ## §6 段 2 全体 -/

/-- ★★★★★**Yoshida 2008 Theorem 6.11 (Hasse-Arf) の段 2 —— 全体**。

原文 (Yoshida08 p.16):
> Theorem 6.11 (Hasse-Arf). If G is abelian, n ∈ Z[bb]_≥0 and G_n = G_n+1, then φ_G(n) ∈ Z[bb]_≥0.

★★逐語の `G_n = G_n+1` は `pdftotext` が `≠` の斜線を落とした形である。

原文の「First assume G = G_1」の場合(`h1`)を、原典 §6.1 の底の設定だけから閉じた形:

* `habel` : ★★**`G` 可換**(命題で渡す)。落とすと偽。
* `h1` : ★**`G = G_1`**。落とすと段 2 は成り立たない(段 3 の場合になる)。
* `hne` : ★★**`G_n ≠ G_{n+1}`**。落とすと空虚。
* `hπ'` / `hresA` / `hadj` / `hAinj` : 原典 §6.1 の設定
  (`π′` は素元、`K′/K` 完全分岐、`𝒪_{K′} = 𝒪_K[π′]`、`𝒪_K → 𝒪_{K′}` 単射)。
* `p` / `hp` / `CharP` : 剰余体の標数。

★段取り: Y15b `exists_cyclic_quotient_of_lowerRamificationGroup_ne` で
`H`・生成元 `σ`・跳びの証人 `x` を作り、`H.Normal`(可換性)と
`Fintype`・素元 `ϖ`(固定環が DVR)を補って §5 に渡す。
★★**`hind` は §3 が供給するので、ここに帰納法は無い。** -/
theorem exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top
    {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (habel : ∀ x y : G, x * y = y * x)
    (h1 : lowerRamificationGroup B G 1 = ⊤)
    {n : ℕ} (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    ∃ j : ℕ, herbrandPhiGroup G π' (n : ℝ) = (j : ℝ) := by
  classical
  obtain ⟨H, σ, x, hle, hxn, hxH, hgen⟩ :=
    exists_cyclic_quotient_of_lowerRamificationGroup_ne (B := B) habel hne
  haveI hnorm : H.Normal :=
    ⟨fun a ha g => by rw [habel g a, mul_assoc, mul_inv_cancel, mul_one]; exact ha⟩
  haveI : Fintype ↥H := Fintype.ofFinite _
  haveI : Fintype (G ⧸ H) := Fintype.ofFinite _
  obtain ⟨ϖ, hϖ⟩ := exists_irreducible ↥(fixedRing B H)
  exact exists_natCast_herbrandPhiGroup_of_cyclic_quotient_of_jump (A := A) p hp hπ' hϖ hresA
    hadj hAinj h1 hgen hxn hxH hle

/-- ★段 2 の結論の非負性(原文の `φ_G(n) ∈ Z≥0` の `≥0`)。 -/
theorem herbrandPhiGroup_nonneg_of_lowerRamificationGroup_one_eq_top
    {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (habel : ∀ x y : G, x * y = y * x)
    (h1 : lowerRamificationGroup B G 1 = ⊤)
    {n : ℕ} (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    0 ≤ herbrandPhiGroup G π' (n : ℝ) :=
  nonneg_of_exists_natCast (exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top
    (A := A) p hp hπ' hresA hadj hAinj habel h1 hne)

end ABC3.Found.PGC
