import ABC3.Found.PGC.LubinTateRamificationBreak

/-!
# `hρ` を外す —— Yoshida 2008 Proposition 6.14(`n = 1`)の部分群の帳簿

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Proposition 6.14(物理 p.17)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-6.html` の `#prop-6-14`。

原文 (Yoshida08 p.17):
> Proposition 6.14. Let x ∈ K^× with v(x) = n > 0. Let L = K_n and K^m_x as in Definition 5.3. Then we have Gal(K^m_x/L)^m = {id} for all m ≥ 1 (see Definition 6.12).

## 位置づけ —— 新しい数学は無い

道 B の残りは 2 つで、本ファイルはそのうち片方である。

* `Found/PGC/LubinTateUpperRamificationVanish.lean` が Prop 6.14(`n = 1` 形)を
  **段 1 だけを仮定 `hρ` に残して**閉じている。
* `Found/PGC/LubinTateRamificationBreak.lean` がその段 1 の数学
  (`i(σ) = q^{v_K(u−1)}`)を閉じている。

本ファイルは**その 2 つを繋ぐ部分群の帳簿だけ**を書く。原典で言えば

    Thus for G = Gal(K^m_x/L) and 1 ≤ i ≤ m, we have |G_n| = |ρ^{-1}_{f,m}(1+p^i)| = q^{m−i}
    for q^{i−1} − 1 < n ≤ q^i − 1.

の、**前半の集合の等式** `G_n = ρ^{-1}_{f,m}(1+p^i)` にあたる 1 行である
(後半の位数 `= q^{m−i}` は `LubinTateUpperRamificationVanish.lean` の
`natCard_map_principalUnits` が既に持っている)。

★添字は 0 始まりにずらしてある: 原典の `1 ≤ i ≤ m` はここでの `i+1`
(`0 ≤ i ≤ M`, `m = M+1`)。原典の範囲 `q^{i−1} − 1 < n ≤ q^i − 1` は
ℕ の切り詰め引き算を出さない形 `q^i ≤ n < q^{i+1}` に書き換えてある。

★**本ファイルは下付き分岐群 `G_n` の帳簿である。** 上付き `G^t` が出るのは
最後の 2 本(消費側 `upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap`
にそのまま流すところ)だけで、そこは既存の補題を呼ぶだけである。

## 抽象核(§0。分岐・付値・Lubin-Tate・Galois の語彙が 1 つも出てこない)

* `le_of_antitone_natPred` —— `P (i+1) → P i` なる ℕ 上の述語は下向きに閉じている。
* `exists_le_iff_of_antitone_natPred` —— そのような `P` が `P 0` を満たし `P m` を
  満たさないなら、**`P j ↔ j ≤ i₀` となる `i₀ < m` が存在する**。
  ★これが「フィルターの段が一意に決まる」の中身である。`∃ i, P i ∧ ¬ P (i+1)`
  ではなく**同値の形**で出すのが要点で、こうしておくと
  「`u ∈ 1+𝔭^{i+1}` ⟺ `i+1 ≤ i₀`」がそのまま取れる(下の `hV`)。
* `mem_map_mk'_iff_of_le` —— `N ≤ H` なら `ū ∈ H.map (mk' N) ↔ u ∈ H`。
  `Subgroup.comap_map_eq_self` + `QuotientGroup.ker_mk'` の 1 行。
  ★`N ≤ H` は落とせない(`H = 1`, `N ≠ 1` で反例)。だから下の帳簿では
  `j ≤ M+1` を付けている。

★3 本とも群論・ℕ の算術だけで、原典の設定に一切依らない。

## 段取り(3 手。前の実装者が残した分解のとおり)

1. `ramIndex_torsionGen_dichotomy` —— `σ` ごとの二分律。
   `ρ_{f,m}(σ)` の代表 `u ∈ 𝒪_K^×` を取り、
   * `u ∈ 1+𝔭^m` なら `ρ(σ) = 1`、`ρ` は同型なので `σ = 1`;
   * そうでなければ**一意な `i₀ ≤ M`** があって `i(σ) = q^{i₀}` かつ
     `ρ(σ) ∈ (1+𝔭^j)/(1+𝔭^m) ⟺ j ≤ i₀`(`j ≤ M+1`)。
2. `lowerRamificationGroupAdjoin_eq_comap_map_principalUnits` ——
   `G_n = ρ^{-1}((1+𝔭^{i+1})/(1+𝔭^{M+1}))`(`q^i ≤ n < q^{i+1}`)。
   B4 の抽象核 `lowerRamificationGroup_eq_of_ramIndex_pow` に
   `hpow`(二分律の前半)と `hV`(後半)を流すだけ。
3. `upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_uniformizer` ——
   2 を `hρ` として消費側に差し込む。★**`hρ` は消えた。**

## 「一意な `i₀`」はどこから出るか(退化の自己検査)

`P j := (u ∈ principalUnits K π j)` と置くと

* `P` は**真に降下する列の指示関数**である: `principalUnits_antitone` が
  `j ≤ j' → 1+𝔭^{j'} ⊆ 1+𝔭^j` を与えるので `P (j+1) → P j`;
* `P 0` は**常に真**(`𝔭^0 = 𝒪_K` なので `u = 1 + (u−1)·1`。`mem_principalUnits_zero`);
* `P (M+1)` が偽なのが「`σ ≠ 1`」の場合である。

この 3 つから `exists_le_iff_of_antitone_natPred` が `i₀` を出す。★`i₀` は
**`P j ↔ j ≤ i₀` を満たすものとして一意**である(2 つあれば互いに `≤`)。
これが原典の `v_K(u−1) = i` の中身であり、原典が黙って使っている
「`v_K` の値が一意に決まる」に対応する。

## `m = 0` の除外

本ファイルの `m` は常に `M + 1 ≥ 1` である(`M : ℕ`)。`m = 0` では
`µ_{f,0} = {0}` で拡大が自明になり主張が空虚になるので、段 1 と同じく除いてある。
`0 ≤ i₀ < m` は `M = i₀ + k` と書くことで埋め込んであり、
★**ℕ / ℕ∞ の切り詰め引き算・除算は 1 度も書いていない。**

## 逸脱の記録

1. ★**`hρ` は仮定ではなくなった。** 消費側
   `upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap` の最後の仮定
   `hρ` を本ファイルが供給する。残る仮定は原典の設定と
   「`α` が `𝒪_{K^m_x}` の素元」(`huni`)だけであり、後者すら
   `upperRamificationGroup_torsionGen_eq_bot` では外してある。
2. ★**二分律の第 2 主張は `j ≤ M + 1` に制限してある。**
   `mem_map_mk'_iff_of_le` が `1+𝔭^{M+1} ≤ 1+𝔭^j` を要求するためで、
   消費側が使うのは `j = i+1 ≤ M+1` だけなので影響しない。
   (`j > M+1` でも `ρ(σ) ∈ (1+𝔭^j)/(1+𝔭^{M+1}) ↔ u ∈ 1+𝔭^{M+1}` から
   両辺偽で成り立つが、`H ⊔ N` を経由するので書かない。)
3. ★**帳簿の仮定は `i ≤ M` である。** 消費側の `hρ` は `i + j = M` の形で
   書かれているが `j` は使わないので、より弱い `i ≤ M` で立てた
   (消費側には `(by omega)` で渡している)。
4. ★**`α` の取り替えは問題にならない。** `lowerRamificationGroupAdjoin` は
   `α` を含まない定義(`(𝔪^{n+1}).inertia`)であり、`α` は
   `lowerRamificationGroup_eq_of_ramIndex_pow` の中で
   `mem_pow_maximalIdeal_iff_of_adjoin_eq_top` を通してしか使われないので、
   素元の選び方に依らない。本ファイルは `α = torsionGen`(= `x` 自身)で帳簿を作り、
   消費側には任意の素元 `α` で渡している。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ### §0 抽象核

★分岐・付値・Lubin-Tate・Galois の語彙が 1 つも出てこない 3 本。 -/

/-- **ℕ 上の「1 段ずつ降りる」述語は下向きに閉じている**。

`P (i+1) → P i` を仮定すると `j ≤ k → P k → P j`。★`Nat.le_induction` を使わずに
`k` の素朴な帰納法で書いてある(`j` を固定しないので `omega` で場合分けできる)。 -/
theorem le_of_antitone_natPred (P : ℕ → Prop) (hanti : ∀ i : ℕ, P (i + 1) → P i) :
    ∀ j k : ℕ, j ≤ k → P k → P j := by
  intro j k
  induction k with
  | zero => intro h hP; exact (show j = 0 by omega) ▸ hP
  | succ n ih =>
      intro h hP
      rcases Nat.lt_or_ge j (n + 1) with hj | hj
      · exact ih (by omega) (hanti n hP)
      · exact (show j = n + 1 by omega) ▸ hP

/-- ★★**「フィルターの段が一意に決まる」の抽象核**。

`P` が 1 段ずつ降りる述語で `P 0` が真、`P m` が偽ならば、
**`P j ↔ j ≤ i₀` となる `i₀ < m` が存在する**。

★`∃ i, P i ∧ ¬ P (i+1)`(切れ目の存在)ではなく**同値の形**で出しているのが要点。
消費側は「`u ∈ 1+𝔭^{j}` ⟺ `j ≤ i₀`」をそのまま使う。
★`i₀` はこの性質で一意に決まる(2 つあれば互いに `≤` になるので等しい)。 -/
theorem exists_le_iff_of_antitone_natPred (P : ℕ → Prop)
    (hanti : ∀ i : ℕ, P (i + 1) → P i) (h0 : P 0) (m : ℕ) (hm : ¬ P m) :
    ∃ i₀ : ℕ, i₀ < m ∧ ∀ j : ℕ, (P j ↔ j ≤ i₀) := by
  have hmono := le_of_antitone_natPred P hanti
  induction m with
  | zero => exact absurd h0 hm
  | succ n ih =>
      by_cases hn : P n
      · refine ⟨n, Nat.lt_succ_self n, fun j => ⟨fun hPj => ?_, fun hj => hmono j n hj hn⟩⟩
        by_contra hc
        exact hm (hmono (n + 1) j (by omega) hPj)
      · obtain ⟨i₀, hi₀, hiff⟩ := ih hn
        exact ⟨i₀, by omega, hiff⟩

/-- ★**商群での所属は `N ≤ H` のとき元の所属と同値**: `ū ∈ H.map (mk' N) ↔ u ∈ H`。

`Subgroup.comap_map_eq_self`(`ker f ≤ H` なら `comap f (map f H) = H`)と
`QuotientGroup.ker_mk'` だけ。

★`N ≤ H` は落とせない(`H = ⊥`, `N ≠ ⊥` なら左辺は真、右辺は偽)。
一般形は mathlib に在る —— `QuotientGroup.comap_map_mk' : comap (mk' N) (map (mk' N) H) = N ⊔ H`
(`GroupTheory/QuotientGroup/Basic.lean:364`)。★本補題はその `N ≤ H` の場合の
所属版で、`exact?` は当てられなかったので 1 行で書いてある。 -/
theorem mem_map_mk'_iff_of_le {G : Type*} [Group G] (N : Subgroup G) [N.Normal]
    (H : Subgroup G) (hle : N ≤ H) (u : G) :
    ((u : G ⧸ N)) ∈ H.map (QuotientGroup.mk' N) ↔ u ∈ H := by
  rw [show ((u : G ⧸ N)) = QuotientGroup.mk' N u from rfl, ← Subgroup.mem_comap,
    Subgroup.comap_map_eq_self (f := (QuotientGroup.mk' N)) (by rwa [QuotientGroup.ker_mk'])]

/-! ### §1 主単数のフィルターの底 -/

/-- `1 + 𝔭^0 = 𝒪_K^×`——`𝔭^0 = 𝒪_K` なので `u = 1 + (u−1)·1`。

★これが上の抽象核の `P 0` にあたる。★`u` が単数であることしか使っていない。 -/
theorem mem_principalUnits_zero {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) (u : (𝒪[K.carrier])ˣ) : u ∈ principalUnits K π 0 := by
  rw [mem_principalUnits_iff]
  exact ⟨(u : 𝒪[K.carrier]) - 1, by ring⟩

/-! ### §2 具体層 —— 段 1 の二分律と、下付き分岐群の帳簿 -/

section Bookkeeping

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-- ★★★★**`σ` ごとの二分律**(原典の「For σ ≠ id, ... If v_K(u−1) = i for 0 ≤ i < m」)。

`ρ_{f,m}(σ)` の代表 `u ∈ 𝒪_K^×` を取ると、次のいずれかが成り立つ:

* `σ = 1`(`u ∈ 1+𝔭^{M+1}`、すなわち `ρ(σ) = 1` の場合。`ρ` が同型だから);
* **一意な `i₀ ≤ M`** があって
  * `i(σ) = q^{i₀}`(段 1、`ramIndex_torsionGen_eq_pow_of_principalUnits`)、かつ
  * `ρ(σ) ∈ (1+𝔭^j)/(1+𝔭^{M+1}) ⟺ j ≤ i₀`(`j ≤ M+1`)。

★`i₀` の存在と一意性は抽象核 `exists_le_iff_of_antitone_natPred` から出る
(`P j := u ∈ 1+𝔭^j` は `principalUnits_antitone` で 1 段ずつ降り、
`P 0` は `mem_principalUnits_zero`、`P (M+1)` は場合分けの仮定で偽)。
★原典の `v_K(u−1) = i` に対応するが、**ℕ の引き算も除算も出てこない**。

★逸脱: 第 2 主張は `j ≤ M+1` に制限してある(冒頭「逸脱の記録 2」)。 -/
theorem ramIndex_torsionGen_dichotomy
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (M : ℕ) (hn : 1 ≤ M + 1) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (M + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    σ = 1 ∨ ∃ i₀ : ℕ, i₀ ≤ M ∧
      ramIndex (torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1) x hxn hmem) σ
        = (((pp ^ ff) ^ i₀ : ℕ) : ℕ∞) ∧
      ∀ j : ℕ, j ≤ M + 1 →
        (galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x hxψ hxn hmem σ
          ∈ (principalUnits K π j).map (QuotientGroup.mk' (principalUnits K π (M + 1)))
        ↔ j ≤ i₀) := by
  obtain ⟨u, hu⟩ := QuotientGroup.mk_surjective (s := principalUnits K π (M + 1))
    (galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x hxψ hxn hmem σ)
  by_cases hmM : u ∈ principalUnits K π (M + 1)
  · left
    have h1 : galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
        hxψ hxn hmem σ = 1 := by
      rw [← hu]; exact (QuotientGroup.eq_one_iff u).mpr hmM
    exact (MulEquiv.map_eq_one_iff _).mp h1
  · right
    obtain ⟨i₀, hi₀, hiff⟩ := exists_le_iff_of_antitone_natPred
      (fun j => u ∈ principalUnits K π j)
      (fun i hi => principalUnits_antitone K π (Nat.le_succ i) hi)
      (mem_principalUnits_zero K π u) (M + 1) hmM
    refine ⟨i₀, by omega, ?_, ?_⟩
    · obtain ⟨k, hk⟩ : ∃ k, M = i₀ + k := ⟨M - i₀, by omega⟩
      subst hk
      exact ramIndex_torsionGen_eq_pow_of_principalUnits K hq hπmax hπne0 f hf0 hf1 hf i₀ k x
        hxψ hxn hmem σ u hu.symm ((hiff i₀).mpr le_rfl)
        (fun hc => absurd ((hiff (i₀ + 1)).mp hc) (by omega))
    · intro j hj
      rw [← hu, mem_map_mk'_iff_of_le _ _ (principalUnits_antitone K π hj) u]
      exact hiff j

/-- ★★★★★★★★**原典の `G_n = ρ^{-1}_{f,m}(1 + p^i)`**(添字を 0 始まりにずらした形)。

    `q^i ≤ n < q^{i+1}`,  `i ≤ M`,  `m = M+1`  ⟹
      `G_n = ρ^{-1}_{f,m}((1+𝔭^{i+1})/(1+𝔭^{m}))`.

原典の該当箇所(1 行):

    Thus for G = Gal(K^m_x/L) and 1 ≤ i ≤ m, we have |G_n| = |ρ^{-1}_{f,m}(1+p^i)| = q^{m−i}
    for q^{i−1} − 1 < n ≤ q^i − 1.

★段取りは B4 の抽象核 `lowerRamificationGroup_eq_of_ramIndex_pow` に

* `hpow`: `i(σ)` は `⊤` か `q` 冪しか取らない —— 二分律の前半;
* `hV`: `σ ∈ ρ^{-1}(...) ⟺ q^{i+1} ≤ i(σ)` —— 二分律の後半 +
  `q^{i+1} ≤ q^{i₀} ⟺ i+1 ≤ i₀`(`Nat.pow_le_pow_iff_right`、`2 ≤ q` が要る)

を流すだけ。★素元 `α` には `x` 自身(`torsionGen`)を取り、それが素元であることは
`irreducible_torsionGen` から出す(仮定ではない)。
★`Algebra.adjoin 𝒪_K {α} = ⊤` は `adjoin_uniformizer_eq_top_adjoinIntegers`
(＝ #59 の定型 (a): **整数側 `adjoinIntegers K x` に寄せる**)。 -/
theorem lowerRamificationGroupAdjoin_eq_comap_map_principalUnits
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (M : ℕ) (hn : 1 ≤ M + 1) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (M + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (i n : ℕ) (hiM : i ≤ M) (hle : (pp ^ ff) ^ i ≤ n) (hlt : n < (pp ^ ff) ^ (i + 1)) :
    lowerRamificationGroupAdjoin K x n
      = Subgroup.comap (galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
          hxψ hxn hmem).toMonoidHom
        ((principalUnits K π (i + 1)).map
          (QuotientGroup.mk' (principalUnits K π (M + 1)))) := by
  have hq2 : 2 ≤ pp ^ ff := by rw [← hq]; exact Fintype.one_lt_card
  have hirr : Irreducible (torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1) x hxn hmem) :=
    irreducible_torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x hxψ hxn hmem
  have huni : IsLocalRing.maximalIdeal (adjoinIntegers K x)
      = Ideal.span {torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1) x hxn hmem} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer _).mp hirr
  have ht : IsTotallyRamifiedAdjoin K x :=
    isTotallyRamifiedAdjoin_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
      hxψ hxn hmem
  have hadj : Algebra.adjoin 𝒪[K.carrier]
      ({torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1) x hxn hmem} :
        Set (adjoinIntegers K x)) = ⊤ :=
    adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni
  have hdich := ramIndex_torsionGen_dichotomy K hq hπmax hπne0 f hf0 hf1 hf M hn x hxψ hxn hmem
  show lowerRamificationGroup (adjoinIntegers K x) _ n = _
  refine lowerRamificationGroup_eq_of_ramIndex_pow (A := 𝒪[K.carrier]) huni hadj (pp ^ ff) hq2
    ?_ _ hle hlt ?_
  · intro σ
    rcases hdich σ with h | ⟨i₀, -, hram, -⟩
    · left; rw [h]; exact ramIndex_one _
    · right; exact ⟨i₀, hram⟩
  · intro σ
    rw [Subgroup.mem_comap]
    rcases hdich σ with h | ⟨i₀, hi₀, hram, hmemiff⟩
    · rw [h, ramIndex_one]
      simp
    · rw [hram, Nat.cast_le, Nat.pow_le_pow_iff_right hq2, MulEquiv.coe_toMonoidHom]
      exact hmemiff (i + 1) (by omega)

/-! ### §3 `hρ` を外す -/

def upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_uniformizer.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Proposition 6.14", sectionId := "prop-6-14" }

/-- ★★★★★★★★★★★★★★★★**Yoshida 2008 Proposition 6.14(`n = 1` の場合)——
`hρ` を外した形**。

原文 (Yoshida08 p.17):
> Proposition 6.14. Let x ∈ K^× with v(x) = n > 0. Let L = K_n and K^m_x as in Definition 5.3. Then we have Gal(K^m_x/L)^m = {id} for all m ≥ 1 (see Definition 6.12).

`L = K_1 = K`、`K^m_x = K(x)`(`x` は `ψ_{M+1}` の根)、`m = M+1 ≥ 1`。

★★`Found/PGC/LubinTateUpperRamificationVanish.lean` の
`upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap` が残していた仮定 `hρ`
(＝原典の段 1)を、本ファイルの
`lowerRamificationGroupAdjoin_eq_comap_map_principalUnits` が供給する。
★**残る仮定は原典の設定と `huni`(`α` が `𝒪_{K(x)}` の素元)だけ**であり、
段 1・段 2 に関する仮定は 1 つも無い。

★逸脱: `n = 1` の場合のみ(決定 D29)。★消費側の `hρ` は `i + j = M` の形で
`i` を取るが、本ファイルの帳簿は `i ≤ M` しか要らないので `(by omega)` で渡している
(冒頭「逸脱の記録 3」)。 -/
theorem upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_uniformizer
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])] {ff : ℕ}
    (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (M : ℕ) (hn : 1 ≤ M + 1) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (M + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Fintype ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))]
    {α : adjoinIntegers K x}
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α}) :
    upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      α ((M + 1 : ℕ) : ℝ) = ⊥ :=
  upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap K hq hπmax hπne0 f hf0 hf1 hf M hn x
    hxψ hxn hmem huni
    (fun i j n hij hle hlt =>
      lowerRamificationGroupAdjoin_eq_comap_map_principalUnits K hq hπmax hπne0 f hf0 hf1 hf M hn x
        hxψ hxn hmem i n (by omega) hle hlt)

/-- ★★★★**Prop 6.14(`n = 1`)——仮定を 1 つも足さない形**。

素元 `α` として捩れ点 `x` 自身(`torsionGen`)を取れば、`huni` すら要らない
(`irreducible_torsionGen` が `x` を素元だと言う)。

★★**残っているのは原典の設定(`f` が Lubin-Tate 級数、`x` が `ψ_{M+1}` の原始根)
だけである。** -/
theorem upperRamificationGroup_torsionGen_eq_bot
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])] {ff : ℕ}
    (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (M : ℕ) (hn : 1 ≤ M + 1) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (M + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Fintype ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))] :
    upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      (torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1) x hxn hmem) ((M + 1 : ℕ) : ℝ) = ⊥ :=
  upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_uniformizer K hq hπmax hπne0 f hf0 hf1 hf
    M hn x hxψ hxn hmem
    ((IsDiscreteValuationRing.irreducible_iff_uniformizer _).mp
      (irreducible_torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x hxψ hxn hmem))

end Bookkeeping

end ABC3.Found.PGC
