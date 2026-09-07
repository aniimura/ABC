import ABC3.Found.PGC.UpperRamificationGroup
import ABC3.Found.PGC.TorsionPointCriterion
import ABC3.Found.PGC.LubinTateTotallyRamified

/-!
# `Gal(K(µ_{f,m})/K)^m = {id}` —— Yoshida 2008 Proposition 6.14 の `n = 1` 形

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Proposition 6.14(物理 p.17)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-6.html` の `#prop-6-14`。

原文 (Yoshida08 p.17):
> Proposition 6.14. Let x ∈ K^× with v(x) = n > 0. Let L = K_n and K^m_x as in Definition 5.3. Then we have Gal(K^m_x/L)^m = {id} for all m ≥ 1 (see Definition 6.12).

★★★**本ノードは Prop 6.14 の `n = 1` の場合だけを扱う。**
決定 D29 により、Local Kronecker-Weber(Theorem 6.15)にはこれで足りる ——
Thm 6.15 の逐語は「Take a σ ∈ W(K^LT/K) with v(σ) = n > 0」で始まるので `σ` は
こちらで選んでよく、`n = 1` を取れる。`n = 1` では `L = K_1 = K` となって相対
Lubin-Tate が絶対 Lubin-Tate に潰れ、`K^m_x = K(µ_{f,m})` が木の
`IntermediateField.adjoin K.carrier {x}`(`x` は `ψ_m` の根、すなわち原始的な
`π^m`-捩れ点)そのものになる。★**一般の `n` については何も主張していない。**

## 原典の証明(全 8 行。★`>` を付けない —— 逐語引用は上の 1 本だけである)

    Proof. Let K^m_x = L^m_f and α ∈ µ^×_{f,m}. For σ ∈ Gal(K^m_x/L)\{id}, we have
    i(σ) = v(σ(α)−α) by Proposition 4.4(ii), where v = v_{K^m_x}. If ρ_{f,m}(σ) = u mod p^m
    ∈ (O/p^m)^× (see Proposition 4.4(iii)), then σ(α) = [u]_f(α). For σ ≠ id, set
    β := [u−1]_f(α). If v_K(u−1) = i for 0 ≤ i < m, then β ∈ µ^×_{f,m−i} by Lemma 4.3(ii).
    Hence β is a uniformizer of K^{m−i}_x by Proposition 4.4(ii), which shows v(β) = q^i.
    Now σ(α) = [u]_f(α) = α +_f β ≡ α+β (mod αβ), hence i(σ) = v(σ(α)−α) = v(β) = q^i.
    Thus for G = Gal(K^m_x/L) and 1 ≤ i ≤ m, we have |G_n| = |ρ^{-1}_{f,m}(1+p^i)| = q^{m−i}
    for q^{i−1}−1 < n ≤ q^i−1. Thus φ_G(q^m−1) = (1/|G|) Σ_{i=1}^{q^m−1} |G_i|
    = (1/((q−1)q^{m−1})) (Σ_{i=1}^{m}(q^i−q^{i−1})·q^{m−i}) = m and G^m = G_{q^m−1} = {id}.

## 本ファイルが入れたもの(4 段のうち 段 2・段 3・段 4)

原典の証明は 4 段に分かれる。本ファイルは **段 2 の数え上げ・段 3 の算術・段 4 の配管**を
入れ、**段 1(`i(σ) = q^i`)だけを仮定として残した**。

| 段 | 内容 | 本ファイル |
|---|---|---|
| 1 | `i(σ) = q^{v_K(u−1)}`(`ρ` で `σ` を単数に翻訳し、`v(σα−α)` を割り算に) | ★**仮定 `hρ`。未着手** |
| 2 | `|ρ^{-1}(1+p^i)| = q^{m−i}` の数え上げ | ✓ `natCard_map_principalUnits` |
| 3 | `φ_G(q^m−1) = m` の算術 | ✓ `sum_Ico_pow_blocks` / `sum_Icc_pow_blocks` |
| 4 | `G^m = G_{q^m−1} = {id}`(Def 6.12 で上付きに移す) | ✓ `upperRamificationGroup_eq_bot_of_natCard_lowerRamificationGroup` |

主定理は 2 本ある。

* `upperRamificationGroup_iteratedLubinTatePsi_eq_bot`
  —— 段 1・段 2 を **`|G_n|` の値そのもの**(`hlower`)として受け取る形。
* `upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap`
  —— 段 1 だけを **`G_n = ρ^{-1}(1+p^{i+1})`**(`hρ`)として受け取る形。
  段 2 は本ファイルが `natCard_map_principalUnits` で消費している。★**こちらが本命**で、
  Prop 6.14(`n = 1`)に残っているのは `hρ` ただ 1 本である。

## 抽象核(§0–§2。分岐・付値・Lubin-Tate の語彙が 1 つも出てこない)

★原典の段 3 は「`q` 進のブロック分割の和」という**純算術**である:

    `∑_{n ∈ [1, q^{m}−1]} c(n) = m·(q^m − q^{m−1})`   ただし `c` は
    ブロック `[q^i, q^{i+1})` の上で定数 `q^{m−1−i}`

これが `sum_Ico_pow_blocks`(`K` 倍の重み付きで帰納法を回す形)である。
★**重み `K` を入れたのが鍵**: 素朴に `M` について帰納法を回すと、下半分
`[1, q^M)` に現れる値が `q^{j+1}` であって帰納法の仮定の `q^j` と食い違う。
`K := K·q` と持ち替えることで食い違いが消える(切り詰め引き算も除算も出てこない)。

★段 2 の数え上げも**純群論**に落ちる: `natCard_subgroup_of_natCard_quotient`
(`|Q| = a·b` かつ `|Q/V| = b` なら `|V| = a`、Lagrange だけ)。

★段 1 と段 2 のあいだの「ブロックの読み替え」も抽象核に切り出した:
`lowerRamificationGroup_eq_of_ramIndex_pow`。`i(σ)` が `q` 冪(または `⊤`)しか
取らないなら、`q^i ≤ n < q^{i+1}` の範囲で `G_n` は `n` に依らず
`{σ | q^{i+1} ≤ i(σ)}` に等しい、という主張である。

## 逸脱の記録

1. ★★**`n = 1` に固定した**(上記 D29)。原典の `L = K_n`・`K^m_x` は本ファイルでは
   `L = K`・`K^m_x = K(x)`(`x` は `ψ_m` の根)である。
2. ★**段 1 を仮定に置いた**(`hlower` / `hρ`)。原典は証明している。
   `hρ` を落とすと主張は**空虚ではなく偽になり得る**ので、消費側は必ず段 1 を供給すること。
3. ★`m ≥ 1` は `m = M + 1` と書く形で埋め込んである(★`m = 0` では `µ_{f,0} = {0}` で
   拡大が自明になり、主張が空虚になる)。
4. ★上付きと下付きを取り違えないこと。★**上付きは番号が大きいほど小さい**。
   `G^m = ⊥` は「`m` で消える」であって `G_m = ⊥` ではない
   (実際 `G_m ≠ ⊥` である ―― 消えるのは `G_{q^m−1}`)。
5. ★`Fintype (Gal(K(x)/K))` と `IsDiscreteValuationRing (adjoinIntegers K x)` は
   木ではインスタンスになっていないので、前者は束縛子で受け取り、後者は
   `attribute [local instance]`(`ConjugateSumValuation.lean` と同じ扱い)で入れている。

## 測ったこと(在庫調査係が「測れなかった」とした 3 点)

1. **Lemma 6.10(i) の在庫**: `herbrandPhiGroup_natCast` は**そのまま消費できた**
   (`(A := 𝒪[K.carrier])` を明示するだけ)。
2. **`ρ^{-1}(1+p^i)` の位数**: ★**素直に出た**。第三同型定理
   (`QuotientGroup.quotientQuotientEquivQuotient`)と既存の
   `card_principalUnitsQuotient`(`|𝒪^×/(1+p^n)| = q^n − q^{n−1}`)を Lagrange で
   割るだけで `q^{m−i}` になる(`natCard_map_principalUnits`)。★新しい数学は要らなかった。
3. **`µ^×_{f,m−i}` の非退化側**: 木には**原始捩れ点の集合そのもの**
   (`iteratedLubinTatePsiTorsionPoints` ＝ `ψ_n` の根)があり、
   `lubinTateActionAtTorsionPoint_eq_zero_iff_dvd_of_mem_iteratedLubinTatePsiTorsionPoints`
   (`[a](α) = 0 ↔ π^n ∣ a`、★**両向き**)と
   `iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints`
   (`Λ_n \ Λ_{n−1} = ψ_n の根`)がある。★**「`v_K(a) = i` なら `[a](α) ∈ µ^×_{f,m−i}`」
   そのものは無い**が、上の 2 本から短く出るはずである(未着手)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open IsLocalRing IsDiscreteValuationRing
open scoped NormedField Valued Classical

/-! ## §0 抽象核(純算術)—— `q` 進ブロックの和

★分岐・付値・Galois の語彙が 1 つも出てこない。 -/

/-- ★★**段 3 の抽象核** —— `c` が `[q^i, q^{i+1})` の上で定数 `K·q^j`(`i + j = M`)なら

    `∑_{n ∈ [1, q^{M+1})} c(n) = K·(M+1)·(q^{M+1} − q^M)`

(`q = r+1`、`q^{M+1} − q^M = r·q^M` と書いて**切り詰め引き算を消してある**)。

★**重み `K` を入れたのが鍵である。** `K` 無しで `M` について帰納法を回すと、
下半分 `[1, q^M)` に現れる値が帰納法の仮定の `q^j` ではなく `q^{j+1}` になって
食い違う。`K := K·q` と持ち替えれば食い違いが消える。 -/
theorem sum_Ico_pow_blocks (r : ℕ) :
    ∀ (M K : ℕ) (c : ℕ → ℕ),
      (∀ i j : ℕ, i + j = M → ∀ n ∈ Finset.Ico ((r + 1) ^ i) ((r + 1) ^ (i + 1)),
          c n = K * (r + 1) ^ j) →
      ∑ n ∈ Finset.Ico 1 ((r + 1) ^ (M + 1)), c n = K * ((M + 1) * (r * (r + 1) ^ M)) := by
  intro M
  induction M with
  | zero =>
      intro K c hc
      have hval : ∀ n ∈ Finset.Ico 1 ((r + 1) ^ 1), c n = K := by
        intro n hn
        have := hc 0 0 rfl n (by simpa using hn)
        simpa using this
      rw [Finset.sum_congr rfl hval]
      simp [Nat.mul_comm]
  | succ M ih =>
      intro K c hc
      have h1 : 1 ≤ (r + 1) ^ (M + 1) := Nat.one_le_pow _ _ (Nat.succ_pos r)
      have h2 : (r + 1) ^ (M + 1) ≤ (r + 1) ^ (M + 1 + 1) :=
        Nat.pow_le_pow_right (Nat.succ_pos r) (by omega)
      have hsplit := Finset.sum_Ico_consecutive c h1 h2
      have hbot : ∑ n ∈ Finset.Ico 1 ((r + 1) ^ (M + 1)), c n
          = (K * (r + 1)) * ((M + 1) * (r * (r + 1) ^ M)) := by
        refine ih (K * (r + 1)) c ?_
        intro i j hij n hn
        have := hc i (j + 1) (by omega) n hn
        rw [this]; ring
      have htopval : ∀ n ∈ Finset.Ico ((r + 1) ^ (M + 1)) ((r + 1) ^ (M + 1 + 1)), c n = K := by
        intro n hn
        have := hc (M + 1) 0 (by omega) n hn
        simpa using this
      have hcard : (r + 1) ^ (M + 1 + 1) - (r + 1) ^ (M + 1) = r * (r + 1) ^ (M + 1) := by
        have hpow : (r + 1) ^ (M + 1 + 1) = r * (r + 1) ^ (M + 1) + (r + 1) ^ (M + 1) := by ring
        omega
      have htop : ∑ n ∈ Finset.Ico ((r + 1) ^ (M + 1)) ((r + 1) ^ (M + 1 + 1)), c n
          = (r * (r + 1) ^ (M + 1)) * K := by
        rw [Finset.sum_congr rfl htopval, Finset.sum_const, Nat.card_Ico, smul_eq_mul, hcard]
      rw [← hsplit, hbot, htop]
      ring

/-- ★**原典の `Σ_{i=1}^{q^m−1}` の形**(`K = 1`、区間が `Icc 1 (q^m − 1)`)。

`Finset.Icc 1 (q^{M+1} − 1) = Finset.Ico 1 (q^{M+1})` と読み替えて
`sum_Ico_pow_blocks` に渡すだけ。 -/
theorem sum_Icc_pow_blocks (r M : ℕ) (c : ℕ → ℕ)
    (hc : ∀ i j n : ℕ, i + j = M → (r + 1) ^ i ≤ n → n < (r + 1) ^ (i + 1) →
      c n = (r + 1) ^ j) :
    ∑ n ∈ Finset.Icc 1 ((r + 1) ^ (M + 1) - 1), c n = (M + 1) * (r * (r + 1) ^ M) := by
  have h1 : 1 ≤ (r + 1) ^ (M + 1) := Nat.one_le_pow _ _ (Nat.succ_pos r)
  have hset : Finset.Icc 1 ((r + 1) ^ (M + 1) - 1) = Finset.Ico 1 ((r + 1) ^ (M + 1)) := by
    ext n
    simp only [Finset.mem_Icc, Finset.mem_Ico]
    omega
  rw [hset]
  have := sum_Ico_pow_blocks r M 1 c (by
    intro i j hij n hn
    rw [Finset.mem_Ico] at hn
    rw [hc i j n hij hn.1 hn.2, one_mul])
  simpa using this

/-! ## §1 抽象核(純群論)

★段 2 の「`|ρ^{-1}(1+p^i)| = q^{m−i}`」は、群同型で移した先で
**Lagrange だけ**になる。 -/

/-- ★**段 2 の抽象核** —— `|Q| = a·b` かつ `|Q/V| = b`(`b > 0`)なら `|V| = a`。 -/
theorem natCard_subgroup_of_natCard_quotient {Q : Type*} [Group Q] (V : Subgroup Q)
    {a b : ℕ} (hQ : Nat.card Q = a * b) (hquot : Nat.card (Q ⧸ V) = b) (hb : 0 < b) :
    Nat.card V = a := by
  have h := Subgroup.card_eq_card_quotient_mul_card_subgroup V
  rw [hQ, hquot] at h
  exact (Nat.eq_of_mul_eq_mul_left hb (by rw [← h]; ring)).symm

/-- ★**群同型による引き戻しは位数を変えない**(原典の `ρ^{-1}` の位数を
値域側で数えるための橋)。 -/
theorem natCard_comap_mulEquiv {G Q : Type*} [Group G] [Group Q] (e : G ≃* Q) (V : Subgroup Q) :
    Nat.card (V.comap e.toMonoidHom) = Nat.card V :=
  Nat.card_congr (Equiv.subtypeEquiv e.toEquiv (by intro a; simp))

/-! ## §2 抽象核(分岐群)—— 段 3・段 4

★Lubin-Tate も付値も出てこない。任意の離散付値環 `B` と有限群 `G` の話である。 -/

section AbstractCore

variable {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
  [IsDiscreteValuationRing B]
variable {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]

omit [Fintype G] in
/-- ★★**段 1 と段 2 のあいだの「ブロックの読み替え」**。

`i(σ)` が `q` 冪(または `⊤`)しか取らないとき、`q^i ≤ n < q^{i+1}` の範囲で
下付き分岐群 `G_n` は `n` に依らず `V = {σ | q^{i+1} ≤ i(σ)}` に等しい。

★原典の「`q^{i−1}−1 < n ≤ q^i−1` の範囲で `|G_n| = |ρ^{-1}(1+p^i)|`」の、
**分岐の語彙を使わない中身**である(添字は 0 始まりにずらしてある: 原典の `i` は
ここでの `i+1`)。 -/
theorem lowerRamificationGroup_eq_of_ramIndex_pow
    {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (q : ℕ) (hq : 2 ≤ q)
    (hpow : ∀ σ : G, ramIndex α σ = ⊤ ∨ ∃ k : ℕ, ramIndex α σ = ((q ^ k : ℕ) : ℕ∞))
    (V : Subgroup G) {i n : ℕ} (hle : q ^ i ≤ n) (hlt : n < q ^ (i + 1))
    (hV : ∀ σ : G, σ ∈ V ↔ ((q ^ (i + 1) : ℕ) : ℕ∞) ≤ ramIndex α σ) :
    lowerRamificationGroup B G n = V := by
  ext σ
  rw [mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj, hV σ]
  rcases hpow σ with h | ⟨k, h⟩
  · rw [h]
    exact iff_of_true (lt_of_le_of_ne le_top (ENat.coe_ne_top n)) le_top
  · rw [h]
    simp only [Nat.cast_lt, Nat.cast_le]
    constructor
    · intro hn
      have hik : i + 1 ≤ k := by
        by_contra hc
        have hki : k ≤ i := by omega
        have hkle : q ^ k ≤ q ^ i := Nat.pow_le_pow_right (by omega) hki
        omega
      exact Nat.pow_le_pow_right (by omega) hik
    · intro hk
      omega

/-- ★★★**段 3 + 段 4 の抽象核** —— 原典の最後の 2 行そのもの。

`|G| = q^{m} − q^{m−1}`(`q = r+1`, `m = M+1` と書いて `r·q^M`)で、
下付き分岐群の位数がブロックごとに `|G_n| = q^{j}`(`q^i ≤ n < q^{i+1}`, `i + j = M`)
なら、**上付き分岐群 `G^{M+1}` は自明**である。

★段取りは原典どおり:

1. `φ_G(q^{M+1}−1) = (1/|G|)·Σ_{n=1}^{q^{M+1}−1} |G_n|`(Lemma 6.10(i)、在庫
   `herbrandPhiGroup_natCast`)。
2. 分子は `sum_Icc_pow_blocks` で `(M+1)·(q^{M+1} − q^M)`、分母は `|G|` で同じ値
   ——★**割ってちょうど `M+1`**。
3. `G^{φ_G(n)} = G_n`(在庫 `upperRamificationGroup_herbrandPhiGroup`)で上付きに移し、
   `n = q^{M+1}−1` は最上段のブロック(`i = M`, `j = 0`)なので `|G_n| = q^0 = 1`、
   すなわち `G_n = ⊥`。

★★**`r > 0`(＝ `q ≥ 2`)を仮定に置いていない** —— `Nat.card G > 0` から出るので、
原典が黙って使っている「剰余体は 2 元以上」を持ち込まずに済んでいる。 -/
theorem upperRamificationGroup_eq_bot_of_natCard_lowerRamificationGroup
    {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (r M : ℕ)
    (hGcard : Nat.card G = r * (r + 1) ^ M)
    (hlower : ∀ i j n : ℕ, i + j = M → (r + 1) ^ i ≤ n → n < (r + 1) ^ (i + 1) →
      Nat.card (lowerRamificationGroup B G n) = (r + 1) ^ j) :
    upperRamificationGroup G α ((M + 1 : ℕ) : ℝ) = ⊥ := by
  have hpow1 : 1 ≤ (r + 1) ^ (M + 1) := Nat.one_le_pow _ _ (Nat.succ_pos r)
  have hpowM : 1 ≤ (r + 1) ^ M := Nat.one_le_pow _ _ (Nat.succ_pos r)
  have hGpos : 0 < Nat.card G := Nat.card_pos
  have hrpos : 0 < r := by
    rcases Nat.eq_zero_or_pos r with h | h
    · rw [h, Nat.zero_mul] at hGcard; omega
    · exact h
  have hden : 0 < r * (r + 1) ^ M := Nat.mul_pos hrpos hpowM
  -- Lemma 6.10(i) の分子
  have hsum : ∑ n ∈ Finset.Icc 1 ((r + 1) ^ (M + 1) - 1),
      Nat.card (lowerRamificationGroup B G n) = (M + 1) * (r * (r + 1) ^ M) :=
    sum_Icc_pow_blocks r M _ hlower
  -- `φ_G(q^m − 1) = m`
  have hphi : herbrandPhiGroup G α (((r + 1) ^ (M + 1) - 1 : ℕ) : ℝ) = ((M + 1 : ℕ) : ℝ) := by
    rw [herbrandPhiGroup_natCast (A := A) huni hadj]
    have hcast : (∑ n ∈ Finset.Icc 1 ((r + 1) ^ (M + 1) - 1),
        (Nat.card (lowerRamificationGroup B G n) : ℝ))
        = (((M + 1) * (r * (r + 1) ^ M) : ℕ) : ℝ) := by
      rw [← hsum]
      push_cast
      rfl
    rw [hcast, hGcard]
    push_cast
    field_simp
  -- 最上段のブロックでは `|G_n| = q^0 = 1`
  have hbot : lowerRamificationGroup B G ((r + 1) ^ (M + 1) - 1) = ⊥ := by
    have hle : (r + 1) ^ M ≤ (r + 1) ^ (M + 1) - 1 := by
      have hh : (r + 1) ^ (M + 1) = r * (r + 1) ^ M + (r + 1) ^ M := by ring
      have : 1 ≤ r * (r + 1) ^ M := hden
      omega
    have hlt : (r + 1) ^ (M + 1) - 1 < (r + 1) ^ (M + 1) := by omega
    have := hlower M 0 ((r + 1) ^ (M + 1) - 1) (by omega) hle hlt
    simp only [pow_zero] at this
    exact Subgroup.card_eq_one.mp this
  calc upperRamificationGroup G α ((M + 1 : ℕ) : ℝ)
      = upperRamificationGroup G α
          (herbrandPhiGroup G α (((r + 1) ^ (M + 1) - 1 : ℕ) : ℝ)) := by rw [hphi]
    _ = ramificationGroupReal α (((r + 1) ^ (M + 1) - 1 : ℕ) : ℝ) :=
        upperRamificationGroup_herbrandPhiGroup _ _
    _ = lowerRamificationGroup B G ((r + 1) ^ (M + 1) - 1) :=
        ramificationGroupReal_natCast (A := A) huni hadj _
    _ = ⊥ := hbot

end AbstractCore

/-! ## §3 段 2 —— `|ρ^{-1}(1+p^i)| = q^{m−i}`

★原典が「`|G_n| = |ρ^{-1}_{f,m}(1+p^i)| = q^{m−i}`」と 1 行で済ませているところ。
`ρ` の値域 `(𝒪_K)^×/(1+p^m)` の側で数える。 -/

/-- 主単数の列は減少列: `i ≤ m ⟹ 1+p^m ⊆ 1+p^i`。 -/
theorem principalUnits_antitone {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (π : 𝒪[K.carrier]) {i m : ℕ} (h : i ≤ m) :
    principalUnits K π m ≤ principalUnits K π i := by
  intro v hv
  have hsub : Ideal.span ({π ^ m} : Set 𝒪[K.carrier])
      ≤ Ideal.span ({π ^ i} : Set 𝒪[K.carrier]) :=
    (Ideal.span_singleton_le_span_singleton).mpr (pow_dvd_pow π h)
  exact hsub hv

/-- `((𝒪_K)^×/(1+p^{i+1+j})) / ((1+p^{i+1})/(1+p^{i+1+j})) ≅ (𝒪_K)^×/(1+p^{i+1})`
の位数 —— 第三同型定理と既存の `card_principalUnitsQuotient` だけ。 -/
theorem natCard_quotient_map_principalUnits {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])] {ff : ℕ}
    (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0) (I j : ℕ) :
    Nat.card ((𝒪[K.carrier]ˣ ⧸ principalUnits K π (I + 1 + j)) ⧸
        (principalUnits K π (I + 1)).map
          (QuotientGroup.mk' (principalUnits K π (I + 1 + j))))
      = (pp ^ ff) ^ (I + 1) - (pp ^ ff) ^ I := by
  rw [Nat.card_congr (QuotientGroup.quotientQuotientEquivQuotient
    (principalUnits K π (I + 1 + j)) (principalUnits K π (I + 1))
    (principalUnits_antitone K π (by omega))).toEquiv]
  have h := card_principalUnitsQuotient K hq hπmax hπne0 (I + 1) (by omega)
  rw [Nat.add_sub_cancel] at h
  exact h

/-- ★★★**段 2 の本体** —— `(𝒪_K/p^{m})^×` のなかの `1+p^{i}` の像の位数は `q^{m−i}`。

`m = I+1+j`・`i = I+1` と書いてあるので、結論は `q^j`(＝ `q^{m−i}`)。

★段取り: `|Q| = q^m − q^{m−1}`(在庫 `card_principalUnitsQuotient`)、
`|Q/V| = q^i − q^{i−1}`(第三同型定理、`natCard_quotient_map_principalUnits`)、
両者を Lagrange(抽象核 `natCard_subgroup_of_natCard_quotient`)で割る。
★**新しい数学は要らなかった。** -/
theorem natCard_map_principalUnits {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])] {ff : ℕ}
    (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0) (I j : ℕ) :
    Nat.card ((principalUnits K π (I + 1)).map
        (QuotientGroup.mk' (principalUnits K π (I + 1 + j))))
      = (pp ^ ff) ^ j := by
  have hq2 : 2 ≤ pp ^ ff := by rw [← hq]; exact Fintype.one_lt_card
  obtain ⟨r, hr⟩ : ∃ r : ℕ, pp ^ ff = r + 1 := ⟨pp ^ ff - 1, by omega⟩
  have hrpos : 0 < r := by omega
  have hquot := natCard_quotient_map_principalUnits K hq hπmax hπne0 I j
  have hQ := card_principalUnitsQuotient K hq hπmax hπne0 (I + 1 + j) (by omega)
  rw [show I + 1 + j - 1 = I + j from by omega] at hQ
  rw [hr] at hquot hQ ⊢
  have hb : (r + 1) ^ (I + 1) - (r + 1) ^ I = r * (r + 1) ^ I := by
    have hpow : (r + 1) ^ (I + 1) = r * (r + 1) ^ I + (r + 1) ^ I := by ring
    omega
  have hbpos : 0 < r * (r + 1) ^ I :=
    Nat.mul_pos hrpos (Nat.pow_pos (Nat.succ_pos r))
  refine natCard_subgroup_of_natCard_quotient (a := (r + 1) ^ j) (b := r * (r + 1) ^ I) _
    ?_ ?_ hbpos
  · rw [hQ]
    have hpow : (r + 1) ^ (I + 1 + j) = r * (r + 1) ^ (I + j) + (r + 1) ^ (I + j) := by ring
    have hmul : (r + 1) ^ j * (r * (r + 1) ^ I) = r * (r + 1) ^ (I + j) := by ring
    omega
  · rw [hquot]
    exact hb

/-! ## §4 具体層 —— `K(µ_{f,m})` への代入

★`IsDiscreteValuationRing (adjoinIntegers K x)` は木でインスタンスになっていないので、
`ConjugateSumValuation.lean` と同じく局所インスタンスとして入れる。 -/

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-- ★★★★**Yoshida 2008 Proposition 6.14(`n = 1` の場合)** ——
`|G_n|` をブロックごとに与える形。

`x` を `ψ_{M+1}` の根(＝原始的な `π^{M+1}`-捩れ点、原典の `α ∈ µ^×_{f,m}`)、
`G = Gal(K(x)/K)`、`α` を `𝒪_{K(x)}` の一意化元とするとき、
下付き分岐群の位数が原典の値
`|G_n| = q^{j}`(`q^i ≤ n < q^{i+1}`, `i + j = M`)
であれば **`G^{M+1} = {id}`**。

★`|G| = q^{M+1} − q^M` は在庫から出る(`galoisReciprocityEquiv` +
`card_units_quotient_span_pi_pow`)ので**仮定に置いていない**。
★★原典が `[K^m_x : L] = (q−1)q^{m−1}` を Prop 4.4(ii) から引くところは、
木では `finrank_adjoin_iteratedLubinTatePsi` を経由せず
**Galois 群の位数そのもの**として取れる(`IsGalois` の instance 群を触らずに済む
——★これが「原典より安い道」)。

★逸脱: `hlower`(段 1・段 2)は仮定である。冒頭「逸脱の記録 2」。 -/
theorem upperRamificationGroup_iteratedLubinTatePsi_eq_bot
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
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α})
    (hlower : ∀ i j n : ℕ, i + j = M → (pp ^ ff) ^ i ≤ n → n < (pp ^ ff) ^ (i + 1) →
      Nat.card (lowerRamificationGroupAdjoin K x n) = (pp ^ ff) ^ j) :
    upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      α ((M + 1 : ℕ) : ℝ) = ⊥ := by
  have hqpos : 0 < pp ^ ff := by rw [← hq]; exact Fintype.card_pos
  obtain ⟨r, hr⟩ : ∃ r : ℕ, pp ^ ff = r + 1 := ⟨pp ^ ff - 1, by omega⟩
  have ht : IsTotallyRamifiedAdjoin K x :=
    isTotallyRamifiedAdjoin_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
      hxψ hxn hmem
  have hadj : Algebra.adjoin 𝒪[K.carrier] ({α} : Set (adjoinIntegers K x)) = ⊤ :=
    adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni
  have hGcard : Nat.card ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      = r * (r + 1) ^ M := by
    rw [Nat.card_congr (galoisReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
      hxψ hxn hmem).toEquiv]
    have h := card_units_quotient_span_pi_pow K hq hπmax hπne0 (M + 1) hn
    rw [Nat.add_sub_cancel, hr] at h
    rw [h]
    have hpow : (r + 1) ^ (M + 1) = r * (r + 1) ^ M + (r + 1) ^ M := by ring
    omega
  refine upperRamificationGroup_eq_bot_of_natCard_lowerRamificationGroup
    (A := 𝒪[K.carrier]) huni hadj r M hGcard ?_
  intro i j n hij hle hlt
  have := hlower i j n hij (by rw [hr]; exact hle) (by rw [hr]; exact hlt)
  rw [hr] at this
  exact this

/-! ## §5 組み立て —— 残っているのは段 1 ただ 1 本 -/

def upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Proposition 6.14", sectionId := "prop-6-14" }

/-- ★★★★★★★★**Yoshida 2008 Proposition 6.14(`n = 1` の場合)** ——
残る仮定は原典の段 1 だけ。

原文 (Yoshida08 p.17):
> Proposition 6.14. Let x ∈ K^× with v(x) = n > 0. Let L = K_n and K^m_x as in Definition 5.3. Then we have Gal(K^m_x/L)^m = {id} for all m ≥ 1 (see Definition 6.12).

★★**`n = 1` の場合のみ**(冒頭「逸脱の記録 1」、決定 D29)。
`L = K_1 = K`、`K^m_x = K(x)`(`x` は `ψ_{M+1}` の根)、`m = M+1 ≥ 1`。

仮定 `hρ` が原典の段 1、すなわち

    `G_n = ρ^{-1}_{f,m}(1 + p^{i+1})`   (`q^i ≤ n < q^{i+1}`, `i + j = M`)

である(原典は `i(σ) = v(σ(α)−α) = q^{v_K(u−1)}` を経由してこれを出す)。
★添字は 0 始まりにずらしてある: 原典の `1 ≤ i ≤ m` はここでの `i+1`。

★★**段 2(`|ρ^{-1}(1+p^{i+1})| = q^{j}`)は本ファイルが消費している**
(`natCard_map_principalUnits`)ので、`hρ` を供給すれば Prop 6.14(`n = 1`)は閉じる。

★逸脱: `hρ` は仮定である。冒頭「逸脱の記録 2」。 -/
theorem upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap
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
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α})
    (hρ : ∀ i j n : ℕ, i + j = M → (pp ^ ff) ^ i ≤ n → n < (pp ^ ff) ^ (i + 1) →
      lowerRamificationGroupAdjoin K x n
        = Subgroup.comap (galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
            hxψ hxn hmem).toMonoidHom
          ((principalUnits K π (i + 1)).map
            (QuotientGroup.mk' (principalUnits K π (M + 1))))) :
    upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      α ((M + 1 : ℕ) : ℝ) = ⊥ := by
  refine upperRamificationGroup_iteratedLubinTatePsi_eq_bot K hq hπmax hπne0 f hf0 hf1 hf M hn x
    hxψ hxn hmem huni ?_
  intro i j n hij hle hlt
  rw [hρ i j n hij hle hlt, natCard_comap_mulEquiv]
  have hidx : i + 1 + j = M + 1 := by omega
  have hcard := natCard_map_principalUnits K hq hπmax hπne0 i j
  rw [hidx] at hcard
  exact hcard

end ABC3.Found.PGC
