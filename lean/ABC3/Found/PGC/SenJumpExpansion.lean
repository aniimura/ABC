import ABC3.Found.PGC.SenValuationCongruence

/-!
# 跳び位置の `p` 進展開(Yoshida 2008 Corollary 6.7)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Corollary 6.7**(物理 p.15)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
`section-6.html` の `id="cor-6-7"`(`data-pdf-page="15"`, `data-item="Corollary 6.7"`)。

原文 (Yoshida08 p.15。`0_Source` の `.txt` 1022–1023 行がこの主張の**全文**である):
> Corollary 6.7. Assume G ∼= Z/pmZ. Then there exist n0, n1, . . . , nm−1 ∈Z≥1 such that,
> for 1 ≤j ≤m −1, we have |Gn| = pm−j if and only if Pj−1
> i=0 nipi < n ≤Pj
> i=0 nipi.

読み下すと(上付き・下付きを `^` / `_` に開き、`P` と潰れている総和記号を `Σ` に戻す):
`G ≅ ℤ/p^mℤ` のとき `n_0, …, n_{m−1} ∈ ℤ_{≥1}` があって、`1 ≤ j ≤ m−1` について
`|G_n| = p^{m−j} ⟺ Σ_{i=0}^{j−1} n_i p^i < n ≤ Σ_{i=0}^{j} n_i p^i`。
★`pdftotext` は `≅` を `∼=` に割り、`Σ` を `P` に潰す。

## ★★原典は Corollary 6.7 に証明を与えていない

`.txt` の 1023 行の直後(1024 行)は「6.2. The Hasse-Arf theorem.」で、**Proof 段落が無い**。
原典が読者に投げた箇所(省略の合図)であり、
★**以下の証明は我々の構成である**(原典に段取りは存在しない)。
隣の Proposition 6.6 には完全な証明があり(`.txt` 986–1021 行)、本ファイルはその (i)(iii) を
`Found/PGC/SenJumpFiltration.lean` / `Found/PGC/SenValuationCongruence.lean` 経由で使う。
★したがって「原典の証明と一致するか」は問えない。★**問えるのは
「原典が要求している `n_i ∈ ℤ_{≥1}` を満たす取り方になっているか」だけ**であり、
それは結論の `∀ i, 1 ≤ nn i` として Lean で保証してある。

## 本ファイルの範囲

Corollary 6.7 を**両向き**(`⟹` と `⟸`)で埋めた。入力は 3 本とも既に在庫にある:

* **(i) 前半**(狭義増加)`ramIndex_pow_pow_lt_succ`(`Found/PGC/SenJumpFiltration.lean`)
* **(i) 後半**(区間判定)`lowerRamificationGroup_inf_zpowers_eq_iff`(同上)
* **(iii)**(合同)`ramIndex_pow_pow_congr`(`Found/PGC/SenValuationCongruence.lean`)

## ★★原文の `< n ≤` と本木の `≤ n <` の食い違いは erratum ではない

★★**原典に証明が無い以上、この読み方の根拠は「Lean で両向きが閉じたこと」だけである。**
下の等式 `i_j = 1 + Σ_{i≤j} n_i p^i` は
`exists_padicDigits_card_lowerRamificationGroup` の証明の中(`e1` / `e2`)で実際に使われ、
そこから原文の `< n ≤` が `≤ n <` に化けて両向きが閉じている。

`i_j := i(σ^{p^j})` と置くと (i) 後半は `i_{j−1} ≤ n < i_j` である。一方、原文の展開は
`Σ_{i<j} n_i p^i < n ≤ Σ_{i≤j} n_i p^i` と**両端が逆**に見える。これは

  ★**`n_0 := i_0 − 1`(定数項だけ 1 を差し引く)**

と取れば一致する。実際、本ファイルの `exists_padicDigits_of_lt_of_dvd_sub` が与える等式は

  `i_j = 1 + Σ_{i≤j} n_i p^i`   (`j < m`)

なので `Σ_{i≤j} n_i p^i = i_j − 1` であり、

  `Σ_{i<j} n_i p^i < n ≤ Σ_{i≤j} n_i p^i`
    ⟺ `i_{j−1} − 1 < n ≤ i_j − 1` ⟺ `i_{j−1} ≤ n < i_j`.

★そして `n_0 = i_0 − 1 ≥ 1` は `σ ∈ G_1`(⟺ `1 < i(σ)`、すなわち `i_0 ≥ 2`)から**ちょうど**出る。
★★**原文が `n_i ∈ Z_{≥1}` と書いていること自体が、この読み方を強制している**
——`n_0 := i_0` と読むと `Σ_{i≤j} n_i p^i = i_j` になり、`i_0 ≥ 1` しか出ないうえ
両端の向きが合わない。したがって `−1` は逸脱ではなく**原典の意図の復元**である。
★本ファイルの証明はこの読み方で実際に両向きが閉じることを Lean で確かめている
(`exists_padicDigits_card_lowerRamificationGroup` の最後の 6 行)。

## 何を証明したか

**§1 抽象核 —— 純粋な `ℕ`/`ℤ`(分岐・付値・Galois の語彙が 1 つも出てこない)**

* `exists_index_boundary` : `P 0` かつ `¬ P M` なら `k < M` で `P k ∧ ¬ P (k+1)` が取れる。
  ★「`n` がどの区間に落ちるか」の存在(逆向きで使う)。単調性も順序も要らない。
* ★`exists_padicDigits_of_lt_of_dvd_sub` : `f : ℕ → ℕ` が
  (a) `j+1 < M` で狭義増加 (b) `p^{j+1} ∣ f(j+1) − f j` (c) `2 ≤ f 0` を満たすなら
  `n : ℕ → ℕ` があって `∀ i, 1 ≤ n i` かつ `j < M` で `f j = 1 + Σ_{i≤j} n_i p^i`。
  ★**本体(望遠鏡和)**。`n_0 := f 0 − 1`、`n_{k+1} := (f(k+1) − f k)/p^{k+1}`。
  `n_{k+1} ≥ 1` は「正で `p^{k+1}` で割れる ⟹ `p^{k+1}` 以上」から出る。

**§2 抽象核 —— 巡回 `p` 群(純群論)**

* `orderOf_pow_pow_eq` : `orderOf σ = p^m`、`j ≤ m` ⟹ `orderOf (σ^{p^j}) = p^{m−j}`。
* `exists_generator_of_isCyclic` : `IsCyclic G` と `|G| = p^m` から
  `⟨σ⟩ = ⊤ ∧ orderOf σ = p^m` なる `σ` が取れる。★原文の仮定「`G ≅ ℤ/p^mℤ`」を
  本ファイルの仮定の対 `(htop, hord)` に翻訳する 1 本(逸脱 2)。

**§3 具体層 —— `|G_n|` と跳び位置**

* `card_lowerRamificationGroup_of_mem_Ico` : `i_j ≤ n < i_{j+1}` ⟹ `|G_n| = p^{m−(j+1)}`。
  ★(i) 後半 + `⟨σ⟩ = ⊤`(∴ `G_n ⊓ ⟨σ⟩ = G_n`) + `|⟨σ^{p^{j+1}}⟩| = p^{m−(j+1)}`。
* ★`card_lowerRamificationGroup_eq_pow_iff` : **逆も成り立つ**。
  `j + 1 ≤ m` のとき `|G_n| = p^{m−(j+1)} ⟺ i_j ≤ n < i_{j+1}`。
  ★逆向きは「位数から部分群を復元する」のではなく、**区間の存在**(§1)で押す:
  `|G_n| < |G|` から `σ ∉ G_n`(∴ `i_0 ≤ n`)、`i_m = ⊤`(∴ `¬ i_m ≤ n`)なので
  `i_k ≤ n < i_{k+1}` なる `k < m` が取れ、順方向から `|G_n| = p^{m−(k+1)}`、
  `p` 冪の単射性で `k = j`。★**巡回群の部分群の分類は使っていない**。
* ★★★`exists_padicDigits_card_lowerRamificationGroup` : **Corollary 6.7 そのもの**。

**§4 `PAdicLocalField` への具体化**(`_adjoin` 版 2 本)

## ★逸脱の記録

1. ★**`hσ : σ ∈ lowerRamificationGroup B G 1` を仮定に置いた**。原典は `G ≅ ℤ/p^mℤ` から
   `G = G_1` を導けるが、それには **Proposition 6.2**(`G_0/G_1` の位数が `p` と素)が要る。
   ★**Prop 6.2 はこの木にまだ無い**ので、導かずに仮定として置いた。
   ★**「`G ≅ ℤ/p^mℤ` から `σ ∈ G_1` を導いた」のではない。** 下流(Thm 6.11)が供給する。
   なお入力 `lowerRamificationGroup_inf_zpowers_eq_iff` 自身が `hσ` を要求している。
2. **「`G ≅ ℤ/p^mℤ`」を `(htop : Subgroup.zpowers σ = ⊤)` と `(hord : orderOf σ = p^m)` の
   対で表した**。同型そのものを持ち回らず、生成元 `σ` を明示に取る形である。
   ★両者が同値であることは §2 `exists_generator_of_isCyclic` に置いた
   (`IsCyclic G` + `Nat.card G = p^m` ⟹ そのような `σ` が存在)。
   ★★この形にする実質的な理由は逸脱 1 である —— `σ ∈ G_1` は**生成元ごとの仮定**なので、
   同型だけを仮定すると「どの生成元について `G_1` に入るのか」が言えない。
3. **`1 ≤ j ≤ m − 1` を `1 ≤ j ∧ j + 1 ≤ m` と書いた**。`ℕ` の切り詰め引き算 `m − 1` を
   仮定に出さないため(`m = 0` で空虚に真になる事故の回避)。数学的内容は同じ。
   ★**端点は 1 つも落としていない**: `m ≥ 1` のとき `j ≤ m − 1 ⟺ j + 1 ≤ m` で、
   `j = 1` も `j = m − 1` も入っている。`m = 0` では原文の範囲自体が空なので、
   そのとき本定理は `nn := fun _ => 1` を返す(結論の `∀` が空虚に真)。
4. **`n : ℕ → ℕ` を全域で `1 ≤ n i` とした**。原文は `n_0, …, n_{m−1}` の `m` 個だけを
   要求する。`i ≥ m` はジャンク値 `1` を置いてある。★仮定ではなく結論を強めた方向。
5. **中間補題 `card_lowerRamificationGroup_eq_pow_iff` の `j` の範囲が原文より広い**。
   原文は `1 ≤ j ≤ m − 1`(指数 `m − j ≥ 1`)だが、この補題は `j + 1 ≤ m`、すなわち
   指数 `0`(`G_n` が自明)まで含む。★`i_m = ⊤` なので `n ≥ i_{m−1} ⟺ G_n = 1` が
   同じ式で書ける。原文の範囲に絞ったのが Corollary 6.7 本体の方である。
6. **添字を 1 つずらしてある**(`SenJumpFiltration.lean` の逸脱 2 を引き継ぐ)。
   §3 の中間補題は `i_j ≤ n < i_{j+1} ⟺ |G_n| = p^{m−(j+1)}` と書く。
   Corollary 6.7 本体(`exists_padicDigits_card_lowerRamificationGroup`)は原文どおり
   `1 ≤ j` で書いてあり、証明の中で `j = j' + 1` と割ってこの補題に渡している。

## ★退化の自己検査

* **(D1) `j + 1 ≤ m` を落とすと `n_j` が定義できない**。`j ≥ m` では `i_j = ⊤` で、
  `f j := (i_j).toNat = 0` に潰れるので狭義増加(a)が破れる。★`ℕ∞` のまま
  `i_{j+1} − i_j` と書くと**切り詰め引き算で空虚に真**になる(`lean-idioms.md` #102)。
  本ファイルは (iii) の入力 `ramIndex_pow_pow_congr` がすでに有限代表 `c d : ℕ` を
  取り出した形なので、この罠には触れていない —— `f := fun k => (i_k).toNat` を作る際に
  `hfin : k < m → i_k = ((i_k).toNat : ℕ∞)` を毎回経由している。
* **(D2) `σ ∈ G_1` を落とすと `n_0 ≥ 1` が出ない**。`n_0 = i_0 − 1` なので、
  `i_0 = 1`(= `σ ∈ G_0 \ G_1`、tame な跳び)のとき `n_0 = 0` になり原文の
  `n_i ∈ Z_{≥1}` が破れる。★`hσ` は `hf0 : 2 ≤ f 0` でちょうど 1 度使われる。
* **(D3) `hadj`(`π` が `B` を生成する)を落とすと `i(σ)` が分岐を測らない**。
  `mem_lowerRamificationGroup_iff_lt_ramIndex` が成り立たなくなる。
* **(D4) `htop`(`G` が `σ` で生成される)を落とすと両向きとも壊れる**。
  (i) 後半が測るのは `G_n ⊓ ⟨σ⟩` であって `G_n` ではない。`⟨σ⟩ ≠ ⊤` なら
  `|G_n|` は `⟨σ⟩` の外まで数えるので `p^{m−j}` にならない。
-/

namespace ABC3.Found.PGC

open IsLocalRing ABC3.Skeleton.PGC IsDiscreteValuationRing
open scoped NNReal Valued

def exists_padicDigits_card_lowerRamificationGroup.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 15, item := "Corollary 6.7", sectionId := "cor-6-7" }

/-! ## §1 抽象核 —— 純粋な `ℕ`/`ℤ`

★分岐・付値・Galois の語彙は 1 つも現れない。 -/

/-- **境界の存在** —— `P 0` かつ `¬ P M` なら、`P k ∧ ¬ P (k+1)` なる `k < M` が取れる。

★逆向き(`|G_n| = p^{m−j} ⟹ 区間`)で「`n` がどの区間に落ちるか」を出すのに使う。
単調性も順序も要らない、純粋な有限探索である。 -/
theorem exists_index_boundary {P : ℕ → Prop} (h0 : P 0) :
    ∀ M : ℕ, ¬ P M → ∃ k, k < M ∧ P k ∧ ¬ P (k + 1) := by
  intro M
  induction M with
  | zero => intro h; exact absurd h0 h
  | succ M ih =>
    intro hM
    by_cases hPM : P M
    · exact ⟨M, Nat.lt_succ_self M, hPM, hM⟩
    · obtain ⟨k, hk, hk1, hk2⟩ := ih hPM
      exact ⟨k, by omega, hk1, hk2⟩

/-- ★★**Corollary 6.7 の抽象核(望遠鏡和)** —— `f : ℕ → ℕ` が

* (a) `j + 1 < M` の範囲で**狭義増加**
* (b) `j + 1 < M` の範囲で `p^{j+1} ∣ f(j+1) − f j`(`ℤ` の中の引き算)
* (c) `2 ≤ f 0`

を満たすなら、`n : ℕ → ℕ` があって `∀ i, 1 ≤ n i` かつ

`f j = 1 + Σ_{i ≤ j} n_i p^i`  (`j < M`).

★`n_0 := f 0 − 1`、`n_{k+1} := (f(k+1) − f k)/p^{k+1}`。
`n_{k+1} ≥ 1` は「正の数が `p^{k+1}` で割れる」から、`n_0 ≥ 1` は (c) から出る。
★`M ≤ i` の側はジャンク値 `1` を置く(逸脱 4)。 -/
theorem exists_padicDigits_of_lt_of_dvd_sub {p M : ℕ} {f : ℕ → ℕ}
    (hmono : ∀ j, j + 1 < M → f j < f (j + 1))
    (hdvd : ∀ j, j + 1 < M → (p : ℤ) ^ (j + 1) ∣ (f (j + 1) : ℤ) - (f j : ℤ))
    (hf0 : 2 ≤ f 0) :
    ∃ n : ℕ → ℕ, (∀ i, 1 ≤ n i) ∧
      ∀ j, j < M → f j = 1 + ∑ i ∈ Finset.range (j + 1), n i * p ^ i := by
  have key : ∀ i : ℕ, ∃ c : ℕ, 1 ≤ c ∧
      (i < M → f i = (if i = 0 then 1 else f (i - 1)) + c * p ^ i) := by
    intro i
    match i with
    | 0 => exact ⟨f 0 - 1, by omega, fun _ => by simp; omega⟩
    | (k + 1) =>
      by_cases hk : k + 1 < M
      · have hlt := hmono k hk
        have hd := hdvd k hk
        have hd' : p ^ (k + 1) ∣ f (k + 1) - f k := by
          have hcast : ((f (k + 1) - f k : ℕ) : ℤ) = (f (k + 1) : ℤ) - (f k : ℤ) := by
            push_cast [Nat.cast_sub (le_of_lt hlt)]; ring
          rw [← Int.natCast_dvd_natCast, hcast]
          push_cast
          exact hd
        obtain ⟨c, hc⟩ := hd'
        refine ⟨c, ?_, fun _ => ?_⟩
        · rcases Nat.eq_zero_or_pos c with rfl | h
          · simp at hc; omega
          · exact h
        · rw [if_neg (by omega : ¬ (k + 1 = 0))]
          have hc' : f (k + 1) - f k = c * p ^ (k + 1) := by rw [hc, mul_comm]
          simp only [Nat.add_sub_cancel]
          omega
      · exact ⟨1, le_refl 1, fun h => absurd h hk⟩
  choose n hn1 hn2 using key
  refine ⟨n, hn1, ?_⟩
  intro j
  induction j with
  | zero => intro h; have := hn2 0 h; simpa using this
  | succ j ih =>
    intro h
    have h1 := hn2 (j + 1) h
    have h2 := ih (by omega)
    rw [if_neg (by omega : ¬ (j + 1 = 0))] at h1
    simp only [Nat.add_sub_cancel] at h1
    rw [Finset.sum_range_succ, ← Nat.add_assoc, ← h2, h1]

/-! ## §2 抽象核 —— 巡回 `p` 群(純群論)

★ここも分岐・付値・Galois の語彙は現れない。 -/

/-- `orderOf σ = p^m`、`j ≤ m` なら `orderOf (σ^{p^j}) = p^{m−j}`。 -/
theorem orderOf_pow_pow_eq {G : Type*} [Group G] {σ : G} {p m : ℕ} (hp : p.Prime)
    (hord : orderOf σ = p ^ m) {j : ℕ} (hj : j ≤ m) :
    orderOf (σ ^ p ^ j) = p ^ (m - j) := by
  rw [orderOf_pow' _ (pow_ne_zero j hp.pos.ne'), hord,
    Nat.gcd_eq_right (pow_dvd_pow p hj), Nat.pow_div hj hp.pos]

/-- ★**原文の仮定「`G ≅ ℤ/p^mℤ`」の翻訳**(逸脱 2) —— `IsCyclic G` と `|G| = p^m` から

`⟨σ⟩ = ⊤` かつ `orderOf σ = p^m` なる生成元 `σ` が取れる。

★本ファイルの主定理はこの `(htop, hord)` の対を仮定に取る。 -/
theorem exists_generator_of_isCyclic {G : Type*} [Group G] [IsCyclic G] {p m : ℕ}
    (hcard : Nat.card G = p ^ m) :
    ∃ σ : G, Subgroup.zpowers σ = ⊤ ∧ orderOf σ = p ^ m := by
  obtain ⟨g, hg⟩ := IsCyclic.exists_generator (α := G)
  have htop : Subgroup.zpowers g = ⊤ := (Subgroup.eq_top_iff' _).2 hg
  refine ⟨g, htop, ?_⟩
  rw [← Nat.card_zpowers, htop, Subgroup.card_top, hcard]

/-! ## §3 具体層 —— `|G_n|` と跳び位置 `i_j`

`i_j := ramIndex π (σ^{p^j})`。★添字は `SenJumpFiltration.lean` の流儀(1 つずらし)に従う。 -/

/-- `i_j ≤ n < i_{j+1}` なら `|G_n| = p^{m−(j+1)}`(Corollary 6.7 の `⟸` 側)。

★(i) 後半 `lowerRamificationGroup_inf_zpowers_eq_iff` で `G_n ⊓ ⟨σ⟩ = ⟨σ^{p^{j+1}}⟩` を出し、
`⟨σ⟩ = ⊤` で左辺を `G_n` に潰し、`|⟨σ^{p^{j+1}}⟩| = orderOf (σ^{p^{j+1}}) = p^{m−(j+1)}`。 -/
theorem card_lowerRamificationGroup_of_mem_Ico {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (huni : maximalIdeal B = Ideal.span {π})
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    (htop : Subgroup.zpowers σ = ⊤) {j : ℕ} (hj : j + 1 ≤ m) {n : ℕ}
    (h1 : ramIndex π (σ ^ p ^ j) ≤ (n : ℕ∞)) (h2 : (n : ℕ∞) < ramIndex π (σ ^ p ^ (j + 1))) :
    Nat.card (lowerRamificationGroup B G n) = p ^ (m - (j + 1)) := by
  have heq : lowerRamificationGroup B G n = Subgroup.zpowers (σ ^ p ^ (j + 1)) := by
    have h := (lowerRamificationGroup_inf_zpowers_eq_iff (A := A) p hp huni hadj hσ hord hj n).2
      ⟨h1, h2⟩
    rwa [htop, inf_top_eq] at h
  rw [heq, Nat.card_zpowers, orderOf_pow_pow_eq hp hord hj]

/-- ★★**`|G_n| = p^{m−(j+1)} ⟺ i_j ≤ n < i_{j+1}`**(`j + 1 ≤ m`)。

★`⟸` は上の補題。`⟹` が本題で、**巡回群の部分群の分類は使わない**:
`|G_n| = p^{m−(j+1)} < p^m = |G|` から `G_n ≠ ⊤`、`⟨σ⟩ = ⊤` なので `σ ∉ G_n`、
すなわち `i_0 ≤ n`。一方 `i_m = ⊤` だから `¬ (i_m ≤ n)`。よって §1 の
`exists_index_boundary` で `i_k ≤ n < i_{k+1}` なる `k < m` が取れる。`⟸` を `k` に当てると
`|G_n| = p^{m−(k+1)}` で、`p` 冪の単射性(`Nat.pow_right_injective`)から `k = j`。

★★`j + 1 ≤ m` は落とせない(退化検査 D1)。★`htop` も落とせない(D4)。 -/
theorem card_lowerRamificationGroup_eq_pow_iff {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (huni : maximalIdeal B = Ideal.span {π})
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    (htop : Subgroup.zpowers σ = ⊤) {j : ℕ} (hj : j + 1 ≤ m) (n : ℕ) :
    Nat.card (lowerRamificationGroup B G n) = p ^ (m - (j + 1)) ↔
      ramIndex π (σ ^ p ^ j) ≤ (n : ℕ∞) ∧ (n : ℕ∞) < ramIndex π (σ ^ p ^ (j + 1)) := by
  refine ⟨fun hcard => ?_, fun h =>
    card_lowerRamificationGroup_of_mem_Ico (A := A) p hp huni hadj hσ hord htop hj h.1 h.2⟩
  have hnotmem : σ ∉ lowerRamificationGroup B G n := by
    intro hmem
    have htopn : lowerRamificationGroup B G n = ⊤ :=
      eq_top_iff.2 (le_of_eq_of_le htop.symm (Subgroup.zpowers_le.2 hmem))
    rw [htopn, ← htop, Nat.card_zpowers, hord] at hcard
    have := Nat.pow_right_injective hp.two_le hcard
    omega
  have h0 : ramIndex π (σ ^ p ^ 0) ≤ (n : ℕ∞) := by
    rw [pow_zero, pow_one]
    exact not_lt.1 fun hc =>
      hnotmem ((mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj n σ).2 hc)
  have hm : ¬ (ramIndex π (σ ^ p ^ m) ≤ (n : ℕ∞)) := by
    rw [(ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord m).2 (le_refl m)]
    exact fun hc => (ENat.coe_ne_top n) (top_le_iff.1 hc)
  obtain ⟨k, hk, hk1, hk2⟩ :=
    exists_index_boundary (P := fun k => ramIndex π (σ ^ p ^ k) ≤ (n : ℕ∞)) h0 m hm
  have hkcard := card_lowerRamificationGroup_of_mem_Ico (A := A) p hp huni hadj hσ hord htop
    (show k + 1 ≤ m by omega) hk1 (not_le.1 hk2)
  have hkj : k = j := by
    rw [hcard] at hkcard
    have := Nat.pow_right_injective hp.two_le hkcard
    omega
  subst hkj
  exact ⟨hk1, not_le.1 hk2⟩

/-- ★★★**Yoshida 2008 Corollary 6.7** —— `G = ⟨σ⟩`、`orderOf σ = p^m`、`σ ∈ G_1` のとき

`n_0, n_1, … ∈ ℤ_{≥1}` があって、`1 ≤ j`・`j + 1 ≤ m` なる `j` について

`|G_n| = p^{m−j} ⟺ Σ_{i<j} n_i p^i < n ≤ Σ_{i≤j} n_i p^i`.

★`n_0 = i_0 − 1` である(ファイル冒頭「原文の `< n ≤` と本木の `≤ n <`」を見よ)。
★`hσ` は仮定である(逸脱 1、Prop 6.2 未形式化)。★`(htop, hord)` が
「`G ≅ ℤ/p^mℤ`」の翻訳(逸脱 2)。★`1 ≤ j ≤ m−1` は `1 ≤ j ∧ j+1 ≤ m` と書いた(逸脱 3)。 -/
theorem exists_padicDigits_card_lowerRamificationGroup {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (hπ : Irreducible π)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    (htop : Subgroup.zpowers σ = ⊤) :
    ∃ nn : ℕ → ℕ, (∀ i, 1 ≤ nn i) ∧
      ∀ j : ℕ, 1 ≤ j → j + 1 ≤ m → ∀ n : ℕ,
        Nat.card (lowerRamificationGroup B G n) = p ^ (m - j) ↔
          (∑ i ∈ Finset.range j, nn i * p ^ i < n ∧
            n ≤ ∑ i ∈ Finset.range (j + 1), nn i * p ^ i) := by
  have huni := (irreducible_iff_uniformizer π).mp hπ
  rcases Nat.eq_zero_or_pos m with rfl | hm1
  · exact ⟨fun _ => 1, fun _ => le_refl 1, by intro j hj1 hj2; omega⟩
  -- `f k := (i_k).toNat`。`k < m` の範囲で `i_k` は有限(#102 の回避)。
  have hfin : ∀ k, k < m → ramIndex π (σ ^ p ^ k) = ((ramIndex π (σ ^ p ^ k)).toNat : ℕ∞) := by
    intro k hk
    refine (ENat.coe_toNat ?_).symm
    intro hc
    have := (ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord k).1 hc
    omega
  -- (a) 狭義増加(Prop 6.6 (i) 前半)
  have hmono : ∀ j, j + 1 < m →
      (ramIndex π (σ ^ p ^ j)).toNat < (ramIndex π (σ ^ p ^ (j + 1))).toNat := by
    intro j hj
    have hlt := ramIndex_pow_pow_lt_succ (A := A) p hp huni hadj hσ hord (show j < m by omega)
    rw [hfin j (by omega), hfin (j + 1) (by omega)] at hlt
    exact_mod_cast hlt
  -- (b) 合同(Prop 6.6 (iii))
  have hdvd : ∀ j, j + 1 < m → (p : ℤ) ^ (j + 1) ∣
      ((ramIndex π (σ ^ p ^ (j + 1))).toNat : ℤ) - ((ramIndex π (σ ^ p ^ j)).toNat : ℤ) := by
    intro j hj
    exact ramIndex_pow_pow_congr (A := A) p hp hπ hadj hres hσ hord (hfin j (by omega))
      (hfin (j + 1) (by omega))
  -- (c) `2 ≤ i_0` —— ★`σ ∈ G_1` はここで 1 度だけ使う(退化検査 D2)
  have hf0 : 2 ≤ (ramIndex π (σ ^ p ^ 0)).toNat := by
    have h1 : (1 : ℕ∞) < ramIndex π (σ ^ p ^ 0) := by
      rw [pow_zero, pow_one]
      exact (mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj 1 σ).1 hσ
    rw [hfin 0 hm1] at h1
    have h2 : (1 : ℕ) < (ramIndex π (σ ^ p ^ 0)).toNat := by exact_mod_cast h1
    omega
  obtain ⟨nn, hnn1, hnn2⟩ := exists_padicDigits_of_lt_of_dvd_sub (p := p) (M := m)
    (f := fun k => (ramIndex π (σ ^ p ^ k)).toNat) hmono hdvd hf0
  refine ⟨nn, hnn1, ?_⟩
  intro j hj1 hjm n
  obtain ⟨j', rfl⟩ : ∃ j', j = j' + 1 := ⟨j - 1, by omega⟩
  rw [card_lowerRamificationGroup_eq_pow_iff (A := A) p hp huni hadj hσ hord htop
    (show j' + 1 ≤ m by omega) n, hfin j' (by omega), hfin (j' + 1) (by omega)]
  -- ★ここが `n_0 = i_0 − 1` の検算: `i_{j'} = 1 + Σ_{i≤j'} n_i p^i`。
  have e1 : (ramIndex π (σ ^ p ^ j')).toNat
      = 1 + ∑ i ∈ Finset.range (j' + 1), nn i * p ^ i := hnn2 j' (by omega)
  have e2 : (ramIndex π (σ ^ p ^ (j' + 1))).toNat
      = 1 + ∑ i ∈ Finset.range (j' + 1 + 1), nn i * p ^ i := hnn2 (j' + 1) (by omega)
  constructor
  · rintro ⟨ha, hb⟩
    have ha' : (ramIndex π (σ ^ p ^ j')).toNat ≤ n := by exact_mod_cast ha
    have hb' : n < (ramIndex π (σ ^ p ^ (j' + 1))).toNat := by exact_mod_cast hb
    omega
  · rintro ⟨ha, hb⟩
    exact ⟨by exact_mod_cast (show (ramIndex π (σ ^ p ^ j')).toNat ≤ n by omega),
      by exact_mod_cast (show n < (ramIndex π (σ ^ p ^ (j' + 1))).toNat by omega)⟩

/-! ## §4 `PAdicLocalField` への具体化

`B := adjoinIntegers K x`、`A := 𝒪[K.carrier]`、`G := Gal(K(x)/K)`。
★Y7a §7・Y7b §4 と同型で、仮定の束はそのまま流用できる。 -/

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-- **`|G_n| = p^{m−(j+1)} ⟺ i_j ≤ n < i_{j+1}`(`PAdicLocalField` 版)**。 -/
theorem card_lowerRamificationGroupAdjoin_eq_pow_iff (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {π : adjoinIntegers K x} (hπ : Irreducible π)
    {σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))}
    (hσ : σ ∈ lowerRamificationGroupAdjoin K x 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    (htop : Subgroup.zpowers σ = ⊤) {j : ℕ} (hj : j + 1 ≤ m) (n : ℕ) :
    Nat.card (lowerRamificationGroupAdjoin K x n) = p ^ (m - (j + 1)) ↔
      ramIndex π (σ ^ p ^ j) ≤ (n : ℕ∞) ∧ (n : ℕ∞) < ramIndex π (σ ^ p ^ (j + 1)) := by
  haveI := charP_residueField_adjoinIntegers K x
  have huni := (irreducible_iff_uniformizer π).mp hπ
  have hadj := adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni
  exact card_lowerRamificationGroup_eq_pow_iff (A := 𝒪[K.carrier]) p Fact.out huni hadj
    hσ hord htop hj n

/-- ★★★**Yoshida 2008 Corollary 6.7(`PAdicLocalField` 版)**。

`G = ⟨σ⟩ ≅ ℤ/p^mℤ`(= `htop` + `hord`)と `σ ∈ G_1`(= `hσ`、逸脱 1)の下で、
跳びの位置が `p` 進展開 `n_0 + n_1 p + n_2 p^2 + ⋯`(各 `n_i ≥ 1`)の形をとる。 -/
theorem exists_padicDigits_card_lowerRamificationGroupAdjoin (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {π : adjoinIntegers K x} (hπ : Irreducible π)
    {σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))}
    (hσ : σ ∈ lowerRamificationGroupAdjoin K x 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    (htop : Subgroup.zpowers σ = ⊤) :
    ∃ nn : ℕ → ℕ, (∀ i, 1 ≤ nn i) ∧
      ∀ j : ℕ, 1 ≤ j → j + 1 ≤ m → ∀ n : ℕ,
        Nat.card (lowerRamificationGroupAdjoin K x n) = p ^ (m - j) ↔
          (∑ i ∈ Finset.range j, nn i * p ^ i < n ∧
            n ≤ ∑ i ∈ Finset.range (j + 1), nn i * p ^ i) := by
  haveI := charP_residueField_adjoinIntegers K x
  have huni := (irreducible_iff_uniformizer π).mp hπ
  have hadj := adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni
  exact exists_padicDigits_card_lowerRamificationGroup (A := 𝒪[K.carrier]) p Fact.out hπ hadj
    (exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin K x ht) hσ hord htop

end ABC3.Found.PGC
