import ABC3.Found.PGC.AbelianJumpDivisibility
import ABC3.Found.PGC.HerbrandComposition
import ABC3.Found.PGC.SenJumpExpansion

/-!
# Hasse-Arf の定理(Yoshida 2008 Theorem 6.11)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Theorem 6.11 (Hasse-Arf)**
(物理 p.16)。構造化済み原文は
`ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/section-6.html`
の `id="thm-6-11"`(`data-pdf-page="16"`, `data-item="Theorem 6.11 (Hasse-Arf)"`)。

原文 (Yoshida08 p.16):
> Theorem 6.11 (Hasse-Arf). If G is abelian, n ∈ Z[bb]_≥0 and G_n = G_n+1, then φ_G(n) ∈ Z[bb]_≥0.

★★上の逐語は `pdftotext` に見えている形であり、原典ではない。原典は `G_n ≠ G_{n+1}` である。
`≠` の斜線はベクター描画なので `pdftotext` が落とし、出力は `=` だけになる。
テキストだけを読むと主張が反転し、空虚に真な statement ができる。
構造化済み HTML は `data-txt="="` でこの脱落を明示している。
本ファイルの形式化はすべて `≠`(`hne : G_n ≠ G_{n+1}`)側で書いてある
(先行ノード `AbelianJumpDivisibility.lean` / `RamificationJumpDivisibility.lean` と同じ扱い)。

## 原典の証明の 3 つの段と、本ファイルの被覆

原文 (Yoshida08 p.16-17) の Proof は次の 3 段からなる。

1. **`G = G_1` かつ `j = 1`(`G ≅ ℤ/p^mℤ`、巡回)**:
   `G_n ≠ G_{n+1}` なら `n = Σ_{i=0}^{j} n_i p^i`(Corollary 6.7)で、
   `φ_G(n) = (1/p^m)(n_0 p^m + n_1 p·p^{m−1} + ⋯ + n_j p^j·p^{m−j}) ∈ ℤ≥0`(Lemma 6.10(i))。
   → ★**本ファイル §2 で完全に埋めた**(`exists_natCast_herbrandPhiGroup_of_isCyclic`)。
2. **`G = G_1` かつ `j > 1`**: 巡回商 `G/H ≅ ℤ/p^{m_i}ℤ` を取り、Proposition 6.9 と
   Lemma 6.10(ii) で `j` に関する帰納。
   → ★**本ファイルには無い**(理由は下の「残した段」)。
3. **`G ≠ G_1`**: `H = G_1`、`e_0 = |G/H|` として `e_0 ∣ φ_H(n)` を示す(Corollary 6.3)。
   → ★**本ファイル §3 で埋めた**。実質部分(`e_0 ∣ φ_H(n)`)は
   `natCard_quot_dvd_of_herbrandPhi_natCast` として**仮定なしで**証明した。
   組み立て(`φ_G(n) = φ_H(n)/e_0`)は Lemma 6.10(ii) を仮定として受け取る
   (`exists_natCast_herbrandPhiGroup_of_tame_quotient`。下の「逸脱の記録 (4)」)。

## ★★残した段(段 2)と、その理由

段 2 は**配管が越えられない**。原文の帰納は各段で

* 商群 `G/H` が「別の環 `C`(= `H` の固定環)に忠実に作用する」ことを使い、
* その `C` 上で Proposition 6.9(`herbrand_mem_iff`)と Lemma 6.10(ii)
  (`herbrandPhiGroup_comp`)を回す。

木の Proposition 6.9 / Lemma 6.10(ii) は、この固定環 `C` と 7 本の適合条件
(`hcomp` `hHtriv` `hinj` `hfixC` `hres` `hAC` `hfix`)を**仮定として**受け取る形で
立っている。それらを `G/H` について実際に構成するのが `Found/PGC/FixedRingTower.lean`
(同時並行で形式化中)であり、本ファイルの時点では存在しない。
★したがって段 2 は**数学が足りない**のではなく**配管が届いていない**。
新ノードとして報告してある。

★**有限アーベル群の構造定理は mathlib にある**:
`AddCommGroup.equiv_directSum_zmod_of_finite`(`Mathlib/GroupTheory/FiniteAbelian/Basic.lean`)。
双対(点を分離する巡回商)は `CommGroup.exists_apply_ne_one_of_hasEnoughRootsOfUnity`
(`Mathlib/GroupTheory/FiniteAbelian/Duality.lean`)。★どちらも本木の import には**入っていない**
(`Unknown constant` で返る)ので、段 2 に着手するノードは import を 1 行足す必要がある。

## 抽象核と具体層(設計)

★分岐・付値・Galois の語彙が 1 つも出てこないものを §1 に集めた。§2 §3 はそこへの代入である。

| 抽象核(§1) | 内容 |
|---|---|
| `le_of_forall_le_succ` | 隣接不等式の連鎖 |
| `dvd_sum_Ico_of_blockwise` | ★★段 1 の心臓。区間 `[bd 0, bd L)` の上でブロックごとに定数な `c` の和は、各ブロックの「長さ × 値」の和である。各項が `d` で割れれば和も割れる |
| `exists_boundary_of_ne` | ★★跳びは境界でしか起きない。ブロックごとに定数な `F : ℕ → X` が `F n ≠ F (n+1)` を満たすなら `n+1` は境界の 1 つ |
| `truncENat_one_add` | `min{1, x+1} = 1`(`x ≥ 0`) |
| `phiOf_of_two_valued` | ★★「順分岐なら φ は線型」。`f` が `T` 上で `⊤`、`T` の外で `1` なら `φ_f(x) = |T| x / |S|` |

★`dvd_sum_Ico_of_blockwise` と `exists_boundary_of_ne` は原文が**一言も書いていない**部分である。
原文は「`n = Σ n_i p^i` なら `φ_G(n) = (1/p^m)(n_0 p^m + ⋯)`」と**答だけ**を書いており、
「なぜ `[1,n]` 上の和がブロックに割れるか」「なぜ跳びが `p` 進展開の切れ目でしか起きないか」は
畳まれている。★この 2 本がその中身である。

## ★★段 1 の証明は原文と経路が違う(逸脱の記録 (1))

原文は Corollary 6.7 の `p` 進展開 `n_i` を経由して
`φ_G(n) = n_0 + n_1 + ⋯ + n_j` という**値**を出す。本ファイルは
**値を出さず整除性だけを出す**:

`i_j := i(σ^{p^j})` と置くと、ブロック `[i_{j−1}, i_j)` の上で `|G_k| = p^{m−j}` であり、
ブロックの長さは `i_j − i_{j−1}`、これは Proposition 6.6 (iii)
(`ramIndex_pow_pow_congr`)により `p^j` で割れる。よって
「長さ × 値」は `p^j · p^{m−j} = p^m` で割れる。全ブロックを足して `p^m ∣ Σ_{k=1}^n |G_k|`、
Lemma 6.10(i) で `φ_G(n) = (Σ_{k=1}^n |G_k|)/p^m ∈ ℤ≥0`。

★原文の `n_i` は「`i_j = 1 + Σ_{i≤j} n_i p^i`」の差分にほかならない
(`SenJumpExpansion.lean` 冒頭の検算)ので、**同じ計算を差分の形で書いただけ**である。
Corollary 6.7 の展開そのものを経由しないので、原文が要求する `n_i ≥ 1` も
`1 ≤ j ≤ m−1` の添字条件も要らなくなり、`m = 0` `m = 1` の退化も自動で通る。
★結論(`φ_G(n) ∈ ℤ≥0`)は原文と同じである。

## 逸脱の記録

1. **段 1 は `p` 進展開(Corollary 6.7)を経由しない**(上記)。
   使うのは Proposition 6.6 の (i)(iii)、すなわち Corollary 6.7 の**入力側**である。
2. **`G` の可換性は `∀ x y : G, x * y = y * x` という命題で渡す**(`CommGroup` 構造を
   要求しない)。先行ノード `AbelianJumpDivisibility.lean` の逸脱 1 と同じ理由:
   具体層の `G = Gal(K(x)/K)` は `Group` インスタンスしか持たない。
   ★段 1(巡回)では可換性そのものは使わない(`Subgroup.zpowers σ = ⊤` が強い)。
   ★段 3 は Corollary 6.3 経由で可換性を使う。
3. **`G ≅ ℤ/p^mℤ`(段 1)は `Subgroup.zpowers σ = ⊤ ∧ orderOf σ = p^m` で表す**
   (`SenJumpExpansion.lean` の逸脱 2 と同じ)。さらに `σ ∈ G_1` を仮定する
   (完全分岐 `p` 群であること。`SenJumpExpansion.lean` の逸脱 1 と同じ)。
4. **段 3 の「`φ_G(n) = φ_{G/H}(φ_H(n))`」(Lemma 6.10(ii))は仮定として受け取る**。
   原文が `by Lemma 6.10(ii)` で済ませている一歩であり、本木ではその Lemma が
   固定環 `C` と 7 本の適合条件を要求する。★実質部分(`e_0 ∣ φ_H(n)`)は
   `natCard_quot_dvd_of_herbrandPhi_natCast` として**仮定なしで**証明してあるので、
   固定環が立った時点で `hcomp` を差し込めば段 3 が閉じる。
   ★先行ノード `AbelianJumpDivisibility.lean` の逸脱 4(Lemma 6.10(i) を仮定で受ける)と
   同じ流儀である。
5. **`φ_H(n) ∈ ℤ≥0` は仮定として受け取る**(段 3)。原文も
   「(we know φ_H(n) ∈ Z≥0)」と括弧で断っており、これは段 1-2(`H = G_1` は `p` 群)の
   帰結である。段 2 が無い以上ここは仮定にせざるを得ない。

## 退化の自己検査

* ★★**`G_n ≠ G_{n+1}` を落とすと空虚**。落とすと `n` が任意になり、
  `φ_G(n)` は一般に整数でない。
  `exists_boundary_of_ne` の `hne` がこの仮定を 1 度だけ、しかし決定的に使う場所である。
* ★★**段 3 で `e_0` と `|H|` が互いに素であることを落とすと出ない**。これは
  `AbelianJumpDivisibility.lean` の `coprime_natCard_quot_natCard_lowerRamificationGroup_one`
  が担っており、本ファイルは `dvd_of_sum_eq_natCard_mul` を通してそれを使う。
* ★**`G` 可換を落とすと段 3 が偽**(Hasse-Arf はアーベルでのみ成り立つ)。
  段 3 の `habel` は `dvd_of_sum_eq_natCard_mul` にそのまま渡っている。
* ★**`φ_G(n) ∈ ℤ≥0` の `≥0`**: 結論は `∃ k : ℕ, φ_G(n) = k` の形にしてあるので
  非負性は型に入っている(`herbrandPhiGroup_nonneg_of_isCyclic` で明示的にも出せる)。
* ★**`ℕ∞` の切り詰め引き算・除算は書いていない**(`lean-idioms.md` #102)。
  `i_j` の差は `toNat` を取ってから `ℕ` の中で引いており、有限性(`hfin`)を先に確保している。
  除算は `phiOf` の定義の中にしかない。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

def exists_natCast_herbrandPhiGroup_of_isCyclic.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

def exists_natCast_herbrandPhiGroup_of_tame_quotient.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

/-! ## §1 抽象核

★分岐・付値・Galois の語彙が 1 つも出てこない。純粋な算術・順序論・有限和である。 -/

/-- 隣接不等式 `bd j ≤ bd (j+1)`(`j < L`)から端点の不等式 `bd 0 ≤ bd L` を出す。 -/
theorem le_of_forall_le_succ {bd : ℕ → ℕ} {L : ℕ}
    (h : ∀ j, j < L → bd j ≤ bd (j + 1)) : bd 0 ≤ bd L := by
  induction L with
  | zero => exact le_rfl
  | succ L ih => exact le_trans (ih fun j hj => h j (by omega)) (h L (by omega))

/-- ★★★**段 1 の心臓**(抽象核)——ブロックごとに定数な列の区間和の整除性。

`bd 0 ≤ bd 1 ≤ ⋯ ≤ bd L` を境界とし、`c` が第 `j` ブロック `[bd j, bd (j+1))` の上で
値 `v j` を取るとする。このとき区間和は各ブロックの「長さ × 値」の和であり、
各項が `d` で割れれば全体も `d` で割れる。★ここでは結論だけ(整除性)を述べる。

★★**原文はこの補題を一言も書いていない。** 原文は
`φ_G(n) = (1/p^m)(n_0·p^m + n_1 p·p^{m−1} + ⋯)` と**答だけ**を書いており、
「`[1,n]` 上の和がブロックに割れる」という中身がここである。 -/
theorem dvd_sum_Ico_of_blockwise {d : ℕ} {c v bd : ℕ → ℕ} {L : ℕ}
    (hmono : ∀ j, j < L → bd j ≤ bd (j + 1))
    (hconst : ∀ j, j < L → ∀ k, bd j ≤ k → k < bd (j + 1) → c k = v j)
    (hdvd : ∀ j, j < L → d ∣ (bd (j + 1) - bd j) * v j) :
    d ∣ ∑ k ∈ Finset.Ico (bd 0) (bd L), c k := by
  induction L with
  | zero => simp
  | succ L ih =>
    have h0L : bd 0 ≤ bd L := le_of_forall_le_succ (fun j hj => hmono j (by omega))
    have hLL : bd L ≤ bd (L + 1) := hmono L (by omega)
    rw [← Finset.sum_Ico_consecutive c h0L hLL]
    refine dvd_add (ih (fun j hj => hmono j (by omega)) (fun j hj => hconst j (by omega))
      (fun j hj => hdvd j (by omega))) ?_
    have hb : ∑ k ∈ Finset.Ico (bd L) (bd (L + 1)), c k = (bd (L + 1) - bd L) * v L := by
      rw [Finset.sum_congr rfl (fun k hk => hconst L (by omega) k (Finset.mem_Ico.1 hk).1
        (Finset.mem_Ico.1 hk).2), Finset.sum_const, Nat.card_Ico, smul_eq_mul]
    rw [hb]
    exact hdvd L (by omega)

/-- ★★★**跳びは境界でしか起きない**(抽象核)。

`F : ℕ → X` が「`b 0` 未満では一定」「各ブロック `[b j, b (j+1))`(`j < m`)の上で一定」
であるとき、`F n ≠ F (n+1)` なら `n + 1` は境界のどれか、すなわち
`b l = n + 1` なる `l < m` が存在する。

★★**この主張が `G_n ≠ G_{n+1}` という仮定を使い切る唯一の場所**である。
★仮定を落とすと `n` が任意になり結論が偽になる(退化検査)。
★`b` の値域を `ℕ∞` にしてあるのは、具体層で `b m = ⊤`(`σ^{p^m} = 1`)だからである
(`htop` がその形)。 -/
theorem exists_boundary_of_ne {X : Type*} (F : ℕ → X) (b : ℕ → ℕ∞) (m n : ℕ)
    (hbot : ∀ x y : ℕ, (x : ℕ∞) < b 0 → (y : ℕ∞) < b 0 → F x = F y)
    (hblk : ∀ j, j < m → ∀ x y : ℕ, b j ≤ (x : ℕ∞) → (x : ℕ∞) < b (j + 1) →
      b j ≤ (y : ℕ∞) → (y : ℕ∞) < b (j + 1) → F x = F y)
    (htop : ¬ (b m ≤ ((n + 1 : ℕ) : ℕ∞)))
    (hne : F n ≠ F (n + 1)) :
    ∃ l, l < m ∧ b l = ((n + 1 : ℕ) : ℕ∞) := by
  have hnn : ((n : ℕ∞)) < ((n + 1 : ℕ) : ℕ∞) := by exact_mod_cast Nat.lt_succ_self n
  by_cases h0 : b 0 ≤ ((n + 1 : ℕ) : ℕ∞)
  · obtain ⟨l, hl, hl1, hl2⟩ :=
      exists_index_boundary (P := fun j => b j ≤ ((n + 1 : ℕ) : ℕ∞)) h0 m htop
    refine ⟨l, hl, le_antisymm hl1 ?_⟩
    by_contra hcon
    have hlt : b l < ((n + 1 : ℕ) : ℕ∞) := lt_of_le_of_ne hl1 (fun h => hcon (le_of_eq h.symm))
    have hle : b l ≤ (n : ℕ∞) :=
      Order.le_of_lt_add_one (x := b l) (y := (n : ℕ∞)) (by simpa using hlt)
    have hup : ((n + 1 : ℕ) : ℕ∞) < b (l + 1) := not_le.1 hl2
    exact hne (hblk l hl n (n + 1) hle (lt_trans hnn hup) hl1 hup)
  · exact absurd (hbot n (n + 1) (lt_trans hnn (not_le.1 h0)) (not_le.1 h0)) hne

/-- `min{1, x+1} = 1`(`x ≥ 0`)。 -/
theorem truncENat_one_add {x : ℝ} (hx : 0 ≤ x) : truncENat (1 : ℕ∞) (x + 1) = 1 := by
  rw [show (1 : ℕ∞) = ((1 : ℕ) : ℕ∞) by rfl, truncENat_coe, Nat.cast_one]
  exact min_eq_left (by linarith)

/-- ★★★**「順分岐なら φ は線型」**(抽象核)——`f` が `T` の上で `⊤`、`T` の外で `1` なら

`φ_f(x) = |T| · x / |S|`  (`x ≥ 0`).

★具体層では `S = G`、`T = G_1`、`f = i_ϖ`(固定環 `C` の素元による跳び指数)であり、
結論は原文の「φ_{G/H}(n) = n/e_0 for n ∈ R≥0 by definition」になる。
★★原文が `by definition` で畳んでいるのは**この計算**である(「定義から」ではない:
`i_ϖ` が `H` の外で `1` であること = `(G/H)_1 = 1` は Proposition 6.9 が要る)。 -/
theorem phiOf_of_two_valued {S : Type*} [Fintype S] [Nonempty S] [DecidableEq S] (f : S → ℕ∞)
    (T : Finset S) (h1 : ∀ s ∈ T, f s = ⊤) (h2 : ∀ s ∉ T, f s = 1) {x : ℝ} (hx : 0 ≤ x) :
    phiOf f x = (T.card : ℝ) * x / (Nat.card S : ℝ) := by
  have hcard : (Nat.card S : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Nat.card_pos.ne'
  have hT : ∑ s ∈ T, truncENat (f s) (x + 1) = (T.card : ℝ) * (x + 1) := by
    rw [Finset.sum_congr rfl (fun s hs => by rw [h1 s hs, truncENat_top]), Finset.sum_const,
      nsmul_eq_mul]
  have hTc : ∑ s ∈ Tᶜ, truncENat (f s) (x + 1) = ((Tᶜ.card : ℕ) : ℝ) := by
    rw [Finset.sum_congr rfl (fun s hs => by
        rw [h2 s (Finset.mem_compl.1 hs), truncENat_one_add hx]),
      Finset.sum_const, nsmul_eq_mul, mul_one]
  have hcompl : ((Tᶜ.card : ℕ) : ℝ) = (Nat.card S : ℝ) - (T.card : ℝ) := by
    rw [Finset.card_compl, Nat.cast_sub (Finset.card_le_univ T), Nat.card_eq_fintype_card]
  have hsum : ∑ s : S, truncENat (f s) (x + 1)
      = (T.card : ℝ) * (x + 1) + ((Nat.card S : ℝ) - (T.card : ℝ)) := by
    rw [← Finset.sum_add_sum_compl T (fun s => truncENat (f s) (x + 1)), hT, hTc, hcompl]
  rw [phiOf, hsum]
  field_simp
  ring

/-- `∃ k : ℕ, x = k` から `0 ≤ x`。★原典の結論 `φ_G(n) ∈ ℤ≥0` の `≥0` を明示するための系。 -/
theorem nonneg_of_exists_natCast {x : ℝ} (h : ∃ k : ℕ, x = (k : ℝ)) : 0 ≤ x := by
  obtain ⟨k, hk⟩ := h
  rw [hk]
  exact Nat.cast_nonneg k

/-! ## §2 段 1 —— `G ≅ ℤ/p^mℤ`(巡回、完全分岐 `p` 群)の場合

原文:
> When j = 1, i.e. G ∼= Z/pmZ, if Gn ̸= Gn+1 then n = Σ_{i=0}^{j} n_i p^i for some
> 0 ≤ j ≤ m−1 by Corollary 6.7, in which case φG(n) = (1/p^m)(n0·p^m + n1p·p^{m−1} +
> ⋯ + njp^j·p^{m−j}) ∈ Z≥0 by Lemma 6.10(i).

★`G ≅ ℤ/p^mℤ` は `Subgroup.zpowers σ = ⊤ ∧ orderOf σ = p^m` で表す(逸脱 3)。 -/

/-- `x < i(σ)` かつ `⟨σ⟩ = G` なら `G_x = ⊤`。
★ブロック `[1, i_0)` の上での `|G_k| = p^m` の根拠。 -/
theorem lowerRamificationGroup_eq_top_of_lt_ramIndex {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] {π : B}
    (huni : maximalIdeal B = Ideal.span {π}) (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    {σ : G} (hgen : Subgroup.zpowers σ = ⊤) {x : ℕ} (h : (x : ℕ∞) < ramIndex π σ) :
    lowerRamificationGroup B G x = ⊤ :=
  eq_top_iff.2 (le_of_eq_of_le hgen.symm (Subgroup.zpowers_le.2
    ((mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj x σ).2 h)))

/-- ★★**跳びの位置は `i_l := i(σ^{p^l})` のどれか**。

`G = ⟨σ⟩`、`orderOf σ = p^m`、`σ ∈ G_1` のとき、`G_n ≠ G_{n+1}` なら
`l < m` があって `i_l = n + 1`。

★段取り: §1 の `exists_boundary_of_ne` に代入するだけ。3 つの入力は
(a) `x < i_0` なら `G_x = ⊤`(上の補題)、
(b) `i_j ≤ x < i_{j+1}` なら `G_x = ⟨σ^{p^{j+1}}⟩`(Proposition 6.6 (i) 後半)、
(c) `i_m = ⊤`。 -/
theorem exists_ramIndex_eq_of_lowerRamificationGroup_ne {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (huni : maximalIdeal B = Ideal.span {π})
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    (hgen : Subgroup.zpowers σ = ⊤) {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    ∃ l, l < m ∧ ramIndex π (σ ^ p ^ l) = ((n + 1 : ℕ) : ℕ∞) := by
  refine exists_boundary_of_ne (fun k => lowerRamificationGroup B G k)
    (fun j => ramIndex π (σ ^ p ^ j)) m n ?_ ?_ ?_ hne
  · intro x y hx hy
    simp only [pow_zero, pow_one] at hx hy
    rw [lowerRamificationGroup_eq_top_of_lt_ramIndex (A := A) huni hadj hgen hx,
      lowerRamificationGroup_eq_top_of_lt_ramIndex (A := A) huni hadj hgen hy]
  · intro j hj x y hx1 hx2 hy1 hy2
    have hx := (lowerRamificationGroup_inf_zpowers_eq_iff (A := A) p hp huni hadj hσ hord
      (show j + 1 ≤ m by omega) x).2 ⟨hx1, hx2⟩
    have hy := (lowerRamificationGroup_inf_zpowers_eq_iff (A := A) p hp huni hadj hσ hord
      (show j + 1 ≤ m by omega) y).2 ⟨hy1, hy2⟩
    rw [hgen, inf_top_eq] at hx hy
    rw [hx, hy]
  · rw [(ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord m).2 le_rfl]
    exact fun hc => (ENat.coe_ne_top (n + 1)) (top_le_iff.1 hc)

/-- ★★★**段 1 の本体** —— `G = ⟨σ⟩ ≅ ℤ/p^mℤ` かつ `G_n ≠ G_{n+1}` なら

`p^m ∣ Σ_{k=1}^{n} |G_k|`.

★これが原文の `φ_G(n) = (1/p^m)(n_0·p^m + n_1 p·p^{m−1} + ⋯) ∈ Z≥0` の中身である
(逸脱 1: 原文は `p` 進展開の値を出すが、ここは差分の形で整除性だけを出す)。

★段取り(§1 への代入):
* 境界 `bd 0 = 1`、`bd (j+1) = i_j`(= `(ramIndex π (σ^{p^j})).toNat`)、
* 値 `v j = p^{m−j}`、
* 第 `0` ブロック `[1, i_0)` は値 `p^m` なので「長さ × 値」は自明に `p^m` で割れる、
* 第 `j+1` ブロック `[i_j, i_{j+1})` は長さ `i_{j+1} − i_j` が `p^{j+1}` で割れ
  (Proposition 6.6 (iii) = `ramIndex_pow_pow_congr`)、値が `p^{m−(j+1)}` なので
  積は `p^m` で割れる。

★`ℕ∞` の引き算は使わない: `hfin` で `j ≤ l` の範囲の有限性を先に確保し、
`toNat` を取ってから `ℕ` の中で引いている(`lean-idioms.md` #102)。 -/
theorem pow_dvd_sum_natCard_lowerRamificationGroup {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (hπ : Irreducible π)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    (hgen : Subgroup.zpowers σ = ⊤) {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    p ^ m ∣ ∑ k ∈ Finset.Icc 1 n, Nat.card (lowerRamificationGroup B G k) := by
  have huni := (irreducible_iff_uniformizer π).mp hπ
  obtain ⟨l, hlm, hl⟩ := exists_ramIndex_eq_of_lowerRamificationGroup_ne (A := A) p hp huni hadj
    hσ hord hgen hne
  have hmonoE := monotone_ramIndex_pow_pow (A := A) p hp huni hadj hσ hord
  set f : ℕ → ℕ := fun j => (ramIndex π (σ ^ p ^ j)).toNat with hfdef
  have hfin : ∀ j, j ≤ l → ramIndex π (σ ^ p ^ j) = ((f j : ℕ) : ℕ∞) := by
    intro j hj
    refine (ENat.coe_toNat ?_).symm
    intro hc
    have hle := hmonoE hj
    simp only at hle
    rw [hc, hl] at hle
    exact (ENat.coe_ne_top (n + 1)) (top_le_iff.1 hle)
  have hfl : f l = n + 1 := by
    have h := hfin l le_rfl
    rw [hl] at h
    exact_mod_cast h.symm
  have hfmono : ∀ j, j + 1 ≤ l → f j ≤ f (j + 1) := by
    intro j hj
    have hle := hmonoE (show j ≤ j + 1 by omega)
    simp only at hle
    rw [hfin j (by omega), hfin (j + 1) hj] at hle
    exact_mod_cast hle
  have hf0 : 1 < f 0 := by
    have h1 : (1 : ℕ∞) < ramIndex π σ :=
      (mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj 1 σ).1 hσ
    rw [show σ = σ ^ p ^ 0 by rw [pow_zero, pow_one], hfin 0 (by omega)] at h1
    exact_mod_cast h1
  set bd : ℕ → ℕ := fun j => if j = 0 then 1 else f (j - 1) with hbddef
  have hbdsucc : ∀ j, bd (j + 1) = f j := by intro j; simp [hbddef]
  have hbd0 : bd 0 = 1 := by simp [hbddef]
  have key : p ^ m ∣ ∑ k ∈ Finset.Ico (bd 0) (bd (l + 1)),
      Nat.card (lowerRamificationGroup B G k) := by
    refine dvd_sum_Ico_of_blockwise (v := fun j => p ^ (m - j)) ?_ ?_ ?_
    · rintro (_ | j) hj
      · rw [hbd0, hbdsucc]; omega
      · rw [hbdsucc, hbdsucc]; exact hfmono j (by omega)
    · rintro (_ | j) hj k hk1 hk2
      · rw [hbd0] at hk1
        rw [hbdsucc] at hk2
        have hlt : (k : ℕ∞) < ramIndex π σ := by
          rw [show σ = σ ^ p ^ 0 by rw [pow_zero, pow_one], hfin 0 (by omega)]
          exact_mod_cast hk2
        rw [lowerRamificationGroup_eq_top_of_lt_ramIndex (A := A) huni hadj hgen hlt,
          Subgroup.card_top, ← Subgroup.card_top (G := G), ← hgen, Nat.card_zpowers, hord]
        simp
      · rw [hbdsucc] at hk1
        rw [hbdsucc] at hk2
        refine card_lowerRamificationGroup_of_mem_Ico (A := A) p hp huni hadj hσ hord hgen
          (show j + 1 ≤ m by omega) ?_ ?_
        · rw [hfin j (by omega)]; exact_mod_cast hk1
        · rw [hfin (j + 1) (by omega)]; exact_mod_cast hk2
    · rintro (_ | j) hj
      · rw [hbd0, hbdsucc]
        simp
      · rw [hbdsucc, hbdsucc]
        have hcong := ramIndex_pow_pow_congr (A := A) p hp hπ hadj hres hσ hord
          (hfin j (by omega)) (hfin (j + 1) (by omega))
        have hle := hfmono j (by omega)
        have hnat : p ^ (j + 1) ∣ f (j + 1) - f j := by
          have hz : ((p ^ (j + 1) : ℕ) : ℤ) ∣ ((f (j + 1) - f j : ℕ) : ℤ) := by
            push_cast [Nat.cast_sub hle]
            exact hcong
          exact_mod_cast hz
        have hsplit : p ^ m = p ^ (j + 1) * p ^ (m - (j + 1)) := by
          rw [← pow_add]; congr 1; omega
        rw [hsplit]
        exact mul_dvd_mul hnat dvd_rfl
  rw [hbd0, hbdsucc, hfl, Finset.Ico_add_one_right_eq_Icc] at key
  exact key

/-- ★★★★**Yoshida 2008 Theorem 6.11 (Hasse-Arf) の段 1** ——
`G ≅ ℤ/p^mℤ`(巡回)かつ `G_n ≠ G_{n+1}` なら `φ_G(n) ∈ ℤ≥0`。

原文 (Yoshida08 p.16):
> Theorem 6.11 (Hasse-Arf). If G is abelian, n ∈ Z[bb]_≥0 and G_n = G_n+1, then φ_G(n) ∈ Z[bb]_≥0.

★★逐語の `G_n = G_n+1` は `pdftotext` が `≠` の斜線を落とした形である。ここは `≠` で書いてある。

★`∈ ℤ≥0` は `∃ k : ℕ, φ_G(n) = (k : ℝ)` と書いた(非負性が型に入る)。
★段取り: 上の `p^m ∣ Σ_{k=1}^n |G_k|` に Lemma 6.10(i)(`herbrandPhiGroup_natCast`)を当て、
`|G| = p^m` で割るだけ。 -/
theorem exists_natCast_herbrandPhiGroup_of_isCyclic {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B] [FaithfulSMul G B]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p] {π : B} (hπ : Irreducible π)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    (hgen : Subgroup.zpowers σ = ⊤) {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    ∃ k : ℕ, herbrandPhiGroup G π (n : ℝ) = (k : ℝ) := by
  have huni := (irreducible_iff_uniformizer π).mp hπ
  obtain ⟨k, hk⟩ := pow_dvd_sum_natCard_lowerRamificationGroup (A := A) p hp hπ hadj hres hσ
    hord hgen hne
  refine ⟨k, ?_⟩
  have hcardG : Nat.card G = p ^ m := by
    rw [← Subgroup.card_top (G := G), ← hgen, Nat.card_zpowers, hord]
  have hpos : (0 : ℝ) < ((p ^ m : ℕ) : ℝ) := by exact_mod_cast Nat.pow_pos hp.pos
  rw [herbrandPhiGroup_natCast (A := A) huni hadj n, ← Nat.cast_sum, hk, hcardG,
    Nat.cast_mul, mul_comm, mul_div_assoc, div_self hpos.ne', mul_one]

/-! ## §3 段 3 —— `G ≠ G_1` の場合

原文 (Yoshida08 p.16-17):
> Now when G ̸= G1, set H = G1 and |G/H| = e0. As φG/H(n) = n/e0 for n ∈R≥0 by definition,
> by Lemma 6.10(ii) it suffices to show e0 | φH(n) when n ∈Z≥0 and Gn ̸= Gn+1
> (we know φH(n) ∈Z≥0). If n = 0 then φH(0) = 0. Let n > 0. For any i ∈Z≥1 (where Hi = Gi)
> with Hi ̸= Hi+1, we have e0 | i by Corollary 6.3, hence e0 | Σ_{i=1}^n |Hi|.
> As e0 and |H| are coprime, we have e0 | φH(n) by Lemma 6.10(i).

★原文の `hence e0 | Σ|Hi|` と `As e0 and |H| are coprime` は先行ノード
`AbelianJumpDivisibility.lean`(Y13)が `dvd_sum_natCard_lowerRamificationGroup` /
`coprime_natCard_quot_natCard_lowerRamificationGroup_one` /
`dvd_of_sum_eq_natCard_mul` として埋めてある。本節はそこへ Lemma 6.10(i) を差し込む。
★`n = 0` の場合(`φ_H(0) = 0`)は `Finset.Icc 1 0 = ∅` として自動的に含まれている。 -/

/-- `i ≥ 1` なら `G_i ≤ G_1` なので `|G_i ∩ G_1| = |G_i|`。

★Lemma 6.10(i) の部分群版 `herbrandPhi_natCast` は
`|(G_i).subgroupOf H|` の形で和を返すので、`H = G_1` のときこれで `G` 側の添字に戻す。
原文が `(where H_i = G_i)` と括弧で断っている一歩である。 -/
theorem natCard_subgroupOf_lowerRamificationGroup_one {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {i : ℕ} (hi : 1 ≤ i) :
    Nat.card ((lowerRamificationGroup B G i).subgroupOf (lowerRamificationGroup B G 1))
      = Nat.card (lowerRamificationGroup B G i) :=
  Nat.card_congr (Subgroup.subgroupOfEquivOfLe (lowerRamificationGroup_antitone B G hi)).toEquiv

/-- ★★**Lemma 6.10(i) を `Σ_{i=1}^n |G_i| = |G_1| · φ_{G_1}(n)` の形に読み替える**。

`φ_H(n) = k ∈ ℕ`(`H = G_1`)なら `Σ_{i=1}^n |G_i| = |G_1| · k`。
★これが Y13 の `dvd_of_sum_eq_natCard_mul` の受け渡し口 `hk` にそのまま入る。 -/
theorem sum_natCard_lowerRamificationGroup_eq_of_herbrandPhi {A B : Type*} [CommRing A]
    [CommRing B] [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [Fintype (lowerRamificationGroup B G 1)]
    {π : B} (huni : maximalIdeal B = Ideal.span {π})
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤) {n k : ℕ}
    (hk : herbrandPhi π (lowerRamificationGroup B G 1) (n : ℝ) = (k : ℝ)) :
    ∑ i ∈ Finset.Icc 1 n, Nat.card (lowerRamificationGroup B G i)
      = Nat.card (lowerRamificationGroup B G 1) * k := by
  have hcard : (0 : ℝ) < (Nat.card (lowerRamificationGroup B G 1) : ℝ) := by
    exact_mod_cast Nat.card_pos
  rw [herbrandPhi_natCast (A := A) huni hadj _ n,
    Finset.sum_congr rfl (fun i hi => by
      rw [natCard_subgroupOf_lowerRamificationGroup_one (B := B) (G := G)
        (Finset.mem_Icc.1 hi).1]),
    div_eq_iff hcard.ne'] at hk
  have hR : ((∑ i ∈ Finset.Icc 1 n, Nat.card (lowerRamificationGroup B G i) : ℕ) : ℝ)
      = ((Nat.card (lowerRamificationGroup B G 1) * k : ℕ) : ℝ) := by
    push_cast
    rw [hk]; ring
  exact_mod_cast hR

/-- ★★★**原文「it suffices to show e_0 | φ_H(n)」の "show" の部分**(段 3 の実質)。

`G` 可換・完全分岐(`G_0 = ⊤`)・`G_n ≠ G_{n+1}` で、`φ_{G_1}(n) = k ∈ ℕ` なら
`e_0 = |G/G_1|` は `k` を割る。★**仮定は原文が置いているものだけ**である
(Lemma 6.10(ii) も固定環も要らない)。

★段取り: Lemma 6.10(i)(上の補題)で `Σ_{i=1}^n |G_i| = |G_1| · k` にし、
Y13 の `dvd_of_sum_eq_natCard_mul`(= Corollary 6.3 + 互いに素)へ差し込むだけ。 -/
theorem natCard_quot_dvd_of_herbrandPhi_natCast {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [Finite G] [FaithfulSMul G B]
    [Fintype (lowerRamificationGroup B G 1)]
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p] {π : B} (hπ : Irreducible π)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x) {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) {k : ℕ}
    (hk : herbrandPhi π (lowerRamificationGroup B G 1) (n : ℝ) = (k : ℝ)) :
    Nat.card (G ⧸ lowerRamificationGroup B G 1) ∣ k :=
  dvd_of_sum_eq_natCard_mul (A := A) p ((irreducible_iff_uniformizer π).mp hπ) hπ.ne_zero hadj
    h0 habel hne
    (sum_natCard_lowerRamificationGroup_eq_of_herbrandPhi (A := A)
      ((irreducible_iff_uniformizer π).mp hπ) hadj hk)

/-- ★★**原文「As φ_{G/H}(n) = n/e_0 for n ∈ ℝ≥0 by definition」**(具体層)。

固定環 `C` の素元 `ϖ` による跳び指数 `i_ϖ` が `H` の上で `⊤`、`H` の外で `1` なら
`φ_{G/H}(x) = x/e_0`。★§1 の `phiOf_of_two_valued` に代入するだけ。

★★**原文の `by definition` は正確ではない**: `i_ϖ` が `H` の外で `1` であること
(= `(G/H)_1 = 1`、順分岐)は Proposition 6.9 の帰結であって定義ではない。
本補題はその 2 つの入力を仮定として明示している。 -/
theorem herbrandPhiGroup_eq_div_natCard_quot {G : Type*} [Group G] [Fintype G] {C : Type*}
    [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    (H : Subgroup G) {ϖ : C} (htop : ∀ σ : G, σ ∈ H → ramIndex ϖ σ = ⊤)
    (hone : ∀ σ : G, σ ∉ H → ramIndex ϖ σ = 1) {x : ℝ} (hx : 0 ≤ x) :
    herbrandPhiGroup G ϖ x = x / (Nat.card (G ⧸ H) : ℝ) := by
  classical
  have hTcard : (Finset.univ.filter (fun σ : G => σ ∈ H)).card = Nat.card H := by
    rw [Nat.card_eq_fintype_card, Fintype.card_subtype]
  have hGH : (Nat.card G : ℝ) = (Nat.card (G ⧸ H) : ℝ) * (Nat.card H : ℝ) := by
    rw [← Nat.cast_mul, ← Subgroup.card_eq_card_quotient_mul_card_subgroup]
  have hHpos : (Nat.card H : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Nat.card_pos.ne'
  have hQpos : (Nat.card (G ⧸ H) : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Nat.card_pos.ne'
  rw [herbrandPhiGroup, phiOf_of_two_valued _ (Finset.univ.filter (fun σ : G => σ ∈ H))
    (fun s hs => htop s (Finset.mem_filter.1 hs).2)
    (fun s hs => hone s (fun hc => hs (Finset.mem_filter.2 ⟨Finset.mem_univ s, hc⟩))) hx,
    hTcard, hGH]
  field_simp

/-- ★★★**段 3 の組み立て** —— `φ_G(n) = φ_{G_1}(n)/e_0` が分かっていれば
`φ_G(n) ∈ ℤ≥0`。

★`hquot` が原文の「by Lemma 6.10(ii)」に当たる一歩(逸脱 4)。
実質部分(`e_0 ∣ φ_{G_1}(n)`)は上の `natCard_quot_dvd_of_herbrandPhi_natCast` が
仮定なしで出している。 -/
theorem exists_natCast_herbrandPhiGroup_of_dvd_quot {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [Fintype G] [FaithfulSMul G B]
    [Fintype (lowerRamificationGroup B G 1)]
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p] {π : B} (hπ : Irreducible π)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x) {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) {k : ℕ}
    (hk : herbrandPhi π (lowerRamificationGroup B G 1) (n : ℝ) = (k : ℝ))
    (hquot : herbrandPhiGroup G π (n : ℝ)
      = herbrandPhi π (lowerRamificationGroup B G 1) (n : ℝ)
          / (Nat.card (G ⧸ lowerRamificationGroup B G 1) : ℝ)) :
    ∃ k' : ℕ, herbrandPhiGroup G π (n : ℝ) = (k' : ℝ) := by
  obtain ⟨k', hk'⟩ := natCard_quot_dvd_of_herbrandPhi_natCast (A := A) p hπ hadj h0 habel hne hk
  refine ⟨k', ?_⟩
  have he0 : (0 : ℝ) < (Nat.card (G ⧸ lowerRamificationGroup B G 1) : ℝ) := by
    exact_mod_cast Nat.card_pos
  rw [hquot, hk, hk']
  push_cast
  field_simp

/-- ★★★★**Yoshida 2008 Theorem 6.11 (Hasse-Arf) の段 3** ——
`G ≠ G_1` の場合。`H = G_1` として `φ_G(n) ∈ ℤ≥0`。

原文 (Yoshida08 p.16):
> Theorem 6.11 (Hasse-Arf). If G is abelian, n ∈ Z[bb]_≥0 and G_n = G_n+1, then φ_G(n) ∈ Z[bb]_≥0.

★★逐語の `G_n = G_n+1` は `pdftotext` が `≠` の斜線を落とした形である。ここは `≠` で書いてある。

仮定は原文の 3 つ:
* `hk` : `φ_H(n) ∈ ℤ≥0`(原文の「we know φ_H(n) ∈ Z≥0」。段 1-2 の帰結、逸脱 5)、
* `hcomp` : `φ_G = φ_{G/H} ∘ φ_H`(Lemma 6.10(ii)、逸脱 4)、
* `htopC` / `honeC` : `i_ϖ` が `H` の上で `⊤`、外で `1`(原文の「φ_{G/H}(n) = n/e_0」)。

★★可換性 `habel` を落とすと偽である(Hasse-Arf はアーベルでのみ成り立つ)。
`habel` は Corollary 6.3 経由で `e_0 ∣ i`(跳びの位置)に効いている。 -/
theorem exists_natCast_herbrandPhiGroup_of_tame_quotient {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [Fintype G] [FaithfulSMul G B]
    [Fintype (lowerRamificationGroup B G 1)]
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p] {π : B} (hπ : Irreducible π)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x) {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) {k : ℕ}
    (hk : herbrandPhi π (lowerRamificationGroup B G 1) (n : ℝ) = (k : ℝ)) {ϖ : C}
    (hcomp : herbrandPhiGroup G π (n : ℝ)
      = herbrandPhiGroup G ϖ (herbrandPhi π (lowerRamificationGroup B G 1) (n : ℝ)))
    (htopC : ∀ σ : G, σ ∈ lowerRamificationGroup B G 1 → ramIndex ϖ σ = ⊤)
    (honeC : ∀ σ : G, σ ∉ lowerRamificationGroup B G 1 → ramIndex ϖ σ = 1) :
    ∃ k' : ℕ, herbrandPhiGroup G π (n : ℝ) = (k' : ℝ) := by
  refine exists_natCast_herbrandPhiGroup_of_dvd_quot (A := A) p hπ hadj h0 habel hne hk ?_
  rw [hcomp, herbrandPhiGroup_eq_div_natCard_quot _ htopC honeC
    (by rw [hk]; exact Nat.cast_nonneg k)]

/-! ## §4 原典の `∈ ℤ≥0` の `≥0` を明示する系 -/

/-- ★段 1 の結論の非負性(原文の `φ_G(n) ∈ Z≥0` の `≥0`)。 -/
theorem herbrandPhiGroup_nonneg_of_isCyclic {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B] [FaithfulSMul G B]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p] {π : B} (hπ : Irreducible π)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    (hgen : Subgroup.zpowers σ = ⊤) {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    0 ≤ herbrandPhiGroup G π (n : ℝ) :=
  nonneg_of_exists_natCast (exists_natCast_herbrandPhiGroup_of_isCyclic (A := A) p hp hπ hadj
    hres hσ hord hgen hne)

end ABC3.Found.PGC
