import ABC3.Found.PGC.LubinTateRamificationBookkeeping
import ABC3.Found.PGC.RamificationImageUnits

/-!
# Yoshida 2008 Proposition 6.14 の「鋭い形」——  `G^i = ρ_(f,m)^(-1)(1 + p^i)`

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Proposition 6.14(物理 p.17)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-6.html` の `#prop-6-14`。

原文 (Yoshida08 p.17):
> Proposition 6.14. Let x ∈ K^× with v(x) = n > 0. Let L = K_n and K^m_x as in Definition 5.3. Then we have Gal(K^m_x/L)^m = {id} for all m ≥ 1 (see Definition 6.12).

## 本ノードが埋めたもの ——  命題文ではなく「証明の中身」

命題文そのもの(上付き分岐群が `m` で**消える**)は木に既にある
(`Found/PGC/LubinTateRamificationBookkeeping.lean` の
`upperRamificationGroup_torsionGen_eq_bot`)。本ファイルが埋めたのは
**その証明の途中で原典が立てている、より強い等式**である。原典の逐語:

    Thus for G = Gal(K^m_x/L) and 1 ≤ i ≤ m, we have |G_n| = |ρ^{-1}_{f,m}(1+p^i)| = q^{m−i}
    for q^{i−1}−1 < n ≤ q^i−1. Thus φ_G(q^m−1) = (1/|G|) Σ_{i=1}^{q^m−1} |G_i|
    = (1/((q−1)q^{m−1})) (Σ_{i=1}^{m}(q^i−q^{i−1})·q^{m−i}) = m and G^m = G_{q^m−1} = {id}.

原典は `φ_G` を `q^m − 1` **1 点でしか**評価していないが、同じ計算がそのまま
`q^I − 1`(`0 ≤ I ≤ m`)で走る。すなわち

    φ_G(q^I − 1) = I           (`herbrandPhiGroup_pow_sub_one`)
    G^I = G_{q^I − 1}          (`upperRamificationGroup_eq_lowerRamificationGroup_pow_sub_one`)
    G^I = ρ_{f,m}^{-1}(1 + p^I) (`upperRamificationGroup_..._eq_comap_map_principalUnits`)

であり、`I = m` が原典の `G^m = {id}` である。★**これが「鋭い形」**で、
下付きの `q` 冪刻みの添字が Herbrand を通ると上付きの **1 刻み**に化ける。

## 埋めた宣言(★正直な線引き)

**埋めたもの**

* §0 **抽象核(純算術)**: `sum_Ico_pow_blocks_le` / `sum_Icc_pow_blocks_le`
  —— `q` 進ブロックの**部分和**。分岐・付値・群の語彙が 1 語も出てこない。
* §1 **抽象核(抽象分岐)**: `herbrandPhiGroup_pow_sub_one` /
  `upperRamificationGroup_eq_lowerRamificationGroup_pow_sub_one` /
  `upperRamificationGroup_eq_of_lowerRamificationGroup_blocks`
  —— 任意の離散付値環 `B` と有限群 `G`。Lubin-Tate も `p` 進体も出てこない。
* §2 **具体層**: `upperRamificationGroup_iteratedLubinTatePsi_eq_comap_map_principalUnits`
  (`1 ≤ I ≤ M+1`)と、`I = 0` を含む `_of_le` 版、素元を `torsionGen` に固定した版。
  ★**追加仮定はゼロ**(段 1 は `LubinTateRamificationBookkeeping.lean` が供給済み)。
* §3 **相互律の両立性**: `reciprocityMap_psiGenSeq_eq_mk_reciprocityUnits`
  —— 有限段の相互律 `ρ_{f,m}` が、極限の相互律 `Art : Γ_K → 𝒪_K^×` の
  `U^m` を法とする還元であること。★これが「無限次への持ち上げ」の要である。
  併せて `algEquivRestrictSelfHom`(制限写像を群準同型として束ねたもの)。
* §4 **`Γ_K` への持ち上げ**:
  `comap_upperRamificationGroup_eq_absGalPrincipalLevel` ——

      `res_{m}^{-1}(Gal(K_{f,m}/K)^I) = Art^{-1}(U^I_K)`   (`I ≤ m`)

  ★これが消費先 `Found/PGC/RamificationImageUnits.lean` の
  `F.S (absGalPrincipalLevel (n+i)) v = absGalPrincipalLevel n` の**中身**である。

**埋めていないもの(★次のノード)**

1. ★★**段データ `F : StageFiltration K.absGal` 自身は本ファイルでは作らない。**
   `RamificationFiltrationBuild.lean` の `compat` が残っているためで、
   本ファイルは `F.S N v` の**中身にあたる部分群の等式**を、`F` を経由せずに
   `Subgroup.comap (algEquivRestrictSelfHom …)` の形で出している。
   `F` が立った時点で、`F.S (absGalPrincipalLevel (n+i)) (n:ℝ)` がこの
   `comap` に等しいことを言えば消費先の仮定 2 本が閉じる。
2. ★**実数の添字 `v`。** 本ファイルの `I` は自然数である。実数 `v` について
   `G^v = G^{⌈v⌉}` を言うには跳びの整数性(Hasse-Arf)が要る。★別ノード。
3. ★**`n = 1` の場合のみ。** 原典の `L = K_n` は本ファイルでは `L = K` である
   (決定 D29。`LubinTateUpperRamificationVanish.lean` の逸脱 1 を引き継ぐ)。

## ★添字をどう決めたか(★自分で確かめた)

原典の下付きの範囲は `q^{i−1} − 1 < n ≤ q^i − 1`、すなわち ℕ の切り詰め引き算を
出さない形で `q^{i−1} ≤ n < q^i` である。木の在庫
`lowerRamificationGroupAdjoin_eq_comap_map_principalUnits` は原典の `i` を `i+1` に
ずらして `q^i ≤ n < q^{i+1}` で `G_n = ρ^{-1}(U^{i+1})` と書いてある。本ファイルの
`upperRamificationGroup_eq_of_lowerRamificationGroup_blocks` はこの「0 始まりの `i`」を
受け取り、`I = i+1`(＝原典の `i`)で結論する。したがって
**上付きの添字 `I` は原典の `i` と一致し、ずれは入らない。**

★★**独立な検算を 1 本入れた**: `ker_algEquivRestrictSelfHom_eq_absGalPrincipalLevel`。
`I = M+1` の場合に `G^{M+1} = {id}`(在庫)を代入すると

    `ker(Γ_K → Gal(K_{f,M+1}/K)) = Art^{-1}(U^{M+1}_K)`

が出る。左辺は Galois 対応から `Gal(K̄/K_{f,M+1})` であり、右辺は消費先の
`absGalPrincipalLevel (M+1)` である。★**両者が一致した**ので、消費先
`RamificationImageUnits.lean` が「レベル `m+1` と段 `m+1` はずれない」と
書いていることは正しい。★見立てどおりだった(前の波の判断を追試した)。

## ★原典より短い道(名指し)

原典は `φ_G(q^m − 1) = m` を**1 点**で計算し、そこから `G^m = {id}` だけを取り出す。
本ファイルは同じ計算を `q^I − 1` で回すだけで **`0 ≤ I ≤ m` のすべての段**を得た。
★**新しい数学は 1 行も要らなかった** —— 必要だったのは、既存の
`sum_Ico_pow_blocks`(全区間)を**部分区間**に一般化することだけである
(`sum_Ico_pow_blocks_le`。重み `K` を持ち歩く代わりに、上端の `I` を動かして
`M` を固定する形にした方が短い)。

★もう 1 つ: §3 の両立性は Serre XV §2 の Herbrand 経由を通らず、
`unitsEquivCompatibleUnits` の `apply_symm_apply` 1 本と
`principalUnitsQuotientEquiv_apply_mk` だけで出る(4 行)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★**`n = 1` の場合のみ**(上記、決定 D29)。原典の `K^m_x` は本ファイルでは
   `K(x)`(`x` は `ψ_{M+1}` の原始根)である。
2. ★**`m = M + 1 ≥ 1`** と書く形で `m ≥ 1` を埋め込んである(`m = 0` では
   拡大が自明になり主張が空虚になる)。
3. ★**`I = 0` も認めた**(`_of_le` 版)。原典の範囲は `1 ≤ i ≤ m` だが、
   `I = 0` では両辺とも `⊤` になって成り立つ(`upperRamificationGroup_zero` と
   `principalUnits_zero_eq_top`)。★**原典の主張を弱めても強めてもいない。**
4. ★**既存の `Found/PGC/*.lean` は 1 行も書き換えていない**(import のみ。D28)。
5. ★本ファイルは Prop 6.14 の**証明の中身**を扱うので、主定理には `.src` を
   付けて Prop 6.14 に紐づけてある。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open IsLocalRing
open scoped NormedField Valued Classical

/-! ## §0 抽象核(純算術)—— `q` 進ブロックの**部分和**

★分岐・付値・群・Galois の語彙が 1 つも出てこない。 -/

/-- ★★**抽象核(純算術)** —— `c` が `[q^i, q^{i+1})` の上で定数 `q^j`(`i + j = M`)なら、
**上端を `q^I` で切った部分和**は

    `∑_{n ∈ [1, q^I)} c(n) = I·(q^{M+1} − q^M)`   (`0 ≤ I ≤ M+1`)

(`q = r+1`、`q^{M+1} − q^M = r·q^M` と書いて**切り詰め引き算を消してある**)。

★在庫 `sum_Ico_pow_blocks`(`Found/PGC/LubinTateUpperRamificationVanish.lean`)は
`I = M+1` の場合である。★**重み `K` を持ち歩く代わりに上端 `I` を動かす**方が短い:
帰納法の 1 段で切り落とすブロックが `[q^I, q^{I+1})` ただ 1 つになり、
そこでの `c` の値 `q^j`(`I + j = M`)がそのまま `q^I·q^j = q^M` に潰れる。 -/
theorem sum_Ico_pow_blocks_le (r M : ℕ) :
    ∀ (I : ℕ) (c : ℕ → ℕ), I ≤ M + 1 →
      (∀ i j : ℕ, i + j = M → ∀ n ∈ Finset.Ico ((r + 1) ^ i) ((r + 1) ^ (i + 1)),
          c n = (r + 1) ^ j) →
      ∑ n ∈ Finset.Ico 1 ((r + 1) ^ I), c n = I * (r * (r + 1) ^ M) := by
  intro I
  induction I with
  | zero => intro c _ _; simp
  | succ I ih =>
      intro c hI hc
      obtain ⟨j, hj⟩ : ∃ j, I + j = M := ⟨M - I, by omega⟩
      have h1 : 1 ≤ (r + 1) ^ I := Nat.one_le_pow _ _ (Nat.succ_pos r)
      have h2 : (r + 1) ^ I ≤ (r + 1) ^ (I + 1) :=
        Nat.pow_le_pow_right (Nat.succ_pos r) (by omega)
      have hsplit := Finset.sum_Ico_consecutive c h1 h2
      have hbot := ih c (by omega) hc
      have hcard : (r + 1) ^ (I + 1) - (r + 1) ^ I = r * (r + 1) ^ I := by
        have hpow : (r + 1) ^ (I + 1) = r * (r + 1) ^ I + (r + 1) ^ I := by ring
        omega
      have hpowM : (r + 1) ^ I * (r + 1) ^ j = (r + 1) ^ M := by rw [← pow_add, hj]
      have htop : ∑ n ∈ Finset.Ico ((r + 1) ^ I) ((r + 1) ^ (I + 1)), c n = r * (r + 1) ^ M := by
        rw [Finset.sum_congr rfl (fun n hn => hc I j hj n hn), Finset.sum_const, Nat.card_Ico,
          smul_eq_mul, hcard, mul_assoc, hpowM]
      rw [← hsplit, hbot, htop]
      ring

/-- ★**原典の `Σ_{i=1}^{q^I−1}` の形**(区間が `Icc 1 (q^I − 1)`)。

`Finset.Icc 1 (q^I − 1) = Finset.Ico 1 (q^I)` と読み替えて `sum_Ico_pow_blocks_le` に
渡すだけ。★`I = 0` でも成り立つ(両辺 `0`)。 -/
theorem sum_Icc_pow_blocks_le (r M I : ℕ) (hI : I ≤ M + 1) (c : ℕ → ℕ)
    (hc : ∀ i j n : ℕ, i + j = M → (r + 1) ^ i ≤ n → n < (r + 1) ^ (i + 1) →
      c n = (r + 1) ^ j) :
    ∑ n ∈ Finset.Icc 1 ((r + 1) ^ I - 1), c n = I * (r * (r + 1) ^ M) := by
  have h1 : 1 ≤ (r + 1) ^ I := Nat.one_le_pow _ _ (Nat.succ_pos r)
  have hset : Finset.Icc 1 ((r + 1) ^ I - 1) = Finset.Ico 1 ((r + 1) ^ I) := by
    ext n; simp only [Finset.mem_Icc, Finset.mem_Ico]; omega
  rw [hset]
  exact sum_Ico_pow_blocks_le r M I c hI (by
    intro i j hij n hn
    rw [Finset.mem_Ico] at hn
    exact hc i j n hij hn.1 hn.2)

/-! ## §1 抽象核(抽象分岐)—— Herbrand を `q^I − 1` で評価する

★Lubin-Tate も `p` 進体も出てこない。任意の離散付値環 `B` と有限群 `G` の話である。 -/

section AbstractCore

variable {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
  [IsDiscreteValuationRing B]
variable {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]

/-- ★★★**抽象核** —— 原典の `φ_G(q^m − 1) = m` を**すべての `I ≤ m` で**。

`|G| = r·q^M`(`q = r+1`)で、下付き分岐群の位数がブロックごとに
`|G_n| = q^j`(`q^i ≤ n < q^{i+1}`, `i + j = M`)なら

    `φ_G(q^I − 1) = I`   (`0 ≤ I ≤ M+1`)

★段取りは原典どおり: Lemma 6.10(i)(在庫 `herbrandPhiGroup_natCast`)で
`φ_G(N) = (1/|G|)·Σ_{n=1}^{N} |G_n|` に直し、分子を §0 の `sum_Icc_pow_blocks_le` で
`I·(r·q^M)` に潰し、分母 `|G| = r·q^M` で割る。

★★**`r > 0`(＝剰余体が 2 元以上)を仮定に置いていない** —— `Nat.card G > 0` から出る。 -/
theorem herbrandPhiGroup_pow_sub_one
    {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (r M : ℕ)
    (hGcard : Nat.card G = r * (r + 1) ^ M)
    (hcard : ∀ i j n : ℕ, i + j = M → (r + 1) ^ i ≤ n → n < (r + 1) ^ (i + 1) →
      Nat.card (lowerRamificationGroup B G n) = (r + 1) ^ j)
    (I : ℕ) (hI : I ≤ M + 1) :
    herbrandPhiGroup G α ((((r + 1) ^ I - 1 : ℕ) : ℝ)) = (I : ℝ) := by
  have hGpos : 0 < Nat.card G := Nat.card_pos
  have hpowM : 1 ≤ (r + 1) ^ M := Nat.one_le_pow _ _ (Nat.succ_pos r)
  have hrpos : 0 < r := by
    rcases Nat.eq_zero_or_pos r with h | h
    · rw [h, Nat.zero_mul] at hGcard; omega
    · exact h
  have hden : 0 < r * (r + 1) ^ M := Nat.mul_pos hrpos hpowM
  rw [herbrandPhiGroup_natCast (A := A) huni hadj]
  have hsum := sum_Icc_pow_blocks_le r M I hI _ hcard
  have hcast : (∑ n ∈ Finset.Icc 1 ((r + 1) ^ I - 1),
      (Nat.card (lowerRamificationGroup B G n) : ℝ))
      = ((I * (r * (r + 1) ^ M) : ℕ) : ℝ) := by
    rw [← hsum]; push_cast; rfl
  rw [hcast, hGcard]
  have hne : ((r * (r + 1) ^ M : ℕ) : ℝ) ≠ 0 := Nat.cast_ne_zero.2 (by omega)
  push_cast
  push_cast at hne
  field_simp

/-- ★★**抽象核** —— `G^I = G_{q^I − 1}`(`0 ≤ I ≤ M+1`)。

`φ_G(q^I − 1) = I` を `G^{φ_G(n)} = G_n`(在庫 `upperRamificationGroup_herbrandPhiGroup`)に
差し込むだけ。★原典が `I = M+1` でしか書いていない段である。 -/
theorem upperRamificationGroup_eq_lowerRamificationGroup_pow_sub_one
    {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (r M : ℕ)
    (hGcard : Nat.card G = r * (r + 1) ^ M)
    (hcard : ∀ i j n : ℕ, i + j = M → (r + 1) ^ i ≤ n → n < (r + 1) ^ (i + 1) →
      Nat.card (lowerRamificationGroup B G n) = (r + 1) ^ j)
    (I : ℕ) (hI : I ≤ M + 1) :
    upperRamificationGroup G α ((I : ℕ) : ℝ) = lowerRamificationGroup B G ((r + 1) ^ I - 1) := by
  have hphi := herbrandPhiGroup_pow_sub_one (A := A) huni hadj r M hGcard hcard I hI
  calc upperRamificationGroup G α ((I : ℕ) : ℝ)
      = upperRamificationGroup G α
          (herbrandPhiGroup G α ((((r + 1) ^ I - 1 : ℕ) : ℝ))) := by rw [hphi]
    _ = ramificationGroupReal α ((((r + 1) ^ I - 1 : ℕ) : ℝ)) :=
        upperRamificationGroup_herbrandPhiGroup _ _
    _ = lowerRamificationGroup B G ((r + 1) ^ I - 1) :=
        ramificationGroupReal_natCast (A := A) huni hadj _

/-- ★★★★**抽象核 —— 鋭い形の心臓部**。

下付き分岐群がブロックごとに部分群 `V (i+1)` に等しい(`q^i ≤ n < q^{i+1}`)なら、
**上付き分岐群はその `V` を 1 刻みで拾う**:

    `G^I = V I`   (`1 ≤ I ≤ M+1`)

★`q^I − 1` が属するブロックは `i = I−1` であり(`q^{I−1} ≤ q^I − 1 < q^I`、
`q ≥ 2` が要る)、そこでの値が `V ((I−1)+1) = V I` である。
★`q ≥ 2` は `Nat.card G = r·q^M > 0` から出るので**仮定に置いていない**。
★分岐の語彙以外(Lubin-Tate・付値・体)は 1 つも出てこない。 -/
theorem upperRamificationGroup_eq_of_lowerRamificationGroup_blocks
    {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (r M : ℕ)
    (hGcard : Nat.card G = r * (r + 1) ^ M)
    (hcard : ∀ i j n : ℕ, i + j = M → (r + 1) ^ i ≤ n → n < (r + 1) ^ (i + 1) →
      Nat.card (lowerRamificationGroup B G n) = (r + 1) ^ j)
    (V : ℕ → Subgroup G)
    (hblock : ∀ i n : ℕ, i ≤ M → (r + 1) ^ i ≤ n → n < (r + 1) ^ (i + 1) →
      lowerRamificationGroup B G n = V (i + 1))
    (I : ℕ) (h1 : 1 ≤ I) (h2 : I ≤ M + 1) :
    upperRamificationGroup G α ((I : ℕ) : ℝ) = V I := by
  obtain ⟨I', rfl⟩ : ∃ I', I = I' + 1 := ⟨I - 1, by omega⟩
  have hGpos : 0 < Nat.card G := Nat.card_pos
  have hpowM : 1 ≤ (r + 1) ^ M := Nat.one_le_pow _ _ (Nat.succ_pos r)
  have hrpos : 0 < r := by
    rcases Nat.eq_zero_or_pos r with h | h
    · rw [h, Nat.zero_mul] at hGcard; omega
    · exact h
  have hp : (1 : ℕ) ≤ (r + 1) ^ I' := Nat.one_le_pow _ _ (Nat.succ_pos r)
  have hexp : (r + 1) ^ (I' + 1) = r * (r + 1) ^ I' + (r + 1) ^ I' := by ring
  have hrp : 0 < r * (r + 1) ^ I' := Nat.mul_pos hrpos hp
  rw [upperRamificationGroup_eq_lowerRamificationGroup_pow_sub_one (A := A) huni hadj r M
    hGcard hcard (I' + 1) h2]
  exact hblock I' _ (by omega) (by omega) (by omega)

end AbstractCore

/-! ## §2 具体層 —— `K(µ_{f,m})` への代入(★追加仮定ゼロ) -/

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

def upperRamificationGroup_iteratedLubinTatePsi_eq_comap_map_principalUnits.src :
    ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Proposition 6.14", sectionId := "prop-6-14" }

/-- ★★★★★★★★**Yoshida 2008 Proposition 6.14 の「鋭い形」(`n = 1` の場合)**。

原文 (Yoshida08 p.17):
> Proposition 6.14. Let x ∈ K^× with v(x) = n > 0. Let L = K_n and K^m_x as in Definition 5.3. Then we have Gal(K^m_x/L)^m = {id} for all m ≥ 1 (see Definition 6.12).

原典の証明の途中式(逐語):

    Thus for G = Gal(K^m_x/L) and 1 ≤ i ≤ m, we have |G_n| = |ρ^{-1}_{f,m}(1+p^i)| = q^{m−i}
    for q^{i−1}−1 < n ≤ q^i−1.

を上付きに直した形である。`G = Gal(K(x)/K)`(`x` は `ψ_{M+1}` の原始根、`m = M+1`)、
`ρ_{f,m}` は Lubin-Tate 相互律 `galoisUnitReciprocityEquiv`、`1 + p^I` は
`principalUnits K π I` の `U^{M+1}` を法とする像であり、結論は

    `G^I = ρ_{f,m}^{-1}(U^I / U^{M+1})`   (`1 ≤ I ≤ M+1`)。

★`I = M+1` を代入すると右辺は `ρ^{-1}(1) = {id}` になり、原典の命題文
`upperRamificationGroup_torsionGen_eq_bot` に一致する。

★★**追加仮定はゼロ**である。段 1(`G_n = ρ^{-1}(U^{i+1})`)は
`Found/PGC/LubinTateRamificationBookkeeping.lean` の
`lowerRamificationGroupAdjoin_eq_comap_map_principalUnits` が供給し、
段 2(位数 `q^{m−i}`)は `Found/PGC/LubinTateUpperRamificationVanish.lean` の
`natCard_map_principalUnits` が供給する。本ファイルは段 3(Herbrand)を
`I` すべてに一般化しただけである。

★逸脱: `n = 1` の場合のみ(冒頭「逸脱の記録 1」)。 -/
theorem upperRamificationGroup_iteratedLubinTatePsi_eq_comap_map_principalUnits
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
    (I : ℕ) (h1 : 1 ≤ I) (h2 : I ≤ M + 1) :
    upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      α (I : ℝ)
      = Subgroup.comap (galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
            hxψ hxn hmem).toMonoidHom
          ((principalUnits K π I).map
            (QuotientGroup.mk' (principalUnits K π (M + 1)))) := by
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
  have hblock : ∀ i n : ℕ, i ≤ M → (r + 1) ^ i ≤ n → n < (r + 1) ^ (i + 1) →
      lowerRamificationGroup (adjoinIntegers K x)
          ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) n
        = Subgroup.comap (galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
              hxψ hxn hmem).toMonoidHom
            ((principalUnits K π (i + 1)).map
              (QuotientGroup.mk' (principalUnits K π (M + 1)))) := by
    intro i n hi hle hlt
    exact lowerRamificationGroupAdjoin_eq_comap_map_principalUnits K hq hπmax hπne0 f hf0 hf1 hf
      M hn x hxψ hxn hmem i n hi (by rw [hr]; exact hle) (by rw [hr]; exact hlt)
  refine upperRamificationGroup_eq_of_lowerRamificationGroup_blocks (A := 𝒪[K.carrier]) huni hadj
    r M hGcard ?_
    (fun J => Subgroup.comap (galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
        hxψ hxn hmem).toMonoidHom
      ((principalUnits K π J).map (QuotientGroup.mk' (principalUnits K π (M + 1)))))
    hblock I h1 h2
  intro i j n hij hle hlt
  rw [hblock i n (by omega) hle hlt, natCard_comap_mulEquiv]
  have hidx : i + 1 + j = M + 1 := by omega
  have hc := natCard_map_principalUnits K hq hπmax hπne0 i j
  rw [hidx, hr] at hc
  exact hc

/-- ★`I = 0` を含めた形。`I = 0` では両辺とも `⊤` である
(`upperRamificationGroup_zero` と `principalUnits_zero_eq_top`)。

★逸脱: 原典の範囲は `1 ≤ i ≤ m` である(冒頭「逸脱の記録 3」)。 -/
theorem upperRamificationGroup_iteratedLubinTatePsi_eq_comap_map_principalUnits_of_le
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
    (I : ℕ) (h2 : I ≤ M + 1) :
    upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      α (I : ℝ)
      = Subgroup.comap (galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
            hxψ hxn hmem).toMonoidHom
          ((principalUnits K π I).map
            (QuotientGroup.mk' (principalUnits K π (M + 1)))) := by
  rcases Nat.eq_zero_or_pos I with rfl | h1
  · rw [principalUnits_zero_eq_top K π,
      Subgroup.map_top_of_surjective _ (QuotientGroup.mk'_surjective _),
      Subgroup.comap_top, Nat.cast_zero, upperRamificationGroup_zero huni]
  · exact upperRamificationGroup_iteratedLubinTatePsi_eq_comap_map_principalUnits K hq hπmax hπne0
      f hf0 hf1 hf M hn x hxψ hxn hmem huni I h1 h2

/-- ★**素元を捩れ点自身に取った形**(仮定 `huni` が消える)。

`irreducible_torsionGen`(在庫)が `x` 自身を `𝒪_{K(x)}` の素元だと言う。
★★**残っているのは原典の設定だけである。** -/
theorem upperRamificationGroup_torsionGen_eq_comap_map_principalUnits
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
    (I : ℕ) (h2 : I ≤ M + 1) :
    upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      (torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1) x hxn hmem) (I : ℝ)
      = Subgroup.comap (galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x
            hxψ hxn hmem).toMonoidHom
          ((principalUnits K π I).map
            (QuotientGroup.mk' (principalUnits K π (M + 1)))) :=
  upperRamificationGroup_iteratedLubinTatePsi_eq_comap_map_principalUnits_of_le K hq hπmax hπne0
    f hf0 hf1 hf M hn x hxψ hxn hmem
    ((IsDiscreteValuationRing.irreducible_iff_uniformizer _).mp
      (irreducible_torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1) hn x hxψ hxn hmem)) I h2

/-! ## §3 有限段の相互律は極限の相互律の還元である

★これが「有限段の鋭い形」を `Γ_K` へ持ち上げる唯一の道具である。
★Serre XV §2 の Herbrand 経由は通らない —— 4 行で済む。 -/

section Compat

variable {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))

/-- ★★★**有限段の相互律 = 極限の相互律の還元**。

`Art : Γ_K → 𝒪_K^×`(在庫 `reciprocityUnits`)を `U^{m}` で割ったものが、
第 `m` 段の相互律 `ρ_{f,m}` に一致する(`m = M+1`、生成元は `psiGenSeq M`)。

★段取り: `unitsEquivCompatibleUnits` の `apply_symm_apply` で
「`Art σ` の第 `m` 成分 = `reciprocityMapLimitFamily σ m`」を取り、
`principalUnitsQuotientEquiv_apply_mk` で商群側に移して単射性で割る。
★**Herbrand も分岐も出てこない。** -/
theorem reciprocityMap_psiGenSeq_eq_mk_reciprocityUnits
    (m : ℕ) (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ
      = QuotientGroup.mk (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ) := by
  have h : unitsEquivCompatibleUnits K hπmax
      (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ)
      = reciprocityMapLimit K hq hπmax hπne0 f hf0 hf1 hf σ := MulEquiv.apply_symm_apply _ _
  have hkey : (Units.map ((Ideal.Quotient.mk
        (Ideal.span ({π ^ (m + 1)} : Set (𝒪[K.carrier])))).toMonoidHom))
      (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ)
      = reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ (m + 1) :=
    congrFun (congrArg Subtype.val h) (m + 1)
  apply (principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)).injective
  rw [principalUnitsQuotientEquiv_apply_mk]
  show _ = (Units.map ((Ideal.Quotient.mk
        (Ideal.span ({π ^ (m + 1)} : Set (𝒪[K.carrier])))).toMonoidHom))
      (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ)
  rw [hkey]
  rfl

variable (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]

/-- `galoisUnitReciprocityEquiv` の値は `galoisUnitReciprocityMap` そのもの
(`MulEquiv.ofBijective` の外し方)。 -/
theorem galoisUnitReciprocityEquiv_apply
    (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem y
      = galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem y := rfl

/-- `ρ_{f,m}(σ|_{K(x)}) = reciprocityMap σ` —— 在庫
`galoisUnitReciprocityMap_eq_reciprocityMap` の同型版。 -/
theorem galoisUnitReciprocityEquiv_algEquivRestrictSelf
    (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
        (algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ)
      = reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ :=
  galoisUnitReciprocityMap_eq_reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ

/-- ★**制限写像 `Γ_K → Gal(K(x)/K)` を群準同型として束ねたもの**。

在庫 `algEquivRestrictSelf` は関数としてしか無かったので、`coe_algEquivRestrictSelf`
(`K.closure` へ戻せば `σ` そのもの、`rfl`)から乗法性を出して `MonoidHom` にした。
★これで `Subgroup.comap` が使えるようになる。 -/
noncomputable def algEquivRestrictSelfHom :
    (K.closure ≃ₐ[K.carrier] K.closure) →*
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) where
  toFun := algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
  map_one' := AlgEquiv.ext fun _ => Subtype.ext rfl
  map_mul' _ _ := AlgEquiv.ext fun _ => Subtype.ext rfl

@[simp] theorem algEquivRestrictSelfHom_apply (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    algEquivRestrictSelfHom K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
      = algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ := rfl

end Compat

/-! ## §4 `Γ_K` への持ち上げ —— 消費先 `RamificationImageUnits.lean` の入力 -/

/-- ★★★★**鋭い形を `Γ_K` の言葉で(所属版)**。

`σ ∈ Γ_K` の第 `M+1` 段への制限が `Gal(K_{f,M+1}/K)^I` に属することと、
`Art(σ) ∈ U^I_K` であることは同値(`I ≤ M+1`)。

★★これが消費先 `Found/PGC/RamificationImageUnits.lean` の

    `F.S (absGalPrincipalLevel (n+i)) v = absGalPrincipalLevel n`

の**中身**である(`M+1 = n+i`, `I = n`)。★`F` そのものはまだ木に無いので、
本ファイルは `F` を経由せずに部分群の等式として出している。 -/
theorem mem_upperRamificationGroup_algEquivRestrictSelf_iff
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
    (M I : ℕ) (h2 : I ≤ M + 1) (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf (M + 1) (by omega)
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hψ
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hn
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hmem σ
      ∈ upperRamificationGroup
          ((IntermediateField.adjoin K.carrier
              ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier
              ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt} : Set K.closure)))
          (torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1)
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hn
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hmem) (I : ℝ)
      ↔ σ ∈ absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf I := by
  rw [upperRamificationGroup_torsionGen_eq_comap_map_principalUnits K hq hπmax hπne0 f hf0 hf1 hf
      M (by omega) (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hψ
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hn
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hmem I h2,
    Subgroup.mem_comap, MulEquiv.coe_toMonoidHom,
    galoisUnitReciprocityEquiv_algEquivRestrictSelf,
    reciprocityMap_psiGenSeq_eq_mk_reciprocityUnits,
    mem_map_mk'_iff_of_le _ _ (principalUnits_antitone K π h2)]
  exact Iff.rfl

/-- ★★★★★**鋭い形を `Γ_K` の言葉で(部分群版)** ——

    `res_{M+1}^{-1}(Gal(K_{f,M+1}/K)^I) = Art^{-1}(U^I_K)`   (`I ≤ M+1`)。

★★**消費先 `Found/PGC/RamificationImageUnits.lean` の
`F.S (absGalPrincipalLevel (n+i)) v = absGalPrincipalLevel n` はこれである**
(`M+1 = n+i`, `I = n`)。段データ `F` が立てば、そのまま仮定 2 本が閉じる。 -/
theorem comap_upperRamificationGroup_eq_absGalPrincipalLevel
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
    (M I : ℕ) (h2 : I ≤ M + 1) :
    Subgroup.comap
        (algEquivRestrictSelfHom K hq hπmax hπne0 f hf0 hf1 hf (M + 1) (by omega)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hmem)
        (upperRamificationGroup
          ((IntermediateField.adjoin K.carrier
              ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier
              ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt} : Set K.closure)))
          (torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1)
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hn
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hmem) (I : ℝ))
      = absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf I := by
  ext σ
  exact mem_upperRamificationGroup_algEquivRestrictSelf_iff K hq hπmax hπne0 f hf0 hf1 hf M I h2 σ

/-- ★★★**添字の独立な検算**:
`ker(Γ_K → Gal(K_{f,M+1}/K)) = Art^{-1}(U^{M+1}_K)`。

`I = M+1` を上の定理に代入し、`upperRamificationGroup_torsionGen_eq_bot`(在庫、
＝ Prop 6.14 の命題文)で左辺の上付き分岐群を `⊥` に潰しただけ。

★★左辺は Galois 対応で `Gal(K̄/K_{f,M+1})` であり、右辺は消費先の
`absGalPrincipalLevel (M+1)` である。両者が一致したので、
**「レベル `m` と塔の段 `m` はずれない」という消費先の主張は正しい。**
(消費先 `absGalPrincipalLevel_eq_fixingSubgroup` と独立な経路での確認になっている。) -/
theorem ker_algEquivRestrictSelfHom_eq_absGalPrincipalLevel
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
    (M : ℕ) :
    (algEquivRestrictSelfHom K hq hπmax hπne0 f hf0 hf1 hf (M + 1) (by omega)
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hψ
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hn
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hmem).ker
      = absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (M + 1) := by
  rw [← comap_upperRamificationGroup_eq_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf
      M (M + 1) le_rfl,
    upperRamificationGroup_torsionGen_eq_bot K hq hπmax hπne0 f hf0 hf1 hf M (by omega)
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hψ
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hn
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hmem]
  exact (MonoidHom.comap_bot _).symm

/-- ★★**消費先の形そのもの**: 段を `n + i`(`i` を動かす)、上付きの番号を `n` に取る。

    `∀ i, res_{n+i}^{-1}(Gal(K_{f,n+i}/K)^n) = Art^{-1}(U^n_K)`

★消費先 `map_limit_reciprocityUnits_eq_principalUnits` は
「`F.S (absGalPrincipalLevel (n+i)) v = absGalPrincipalLevel n` が**すべての `i` で**」
を要求する。★**その `i` すべてを本補題が与える。** -/
theorem forall_comap_upperRamificationGroup_eq_absGalPrincipalLevel
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
    (n : ℕ) (i : ℕ) :
    Subgroup.comap
        (algEquivRestrictSelfHom K hq hπmax hπne0 f hf0 hf1 hf (n + i + 1) (by omega)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (n + i)).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (n + i)).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (n + i)).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (n + i)).hmem)
        (upperRamificationGroup
          ((IntermediateField.adjoin K.carrier
              ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (n + i)).pt} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier
              ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (n + i)).pt} : Set K.closure)))
          (torsionGen K hq hπmax hπne0 f hf0 hf1 hf (n + i + 1)
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (n + i)).pt
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (n + i)).hn
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (n + i)).hmem) ((n + 1 : ℕ) : ℝ))
      = absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + 1) :=
  comap_upperRamificationGroup_eq_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf
    (n + i) (n + 1) (by omega)

end ABC3.Found.PGC
