import ABC3.Found.PGC.HasseArf
import ABC3.Found.PGC.FixedRingAction
import Mathlib.GroupTheory.FiniteAbelian.Basic

/-!
# Hasse-Arf の定理・段 2(`j > 1` の帰納)(Yoshida 2008 Theorem 6.11)

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
本ファイルの形式化はすべて `≠`(`hne : G_n ≠ G_{n+1}`)側で書いてある。

## ★本ノードの守備範囲

★★**本ノードは Theorem 6.11 の段 2(`G = G_1` かつ `j > 1`)のみ**である。
**段 1(`j = 1`、巡回)と段 3(`G ≠ G_1`)は `Found/PGC/HasseArf.lean`(Y15)**が
`sorry` 0 で埋めてある(`exists_natCast_herbrandPhiGroup_of_isCyclic` /
`exists_natCast_herbrandPhiGroup_of_tame_quotient`)。★本ファイルは既存の
`Found/PGC/*.lean` を 1 行も書き換えていない(import のみ)。

原文 (Yoshida08 p.16) の段 2:
> For j > 1, if Gn ̸= Gn+1 we can find H with G/H ∼= Z/pmiZ, and GnH/H ̸= Gn+1H/H.
> We have φH(n) ∈Z≥0 by inductive hypothesis, and
> (G/H)φH(n) ̸= (G/H)φH(n+1) = (G/H)φH(n)+1 by Proposition 6.9.
> As G/H is cyclic, we see φG/H(φH(n)) ∈Z≥0, which is φG(n) by Lemma 6.10(ii).

この 1 文は 4 つの主張を畳んでいる。本ファイルの被覆は次のとおり。

| 原文の句 | 本ファイル |
|---|---|
| `we can find H with G/H ≅ ℤ/p^{m_i}ℤ, and G_nH/H ≠ G_{n+1}H/H` | §2 `exists_cyclic_quotient_of_lowerRamificationGroup_ne`(★**純群論の抽象核** §1 への代入) |
| `(G/H)_{φ_H(n)} ≠ (G/H)_{φ_H(n+1)}` by Prop 6.9 | §5 の証明前半(`herbrand_coe_mul_coe_eq` を 2 回) |
| `(G/H)_{φ_H(n+1)} = (G/H)_{φ_H(n)+1}` | §3 `ramificationGroupReal_eq_succ_of_mem_Ioc` + `herbrandPhi_succ_le` |
| `As G/H is cyclic, φ_{G/H}(φ_H(n)) ∈ ℤ≥0` | Y15 の段 1 を商 `G ⧸ H` に代入(§4 `quotientMulSemiringActionOfTrivial` + §6 `faithfulSMul_quotient_fixedRing`) |
| `which is φ_G(n) by Lemma 6.10(ii)` | Y12 `herbrandPhiGroup_comp` + §4 `herbrandPhiGroup_quotient_eq` |

★**`by inductive hypothesis`(`φ_H(n) ∈ ℤ≥0` の再帰呼び出し)だけは仮定 `hind` で
受けている**(逸脱 3)。理由は数学ではなく配管で、下の「新しく必要になったノード」にある。

## ★★抽象核(§1・§3・§4)と具体層(§2・§5)

★原文が「we can find H」と 1 語で畳んだ部分は、**分岐も付値も Galois も出てこない
純群論**である。§1 に切り出した。

| §1 抽象核(純群論) | 内容 |
|---|---|
| `exists_cyclic_quotient_of_hom` | 巡回群への準同型 1 本から `H = ker ψ` と生成元 `σ` を取り出す |
| `exists_cyclic_quotient_not_mem` | ★★**有限アーベル群 `G`、`N ≤ G`、`x ∉ N` なら `N ≤ H`・`x ∉ H`・`G/H` 巡回 なる `H` がある**(有限アーベル群の構造定理) |
| `coe_mul_ne_coe_mul_of_not_mem` | `N ≤ H`・`x ∈ M`・`x ∉ H` なら `MH ≠ NH`(原文の `G_nH/H ≠ G_{n+1}H/H`) |
| `coe_mul_coe_eq_coe` | `N ≤ H` なら `NH = H` |

| §3 抽象核(順序・実数) | 内容 |
|---|---|
| `ramificationGroupReal_eq_succ_of_mem_Ioc` | ★`k < x ≤ k+1`(`k : ℕ`)なら `G_x = G_{k+1}`。★**`hadj` を使わない**(`RealLeENat` の整数性だけ) |
| `truncENat_le_add_one` / `herbrandSum_succ_le` / `herbrandPhi_succ_le` | ★`φ_H(n+1) ≤ φ_H(n) + 1`。原文が `(G/H)_{φ_H(n+1)} = (G/H)_{φ_H(n)+1}` と書くとき暗黙に使っている |
| `truncENat_eq_of_realLeENat` / `herbrandPhi_eq_self_of_forall` | ★`i(τ) ≥ n+1` が `H` 上で全部成り立てば `φ_H(n) = n`。`G = G_1` から `φ_H(1) = 1` を出すのに使う |

| §4 抽象核(群作用の商) | 内容 |
|---|---|
| `quotientMulSemiringActionOfTrivial` | ★`H` が `C` に自明に作用するなら `G ⧸ H` が `C` に環作用する(`MulSemiringAction.compHom` + `QuotientGroup.lift`) |
| `quotientSMul_mk_fixedRing` | ★それを Y16 の固定環に代入した形(`hq` は `rfl`) |
| `ramIndex_quotient_mk` / `mem_ramificationGroupReal_quotient_mk_iff` | `i_ϖ(σ̄) = i_ϖ(σ)` とその帰結 |
| `herbrandPhiGroup_quotient_eq` | ★★`φ` を `G ⧸ H` 上で計算しても `G` 上で計算しても同じ(Y12 の `phiOf_quotient` に代入) |

## ★★Y15 が「配管が越えられない」とした 7 本の現況(実機で確かめた)

Y15 は段 2 を落とし、理由を「Prop 6.9 / Lemma 6.10(ii) が固定環 `C` と 7 本の適合条件
(`hcomp` `hHtriv` `hinj` `hfixC` `hres` `hAC` `hfix`)を仮定で受ける形になっている」と
書いた。★**そのうち 6 本は Y14 + Y16 で本当に埋まった**:

| 債務 | 供給 |
|---|---|
| `[MulSemiringAction G C]` | Y16 `fixedRingMulSemiringAction`(★`H ⊴ G` が要る) |
| `hcomp` | Y16 `algebraMap_smul_fixedRing`(`rfl`) |
| `hHtriv` | Y16 `smul_fixedRing_eq_self` |
| `hinj` | Y14 `fixedRing_injective` |
| `hfixC` | Y14 `exists_algebraMap_fixedRing` |
| `hres` | Y14 `exists_sub_mem_fixedRing` |
| `hAC` | Y14 `exists_algebraMap_fixedRing_eq` |
| `[IsDiscreteValuationRing C]` | Y14 `fixedRing_isDiscreteValuationRing`(★`[Fintype ↥H]`) |

★**残る 1 本は `hfix`(`𝒪_{K′}^H ⊆ 𝒪[π″]`、Y17 の持ち場)である。** 本ファイルは
これを仮定として受け取る(Y17 が着地したら差し込める形にしてある)。

★★**しかしそれだけでは足りなかった**(★これが本ノードの実測である)。段 2 は
Y15 の段 1 を**商 `G ⧸ H` が固定環 `C` に作用する形で**使うので、Y10/Y12 の 7 本とは
別に次が要る。★本ファイルはこれらを**仮定・インスタンス引数として明示**している。

1. ★`[MulSemiringAction (G ⧸ H) C]` —— ★§4 `quotientMulSemiringActionOfTrivial` が
   **供給する**(Y16 の `smul_fixedRing_eq_self` を渡すだけ)。新しい債務ではない。
2. ★★`[Algebra A C]` —— **まだ無い**。`C = B^H` は `Subring B` なので
   `Algebra C B` は付くが、`Algebra A C`(底 `𝒪 ⊆ 𝒪_{K″}`)は木にインスタンスが無い。
   Y14 の `exists_algebraMap_fixedRing_eq` は「像が固定環に入る」ことしか言っていない。
3. ★`[FaithfulSMul (G ⧸ H) C]` —— ★**§6 `faithfulSMul_quotient_fixedRing` が
   `C = B^H` に対して証明した**(Y16 の Lemma 6.8 + `ENat` の有限和が `⊤` なら
   項の 1 つが `⊤`)。★`hfix` の下でのみ言えるので、`hfix` が要る点は変わらない。
4. `hadjC`(`𝒪_{K″} = 𝒪[ϖ]`)・`hresC`(`K″/K` が完全分岐)・`[CharP (ResidueField C) p]`・
   `[Algebra A C]` —— **まだ無い**。★`hadjC` は Y17 の `hfix` の `K″` 側であり、
   同じ債務である。

★★**したがって「配管が越えられない」という Y15 の判定は、いまは「あと 2 本」
(`[Algebra A C]` と `[CharP (ResidueField C) p]`、および Y17 の `hfix`/`hadjC`)まで
縮んだ**。★新ノードとして報告してある(下の「新しく必要になったノード」)。

## 逸脱の記録

1. **`G` の可換性は `∀ x y : G, x * y = y * x` という命題で渡す**(`CommGroup` 構造を
   要求しない)。先行ノード `AbelianJumpDivisibility.lean`(Y13)・`HasseArf.lean`(Y15)の
   逸脱と同じ理由: 具体層の `G = Gal(K′/K)` は `Group` インスタンスしか持たない。
   ★§1 の証明の中でだけ `letI : CommGroup (G ⧸ N)` を立てている(`G` 自身の
   インスタンスは触らない)。
2. ★★**原文の `G/H ≅ ℤ/p^{m_i}ℤ` は「`G/H` が巡回」とだけ読み替えた。**
   原文は有限アーベル `p` 群の直和分解の第 `i` 成分を名指ししているが、段 2 が使うのは
   「`G/H` が巡回」だけである(位数が `p` べきであることは `G` が `p` 群であることから
   自動で出る)。★**構造定理そのものではなく「点を分離する巡回商」を使った。**
   理由は下の「構造定理か巡回商か」。
3. **`j` に関する帰納そのもの(`φ_H(n) ∈ ℤ≥0` の再帰呼び出し)は仮定 `hind` で受ける。**
   原文の「by inductive hypothesis」に当たる。★理由は数学ではなく配管である:
   帰納の各段で塔 `K′/K^H/K` を立て直す必要があり、上記 2・3 が未供給だからである。
   ★**`hind` を除いた残りはすべて証明してある。**
4. **`G/H` の巡回性は `∀ g : G, ∃ k : ℤ, (σ ^ k)⁻¹ * g ∈ H` の形で持ち回る。**
   `Subgroup.zpowers (σ : G ⧸ H) = ⊤` と同値(§5 で変換している)。★商群を作らずに
   §1 を述べられるので、`[H.Normal]` を §1 が要求しなくて済む。
5. **`G = G_1`(原文の「First assume G = G_1」)は `h1 : G_1 = ⊤` で受ける。**
   ★ここから `σ̄ ∈ (G/H)_1`(段 1 の入力)と `orderOf σ̄ = p^m` の両方を導いている
   (仮定として置いていない)。

## ★構造定理か「非自明な巡回商」か

★**巡回商を選んだ。** 理由:

* 段 2 が実際に使うのは「`G_n H/H ≠ G_{n+1}H/H` かつ `G/H` 巡回」だけであり、
  直和分解の成分の位数 `p^{m_i}` は 1 度も使わない(位数が `p` べきであることは
  `IsPGroup.to_quotient` で別に出る)。
* 分解 `G ≅ ⊕ ℤ/p^{m_i}ℤ` を経由すると、`H` を「第 `i` 成分の核」として構成し直す
  手間がかかる。**点を分離する準同型 1 本**なら `H := ker ψ` で終わる。
* ★ただし**その準同型を作るのに構造定理は使っている**
  (`CommGroup.equiv_prod_multiplicative_zmod_of_finite`)。乗法版があるので
  `Additive`/`Multiplicative` の往復が要らない。
  ★import を 1 行足した: `Mathlib.GroupTheory.FiniteAbelian.Basic`
  (★`Unknown constant` は「無い」ではなく「import していない」——`lean-idioms.md` #68)。

## 退化の自己検査

* ★★**`G` が可換であることを落とすと §1 は偽**。非可換有限群(例: `A_5`)は
  非自明な巡回商を持たないことがある。§1 は `habel` を構造定理に渡す 1 箇所で使う。
* ★★**`G_n ≠ G_{n+1}` を落とすと空虚**(★**原典のテキストは `=` にしか見えない**)。
  §2 の `hne` がその仮定であり、`SetLike.exists_of_lt` で「跳びの証人 `x`」を作る
  唯一の入口である。これを落とすと `x` が取れず、`H` の条件 `x ∉ H` が空になる。
* ★**`G = G_1`(`h1`)を落とすと段 2 は成り立たない**(段 3 が `G ≠ G_1` を扱う)。
  §5 では `h1` を 2 度使う: `φ_H(1) = 1`(⇒ `σ̄ ∈ (G/H)_1`)と `IsPGroup p G`
  (⇒ `orderOf σ̄ = p^m`)。
* ★**`H ⊴ G`** を落とすと `G ⧸ H` が群にならない(§4・§5 のインスタンス引数)。
  ★§1 は `H.Normal` を**要求しない**(商を作らない形で述べたから。逸脱 4)。
* ★**`[Fintype ↥H]`** を落とすと `φ_H` が定義できない(`herbrandSum` の有限和)。
  ★Y14 の `fixedRing_isDiscreteValuationRing` も同じ理由でこれを要求する。
* ★**`hind`(`φ_H(n) ∈ ℤ≥0`)を落とすと結論は出ない**。原文も
  「We have φH(n) ∈Z≥0 by inductive hypothesis」と明示している。
* ★`ℕ∞` の切り詰め引き算・除算は 1 箇所も書いていない(`lean-idioms.md` #102)。
  除算は `herbrandPhi` / `phiOf` の定義の中にしかない。

## 新しく必要になったノード

1. ★★`Algebra A ↥(fixedRing B H)`(底の環が固定環に入ることのインスタンス化)と、
   それに伴う `SMulCommClass (G ⧸ H) A ↥(fixedRing B H)` / `hresC`。
2. `CharP (ResidueField ↥(fixedRing B H)) p`(完全分岐なので `K′` と同じ剰余体)。
3. ★Y17 の `hfix` / `hadjC`(`𝒪_{K″} = 𝒪[ϖ]`。本ファイルは仮定として受け取っている)。

★★`FaithfulSMul (G ⧸ H) ↥(fixedRing B H)` は**本ノードが払った**
(§6 `faithfulSMul_quotient_fixedRing`)。
★`MulSemiringAction (G ⧸ H) ↥(fixedRing B H)` と `hq` も**本ノードが払った**
(§4 `quotientMulSemiringActionOfTrivial` / `quotientSMul_mk_fixedRing`、★`rfl`)。

★上の 3 本が入ると、本ファイルの §6 から対応する仮定が消え、`hind`(帰納法の仮定)を
除いて段 2 が底の設定だけで閉じる。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing Pointwise

def exists_cyclic_quotient_of_lowerRamificationGroup_ne.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

def exists_natCast_herbrandPhiGroup_of_cyclic_quotient.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

/-! ## §1 抽象核 —— 「非自明な巡回商が取れる」(純群論)

★分岐・付値・Galois の語彙が 1 つも出てこない。商群すら作らない(逸脱 4)。 -/

/-- ★★**抽象核** —— 巡回群への準同型 `ψ` が `x` を殺さないなら、`H := ker ψ` は
`N` を含み `x` を含まず、`G/H` は 1 元 `σ` で生成される。

「`G/H` が巡回」は `∀ g, ∃ k : ℤ, (σ^k)⁻¹ * g ∈ H` と書いてある(商群を作らない)。 -/
theorem exists_cyclic_quotient_of_hom {G Z : Type*} [Group G] [Group Z] [IsCyclic Z]
    (ψ : G →* Z) {N : Subgroup G} (hN : N ≤ ψ.ker) {x : G} (hx : ψ x ≠ 1) :
    ∃ (H : Subgroup G) (σ : G), N ≤ H ∧ x ∉ H ∧ ∀ g : G, ∃ k : ℤ, (σ ^ k)⁻¹ * g ∈ H := by
  obtain ⟨y, hy⟩ := IsCyclic.exists_generator (α := ψ.range)
  obtain ⟨σ, hσ⟩ := y.2
  refine ⟨ψ.ker, σ, hN, fun hc => hx (MonoidHom.mem_ker.1 hc), fun g => ?_⟩
  obtain ⟨k, hk⟩ := hy ⟨ψ g, ⟨g, rfl⟩⟩
  refine ⟨k, ?_⟩
  have hk' : ((y : Z)) ^ k = ψ g := congrArg Subtype.val hk
  rw [MonoidHom.mem_ker, map_mul, map_inv, map_zpow, hσ, hk', inv_mul_cancel]

/-- ★★★**抽象核(本ノードの心臓)** —— 有限アーベル群 `G` と部分群 `N`、`x ∉ N` に対し、
`N ≤ H`・`x ∉ H`・`G/H` 巡回 なる部分群 `H` が取れる。

原文 (Yoshida08 p.16) の「we can find H with G/H ∼= Z/pmiZ」に当たる。
★**原文は直和分解の第 `i` 成分を名指ししているが、段 2 が使うのは「巡回」だけである**
(逸脱 2)。位数が `p` べきであることは `G` が `p` 群であることから別に出る。

★★**可換性を落とすと偽**(非可換単純群は非自明な巡回商を持たない)。
★可換性は `∀ x y, x * y = y * x` の命題で受ける(逸脱 1)。`CommGroup` 構造を立てるのは
証明の中の `G ⧸ N` に対してだけである。

段取り: `G ⧸ N` に有限アーベル群の構造定理
(`CommGroup.equiv_prod_multiplicative_zmod_of_finite`)を当て、`x̄ ≠ 1` を潰さない
成分 `i` への射影 `ψ` を取り、`exists_cyclic_quotient_of_hom` に渡す。 -/
theorem exists_cyclic_quotient_not_mem {G : Type*} [Group G] [Finite G]
    (habel : ∀ x y : G, x * y = y * x) {N : Subgroup G} {x : G} (hx : x ∉ N) :
    ∃ (H : Subgroup G) (σ : G), N ≤ H ∧ x ∉ H ∧ ∀ g : G, ∃ k : ℤ, (σ ^ k)⁻¹ * g ∈ H := by
  haveI hNn : N.Normal := ⟨fun a ha g => by rwa [habel g a, mul_assoc, mul_inv_cancel, mul_one]⟩
  letI : CommGroup (G ⧸ N) :=
    { (inferInstance : Group (G ⧸ N)) with
      mul_comm := by
        intro a b
        induction a using QuotientGroup.induction_on with
        | _ a => induction b using QuotientGroup.induction_on with
          | _ b => exact congrArg (QuotientGroup.mk : G → G ⧸ N) (habel a b) }
  obtain ⟨ι, _, n, _, ⟨e⟩⟩ := CommGroup.equiv_prod_multiplicative_zmod_of_finite (G ⧸ N)
  have hxne : (QuotientGroup.mk x : G ⧸ N) ≠ 1 := fun hc => hx ((QuotientGroup.eq_one_iff x).1 hc)
  obtain ⟨i, hi⟩ : ∃ i : ι, e (QuotientGroup.mk x) i ≠ 1 := by
    contrapose! hxne
    exact (MulEquiv.map_eq_one_iff e).mp (funext hxne)
  refine exists_cyclic_quotient_of_hom
    (((Pi.evalMonoidHom (fun i => Multiplicative (ZMod (n i))) i).comp e.toMonoidHom).comp
      (QuotientGroup.mk' N)) (fun a ha => ?_) hi
  simp only [MonoidHom.mem_ker, MonoidHom.coe_comp, Function.comp_apply,
    QuotientGroup.mk'_apply, (QuotientGroup.eq_one_iff a).2 ha, map_one]

/-- ★**抽象核** —— `N ≤ H`、`x ∈ M`、`x ∉ H` なら `M·H ≠ N·H`。

原文の `G_nH/H ≠ G_{n+1}H/H` を、商群を作らずに `G` の中の集合の等式で書いた形
(Y11 の `herbrand_coe_mul_coe_eq` がこの形で Proposition 6.9 を述べている)。 -/
theorem coe_mul_ne_coe_mul_of_not_mem {G : Type*} [Group G] {M N H : Subgroup G}
    (hNH : N ≤ H) {x : G} (hxM : x ∈ M) (hxH : x ∉ H) :
    (M : Set G) * (H : Set G) ≠ (N : Set G) * (H : Set G) := by
  intro hc
  have hx : x ∈ (M : Set G) * (H : Set G) := ⟨x, hxM, 1, H.one_mem, mul_one x⟩
  rw [hc] at hx
  obtain ⟨a, ha, b, hb, hab⟩ := hx
  exact hxH (hab ▸ H.mul_mem (hNH ha) hb)

/-- ★**抽象核** —— `N ≤ H` なら `N·H = H`。 -/
theorem coe_mul_coe_eq_coe {G : Type*} [Group G] {N H : Subgroup G} (hNH : N ≤ H) :
    (N : Set G) * (H : Set G) = (H : Set G) := by
  ext g
  constructor
  · rintro ⟨a, ha, b, hb, rfl⟩
    exact H.mul_mem (hNH ha) hb
  · intro hg
    exact ⟨1, N.one_mem, g, hg, one_mul g⟩

/-! ## §2 具体層 —— 跳びを保つ巡回商が取れる -/

/-- ★★★**原文「we can find H with G/H ∼= Z/pmiZ, and GnH/H ̸= Gn+1H/H」**。

`G` 可換・`G_n ≠ G_{n+1}` なら、`G_{n+1} ≤ H`・`G/H` 巡回 で、しかも
跳びの証人 `x ∈ G_n ∖ H` を持つ `H` が取れる。

★`G_n H ≠ G_{n+1} H` は結論の `x` から `coe_mul_ne_coe_mul_of_not_mem` で出る形にした
(§5 はその形では使わず、`x` を直接使う)。

★★**`hne` を落とすと空虚**: `x` が取れなくなる。★★**`habel` を落とすと偽**。 -/
theorem exists_cyclic_quotient_of_lowerRamificationGroup_ne {B : Type*} [CommRing B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [Finite G]
    (habel : ∀ x y : G, x * y = y * x) {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    ∃ (H : Subgroup G) (σ : G) (x : G),
      lowerRamificationGroup B G (n + 1) ≤ H ∧ x ∈ lowerRamificationGroup B G n ∧ x ∉ H ∧
        ∀ g : G, ∃ k : ℤ, (σ ^ k)⁻¹ * g ∈ H := by
  obtain ⟨x, hxn, hxn1⟩ :=
    SetLike.exists_of_lt (lt_of_le_of_ne (lowerRamificationGroup_antitone B G (Nat.le_succ n))
      (Ne.symm hne))
  obtain ⟨H, σ, hNH, hxH, hgen⟩ := exists_cyclic_quotient_not_mem habel hxn1
  exact ⟨H, σ, x, hNH, hxn, hxH, hgen⟩

/-! ## §3 抽象核 —— 実数添字の分岐群と `φ_H` の 1 段の増分 -/

/-- ★★**抽象核** —— `k : ℕ` と `k < x ≤ k+1` なら `G_x = G_{k+1}`。

原文の `(G/H)_{φ_H(n+1)} = (G/H)_{φ_H(n)+1}` はこれである。

★★**`hadj`(`𝒪 = 𝒪[π]`)を使わない。** 使うのは `RealLeENat` と `i(σ)` が
`ℕ∞` 値であること(整数性)だけである。Y11 の `ramificationGroupReal_eq_of_mem_Ioc` は
`ℕ` 添字の `G_m` に落とすので `hadj` を要求するが、本補題は実数添字どうしを比べるので
要らない。 -/
theorem ramificationGroupReal_eq_succ_of_mem_Ioc {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] (α : B) {k : ℕ}
    {x : ℝ} (h1 : (k : ℝ) < x) (h2 : x ≤ (k : ℝ) + 1) :
    ramificationGroupReal (G := G) α x = ramificationGroupReal (G := G) α ((k + 1 : ℕ) : ℝ) := by
  ext σ
  simp only [mem_ramificationGroupReal]
  rcases eq_or_ne (ramIndex α σ) ⊤ with h | h
  · rw [h]; exact iff_of_true RealLeENat.top RealLeENat.top
  · rw [← ENat.coe_toNat h, realLeENat_coe, realLeENat_coe]
    push_cast
    constructor
    · intro hh
      have h3 : (k : ℝ) + 1 < ((ramIndex α σ).toNat : ℝ) := by linarith
      have h4 : k + 1 < (ramIndex α σ).toNat := by exact_mod_cast h3
      have h5 : ((k + 1 : ℕ) : ℝ) + 1 ≤ ((ramIndex α σ).toNat : ℝ) := by exact_mod_cast h4
      push_cast at h5
      linarith
    · intro hh; linarith

/-- ★**抽象核** —— `min{x, r+1} ≤ min{x, r} + 1`。 -/
theorem truncENat_le_add_one (x : ℕ∞) (r : ℝ) : truncENat x (r + 1) ≤ truncENat x r + 1 := by
  unfold truncENat; split
  · linarith
  · rcases le_total ((x.toNat : ℝ)) r with h | h
    · rw [min_eq_left h, min_eq_left (by linarith : (x.toNat : ℝ) ≤ r + 1)]; linarith
    · rw [min_eq_right h]
      exact (min_le_right _ _).trans (by linarith)

/-- `Σ_{τ∈H} min{i(τ), n+2} ≤ Σ_{τ∈H} min{i(τ), n+1} + |H|`。 -/
theorem herbrandSum_succ_le {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] (α : B) (H : Subgroup G) [Fintype H] (n : ℝ) :
    herbrandSum α H (n + 1) ≤ herbrandSum α H n + (Nat.card H : ℝ) := by
  unfold herbrandSum
  calc ∑ τ : H, truncENat (ramIndex α (τ : G)) (n + 1 + 1)
      ≤ ∑ τ : H, (truncENat (ramIndex α (τ : G)) (n + 1) + 1) :=
        Finset.sum_le_sum fun τ _ => truncENat_le_add_one _ _
    _ = (∑ τ : H, truncENat (ramIndex α (τ : G)) (n + 1)) + (Nat.card H : ℝ) := by
        rw [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ, nsmul_eq_mul, mul_one,
          Nat.card_eq_fintype_card]

/-- ★★**`φ_H(n+1) ≤ φ_H(n) + 1`** —— 原文が `(G/H)_{φ_H(n+1)} = (G/H)_{φ_H(n)+1}` と
書くときに(`φ_H` の狭義単調性と合わせて)暗黙に使っている評価。 -/
theorem herbrandPhi_succ_le {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] (α : B) (H : Subgroup G) [Fintype H] (n : ℝ) :
    herbrandPhi α H (n + 1) ≤ herbrandPhi α H n + 1 := by
  have hpos : (0 : ℝ) < (Nat.card H : ℝ) := by exact_mod_cast Nat.card_pos (α := H)
  have h := herbrandSum_succ_le α H n
  rw [herbrandPhi, herbrandPhi]
  have h2 : herbrandSum α H (n + 1) / (Nat.card H : ℝ)
      ≤ herbrandSum α H n / (Nat.card H : ℝ) + 1 := by
    rw [div_add' _ _ _ hpos.ne', div_le_div_iff_of_pos_right hpos, one_mul]
    exact h
  linarith

/-- `i(σ) ≥ r` なら `min{i(σ), r} = r`。 -/
theorem truncENat_eq_of_realLeENat {x : ℕ∞} {r : ℝ} (h : RealLeENat r x) : truncENat x r = r := by
  rw [truncENat]
  split
  · rfl
  · rename_i hx
    rcases h with h | h
    · exact absurd h hx
    · exact min_eq_right h

/-- ★★**抽象核** —— `H` のすべての元が `i(τ) ≥ n+1` を満たすなら `φ_H(n) = n`。

★`G = G_1`(`h1`)から `φ_H(1) = 1` を出すのに使う。 -/
theorem herbrandPhi_eq_self_of_forall {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] (α : B)
    (H : Subgroup G) [Fintype H] (n : ℝ)
    (h : ∀ τ : H, RealLeENat (n + 1) (ramIndex α (τ : G))) :
    herbrandPhi α H n = n := by
  have hpos : (0 : ℝ) < (Nat.card H : ℝ) := by exact_mod_cast Nat.card_pos (α := H)
  have hsum : herbrandSum α H n = (Nat.card H : ℝ) * (n + 1) := by
    unfold herbrandSum
    rw [Finset.sum_congr rfl (fun τ _ => truncENat_eq_of_realLeENat (h τ)), Finset.sum_const,
      Finset.card_univ, nsmul_eq_mul, Nat.card_eq_fintype_card]
  rw [herbrandPhi, hsum]
  field_simp
  ring

/-! ## §4 抽象核 —— `H` が自明に作用するなら商群 `G ⧸ H` が作用する -/

/-- ★★**抽象核** —— `H ⊴ G` が `C` に自明に作用するなら、`G ⧸ H` が `C` に環作用する。

★Y16 `fixedRingMulSemiringAction` は `G` そのものの作用を与える。段 2 は
Y15 の段 1(`[FaithfulSMul G B]` を要求する)を商に代入するので、
**忠実になり得る `G ⧸ H` の作用**が要る。それを供給するのが本 `def` である。

★`MulSemiringAction.compHom` と `QuotientGroup.lift` の合成なので、
`(σ̄) • c = σ • c` は `rfl` で出る(`quotientSMul_mk`)。
★インスタンスにはしない(`G ⧸ H` の作用は文脈ごとに選ぶものなので、
`letI` で入れる)。 -/
@[reducible]
def quotientMulSemiringActionOfTrivial {G : Type*} [Group G] {H : Subgroup G} [H.Normal]
    (C : Type*) [Ring C] [MulSemiringAction G C] (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c) :
    MulSemiringAction (G ⧸ H) C :=
  MulSemiringAction.compHom C
    (QuotientGroup.lift H (MulSemiringAction.toRingAut G C)
      (fun ρ hρ => by ext c; exact hHtriv ρ hρ c))

/-- 商の作用は代表元の作用である。★`quotientMulSemiringActionOfTrivial` の定義から `rfl`。 -/
theorem quotientSMul_mk {G : Type*} [Group G] {H : Subgroup G} [H.Normal]
    (C : Type*) [Ring C] [MulSemiringAction G C] (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (σ : G) (c : C) :
    letI := quotientMulSemiringActionOfTrivial C hHtriv
    (QuotientGroup.mk σ : G ⧸ H) • c = σ • c := rfl

/-- ★§4 の作用を Y16 の固定環 `C = B^H` に代入した形。★`hq` は **`rfl` で出る**。

★これで §5 / §6 の仮定 `[MulSemiringAction (G ⧸ H) C]` と `hq` はどちらも
供給済みになる(残る債務はファイル冒頭の 3 本である)。 -/
theorem quotientSMul_mk_fixedRing {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {H : Subgroup G} [H.Normal] (σ : G) (c : ↥(fixedRing B H)) :
    letI := quotientMulSemiringActionOfTrivial (↥(fixedRing B H))
      (smul_fixedRing_eq_self (B := B) (H := H))
    (QuotientGroup.mk σ : G ⧸ H) • c = σ • c := rfl

/-- `i_ϖ` は代表元に依らない。 -/
theorem ramIndex_quotient_mk {G : Type*} [Group G] {H : Subgroup G} [H.Normal] {C : Type*}
    [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    [MulSemiringAction (G ⧸ H) C]
    (hq : ∀ (σ : G) (c : C), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c) (ϖ : C) (σ : G) :
    ramIndex ϖ (QuotientGroup.mk σ : G ⧸ H) = ramIndex ϖ σ := by
  simp only [ramIndex, hq]

/-- ★★**`φ` を `G ⧸ H` の上で計算しても `G` の上で計算しても同じ**。

Y12 の `phiOf_quotient`(ファイバーの個数が `|H|`)に代入するだけ。
★段 2 は Y12 の Lemma 6.10(ii) を `herbrandPhiGroup G ϖ`(`G` 上の和)の形で受け取り、
Y15 の段 1 は `herbrandPhiGroup (G ⧸ H) ϖ`(商の上の和)を返すので、この橋が要る。 -/
theorem herbrandPhiGroup_quotient_eq {G : Type*} [Group G] [Fintype G] {H : Subgroup G}
    [H.Normal] [Fintype (G ⧸ H)] {C : Type*} [CommRing C] [IsDomain C]
    [IsDiscreteValuationRing C] [MulSemiringAction G C] [MulSemiringAction (G ⧸ H) C]
    (hq : ∀ (σ : G) (c : C), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c) (ϖ : C) (n : ℝ) :
    herbrandPhiGroup (G ⧸ H) ϖ n = herbrandPhiGroup G ϖ n := by
  rw [herbrandPhiGroup, herbrandPhiGroup, ← phiOf_quotient H (fun q : G ⧸ H => ramIndex ϖ q) n]
  simp only [ramIndex_quotient_mk hq]

/-- 実数添字の分岐群は商でも同じ元を拾う。 -/
theorem mem_ramificationGroupReal_quotient_mk_iff {G : Type*} [Group G] {H : Subgroup G}
    [H.Normal] {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C]
    [MulSemiringAction G C] [MulSemiringAction (G ⧸ H) C]
    (hq : ∀ (σ : G) (c : C), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c) (ϖ : C) (x : ℝ) (σ : G) :
    (QuotientGroup.mk σ : G ⧸ H) ∈ ramificationGroupReal (G := G ⧸ H) ϖ x
      ↔ σ ∈ ramificationGroupReal (G := G) ϖ x := by
  simp only [mem_ramificationGroupReal, ramIndex_quotient_mk hq]

/-! ## §5 具体層 —— 段 2 の組み立て -/

/-- ★★★★**Yoshida 2008 Theorem 6.11 (Hasse-Arf) の段 2**(`G = G_1`、`j > 1`)。

原文 (Yoshida08 p.16):
> Theorem 6.11 (Hasse-Arf). If G is abelian, n ∈ Z[bb]_≥0 and G_n = G_n+1, then φ_G(n) ∈ Z[bb]_≥0.

★★逐語の `G_n = G_n+1` は `pdftotext` が `≠` の斜線を落とした形である。ここは
「跳びの証人 `x ∈ G_n ∖ H`」の形で `≠` 側を使っている。

原文の段 2 の 1 文:
> For j > 1, if Gn ̸= Gn+1 we can find H with G/H ∼= Z/pmiZ, and GnH/H ̸= Gn+1H/H.
> We have φH(n) ∈Z≥0 by inductive hypothesis, and
> (G/H)φH(n) ̸= (G/H)φH(n+1) = (G/H)φH(n)+1 by Proposition 6.9.
> As G/H is cyclic, we see φG/H(φH(n)) ∈Z≥0, which is φG(n) by Lemma 6.10(ii).

仮定は原文が置いているものと、木の Prop 6.9 / Lemma 6.10(ii) が要求する塔の適合条件である:

* `h1` : ★**`G = G_1`**(原文「First assume G = G_1」)。★落とすと段 2 は偽。
* `hind` : ★**`φ_H(n) = k ∈ ℤ≥0`**(原文「by inductive hypothesis」。逸脱 3)。
* `hgen` : `G/H` が `σ̄` で生成される(原文「G/H ∼= Z/pmiZ」の使う部分。逸脱 2)。
* `hxn` / `hxH` / `hle` : ★**跳びが商に降りる**(原文「GnH/H ̸= Gn+1H/H」)。
  §2 `exists_cyclic_quotient_of_lowerRamificationGroup_ne` がこの 3 つを供給する。
* `hcomp` `hHtriv` `hinj` `hfixC` `hresB` `hAC` `hadj` `hfix` : Y10/Y11/Y12 の塔の適合条件。
  ★このうち `hfix` 以外は `C = B^H` に対し Y14 + Y16 が供給する(ファイル冒頭の表)。
* `hq` : 商 `G ⧸ H` の作用が代表元の作用に一致する(§4 が供給する)。
* `hadjC` `hresC` `[Algebra A C]` `[FaithfulSMul (G ⧸ H) C]` `[CharP (ResidueField C) p]` :
  ★★**`K″/K` 側の設定。まだ木に無い(ファイル冒頭「新しく必要になったノード」)。**

★★`orderOf σ̄ = p^m` と `σ̄ ∈ (G/H)_1`(Y15 の段 1 の入力)は**仮定していない**。
`h1` から `IsPGroup p G` と `φ_H(1) = 1` を経由して導いている。 -/
theorem exists_natCast_herbrandPhiGroup_of_cyclic_quotient
    {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B]
    {H : Subgroup G} [H.Normal] [Fintype H] [Fintype (G ⧸ H)]
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [Algebra C B] [Algebra A C]
    [MulSemiringAction G C] [MulSemiringAction (G ⧸ H) C] [SMulCommClass (G ⧸ H) A C]
    [FaithfulSMul (G ⧸ H) C]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p] [CharP (ResidueField C) p]
    {π' : B} (hπ' : Irreducible π') {ϖ : C} (hϖ : Irreducible ϖ)
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hinj : Function.Injective (algebraMap C B))
    (hfixC : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hresB : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) →
      y ∈ Algebra.adjoin A ({algebraMap C B ϖ} : Set B))
    (hq : ∀ (σ : G) (c : C), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c)
    (hadjC : Algebra.adjoin A ({ϖ} : Set C) = ⊤)
    (hresC : ∀ c : C, ∃ a : A, c - algebraMap A C a ∈ maximalIdeal C)
    (h1 : lowerRamificationGroup B G 1 = ⊤)
    {σ : G} (hgen : ∀ g : G, ∃ t : ℤ, (σ ^ t)⁻¹ * g ∈ H)
    {n k : ℕ} (hind : herbrandPhi π' H (n : ℝ) = (k : ℝ))
    {x : G} (hxn : x ∈ lowerRamificationGroup B G n) (hxH : x ∉ H)
    (hle : lowerRamificationGroup B G (n + 1) ≤ H) :
    ∃ j : ℕ, herbrandPhiGroup G π' (n : ℝ) = (j : ℝ) := by
  haveI : Fact p.Prime := ⟨hp⟩
  have huniB : maximalIdeal B = Ideal.span {π'} := (irreducible_iff_uniformizer π').mp hπ'
  have huniC : maximalIdeal C = Ideal.span {ϖ} := (irreducible_iff_uniformizer ϖ).mp hϖ
  -- `σ̄` が `G ⧸ H` を生成する(逸脱 4 の形を `Subgroup.zpowers` に直す)。
  have hgenQ : Subgroup.zpowers (QuotientGroup.mk σ : G ⧸ H) = ⊤ := by
    rw [Subgroup.eq_top_iff']
    intro q
    induction q using QuotientGroup.induction_on with
    | _ g =>
      obtain ⟨t, ht⟩ := hgen g
      refine ⟨t, ?_⟩
      show (QuotientGroup.mk σ : G ⧸ H) ^ t = QuotientGroup.mk g
      rw [← QuotientGroup.mk_zpow]
      exact QuotientGroup.eq.2 ht
  -- 原文「First assume G = G1」から `G` が `p` 群であること、したがって
  -- `orderOf σ̄` が `p` べきであることが出る。
  have hpG : IsPGroup p G :=
    (h1 ▸ isPGroup_lowerRamificationGroup_one (A := A) p huniB hπ'.ne_zero hadj).of_equiv
      Subgroup.topEquiv
  obtain ⟨m, hord⟩ := (IsPGroup.iff_orderOf (p := p)).1 (hpG.to_quotient H)
    (QuotientGroup.mk σ : G ⧸ H)
  -- Proposition 6.9 を `n`、`n+1`、`1` の 3 点で使う。
  have h69 : ∀ r : ℝ,
      ((ramificationGroupReal (G := G) π' r : Subgroup G) : Set G) * (H : Set G)
        = ((ramificationGroupReal (G := G) ϖ (herbrandPhi π' H r) : Subgroup G) : Set G) :=
    fun r => herbrand_coe_mul_coe_eq (A := A) hcomp hHtriv hπ' hinj hfixC hresB hAC rfl hadj hfix r
  -- 原文「First assume G = G1」⇒ `φ_H(1) = 1` ⇒ `σ̄ ∈ (G/H)_1`。
  have hphi1 : herbrandPhi π' H (1 : ℝ) = 1 := by
    refine herbrandPhi_eq_self_of_forall π' H 1 (fun τ => ?_)
    have hmem : (τ : G) ∈ ramificationGroupReal (G := G) π' (1 : ℝ) := by
      have hcast := ramificationGroupReal_natCast (A := A) (G := G) huniB hadj 1
      rw [Nat.cast_one] at hcast
      rw [hcast, h1]
      trivial
    exact hmem
  have hσ1 : (QuotientGroup.mk σ : G ⧸ H) ∈ lowerRamificationGroup C (G ⧸ H) 1 := by
    have hmemG : σ ∈ ramificationGroupReal (G := G) ϖ (1 : ℝ) := by
      have hx1 : σ ∈ ((ramificationGroupReal (G := G) π' (1 : ℝ) : Subgroup G) : Set G)
          * (H : Set G) := by
        refine ⟨σ, ?_, 1, H.one_mem, mul_one σ⟩
        have hcast := ramificationGroupReal_natCast (A := A) (G := G) huniB hadj 1
        rw [Nat.cast_one] at hcast
        rw [SetLike.mem_coe, hcast, h1]
        trivial
      rw [h69 1, hphi1] at hx1
      exact hx1
    have hcastC := ramificationGroupReal_natCast (A := A) (G := G ⧸ H) huniC hadjC 1
    rw [Nat.cast_one] at hcastC
    rw [← hcastC]
    exact (mem_ramificationGroupReal_quotient_mk_iff hq ϖ 1 σ).2 hmemG
  -- 原文「GnH/H ̸= Gn+1H/H」⇒「(G/H)_{φ_H(n)} ≠ (G/H)_{φ_H(n)+1}」。
  have hcast : ((n + 1 : ℕ) : ℝ) = (n : ℝ) + 1 := by push_cast; ring
  have hlow : (k : ℝ) < herbrandPhi π' H (((n + 1 : ℕ)) : ℝ) := by
    rw [← hind]
    exact (strictMono_herbrandPhi π' H) (by rw [hcast]; linarith)
  have hup : herbrandPhi π' H (((n + 1 : ℕ)) : ℝ) ≤ (k : ℝ) + 1 := by
    rw [hcast, ← hind]
    exact herbrandPhi_succ_le π' H (n : ℝ)
  have hxk : x ∈ ramificationGroupReal (G := G) ϖ (k : ℝ) := by
    have hx1 : x ∈ ((ramificationGroupReal (G := G) π' ((n : ℕ) : ℝ) : Subgroup G) : Set G)
        * (H : Set G) := by
      refine ⟨x, ?_, 1, H.one_mem, mul_one x⟩
      rw [SetLike.mem_coe, ramificationGroupReal_natCast (A := A) huniB hadj n]
      exact hxn
    rw [h69 ((n : ℕ) : ℝ), hind] at hx1
    exact hx1
  have hxk1 : x ∉ ramificationGroupReal (G := G) ϖ (((k + 1 : ℕ)) : ℝ) := by
    have heq := h69 (((n + 1 : ℕ)) : ℝ)
    rw [ramificationGroupReal_natCast (A := A) huniB hadj (n + 1), coe_mul_coe_eq_coe hle,
      ramificationGroupReal_eq_succ_of_mem_Ioc ϖ hlow hup] at heq
    intro hc
    refine hxH ?_
    have hmem : x ∈ (H : Set G) := by rw [heq]; exact hc
    exact hmem
  have hneQ : lowerRamificationGroup C (G ⧸ H) k ≠ lowerRamificationGroup C (G ⧸ H) (k + 1) := by
    rw [← ramificationGroupReal_natCast (A := A) (G := G ⧸ H) huniC hadjC k,
      ← ramificationGroupReal_natCast (A := A) (G := G ⧸ H) huniC hadjC (k + 1)]
    intro hc
    refine hxk1 ?_
    have := (mem_ramificationGroupReal_quotient_mk_iff hq ϖ ((k : ℕ) : ℝ) x).2 hxk
    rw [hc] at this
    exact (mem_ramificationGroupReal_quotient_mk_iff hq ϖ (((k + 1 : ℕ)) : ℝ) x).1 this
  -- 原文「As G/H is cyclic, we see φG/H(φH(n)) ∈Z≥0」—— Y15 の段 1。
  obtain ⟨j, hj⟩ := exists_natCast_herbrandPhiGroup_of_isCyclic (A := A) (B := C) (G := G ⧸ H)
    p hp hϖ hadjC hresC hσ1 hord hgenQ hneQ
  -- 原文「which is φG(n) by Lemma 6.10(ii)」—— Y12。
  refine ⟨j, ?_⟩
  rw [herbrandPhiGroup_comp (A := A) hcomp hHtriv hπ' hinj hfixC hresB hAC rfl hadj hfix
    ((n : ℕ) : ℝ), hind, ← herbrandPhiGroup_quotient_eq hq ϖ ((k : ℕ) : ℝ), hj]

/-! ## §6 `C = B^H` に代入した形 —— Y14 + Y16 が払える債務を実際に払う -/

/-- ★★★**`G ⧸ H = Gal(K″/K)` は固定環 `B^H = 𝒪_{K″}` に忠実に作用する。**

★ファイル冒頭「新しく必要になったノード」の 2 番目をここで払う(★`hfix` の下で)。

★★**原文にはない一歩である。** Yoshida は `G/H = Gal(K″/K)` と書いた時点で
忠実性を暗黙に使っているが(`Gal` の定義に入っている)、本木の `G` は抽象的な
`MulSemiringAction` なので示す必要がある。段取りは Lemma 6.8 の帰結:

`σ̄` が `C` に自明に作用するなら `i_ϖ(σ) = ⊤` なので、Y16 の Lemma 6.8
`|H| · i_ϖ(σ) = Σ_{τ∈H} i_{π′}(στ)` の左辺は `⊤`。有限和が `⊤` なら項の 1 つが `⊤`、
すなわち `στ = 1` となる `τ ∈ H` があるので `σ ∈ H`、つまり `σ̄ = 1`。

★`[FaithfulSMul G B]` を落とすと `στ = 1` の一歩(`eq_one_of_smul_eq_of_adjoin_eq_top`)が
崩れる。★`[Fintype ↥H]` を落とすと `|H| = 0` で左辺が `0` になり空虚。 -/
theorem faithfulSMul_quotient_fixedRing {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype H]
    [MulSemiringAction (G ⧸ H) ↥(fixedRing B H)]
    (hq : ∀ (σ : G) (c : ↥(fixedRing B H)), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c)
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadjA : Algebra.adjoin A ({π'} : Set B) = ⊤) {ϖ : ↥(fixedRing B H)}
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) →
      y ∈ Algebra.adjoin A ({algebraMap (↥(fixedRing B H)) B ϖ} : Set B)) :
    FaithfulSMul (G ⧸ H) ↥(fixedRing B H) := by
  have key : ∀ q : G ⧸ H, (∀ c : ↥(fixedRing B H), q • c = c) → q = 1 := by
    intro q hqc
    induction q using QuotientGroup.induction_on with
    | _ σ =>
      have hσϖ : σ • ϖ = ϖ := by rw [← hq σ ϖ]; exact hqc ϖ
      have htop : (Nat.card H : ℕ∞) * ramIndex ϖ σ = ⊤ := by
        rw [(ramIndex_eq_top_iff ϖ σ).2 hσϖ]
        exact ENat.mul_top (by exact_mod_cast (Nat.card_pos (α := H)).ne')
      rw [card_mul_ramIndex_eq_sum_ramIndex_fixedRing (A := A) hπ' hresA hadjA ϖ hfix σ] at htop
      obtain ⟨τ, -, hτ⟩ := (WithTop.sum_eq_top (M := ℕ)).1 htop
      have h1 : σ * (τ : G) = 1 :=
        eq_one_of_smul_eq_of_adjoin_eq_top hadjA ((ramIndex_eq_top_iff π' _).1 hτ)
      refine (QuotientGroup.eq_one_iff σ).2 ?_
      have : σ = (τ : G)⁻¹ := by
        rw [← mul_one σ, ← mul_inv_cancel (τ : G), ← mul_assoc, h1, one_mul]
      rw [this]
      exact H.inv_mem τ.2
  refine ⟨fun {q₁ q₂} h => ?_⟩
  have hone : q₂⁻¹ * q₁ = 1 := by
    refine key _ (fun c => ?_)
    rw [mul_smul, h c, ← mul_smul, inv_mul_cancel, one_smul]
  exact (inv_mul_eq_one.1 hone).symm

/-- ★★★**§5 を固定環 `C = 𝒪_{K′}^H` に代入した形**。

★★**Y15 が「越えられない」とした 7 本のうち 6 本がここで消えている**:
`hcomp`(Y16 `algebraMap_smul_fixedRing`)・`hHtriv`(Y16 `smul_fixedRing_eq_self`)・
`hinj`(Y14 `fixedRing_injective`)・`hfixC`(Y14 `exists_algebraMap_fixedRing`)・
`hres`(Y14 `exists_sub_mem_fixedRing`)・`hAC`(Y14 `exists_algebraMap_fixedRing_eq`)、
および `[MulSemiringAction G C]`(Y16 のインスタンス)と
`[IsDiscreteValuationRing C]`(Y14 のインスタンス)。

★**残る 1 本は `hfix`(`𝒪_{K′}^H ⊆ 𝒪[π″]`、Y17 の持ち場)である。**

★★**そして段 2 に固有の 4 本が残る**(ファイル冒頭「新しく必要になったノード」):
`[Algebra A ↥(fixedRing B H)]`・`[FaithfulSMul (G ⧸ H) ↥(fixedRing B H)]`・
`[CharP (ResidueField ↥(fixedRing B H)) p]`・`hadjC`(= `hfix` の `K″` 側)。
★`[MulSemiringAction (G ⧸ H) ↥(fixedRing B H)]` と `hq` は §4
(`quotientMulSemiringActionOfTrivial` に Y16 の `smul_fixedRing_eq_self` を渡す)が供給する。 -/
theorem exists_natCast_herbrandPhiGroup_of_cyclic_quotient_fixedRing
    {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B]
    {H : Subgroup G} [H.Normal] [Fintype H] [Fintype (G ⧸ H)]
    [Algebra A ↥(fixedRing B H)] [MulSemiringAction (G ⧸ H) ↥(fixedRing B H)]
    [SMulCommClass (G ⧸ H) A ↥(fixedRing B H)]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p]
    [CharP (ResidueField ↥(fixedRing B H)) p]
    {π' : B} (hπ' : Irreducible π') {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) →
      y ∈ Algebra.adjoin A ({algebraMap (↥(fixedRing B H)) B ϖ} : Set B))
    (hq : ∀ (σ : G) (c : ↥(fixedRing B H)), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c)
    (hadjC : Algebra.adjoin A ({ϖ} : Set ↥(fixedRing B H)) = ⊤)
    (hresC : ∀ c : ↥(fixedRing B H),
      ∃ a : A, c - algebraMap A (↥(fixedRing B H)) a ∈ maximalIdeal (↥(fixedRing B H)))
    (h1 : lowerRamificationGroup B G 1 = ⊤)
    {σ : G} (hgen : ∀ g : G, ∃ t : ℤ, (σ ^ t)⁻¹ * g ∈ H)
    {n k : ℕ} (hind : herbrandPhi π' H (n : ℝ) = (k : ℝ))
    {x : G} (hxn : x ∈ lowerRamificationGroup B G n) (hxH : x ∉ H)
    (hle : lowerRamificationGroup B G (n + 1) ≤ H) :
    ∃ j : ℕ, herbrandPhiGroup G π' (n : ℝ) = (j : ℝ) := by
  haveI := faithfulSMul_quotient_fixedRing (A := A) hq hπ' hresA hadj hfix
  exact exists_natCast_herbrandPhiGroup_of_cyclic_quotient (A := A) (C := ↥(fixedRing B H))
    p hp hπ' hϖ algebraMap_smul_fixedRing smul_fixedRing_eq_self fixedRing_injective
    exists_algebraMap_fixedRing (exists_sub_mem_fixedRing hresA)
    (fun a => exists_algebraMap_fixedRing_eq a) hadj hfix hq hadjC hresC h1 hgen hind hxn hxH hle

end ABC3.Found.PGC
