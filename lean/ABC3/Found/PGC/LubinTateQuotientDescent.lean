import ABC3.Found.PGC.AbelianSubfieldInLubinTate
import ABC3.Found.PGC.LubinTateDegree

/-!
# 上付き分岐群の商への降下 —— B5 の `htriv` / `hdeg` を外す

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Corollary 6.13 (i)** と
**Theorem 6.15**(どちらも物理 p.17)。構造化済み原文は
`ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/section-6.html` の
`id="cor-6-13"`(`data-pdf-page="17"`, `data-item="Corollary 6.13"`)と
`id="thm-6-15"`(`data-pdf-page="17"`, `data-item="Theorem 6.15 (Local Kronecker-Weber)"`)。

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

## 本ファイルが担当するもの

直前のノード `Found/PGC/AbelianSubfieldInLubinTate.lean`(道 B の B5)は
Theorem 6.15 の最後の 4 文を組み立てたが、そこに 2 つの仮定を残した:

| B5 の仮定 | 原典での意味 | 本ファイル |
|---|---|---|
| `htriv` : `upperRamificationGroup G ϖ m = H` | `Gal(K^m_x/L)^m = {id}` | §1 の降下で「商の上付き分岐群が `⊥`」から出す |
| `hdeg` : `Nat.card (G ⧸ H) = (q−1)q^k` | `[K^m_x : L] = (q−1)q^{m−1}` | §3・§4 で**体の次数**に直し、§3 で Lubin-Tate から供給する |

★★**`hdeg` は本ファイルで完全に外れた。** `le_of_two_quotients_upperRamification_lubinTate`
は `hdeg` を 1 つも持たない(`finrank_adjoin_iteratedLubinTatePsi_succ` が供給する)。
★**`htriv` は「商の上付き分岐群が `⊥`」(`hbot`)に置き換わった。**
`hbot` を B4(Proposition 6.14)の結論から出すには、まだ環同型
`fixedRing 𝒪_{K(x)} H ≅ 𝒪_{K^m_x}` が要る(末尾「残っている穴」)。
★**その環同型さえ来れば §2 の移送核がそのまま使える形にしてある。**

## §1 抽象核 —— 上付き分岐群の商への降下(Corollary 6.13 (i))

`H ⊴ G` が `C` に自明に作用する(つまり `ϖ` が固定環の元である)とき

    `G^m = (QuotientGroup.mk' H)⁻¹((G/H)^m)`          (`upperRamificationGroup_comap_quotient`)
    `(mk' H)(G^m) = (G/H)^m`                          (`map_upperRamificationGroup_quotient`)

★**分岐・付値・Galois の語彙は出るが、Lubin-Tate も局所体も出てこない。**
`C` は「離散付値環 + `G` 作用」だけである。

★★**原典の証明とは道が違う。** 原典は Proposition 6.9 と Lemma 6.10(ii)
(Herbrand の合成公式)で `G^mH/H = (G/H)^{φ_H(φ^{−1}_G(m))} = (G/H)^m` と計算するが、
木の `upperRamificationGroup G ϖ`(`ϖ` は**固定環の素元**)は定義からして
`φ_{G/H}` で番号付けされている(在庫 `herbrandPhiGroup_quotient_eq`)ので、
**Herbrand の合成公式を一度も使わずに** `φ` の一致 → `ψ` の一致 → 群の一致 と
3 段で降りる。★これが「原典より短い道」である。

## §2 抽象核 —— 同型による移送

群同型 `e : G ≃* G′` と環同型 `ψ : C ≃+* C′` が作用を保つ(`hcompat`)なら

    `addVal C′ (ψ c) = addVal C c`      (`addVal_ringEquiv`、★純 mathlib)
    `i_{ψϖ}(e σ) = i_ϖ(σ)`             (`ramIndex_congr_equiv`)
    `φ_{G′} = φ_G`, `ψ_{G′} = ψ_G`, `G^m = e⁻¹(G′^m)`

★これは**「残っている穴」を消費するための道具**である。`fixedRing 𝒪_{K(x)} H ≅ 𝒪_{K^m_x}`
と `G ⧸ H ≃* Gal(K^m_x/K)` が来れば、B4 の結論
`upperRamificationGroup Gal(K^m_x/K) α₀ m = ⊥` が
`upperRamificationGroup (G ⧸ H) ϖ m = ⊥`(= 本ファイルの `hbot`)に移る。
★★**環同型が作用を保つこと(`hcompat`)を落とすと移送は偽になる**(作用を保たない
環同型では `i` が変わる)ので、`hcompat` は仮定として明示的に受け取っている。

## §3 橋渡し(具体層)

| 補題 | 内容 |
|---|---|
| `natCard_residueField_adjoinIntegers_eq` | 完全分岐 `ht` から `q_{K(x)} = q_K`(在庫 `isTotallyRamifiedAdjoin_iff_residueDegree`) |
| `natCard_quotient_fixingSubgroup_restrict` | ★`Nat.card (G ⧸ Gal(K(x)/E)) = [E : K]`(mathlib `finrank_eq_fixingSubgroup_index` + `restrict_algEquiv`) |
| `finrank_adjoin_iteratedLubinTatePsi_succ` | ★`[K^{k+1}_x : K] = (q−1)q^k`(在庫 `finrank_adjoin_iteratedLubinTatePsi` の `n = k+1` 形) |

★★**`restrictNormalHom` を 1 度も使っていない。** 素朴には
`G ⧸ Gal(K(x)/E) ≃* Gal(E/K)`(`restrictNormalHom` の全射性 + 核)を作って
`IsGalois.card_aut_eq_finrank` を当てるところだが、mathlib の
`IntermediateField.finrank_eq_fixingSubgroup_index` は
**`[Normal K E]` を要求せずに** `[IsGalois K L]` だけで `[E : K] = index` を与える。
★`lean-idioms.md` #154(`restrictNormalHom` に `exact` を当てると heartbeats で落ちる)も
これで丸ごと回避できた。

## §4 具体層 —— B5 の 3 本を `hbot` / 体の次数の形に書き直す

`upperRamificationGroupAdjoin_le_of_quotient_upper_eq_bot` /
`le_of_two_quotients_upperRamification_of_finrank` /
`le_of_two_quotients_upperRamification_lubinTate`。
★**B5 の宣言は 1 つも書き換えていない**(本ファイルは新しい名前で言い直すだけである)。

## ★在庫調査の記録(`lean-idioms.md` #158 —— 同じ名前で書いて `already been declared` を出す)

★**この手で 6 本が「既に木にある」と分かった。**
最初の版では `quotientMulSemiringAction` / `ramIndex_quotient_mk` /
`herbrandPhiGroup_quotient` を自分で書いていたが、`Found/PGC/HasseArfInduction.lean` の §4 に

* `quotientMulSemiringActionOfTrivial`(`MulSemiringAction.compHom` + `QuotientGroup.lift`)
* `quotientSMul_mk` / `quotientSMul_mk_fixedRing`
* `ramIndex_quotient_mk` / `mem_ramificationGroupReal_quotient_mk_iff`
* `herbrandPhiGroup_quotient_eq`

が**そのままの形で在った**(`ramIndex_quotient_mk` は仮定の形まで一致していた)。
★本ファイルはそれらを消費して、**`ψ` と `G^m` の 2 段だけ**を足している。

★★**`MulSemiringAction (G ⧸ H) ↥(fixedRing B H)` を `instance` にはしていない。**
`HasseArfInduction` が「`G ⧸ H` の作用は文脈ごとに選ぶものなので `letI` で入れる」と
明記しているので、木の流儀に合わせて `[MulSemiringAction (G ⧸ H) C]` + `hq` を
**仮定として受け取る**形にした(在庫 `ramIndex_quotient_mk` と同じ形)。
消費側は `letI := quotientMulSemiringActionOfTrivial _ smul_fixedRing_eq_self` と
`hq := quotientSMul_mk_fixedRing` を渡せばよい。

## 逸脱の記録

1. **`L = K`(`n = 1`)に固定**。決定 D29 と B4・B5 に合わせた。原典の `L = K_n` 一般は
   扱っていない。
2. **`m` は `k + 1` の形に固定**(`ℕ` の切り詰め引き算・`ℕ∞` の除算を書かないため。
   `lean-idioms.md` #102、B4・B5 と同じ逃げ方)。原典の `m ≥ 1` はこれで尽くされる。
3. **§1 は原典 Corollary 6.13 (i) そのものではない。** 原典の `G^m` は大きい体
   `K′` の素元で番号付けした群だが、本ファイルの `upperRamificationGroup G ϖ m` は
   **固定環 `C` の素元 `ϖ`** で番号付けした群である。両者が一致するのが
   Proposition 6.9 + Lemma 6.10(ii)(在庫 `herbrandPhiGroup_comp`)であり、
   B5 の `upperRamificationGroupAdjoin_le_of_quotient_eq_bot` がそこを担っている。
   ★本ファイルが示すのは「`ϖ` で番号付けした群が本当に `(G/H)^m` の引き戻しである」
   という部分である。★**Y18 が置いた「引き戻しの形」の正当化**にあたる。
4. **`hbot` は仮定である。** 原典の `Gal(K^m_x/L)^m = {id}`(Proposition 6.14 = B4)を
   商群 `G ⧸ H` と固定環 `C` の言葉で述べたもの。B4 の結論とは群も環も違うので、
   末尾「残っている穴」の環同型が要る。

## 退化の自己検査

* ★**`H ⊴ G` はどこから来るか**: §1 では仮定 `[H.Normal]` である。§4 では
  `H = Gal(K(x)/E₁)` が正規であること、すなわち `E₁/K` が正規であることに対応し、
  B5 と同じく `[((IntermediateField.restrict h₁).fixingSubgroup).Normal]` として
  受け取っている。★正規性を落とすと `G ⧸ H` が群にならず、主張が**書けない**。
* ★**上付きと下付きを取り違えていない**: 本ファイルが商へ降ろすのは
  `upperRamificationGroup`(= `G^m`)である。下付き `G_n` は商で素直に降りない
  (それが Herbrand の要点)ので、本ファイルには下付きの降下は 1 本も無い。
  ★§1 の証明が使うのは `φ` の一致(`herbrandPhiGroup_quotient_eq`)であって、
  下付きの降下ではない。
* ★**`H` が `C` に自明に作用することを落とすと §1 は偽**。`hq`(および在庫
  `ramIndex_quotient_mk` の同名の仮定)がそれである。`ϖ` が固定環の元でなければ
  `i_ϖ(σ)` は剰余類上で定数にならない。
* ★**§2 の `hcompat` を落とすと移送は偽**(上に書いた)。
* ★**`le_upperRamificationGroup_of_smul_eq_self` は `H ≤ G^m` を言う。**
  これが無いと Corollary 6.13 (i) の左辺 `G^m H/H` が `G^m/H` に潰れない。
  ★`i_ϖ(τ) = ⊤`(`τ ∈ H`)から出る。`m` に依らないことに注意。
* ★**`ℕ∞` の切り詰め引き算・除算は 1 つも書いていない。** §3 の
  `(q−1)q^k = q^{k+1} − q^k` は `ℕ` の等式であって `ℕ∞` ではない。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing
open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 抽象核 —— 上付き分岐群の商への降下

★局所体も Lubin-Tate も出てこない。`C` は離散付値環、`G` はその上に作用する有限群、
`H ⊴ G` は `C` に自明に作用する(在庫 `HasseArfInduction` の §4 と同じ設定である)。 -/

section Descent

variable {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C]
variable {G : Type*} [Group G] [MulSemiringAction G C]
variable {H : Subgroup G} [H.Normal] [MulSemiringAction (G ⧸ H) C]

/-- `ψ_G = φ_G^{−1}` も商で変わらない。

在庫 `herbrandPhiGroup_quotient_eq`(`φ` の一致)を**関数の等式**に持ち上げてから
`Function.invFun` に入れるだけ。★`φ` が全単射であることは使わない
(`herbrandPsiGroup` は `invFun` で定義されているので、`φ` が等しければ `ψ` も等しい)。 -/
theorem herbrandPsiGroup_quotient_eq [Fintype G] [Fintype (G ⧸ H)]
    (hq : ∀ (σ : G) (c : C), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c)
    (ϖ : C) (m : ℝ) : herbrandPsiGroup (G ⧸ H) ϖ m = herbrandPsiGroup G ϖ m := by
  have hfun : herbrandPhiGroup (G ⧸ H) ϖ = herbrandPhiGroup G ϖ :=
    funext fun n => herbrandPhiGroup_quotient_eq hq ϖ n
  unfold herbrandPsiGroup
  rw [hfun]

/-- ★★★★**上付き分岐群の商への降下**(Corollary 6.13 (i) の引き戻し形)——

    `G^m = (G → G/H)⁻¹((G/H)^m)`

`ϖ` が固定環の素元(= `H` が `C` に自明に作用する)ときの `G^m` は、
定義からして商 `G/H` の上付き分岐群の引き戻しである。

★段取りは 2 行: `ψ` が一致すること(上の `herbrandPsiGroup_quotient_eq`)と
`i_ϖ` が代表元に依らないこと(在庫 `ramIndex_quotient_mk`)。
★★**Herbrand の合成公式(Prop 6.9 / Lemma 6.10(ii))を使っていない。**
原典 Corollary 6.13 (i) の証明はそれを使うが、`ϖ` で番号付けした群については
`φ` の一致が直接出るので不要である(冒頭「逸脱の記録 3」)。 -/
theorem upperRamificationGroup_comap_quotient [Fintype G] [Fintype (G ⧸ H)]
    (hq : ∀ (σ : G) (c : C), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c)
    (ϖ : C) (m : ℝ) :
    upperRamificationGroup G ϖ m
      = Subgroup.comap (QuotientGroup.mk' H) (upperRamificationGroup (G ⧸ H) ϖ m) := by
  ext σ
  rw [Subgroup.mem_comap, mem_upperRamificationGroup, mem_upperRamificationGroup,
    QuotientGroup.mk'_apply, ramIndex_quotient_mk hq, herbrandPsiGroup_quotient_eq hq]

omit [H.Normal] [MulSemiringAction (G ⧸ H) C] in
/-- ★`H` が `C` に自明に作用するなら `H ≤ G^m`(`m` に依らない)。

`τ ∈ H` は `ϖ` を動かさないので `i_ϖ(τ) = ⊤` であり、`RealLeENat` は常に真。
★これが Corollary 6.13 (i) の左辺 `G^m H/H` が `G^m/H` に潰れる理由である。 -/
theorem le_upperRamificationGroup_of_smul_eq_self [Fintype G]
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c) (ϖ : C) (m : ℝ) :
    H ≤ upperRamificationGroup G ϖ m := by
  intro σ hσ
  have htop : ramIndex ϖ σ = ⊤ := (ramIndex_eq_top_iff ϖ σ).2 (hHtriv σ hσ ϖ)
  rw [mem_upperRamificationGroup, htop]
  exact RealLeENat.top

def map_upperRamificationGroup_quotient.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Corollary 6.13", sectionId := "cor-6-13" }

/-- ★★★★★**Yoshida 2008 Corollary 6.13 (i) の字面の形** —— `G^m H/H = (G/H)^m`。

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

原文の `G^mH/H` は、`H ≤ G^m`(上の `le_upperRamificationGroup_of_smul_eq_self`)なので
`G^m` の商への像 `Subgroup.map (QuotientGroup.mk' H) (G^m)` に他ならない。
上の引き戻し形に `map ∘ comap = id`(全射)を当てるだけ。

★逸脱: ここでの `G^m` は**固定環の素元 `ϖ` で番号付けした群**である
(冒頭「逸脱の記録 3」)。原典の `G^m`(大きい体の素元で番号付け)と一致することは
Proposition 6.9 + Lemma 6.10(ii) であり、本ファイルの範囲外である。 -/
theorem map_upperRamificationGroup_quotient [Fintype G] [Fintype (G ⧸ H)]
    (hq : ∀ (σ : G) (c : C), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c)
    (ϖ : C) (m : ℝ) :
    Subgroup.map (QuotientGroup.mk' H) (upperRamificationGroup G ϖ m)
      = upperRamificationGroup (G ⧸ H) ϖ m := by
  rw [upperRamificationGroup_comap_quotient hq,
    Subgroup.map_comap_eq_self_of_surjective (QuotientGroup.mk'_surjective H)]

/-- ★★★**B5 の `htriv` の供給元** —— 商の上付き分岐群が `⊥` なら `G^m = H`。

引き戻し形に `h` を代入し、`comap f ⊥ = ker f`(`MonoidHom.comap_bot`)と
`ker (mk' H) = H`(`QuotientGroup.ker_mk'`)を当てるだけ。

★★これが原典の「`Gal(K^m_x/L)^m = {id}` ならば `Gal(K′K^m_x/L)^m ⊆ Gal(K′K^m_x/K^m_x)`」
の第 1 段である(Corollary 6.13 (ii) の証明の最初の 1 文)。 -/
theorem upperRamificationGroup_eq_of_quotient_eq_bot [Fintype G] [Fintype (G ⧸ H)]
    (hq : ∀ (σ : G) (c : C), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c)
    (ϖ : C) (m : ℝ) (h : upperRamificationGroup (G ⧸ H) ϖ m = ⊥) :
    upperRamificationGroup G ϖ m = H := by
  rw [upperRamificationGroup_comap_quotient hq, h, MonoidHom.comap_bot, QuotientGroup.ker_mk']

end Descent

/-! ## §2 抽象核 —— 同型による移送

★「残っている穴」(`fixedRing 𝒪_{K(x)} H ≅ 𝒪_{K^m_x}`)が来たときに、B4 の結論を
本ファイルの `hbot` に移すための道具。★**分岐の語彙すら出ない補題が 1 本ある**
(`addVal_ringEquiv` は純 mathlib)。 -/

section Transport

/-- ★★**離散付値環の間の環同型は `addVal` を保つ**。★純 mathlib(0.14 秒・一発)。

`c = u·ϖ^n` と書けば `ψ c = (ψ u)·(ψ ϖ)^n` で、`ψ u` は単元・`ψ ϖ` は既約
(`MulEquiv.irreducible_iff`)だから、両辺とも `addVal_def'` で `n` になる。
★`c = 0` の場合だけ別扱い(両辺 `⊤`)。 -/
theorem addVal_ringEquiv {C C' : Type*} [CommRing C] [IsDomain C]
    [IsDiscreteValuationRing C] [CommRing C'] [IsDomain C'] [IsDiscreteValuationRing C']
    (ψ : C ≃+* C') (c : C) : addVal C' (ψ c) = addVal C c := by
  rcases eq_or_ne c 0 with rfl | hc
  · simp
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible C
  obtain ⟨n, u, rfl⟩ := IsDiscreteValuationRing.eq_unit_mul_pow_irreducible hc hϖ
  have hψϖ : Irreducible (ψ ϖ) := (MulEquiv.irreducible_iff (f := (ψ : C ≃* C'))).2 hϖ
  have hmap : ψ ((u : C) * ϖ ^ n)
      = ((Units.map (ψ : C →* C') u : C'ˣ) : C') * (ψ ϖ) ^ n := by simp
  rw [addVal_def' u hϖ n, hmap, addVal_def' _ hψϖ n]

variable {C C' : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C]
variable [CommRing C'] [IsDomain C'] [IsDiscreteValuationRing C']
variable {G G' : Type*} [Group G] [MulSemiringAction G C] [Group G'] [MulSemiringAction G' C']

/-- `i` は(作用を保つ)同型で移る。★`hcompat` を落とすと偽。 -/
theorem ramIndex_congr_equiv (e : G ≃* G') (ψ : C ≃+* C')
    (hcompat : ∀ (σ : G) (c : C), ψ (σ • c) = e σ • ψ c) (ϖ : C) (σ : G) :
    ramIndex (ψ ϖ) (e σ) = ramIndex ϖ σ := by
  rw [ramIndex, ramIndex, ← hcompat, ← map_sub, addVal_ringEquiv]

/-- `φ` は(作用を保つ)同型で移る。★`|G| = |G′|` と和の並べ替えだけ。 -/
theorem herbrandPhiGroup_congr_equiv [Fintype G] [Fintype G'] (e : G ≃* G') (ψ : C ≃+* C')
    (hcompat : ∀ (σ : G) (c : C), ψ (σ • c) = e σ • ψ c) (ϖ : C) (n : ℝ) :
    herbrandPhiGroup G' (ψ ϖ) n = herbrandPhiGroup G ϖ n := by
  have hsum : ∑ τ : G', truncENat (ramIndex (ψ ϖ) τ) (n + 1)
      = ∑ σ : G, truncENat (ramIndex ϖ σ) (n + 1) :=
    (Fintype.sum_equiv e.toEquiv _ _
      (fun σ => (congrArg (fun t => truncENat t (n + 1))
        (ramIndex_congr_equiv e ψ hcompat ϖ σ)).symm)).symm
  have hcard : (Nat.card G' : ℝ) = (Nat.card G : ℝ) := by
    exact_mod_cast (Nat.card_congr e.toEquiv).symm
  simp only [herbrandPhiGroup, phiOf]
  rw [hsum, hcard]

/-- `ψ = φ^{−1}` も同型で移る。 -/
theorem herbrandPsiGroup_congr_equiv [Fintype G] [Fintype G'] (e : G ≃* G') (ψ : C ≃+* C')
    (hcompat : ∀ (σ : G) (c : C), ψ (σ • c) = e σ • ψ c) (ϖ : C) (m : ℝ) :
    herbrandPsiGroup G' (ψ ϖ) m = herbrandPsiGroup G ϖ m := by
  have hfun : herbrandPhiGroup G' (ψ ϖ) = herbrandPhiGroup G ϖ :=
    funext fun n => herbrandPhiGroup_congr_equiv e ψ hcompat ϖ n
  unfold herbrandPsiGroup
  rw [hfun]

/-- ★★**上付き分岐群は(作用を保つ)同型で移る**。 -/
theorem upperRamificationGroup_comap_equiv [Fintype G] [Fintype G'] (e : G ≃* G')
    (ψ : C ≃+* C') (hcompat : ∀ (σ : G) (c : C), ψ (σ • c) = e σ • ψ c) (ϖ : C) (m : ℝ) :
    upperRamificationGroup G ϖ m
      = Subgroup.comap (e : G →* G') (upperRamificationGroup G' (ψ ϖ) m) := by
  ext σ
  rw [Subgroup.mem_comap, mem_upperRamificationGroup, mem_upperRamificationGroup,
    MonoidHom.coe_coe, ramIndex_congr_equiv e ψ hcompat, herbrandPsiGroup_congr_equiv e ψ hcompat]

/-- ★★★**「残っている穴」の消費口** —— 移った先で `⊥` なら元でも `⊥`。

★これに `e : G ⧸ H ≃* Gal(K^m_x/K)` と `ψ : fixedRing 𝒪_{K(x)} H ≃+* 𝒪_{K^m_x}` を
渡せば、B4 の結論が §4 の `hbot` になる。 -/
theorem upperRamificationGroup_eq_bot_of_equiv [Fintype G] [Fintype G'] (e : G ≃* G')
    (ψ : C ≃+* C') (hcompat : ∀ (σ : G) (c : C), ψ (σ • c) = e σ • ψ c) (ϖ : C) (m : ℝ)
    (h : upperRamificationGroup G' (ψ ϖ) m = ⊥) : upperRamificationGroup G ϖ m = ⊥ := by
  rw [upperRamificationGroup_comap_equiv e ψ hcompat, h, MonoidHom.comap_bot,
    (MonoidHom.ker_eq_bot_iff (e : G →* G')).2 e.injective]

end Transport

/-! ## §3・§4 具体層 -/

section Compositum

variable (K : PAdicLocalField p) (x : K.closure)
  [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-- 完全分岐なら剰余体の濃度は底のまま —— `q_{K(x)} = q_K`。

★在庫 `isTotallyRamifiedAdjoin_iff_residueDegree` の左から右。
これで B5 の `hdeg` に出てくる `q` を、原典どおり**底の体 `K` の剰余体の濃度**で
書けるようになる。 -/
theorem natCard_residueField_adjoinIntegers_eq (ht : IsTotallyRamifiedAdjoin K x) :
    Nat.card (IsLocalRing.ResidueField (adjoinIntegers K x)) = Nat.card 𝓀[K.carrier] :=
  (isTotallyRamifiedAdjoin_iff_residueDegree K x).1 ht

omit [FiniteDimensional K.carrier
  (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] in
/-- ★★**`|Gal(K(x)/K) ⧸ Gal(K(x)/E)| = [E : K]`**。

★★**mathlib の `IntermediateField.finrank_eq_fixingSubgroup_index` は
`[Normal K E]` を要求しない**(`[IsGalois K L]` だけでよい)。そのため
`AlgEquiv.restrictNormalHom` の全射性・核を経由する必要がなく、
`lean-idioms.md` #154(`restrictNormalHom` に `exact` を当てると heartbeats で落ちる)にも
触れずに済む。★2 層の中間体は `IntermediateField.restrict` に閉じ込め、
mathlib の `restrict_algEquiv` だけで `[E : K]` に戻している(#59 の定型 (c))。 -/
theorem natCard_quotient_fixingSubgroup_restrict
    [IsGalois K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {E : IntermediateField K.carrier K.closure}
    (h : E ≤ IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    Nat.card ((((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
        ⧸ (IntermediateField.restrict h).fixingSubgroup)
      = Module.finrank K.carrier E := by
  rw [← Subgroup.index_eq_card, ← IntermediateField.finrank_eq_fixingSubgroup_index]
  exact (LinearEquiv.finrank_eq (IntermediateField.restrict_algEquiv h).toLinearEquiv).symm

def upperRamificationGroupAdjoin_le_of_quotient_upper_eq_bot.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Corollary 6.13", sectionId := "cor-6-13" }

/-- ★★★★**Corollary 6.13 (ii) の前半、`htriv` を外した形**。

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

B5 の `upperRamificationGroupAdjoin_le_of_quotient_eq_bot` は
`htriv : upperRamificationGroup G ϖ m = H` を仮定として受け取っていたが、
本補題はそれを **`hbot : (G ⧸ H)^m = ⊥`**(原典の `Gal(K′′/K)^m = {id}` そのもの)に
置き換える。★橋は §1 の `upperRamificationGroup_eq_of_quotient_eq_bot`。

★`[MulSemiringAction (G ⧸ H) C]` と `hq` は在庫 `ramIndex_quotient_mk` と同じ流儀で
仮定として受け取っている(冒頭「在庫調査の記録」)。消費側は
`letI := quotientMulSemiringActionOfTrivial _ smul_fixedRing_eq_self` と
`hq := quotientSMul_mk_fixedRing` を渡す。 -/
theorem upperRamificationGroupAdjoin_le_of_quotient_upper_eq_bot
    (ht : IsTotallyRamifiedAdjoin K x)
    {α : adjoinIntegers K x}
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α})
    [Fintype ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))]
    {H : Subgroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))}
    [H.Normal] [Fintype ↥H]
    [Fintype (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ⧸ H)]
    [MulSemiringAction ((((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))) ⧸ H)
      ↥(fixedRing (adjoinIntegers K x) H)]
    (hq : ∀ (σ : ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
      (c : ↥(fixedRing (adjoinIntegers K x) H)),
      (QuotientGroup.mk σ : (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))) ⧸ H)
        • c = σ • c)
    {ϖ : ↥(fixedRing (adjoinIntegers K x) H)} (hϖ : Irreducible ϖ) {m : ℝ}
    (hbot : upperRamificationGroup ((((IntermediateField.adjoin K.carrier
      ({x} : Set K.closure)) ≃ₐ[K.carrier]
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))) ⧸ H) ϖ m = ⊥) :
    upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) α m ≤ H :=
  upperRamificationGroupAdjoin_le_of_quotient_eq_bot K x ht huni hϖ
    (upperRamificationGroup_eq_of_quotient_eq_bot hq ϖ m hbot)

def le_of_two_quotients_upperRamification_of_finrank.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★★★**Yoshida 2008 Theorem 6.15 の "thus K′ ⊂ Km_x"、`htriv` と `hdeg` を
原典の字面に戻した形**。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

B5 の `le_of_two_quotients_upperRamification` から

* `htriv` / `htriv′`(`upperRamificationGroup G ϖ m = H` という**引き戻しの形**)は
  `hbot` / `hbot′`(**商の上付き分岐群が `⊥`** = 原典の `Gal(K′/L)^m = Gal(K^m_x/L)^m = {id}`)に、
* `hdeg`(`Nat.card (G ⧸ H) = (q−1)q^k` という**群の位数**の形)は
  `hdeg`(`[E₁ : K] = (q−1)q^k` という**体の次数**の形、`q` は**底の `K` の剰余体の濃度**
  = 原典の `[K^m_x : L] = (q−1)q^{m−1}`)に

置き換わっている。★★**どちらも群の言葉から原典の字面に戻っている。**

★次数の橋は §3 の `natCard_quotient_fixingSubgroup_restrict`、
`q` の橋は §3 の `natCard_residueField_adjoinIntegers_eq`(完全分岐 `ht` を使う)。 -/
theorem le_of_two_quotients_upperRamification_of_finrank
    (ht : IsTotallyRamifiedAdjoin K x)
    [IsGalois K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {α : adjoinIntegers K x}
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α})
    [Fintype ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))]
    (habel : ∀ σ τ : ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))),
      σ * τ = τ * σ)
    {E₁ E₂ : IntermediateField K.carrier K.closure}
    (h₁ : E₁ ≤ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (h₂ : E₂ ≤ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (hsup : E₁ ⊔ E₂ = IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [((IntermediateField.restrict h₁).fixingSubgroup).Normal]
    [((IntermediateField.restrict h₂).fixingSubgroup).Normal]
    [Fintype (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      ⧸ (IntermediateField.restrict h₁).fixingSubgroup)]
    [Fintype (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      ⧸ (IntermediateField.restrict h₂).fixingSubgroup)]
    [MulSemiringAction ((((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
      ⧸ (IntermediateField.restrict h₁).fixingSubgroup)
      ↥(fixedRing (adjoinIntegers K x) (IntermediateField.restrict h₁).fixingSubgroup)]
    [MulSemiringAction ((((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
      ⧸ (IntermediateField.restrict h₂).fixingSubgroup)
      ↥(fixedRing (adjoinIntegers K x) (IntermediateField.restrict h₂).fixingSubgroup)]
    (hq : ∀ (σ : ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
      (c : ↥(fixedRing (adjoinIntegers K x) (IntermediateField.restrict h₁).fixingSubgroup)),
      (QuotientGroup.mk σ : (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
        ⧸ (IntermediateField.restrict h₁).fixingSubgroup) • c = σ • c)
    (hq' : ∀ (σ : ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
      (c : ↥(fixedRing (adjoinIntegers K x) (IntermediateField.restrict h₂).fixingSubgroup)),
      (QuotientGroup.mk σ : (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
        ⧸ (IntermediateField.restrict h₂).fixingSubgroup) • c = σ • c)
    {ϖ : ↥(fixedRing (adjoinIntegers K x) (IntermediateField.restrict h₁).fixingSubgroup)}
    (hϖ : Irreducible ϖ)
    {ϖ' : ↥(fixedRing (adjoinIntegers K x) (IntermediateField.restrict h₂).fixingSubgroup)}
    (hϖ' : Irreducible ϖ') (k : ℕ)
    (hbot : upperRamificationGroup ((((IntermediateField.adjoin K.carrier
        ({x} : Set K.closure)) ≃ₐ[K.carrier]
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
        ⧸ (IntermediateField.restrict h₁).fixingSubgroup) ϖ ((k + 1 : ℕ) : ℝ) = ⊥)
    (hbot' : upperRamificationGroup ((((IntermediateField.adjoin K.carrier
        ({x} : Set K.closure)) ≃ₐ[K.carrier]
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))
        ⧸ (IntermediateField.restrict h₂).fixingSubgroup) ϖ' ((k + 1 : ℕ) : ℝ) = ⊥)
    (hdeg : Module.finrank K.carrier E₁
      = (Nat.card 𝓀[K.carrier] - 1) * Nat.card 𝓀[K.carrier] ^ k) :
    E₂ ≤ E₁ := by
  haveI : Fintype ↥((IntermediateField.restrict h₁).fixingSubgroup) := Fintype.ofFinite _
  haveI : Fintype ↥((IntermediateField.restrict h₂).fixingSubgroup) := Fintype.ofFinite _
  refine le_of_two_quotients_upperRamification K x ht huni habel h₁ h₂ hsup hϖ hϖ' k
    (upperRamificationGroup_eq_of_quotient_eq_bot hq ϖ _ hbot)
    (upperRamificationGroup_eq_of_quotient_eq_bot hq' ϖ' _ hbot') ?_
  rw [natCard_quotient_fixingSubgroup_restrict K x h₁,
    natCard_residueField_adjoinIntegers_eq K x ht]
  exact hdeg

end Compositum

def finrank_adjoin_iteratedLubinTatePsi_succ.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★**原典の `(q −1)qm−1 = [Km_x : L]`(`n = m = k+1`、`L = K`)**。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

在庫 `finrank_adjoin_iteratedLubinTatePsi`(`[K(Λ_n) : K] = q^n − q^{n−1}`)を
`n = k+1` に固定し、`q^{k+1} − q^k = (q−1)q^k` と書き直しただけ。
★★これが上の `le_of_two_quotients_upperRamification_of_finrank` の `hdeg` を
**そのまま埋める**。★`ℕ` の引き算は `q^{k+1} − q^k`(`q ≥ 1` なので安全)であって
`ℕ∞` ではない(`lean-idioms.md` #102)。 -/
theorem finrank_adjoin_iteratedLubinTatePsi_succ
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hqcard : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (k : ℕ) (x₀ : K.closure)
    (hxψ : x₀ ∈ iteratedLubinTatePsiTorsionPoints K hqcard hπmax hπne0 f hf0 hf1 hf (k + 1)
      (Nat.le_add_left 1 k))
    (hxn : x₀ ∈ iteratedLubinTateTorsionPoints K hqcard hπmax hπne0 f hf0 hf1 hf (k + 1))
    (hmem : x₀ ∈ IntermediateField.adjoin K.carrier ({x₀} : Set K.closure))
    [FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({x₀} : Set K.closure))] :
    Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x₀} : Set K.closure))
      = (Nat.card 𝓀[K.carrier] - 1) * Nat.card 𝓀[K.carrier] ^ k := by
  have h := finrank_adjoin_iteratedLubinTatePsi K hqcard hπmax hπne0 f hf0 hf1 hf (k + 1)
    (Nat.le_add_left 1 k) x₀ hxψ hxn hmem
  have hcard : Nat.card 𝓀[K.carrier] = pp ^ ff := by
    rw [Nat.card_eq_fintype_card]; exact hqcard
  rw [h, hcard, Nat.add_sub_cancel, Nat.sub_mul, one_mul, Nat.mul_comm (pp ^ ff)
    ((pp ^ ff) ^ k), pow_succ]

/-! ## 残っている穴(★次のノードになるもの。1 本だけ)

★★**`fixedRing (adjoinIntegers K x) H ≅ adjoinIntegers K x₀` の同一視**
(`H = Gal(K(x)/K(x₀))`)。作用と `Algebra 𝒪_K` を保つ環同型と、
群同型 `G ⧸ H ≃* Gal(K(x₀)/K)`(mathlib の
`AlgEquiv.restrictNormalHom_surjective` + `IntermediateField.restrictNormalHom_ker` +
`QuotientGroup.quotientKerEquivOfSurjective`)が要る。

★★**それが来たときの受け口は本ファイルの §2 に用意してある**:
`upperRamificationGroup_eq_bot_of_equiv` に `e` と `ψ` と `hcompat` を渡せば、
B4 の結論 `upperRamificationGroup Gal(K^m_x/K) α₀ m = ⊥`
(`Found/PGC/LubinTateUpperRamificationVanish.lean` の
`upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap`)が
§4 の `hbot` になり、`htriv` は完全に外れる。

★環同型の中身は「`𝒪_{K(x)}` の `H` 不変元がちょうど `𝒪_{K(x₀)}`」という主張であり、
`adjoinField` / `adjoinIntegers` の境界(`lean-idioms.md` #69)に触れる。
★**不変部分環の側に寄せて書く**のが安全である(B5 の判断と同じ)。 -/

end ABC3.Found.PGC
