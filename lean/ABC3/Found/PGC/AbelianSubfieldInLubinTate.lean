import ABC3.Found.PGC.LubinTateUpperRamificationVanish
import ABC3.Found.PGC.UpperRamificationIndex
import ABC3.Found.PGC.AbelianJumpDivisibility
import ABC3.Found.PGC.FixedRingMonogenic
import ABC3.Found.PGC.AbelianClosureSplit

/-!
# `E_σ ⊂ K^m_x` —— [Yoshida08] Theorem 6.15 の最後の数学(道 B の B5)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Theorem 6.15**(物理 p.17)。
構造化済み原文は
`ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/section-6.html` の
`id="thm-6-15"`(`data-pdf-page="17"`, `data-item="Theorem 6.15 (Local Kronecker-Weber)"`)と
`id="cor-6-13"`(`data-pdf-page="17"`, `data-item="Corollary 6.13"`)。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

## 本ファイルが担当する原典の段落

Theorem 6.15 の Proof の末尾 4 文(`0_Source` の `.txt` 1195–1198 行)である
(◆抽出器の字形の潰れは直していない):

    Let K′/L be any finite Galois extension contained in Eσ. It is totally ramified, and
    Gal(K′/L)m = {id} for a large m. Then we have Gal(K′Km_x/L)m = {id} by Proposition 6.14
    and Corollary 6.13(ii), hence [K′Km_x : L] | (q −1)qm−1 = [Km_x : L] by Corollary 6.13(iii),
    thus K′ ⊂ Km_x.

★決定 D29 により `n = 1`、すなわち `L = K_1 = K` である。以下 `L = K` と読む。

## 段取り(原典の 4 文がそのまま §1–§6 に対応する)

| 原典の文 | 本ファイル |
|---|---|
| `E_σ ∩ K^ur = L` | §2 `fixedField_zpowers_inf_unramifiedClosure_eq_bot` |
| It is totally ramified | 在庫 `isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot` |
| Cor 6.13(ii) を合成体に当てる | §3 `upperRamificationGroupAdjoin_eq_bot_of_two_quotients` |
| Cor 6.13(iii) を合成体に当てる | §4 `natCard_quot_upperRamificationGroupAdjoin_succ_dvd` |
| 次数比較 ⇒ `K′ ⊂ K^m_x` | §6 `subgroup_eq_bot_of_two_quotients` / §7 `le_of_two_quotients_upperRamification` |

## ★★本ファイルの重心 —— Corollary 6.13(ii)(iii) の `C`/`C′` 供給層

在庫の `upperRamificationGroup_eq_bot_of_two_quotients`(Y18)は
**仮定 14 本**を持ち、そのうち `C`(= `𝒪_{K′′}`)と `C′`(= `𝒪_{K′}`)を
`Algebra C B` + `MulSemiringAction G C` + `hfixC` + `hfix` + `hAC` + `hϖ` で供給する層を
誰も書いていなかった。本ファイルの §3 がそれである。

★★**供給の決め手は「`C` を `fixedRing B H`(不変部分環)として取る」ことである。**
体の側で `𝒪_{K′′}` を作って `Algebra 𝒪_{K′′} 𝒪_{K′K′′}` を張ろうとすると
`adjoinField` / `adjoinIntegers` の境界(`lean-idioms.md` #69)に当たるが、
不変部分環なら `Algebra` も `MulSemiringAction` も既に木にある
(`FixedRingTower` / `FixedRingAction` / `FixedRingMonogenic`)。
14 本の仮定の供給元は §3 の docstring に一覧してある。

★★**もう 1 つの決め手は「底を `A := 𝒪_K` に取れること」である。**
`RamificationFiltrationBuild` の `stage*` 群は底を `𝒪_{L₀}`(惰性群の不変環)に取るが、
本ファイルの設定では `K′K^m_x/K` が**完全分岐**なので `L₀ = K` であり、
`A := 𝒪[K.carrier]`・`G := Gal(K(x)/K)` のまま `Algebra.adjoin A {α} = ⊤` が
`adjoin_uniformizer_eq_top_adjoinIntegers` で出る。★これは B4 と同じ流儀である。

## ★★★持ち込んではいけない近道(反例)

Corollary 6.13(ii) は `K′K′′/K` が**完全分岐**であることを仮定する。
`E_σ` を捨てて「2 つの完全分岐アーベル拡大の合成」で代用すると**壊れる**:
`K = ℚ_p` で `ℚ_p(√p)` と `ℚ_p(√(up))`(`u` は非平方単数)はどちらも完全分岐だが、
その合成は `ℚ_p(√u)`(不分岐)を含む。★**完全分岐性は合成で保たれない。**
本ファイルが壊れないのは、原典と同じく `K′` も `K^m_x` も**同じ `E_σ` の中**に取り、
その合成体 `K(x)` の完全分岐性 `ht : IsTotallyRamifiedAdjoin K x` を
**仮定として明示的に受け取る**からである(§3・§4・§6・§7 の `ht`)。

## ★`hHH : H ⊓ H′ = ⊥` をどこから出したか

木での読み替えは `upperRamificationGroup_eq_bot_of_two_quotients` の
`hHH : H ⊓ H′ = ⊥` である。本ファイルは §5 の
`inf_fixingSubgroup_restrict_eq_bot` で供給する。すなわち

* `H := Gal(K(x)/E₁)`、`H′ := Gal(K(x)/E₂)`(どちらも `IntermediateField.fixingSubgroup`)、
* `E₁ ⊔ E₂ = K(x)`(合成体である、という仮定 `hsup`)

から、在庫 `fixingSubgroup_sup`(`Found/PGC/CompositumSurjection.lean`)と
mathlib の `IntermediateField.fixingSubgroup_top` で出す。
★★**環論版の `inf_eq_bot_of_closure_eq_top`(`𝒪_{K′K′′} = 𝒪_{K′}·𝒪_{K′′}`)は使っていない** ——
整数環が合成で生成されるとは限らないので、体の側の Galois 対応で出すほうが安全である。

## 逸脱の記録

1. **`L = K`(`n = 1`)に固定した。** 決定 D29 と B4(`LubinTateUpperRamificationVanish`)に
   合わせた。原典の `L = K_n` 一般は扱っていない。
2. **`E_σ` そのものではなく、`E_σ` の中の有限次部分拡大を扱う。** これは逸脱ではなく
   原典の段取りそのもの(「Let K′/L be any finite Galois extension contained in E_σ」)。
   無限次の `E_σ` に直接当てる必要はない。
3. **`htriv` / `htriv′`(原典の `Gal(K′/L)^m = Gal(K^m_x/L)^m = {id}`)は仮定である。**
   これは Y18 の `upperRamificationGroup_eq_bot_of_two_quotients` が要求する
   「引き戻しの形」`upperRamificationGroup G ϖ m = H` そのままである。
   ★B4 の結論 `upperRamificationGroup Gal(K(x₀)/K) α₀ m = ⊥` は**群も環も違う**
   (`Gal(K^m_x/K)` が `𝒪_{K^m_x}` に作用する形)ので、そこから `htriv` を出すには
   「`G` が `C` に `G/H` を経由して作用するとき `G^m = 引き戻し((G/H)^m)`」という
   降下補題と、`fixedRing B H ≅ 𝒪_{K^m_x}` の同一視が要る。★**本ファイルには入っていない**
   (末尾「残っている穴」を見よ)。
4. **`hdeg`(原典の `[K^m_x : L] = (q−1)q^{m−1}`)は仮定である。** 木では
   `finrank_adjoin_iteratedLubinTatePsi` あるいは
   `galoisReciprocityEquiv` + `card_units_quotient_span_pi_pow` から出るが、
   どちらも `K^m_x = K(x₀)` の言葉であり、本ファイルの `H = Gal(K(x)/K^m_x)` と
   結び付けるには逸脱 3 と同じ同一視が要る。★`Nat.card (G ⧸ H)` の形で受け取っている。
5. **`m` は `k + 1` の形に固定した**(`ℕ` の切り詰め引き算・`ℕ∞` の除算を書かないため。
   `lean-idioms.md` #102、B4 と同じ逃げ方)。原典の `m ≥ 1` はこれで尽くされる。

## 退化の自己検査

* ★**上付きと下付きを取り違えていない**: 本ファイルが `⊥` にするのは
  `upperRamificationGroup … α ((k+1 : ℕ) : ℝ)`(= `G^{k+1}`)であって
  `lowerRamificationGroup … (k+1)`(= `G_{k+1}`)ではない。
  下付きが出てくるのは `lowerRamificationGroupAdjoin_zero_eq_top`(`G_0 = ⊤`、完全分岐)
  を Corollary 6.13(iii) の `h0` に渡すところだけである。
* ★**結論は包含 `E₂ ≤ E₁` であって等号ではない**(原典の `K′ ⊂ K^m_x`)。
  等号 `K^LT = K^ab` は B6 で `K^LT ≤ K^ab` と併せて初めて出る。
* ★**完全分岐性 `ht` を落とすと偽**(上の `ℚ_p(√p)`・`ℚ_p(√(up))` の反例)。
  `ht` は `hadj : Algebra.adjoin 𝒪_K {α} = ⊤` と `hresA`(剰余体が伸びない)の
  両方に使われており、どちらも Corollary 6.13(ii)(iii) の必須仮定である。
* ★**`hHH` を落とすと偽**: `H = H′` と取れば `H ⊓ H′ = H` で、
  `G^m ≤ H` だけからは `G^m = ⊥` は出ない。
* ★**`habel`(可換性)を落とすと Corollary 6.13(iii) が偽**(Hasse-Arf が使えない)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open IsLocalRing IsDiscreteValuationRing
open scoped NormedField Valued Classical

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 抽象核 —— 分岐・付値・Lubin-Tate の語彙が 1 つも出てこない 4 本

★いずれも mathlib だけで書けて、`lean_check` が 0.10–0.23 秒で返った。 -/

/-- ★★**抽象核 1** —— **稠密な部分群の固定体は底体**。

`Dense (H : Set Gal(L/k))` なら `L^H = k`。原典が
「`Gal(K^ab/E_σ) ≅ Ẑ` with `σ ↦ 1` by the definition of `E_σ`」と書いて
`Ẑ` の Hopf 性を経由するところを、**`Ẑ` を一切使わずに**位相だけで済ませる。

段取り: `x ∈ L^H` を取ると `{g | g x = x} = k(x).fixingSubgroup` は
Krull 位相で閉(`IntermediateField.fixingSubgroup_isClosed`、`k(x)/k` は有限次)。
これが稠密な `H` を含むので全体に等しく、`x ∈ L^⊤ = k`
(`InfiniteGalois.fixedField_bot`)。 -/
theorem fixedField_eq_bot_of_dense {k L : Type*} [Field k] [Field L] [Algebra k L]
    [IsGalois k L] (H : Subgroup (L ≃ₐ[k] L))
    (hdense : Dense (H : Set (L ≃ₐ[k] L))) :
    IntermediateField.fixedField H = ⊥ := by
  rw [← InfiniteGalois.fixedField_bot (k := k) (K := L)]
  refine le_antisymm ?_ (IntermediateField.fixedField_le le_top)
  intro y hy
  have hint : IsIntegral k y := Algebra.IsIntegral.isIntegral y
  have hfin : FiniteDimensional k (IntermediateField.adjoin k ({y} : Set L)) :=
    IntermediateField.adjoin.finiteDimensional hint
  have hadjle : IntermediateField.adjoin k ({y} : Set L) ≤ IntermediateField.fixedField H :=
    IntermediateField.adjoin_le_iff.2 (Set.singleton_subset_iff.2 hy)
  have hclosed : IsClosed (((IntermediateField.adjoin k ({y} : Set L)).fixingSubgroup :
      Set (L ≃ₐ[k] L))) :=
    IntermediateField.fixingSubgroup_isClosed _
  have hsub : (H : Set (L ≃ₐ[k] L))
      ⊆ (IntermediateField.adjoin k ({y} : Set L)).fixingSubgroup := by
    intro g hg
    simp only [SetLike.mem_coe, IntermediateField.mem_fixingSubgroup_iff]
    intro z hz
    exact hadjle hz ⟨g, hg⟩
  have huniv : (Set.univ : Set (L ≃ₐ[k] L)) ⊆
      (IntermediateField.adjoin k ({y} : Set L)).fixingSubgroup := by
    rw [← hdense.closure_eq]
    exact hclosed.closure_subset_iff.2 hsub
  intro g
  have hg := huniv (Set.mem_univ g.1)
  simp only [SetLike.mem_coe, IntermediateField.mem_fixingSubgroup_iff] at hg
  exact hg y (IntermediateField.subset_adjoin k _ rfl)

/-- ★★**抽象核 2** —— 原典の

    [K′K^m_x : L] | (q−1)q^{m−1} = [K^m_x : L], thus K′ ⊂ K^m_x

の**次数の部分だけ**を取り出したもの。`M ⊔ N` の次数が `N` の次数を割るなら `M ≤ N`。

★本ファイルの主定理は次数ではなく**群の位数**でこの段を通しているが
(§6、そちらのほうが `FiniteDimensional` の instance を並べずに済む)、
原典の字面どおりの道もここに残しておく。★どちらでも結論は同じである。 -/
theorem le_of_finrank_sup_dvd {F L : Type*} [Field F] [Field L] [Algebra F L]
    (M N : IntermediateField F L) [FiniteDimensional F ↥(M ⊔ N)]
    (h : Module.finrank F ↥(M ⊔ N) ∣ Module.finrank F N) : M ≤ N := by
  haveI : FiniteDimensional F N := by
    refine Module.Finite.of_injective
      (IntermediateField.inclusion (le_sup_right : N ≤ M ⊔ N)).toLinearMap ?_
    intro a b hab
    exact (IntermediateField.inclusion (le_sup_right : N ≤ M ⊔ N)).injective hab
  have hpos : 0 < Module.finrank F N := Module.finrank_pos
  have heq : N = M ⊔ N :=
    IntermediateField.eq_of_le_of_finrank_le le_sup_right (Nat.le_of_dvd hpos h)
  exact heq ▸ le_sup_left

/-- ★★**抽象核 3(純群論)** —— 原典の次数比較を**群の位数**で行う形。

`|G|` が `d` を割り、`|G/H| = d` なら `H = ⊥`。
`|H| · |G/H| = |G|` から `|H| · d ∣ d`、`d > 0` なので `|H| = 1`。

★原典が `[K′K^m_x : L] | (q−1)q^{m−1} = [K^m_x : L]` から `K′ ⊂ K^m_x` を出すところ。
`(q−1)q^{m−1}` が具体的に何であるかは使っていない(`d` は任意)。 -/
theorem eq_bot_of_natCard_dvd_of_natCard_quotient {G : Type*} [Group G] [Finite G]
    (H : Subgroup G) (d : ℕ) (hd : Nat.card G ∣ d) (hq : Nat.card (G ⧸ H) = d) : H = ⊥ := by
  have hdpos : 0 < d := by rw [← hq]; exact Nat.card_pos
  have hmul : Nat.card H * H.index = Nat.card G := Subgroup.card_mul_index H
  rw [Subgroup.index_eq_card, hq] at hmul
  have hdvd2 : Nat.card H * d ∣ d := hmul ▸ hd
  have heq : Nat.card H * d = d := Nat.dvd_antisymm hdvd2 (dvd_mul_left d (Nat.card H))
  have h1 : Nat.card H = 1 := by
    refine Nat.eq_of_mul_eq_mul_right hdpos ?_
    rw [heq, one_mul]
  exact Subgroup.eq_bot_of_card_eq _ h1

/-- ★**抽象核 4** —— `E ⊔ F = ⊤` なら `Gal(L/E) ∩ Gal(L/F) = 1`。

★これが Y18 の `hHH : H ⊓ H′ = ⊥` の出所である(本ファイル冒頭を見よ)。
在庫 `fixingSubgroup_sup`(`Found/PGC/CompositumSurjection.lean`)と
mathlib の `IntermediateField.fixingSubgroup_top` の合成。 -/
theorem inf_fixingSubgroup_eq_bot_of_sup_eq_top {k L : Type*} [Field k] [Field L] [Algebra k L]
    {E F : IntermediateField k L} (h : E ⊔ F = ⊤) :
    E.fixingSubgroup ⊓ F.fixingSubgroup = ⊥ := by
  rw [← fixingSubgroup_sup, h, IntermediateField.fixingSubgroup_top]

/-! ## §2 B5-0 —— `E_σ ∩ K^ur = L` -/

def fixedField_zpowers_inf_unramifiedClosure_eq_bot.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★**Yoshida 2008 Theorem 6.15 の "Then Eσ ∩ Kur = L"**(`n = 1` なので `L = K`)。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

`σ|_{K^ur}` が算術 Frobenius なら、`σ` の固定体は `K^ur` と自明にしか交わらない。

段取りは 3 行:
1. `InfiniteGalois.restrict_fixedField` で `E_σ ⊓ K^ur = lift (fixedField (⟨σ⟩ の像))`。
2. `MonoidHom.map_zpowers` で像は `⟨算術 Frobenius⟩`。
3. 在庫 `dense_zpowers_frobenius` が稠密性を与えるので §1 の抽象核 1 で `⊥`。

★★**`Ẑ` を経由していない**(原典は `Gal(K^ab/E_σ) ≅ Ẑ` を使う)。
★これで在庫 `isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot` が使えるようになり、
原典の "It is totally ramified" が木でも出る。 -/
theorem fixedField_zpowers_inf_unramifiedClosure_eq_bot (K : PAdicLocalField p) {σ : K.absGal}
    (hσ : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K) :
    IntermediateField.fixedField (Subgroup.zpowers σ) ⊓ unramifiedClosure K = ⊥ := by
  haveI := isGalois_closure K
  haveI := normal_unramifiedClosure K
  haveI := isGalois_unramifiedClosure K
  have hmap : Subgroup.map
      (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) ↥(unramifiedClosure K))
      (Subgroup.zpowers σ) = Subgroup.zpowers (arithFrobenius K) := by
    rw [MonoidHom.map_zpowers, hσ]
  have hdense : Dense ((Subgroup.zpowers (arithFrobenius K) : Subgroup (unramGal K)) :
      Set (unramGal K)) := by
    rw [Subgroup.coe_zpowers]
    exact dense_zpowers_frobenius K (fun N hN => arithFrobenius_mem_unramLevelGeneratorSet K N hN)
  rw [InfiniteGalois.restrict_fixedField (Subgroup.zpowers σ) (unramifiedClosure K), hmap,
    fixedField_eq_bot_of_dense _ hdense]
  simp

section Compositum

variable (K : PAdicLocalField p) (x : K.closure)
  [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-! ## §3 B5-1 —— Corollary 6.13(ii) の `C`/`C′` 供給層

★★**ここが本持ち場の工数の全部である。** 在庫の Y18
`upperRamificationGroup_le_of_quotient_eq_bot` / `upperRamificationGroup_eq_bot_of_two_quotients`
は仮定を 11 本 / 14 本持っており、`C`(不変部分環)を渡す層が無かった。 -/

/-- ★★★★**Corollary 6.13(ii) の前半の供給層** ——
`(G/H)^m = {id}`(引き戻しの形)なら `G^m ⊆ H`。

`B := 𝒪_{K(x)}`、`A := 𝒪_K`、`G := Gal(K(x)/K)`、`C := 𝒪_{K(x)}^H`(不変部分環)。

★★**Y18 の 11 本の仮定をどう供給したか**(これが本補題の中身である):

| Y18 の仮定 | 供給元 |
|---|---|
| `hcomp` : `ι(ρ • c) = ρ • ι c` | `algebraMap_smul_fixedRing`(`FixedRingAction`) |
| `hHtriv` : `H` は `C` に自明に作用 | `smul_fixedRing_eq_self`(同上) |
| `hπ′` : `α` は `B` の素元 | `irreducible_iff_uniformizer` を `huni` に |
| `hinj` : `C → B` が単射 | `fixedRing_injective`(`FixedRingTower`) |
| `hfixC` : `H` 不変元は `C` から来る | `exists_algebraMap_fixedRing`(同上) |
| `hres` : `B` の元は `C` の元と `𝔪_B` を法として合同 | `exists_sub_mem_fixedRing` + `hresA` |
| `hAC` : `A` の像は `C` の像に入る | `exists_algebraMap_fixedRing_eq`(同上) |
| `hϖ` : `ι ϖ = π″` | `rfl`(`π″` を `ι ϖ` と取る) |
| `hadj` : `𝒪_K[α] = 𝒪_{K(x)}` | `adjoin_uniformizer_eq_top_adjoinIntegers`(★`ht` が要る) |
| `hfix` : `H` 不変元は `A[ι ϖ]` に入る | `fixedRing_mem_adjoin_uniformizer`(`FixedRingMonogenic`) |
| `[Fintype G]` / `[Fintype H]` | 仮定(`K(x)/K` は有限次) |

`hresA`(剰余体が伸びない)は完全分岐 `ht` から
`exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin` で出る。
★★**`ht` を落とすと `hadj` も `hresA` も出ない**(冒頭の反例)。 -/
theorem upperRamificationGroupAdjoin_le_of_quotient_eq_bot
    (ht : IsTotallyRamifiedAdjoin K x)
    {α : adjoinIntegers K x}
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α})
    [Fintype ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))]
    {H : Subgroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))}
    [H.Normal] [Fintype ↥H]
    {ϖ : ↥(fixedRing (adjoinIntegers K x) H)} (hϖ : Irreducible ϖ) {m : ℝ}
    (htriv : upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ϖ m = H) :
    upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) α m ≤ H := by
  haveI := isDiscreteValuationRing_carrierIntegers K
  haveI := module_finite_adjoinIntegers K x
  haveI : IsNoetherian 𝒪[K.carrier] (adjoinIntegers K x) :=
    isNoetherian_of_isNoetherianRing_of_finite _ _
  have hresA := exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin K x ht
  have hAne : ∃ a ∈ IsLocalRing.maximalIdeal 𝒪[K.carrier],
      algebraMap 𝒪[K.carrier] (adjoinIntegers K x) a ≠ 0 :=
    exists_mem_maximalIdeal_map_ne_zero (injective_algebraMap_adjoinIntegers K x)
  exact upperRamificationGroup_le_of_quotient_eq_bot
    (A := 𝒪[K.carrier])
    (fun ρ c => algebraMap_smul_fixedRing ρ c) smul_fixedRing_eq_self
    ((IsDiscreteValuationRing.irreducible_iff_uniformizer _).2 huni)
    fixedRing_injective exists_algebraMap_fixedRing
    (fun b => exists_sub_mem_fixedRing hresA b) exists_algebraMap_fixedRing_eq rfl
    (adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni)
    (fixedRing_mem_adjoin_uniformizer hresA hAne hϖ) htriv

def upperRamificationGroupAdjoin_eq_bot_of_two_quotients.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Corollary 6.13", sectionId := "cor-6-13" }

/-- ★★★★★**Yoshida 2008 Corollary 6.13 (ii) の `PAdicLocalField` 版**。

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

`K(x)` が原文の合成体 `K′K′′`(本ファイルの用途では `K′K^m_x`)、
`H = Gal(K(x)/K′′)`、`H′ = Gal(K(x)/K′)`。
`ht` が原文の「with `K′K′′/K` totally ramified」である。

★段取りは原文どおり: 上の供給層を 2 回使って `G^m ≤ H` と `G^m ≤ H′` を出し、
`hHH : H ⊓ H′ = ⊥` から `G^m = ⊥`(在庫 `eq_bot_of_le_of_inf_eq_bot`)。
★`hHH` は §5 の `inf_fixingSubgroup_restrict_eq_bot` が供給する。 -/
theorem upperRamificationGroupAdjoin_eq_bot_of_two_quotients
    (ht : IsTotallyRamifiedAdjoin K x)
    {α : adjoinIntegers K x}
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α})
    [Fintype ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))]
    {H H' : Subgroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))}
    [H.Normal] [H'.Normal] [Fintype ↥H] [Fintype ↥H']
    {ϖ : ↥(fixedRing (adjoinIntegers K x) H)} (hϖ : Irreducible ϖ)
    {ϖ' : ↥(fixedRing (adjoinIntegers K x) H')} (hϖ' : Irreducible ϖ')
    (hHH : H ⊓ H' = ⊥) {m : ℝ}
    (htriv : upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ϖ m = H)
    (htriv' : upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ϖ' m = H') :
    upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) α m = ⊥ :=
  eq_bot_of_le_of_inf_eq_bot
    (upperRamificationGroupAdjoin_le_of_quotient_eq_bot K x ht huni hϖ htriv)
    (upperRamificationGroupAdjoin_le_of_quotient_eq_bot K x ht huni hϖ' htriv') hHH

/-! ## §4 B5-2 —— Corollary 6.13(iii) の `Adjoin` 版 -/

def natCard_quot_upperRamificationGroupAdjoin_succ_dvd.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Corollary 6.13", sectionId := "cor-6-13" }

/-- ★★★★**Yoshida 2008 Corollary 6.13 (iii) の `PAdicLocalField` 版**。

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

`G := Gal(K(x)/K)` が可換で `K(x)/K` が完全分岐なら
`|G/G^{k+1}|` は `(q−1)q^k` を割る(`q := |𝓀_{K(x)}|`)。

★Y18 の抽象版 `natCard_quot_upperRamificationGroup_succ_dvd` の 6 本の仮定の供給元:

* `hπ′` = `irreducible_iff_uniformizer` を `huni` に、
* `hresA` = `exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin`(★`ht`)、
* `hadj` = `adjoin_uniformizer_eq_top_adjoinIntegers`(★`ht`)、
* `hAinj` = `injective_algebraMap_adjoinIntegers`、
* `h0` = `lowerRamificationGroupAdjoin_zero_eq_top`(★`ht`。これだけが**下付き**)、
* `[CharP (ResidueField B) p]` = `charP_residueField_adjoinIntegers`。

★`habel` は原文の "Let G be abelian" で、落とすと Hasse-Arf が使えず偽になる。
★`m = k + 1` の形にしたのは `ℕ` の切り詰め引き算を避けるため(逸脱 5)。 -/
theorem natCard_quot_upperRamificationGroupAdjoin_succ_dvd
    (ht : IsTotallyRamifiedAdjoin K x)
    {α : adjoinIntegers K x}
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α})
    [Fintype ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))]
    (habel : ∀ σ τ : ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))),
      σ * τ = τ * σ)
    (k : ℕ) :
    Nat.card (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ⧸
        upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
          ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) α
          ((k + 1 : ℕ) : ℝ))
      ∣ (Nat.card (IsLocalRing.ResidueField (adjoinIntegers K x)) - 1)
          * Nat.card (IsLocalRing.ResidueField (adjoinIntegers K x)) ^ k := by
  haveI := isDiscreteValuationRing_carrierIntegers K
  haveI := module_finite_adjoinIntegers K x
  haveI : IsNoetherian 𝒪[K.carrier] (adjoinIntegers K x) :=
    isNoetherian_of_isNoetherianRing_of_finite _ _
  haveI := charP_residueField_adjoinIntegers K x
  haveI : Fintype ↥(lowerRamificationGroup (adjoinIntegers K x)
      ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) 1) :=
    Fintype.ofFinite _
  exact natCard_quot_upperRamificationGroup_succ_dvd (A := 𝒪[K.carrier]) p
    ((IsDiscreteValuationRing.irreducible_iff_uniformizer _).2 huni)
    (exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin K x ht)
    (adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni)
    (injective_algebraMap_adjoinIntegers K x)
    (lowerRamificationGroupAdjoin_zero_eq_top K x ht) habel k

/-! ## §5 抽象核(続き)—— 「`Gal(N/E₁) = 1` なら `E₂ ≤ E₁`」

★★`IntermediateField.restrict` / `lift` で 1 層に閉じてある(`lean-idioms.md` #59 回避)。
2 層をまたぐ `rfl` は 1 つも書いていない。 -/

/-- ★★**抽象核 5** —— `E₁ ≤ N` の `N` の中での固定部分群が自明なら `E₁ = N`。
したがって `N` に入るどんな `E₂` も `E₁` に入る。

`fixedField (fixingSubgroup F) = F`(`IsGalois.fixedField_fixingSubgroup`)と
`fixedField ⊥ = ⊤` から `F = ⊤`、`lift_restrict` / `lift_top` で `E₁ = N`。 -/
theorem le_of_fixingSubgroup_restrict_eq_bot {k L : Type*} [Field k] [Field L] [Algebra k L]
    {N : IntermediateField k L} [FiniteDimensional k N] [IsGalois k ↥N]
    {E₁ E₂ : IntermediateField k L} (h₁ : E₁ ≤ N) (h₂ : E₂ ≤ N)
    (hbot : (IntermediateField.restrict h₁).fixingSubgroup = ⊥) : E₂ ≤ E₁ := by
  have htop : IntermediateField.restrict h₁ = ⊤ := by
    have h := IsGalois.fixedField_fixingSubgroup (K := IntermediateField.restrict h₁)
    rw [hbot, IntermediateField.fixedField_bot] at h
    exact h.symm
  have hE : E₁ = N := by
    rw [← IntermediateField.lift_restrict (h := h₁), htop, IntermediateField.lift_top]
  rw [hE]; exact h₂

/-- ★**抽象核 6** —— `E₁ ⊔ E₂ = N` を `N` の中に降ろすと `⊤`。`lift` の単射性で出る。 -/
theorem sup_restrict_eq_top {k L : Type*} [Field k] [Field L] [Algebra k L]
    {N : IntermediateField k L} {E₁ E₂ : IntermediateField k L} (h₁ : E₁ ≤ N) (h₂ : E₂ ≤ N)
    (hsup : E₁ ⊔ E₂ = N) :
    IntermediateField.restrict h₁ ⊔ IntermediateField.restrict h₂ = ⊤ := by
  refine IntermediateField.lift_injective N ?_
  rw [IntermediateField.lift_sup, IntermediateField.lift_restrict, IntermediateField.lift_restrict,
    IntermediateField.lift_top, hsup]

/-- ★★★**`hHH : H ⊓ H′ = ⊥` の供給元**(本ファイル冒頭の約束)。

`N` が `E₁` と `E₂` の合成体なら `Gal(N/E₁) ∩ Gal(N/E₂) = 1`。 -/
theorem inf_fixingSubgroup_restrict_eq_bot {k L : Type*} [Field k] [Field L] [Algebra k L]
    {N : IntermediateField k L} {E₁ E₂ : IntermediateField k L} (h₁ : E₁ ≤ N) (h₂ : E₂ ≤ N)
    (hsup : E₁ ⊔ E₂ = N) :
    (IntermediateField.restrict h₁).fixingSubgroup
        ⊓ (IntermediateField.restrict h₂).fixingSubgroup = ⊥ :=
  inf_fixingSubgroup_eq_bot_of_sup_eq_top (sup_restrict_eq_top h₁ h₂ hsup)

/-! ## §6 主定理(群の言葉)—— 原典の 3 文をつなぐ -/

def subgroup_eq_bot_of_two_quotients.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★**Yoshida 2008 Theorem 6.15 の Proof の末尾 3 文**(群の言葉)。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

原典の該当箇所(`.txt` 1196–1198 行、抽出器の字形の潰れは直していない):

    Then we have Gal(K′Km_x/L)m = {id} by Proposition 6.14 and Corollary 6.13(ii),
    hence [K′Km_x : L] | (q −1)qm−1 = [Km_x : L] by Corollary 6.13(iii), thus K′ ⊂ Km_x.

`G := Gal(K(x)/K)`(`K(x)` が合成体 `K′K^m_x`)、`H := Gal(K(x)/K^m_x)`、
`H′ := Gal(K(x)/K′)`。結論 `H = ⊥` が原典の `K′ ⊂ K^m_x` である。

段取りは 3 行:
1. §3(Cor 6.13(ii))で `G^{k+1} = ⊥`。
2. §4(Cor 6.13(iii))で `|G/G^{k+1}|` が `(q−1)q^k` を割る。`G^{k+1} = ⊥` なので
   左辺は `|G|` そのもの(`QuotientGroup.quotientBot`)。
3. §1 の抽象核 3 に `d := (q−1)q^k` と `hdeg : |G/H| = (q−1)q^k` を入れて `H = ⊥`。

★★逸脱 3・4: `htriv` / `htriv′` / `hdeg` は仮定である(冒頭「逸脱の記録」)。
`htriv` は原典の `Gal(K^m_x/L)^m = {id}`(Prop 6.14 = B4)、
`hdeg` は原典の `[K^m_x : L] = (q−1)q^{m−1}` にあたる。 -/
theorem subgroup_eq_bot_of_two_quotients
    (ht : IsTotallyRamifiedAdjoin K x)
    {α : adjoinIntegers K x}
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α})
    [Fintype ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))]
    (habel : ∀ σ τ : ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))),
      σ * τ = τ * σ)
    {H H' : Subgroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))}
    [H.Normal] [H'.Normal] [Fintype ↥H] [Fintype ↥H']
    {ϖ : ↥(fixedRing (adjoinIntegers K x) H)} (hϖ : Irreducible ϖ)
    {ϖ' : ↥(fixedRing (adjoinIntegers K x) H')} (hϖ' : Irreducible ϖ')
    (hHH : H ⊓ H' = ⊥) (k : ℕ)
    (htriv : upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ϖ
      ((k + 1 : ℕ) : ℝ) = H)
    (htriv' : upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ϖ'
      ((k + 1 : ℕ) : ℝ) = H')
    (hdeg : Nat.card (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ⧸ H)
      = (Nat.card (IsLocalRing.ResidueField (adjoinIntegers K x)) - 1)
          * Nat.card (IsLocalRing.ResidueField (adjoinIntegers K x)) ^ k) :
    H = ⊥ := by
  have hbot := upperRamificationGroupAdjoin_eq_bot_of_two_quotients K x ht huni hϖ hϖ' hHH
    htriv htriv'
  have hdvd := natCard_quot_upperRamificationGroupAdjoin_succ_dvd K x ht huni habel k
  rw [hbot] at hdvd
  rw [Nat.card_congr (QuotientGroup.quotientBot (G := ((IntermediateField.adjoin K.carrier
    ({x} : Set K.closure)) ≃ₐ[K.carrier]
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure))))).toEquiv] at hdvd
  exact eq_bot_of_natCard_dvd_of_natCard_quotient H _ hdvd hdeg

/-! ## §7 主定理(体の言葉)—— `K′ ⊂ K^m_x` -/

def le_of_two_quotients_upperRamification.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★★★**Yoshida 2008 Theorem 6.15 の "thus K′ ⊂ Km_x"**(体の言葉)。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

`E₁ := K^m_x`、`E₂ := K′`、`K(x) := E₁ ⊔ E₂`(合成体、仮定 `hsup`)としたとき
`E₂ ≤ E₁`、すなわち `K′ ⊂ K^m_x`。

★★**`E₂` は `E_σ` に含まれる任意の有限次 Galois 拡大でよい**(原典の
「Let K′/L be any finite Galois extension contained in Eσ」)。
`E_σ` そのものは出てこない —— `E_σ` が効くのは
「合成体 `K(x)` が完全分岐」(`ht`)という 1 点だけであり、それは
`E_σ ⊓ K^ur = ⊥`(§2)と在庫 `isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot`
から出る。★冒頭の反例のとおり、`E_σ` を捨てて `ht` を落とすと**偽**になる。

★段取り: §5 で `hHH` を作り、§6 で `Gal(K(x)/E₁) = ⊥` を出し、
§5 の抽象核 5 で `E₁ = K(x)` に上げて `E₂ ≤ E₁`。

★逸脱 3・4: `htriv`(Prop 6.14 = B4 の結論)と `hdeg`(`[K^m_x : L] = (q−1)q^{m−1}`)は
仮定である。冒頭「逸脱の記録」と「残っている穴」を見よ。 -/
theorem le_of_two_quotients_upperRamification
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
    {ϖ : ↥(fixedRing (adjoinIntegers K x) (IntermediateField.restrict h₁).fixingSubgroup)}
    (hϖ : Irreducible ϖ)
    {ϖ' : ↥(fixedRing (adjoinIntegers K x) (IntermediateField.restrict h₂).fixingSubgroup)}
    (hϖ' : Irreducible ϖ') (k : ℕ)
    (htriv : upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ϖ
      ((k + 1 : ℕ) : ℝ) = (IntermediateField.restrict h₁).fixingSubgroup)
    (htriv' : upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ϖ'
      ((k + 1 : ℕ) : ℝ) = (IntermediateField.restrict h₂).fixingSubgroup)
    (hdeg : Nat.card (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ⧸
        (IntermediateField.restrict h₁).fixingSubgroup)
      = (Nat.card (IsLocalRing.ResidueField (adjoinIntegers K x)) - 1)
          * Nat.card (IsLocalRing.ResidueField (adjoinIntegers K x)) ^ k) :
    E₂ ≤ E₁ := by
  haveI : Fintype ↥((IntermediateField.restrict h₁).fixingSubgroup) := Fintype.ofFinite _
  haveI : Fintype ↥((IntermediateField.restrict h₂).fixingSubgroup) := Fintype.ofFinite _
  refine le_of_fixingSubgroup_restrict_eq_bot h₁ h₂ ?_
  exact subgroup_eq_bot_of_two_quotients K x ht huni habel hϖ hϖ'
    (inf_fixingSubgroup_restrict_eq_bot h₁ h₂ hsup) k htriv htriv' hdeg

end Compositum

/-! ## 残っている穴(★次のノードになるもの。名指しで 2 本)

1. ★★**上付き分岐群の商への降下**。
   `G` が `C` に `G ⧸ H` を経由して作用するとき

       `upperRamificationGroup G ϖ m = Subgroup.comap (QuotientGroup.mk' H)
          (upperRamificationGroup (G ⧸ H) ϖ m)`

   `i_G(σ)` が `σ` の像だけで決まるので `φ` が一致する、というだけの主張で、
   材料は `HerbrandComposition` の `phiOf_quotient` / `exists_quotient_ramIndex` /
   `herbrandPhiGroup_eq_phiOf_quotient` に揃っている。★要るのは
   `MulSemiringAction (G ⧸ H) C` の構成(`H` が核に入るので `QuotientGroup.lift` で作れる)。
   ★これが入ると本ファイルの `htriv` が「商の上付き分岐群が自明」から出せる。

2. ★★**`fixedRing (adjoinIntegers K x) H ≅ adjoinIntegers K x₀` の同一視**
   (`H = Gal(K(x)/K(x₀))`)。作用と `Algebra 𝒪_K` を保つ環同型が要る。
   ★これが入ると B4 の結論 `upperRamificationGroup Gal(K(x₀)/K) α₀ m = ⊥` が
   本ファイルの `htriv` に翻訳でき、同時に `hdeg` も
   `galoisReciprocityEquiv` + `card_units_quotient_span_pi_pow` から出る。
   ★`adjoinField` / `adjoinIntegers` の境界(`lean-idioms.md` #69)に触れるので、
   **不変部分環の側に寄せて書く**のが安全だと思われる。 -/

end ABC3.Found.PGC
