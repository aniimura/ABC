import ABC3.Found.PGC.LubinTateTowerFIndependent

/-!
# `[K̂^m_f : K̂^ur] = q^n − q^{n−1}` —— Yoshida 2008 Proposition 4.4(ii) の次数の部分

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Proposition 4.4(物理 p.7)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-4.html` の `#prop-4-4`。

原文 (Yoshida08 p.7):
> Proposition 4.4. Let m ≥ 1 and f ∈ O[scr]_L[X] as above, with the linear coefficient π. (i) The set µ_f,m is an O[scr]-module by +_F_f and [·]_f. For any α ∈ µ^×_f,m := µ_f,m \ µ_f,m−1, the following is an isomorphism of O[scr]-modules: O[scr]/p[frak]^m ∋ a mod p[frak]^m −→ [a]_f(α) ∈ µ_f,m. (ii) If α ∈ µ^×_f,m, then L^m_f = L(α), N_L^m_f/L(−α) = π^ϕ^m−1 and α is a uniformizer of L^m_f. The L^m_f/L is totally ramified Galois extension of degree |µ^×_f,m| = q^m−1(q − 1). (iii) We have canonical isomorphisms of abelian groups: ρ_f,m : Gal(L^m_f/L) ∼ =−→ Aut_O[scr](µ_f,m) ∼ =−→ (O[scr]/p[frak]^m)^×. (α → [u]_f(α), ∀α ∈ µ_f,m) −→ u mod p[frak]^m

★★本ノードは Proposition 4.4(ii) の**次数の部分のみ**、しかも **`L = K̂`(＝
`L = K̂^{ur}`)の場合**である。原典の残り(`N_{L^m_f/L}(−α) = π^{ϕ^{m−1}}`・
`α` が素元・完全分岐・Galois・(i)・(iii))は本ファイルには**含まれない**。
ただし `L^m_f = L(α)`(原典 (ii) の最初の等式)は次数を出すのに必要なので、
`L = K` の場合(`iteratedLubinTateTorsionPoints_subset_adjoin`)と
`L = K̂^{ur}` の場合(`lubinTateCompletionField_eq_adjoin_simple`)の両方を
本ファイルで示している。

## 原典の証明のうち本ファイルが使う部分(p.7 の Proof、`.txt` の 364–386 行)

> (ii): We have µ_{f,m} ⊂ L(α) by (i), hence L′ = L(α) and L′/L is Galois.
> Now the constant term of f_m/f_{m−1} reads π^{ϕ^{m−1}} = ∏_{α∈µ^×_{f,m}}(−α),
> and taking the v_{L′} of both sides shows e(L′/L) = Σ v_{L′}(−α) ≥ |µ^×_{f,m}| by
> Lemma 4.3(i). But |µ^×_{f,m}| = deg(f_m/f_{m−1}) ≥ [L′ : L] ≥ e(L′/L),
> hence all are equalities and f_m/f_{m−1} is irreducible.

★★**本ファイルは原典より短い道を通る(逸脱ではなく手筋の差)。**
原典は「付値の和 ≥ 根の個数」と「次数 ≥ 拡大次数 ≥ 分岐指数」を挟み撃ちして
既約性を出すが、木にはすでに **`ψ_n = f_n/f_{n−1}` が `𝒪_K` 上 Eisenstein**
(`isEisensteinAt_iteratedLubinTatePsi`、`LubinTateActionPsi.lean`)が在る。
Eisenstein は**係数環を変えても Eisenstein のまま**(抽象核
`isEisensteinAt_map_of_maximalIdeal_eq_span`)なので、`𝒪_{K̂^{ur}}` 上で
そのまま Gauss の補題に載る。★**付値も分岐も一度も出てこない。**

## 何を足したか

### 抽象核(分岐・付値・Lubin-Tate の語彙が 1 つも出てこない)

| 宣言 | 内容 |
|---|---|
| `finrank_adjoin_of_isEisensteinAt` | `A` 整閉整域・`F = Frac A`・`g ∈ A[X]` がモニック Eisenstein なら、`M ⊇ F` の中の根 `α` について `[F(α) : F] = deg g`。Eisenstein の既約性 + Gauss の補題 + `minpoly` |
| `isEisensteinAt_map_of_maximalIdeal_eq_span` | ★★**Eisenstein は係数環の変更で保たれる**: `𝔪_A = (π)` で `φ π` が `B` の非単元・非零なら、`g.map φ` は `(φ π)` で Eisenstein。★局所準同型も DVR も要らない(単元は単元へ写るだけ) |
| `adjoin_image_eq_adjoin_simple` | `S ⊆ F⟮x⟯` かつ `x ∈ S` なら、`F`-代数射 `φ` の像について `F'⟮φ '' S⟯ = F'⟮φ x⟯`(`F'` は `F` 上の任意の底) |

★**`isEisensteinAt_map_of_maximalIdeal_eq_span` が「完備化を渡る段」の実体**である。
在庫調査係は「`K^ur → K̂^ur` の完備化を渡る段が無い」と名指ししたが、実機で測ると
**整数環のレベルではすでに木に在った**(`maximalIdeal_unramifiedCompletionInt_eq_span`、
`UnramifiedCompletionDVR.lean:231`)。足りなかったのは
「その等式を Eisenstein 多項式に載せる」多項式側の 1 本だけである。

### 具体層

| 宣言 | 内容 |
|---|---|
| `algebraMap_comp_baseIntHom_eq` | `𝒪_K → 𝒪_{K̂^{ur}} → K̂^{ur} → ℂ_K` は `𝒪_K → K^al → ℂ_K` に等しい |
| `isEisensteinAt_map_baseIntHom_iteratedLubinTatePsi` | `ψ_n` を `𝒪_{K̂^{ur}}` へ写しても Eisenstein |
| `aeval_closureCompletionCoe_map_iteratedLubinTatePsi` | 原始点 `α` の `ℂ_K` での像は、`K̂^{ur}` 上へ写した `ψ_n` の根 |
| `finrank_adjoin_closureCompletionCoe_iteratedLubinTatePsi` | ★**`[K̂^{ur}(α) : K̂^{ur}] = q^n − q^{n−1}`** |
| `iteratedLubinTateTorsionPoints_subset_adjoin` | ★★**`Λ_n ⊆ K(α)`**(原典 (ii) の `µ_{f,m} ⊂ L(α)`、`L = K` 版)。`n` の帰納法 |
| `lubinTateCompletionField_eq_adjoin_simple` | ★**`K̂^m_f = K̂^{ur}(α)`**(原典 (ii) の `L^m_f = L(α)`) |
| `finrank_lubinTateCompletionField` | ★★★**`[K̂^m_f : K̂^{ur}] = q^n − q^{n−1}`** |

## 退化の自己検査

* ★★**`α ∈ µ^×_{f,m}`(位数がちょうど `π^n`)を落とすと偽**——低い層の点 `α ∈ Λ_{n−1}`
  では `K̂^{ur}(α)` の次数は下がる。本ファイルはどの宣言でも
  `hxψ : x ∈ iteratedLubinTatePsiTorsionPoints … n hn`(＝`ψ_n` の根)を要求する。
* ★**`m ≥ 1` を落とすと `q^{m−1}(q−1)` が意味を失う**——`hn : 1 ≤ n` は必須。
  `n = 0` では `Λ_0 = {0}` で `ψ_0` が定義されない。
* ★★**`𝒪_{K̂^{ur}}` が DVR であることを落とすと Gauss の補題が使えない**——
  抽象核 `finrank_adjoin_of_isEisensteinAt` は `[IsIntegrallyClosed A]` を要求し、
  具体層はそれを `isDiscreteValuationRing_unramifiedCompletionInt` から得ている。
  これは `K̂^{ur}` の**完備性**に依る事実である。
* ★★**`𝔪_{𝒪_{K̂^{ur}}} = (π)`(＝ `K̂^{ur}/K` が不分岐)を落とすと偽**——
  もし `π` が `𝒪_{K̂^{ur}}` で `ϖ^e`(`e ≥ 2`)になれば `ψ_n` の定数項が `𝔪²` に落ちて
  Eisenstein でなくなり、次数は下がりうる。本ファイルは
  `maximalIdeal_unramifiedCompletionInt_eq_span` でこれを供給している。
* ★★**剰余体の有限性は `𝒪_K` 側にしか使っていない。** `𝒪_{K̂^{ur}}` の剰余体は代数閉
  (＝無限)なので `[Fintype (ResidueField 𝒪_{K̂^{ur}})]` を要求する在庫は使えないが、
  本ファイルはそれを 1 つも使っていない(決定 D28 の予測どおり)。
  `[Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]` は `𝒪_K` の剰余体についてであり、
  こちらは有限で正しい。

## 逸脱(記録)

1. **`L = K̂^{ur}` に固定した。** 原典 §4.1 は「a complete unramified extension `L` of `K`
   を固定する」と一般に書くが、本ファイルは `L = K̂^{ur}`(最大不分岐拡大の完備化)の
   場合だけを扱う。Proposition 4.7(ii) 以降で実際に使われるのはこの場合である。
2. **`K̂^m_f` の定義は #9(`LubinTateTowerFIndependent.lean`)のものを使う。**
   すなわち `ℂ_K` の中の `IntermediateField.adjoin (unramifiedCompletion K) (ι_al '' Λ_{f,n})`。
   原典の `L^m_f := L(µ_{f,m})` との対応は #9 の docstring に記録済み。
3. **次数の書き方**。原典は `q^{m−1}(q − 1)` と書くが、木の在庫
   (`finrank_adjoin_iteratedLubinTatePsi`・`natDegree_iteratedLubinTatePsi`)に合わせて
   `q^n − q^{n−1}` と書く(`ℕ` の引き算だが `q^{n−1} ≤ q^n` なので切り詰めは起きない)。
4. **原典の証明の手筋を使っていない**(上の「原典より短い道」を参照)。結論は同じ。

## `L = K` 版は流用できたか

★**次数そのもの(`finrank_adjoin_iteratedLubinTatePsi`、`LubinTateDegree.lean:31`)は
流用できなかった。** あちらは `Gal(K(α)/K) ≃ (𝒪_K/π^n)^×`(Lubin-Tate の主定理)から
次数を読む道で、底を `K̂^{ur}` に取り替えると Galois 群の同定からやり直しになる。
本ファイルは**多項式の既約性から直接**次数を読む(Eisenstein → Gauss → `minpoly`)。
★一方 `Λ_n ⊆ K(α)` の側は、木の `iteratedLubinTatePsiTorsionPoints_subset_adjoin`
(同じ層の原始点)と `lubinTateActionAtTorsionPoint_pi_mem_iteratedLubinTatePsiTorsionPoints`
(`[π]α` は 1 つ下の層の原始点)を帰納法で繋ぐだけで済んだ。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal Valued Classical

/-! ## 0. 抽象核

分岐・付値・Lubin-Tate の語彙が 1 つも出てこない部分。原典の設定に依らない。 -/

section AbstractCore

open Polynomial in
/-- ★★★★★★★★★★★★★★★★**Eisenstein 多項式の根が生成する拡大の次数は
多項式の次数**。

`A` を整閉整域、`F = Frac A`、`M ⊇ F` を体とする。`g ∈ A[X]` がモニックで
素イデアル `𝓟` について Eisenstein なら、`M` の中の `g` の根 `α` について
`[F(α) : F] = deg g`。

段取り: Eisenstein の判定法(`Polynomial.IsEisensteinAt.irreducible`)で `A[X]` の中の
既約性を出し、Gauss の補題(`Monic.irreducible_iff_irreducible_map_fraction_map`)で
`F[X]` へ移し、`minpoly.eq_of_irreducible_of_monic` で最小多項式と同定して
`IntermediateField.adjoin.finrank` を読むだけ。

★`Algebra A M` も `IsScalarTower A F M` も**要らない**——根の条件を
`aeval α (g.map (algebraMap A F)) = 0` の形で受け取れば `F` から下は見えない。 -/
theorem finrank_adjoin_of_isEisensteinAt
    {A F M : Type*} [CommRing A] [IsDomain A] [IsIntegrallyClosed A]
    [Field F] [Algebra A F] [IsFractionRing A F] [Field M] [Algebra F M]
    {𝓟 : Ideal A} (h𝓟 : 𝓟.IsPrime)
    {g : Polynomial A} (hmonic : g.Monic) (heis : g.IsEisensteinAt 𝓟)
    (hdegpos : 0 < g.natDegree)
    (α : M) (hroot : Polynomial.aeval α (g.map (algebraMap A F)) = 0) :
    Module.finrank F ↥(IntermediateField.adjoin F ({α} : Set M)) = g.natDegree := by
  have hirrA : Irreducible g := heis.irreducible h𝓟 hmonic.isPrimitive hdegpos
  have hmonicF : (g.map (algebraMap A F)).Monic := hmonic.map _
  have hirrF : Irreducible (g.map (algebraMap A F)) :=
    (hmonic.irreducible_iff_irreducible_map_fraction_map (K := F)).mp hirrA
  have hmin : g.map (algebraMap A F) = minpoly F α :=
    minpoly.eq_of_irreducible_of_monic hirrF hroot hmonicF
  have hint : IsIntegral F α := ⟨_, hmonicF, hroot⟩
  rw [IntermediateField.adjoin.finrank hint, ← hmin, hmonic.natDegree_map]

open Polynomial in
/-- ★★★★★★★★★★★★★★★★★★★★**Eisenstein は係数環の変更で保たれる**。

`A` は局所整域で `𝔪_A = (π)`、`φ : A →+* B` は任意の環準同型、`B` は整域。
`φ π` が `B` で**非零**かつ**非単元**なら、`A` 上 `𝔪_A`-Eisenstein なモニック多項式 `g`
の像 `g.map φ` は `(φ π)` について Eisenstein。

★段取りは 3 行に尽きる: `g.coeff 0 = u·π`(`u` は `A` の単元、`𝔪_A²` に入らないから)
⇒ `φ (g.coeff 0) = φ u · φ π` で `φ u` は**単元**(単元は単元へ写る)
⇒ `φ π · φ u ∈ ((φ π)²)` は `φ π` を約せて `φ u ∈ (φ π)`、すなわち `φ π` が単元となり矛盾。

★★**`φ` が局所準同型であることも `B` が DVR であることも要らない。**
必要なのは「`φ π` が `B` の非単元」だけであり、これは
`𝔪_B = (φ π)`(不分岐)から自動的に出る。 -/
theorem isEisensteinAt_map_of_maximalIdeal_eq_span
    {A B : Type*} [CommRing A] [IsLocalRing A] [IsDomain A] [CommRing B] [IsDomain B]
    (φ : A →+* B) {π : A}
    (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (hϖ0 : φ π ≠ 0) (hϖu : ¬ IsUnit (φ π))
    {g : Polynomial A} (hmonic : g.Monic)
    (heis : g.IsEisensteinAt (IsLocalRing.maximalIdeal A)) (hdegpos : 0 < g.natDegree) :
    (g.map φ).IsEisensteinAt (Ideal.span {φ π}) := by
  have hntop : (Ideal.span {φ π}) ≠ ⊤ := by
    simpa [Ideal.span_singleton_eq_top] using hϖu
  have hdegmap : (g.map φ).natDegree = g.natDegree := hmonic.natDegree_map φ
  refine ⟨?_, ?_, ?_⟩
  · show (g.map φ).leadingCoeff ∉ _
    rw [(hmonic.map φ).leadingCoeff]
    intro h
    exact hntop (Ideal.eq_top_of_isUnit_mem _ h isUnit_one)
  · intro k hk
    rw [hdegmap] at hk
    rw [Polynomial.coeff_map]
    have hmem := heis.mem hk
    rw [hπmax, Ideal.mem_span_singleton'] at hmem
    obtain ⟨c, hc⟩ := hmem
    rw [← hc, map_mul, Ideal.mem_span_singleton']
    exact ⟨φ c, rfl⟩
  · rw [Polynomial.coeff_map, Ideal.span_singleton_pow]
    have h0 : g.coeff 0 ∈ IsLocalRing.maximalIdeal A := heis.mem hdegpos
    rw [hπmax, Ideal.mem_span_singleton'] at h0
    obtain ⟨u, hu⟩ := h0
    have hune : u ∉ IsLocalRing.maximalIdeal A := by
      intro humem
      apply heis.notMem
      rw [hπmax, Ideal.span_singleton_pow, ← hu]
      rw [hπmax, Ideal.mem_span_singleton'] at humem
      obtain ⟨d, hd⟩ := humem
      exact Ideal.mem_span_singleton'.mpr ⟨d, by rw [← hd]; ring⟩
    have huunit : IsUnit u := by
      rwa [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff, not_not] at hune
    intro hmem
    rw [← hu, map_mul, Ideal.mem_span_singleton'] at hmem
    obtain ⟨w, hw⟩ := hmem
    have hkey : φ u = φ π * w := by
      apply mul_left_cancel₀ hϖ0
      linear_combination -hw
    exact hϖu (isUnit_of_mul_isUnit_left (by rw [← hkey]; exact huunit.map φ))

/-- ★★★★★★★★★★★★**単一生成の像への還元**。

`φ : L →ₐ[F] M`、`S ⊆ F⟮x⟯`、`x ∈ S` とすると、`F` 上の任意の底 `F'`(`F ⊆ F' ⊆ M`)に
ついて `F'⟮φ '' S⟯ = F'⟮φ x⟯`。

★原典 (ii) の `L^m_f = L(α)` を `L = K` から `L = K̂^{ur}` へ持ち上げるのがこれ。
「`K` 上で 1 元生成なら、底を大きくしても同じ 1 元で生成される」という当たり前の事実だが、
中間体の `map` と `restrictScalars` を挟むので 1 本にしておく。 -/
theorem adjoin_image_eq_adjoin_simple
    {F L M F' : Type*} [Field F] [Field L] [Field M] [Field F']
    [Algebra F L] [Algebra F M] [Algebra F' M] [Algebra F F'] [IsScalarTower F F' M]
    (φ : L →ₐ[F] M) {x : L} {S : Set L} (hx : x ∈ S)
    (hS : S ⊆ (IntermediateField.adjoin F ({x} : Set L) : IntermediateField F L)) :
    IntermediateField.adjoin F' (φ '' S) = IntermediateField.adjoin F' ({φ x} : Set M) := by
  refine le_antisymm (IntermediateField.adjoin_le_iff.mpr ?_) ?_
  · rintro _ ⟨y, hy, rfl⟩
    have h1 : φ y ∈ (IntermediateField.adjoin F ({x} : Set L)).map φ := ⟨y, hS hy, rfl⟩
    rw [IntermediateField.adjoin_map, Set.image_singleton] at h1
    have h2 : IntermediateField.adjoin F ({φ x} : Set M)
        ≤ (IntermediateField.adjoin F' ({φ x} : Set M)).restrictScalars F :=
      IntermediateField.adjoin_le_iff.mpr (by
        rintro _ rfl
        exact IntermediateField.subset_adjoin F' _ rfl)
    exact h2 h1
  · exact IntermediateField.adjoin_simple_le_iff.mpr
      (IntermediateField.subset_adjoin F' _ ⟨x, hx, rfl⟩)

end AbstractCore

variable {p : ℕ} [Fact p.Prime]

/-! ## 1. 二つの経路が一致する —— `𝒪_K → ℂ_K`

`ψ_n` は `𝒪_K` 係数の多項式である。それを `ℂ_K` の中で評価する道は 2 つある:

* `𝒪_K → 𝒪_{K̂^{ur}} → K̂^{ur} → ℂ_K`(本ファイルが Eisenstein を載せる道)
* `𝒪_K → K^al → ℂ_K`(木が捩れ点を作っている道)

この 2 つが同じ環準同型であることを先に確かめておく。 -/

/-- `𝒪_K → 𝒪_{K̂^{ur}} → K̂^{ur} → ℂ_K` は `𝒪_K → K^al → ℂ_K` に等しい。 -/
theorem algebraMap_comp_baseIntHom_eq (K : PAdicLocalField p) :
    ((algebraMap (unramifiedCompletion K) (closureCompletion K)).comp
        (algebraMap ↥(unramifiedCompletionInt K) (unramifiedCompletion K))).comp (baseIntHom K)
      = (closureCompletionCoe K).comp (algebraMap 𝒪[K.carrier] K.closure) := by
  ext a
  simp only [RingHom.coe_comp, Function.comp_apply]
  show algebraMap (unramifiedCompletion K) (closureCompletion K)
      ((baseIntHom K a : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K) = _
  rw [baseIntHom_coe, algebraMap_unramifiedCompletion_closureCompletion_eq,
    unramifiedToClosureCompletion_algebraMap,
    IsScalarTower.algebraMap_apply 𝒪[K.carrier] K.carrier K.closure,
    UniformSpace.Completion.algebraMap_def]
  rfl

/-! ## 2. `ψ_n` を `𝒪_{K̂^{ur}}` へ渡す —— 「完備化を渡る段」 -/

/-- `deg ψ_n = q^n − q^{n−1} > 0`(`q > 1` と `n ≥ 1` から)。 -/
theorem sub_pow_pred_pos (K : PAdicLocalField p)
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {pp ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    (n : ℕ) (hn : 1 ≤ n) : 0 < (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  have h2 : 1 < pp ^ ff := hq ▸ Fintype.one_lt_card
  have hlt : (pp ^ ff) ^ (n - 1) < (pp ^ ff) ^ n := Nat.pow_lt_pow_right h2 (by omega)
  omega

/-- ★★★★★★★★★★★★★★★★**`ψ_n` を `𝒪_{K̂^{ur}}` へ写しても Eisenstein**。

これが「`K^{ur} → K̂^{ur}` の完備化を渡る段」の実体である。渡すのに要るのは
**`𝔪_{𝒪_{K̂^{ur}}} = (π)`**(`maximalIdeal_unramifiedCompletionInt_eq_span`、木に在る)
だけで、あとは抽象核 `isEisensteinAt_map_of_maximalIdeal_eq_span` に代入する。

★★`𝒪_{K̂^{ur}}` の剰余体は代数閉(＝無限)なので、`ψ_n` を「`A = 𝒪_{K̂^{ur}}` で
もう一度作り直す」道(`[Fintype (ResidueField A)]` を要求する Lubin-Tate 在庫)は
**論理的に取れない**(決定 D28)。本ファイルは `ψ_n` を `𝒪_K` 上で作ったまま**写す**。 -/
theorem isEisensteinAt_map_baseIntHom_iteratedLubinTatePsi (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map (baseIntHom K)).IsEisensteinAt
      (IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K)) := by
  have hπ0 : (π : K.carrier) ≠ 0 := by simpa using hπne0
  have hspan : IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K)
      = Ideal.span {baseIntHom K π} := by
    rw [baseIntHom_eq_uniformizerCompletionInt]
    exact maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax
  rw [hspan]
  refine isEisensteinAt_map_of_maximalIdeal_eq_span (baseIntHom K) hπmax ?_ ?_
    (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
    (isEisensteinAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn) ?_
  · rw [baseIntHom_eq_uniformizerCompletionInt]
    exact uniformizerCompletionInt_ne_zero K hπ0
  · intro hu
    exact (IsLocalRing.maximalIdeal.isMaximal ↥(unramifiedCompletionInt K)).ne_top
      (by rw [hspan]; exact Ideal.span_singleton_eq_top.mpr hu)
  · rw [natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]
    exact sub_pow_pred_pos K hq n hn

/-- 原始点 `α ∈ Λ^×_{f,n}` の `ℂ_K` での像は、`K̂^{ur}` へ写した `ψ_n` の根。 -/
theorem aeval_closureCompletionCoe_map_iteratedLubinTatePsi (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hx : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    Polynomial.aeval (closureCompletionCoe K x)
        (((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map (baseIntHom K)).map
          (algebraMap ↥(unramifiedCompletionInt K) (unramifiedCompletion K))) = 0 := by
  have hroot0 : Polynomial.eval x
      ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map
        (algebraMap 𝒪[K.carrier] K.closure)) = 0 :=
    (Polynomial.mem_roots'.mp (Multiset.mem_toFinset.mp hx)).2
  rw [Polynomial.aeval_def, Polynomial.eval₂_map, Polynomial.eval₂_map,
    algebraMap_comp_baseIntHom_eq, ← Polynomial.hom_eval₂, ← Polynomial.eval_map, hroot0, map_zero]

/-- ★★★★★★★★★★★★★★★★★★★★**`[K̂^{ur}(α) : K̂^{ur}] = q^n − q^{n−1}`**
(`α ∈ Λ^×_{f,n}`)。

原典 Proposition 4.4(ii) の「degree `|µ^×_{f,m}| = q^{m−1}(q − 1)`」の部分を、
`L = K̂^{ur}` の場合に、`L^m_f = L(α)` を経由せずに**単純拡大について**述べたもの。
`Λ_n` 全体を添加した `K̂^m_f` についての形は `finrank_lubinTateCompletionField`。 -/
theorem finrank_adjoin_closureCompletionCoe_iteratedLubinTatePsi (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hx : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    Module.finrank (unramifiedCompletion K)
        ↥(IntermediateField.adjoin (unramifiedCompletion K)
          ({closureCompletionCoe K x} : Set (closureCompletion K)))
      = (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  haveI := isDiscreteValuationRing_unramifiedCompletionInt K
  have hmonic : ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map
      (baseIntHom K)).Monic :=
    (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic.map _
  have hdeg : ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map
      (baseIntHom K)).natDegree = (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
    rw [(isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n
        hn).monic.natDegree_map,
      natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]
  rw [← hdeg]
  refine finrank_adjoin_of_isEisensteinAt
    (IsLocalRing.maximalIdeal.isMaximal ↥(unramifiedCompletionInt K)).isPrime
    hmonic (isEisensteinAt_map_baseIntHom_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn)
    ?_ _
    (aeval_closureCompletionCoe_map_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x hx)
  rw [hdeg]
  exact sub_pow_pred_pos K hq n hn

/-! ## 3. `L^m_f = L(α)` —— まず `L = K`、つぎに `L = K̂^{ur}` -/

/-- ★★★★★★★★★★★★★★★★★★**`Λ_n ⊆ K(α)`**(`α ∈ Λ^×_{f,n}`)——
原典 Proposition 4.4(ii) の `µ_{f,m} ⊂ L(α)` の `L = K` 版。

原典は (i)(＝`𝒪/𝔭^m ≅ µ_{f,m}`、`a ↦ [a]_f(α)`)から一息に出すが、本ファイルは
`n` の帰納法で回る:

* `Λ_n = Λ_{n−1} ∪ Λ^×_n`(`iteratedLubinTateTorsionPoints_eq_union`)
* `Λ^×_n ⊆ K(α)`(`iteratedLubinTatePsiTorsionPoints_subset_adjoin`、単数の作用の全射性)
* `[π]α` は 1 つ下の層の原始点で、しかも `K(α)` の元
  (`lubinTateActionAtTorsionPoint_pi_mem_iteratedLubinTatePsiTorsionPoints`。
  ★`lubinTateActionAtTorsionPoint` の値がそもそも `adjoinIntegers K α ⊆ K(α)` にある)

★★これで (i) 全体(`𝒪`-加群としての同型)を作らずに済んだ。
★`n = 1` の底は `Λ_0 = {0}`(`iteratedLubinTateTorsionPoints_zero`)。 -/
theorem iteratedLubinTateTorsionPoints_subset_adjoin (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff)) :
    ∀ (n : ℕ) (hn : 1 ≤ n) (x : K.closure),
      x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn →
      x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n →
      ∀ y ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n,
        y ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
  intro n
  induction n with
  | zero => intro hn; omega
  | succ m ih =>
    intro _ x hxψ hxn y hy
    haveI := finiteDimensional_adjoin_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0
      f hf0 hf1 hf (m + 1) x hxn
    have hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure) :=
      IntermediateField.mem_adjoin_simple_self _ _
    rw [iteratedLubinTateTorsionPoints_eq_union K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega),
      Finset.mem_union] at hy
    rcases hy with hy | hy
    · rcases Nat.eq_zero_or_pos m with rfl | hm
      · simp only [Nat.add_sub_cancel] at hy
        rw [iteratedLubinTateTorsionPoints_zero K hq hπmax hπne0 f hf0 hf1 hf,
          Finset.mem_singleton] at hy
        subst hy
        exact zero_mem _
      · simp only [Nat.add_sub_cancel] at hy
        set w := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1) x hxn hmem π
          with hw_def
        have hzψ : ((w : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)
            ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf m hm :=
          lubinTateActionAtTorsionPoint_pi_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0
            f hf0 hf1 hf m hm x hxψ hxn hmem
        have hzn : ((w : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)
            ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf m :=
          lubinTateActionAtTorsionPoint_pi_mem_pred K hq hπmax hπne0 f hf0 hf1 hf m x hxn hmem
        have hle : IntermediateField.adjoin K.carrier
            ({((w : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)}
              : Set K.closure)
            ≤ IntermediateField.adjoin K.carrier ({x} : Set K.closure) :=
          IntermediateField.adjoin_simple_le_iff.mpr
            (w : IntermediateField.adjoin K.carrier ({x} : Set K.closure)).2
        exact hle (ih hm _ hzψ hzn y hy)
    · exact iteratedLubinTatePsiTorsionPoints_subset_adjoin K hq hπmax hπne0 f hf0 hf1 hf
        (m + 1) (by omega) x hxψ hxn hmem y hy

/-- ★★★★★★★★★★★★★★★★★★★★**`K̂^m_f = K̂^{ur}(α)`**(`α ∈ Λ^×_{f,n}`)——
原典 Proposition 4.4(ii) の `L^m_f = L(α)`、`L = K̂^{ur}` 版。

`Λ_n ⊆ K(α)`(上)を `ι_al : K^al →ₐ[K] ℂ_K` で送り、抽象核
`adjoin_image_eq_adjoin_simple` に渡すだけ。 -/
theorem lubinTateCompletionField_eq_adjoin_simple (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n
      = IntermediateField.adjoin (unramifiedCompletion K)
          ({closureCompletionCoe K x} : Set (closureCompletion K)) :=
  adjoin_image_eq_adjoin_simple (F' := unramifiedCompletion K) (closureCompletionAlgHom K)
    (S := ((iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n : Finset K.closure) :
      Set K.closure))
    hxn
    (fun y hy => iteratedLubinTateTorsionPoints_subset_adjoin K hq hπmax hπne0 f hf0 hf1 hf
      n hn x hxψ hxn y hy)

/-! ## 4. 主結果 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`[K̂^m_f : K̂^{ur}] = q^n − q^{n−1}`**(Yoshida 2008 Proposition 4.4(ii)、`L = K̂^{ur}`)。

`K̂^m_f = K̂^{ur}(α)`(`lubinTateCompletionField_eq_adjoin_simple`)と
`[K̂^{ur}(α) : K̂^{ur}] = q^n − q^{n−1}`
(`finrank_adjoin_closureCompletionCoe_iteratedLubinTatePsi`)を繋ぐだけ。
原始点 `α` の存在は `iteratedLubinTatePsiTorsionPoints_nonempty`(木に在る)から。 -/
theorem finrank_lubinTateCompletionField (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    Module.finrank (unramifiedCompletion K)
        ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n)
      = (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  obtain ⟨x, hxψ⟩ :=
    iteratedLubinTatePsiTorsionPoints_nonempty K hq hπmax hπne0 f hf0 hf1 hf n hn
  have hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n := by
    rw [iteratedLubinTateTorsionPoints_eq_union K hq hπmax hπne0 f hf0 hf1 hf n hn,
      Finset.mem_union]
    exact Or.inr hxψ
  rw [lubinTateCompletionField_eq_adjoin_simple K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn]
  exact finrank_adjoin_closureCompletionCoe_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf
    n hn x hxψ

def finrank_lubinTateCompletionField.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 7, item := "Proposition 4.4", sectionId := "prop-4-4" }

end ABC3.Found.PGC
