import ABC3.Found.PGC.RamificationFiltrationZero
import ABC3.Found.PGC.LubinTateSharpRamificationImage

/-!
# `Art(Γ_K^n) = U^n_K` —— ★★段データの値を Lubin-Tate 塔で計算する

`Found/PGC/RamificationImageUnits.lean`(Y21)は pGC p.4 の

> Then it is well-known (Theorem 1 of [3], p. 155) that the image of Γ^v_K in Γ^ab_K is equal to U^v_K ⊆ U_K.

を、段データ `F : StageFiltration K.absGal` を**仮定として受け取り**、`F` に課される
条件を 2 本(`hbase` / `hstage`)に絞った形で証明した。
`Found/PGC/UnramifiedBaseChangeInvariance.lean`(Y19e)は `compat` を無条件で埋めて
**本物の `ramificationFiltration p`** を作った。

★★**本ファイルはその 2 本を落とす。** 到達点は

```
map_ramificationFiltration_reciprocityUnits_eq_principalUnits :
    Subgroup.map (reciprocityUnits K ...) ((ramificationFiltration p).Gv K (n : ℝ))
      = principalUnits K π n
```

であり、**仮定は 1 つも無い**(段データは Y19e の本物、退化 witness ではない)。

## ★★★配られた持ち場との差(隠さない)

持ち場は「`StageFiltration K.absGal` の `compat` を埋めよ」だった。
★**着手時に測ったところ `compat` は既に埋まっていた**(`UnramifiedBaseChangeInvariance.lean`
の `absGalStage_compat` / `unramifiedBaseChange`、commit `916b7082`)。
持ち場が挙げていた 2 つの残り((1) 固定環と中間体の整数環の同定、(2) 不分岐底変換で
上付き番号付けが変わらないこと)も、そこで `stageBaseHom` / `ramIndex_ringHom_eq` として
入っている。★**本ファイルは持ち場の「ゴール」(消費先の仮定 2 本が落ちること)の方を埋めた。**

## 何が入ったか

### §1 抽象核 —— ★分岐・付値・Galois の語彙が出てこない(2 本は 1 語も出ない)

* `coe_mul_coe_eq_self` / `eq_of_coe_mul_coe_eq_coe` —— **純群論**。
  `N ≤ S` なら `S · N = S`。したがって `S · N = T` から `S = T` が出る。
  ★★**これが「段データの選択独立性」の全部**である(下記 §2)。
* `ramIndex_eq_of_adjoin_eq_top` —— ★**`i(σ) = v(σα − α)` は
  「`𝒪_L` を 1 元生成する素元」の取り方に依らない**。`A` と `A'` は**別の底環でよい**。
  証明は `σ ∈ G_n ⟺ n < i(σ)` を両方の素元で書いて `ℕ∞` の外延性を使うだけ(3 行)。
* `upperRamificationGroup_congr_of_ramIndex` / `upperRamificationGroup_congr_adjoin`
  —— したがって `φ_G`・`ψ_G`・`G^m` もすべて素元の取り方に依らない。
  ★Y19e §2 の `upperRamificationGroup_eq_comap` に `f = id` を入れるだけ。

### §2 ★★★段データの選択独立性(Y19d 逸脱 1・2 の解消)

* `stage_eq_stage` —— ★★**同じ `N` の 2 つの `StageGenerator` は同じ段を与える**。
  ★**`compat` を `M = N` で使うだけ**である: `stage_mul_coe_eq h K g g' le_rfl v` は
  `(g.stage v) · N = g'.stage v` を与え、左辺は `N ≤ g.stage v`(`le_stage`)から
  `g.stage v` に潰れる。★Y19d の逸脱 2(「異なる `x` が同じ `S N v` を与えることは
  証明していない」)と Y19e の「証明していないこと」がこれで消える。
* `absGalStage_eq_stage_any` —— したがって `absGalStage K N v` は**どの生成元でも計算できる**。
  ★★**これが本ファイルを可能にした鍵**である。これが無いと `absGalStage` の中の
  `Classical.choice` が Lubin-Tate の生成元と噛み合わない。
* `adjoin_uniformizer_inertiaFixedRing_eq_top` /
  `stageUpperRamification_eq_map_upperRamificationGroup` —— ★段の一意化元も任意でよい
  (Y19d 逸脱 1 の解消)。

### §3 完全分岐段では `Gal(L/L₀)^v` が `Gal(L/K)^v` になる

* `inertiaGalAdjoin_eq_top_of_isTotallyRamified` —— `K(x)/K` 完全分岐 ⇒ `I(L/K) = ⊤`。
* `stageUpperRamification_eq_upperRamificationGroup` —— ★★★したがって段の上付き分岐群は
  **`Gal(K(x)/K)` 全体の上付き分岐群**であり、しかも `𝒪_{K(x)}` の**どの素元**で
  計算してもよい。★決め手は §1 の `ramIndex_eq_of_adjoin_eq_top` で、
  底環を `𝒪_L^{I}`(段側)と `𝒪_K`(Lubin-Tate 側)で**取り替える**ところである。

### §4 Lubin-Tate 塔での計算

* `absGalPrincipalLevel_zero`(`Art^{-1}(U^0) = Γ_K`)、
  `isGalois_adjoin_psiGenSeq` / `normal_adjoin_psiGenSeq`
  —— ★**塔の第 `m+1` 段が `K` 上正規**であることを、
  `absGalPrincipalLevel` が正規部分群であること(Y21)から `InfiniteGalois.normal_iff_isGalois`
  で降ろした。★**新しい数学は要らなかった。**
* `psiStageGenerator` —— `absGalPrincipalLevel (m+1)` の `StageGenerator`。
* `algEquivRestrictSelfHom_eq_restrictNormalHom` —— Y22 の `algEquivRestrictSelfHom` は
  `AlgEquiv.restrictNormalHom` そのもの(`coe_algEquivRestrictSelf` は `rfl` なので 2 行)。
* ★★★★`absGalStage_absGalPrincipalLevel` ——

      `absGalStage K (Art^{-1}(U^m)) (I : ℝ) = Art^{-1}(U^I)`   (`I ≤ m`)

  ★★**これが消費先の `hstage` そのもの**である。

### §5 ★★★★★消費先の仮定 2 本が落ちる

* `map_ramificationFiltration_reciprocityUnits_eq_principalUnits` /
  `..._eq_upperPrincipalUnits` —— **`Art(Γ_K^n) = U^n_K`(仮定ゼロ)**。
* `ramificationFiltration_Gv_le_absGalPrincipalLevel` /
  `ramificationFiltration_Gv_le_fixingSubgroup` —— `Γ_K^n` は
  Lubin-Tate 塔の第 `n` 段 `K_{f,n}` を各点固定する。
* ★★`map_absInertia_reciprocityUnits_eq_top` —— **退化の自己検査**。
  `n = 0` を入れると `Art(I_K) = 𝒪_K^×`(`Γ_K^0 = I_K` は Y19f)。
  ★古典的に正しい主張が出るので、計算が退化していない。

## ★何が入っていないか(正直な線引き)

1. ★★**`v` は自然数である。** 実数 `v` について `Γ_K^v = Γ_K^{⌈v⌉}` を言うには
   跳びの整数性(Hasse-Arf)が要る。★Y21 の「逸脱の記録 3」をそのまま引き継ぐ。
   ★**別ノード**。
2. ★★**`Γ^ab_K` ではなく `𝒪_K^×` 成分で述べている。** Y21 の「逸脱の記録 1」を引き継ぐ。
   原典の `Γ^v_K → Γ^ab_K` の像を言うには「`v ≥ 0` で `Γ_K^v ⊆ I_K`」が要り、
   それは `Found/PGC/RamificationFiltrationZero.lean` の
   `ramificationFiltration_Gv_le_absInertia` に在る。★**両者を繋ぐノードは書いていない。**
3. ★**`Γ_K^n = Art^{-1}(U^n_K)` は成り立たない**(右辺は `ker(Art)` を丸ごと含む)。
   本ファイルが言うのは包含 `⊆`(`ramificationFiltration_Gv_le_absGalPrincipalLevel`)と
   **像の等号**だけである。

## ★原典より短い道(名指し)

1. ★★★**選択独立性は `compat` の系である。** 原典(および普通の教科書)は上付き分岐群の
   well-defined 性を「生成元・素元に依らない」ことから示すが、本ファイルは逆向きに、
   **`compat` を `M = N` に当てて**それを出した(`stage_eq_stage`、証明 2 行)。
   ★Y19e が「生成元の独立性は要らなかった」と書いた場所が、実は
   **独立性そのものを与えていた**。
2. ★**素元の独立性に微分も Herbrand も要らない。** `ramIndex_eq_of_adjoin_eq_top` は
   `σ ∈ G_n ⟺ n < i(σ)` を 2 通りに書くだけで、しかも**底環を取り替えてよい**。
   これが無いと「段側の底環 `𝒪_L^{I}`」と「Lubin-Tate 側の底環 `𝒪_K`」が繋がらない。
3. ★**塔の正規性は群の正規性から降ろせる。** `Gal(K̄/K_{f,m})` が正規(Y21、
   可換群への準同型の引き戻しだから)⇒ `K_{f,m}/K` が Galois、という向きに使った。
   Lubin-Tate 理論の「`K_{f,m}/K` はアーベル」を再証明していない。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★**`v` は自然数**(上記「入っていないもの」1。Y21 逸脱 3 の引き継ぎ)。
2. ★**`Γ^ab_K` ではなく `𝒪_K^×`**(Y21 逸脱 1 の引き継ぎ)。
3. ★**Lubin-Tate の底 `L = K`**(Y22 逸脱 1・決定 D29 の引き継ぎ)。
   原典 Yoshida Prop 6.14 の `L = K_n` は本ファイルでも `K` である。
4. ★**既存の `Found/PGC/*.lean` は 1 行も書き換えていない**(import のみ。D28)。
   `Found/PGC/LubinTate*.lean` は読んだだけである。
5. ★本ファイルは原典が名前を付けていない組み立てノードなので `.src` を持たない。
   使っている原典項目は pGC p.4 の "well-known" の段落と Yoshida Prop 6.14 であり、
   `.src` はそれぞれ `RamificationImageUnits.lean` / `LubinTateSharpRamificationImage.lean`
   にある。

## 配管の記録(`lean-idioms.md` 行き)

* ★**`haveI : Fintype ↥(inertiaGalAdjoin K x) := Fintype.ofFinite _` を書くと
  `rw` が当たらなくなる**。木には既に `fintypeInertiaGalAdjoin` が instance で在るので、
  ★**書かないほうが通る**。
  エラーは `Tactic `rewrite` failed: Did not find an occurrence of the pattern
  upperRamificationGroup (↥(inertiaGalAdjoin K x)) (stageUniformizer K x) v` である。
* ★`IsGalois` から `Normal` を取り出すのに `.toNormal` は使えない
  (``Invalid field `toNormal`: The environment does not contain `IsGalois.toNormal```)。
  `haveI := hgal; infer_instance` と書く。
* ★中間体の 2 層の塔は 1 つも作っていない(#59 回避)。
-/

namespace ABC3.Found.PGC

open IsLocalRing ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped Pointwise NNReal Valued

/-! ## §1 抽象核

★`coe_mul_coe_eq_self` / `eq_of_coe_mul_coe_eq_coe` は**純群論**
(分岐・付値・Galois の語彙が 1 語も出てこない)。 -/

/-- ★★**抽象核(純群論)** —— `N ≤ S` なら `S · N = S`。 -/
theorem coe_mul_coe_eq_self {Γ : Type*} [Group Γ] {S N : Subgroup Γ} (h : N ≤ S) :
    ((S : Set Γ) * (N : Set Γ)) = (S : Set Γ) :=
  Set.Subset.antisymm (by rintro _ ⟨a, ha, b, hb, rfl⟩; exact mul_mem ha (h hb))
    (fun x hx => ⟨x, hx, 1, one_mem _, mul_one x⟩)

/-- ★★★**抽象核(純群論)** —— `N ≤ S` かつ `S · N = T` なら `S = T`。

★★**これが「段データが選択に依らない」ことの全部**である:
`StageFiltration.compat` を `M = N` に当てると `S N v · N = S N v`(別の選択で作った側)
になり、左辺は本補題で潰れる。 -/
theorem eq_of_coe_mul_coe_eq_coe {Γ : Type*} [Group Γ] {S T N : Subgroup Γ} (hS : N ≤ S)
    (h : ((S : Set Γ) * (N : Set Γ)) = (T : Set Γ)) : S = T :=
  SetLike.coe_injective (by rw [← h, coe_mul_coe_eq_self hS])

/-- ★★★**抽象核** —— `i(σ) := v(σα − α)` は「`B` を 1 元生成する素元 `α`」の
取り方に依らない。★★**底環 `A` も取り替えてよい**(`A` と `A'` は無関係でよい)。

証明は `σ ∈ G_n ⟺ n < i(σ)`(`mem_lowerRamificationGroup_iff_lt_ramIndex`)を
`α` と `β` の両方で書き、`ℕ∞` の外延性(`enat_eq_of_forall_natCast_lt_iff`)を使うだけ。
★★**`G_n` は `α` にも `A` にも依らない**(`𝔪^{n+1}` の inertia)という一点が効いている。

★分岐の語彙は `ramIndex` / `lowerRamificationGroup` だけで、付値も Galois も出てこない。 -/
theorem ramIndex_eq_of_adjoin_eq_top {A A' B : Type*} [CommRing A] [CommRing A'] [CommRing B]
    [Algebra A B] [Algebra A' B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] [SMulCommClass G A' B]
    {α β : B} (huα : maximalIdeal B = Ideal.span {α}) (haα : Algebra.adjoin A ({α} : Set B) = ⊤)
    (huβ : maximalIdeal B = Ideal.span {β}) (haβ : Algebra.adjoin A' ({β} : Set B) = ⊤)
    (σ : G) : ramIndex α σ = ramIndex β σ := by
  refine enat_eq_of_forall_natCast_lt_iff fun n => ?_
  rw [← mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huα haα,
    ← mem_lowerRamificationGroup_iff_lt_ramIndex (A := A') huβ haβ]

/-- ★**抽象核** —— `i` が一致すれば `G^m` も一致する(Y19e §2 に `f = id` を入れるだけ)。 -/
theorem upperRamificationGroup_congr_of_ramIndex {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [Fintype G] [MulSemiringAction G B]
    {α β : B} (h : ∀ σ : G, ramIndex α σ = ramIndex β σ) (m : ℝ) :
    upperRamificationGroup G α m = upperRamificationGroup G β m := by
  simpa using upperRamificationGroup_eq_comap (f := MonoidHom.id G) Function.surjective_id h m

/-- ★★**抽象核の到達点** —— **`G^m` は 1 元生成する素元の取り方に依らない**。
★底環も取り替えてよい。 -/
theorem upperRamificationGroup_congr_adjoin {A A' B : Type*} [CommRing A] [CommRing A']
    [CommRing B] [Algebra A B] [Algebra A' B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [Fintype G] [MulSemiringAction G B]
    [SMulCommClass G A B] [SMulCommClass G A' B] {α β : B}
    (huα : maximalIdeal B = Ideal.span {α}) (haα : Algebra.adjoin A ({α} : Set B) = ⊤)
    (huβ : maximalIdeal B = Ideal.span {β}) (haβ : Algebra.adjoin A' ({β} : Set B) = ⊤) (m : ℝ) :
    upperRamificationGroup G α m = upperRamificationGroup G β m :=
  upperRamificationGroup_congr_of_ramIndex
    (fun σ => ramIndex_eq_of_adjoin_eq_top huα haα huβ haβ σ) m

/-! ## §2 段データの選択独立性(Y19d 逸脱 1・2 の解消) -/

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**同じ `N` の 2 つの段生成元は同じ段を与える**。

★★**証明は `compat` を `M = N` に当てるだけ**である(2 行)。
Y19d の逸脱 2(「異なる `x` が同じ `S N v` を与えることは証明していない」)と、
Y19e の「★証明していないこと: `G^m` が生成元の選択に依らないこと」がこれで消える。 -/
theorem stage_eq_stage (K : PAdicLocalField p) {N : Subgroup K.absGal}
    (g g' : StageGenerator K N) (v : ℝ) : g.stage v = g'.stage v :=
  eq_of_coe_mul_coe_eq_coe (le_stage g v)
    (stage_mul_coe_eq unramifiedBaseChange K g g' le_rfl v)

/-- ★★★★**`absGalStage K N v` はどの段生成元でも計算できる**。

★★`absGalStage` の中の `Classical.choice` が外れる。**これが本ファイルの鍵**であり、
これが無いと段データの値を Lubin-Tate の生成元で計算できない。 -/
theorem absGalStage_eq_stage_any (K : PAdicLocalField p) {N : Subgroup K.absGal}
    (g : StageGenerator K N) (v : ℝ) : absGalStage K N v = g.stage v := by
  rw [absGalStage_eq_stage ⟨g⟩ v]
  exact stage_eq_stage K _ g v

/-- ★段の一意化元はどれを取ってもよい(`𝒪_L = 𝒪_{L₀}[β]`、Yoshida Lemma 5.11)。 -/
theorem adjoin_uniformizer_inertiaFixedRing_eq_top (K : PAdicLocalField p) (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {β : adjoinIntegers K x} (hβ : maximalIdeal (adjoinIntegers K x) = Ideal.span {β}) :
    Algebra.adjoin ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x))
      ({β} : Set (adjoinIntegers K x)) = ⊤ := by
  haveI := module_finite_inertiaFixedRing K x
  exact adjoin_uniformizer_eq_top hβ (exists_sub_mem_maximalIdeal_inertiaFixedRing K x)
    (exists_pow_maximalIdeal_le_map hβ
      (map_maximalIdeal_ne_bot_of_injective fixedRing_injective))

/-- ★★**段の上付き分岐群は一意化元の選び方に依らない**(Y19d 逸脱 1 の解消)。 -/
theorem stageUpperRamification_eq_map_upperRamificationGroup (K : PAdicLocalField p)
    (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {β : adjoinIntegers K x} (hβ : maximalIdeal (adjoinIntegers K x) = Ideal.span {β}) (v : ℝ) :
    stageUpperRamification K x v
      = Subgroup.map (inertiaGalAdjoin K x).subtype
          (upperRamificationGroup ↥(inertiaGalAdjoin K x) β v) := by
  rw [stageUpperRamification, upperRamificationGroup_congr_adjoin
    (A := ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x)))
    (A' := ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x)))
    (maximalIdeal_eq_span_stageUniformizer K x) (adjoin_stageUniformizer_eq_top K x)
    hβ (adjoin_uniformizer_inertiaFixedRing_eq_top K x hβ) v]

/-! ## §3 完全分岐段では `Gal(L/L₀)^v = Gal(L/K)^v` -/

/-- `K(x)/K` が完全分岐なら惰性群は `Gal(L/K)` 全体。 -/
theorem inertiaGalAdjoin_eq_top_of_isTotallyRamified (K : PAdicLocalField p) (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) : inertiaGalAdjoin K x = ⊤ :=
  inertiaGal_eq_top_of_inf_eq_bot K _
    (by rw [inf_comm]; exact totallyRamifiedAdjoin_inf_unramifiedClosure K ht)

/-- ★★★**完全分岐段の上付き分岐群は `Gal(K(x)/K)` 全体のもの**であり、
`𝒪_{K(x)}` の**どの素元 `α`** で計算してもよい。

★決め手は §1 の `ramIndex_eq_of_adjoin_eq_top` で、底環を
`𝒪_L^{I(L/K)}`(段側 = Yoshida §6.1 の設定)と `𝒪_K`(Lubin-Tate 側)で
**取り替える**ところである。両者はどちらも `𝒪_{K(x)}` を 1 元生成する
(前者は `adjoin_stageUniformizer_eq_top`、後者は
`adjoin_uniformizer_eq_top_adjoinIntegers`)。 -/
theorem stageUpperRamification_eq_upperRamificationGroup (K : PAdicLocalField p) (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {α : adjoinIntegers K x}
    (hα : maximalIdeal (adjoinIntegers K x) = Ideal.span {α}) (v : ℝ) :
    stageUpperRamification K x v
      = upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
          ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) α v := by
  have htop : inertiaGalAdjoin K x = ⊤ := inertiaGalAdjoin_eq_top_of_isTotallyRamified K x ht
  have hsurj : Function.Surjective (inertiaGalAdjoin K x).subtype := fun g =>
    ⟨⟨g, by rw [htop]; exact Subgroup.mem_top g⟩, rfl⟩
  have hram : ∀ σ : ↥(inertiaGalAdjoin K x),
      ramIndex (stageUniformizer K x) σ = ramIndex α ((inertiaGalAdjoin K x).subtype σ) :=
    fun σ => ramIndex_eq_of_adjoin_eq_top
      (A := ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x))) (A' := 𝒪[K.carrier])
      (maximalIdeal_eq_span_stageUniformizer K x) (adjoin_stageUniformizer_eq_top K x)
      hα (adjoin_uniformizer_eq_top_adjoinIntegers K x ht hα) σ
  rw [stageUpperRamification, upperRamificationGroup_eq_comap hsurj hram v,
    Subgroup.map_comap_eq, Subgroup.range_subtype, htop, top_inf_eq]

/-! ## §4 Lubin-Tate 塔での計算 -/

section LubinTate

variable (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))

/-- `Art^{-1}(U^0_K) = Γ_K`(`U^0_K = 𝒪_K^×`)。 -/
theorem absGalPrincipalLevel_zero :
    absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf 0 = ⊤ := by
  ext σ
  simp [absGalPrincipalLevel, principalUnits_zero_eq_top K]

/-- ★★**Lubin-Tate 塔の第 `m+1` 段は `K` 上 Galois**。

★★**Lubin-Tate 理論を再証明していない**: Y21 が示した「`Art^{-1}(U^{m+1}_K)` は
`Γ_K` の正規部分群」(可換群への準同型の引き戻しだから)を、
`InfiniteGalois.normal_iff_isGalois` で体の側へ降ろしただけである。 -/
theorem isGalois_adjoin_psiGenSeq (m : ℕ) :
    IsGalois K.carrier
      (IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) := by
  haveI := isGalois_closure K
  have hnorm : (IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} :
        Set K.closure)).fixingSubgroup.Normal := by
    rw [← absGalPrincipalLevel_eq_fixingSubgroup K hq hπmax hπne0 f hf0 hf1 hf m]
    exact normal_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (m + 1)
  exact (InfiniteGalois.normal_iff_isGalois _).mp hnorm

/-- 上の正規性の部分。★`IsGalois` から `Normal` を取り出すのに `.toNormal` は使えない。 -/
theorem normal_adjoin_psiGenSeq (m : ℕ) :
    Normal K.carrier
      (IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) := by
  haveI := isGalois_adjoin_psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m
  infer_instance

/-- ★★**Lubin-Tate 塔の第 `m+1` 段は `Art^{-1}(U^{m+1}_K)` の段生成元**。 -/
noncomputable def psiStageGenerator (m : ℕ) :
    StageGenerator K (absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (m + 1)) where
  gen := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
  finiteDimensional := inferInstance
  normal := normal_adjoin_psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m
  fixingSubgroup_eq := (absGalPrincipalLevel_eq_fixingSubgroup K hq hπmax hπne0 f hf0 hf1 hf m).symm

/-- ★Y22 の `algEquivRestrictSelfHom` は `AlgEquiv.restrictNormalHom` そのもの。

★Y22 は `Normal` インスタンスを立てずに済ませるために自作していた(#59 回避)。
本ファイルは §4 で正規性を得たので、両者を繋げる。 -/
theorem algEquivRestrictSelfHom_eq_restrictNormalHom
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    algEquivRestrictSelfHom K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
      = AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
          ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : Type _) := by
  refine MonoidHom.ext fun σ => AlgEquiv.ext fun z => Subtype.ext ?_
  rw [algEquivRestrictSelfHom_apply, coe_algEquivRestrictSelf, coe_restrictNormalHom_apply]

/-- ★★★★★**段データの値が Lubin-Tate 塔で計算できる** ——

    `absGalStage K (Art^{-1}(U^m_K)) (I : ℝ) = Art^{-1}(U^I_K)`   (`I ≤ m`)。

★★**これが `Found/PGC/RamificationImageUnits.lean` の `hstage` そのもの**である。

段取り:
* §2 `absGalStage_eq_stage_any` で `Classical.choice` を Lubin-Tate の生成元に取り替える。
* §3 `stageUpperRamification_eq_upperRamificationGroup` で
  `Gal(L/L₀)^I` を `Gal(K_{f,m}/K)^I`(素元は `torsionGen`)に書き換える。
* §4 `algEquivRestrictSelfHom_eq_restrictNormalHom` で制限射を揃え、
  Y22 の `comap_upperRamificationGroup_eq_absGalPrincipalLevel`(Yoshida Prop 6.14 の
  鋭い形)を当てる。

★`m = 0` は `Art^{-1}(U^0) = Γ_K` なので `le_S` だけで閉じる。 -/
theorem absGalStage_absGalPrincipalLevel (m I : ℕ) (h2 : I ≤ m) :
    absGalStage K (absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf m) ((I : ℕ) : ℝ)
      = absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf I := by
  match m, h2 with
  | 0, h2 =>
    have hI : I = 0 := Nat.le_zero.1 h2
    subst hI
    rw [absGalPrincipalLevel_zero]
    exact eq_top_iff.2 (le_absGalStage K
      (by rw [← absGalPrincipalLevel_zero K hq hπmax hπne0 f hf0 hf1 hf]
          exact absGalPrincipalLevel_mem_openNormalBase K hq hπmax hπne0 f hf0 hf1 hf 0) _)
  | (M + 1), h2 =>
    haveI := normal_adjoin_psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M
    have hα : maximalIdeal ↥(adjoinIntegers K (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt)
        = Ideal.span {torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1)
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hn
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hmem} :=
      (IsDiscreteValuationRing.irreducible_iff_uniformizer _).1
        (irreducible_torsionGen K hq hπmax hπne0 f hf0 hf1 hf (M + 1) (by omega) _
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hmem)
    rw [absGalStage_eq_stage_any K (psiStageGenerator K hq hπmax hπne0 f hf0 hf1 hf M)]
    show Subgroup.comap (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        ((IntermediateField.adjoin K.carrier
          ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt} : Set K.closure)) : Type _))
        (stageUpperRamification K (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt
          ((I : ℕ) : ℝ)) = _
    rw [stageUpperRamification_eq_upperRamificationGroup K _
        (isTotallyRamifiedAdjoin_psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M) hα,
      ← algEquivRestrictSelfHom_eq_restrictNormalHom K hq hπmax hπne0 f hf0 hf1 hf (M + 1)
        (by omega) _ (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hψ
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hn
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).hmem]
    exact comap_upperRamificationGroup_eq_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf M I h2

/-! ## §5 ★★★★★消費先の仮定 2 本が落ちる -/

/-- ★★★★★**pGC p.4 の "well-known" が仮定ゼロで立った** ——

    `Art(Γ_K^n) = U^n_K`   (`n : ℕ`)。

原文 (pGC p.4):
> Then it is well-known (Theorem 1 of [3], p. 155) that the image of Γ^v_K in Γ^ab_K is equal to U^v_K ⊆ U_K.

★`Γ_K^v` は Y19e の**本物**(`ramificationFiltration p`)であって、
退化 witness(`trivialStageFiltration` / `inertiaStageFiltration` /
`lubinTateStageFiltration`)ではない。
★逸脱: `Γ^ab_K` ではなく `𝒪_K^×` 成分(冒頭「逸脱の記録 2」)、`v` は自然数(同 1)。 -/
theorem map_ramificationFiltration_reciprocityUnits_eq_principalUnits (n : ℕ) :
    Subgroup.map (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf)
        ((ramificationFiltration p).Gv K ((n : ℕ) : ℝ))
      = principalUnits K π n :=
  map_limit_reciprocityUnits_eq_principalUnits K hq hπmax hπne0 f hf0 hf1 hf
    (absGalStageFiltration K (absGalStage_compat unramifiedBaseChange K)) n
    (fun i => absGalPrincipalLevel_mem_openNormalBase K hq hπmax hπne0 f hf0 hf1 hf (n + i))
    (fun i => absGalStage_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i) n
      (Nat.le_add_right n i))

/-- ★同じことを原典の記号 `U^v_K` で。 -/
theorem map_ramificationFiltration_reciprocityUnits_eq_upperPrincipalUnits (n : ℕ) :
    Subgroup.map (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf)
        ((ramificationFiltration p).Gv K ((n : ℕ) : ℝ))
      = upperPrincipalUnits K π ((n : ℕ) : ℝ) :=
  map_limit_reciprocityUnits_eq_upperPrincipalUnits K hq hπmax hπne0 f hf0 hf1 hf
    (absGalStageFiltration K (absGalStage_compat unramifiedBaseChange K)) n
    (fun i => absGalPrincipalLevel_mem_openNormalBase K hq hπmax hπne0 f hf0 hf1 hf (n + i))
    (fun i => absGalStage_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i) n
      (Nat.le_add_right n i))

/-- ★**包含の側**: `Γ_K^n ⊆ Art^{-1}(U^n_K)`。★等号ではない(右辺は `ker(Art)` を含む)。 -/
theorem ramificationFiltration_Gv_le_absGalPrincipalLevel (n : ℕ) :
    (ramificationFiltration p).Gv K ((n : ℕ) : ℝ)
      ≤ absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n :=
  limit_le_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf
    (absGalStageFiltration K (absGalStage_compat unramifiedBaseChange K)) n
    (absGalPrincipalLevel_mem_openNormalBase K hq hπmax hπne0 f hf0 hf1 hf n)
    (le_of_eq (absGalStage_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n n le_rfl))

/-- ★★★**退化の自己検査** —— `n = 0` を主定理に入れると

    `Art(I_K) = 𝒪_K^×`

が出る(`Γ_K^0 = I_K` は `Found/PGC/RamificationFiltrationZero.lean`、
`U^0_K = 𝒪_K^×` は `principalUnits_zero_eq_top`)。

★★**これは古典的に正しい主張**(局所類体論で相互写像は惰性群を単数群へ写す)であり、
本ファイルの計算が退化していないことの独立な確認になっている。
★退化 witness `trivialStageFiltration` では `Γ^0 = ⊤` なので、この式は
「`Art(Γ_K) = 𝒪_K^×`」という**別の(より弱い)主張**に化けてしまう。 -/
theorem map_absInertia_reciprocityUnits_eq_top :
    Subgroup.map (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf) (absInertia K) = ⊤ := by
  have h := map_ramificationFiltration_reciprocityUnits_eq_principalUnits K hq hπmax hπne0
    f hf0 hf1 hf 0
  rw [Nat.cast_zero, ramificationFiltration_Gv_zero K, principalUnits_zero_eq_top K] at h
  exact h

/-- ★体の言葉で: **`Γ_K^{m+1}` は Lubin-Tate 塔の第 `m+1` 段 `K_{f,m+1}` を各点固定する**。 -/
theorem ramificationFiltration_Gv_le_fixingSubgroup (m : ℕ) :
    (ramificationFiltration p).Gv K (((m + 1 : ℕ) : ℕ) : ℝ)
      ≤ (IntermediateField.adjoin K.carrier
          ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)).fixingSubgroup := by
  rw [← absGalPrincipalLevel_eq_fixingSubgroup K hq hπmax hπne0 f hf0 hf1 hf m]
  exact ramificationFiltration_Gv_le_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (m + 1)

end LubinTate

end ABC3.Found.PGC
