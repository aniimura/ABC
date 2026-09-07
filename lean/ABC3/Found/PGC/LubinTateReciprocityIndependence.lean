import ABC3.Found.PGC.LubinTateUniformizerIndependence
import ABC3.Found.PGC.WeilReciprocityExtension

/-!
# `ρf,m = ρf′,m` —— 相互写像は Lubin-Tate 系の取り方に依らない(Yoshida 2008 Corollary 4.9 後半)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Corollary 4.9(物理 p.9)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-4.html` の `#cor-4-9`。

原文 (Yoshida08 p.9):
> Corollary 4.9. The K^m_f and ρ_f,m, hence also K^LT_f and ρ_f, of Proposition 4.7(ii)
> do not depend on f. (We will drop the subscript f and write K^m, ρ_m, K^LT and ρ.)

原典の証明(p.9、`0_Source` の `.txt` の 516–529 行を直読した全文):

> Proof. For f, f′ with linear coefficients π, π′, take θ ∈ Θ^{K̂,×}_{π,π′} and
> [θ] : µ_{f,m} ≅ µ_{f′,m} by Proposition 4.8. Lemma 4.6 shows K̂^m_f = K̂^m_{f′}.
> If σ(α) = [xπj](α) for σ ∈ W(Kmf/K), then
> σ([θ](α)) = [θ](j)[xπj](α) = [xπ′j][θ](α) by Lemma 4.5, hence ρf,m = ρf′,m. □

前半(体の一致 `K̂^m_f = K̂^m_{f′}`)は `Found/PGC/LubinTateTowerFIndependent.lean` の
`lubinTateCompletionField_eq` と `Found/PGC/LubinTateUniformizerIndependence.lean` の
`lubinTateLevelField_sup_unramifiedClosure_eq_of_isUnit_mul` が既に出している。
本ファイルが担当するのは残りの**写像の一致**である。

## ★★★ 射程 —— `ρ` の定義域には `K^ur` が最初から入っている

直前の節点(`LubinTateUniformizerIndependence.lean`)は、体の側では
`⊔ K^ur` を落とすと主張が偽になることを見つけた。★では `ρ` の側はどうか。
`.txt` を直読して**射程を決める文**を特定した:

* `.txt` 433 行(Proposition 4.7(ii) の冒頭)—— 逐語:
  > (ii) Let L = K̂. The ρf,m of Proposition 4.4(iii) extend to isomorphisms:
  > ρf,m : W(Kmf/K) ≅ K^×/(1 + p^m).
  > (ϕ^j on K̂, α ↦ [xπj](α), ∀α ∈ µ_{f,m}) ↦ x mod 1 + p^m  (v(x) = −j)
* `.txt` 530 行(Definition 4.10 の冒頭)—— 逐語:
  > Definition 4.10. For any f ∈ 𝒪_L[X] with L/K finite, set K^m := K^ur L^m_f.

★したがって `ρf,m` の定義域は最初から `W(Kmf/K)`、すなわち
`K̂^m_f = K̂ · L^m_f ⊇ K̂`(= `K^ur` の完備化)の上の Weil 群である。
代数側の対応物は Definition 4.10 の `K^m = K^ur L^m_f` の上の `W(K^m/K)`。

★★**結論: `ρ` の側では `⊔ K^ur` を「足す」必要が無い。最初から入っている。**
体の側の反例(`K = ℚ_p`、`π = p`、`π′ = −p` で `K_p ≠ K_{−p}`)は `ρ` を脅かさない
——`Art(−1)` が `µ_{p^∞}` 上で反転として働くという現象は、`K^ur` を含む体の上で
見れば `Frobenius` の寄与に吸収される。★逆に、`K^ur` を落として
「`Gal(K_π/K)` の上の相互写像が `π` に依らない」と述べると、そもそも
`K_π ≠ K_{π′}` なので**共通の定義域が無く、主張が型として立たない**。

★★退化検査(D1): `j = 0`(慣性部分 `Gal(K̂^m_f/K̂) ≅ (𝒪/𝔭^m)^×`)では
`π_0 = π′_0 = 1` なので Lemma 4.5 は自明になり、独立性は `[θ]` の
`𝒪`-線型性だけから出る(`smul_spec_transfer_of_one`)。★本 Corollary の
内容はすべて `j ≠ 0`、すなわち `K^ur` の側にある。

## ★ 原典が実際に使っている等式は 1 本だけ

証明の心臓は `[θ](j) ∘ [xπj] = [xπ′j] ∘ [θ]`。冪級数の合成則
`[a] ∘ [b] = [ab]`(Proposition 3.5、木では `LubinTateAction_comp`)で潰すと

  `θ^{(j)} · (x π_j) = (x π′_j) · θ`

になり、これは Lemma 4.5 第 1 主張 `θ^{(j)} π_j = θ π′_j`
(`Found/PGC/UniformizerCocycle.lean` の `zpow_apply_mul_uniformizerZ`)に
`x` を掛けただけである。★`x` の側には `ϕ`-不変性が要らない
(`x` は `ϕ^j` の外側にいる)——`WeilReciprocityExtension.lean` の
`zpow_apply_mul_uniformizerZ_mul` が準同型性の側で見つけたのと同じ現象。

## 抽象核と具体層

| 段 | 宣言 | 分岐・付値・Galois・Lubin-Tate の語彙 |
|---|---|---|
| 抽象核 A | `smul_spec_transfer` | 出てこない(`CommMonoid` の作用と 1 本の関数だけ) |
| 抽象核 B | `existsUnique_of_generator` | 出てこない |
| 抽象核 C | `eq_of_spec_of_existsUnique` | 出てこない(型 2 つと述語 1 つ) |
| 抽象核 D | `reciprocity_eq_of_intertwiner` | 出てこない(★これが Corollary 4.9 後半の全内容) |
| 抽象核 E | `reciprocity_eq_of_intertwiner_of_generator` | 出てこない |
| 橋 | `weilGroup_reciprocity_eq_of_intertwiner` | `weilGroup`(抽象 Weil 群)だけ |
| 具体層 1 | `zpow_mul_uniformizerZ_units` / `map_zpow_mul_uniformizerZ` | `π_j`・`Θ`(Lemma 4.5 の消費) |
| 具体層 2 | `reciprocity_eq_of_thetaSet` | 同上(★原典の字面) |
| 具体層 3 | `fixingSubgroup_lubinTateLevelField_sup_unramifiedClosure_eq_of_isUnit_mul` | 木の Lubin-Tate 機械 |

## 逸脱の記録

1. ★**`µ_{f,m}` と `[·]` は作らず、それらが供給する型・作用・全射性を仮定に取った。**
   `Found/PGC/WeilReciprocityExtension.lean` の `prop_4_7_i_equiv`(Proposition 4.7(i))が
   同じ設計を取っており(そちらの逸脱 2)、本ファイルはその続きなので合わせた。
   具体的には `MulAction S N`(`S` = 係数のなす可換モノイド、`N` = 捩れ点の住む型)と
   `ι : Lˣ →* S`(`u ↦ [u]`)を仮定に置く。
2. ★**`[u]` は原典では `v(u) = 0` のときにしか定義されない**が、上の `ι` は `Lˣ` 全体で
   モノイド準同型であることを要求している。消費側は `S := Function.End N`
   (`N` = `𝒪_{ℂ_K}` の極大イデアル)と取れば、`[u]` は `u` が単元でなくとも
   冪級数として `N` 上で評価できるので、この形で使える。★この逸脱は
   「仮定を強めた」向きなので、消費側が供給できなければ弱い版
   (`smul_spec_transfer` を直に使う)へ落とせる。
3. ★**`ρ` の値域 `K^×/(1+𝔭^m)` は抽象型 `Q` と代表元写像 `rep : Q → S` で表した。**
   原典の `x mod 1+𝔭^m` は `Q` の元、`[xπj]` に現れる `x` は `rep` の像である。
   `rep` に準同型性も単射性も要求していない(要らない)。
4. ★**`W(Kmf/K)` は抽象群 `G` と次数写像 `deg : G → ℤ` で表した。**
   `G` が群であることすら使っていない(`deg` は単なる関数でよい)。
   Weil 群の形にした版は `weilGroup_reciprocity_eq_of_intertwiner`。
5. 原典 Definition 4.10 は `f ∈ 𝒪_L[X]`(`L/K` 有限)を許すが、具体層 3 は
   直前の節点に合わせて `L = K` の場合だけを扱う。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal Valued Classical

def reciprocity_eq_of_intertwiner.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

def reciprocity_eq_of_thetaSet.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

def fixingSubgroup_lubinTateLevelField_sup_unramifiedClosure_eq_of_isUnit_mul.src :
    ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

def reciprocity_eq_of_bracket.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

def reciprocity_eq_of_thetaSet_bracket.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

/-! ## 1. 抽象核 —— ねじれた同型に沿った「特徴づけ」の輸送

★分岐・付値・Galois・Lubin-Tate の語彙が 1 つも出てこない。
可換モノイド `S` の作用と、1 本の関数 `σ : N → N` だけである。 -/

section AbstractCore

/-- ★★★**抽象核 A —— 原典の証明の 1 行そのもの**。

> σ([θ](α)) = [θ](j)[xπj](α) = [xπ′j][θ](α)

`σ` が `θ` を「`c` へずらしながら」通り抜ける(`hsemi`)とき、
`M` の上で `σ = a •` なら、`θ • M` の上では `σ = a′ •` になる。
必要なのは `c * a = a′ * θ` という**係数の 1 本の等式**だけ
(原典ではこれが Lemma 4.5 である)。

★`S` は可換であることしか使わない。`N` にも何も要らない。
★`σ` が全単射である必要も、`M` が部分加群である必要も無い。 -/
theorem smul_spec_transfer {S : Type*} [CommMonoid S] {N : Type*} [MulAction S N]
    (σ : N → N) (θ c a a' : S) (M : Set N)
    (hsemi : ∀ α : N, σ (θ • α) = c • σ α)
    (hθ : c * a = a' * θ)
    (hσ : ∀ α ∈ M, σ α = a • α) :
    ∀ β ∈ (fun α => θ • α) '' M, σ β = a' • β := by
  rintro _ ⟨α, hα, rfl⟩
  rw [hsemi, hσ α hα, ← mul_smul, ← mul_smul, hθ]

/-- ★★**退化検査 (D1) —— `j = 0`(慣性部分)では Lemma 4.5 は要らない**。

`π_0 = π′_0 = 1` なので `c = θ`(すなわち `θ^{(0)} = θ`)と `a′ = a` が取れ、
係数の等式 `θ * a = a * θ` は可換性だけから出る。
★本 Corollary の内容がすべて `j ≠ 0`(= `K^ur` の側)にあることの形式的な確認。 -/
theorem smul_spec_transfer_of_one {S : Type*} [CommMonoid S] {N : Type*} [MulAction S N]
    (σ : N → N) (θ a : S) (M : Set N)
    (hsemi : ∀ α : N, σ (θ • α) = θ • σ α)
    (hσ : ∀ α ∈ M, σ α = a • α) :
    ∀ β ∈ (fun α => θ • α) '' M, σ β = a • β :=
  smul_spec_transfer σ θ θ a a M hsemi (mul_comm _ _) hσ

/-- ★★**抽象核 B —— 一意性は「生成元 1 個」から出る**。

原典 Proposition 4.7(ii) の
> σ(α) = [xπj](α) for a unique x mod 1 + p^m by (i)
に対応する。`(i)` が言っているのは「`q ↦ rep q • β₀` が単射」であり、
それ以上のことは使わない(★`M'` の全体を見る必要が無い)。 -/
theorem existsUnique_of_generator {Q S N : Type*} [Monoid S] [MulAction S N]
    (rep : Q → S) (σ : N → N) (M' : Set N) (β₀ : N) (hβ₀ : β₀ ∈ M')
    (hinj : Function.Injective fun q : Q => rep q • β₀)
    (hex : ∃ q : Q, ∀ β ∈ M', σ β = rep q • β) :
    ∃! q : Q, ∀ β ∈ M', σ β = rep q • β := by
  obtain ⟨q, hq⟩ := hex
  refine ⟨q, hq, fun q' hq' => hinj ?_⟩
  show rep q' • β₀ = rep q • β₀
  rw [← hq' β₀ hβ₀, hq β₀ hβ₀]

/-- ★**抽象核 C —— 同じ特徴づけを満たす 2 本の写像は一致する**。

型 2 つと述語 1 つしか出てこない。★`ρ = ρ′` という結論の形は
これが供給する(原典の "hence ρf,m = ρf′,m")。 -/
theorem eq_of_spec_of_existsUnique {G Q : Type*} {Spec : G → Q → Prop}
    (ρ ρ' : G → Q) (hρ : ∀ g, Spec g (ρ g)) (hρ' : ∀ g, Spec g (ρ' g))
    (huniq : ∀ g, ∃! q, Spec g q) : ρ = ρ' :=
  funext fun g => (huniq g).unique (hρ g) (hρ' g)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**抽象核 D —— Corollary 4.9 後半の全内容**。

原文 (Yoshida08 p.9):
> If σ(α) = [xπj](α) for σ ∈ W(Kmf/K), then
> σ([θ](α)) = [θ](j)[xπj](α) = [xπ′j][θ](α) by Lemma 4.5, hence ρf,m = ρf′,m.

記号の対応:

| 原典 | ここ |
|---|---|
| `W(Kmf/K)` | `G`(★群である必要すら無い) |
| `σ ↦ (σ\|_{K̂} = ϕ^j)` | `deg : G → ℤ` |
| `K^×/(1+𝔭^m)` | `Q`、`x` の取り出しは `rep : Q → S` |
| `µ_{f,m}` / `µ_{f′,m}` | `M` / `M′`(`N` の部分集合) |
| `[a]` | `S` の `N` への作用 |
| `[θ]` | `θ • ·` |
| `[θ](j)` | `c j • ·` |
| `π_j` / `π′_j` | `t j` / `t′ j` |
| Lemma 4.5 | `hθ : ∀ j, c j * t j = t′ j * θ` |
| Lemma 4.6(`[θ]` が同型) | `hMM′ : M′ ⊆ θ • M` |
| Proposition 4.7(ii) の一意性 | `huniq` |

★`hMM′` は「`[θ]` が全射」しか使わない(単射性も `𝒪`-線型性も不要)。
★`deg` が準同型であることも使わない。 -/
theorem reciprocity_eq_of_intertwiner {G Q S N : Type*} [CommMonoid S] [MulAction S N]
    (deg : G → ℤ) (act : G → N → N) (rep : Q → S)
    (θ : S) (c t t' : ℤ → S) (M M' : Set N)
    (hsemi : ∀ (g : G) (α : N), act g (θ • α) = c (deg g) • act g α)
    (hθ : ∀ j : ℤ, c j * t j = t' j * θ)
    (hMM' : M' ⊆ (fun α => θ • α) '' M)
    (ρ ρ' : G → Q)
    (hρ : ∀ g, ∀ α ∈ M, act g α = (rep (ρ g) * t (deg g)) • α)
    (hρ' : ∀ g, ∀ β ∈ M', act g β = (rep (ρ' g) * t' (deg g)) • β)
    (huniq : ∀ g, ∃! q : Q, ∀ β ∈ M', act g β = (rep q * t' (deg g)) • β) :
    ρ = ρ' := by
  refine eq_of_spec_of_existsUnique
    (Spec := fun g q => ∀ β ∈ M', act g β = (rep q * t' (deg g)) • β) ρ ρ' ?_ hρ' huniq
  intro g β hβ
  refine smul_spec_transfer (act g) θ (c (deg g)) (rep (ρ g) * t (deg g))
    (rep (ρ g) * t' (deg g)) M (hsemi g) ?_ (hρ g) β (hMM' hβ)
  calc c (deg g) * (rep (ρ g) * t (deg g))
      = rep (ρ g) * (c (deg g) * t (deg g)) := mul_left_comm _ _ _
    _ = rep (ρ g) * (t' (deg g) * θ) := by rw [hθ]
    _ = rep (ρ g) * t' (deg g) * θ := (mul_assoc _ _ _).symm

/-- ★★**抽象核 E** —— 抽象核 D の `huniq` を抽象核 B(生成元 1 個)で置き換えた形。

原典が Proposition 4.7(i) を「`α ∈ µ^×_{f,m}` を 1 つ取る」形で使うのに合わせた。
`β₀` は `[θ](α)`(`α` は原典の生成元)に当たる。 -/
theorem reciprocity_eq_of_intertwiner_of_generator
    {G Q S N : Type*} [CommMonoid S] [MulAction S N]
    (deg : G → ℤ) (act : G → N → N) (rep : Q → S)
    (θ : S) (c t t' : ℤ → S) (M M' : Set N) (β₀ : N) (hβ₀ : β₀ ∈ M')
    (hsemi : ∀ (g : G) (α : N), act g (θ • α) = c (deg g) • act g α)
    (hθ : ∀ j : ℤ, c j * t j = t' j * θ)
    (hMM' : M' ⊆ (fun α => θ • α) '' M)
    (hinj : ∀ g : G, Function.Injective fun q : Q => (rep q * t' (deg g)) • β₀)
    (ρ ρ' : G → Q)
    (hρ : ∀ g, ∀ α ∈ M, act g α = (rep (ρ g) * t (deg g)) • α)
    (hρ' : ∀ g, ∀ β ∈ M', act g β = (rep (ρ' g) * t' (deg g)) • β) :
    ρ = ρ' :=
  reciprocity_eq_of_intertwiner deg act rep θ c t t' M M' hsemi hθ hMM' ρ ρ' hρ hρ'
    (fun g => existsUnique_of_generator (fun q => rep q * t' (deg g)) (act g) M' β₀ hβ₀
      (hinj g) ⟨ρ' g, hρ' g⟩)


/-! ## 1b. 抽象核 A′/D′ —— `MulAction` インスタンスを要求しない形

★消費側(冪級数)にとって `[a]` は**ただの関数** `S → N → N` であり、
成り立つのは合成則 `[a] ∘ [b] = [ab]`(原典 Proposition 3.5、木では
`LubinTateAction_comp`)だけである。`MulAction S N` のインスタンスを
探しに行かせないため、合成則を仮定に書いた版を用意する。
★中身は 1 節と同じ 1 行の計算である。 -/

/-- ★★★**抽象核 A′** —— 抽象核 A の「作用を関数で書いた」版。

`hbr` は原典 Proposition 3.5 の `[a] ∘ [b] = [ab]`。
★`br 1 = id` も、`br` が単射であることも使わない。 -/
theorem map_spec_transfer {S : Type*} [CommMonoid S] {N : Type*}
    (br : S → N → N) (hbr : ∀ (a b : S) (n : N), br a (br b n) = br (a * b) n)
    (σ : N → N) (θ c a a' : S) (M : Set N)
    (hsemi : ∀ α : N, σ (br θ α) = br c (σ α))
    (hθ : c * a = a' * θ)
    (hσ : ∀ α ∈ M, σ α = br a α) :
    ∀ β ∈ (fun α => br θ α) '' M, σ β = br a' β := by
  rintro _ ⟨α, hα, rfl⟩
  rw [hsemi, hσ α hα, hbr, hθ, ← hbr]

/-- ★**抽象核 B′** —— 抽象核 B の「作用を関数で書いた」版。 -/
theorem existsUnique_of_generator_bracket {Q S N : Type*} (br : S → N → N)
    (rep : Q → S) (σ : N → N) (M' : Set N) (β₀ : N) (hβ₀ : β₀ ∈ M')
    (hinj : Function.Injective fun q : Q => br (rep q) β₀)
    (hex : ∃ q : Q, ∀ β ∈ M', σ β = br (rep q) β) :
    ∃! q : Q, ∀ β ∈ M', σ β = br (rep q) β := by
  obtain ⟨q, hq⟩ := hex
  refine ⟨q, hq, fun q' hq' => hinj ?_⟩
  show br (rep q') β₀ = br (rep q) β₀
  rw [← hq' β₀ hβ₀, hq β₀ hβ₀]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**抽象核 D′ —— Corollary 4.9 後半(合成則を仮定に書いた版)**。

抽象核 D と同じ内容だが、`MulAction S N` の代わりに
`br : S → N → N` と合成則 `hbr` を取る。★冪級数の側から使うときはこちらである。 -/
theorem reciprocity_eq_of_bracket {G Q S N : Type*} [CommMonoid S]
    (br : S → N → N) (hbr : ∀ (a b : S) (n : N), br a (br b n) = br (a * b) n)
    (deg : G → ℤ) (act : G → N → N) (rep : Q → S)
    (θ : S) (c t t' : ℤ → S) (M M' : Set N)
    (hsemi : ∀ (g : G) (α : N), act g (br θ α) = br (c (deg g)) (act g α))
    (hθ : ∀ j : ℤ, c j * t j = t' j * θ)
    (hMM' : M' ⊆ (fun α => br θ α) '' M)
    (ρ ρ' : G → Q)
    (hρ : ∀ g, ∀ α ∈ M, act g α = br (rep (ρ g) * t (deg g)) α)
    (hρ' : ∀ g, ∀ β ∈ M', act g β = br (rep (ρ' g) * t' (deg g)) β)
    (huniq : ∀ g, ∃! q : Q, ∀ β ∈ M', act g β = br (rep q * t' (deg g)) β) :
    ρ = ρ' := by
  refine eq_of_spec_of_existsUnique
    (Spec := fun g q => ∀ β ∈ M', act g β = br (rep q * t' (deg g)) β) ρ ρ' ?_ hρ' huniq
  intro g β hβ
  refine map_spec_transfer br hbr (act g) θ (c (deg g)) (rep (ρ g) * t (deg g))
    (rep (ρ g) * t' (deg g)) M (hsemi g) ?_ (hρ g) β (hMM' hβ)
  calc c (deg g) * (rep (ρ g) * t (deg g))
      = rep (ρ g) * (c (deg g) * t (deg g)) := mul_left_comm _ _ _
    _ = rep (ρ g) * (t' (deg g) * θ) := by rw [hθ]
    _ = rep (ρ g) * t' (deg g) * θ := (mul_assoc _ _ _).symm

/-! ### 群準同型として述べる版と、段数の極限へ渡す段 -/

/-- ★**`ρ` を群準同型として述べた版**。木の `reciprocityMapLimit` は `MonoidHom` なので、
消費側はこちらを使う。★`DFunLike.coe_injective` で包むだけ。 -/
theorem reciprocityHom_eq_of_intertwiner {G Q S N : Type*} [Monoid G] [Monoid Q]
    [CommMonoid S] [MulAction S N]
    (deg : G → ℤ) (act : G → N → N) (rep : Q → S)
    (θ : S) (c t t' : ℤ → S) (M M' : Set N)
    (hsemi : ∀ (g : G) (α : N), act g (θ • α) = c (deg g) • act g α)
    (hθ : ∀ j : ℤ, c j * t j = t' j * θ)
    (hMM' : M' ⊆ (fun α => θ • α) '' M)
    (ρ ρ' : G →* Q)
    (hρ : ∀ g, ∀ α ∈ M, act g α = (rep (ρ g) * t (deg g)) • α)
    (hρ' : ∀ g, ∀ β ∈ M', act g β = (rep (ρ' g) * t' (deg g)) • β)
    (huniq : ∀ g, ∃! q : Q, ∀ β ∈ M', act g β = (rep q * t' (deg g)) • β) :
    ρ = ρ' :=
  DFunLike.coe_injective
    (reciprocity_eq_of_intertwiner deg act rep θ c t t' M M' hsemi hθ hMM'
      (fun g => ρ g) (fun g => ρ' g) hρ hρ' huniq)

/-- ★★**原典の "hence also K^LT_f and ρ_f" —— 段数の極限へ渡す段**。

各段 `m` で `ρ_m = ρ′_m` なら、それらを束ねた `ρ = ρ′`。
★必要なのは「成分で分離する」(`hsep`)だけで、射影系であることも
`Q` が位相群であることも使わない。 -/
theorem eq_of_forall_component {G Q : Type*} {ι : Type*} {Qm : ι → Type*}
    (comp : ∀ i, Q → Qm i) (hsep : ∀ q q' : Q, (∀ i, comp i q = comp i q') → q = q')
    (ρ ρ' : G → Q) (h : ∀ i, (fun g => comp i (ρ g)) = fun g => comp i (ρ' g)) :
    ρ = ρ' :=
  funext fun g => hsep _ _ fun i => congrFun (h i) g

end AbstractCore

/-! ## 2. 橋 —— 抽象 Weil 群の形に直す

`Found/PGC/WeilReciprocityExtension.lean` の `weilGroup` / `weilDeg`
(Proposition 4.7 の抽象化)に載せ替えるだけ。★抽象核 D をそのまま流用する。 -/

section WeilBridge

variable {Γ Δ : Type*} [Group Γ] [Group Δ]

/-- ★★**Corollary 4.9 後半、抽象 Weil 群の上の形**。

`G := W(Kmf/K)` を `weilGroup res ϕ` に、`deg` を `weilDeg`(の加法表示)に
取っただけ。★★原典の `σ|_{K̂} = ϕ^j` はここでは `weilDeg` の値そのものである。 -/
theorem weilGroup_reciprocity_eq_of_intertwiner
    (res : Γ →* Δ) (ϕ : Δ) {Q S N : Type*} [CommMonoid S] [MulAction S N]
    (act : weilGroup res ϕ → N → N) (rep : Q → S)
    (θ : S) (c t t' : ℤ → S) (M M' : Set N)
    (hsemi : ∀ (g : weilGroup res ϕ) (α : N),
      act g (θ • α) = c (Multiplicative.toAdd (weilDeg res ϕ g)) • act g α)
    (hθ : ∀ j : ℤ, c j * t j = t' j * θ)
    (hMM' : M' ⊆ (fun α => θ • α) '' M)
    (ρ ρ' : weilGroup res ϕ → Q)
    (hρ : ∀ g, ∀ α ∈ M,
      act g α = (rep (ρ g) * t (Multiplicative.toAdd (weilDeg res ϕ g))) • α)
    (hρ' : ∀ g, ∀ β ∈ M',
      act g β = (rep (ρ' g) * t' (Multiplicative.toAdd (weilDeg res ϕ g))) • β)
    (huniq : ∀ g, ∃! q : Q, ∀ β ∈ M',
      act g β = (rep q * t' (Multiplicative.toAdd (weilDeg res ϕ g))) • β) :
    ρ = ρ' :=
  reciprocity_eq_of_intertwiner
    (fun g => Multiplicative.toAdd (weilDeg res ϕ g)) act rep θ c t t' M M'
    hsemi hθ hMM' ρ ρ' hρ hρ' huniq

end WeilBridge

/-! ## 3. 具体層 1 —— Lemma 4.5 が係数の等式 `hθ` を供給する

`Found/PGC/UniformizerCocycle.lean` の `uniformizerZ`(= `π_j`)と
`ThetaSet`(= `Θ^L_{π,π′}`)の上で、抽象核 D の `hθ` を作る。 -/

section Concrete

variable {L : Type*} [Field L]

/-- ★★★**Lemma 4.5 の単元群での形** `θ^{(j)} π_j = π′_j θ`(`Lˣ` の等式)。

`Found/PGC/UniformizerCocycle.lean` の `zpow_apply_mul_uniformizerZ`
(`L` の等式)を `Units.ext` で持ち上げただけ。★★これが本ファイルが
Lemma 4.5 を消費する唯一の場所である。

★`(ϕ^j)` の側は `unitsRingAutHom`(`RingAut L →* MulAut Lˣ`)を通す。
群準同型として作ってあるので `(unitsRingAutHom L ϕ)^j` と `ϕ^j` が一致する
(`coe_unitsRingAutHom_zpow`)。 -/
theorem zpow_mul_uniformizerZ_units (ϕ : RingAut L) (π π' θ : Lˣ)
    (hθ : ((θ : L)) ∈ ThetaSet ϕ (π : L) (π' : L)) (j : ℤ) :
    ((unitsRingAutHom L ϕ) ^ j) θ * uniformizerZ ϕ π j = uniformizerZ ϕ π' j * θ := by
  apply Units.ext
  rw [Units.val_mul, Units.val_mul, coe_unitsRingAutHom_zpow,
    zpow_apply_mul_uniformizerZ ϕ π π' hθ j, mul_comm]

/-- ★★**係数のモノイドへ運んだ形** —— 抽象核 D の `hθ` そのもの。

`ι : Lˣ →* S` は原典の `u ↦ [u]`(冪級数の合成をモノイド演算とする)に当たる。
★逸脱 2: 原典の `[u]` は `v(u) = 0` のときにしか整係数でないが、`ι` は
`Lˣ` 全体で定義された準同型として要求している。消費側は
`S := Function.End N`(`N` は `𝒪_{ℂ_K}` の極大イデアル)と取れる。 -/
theorem map_zpow_mul_uniformizerZ {S : Type*} [CommMonoid S] (ι : Lˣ →* S)
    (ϕ : RingAut L) (π π' θ : Lˣ)
    (hθ : ((θ : L)) ∈ ThetaSet ϕ (π : L) (π' : L)) (j : ℤ) :
    ι (((unitsRingAutHom L ϕ) ^ j) θ) * ι (uniformizerZ ϕ π j)
      = ι (uniformizerZ ϕ π' j) * ι θ := by
  rw [← map_mul, ← map_mul, zpow_mul_uniformizerZ_units ϕ π π' θ hθ j]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Corollary 4.9 後半、原典の字面の形** `ρf,m = ρf′,m`。

原文 (Yoshida08 p.9):
> If σ(α) = [xπj](α) for σ ∈ W(Kmf/K), then
> σ([θ](α)) = [θ](j)[xπj](α) = [xπ′j][θ](α) by Lemma 4.5, hence ρf,m = ρf′,m.

抽象核 D に、`π_j` と `Θ^L_{π,π′}` を代入しただけ。
`hθ`(= Lemma 4.5)は `map_zpow_mul_uniformizerZ` が供給する。

★★射程: `deg` は原典の `σ|_{K̂} = ϕ^j` の `j` であり、`j ≠ 0` を許す。
すなわち `G` は**慣性部分ではなく Weil 群全体**である
(= `K^ur` が定義域に入っている)。ファイル冒頭の射程の議論を見よ。 -/
theorem reciprocity_eq_of_thetaSet {G Q S N : Type*} [CommMonoid S] [MulAction S N]
    (ϕ : RingAut L) (π π' θ : Lˣ)
    (hθ : ((θ : L)) ∈ ThetaSet ϕ (π : L) (π' : L))
    (ι : Lˣ →* S)
    (deg : G → ℤ) (act : G → N → N) (rep : Q → S) (M M' : Set N)
    (hsemi : ∀ (g : G) (α : N),
      act g (ι θ • α) = ι (((unitsRingAutHom L ϕ) ^ (deg g)) θ) • act g α)
    (hMM' : M' ⊆ (fun α => ι θ • α) '' M)
    (ρ ρ' : G → Q)
    (hρ : ∀ g, ∀ α ∈ M, act g α = (rep (ρ g) * ι (uniformizerZ ϕ π (deg g))) • α)
    (hρ' : ∀ g, ∀ β ∈ M', act g β = (rep (ρ' g) * ι (uniformizerZ ϕ π' (deg g))) • β)
    (huniq : ∀ g, ∃! q : Q, ∀ β ∈ M',
      act g β = (rep q * ι (uniformizerZ ϕ π' (deg g))) • β) :
    ρ = ρ' :=
  reciprocity_eq_of_intertwiner deg act rep (ι θ)
    (fun j => ι (((unitsRingAutHom L ϕ) ^ j) θ))
    (fun j => ι (uniformizerZ ϕ π j)) (fun j => ι (uniformizerZ ϕ π' j)) M M'
    hsemi (map_zpow_mul_uniformizerZ ι ϕ π π' θ hθ) hMM' ρ ρ' hρ hρ' huniq

/-- ★★**退化検査 (D2) —— `deg g = 0`(慣性部分)では `π_j` が消える**。

`uniformizerZ ϕ π 0 = 1`(`uniformizerZ_zero`)なので、
`ρ` の特徴づけは `act g α = rep (ρ g) • α` になり、`π` も `π′` も現れない。
★★これが「本 Corollary の内容はすべて `j ≠ 0` にある」ことの形式的な確認である。 -/
theorem uniformizerZ_deg_zero_eq_one (ϕ : RingAut L) (π : Lˣ)
    {S : Type*} [CommMonoid S] (ι : Lˣ →* S) :
    ι (uniformizerZ ϕ π 0) = 1 := by
  rw [uniformizerZ_zero, map_one]

/-- ★退化検査 (D3) —— `π′ = π` なら `θ` は `ϕ`-不変で、
Lemma 4.5 は「`θ^{(j)} = θ`」に退化する。 -/
theorem zpow_eq_self_of_mem_thetaSet_self (ϕ : RingAut L) (π θ : Lˣ)
    (hθ : ((θ : L)) ∈ ThetaSet ϕ (π : L) (π : L)) (j : ℤ) :
    ((unitsRingAutHom L ϕ) ^ j) θ = θ := by
  have h := zpow_mul_uniformizerZ_units ϕ π π θ hθ j
  exact mul_right_cancel (h.trans (mul_comm (uniformizerZ ϕ π j) θ))


/-- ★★★**Corollary 4.9 後半、合成則版**(消費側が実際に使う形)。

`reciprocity_eq_of_thetaSet` の `MulAction` を `br` + 合成則に置き換えただけ。
★冪級数 `[·]` は `MulAction` インスタンスを持たないので、木の Lubin-Tate 機械から
使うときはこちらになる。 -/
theorem reciprocity_eq_of_thetaSet_bracket {G Q S N : Type*} [CommMonoid S]
    (br : S → N → N) (hbr : ∀ (a b : S) (n : N), br a (br b n) = br (a * b) n)
    (ϕ : RingAut L) (π π' θ : Lˣ)
    (hθ : ((θ : L)) ∈ ThetaSet ϕ (π : L) (π' : L))
    (ι : Lˣ →* S)
    (deg : G → ℤ) (act : G → N → N) (rep : Q → S) (M M' : Set N)
    (hsemi : ∀ (g : G) (α : N),
      act g (br (ι θ) α) = br (ι (((unitsRingAutHom L ϕ) ^ (deg g)) θ)) (act g α))
    (hMM' : M' ⊆ (fun α => br (ι θ) α) '' M)
    (ρ ρ' : G → Q)
    (hρ : ∀ g, ∀ α ∈ M, act g α = br (rep (ρ g) * ι (uniformizerZ ϕ π (deg g))) α)
    (hρ' : ∀ g, ∀ β ∈ M', act g β = br (rep (ρ' g) * ι (uniformizerZ ϕ π' (deg g))) β)
    (huniq : ∀ g, ∃! q : Q, ∀ β ∈ M',
      act g β = br (rep q * ι (uniformizerZ ϕ π' (deg g))) β) :
    ρ = ρ' :=
  reciprocity_eq_of_bracket br hbr deg act rep (ι θ)
    (fun j => ι (((unitsRingAutHom L ϕ) ^ j) θ))
    (fun j => ι (uniformizerZ ϕ π j)) (fun j => ι (uniformizerZ ϕ π' j)) M M'
    hsemi (map_zpow_mul_uniformizerZ ι ϕ π π' θ hθ) hMM' ρ ρ' hρ hρ' huniq

end Concrete

/-! ## 4. 具体層 2 —— `ρ_m` の定義域が `f` に依らない(木の Lubin-Tate 機械)

`ρ_m` は `Gal(K^m/K)`(`K^m = K^ur · L^m_f`、Definition 4.10)の上の写像なので、
「`ρf,m = ρf′,m` が型として立つ」ためには**定義域が一致する**ことが要る。
直前の節点 `LubinTateUniformizerIndependence.lean` が体の一致を出しているので、
固定部分群のレベルへ移すだけである。

★★ここが射程の要である: 定義域は `K^ur ⊔ K(Λ_{f,n})` であって
`K(Λ_{f,n})` **ではない**。後者では `π ≠ π′` のとき体が一致しないので、
`ρf,m` と `ρf′,m` は共通の定義域を持たない。 -/

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**`ρ_m` の定義域は `f`(と素元)に依らない** —— `Gal(K̄/K)` の中で
`K^m = K^ur · K(Λ_{f,n})` を固定する部分群が一致する。

★これがあって初めて `ρf,m = ρf′,m` が型として立つ。
`lubinTateLevelField_sup_unramifiedClosure_eq_of_isUnit_mul`(直前の節点)を
`IntermediateField.fixingSubgroup` で送っただけ。 -/
theorem fixingSubgroup_lubinTateLevelField_sup_unramifiedClosure_eq_of_isUnit_mul
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u)
    (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {u * π})
    (hϖne0 : u * π ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = u * π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    (lubinTateLevelField K hq hπmax hπne0 f hf0 hf1 hf n ⊔ unramifiedClosure K).fixingSubgroup
      = (lubinTateLevelField K hq hϖmax hϖne0 g hg0 hg1 hg n
          ⊔ unramifiedClosure K).fixingSubgroup :=
  congrArg IntermediateField.fixingSubgroup
    (lubinTateLevelField_sup_unramifiedClosure_eq_of_isUnit_mul K hq hπmax hπne0 f hf0 hf1 hf
      hu hϖmax hϖne0 g hg0 hg1 hg n)

/-- ★★**`ρ`(段数を渡った極限)の定義域も `f` に依らない** ——
`K^LT = K^ur · K_π` を固定する部分群が一致する。

原典 Corollary 4.9 の "hence also K^LT_f and ρ_f" に当たる。 -/
theorem fixingSubgroup_lubinTateClosure_sup_unramifiedClosure_eq_of_isUnit_mul
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u)
    (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {u * π})
    (hϖne0 : u * π ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = u * π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff)) :
    (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K).fixingSubgroup
      = (lubinTateClosure K hq hϖmax hϖne0 g hg0 hg1 hg
          ⊔ unramifiedClosure K).fixingSubgroup :=
  congrArg IntermediateField.fixingSubgroup
    (lubinTateClosure_sup_unramifiedClosure_eq_of_isUnit_mul K hq hπmax hπne0 f hf0 hf1 hf
      hu hϖmax hϖne0 g hg0 hg1 hg)

end ABC3.Found.PGC
