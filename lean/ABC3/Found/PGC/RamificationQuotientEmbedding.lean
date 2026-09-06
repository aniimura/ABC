import ABC3.Found.PGC.LowerRamificationGroup

/-!
# 分岐群の商の埋め込み `G_n/G_{n+1} ↪ 𝓀ˣ` / `𝓀⁺`(Yoshida 2008 Proposition 6.2)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Proposition 6.2**(物理 p.14)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
`section-6.html` の `id="prop-6-2"`。原文の主張は

* `θ_0 : G_0/G_1 ∋ σ ↦ σ(π)/π mod 𝔭 ∈ (𝒪/𝔭)^× ≅ 𝔽_q^×`
* `θ_n : G_n/G_{n+1} ∋ σ ↦ (σ(π)/π) − 1 mod 𝔭^{n+1} ∈ 𝔭^n/𝔭^{n+1} ≅ 𝔽_q`(`n ≥ 1`)

の 2 本が**単射群準同型**であり、しかも**π の取り方に依らない**というもの。

★★**添字の場合分けは本質**である。`n = 0` の行き先は**乗法群** `𝓀ˣ`、`n ≥ 1` の行き先は
**加法群** `𝓀⁺` であり、一様に `↪ 𝓀` と書くと偽になる。Lean の証明で場所が特定できる:
`ramCoeff_mul` が出す補正因子 `(1 + α^n c_σ)^{n+1}` は `n ≥ 1` のときだけ `≡ 1 (mod 𝔪)` で
(`residue_ramCoeff_mul` の `zero_pow` がそこ)、`n = 0` ではこの因子が `σ(π)/π` そのもの
——つまり加法性が壊れる分だけがちょうど乗法群の構造になっている。

## 構成

`α` を素元、`c_σ := (σα − α)/α^{n+1}`(`ramCoeff`)と置くと

* `n = 0`: `σα = α · (1 + c_σ)` で `1 + c_σ` は単元。`θ_0 σ := (1 + c_σ) mod 𝔪 ∈ 𝓀ˣ`。
* `n ≥ 1`: `θ_n σ := c_σ mod 𝔪 ∈ 𝓀⁺`。
* Yoshida 正規化: `Θ_n σ := α^n c_σ mod 𝔪^{n+1} ∈ 𝔭^n/𝔭^{n+1}`(`n ≥ 1`)。

核はどちらも `G_{n+1}` であり(`thetaMul_ker` / `thetaAdd_ker` / `thetaYoshida_ker`)、
`QuotientGroup.lift` で商からの単射準同型になる。

## ★★逸脱の記録(1): 「π の取り方に依らない」が成り立つのは Yoshida 正規化の方だけ

原文は θ_0 と θ_n を並べて "independent of the choice of π" と書く。形式化して分かったのは
**この一言が 2 つの正規化で意味が違う**ことである。`β = α w`(`w` 単元)と取り替えると

* `ramCoeff_independent`: `α^n c^α_σ − β^n c^β_σ ∈ 𝔪^{n+1}`(**全ての `n ≥ 0`**)。
  すなわち原文どおり `𝔭^n/𝔭^{n+1}` に値を取る `Θ_n` は**厳密に π 非依存**
  (`thetaYoshida_independent`)。`n = 0` では `θ_0` の π 非依存性そのもの
  (`thetaMul_independent`)。
* ところが `𝓀` に値を取る `θ_n`(= `α^{n+1}` で割る Serre 型の正規化、本ファイルの `thetaAdd`)は
  `n ≥ 1` では**一般に π 非依存ではない**: `residue_ramCoeff_change_uniformizer` が
  `c̄^α_σ = w̄^n · c̄^β_σ` を与える。`𝔭^n/𝔭^{n+1} ≅ 𝔽_q` という**同一視自身が π 依存**
  (`π^n` を基底に取る)なので、原文の `≅` を抜けた先では独立性は `n = 0` でしか残らない。

⇒ 本ファイルは**両方**を作り、どちらがどこまで π 非依存かを定理として記録する。
`θ_n` を使う後続(Corollary 6.3、Proposition 6.6 (Sen))は「単射」と「加法性」しか
使わないので、この違いは消費側に影響しない。

## 逸脱の記録(2): 群の書き方

* `𝓀⁺` を行き先とする**群**準同型は Lean では `→* Multiplicative 𝓀` と書く。
  `Multiplicative (ResidueField B)` は剰余体の加法群を乗法的に見たもので、
  `Multiplicative.ofAdd` は恒等写像である。
* 商 `G_n/G_{n+1}` は `↥G_n ⧸ (G_{n+1}).subgroupOf G_n` と書く。正規性は
  `LowerRamificationGroup.lean` の `lowerRamificationGroup_normal` と mathlib の
  `Subgroup.normal_subgroupOf` から自動で出る(`subgroupOf_lowerRamificationGroup` が
  `rfl` なので、この商は「`G_n` の中で計算した `(G_n)_{n+1}` による商」と同じ対象)。
* `𝔭^n/𝔭^{n+1}` は `B ⧸ 𝔪^{n+1}` の**中**で実現した(部分商を別の型にしない)。
  像が `𝔭^n/𝔭^{n+1}` に入ることは `thetaYoshida_mem_pow` が述べる。
* `(𝒪/𝔭)^× ≅ 𝔽_q^×` と `𝔭^n/𝔭^{n+1} ≅ 𝔽_q` の**同一視は行っていない**
  (剰余体は `IsLocalRing.ResidueField` のまま)。原文の `𝔽_q` は剰余体の別名にすぎず、
  後続で必要になるのは位数だけである。

## 逸脱の記録(3): 仮定

抽象核は「`B` は局所整域、`𝔪_B = (α)`、`α ≠ 0`、`B = A[α]`」だけを使う
(原典の `B` は完備離散付値環)。`α ≠ 0` は DVR なら自動(`Irreducible.ne_zero`)で、
具体層 §7 ではそう供給している。`B = A[α]` は Lemma 5.11
(`adjoin_uniformizer_eq_top_adjoinIntegers`)の**結論**として渡すので、
完全分岐性は仮定に残らない。

★**結論に自由なパラメータ `α` を残さない**: 抽象層では `α` は写像を書くために必要な
パラメータだが、(a) 具体層 §7 の定理はすべて `α` を `∃` の内側に閉じ込め、
(b) `thetaMul_independent` / `thetaYoshida_independent` が「別の素元を取っても
同じ写像」を保証する。

## ★退化の自己検査

* **`i = 0` を加法群 `𝓀⁺` に入れると偽**。`G_0/G_1` は位数が `q − 1` を割る巡回群であり、
  `𝓀⁺` は指数 `p` の初等アーベル群だから、`q > 2` のときは埋め込めない。
  本ファイルは `n = 0` を `thetaMul`(行き先 `𝓀ˣ`)、`n ≥ 1` を `thetaAdd`(行き先
  `Multiplicative 𝓀`)に**分けており**、`thetaAdd` は `hn : 1 ≤ n` を要求する。
  `hn` を落とすと `residue_ramCoeff_mul` の `zero_pow` が壊れる(証明が通らない)。
* **完全分岐を落とすと `G_0 = ⊤` が偽**になる(`G_0` は惰性群になる)。`thetaMulTop`
  (`G/G_1 ↪ 𝓀ˣ` の形)はこれを `h0 : G_0 = ⊤` として明示的に要求し、具体層では
  `lowerRamificationGroupAdjoin_zero_eq_top`(`IsTotallyRamifiedAdjoin` が要る)から供給する。
  `thetaMul` 自身は `G_0` の上でしか定義していないので、完全分岐でなくても正しい。
* **単射性を落とすと `θ ≡ 1` が通ってしまう**(★退化検査の錨)。これを排除しているのは
  核の計算 `thetaMul_ker` / `thetaAdd_ker` / `thetaYoshida_ker` であり、その `←` 向きは
  `mem_lowerRamificationGroup_iff` を通じて **Lemma 5.11(`B = A[α]`)を使う**。
  すなわち `hadj` を落とすと単射性だけが落ちる。
* **`hα0 : α ≠ 0` を落とすと `divOfDvd` の一意性が壊れる**(`B` が体のとき `𝔪 = (0)`)。
  そのとき `G_n` は全て自明なので主張自体は無害だが、証明が通らない。

## 実測(抽象核と具体層の分離)

前ノードと同じく**中核を「一般の可換環・一般の群作用」に切り出す**方針が効いた。
§1–§6(抽象核)の `lean_check` は 0.05–0.52 秒、§7(`PAdicLocalField`)は 4 本まとめて
2.5 秒。`PAdicLocalField` のインスタンス探索は §7 の 5 宣言だけに閉じている。
-/

namespace ABC3.Found.PGC

open IsLocalRing ABC3.Skeleton.PGC
open scoped NNReal Valued

def thetaMul.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 14, item := "Proposition 6.2", sectionId := "prop-6-2" }

/-! ## §1 整域での「割り算」

`σα − α ∈ 𝔪^{n+1} = (α^{n+1})` から商 `c_σ` を取り出すための最小限の道具。
整域なので商は一意で、`divOfDvd_mul_cancel` がそれを述べる。 -/

open scoped Classical in
/-- 「割り切れるときの商」を取る全域関数。`a ∣ y` のときだけ意味を持ち、
整域では `a ≠ 0` の下で一意(`divOfDvd_mul_cancel`)。 -/
noncomputable def divOfDvd {B : Type*} [CommRing B] (a y : B) : B :=
  if h : a ∣ y then h.choose else 0

theorem mul_divOfDvd {B : Type*} [CommRing B] {a y : B} (h : a ∣ y) :
    a * divOfDvd a y = y := by
  rw [divOfDvd, dif_pos h]
  exact h.choose_spec.symm

theorem divOfDvd_mul_cancel {B : Type*} [CommRing B] [IsDomain B] {a : B} (ha : a ≠ 0) (c : B) :
    divOfDvd a (a * c) = c :=
  mul_left_cancel₀ ha (mul_divOfDvd ⟨c, rfl⟩)

theorem divOfDvd_zero {B : Type*} [CommRing B] [IsDomain B] {a : B} (ha : a ≠ 0) :
    divOfDvd a 0 = 0 := by
  simpa using divOfDvd_mul_cancel ha 0

/-! ## §2 分岐係数 `c_σ = (σα − α)/α^{n+1}`

原典 `Prop 6.2` の写像の中身。`σ ∈ G_n` の定義そのもの(`σα − α ∈ 𝔪^{n+1}`)から
割り算ができ、`ramCoeff_mul` が「準同型からのずれ」を完全に記述する。 -/

/-- **分岐係数** `c_σ := (σα − α)/α^{n+1}`。`σ ∈ G_n` のときに意味を持つ
(`pow_mul_ramCoeff`)。 -/
noncomputable def ramCoeff {B : Type*} [CommRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (α : B) (n : ℕ) (σ : G) : B := divOfDvd (α ^ (n + 1)) (σ • α - α)

/-- `σ ∈ G_n` なら `α^{n+1} c_σ = σα − α`。★ここで使うのは `G_n` の定義(`Ideal.inertia`)を
`x := α` に当てるだけで、Lemma 5.11 は要らない。 -/
theorem pow_mul_ramCoeff {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) {n : ℕ} {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G n) :
    α ^ (n + 1) * ramCoeff α n σ = σ • α - α := by
  apply mul_divOfDvd
  rw [← Ideal.mem_span_singleton, ← Ideal.span_singleton_pow, ← hα]
  exact mem_lowerRamificationGroup_iff_forall.mp hσ α

/-- `σα = α + α^{n+1} c_σ`(`pow_mul_ramCoeff` の移項)。 -/
theorem smul_eq_add_pow_mul_ramCoeff {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) {n : ℕ} {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G n) :
    σ • α = α + α ^ (n + 1) * ramCoeff α n σ := by
  rw [pow_mul_ramCoeff hα hσ]; ring

theorem ramCoeff_one {B : Type*} [CommRing B] [IsDomain B] {G : Type*} [Group G]
    [MulSemiringAction G B] {α : B} (hα0 : α ≠ 0) (n : ℕ) :
    ramCoeff α n (1 : G) = 0 := by
  rw [ramCoeff, one_smul, sub_self]
  exact divOfDvd_zero (pow_ne_zero _ hα0)

/-- **合成則(厳密形)**——`c_{στ} = c_σ + (1 + α^n c_σ)^{n+1} · σ(c_τ)`。

★この 1 本が `n = 0` と `n ≥ 1` の分かれ目を全部持っている。補正因子
`(1 + α^n c_σ)^{n+1}` は `n ≥ 1` なら `α^n ∈ 𝔪` なので法 `𝔪` で `1` になり
`θ_n` が加法的になる(`residue_ramCoeff_mul`)。`n = 0` ではこの因子は
`(1 + c_σ)^1 = σα/α` そのもので、消えずに残る——そこが乗法群の構造になる
(`one_add_ramCoeff_zero_mul`)。 -/
theorem ramCoeff_mul {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    {n : ℕ} {σ τ : G} (hσ : σ ∈ lowerRamificationGroup B G n)
    (hτ : τ ∈ lowerRamificationGroup B G n) :
    ramCoeff α n (σ * τ)
      = ramCoeff α n σ + (1 + α ^ n * ramCoeff α n σ) ^ (n + 1) * (σ • ramCoeff α n τ) := by
  have hστ := (lowerRamificationGroup B G n).mul_mem hσ hτ
  apply mul_left_cancel₀ (pow_ne_zero (n + 1) hα0)
  rw [pow_mul_ramCoeff hα hστ, mul_smul, smul_eq_add_pow_mul_ramCoeff hα hτ,
    smul_add, smul_mul', smul_pow', smul_eq_add_pow_mul_ramCoeff hα hσ]
  have h2 : α + α ^ (n + 1) * ramCoeff α n σ = α * (1 + α ^ n * ramCoeff α n σ) := by ring
  rw [h2, mul_pow]
  ring

/-! ### 剰余体への降下に要る小道具 -/

theorem sub_mem_maximalIdeal_of_mem_zero {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {σ : G} (hσ : σ ∈ lowerRamificationGroup B G 0) (x : B) :
    σ • x - x ∈ maximalIdeal B := by
  simpa using mem_lowerRamificationGroup_iff_forall.mp hσ x

/-- `G_0` の元は剰余体に自明に作用する(`G_0` が惰性群であることの言い換え)。 -/
theorem residue_smul_of_mem_zero {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {σ : G} (hσ : σ ∈ lowerRamificationGroup B G 0) (x : B) :
    residue B (σ • x) = residue B x := by
  rw [← sub_eq_zero, ← map_sub]
  exact (residue_eq_zero_iff _).mpr (sub_mem_maximalIdeal_of_mem_zero hσ x)

theorem mem_zero_of_mem {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {n : ℕ} {σ : G} (hσ : σ ∈ lowerRamificationGroup B G n) :
    σ ∈ lowerRamificationGroup B G 0 :=
  lowerRamificationGroup_antitone B G (Nat.zero_le n) hσ

theorem residue_uniformizer {B : Type*} [CommRing B] [IsLocalRing B] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) : residue B α = 0 :=
  (residue_eq_zero_iff _).mpr (hα ▸ Ideal.mem_span_singleton_self α)

/-- `α^k` を掛けても「`(α)` に入るか」は変わらない(整域での約分)。 -/
theorem mem_span_singleton_iff_pow_mul {B : Type*} [CommRing B] [IsDomain B] {α : B}
    (hα0 : α ≠ 0) (k : ℕ) (c : B) :
    c ∈ Ideal.span ({α} : Set B) ↔ α ^ k * c ∈ Ideal.span ({α ^ (k + 1)} : Set B) := by
  simp only [Ideal.mem_span_singleton']
  constructor
  · rintro ⟨d, rfl⟩
    exact ⟨d, by ring⟩
  · rintro ⟨d, hd⟩
    refine ⟨d, ?_⟩
    apply mul_left_cancel₀ (pow_ne_zero k hα0)
    rw [← hd]; ring

theorem pow_mul_mem_pow_succ {B : Type*} [CommRing B] [IsLocalRing B] {α : B}
    (hαm : α ∈ maximalIdeal B) (n : ℕ) {z : B} (hz : z ∈ maximalIdeal B) :
    α ^ n * z ∈ (maximalIdeal B) ^ (n + 1) := by
  rw [pow_succ]
  exact Ideal.mul_mem_mul (Ideal.pow_mem_pow hαm n) hz

theorem pow_mul_mem_pow_succ_iff {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0) (n : ℕ) (c : B) :
    α ^ n * c ∈ (maximalIdeal B) ^ (n + 1) ↔ c ∈ maximalIdeal B := by
  rw [hα, Ideal.span_singleton_pow, ← mem_span_singleton_iff_pow_mul hα0 n]

/-- **核の判定**——`c_σ ≡ 0 (mod 𝔪) ⟺ σ ∈ G_{n+1}`。

★`←` 向きに **Lemma 5.11(`B = A[α]`)**が要る(`mem_lowerRamificationGroup_iff`)。
これが「単射性」を支えており、`hadj` を落とすと退化写像 `θ ≡ 1` が排除できなくなる。 -/
theorem residue_ramCoeff_eq_zero_iff {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℕ} {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G n) :
    residue B (ramCoeff α n σ) = 0 ↔ σ ∈ lowerRamificationGroup B G (n + 1) := by
  rw [residue_eq_zero_iff, hα, mem_span_singleton_iff_pow_mul hα0 (n + 1),
    pow_mul_ramCoeff hα hσ, mem_lowerRamificationGroup_iff hadj, hα, Ideal.span_singleton_pow]

/-! ## §3 `θ_0 : G_0/G_1 ↪ 𝓀ˣ`(Prop 6.2 の前半)

原文 `θ_0 : G_0/G_1 ∋ σ ⟼ σ(π)/π mod 𝔭_{K'} ∈ (𝒪_{K'}/𝔭_{K'})^×`。
`σ(π)/π = 1 + c_σ` であり(`mul_one_add_ramCoeff_zero`)、これが単元であることは
`σ` と `σ⁻¹` の関係式から出る。 -/

/-- 原文の `σ(π)/π` そのもの: `α · (1 + c_σ) = σα`。 -/
theorem mul_one_add_ramCoeff_zero {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 0) :
    α * (1 + ramCoeff α 0 σ) = σ • α := by
  linear_combination pow_mul_ramCoeff hα hσ

/-- `σ(π)/π` は単元。`σ` を `α·(1 + c_{σ⁻¹}) = σ⁻¹α` に施して `α` で約分すると
`(1 + c_σ) · σ(1 + c_{σ⁻¹}) = 1`。 -/
theorem isUnit_one_add_ramCoeff_zero {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 0) : IsUnit (1 + ramCoeff α 0 σ) := by
  have hinv : σ⁻¹ ∈ lowerRamificationGroup B G 0 := (lowerRamificationGroup B G 0).inv_mem hσ
  have h1 : α * (1 + ramCoeff α 0 σ) = σ • α := mul_one_add_ramCoeff_zero hα hσ
  have h2 : α * (1 + ramCoeff α 0 σ⁻¹) = σ⁻¹ • α := mul_one_add_ramCoeff_zero hα hinv
  have h3 : (σ • α) * (σ • (1 + ramCoeff α 0 σ⁻¹)) = α := by
    have h := congrArg (fun z : B => σ • z) h2
    simpa [smul_mul', smul_smul] using h
  refine IsUnit.of_mul_eq_one (σ • (1 + ramCoeff α 0 σ⁻¹)) ?_
  apply mul_left_cancel₀ hα0
  rw [mul_one, ← mul_assoc, h1, h3]

/-- **`n = 0` での合成則**——`1 + c_{στ} = (1 + c_σ) · σ(1 + c_τ)`。
`ramCoeff_mul` の補正因子 `(1 + α^0 c_σ)^{0+1}` がちょうど `1 + c_σ` になる。 -/
theorem one_add_ramCoeff_zero_mul {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0) {σ τ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 0) (hτ : τ ∈ lowerRamificationGroup B G 0) :
    1 + ramCoeff α 0 (σ * τ) = (1 + ramCoeff α 0 σ) * (σ • (1 + ramCoeff α 0 τ)) := by
  rw [ramCoeff_mul hα hα0 hσ hτ, smul_add, smul_one]
  ring

/-- **Yoshida Prop 6.2 (θ_0)**——`G_0 → 𝓀ˣ`、`σ ↦ σ(π)/π mod 𝔭`。

★行き先は**乗法群** `𝓀ˣ`。加法群にすると偽(ファイル冒頭「退化の自己検査」)。
準同型性は `one_add_ramCoeff_zero_mul` と「`G_0` は剰余体に自明に作用する」の 2 つから。 -/
noncomputable def thetaMul {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α})
    (hα0 : α ≠ 0) :
    (lowerRamificationGroup B G 0) →* (ResidueField B)ˣ where
  toFun σ := ((isUnit_one_add_ramCoeff_zero hα hα0 σ.2).map (residue B)).unit
  map_one' := by
    apply Units.ext
    simp [ramCoeff_one hα0]
  map_mul' σ τ := by
    apply Units.ext
    push_cast [IsUnit.unit_spec]
    rw [one_add_ramCoeff_zero_mul hα hα0 σ.2 τ.2, map_mul, residue_smul_of_mem_zero σ.2]

@[simp] theorem thetaMul_apply_coe {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (σ : lowerRamificationGroup B G 0) :
    ((thetaMul hα hα0 σ : (ResidueField B)ˣ) : ResidueField B)
      = residue B (1 + ramCoeff α 0 (σ : G)) :=
  IsUnit.unit_spec _

/-- **`ker θ_0 = G_1`**。★これが単射性の全部であり、Lemma 5.11 を使う唯一の場所。 -/
theorem thetaMul_ker {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) :
    (thetaMul (G := G) hα hα0).ker
      = (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0) := by
  ext σ
  rw [MonoidHom.mem_ker, Units.ext_iff, thetaMul_apply_coe, Units.val_one, map_add, map_one,
    add_eq_left, Subgroup.mem_subgroupOf]
  exact residue_ramCoeff_eq_zero_iff hα hα0 hadj σ.2

/-- **Yoshida Prop 6.2 (θ_0、商の形)**: `G_0/G_1 → 𝓀ˣ`。 -/
noncomputable def thetaMulQuot {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) :
    (lowerRamificationGroup B G 0 ⧸
        (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0))
      →* (ResidueField B)ˣ :=
  QuotientGroup.lift _ (thetaMul hα hα0)
    (fun x hx => by rw [← thetaMul_ker hα hα0 hadj] at hx; exact hx)

/-- **`θ_0 : G_0/G_1 ↪ 𝓀ˣ` は単射**(Prop 6.2 前半)。 -/
theorem thetaMulQuot_injective {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) :
    Function.Injective (thetaMulQuot (G := G) hα hα0 hadj) := by
  rw [injective_iff_map_eq_one]
  intro q hq
  obtain ⟨σ, rfl⟩ := QuotientGroup.mk_surjective q
  rw [QuotientGroup.eq_one_iff, ← thetaMul_ker hα hα0 hadj, MonoidHom.mem_ker]
  exact hq

/-! ### 完全分岐のとき(`G_0 = ⊤`)は `G/G_1 ↪ 𝓀ˣ`

原典は `Def 6.1` の直後に「`K'/K` は完全分岐だから `G = G_0`」と書いており、
Prop 6.2 の θ_0 は事実上 `G/G_1 ↪ 𝓀ˣ` として使われる。 -/

/-- `G_0 = ⊤` のときの包含 `G →* G_0`。 -/
def toLowerRamificationGroupZero {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] (h0 : lowerRamificationGroup B G 0 = ⊤) :
    G →* (lowerRamificationGroup B G 0) where
  toFun σ := ⟨σ, by rw [h0]; exact Subgroup.mem_top σ⟩
  map_one' := rfl
  map_mul' _ _ := rfl

theorem ker_thetaMul_comp {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (h0 : lowerRamificationGroup B G 0 = ⊤) :
    ((thetaMul hα hα0).comp (toLowerRamificationGroupZero h0)).ker
      = lowerRamificationGroup B G 1 := by
  ext σ
  rw [MonoidHom.mem_ker, MonoidHom.comp_apply, ← MonoidHom.mem_ker, thetaMul_ker hα hα0 hadj,
    Subgroup.mem_subgroupOf]
  exact Iff.rfl

/-- **完全分岐版の θ_0**: `G/G_1 → 𝓀ˣ`。 -/
noncomputable def thetaMulTop {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (h0 : lowerRamificationGroup B G 0 = ⊤) :
    (G ⧸ lowerRamificationGroup B G 1) →* (ResidueField B)ˣ :=
  QuotientGroup.lift _ ((thetaMul hα hα0).comp (toLowerRamificationGroupZero h0))
    (le_of_eq (ker_thetaMul_comp hα hα0 hadj h0).symm)

theorem thetaMulTop_injective {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (h0 : lowerRamificationGroup B G 0 = ⊤) :
    Function.Injective (thetaMulTop hα hα0 hadj h0) := by
  rw [injective_iff_map_eq_one]
  intro q hq
  obtain ⟨σ, rfl⟩ := QuotientGroup.mk_surjective q
  rw [QuotientGroup.eq_one_iff, ← ker_thetaMul_comp hα hα0 hadj h0, MonoidHom.mem_ker]
  exact hq

/-! ## §4 `θ_n : G_n/G_{n+1} ↪ 𝓀⁺`(`n ≥ 1`、Prop 6.2 の後半)

`𝔭^n/𝔭^{n+1} ≅ 𝔽_q` を `π^n` で割って実現した版(Serre 型の正規化)。
★この正規化は π 非依存**ではない**(§5 の `residue_ramCoeff_change_uniformizer`)。 -/

/-- **加法性**——`n ≥ 1` なら `c̄_{στ} = c̄_σ + c̄_τ`。

★`hn : 1 ≤ n` を使うのはただ 1 箇所、`residue (α^n) = 0` を出す `zero_pow` である。
`n = 0` だと `α^0 = 1` で `residue = 1 ≠ 0` になり、補正因子が消えない。 -/
theorem residue_ramCoeff_mul {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α})
    (hα0 : α ≠ 0) {n : ℕ} (hn : 1 ≤ n) {σ τ : G} (hσ : σ ∈ lowerRamificationGroup B G n)
    (hτ : τ ∈ lowerRamificationGroup B G n) :
    residue B (ramCoeff α n (σ * τ))
      = residue B (ramCoeff α n σ) + residue B (ramCoeff α n τ) := by
  have hz : residue B (α ^ n) = 0 := by
    rw [map_pow, residue_uniformizer hα, zero_pow (by omega)]
  rw [ramCoeff_mul hα hα0 hσ hτ, map_add, map_mul, map_pow, map_add, map_one, map_mul, hz,
    zero_mul, add_zero, one_pow, one_mul, residue_smul_of_mem_zero (mem_zero_of_mem hσ)]

/-- **Yoshida Prop 6.2 (θ_n、`n ≥ 1`)**——`G_n → 𝓀⁺`、`σ ↦ (σα − α)/α^{n+1} mod 𝔪`。

★行き先の `Multiplicative (ResidueField B)` は剰余体の**加法群**を乗法的に見たもの
(Lean で「加法群への群準同型」を書く標準の書き方)。 -/
noncomputable def thetaAdd {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α})
    (hα0 : α ≠ 0) {n : ℕ} (hn : 1 ≤ n) :
    (lowerRamificationGroup B G n) →* Multiplicative (ResidueField B) where
  toFun σ := Multiplicative.ofAdd (residue B (ramCoeff α n (σ : G)))
  map_one' := by simp [ramCoeff_one hα0]
  map_mul' σ τ := by
    simp only [Subgroup.coe_mul]
    rw [residue_ramCoeff_mul hα hα0 hn σ.2 τ.2]
    rfl

@[simp] theorem thetaAdd_apply {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α})
    (hα0 : α ≠ 0) {n : ℕ} (hn : 1 ≤ n) (σ : lowerRamificationGroup B G n) :
    thetaAdd hα hα0 hn σ = Multiplicative.ofAdd (residue B (ramCoeff α n (σ : G))) := rfl

/-- **`ker θ_n = G_{n+1}`**(`n ≥ 1`)。 -/
theorem thetaAdd_ker {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℕ} (hn : 1 ≤ n) :
    (thetaAdd (G := G) hα hα0 hn).ker
      = (lowerRamificationGroup B G (n + 1)).subgroupOf (lowerRamificationGroup B G n) := by
  ext σ
  simp only [MonoidHom.mem_ker, Subgroup.mem_subgroupOf, thetaAdd_apply, ofAdd_eq_one]
  exact residue_ramCoeff_eq_zero_iff hα hα0 hadj σ.2

/-- **Yoshida Prop 6.2 (θ_n、商の形)**: `G_n/G_{n+1} → 𝓀⁺`(`n ≥ 1`)。 -/
noncomputable def thetaAddQuot {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℕ} (hn : 1 ≤ n) :
    (lowerRamificationGroup B G n ⧸
        (lowerRamificationGroup B G (n + 1)).subgroupOf (lowerRamificationGroup B G n))
      →* Multiplicative (ResidueField B) :=
  QuotientGroup.lift _ (thetaAdd hα hα0 hn)
    (fun x hx => by rw [← thetaAdd_ker hα hα0 hadj hn] at hx; exact hx)

/-- **`θ_n : G_n/G_{n+1} ↪ 𝓀⁺` は単射**(`n ≥ 1`、Prop 6.2 後半)。 -/
theorem thetaAddQuot_injective {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℕ} (hn : 1 ≤ n) :
    Function.Injective (thetaAddQuot (G := G) hα hα0 hadj hn) := by
  rw [injective_iff_map_eq_one]
  intro q hq
  obtain ⟨σ, rfl⟩ := QuotientGroup.mk_surjective q
  rw [QuotientGroup.eq_one_iff, ← thetaAdd_ker hα hα0 hadj hn, MonoidHom.mem_ker]
  exact hq

/-! ## §5 素元の取り替え——原文の "independent of the choice of π"

`β = α w`(`w` は単元)と取り替えたときの比較。★中心になるのは
`ramCoeff_independent`(`α^n c^α_σ ≡ β^n c^β_σ mod 𝔪^{n+1}`)で、これは
**原文の正規化** `(σπ/π) − 1 mod 𝔭^{n+1} ∈ 𝔭^n/𝔭^{n+1}` がそのまま π 非依存だという主張。
`𝓀` に落とした `θ_n` の方は `w̄^n` だけねじれる(冒頭「逸脱の記録(1)」)。 -/

/-- **素元非依存性(原文の正規化)**——`α^n c^α_σ − β^n c^β_σ ∈ 𝔪^{n+1}`。全ての `n ≥ 0`。

証明は `σβ = (σα)(σw)` を展開して `α^{n+1}` で約分し、`σw − w ∈ 𝔪^{n+1}` を使うだけ。
★`G_n` の定義が `Ideal.inertia`(π 非依存)で与えられているおかげで、
「`σ ∈ G_n` を仮定してよい」という部分が**無料**になっている
(原典は `G_n` の π 非依存性をここで別途処理する必要がある)。 -/
theorem ramCoeff_independent {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α β : B} (hα : maximalIdeal B = Ideal.span {α})
    (hα0 : α ≠ 0) (hβ : maximalIdeal B = Ideal.span {β}) {n : ℕ} {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G n) :
    α ^ n * ramCoeff α n σ - β ^ n * ramCoeff β n σ ∈ (maximalIdeal B) ^ (n + 1) := by
  obtain ⟨w, hw⟩ : Associated α β := Ideal.span_singleton_eq_span_singleton.mp (hα ▸ hβ)
  obtain ⟨d, hd⟩ : ∃ d : B, σ • (w : B) - (w : B) = α ^ (n + 1) * d := by
    have h := mem_lowerRamificationGroup_iff_forall.mp hσ (w : B)
    rw [hα, Ideal.span_singleton_pow, Ideal.mem_span_singleton] at h
    exact h
  have hX : σ • α = α + α ^ (n + 1) * ramCoeff α n σ := smul_eq_add_pow_mul_ramCoeff hα hσ
  have hZ : σ • (w : B) = (w : B) + α ^ (n + 1) * d := by rw [← hd]; ring
  have hsβ : σ • β = (σ • α) * (σ • (w : B)) := by rw [← hw, smul_mul']
  have hpow : β ^ (n + 1) = α ^ (n + 1) * (w : B) ^ (n + 1) := by rw [← hw, mul_pow]
  have hpowN : β ^ n = α ^ n * (w : B) ^ n := by rw [← hw, mul_pow]
  have hstar : (w : B) ^ (n + 1) * ramCoeff β n σ
      = α * d + ramCoeff α n σ * (w : B) + α ^ (n + 1) * ramCoeff α n σ * d := by
    apply mul_left_cancel₀ (pow_ne_zero (n + 1) hα0)
    have e1 : α ^ (n + 1) * ((w : B) ^ (n + 1) * ramCoeff β n σ)
        = β ^ (n + 1) * ramCoeff β n σ := by rw [hpow]; ring
    rw [e1, pow_mul_ramCoeff hβ hσ, hsβ, hX, hZ, ← hw]
    ring
  have hgoal : (w : B) * (α ^ n * ramCoeff α n σ - β ^ n * ramCoeff β n σ)
      = α ^ (n + 1) * (-(d + α ^ n * ramCoeff α n σ * d)) := by
    rw [hpowN]
    linear_combination (-(α ^ n)) * hstar
  rw [hα, Ideal.span_singleton_pow, ← Ideal.unit_mul_mem_iff_mem _ w.isUnit, hgoal,
    Ideal.mem_span_singleton]
  exact ⟨_, rfl⟩

/-- **θ_0 は素元の取り方に依らない**(`n = 0` では `𝓀ˣ` への写像そのものが不変)。 -/
theorem thetaMul_independent {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α β : B} (hα : maximalIdeal B = Ideal.span {α})
    (hα0 : α ≠ 0) (hβ : maximalIdeal B = Ideal.span {β}) (hβ0 : β ≠ 0) :
    thetaMul (G := G) hα hα0 = thetaMul (G := G) hβ hβ0 := by
  refine MonoidHom.ext fun σ => Units.ext ?_
  rw [thetaMul_apply_coe, thetaMul_apply_coe, ← sub_eq_zero, ← map_sub, residue_eq_zero_iff]
  have he : (1 + ramCoeff α 0 (σ : G)) - (1 + ramCoeff β 0 (σ : G))
      = α ^ 0 * ramCoeff α 0 (σ : G) - β ^ 0 * ramCoeff β 0 (σ : G) := by ring
  rw [he]
  simpa using ramCoeff_independent hα hα0 hβ (n := 0) σ.2

/-- ★★**`𝓀` に値を取る `θ_n` は `n ≥ 1` では一般に素元非依存ではない**——
`β = α w` と取り替えると `c̄^α_σ = w̄^n · c̄^β_σ` とねじれる。

これは原文の `𝔭^n/𝔭^{n+1} ≅ 𝔽_q` という同一視が `π^n` の選択を含むためで、
矛盾ではない(冒頭「逸脱の記録(1)」)。`n = 0` では `w̄^0 = 1` で
`thetaMul_independent` に一致する。 -/
theorem residue_ramCoeff_change_uniformizer {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {α β : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hβ : maximalIdeal B = Ideal.span {β}) {n : ℕ} {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G n) {w : Bˣ} (hw : α * (w : B) = β) :
    residue B (ramCoeff α n σ)
      = (residue B (w : B)) ^ n * residue B (ramCoeff β n σ) := by
  have hpowN : β ^ n = α ^ n * (w : B) ^ n := by rw [← hw, mul_pow]
  have h2 : α ^ n * (ramCoeff α n σ - (w : B) ^ n * ramCoeff β n σ)
      ∈ (maximalIdeal B) ^ (n + 1) := by
    have he : α ^ n * (ramCoeff α n σ - (w : B) ^ n * ramCoeff β n σ)
        = α ^ n * ramCoeff α n σ - β ^ n * ramCoeff β n σ := by rw [hpowN]; ring
    rw [he]
    exact ramCoeff_independent hα hα0 hβ hσ
  have h3 := (residue_eq_zero_iff _).mpr ((pow_mul_mem_pow_succ_iff hα hα0 n _).mp h2)
  rw [map_sub, sub_eq_zero, map_mul, map_pow] at h3
  exact h3

/-! ### Yoshida 正規化 `Θ_n : G_n → 𝔭^n/𝔭^{n+1}`(`n ≥ 1`)

原文の字面どおり `σ ↦ (σπ/π) − 1 mod 𝔭^{n+1}`。`(σπ/π) − 1 = α^n c_σ` である。
こちらが**厳密に π 非依存**な方(`thetaYoshida_independent`)。 -/

theorem ramCoeff_mul_sub_mem {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α})
    (hα0 : α ≠ 0) {n : ℕ} (hn : 1 ≤ n) {σ τ : G} (hσ : σ ∈ lowerRamificationGroup B G n)
    (hτ : τ ∈ lowerRamificationGroup B G n) :
    α ^ n * ramCoeff α n (σ * τ) - (α ^ n * ramCoeff α n σ + α ^ n * ramCoeff α n τ)
      ∈ (maximalIdeal B) ^ (n + 1) := by
  have h : ramCoeff α n (σ * τ) - (ramCoeff α n σ + ramCoeff α n τ) ∈ maximalIdeal B := by
    rw [← residue_eq_zero_iff, map_sub, map_add, residue_ramCoeff_mul hα hα0 hn hσ hτ, sub_self]
  have he : α ^ n * ramCoeff α n (σ * τ) - (α ^ n * ramCoeff α n σ + α ^ n * ramCoeff α n τ)
      = α ^ n * (ramCoeff α n (σ * τ) - (ramCoeff α n σ + ramCoeff α n τ)) := by ring
  rw [he]
  exact pow_mul_mem_pow_succ (hα ▸ Ideal.mem_span_singleton_self α) n h

/-- **Yoshida Prop 6.2 の正規化そのまま**——`Θ_n : G_n → B/𝔭^{n+1}`、
`σ ↦ (σπ/π) − 1 mod 𝔭^{n+1}`(値は `𝔭^n/𝔭^{n+1}` に入る、`thetaYoshida_mem_pow`)。 -/
noncomputable def thetaYoshida {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α})
    (hα0 : α ≠ 0) {n : ℕ} (hn : 1 ≤ n) :
    (lowerRamificationGroup B G n) →* Multiplicative (B ⧸ (maximalIdeal B) ^ (n + 1)) where
  toFun σ := Multiplicative.ofAdd
    (Ideal.Quotient.mk ((maximalIdeal B) ^ (n + 1)) (α ^ n * ramCoeff α n (σ : G)))
  map_one' := by simp [ramCoeff_one hα0]
  map_mul' σ τ := by
    have h : Ideal.Quotient.mk ((maximalIdeal B) ^ (n + 1))
          (α ^ n * ramCoeff α n ((σ : G) * (τ : G)))
        = Ideal.Quotient.mk ((maximalIdeal B) ^ (n + 1)) (α ^ n * ramCoeff α n (σ : G))
          + Ideal.Quotient.mk ((maximalIdeal B) ^ (n + 1)) (α ^ n * ramCoeff α n (τ : G)) := by
      rw [← map_add, Ideal.Quotient.eq]
      exact ramCoeff_mul_sub_mem hα hα0 hn σ.2 τ.2
    exact congrArg Multiplicative.ofAdd h

@[simp] theorem thetaYoshida_apply {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0) {n : ℕ} (hn : 1 ≤ n)
    (σ : lowerRamificationGroup B G n) :
    thetaYoshida hα hα0 hn σ = Multiplicative.ofAdd
      (Ideal.Quotient.mk ((maximalIdeal B) ^ (n + 1)) (α ^ n * ramCoeff α n (σ : G))) := rfl

/-- **像は `𝔭^n/𝔭^{n+1}` に入る**(原文の行き先)。 -/
theorem thetaYoshida_mem_pow {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0) {n : ℕ} (hn : 1 ≤ n)
    (σ : lowerRamificationGroup B G n) :
    Multiplicative.toAdd (thetaYoshida hα hα0 hn σ)
      ∈ Ideal.map (Ideal.Quotient.mk ((maximalIdeal B) ^ (n + 1))) ((maximalIdeal B) ^ n) :=
  Ideal.mem_map_of_mem _
    (Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow (hα ▸ Ideal.mem_span_singleton_self α) n))

theorem thetaYoshida_ker {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℕ} (hn : 1 ≤ n) :
    (thetaYoshida (G := G) hα hα0 hn).ker
      = (lowerRamificationGroup B G (n + 1)).subgroupOf (lowerRamificationGroup B G n) := by
  ext σ
  rw [MonoidHom.mem_ker, Subgroup.mem_subgroupOf, thetaYoshida_apply, ofAdd_eq_one,
    Ideal.Quotient.eq_zero_iff_mem, pow_mul_mem_pow_succ_iff hα hα0, ← residue_eq_zero_iff]
  exact residue_ramCoeff_eq_zero_iff hα hα0 hadj σ.2

/-- ★★**原文の "independent of the choice of π" の厳密形**——
`Θ_n` は素元の取り方に**完全に**依らない(写像として等しい)。 -/
theorem thetaYoshida_independent {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {α β : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hβ : maximalIdeal B = Ideal.span {β}) (hβ0 : β ≠ 0) {n : ℕ} (hn : 1 ≤ n) :
    thetaYoshida (G := G) hα hα0 hn = thetaYoshida (G := G) hβ hβ0 hn := by
  ext σ
  have h : Ideal.Quotient.mk ((maximalIdeal B) ^ (n + 1)) (α ^ n * ramCoeff α n (σ : G))
      = Ideal.Quotient.mk ((maximalIdeal B) ^ (n + 1)) (β ^ n * ramCoeff β n (σ : G)) := by
    rw [Ideal.Quotient.eq]
    exact ramCoeff_independent hα hα0 hβ σ.2
  exact congrArg Multiplicative.ofAdd h

/-- **Yoshida 正規化(商の形)**: `G_n/G_{n+1} → 𝔭^n/𝔭^{n+1} ⊆ B/𝔭^{n+1}`。 -/
noncomputable def thetaYoshidaQuot {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℕ} (hn : 1 ≤ n) :
    (lowerRamificationGroup B G n ⧸
        (lowerRamificationGroup B G (n + 1)).subgroupOf (lowerRamificationGroup B G n))
      →* Multiplicative (B ⧸ (maximalIdeal B) ^ (n + 1)) :=
  QuotientGroup.lift _ (thetaYoshida hα hα0 hn)
    (fun x hx => by rw [← thetaYoshida_ker hα hα0 hadj hn] at hx; exact hx)

theorem thetaYoshidaQuot_injective {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℕ} (hn : 1 ≤ n) :
    Function.Injective (thetaYoshidaQuot (G := G) hα hα0 hadj hn) := by
  rw [injective_iff_map_eq_one]
  intro q hq
  obtain ⟨σ, rfl⟩ := QuotientGroup.mk_surjective q
  rw [QuotientGroup.eq_one_iff, ← thetaYoshida_ker hα hα0 hadj hn, MonoidHom.mem_ker]
  exact hq

/-! ## §6 系: 各段はアーベル(原文の "they show that G is supersoluble")

単射準同型の行き先が可換群なので、`G_n/G_{n+1}` は可換。原文が Prop 6.2 の括弧で
述べている超可解性の内訳はこれ(と `G_n` が正規列をなすこと)である。 -/

theorem mul_comm_of_quot_lowerRamificationGroup {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℕ} (hn : 1 ≤ n)
    (x y : lowerRamificationGroup B G n ⧸
      (lowerRamificationGroup B G (n + 1)).subgroupOf (lowerRamificationGroup B G n)) :
    x * y = y * x :=
  thetaAddQuot_injective hα hα0 hadj hn (by rw [map_mul, map_mul, mul_comm])

theorem mul_comm_of_quot_lowerRamificationGroup_zero {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤)
    (x y : lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0)) :
    x * y = y * x :=
  thetaMulQuot_injective hα hα0 hadj (by rw [map_mul, map_mul, mul_comm])

/-! ## §7 `PAdicLocalField` への具体化

`A := 𝒪[K.carrier]`、`B := adjoinIntegers K x`(`= 𝒪_{K(x)}`)、`G := Gal(K(x)/K)`。
素元 `α` は `∃` の内側に閉じ込め、結論に自由なパラメータを残さない
(取り方に依らないことは §5 が保証する)。 -/

variable {p : ℕ} [Fact p.Prime]

instance lowerRamificationGroupAdjoin_normal (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (n : ℕ) : (lowerRamificationGroupAdjoin K x n).Normal :=
  lowerRamificationGroup_normal _ _ n

/-- 具体層で 4 回使う「素元 + 非零 + 単項生成」の取得(Lemma 5.11)。 -/
theorem exists_uniformizer_ne_zero_adjoin_eq (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) :
    ∃ α : adjoinIntegers K x,
      IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α} ∧ α ≠ 0 ∧
      Algebra.adjoin 𝒪[K.carrier] ({α} : Set (adjoinIntegers K x)) = ⊤ := by
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  obtain ⟨α, hα⟩ := IsDiscreteValuationRing.exists_irreducible (adjoinIntegers K x)
  have hspan : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer α).mp hα
  exact ⟨α, hspan, hα.ne_zero, adjoin_uniformizer_eq_top_adjoinIntegers K x ht hspan⟩

/-- **Yoshida Prop 6.2 前半(`PAdicLocalField` 版)**: `G_0/G_1 ↪ 𝓀ˣ`。 -/
theorem exists_injective_thetaMul_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) :
    ∃ θ : (lowerRamificationGroupAdjoin K x 0 ⧸
        (lowerRamificationGroupAdjoin K x 1).subgroupOf (lowerRamificationGroupAdjoin K x 0))
        →* (IsLocalRing.ResidueField (adjoinIntegers K x))ˣ,
      Function.Injective θ := by
  obtain ⟨α, hspan, hα0, hadj⟩ := exists_uniformizer_ne_zero_adjoin_eq K x ht
  exact ⟨thetaMulQuot hspan hα0 hadj, thetaMulQuot_injective _ _ _⟩

/-- **完全分岐版(`G = G_0`)**: `G/G_1 ↪ 𝓀ˣ`。 -/
theorem exists_injective_thetaMulTop_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) :
    ∃ θ : ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
          ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ⧸
        lowerRamificationGroupAdjoin K x 1
        →* (IsLocalRing.ResidueField (adjoinIntegers K x))ˣ,
      Function.Injective θ := by
  obtain ⟨α, hspan, hα0, hadj⟩ := exists_uniformizer_ne_zero_adjoin_eq K x ht
  exact ⟨thetaMulTop hspan hα0 hadj (lowerRamificationGroupAdjoin_zero_eq_top K x ht),
    thetaMulTop_injective _ _ _ _⟩

/-- **Yoshida Prop 6.2 後半(`PAdicLocalField` 版、`n ≥ 1`)**: `G_n/G_{n+1} ↪ 𝓀⁺`。 -/
theorem exists_injective_thetaAdd_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {n : ℕ} (hn : 1 ≤ n) :
    ∃ θ : (lowerRamificationGroupAdjoin K x n ⧸
        (lowerRamificationGroupAdjoin K x (n + 1)).subgroupOf
          (lowerRamificationGroupAdjoin K x n))
        →* Multiplicative (IsLocalRing.ResidueField (adjoinIntegers K x)),
      Function.Injective θ := by
  obtain ⟨α, hspan, hα0, hadj⟩ := exists_uniformizer_ne_zero_adjoin_eq K x ht
  exact ⟨thetaAddQuot hspan hα0 hadj hn, thetaAddQuot_injective _ _ _ _⟩

/-- **Yoshida 正規化(`PAdicLocalField` 版、`n ≥ 1`)**: `G_n/G_{n+1} ↪ 𝔭^n/𝔭^{n+1}`。
★この版は素元の取り方に依らない(`thetaYoshida_independent`)。 -/
theorem exists_injective_thetaYoshida_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {n : ℕ} (hn : 1 ≤ n) :
    ∃ θ : (lowerRamificationGroupAdjoin K x n ⧸
        (lowerRamificationGroupAdjoin K x (n + 1)).subgroupOf
          (lowerRamificationGroupAdjoin K x n))
        →* Multiplicative (adjoinIntegers K x ⧸
          (IsLocalRing.maximalIdeal (adjoinIntegers K x)) ^ (n + 1)),
      Function.Injective θ := by
  obtain ⟨α, hspan, hα0, hadj⟩ := exists_uniformizer_ne_zero_adjoin_eq K x ht
  exact ⟨thetaYoshidaQuot hspan hα0 hadj hn, thetaYoshidaQuot_injective _ _ _ _⟩

end ABC3.Found.PGC
