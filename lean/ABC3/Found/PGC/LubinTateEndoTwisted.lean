import ABC3.Found.PGC.LubinTateThetaLinear

/-!
# `𝒪_L` 係数の `[θ]_{f,f′}`(Yoshida 2008 Proposition 3.5 の Frobenius ねじれ版)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Proposition 3.5(物理 p.5)
および Corollary 3.7(ii)(物理 p.6)。構造化済み原文は
`ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-3.html` の `#prop-3-5` / `#cor-3-7`。

原文 (Yoshida08 p.5):
> Proposition 3.5. Let f, f′ ∈ O[scr]_L[[X]] be as above, with linear coefficients π, π′ respectively. (i) There exists a unique formal group F_f over O[scr]_L such that f ∈ Hom_O[scr]_L(F_f, F^ϕ_f). (We call F_f the Lubin-Tate group associated to f.) (ii) There is a unique map [·]_f,f′ : Θ^L_π,π′ → (X) ⊂ O[scr]_L[[X]] such that: [θ]_f,f′(X) ≡ θX (mod deg 2), f′ ◦ [θ]_f,f′ = [θ]^ϕ_f,f′ ◦ f. It satisfies [θ]_f,f′ +_F_f′ [θ′]_f,f′ = [θ + θ′]_f,f′, [θ′]_f′,f′′ ◦ [θ]_f,f′ = [θθ′]_f,f′′. (iii) We have [θ]_f,f′ ∈ Hom_O[scr]_L(F_f, F_f′) for all θ ∈ Θ^L_π,π′.

原文 (Yoshida08 p.6, Corollary 3.7(ii)):
> (ii) If θ ∈ ΘL,× π,π′ := ΘL π,π′ ∩ O× L , then [θ]f,f′ is an isomorphism
> with the inverse [θ−1]f′,f.

原文 (Yoshida08 p.5, Definition 3.3):
> Definition 3.3. For uniformizers π, π′of L, set ΘL π,π′ := {θ ∈OL | θϕ/θ = π′/π}.
> It is an additive group. If θ ∈ΘL π,π′ and θ′ ∈ΘL π′,π′′, then θθ′ ∈ΘL π,π′′.

## 直前の波との関係(★ここが本ファイルの存在理由)

`Found/PGC/LubinTateThetaLinear.lean` は Proposition 3.5(ii)(iii) の
「対版」(`ϕ = id` かつ `π = π′`)を閉じたが、その逸脱記録に
「Frobenius ねじれ `[θ]^ϕ` を落とした」「`π = π′` に限る」と 2 点が残っていた。
★本ファイルは**その 2 点をどちらも外す**:

* `Θ^L_{π,π′}`(`π′θ = π ϕ(θ)`)の一般の `θ` を扱い、
* 関数等式は `f′ ∘ [θ] = [θ]^ϕ ∘ f`(係数に `ϕ` を当てる)であり、
* `π ≠ π′` を許す(`f` の 1 次係数 `π` は `𝔪` の元でありさえすればよい)。

★木に `𝒪_L` 係数の `LubinTateEndo` は**無かった**(`#158` の衝突検査で
`LubinTateEndoTwisted` ほか 10 名すべて衝突 0)。木の `LubinTateEndo`
(`Found/PGC/LubinTateEndoLimit.lean:126`)は
`hq : Fintype.card (ResidueField A) = q` を要求する——★すなわち剰余体が
**ちょうど `𝔽_q`** の場合しか扱えず、`k_L = 𝔽_{q^d}` では可除性
(`residue_divides_R_endo`)の Frobenius 恒等式が成り立たない。
★ねじれはまさにこの穴を埋める:剰余体が大きいときは `a ↦ a^q` が恒等でなく、
`ϕ` がその持ち上げになる。

## 何を足したか

### 純抽象核(分岐・付値・Galois・冪級数のどれも出てこない)

| 宣言 | 内容 |
|---|---|
| `comp_twist_intertwine_comp` | 結合律と「乗法的なねじれ `t`」だけ: `v∘a = t(a)∘u`・`w∘b = t(b)∘v` ならば `w∘(b∘a) = t(b∘a)∘u` |
| `exists_sub_linearMap_eq` | 中山: 有限加群上 `T(M) ⊆ I·M`(`I ≤ jacobson ⊥`)ならば `1 − T` は全射 |
| `mul_mem_smul_top` | `I • ⊤` は `A` の乗法で閉じる |

★1 本目は `M` と `comp : M → M → M` と `t : M → M` しか使わない
(群でもモノイドでもなくてよい)。2 本目は環論だが**分岐の語彙が 1 つも出ない**。

### 半線型方程式(ねじれの本体)

`ϕ`-半線型な方程式 `c − u·ϕ(c) = b`(`u ∈ 𝔪`)の**可解性**が、ねじれ版の
次数ごとの再帰を回す唯一の新しい入力である。★本ファイルはこれを舞台
(`TwistedLT`)の仮定 `hsolve` として括り出し、十分条件を 2 つ与える:

| 宣言 | 十分条件 |
|---|---|
| `exists_semilinear_solution_of_moduleFinite` | `A` が `R` 上**有限**・`ϕ` が `R` 上恒等・`u ∈ I·A`(`I ≤ jacobson ⊥`) |
| `hsolve_of_moduleFinite` | 上を `𝔪_A ⊆ I·A` と合わせて `hsolve` の形に整えたもの |
| (`ϕ = id` の場合) | `1 − u` が局所環の単数(`TwistedLT.ofUntwisted` の中) |

★**`L/K` が不分岐であることは `𝔪_L ⊆ 𝔪_K·𝒪_L` としてここに 1 度だけ入る**
(`hsolve_of_moduleFinite` の `hunram`)。有限性は `Module.Finite R A`。

同様に、一意性側の入力は「半線型の消去」`π′d = ϕ(d)·π^m (m ≥ 2) ⇒ d = 0`
であり、`semilinear_cancel_of_isNoetherian` が Noether 局所整域で証明する
(Krull の交叉定理 `Ideal.iInf_pow_eq_bot_of_isLocalRing`)。

### 抽象層(一般の可換環上の冪級数。分岐・Galois の語彙は出ない)

| 宣言 | 内容 |
|---|---|
| `map_twisted_obstruction_eq_zero` | 可除性: `f′∘φ − φ^ϕ∘f` の剰余体への還元は 0 |
| `base_case_twisted` | 出発点 `θ·X` は次数 `≤1` で障害を消す |
| `powerSeries_uniqueness_twisted` | ★ねじれ版一意性(原典 Lemma 3.4 の `t = 1`)。★`f ≠ f′` を**直接**扱う |
| `subst_twisted_intertwine_comp` | ねじれ絡み作用素の合成 |
| `subst_comp_eq_of_twisted_intertwine_pair` | ねじれ合成則の骨 |
| `theta_mul_mem` | 原典 Definition 3.3 の `θ ∈ Θ_{π,π′}`, `θ′ ∈ Θ_{π′,π″}` ⇒ `θθ′ ∈ Θ_{π,π″}` |

### 具体層

| 宣言 | 内容 |
|---|---|
| `TwistedLT` | 舞台一式(`ϕ`・`π`・`π′`・`f`・`f′`・`hϕres`・`hsolve`) |
| `LubinTateEndoTwisted` | ★`[θ]_{f,f′}` そのもの(近似列 `TwistedLT.φSeq` の極限) |
| `coeff_one_LubinTateEndoTwisted` | `[θ](X) ≡ θX (mod deg 2)` |
| `LubinTateEndoTwisted_functional_equation` | ★`f′ ∘ [θ] = [θ]^ϕ ∘ f` |
| `eq_LubinTateEndoTwisted` | ★原典 (ii) の**一意性**(条件を満たす級数は `[θ]` に限る) |
| `subst_LubinTateEndoTwisted_LubinTateEndoTwisted` | ★原典 (ii) の**合成則** `[θ′]∘[θ] = [θθ′]` |
| `subst_LubinTateEndoTwisted_eq_X` | ★原典 Corollary 3.7(ii) のねじれ版 |

## ★段取りとの差分

1. ★**ねじれ版の一意性は「二面版」を作る手間が要らない。** 直前の波は木の
   単一 `h` 版 `powerSeries_uniqueness` から逆向き絡み作用素 `η` と左簡約で
   二面版を組み立てた。本ファイルの `powerSeries_uniqueness_twisted` は
   `α − β` の次数ごとの帰納法を**直接**書いており、最初から `f ≠ f′` である
   (`η` も左簡約も出てこない)。★半線型の消去 `π′d = ϕ(d)π^m ⇒ d = 0` が
   `π^{n+1} − π` の因数分解の役を果たす。
2. ★**存在の構成に `IsDomain` も `Fintype (ResidueField A)` も要らなかった。**
   木の 1 変数構成は「`1 − π^n` が単数」を使うために局所環を、可除性のために
   `card k = q` を使う。ねじれ版では前者が `hsolve` に、後者が `hϕres`
   (`ϕ` が剰余体の `q` 乗写像を持ち上げる)に置き換わるので、
   ★**`[CommRing A] [IsLocalRing A]` と `[ExpChar (ResidueField A) pp]` だけで回る**。
   (一意性側では `IsDomain`・`IsNoetherianRing` を使う。)
3. ★**思ったより安かった箇所**: mathlib の `PowerSeries.map_expand` と
   `MvPowerSeries.map_iterateFrobenius_expand` が揃っていたので、可除性の
   ねじれ版(`map (r∘ϕ) φ = map (iterateFrobenius) (map r φ)` 経由)は
   木の `residue_divides_R_endo` とほぼ同じ長さで済んだ。

## 逸脱の記録

1. ★**`hsolve`(半線型方程式の可解性)を舞台の仮定に置いた。** 原典は
   `L` が完備であることから `α = −β − Σ_{i≥1}(π^{m+1}/π′)^{1+ϕ+⋯}β^{ϕ^i}` と
   **無限級数**で解く。本ファイルはその 1 行を仮定として括り出し、
   十分条件を「有限性 + 中山」(`exists_semilinear_solution_of_moduleFinite`)
   および「`ϕ = id`」(`TwistedLT.ofUntwisted`)で与えた。
   ★**完備(非有限)な `L`(= `K̂^ur`)に対する `hsolve` は本ファイルに無い**
   ——`IsAdicComplete` からの構成は新しいノード(下記)。
2. ★同様に一意性の入力 `hcancel` も仮定として括り出し、
   `semilinear_cancel_of_isNoetherian`(Noether 局所整域)で供給した。
3. ★原典は `π, π′` を**どちらも `L` の素元**とするが、本ファイルは
   `π ∈ 𝔪`(`hπmem`)と `𝔪 = (π′)`(`hπ'max`)しか要求しない
   ——原典より弱い仮定(一般化)である。一意性側だけは `π′ ≠ 0` 等を使う。
4. 原典 Proposition 3.5(ii) の**加法性** `[θ] +_{F_{f′}} [θ′] = [θ+θ′]` と
   (i)(形式群 `F_f` の存在)と (iii)(`[θ]` が形式群の準同型)は入っていない
   (2 変数の一意性のねじれ版が要る)。

## 退化の自己検査

* ★★`ϕ = id`(かつ `π = π′`, `f = f′`)に潰すと、
  `LubinTateEndoTwisted_ofUntwisted_eq_LubinTateEndo` が
  **木の `LubinTateEndo` と同じ冪級数になる**ことを証明している
  (等式そのもの。木の `powerSeries_uniqueness` で照合した)。
  ★このとき `hϕres` は `FiniteField.pow_card`(`x^q = x`)から、
  `hsolve` は `1 − u` が単数であることから出る——★**舞台が空でない**。
* ★`π ≠ π′` を許している(`TwistedLT` は `π` と `π'` を別々の場に持ち、
  `hπmem`・`hπ'max` しか課さない)。合成則
  `subst_LubinTateEndoTwisted_LubinTateEndoTwisted` は
  `f → f′ → f″` の 3 つの異なる級数の間の合成である。
* ★`L/K` の不分岐性を使った場所は 2 箇所に閉じ込めてある:
  (a) `hϕres`(`ϕ` が剰余体上 `q` 乗写像を誘導する = 算術 Frobenius であること)、
  (b) `hsolve_of_moduleFinite` の `hunram : 𝔪_A ⊆ I·A`(`𝔪_L = 𝔪_K𝒪_L`)。
* ★`ℕ∞` の切り詰め引き算・除算は 1 度も書いていない。`m − 1` の類は
  `obtain ⟨t, rfl⟩ : ∃ t, m = t + 1` で ℕ の引き算ごと避けた。

## 新しく必要になったノード

* ★`hsolve` の**完備版**: `IsAdicComplete 𝔪 A` と `ϕ(𝔪) ⊆ 𝔪` から
  `c − u·ϕ(c) = b` を解く(`Σ T^i b` の収束)。`L = K̂^ur` に要る。
* ★`TwistedLT` を木の実物(`unramifiedCompletionInt` / `arithFrobenius`)へ
  当てはめる配管。
* Corollary 4.9 具体化の項目 2(完備化レベルの半線型性
  `σ([θ](α)) = [θ^{(j)}](σα)`)——本ファイルは冪級数の等式までで、
  点で評価する層は入っていない。
* Proposition 3.5 の (i)(形式群 `F_f` の存在)・(iii)・加法性のねじれ版
  (2 変数の一意性補題のねじれ版が要る)。
-/

namespace ABC3.Found.PGC

def LubinTateEndoTwisted.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Proposition 3.5", sectionId := "prop-3-5" }

def LubinTateEndoTwisted_functional_equation.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Proposition 3.5", sectionId := "prop-3-5" }

def eq_LubinTateEndoTwisted.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Proposition 3.5", sectionId := "prop-3-5" }

def subst_LubinTateEndoTwisted_LubinTateEndoTwisted.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Proposition 3.5", sectionId := "prop-3-5" }

def subst_LubinTateEndoTwisted_eq_X.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 6, item := "Corollary 3.7", sectionId := "cor-3-7" }

/-! ## 0. 純抽象核

★型 `M` と二項演算 `comp` と「乗法的なねじれ `t`」しか使わない。冪級数も局所環も
分岐も出てこない。`comp x y` は「`x` のあとに `y` を代入する」= `x ∘ y` と読む。 -/

/-- ★★★★★★★★**純抽象核(ねじれ版)** —— ねじれ絡み作用素の合成はまた
ねじれ絡み作用素。

`v ∘ a = t(a) ∘ u`(`a : u → v`)と `w ∘ b = t(b) ∘ v`(`b : v → w`)から
`w ∘ (b ∘ a) = t(b ∘ a) ∘ u`(`b ∘ a : u → w`)。

仮定は結合律 `hassoc` と `t` の乗法性 `ht` だけ。群でもモノイドでもなくてよい。
★`t = id` と置くと直前の波の `comp_intertwine_comp` に一致する。 -/
theorem comp_twist_intertwine_comp {M : Type*} (comp : M → M → M)
    (hassoc : ∀ x y z : M, comp (comp x y) z = comp x (comp y z))
    (t : M → M) (ht : ∀ x y : M, t (comp x y) = comp (t x) (t y))
    (u v w a b : M) (ha : comp v a = comp (t a) u) (hb : comp w b = comp (t b) v) :
    comp w (comp b a) = comp (t (comp b a)) u := by
  calc comp w (comp b a) = comp (comp w b) a := (hassoc w b a).symm
    _ = comp (comp (t b) v) a := by rw [hb]
    _ = comp (t b) (comp v a) := hassoc (t b) v a
    _ = comp (t b) (comp (t a) u) := by rw [ha]
    _ = comp (comp (t b) (t a)) u := (hassoc (t b) (t a) u).symm
    _ = comp (t (comp b a)) u := by rw [ht]

/-- `I • ⊤` は `A` の乗法で閉じている。 -/
theorem mul_mem_smul_top {R A : Type*} [CommRing R] [CommRing A] [Algebra R A] {I : Ideal R}
    {u : A} (hu : u ∈ I • (⊤ : Submodule R A)) (y : A) :
    u * y ∈ I • (⊤ : Submodule R A) := by
  refine Submodule.smul_induction_on hu ?_ ?_
  · intro r hr x _
    rw [smul_mul_assoc]
    exact Submodule.smul_mem_smul hr Submodule.mem_top
  · intro x z hx hz
    rw [add_mul]
    exact Submodule.add_mem _ hx hz

/-- ★★★★★★★★**純抽象核(中山)** —— 有限加群上、像がジャコブソン根基に
含まれるイデアル倍に入る線型作用素 `T` に対し `1 − T` は全射。

分岐も付値も冪級数も出てこない。ねじれ版の再帰を回す唯一の新しい入力
(半線型方程式 `c − u·ϕ(c) = b` の可解性)はこれに帰着する。 -/
theorem exists_sub_linearMap_eq {R M : Type*} [CommRing R] [AddCommGroup M] [Module R M]
    [Module.Finite R M] {I : Ideal R} (hI : I ≤ Ideal.jacobson ⊥) (T : M →ₗ[R] M)
    (hT : ∀ x : M, T x ∈ I • (⊤ : Submodule R M)) (b : M) : ∃ c : M, c - T c = b := by
  have hle : (⊤ : Submodule R M) ≤ LinearMap.range (LinearMap.id - T) ⊔ I • ⊤ := by
    intro x _
    have hx : x = ((LinearMap.id - T : M →ₗ[R] M)) x + T x := by
      simp only [LinearMap.sub_apply, LinearMap.id_apply]
      abel
    rw [hx]
    exact Submodule.add_mem_sup (LinearMap.mem_range_self _ x) (hT x)
  have htop := Submodule.le_of_le_smul_of_le_jacobson_bot (Module.Finite.fg_top) hI hle
  obtain ⟨c, hc⟩ := htop (Submodule.mem_top : b ∈ ⊤)
  refine ⟨c, ?_⟩
  rw [← hc]
  simp only [LinearMap.sub_apply, LinearMap.id_apply]

/-- ★★★★**半線型方程式の可解性(有限の場合)**: `A` が `R` 上有限で `ϕ` が `R` 上
恒等、`u ∈ I·A`(`I ≤ jacobson ⊥`)ならば `c − u·ϕ(c) = b` は解ける。 -/
theorem exists_semilinear_solution_of_moduleFinite {R A : Type*} [CommRing R] [CommRing A]
    [Algebra R A] [Module.Finite R A] {I : Ideal R} (hI : I ≤ Ideal.jacobson ⊥) (ϕ : A →+* A)
    (hϕR : ∀ r : R, ϕ (algebraMap R A r) = algebraMap R A r)
    {u : A} (hu : u ∈ I • (⊤ : Submodule R A)) (b : A) : ∃ c : A, c - u * ϕ c = b := by
  let T : A →ₗ[R] A :=
    { toFun := fun x => u * ϕ x
      map_add' := by intro x y; show u * ϕ (x + y) = u * ϕ x + u * ϕ y; rw [map_add, mul_add]
      map_smul' := by
        intro r x
        show u * ϕ (r • x) = r • (u * ϕ x)
        rw [Algebra.smul_def, map_mul, hϕR, Algebra.smul_def]
        ring }
  exact exists_sub_linearMap_eq hI T (fun x => mul_mem_smul_top hu (ϕ x)) b

/-- ★★★★`TwistedLT` の `hsolve` の十分条件。★**`L/K` の不分岐性はここに
`hunram : 𝔪_A ⊆ I·A`(`𝔪_L = 𝔪_K𝒪_L`)として 1 度だけ入る。** -/
theorem hsolve_of_moduleFinite {R A : Type*} [CommRing R] [CommRing A] [IsLocalRing A]
    [Algebra R A] [Module.Finite R A] {I : Ideal R} (hI : I ≤ Ideal.jacobson ⊥) (ϕ : A →+* A)
    (hϕR : ∀ r : R, ϕ (algebraMap R A r) = algebraMap R A r)
    (hunram : ∀ x ∈ IsLocalRing.maximalIdeal A, x ∈ I • (⊤ : Submodule R A)) :
    ∀ u ∈ IsLocalRing.maximalIdeal A, ∀ b : A, ∃ c : A, c - u * ϕ c = b :=
  fun u hu b => exists_semilinear_solution_of_moduleFinite hI ϕ hϕR (hunram u hu) b

/-! ## 1. 代入モノイド上の Frobenius ねじれ

直前の波が作った `SubstNilp`(定数項 `0` の冪級数)と `substComp`(合成)の上に、
係数への `ϕ` の作用を乗せる。★結合律は `substComp_assoc` を使い回すだけ。 -/

/-- 代入モノイド `SubstNilp` 上の Frobenius ねじれ `x ↦ x^ϕ`。 -/
noncomputable def substTwist {A : Type*} [CommRing A] (ϕ : A →+* A) (x : SubstNilp A) :
    SubstNilp A :=
  ⟨PowerSeries.map ϕ x.1, by
    show MvPowerSeries.constantCoeff (MvPowerSeries.map ϕ x.1) = 0
    rw [MvPowerSeries.constantCoeff_map]
    show ϕ (PowerSeries.constantCoeff x.1) = 0
    rw [x.2, map_zero]⟩

/-- ★ねじれは合成と可換(`(x ∘ y)^ϕ = x^ϕ ∘ y^ϕ`)——純抽象核の `ht` を供給する。 -/
theorem substTwist_substComp {A : Type*} [CommRing A] (ϕ : A →+* A) (x y : SubstNilp A) :
    substTwist ϕ (substComp x y) = substComp (substTwist ϕ x) (substTwist ϕ y) :=
  Subtype.ext (powerSeries_map_subst (hasSubst_of_constantCoeff_zero y.2) ϕ x.1)

/-- 純抽象核(ねじれ版)の冪級数版: `a : u → v` と `b : v → w` の合成 `b ∘ a` は `u → w`。 -/
theorem subst_twisted_intertwine_comp {A : Type*} [CommRing A] (ϕ : A →+* A)
    {u v w a b : PowerSeries A}
    (hu0 : PowerSeries.constantCoeff u = 0) (hv0 : PowerSeries.constantCoeff v = 0)
    (hw0 : PowerSeries.constantCoeff w = 0) (ha0 : PowerSeries.constantCoeff a = 0)
    (hb0 : PowerSeries.constantCoeff b = 0)
    (ha : PowerSeries.subst a v = PowerSeries.subst u (PowerSeries.map ϕ a))
    (hb : PowerSeries.subst b w = PowerSeries.subst v (PowerSeries.map ϕ b)) :
    PowerSeries.subst (PowerSeries.subst a b) w =
      PowerSeries.subst u (PowerSeries.map ϕ (PowerSeries.subst a b)) :=
  congrArg Subtype.val
    (comp_twist_intertwine_comp (M := SubstNilp A) substComp substComp_assoc
      (substTwist ϕ) (substTwist_substComp ϕ)
      ⟨u, hu0⟩ ⟨v, hv0⟩ ⟨w, hw0⟩ ⟨a, ha0⟩ ⟨b, hb0⟩ (Subtype.ext ha) (Subtype.ext hb))

/-! ## 2. 冪級数の小道具(一般の可換環) -/

/-- `map ϕ (a • X^m) = ϕ(a) • X^m`。 -/
theorem map_smul_X_pow {A : Type*} [CommRing A] (ϕ : A →+* A) (a : A) (m : ℕ) :
    PowerSeries.map ϕ (a • (PowerSeries.X : PowerSeries A) ^ m) =
      ϕ a • (PowerSeries.X : PowerSeries A) ^ m := by
  ext n
  simp only [PowerSeries.coeff_map, PowerSeries.coeff_smul, smul_eq_mul, map_mul,
    PowerSeries.coeff_X_pow]
  split_ifs <;> simp

theorem constantCoeff_smul_X_pow {A : Type*} [CommRing A] (a : A) {m : ℕ} (hm : m ≠ 0) :
    PowerSeries.constantCoeff (a • (PowerSeries.X : PowerSeries A) ^ m) = 0 := by
  rw [PowerSeries.constantCoeff_smul]
  rw [show PowerSeries.constantCoeff ((PowerSeries.X : PowerSeries A) ^ m) = 0 by
    rw [map_pow, PowerSeries.constantCoeff_X, zero_pow hm]]
  rw [smul_zero]

theorem order_X_pow_eq {A : Type*} [CommRing A] [Nontrivial A] (m : ℕ) :
    MvPowerSeries.order ((PowerSeries.X : PowerSeries A) ^ m) = ((m : ℕ) : ℕ∞) := by
  rw [← PowerSeries.order_eq_order]
  exact PowerSeries.order_eq.mpr ⟨fun i hi => by
      have hcx := PowerSeries.coeff_X_pow (R := A) i m
      rw [if_pos (by exact_mod_cast hi)] at hcx
      rw [hcx]; exact one_ne_zero, fun i hi => by
      rw [PowerSeries.coeff_X_pow]
      rw [if_neg (by intro heq; rw [heq] at hi; exact absurd hi (lt_irrefl _))]⟩

/-- 係数環の写像は次数を下げない。 -/
theorem nat_le_order_map {A B : Type*} [CommRing A] [CommRing B] (ϕ : A →+* B)
    {δ : PowerSeries A} {m : ℕ} (h : ((m : ℕ) : ℕ∞) ≤ MvPowerSeries.order δ) :
    ((m : ℕ) : ℕ∞) ≤ MvPowerSeries.order (PowerSeries.map ϕ δ) := by
  rw [← PowerSeries.order_eq_order]
  apply PowerSeries.nat_le_order
  intro i hi
  rw [PowerSeries.coeff_map]
  have hz : PowerSeries.coeff i δ = 0 := by
    apply PowerSeries.coeff_of_lt_order
    rw [PowerSeries.order_eq_order]
    exact lt_of_lt_of_le (by exact_mod_cast hi) h
  rw [hz, map_zero]

/-! ## 3. 可除性(ねじれ版)

★木の `residue_divides_R_endo`(`Found/PGC/LubinTateEndoDivisibility.lean`)は
剰余体の位数がちょうど `q` であることを使う。ねじれ版では代わりに
「`ϕ` が剰余体上 `q` 乗写像を誘導する」(`hϕ`)を使うので、
★**剰余体が `𝔽_q` より大きくてよい**(`Fintype` すら要らない)。 -/

/-- ★★★★**ねじれ版の可除性(剰余体版)**: `f`・`f′` がどちらも剰余体上で `X^q` に
還元され、`ϕ` が剰余体上 `q` 乗写像を誘導するなら、ねじれ障害
`f′ ∘ φ − φ^ϕ ∘ f` の剰余体への還元は恒等的に 0。

証明: `f′∘φ` の還元は `φ̄^q`。`φ^ϕ∘f` の還元は
`expand q (map (frobenius^ff) φ̄) = map (frobenius^ff) (expand q φ̄) = φ̄^q`
(mathlib の `PowerSeries.map_expand` と `MvPowerSeries.map_iterateFrobenius_expand`)。
★`map r ∘ map ϕ = map (frobenius^ff) ∘ map r` は木の在庫
`Found/PGC/DworkThetaStep2.lean::map_map_eq_map_iterateFrobenius` をそのまま使う。 -/
theorem map_twisted_obstruction_eq_zero {A κ : Type*} [CommRing A] [CommRing κ]
    {pp ff : ℕ} [ExpChar κ pp] (r : A →+* κ) (ϕ : A →+* A)
    (hϕ : ∀ a : A, r (ϕ a) = r a ^ pp ^ ff)
    (f f' : PowerSeries A) (hf0 : PowerSeries.constantCoeff f = 0)
    (hf : PowerSeries.map r f = PowerSeries.X ^ pp ^ ff)
    (hf' : PowerSeries.map r f' = PowerSeries.X ^ pp ^ ff)
    (φ : PowerSeries A) (hφ0 : PowerSeries.constantCoeff φ = 0) :
    PowerSeries.map r (PowerSeries.subst φ f' - PowerSeries.subst f (PowerSeries.map ϕ φ)) = 0 := by
  have hq0 : pp ^ ff ≠ 0 := pow_ne_zero ff (expChar_ne_zero κ pp)
  have hφHS : PowerSeries.HasSubst φ := by
    show IsNilpotent (PowerSeries.constantCoeff φ); rw [hφ0]; exact IsNilpotent.zero
  have hfHS : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0]; exact IsNilpotent.zero
  set φbar := PowerSeries.map r φ with hφbar
  have hφbarHS : PowerSeries.HasSubst φbar := by
    show IsNilpotent (MvPowerSeries.constantCoeff (MvPowerSeries.map r φ))
    rw [MvPowerSeries.constantCoeff_map]
    show IsNilpotent (r (PowerSeries.constantCoeff φ))
    rw [hφ0, map_zero]; exact IsNilpotent.zero
  rw [map_sub]
  have h1 : PowerSeries.map r (PowerSeries.subst φ f') = φbar ^ pp ^ ff := by
    show MvPowerSeries.map r (PowerSeries.subst φ f') = _
    rw [powerSeries_map_subst hφHS r f']
    show PowerSeries.subst φbar (PowerSeries.map r f') = _
    rw [hf']
    exact substXpow_eq_pow hφbarHS
  have h2 : PowerSeries.map r (PowerSeries.subst f (PowerSeries.map ϕ φ)) = φbar ^ pp ^ ff := by
    show MvPowerSeries.map r (PowerSeries.subst f (PowerSeries.map ϕ φ)) = _
    rw [powerSeries_map_subst hfHS r (PowerSeries.map ϕ φ)]
    show PowerSeries.subst (PowerSeries.map r f) (PowerSeries.map r (PowerSeries.map ϕ φ)) = _
    -- ★木の在庫: `Found/PGC/DworkThetaStep2.lean::map_map_eq_map_iterateFrobenius`
    rw [hf, map_map_eq_map_iterateFrobenius r ϕ hϕ φ, ← hφbar]
    calc PowerSeries.subst ((PowerSeries.X : PowerSeries κ) ^ pp ^ ff)
          (PowerSeries.map (iterateFrobenius κ pp ff) φbar)
        = MvPowerSeries.expand (pp ^ ff) hq0 (PowerSeries.map (iterateFrobenius κ pp ff) φbar) :=
          (expand_eq_subst_pow hq0 _).symm
      _ = PowerSeries.map (iterateFrobenius κ pp ff) (PowerSeries.expand (pp ^ ff) hq0 φbar) :=
          (PowerSeries.map_expand _ hq0 _ _).symm
      _ = φbar ^ pp ^ ff :=
          MvPowerSeries.map_iterateFrobenius_expand pp (expChar_ne_zero κ pp) φbar ff
  rw [h1, h2, sub_self]

/-! ## 4. 出発点 -/

/-- ★出発点 `θ • X` はねじれ版の障害を次数 `≤1` の範囲で消す。

★ここが `θ ∈ Θ^L_{π,π′}`(`π′θ = π ϕ(θ)`、原典 Definition 3.3)を使う唯一の場所である
——1 次の係数は `f′∘(θX)` 側が `π′θ`、`(θX)^ϕ∘f` 側が `ϕ(θ)π` だから。 -/
theorem base_case_twisted {A : Type*} [CommRing A] {π π' θ : A} (ϕ : A →+* A)
    (hθ : π' * θ = π * ϕ θ)
    (f : PowerSeries A) (hf0 : PowerSeries.constantCoeff f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (f' : PowerSeries A) (hf'0 : PowerSeries.coeff 0 f' = 0) (hf'1 : PowerSeries.coeff 1 f' = π')
    (k : ℕ) (hk : k ≤ 1) :
    PowerSeries.coeff k
      (PowerSeries.subst (θ • (PowerSeries.X : PowerSeries A)) f' -
        PowerSeries.subst f (PowerSeries.map ϕ (θ • (PowerSeries.X : PowerSeries A)))) = 0 := by
  have hφ₁cc : PowerSeries.constantCoeff (θ • (PowerSeries.X : PowerSeries A)) = 0 := by
    have h := constantCoeff_smul_X_pow (A := A) θ (m := 1) one_ne_zero
    rwa [pow_one] at h
  have hφ₁order : (1 : ℕ∞) ≤ MvPowerSeries.order (θ • (PowerSeries.X : PowerSeries A)) :=
    MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hφ₁cc
  have h0order : (1 : ℕ∞) ≤ MvPowerSeries.order (0 : PowerSeries A) := by
    simp [MvPowerSeries.order_zero]
  have h00 : PowerSeries.constantCoeff (0 : PowerSeries A) = 0 := map_zero _
  have hFside := coeff_subst_linearize_1var h00 hφ₁cc h0order hφ₁order le_rfl f' π' hf'0 hf'1 k
    (by exact_mod_cast hk)
  rw [zero_add, coeff_subst_zero_eq_zero_1var f' hf'0 k, sub_zero] at hFside
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0]; exact IsNilpotent.zero
  have hGeq : PowerSeries.subst f (PowerSeries.map ϕ (θ • (PowerSeries.X : PowerSeries A))) =
      ϕ θ • f := by
    have h1 := map_smul_X_pow ϕ θ 1
    rw [pow_one] at h1
    rw [h1, PowerSeries.subst_smul hHSf, PowerSeries.subst_X hHSf]
  rw [map_sub, hFside, hGeq, PowerSeries.coeff_smul, PowerSeries.coeff_X, PowerSeries.coeff_smul]
  split_ifs with h
  · subst h
    rw [hf1, smul_eq_mul, smul_eq_mul, mul_one, hθ]
    ring
  · have hk0 : k = 0 := by omega
    subst hk0
    rw [smul_zero, mul_zero, PowerSeries.coeff_zero_eq_constantCoeff_apply, hf0, smul_zero, sub_zero]

/-! ## 5. 一意性(ねじれ版)

★直前の波は木の単一 `h` 版から「二面版」を組み立てたが、ねじれ版は
`α − β` の次数ごとの帰納法を**直接**書けば最初から `f ≠ f′` を扱える。 -/

/-- 原典 Definition 3.3 の「`θ ∈ Θ_{π,π′}`・`θ′ ∈ Θ_{π′,π″}` ならば `θθ′ ∈ Θ_{π,π″}`」。 -/
theorem theta_mul_mem {A : Type*} [CommRing A] (ϕ : A →+* A) {π π' π'' θ θ' : A}
    (hθ : π' * θ = π * ϕ θ) (hθ' : π'' * θ' = π' * ϕ θ') :
    π'' * (θ * θ') = π * ϕ (θ * θ') := by
  rw [map_mul]
  calc π'' * (θ * θ') = (π'' * θ') * θ := by ring
    _ = (π' * ϕ θ') * θ := by rw [hθ']
    _ = (π' * θ) * ϕ θ' := by ring
    _ = (π * ϕ θ) * ϕ θ' := by rw [hθ]
    _ = π * (ϕ θ * ϕ θ') := by ring

/-- ★★★★★★★★**ねじれ版の一意性補題**(原典 Lemma 3.4 の `t = 1`・Frobenius ねじれ版)。

`f′ ∘ α = α^ϕ ∘ f` と `f′ ∘ β = β^ϕ ∘ f` を満たし、定数項が `0` で 1 次係数が
一致する 2 つの冪級数は等しい。

証明: `δ := α − β` の次数が `≥ n+1`(`n ≥ 1`)なら、次数 `n+1` の係数の比較で
`π′δ_{n+1} = ϕ(δ_{n+1})π^{n+1}` が出る。`hcancel` でこれが `δ_{n+1} = 0` を与え、
次数が 1 つ上がる。★`hcancel` が木の untwisted 版における
「`π^{n+1} − π = π(1 − π^n)` が非零」の役を果たす。 -/
theorem powerSeries_uniqueness_twisted {A : Type*} [CommRing A] {π π' : A} (ϕ : A →+* A)
    (hcancel : ∀ m : ℕ, 2 ≤ m → ∀ d : A, π' * d = ϕ d * π ^ m → d = 0)
    {f f' : PowerSeries A}
    (hf0 : PowerSeries.constantCoeff f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf'0 : PowerSeries.coeff 0 f' = 0) (hf'1 : PowerSeries.coeff 1 f' = π')
    {α β : PowerSeries A} (hα0 : PowerSeries.constantCoeff α = 0)
    (hβ0 : PowerSeries.constantCoeff β = 0)
    (hlead : PowerSeries.coeff 1 α = PowerSeries.coeff 1 β)
    (hα : PowerSeries.subst α f' = PowerSeries.subst f (PowerSeries.map ϕ α))
    (hβ : PowerSeries.subst β f' = PowerSeries.subst f (PowerSeries.map ϕ β)) :
    α = β := by
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0]; exact IsNilpotent.zero
  set δ := α - β with hδ_def
  have hδ0 : PowerSeries.constantCoeff δ = 0 := by rw [hδ_def, map_sub, hα0, hβ0, sub_zero]
  have hδ1 : PowerSeries.coeff 1 δ = 0 := by rw [hδ_def, map_sub, hlead, sub_self]
  have hbase : ((2 : ℕ) : ℕ∞) ≤ δ.order := by
    apply PowerSeries.nat_le_order
    intro i hi
    interval_cases i
    · rw [PowerSeries.coeff_zero_eq_constantCoeff]; exact hδ0
    · exact hδ1
  have hstep : ∀ n : ℕ, 1 ≤ n → ((n + 1 : ℕ) : ℕ∞) ≤ δ.order →
      ((n + 2 : ℕ) : ℕ∞) ≤ δ.order := by
    intro n hn1 hδorder
    apply PowerSeries.nat_le_order
    intro i hi
    rcases lt_or_eq_of_le (by omega : i ≤ n + 1) with hilt | hieq
    · exact PowerSeries.coeff_of_lt_order i (lt_of_lt_of_le (by exact_mod_cast hilt) hδorder)
    · subst hieq
      have hδorderMv : ((n + 1 : ℕ) : ℕ∞) ≤ MvPowerSeries.order (δ : PowerSeries A) := by
        rw [← PowerSeries.order_eq_order]; exact hδorder
      have hβorder : (1 : ℕ∞) ≤ MvPowerSeries.order (β : PowerSeries A) :=
        MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hβ0
      have hlin := coeff_subst_linearize_1var hβ0 hδ0 hβorder hδorderMv (by omega : 1 ≤ n + 1)
        f' π' hf'0 hf'1 (n + 1) le_rfl
      rw [show β + δ = α from by rw [hδ_def]; ring] at hlin
      have hmapsub : PowerSeries.map ϕ α - PowerSeries.map ϕ β = PowerSeries.map ϕ δ := by
        rw [hδ_def, map_sub]
      have hsplit : PowerSeries.subst f (PowerSeries.map ϕ α) -
          PowerSeries.subst f (PowerSeries.map ϕ β) =
          PowerSeries.subst f (PowerSeries.map ϕ δ) := by
        rw [← hmapsub, PowerSeries.subst_sub hHSf]
      rw [hα, hβ] at hlin
      have hcombine : PowerSeries.coeff (n + 1) (PowerSeries.subst f (PowerSeries.map ϕ δ)) =
          π' * PowerSeries.coeff (n + 1) δ := by
        rw [← hsplit, map_sub]; exact hlin
      have hmapδorder : ((n + 1 : ℕ) : ℕ∞) ≤
          PowerSeries.order (PowerSeries.map ϕ δ : PowerSeries A) := by
        rw [PowerSeries.order_eq_order]
        exact nat_le_order_map ϕ hδorderMv
      have hkey := coeff_subst_eq_of_order_ge (δ := PowerSeries.map ϕ δ) (h := f) (π := π) (n := n)
        hf0 hf1 hmapδorder hHSf
      rw [hkey, PowerSeries.coeff_map] at hcombine
      exact hcancel (n + 1) (by omega) _ hcombine.symm
  have hall : ∀ n : ℕ, 1 ≤ n → ((n + 1 : ℕ) : ℕ∞) ≤ δ.order := by
    intro n hn
    induction n, hn using Nat.le_induction with
    | base => exact hbase
    | succ n hn1 ih => exact hstep n hn1 ih
  have horder : δ.order = ⊤ := by
    by_contra hne
    obtain ⟨m, hm⟩ := WithTop.ne_top_iff_exists.mp hne
    have hle := hall (m + 1) (by omega)
    rw [← hm] at hle
    have hcontra : ((m : ℕ) : ℕ∞) < ((m + 1 + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega)
    exact absurd hle (not_le.mpr hcontra)
  have hδzero : δ = 0 := PowerSeries.order_eq_top.mp horder
  rw [hδ_def, sub_eq_zero] at hδzero
  exact hδzero

/-- ★★★★**半線型の消去**(ねじれ版一意性の入力 `hcancel` の十分条件)。

Noether 局所整域で `𝔪 = (π′)`・`π ∈ 𝔪` なら、`π′d = ϕ(d)π^m`(`m ≥ 2`)から
`d ∈ ⋂_k 𝔪^k = 0`(Krull の交叉定理)。★`ϕ` が `𝔪` を保つことだけ使う。 -/
theorem semilinear_cancel_of_isNoetherian {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsNoetherianRing A] {π π' : A}
    (hπ'max : IsLocalRing.maximalIdeal A = Ideal.span {π'}) (hπ'ne0 : π' ≠ 0)
    (hπmem : π ∈ IsLocalRing.maximalIdeal A) (ϕ : A →+* A)
    (hϕloc : ∀ x ∈ IsLocalRing.maximalIdeal A, ϕ x ∈ IsLocalRing.maximalIdeal A)
    (m : ℕ) (hm : 2 ≤ m) (d : A) (h : π' * d = ϕ d * π ^ m) : d = 0 := by
  obtain ⟨t, rfl⟩ : ∃ t, m = t + 1 := ⟨m - 1, by omega⟩
  have ht1 : 1 ≤ t := by omega
  have hπ'mem : π' ∈ IsLocalRing.maximalIdeal A := hπ'max ▸ Ideal.mem_span_singleton_self π'
  obtain ⟨s, hs⟩ : ∃ s : A, π = π' * s := by
    have hmem : π ∈ Ideal.span ({π'} : Set A) := hπ'max ▸ hπmem
    obtain ⟨s, hs⟩ := Ideal.mem_span_singleton'.mp hmem
    exact ⟨s, by rw [← hs]; ring⟩
  have hmapm : Ideal.map ϕ (IsLocalRing.maximalIdeal A) ≤ IsLocalRing.maximalIdeal A := by
    rw [Ideal.map_le_iff_le_comap]
    intro y hy
    exact hϕloc y hy
  have hϕpow : ∀ (k : ℕ) (x : A), x ∈ IsLocalRing.maximalIdeal A ^ k →
      ϕ x ∈ IsLocalRing.maximalIdeal A ^ k := by
    intro k x hx
    have hmap : Ideal.map ϕ (IsLocalRing.maximalIdeal A ^ k) ≤ IsLocalRing.maximalIdeal A ^ k := by
      rw [Ideal.map_pow]
      exact Ideal.pow_right_mono hmapm k
    exact hmap (Ideal.mem_map_of_mem ϕ hx)
  have hcoefmem : π' ^ t * s ^ (t + 1) ∈ IsLocalRing.maximalIdeal A :=
    Ideal.mul_mem_right _ _ (Ideal.pow_mem_of_mem _ hπ'mem t ht1)
  have hd_eq : d = ϕ d * (π' ^ t * s ^ (t + 1)) := by
    refine (mul_left_cancel₀ hπ'ne0 ?_).symm
    rw [h, hs]; ring
  have hdk : ∀ k : ℕ, d ∈ IsLocalRing.maximalIdeal A ^ k := by
    intro k
    induction k with
    | zero => simp
    | succ k ih =>
      rw [hd_eq, pow_succ]
      exact Ideal.mul_mem_mul (hϕpow k d ih) hcoefmem
  have hmemInf : d ∈ ⨅ k : ℕ, IsLocalRing.maximalIdeal A ^ k := Ideal.mem_iInf.mpr hdk
  rw [Ideal.iInf_pow_eq_bot_of_isLocalRing _ (IsLocalRing.maximalIdeal.isMaximal A).ne_top]
    at hmemInf
  exact hmemInf

/-- ★★★★★★★★**抽象層のねじれ合成則**(原典 Proposition 3.5(ii) 第2式の骨)。

`a : F_u → F_v`・`b : F_v → F_w`・`c : F_u → F_w` がいずれもねじれ絡み作用素で、
`coeff 1 c = coeff 1 b · coeff 1 a` なら `b ∘ a = c`。 -/
theorem subst_comp_eq_of_twisted_intertwine_pair {A : Type*} [CommRing A] {π π'' : A}
    (ϕ : A →+* A) (hcancel : ∀ m : ℕ, 2 ≤ m → ∀ d : A, π'' * d = ϕ d * π ^ m → d = 0)
    {u v w a b c : PowerSeries A}
    (hu0 : PowerSeries.constantCoeff u = 0) (hu1 : PowerSeries.coeff 1 u = π)
    (hv0 : PowerSeries.constantCoeff v = 0)
    (hw0 : PowerSeries.constantCoeff w = 0) (hw1 : PowerSeries.coeff 1 w = π'')
    (ha0 : PowerSeries.constantCoeff a = 0) (hb0 : PowerSeries.constantCoeff b = 0)
    (hc0 : PowerSeries.constantCoeff c = 0)
    (ha : PowerSeries.subst a v = PowerSeries.subst u (PowerSeries.map ϕ a))
    (hb : PowerSeries.subst b w = PowerSeries.subst v (PowerSeries.map ϕ b))
    (hc : PowerSeries.subst c w = PowerSeries.subst u (PowerSeries.map ϕ c))
    (hlead : PowerSeries.coeff 1 c = PowerSeries.coeff 1 b * PowerSeries.coeff 1 a) :
    PowerSeries.subst a b = c := by
  refine powerSeries_uniqueness_twisted ϕ hcancel hu0 hu1
    ((PowerSeries.coeff_zero_eq_constantCoeff_apply w).trans hw0) hw1
    (PowerSeries.constantCoeff_subst_eq_zero ha0 b hb0) hc0 ?_
    (subst_twisted_intertwine_comp ϕ hu0 hv0 hw0 ha0 hb0 ha hb) hc
  rw [coeff_one_subst_1var ha0, hlead]

/-! ## 6. 舞台 -/

/-- ねじれ版 Lubin-Tate の舞台一式。

`ϕ` は原典の算術 Frobenius(`hϕres` が「剰余体上 `q` 乗写像を誘導する」という
唯一の要請)、`f`(1 次係数 `π`)から `f′`(1 次係数 `π′`)への
`[·]_{f,f′}` を作るための道具立てである。★`hsolve` は半線型方程式の可解性で、
`hsolve_of_moduleFinite` などが供給する。 -/
structure TwistedLT (A : Type*) [CommRing A] [IsLocalRing A] (pp ff : ℕ) where
  /-- 係数環の自己準同型(原典の Frobenius `ϕ`)。 -/
  ϕ : A →+* A
  /-- `f` の 1 次係数。 -/
  π : A
  /-- `f′` の 1 次係数(`𝔪` の生成元)。 -/
  π' : A
  hπmem : π ∈ IsLocalRing.maximalIdeal A
  hπ'max : IsLocalRing.maximalIdeal A = Ideal.span {π'}
  hϕres : ∀ a : A, IsLocalRing.residue A (ϕ a) = IsLocalRing.residue A a ^ pp ^ ff
  hsolve : ∀ u ∈ IsLocalRing.maximalIdeal A, ∀ b : A, ∃ c : A, c - u * ϕ c = b
  f : PowerSeries A
  hf0 : PowerSeries.constantCoeff f = 0
  hf1 : PowerSeries.coeff 1 f = π
  hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ pp ^ ff
  f' : PowerSeries A
  hf'0 : PowerSeries.constantCoeff f' = 0
  hf'1 : PowerSeries.coeff 1 f' = π'
  hf'res : PowerSeries.map (IsLocalRing.residue A) f' = PowerSeries.X ^ pp ^ ff

namespace TwistedLT

variable {A : Type*} [CommRing A] [IsLocalRing A] {pp ff : ℕ} (S : TwistedLT A pp ff)

theorem coeff_zero_f : PowerSeries.coeff 0 S.f = 0 :=
  (PowerSeries.coeff_zero_eq_constantCoeff_apply S.f).trans S.hf0

theorem coeff_zero_f' : PowerSeries.coeff 0 S.f' = 0 :=
  (PowerSeries.coeff_zero_eq_constantCoeff_apply S.f').trans S.hf'0

theorem hasSubst_f : PowerSeries.HasSubst S.f := by
  show IsNilpotent (PowerSeries.constantCoeff S.f); rw [S.hf0]; exact IsNilpotent.zero

/-- ねじれ版の障害 `f′ ∘ φ − φ^ϕ ∘ f`。 -/
noncomputable def obstruction (φ : PowerSeries A) : PowerSeries A :=
  PowerSeries.subst φ S.f' - PowerSeries.subst S.f (PowerSeries.map S.ϕ φ)

/-- 次数ごとの再帰の不変量。 -/
def SeqInvariant (φ : PowerSeries A) (k : ℕ) : Prop :=
  (∀ j ≤ k + 1, PowerSeries.coeff j (S.obstruction φ) = 0) ∧ PowerSeries.constantCoeff φ = 0

/-- 舞台の合成 `f → f′ → f″`(`S₁.trans S₂` は `f` から `f″` への舞台)。 -/
def trans (S₁ S₂ : TwistedLT A pp ff) : TwistedLT A pp ff where
  ϕ := S₁.ϕ
  π := S₁.π
  π' := S₂.π'
  hπmem := S₁.hπmem
  hπ'max := S₂.hπ'max
  hϕres := S₁.hϕres
  hsolve := S₁.hsolve
  f := S₁.f
  hf0 := S₁.hf0
  hf1 := S₁.hf1
  hfres := S₁.hfres
  f' := S₂.f'
  hf'0 := S₂.hf'0
  hf'1 := S₂.hf'1
  hf'res := S₂.hf'res

end TwistedLT

/-! ## 7. 次数ごとの 1 ステップ -/

/-- ★★★★★**ねじれ版の 1 ステップ**。障害が次数 `≤n`(`n ≠ 0`)で消えているとき、
スカラー `c` を足して次数 `≤n+1` まで消せる。

`c` は次の 3 つで決まる:
可除性(`map_twisted_obstruction_eq_zero` + `exists_scalar_dvd_of_map_residue_eq_zero`)で
障害の `n+1` 次係数を `π′δ` と書き、`π^{n+1} = π′u`(`u ∈ 𝔪`、`n ≥ 1` だから)と
分解して、半線型方程式 `c − u·ϕ(c) = −δ` を `hsolve` で解く。
★木の untwisted 版が「`1 − π^n` が単数」で割っていたところが、ここでは `hsolve`。 -/
theorem exists_next_step_twisted {A : Type*} [CommRing A] [IsLocalRing A]
    {pp ff : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp]
    (ϕ : A →+* A) {π π' : A}
    (hπ'max : IsLocalRing.maximalIdeal A = Ideal.span {π'})
    (hπmem : π ∈ IsLocalRing.maximalIdeal A)
    (hϕres : ∀ a : A, IsLocalRing.residue A (ϕ a) = IsLocalRing.residue A a ^ pp ^ ff)
    (hsolve : ∀ u ∈ IsLocalRing.maximalIdeal A, ∀ b : A, ∃ c : A, c - u * ϕ c = b)
    (f : PowerSeries A) (hf0 : PowerSeries.constantCoeff f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ pp ^ ff)
    (f' : PowerSeries A) (hf'0 : PowerSeries.coeff 0 f' = 0) (hf'1 : PowerSeries.coeff 1 f' = π')
    (hf'res : PowerSeries.map (IsLocalRing.residue A) f' = PowerSeries.X ^ pp ^ ff)
    {φ : PowerSeries A} (hφ0 : PowerSeries.constantCoeff φ = 0)
    {n : ℕ} (hn : n ≠ 0)
    (hinv : ∀ k ≤ n, PowerSeries.coeff k
      (PowerSeries.subst φ f' - PowerSeries.subst f (PowerSeries.map ϕ φ)) = 0) :
    ∃ c : A,
      PowerSeries.constantCoeff (φ + c • PowerSeries.X ^ (n + 1)) = 0 ∧
      ∀ k ≤ n + 1, PowerSeries.coeff k
        (PowerSeries.subst (φ + c • PowerSeries.X ^ (n + 1)) f' -
          PowerSeries.subst f (PowerSeries.map ϕ (φ + c • PowerSeries.X ^ (n + 1)))) = 0 := by
  have hRdvd := map_twisted_obstruction_eq_zero (IsLocalRing.residue A) ϕ hϕres f f' hf0 hfres
    hf'res φ hφ0
  obtain ⟨c₀, hc₀⟩ := exists_scalar_dvd_of_map_residue_eq_zero hπ'max hRdvd
  set δ := PowerSeries.coeff (n + 1) c₀ with hδ_def
  have hπ'mem : π' ∈ IsLocalRing.maximalIdeal A := hπ'max ▸ Ideal.mem_span_singleton_self π'
  obtain ⟨s, hs⟩ : ∃ s : A, π = π' * s := by
    have hmem : π ∈ Ideal.span ({π'} : Set A) := hπ'max ▸ hπmem
    obtain ⟨s, hs⟩ := Ideal.mem_span_singleton'.mp hmem
    exact ⟨s, by rw [← hs]; ring⟩
  set u : A := π' ^ n * s ^ (n + 1) with hu_def
  have humem : u ∈ IsLocalRing.maximalIdeal A := by
    have hpn : π' ^ n ∈ IsLocalRing.maximalIdeal A := by
      have hn1 : n = (n - 1) + 1 := by omega
      rw [hn1, pow_succ]
      exact Ideal.mul_mem_left _ _ hπ'mem
    exact Ideal.mul_mem_right _ _ hpn
  have hπpow : π ^ (n + 1) = π' * u := by
    rw [hs, hu_def, mul_pow]; ring
  obtain ⟨c, hc⟩ := hsolve u humem (-δ)
  have hnew0 : PowerSeries.constantCoeff (c • (PowerSeries.X : PowerSeries A) ^ (n + 1)) = 0 :=
    constantCoeff_smul_X_pow c (Nat.succ_ne_zero n)
  refine ⟨c, ?_, ?_⟩
  · rw [map_add, hφ0, zero_add]; exact hnew0
  · intro k hk
    have hforder : (1 : ℕ∞) ≤ MvPowerSeries.order (φ : PowerSeries A) :=
      MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hφ0
    have hneworder : ((n + 1 : ℕ) : ℕ∞) ≤
        MvPowerSeries.order (c • (PowerSeries.X : PowerSeries A) ^ (n + 1)) :=
      (order_X_pow_eq (A := A) (n + 1)) ▸ MvPowerSeries.le_order_smul
    have hflin := coeff_subst_linearize_1var hφ0 hnew0 hforder hneworder (Nat.succ_pos n) f' π'
      hf'0 hf'1 k (by exact_mod_cast hk)
    have hfnew : PowerSeries.coeff k (c • (PowerSeries.X : PowerSeries A) ^ (n + 1)) =
        if k = n + 1 then c else 0 := by
      rw [PowerSeries.coeff_smul, PowerSeries.coeff_X_pow]
      split_ifs with h <;> simp
    have hHSf : PowerSeries.HasSubst f := by
      show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0]; exact IsNilpotent.zero
    have hgsplit : PowerSeries.subst f (PowerSeries.map ϕ (φ + c • PowerSeries.X ^ (n + 1))) =
        PowerSeries.subst f (PowerSeries.map ϕ φ) + ϕ c • f ^ (n + 1) := by
      rw [map_add, map_smul_X_pow, PowerSeries.subst_add hHSf, PowerSeries.subst_smul hHSf,
        PowerSeries.subst_pow hHSf, PowerSeries.subst_X hHSf]
    have hgnew : PowerSeries.coeff k (ϕ c • f ^ (n + 1)) =
        if k = n + 1 then ϕ c * π ^ (n + 1) else 0 := by
      rw [PowerSeries.coeff_smul]
      split_ifs with h
      · rw [h, coeff_pow_self_eq_pow hf0 hf1 n, smul_eq_mul]
      · rw [coeff_lt_pow_self_eq_zero hf0 n (by omega), smul_zero]
    have hflin' : PowerSeries.coeff k (PowerSeries.subst (φ + c • PowerSeries.X ^ (n + 1)) f') =
        PowerSeries.coeff k (PowerSeries.subst φ f') +
          π' * PowerSeries.coeff k (c • PowerSeries.X ^ (n + 1)) := by
      linear_combination hflin
    rw [map_sub, hgsplit, map_add, hflin', hfnew, hgnew]
    by_cases hkn1 : k = n + 1
    · subst hkn1
      simp only [if_true]
      have hc0eq : PowerSeries.coeff (n + 1)
          (PowerSeries.subst φ f' - PowerSeries.subst f (PowerSeries.map ϕ φ))
          = PowerSeries.coeff (n + 1) (π' • c₀) := by rw [hc₀]
      rw [map_sub, PowerSeries.coeff_smul] at hc0eq
      have hobs' : PowerSeries.coeff (n + 1) (PowerSeries.subst φ f') -
          PowerSeries.coeff (n + 1) (PowerSeries.subst f (PowerSeries.map ϕ φ)) = π' * δ := by
        rw [hc0eq, hδ_def, smul_eq_mul]
      linear_combination hobs' + π' * hc - ϕ c * hπpow
    · simp only [if_neg hkn1]
      have hobs0 := hinv k (by omega)
      rw [map_sub] at hobs0
      linear_combination hobs0

/-! ## 8. 近似列 -/

namespace TwistedLT

variable {A : Type*} [CommRing A] [IsLocalRing A] {pp ff : ℕ}
  [ExpChar (IsLocalRing.ResidueField A) pp]

/-- ★★**ねじれ版の次数ごとの近似列**。`φSeq k` は障害が次数 `≤k+1` まで消えている
定数項 `0` の冪級数。出発点は `θ • X`。 -/
noncomputable def φSeq (S : TwistedLT A pp ff) (θ : A) (hθ : S.π' * θ = S.π * S.ϕ θ) :
    (k : ℕ) → {φ : PowerSeries A // S.SeqInvariant φ k} :=
  fun k => Nat.rec
    (⟨θ • (PowerSeries.X : PowerSeries A),
        fun j hj => base_case_twisted S.ϕ hθ S.f S.hf0 S.hf1 S.f' S.coeff_zero_f' S.hf'1 j hj, by
          have h := constantCoeff_smul_X_pow (A := A) θ (m := 1) one_ne_zero
          rwa [pow_one] at h⟩ :
      {φ : PowerSeries A // S.SeqInvariant φ 0})
    (fun k prev =>
      let cex := exists_next_step_twisted S.ϕ S.hπ'max S.hπmem S.hϕres S.hsolve
        S.f S.hf0 S.hf1 S.hfres S.f' S.coeff_zero_f' S.hf'1 S.hf'res prev.2.2
        (n := k + 1) (Nat.succ_ne_zero k) prev.2.1
      (⟨prev.1 + cex.choose • (PowerSeries.X : PowerSeries A) ^ (k + 1 + 1),
          cex.choose_spec.2, cex.choose_spec.1⟩ :
        {φ : PowerSeries A // S.SeqInvariant φ (k + 1)}))
    k

/-- `φSeq (k+1)` は `φSeq k` にスカラー倍の `X^{k+2}` を足したもの。 -/
theorem φSeq_succ_eq (S : TwistedLT A pp ff) (θ : A) (hθ : S.π' * θ = S.π * S.ϕ θ) (k : ℕ) :
    ∃ c : A, (S.φSeq θ hθ (k + 1)).1 =
      (S.φSeq θ hθ k).1 + c • (PowerSeries.X : PowerSeries A) ^ (k + 2) := by
  set prev := S.φSeq θ hθ k with hprev
  set cex := exists_next_step_twisted S.ϕ S.hπ'max S.hπmem S.hϕres S.hsolve
    S.f S.hf0 S.hf1 S.hfres S.f' S.coeff_zero_f' S.hf'1 S.hf'res prev.2.2
    (n := k + 1) (Nat.succ_ne_zero k) prev.2.1 with hcex
  exact ⟨cex.choose, rfl⟩

/-- 近似列の差の次数は離れているほど高い。 -/
theorem order_diff_φSeq_ge (S : TwistedLT A pp ff) (θ : A) (hθ : S.π' * θ = S.π * S.ϕ θ)
    (k : ℕ) : ∀ m, k ≤ m →
    ((k + 2 : ℕ) : ℕ∞) ≤
      MvPowerSeries.order ((S.φSeq θ hθ m).1 - (S.φSeq θ hθ k).1 : PowerSeries A) := by
  intro m hkm
  induction m, hkm using Nat.le_induction with
  | base => simp [MvPowerSeries.order_zero]
  | succ n hn ih =>
    obtain ⟨c, hceq⟩ := S.φSeq_succ_eq θ hθ n
    have hstep : (S.φSeq θ hθ (n + 1)).1 - (S.φSeq θ hθ k).1 =
        ((S.φSeq θ hθ n).1 - (S.φSeq θ hθ k).1) +
          c • (PowerSeries.X : PowerSeries A) ^ (n + 2) := by
      rw [hceq]; ring
    rw [hstep]
    have hcorder : ((n + 2 : ℕ) : ℕ∞) ≤
        MvPowerSeries.order (c • (PowerSeries.X : PowerSeries A) ^ (n + 2)) :=
      (order_X_pow_eq (A := A) (n + 2)) ▸ MvPowerSeries.le_order_smul
    have hnk : ((k + 2 : ℕ) : ℕ∞) ≤ ((n + 2 : ℕ) : ℕ∞) := by
      exact_mod_cast (by omega : k + 2 ≤ n + 2)
    have hmin : ((k + 2 : ℕ) : ℕ∞) ≤
        min (MvPowerSeries.order ((S.φSeq θ hθ n).1 - (S.φSeq θ hθ k).1 : PowerSeries A))
          (MvPowerSeries.order (c • (PowerSeries.X : PowerSeries A) ^ (n + 2))) := by
      rw [le_min_iff]; exact ⟨ih, le_trans hnk hcorder⟩
    exact le_trans hmin MvPowerSeries.min_order_le_add

/-- 係数の安定性。 -/
theorem coeff_φSeq_stable (S : TwistedLT A pp ff) (θ : A) (hθ : S.π' * θ = S.π * S.ϕ θ)
    (e k m : ℕ) (hkm : k ≤ m) (he : e ≤ k + 1) :
    PowerSeries.coeff e (S.φSeq θ hθ m).1 = PowerSeries.coeff e (S.φSeq θ hθ k).1 := by
  have hord := S.order_diff_φSeq_ge θ hθ k m hkm
  have hlt : ((e : ℕ) : ℕ∞) <
      PowerSeries.order ((S.φSeq θ hθ m).1 - (S.φSeq θ hθ k).1 : PowerSeries A) := by
    rw [PowerSeries.order_eq_order]
    calc ((e : ℕ) : ℕ∞) ≤ ((k + 1 : ℕ) : ℕ∞) := by exact_mod_cast he
      _ < ((k + 2 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : k + 1 < k + 2)
      _ ≤ _ := hord
  have hz := PowerSeries.coeff_of_lt_order e hlt
  rw [map_sub] at hz
  exact sub_eq_zero.mp hz

end TwistedLT

/-! ## 9. 極限と関数等式 -/

section Limit

variable {A : Type*} [CommRing A] [IsLocalRing A] {pp ff : ℕ}
  [ExpChar (IsLocalRing.ResidueField A) pp]

/-- ★★★★★★★★★**`𝒪_L` 係数の `[θ]_{f,f′}`**(原典 Proposition 3.5(ii) の冪級数)。

`TwistedLT.φSeq` の極限——次数 `n` の係数は `n` 次まで進んだ近似から読み取る。 -/
noncomputable def LubinTateEndoTwisted (S : TwistedLT A pp ff) (θ : A)
    (hθ : S.π' * θ = S.π * S.ϕ θ) : PowerSeries A :=
  PowerSeries.mk (fun n => PowerSeries.coeff n (S.φSeq θ hθ n).1)

theorem coeff_LubinTateEndoTwisted (S : TwistedLT A pp ff) (θ : A)
    (hθ : S.π' * θ = S.π * S.ϕ θ) (n : ℕ) :
    PowerSeries.coeff n (LubinTateEndoTwisted S θ hθ) =
      PowerSeries.coeff n (S.φSeq θ hθ n).1 :=
  PowerSeries.coeff_mk n _

theorem coeff_LubinTateEndoTwisted_eq_φSeq (S : TwistedLT A pp ff) (θ : A)
    (hθ : S.π' * θ = S.π * S.ϕ θ) (n m : ℕ) (hn : n ≤ m) :
    PowerSeries.coeff n (LubinTateEndoTwisted S θ hθ) =
      PowerSeries.coeff n (S.φSeq θ hθ m).1 := by
  rw [coeff_LubinTateEndoTwisted]
  exact (S.coeff_φSeq_stable θ hθ n n m hn (by omega)).symm

theorem order_diff_LubinTateEndoTwisted_φSeq_ge (S : TwistedLT A pp ff) (θ : A)
    (hθ : S.π' * θ = S.π * S.ϕ θ) (m : ℕ) :
    ((m + 1 : ℕ) : ℕ∞) ≤
      MvPowerSeries.order (LubinTateEndoTwisted S θ hθ - (S.φSeq θ hθ m).1 : PowerSeries A) := by
  rw [← PowerSeries.order_eq_order]
  apply PowerSeries.nat_le_order
  intro n hn
  have hnm : n ≤ m := by omega
  have heq := coeff_LubinTateEndoTwisted_eq_φSeq S θ hθ n m hnm
  rw [map_sub, heq, sub_self]

theorem constantCoeff_LubinTateEndoTwisted (S : TwistedLT A pp ff) (θ : A)
    (hθ : S.π' * θ = S.π * S.ϕ θ) :
    PowerSeries.constantCoeff (LubinTateEndoTwisted S θ hθ) = 0 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply, coeff_LubinTateEndoTwisted,
    PowerSeries.coeff_zero_eq_constantCoeff_apply]
  exact (S.φSeq θ hθ 0).2.2

/-- ★原典 (ii) の第 1 条件 `[θ]_{f,f′}(X) ≡ θX (mod deg 2)`。 -/
theorem coeff_one_LubinTateEndoTwisted (S : TwistedLT A pp ff) (θ : A)
    (hθ : S.π' * θ = S.π * S.ϕ θ) :
    PowerSeries.coeff 1 (LubinTateEndoTwisted S θ hθ) = θ := by
  rw [coeff_LubinTateEndoTwisted]
  obtain ⟨c, hceq⟩ := S.φSeq_succ_eq θ hθ 0
  show PowerSeries.coeff 1 (S.φSeq θ hθ 1).1 = θ
  rw [hceq]
  have hbase : (S.φSeq θ hθ 0).1 = θ • (PowerSeries.X : PowerSeries A) := rfl
  rw [hbase, map_add, PowerSeries.coeff_smul, PowerSeries.coeff_X, if_pos rfl, smul_eq_mul, mul_one]
  rw [PowerSeries.coeff_smul, PowerSeries.coeff_X_pow, if_neg (by omega), smul_zero, add_zero]

/-- `f′` 側の安定性。 -/
theorem coeff_subst_LubinTateEndoTwisted_eq_φSeq (S : TwistedLT A pp ff) (θ : A)
    (hθ : S.π' * θ = S.π * S.ϕ θ) (n m : ℕ) (hn : n ≤ m) :
    PowerSeries.coeff n (PowerSeries.subst (LubinTateEndoTwisted S θ hθ) S.f') =
      PowerSeries.coeff n (PowerSeries.subst (S.φSeq θ hθ m).1 S.f') := by
  set φm := (S.φSeq θ hθ m).1 with hφm_def
  set δ := LubinTateEndoTwisted S θ hθ - φm with hδ_def
  have hφmadd : φm + δ = LubinTateEndoTwisted S θ hθ := by rw [hδ_def]; ring
  have hφm0 : PowerSeries.constantCoeff φm = 0 := (S.φSeq θ hθ m).2.2
  have hφmorder : (1 : ℕ∞) ≤ MvPowerSeries.order (φm : PowerSeries A) :=
    MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hφm0
  have hδorder : ((m + 1 : ℕ) : ℕ∞) ≤ MvPowerSeries.order (δ : PowerSeries A) :=
    order_diff_LubinTateEndoTwisted_φSeq_ge S θ hθ m
  have hδ0 : PowerSeries.constantCoeff δ = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]
    apply PowerSeries.coeff_of_lt_order
    rw [PowerSeries.order_eq_order]
    calc ((0 : ℕ) : ℕ∞) < ((m + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : 0 < m + 1)
      _ ≤ _ := hδorder
  have hlin := coeff_subst_linearize_1var hφm0 hδ0 hφmorder hδorder (by omega : 1 ≤ m + 1)
    S.f' S.π' S.coeff_zero_f' S.hf'1 n (by exact_mod_cast (by omega : n ≤ m + 1))
  rw [hφmadd] at hlin
  have hδn0 : PowerSeries.coeff n δ = 0 := by
    apply PowerSeries.coeff_of_lt_order
    rw [PowerSeries.order_eq_order]
    calc ((n : ℕ) : ℕ∞) ≤ ((m : ℕ) : ℕ∞) := by exact_mod_cast hn
      _ < ((m + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : m < m + 1)
      _ ≤ _ := hδorder
  rw [hδn0, mul_zero] at hlin
  exact sub_eq_zero.mp hlin

/-- ねじれ側(`φ^ϕ ∘ f`)の安定性。 -/
theorem coeff_twistSubst_LubinTateEndoTwisted_eq_φSeq (S : TwistedLT A pp ff) (θ : A)
    (hθ : S.π' * θ = S.π * S.ϕ θ) (n m : ℕ) (hn : n ≤ m) :
    PowerSeries.coeff n
        (PowerSeries.subst S.f (PowerSeries.map S.ϕ (LubinTateEndoTwisted S θ hθ))) =
      PowerSeries.coeff n (PowerSeries.subst S.f (PowerSeries.map S.ϕ (S.φSeq θ hθ m).1)) := by
  set φm := (S.φSeq θ hθ m).1 with hφm_def
  set δ := LubinTateEndoTwisted S θ hθ - φm with hδ_def
  have hδorder : ((m + 1 : ℕ) : ℕ∞) ≤ MvPowerSeries.order (δ : PowerSeries A) :=
    order_diff_LubinTateEndoTwisted_φSeq_ge S θ hθ m
  have hmapδorder : ((m + 1 : ℕ) : ℕ∞) ≤
      MvPowerSeries.order (PowerSeries.map S.ϕ δ : PowerSeries A) :=
    nat_le_order_map S.ϕ hδorder
  have hmapsub : PowerSeries.map S.ϕ (LubinTateEndoTwisted S θ hθ) - PowerSeries.map S.ϕ φm =
      PowerSeries.map S.ϕ δ := by
    rw [hδ_def, map_sub]
  have hsub : PowerSeries.subst S.f (PowerSeries.map S.ϕ (LubinTateEndoTwisted S θ hθ)) -
      PowerSeries.subst S.f (PowerSeries.map S.ϕ φm) =
      PowerSeries.subst S.f (PowerSeries.map S.ϕ δ) := by
    rw [← hmapsub, PowerSeries.subst_sub S.hasSubst_f]
  have horder : ((m + 1 : ℕ) : ℕ∞) ≤
      MvPowerSeries.order (PowerSeries.subst S.f (PowerSeries.map S.ϕ δ) : PowerSeries A) := by
    have hle := PowerSeries.le_order_subst (S.f : MvPowerSeries Unit A) S.hasSubst_f
      (PowerSeries.map S.ϕ δ)
    have hforder : (1 : ℕ∞) ≤ MvPowerSeries.order (S.f : PowerSeries A) :=
      MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr S.hf0
    have hmapδorder' : ((m + 1 : ℕ) : ℕ∞) ≤
        PowerSeries.order (PowerSeries.map S.ϕ δ : PowerSeries A) := by
      rw [PowerSeries.order_eq_order]; exact hmapδorder
    calc ((m + 1 : ℕ) : ℕ∞) = 1 * ((m + 1 : ℕ) : ℕ∞) := by ring
      _ ≤ MvPowerSeries.order (S.f : PowerSeries A) * ((m + 1 : ℕ) : ℕ∞) := by gcongr
      _ ≤ MvPowerSeries.order (S.f : PowerSeries A) *
            PowerSeries.order (PowerSeries.map S.ϕ δ : PowerSeries A) := by gcongr
      _ ≤ _ := hle
  have hz : PowerSeries.coeff n (PowerSeries.subst S.f (PowerSeries.map S.ϕ δ)) = 0 := by
    apply PowerSeries.coeff_of_lt_order
    rw [PowerSeries.order_eq_order]
    calc ((n : ℕ) : ℕ∞) ≤ ((m : ℕ) : ℕ∞) := by exact_mod_cast hn
      _ < ((m + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : m < m + 1)
      _ ≤ _ := horder
  rw [← hsub, map_sub] at hz
  exact sub_eq_zero.mp hz

/-- ★★★★★★★★★**ねじれ版の関数等式** `f′ ∘ [θ]_{f,f′} = [θ]^ϕ_{f,f′} ∘ f`
(原典 Proposition 3.5(ii) の第 2 条件)。

★代入の記法では `PowerSeries.subst a b = b ∘ a` なので、左辺は `subst [θ] f′`、
右辺は `subst f (map ϕ [θ])` である。 -/
theorem LubinTateEndoTwisted_functional_equation (S : TwistedLT A pp ff) (θ : A)
    (hθ : S.π' * θ = S.π * S.ϕ θ) :
    PowerSeries.subst (LubinTateEndoTwisted S θ hθ) S.f' =
      PowerSeries.subst S.f (PowerSeries.map S.ϕ (LubinTateEndoTwisted S θ hθ)) := by
  apply PowerSeries.ext
  intro n
  have h1 := coeff_subst_LubinTateEndoTwisted_eq_φSeq S θ hθ n n le_rfl
  have h2 := coeff_twistSubst_LubinTateEndoTwisted_eq_φSeq S θ hθ n n le_rfl
  have h3 : PowerSeries.coeff n
      (PowerSeries.subst (S.φSeq θ hθ n).1 S.f' -
        PowerSeries.subst S.f (PowerSeries.map S.ϕ (S.φSeq θ hθ n).1)) = 0 :=
    (S.φSeq θ hθ n).2.1 n (by omega)
  rw [map_sub] at h3
  rw [h1, h2]
  exact sub_eq_zero.mp h3

/-! ## 10. 原典 Proposition 3.5(ii) の残り(一意性・合成則)と Corollary 3.7(ii) -/

/-- ★★★★★★★★**原典 Proposition 3.5(ii) の一意性**: 条件
(定数項 `0`・1 次係数 `θ`・関数等式)を満たす冪級数は `[θ]_{f,f′}` に限る。 -/
theorem eq_LubinTateEndoTwisted (S : TwistedLT A pp ff) (θ : A)
    (hθ : S.π' * θ = S.π * S.ϕ θ)
    (hcancel : ∀ m : ℕ, 2 ≤ m → ∀ d : A, S.π' * d = S.ϕ d * S.π ^ m → d = 0)
    {α : PowerSeries A} (hα0 : PowerSeries.constantCoeff α = 0)
    (hα1 : PowerSeries.coeff 1 α = θ)
    (hα : PowerSeries.subst α S.f' = PowerSeries.subst S.f (PowerSeries.map S.ϕ α)) :
    α = LubinTateEndoTwisted S θ hθ :=
  powerSeries_uniqueness_twisted S.ϕ hcancel S.hf0 S.hf1 S.coeff_zero_f' S.hf'1 hα0
    (constantCoeff_LubinTateEndoTwisted S θ hθ)
    (by rw [hα1, coeff_one_LubinTateEndoTwisted]) hα
    (LubinTateEndoTwisted_functional_equation S θ hθ)

/-- ★★★★★★★★★**原典 Proposition 3.5(ii) の合成則(Frobenius ねじれ版)**
`[θ′]_{f′,f″} ∘ [θ]_{f,f′} = [θθ′]_{f,f″}`。

★`π`・`π′`・`π″` は 3 つとも異なってよい。`hθ''` は原典 Definition 3.3 の
`θθ′ ∈ Θ_{π,π″}` であり、`theta_mul_mem` が `hθ`・`hθ′` から供給する
(`S₂.π = S₁.π′` のとき)。 -/
theorem subst_LubinTateEndoTwisted_LubinTateEndoTwisted (S₁ S₂ : TwistedLT A pp ff)
    (hϕ : S₂.ϕ = S₁.ϕ) (hff : S₂.f = S₁.f')
    (hcancel : ∀ m : ℕ, 2 ≤ m → ∀ d : A, S₂.π' * d = S₁.ϕ d * S₁.π ^ m → d = 0)
    (θ θ' : A) (hθ : S₁.π' * θ = S₁.π * S₁.ϕ θ) (hθ' : S₂.π' * θ' = S₂.π * S₂.ϕ θ')
    (hθ'' : S₂.π' * (θ * θ') = S₁.π * S₁.ϕ (θ * θ')) :
    PowerSeries.subst (LubinTateEndoTwisted S₁ θ hθ) (LubinTateEndoTwisted S₂ θ' hθ') =
      LubinTateEndoTwisted (S₁.trans S₂) (θ * θ') hθ'' := by
  have ha := LubinTateEndoTwisted_functional_equation S₁ θ hθ
  have hb := LubinTateEndoTwisted_functional_equation S₂ θ' hθ'
  have hc := LubinTateEndoTwisted_functional_equation (S₁.trans S₂) (θ * θ') hθ''
  rw [hϕ, hff] at hb
  refine subst_comp_eq_of_twisted_intertwine_pair S₁.ϕ hcancel S₁.hf0 S₁.hf1 S₁.hf'0
    S₂.hf'0 S₂.hf'1
    (constantCoeff_LubinTateEndoTwisted S₁ θ hθ)
    (constantCoeff_LubinTateEndoTwisted S₂ θ' hθ')
    (constantCoeff_LubinTateEndoTwisted (S₁.trans S₂) (θ * θ') hθ'')
    ha hb hc ?_
  rw [coeff_one_LubinTateEndoTwisted (S₁.trans S₂) (θ * θ') hθ'',
    coeff_one_LubinTateEndoTwisted S₁ θ hθ, coeff_one_LubinTateEndoTwisted S₂ θ' hθ', mul_comm]

/-- `f = f′` のとき `[1]_{f,f} = X`。 -/
theorem LubinTateEndoTwisted_one_eq_X (S : TwistedLT A pp ff) (hff : S.f = S.f')
    (hcancel : ∀ m : ℕ, 2 ≤ m → ∀ d : A, S.π' * d = S.ϕ d * S.π ^ m → d = 0)
    (h1 : S.π' * 1 = S.π * S.ϕ 1) :
    LubinTateEndoTwisted S 1 h1 = PowerSeries.X := by
  refine (eq_LubinTateEndoTwisted S 1 h1 hcancel PowerSeries.constantCoeff_X
    PowerSeries.coeff_one_X ?_).symm
  rw [PowerSeries.map_X, PowerSeries.X_subst, PowerSeries.subst_X S.hasSubst_f, hff]

/-- ★★★★★★★★**原典 Corollary 3.7(ii) の Frobenius ねじれ版**:
`θ′θ = 1` ならば `[θ′]_{f′,f} ∘ [θ]_{f,f′} = X`(すなわち `[θ]_{f,f′}` は同型で
逆は `[θ^{-1}]_{f′,f}`)。 -/
theorem subst_LubinTateEndoTwisted_eq_X (S₁ S₂ : TwistedLT A pp ff)
    (hϕ : S₂.ϕ = S₁.ϕ) (hff : S₂.f = S₁.f') (hff2 : S₂.f' = S₁.f)
    (hcancel : ∀ m : ℕ, 2 ≤ m → ∀ d : A, S₂.π' * d = S₁.ϕ d * S₁.π ^ m → d = 0)
    (θ θ' : A) (hθ : S₁.π' * θ = S₁.π * S₁.ϕ θ) (hθ' : S₂.π' * θ' = S₂.π * S₂.ϕ θ')
    (hθθ' : θ' * θ = 1) :
    PowerSeries.subst (LubinTateEndoTwisted S₁ θ hθ) (LubinTateEndoTwisted S₂ θ' hθ') =
      PowerSeries.X := by
  have ha := LubinTateEndoTwisted_functional_equation S₁ θ hθ
  have hb := LubinTateEndoTwisted_functional_equation S₂ θ' hθ'
  rw [hϕ, hff] at hb
  refine subst_comp_eq_of_twisted_intertwine_pair S₁.ϕ hcancel S₁.hf0 S₁.hf1 S₁.hf'0
    S₂.hf'0 S₂.hf'1
    (constantCoeff_LubinTateEndoTwisted S₁ θ hθ)
    (constantCoeff_LubinTateEndoTwisted S₂ θ' hθ')
    PowerSeries.constantCoeff_X ha hb ?_ ?_
  · rw [PowerSeries.map_X, PowerSeries.X_subst, PowerSeries.subst_X S₁.hasSubst_f, hff2]
  · rw [PowerSeries.coeff_one_X, coeff_one_LubinTateEndoTwisted S₁ θ hθ,
      coeff_one_LubinTateEndoTwisted S₂ θ' hθ', hθθ']

end Limit

/-! ## 11. 退化の自己検査(`ϕ = id`)

★舞台が空でないこと、および本ファイルが木の既存の構成の**真の一般化**である
ことを、`ϕ = id`(`π = π′`, `f = f′`)への退化で確かめる。 -/

/-- ★**退化の自己検査(その1)**: `ϕ = id`・`π = π′`・`f = f′` は本当に舞台の一例になる。

`hϕres` は `FiniteField.pow_card`(剰余体の位数が `q` なら `x^q = x`)から、
`hsolve` は「`u ∈ 𝔪` なら `1 − u` が単数」から出る。 -/
noncomputable def TwistedLT.ofUntwisted {A : Type*} [CommRing A] [IsLocalRing A]
    [Fintype (IsLocalRing.ResidueField A)] {pp ff : ℕ}
    [ExpChar (IsLocalRing.ResidueField A) pp]
    (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.constantCoeff f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ pp ^ ff) :
    TwistedLT A pp ff where
  ϕ := RingHom.id A
  π := π
  π' := π
  hπmem := hπmax ▸ Ideal.mem_span_singleton_self π
  hπ'max := hπmax
  hϕres := fun a => by
    show IsLocalRing.residue A a = IsLocalRing.residue A a ^ pp ^ ff
    rw [← hq]
    exact (FiniteField.pow_card _).symm
  hsolve := fun u hu b => by
    obtain ⟨v, hv⟩ := IsLocalRing.isUnit_one_sub_self_of_mem_nonunits u
      (IsLocalRing.mem_maximalIdeal u |>.mp hu)
    refine ⟨((v⁻¹ : Aˣ) : A) * b, ?_⟩
    show ((v⁻¹ : Aˣ) : A) * b - u * (((v⁻¹ : Aˣ) : A) * b) = b
    have hkey : (1 - u) * (((v⁻¹ : Aˣ) : A) * b) = b := by
      rw [← hv, ← mul_assoc, Units.mul_inv, one_mul]
    linear_combination hkey
  f := f
  hf0 := hf0
  hf1 := hf1
  hfres := hfres
  f' := f
  hf'0 := hf0
  hf'1 := hf1
  hf'res := hfres

/-- ★★**退化の自己検査(その2)**: `ϕ = id`(かつ `π = π′`, `f = f′`)に潰すと、
本ファイルの `LubinTateEndoTwisted` は木の既存の `LubinTateEndo`
(`Found/PGC/LubinTateEndoLimit.lean`)と**同じ冪級数**になる。

★等式そのものを木の一意性補題 `powerSeries_uniqueness` で照合している
——「退化していない」ことの証拠であり、同時に本ファイルが木の構成の
真の一般化であることの証拠でもある。 -/
theorem LubinTateEndoTwisted_ofUntwisted_eq_LubinTateEndo {A : Type*} [CommRing A] [IsLocalRing A]
    [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] {pp ff : ℕ}
    [ExpChar (IsLocalRing.ResidueField A) pp]
    (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.constantCoeff f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ pp ^ ff) (a : A)
    (hθ : (TwistedLT.ofUntwisted hq hπmax f hf0 hf1 hfres).π' * a =
      (TwistedLT.ofUntwisted hq hπmax f hf0 hf1 hfres).π *
        (TwistedLT.ofUntwisted hq hπmax f hf0 hf1 hfres).ϕ a) :
    LubinTateEndoTwisted (TwistedLT.ofUntwisted hq hπmax f hf0 hf1 hfres) a hθ =
      LubinTateEndo hq hπmax f hf0 hf1 hfres f
        ((PowerSeries.coeff_zero_eq_constantCoeff_apply f).trans hf0) hf1 hfres a := by
  set S := TwistedLT.ofUntwisted hq hπmax f hf0 hf1 hfres with hS
  have hf0' : PowerSeries.coeff 0 f = 0 :=
    (PowerSeries.coeff_zero_eq_constantCoeff_apply f).trans hf0
  set α := LubinTateEndoTwisted S a hθ with hα
  set β := LubinTateEndo hq hπmax f hf0 hf1 hfres f hf0' hf1 hfres a with hβ
  have hαcomm : PowerSeries.subst f α = PowerSeries.subst α f := by
    have h := LubinTateEndoTwisted_functional_equation S a hθ
    have hSf' : S.f' = f := rfl
    have hSf : S.f = f := rfl
    have hSϕ : PowerSeries.map S.ϕ α = α := by
      show PowerSeries.map (RingHom.id A) α = α
      simp
    rw [hSf', hSf, hSϕ] at h
    exact h.symm
  have hβcomm : PowerSeries.subst f β = PowerSeries.subst β f :=
    (LubinTateEndo_functional_equation hq hπmax f hf0 hf1 hfres f hf0' hf1 hfres a).symm
  refine powerSeries_uniqueness hπmax hπne0 hf0 hf1
    (constantCoeff_LubinTateEndoTwisted S a hθ)
    (constantCoeff_LubinTateEndo hq hπmax f hf0 hf1 hfres f hf0' hf1 hfres a) ?_ hαcomm hβcomm
  rw [coeff_one_LubinTateEndoTwisted S a hθ,
    coeff_one_LubinTateEndo hq hπmax f hf0 hf1 hfres f hf0' hf1 hfres a]

end ABC3.Found.PGC
