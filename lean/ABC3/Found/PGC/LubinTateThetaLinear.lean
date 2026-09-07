import ABC3.Found.PGC.LubinTateReciprocityIndependence

/-!
# `[θ]_{f,f′}` の `𝒪`-線型性と合成則(Yoshida 2008 Proposition 3.5(ii)(iii) の対版)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Proposition 3.5(物理 p.5)
および Corollary 3.7(ii)(物理 p.6)。構造化済み原文は
`ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-3.html` の `#prop-3-5` / `#cor-3-7`。

原文 (Yoshida08 p.5, Proposition 3.5(ii)(iii)) の逐語:

> (ii) There is a unique map [·]f,f′ : ΘL π,π′ → (X) ⊂ OL[[X]] such that:
> [θ]f,f′(X) ≡ θX (mod deg 2), f ′ ◦ [θ]f,f′ = [θ]ϕ f,f′ ◦ f.
> It satisfies [θ]f,f′ +Ff′ [θ′]f,f′ = [θ + θ′]f,f′, [θ′]f′,f′′ ◦ [θ]f,f′ = [θθ′]f,f′′.
> (iii) We have [θ]f,f′ ∈ HomOL(Ff, Ff′) for all θ ∈ ΘL π,π′.

原文 (Yoshida08 p.6, Corollary 3.7(ii)) の逐語:

> (ii) If θ ∈ ΘL,× π,π′ := ΘL π,π′ ∩ O× L , then [θ]f,f′ is an isomorphism
> with the inverse [θ−1]f′,f.

本ファイルが閉じるのは (ii) の**合成則** `[θ′] ∘ [θ] = [θθ′]`、そこから出る
(iii) の `𝒪`-線型性(`[θ] ∘ [a]_f = [a]_{f′} ∘ [θ]`)、および Corollary 3.7(ii)
(`[θ]_{f,f′}` の可逆性)である。(ii) の**加法性** `[θ] +_{F_{f′}} [θ′] = [θ+θ′]`
は 2 変数の一意性を要するので本ファイルの担当ではない(下記「新しく必要なノード」)。

## 消費先

`Found/PGC/LubinTateReciprocityIndependence.lean`(Corollary 4.9 後半)が
具体層を「仮定」から「木の実物」へ置き換えるために名指しした 3 つのうちの
**1 番目**である。残り 2 つ(完備化レベルでの半線型性 `σ([θ](α)) = [θ^{(j)}](σα)`、
`ρ_{f,m}` そのものの完備化レベルでの定義)は本ファイルには入っていない。

## 何を足したか

### 純抽象核(冪級数も局所環も出てこない —— 型 1 つと二項演算 1 つだけ)

| 宣言 | 内容 |
|---|---|
| `comp_intertwine_left` | `η∘f = g∘η`・`f∘α = α∘g` ならば `g∘(η∘α) = (η∘α)∘g` |
| `comp_cancel_of_left_inv` | `θ∘η = e`・`e∘x = x`・`η∘α = η∘β` ならば `α = β` |
| `comp_intertwine_comp` | `v∘a = a∘u`・`w∘b = b∘v` ならば `w∘(b∘a) = (b∘a)∘u` |

★どれも仮定は「結合律 `hassoc`」だけである。群でもモノイドでもなく、
`comp : M → M → M` と `hassoc : ∀ x y z, comp (comp x y) z = comp x (comp y z)`
しか要らない(`e` を使う 2 本目だけ `hid : ∀ x, comp e x = x` を足す)。
★`#print axioms` は 3 本とも `does not depend on any axioms`。

### 代入モノイド(一般の可換環 `A` 上、定数項 `0` の冪級数)

| 宣言 | 内容 |
|---|---|
| `SubstNilp` | `{p : PowerSeries A // constantCoeff p = 0}` |
| `substComp` | `x ∘ y`(= `PowerSeries.subst y.1 x.1`)。この部分型で閉じる |
| `substComp_assoc` | 結合律(mathlib の `PowerSeries.subst_comp_subst_apply` 1 本) |
| `substOne` / `substOne_comp` | 単位元 `X` と `X ∘ x = x` |

★これで純抽象核の `hassoc` が**無条件**に供給できる。`PowerSeries A` のままだと
`PowerSeries.subst_comp_subst_apply` が `HasSubst` を要求するので `hassoc` が立たない。

### 抽象層(一般の可換環上の冪級数。分岐・付値・Galois・Lubin-Tate の語彙は出ない)

| 宣言 | 内容 |
|---|---|
| `hasSubst_of_constantCoeff_zero` | 定数項 `0` ならば `HasSubst` |
| `subst_comp_intertwine_left` | 純抽象核 1 の冪級数版 |
| `subst_cancel_left` | 純抽象核 2 の冪級数版 |
| `subst_intertwine_comp` | 純抽象核 3 の冪級数版 |

### 二面版一意性(局所環・整域・`𝔪 = (π)` を要する。それでも分岐の語彙は出ない)

| 宣言 | 内容 |
|---|---|
| `powerSeries_uniqueness_pair` | `f ≠ g` を許す一意性: `f∘α = α∘g`・`f∘β = β∘g`・`α′(0) = β′(0)` ならば `α = β` |
| `subst_comm_of_intertwine_pair` | `θ∘ψ = φ∘θ`(`𝒪`-線型性そのもの) |
| `subst_comp_eq_of_intertwine_pair` | `b∘a = c`(合成則そのもの) |

### 具体層(木の Lubin-Tate 機械へ代入するだけ)

| 宣言 | 内容 |
|---|---|
| `subst_lubinTateEndo_lubinTateEndo` | `[b]_{v,w} ∘ [a]_{u,v} = [ab]_{u,w}`(Prop 3.5(ii) 合成則) |
| `subst_lubinTateEndo_lubinTateEndo_eq_X` | `ab = 1` ならば `[a]_{u,v} ∘ [b]_{v,u} = X`(Cor 3.7(ii)) |
| `subst_lubinTateAction_lubinTateEndo` | `[θ]_{g,f} ∘ [a]_g = [a]_f ∘ [θ]_{g,f}`(Prop 3.5(iii) の `𝒪`-線型性) |

## ★段取りとの差分 —— 「二面版は要らない」は半分だけ当たっていた

持ち場の見立ては

> `γ := η ∘ [a]_f ∘ θ` が `g` と可換で `coeff 1 γ = a` を示せば、
> 既存の単一 `h` 版で `γ = [a]_g` が出る。二面版は要らない。

だった。★**「単一 `h` 版だけで足りる」は本当**である。しかし `γ` を 1 回こしらえる
代わりに、**同じ左簡約を 1 回だけ切り出して二面版一意性そのものを作る**方が安い:

* `η ∘ α` と `η ∘ β` はどちらも `g` と可換(`comp_intertwine_left`)で 1 次係数が等しい
  ので、単一 `h` 版で `η ∘ α = η ∘ β`。
* `θ ∘ η = X` で左から簡約して `α = β`(`comp_cancel_of_left_inv`)。

★これが `powerSeries_uniqueness_pair` である。**二面版は「要らない」のではなく
「単一 `h` 版から 5 行で出る」**。しかも一度出してしまえば、`𝒪`-線型性だけでなく
**合成則 `[θ′]∘[θ] = [θθ′]` と Corollary 3.7(ii) の可逆性も同じ 1 本から落ちる**
——`γ` を作る道では合成則までは届かない。★原典が Lemma 3.4 を二面版で述べているのは
正しく、木に単一 `h` 版しか無かったのが差分だった。

★`aeval_subst_eq_aeval_aeval`(持ち場が「取り回しが本体」と見立てた補題)は
**1 度も使っていない**。それは木の `Found/PGC/AdjoinIntegers.lean` にある
「点で評価するときの連鎖律」であって、冪級数どうしの合成の結合律ではない。
後者は mathlib の `PowerSeries.subst_comp_subst_apply` であり、しかも
`substComp_assoc` に 1 度閉じ込めたので、以後の証明に結合律は 1 度も現れない。

## 逸脱の記録

1. ★**Frobenius ねじれ `[θ]^ϕ` を落とした**。原典は `𝒪_L` 係数(`L/K` 不分岐)で
   `f′ ◦ [θ]_{f,f′} = [θ]^ϕ_{f,f′} ◦ f` と述べる。本ファイルは木の
   `LubinTateEndo`(`Found/PGC/LubinTateEndoLimit.lean`)に合わせて `ϕ = id` の場合、
   すなわち `f′ ◦ [θ] = [θ] ◦ f` だけを扱う。★これが持ち場の言う「対版」である。
   完備化レベルの半線型版(`σ([θ](α)) = [θ^{(j)}](σα)`)は本ファイルに入っていない。
2. ★**`π = π′` に限る**。原典は `Θ^L_{π,π′}`(`θπ = π′θ^ϕ`)で `π ≠ π′` を許すが、
   木の `LubinTateEndo` は `f` と `g` が同じ `π` を持つ場合にしか定義されていない
   (`hf1 : coeff 1 f = π` と `hg1 : coeff 1 g = π` の両方を取る)。
   したがって本ファイルの `[a]_{u,v}` は `a ∈ 𝒪 = Θ^L_{π,π}` の場合である。
3. ★**`θ` の可逆性の出どころ**: `subst_lubinTateEndo_lubinTateEndo_eq_X` は
   合成則 + `[1]_f = X`(木の `LubinTateAction_one_eq_X`)から出している。
   一方、合成則の証明の中で使う `θ ∘ η = X`(`θ = [1]_{u,w}`・`η = [1]_{w,u}`)は
   木の `Found/PGC/LubinTateFieldFIndependent.lean::subst_eq_X_of_intertwine_pair`
   から取っている(そちらは単一 `h` 版の一意性を `h := w` で使う)。★循環はしていない
   ——前者は `a`・`b` が一般、後者は `a = b = 1` の場合だけで、後者が先に立つ。
4. 原典 Prop 3.5(ii) の**加法性** `[θ] +_{F_{f′}} [θ′] = [θ+θ′]_{f,f′}` は入れていない。
   `F_{f′}` は 2 変数(`MvPowerSeries (Fin 2) A`)なので、二面版一意性の 2 変数版
   (木の `mvPowerSeries_uniqueness_general` の二面版)が要る。

## 退化の自己検査

* ★`[a]_f` と `[a]_g` は**別の形式群の間の準同型**である。本ファイルの
  `subst_lubinTateAction_lubinTateEndo` は `LubinTateAction hq hπmax g … a`(= `[a]_g`)と
  `LubinTateAction hq hπmax f … a`(= `[a]_f`)を**別の項として**書いており、
  `f = g` に退化させてはいない。`f = g` と置くと木の `LubinTateAction_comp`
  (単一 `f` 版)の可換性に落ちる。
* ★`η = [1]_{w,u}` の可逆性(`θ ∘ η = X`)は
  `subst_eq_X_of_intertwine_pair`(木、`LubinTateFieldFIndependent.lean:164`)から
  取っている。★仮定として置いたのではない。
* ★1 次係数の計算(`coeff 1 (b ∘ a) = coeff 1 b · coeff 1 a`)は
  `coeff_one_subst_1var`(木、`LubinTateActionMul.lean:37`)で明示的に計算しており、
  これが一意性の入力である(`hlead`)。`a * b` と `b * a` の入れ替えは
  `mul_comm` を明示的に書いている。
* ★`ℕ∞` の切り詰め引き算・除算は 1 度も書いていない(`order` の比較は
  木の既存補題の内部にしか現れない)。

## 新しく必要になったノード

* Prop 3.5(ii) の**加法性**(`[θ] +_{F_{f′}} [θ′] = [θ+θ′]`)—— 二面版の
  2 変数一意性が要る。
* Prop 3.5 の**Frobenius ねじれ版**(`f′ ◦ [θ] = [θ]^ϕ ◦ f`、`π ≠ π′`)——
  `LubinTateEndo` の `𝒪_L` 係数版そのものが木に無い。
* Corollary 4.9 具体化の残り 2 つ(半線型性・`ρ_{f,m}` の完備化レベル定義)。
-/

namespace ABC3.Found.PGC

def subst_lubinTateEndo_lubinTateEndo.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Proposition 3.5", sectionId := "prop-3-5" }

def subst_lubinTateAction_lubinTateEndo.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Proposition 3.5", sectionId := "prop-3-5" }

def subst_lubinTateEndo_lubinTateEndo_eq_X.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 6, item := "Corollary 3.7", sectionId := "cor-3-7" }

/-! ## 0. 純抽象核

★型 `M` と二項演算 `comp` と結合律しか使わない。冪級数も局所環も分岐も出てこない。
`comp x y` は「`x` のあとに `y` を代入する」= `x ∘ y` と読む。 -/

/-- ★★★★★★★★**純抽象核 1** —— 絡み作用素は左から掛けると可換元を作る。

`η ∘ f = g ∘ η`(`η` は `f` を `g` へ移す)と `f ∘ α = α ∘ g`(`α` は `g` を `f` へ移す)
から `g ∘ (η ∘ α) = (η ∘ α) ∘ g`、すなわち `η ∘ α` は `g` と可換。

仮定は結合律だけ。群でもモノイドでもなくてよい。 -/
theorem comp_intertwine_left {M : Type*} (comp : M → M → M)
    (hassoc : ∀ x y z : M, comp (comp x y) z = comp x (comp y z))
    (f g η α : M) (hη : comp η f = comp g η) (hα : comp f α = comp α g) :
    comp g (comp η α) = comp (comp η α) g := by
  calc comp g (comp η α) = comp (comp g η) α := (hassoc g η α).symm
    _ = comp (comp η f) α := by rw [hη]
    _ = comp η (comp f α) := hassoc η f α
    _ = comp η (comp α g) := by rw [hα]
    _ = comp (comp η α) g := (hassoc η α g).symm

/-- ★★★★★★★★**純抽象核 2** —— 左逆元による簡約。

`θ ∘ η = e` かつ `e` が左単位元なら、`η ∘ α = η ∘ β` から `α = β`。
`θ` が両側逆元である必要も、`e` が右単位元である必要も無い。 -/
theorem comp_cancel_of_left_inv {M : Type*} (comp : M → M → M)
    (hassoc : ∀ x y z : M, comp (comp x y) z = comp x (comp y z))
    (e θ η : M) (hθη : comp θ η = e) (hid : ∀ x : M, comp e x = x)
    (α β : M) (h : comp η α = comp η β) : α = β := by
  calc α = comp e α := (hid α).symm
    _ = comp (comp θ η) α := by rw [hθη]
    _ = comp θ (comp η α) := hassoc θ η α
    _ = comp θ (comp η β) := by rw [h]
    _ = comp (comp θ η) β := (hassoc θ η β).symm
    _ = comp e β := by rw [hθη]
    _ = β := hid β

/-- ★★★★★★★★**純抽象核 3** —— 絡み作用素の合成はまた絡み作用素。

`v ∘ a = a ∘ u`(`a : u → v`)と `w ∘ b = b ∘ v`(`b : v → w`)から
`w ∘ (b ∘ a) = (b ∘ a) ∘ u`(`b ∘ a : u → w`)。仮定は結合律だけ。 -/
theorem comp_intertwine_comp {M : Type*} (comp : M → M → M)
    (hassoc : ∀ x y z : M, comp (comp x y) z = comp x (comp y z))
    (u v w a b : M) (ha : comp v a = comp a u) (hb : comp w b = comp b v) :
    comp w (comp b a) = comp (comp b a) u := by
  calc comp w (comp b a) = comp (comp w b) a := (hassoc w b a).symm
    _ = comp (comp b v) a := by rw [hb]
    _ = comp b (comp v a) := hassoc b v a
    _ = comp b (comp a u) := by rw [ha]
    _ = comp (comp b a) u := (hassoc b a u).symm

/-! ## 1. 代入モノイド —— 定数項 `0` の冪級数

★純抽象核の `hassoc` を**無条件に**供給するための舞台。`PowerSeries A` のままだと
mathlib の `PowerSeries.subst_comp_subst_apply` が `HasSubst` を要求するので
`∀ x y z, …` の形の結合律が立たない。定数項 `0` に制限すれば `HasSubst` は自動である。 -/

/-- 定数項が `0` の冪級数は代入できる(`HasSubst`)。 -/
theorem hasSubst_of_constantCoeff_zero {A : Type*} [CommRing A] {p : PowerSeries A}
    (hp : PowerSeries.constantCoeff p = 0) : PowerSeries.HasSubst p := by
  show IsNilpotent (PowerSeries.constantCoeff p)
  rw [hp]
  exact IsNilpotent.zero

/-- 定数項が `0` の冪級数のなす部分型。代入(合成)で閉じている。 -/
def SubstNilp (A : Type*) [CommRing A] :=
  {p : PowerSeries A // PowerSeries.constantCoeff p = 0}

/-- `SubstNilp` の上の合成 `x ∘ y`(`y` を `x` に代入する)。 -/
noncomputable def substComp {A : Type*} [CommRing A] (x y : SubstNilp A) : SubstNilp A :=
  ⟨PowerSeries.subst y.1 x.1, PowerSeries.constantCoeff_subst_eq_zero y.2 x.1 x.2⟩

theorem substComp_val {A : Type*} [CommRing A] (x y : SubstNilp A) :
    (substComp x y).1 = PowerSeries.subst y.1 x.1 := rfl

/-- ★合成の結合律。mathlib の `PowerSeries.subst_comp_subst_apply` 1 本。
★本ファイルで結合律を使うのはここだけである。 -/
theorem substComp_assoc {A : Type*} [CommRing A] (x y z : SubstNilp A) :
    substComp (substComp x y) z = substComp x (substComp y z) :=
  Subtype.ext (PowerSeries.subst_comp_subst_apply (hasSubst_of_constantCoeff_zero y.2)
    (hasSubst_of_constantCoeff_zero z.2) x.1)

/-- 合成の単位元 `X`。 -/
noncomputable def substOne {A : Type*} [CommRing A] : SubstNilp A :=
  ⟨PowerSeries.X, PowerSeries.constantCoeff_X⟩

theorem substOne_comp {A : Type*} [CommRing A] (x : SubstNilp A) : substComp substOne x = x :=
  Subtype.ext (PowerSeries.subst_X (hasSubst_of_constantCoeff_zero x.2))

/-! ## 2. 抽象層 —— 一般の可換環上の冪級数

★分岐・付値・Galois・Lubin-Tate の語彙は 1 つも出てこない。
`PowerSeries.subst h α = α(h)` の規約で、`subst α f = subst g α` は
`f(α(X)) = α(g(X))`、すなわち `f ∘ α = α ∘ g` を意味する。 -/

/-- 純抽象核 1 の冪級数版: `g ∘ η = η ∘ f`・`f ∘ α = α ∘ g` ならば
`η ∘ α` は `g` と可換。 -/
theorem subst_comp_intertwine_left {A : Type*} [CommRing A] {f g η α : PowerSeries A}
    (hf0 : PowerSeries.constantCoeff f = 0) (hg0 : PowerSeries.constantCoeff g = 0)
    (hη0 : PowerSeries.constantCoeff η = 0) (hα0 : PowerSeries.constantCoeff α = 0)
    (hη : PowerSeries.subst η g = PowerSeries.subst f η)
    (hα : PowerSeries.subst α f = PowerSeries.subst g α) :
    PowerSeries.subst (PowerSeries.subst α η) g =
      PowerSeries.subst g (PowerSeries.subst α η) :=
  congrArg Subtype.val
    (comp_intertwine_left (M := SubstNilp A) substComp substComp_assoc
      ⟨f, hf0⟩ ⟨g, hg0⟩ ⟨η, hη0⟩ ⟨α, hα0⟩ (Subtype.ext hη.symm) (Subtype.ext hα))

/-- 純抽象核 2 の冪級数版: `θ ∘ η = X` ならば `η` は左から簡約できる。 -/
theorem subst_cancel_left {A : Type*} [CommRing A] {θ η α β : PowerSeries A}
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hη0 : PowerSeries.constantCoeff η = 0)
    (hα0 : PowerSeries.constantCoeff α = 0) (hβ0 : PowerSeries.constantCoeff β = 0)
    (hθη : PowerSeries.subst η θ = PowerSeries.X)
    (h : PowerSeries.subst α η = PowerSeries.subst β η) : α = β :=
  congrArg Subtype.val
    (comp_cancel_of_left_inv (M := SubstNilp A) substComp substComp_assoc
      substOne ⟨θ, hθ0⟩ ⟨η, hη0⟩ (Subtype.ext hθη) substOne_comp
      ⟨α, hα0⟩ ⟨β, hβ0⟩ (Subtype.ext h))

/-- 純抽象核 3 の冪級数版: `a : u → v` と `b : v → w` の合成 `b ∘ a` は `u → w`。 -/
theorem subst_intertwine_comp {A : Type*} [CommRing A] {u v w a b : PowerSeries A}
    (hu0 : PowerSeries.constantCoeff u = 0) (hv0 : PowerSeries.constantCoeff v = 0)
    (hw0 : PowerSeries.constantCoeff w = 0) (ha0 : PowerSeries.constantCoeff a = 0)
    (hb0 : PowerSeries.constantCoeff b = 0)
    (ha : PowerSeries.subst a v = PowerSeries.subst u a)
    (hb : PowerSeries.subst b w = PowerSeries.subst v b) :
    PowerSeries.subst (PowerSeries.subst a b) w =
      PowerSeries.subst u (PowerSeries.subst a b) :=
  congrArg Subtype.val
    (comp_intertwine_comp (M := SubstNilp A) substComp substComp_assoc
      ⟨u, hu0⟩ ⟨v, hv0⟩ ⟨w, hw0⟩ ⟨a, ha0⟩ ⟨b, hb0⟩ (Subtype.ext ha) (Subtype.ext hb))

/-! ## 3. 二面版一意性 —— 原典 Lemma 3.4(`t = 1`)の `f ≠ f′` 版

★木にあるのは単一 `h` 版 `powerSeries_uniqueness`
(`Found/PGC/LubinTateUniqueness.lean:174`、`h ∘ α = α ∘ h` の形)だけだった。
`η` による左簡約でここから二面版を出す。 -/

/-- ★★★★★★★★★**二面版の一意性補題**(原典 Lemma 3.4 の `t = 1`・`f ≠ f′` の場合)。

`f ∘ α = α ∘ g`・`f ∘ β = β ∘ g` を満たし、定数項が `0` で 1 次係数が一致する
2 つの冪級数 `α, β` は等しい。ただし逆向きの絡み作用素 `η`(`g ∘ η = η ∘ f`)と
その左逆元 `θ`(`θ ∘ η = X`)が与えられているとする。

証明: `η ∘ α` と `η ∘ β` はどちらも `g` と可換(`subst_comp_intertwine_left`)で
1 次係数が等しい(`coeff_one_subst_1var`)ので、単一 `h` 版の一意性
`powerSeries_uniqueness`(`h := g`)から `η ∘ α = η ∘ β`。あとは `θ` で左簡約する
(`subst_cancel_left`)。

★分岐・Galois の語彙は 1 つも出てこないが、`powerSeries_uniqueness` が
`π^{n+1} - π` で割るために `IsLocalRing`・`IsDomain`・`𝔪 = (π)`・`π ≠ 0` を要求する。
★`f` 側には 1 次係数の仮定が要らない(`hg1` だけでよい)。 -/
theorem powerSeries_uniqueness_pair {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    {f g θ η α β : PowerSeries A}
    (hf0 : PowerSeries.constantCoeff f = 0) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hη0 : PowerSeries.constantCoeff η = 0)
    (hθη : PowerSeries.subst η θ = PowerSeries.X)
    (hηint : PowerSeries.subst η g = PowerSeries.subst f η)
    (hα0 : PowerSeries.constantCoeff α = 0) (hβ0 : PowerSeries.constantCoeff β = 0)
    (hlead : PowerSeries.coeff 1 α = PowerSeries.coeff 1 β)
    (hαint : PowerSeries.subst α f = PowerSeries.subst g α)
    (hβint : PowerSeries.subst β f = PowerSeries.subst g β) :
    α = β := by
  refine subst_cancel_left hθ0 hη0 hα0 hβ0 hθη ?_
  refine powerSeries_uniqueness hπmax hπne0 hg0 hg1
    (PowerSeries.constantCoeff_subst_eq_zero hα0 η hη0)
    (PowerSeries.constantCoeff_subst_eq_zero hβ0 η hη0) ?_ ?_ ?_
  · rw [coeff_one_subst_1var hα0, coeff_one_subst_1var hβ0, hlead]
  · exact (subst_comp_intertwine_left hf0 hg0 hη0 hα0 hηint hαint).symm
  · exact (subst_comp_intertwine_left hf0 hg0 hη0 hβ0 hηint hβint).symm

/-- ★★★★★★★★**抽象層の `𝒪`-線型性**(原典 Proposition 3.5(iii) の骨)。

`θ : F_g → F_f` が絡み作用素(`f ∘ θ = θ ∘ g`)で、`φ` が `f` と可換・`ψ` が `g` と可換で
1 次係数が等しいなら、`θ ∘ ψ = φ ∘ θ`。★`φ = [a]_f`・`ψ = [a]_g` と読む。

両辺はどちらも `F_g → F_f` の絡み作用素(`subst_intertwine_comp` を 2 回)で、
1 次係数が `coeff 1 θ · coeff 1 ψ = coeff 1 φ · coeff 1 θ` と一致するので、
二面版一意性で等しい。★`coeff 1 θ = 1` は要らない。 -/
theorem subst_comm_of_intertwine_pair {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    {f g θ η φ ψ : PowerSeries A}
    (hf0 : PowerSeries.constantCoeff f = 0) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hη0 : PowerSeries.constantCoeff η = 0)
    (hθη : PowerSeries.subst η θ = PowerSeries.X)
    (hηint : PowerSeries.subst η g = PowerSeries.subst f η)
    (hθint : PowerSeries.subst θ f = PowerSeries.subst g θ)
    (hφ0 : PowerSeries.constantCoeff φ = 0) (hψ0 : PowerSeries.constantCoeff ψ = 0)
    (hφ : PowerSeries.subst φ f = PowerSeries.subst f φ)
    (hψ : PowerSeries.subst ψ g = PowerSeries.subst g ψ)
    (hlead : PowerSeries.coeff 1 φ = PowerSeries.coeff 1 ψ) :
    PowerSeries.subst ψ θ = PowerSeries.subst θ φ := by
  refine powerSeries_uniqueness_pair hπmax hπne0 hf0 hg0 hg1 hθ0 hη0 hθη hηint
    (PowerSeries.constantCoeff_subst_eq_zero hψ0 θ hθ0)
    (PowerSeries.constantCoeff_subst_eq_zero hθ0 φ hφ0) ?_ ?_ ?_
  · rw [coeff_one_subst_1var hψ0, coeff_one_subst_1var hθ0, hlead, mul_comm]
  · exact subst_intertwine_comp hg0 hg0 hf0 hψ0 hθ0 hψ hθint
  · exact subst_intertwine_comp hg0 hf0 hf0 hθ0 hφ0 hθint hφ

/-- ★★★★★★★★**抽象層の合成則**(原典 Proposition 3.5(ii) の第 2 式の骨)。

`a : F_u → F_v`・`b : F_v → F_w`・`c : F_u → F_w` がいずれも絡み作用素で、
`coeff 1 c = coeff 1 b · coeff 1 a` なら `b ∘ a = c`。 -/
theorem subst_comp_eq_of_intertwine_pair {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    {u v w θ η a b c : PowerSeries A}
    (hu0 : PowerSeries.constantCoeff u = 0) (hu1 : PowerSeries.coeff 1 u = π)
    (hv0 : PowerSeries.constantCoeff v = 0) (hw0 : PowerSeries.constantCoeff w = 0)
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hη0 : PowerSeries.constantCoeff η = 0)
    (hθη : PowerSeries.subst η θ = PowerSeries.X)
    (hηint : PowerSeries.subst η u = PowerSeries.subst w η)
    (ha0 : PowerSeries.constantCoeff a = 0) (hb0 : PowerSeries.constantCoeff b = 0)
    (hc0 : PowerSeries.constantCoeff c = 0)
    (haint : PowerSeries.subst a v = PowerSeries.subst u a)
    (hbint : PowerSeries.subst b w = PowerSeries.subst v b)
    (hcint : PowerSeries.subst c w = PowerSeries.subst u c)
    (hlead : PowerSeries.coeff 1 c = PowerSeries.coeff 1 b * PowerSeries.coeff 1 a) :
    PowerSeries.subst a b = c := by
  refine powerSeries_uniqueness_pair hπmax hπne0 hw0 hu0 hu1 hθ0 hη0 hθη hηint
    (PowerSeries.constantCoeff_subst_eq_zero ha0 b hb0) hc0 ?_ ?_ hcint
  · rw [coeff_one_subst_1var ha0, hlead]
  · exact subst_intertwine_comp hu0 hv0 hw0 ha0 hb0 haint hbint

/-! ## 4. 具体層 —— 木の Lubin-Tate 機械への代入

木の `LubinTateEndo hq hπmax u hu0 hu1 hu v hv0 hv1 hv a` は原典の `[a]_{u,v}`
(始域 `F_u`・終域 `F_v`)であり、関数等式は
`PowerSeries.subst (LubinTateEndo …) v = PowerSeries.subst u (LubinTateEndo …)`、
すなわち `v ∘ [a]_{u,v} = [a]_{u,v} ∘ u` である。 -/

/-- ★★★★★★★★★**原典 Proposition 3.5(ii) の合成則** `[θ′]_{f′,f″} ∘ [θ]_{f,f′} = [θθ′]_{f,f″}`。

3 つの Lubin-Tate 級数 `u, v, w`(同じ素元 `π`・同じ剰余還元 `X^q`)に対し
`[b]_{v,w} ∘ [a]_{u,v} = [ab]_{u,w}`。★代入の記法では
`PowerSeries.subst [a]_{u,v} [b]_{v,w}` が `[b] ∘ [a]` である。

証明は二面版一意性 `powerSeries_uniqueness_pair` に代入するだけ。逆向きの
絡み作用素 `η := [1]_{w,u}` とその左逆元 `θ := [1]_{u,w}` は木の
`subst_eq_X_of_intertwine_pair` が供給する。

★木の `LubinTateAction_comp`(`Found/PGC/LubinTateActionMul.lean:77`)は
`u = v = w` の場合である。本補題はそれを 3 つの異なる級数へ広げたもの。 -/
theorem subst_lubinTateEndo_lubinTateEndo {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (u : PowerSeries A) (hu0 : PowerSeries.constantCoeff u = 0)
    (hu1 : PowerSeries.coeff 1 u = π)
    (hu : PowerSeries.map (IsLocalRing.residue A) u = PowerSeries.X ^ (pp ^ ff))
    (v : PowerSeries A) (hv0 : PowerSeries.constantCoeff v = 0)
    (hv1 : PowerSeries.coeff 1 v = π)
    (hv : PowerSeries.map (IsLocalRing.residue A) v = PowerSeries.X ^ (pp ^ ff))
    (w : PowerSeries A) (hw0 : PowerSeries.constantCoeff w = 0)
    (hw1 : PowerSeries.coeff 1 w = π)
    (hw : PowerSeries.map (IsLocalRing.residue A) w = PowerSeries.X ^ (pp ^ ff))
    (a b : A) :
    PowerSeries.subst
        (LubinTateEndo hq hπmax u hu0 hu1 hu v
          ((PowerSeries.coeff_zero_eq_constantCoeff_apply v).trans hv0) hv1 hv a)
        (LubinTateEndo hq hπmax v hv0 hv1 hv w
          ((PowerSeries.coeff_zero_eq_constantCoeff_apply w).trans hw0) hw1 hw b) =
      LubinTateEndo hq hπmax u hu0 hu1 hu w
        ((PowerSeries.coeff_zero_eq_constantCoeff_apply w).trans hw0) hw1 hw (a * b) := by
  have hu0' : PowerSeries.coeff 0 u = 0 :=
    (PowerSeries.coeff_zero_eq_constantCoeff_apply u).trans hu0
  have hv0' : PowerSeries.coeff 0 v = 0 :=
    (PowerSeries.coeff_zero_eq_constantCoeff_apply v).trans hv0
  have hw0' : PowerSeries.coeff 0 w = 0 :=
    (PowerSeries.coeff_zero_eq_constantCoeff_apply w).trans hw0
  -- `θ := [1]_{u,w}` と `η := [1]_{w,u}`、および `θ ∘ η = X`
  have hθ0 := constantCoeff_LubinTateEndo hq hπmax u hu0 hu1 hu w hw0' hw1 hw 1
  have hθ1 := coeff_one_LubinTateEndo hq hπmax u hu0 hu1 hu w hw0' hw1 hw 1
  have hθint := LubinTateEndo_functional_equation hq hπmax u hu0 hu1 hu w hw0' hw1 hw 1
  have hη0 := constantCoeff_LubinTateEndo hq hπmax w hw0 hw1 hw u hu0' hu1 hu 1
  have hη1 := coeff_one_LubinTateEndo hq hπmax w hw0 hw1 hw u hu0' hu1 hu 1
  have hηint := LubinTateEndo_functional_equation hq hπmax w hw0 hw1 hw u hu0' hu1 hu 1
  have hθη := subst_eq_X_of_intertwine_pair hπmax hπne0 hw0 hw1 hu0 hθ0 hθ1 hη0 hη1 hθint hηint
  -- `a`・`b`・`ab` の 3 本
  have ha0 := constantCoeff_LubinTateEndo hq hπmax u hu0 hu1 hu v hv0' hv1 hv a
  have ha1 := coeff_one_LubinTateEndo hq hπmax u hu0 hu1 hu v hv0' hv1 hv a
  have haint := LubinTateEndo_functional_equation hq hπmax u hu0 hu1 hu v hv0' hv1 hv a
  have hb0 := constantCoeff_LubinTateEndo hq hπmax v hv0 hv1 hv w hw0' hw1 hw b
  have hb1 := coeff_one_LubinTateEndo hq hπmax v hv0 hv1 hv w hw0' hw1 hw b
  have hbint := LubinTateEndo_functional_equation hq hπmax v hv0 hv1 hv w hw0' hw1 hw b
  have hc0 := constantCoeff_LubinTateEndo hq hπmax u hu0 hu1 hu w hw0' hw1 hw (a * b)
  have hc1 := coeff_one_LubinTateEndo hq hπmax u hu0 hu1 hu w hw0' hw1 hw (a * b)
  have hcint := LubinTateEndo_functional_equation hq hπmax u hu0 hu1 hu w hw0' hw1 hw (a * b)
  refine powerSeries_uniqueness_pair hπmax hπne0 hw0 hu0 hu1 hθ0 hη0 hθη hηint
    (PowerSeries.constantCoeff_subst_eq_zero ha0 _ hb0) hc0 ?_ ?_ hcint
  · rw [coeff_one_subst_1var ha0, ha1, hb1, hc1, mul_comm]
  · exact subst_intertwine_comp hu0 hv0 hw0 ha0 hb0 haint hbint

/-- ★★★★★★★★**原典 Corollary 3.7(ii)** —— `θ` が単元なら `[θ]_{f,f′}` は同型で
逆は `[θ^{-1}]_{f′,f}`。

`a * b = 1` のとき `[a]_{u,v} ∘ [b]_{v,u} = X`。合成則 `subst_lubinTateEndo_lubinTateEndo`
と `[1]_u = X`(木の `LubinTateAction_one_eq_X`)から出る。★`u` と `v`・`a` と `b` を
入れ替えて 2 度当てれば両側の逆であることが出る。 -/
theorem subst_lubinTateEndo_lubinTateEndo_eq_X {A : Type*} [CommRing A] [IsLocalRing A]
    [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (u : PowerSeries A) (hu0 : PowerSeries.constantCoeff u = 0)
    (hu1 : PowerSeries.coeff 1 u = π)
    (hu : PowerSeries.map (IsLocalRing.residue A) u = PowerSeries.X ^ (pp ^ ff))
    (v : PowerSeries A) (hv0 : PowerSeries.constantCoeff v = 0)
    (hv1 : PowerSeries.coeff 1 v = π)
    (hv : PowerSeries.map (IsLocalRing.residue A) v = PowerSeries.X ^ (pp ^ ff))
    (a b : A) (hab : a * b = 1) :
    PowerSeries.subst
        (LubinTateEndo hq hπmax u hu0 hu1 hu v
          ((PowerSeries.coeff_zero_eq_constantCoeff_apply v).trans hv0) hv1 hv a)
        (LubinTateEndo hq hπmax v hv0 hv1 hv u
          ((PowerSeries.coeff_zero_eq_constantCoeff_apply u).trans hu0) hu1 hu b) =
      PowerSeries.X := by
  have hu0' : PowerSeries.coeff 0 u = 0 :=
    (PowerSeries.coeff_zero_eq_constantCoeff_apply u).trans hu0
  rw [subst_lubinTateEndo_lubinTateEndo hq hπmax hπne0 u hu0 hu1 hu v hv0 hv1 hv u hu0 hu1 hu a b,
    hab]
  exact LubinTateAction_one_eq_X hq hπmax hπne0 u hu0' hu1 hu

/-- ★★★★★★★★★**原典 Proposition 3.5(iii) の `𝒪`-線型性**:
`[θ]_{g,f} ∘ [a]_g = [a]_f ∘ [θ]_{g,f}`。

★これが持ち場の 1 本目 —— 消費側(`Found/PGC/LubinTateReciprocityIndependence.lean` の
Corollary 4.9 後半)が名指しした「`subst ([a]_g) θ = subst θ ([a]_f)`」そのものである。

★`[a]_g` と `[a]_f` は**別の形式群 `F_g`・`F_f` の自己準同型**であり、
`[θ]_{g,f}` はその 2 つの間の準同型である。どちらの側も
`LubinTateAction`(木の `Found/PGC/LubinTateActionEndomorphism.lean:55`)で書いてある。

証明は合成則を 2 回当てるだけ: 左辺 `= [aθ]_{g,f}`、右辺 `= [θa]_{g,f}`、
あとは `mul_comm`。 -/
theorem subst_lubinTateAction_lubinTateEndo {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (g : PowerSeries A) (hg0 : PowerSeries.coeff 0 g = 0) (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))
    (c a : A) :
    PowerSeries.subst
        (LubinTateAction hq hπmax g hg0 hg1 hg a)
        (LubinTateEndo hq hπmax g
          ((PowerSeries.coeff_zero_eq_constantCoeff_apply g).symm.trans hg0) hg1 hg
          f hf0 hf1 hf c) =
      PowerSeries.subst
        (LubinTateEndo hq hπmax g
          ((PowerSeries.coeff_zero_eq_constantCoeff_apply g).symm.trans hg0) hg1 hg
          f hf0 hf1 hf c)
        (LubinTateAction hq hπmax f hf0 hf1 hf a) := by
  have hg0c : PowerSeries.constantCoeff g = 0 :=
    (PowerSeries.coeff_zero_eq_constantCoeff_apply g).symm.trans hg0
  have hf0c : PowerSeries.constantCoeff f = 0 :=
    (PowerSeries.coeff_zero_eq_constantCoeff_apply f).symm.trans hf0
  have hLg : LubinTateAction hq hπmax g hg0 hg1 hg a =
      LubinTateEndo hq hπmax g hg0c hg1 hg g hg0 hg1 hg a := rfl
  have hLf : LubinTateAction hq hπmax f hf0 hf1 hf a =
      LubinTateEndo hq hπmax f hf0c hf1 hf f hf0 hf1 hf a := rfl
  rw [hLg, hLf,
    subst_lubinTateEndo_lubinTateEndo hq hπmax hπne0 g hg0c hg1 hg g hg0c hg1 hg f hf0c hf1 hf a c,
    subst_lubinTateEndo_lubinTateEndo hq hπmax hπne0 g hg0c hg1 hg f hf0c hf1 hf f hf0c hf1 hf c a,
    mul_comm]

end ABC3.Found.PGC
