import ABC3.Found.PGC.UniformizerCocycle

/-!
# `ρ_{f,m} : W(K̂^m_f/K) ≅ K^×/(1+𝔭^m)` —— Weil 群への延長(Yoshida 2008 Proposition 4.7)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Proposition 4.7(物理 p.8)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-4.html` の `#prop-4-7`。

原文 (Yoshida08 p.8):
> Proposition 4.7. Let m ≥ 1 and f ∈ O[scr]_L[X] as above, with the linear coefficient π. (i) The L^m_f is Galois over K, and the following map is bijective for any α ∈ µ^×_f,m: K^×/(1 + p[frak]^m) ∋ x mod 1 + p[frak]^m −→ [xπ_j]_f,f^(j)(α) ∈ _j∈Z[bb] µ^(j),×_f,m (v(x) = −j). (ii) Let L = K. The ρ_f,m of Proposition 4.4(iii) extend to isomorphisms: ρ_f,m : W(K^m_f/K) ∼ =−→ K^×/(1 + p[frak]^m). (ϕ^j on K, α → [xπ_j](α), ∀α ∈ µ_f,m) −→ x mod 1 + p[frak]^m (v(x) = −j) Setting K^LT_f := _m≥1 K^m_f, we get ρ_f : W(K^LT_f/K) ∼ =−→ K^× by passing to the limit.

★`pdftotext` は `∐`(非交和)と `⋃` を落とす。上の逐語で `_j∈Z[bb] µ^(j),×_f,m` と
見えている箇所の頭には `∐` があり、`_m≥1 K^m_f` の頭には `⋃` がある(原典 PDF 直読)。
★本ファイルは `∐` を **`Σ j : ℤ, B j`(依存和)** で書いた(下記「`∐` の書き方」)。

原典の証明(物理 p.8、`0_Source` の `.txt` の 450–484 行を直読した全文):

> Proof. (i): If v(x) = −j, then xπ_j ∈ Θ^{L,×}_{π,π^{(j)}} by Lemma 4.5, hence
> [xπ_j] : µ_{f,m} ≅→ µ^{(j)}_{f,m} by Lemma 4.6. As [xπ_j] is O-linear,
> v^{−1}(−j)/(1 + p^m) ∋ x ↦ [xπ_j](α) ∈ µ^{(j),×}_{f,m} is bijective for each j.
> As L(α) = L^m_f = L^m_{f^{(j)}} by Proposition 4.4(ii) and Lemma 4.6, the ϕ^j ∈ Aut(L/K)
> extends to L^m_f by α ↦ α′ for each α′ ∈ µ^{(j),×}_{f,m}, hence L^m_f is Galois over K.
> (ii): Let σ ∈ W(K̂^m_f/K) with σ|_{K̂} = ϕ^j. If α ∈ µ^×_{f,m}, then σ(α) ∈ µ^{(j),×}_{f,m},
> hence σ(α) = [xπ_j](α) for a unique x mod 1+p^m by (i). This holds for all α ∈ µ_{f,m}
> because σ([a]_f(α)) = [a]^{(j)}_f(σ(α)) = [a]_{f^{(j)}}[xπ_j](α) = [xπ_j][a]_f(α) for all a ∈ O
> (this shows the compatibility of ρ_{f,m} for varying m). The map ρ_{f,m} is a group
> homomorphism because if τ(α) = [yπ_{j′}](α), then στ(α) = σ([yπ_{j′}](α))
> = [yπ_{j′}]^{(j)}[xπ_j](α) = [y π^{(j)}_{j′} · xπ_j](α) = [xy · π_{j+j′}](α).
> It is bijective because it restricts to Gal(K̂^m_f/K̂) ≅ (O/p^m)^× = O^×/(1+p^m)
> by Proposition 4.4(iii) and the quotient W(K̂/K) = Frob^ℤ_K is mapped onto
> K^×/O^× ≅ ℤ, i.e. v ∘ ρ_{f,m} = v. □

記号は §4.2(物理 p.7)による(`π_j ∈ L^×`、`v_L(π_j) = j`、`( )^{(i)} := ( )^{ϕ^i}`)。
`π_j` そのものは Lemma 4.5(`Found/PGC/UniformizerCocycle.lean::uniformizerZ`)にある。

## 何を足したか

### 抽象核 A —— ℤ-次数による分解(原典の `∐_{j∈ℤ}`)

★分岐・付値・Galois・Lubin-Tate の語彙が 1 つも現れない。ただの `Equiv` の話である。

| 宣言 | 内容 |
|---|---|
| `negFiberSigmaEquiv` | ★★`w : X → ℤ` と各次数のファイバー同型 `e j` から `X ≃ Σ j : ℤ, B j`。写像は `x ↦ ⟨-w x, e _ ⟨x, _⟩⟩`(★符号は原典の `v(x) = −j`) |
| `negFiberSigmaEquiv_fst` / `_fst_eq` | ★`j` は `x` から `-w x` として決まる —— 非交和であることの中身 |
| `bijective_of_fiberwise_bijective` | ★原典の「the following map is bijective」の字面そのもの |
| `valFiberEquivKer` | 付値のファイバーは核の平行移動 `{y | w y = c} ≃ ker w` |
| `gradedKerSigmaEquiv` | ★A の組み立て —— `ker w ≃ B j`(Prop 4.4(iii))から `Q ≃ Σ j, B j` |

### 抽象核 B —— 「`×`(位数ちょうど `m`)」の正体は「生成元」

原典の `µ^×_{f,m}` は「`µ_{f,m}` を `𝒪/𝔭^m`-加群として生成する元」である。
自由階数 1 の加群では生成元の集合はちょうど `Rˣ` と全単射になる。

| 宣言 | 内容 |
|---|---|
| `unitsEquivGenerators` | ★★`Rˣ ≃ {β | r ↦ r • β が全射}`(`α` は自由階数 1 の基底) |
| `not_mem_range_of_not_generator` | ★退化 —— 生成元でない `β`(低い層の点)は像に入らない |

### 抽象核 C —— 短 5 補題(次数つき)

原典の「It is bijective because it restricts to … and the quotient … is mapped onto ℤ」。

| 宣言 | 内容 |
|---|---|
| `bijective_of_bijOn_ker` | `dH ∘ ρ = dG`、`dG` 全射、`ρ` が核の間で全単射 ⟹ `ρ` は全単射 |
| `bijective_of_bijOn_ker_neg` | ★符号つき版 `dH ∘ ρ = dG⁻¹`(原典の `v(x) = −j`) |
| `bijOn_of_subtype_bijective` | 部分群の間の全単射性から `Set.BijOn` へ |
| `valQuot` / `ker_valQuot` | `v` の商 `A ⧸ N` への降下(`N = 1+𝔭^m ⊆ 𝒪^× = ker v`) |

### 抽象核 D —— Weil 群を「対」として持つ

★★原典は `W(K̂^m_f/K)` の元を「`K̂` の上で `ϕ^j`、`µ_{f,m}` の上で `α ↦ [xπ_j](α)`」という
**対**として書いている。そこで Weil 群を新しく定義せず、
`Γ × Multiplicative ℤ` の部分群 `{(σ, j) | res σ = ϕ^j}` として持った。
群構造は積の部分群として自動で付く(`Subgroup` の構造を 1 つ書くだけで済む)。

| 宣言 | 内容 |
|---|---|
| `weilGroup` | ★★`W = {(σ, j) ∈ Γ × ℤ | res σ = ϕ^j}` |
| `weilDeg` | Frobenius 次数 `j`(第 2 成分)。★`ϕ` の位数無限を仮定しなくてよい |
| `weilDeg_surjective` | `ϕ ∈ res.range` なら `deg` は全射(原典の `W(K̂/K) = Frob^ℤ`) |
| `weilDegKerEquiv` | ★`ker(deg) ≅ ker(res)` —— 原典の `Gal(K̂^m_f/K̂)` |
| `weilReciprocity_bijective` | ★★★**Proposition 4.7(ii)** の中身 |
| `weilReciprocity_bijective_of_kerEquiv` | ★Prop 4.4(iii) を `res.ker ≃* v.ker` の形で受け取る版 |
| `weilReciprocityEquiv` | ★★★`W ≃* K^×/(1+𝔭^m)` |
| `weilReciprocity_toAdd_val` | ★符号 `v(ρ σ) = −deg σ`(退化検査 D2) |

### 具体層 —— Lemma 4.5 の `π_j` との接続

| 宣言 | 内容 |
|---|---|
| `fixedUnits` | `K^× = {x ∈ L^× | ϕ x = x}`(`ϕ` 不変な単元) |
| `zpow_apply_eq_self` | `ϕ x = x ⟹ (ϕ^j) x = x`(★`j : ℤ`。`Int.induction_on` の 3 場合) |
| `mul_uniformizerZ_mem_thetaSet` | ★★**証明 (i) の第 1 文** `xπ_j ∈ Θ^L_{π,π^{(j)}}` |
| `map_mul_uniformizerZ` | `v(xπ_j) = v(x)·v(π)^j` |
| `map_mul_uniformizerZ_eq_one` | ★★**`Θ^{L,×}` の `×`** —— `v(x) = v(π)^{−j}` なら `v(xπ_j) = 1` |
| `zpow_apply_mul_uniformizerZ_mul` | ★★★**証明 (ii) の準同型性の核** `(yπ_{j′})^{(j)}·(xπ_j) = xy·π_{j+j′}` |

### 組み立て

| 宣言 | 内容 |
|---|---|
| `prop_4_7_i_equiv` | ★★★**Proposition 4.7(i)** —— `Q ≃ Σ j : ℤ, µ^{(j),×}` |
| `prop_4_7_i_equiv_fst` | ★`j = −v(x)`(符号) |

## ★`∐` をどう書いたか

`Σ j : ℤ, B j`(依存和)で書いた。★**両方の読みを同時に満たしている**:

* 型としては依存和 = 非交和である(`⟨j, b⟩ = ⟨j′, b′⟩` は `j = j′` を含む)。
* 同時に、`negFiberSigmaEquiv_fst` が `(Θ x).1 = -w x` を与えるので、
  **`j` は `x` から `v(x) = −j` によって決まる**。原典の丸括弧の中身そのもの。

★構造化係が申告した「170dpi では `∐` と `⋃` が紛らわしい」は、
本ファイルの書き方では**どちらでも同じ**にはならない: `⋃` と読むと
`j` の一意性が落ちて `Σ` ではなく像の合併になり、全単射が主張できない。
原典の丸括弧 `(v(x) = −j)` が `j` を一意に決めているので `∐` が正しい。

## 逸脱の記録

1. **(i) の「`L^m_f` is Galois over K」を落とした。**
   原典 (i) は前半(Galois 性)と後半(全単射)の 2 つの主張を含むが、
   本ファイルは**後半だけ**を扱う。理由: Galois 性の証明は
   「`L(α) = L^m_f = L^m_{f^{(j)}}` by Proposition 4.4(ii) and Lemma 4.6」に依存し、
   Proposition 4.4 も Lemma 4.6 もまだ木に無い。
   ★Galois 性は新ノード(下記「新しく必要になったノード」)。
2. **`µ_{f,m}`・`[θ]` を作らず、それらが供給する型と全単射を仮定として受け取った。**
   `prop_4_7_i_equiv` は `M : ℤ → Type*`(= `µ^{(j)}_{f,m}`)と
   `u : ker w ≃* Rˣ`(= Proposition 4.4(iii))と
   `α j`(= 基底、原典の `α ∈ µ^×_{f,m}` とその像)を仮定に取る。
   ★これは弱めではなく一般化である。木に Lubin-Tate 加群 `µ_{f,m}` が
   `𝒪/𝔭^m`-加群として立った時点で代入すれば原典の主張になる。
3. **(ii) の `ρ_{f,m}` の構成((i) からの一意性)を落とし、`ρ` を仮定として受け取った。**
   原典 (ii) の前半(「σ(α) = [xπ_j](α) for a unique x」から `ρ` を作る)は
   `µ_{f,m}` 上の `W` の作用を要する。本ファイルが与えるのは
   **`ρ` が同型であること**(原典の証明の最後の 2 文)である。
   ★準同型性の計算((ii) の 3 文目)は具体層の
   `zpow_apply_mul_uniformizerZ_mul` として原典どおりの形で入れた。
4. **`ρ_f`(`m → ∞` の極限)は入れていない。** 新ノード。
5. `L^m_f` の付値を `Lˣ →* N`(`N` は任意の可換群)の形で受け取った。
   `N = Multiplicative ℤ` を代入すれば原典の `v_L` になる。
   ★Lemma 4.5 の `map_zpowProd` がこの形なので、そこに合わせた。
6. `weilGroup` は `ϕ` の位数が無限であることを仮定しない。
   ★位数有限なら `deg` が「見かけの次数」を覚えているだけの被覆になるが、
   原典の設定(`ϕ` = 算術 Frobenius、`Gal(K̂/K) = Ẑ`)では `ϕ` は無限位数なので
   `weilGroup res ϕ` は原典の `W(K̂^m_f/K)` と一致する。
   ★仮定を落としたのは、`deg` を作るのに `Classical.choose` が要らなくなるからである。

## ★退化の自己検査

* (D1) **`α ∈ µ^×_{f,m}`(`×` = 位数ちょうど `m`)を落とすと全単射でない。**
  抽象核 B がこれを正面から扱う: 標的は `{β | r ↦ r • β が全射}`(= 生成元)であり、
  `not_mem_range_of_not_generator` が「生成元でない `β` は像に入らない」を言う。
  低い層の点(`𝔭^{m−1}`-捩れ)は生成元でないので、`×` を落とすと単射性が壊れる。
* (D2) **`v(x) = −j` の符号を落とすと写像が逆向きになる。**
  `negFiberSigmaEquiv` の第 1 成分は `-w x`(`negFiberSigmaEquiv_fst`)。
  (ii) 側は `bijective_of_bijOn_ker_neg` が `v ∘ ρ = deg⁻¹` を要求し、
  `weilReciprocity_toAdd_val` が `toAdd (v (ρ σ)) = -(toAdd (deg σ))` を出す。
  ★符号を落とした `bijective_of_bijOn_ker`(符号なし)も別に置いてあるので、
  どちらを使ったかが型で見える。
* (D3) **`m ≥ 1` を落とすと空虚。** `1 + 𝔭^m ⊆ 𝒪^×` は `m ≥ 1` でしか成り立たない
  (`m = 0` なら `1 + 𝔭^0 = 𝒪` で、そもそも `K^×` の部分群でない)。
  本ファイルではこれが `hN : N ≤ v.ker`(`valQuot`)として現れる。
  ★`N ≤ ker v` を落とすと `v` が商に降りず `w : Q →* Multiplicative ℤ` が作れない。
* (D4) **(ii) の `L = K̂`。** 原典が明示している限定であり、落とすと `W(K̂^m_f/K)` が
  意味を持たない(`K̂` の上の Frobenius `ϕ` が無いと `weilGroup` の `ϕ` が取れない)。
  本ファイルでは `res : Γ →* Δ` と `ϕ : Δ` を仮定に取ることでこれを型に載せた
  (`Γ = Gal(K̂^m_f/K)`、`Δ = Gal(K̂/K)`、`ϕ` = 算術 Frobenius)。
  ★`ϕ ∈ res.range` を落とすと `deg` が全射でなくなり、`ρ` の全射性が壊れる。
* (D5) `j = 0` で `π_0 = 1`、`weilDeg = 1`、`ker(deg) = Gal(K̂^m_f/K̂)`。
  この層だけを見れば Proposition 4.4(iii) そのものに退化する
  (`weilDegKerEquiv` + `weilReciprocity_bijective_of_kerEquiv` の `e`)。
* (D6) `ℕ∞` の切り詰め引き算は使っていない。次数は `ℤ` のままである。
-/

namespace ABC3.Found.PGC

def prop_4_7_i_equiv.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 8, item := "Proposition 4.7", sectionId := "prop-4-7" }

def weilReciprocityEquiv.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 8, item := "Proposition 4.7", sectionId := "prop-4-7" }

/-! ## §1 抽象核 A —— ℤ-次数による分解(原典の `∐_{j∈ℤ}`)

★ここには分岐・付値・Galois・Lubin-Tate の語彙が 1 つも現れない。
`w : X → ℤ` と、各次数のファイバーの同型があるだけの `Equiv` の話である。 -/

section GradedDecomposition

variable {X : Type*} {B : ℤ → Type*}

/-- ★★**`X ≃ ∐_{j∈ℤ} B j`**。

`w : X → ℤ` が付値(原典の `v`)、`B j` が原典の `µ^{(j),×}_{f,m}`。
各次数のファイバー `{x | w x = -j}` が `B j` と同型なら、全体が非交和と同型になる。
★写像は `x ↦ ⟨-w x, e _ ⟨x, _⟩⟩` で、**第 1 成分の符号が原典の `v(x) = −j`** である。 -/
def negFiberSigmaEquiv (w : X → ℤ) (e : ∀ j : ℤ, {x : X // w x = -j} ≃ B j) :
    X ≃ Σ j : ℤ, B j where
  toFun x := ⟨-w x, e (-w x) ⟨x, by ring⟩⟩
  invFun p := ((e p.1).symm p.2 : X)
  left_inv x := by simp
  right_inv := by
    rintro ⟨j, b⟩
    have key : ∀ (j : ℤ) (y : {x : X // w x = -j}),
        (⟨-w (y : X), e (-w (y : X)) ⟨(y : X), by ring⟩⟩ : Σ j : ℤ, B j) = ⟨j, e j y⟩ := by
      rintro j ⟨x, hx⟩
      have hj : j = -w x := by omega
      subst hj
      rfl
    simpa using key j ((e j).symm b)

theorem negFiberSigmaEquiv_apply (w : X → ℤ) (e : ∀ j : ℤ, {x : X // w x = -j} ≃ B j) (x : X) :
    negFiberSigmaEquiv w e x = ⟨-w x, e (-w x) ⟨x, by ring⟩⟩ := rfl

/-- ★**次数は `x` から決まる** —— これが「非交和」であることの中身。 -/
theorem negFiberSigmaEquiv_fst (w : X → ℤ) (e : ∀ j : ℤ, {x : X // w x = -j} ≃ B j) (x : X) :
    (negFiberSigmaEquiv w e x).1 = -w x := rfl

/-- ★**原典の丸括弧 `(v(x) = −j)`** —— `w x = -j` なら次数成分はちょうど `j`。 -/
theorem negFiberSigmaEquiv_fst_eq (w : X → ℤ) (e : ∀ j : ℤ, {x : X // w x = -j} ≃ B j)
    {x : X} {j : ℤ} (hx : w x = -j) : (negFiberSigmaEquiv w e x).1 = j := by
  rw [negFiberSigmaEquiv_fst, hx, neg_neg]

theorem negFiberSigmaEquiv_symm_apply (w : X → ℤ) (e : ∀ j : ℤ, {x : X // w x = -j} ≃ B j)
    (p : Σ j : ℤ, B j) : (negFiberSigmaEquiv w e).symm p = ((e p.1).symm p.2 : X) := rfl

/-- ★**原典の「the following map is bijective」の字面そのもの**。

各次数で全単射なら、`x ↦ ⟨-w x, Φ (-w x) x⟩` は全単射。 -/
theorem bijective_of_fiberwise_bijective (w : X → ℤ)
    (Φ : ∀ j : ℤ, {x : X // w x = -j} → B j) (hΦ : ∀ j, Function.Bijective (Φ j)) :
    Function.Bijective (fun x : X => (⟨-w x, Φ (-w x) ⟨x, by ring⟩⟩ : Σ j : ℤ, B j)) :=
  (negFiberSigmaEquiv w (fun j => Equiv.ofBijective _ (hΦ j))).bijective

/-- 付値 `w` のファイバーは核の平行移動である。

原典の `v^{−1}(−j)/(1 + 𝔭^m) ≃ 𝒪^×/(1 + 𝔭^m)` に当たる部分。 -/
def valFiberEquivKer {Q A : Type*} [Group Q] [Group A] (w : Q →* A) {c : A}
    (t : Q) (ht : w t = c) : {y : Q // w y = c} ≃ w.ker where
  toFun y := ⟨t⁻¹ * y, by simp [MonoidHom.mem_ker, ht, y.2]⟩
  invFun k := ⟨t * k, by simpa [ht] using (MonoidHom.mem_ker.mp k.2)⟩
  left_inv y := by ext; simp
  right_inv k := by ext; simp

theorem valFiberEquivKer_apply {Q A : Type*} [Group Q] [Group A] (w : Q →* A) {c : A}
    (t : Q) (ht : w t = c) (y : {y : Q // w y = c}) :
    ((valFiberEquivKer w t ht y : w.ker) : Q) = t⁻¹ * y := rfl

/-- ★**抽象核 A の組み立て** —— 核の同型(原典では Proposition 4.4(iii))から
全体の非交和分解を作る。

`t j` は次数 `-j` の元(原典では `π_{-j}` の類。`v` の全射性から取れる)。 -/
def gradedKerSigmaEquiv {Q : Type*} [Group Q] (w : Q →* Multiplicative ℤ)
    (t : ℤ → Q) (ht : ∀ j : ℤ, w (t j) = Multiplicative.ofAdd (-j))
    (e : ∀ j : ℤ, w.ker ≃ B j) : Q ≃ Σ j : ℤ, B j :=
  negFiberSigmaEquiv (fun y => Multiplicative.toAdd (w y))
    (fun j => (valFiberEquivKer w (t j) (ht j)).trans (e j))

theorem gradedKerSigmaEquiv_fst {Q : Type*} [Group Q] (w : Q →* Multiplicative ℤ)
    (t : ℤ → Q) (ht : ∀ j : ℤ, w (t j) = Multiplicative.ofAdd (-j))
    (e : ∀ j : ℤ, w.ker ≃ B j) (y : Q) :
    (gradedKerSigmaEquiv w t ht e y).1 = -Multiplicative.toAdd (w y) := rfl

end GradedDecomposition

/-! ## §2 抽象核 B —— 「`×`(位数ちょうど `m`)」の正体は「生成元」

原典の `µ^×_{f,m}` は `µ_{f,m}` の `𝒪/𝔭^m`-加群としての生成元の集合である。
`µ_{f,m}` は自由階数 1 なので、生成元の集合は `(𝒪/𝔭^m)^×` と全単射になる。
★ここにも分岐・付値・Lubin-Tate は現れない。可換環上の自由階数 1 加群の話である。 -/

section FreeRankOne

variable {R M : Type*} [CommRing R] [AddCommGroup M] [Module R M]

/-- ★★**`Rˣ ≃ {生成元}`**。

`α` が `M` の基底(`r ↦ r • α` が全単射)のとき、`M` の生成元
(`r ↦ r • β` が全射になる `β`)はちょうど `u • α`(`u ∈ Rˣ`)である。

★原典の `µ^×_{f,m}`(位数がちょうど `𝔭^m` の捩れ点)がこの「生成元」であり、
Proposition 4.4(iii) の `(𝒪/𝔭^m)^×` がこの `Rˣ` である。 -/
noncomputable def unitsEquivGenerators (α : M) (hα : Function.Bijective (fun r : R => r • α)) :
    Rˣ ≃ {β : M // Function.Surjective (fun r : R => r • β)} where
  toFun u := ⟨(u : R) • α, by
    intro m
    obtain ⟨r, hr⟩ := hα.2 m
    refine ⟨r * (u⁻¹ : Rˣ), ?_⟩
    show (r * (u⁻¹ : Rˣ)) • ((u : R) • α) = m
    rw [← mul_smul, mul_assoc, Units.inv_mul, mul_one]
    exact hr⟩
  invFun β :=
    { val := (Equiv.ofBijective _ hα).symm β.1
      inv := Classical.choose (β.2 α)
      val_inv := hα.1 (by
        have hr : ((Equiv.ofBijective _ hα).symm β.1) • α = (β : M) :=
          (Equiv.ofBijective _ hα).apply_symm_apply β.1
        have hs : (Classical.choose (β.2 α)) • (β : M) = α := Classical.choose_spec (β.2 α)
        show (((Equiv.ofBijective _ hα).symm β.1) * Classical.choose (β.2 α)) • α = (1 : R) • α
        rw [mul_comm, mul_smul, hr, hs, one_smul])
      inv_val := hα.1 (by
        have hr : ((Equiv.ofBijective _ hα).symm β.1) • α = (β : M) :=
          (Equiv.ofBijective _ hα).apply_symm_apply β.1
        have hs : (Classical.choose (β.2 α)) • (β : M) = α := Classical.choose_spec (β.2 α)
        show ((Classical.choose (β.2 α)) * ((Equiv.ofBijective _ hα).symm β.1)) • α = (1 : R) • α
        rw [mul_smul, hr, hs, one_smul]) }
  left_inv u := by ext; exact (Equiv.ofBijective _ hα).symm_apply_apply (u : R)
  right_inv β := by ext; exact (Equiv.ofBijective _ hα).apply_symm_apply β.1

theorem unitsEquivGenerators_apply (α : M) (hα : Function.Bijective (fun r : R => r • α))
    (u : Rˣ) :
    ((unitsEquivGenerators α hα u : {β : M // Function.Surjective (fun r : R => r • β)}) : M)
      = (u : R) • α := rfl

/-- ★**退化検査 (D1)** —— 生成元でない `β`(低い層の捩れ点)は像に入らない。
`µ^×_{f,m}` の `×` を落とすと全単射でなくなる、の正確な形。 -/
theorem not_mem_range_of_not_generator (α : M) (hα : Function.Bijective (fun r : R => r • α))
    (β : M) (hβ : ¬ Function.Surjective (fun r : R => r • β)) :
    ∀ u : Rˣ, (u : R) • α ≠ β := fun u h => hβ (h ▸ (unitsEquivGenerators α hα u).2)

end FreeRankOne

/-! ## §3 抽象核 C —— 商への降下と短 5 補題

原典 (ii) の最後の 2 文
「It is bijective because it restricts to `Gal(K̂^m_f/K̂) ≅ (𝒪/𝔭^m)^×` …
and the quotient `W(K̂/K) = Frob^ℤ` is mapped onto `K^×/𝒪^× ≅ ℤ`」に当たる。 -/

section Descent

variable {A : Type*} [CommGroup A]

/-- 付値 `v` の商 `A ⧸ N` への降下。

★原典では `A = K^×`、`N = 1 + 𝔭^m`。**`m ≥ 1` が `N ≤ ker v` を保証する**
(退化検査 D3): `m = 0` なら `1 + 𝔭^0 = 𝒪` で `K^×` の部分群にすらならない。 -/
def valQuot (N : Subgroup A) (v : A →* Multiplicative ℤ) (hN : N ≤ v.ker) :
    A ⧸ N →* Multiplicative ℤ := QuotientGroup.lift N v hN

@[simp] theorem valQuot_mk (N : Subgroup A) (v : A →* Multiplicative ℤ) (hN : N ≤ v.ker) (a : A) :
    valQuot N v hN (QuotientGroup.mk a) = v a := rfl

theorem valQuot_surjective (N : Subgroup A) (v : A →* Multiplicative ℤ) (hN : N ≤ v.ker)
    (hv : Function.Surjective v) : Function.Surjective (valQuot N v hN) := fun c => by
  obtain ⟨a, ha⟩ := hv c
  exact ⟨QuotientGroup.mk a, ha⟩

/-- ★`ker(v̄) = 𝒪^×/(1 + 𝔭^m)` —— 商の核はもとの核の像。 -/
theorem ker_valQuot (N : Subgroup A) (v : A →* Multiplicative ℤ) (hN : N ≤ v.ker) :
    (valQuot N v hN).ker = v.ker.map (QuotientGroup.mk' N) := by
  ext x
  induction x using QuotientGroup.induction_on with
  | H a =>
    refine ⟨fun h => ⟨a, h, rfl⟩, ?_⟩
    rintro ⟨b, hb, hbe⟩
    have hba : (QuotientGroup.mk b : A ⧸ N) = QuotientGroup.mk a := hbe
    rw [MonoidHom.mem_ker, ← hba, valQuot_mk]
    exact hb

end Descent

section ShortFive

/-- ★**短 5 補題(次数つき)** —— 核の上で全単射、かつ次数を保つなら全単射。

原典: 「it restricts to `Gal(K̂^m_f/K̂) ≅ (𝒪/𝔭^m)^×` … and the quotient
`W(K̂/K) = Frob^ℤ_K` is mapped onto `K^×/𝒪^× ≅ ℤ`」。 -/
theorem bijective_of_bijOn_ker {G H A : Type*} [Group G] [Group H] [Group A]
    (ρ : G →* H) (dG : G →* A) (dH : H →* A) (hcomm : ∀ g, dH (ρ g) = dG g)
    (hd : Function.Surjective dG)
    (hker : Set.BijOn ρ (dG.ker : Set G) (dH.ker : Set H)) :
    Function.Bijective ρ := by
  constructor
  · rw [← MonoidHom.ker_eq_bot_iff, eq_bot_iff]
    intro g hg
    have hg1 : ρ g = 1 := hg
    have hgk : g ∈ (dG.ker : Set G) := by
      simp only [SetLike.mem_coe, MonoidHom.mem_ker, ← hcomm, hg1, map_one]
    have h1k : (1 : G) ∈ (dG.ker : Set G) := by simp
    simpa using hker.injOn hgk h1k (by simp [hg1])
  · intro h
    obtain ⟨g0, hg0⟩ := hd (dH h)
    have hmem : h * (ρ g0)⁻¹ ∈ (dH.ker : Set H) := by
      simp only [SetLike.mem_coe, MonoidHom.mem_ker, map_mul, map_inv, hcomm, hg0]
      simp
    obtain ⟨k, hk, hkeq⟩ := hker.surjOn hmem
    exact ⟨k * g0, by rw [map_mul, hkeq]; group⟩

/-- ★**符号つき版** —— 原典の `v(x) = −j`。

`dH ∘ ρ = dG⁻¹`(次数を**反転**して保つ)でも結論は同じ。
★退化検査 (D2): 符号なしの `bijective_of_bijOn_ker` と別宣言にしてあるので、
どちらを使ったかが型で見える。 -/
theorem bijective_of_bijOn_ker_neg {G H A : Type*} [Group G] [Group H] [CommGroup A]
    (ρ : G →* H) (dG : G →* A) (dH : H →* A) (hcomm : ∀ g, dH (ρ g) = (dG g)⁻¹)
    (hd : Function.Surjective dG)
    (hker : Set.BijOn ρ (dG.ker : Set G) (dH.ker : Set H)) :
    Function.Bijective ρ := by
  refine bijective_of_bijOn_ker ρ dG⁻¹ dH hcomm ?_ ?_
  · intro a
    obtain ⟨g, hg⟩ := hd a⁻¹
    exact ⟨g, by simp [show dG⁻¹ g = (dG g)⁻¹ from rfl, hg]⟩
  · have hk : (dG⁻¹).ker = dG.ker := by
      ext x
      simp only [MonoidHom.mem_ker, show dG⁻¹ x = (dG x)⁻¹ from rfl, inv_eq_one]
    rw [hk]
    exact hker

/-- 部分群の間の全単射性から `Set.BijOn` へ。 -/
theorem bijOn_of_subtype_bijective {G H : Type*} [Group G] [Group H] (ρ : G →* H)
    (S : Subgroup G) (T : Subgroup H) (hmap : ∀ g ∈ S, ρ g ∈ T)
    (hbij : Function.Bijective (fun g : S => (⟨ρ g, hmap g g.2⟩ : T))) :
    Set.BijOn ρ (S : Set G) (T : Set H) := by
  refine ⟨fun g hg => hmap g hg, ?_, ?_⟩
  · intro a ha b hb hab
    have h' : (fun g : S => (⟨ρ g, hmap g g.2⟩ : T)) ⟨a, ha⟩
        = (fun g : S => (⟨ρ g, hmap g g.2⟩ : T)) ⟨b, hb⟩ := Subtype.ext hab
    exact congrArg Subtype.val (hbij.1 h')
  · intro h hh
    obtain ⟨g, hg⟩ := hbij.2 ⟨h, hh⟩
    exact ⟨(g : G), g.2, congrArg Subtype.val hg⟩

end ShortFive

/-! ## §4 抽象核 D —— Weil 群を「対」として持つ

★★原典は `W(K̂^m_f/K)` の元を
「`K̂` の上で `ϕ^j`、`µ_{f,m}` の上で `α ↦ [xπ_j](α)`」という**対**として書いている。
そこで Weil 群を新しく定義せず、積 `Γ × Multiplicative ℤ` の部分群
`{(σ, j) | res σ = ϕ^j}` として持つ。
★`Γ = Gal(K̂^m_f/K)`、`Δ = Gal(K̂/K)`、`res` = 制限、`ϕ` = 算術 Frobenius。
★`ϕ` の位数が無限であることを仮定しないので、`deg` を作るのに選択公理が要らない。 -/

section Weil

variable {Γ Δ : Type*} [Group Γ] [Group Δ]

/-- ★★**Weil 群** `W = {(σ, j) ∈ Γ × ℤ | res σ = ϕ^j}`。

原典の `W(K̂^m_f/K)`(`Γ = Gal(K̂^m_f/K)` の、`Δ = Gal(K̂/K)` への制限が
`ϕ^ℤ` に落ちる元全体)。★群構造は積の部分群として自動で付く。 -/
def weilGroup (res : Γ →* Δ) (ϕ : Δ) : Subgroup (Γ × Multiplicative ℤ) where
  carrier := {p | res p.1 = ϕ ^ (Multiplicative.toAdd p.2)}
  one_mem' := by simp
  mul_mem' {a b} ha hb := by
    have ha' : res a.1 = ϕ ^ (Multiplicative.toAdd a.2) := ha
    have hb' : res b.1 = ϕ ^ (Multiplicative.toAdd b.2) := hb
    show res (a.1 * b.1) = ϕ ^ (Multiplicative.toAdd (a.2 * b.2))
    rw [map_mul, ha', hb', ← zpow_add]
    rfl
  inv_mem' {a} ha := by
    have ha' : res a.1 = ϕ ^ (Multiplicative.toAdd a.2) := ha
    show res a.1⁻¹ = ϕ ^ (Multiplicative.toAdd a.2⁻¹)
    rw [map_inv, ha', show Multiplicative.toAdd a.2⁻¹ = -(Multiplicative.toAdd a.2) from rfl,
      zpow_neg]

theorem mem_weilGroup {res : Γ →* Δ} {ϕ : Δ} {p : Γ × Multiplicative ℤ} :
    p ∈ weilGroup res ϕ ↔ res p.1 = ϕ ^ (Multiplicative.toAdd p.2) := Iff.rfl

/-- Frobenius 次数 `j`(対の第 2 成分)。原典の `σ|_{K̂} = ϕ^j` の `j`。 -/
def weilDeg (res : Γ →* Δ) (ϕ : Δ) : weilGroup res ϕ →* Multiplicative ℤ :=
  (MonoidHom.snd Γ (Multiplicative ℤ)).comp (weilGroup res ϕ).subtype

@[simp] theorem weilDeg_apply (res : Γ →* Δ) (ϕ : Δ) (σ : weilGroup res ϕ) :
    weilDeg res ϕ σ = (σ : Γ × Multiplicative ℤ).2 := rfl

/-- ★原典の `W(K̂/K) = Frob^ℤ_K`(次数写像は全射)。

★退化検査 (D4): `ϕ ∈ res.range` を落とすと全射でなくなり、`ρ` の全射性が壊れる。 -/
theorem weilDeg_surjective {res : Γ →* Δ} {ϕ : Δ} (hϕ : ϕ ∈ res.range) :
    Function.Surjective (weilDeg res ϕ) := by
  obtain ⟨γ, hγ⟩ := hϕ
  intro j
  refine ⟨⟨(γ ^ (Multiplicative.toAdd j), j), ?_⟩, rfl⟩
  show res (γ ^ (Multiplicative.toAdd j)) = ϕ ^ (Multiplicative.toAdd j)
  rw [map_zpow, hγ]

/-- ★**`ker(deg) ≅ ker(res)`** —— 原典の `Gal(K̂^m_f/K̂)`。

次数 `0` の層は、`K̂` の上で恒等な元、すなわち `Gal(K̂^m_f/K̂)` である。 -/
def weilDegKerEquiv (res : Γ →* Δ) (ϕ : Δ) : (weilDeg res ϕ).ker ≃* res.ker where
  toFun x := ⟨(x : weilGroup res ϕ).1.1, by
    have h1 : ((x : weilGroup res ϕ) : Γ × Multiplicative ℤ).2 = 1 := x.2
    have h2 : res ((x : weilGroup res ϕ) : Γ × Multiplicative ℤ).1
        = ϕ ^ (Multiplicative.toAdd ((x : weilGroup res ϕ) : Γ × Multiplicative ℤ).2) :=
      (x : weilGroup res ϕ).2
    rw [MonoidHom.mem_ker, h2, h1]
    simp⟩
  invFun k := ⟨⟨((k : Γ), 1), by
    show res (k : Γ) = ϕ ^ (Multiplicative.toAdd (1 : Multiplicative ℤ))
    rw [MonoidHom.mem_ker.mp k.2]
    simp⟩, rfl⟩
  left_inv x := by
    ext
    · rfl
    · exact (x.2 : ((x : weilGroup res ϕ) : Γ × Multiplicative ℤ).2 = 1).symm
  right_inv k := rfl
  map_mul' x y := rfl

/-- ★★★**Proposition 4.7(ii) の中身** —— `ρ_{f,m}` が全単射であること。

原典の最後の 2 文をそのまま写した:
`ρ` は核 `Gal(K̂^m_f/K̂)` の上で `𝒪^×/(1+𝔭^m)` への全単射に**制限され**(`hker`)、
次数 `Frob^ℤ` は `K^×/𝒪^× ≅ ℤ` の上へ**落ちる**(`hsign` + `weilDeg_surjective`)。
★`hsign` の `⁻¹` が原典の `v(x) = −j` である。 -/
theorem weilReciprocity_bijective {G : Type*} [Group G] {res : Γ →* Δ} {ϕ : Δ}
    (hϕ : ϕ ∈ res.range) (ρ : weilGroup res ϕ →* G) (v : G →* Multiplicative ℤ)
    (hsign : ∀ σ, v (ρ σ) = (weilDeg res ϕ σ)⁻¹)
    (hker : Set.BijOn ρ ((weilDeg res ϕ).ker : Set (weilGroup res ϕ)) (v.ker : Set G)) :
    Function.Bijective ρ :=
  bijective_of_bijOn_ker_neg ρ (weilDeg res ϕ) v hsign (weilDeg_surjective hϕ) hker

/-- ★**Proposition 4.4(iii) を `res.ker ≃* v.ker` の形で受け取る版**。

原典が「by Proposition 4.4(iii)」と言っている入力を、そのままの形
(`Gal(K̂^m_f/K̂) ≅ 𝒪^×/(1+𝔭^m)`)で受け取る。 -/
theorem weilReciprocity_bijective_of_kerEquiv {G : Type*} [Group G] {res : Γ →* Δ} {ϕ : Δ}
    (hϕ : ϕ ∈ res.range) (ρ : weilGroup res ϕ →* G) (v : G →* Multiplicative ℤ)
    (hsign : ∀ σ, v (ρ σ) = (weilDeg res ϕ σ)⁻¹)
    (e : res.ker ≃* v.ker)
    (hcompat : ∀ σ : (weilDeg res ϕ).ker,
      ρ (σ : weilGroup res ϕ) = (e (weilDegKerEquiv res ϕ σ) : G)) :
    Function.Bijective ρ := by
  have hmap : ∀ σ ∈ (weilDeg res ϕ).ker, ρ σ ∈ v.ker := by
    intro σ hσ
    rw [hcompat ⟨σ, hσ⟩]
    exact (e (weilDegKerEquiv res ϕ ⟨σ, hσ⟩)).2
  refine weilReciprocity_bijective hϕ ρ v hsign (bijOn_of_subtype_bijective ρ _ _ hmap ?_)
  have hfun : (fun σ : (weilDeg res ϕ).ker => (⟨ρ σ, hmap σ σ.2⟩ : v.ker))
      = fun σ => e (weilDegKerEquiv res ϕ σ) := funext fun σ => Subtype.ext (hcompat σ)
  rw [hfun]
  exact ((weilDegKerEquiv res ϕ).toEquiv.trans e.toEquiv).bijective

/-- ★★★**Proposition 4.7(ii)** `ρ_{f,m} : W(K̂^m_f/K) ≅ K^×/(1 + 𝔭^m)`。 -/
noncomputable def weilReciprocityEquiv {G : Type*} [Group G] (res : Γ →* Δ) (ϕ : Δ)
    (hϕ : ϕ ∈ res.range) (ρ : weilGroup res ϕ →* G) (v : G →* Multiplicative ℤ)
    (hsign : ∀ σ, v (ρ σ) = (weilDeg res ϕ σ)⁻¹)
    (hker : Set.BijOn ρ ((weilDeg res ϕ).ker : Set (weilGroup res ϕ)) (v.ker : Set G)) :
    weilGroup res ϕ ≃* G :=
  MulEquiv.ofBijective ρ (weilReciprocity_bijective hϕ ρ v hsign hker)

@[simp] theorem weilReciprocityEquiv_apply {G : Type*} [Group G] (res : Γ →* Δ) (ϕ : Δ)
    (hϕ : ϕ ∈ res.range) (ρ : weilGroup res ϕ →* G) (v : G →* Multiplicative ℤ)
    (hsign : ∀ σ, v (ρ σ) = (weilDeg res ϕ σ)⁻¹)
    (hker : Set.BijOn ρ ((weilDeg res ϕ).ker : Set (weilGroup res ϕ)) (v.ker : Set G))
    (σ : weilGroup res ϕ) : weilReciprocityEquiv res ϕ hϕ ρ v hsign hker σ = ρ σ := rfl

/-- ★**退化検査 (D2) の符号** —— `v(ρ σ) = −deg σ`(原典の `v(x) = −j`)。 -/
theorem weilReciprocity_toAdd_val {G : Type*} [Group G] {res : Γ →* Δ} {ϕ : Δ}
    (ρ : weilGroup res ϕ →* G) (v : G →* Multiplicative ℤ)
    (hsign : ∀ σ, v (ρ σ) = (weilDeg res ϕ σ)⁻¹) (σ : weilGroup res ϕ) :
    Multiplicative.toAdd (v (ρ σ)) = -(Multiplicative.toAdd (weilDeg res ϕ σ)) := by
  rw [hsign]
  rfl

end Weil

/-! ## §5 (i) の組み立て —— `K^×/(1+𝔭^m) ≃ ∐_{j∈ℤ} µ^{(j),×}_{f,m}` -/

/-- ★★★**Proposition 4.7(i)**(全単射の部分)。

原文 (Yoshida08 p.8):
> the following map is bijective for any α ∈ µ^×_f,m:
> K^×/(1 + p[frak]^m) ∋ x mod 1 + p[frak]^m −→ [xπ_j]_f,f^(j)(α) ∈ _j∈Z[bb] µ^(j),×_f,m (v(x) = −j).

* `Q` = `K^×/(1 + 𝔭^m)`、`w` = そこへ降りた付値(`valQuot`)。
* `t j` = 次数 `-j` の元(原典では `π_{-j}` の類)。`w` の全射性から取れる。
* `u : w.ker ≃* Rˣ` = **Proposition 4.4(iii)**(`𝒪^×/(1+𝔭^m) ≅ (𝒪/𝔭^m)^×`)。
* `M j` = `µ^{(j)}_{f,m}`、`α j` = その基底(原典の `α ∈ µ^×_{f,m}` とその `[π_j]` 像)。
* 標的 `{β | r ↦ r • β が全射}` = `µ^{(j),×}_{f,m}`(★`×` = 生成元、退化検査 D1)。

★`∐` は `Σ j : ℤ`(依存和)で書いた。次数成分が `-w y` であることは
`prop_4_7_i_equiv_fst` が言う —— これが原典の丸括弧 `(v(x) = −j)` である。

★逸脱 1: 原典 (i) の前半「`L^m_f` は `K` 上 Galois」は落とした(Prop 4.4(ii) と
Lemma 4.6 が木に無い)。
★逸脱 2: `µ_{f,m}` と `[θ]` は作らず、それらが供給する型と全単射を仮定に取った。 -/
noncomputable def prop_4_7_i_equiv {Q : Type*} [Group Q] {R : Type*} [CommRing R]
    {M : ℤ → Type*} [∀ j, AddCommGroup (M j)] [∀ j, Module R (M j)]
    (w : Q →* Multiplicative ℤ)
    (t : ℤ → Q) (ht : ∀ j : ℤ, w (t j) = Multiplicative.ofAdd (-j))
    (u : w.ker ≃* Rˣ)
    (α : ∀ j : ℤ, M j) (hα : ∀ j, Function.Bijective (fun r : R => r • α j)) :
    Q ≃ Σ j : ℤ, {β : M j // Function.Surjective (fun r : R => r • β)} :=
  gradedKerSigmaEquiv w t ht (fun j => (u.toEquiv).trans (unitsEquivGenerators (α j) (hα j)))

/-- ★**原典の丸括弧 `(v(x) = −j)`** —— 次数成分は `−v(x)`。 -/
theorem prop_4_7_i_equiv_fst {Q : Type*} [Group Q] {R : Type*} [CommRing R]
    {M : ℤ → Type*} [∀ j, AddCommGroup (M j)] [∀ j, Module R (M j)]
    (w : Q →* Multiplicative ℤ)
    (t : ℤ → Q) (ht : ∀ j : ℤ, w (t j) = Multiplicative.ofAdd (-j))
    (u : w.ker ≃* Rˣ)
    (α : ∀ j : ℤ, M j) (hα : ∀ j, Function.Bijective (fun r : R => r • α j)) (y : Q) :
    (prop_4_7_i_equiv w t ht u α hα y).1 = -Multiplicative.toAdd (w y) := rfl

/-! ## §6 具体層 —— Lemma 4.5 の `π_j` との接続

★原典の証明が `π_j` について実際に使っている 3 つの事実を、
`Found/PGC/UniformizerCocycle.lean` の `uniformizerZ` / `ThetaSet` の上で出す。
`ϕ : RingAut L` は任意の環自己同型でよい(Lemma 4.5 の逸脱 5 をそのまま引き継ぐ)。 -/

section Concrete

variable {L : Type*} [Field L]

/-- `K^× = {x ∈ L^× | ϕ x = x}` —— `ϕ` 不変な単元のなす部分群。

原典では `L = K̂` の Frobenius `ϕ` の固定体が `K` であり、`x ∈ K^×` がこれに当たる。 -/
def fixedUnits (ϕ : RingAut L) : Subgroup Lˣ where
  carrier := {x : Lˣ | ϕ (x : L) = (x : L)}
  one_mem' := by simp
  mul_mem' {a b} ha hb := by
    have ha' : ϕ (a : L) = (a : L) := ha
    have hb' : ϕ (b : L) = (b : L) := hb
    show ϕ ((a * b : Lˣ) : L) = ((a * b : Lˣ) : L)
    rw [Units.val_mul, map_mul, ha', hb']
  inv_mem' {a} ha := by
    have ha' : ϕ (a : L) = (a : L) := ha
    show ϕ ((a⁻¹ : Lˣ) : L) = ((a⁻¹ : Lˣ) : L)
    rw [Units.val_inv_eq_inv_val, map_inv₀, ha']

theorem mem_fixedUnits {ϕ : RingAut L} {x : Lˣ} : x ∈ fixedUnits ϕ ↔ ϕ (x : L) = (x : L) := Iff.rfl

/-- `ϕ` で固定される元は `ϕ^j`(`j : ℤ`)でも固定される。

★`Int.induction_on` の 3 場合(`zero`/`succ`/`pred`)。負の `j` を落としていない。 -/
theorem zpow_apply_eq_self (ϕ : RingAut L) {x : L} (hx : ϕ x = x) (j : ℤ) :
    (ϕ ^ j : RingAut L) x = x := by
  have hinv : ϕ⁻¹ x = x := by
    conv_lhs => rw [← hx]
    exact ϕ.symm_apply_apply x
  induction j using Int.induction_on with
  | zero => rfl
  | succ i ih =>
      rw [zpow_add_one, show ((ϕ ^ (i : ℤ) * ϕ : RingAut L)) x = (ϕ ^ (i : ℤ)) (ϕ x) from rfl,
        hx, ih]
  | pred i ih =>
      rw [zpow_sub_one,
        show ((ϕ ^ (-i : ℤ) * ϕ⁻¹ : RingAut L)) x = (ϕ ^ (-i : ℤ)) (ϕ⁻¹ x) from rfl, hinv, ih]

/-- ★★**原典の証明 (i) の第 1 文** `xπ_j ∈ Θ^L_{π,π^{(j)}}`。

> If v(x) = −j, then xπ_j ∈ Θ^{L,×}_{π,π^{(j)}} by Lemma 4.5

`Θ` は `ϕ`-不変な元での掛け算で閉じているので、Lemma 4.5 の第 2 主張
(`coe_uniformizerZ_mem_thetaSet`)に `x ∈ K` を掛けるだけで出る。
★`×`(単元であること)は付値の話であり `map_mul_uniformizerZ_eq_one` が担当する。 -/
theorem mul_uniformizerZ_mem_thetaSet (ϕ : RingAut L) (π : Lˣ) {x : L} (hx : ϕ x = x) (j : ℤ) :
    x * ((uniformizerZ ϕ π j : Lˣ) : L)
      ∈ ThetaSet ϕ (π : L) ((ϕ ^ j : RingAut L) (π : L)) := by
  have h := coe_uniformizerZ_mem_thetaSet ϕ π j
  rw [mem_thetaSet_iff] at h ⊢
  rw [map_mul, hx, mul_assoc, h, ← mul_assoc]

/-- `v(xπ_j) = v(x)·v(π)^j` —— 原典 §4.2 の `v_L(π_j) = j` の帰結。

`v` は `ϕ` 不変な群準同型(`N = Multiplicative ℤ` を代入すれば原典の `v_L`)。 -/
theorem map_mul_uniformizerZ {N : Type*} [CommGroup N] (ϕ : RingAut L) (v : Lˣ →* N)
    (hv : ∀ u : Lˣ, v ((unitsRingAutHom L ϕ) u) = v u) (π : Lˣ) (x : Lˣ) (j : ℤ) :
    v (x * uniformizerZ ϕ π j) = v x * (v π) ^ j := by
  rw [map_mul, uniformizerZ, map_zpowProd v hv]

/-- ★★**`Θ^{L,×}` の `×`** —— `v(x) = v(π)^{−j}`(原典の `v(x) = −j`)なら
`xπ_j` は付値 `0`、すなわち単元である。

★退化検査 (D2): 指数の `-j` を `j` にすると結論が `v x ^ 2` 型になって崩れる。 -/
theorem map_mul_uniformizerZ_eq_one {N : Type*} [CommGroup N] (ϕ : RingAut L) (v : Lˣ →* N)
    (hv : ∀ u : Lˣ, v ((unitsRingAutHom L ϕ) u) = v u) (π : Lˣ) (x : Lˣ) (j : ℤ)
    (hx : v x = (v π) ^ (-j)) : v (x * uniformizerZ ϕ π j) = 1 := by
  rw [map_mul_uniformizerZ ϕ v hv π x j, hx, zpow_neg, inv_mul_cancel]

/-- ★★★**原典の証明 (ii) の準同型性の核**。

> στ(α) = σ([yπ_{j′}](α)) = [yπ_{j′}]^{(j)}[xπ_j](α) = [y π^{(j)}_{j′} · xπ_j](α)
> = [xy · π_{j+j′}](α)

すなわち `(yπ_{j′})^{(j)} · (xπ_j) = xy · π_{j+j′}`。
使うのは (a) `y ∈ K`(`ϕ y = y`)と (b) §4.2 のコサイクル関係
`π_{j+j′} = π^{(j)}_{j′} π_j`(`coe_uniformizerZ_add`)だけである。

★`x` については `ϕ x = x` が**要らない**(`x` は `ϕ^j` の外側にある)。
原典は両方を `K^×` に取っているが、この等式に効くのは `y` の側だけである。 -/
theorem zpow_apply_mul_uniformizerZ_mul (ϕ : RingAut L) (π : Lˣ) (x : L) {y : L}
    (hy : ϕ y = y) (j j' : ℤ) :
    (ϕ ^ j : RingAut L) (y * ((uniformizerZ ϕ π j' : Lˣ) : L))
        * (x * ((uniformizerZ ϕ π j : Lˣ) : L))
      = (x * y) * ((uniformizerZ ϕ π (j + j') : Lˣ) : L) := by
  rw [map_mul, zpow_apply_eq_self ϕ hy, coe_uniformizerZ_add]
  ring

end Concrete

end ABC3.Found.PGC
