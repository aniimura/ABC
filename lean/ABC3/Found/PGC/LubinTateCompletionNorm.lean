import ABC3.Found.PGC.LubinTateCompletionGalois

/-!
# `N_{K̂^m_f/K̂^{ur}}(−α)` と「`α` は `K̂^m_f` の素元」—— Yoshida 2008 Proposition 4.4(ii) の残り

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Proposition 4.4(物理 p.7)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-4.html` の `#prop-4-4`。

原文 (Yoshida08 p.7):
> Proposition 4.4. Let m ≥ 1 and f ∈ O[scr]_L[X] as above, with the linear coefficient π. (i) The set µ_f,m is an O[scr]-module by +_F_f and [·]_f. For any α ∈ µ^×_f,m := µ_f,m \ µ_f,m−1, the following is an isomorphism of O[scr]-modules: O[scr]/p[frak]^m ∋ a mod p[frak]^m −→ [a]_f(α) ∈ µ_f,m. (ii) If α ∈ µ^×_f,m, then L^m_f = L(α), N_L^m_f/L(−α) = π^ϕ^m−1 and α is a uniformizer of L^m_f. The L^m_f/L is totally ramified Galois extension of degree |µ^×_f,m| = q^m−1(q − 1). (iii) We have canonical isomorphisms of abelian groups: ρ_f,m : Gal(L^m_f/L) ∼ =−→ Aut_O[scr](µ_f,m) ∼ =−→ (O[scr]/p[frak]^m)^×. (α → [u]_f(α), ∀α ∈ µ_f,m) −→ u mod p[frak]^m

★★本ノードは Proposition 4.4(ii) の **`N_{L^m_f/L}(−α)` の部分と「`α` は `L^m_f` の素元」の
部分**を、`L = K̂^{ur}` の場合について扱う。(ii) の残り(`L^m_f = L(α)`・次数)は
`LubinTateCompletionDegree.lean`(§4-a)が、(iii) は `LubinTateCompletionGalois.lean`(§4-b)が
既に持っている。

## ★★★原典との差 —— `N(−α) = π^{ϕ^{m−1}}` は木の設定では字面のままでは**偽**である

★**原典は `f ∈ 𝒪_L[X]`(多項式)に限っている**(逐語 "Let m ≥ 1 and f ∈ O_L[X] as above"、
および Lemma 4.3(i) の証明中の "as f_m is a monic in O_L[X]")。そのとき
`f_m/f_{m−1}` は多項式の商で、`f_m = f^{ϕ^{m−1}} ∘ f_{m−1}` と `f_{m−1}(0) = 0` から
定数項はちょうど `π^{ϕ^{m−1}}` になる。

★★**ところが木の Lubin-Tate 機械は `f` を冪級数として扱う**
(`f : PowerSeries 𝒪_K`、`coeff 1 f = π`、`map residue f = X^q`)。そこでの `ψ_n` は
`r_n := [π^n]_f/[π^{n−1}]_f` の **Weierstrass の distinguished 部分**であって、
商そのものではない。`r_n = ψ_n · U'_n`(`U'_n` は単元冪級数)の定数項を比べると

  `ψ_n.coeff 0 · U'_n(0) = π`  (`iteratedLubinTatePsi_coeff_zero_mul`、木に在る)

であり、**`U'_n(0) = 1` とは限らない**。したがって `ψ_n.coeff 0 = π` は一般には成り立たない。

★★**反例**(手で確かめたもの、Lean には書いていない): `K = ℚ_2`、`𝒪 = ℤ_2`、`q = 2`、`π = 2`、
`f = 2X + X² + 2X³`。これは木の仮定 `coeff 0 f = 0`・`coeff 1 f = π`・
`map residue f = X²` をすべて満たす。`n = 1` では `r_1 = f/X = 2 + X + 2X²` で、
`𝔪` の中の唯一の根は `α ≡ 22 (mod 32)`。`ψ_1 = X − α`、`ψ_1(0) = −α ≡ 10 (mod 32)`。
`[K̂^1_f : K̂^{ur}] = q − 1 = 1` なので `N(−α) = −α ≡ 10 ≠ 2 = π (mod 32)`。
★`v(−α) = 1` は保たれる(素元ではある)が、**等号は壊れる**。

★★したがって本ファイルは次の 3 段で述べる:

1. `N_{K̂^m_f/K̂^{ur}}(−β) = ψ_n.coeff 0`(**無仮定・厳密**、`norm_neg_eq_coeff_zero_psiCompletion`)。
2. `ψ_n.coeff 0` は `π` と**単元倍だけ違う**(`exists_unit_norm_neg` /
   `span_coeff_zero_iteratedLubinTatePsi` / `norm_norm_neg`)。すなわち
   **`N(−β)` は `𝒪_K` の素元**であり、`(N(−β)) = (π)`。
3. `U'_n(0) = 1` を仮定すれば原典どおり `N(−β) = π`(`norm_neg_eq_uniformizer`)。
   ★この仮定は原典の設定(`f` が多項式)では成り立つはずである(そのとき `r_n` 自身が
   distinguished になる)が、★★**それは本ファイルでは証明していない**
   ——「`f` が多項式なら `U'_n(0) = 1`」は新しい節点である。
   ★**木の一般の `f` に対して `U'_n(0) = 1` を証明することはできない**(上の反例)。

★**`π^{ϕ^{m−1}}` の Frobenius の冪について**: 本ファイルは `L = K̂^{ur}` に固定しており、
`π ∈ 𝒪_K` は `ϕ` で固定されるので `π^{ϕ^{m−1}} = π` である。本体の見当どおり
**指数は消える**(原典が `ϕ` を書くのは `f` の係数が一般の `𝒪_L` にあり、`f_m` が
`f^{ϕ^{m−1}} ∘ … ∘ f` と捻れているからである)。

## 何を示したか

### 抽象核(§0。分岐・付値・Lubin-Tate の語彙が 1 つも出てこない)

| 宣言 | 内容 |
|---|---|
| `norm_neg_powerBasis_gen` | ★**`N(−pb.gen) = minpoly の定数項`**。`Algebra.norm(−1) = (−1)^{dim}` と mathlib の `PowerBasis.norm_gen_eq_coeff_zero_minpoly`(`= (−1)^{dim}·coeff 0`)の符号がちょうど打ち消し合う。**これが原典が `−α` と書く理由である** |
| `norm_neg_adjoinSimple_gen` | 同じことを `F⟮α⟯` の生成元について |
| `norm_neg_eq_coeff_zero_minpoly` | `E = F⟮α⟯` で `β ∈ E` の像が `α` なら `N_{E/F}(−β) = (minpoly_F α).coeff 0` |

### 具体層

| 宣言 | 内容 |
|---|---|
| `monic_psiCompletion` / `coeff_psiCompletion` | `ψ_n` を `K̂^{ur}` 係数で見たものの基本性質 |
| `isIntegral_closureCompletionCoe_torsion` | `ι α` は `K̂^{ur}` 上整 |
| `exists_unit_coeff_zero_iteratedLubinTatePsi` | ★`ψ_n.coeff 0 = π · u`(`u ∈ 𝒪_K^×`) |
| `span_coeff_zero_iteratedLubinTatePsi` | `(ψ_n.coeff 0) = (π)` as ideals |
| ★★`norm_neg_eq_coeff_zero_psiCompletion` | ★★**`N_{K̂^m_f/K̂^{ur}}(−β) = ψ_n.coeff 0`** |
| `norm_neg_mul_constantCoeff_stepU` | `N(−β) · U'_n(0) = π` |
| ★`exists_unit_norm_neg` | ★**`N(−β) = π · u`** |
| `norm_neg_eq_uniformizer` | `U'_n(0) = 1` のとき `N(−β) = π`(原典の字面) |
| `norm_norm_neg` | `‖N(−β)‖ = ‖π‖`——`N(−β)` は `𝒪_{K̂^{ur}}` の素元 |
| `aeval_carrier_iteratedLubinTatePsi` / `norm_pow_torsion_point` | `‖α‖^{q^n−q^{n−1}} = ‖π‖` |
| ★★`norm_pow_finrank_closureCompletionCoe` | ★★**`‖α‖^{[K̂^m_f : K̂^{ur}]} = ‖π‖`**——「`α` は `K̂^m_f` の素元」 |
| ★★★`exists_uniformizer_norm_neg` | ★★★**4 つをまとめた形**(原始点の選択を隠す) |

## 「`α` は `L^m_f` の素元」をどう書いたか(逸脱の記録)

★**`𝒪_{K̂^m_f}` を作っていない。** 中間体 `K̂^m_f ⊆ ℂ_K` の整数環を立てると
`adjoinField`/`adjoinIntegers` の境界に入り(`tools/lean-idioms.md` #69、実測 212 秒 timeout)、
かつ `𝒪_{ℂ_K}` は DVR ではない(値群が稠密。`ClosureCompletion.lean` の「作らなかったもの」)。
★代わりに**ノルムの言葉**で書いた:

> `‖α‖ ^ [K̂^m_f : K̂^{ur}] = ‖π‖`

`K̂^m_f/K̂^{ur}` は完全分岐(`e = [K̂^m_f : K̂^{ur}]`)なので、これは
`v(α) = v(π)/e`、すなわち「`α` が値群を生成する」＝素元、と同値である。
★木の `TotallyRamified.lean` も同じ書き方(`isTotallyRamifiedAdjoin_of_norm_pow_eq`)を採る。

## 退化の自己検査

* ★★**`α ∈ Λ^×_{f,n}`(位数がちょうど `π^n`)を落とすと偽**——`α ∈ Λ_{n−1}` では
  `minpoly_{K̂^{ur}}(ι α) ≠ ψ_n` で、`N(−α)` は `ψ_{n'}` の定数項(`n' < n`)になる。
  本ファイルはどの宣言でも `hxψ : x ∈ iteratedLubinTatePsiTorsionPoints … n hn` を要求する。
* ★★**`−α` の符号は落とせない**——`N(α) = (−1)^{deg} ψ_n(0)` であり、
  `deg = q^n − q^{n−1}` は `q` が偶数(`p = 2`)のとき奇数になりうる。
  そのとき `N(α) = −ψ_n(0) ≠ ψ_n(0)`。★`norm_neg_powerBasis_gen` の 2 つの `(−1)^{dim}` が
  打ち消し合うのは `−α` の側だけである。
* ★**`n ≥ 1`(`hn`)を落とすと `ψ_n` が定義されない**(`iteratedLubinTatePsi` は `1 ≤ n` を要求)。
* ★★**剰余体が代数閉(＝無限)であることを忘れない**——`𝒪_{K̂^{ur}}` の剰余体は無限なので
  `[Fintype (ResidueField 𝒪_{K̂^{ur}})]` を要求する在庫は使えない(決定 D28)。
  ★本ファイルはそれを 1 つも使っていない(§4-a・§4-b と同じ)。
* ★**`U'_n(0) = 1` を落とすと `norm_neg_eq_uniformizer` は偽**(冒頭の反例)。
  ★逆に `norm_neg_eq_coeff_zero_psiCompletion` と `exists_unit_norm_neg` は無仮定で正しい。

## 逸脱(記録)

1. **`L = K̂^{ur}` に固定した。** 原典 §4.1 は「a complete unramified extension `L` of `K`」と
   一般に書くが、本ファイルは `L = K̂^{ur}` の場合だけを扱う(§4-a・§4-b と同じ逸脱)。
2. ★★**`N(−α) = π^{ϕ^{m−1}}` を「`π` と単元倍だけ違う」に弱めた**(上の「原典との差」)。
   原典の字面は `f` が多項式のときだけ正しく、木の `f` は冪級数である。
   ★原典どおりの等号は `norm_neg_eq_uniformizer` に仮定 `U'_n(0) = 1` として残してある。
3. **「素元」をノルムの言葉で書いた**(上の「どう書いたか」)。
4. `exists_uniformizer_norm_neg` は原始点 `α` の選択を隠した形である。原典の主張は
   「任意の `α ∈ µ^×_{f,m}` について」なので、選択を明示した形
   (`norm_neg_eq_coeff_zero_psiCompletion` など)の方が原典に忠実である。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal Valued Classical

/-! ## 0. 抽象核

分岐・付値・Lubin-Tate の語彙が 1 つも出てこない部分。原典の設定に依らない。 -/

section AbstractCore

open Polynomial in
/-- ★★★★★★★★★★★★★★★★**`N(−pb.gen)` は最小多項式の定数項そのもの**。

mathlib の `PowerBasis.norm_gen_eq_coeff_zero_minpoly` は
`N(pb.gen) = (−1)^{dim} · (minpoly).coeff 0` である。一方
`N(−1) = N(algebraMap (−1)) = (−1)^{dim}` なので、`−pb.gen` を取ると
2 つの `(−1)^{dim}` が打ち消し合う。

★**原典が `N(−α)`(`N(α)` ではなく)と書く理由がこれである。** -/
theorem norm_neg_powerBasis_gen {R S : Type*} [CommRing R] [Ring S] [Algebra R S]
    (pb : PowerBasis R S) :
    Algebra.norm R (-pb.gen) = (minpoly R pb.gen).coeff 0 := by
  have hneg : (-pb.gen : S) = algebraMap R S (-1) * pb.gen := by
    rw [map_neg, map_one, neg_one_mul]
  have hsq : ((-1 : R) ^ pb.dim) * ((-1 : R) ^ pb.dim) = 1 := by
    rw [← pow_add, ← two_mul, pow_mul]; simp
  rw [hneg, map_mul, Algebra.norm_algebraMap_of_basis pb.basis,
    Algebra.PowerBasis.norm_gen_eq_coeff_zero_minpoly pb, Fintype.card_fin, ← mul_assoc, hsq,
    one_mul]

open Polynomial in
/-- `F⟮α⟯` の生成元についての `norm_neg_powerBasis_gen`。 -/
theorem norm_neg_adjoinSimple_gen {F M : Type*} [Field F] [Field M] [Algebra F M] {α : M}
    (hα : IsIntegral F α) :
    Algebra.norm F (-(IntermediateField.AdjoinSimple.gen F α)) = (minpoly F α).coeff 0 := by
  have h := norm_neg_powerBasis_gen (IntermediateField.adjoin.powerBasis hα)
  rw [IntermediateField.adjoin.powerBasis_gen hα, IntermediateField.minpoly_gen] at h
  exact h

open Polynomial in
/-- ★★★★★★★★★★★★**中間体 `E` が `F⟮α⟯` と一致するとき、`E` の中の `α` の
ノルムは最小多項式の定数項**。

`E = F⟮α⟯` は**等式**でしか与えられないことが多い(木では
`lubinTateCompletionField_eq_adjoin_simple`)ので、`IntermediateField.equivOfEq` で
`E ≃ₐ[F] F⟮α⟯` を作り、`Algebra.norm_eq_of_algEquiv` で移す。 -/
theorem norm_neg_eq_coeff_zero_minpoly {F M : Type*} [Field F] [Field M] [Algebra F M] {α : M}
    (hα : IsIntegral F α) {E : IntermediateField F M}
    (hE : E = IntermediateField.adjoin F ({α} : Set M)) (β : ↥E) (hβ : (β : M) = α) :
    Algebra.norm F (-β) = (minpoly F α).coeff 0 := by
  have h1 : (IntermediateField.equivOfEq hE) β = IntermediateField.AdjoinSimple.gen F α := by
    apply Subtype.ext
    simpa using hβ
  have h2 := Algebra.norm_eq_of_algEquiv (IntermediateField.equivOfEq hE) (-β)
  rw [map_neg, h1] at h2
  rw [← h2, norm_neg_adjoinSimple_gen hα]

end AbstractCore

variable {p : ℕ} [Fact p.Prime]

section Concrete

variable (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))

/-! ## 1. `ψ_n` の定数項 —— 「`π` と単元倍だけ違う」 -/

/-- `ψ_n` を `K̂^{ur}` 係数で見たものはモニック。 -/
theorem monic_psiCompletion (n : ℕ) (hn : 1 ≤ n) :
    (psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn).Monic :=
  (((isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic.map
    (baseIntHom K)).map _)

/-- `ψ_n` を `K̂^{ur}` 係数で見たものの係数は、`𝒪_K → 𝒪_{K̂^{ur}} → K̂^{ur}` の像。 -/
theorem coeff_psiCompletion (n : ℕ) (hn : 1 ≤ n) (i : ℕ) :
    (psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn).coeff i
      = algebraMap ↥(unramifiedCompletionInt K) (unramifiedCompletion K)
          (baseIntHom K ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff i)) := by
  rw [psiCompletion, Polynomial.coeff_map, Polynomial.coeff_map]

/-- `ι α`(`α ∈ Λ^×_{f,n}`)は `K̂^{ur}` 上整——`ψ_n` を証拠にするだけ。 -/
theorem isIntegral_closureCompletionCoe_torsion (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    IsIntegral (unramifiedCompletion K) (closureCompletionCoe K x) :=
  ⟨psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn,
    monic_psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn,
    aeval_closureCompletionCoe_map_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ⟩

/-- ★★★★★★★★★★★★**`ψ_n.coeff 0 = π · u`**(`u ∈ 𝒪_K^×`)。

★**これが原典 `N(−α) = π^{ϕ^{m−1}}` の、木の設定での正しい形である**(冒頭参照)。
`r_n = ψ_n · U'_n` の定数項を比べた `iteratedLubinTatePsi_coeff_zero_mul`
(`ψ_n.coeff 0 · U'_n(0) = π`、木に在る)と、`U'_n` が単元冪級数であること
(`isUnit_iteratedLubinTateStepU`)から、`U'_n(0)` は `𝒪_K` の単元である。 -/
theorem exists_unit_coeff_zero_iteratedLubinTatePsi (n : ℕ) (hn : 1 ≤ n) :
    ∃ u : (𝒪[K.carrier])ˣ,
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0
        = π * (u : 𝒪[K.carrier]) := by
  obtain ⟨u, hu⟩ := PowerSeries.isUnit_constantCoeff
    (iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn)
    (isUnit_iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn)
  have h := iteratedLubinTatePsi_coeff_zero_mul hq hπmax hπne0 f hf0 hf1 hf n hn
  rw [← hu] at h
  exact ⟨u⁻¹, (Units.eq_mul_inv_iff_mul_eq u).mpr h⟩

/-- `(ψ_n.coeff 0) = (π)` —— イデアルとしては `π` と区別が付かない。 -/
theorem span_coeff_zero_iteratedLubinTatePsi (n : ℕ) (hn : 1 ≤ n) :
    Ideal.span {(iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0}
      = Ideal.span ({π} : Set (𝒪[K.carrier])) := by
  obtain ⟨u, hu⟩ := exists_unit_coeff_zero_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn
  rw [hu]
  exact Ideal.span_singleton_mul_right_unit u.isUnit π

/-! ## 2. `N_{K̂^m_f/K̂^{ur}}(−α)` -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`N_{K̂^m_f/K̂^{ur}}(−β) = ψ_n.coeff 0`**(`β ∈ K̂^m_f` の像が `ι α`、`α ∈ Λ^×_{f,n}`)。

Yoshida 2008 Proposition 4.4(ii) の `N_{L^m_f/L}(−α)` の部分、`L = K̂^{ur}` 版。
★**無仮定で厳密に成り立つ**(原典の `= π^{ϕ^{m−1}}` は木の一般の `f` では偽。冒頭参照)。

材料は 2 つだけ:
* 抽象核 `norm_neg_eq_coeff_zero_minpoly`(`N(−β) = minpoly の定数項`)
* §4-a・§4-b の `lubinTateCompletionField_eq_adjoin_simple`(`K̂^m_f = K̂^{ur}(ι α)`)と
  `minpoly_closureCompletionCoe`(`minpoly_{K̂^{ur}}(ι α) = ψ_n`) -/
theorem norm_neg_eq_coeff_zero_psiCompletion (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (β : ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n))
    (hβ : (β : closureCompletion K) = closureCompletionCoe K x) :
    Algebra.norm (unramifiedCompletion K) (-β)
      = (psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0 := by
  rw [norm_neg_eq_coeff_zero_minpoly
      (isIntegral_closureCompletionCoe_torsion K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ)
      (lubinTateCompletionField_eq_adjoin_simple K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn)
      β hβ,
    minpoly_closureCompletionCoe K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ]

/-- `N(−β) · U'_n(0) = π`——原典の等号が壊れる場所を明示した形。 -/
theorem norm_neg_mul_constantCoeff_stepU (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (β : ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n))
    (hβ : (β : closureCompletion K) = closureCompletionCoe K x) :
    Algebra.norm (unramifiedCompletion K) (-β)
        * algebraMap ↥(unramifiedCompletionInt K) (unramifiedCompletion K)
            (baseIntHom K (PowerSeries.constantCoeff
              (iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn)))
      = algebraMap ↥(unramifiedCompletionInt K) (unramifiedCompletion K) (baseIntHom K π) := by
  rw [norm_neg_eq_coeff_zero_psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn β hβ,
    coeff_psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn 0, ← map_mul, ← map_mul,
    iteratedLubinTatePsi_coeff_zero_mul hq hπmax hπne0 f hf0 hf1 hf n hn]

/-- ★★★★★★★★★★★★★★★★★★★★
**`N_{K̂^m_f/K̂^{ur}}(−β) = π · u`**(`u ∈ 𝒪_K^×`)。

★これが原典 Proposition 4.4(ii) の `N_{L^m_f/L}(−α) = π^{ϕ^{m−1}}` の、
木の設定(`f` が冪級数)での正しい形である。`L = K̂^{ur}` では `π ∈ 𝒪_K` が `ϕ` で
固定されるので **`π^{ϕ^{m−1}} = π`** であり、残る差は単元 `u` だけである。 -/
theorem exists_unit_norm_neg (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (β : ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n))
    (hβ : (β : closureCompletion K) = closureCompletionCoe K x) :
    ∃ u : (𝒪[K.carrier])ˣ,
      Algebra.norm (unramifiedCompletion K) (-β)
        = algebraMap ↥(unramifiedCompletionInt K) (unramifiedCompletion K)
            (baseIntHom K (π * (u : 𝒪[K.carrier]))) := by
  obtain ⟨u, hu⟩ := exists_unit_coeff_zero_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn
  refine ⟨u, ?_⟩
  rw [norm_neg_eq_coeff_zero_psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn β hβ,
    coeff_psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn 0, hu]

/-- ★**原典の字面どおりの `N_{L^m_f/L}(−α) = π^{ϕ^{m−1}}`**(`L = K̂^{ur}` では
右辺は `π`)——ただし **`U'_n(0) = 1` を仮定する**。

★この仮定は原典の設定(`f ∈ 𝒪_L[X]` が多項式)では成り立つ:そのとき
`r_n = [π^n]_f/[π^{n−1}]_f` 自身がモニックで `mod 𝔪` が `X^{deg}` になり、
Weierstrass 分解の単元部分は `1` である。
★★**木の一般の `f`(冪級数)では偽である**(冒頭の反例)。仮定を外した形は
`exists_unit_norm_neg`。 -/
theorem norm_neg_eq_uniformizer (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (β : ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n))
    (hβ : (β : closureCompletion K) = closureCompletionCoe K x)
    (hU : PowerSeries.constantCoeff (iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn) = 1) :
    Algebra.norm (unramifiedCompletion K) (-β)
      = algebraMap ↥(unramifiedCompletionInt K) (unramifiedCompletion K) (baseIntHom K π) := by
  have h := norm_neg_mul_constantCoeff_stepU K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn β hβ
  rwa [hU, map_one, map_one, mul_one] at h

/-- ★★★★★★★★★★★★**`‖N_{K̂^m_f/K̂^{ur}}(−β)‖ = ‖π‖`**——
`N(−β)` は `𝒪_{K̂^{ur}}` の素元である(単元の差はノルムに見えない)。 -/
theorem norm_norm_neg (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (β : ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n))
    (hβ : (β : closureCompletion K) = closureCompletionCoe K x) :
    ‖Algebra.norm (unramifiedCompletion K) (-β)‖ = ‖π‖ := by
  rw [norm_neg_eq_coeff_zero_psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn β hβ,
    coeff_psiCompletion K hq hπmax hπne0 f hf0 hf1 hf n hn 0]
  show ‖((baseIntHom K ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0)
    : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)‖ = ‖π‖
  rw [norm_baseIntHom]
  exact norm_iteratedLubinTatePsi_coeff_zero K hq hπmax hπne0 f hf0 hf1 hf n hn

/-! ## 3. 「`α` は `K̂^m_f` の素元」 -/

/-- `α ∈ Λ^×_{f,n}` は `K` 係数に写した `ψ_n` の根(`aeval` の形)。 -/
theorem aeval_carrier_iteratedLubinTatePsi (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    Polynomial.aeval x (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)) = 0 := by
  have hroot0 : Polynomial.eval x
      ((iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map
        (algebraMap 𝒪[K.carrier] K.closure)) = 0 :=
    (Polynomial.mem_roots'.mp (Multiset.mem_toFinset.mp hxψ)).2
  rw [Polynomial.aeval_def, Polynomial.eval₂_map, ← Polynomial.eval_map,
    ← IsScalarTower.algebraMap_eq]
  exact hroot0

/-- ★★★★★★★★**`‖α‖^{q^n − q^{n−1}} = ‖π‖`**(`α ∈ Λ^×_{f,n}`)。

木の `spectralNorm_root_iteratedLubinTatePsi`(`LubinTatePsiNorm.lean`、
根のスペクトルノルムは `‖π‖^{1/deg}`)と `norm_eq_spectralNorm_closure`
(`K^{al}` のノルムはスペクトルノルム、`rfl`)を繋ぎ、`rpow` を `npow` に戻すだけ。 -/
theorem norm_pow_torsion_point (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    ‖x‖ ^ ((pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1)) = ‖π‖ := by
  have hdpos : 0 < (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := sub_pow_pred_pos K hq n hn
  have hspec := spectralNorm_root_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x
    (aeval_carrier_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ)
  rw [← norm_eq_spectralNorm_closure,
    natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn] at hspec
  have hd : (((pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) : ℕ) : ℝ) ≠ 0 := by
    exact_mod_cast hdpos.ne'
  rw [hspec, ← Real.rpow_natCast (‖π‖ ^ (1 / (((pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) : ℕ) : ℝ))) _,
    ← Real.rpow_mul (norm_nonneg _), one_div_mul_cancel hd, Real.rpow_one]

/-- ★★★★★★★★★★★★★★★★★★★★★★
**`α` は `K̂^m_f` の素元**(Yoshida 2008 Proposition 4.4(ii)、`L = K̂^{ur}` 版):

`‖α‖ ^ [K̂^m_f : K̂^{ur}] = ‖π‖`

★`K̂^m_f/K̂^{ur}` は完全分岐(`e = [K̂^m_f : K̂^{ur}]`)なので、これは
`v(α) = v(π)/e`、すなわち `α` が `K̂^m_f` の値群を生成することを言っている。
★整数環 `𝒪_{K̂^m_f}` を作らない書き方については冒頭の「逸脱」を見よ。 -/
theorem norm_pow_finrank_closureCompletionCoe (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    ‖closureCompletionCoe K x‖ ^ (Module.finrank (unramifiedCompletion K)
        ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n)) = ‖π‖ := by
  rw [finrank_lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n hn,
    norm_closureCompletionCoe]
  exact norm_pow_torsion_point K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ

/-! ## 4. 主結果 -/

/-- `ι α`(`α ∈ Λ^×_{f,n}`)は `K̂^m_f` の元。 -/
theorem closureCompletionCoe_mem_lubinTateCompletionField (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n) :
    closureCompletionCoe K x ∈ lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n := by
  rw [lubinTateCompletionField_eq_adjoin_simple K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn]
  exact IntermediateField.mem_adjoin_simple_self _ _

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Yoshida 2008 Proposition 4.4(ii) の残り**(`L = K̂^{ur}`)を 1 つにまとめた形:

`K̂^m_f` の元 `β` で
1. `K̂^m_f = K̂^{ur}(β)`
2. `‖β‖ ^ [K̂^m_f : K̂^{ur}] = ‖π‖`  (★`β` は `K̂^m_f` の素元)
3. `‖N_{K̂^m_f/K̂^{ur}}(−β)‖ = ‖π‖`  (★`N(−β)` は `𝒪_{K̂^{ur}}` の素元)
4. `N_{K̂^m_f/K̂^{ur}}(−β) = π · u`(`u ∈ 𝒪_K^×`)

を満たすものが存在する。★4 が原典の `N_{L^m_f/L}(−α) = π^{ϕ^{m−1}}` に対応する
(冒頭のとおり、木の `f` は冪級数なので単元 `u` は外せない)。 -/
theorem exists_uniformizer_norm_neg (n : ℕ) (hn : 1 ≤ n) :
    ∃ β : ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n),
      lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n
          = IntermediateField.adjoin (unramifiedCompletion K)
              ({(β : closureCompletion K)} : Set (closureCompletion K))
        ∧ ‖(β : closureCompletion K)‖ ^ (Module.finrank (unramifiedCompletion K)
            ↥(lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n)) = ‖π‖
        ∧ ‖Algebra.norm (unramifiedCompletion K) (-β)‖ = ‖π‖
        ∧ ∃ u : (𝒪[K.carrier])ˣ, Algebra.norm (unramifiedCompletion K) (-β)
            = algebraMap ↥(unramifiedCompletionInt K) (unramifiedCompletion K)
                (baseIntHom K (π * (u : 𝒪[K.carrier]))) := by
  obtain ⟨x, hxψ⟩ :=
    iteratedLubinTatePsiTorsionPoints_nonempty K hq hπmax hπne0 f hf0 hf1 hf n hn
  have hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n := by
    rw [iteratedLubinTateTorsionPoints_eq_union K hq hπmax hπne0 f hf0 hf1 hf n hn,
      Finset.mem_union]
    exact Or.inr hxψ
  refine ⟨⟨closureCompletionCoe K x,
    closureCompletionCoe_mem_lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn⟩,
    lubinTateCompletionField_eq_adjoin_simple K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn,
    norm_pow_finrank_closureCompletionCoe K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ,
    norm_norm_neg K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn _ rfl,
    exists_unit_norm_neg K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn _ rfl⟩

def exists_uniformizer_norm_neg.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 7, item := "Proposition 4.4", sectionId := "prop-4-4" }

end Concrete

end ABC3.Found.PGC
