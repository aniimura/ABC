/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import ABC3.Found.NumberField.RatHeightOne
import Mathlib.RingTheory.Ideal.Norm.AbsNorm
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.NumberTheory.NumberField.DedekindZeta
import Mathlib.RingTheory.DedekindDomain.Ideal.Lemmas
import Mathlib.FieldTheory.Galois.Basic

/-!
# Tchebotarev の密度定理(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> §4, Theorem 10], it follows that [Li : Q] is equal to the maximum of the deg(Li, vi),

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★これは §6 に残る**唯一の在庫不足**であった

`Theorem 6.4, (iv)` は Tchebotarev を **3 箇所**で使う:

| 箇所 | 主張 | 節点 |
|---|---|---|
| (a) | `[L_i:ℚ]` は `deg(L_i, v_i)` の最大値 | `cheb-max-deg` |
| (b) | `p` が `L_i` で完全分解 ⟺ `deg(L_i,v_i) = [L_i:ℚ]` | `cheb-split-count` |
| (c) | 完全分解する素点が包含すれば体が包含する | `cheb-spl-det` |

★★★**迂回路がある** —— (a) は「Galois 拡大で完全分解する素数は無限個」で足り、
これは最小多項式の値の素因子を数える初等的な議論で出る(`cheb-spl-infinite`)。
★**密度を本当に要求するのは (c) だけ**である。

## ★出典と在庫

出典は **[MilneCFT]**(J. S. Milne, *Class Field Theory*、著者が無償公開している講義ノート)の
Chapter V, Theorem 3.23(声明)と Chapter VIII, §7, Theorem 7.4(証明)。
★PDF は `ResearchPaper/0_Source/` に取得済み(**gitignore 済みで公開リポジトリには入らない**)。

★★mathlib の在庫測定(2026-08-20):

| 要るもの | mathlib | 判定 |
|---|---|---|
| Dedekind ζ の `s=1` での留数 | `NumberField.dedekindZeta_residue` | ★★**ある** |
| Dirichlet の算術級数定理 | `Nat.infinite_setOf_prime_and_eq_mod` | ★★**ある**(密度つきではない) |
| `Ideal.absNorm` / `ramificationIdx` / `inertiaDeg` / `primesOver` | 各所 | ★**ある** |
| `Chebotarev` / `Tchebotarev` | 全体を grep して **0 件** | ★**無い** |
| Artin の相互法則 | 類体論は実質不在(`Artin` の 54 件は別物) | ★**無い** |
| Dirichlet 密度 / 射類群 | grep 0 件 | ★**無い** |

★★分解は `ResearchPaper/frdi-decomposition.json` の鎖 `cheb`(22 節点)。

## ★★★★★★**訂正(2026-08-20)—— 律速は Artin の相互法則ではなかった**

当初「(c) は密度を本当に要求し、律速は Artin の相互法則」と書いたが、
それは**共役類つきの一般形**についてであって、
**[FrdI] が使う 3 箇所はいずれも「完全分解」の場合に収まる**。

★完全分解の場合の密度 `1/[L:K]` は **Frobenius/Dedekind の古典的な論法**で出る:

1. `log ζ_L(s) = Σ_{f(𝔭)=1} N𝔭^{-s} + O(1)`（`f ≥ 2` や `m ≥ 2` の項は有界）
2. `L/K` が Galois なら `p` の上の `f` はすべて等しいので、
   `#{𝔭 ∣ p, f=1}` は完全分解なら `[L:K]`、それ以外は `0`
3. `ζ_L` は `s=1` で単純極（**mathlib 在庫** `dedekindZeta_residue`）なので
   `log ζ_L(s) ~ log(1/(s-1))`
4. ゆえに `Σ_{p ∈ Spl(L)} p^{-s} ~ (1/[L:K])·log(1/(s-1))`

★★**使うのは Dedekind ζ の極だけ**で、Artin の相互法則も Hecke L 関数も要らない。
★★★したがって **§6 に類体論は要らない**。

## ★★★★★★★★**判断 D11(2026-09-06)—— 畳んで `Found` へ委譲した**

★**このファイルの数学は `Found/NumberField/`(21 ファイル・4,864 行・sorry 0)で
既に閉じていた**。本ファイルの 3 つの `sorry` は、その結果への**薄い橋**に置き換えた:

| 本ファイルの宣言 | 委譲先(`Found`、sorry 0) |
|---|---|
| `HasDirichletDensity` | `ABC3.Found.NF.HasDirichletDensity`(**重複定義を解消**) |
| `chebotarev_splitsCompletely` | `hasDirichletDensity_splitsCompletely` ← `tendsto_splQ_div_log` |
| `nonempty_algHom_of_splitsCompletely_subset` | `nonempty_algHom_of_splitsCompletely_subset` ← `le_of_SplQ_subset` |
| `infinite_splitsCompletely` | `infinite_splitsCompletely_heightOne` ← `infinite_splitsCompletely_of_isGalois` |

語彙の橋(`HeightOneSpectrum (𝓞 ℚ)` と `Nat.Primes` の同一視)は
`Found/NumberField/RatHeightOne.lean` に置いた
(mathlib の `Rat.HeightOneSpectrum.primesEquiv` 経由)。

### ★★★逸脱の記録(判断 D11、2026-09-06)

★**底を一般の `K` から `ℚ` に固定した**。
これは **statement を変える逸脱**であり、理由は 2 つある:

1. **原典が要求するのは `ℚ` の場合だけ**である ——
   [FrdI] Theorem 6.4, (iv) の `L_i` は `ℚ` 上の数体であり、
   原文の Tchebotarev の 3 つの使い方(上の (a)(b)(c))はいずれも底が `ℚ` である。
2. **一般の底 `K` の版は本プロジェクトの分解に入っていない** ——
   `Found/NumberField/` が閉じたのは底 `ℚ` の場合であり、
   一般の `K` へ上げる作業は原典には要らない。

★**消費者は 0 だったので、この逸脱で壊れる下流は無い**(2026-09-06 実測。
唯一 `Skeleton/FrdI/Thm64Deg.lean` が本ファイルを import しているが、
そこで使うのは `Found.NF.*` であって本ファイルの宣言ではない)。

### ★重複定義の解消(判断 D11 の 3 番目)

`HasDirichletDensity` は `Found/NumberField/DirichletDensity.lean` と**本文まで同一**だった
(片方が `nhdsWithin 1 (Set.Ioi 1)`、片方が `𝓝[>] 1` と書いていただけ)。
★**`abbrev` にして `Found` の 1 本に寄せた**。定義は 1 つしか無い。

## ★本ファイルが置くもの

★**共役類つきの完全な Chebotarev は Frobenius 元を要求する**ので、
ここでは **Frobenius を使わない形**(完全分解する素点の密度 `1/[L:ℚ]`)を置く。
★これは共役類 `C = {1}` の場合であり、`Theorem 6.4, (iv)` が使う (c) を出すのに十分である。
-/

namespace ABC3.Skeleton.Cheb

open NumberField IsDedekindDomain ABC3.Meta
open scoped NumberField

/-! ## ★1. Dirichlet 密度 -/

/-- ★★**Dirichlet 密度** —— `i(T) = lim_{s→1+} (Σ_{𝔭∈T} N𝔭^{-s}) / log(1/(s-1))`。

★★**定義の実体は `Found/NumberField/DirichletDensity.lean` にある**(判断 D11、2026-09-06)。
以前はここに同じ本文の写しが置かれていたが、**重複定義を解消して `Found` に寄せた**。 -/
abbrev HasDirichletDensity (K : Type*) [Field K] [NumberField K]
    (S : Set (HeightOneSpectrum (𝓞 K))) (d : ℝ) : Prop :=
  ABC3.Found.NF.HasDirichletDensity S d

/-! ## ★2. 完全分解 -/

/-- ★**`𝔭` が `L` で完全分解する** —— `𝔭` の上にある素イデアルがちょうど `[L:K]` 個。 -/
def SplitsCompletely (K L : Type*) [Field K] [Field L] [NumberField K] [NumberField L]
    [Algebra K L] (𝔭 : HeightOneSpectrum (𝓞 K)) : Prop :=
  (Ideal.primesOver 𝔭.asIdeal (𝓞 L)).ncard = Module.finrank K L

/-! ## ★3. Chebotarev(Frobenius を使わない形) -/

/-- ★★★**Chebotarev の密度定理**(完全分解の場合)——
`L/ℚ` が Galois なら、`L` で完全分解する素点の Dirichlet 密度は `1/[L:ℚ]`。

原文 (FrdI p.116):
> §4, Theorem 10], it follows that [Li : Q] is equal to the maximum of the deg(Li, vi),

★★[MilneCFT] Chapter VIII, §7, Theorem 7.4 の共役類 `C = {1}` の場合。
★★★**中身は `Found` にある**(判断 D11、2026-09-06)——
`Found.NF.hasDirichletDensity_splitsCompletely`(sorry 0)への薄い橋である。
その中身は `tendsto_splQ_div_log`、すなわち **Dedekind ζ の単純極だけ**を使う証明であり、
Artin の相互法則も Hecke L 関数も使わない。
★**逸脱**: 底を `ℚ` に固定した(上の「逸脱の記録」を見よ)。 -/
theorem chebotarev_splitsCompletely (L : Type) [Field L] [NumberField L] [IsGalois ℚ L] :
    HasDirichletDensity ℚ {𝔭 | SplitsCompletely ℚ L 𝔭}
      (1 / (Module.finrank ℚ L : ℝ)) :=
  ABC3.Found.NF.hasDirichletDensity_splitsCompletely L

/-- ★★★★**完全分解する素点の集合が体を決める**([MilneCFT] Chapter V, Theorem 3.25)。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

★★**`Theorem 6.4, (iv)` が最後に使うのはこれ**である ——
原文「[again by Tchebotarev's density theorem — cf. [NSW], Theorem 12.2.5] that `L₁ ⊆ L₂`」。
★★★**中身は `Found` にある**(判断 D11、2026-09-06)——
`Found.NF.nonempty_algHom_of_splitsCompletely_subset`(sorry 0)への薄い橋である。
その中身は `Spl(LM) = Spl(L) ∩ Spl(M)` と密度の一意性で `[LM:ℚ] = [M:ℚ]` を出す
`le_of_SplQ_subset` である。
★**逸脱**: 底を `ℚ` に固定した(上の「逸脱の記録」を見よ)。 -/
theorem nonempty_algHom_of_splitsCompletely_subset (L M : Type) [Field L] [Field M]
    [NumberField L] [NumberField M] [IsGalois ℚ L] [IsGalois ℚ M]
    (h : {𝔭 | SplitsCompletely ℚ L 𝔭} ⊆ {𝔭 | SplitsCompletely ℚ M 𝔭}) :
    Nonempty (M →ₐ[ℚ] L) :=
  ABC3.Found.NF.nonempty_algHom_of_splitsCompletely_subset L M h

/-! ## ★4. 迂回路 —— 密度を使わない部分 -/

/-- ★★**Galois 拡大で完全分解する素点は無限個**。

★★**これは全 Chebotarev を使わずに出る** —— `L = ℚ(θ)`、`f` を `θ` の最小多項式とすると、
`𝔭 ∤ disc(f)` について「`𝔭` が完全分解 ⟺ `f` が `mod 𝔭` で根を持つ」(`L/ℚ` が Galois のとき)。
「`f` の値の素因子は無限個」は初等的である。
★★★**中身は `Found` にある**(判断 D11、2026-09-06)——
`Found.NF.infinite_splitsCompletely_heightOne`(sorry 0、**仮定なし**)への薄い橋である。
★**逸脱**: 底を `ℚ` に固定した(上の「逸脱の記録」を見よ)。 -/
theorem infinite_splitsCompletely (L : Type) [Field L] [NumberField L] [IsGalois ℚ L] :
    {𝔭 : HeightOneSpectrum (𝓞 ℚ) | SplitsCompletely ℚ L 𝔭}.Infinite :=
  ABC3.Found.NF.infinite_splitsCompletely_heightOne L

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def HasDirichletDensity.src : Source :=
  { paper := "FrdI", pdfPage := 116, item := "Theorem 6.4, (iv) — Dirichlet 密度",
    sectionId := "frdi-thm-6-4" }

def SplitsCompletely.src : Source :=
  { paper := "FrdI", pdfPage := 116, item := "Theorem 6.4, (iv) — 完全分解",
    sectionId := "frdi-thm-6-4" }

def chebotarev_splitsCompletely.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — Tchebotarev の密度定理(完全分解の場合)",
    sectionId := "frdi-thm-6-4" }

/-- ★★分解は `ResearchPaper/frdi-decomposition.json` の鎖 `cheb`(19 節点・8 層)。
★★訂正(2026-08-20): [FrdI] が要求するのは**完全分解の場合だけ**で、それは Dedekind ζ の極（mathlib 在庫）で出る。Artin の相互法則は**共役類つきの一般形**にしか要らない。
★★★判断 D11(2026-09-06): 本宣言は `Found` への薄い橋になった。以下の債務は**すべて `Found` 側で解消済み**である。 -/
def chebotarev_splitsCompletely.needs : List ProofObligation :=
  [ .otherPaper "[MilneCFT]" "Chapter VIII, §7, Theorem 7.4(Chebotarev の密度定理)" 274,
    .otherPaper "[MilneCFT]" "Chapter V, §3, Theorem 3.5(Artin の相互法則)" 178,
    .citation "[mathlib]" "NumberField.dedekindZeta_residue(Dedekind ζ の s=1 での留数)"
      (.inMathlib "NumberField.dedekindZeta_residue") 116,
    .citation "[mathlib]" "Nat.infinite_setOf_prime_and_eq_mod(Dirichlet の算術級数定理)"
      (.inMathlib "Nat.infinite_setOf_prime_and_eq_mod") 116,
    .citation "[mathlib]" "Chebotarev / Artin の相互法則 / Dirichlet 密度 / 射類群"
      (.absent "lean/.lake/packages/mathlib/Mathlib/ 全体を Chebotarev|Tchebotarev|ArtinReciprocity|artinMap|DirichletDensity|RayClassGroup|SplitsCompletely|decompositionGroup|analyticDensity|naturalDensity|primeDensity で grep、いずれも 0 件(2026-09-06 再測。node tools/absent-recheck.mjs の try オプションで再実行できる)。★規約つきに直した(2026-09-06): re:`Chebotarev|Tchebotarev|ArtinReciprocity|artinMap|DirichletDensity|RayClassGroup|SplitsCompletely|decompositionGroup|analyticDensity|naturalDensity|primeDensity`→0") 116,
    .citation "[mathlib]" "Frobenius 元(★2026-08-20 の「不在」は誤り。名前が FrobeniusElement ではないので旧 regex に当たらなかった。流儀は Gal(L/K) ではなく、G が S に作用し R が固定環という RingTheory/Invariant の形)"
      (.inMathlib "IsArithFrobAt") 116,
    .citation "[ABC3]" "hasDirichletDensity_splitsCompletely(本宣言の委譲先、sorry 0)"
      (.inProject "ABC3" "ABC3.Found.NF.hasDirichletDensity_splitsCompletely") 116,
    .citation "[ABC3]" "tendsto_splQ_div_log(類体論を使わない Chebotarev、完全分解の場合、sorry 0)"
      (.inProject "ABC3" "ABC3.Found.NF.tendsto_splQ_div_log") 116,
    .derivation
      "鎖 cheb: 密度の定義 → Frobenius 元 → 射類群 → Artin 写像 → 相互法則 → L(1,χ)≠0 → 射類ごとの等分布 → アーベルの場合 → 巡回への還元 → 2 つの計数((a) 全単射・(b) d/(cf) 対 1)→ 本体" 116,
    .implicitStep
      "★原文は「by Tchebotarev's density theorem [cf., e.g., [Lang2], Chapter VIII, §4, Theorem 10]」と外部文献へ送る" 116,
    .implicitStep
      ("★判断 D11(2026-09-06): 底を一般の K から ℚ に固定した(逸脱)。"
       ++ "原典 Theorem 6.4, (iv) の底は ℚ であり、一般の底の Chebotarev は本プロジェクトの分解に入っていない。"
       ++ "★消費者は 0 だったので下流に影響は無い") 116 ]

def nonempty_algHom_of_splitsCompletely_subset.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解する素点の包含が体の包含を与える",
    sectionId := "frdi-thm-6-4" }

def nonempty_algHom_of_splitsCompletely_subset.needs : List ProofObligation :=
  [ .otherPaper "[MilneCFT]" "Chapter V, §3, Theorem 3.25(Spl が体を決める)" 178,
    .citation "[ABC3]" "chebotarev_splitsCompletely"
      (.inProject "ABC3" "ABC3.Skeleton.Cheb.chebotarev_splitsCompletely") 116,
    .citation "[ABC3]" "nonempty_algHom_of_splitsCompletely_subset(本宣言の委譲先、sorry 0)"
      (.inProject "ABC3" "ABC3.Found.NF.nonempty_algHom_of_splitsCompletely_subset") 116,
    .citation "[ABC3]" "le_of_SplQ_subset(Bauer——Spl の包含から体の包含、sorry 0)"
      (.inProject "ABC3" "ABC3.Found.NF.le_of_SplQ_subset") 116,
    .derivation "Spl(LM/K) = Spl(L/K) ∩ Spl(M/K) と密度の比較で [LM:K] = [M:K]" 116,
    .implicitStep
      "★原文は [NSW] Theorem 12.2.5 へ送る(他論文)" 116,
    .implicitStep
      ("★判断 D11(2026-09-06): 底を一般の K から ℚ に固定した(逸脱)。"
       ++ "Found 側の le_of_SplQ_subset が底 ℚ の版だからであり、原典が要求するのもそれである") 116 ]

def infinite_splitsCompletely.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解する素点は無限個",
    sectionId := "frdi-thm-6-4" }

/-- ★★**これは密度を要求しない** —— 迂回路である。
★★★判断 D11(2026-09-06): `Found` 側で**仮定なし**に閉じている。 -/
def infinite_splitsCompletely.needs : List ProofObligation :=
  [ .otherPaper "[MilneANT]" "Chapter IV(分解群・惰性群・Frobenius の基礎)" 100,
    .citation "[mathlib]" "Ideal.sum_ramification_inertia(Σ e_i f_i = [L:K])"
      (.inMathlib "Ideal.sum_ramification_inertia") 116,
    .citation "[ABC3]" "infinite_splitsCompletely_heightOne(本宣言の委譲先、sorry 0)"
      (.inProject "ABC3" "ABC3.Found.NF.infinite_splitsCompletely_heightOne") 116,
    .citation "[ABC3]" "infinite_splitsCompletely_of_isGalois(Chebotarev を使わない迂回路、仮定なし、sorry 0)"
      (.inProject "ABC3" "ABC3.Found.NF.infinite_splitsCompletely_of_isGalois") 116,
    .derivation
      "L = K(θ)、f を θ の最小多項式として、𝔭 ∤ disc(f) では「完全分解 ⟺ f が mod 𝔭 で根を持つ」。f の値の素因子が無限個であることは初等的(Schur の議論)" 116 ]

end ABC3.Skeleton.Cheb
