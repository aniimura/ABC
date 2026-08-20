/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
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

## ★★これは §6 に残る**唯一の在庫不足**である

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
| Dirichlet 密度 / Frobenius 元 / 射類群 | grep 0 件 | ★**無い** |

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

## ★本ファイルが置くもの

★**共役類つきの完全な Chebotarev は Frobenius 元(`cheb-frob`、mathlib に無い)を要求する**ので、
ここでは **Frobenius を使わない形**(完全分解する素点の密度 `1/[L:K]`)を置く。
★これは共役類 `C = {1}` の場合であり、`Theorem 6.4, (iv)` が使う (c) を出すのに十分である。
-/

namespace ABC3.Skeleton.Cheb

open NumberField IsDedekindDomain ABC3.Meta
open scoped NumberField

/-! ## ★1. Dirichlet 密度 -/

/-- ★★**Dirichlet 密度** —— `i(T) = lim_{s→1+} (Σ_{𝔭∈T} N𝔭^{-s}) / log(1/(s-1))`。

★mathlib に `DirichletDensity` は **0 件**(2026-08-20 実測)。 -/
def HasDirichletDensity (K : Type*) [Field K] [NumberField K]
    (S : Set (HeightOneSpectrum (𝓞 K))) (d : ℝ) : Prop :=
  Filter.Tendsto
    (fun s : ℝ =>
      (∑' 𝔭 : S, ((Ideal.absNorm (𝔭 : HeightOneSpectrum (𝓞 K)).asIdeal : ℝ)) ^ (-s))
        / Real.log (1 / (s - 1)))
    (nhdsWithin 1 (Set.Ioi 1)) (nhds d)

/-! ## ★2. 完全分解 -/

/-- ★**`𝔭` が `L` で完全分解する** —— `𝔭` の上にある素イデアルがちょうど `[L:K]` 個。 -/
def SplitsCompletely (K L : Type*) [Field K] [Field L] [NumberField K] [NumberField L]
    [Algebra K L] (𝔭 : HeightOneSpectrum (𝓞 K)) : Prop :=
  (Ideal.primesOver 𝔭.asIdeal (𝓞 L)).ncard = Module.finrank K L

/-! ## ★3. Chebotarev(Frobenius を使わない形) -/

/-- ★★★**Chebotarev の密度定理**(完全分解の場合)——
`L/K` が Galois なら、`L` で完全分解する素点の Dirichlet 密度は `1/[L:K]`。

原文 (FrdI p.116):
> §4, Theorem 10], it follows that [Li : Q] is equal to the maximum of the deg(Li, vi),

★★[MilneCFT] Chapter VIII, §7, Theorem 7.4 の共役類 `C = {1}` の場合。
★共役類つきの一般形は Frobenius 元(鎖 `cheb` の `cheb-frob`)を要求する。 -/
theorem chebotarev_splitsCompletely (K L : Type) [Field K] [Field L]
    [NumberField K] [NumberField L] [Algebra K L] [IsGalois K L] :
    HasDirichletDensity K {𝔭 | SplitsCompletely K L 𝔭}
      (1 / (Module.finrank K L : ℝ)) := by
  sorry

/-- ★★★★**完全分解する素点の集合が体を決める**([MilneCFT] Chapter V, Theorem 3.25)。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

★★**`Theorem 6.4, (iv)` が最後に使うのはこれ**である ——
原文「[again by Tchebotarev's density theorem — cf. [NSW], Theorem 12.2.5] that `L₁ ⊆ L₂`」。
★証明は `Spl(LM/K) = Spl(L/K) ∩ Spl(M/K)` と Chebotarev(密度から次数が読める)で
`[LM:K] = [M:K]` を出す。 -/
theorem nonempty_algHom_of_splitsCompletely_subset (K L M : Type) [Field K] [Field L] [Field M]
    [NumberField K] [NumberField L] [NumberField M]
    [Algebra K L] [Algebra K M] [IsGalois K L] [IsGalois K M]
    (h : {𝔭 | SplitsCompletely K L 𝔭} ⊆ {𝔭 | SplitsCompletely K M 𝔭}) :
    Nonempty (M →ₐ[K] L) := by
  sorry

/-! ## ★4. 迂回路 —— 密度を使わない部分 -/

/-- ★★**Galois 拡大で完全分解する素点は無限個**。

★★**これは全 Chebotarev を使わずに出る** —— `L = K(θ)`、`f` を `θ` の最小多項式とすると、
`𝔭 ∤ disc(f)` について「`𝔭` が完全分解 ⟺ `f` が `mod 𝔭` で根を持つ」(`L/K` が Galois のとき)。
「`f` の値の素因子は無限個」は初等的である。
★これで原文 p.116 の (a) は足りる見込みである。 -/
theorem infinite_splitsCompletely (K L : Type) [Field K] [Field L]
    [NumberField K] [NumberField L] [Algebra K L] [IsGalois K L] :
    {𝔭 : HeightOneSpectrum (𝓞 K) | SplitsCompletely K L 𝔭}.Infinite := by
  sorry

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
★★訂正(2026-08-20): [FrdI] が要求するのは**完全分解の場合だけ**で、それは Dedekind ζ の極（mathlib 在庫）で出る。Artin の相互法則は**共役類つきの一般形**にしか要らない。 -/
def chebotarev_splitsCompletely.needs : List ProofObligation :=
  [ .otherPaper "[MilneCFT]" "Chapter VIII, §7, Theorem 7.4(Chebotarev の密度定理)" 274,
    .otherPaper "[MilneCFT]" "Chapter V, §3, Theorem 3.5(Artin の相互法則)" 178,
    .citation "[mathlib]" "NumberField.dedekindZeta_residue(Dedekind ζ の s=1 での留数)"
      (.inMathlib "NumberField.dedekindZeta_residue") 116,
    .citation "[mathlib]" "Nat.infinite_setOf_prime_and_eq_mod(Dirichlet の算術級数定理)"
      (.inMathlib "Nat.infinite_setOf_prime_and_eq_mod") 116,
    .citation "[mathlib]" "Chebotarev / Artin の相互法則 / Dirichlet 密度 / Frobenius 元 / 射類群"
      (.absent "lean/.lake/packages/mathlib/Mathlib/ 全体を Chebotarev|Tchebotarev|ArtinReciprocity|DirichletDensity|FrobeniusElement|RayClassGroup で grep、いずれも 0 件(2026-08-20)") 116,
    .derivation
      "鎖 cheb: 密度の定義 → Frobenius 元 → 射類群 → Artin 写像 → 相互法則 → L(1,χ)≠0 → 射類ごとの等分布 → アーベルの場合 → 巡回への還元 → 2 つの計数((a) 全単射・(b) d/(cf) 対 1)→ 本体" 116,
    .implicitStep
      "★原文は「by Tchebotarev's density theorem [cf., e.g., [Lang2], Chapter VIII, §4, Theorem 10]」と外部文献へ送る" 116 ]

def nonempty_algHom_of_splitsCompletely_subset.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解する素点の包含が体の包含を与える",
    sectionId := "frdi-thm-6-4" }

def nonempty_algHom_of_splitsCompletely_subset.needs : List ProofObligation :=
  [ .otherPaper "[MilneCFT]" "Chapter V, §3, Theorem 3.25(Spl が体を決める)" 178,
    .citation "[ABC3]" "chebotarev_splitsCompletely"
      (.inProject "ABC3" "ABC3.Skeleton.Cheb.chebotarev_splitsCompletely") 116,
    .derivation "Spl(LM/K) = Spl(L/K) ∩ Spl(M/K) と密度の比較で [LM:K] = [M:K]" 116,
    .implicitStep
      "★原文は [NSW] Theorem 12.2.5 へ送る(他論文)" 116 ]

def infinite_splitsCompletely.src : Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解する素点は無限個",
    sectionId := "frdi-thm-6-4" }

/-- ★★**これは密度を要求しない** —— 迂回路である。 -/
def infinite_splitsCompletely.needs : List ProofObligation :=
  [ .otherPaper "[MilneANT]" "Chapter IV(分解群・惰性群・Frobenius の基礎)" 100,
    .citation "[mathlib]" "Ideal.sum_ramification_inertia(Σ e_i f_i = [L:K])"
      (.inMathlib "Ideal.sum_ramification_inertia") 116,
    .derivation
      "L = K(θ)、f を θ の最小多項式として、𝔭 ∤ disc(f) では「完全分解 ⟺ f が mod 𝔭 で根を持つ」。f の値の素因子が無限個であることは初等的(Schur の議論)" 116 ]

end ABC3.Skeleton.Cheb
