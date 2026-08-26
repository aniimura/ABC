import ABC3.Meta.Claim
import Mathlib.Data.Rat.Defs
import Mathlib.Data.Real.Basic
import Mathlib.Data.Nat.Prime.Basic

/-!
# [GenEll] §3 局所理論 —— Tate 曲線と局所高さの `Interface`

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.15–p.16。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.15):
> Lemma 3.2. (Local Rank One Subgroups of l-Torsion)

## ★★なぜ `Interface` なのか —— Tate 曲線が mathlib に無い

`Lemma 3.2` / `Definition 3.3` / `Remark 3.3.1` の statement に現れる語:

| statement の語 | mathlib(2026-08-16 実測) |
|---|---|
| **Tate 母数** `q_E ∈ 𝔪_K` | ★**無い**(`EllipticCurve/` 配下の全宣言名を確認) |
| **半安定還元**(特殊ファイバー ≅ `𝔾_m`) | ★**無い**(`Reduction.lean` はあるが Tate 曲線の理論は無い) |
| "mod l" Tate 加群 `M_l(E) = Hom(ℤ/lℤ, E(K̄))` | ★**無い**(`E[n] ≅ (ℤ/n)²` も Galois 作用も無い) |
| 完全列 `0 → 𝔽_l(1) → M_l(E) → 𝔽_l → 0` | ★**無い**(Tate twist が無い) |

★したがって §3 の局所理論は**まるごとこの Interface が受ける**。

## ★★ただし「何を作れば終わりか」は完全に決まる

★★★★★★**2026-08-26 の訂正(欠陥 #7 を部分的に塞ぐ)**。

以前は本 Interface は**公理を 1 つも持たなかった**(データと述語だけ)。
★そのため `Skeleton/GenEll/Section3.lean` の `lemma_3_2` ・ `potLocalHeight_indep` は
**退化した `D` で破れてしまう**——すなわち「まだ証明していない」のではなく
**「この形では偽」**だった(`Check/GenEll/Section3NotProvable.lean` で実測)。

★★そこで `vq_baseChange`(局所高さは分岐指数倍)を足した。
★★★**これは posit ではない**——`Found/GenEll/LocalHeightRamified.lean`(第 359)が
mathlib の `valuation_liesOver` から**定理として**導いている。
★★★★**2026-08-26 の第 2 の訂正**: `Lemma 3.2` の欄を足すとき
`vq_quotMu` を `∀ l : ℕ` で書くと**充足不能**になる
——`l = 0` で `v_K(q_{E/μ_0}) = 0` と `vq_pos` が衝突する。
★原文の `l` は**素数**なので `Nat.Prime l` を仮定に加えた。
★★充足可能であることは `Check/GenEll/Section3NotProvable.lean` の
`tateLocalSatisfiable` で確かめている。
-/

namespace ABC3.Interface.GenEll

open ABC3.Meta

/-- **[GenEll] §3 の局所理論**を受ける `Interface`。

★原文 p.15 の設定「`E → Spec(O_K)` a one-dimensional semi-abelian scheme over `O_K`
such that the generic fiber `E_K` of `E` is proper, while the special fiber `E_k` of `E`
is isomorphic to `(𝔾_m)_k`」を、**点と数値の水準**で受ける。 -/
structure TateLocalData where
  /-- 非アルキメデス局所体 `K`(の添字)。 -/
  LocalField : Type
  /-- `K` 上の、上の条件を満たす 1 次元半アーベル scheme `E`。 -/
  Curve : LocalField → Type
  /-- 原文 `v_K(q_E)` —— **Tate 母数の付値**。`Definition 3.3` の**局所高さ**そのもの。 -/
  vq : (K : LocalField) → Curve K → ℕ
  /-- 原文「the **positive** integer `v_K(q_E) ∈ ℤ_{>0}`」——★正であることは
  `Definition 3.3` が明記しているので、ここで持つ。 -/
  vq_pos : ∀ (K : LocalField) (E : Curve K), 0 < vq K E
  /-- 原文 `deg_∞(E) ≝ log(#(O_K/(q_E)))`。 -/
  degInf : (K : LocalField) → Curve K → ℝ
  /-- 原文「a one-dimensional `𝔽_l`-subspace which is **stabilized** by `G_K`」——
  `M_l(E)` の `G_K` 安定な 1 次元部分空間 `N`。 -/
  StableLine : (K : LocalField) → Curve K → ℕ → Type
  /-- 原文「`N` is equal to the submodule `𝔽_l(1) ⊆ M_l(E)`」。 -/
  IsCyclotomic : {K : LocalField} → {E : Curve K} → {l : ℕ} → StableLine K E l → Prop
  /-- 原文 `E′ ≝ E/μ_l`(`Lemma 3.2, (ii)`)。 -/
  quotMu : (K : LocalField) → (E : Curve K) → ℕ → Curve K
  /-- 原文「a **finite extension** `L` of `K` over which `E_K` has **multiplicative
  reduction**」(`Remark 3.3.1`)。 -/
  MultExt : (K : LocalField) → Curve K → Type
  /-- その拡大体 `L` 自身。 -/
  extField : {K : LocalField} → {E : Curve K} → MultExt K E → LocalField
  /-- 原文 `E_K ×_K L`。 -/
  baseChange : {K : LocalField} → {E : Curve K} → (L : MultExt K E) → Curve (extField L)
  /-- 原文「the **ramification index** of `L/K`」。 -/
  ramIdx : {K : LocalField} → {E : Curve K} → MultExt K E → ℕ
  /-- 分岐指数は正。 -/
  ramIdx_pos : ∀ {K : LocalField} {E : Curve K} (L : MultExt K E), 0 < ramIdx L
  /-- ★★★★★★**局所高さは分岐指数倍になる**——`Remark 3.3.1` の実質。

  原文 (GenEll p.16):
  > that this definition is independent of the choice of L].

  ★原文は「one verifies **immediately**」で畳んでいるが、畳んでいるのは
  **`v_L(q_E) = e(L/K)·v_K(q_E)`** である。

  ★★★★**これは posit ではなく定理である**——
  `Found/GenEll/LocalHeightRamified.lean` の `ordAt_liesOver` が
  mathlib の `HeightOneSpectrum.valuation_liesOver` から導いている(2026-08-26、第 359)。
  ★★ここでは「意図した対象が満たす仕様」として欄に出している。

  ★★★★★★**退化封じとしても効く**——`v_K(q_E)` を定数にしたまま
  分岐指数だけを動かす witness はこれで落ちる。 -/
  vq_baseChange : ∀ {K : LocalField} {E : Curve K} (L L' : MultExt K E),
    vq (extField L) (baseChange L) * ramIdx L'
      = vq (extField L') (baseChange L') * ramIdx L
  /-- 剰余体の位数の対数 `log #(O_K/𝔪_K)`。 -/
  logResidueCard : LocalField → ℝ
  /-- ★剰余体は 2 元以上なので対数は正。 -/
  logResidueCard_pos : ∀ K : LocalField, 0 < logResidueCard K
  /-- ★★★★**`deg_∞(E) = v_K(q_E)·log #(O_K/𝔪_K)`**——原文 `Lemma 3.2, (ii)` の定義。

  原文 (GenEll p.15):
  > Lemma 3.2. (Local Rank One Subgroups of l-Torsion)

  ★これがあるので `deg_∞` を定数にした退化 witness は作れない。 -/
  degInf_eq : ∀ (K : LocalField) (E : Curve K),
    degInf K E = (vq K E : ℝ) * logResidueCard K
  /-- ★★★★★**`q_{E/μ_l} = q_E^l`**——原文 `Lemma 3.2, (ii)` の前半。

  原文 (GenEll p.15):
  > Lemma 3.2. (Local Rank One Subgroups of l-Torsion)

  ★原文の典拠は **[FC] Chapter III, Corollary 7.3** である。
  ★★付値を取ると `v_K(q_{E′}) = l·v_K(q_E)` になる。 -/
  vq_quotMu : ∀ (K : LocalField) (E : Curve K) (l : ℕ), Nat.Prime l →
    vq K (quotMu K E l) = l * vq K E
  /-- ★★★★★★**`Lemma 3.2, (i)`**——`G_K` 安定な 1 次元部分空間は
  `μ_l` であるか、さもなくば局所高さが `l` で割り切れる。

  原文 (GenEll p.15):
  > Lemma 3.2. (Local Rank One Subgroups of l-Torsion)

  ★★原文はこれを『we have the following **well-known** result』として提示し、
  証明を与えていない。典拠は **[FC] Chapter III, Corollary 7.3** の完全列
  `0 → 𝔽_l(1) → M_l(E) → 𝔽_l → 0` である。
  ★★★**完全列からこの主張への段**(拡大類が `q_E` の `l` 乗根の抽出で与えられること)は
  原文に書かれていない——mathlib に `M_l(E)`･Tate twist が無いので、
  ここでは**結論を仕様として受ける**。 -/
  stableLine_dvd_or_cyclotomic : ∀ (K : LocalField) (E : Curve K) (l : ℕ), Nat.Prime l →
    ∀ (N : StableLine K E l), l ∣ vq K E ∨ IsCyclotomic N

/-- ★Track B は何を作らねばならないか。 -/
def TateLocalData.waiting : WaitingFor :=
  { what := "Tate 曲線と Tate 母数 q_E、半安定還元(特殊ファイバー ≅ 𝔾_m)、mod l Tate 加群 M_l(E) = Hom(ℤ/lℤ, E(K̄)) とその G_K 作用、および完全列 0 → 𝔽_l(1) → M_l(E) → 𝔽_l → 0(Tate twist を含む)"
    trackB := "Found/GenEll — ★mathlib は楕円曲線の群構造(`Affine.Point` の `AddCommGroup`)と分点多項式までは持つが、`E[n] ≅ (ℤ/n)²`・Galois 作用・Tate 曲線はいずれも無い(2026-08-16 実測)。★[FC](Faltings–Chai, *Degenerations of Abelian Varieties*)Chapter III, Corollary 7.3 が原文の典拠であり、これも mathlib に無い" }

/-! ## ★出典の紐付け(`.src`) -/

def TateLocalData.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Lemma 3.2",
    sectionId := "genell-lemma-3-2" }

end ABC3.Interface.GenEll
