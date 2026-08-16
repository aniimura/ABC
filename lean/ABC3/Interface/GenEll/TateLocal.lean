import ABC3.Meta.Claim
import Mathlib.Data.Rat.Defs
import Mathlib.Data.Real.Basic

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

★本 Interface は**公理を 1 つも持たない**(データと述語だけ)。
`Lemma 3.2` の主張(i)(ii)、`Remark 3.3.1` の well-defined 性は、
**すべて `Skeleton/GenEll/Section3.lean` 側の `sorry` として残る**。
これが `tools/check.mjs` 冒頭 B5 の穴(条件を posit して `sorry` を消す)を避ける形である。
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

/-- ★Track B は何を作らねばならないか。 -/
def TateLocalData.waiting : WaitingFor :=
  { what := "Tate 曲線と Tate 母数 q_E、半安定還元(特殊ファイバー ≅ 𝔾_m)、mod l Tate 加群 M_l(E) = Hom(ℤ/lℤ, E(K̄)) とその G_K 作用、および完全列 0 → 𝔽_l(1) → M_l(E) → 𝔽_l → 0(Tate twist を含む)"
    trackB := "Found/GenEll — ★mathlib は楕円曲線の群構造(`Affine.Point` の `AddCommGroup`)と分点多項式までは持つが、`E[n] ≅ (ℤ/n)²`・Galois 作用・Tate 曲線はいずれも無い(2026-08-16 実測)。★[FC](Faltings–Chai, *Degenerations of Abelian Varieties*)Chapter III, Corollary 7.3 が原文の典拠であり、これも mathlib に無い" }

/-! ## ★出典の紐付け(`.src`) -/

def TateLocalData.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Lemma 3.2",
    sectionId := "genell-lemma-3-2" }

end ABC3.Interface.GenEll
