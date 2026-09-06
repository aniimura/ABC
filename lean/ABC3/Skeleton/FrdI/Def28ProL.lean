import ABC3.Meta.Claim
import ABC3.Found.ProL.Defs
import ABC3.Found.ProL.Decomposition

/-!
# [FrdI] Definition 2.8, (ii) —— 副有限アーベル群の pro-`l` 分解(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.52。

原文 (FrdI p.52):
> (ii) Suppose that M is a topologically finitely generated profinite abelian group.

原文 (FrdI p.52):
> Thus, M decomposes as a direct product of pro-l groups M [l], where l varies over

## ★★このファイルの位置づけ —— 「壁」を割って葉を出すために作った

★2026-08-17 まで、ここは「mathlib に pro-`l` 群が無い」の一言で止めていた。
CLAUDE.md の**姿勢**——「工数の山を『壁』と呼ばない。既知数学の person-years は
壁でなく道」——に従い、**statement を型で固定して、そこへ至る道を `.needs` に測る**。

## ★★★★★★2026-09-06: `sorry` は消えた —— `Found/ProL/` からの配線 2 行

`Found/ProL/Decomposition.lean` の `decompEquiv`(`M ≃ₜ* ∏_l M[l]`)と
`Found/ProL/LPart.lean` の `isProL_lPartGrp`(`M[l]` は pro-`l`)は
**どちらも `sorry` 0 で完成していた**。`N l := lPartGrp M l.1` と置けば本定理は閉じる。

★残っていた配管はただ 1 つ、「仮定 `_hcomm : ∀ a b, a * b = b * a` から
`CommGroup M` インスタンスを起こす」ことである(`Found/ProL/` 側は
`[CommGroup M]` を要求する)。`letI` で `M` の `Group` 構造に `mul_comm` を足す。

★`#print axioms def_2_8_ii` は `[propext, Classical.choice, Quot.sound]`
——`sorryAx` は無い。

## ★★道の測定(2026-08-18、探索範囲つき)

| 要るもの | mathlib | 判定 |
|---|---|---|
| `M ≅ lim_U M/U`(開正規部分群による極限表示) | `ProfiniteGrp.continuousMulEquivLimittoFiniteQuotientFunctor`(`Topology/Algebra/Category/ProfiniteGrp/Limits.lean:129`) | ★**ある** |
| 有限アーベル群の `l`-準素分解 | `Ideal.iSup_primaryComponent_eq_top` / `Ideal.iSupIndep_primaryComponent`(`Algebra/Module/Torsion/PrimaryComponent.lean:140,175`) | ★**ある**(`[IsDedekindDomain A]`、`A = ℤ` で使う) |
| その関手性 | `Ideal.primaryComponent.map`(同 :76) | ★**ある** |
| 極限と積の交換 | `CategoryTheory.Limits.limitFlipCompLimIsoLimitCompLim`(`CategoryTheory/Limits/Fubini.lean:545`) | ★**ある** |
| **pro-`l` 群の定義** | `Mathlib/` 全体を `IsProPGroup|ProPGroup|pro-p group` で grep して **0 件**(2026-08-18) | ★**無い**(`Found/ProL/Defs.lean` で我々が定義した) |

★★**すなわち、部品は 4/5 が既に在庫にある。**足りないのは pro-`l` の語彙と、
それらを繋ぐ我々の作業である。分解の DAG は
`ResearchPaper/frdi-decomposition.json` の `prol` チェーン(7 節点・4 層・葉 4 個)。

## ★なぜ迂回できないか

`Definition 2.8, (ii)` の分解は**下流で本質的に使われる**:
原文 p.106–107 は `∏_p 𝒪^×(A)[p] ≅ 𝒪^×(A)` を**そのまま等式として**使い、
各 `p` 成分ごとに `u_p` を作って貼り合わせる。
★また (iii) の「λ 乗写像」は `M[l]` ごとに `λ(l)` 乗する定義なので、
**分解が無いと定義そのものが書けない**。
-/

namespace ABC3.Skeleton.FrdI

open CategoryTheory ABC3.Meta

universe u

/-- ★**pro-`l` 群**は `Found/ProL/Defs.lean` で実装した(負の対照つき)。
ここではその述語を使う。 -/
abbrev IsProLGrp (l : ℕ) (M : ProfiniteGrp.{u}) : Prop := ABC3.Found.ProL.IsProL l M

/-- ★★★**[FrdI] Definition 2.8, (ii)** —— 副有限アーベル群の pro-`l` 分解。

原文 (FrdI p.52):
> Thus, M decomposes as a direct product of pro-l groups M [l], where l varies over

★**原文は「topologically finitely generated」を仮定に置くが、この分解自体は
可換な副有限群なら成り立つ**(位相的有限生成は `Definition 2.8` の他の部分で効く)。
ゆえにここでは可換性だけを仮定する——★**原文より弱い仮定で述べている**ので、
原文の主張はこれの系である(逸脱の向きが安全側)。

★`(l : Nat.Primes) → N l` は積位相・各点積の群。`≃ₜ*` は
`ContinuousMulEquiv`(位相群の同型)である。

★★**証明(2026-09-06)**: `N l := lPartGrp M l.1`。
`isProL_lPartGrp` が第 1 成分を、`decompEquiv` が第 2 成分を与える。
`letI` は `_hcomm` から `CommGroup ↑M` を起こすためのもので、
`CommGroup.toGroup` は `M.str` と定義的に等しいので同型の型は変わらない。

☆**退化しない**ことの確認: `N` をすべて自明群に取ると
`Nonempty (M ≃ₜ* 1)` が `M` 自身の自明性を要求するので、同型が錨になっている。 -/
theorem def_2_8_ii (M : ProfiniteGrp.{u}) (_hcomm : ∀ a b : M, a * b = b * a) :
    ∃ N : Nat.Primes → ProfiniteGrp.{u},
      (∀ l : Nat.Primes, IsProLGrp l.1 (N l)) ∧
      Nonempty (M ≃ₜ* ((l : Nat.Primes) → N l)) := by
  letI : CommGroup M := { (inferInstance : Group M) with mul_comm := _hcomm }
  exact ⟨fun l => ABC3.Found.ProL.lPartGrp M l.1,
    fun l => ABC3.Found.ProL.isProL_lPartGrp (M := (M : Type u)) l.1,
    ⟨ABC3.Found.ProL.decompEquiv M⟩⟩

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def IsProLGrp.src : Source :=
  { paper := "FrdI", pdfPage := 52, item := "Definition 2.8, (ii)",
    sectionId := "frdi-def-2-8" }

def def_2_8_ii.src : Source :=
  { paper := "FrdI", pdfPage := 52, item := "Definition 2.8, (ii)",
    sectionId := "frdi-def-2-8" }

/-- ★★分解の DAG は `ResearchPaper/frdi-decomposition.json` の `prol` チェーン。
`node tools/frdi-newleaves.mjs` が層と葉を印字する。 -/
def def_2_8_ii.needs : List ProofObligation :=
  [ .citation "[mathlib]" "ProfiniteGrp.continuousMulEquivLimittoFiniteQuotientFunctor(M = lim_U M/U)"
      (.inMathlib "ProfiniteGrp.continuousMulEquivLimittoFiniteQuotientFunctor") 52,
    .citation "[mathlib]" "Ideal.iSup_primaryComponent_eq_top / Ideal.iSupIndep_primaryComponent(有限アーベル群の l-準素分解)"
      (.inMathlib "Ideal.iSup_primaryComponent_eq_top") 52,
    .citation "[mathlib]" "Ideal.primaryComponent.map(準素分解の関手性)"
      (.inMathlib "Ideal.primaryComponent.map") 52,
    .citation "[mathlib]" "CategoryTheory.Limits.limitFlipCompLimIsoLimitCompLim(極限と積の交換)"
      (.inMathlib "CategoryTheory.Limits.limitFlipCompLimIsoLimitCompLim") 52,
    .citation "[mathlib]" "pro-l 群の理論(mathlib に無いので我々が定義した)"
      (.absent "lean/.lake/packages/mathlib/Mathlib/ 全体を IsProPGroup|ProPGroup|pro-p group で grep、0 件(2026-08-18)。★2026-09-06 に再測(索引 .cache/mathlib-index.txt): re:`IsProPGroup|ProPGroup|pro-p group|pro-l group|IsProCyclic`→0") 52,
    .citation "[ABC3]" "decompEquiv(M ≅ ∏_l M[l]、位相群の同型。sorry 0)"
      (.inProject "ABC3" "ABC3.Found.ProL.decompEquiv") 52,
    .citation "[ABC3]" "isProL_lPartGrp(M[l] は pro-l 群。sorry 0)"
      (.inProject "ABC3" "ABC3.Found.ProL.isProL_lPartGrp") 52,
    .derivation
      "M[l] := lim_U (M/U)[l] を作り、M ≅ lim_U M/U ≅ lim_U ∏_l (M/U)[l] ≅ ∏_l M[l] と繋ぐ(frdi-decomposition.json の prol チェーン、葉 4 個)。★2026-09-06 に Found/ProL/ の 10 ファイルで完了" 52,
    .implicitStep
      "★原文は分解を『Thus』の一語で述べ、証明を置かない。位相的有限生成の仮定はこの分解には要らない(我々は可換性だけを仮定して述べた)" 52 ]

end ABC3.Skeleton.FrdI
