/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Topology.Algebra.Category.ProfiniteGrp.Limits

/-!
# 副有限群の極限表示(使える形) —— チェーン `prol` の葉 `M-limit`

★`ResearchPaper/frdi-decomposition.json` の `prol` チェーンの葉。
最終目標は [FrdI] `Definition 2.8, (ii)`。

## ★測定(2026-08-18)

★mathlib に `ProfiniteGrp.continuousMulEquivLimittoFiniteQuotientFunctor`
(`Topology/Algebra/Category/ProfiniteGrp/Limits.lean:129`)は**ある**。
`toLimit_injective` / `toLimit_surjective` も。

★★**しかし圏論的な極限のままでは使えない。**
下流(準素分解の貼り合わせ)で実際に要るのは次の 2 つの**具体的な**言明である:

| | 宣言 |
|---|---|
| **分離** —— すべての有限商で等しければ等しい | `eq_of_forall_quotient_eq` |
| **持ち上げ** —— 整合な族は `M` の元から来る | `exists_mk_eq_of_compatible` |
| 持ち上げの一意性 | `existsUnique_mk_eq_of_compatible` |

★このファイルはその 2 つを取り出すだけである。**これが「葉」の中身**であり、
圏論的な極限を毎回展開しなくて済むようになる。
-/

namespace ABC3.Found.ProL

open CategoryTheory

universe u

variable (P : ProfiniteGrp.{u})

/-- ★★**分離** —— すべての開正規部分群による商で等しければ、`M` の中で等しい。

★これは「開正規部分群の共通部分が自明」の言い換えである
(`ProfiniteGrp.toLimit_injective`)。 -/
theorem eq_of_forall_quotient_eq {x y : P}
    (h : ∀ U : OpenNormalSubgroup P,
      (QuotientGroup.mk x : P ⧸ U.toSubgroup) = QuotientGroup.mk y) : x = y :=
  P.toLimit_injective (Subtype.ext (funext fun U => h U))

/-- ★★**持ち上げ** —— 遷移射と整合する族は、`M` のある元の像である。

★`ProfiniteGrp.toLimit_surjective`(稠密性 + 像の閉性)の具体形。 -/
theorem exists_mk_eq_of_compatible (f : ∀ U : OpenNormalSubgroup P, P ⧸ U.toSubgroup)
    (hf : ∀ (U V : OpenNormalSubgroup P) (h : U.toSubgroup ≤ V.toSubgroup),
      QuotientGroup.map U.toSubgroup V.toSubgroup (MonoidHom.id P) h (f U) = f V) :
    ∃ x : P, ∀ U, (QuotientGroup.mk x : P ⧸ U.toSubgroup) = f U := by
  obtain ⟨x, hx⟩ := P.toLimit_surjective ⟨f, fun U V π => hf U V (leOfHom π)⟩
  exact ⟨x, fun U => congrFun (congrArg Subtype.val hx) U⟩

/-- ★持ち上げは一意。 -/
theorem existsUnique_mk_eq_of_compatible (f : ∀ U : OpenNormalSubgroup P, P ⧸ U.toSubgroup)
    (hf : ∀ (U V : OpenNormalSubgroup P) (h : U.toSubgroup ≤ V.toSubgroup),
      QuotientGroup.map U.toSubgroup V.toSubgroup (MonoidHom.id P) h (f U) = f V) :
    ∃! x : P, ∀ U, (QuotientGroup.mk x : P ⧸ U.toSubgroup) = f U := by
  obtain ⟨x, hx⟩ := exists_mk_eq_of_compatible P f hf
  refine ⟨x, hx, fun y hy => eq_of_forall_quotient_eq P (fun U => ?_)⟩
  rw [hy U, hx U]

/-- 各有限商への写像は全射(自明だが、下流で `obtain ⟨y, rfl⟩` として使う)。 -/
theorem mk_surjective (U : OpenNormalSubgroup P) :
    Function.Surjective (QuotientGroup.mk : P → P ⧸ U.toSubgroup) :=
  QuotientGroup.mk_surjective

end ABC3.Found.ProL
