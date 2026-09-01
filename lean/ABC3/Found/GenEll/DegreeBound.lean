/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.FieldTheory.PrimitiveElement
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import ABC3.Meta.Claim

/-!
# 第 1214 ブロック —— **共役の個数で次数を上から抑える**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`[M:L] < l` へ向かう段

第 1209 で「残る義務は `[M:L] < l` の 1 つ」になった。
☆原文の道では `M = L(H)` の Galois 群が `𝔽_l^×` に埋まるので
`[M:L] ∣ l−1 < l` である。

★その議論の**数え上げの部分**を、`Galois` を使わない形で取る:

    `M →ₐ[L] L̄` が有限集合 `s` に単射で入る  ⟹  `[M:L] ≤ #s`

☆標数 0（分離的）なら `#(M →ₐ[L] L̄) = [M:L]`（mathlib の `AlgHom.card`）。

★★★あとは「`φ : M →ₐ[L] L̄` は `Q ↦ c·Q` の `c ∈ 𝔽_l^×` で決まる」を示せば
`[M:L] ≤ l−1 < l` が出る——第 1205 と第 1213 がその材料である。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

/-- ★★★★★★★★★★★★
**共役の個数で次数を上から抑える**——★**無条件**（第 1214）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆標数 0（分離的）なら `#(M →ₐ[L] L̄) = [M:L]` なので、
`M →ₐ[L] L̄` が有限集合 `s` に単射で入れば `[M:L] ≤ #s` である。

★★★これが `[M:L] < l`（第 1209 が残した唯一の義務）へ向かう
**数え上げの段**である。 -/
theorem finrank_le_of_algHom_mem_finset
    {L M Lbar : Type*} [Field L] [Field M] [Field Lbar]
    [Algebra L M] [Algebra L Lbar] [IsAlgClosed Lbar]
    [FiniteDimensional L M] [Algebra.IsSeparable L M]
    {α : Type*} [DecidableEq α] (f : (M →ₐ[L] Lbar) → α)
    (hinj : Function.Injective f) (s : Finset α) (hs : ∀ φ, f φ ∈ s) :
    Module.finrank L M ≤ s.card := by
  classical
  have hcard : Fintype.card (M →ₐ[L] Lbar) = Module.finrank L M := AlgHom.card L M Lbar
  rw [← hcard, ← Finset.card_univ]
  exact Finset.card_le_card_of_injOn f (fun a _ => hs a) (fun a _ b _ h => hinj h)

/-! ## ★出典の紐付け(`.src`) -/

def finrank_le_of_algHom_mem_finset.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(共役の個数で次数を上から抑える。★無条件)",
    sectionId := "genell-lemma-3-5" }

def finrank_le_of_algHom_mem_finset.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1214）**——第 1209 が残した唯一の義務 `[M:L] < l` へ" ++
       "向かう**数え上げの段**である。☆標数 0（分離的）なら" ++
       "`#(M →ₐ[L] L̄) = [M:L]`（mathlib の `AlgHom.card`）。" ++
       "★あとは「`φ : M →ₐ[L] L̄` は `Q ↦ c·Q` の `c ∈ 𝔽_l^×` で決まる」を" ++
       "示せば `[M:L] ≤ l−1 < l` が出る——第 1205 と第 1213 がその材料である。") 3 ]

end ABC3.Found.GenEll
