/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.ProL.Decomposition
import ABC3.Found.ProL.ZetaPow

/-!
# DivByP —— `[FrdI] Definition 2.8` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.ProL

universe u
variable {N : Type u} [CommGroup N] [TopologicalSpace N] [IsTopologicalGroup N]
  [CompactSpace N] [TotallyDisconnectedSpace N]

/-- ★★★**pro-`p` 群では `⋂_k N^{p^k} = 1`**。

★各有限商 `N/U` は有限 `p` 群なので `Nat.card (N/U) = p^e`。
`x^{p^e} = y` を落とすと `(mk x)^{p^e} = 1`(Lagrange)なので `mk y = 1`。
★★**位相的有限生成は要らない**。 -/
theorem eq_one_of_forall_pow_pow {p : ℕ} (hp : p.Prime)
    (hpro : ∀ U : OpenNormalSubgroup N, IsPGroup p (N ⧸ U.toSubgroup))
    {y : N} (hy : ∀ k : ℕ, ∃ x : N, x ^ (p ^ k) = y) : y = 1 := by
  haveI : Fact p.Prime := ⟨hp⟩
  refine eq_of_forall_quotient_eq (asProfiniteGrp N) (fun U => ?_)
  show (QuotientGroup.mk y : N ⧸ U.toSubgroup) = QuotientGroup.mk 1
  obtain ⟨e, he⟩ := (IsPGroup.iff_card (p := p) (G := N ⧸ U.toSubgroup)).mp (hpro U)
  obtain ⟨x, hx⟩ := hy e
  have h1 : (QuotientGroup.mk x : N ⧸ U.toSubgroup) ^ (p ^ e) = QuotientGroup.mk y := by
    rw [← QuotientGroup.mk_pow, hx]
  have h2 : (QuotientGroup.mk x : N ⧸ U.toSubgroup) ^ (Nat.card (N ⧸ U.toSubgroup)) = 1 :=
    pow_card_eq_one'
  rw [he] at h2
  rw [← h1, h2]
  exact (QuotientGroup.mk_one _).symm

/-! ### ★出典の紐付け -/

/-- ★locator —— `Definition 2.8, (ii)` の pro-`l` 分解の帰結。 -/
def eq_one_of_forall_pow_pow.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 52, item := "Definition 2.8, (ii) — pro-p 群では ⋂_k N^{p^k} = 1",
    sectionId := "frdi-def-2-8" }

end ABC3.Found.ProL
