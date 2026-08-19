/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24SuppElt

/-!
# [FrdI] Theorem 4.2, (ii) —— 「同じ素点」の**圏論的**判定

原文 (FrdI p.77):
> (ii) There exists a unique isomorphism ΨPrime between the functors

## ★★★★★★測って分かったこと —— `Prime` は「共通の下界」で決まる

`Theorem 4.2, (ii)` は `Prime(Φ(A))` の間の同型 `Ψ_Prime` を作る。
そのためには **`Prime(Φ(A))` が圏論的に復元できる**ことが要る。

★★在庫の `suppElt_disjoint_iff`(`Def24SuppElt.lean:303`)が鍵である:
```
SuppElt ι a ∩ SuppElt ι b = ∅ ↔ ∀ d : M, MLe d a → MLe d b → d = 0
```
★**右辺は `MLe` だけ** —— そして `MLe d x` は Frobenioid では
「pre-step が別の pre-step を経由して分解する」という**圏論的**な条件である
(`prop_4_1_iv` の証明中の `hb2` を見よ)。

★★★したがって primary な 2 元について
「台が一致する ⟺ 0 でない共通の下界を持つ」が言えれば、
**`Prime(Φ(A))` は「primary pre-step の同値類」として圏論的に復元される**。
これが本ファイルの主張である。
-/

namespace ABC3.Found.FrdI

open scoped NNReal

universe w

variable {M : Type w} [AddCommMonoid M] {ι : Prime M → Pf M → ℝ≥0}

/-- ★★★★★★**primary な 2 元の「同素性」は `MLe` だけで書ける**。

★`suppElt_disjoint_iff` の対偶に、`primary ⟹ 台は 1 点` を合わせるだけ。
★★★これが `Ψ_Prime` の存在の数学的な中身である ——
右辺は `Ψ` で移る形（`MLe` は pre-step の分解）なので、
`Prime` の間の対応が**選択なしに**決まる。 -/
theorem suppElt_eq_iff_exists_common (H : IsPerfFactorialWith M ι)
    (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M) {x y : M}
    (hx : IsPrimaryElt x) (hy : IsPrimaryElt y) :
    SuppElt ι x = SuppElt ι y ↔ ∃ d : M, d ≠ 0 ∧ MLe d x ∧ MLe d y := by
  classical
  obtain ⟨p, hp⟩ := suppElt_singleton_of_primary H hperf hdiv hx
  obtain ⟨q, hq⟩ := suppElt_singleton_of_primary H hperf hdiv hy
  constructor
  · intro heq
    -- ★台が一致するなら交わりは空でないので、`suppElt_disjoint_iff` の否定が使える
    have hne : SuppElt ι x ∩ SuppElt ι y ≠ ∅ := by
      rw [heq, Set.inter_self, hq]
      exact fun h => (Set.singleton_ne_empty q) h
    by_contra hc
    refine hne ((suppElt_disjoint_iff H hperf hdiv x y).mpr (fun d hdx hdy => ?_))
    by_contra hd0
    exact hc ⟨d, hd0, hdx, hdy⟩
  · rintro ⟨d, hd0, hdx, hdy⟩
    have hnd : ¬ (SuppElt ι x ∩ SuppElt ι y = ∅) := by
      intro hdisj
      exact hd0 ((suppElt_disjoint_iff H hperf hdiv x y).mp hdisj d hdx hdy)
    -- ★どちらの台も 1 点なので、交わりが空でないことは `p = q` を意味する
    have hpq : p = q := by
      by_contra hne
      refine hnd ?_
      rw [hp, hq]
      exact Set.singleton_inter_eq_empty.mpr (fun h => hne h)
    rw [hp, hq, hpq]

/-- ★★**同じことを「和が primary」でも言える** —— 台の合併が 1 点になること。 -/
theorem isPrimaryElt_add_iff_suppElt_eq (H : IsPerfFactorialWith M ι)
    (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M) {x y : M}
    (hx : IsPrimaryElt x) (hy : IsPrimaryElt y) :
    IsPrimaryElt (x + y) ↔ SuppElt ι x = SuppElt ι y := by
  classical
  obtain ⟨p, hp⟩ := suppElt_singleton_of_primary H hperf hdiv hx
  obtain ⟨q, hq⟩ := suppElt_singleton_of_primary H hperf hdiv hy
  have hadd : SuppElt ι (x + y) = {p} ∪ {q} := by rw [suppElt_add H, hp, hq]
  have hne0 : x + y ≠ 0 := by
    intro h0
    have hz : SuppElt ι (x + y) = ∅ := by rw [h0]; exact suppElt_zero H
    have hmem : p ∈ SuppElt ι (x + y) := by rw [hadd]; exact Or.inl rfl
    rw [hz] at hmem
    exact hmem
  rw [isPrimaryElt_iff_suppElt_singleton H hperf hdiv hne0, hp, hq]
  constructor
  · rintro ⟨r, hr⟩
    rw [hadd] at hr
    have hpr : p = r := by
      have : p ∈ ({r} : Set (Prime M)) := by rw [← hr]; exact Or.inl rfl
      exact this
    have hqr : q = r := by
      have : q ∈ ({r} : Set (Prime M)) := by rw [← hr]; exact Or.inr rfl
      exact this
    rw [hpr, hqr]
  · intro heq
    refine ⟨p, ?_⟩
    rw [hadd]
    have : ({q} : Set (Prime M)) = {p} := heq.symm
    rw [this, Set.union_self]

def suppElt_eq_iff_exists_common.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (ii) — 同素性は共通の下界で決まる",
    sectionId := "frdi-thm-4-2" }

end ABC3.Found.FrdI
