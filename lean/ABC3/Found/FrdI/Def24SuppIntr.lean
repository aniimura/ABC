/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24Supp

/-!
# [FrdI] `Definition 2.4, (i), (d)` —— `Supp` は `ι` を含まない形に書ける

原文 (FrdI p.48):
> Prime(M ), as a ranges over the elements of M pf; if a, b

## ★★なぜ要るか(2026-08-19)

`Definition 4.5, (ii)` の `IsStrictlyRational` は `SuppElt (ι Y) a` を使う。
`ι` は **各 `Y : D` ごとに与えられた族**であり、
`Remark 4.5.1` で `B` と `B^istr` を比べるとき両者の底は
**同型だが等しくない**。★したがって `Supp` が底の同型で移ることが要る。

## ★★★測ったこと —— `Supp` は `ι` を消せる

`factorMap ι a p = boundSup (ι p) (pCarrierPf M p) a = sSup (ι p '' Bound _ _ a)` である。
★`ι p` は `pCarrierPf M p` の上で**単射**(`embedInj`)で `ι p 0 = 0`(`iota_zero`)なので、

```
p ∈ Supp (factorMap ι a)  ⟺  ∃ x ∈ Bound (Pf M) (pCarrierPf M p) a, x ≠ 0
```

★★右辺に `ι` は現れない。したがって `Supp` は**単系だけで決まる**。

★`Def24.lean` は「2 つの埋め込みは `ℝ≥0` の正の実数倍だけ違う」と記録しているが、
**その一意性を経由するより本補題の方が短い**(実測)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe w

variable {M : Type w} [AddCommMonoid M] {ι : Prime M → Pf M → ℝ≥0}

/-! ## ★1. `ι` を含まない特徴づけ -/

/-- ★★★★★**`Supp` は `ι` を含まない形に書ける**。

原文 (FrdI p.48):
> Prime(M ), as a ranges over the elements of M pf; if a, b

★★`ι p` が `pCarrierPf M p` 上で単射で `0 ↦ 0` であることだけを使う。 -/
theorem mem_supp_factorMap_iff (H : IsPerfFactorialWith M ι) (a : Pf M) (p : Prime M) :
    p ∈ Supp (factorMap ι a) ↔ ∃ x ∈ Bound (Pf M) (pCarrierPf M p) a, x ≠ 0 := by
  constructor
  · intro hp
    by_contra hc
    push_neg at hc
    refine hp ?_
    show boundSup (ι p) (pCarrierPf M p) a = 0
    refine le_antisymm ?_ zero_le'
    refine boundSup_le (ι p) (zero_mem_pCarrierPf p) a ?_
    intro x hx
    rw [hc x hx, iota_zero H p]
  · rintro ⟨x, hx, hx0⟩ hcon
    have hcon' : boundSup (ι p) (pCarrierPf M p) a = 0 := hcon
    have h1 : ι p x ≤ boundSup (ι p) (pCarrierPf M p) a :=
      le_boundSup (ι p) _ a (H.bounded a p) hx
    rw [hcon', le_zero_iff] at h1
    exact hx0 (H.embedInj p hx.1 (zero_mem_pCarrierPf p) (h1.trans (iota_zero H p).symm))

def mem_supp_factorMap_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 48,
    item := "Definition 2.4, (i), (d) — Supp は ι を含まない形に書ける",
    sectionId := "frdi-def-2-4" }

end ABC3.Found.FrdI
