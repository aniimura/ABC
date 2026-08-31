/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Interface.GenEll.EllModuli

/-!
# 界面の測定 —— **`mem_lcyclicExc` は `l` の下界を落としている**（`Check`）

**これは原典の主張ではない**（我々の界面についての事実）ので `.src` を持たない。

## ★★★★★★★★2026-08-31 の測定（第 754）

`Interface/GenEll/EllModuli.lean` の

```
lcyclicExc : Set EllClass
galoisFinite_lcyclicExc : GaloisFinite lcyclicExc
mem_lcyclicExc : ∀ (E : Curve) (l : ℕ), Nat.Prime l → SemiStable E →
    HasLCyclic E l → PrimeToLocalHeights E l → cls E ∈ lcyclicExc
```

は、本ファイルの `lcyclic_classes_finite` が示すとおり

> **`l`-巡回部分群をもち `l` が局所高さと素であるような半安定曲線の類は、
> `l` を動かしても全体で Galois-finite な集合に収まる**

を主張している。★★**これは強すぎる。**

## ★★★なぜ強すぎるか

`Lemma 3.5` が与えるのは

    (†)  (l/14)·deg∞(E) ≤ ht^Falt(E) + 2·log(l) + C′

であり、ここから `ht^Falt` の**上界**を出すには `l` が `ht^Falt` に比べて大きいこと
——すなわち原文の**条件 (a)**（`100·d·(ht^Falt + C·d^ϵ) ≤ l`）——が要る。
★`Found/GaloisRep/Lemma37C.lean` の `htFalt_le_of_condA` がまさにその計算である。

☆`l` の下界が無いと、たとえば `l = 2` を取れば
「2-同種をもち、すべての局所高さが奇数であるような半安定曲線」の類が対象になるが、
それは**有限個ではない**（`ℚ` 上ですでに無限個ある）。

## ★★★★どう直すべきか

`Skeleton/GenEll/Section3.lean` の `lemma_3_7` の証明（部分 (c)）を見ると、
`mem_lcyclicExc` を呼ぶ時点で `hor : condA ∨ condB` を持っている。
★★**`mem_lcyclicExc` の仮説にその情報を渡していない**のが原因である。

したがって直し方は

* `lcyclicExc` を定数 `C`（と `KV`）に依存させ、
* `mem_lcyclicExc` に「条件 (a) または条件 (b)」を仮説として渡す

ことである。☆これは `Theorem 3.8`・`Corollary 4.3/4.4` の証明にも波及するので、
**独立した作業として行う**（本ファイルは測定の記録である）。

## ☆同じ形の測定

* `Check/GenEll/EllModuliDegInfPos.lean`（第 745）——界面は `deg∞ > 0` を強制する
* 本ファイル（第 754）——`mem_lcyclicExc` は `l` の下界を落としている

★どちらも **witness を実際に作ろうとして初めて見えた**ものである。
-/

namespace ABC3.Check.GenEll

open ABC3.Interface.GenEll

/-- ★★★★★★★★**界面が主張していること**（`l` の下界なしの有限性）。

★これは `lcyclicExc` の 2 つの欄を並べただけであり、**強すぎることの明示**である。 -/
theorem lcyclic_classes_finite (D : EllModuliData) :
    ∃ S : Set D.EllClass, D.GaloisFinite S ∧
      ∀ (E : D.Curve) (l : ℕ), Nat.Prime l → D.SemiStable E →
        D.HasLCyclic E l → D.PrimeToLocalHeights E l → D.cls E ∈ S :=
  ⟨D.lcyclicExc, D.galoisFinite_lcyclicExc, D.mem_lcyclicExc⟩

/-- ★★★★★**`l` を 1 つ止めても同じ主張が出る**——これが強すぎることの核心。

☆`l = 2` に取れば「2-同種をもち局所高さがすべて奇数の半安定曲線」の類が
全体で Galois-finite だと言っていることになる。 -/
theorem lcyclic_classes_finite_fixed_l (D : EllModuliData) (l : ℕ) (hl : Nat.Prime l) :
    ∃ S : Set D.EllClass, D.GaloisFinite S ∧
      ∀ E : D.Curve, D.SemiStable E → D.HasLCyclic E l → D.PrimeToLocalHeights E l →
        D.cls E ∈ S :=
  ⟨D.lcyclicExc, D.galoisFinite_lcyclicExc,
    fun E hss hcyc hpr => D.mem_lcyclicExc E l hl hss hcyc hpr⟩

end ABC3.Check.GenEll
