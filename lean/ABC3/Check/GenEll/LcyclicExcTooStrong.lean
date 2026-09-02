/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Interface.GenEll.EllModuli

/-!
# 界面の測定と訂正 —— **`mem_lcyclicExc` は `l` の下界を落としていた**（`Check`）

**これは原典の主張ではない**（我々の界面についての事実）ので `.src` を持たない。

## ★★★★★★★★2026-08-31 の測定（第 754）と訂正（第 755）

`Interface/GenEll/EllModuli.lean` の `lcyclicExc` は、以前

```
lcyclicExc : Set EllClass
galoisFinite_lcyclicExc : GaloisFinite lcyclicExc
mem_lcyclicExc : ∀ (E : Curve) (l : ℕ), Nat.Prime l → SemiStable E →
    HasLCyclic E l → PrimeToLocalHeights E l → cls E ∈ lcyclicExc
```

であった。★★この形は

> **`l`-巡回部分群をもち `l` が局所高さと素であるような半安定曲線の類は、
> `l` を動かしても——`l` を 1 つ止めても——全体で Galois-finite な集合に収まる**

を主張しており、**強すぎて witness が作れない**。

## ★★★なぜ強すぎたか

`Lemma 3.5` が与えるのは

    (†)  (l/14)·deg∞(E) ≤ ht^Falt(E) + 2·log(l) + C′

であり、ここから `ht^Falt` の**上界**を出すには `l` が `ht^Falt` に比べて大きいこと
——すなわち原文の**条件 (a)**（`100·d·(ht^Falt + C·d^ϵ) ≤ l`）——が要る。
★`Found/GaloisRep/Lemma37C.lean` の `htFalt_le_of_condA` がまさにその計算である。

☆`l` の下界が無いと、たとえば `l = 2` を取れば
「2-同種をもち、すべての局所高さが奇数であるような半安定曲線」の類が対象になるが、
それは**有限個ではない**（`ℚ` 上ですでに無限個ある）。

## ★★★★訂正（第 755）

`Skeleton/GenEll/Section3.lean` の `lemma_3_7` の証明（部分 (c)）は、
`mem_lcyclicExc` を呼ぶ時点で `hor : condA ∨ condB` を**持っている**。
★そこで欄を

```
lcyclicExc : ℝ → ℝ → Set EllClass → Set EllClass
galoisFinite_lcyclicExc : ∀ C eps KV, CompactlyBounded KV → GaloisFinite (lcyclicExc C eps KV)
mem_lcyclicExc : ∀ C eps KV E l, Prime l → SemiStable E → HasLCyclic E l →
    PrimeToLocalHeights E l →
    ((100·d·(ht^Falt + C·d^eps) ≤ l) ∨ cls E ∈ KV) → cls E ∈ lcyclicExc C eps KV
```

に直した。★★**`lemma_3_7`・`theorem_3_8`・`Corollary 4.3/4.4` の statement は
1 文字も変わっていない**——変わったのは `Exc` の作り方だけである。

## ☆同じ形の測定

* `Check/GenEll/EllModuliDegInfPos.lean`（第 745）——界面は `deg∞ > 0` を強制する
* 本ファイル（第 754-755）——`mem_lcyclicExc` は `l` の下界を落としていた

★どちらも **witness を実際に作ろうとして初めて見えた**ものである。
-/

namespace ABC3.Check.GenEll

open ABC3.Interface.GenEll

/-- ★★★★★★**訂正後に界面が主張していること**。

★`l` の下界（条件 (a)）または `cls E ∈ K_V`（条件 (b)）が**仮説に入っている**。
☆これが無いと `l = 2` で反例が出る（上の docstring）。 -/
theorem lcyclic_classes_finite (D : EllModuliData) (eps : ℝ) (heps : 0 < eps)
    (KV : Set D.EllClass) (hKV : D.CompactlyBounded KV) :
    ∃ C₀ : ℝ, ∀ C : ℝ, C₀ ≤ C →
    ∃ S : Set D.EllClass, D.GaloisFinite S ∧
      ∀ (E : D.Curve) (l : ℕ), Nat.Prime l → D.SemiStable E →
        D.HasLCyclic E l → D.PrimeToLocalHeights E l →
        (((100 * (D.degOfDefinition E : ℝ)
              * (D.faltingsHeight (D.cls E) + C * (D.degOfDefinition E : ℝ) ^ eps) ≤ (l : ℝ))
            ∧ D.HasMultRed E)
          ∨ D.cls E ∈ KV) →
        D.cls E ∈ S :=
  let ⟨C₀, hC₀⟩ := D.galoisFinite_lcyclicExc eps heps KV hKV
  ⟨C₀, fun C hC => ⟨D.lcyclicExc C eps KV, hC₀ C hC, D.mem_lcyclicExc C eps KV⟩⟩

end ABC3.Check.GenEll
