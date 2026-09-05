import Mathlib.GroupTheory.Index
import Mathlib.GroupTheory.QuotientGroup.Basic
import Mathlib.Tactic.Group
import Mathlib.Tactic.Ring

/-!
# 有限指数の部分群では Herbrand 商が変わらない

論文にも我々のモデルにも固有でない、**一般の**群論(`Found/Teichmuller.lean`・
`Found/HenselianSplits.lean` と同じ位置づけ。mathlib へ出せる形で書く)。
mathlib には Herbrand 商が無く(`grep -i herbrand` が 0 件)、ABC3 側にも無かった。

## 何を言っているか

可換群 `A` と、その**有限指数**の部分群 `B` を取る。`n` 乗写像
`(·)^n : A → A` に対して、核 `A[n] := {a | a^n = 1}` と余核 `A/A^n` の
「比」——古典的に **Herbrand 商** `h_n(A) := #A[n] / [A : A^n]` と呼ばれる量
——は `A` と `B` で一致する:

```
[A : A^n] · #B[n]  =  [B : B^n] · #A[n]                          (割り算を使わない形)
```

`#A[n]` や `[A : A^n]` が無限でも(`Nat.card = 0`・`Subgroup.index = 0` の規約で)
両辺が釣り合うので、**有限性の仮定は `B` の指数だけ**でよい。

## 証明——連結準同型を作らずに済ませる

教科書は蛇の補題で 6 項完全列

```
1 → B[n] → A[n] → C[n] → B/B^n → A/A^n → C/C^n → 1     (C := A/B, 有限)
```

を作り、交代積を取る。連結準同型 `C[n] → B/B^n` は代表元の取り方に依らないことを
確かめる必要があり、Lean では選択公理の後始末が面倒になる。ここでは**それを避ける**。

`f := (·)^n : A →* A`(`A` が可換だからこれは準同型)と置き、`B^n := f(B)` と書く。
鍵は**一本の合成準同型**

```
ψ : A →* A ⧸ B^n,     ψ(a) = a^n · B^n
```

だけである。`ψ` の核は `f⁻¹(B^n) = B · A[n]`(`comap_powMap_eq_sup_ker`
——`a^n = b^n` なら `a = b · (b⁻¹a)` で `(b⁻¹a)^n = 1`)、像は `A^n` の
`A ⧸ B^n` での像。`Subgroup.index_ker` を `ψ` と、`A^n ↪ A → A ⧸ B^n` の
二つに使うと

```
[A : B · A[n]]  =  [A^n : B^n]                                   (index_sup_ker)
```

が出る。あとは指数の塔をつなぐだけ:

```
[A^n : B^n] · [A : A^n] = [A : B^n] = [B : B^n] · [A : B]        (relIndex_mul_index 二本)
[A : B]  =  [A : B·A[n]] · [B·A[n] : B]  =  [A^n : B^n] · [A[n] : B[n]]
                                                (第二同型定理 relIndex_sup_left)
```

を代入して `[A^n : B^n]`(`[A : B]` の約数だから `B` の有限指数性より `≠ 0`)を
約せば `[A : A^n] = [B : B^n] · [A[n] : B[n]]`。両辺に `#B[n]` を掛け、
`#B[n] · [A[n] : B[n]] = #A[n]` を使えば主張になる。**割り算は一度も現れない。**

## ★退化の自己検査——仮定を落とすと何が起きるか

* **`hB : B.index ≠ 0`(有限指数)を落とすと主張は偽になる**(自明化ではない)。
  `A = Multiplicative ℤ`、`B = ⊥`、`n = 2` を取る。
  `A^2 = 2ℤ` なので `[A : A^2] = 2`、`A[2] = {0}` なので `#A[2] = 1`、
  `B = ⊥` は自明群だから `[B : B^2] = 1`・`#B[2] = 1`。
  左辺 `2 · 1 = 2`、右辺 `1 · 1 = 1` で食い違う。
  ここで `B.index = [ℤ : 0] = 0`(`Nat.card` の規約で無限指数は `0`)であり、
  ちょうど仮定が破れている。証明のどこで効いたかも見える——上の
  `[A : B]  = [A^n : B^n] · [A[n] : B[n]]` が `0 = 1 · 1` になって壊れる。
* **`n = 0` と `n = 1` では自明に真**。`n = 0` なら `(·)^0` は自明準同型で
  `A^0 = ⊥`・`A[0] = A` だから両辺とも `#A · #B` で、形の上から等しい。
  `n = 1` なら両辺とも `1 · 1`。したがって `n` に関する仮定は要らない代わりに、
  小さい `n` は主張の中身を持たない。
* **可換性を落とすと述べることすらできない。** `powMonoidHom n` は
  `CommMonoid` でしか準同型にならず、`A^n` は非可換群では部分群ですらない。
* **`A` や `A[n]` の有限性は要らない。** 両辺が `0` になって釣り合うのは偶然ではなく、
  `hB` から「`A[n]` 無限 ⟺ `B[n]` 無限」「`[A:A^n]` 無限 ⟺ `[B:B^n]` 無限」が
  従うためである(上の二本の塔の等式が、有限・無限を問わず `Nat.card` の規約の
  ままで成り立つ)。

## どこで使うか

`Found/PGC/UnitsPowP.lean`——`p` 進局所体の単数群 `𝒪_K^×` は
`ℤ_p^{[K:ℚ_p]}` を有限指数部分群として含む。`ℤ_p^d` 側では
`p` 乗指数が `p^d`・`p` 捩れが自明なので、この不変性で `𝒪_K^×` 側へ移すと
`[𝒪_K^× : (𝒪_K^×)^p] = p^{[K:ℚ_p]} · #μ_p(K)` が出る([pGC] 経路 C の (C-d))。
-/

namespace ABC3.Found

open Subgroup

variable {A : Type*} [CommGroup A]

/-! ## `B^n := f(B)` の位置——`B` の中にも `A^n` の中にもある -/

/-- `B^n ≤ B`。 -/
theorem powMap_le_self (B : Subgroup A) (n : ℕ) :
    Subgroup.map (powMonoidHom n : A →* A) B ≤ B := by
  rintro _ ⟨b, hb, rfl⟩
  exact B.pow_mem hb n

/-- `B^n ≤ A^n`。 -/
theorem powMap_le_range (B : Subgroup A) (n : ℕ) :
    Subgroup.map (powMonoidHom n : A →* A) B ≤ (powMonoidHom n : A →* A).range := by
  rintro _ ⟨b, _, rfl⟩
  exact ⟨b, rfl⟩

/-! ## 鍵——`f⁻¹(B^n) = B · A[n]` -/

/-- **`n` 乗写像による `B^n` の逆像は `B · A[n]`**。

`a^n ∈ B^n` は `∃ b ∈ B, b^n = a^n` のこと。可換群なので `(b⁻¹a)^n = 1`、
すなわち `a = b · (b⁻¹a) ∈ B · A[n]`。逆の包含は `b^n ∈ B^n` と `1 ∈ B^n`。
**連結準同型の代わりに働くのがこの等式**である。 -/
theorem comap_powMap_eq_sup_ker (B : Subgroup A) (n : ℕ) :
    Subgroup.comap (powMonoidHom n : A →* A) (Subgroup.map (powMonoidHom n : A →* A) B)
      = B ⊔ (powMonoidHom n : A →* A).ker := by
  apply le_antisymm
  · intro a ha
    obtain ⟨b, hb, hba⟩ := ha
    have hrw : a = b * (b⁻¹ * a) := by group
    rw [hrw]
    refine Subgroup.mul_mem_sup hb ?_
    rw [MonoidHom.mem_ker, map_mul, map_inv, hba]
    simp
  · refine sup_le (fun b hb => Subgroup.mem_map_of_mem _ hb) (fun t ht => ?_)
    rw [MonoidHom.mem_ker] at ht
    simp only [Subgroup.mem_comap, ht]
    exact one_mem _

/-- **★★★★`[A : B · A[n]] = [A^n : B^n]`**——第一同型定理を
`ψ : A →* A ⧸ B^n`(`a ↦ a^n`)と `A^n ↪ A → A ⧸ B^n` の二つに当て、
どちらの像も `A ⧸ B^n` の中の同じ部分群 `A^n·B^n/B^n` であることを使う。
`Subgroup.index_ker`(`ker の指数 = range の濃度`)を二度使うだけで、
商群の濃度を直接計算しない。 -/
theorem index_sup_ker (B : Subgroup A) (n : ℕ) :
    (B ⊔ (powMonoidHom n : A →* A).ker).index
      = (Subgroup.map (powMonoidHom n : A →* A) B).relIndex
          ((powMonoidHom n : A →* A).range) := by
  set f := (powMonoidHom n : A →* A) with hf
  set Bn := Subgroup.map f B with hBn
  have h1 : ((QuotientGroup.mk' Bn).comp f).ker = B ⊔ f.ker := by
    rw [← MonoidHom.comap_ker, QuotientGroup.ker_mk']
    exact comap_powMap_eq_sup_ker B n
  have h3 := Subgroup.index_ker ((QuotientGroup.mk' Bn).comp f)
  rw [h1, MonoidHom.range_comp] at h3
  have h4 := Subgroup.index_ker ((QuotientGroup.mk' Bn).comp (f.range).subtype)
  have h5 : ((QuotientGroup.mk' Bn).comp (f.range).subtype).ker = Bn.subgroupOf f.range := by
    rw [← MonoidHom.comap_ker, QuotientGroup.ker_mk']
    rfl
  rw [h5, MonoidHom.range_comp, Subgroup.range_subtype] at h4
  rw [h3, ← h4]
  rfl

/-! ## 群の同型で移る量

`Found/PGC/UnitsGroupInvariants.lean::index_powRange_of_mulEquiv`(`[G : G^n]` は
同型不変量)の捩れ側の相棒。Herbrand 商の**両方の因子**が同型不変量であることを
言うために要る。 -/

/-- `#G[n]` は群の同型不変量。 -/
theorem card_torsion_of_mulEquiv {G H : Type*} [Group G] [Group H] (e : G ≃* H) (n : ℕ) :
    Nat.card { g : G // g ^ n = 1 } = Nat.card { h : H // h ^ n = 1 } :=
  Nat.card_congr (Equiv.subtypeEquiv e.toEquiv (fun g => by
    constructor
    · intro h
      show e g ^ n = 1
      rw [← map_pow, h, map_one]
    · intro h
      exact e.injective (by rw [map_pow, map_one]; exact h)))

/-! ## 部分群の側の量を `A` の中の量に翻訳する -/

/-- `[B : B^n]`——`B` の中で測った `n` 乗部分群の指数は、`A` の中の `B^n` の
`B` に対する相対指数。 -/
theorem index_powRange_subgroup (B : Subgroup A) (n : ℕ) :
    ((powMonoidHom n : B →* B).range).index
      = (Subgroup.map (powMonoidHom n : A →* A) B).relIndex B := by
  have h : (Subgroup.map (powMonoidHom n : A →* A) B).subgroupOf B
      = (powMonoidHom n : B →* B).range := by
    ext b
    simp only [Subgroup.mem_subgroupOf, MonoidHom.mem_range]
    constructor
    · rintro ⟨a, ha, hab⟩
      exact ⟨⟨a, ha⟩, Subtype.ext hab⟩
    · rintro ⟨c, rfl⟩
      exact ⟨(c : A), c.2, rfl⟩
  rw [Subgroup.relIndex, h]

/-- `B[n] = A[n] ∩ B`——`B` の中の `n` 捩れは、`A[n]` の中の `B` の元。 -/
def torsionSubgroupOfEquiv (B : Subgroup A) (n : ℕ) :
    ↥(B.subgroupOf (powMonoidHom n : A →* A).ker) ≃ { b : B // b ^ n = 1 } where
  toFun x := ⟨⟨(x.1 : A), x.2⟩, Subtype.ext x.1.2⟩
  invFun y := ⟨⟨(y.1 : A), congrArg Subtype.val y.2⟩, y.1.2⟩
  left_inv _ := rfl
  right_inv _ := rfl

/-- `#B[n] · [A[n] : B[n]] = #A[n]`(有限・無限を問わず `Nat.card` の規約で成立)。 -/
theorem card_torsion_mul_relIndex (B : Subgroup A) (n : ℕ) :
    Nat.card { b : B // b ^ n = 1 } * B.relIndex ((powMonoidHom n : A →* A).ker)
      = Nat.card { a : A // a ^ n = 1 } := by
  have h := Subgroup.card_mul_index (B.subgroupOf (powMonoidHom n : A →* A).ker)
  rw [Nat.card_congr (torsionSubgroupOfEquiv B n)] at h
  rw [Subgroup.relIndex, h]
  exact (Nat.card_congr (Equiv.subtypeEquivRight
    (fun _ => by rw [MonoidHom.mem_ker]; exact Iff.rfl))).symm

/-! ## 本体 -/

/-- **[A : A^n] = [B : B^n] · [A[n] : B[n]]**——主定理の割り算のない中核。

`[A : B] ≠ 0` から `[A : B·A[n]] ≠ 0`(`B ≤ B·A[n]` なので指数は約数)が出るので、
二本の塔の等式から共通因子 `[A^n : B^n] = [A : B·A[n]]` を約せる。 -/
theorem index_powRange_eq_mul_relIndex (B : Subgroup A) (hB : B.index ≠ 0) (n : ℕ) :
    ((powMonoidHom n : A →* A).range).index
      = (Subgroup.map (powMonoidHom n : A →* A) B).relIndex B
        * B.relIndex ((powMonoidHom n : A →* A).ker) := by
  set f := (powMonoidHom n : A →* A) with hf
  set Bn := Subgroup.map f B with hBn
  have hne : (B ⊔ f.ker).index ≠ 0 := fun h =>
    hB (Nat.eq_zero_of_zero_dvd (h ▸ Subgroup.index_dvd_of_le (le_sup_left : B ≤ B ⊔ f.ker)))
  have e1 : Bn.relIndex f.range * (f.range).index = Bn.index :=
    Subgroup.relIndex_mul_index (powMap_le_range B n)
  have e2 : Bn.relIndex B * B.index = Bn.index :=
    Subgroup.relIndex_mul_index (powMap_le_self B n)
  have e3 : B.relIndex (B ⊔ f.ker) * (B ⊔ f.ker).index = B.index :=
    Subgroup.relIndex_mul_index le_sup_left
  have e5 : B.relIndex (B ⊔ f.ker) = B.relIndex f.ker := Subgroup.relIndex_sup_left f.ker B
  have e4 : (B ⊔ f.ker).index = Bn.relIndex f.range := index_sup_ker B n
  have h : (B ⊔ f.ker).index * (f.range).index
      = (B ⊔ f.ker).index * (Bn.relIndex B * B.relIndex f.ker) := by
    calc (B ⊔ f.ker).index * (f.range).index
        = Bn.relIndex f.range * (f.range).index := by rw [e4]
      _ = Bn.index := e1
      _ = Bn.relIndex B * B.index := e2.symm
      _ = Bn.relIndex B * (B.relIndex (B ⊔ f.ker) * (B ⊔ f.ker).index) := by rw [e3]
      _ = (B ⊔ f.ker).index * (Bn.relIndex B * B.relIndex f.ker) := by rw [e5]; ring
  exact Nat.eq_of_mul_eq_mul_left (Nat.pos_of_ne_zero hne) h

/-- **★★★★★★★★Herbrand 商は有限指数部分群で不変**:

```
[A : A^n] · #B[n]  =  [B : B^n] · #A[n]
```

`hB : B.index ≠ 0`(= `B` は有限指数)が唯一の仮定。これを落とすと
`A = Multiplicative ℤ`・`B = ⊥`・`n = 2` で `2 · 1 ≠ 1 · 1` となり**偽**になる
(ファイル冒頭の「退化の自己検査」)。 -/
theorem index_pow_mul_card_torsion (B : Subgroup A) (hB : B.index ≠ 0) (n : ℕ) :
    ((powMonoidHom n : A →* A).range).index * Nat.card { b : B // b ^ n = 1 }
      = ((powMonoidHom n : B →* B).range).index * Nat.card { a : A // a ^ n = 1 } := by
  rw [index_powRange_eq_mul_relIndex B hB n, index_powRange_subgroup B n, mul_assoc,
    mul_comm (B.relIndex ((powMonoidHom n : A →* A).ker)), card_torsion_mul_relIndex B n]

end ABC3.Found
