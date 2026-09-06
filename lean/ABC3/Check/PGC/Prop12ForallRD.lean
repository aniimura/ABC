import ABC3.Skeleton.PGC.Section1
import ABC3.Found.PGC.QpResidueField

/-!
# [pGC] Proposition 1.2 の `∀ RD` 版は、原典が偽と述べている命題と**同値**である

原文 (pGC p.3):

> The number q of elements in the residue field of O[scr]_K, and well as the absolute
> degree [K : Q[bb]_p] of K, can be recovered entirely group-theoretically from Γ_K.

`Skeleton/PGC/Section1.lean` の形は

```
∀ RD : ResidueCardinality p, (residueCardAndDegreeObject RD).RecoverableFromAbsGal
```

である。本ファイルはこの形が

```
∀ K K', (Γ_K ≃ₜ* Γ_{K'}) → Nonempty (K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
```

と**同値**であることを証明する(`forall_RD_recoverable_iff_algEquiv`、`sorry` 無し)。

## なぜこれが問題なのか

右辺は「絶対 Galois 群だけから体の同型類が決まる」という主張——原典が Introduction で
次のように述べて否定している、まさにその命題である。

原文 (pGC p.1, Introduction):

> On the one hand, one knows (cf. the Remark in [4] following Theorem 4.2) that the
> Grothendieck Conjecture cannot hold in the naive sense (i.e., if one removes the condition
> of "compatibility with the filtrations" from the outer isomorphisms considered - see, e.g.,
> [8]), so one must put some sort of condition on the outer isomorphisms of Galois groups that
> one considers.

原文 (pGC p.1, Historical Remark):

> I originally set out to prove the naive version of the above Theorem, only to discover
> that this was, in fact, false.

[8] は M. Jarden, J. Ritter, On the Characterization of Local Fields by their Absolute
Galois Groups, J. Number Th. 11 (1979), pp. 1-13——同型でない p 進局所体で絶対 Galois 群が
位相群として同型なものを与える文献である。右辺はその存在によって偽になる。

したがって `∀ RD` 版は Proposition 1.2 ではない。それは

* 原典の主定理(Theorem: `Isom_{Q_p}(K,K') → Out^{Filt}(Γ_K, Γ_{K'})` が全単射)よりも
  **強い**(こちらは分岐フィルターとの両立を課さない)、かつ
* 原典が偽と述べている「naive version」に等しい

命題である。**同型不変性(`card_congr`)を足した修理が行き過ぎた**——
修理前(`Check/PGC/Prop12Degenerate.lean`、第 1012)は安い捻りで偽になっていたが、
修理後は「安い反例が消えた」のであって「正しい主張になった」のではない。

## どこで行き過ぎたのか——`isPrimePow` の f は剰余次数ではない

`ResidueCardinality` は `card K = p^f`(`f > 0`)としか言わないので、`f` と
実際の剰余次数を結ぶ場が無い。`q = p^f`、`[K:ℚ_p] = e·f` という関係で `f` を挟もうとしても
挟めない——`exponent_not_determined` が示すとおり、**同じ体に**別々の `f` を割り当てる
`ResidueCardinality` の項が 2 つ作れる(`p` と `p^2`)。両者とも `card_congr` を満たす。

`card` に課されているのは「ℚ_p-代数同型で不変」だけなので、`card` は
**同型類の任意の関数**でよい。ゆえに `card` を「`K₀` の同型類の指示関数」に取れば、
`∀ RD` 版は「Γ が同型なら同じ同型類」を強制する。これが同値性の中身である。

## 正しい形

量化を外し、実物の `RD` に固定すればよい。それは
`Found/PGC/Prop12Transport.lean::residueCard_and_degree_recoverable_real` として
**無条件に**証明済み(`sorry` 無し)。実物の `card` は `Nat.card 𝓀[K.carrier]` であって
自由なデータではないので、体の同型を経由せずに `Γ` の不変量だけで値が決まる。

★`Skeleton/PGC/Section1.lean` の statement をどう直すかは人の判断待ち(方針書 §2)。
本ファイルは判断の材料を出すだけで、`Skeleton` には触れていない。

## 「落とした条件は主張を偽にするか自明にする」の 10 例目——ただし新種

これまでの 9 例は「条件を落とすと偽/自明になる」だった。9 例目 (D10) で
「素朴な修復が偽の主張を作る」が出た。本件はさらに別の形で、

> **修復が(偽になるほど)強すぎる主張を作った**

である。落ちている条件は「`RD` が実物であること」——すなわち
`card K = Nat.card 𝓀[K.carrier]` を要求する場である。

**これは原典の主張ではない**(我々のモデルと器具についての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Found.PGC ABC3.Interface.PGC
open scoped Classical

variable (p : ℕ) [Fact p.Prime]

/-! ## 1. `ResidueCardinality` の 3 つの場を通る、病的な項 -/

/-- `K₀` の ℚ_p-代数同型類の**指示関数**を `card` に持つ `ResidueCardinality`。

`K₀` と同型なら `p`、そうでなければ `p^2`。`isPrimePow` は f = 1 と f = 2 で通り、
`card_congr` は「同型類の指示関数は同型で不変」だから通る。
すなわち**現在の `Interface` の 3 つの場を全て満たす**。 -/
noncomputable def isoIndicatorRD (K₀ : PAdicLocalField p) : ResidueCardinality p where
  card X := if Nonempty (X.carrier ≃ₐ[ℚ_[p]] K₀.carrier) then p else p ^ 2
  isPrimePow X := by
    by_cases h : Nonempty (X.carrier ≃ₐ[ℚ_[p]] K₀.carrier)
    · exact ⟨1, one_pos, by simp [h]⟩
    · exact ⟨2, two_pos, by simp [h]⟩
  card_congr {X Y} e := by
    have h : Nonempty (X.carrier ≃ₐ[ℚ_[p]] K₀.carrier)
        ↔ Nonempty (Y.carrier ≃ₐ[ℚ_[p]] K₀.carrier) :=
      ⟨fun ⟨f⟩ => ⟨e.symm.trans f⟩, fun ⟨f⟩ => ⟨e.trans f⟩⟩
    simp only [h]

/-- 定数 `p^2` を返す `ResidueCardinality`。3 つの場を自明に満たす。 -/
noncomputable def constRD : ResidueCardinality p where
  card _ := p ^ 2
  isPrimePow _ := ⟨2, two_pos, rfl⟩
  card_congr _ := rfl

/-- **`isPrimePow` の f は体から決まらない**。

同じ体 `ℚ_[p]`(`selfField p`)に、許される 2 つの `ResidueCardinality` が
別の値 `p` と `p^2` を割り当てる。ゆえに「`q = p^f` と `[K:ℚ_p] = e·f` で `f` を挟む」
という道は `Interface` の現在の場からは取れない——`f` は剰余次数と結ばれていない。 -/
theorem exponent_not_determined :
    (isoIndicatorRD p (selfField p)).card (selfField p)
      ≠ (constRD p).card (selfField p) := by
  show (if Nonempty ((selfField p).carrier ≃ₐ[ℚ_[p]] (selfField p).carrier) then p else p ^ 2)
      ≠ p ^ 2
  rw [if_pos ⟨AlgEquiv.refl⟩]
  have hp : 1 < p := (Fact.out : p.Prime).one_lt
  nlinarith [sq_nonneg p]

/-! ## 2. 同値性 -/

/-- **★★★★★★★★★★★★★★★★★★★★`∀ RD` 版は「Γ から体の同型類が決まる」と同値**。

→ は `isoIndicatorRD` を当てるだけ。`K` の側は `AlgEquiv.refl` で `p` になり、
`K'` の側が `p^2` だと `p = p^2` から `p` が素数であることに矛盾する。
ゆえに `K'` は `K` と同型でなければならない。

← は `card_congr`(第 1 成分)と `LinearEquiv.finrank_eq`(第 2 成分)。
★← には次数側の深い定理(`finrank_eq_of_absGal_equiv`)は要らない——体の同型が
あれば次数は直ちに等しい。つまりこの同値は**まるごと `card` 成分の問題**である。 -/
theorem forall_RD_recoverable_iff_algEquiv :
    (∀ RD : ResidueCardinality p, (residueCardAndDegreeObject RD).RecoverableFromAbsGal) ↔
      ∀ {K K' : PAdicLocalField p}, (K.absGal ≃ₜ* K'.absGal) →
        Nonempty (K.carrier ≃ₐ[ℚ_[p]] K'.carrier) := by
  constructor
  · intro hall K K' α
    have h := congrArg Prod.fst (hall (isoIndicatorRD p K) α)
    simp only [residueCardAndDegreeObject, isoIndicatorRD] at h
    rw [if_pos ⟨AlgEquiv.refl⟩] at h
    by_cases hne : Nonempty (K'.carrier ≃ₐ[ℚ_[p]] K.carrier)
    · exact ⟨(Classical.choice hne).symm⟩
    · rw [if_neg hne] at h
      exact absurd h.symm (by
        have hp : 1 < p := (Fact.out : p.Prime).one_lt
        have : p < p ^ 2 := by nlinarith
        omega)
  · intro h RD K K' α
    obtain ⟨e⟩ := h α
    exact Prod.ext (RD.card_congr e) (LinearEquiv.finrank_eq e.toLinearEquiv)

#print axioms forall_RD_recoverable_iff_algEquiv

/-- **反例が 1 組でも出れば `∀ RD` 版は偽**。

[8] (Jarden-Ritter) が与えるのはまさにこの仮説——絶対 Galois 群が位相群として同型なのに
ℚ_p-代数として同型でない 2 つの p 進局所体である。その構成は本リポジトリには無いので、
本補題は「その日が来たら `∀ RD` 版が倒れる」ことだけを述べる。 -/
theorem not_forall_RD_recoverable_of_nonisomorphic {K K' : PAdicLocalField p}
    (α : K.absGal ≃ₜ* K'.absGal) (hne : ¬ Nonempty (K.carrier ≃ₐ[ℚ_[p]] K'.carrier)) :
    ¬ ∀ RD : ResidueCardinality p, (residueCardAndDegreeObject RD).RecoverableFromAbsGal :=
  fun h => hne ((forall_RD_recoverable_iff_algEquiv p).mp h α)

end ABC3.Check.PGC
