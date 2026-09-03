import ABC3.Meta.Claim
import Mathlib.Algebra.CharP.Algebra
import Mathlib.NumberTheory.Padics.PadicIntegers
import Mathlib.NumberTheory.Padics.RingHoms
import Mathlib.RingTheory.AdicCompletion.LocalRing
import Mathlib.RingTheory.DiscreteValuationRing.Basic
import Mathlib.RingTheory.LocalRing.ResidueField.Defs

/-!
# [LocProP] §0 —— `p`-adic field / `p`-adic local field

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、物理 p.14。
**260 dpi 目視確認 2026-09-04**(このページに逐語の食い違いは無い)。

原文 (LocProP p.14):
> Definition 0.1. We shall call K a p-adic field if it is the quotient field of a
> p-adically complete, mixed characteristic discrete valuation ring OK.

## ★設計 —— `𝒪_K` を主対象に取る

原文は「`K` は `p`-adically complete・mixed characteristic な離散付値環 `𝒪_K` の
商体である」と定義する。`K` を主語にすると存在量化(`∃ 𝒪_K, …`)が要るが、
`𝒪_K`(本ファイルでは `O`)を主語に取れば mathlib の道具がそのまま乗る
——`IsDiscreteValuationRing` / `IsAdicComplete` / `CharZero` はどれも環に対する述語。
`K := FractionRing O` は必要なときに取り出せる。

## mathlib 在庫(実測、2026-09-04)

| 要る条件 | mathlib |
|---|---|
| 離散付値環 | `IsDiscreteValuationRing`(`RingTheory/DiscreteValuationRing/Basic.lean`) |
| mixed characteristic(`𝒪_K` は標数 0、剰余体は標数 `p`) | `CharZero O` + `CharP (IsLocalRing.ResidueField O) p` |
| `p`-adically complete | `IsAdicComplete (IsLocalRing.maximalIdeal O) O`(`RingTheory/AdicCompletion/*.lean`) |
| 剰余体が有限(`p`-adic **local** field) | `Finite (IsLocalRing.ResidueField O)` |

★★**`ℤ_[p]` は 4 条件すべてを instance で満たす** —— 特に `p`-adic 完備性は
`PadicInt.instIsAdicCompleteMaximalIdeal` が既にある(2026-09-04 実測、
`exact?` 1 発)。剰余体の同一視は `PadicInt.residueField : ResidueField ℤ_[p] ≃+* ZMod p`。
-/

namespace ABC3.Found.LocProP

/-- **[LocProP] Definition 0.1** の前半 —— `O` が `p`-adically complete・
mixed characteristic な離散付値環であること(`𝒪_K` の側から見た「`K` は `p`-adic field」)。

★`p` は `outParam`——`O` から一意に決まる(剰余標数)ので、使う側で毎回書かなくてよい。 -/
class IsPAdicRing (p : outParam ℕ) [Fact p.Prime] (O : Type*) [CommRing O]
    [IsDomain O] : Prop extends IsDiscreteValuationRing O, CharZero O where
  /-- `𝔪_K`-進位相で完備(原文の "p-adically complete")。 -/
  adicComplete : IsAdicComplete (IsLocalRing.maximalIdeal O) O
  /-- 剰余体 `k` の標数が `p`(原文の "mixed characteristic … discrete valuation ring"、
  `𝒪_K` 自身は標数 0 なので mixed)。 -/
  residueCharP : CharP (IsLocalRing.ResidueField O) p

/-- **[LocProP] Definition 0.1** の後半 —— 剰余体 `k` が有限なら `K` は
"`p`-adic **local** field"。 -/
class IsPAdicLocalRing (p : outParam ℕ) [Fact p.Prime] (O : Type*) [CommRing O]
    [IsDomain O] : Prop extends IsPAdicRing p O where
  /-- 原文: "We shall call K a p-adic local field if k is a finite field." -/
  residueFinite : Finite (IsLocalRing.ResidueField O)

/-- ★★**非退化性の witness** —— `ℤ_[p]` は実際に `IsPAdicLocalRing p` を満たす。 -/
instance instIsPAdicLocalRing (p : ℕ) [Fact p.Prime] : IsPAdicLocalRing p ℤ_[p] where
  adicComplete := inferInstance
  residueCharP := by
    have e := PadicInt.residueField (p := p)
    haveI : CharP (ZMod p) p := ZMod.charP p
    exact CharP.of_ringHom_of_ne_zero e.symm.toRingHom p (Fact.out (p := p.Prime)).pos.ne'
  residueFinite := by
    have e := PadicInt.residueField (p := p)
    exact Finite.of_equiv (ZMod p) e.symm.toEquiv

/-- `K`(`p`-adic field 本体)は `𝒪_K` の商体として取り出す。 -/
abbrev pAdicField (O : Type*) [CommRing O] [IsDomain O] : Type _ := FractionRing O

end ABC3.Found.LocProP
