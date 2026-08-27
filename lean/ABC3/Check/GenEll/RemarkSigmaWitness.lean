import ABC3.Skeleton.GenEll.Section1

/-!
# ★[GenEll] `Remark 1.5.1` の仮定は空虚ではない(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> Remark 1.5.1. In the spirit of Remark 1.4.1, we observe that the log-different

## ★★★なぜこの検査が要るか

2026-08-27(第 420 ブロック)に `remark_1_5_1` を構成へ載せ替えたとき、
原文の証明の**前半(spreading out)を仮定 `hagree` として受けた**。

★**仮定を置いて `sorry` を消すのは、それ自体では進捗ではない**
——満たされない仮定を置けば statement は空虚になる。

★★本ファイルは `hagree` が**実際に満たされる場合がある**ことを機械検証する:
`D = D′`、`ePt = id`、`Σ = ∅` を取れば `hagree` は `rfl` で成り立ち、
`remark_1_5_1` の結論(`BDeq`)が実際に得られる。

★★★これは「非空虚 witness」であって「spreading out を実装した」ではない。
本当に要るのは **`Σ ≠ ∅` の場合に `hagree` を供給する幾何**であり、
それは `.needs` に `[Stacks] 32.22` として登記してある。
-/

namespace ABC3.Check.GenEll

open ABC3.Meta ABC3.Found.GenEll

/-- ★**`hagree` は満たされうる** —— 同じ因子を同じ点で見れば `Σ = ∅` でも一致する。 -/
theorem remark_1_5_1_hagree_satisfiable
    (F : Type) [Field F] [NumberField F] {X : AlgebraicGeometry.Scheme.{0}}
    (D : X.IdealSheafData) (ch : FinitePlace F → ℕ) :
    ∀ xF : specRingOfIntegers F ⟶ X, ∀ v : FinitePlace F, ch v ∉ (∅ : Finset ℕ) →
      (conductorADiv F D xF).fin v = (conductorADiv F D (id xF)).fin v :=
  fun _ _ _ => rfl

/-- ★★**結論が実際に出る** —— `remark_1_5_1` を満たされる仮定で叩くと `BDeq` が得られる。

★これで「仮定が偽だから何でも言える」ではないことが確かめられた。 -/
theorem remark_1_5_1_conclusion_obtained
    (F : Type) [Field F] [NumberField F] {X : AlgebraicGeometry.Scheme.{0}}
    (D : X.IdealSheafData)
    (ch : FinitePlace F → ℕ)
    (hover : ∀ v : FinitePlace F,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)}))
    (hI : ∀ xF, pullbackIdeal F D xF ≠ 0) :
    BDeq (fun xF => logCond F D xF) (fun xF => logCond F D xF) :=
  (ABC3.Skeleton.GenEll.remark_1_5_1 F D D id ∅ (by simp) ch hover hI hI
    (fun _ _ _ => rfl)).2

end ABC3.Check.GenEll
