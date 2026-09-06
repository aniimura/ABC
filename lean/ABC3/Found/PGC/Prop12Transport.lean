import ABC3.Skeleton.PGC.Section1Defs
import ABC3.Found.PGC.ResidueCardTransport
import ABC3.Found.PGC.DegreeTransport
import ABC3.Found.PGC.ResidueCardinality

/-!
# [pGC] Proposition 1.2 の**無条件版**——経路 C の合流点(ノード J)

原文 (pGC p.3):

> The number q of elements in the residue field of O[scr]_K, and well as the absolute
> degree [K : Q[bb]_p] of K, can be recovered entirely group-theoretically from Γ_K.

## ★2026-09-06(D13 の実行): Skeleton 側が本ファイルと同じ形になった

本ファイルを書いた時点では `Skeleton/PGC/Section1.lean::residueCard_and_degree_recoverable` が
`(RD : ResidueCardinality p)` を**仮説として量化した**条件付き形式化で、そこは `sorry` だった
(理由は下の「★量化された形は Prop 1.2 より強い」)。
2026-09-06 に D13 を採り、Skeleton 側も `realResidueCardinality` に固定して**閉じた**ので、
本ファイルの定理は Skeleton 側と同じ内容の別名になっている。
★下の「★量化された形は Prop 1.2 より強い」の分析が D13 の根拠そのものなので、
記録としてこのファイルを残す。

本ファイルはその量化を外し、**実物の `RD`**
(`Found/PGC/ResidueCardinality.lean::realResidueCardinality`、`card = residueCard`)に
限った形を**無条件に**証明する。`sorry` は無い。

## 何と何を貼り合わせたか

| 成分 | 補題 | 置き場所 |
|---|---|---|
| `q` | `residueCard_eq_of_absGal_equiv` | `Found/PGC/ResidueCardTransport.lean`(経路 C ノード G) |
| `[K:ℚ_p]` | `finrank_eq_of_absGal_equiv` | `Found/PGC/DegreeTransport.lean`(経路 C (C-d)) |

両者はともに「`Γ_K ≃ₜ* Γ_{K'}` なら値が等しい」という形をしており、
`residueCardAndDegreeObject` の `transport` は恒等なので、対 `(q, [K:ℚ_p])` の
等号をそのまま組めばよい(`Prod.ext`)。

## 橋は `rfl` だった

`realResidueCardinality.card = residueCard` であり(`realResidueCardinality_card`、`rfl`)、
`residueCard K = Nat.card 𝓀[K.carrier]` も定義そのもの
(`Found/PGC/LocalFieldNorm.lean:118`)。したがって q 側の補題(実際の剰余体の濃度に
ついての主張)は**定義的に**第 1 成分に届く。念のため
`realResidueCardinality_card_eq_natCard` として明示しておく(証明は `rfl`)。

## ★量化された形(`∀ RD`)は Prop 1.2 より強い——なぜここで量化を外すのか

`Skeleton` 側の `∀ RD : ResidueCardinality p, …` は、`ResidueCardinality` の 3 つの場
(`card` / `isPrimePow` / `card_congr`)しか使えないため、
**「Γ_K ≅ Γ_{K'} ならば K ≃ₐ[ℚ_p] K'」と同値**になってしまう
(証明は `Check/PGC/Prop12ForallRD.lean::forall_RD_recoverable_iff_algEquiv`)。
これは原典が Introduction で
「the Grothendieck Conjecture cannot hold in the naive sense」と述べて
[8] (Jarden-Ritter) を挙げている、まさにその主張である。
すなわち `∀ RD` 版は Proposition 1.2 ではなく、**原典が偽と述べている命題**に等しい。

本ファイルの無条件版はその落とし穴を通らない——実物の `card` は
`Nat.card 𝓀[K.carrier]` であって自由なデータではないので、
体の同型を経由せずに `Γ` の不変量(`contHomCard`)だけで値が決まる。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NormedField Valued

variable (p : ℕ) [Fact p.Prime]

/-- 実物の `RD` の `card` は、剰余体の実際の濃度そのもの(`rfl`)。

`realResidueCardinality.card = residueCard`(`realResidueCardinality_card`)であり、
`residueCard K = Nat.card 𝓀[K.carrier]` は定義。橋は 1 本も要らなかった。 -/
theorem realResidueCardinality_card_eq_natCard (K : PAdicLocalField p) :
    (realResidueCardinality p).card K = Nat.card 𝓀[K.carrier] := rfl

/-- **★★★★★★★★★★★★★★★★★★★★★★★★★★[pGC] Proposition 1.2
(実物の q についての無条件版)**。

`q = #𝓀[K]` と `[K:ℚ_p]` の対は、`Γ_K` の位相群としての同型類だけから決まる。

`Skeleton/PGC/Section1.lean::residueCard_and_degree_recoverable` の
`RD` を実物 `realResidueCardinality` に固定した形。仮説は無い。

第 1 成分は `Found/PGC/ResidueCardTransport.lean::residueCard_eq_of_absGal_equiv`
(不変量 `N_n(Γ) = #Hom_cont(Γ, ℤ/n)` の集合 `S(Γ)` の最大元が `q−1`)、
第 2 成分は `Found/PGC/DegreeTransport.lean::finrank_eq_of_absGal_equiv`
(Kummer 双対で `#Hom_cont(Γ_F, ℤ/p) = p^{[F:ℚ_p]+2}`)。

★経路 C は原典の論拠(局所類体論の相互律 + p 進対数)を経由しない。
その逸脱は `ResearchPaper/pgc-goal.md` に記録済み。 -/
theorem residueCard_and_degree_recoverable_real :
    (residueCardAndDegreeObject (realResidueCardinality p)).RecoverableFromAbsGal := by
  intro K K' α
  exact Prod.ext (residueCard_eq_of_absGal_equiv α) (finrank_eq_of_absGal_equiv K K' α)

def residueCard_and_degree_recoverable_real.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.2", sectionId := "prop-1-2" }

#print axioms residueCard_and_degree_recoverable_real

end ABC3.Found.PGC
