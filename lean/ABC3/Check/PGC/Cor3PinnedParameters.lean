import ABC3.Check.PGC.FreeTermFunctionRefutation
import ABC3.Found.PGC.Section3RealParameters

/-!
# [pGC] §3 —— 自由パラメータを実物に固定すると、既知の反例は作れなくなる

## 1. 自由版が退化しているかの判定は既に済んでいた(2026-09-08 の測定)

`Skeleton/PGC/Section3.lean` の `cor_3_1` / `cor_3_3` は、いま在る形のままだと
**偽**である。これは本セッションで新たに証明したものではなく、
`Check/PGC/FreeTermFunctionRefutation.lean`(2026-09-06)が既に `sorry` 無しで
出していた:

* `not_cor_3_1_current_form` —— `isHodgeTate` が項 `K` の自由な述語だから。
* `not_cor_3_3_current_form` —— `toGal` が項 `K` の自由な関数だから。

どちらも `selfField p` と `twistedField p`(ℚ_p-代数として同型だが項としては相異なる
2 つの `PAdicLocalField p`)に別の値を割り当てることで倒れる。
★したがって `Check/PGC/Prop12ForallRD.lean` と同じ形の「同値による測定」を
新しく書く必要は無かった —— こちらは同値より強い**反証**が既に在る。

## 2. 固定すると、その反例は作れない

`Found/PGC/Section3RealParameters.lean::artinToGal` は Artin 写像
(`Found/PGC/ArtinMap.lean::exists_artinMap`、仮定ゼロ)を `Γ_K` へ持ち上げたもので、
**単射**である(`artinToGal_injective`)。ところが反例に使われた `badToGal` は
捻り側で定数 `1` を返す。`U_K` には `1` 以外の元があるから、定数写像は単射でない。
ゆえに

    badToGal p (twistedField p) ≠ artinToGal (twistedField p) D   (どの `D` についても)

が言える(`badToGal_ne_artinToGal`)。★これは「固定した形が真である」ことの証明では
**ない** —— 「既知の反例がその形では作れない」ことだけを言う。

## 3. 固定した形は空虚でもない

`toGal` を実物に固定したうえで、`IsUniformizing` は

* 真になる `ρ` を持ち(`Found/PGC/Section3RealParameters.lean::isUniformizing_artin`)、
* 偽になる `ρ` を持つ(`ρ = 1`、`Check/PGC/Def32Degenerate.lean::not_isUniformizing_one`)。

したがって `Cor33Pinned` の `↔` は「両辺が恒等的に偽なので自明に真」ではない
(`Skeleton/PGC/Section3.lean` の docstring が書いていた板挟みの片側は解けた)。

これは原典の主張ではない(我々のモデルと器具についての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Found.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★反例に使われた `badToGal` は、どの Artin データから来る `toGal` とも一致しない。 -/
theorem badToGal_ne_artinToGal (p : ℕ) [Fact p.Prime] (D : ArtinDatum (twistedField p)) :
    badToGal p (twistedField p) ≠ artinToGal (twistedField p) D := by
  rw [badToGal_twisted]
  exact fun h => artinToGal_ne_const (twistedField p) D 1 h.symm

#print axioms badToGal_ne_artinToGal

/-- ★実物の `toGal` の下でも `ρ = 1` は uniformizing でない —— 固定した `↔` は
「両辺が恒等的に偽」で自明に真になることはない。 -/
theorem not_isUniformizing_artin_one (K : PAdicLocalField p) (D : ArtinDatum K) :
    ¬ IsUniformizing K K.carrier (artinToGal K D) 1 :=
  not_isUniformizing_one K K.carrier (artinToGal K D)

/-- ★実物の `toGal` の下で uniformizing になる `ρ` は在る(再掲)。 -/
theorem exists_isUniformizing_artin (K : PAdicLocalField p) (D : ArtinDatum K) :
    ∃ ρ : K.absGal →* (K.carrier)ˣ, IsUniformizing K K.carrier (artinToGal K D) ρ :=
  ⟨artinUnitChar K D, isUniformizing_artin K D⟩

end ABC3.Check.PGC
