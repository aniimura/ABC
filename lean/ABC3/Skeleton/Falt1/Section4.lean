import ABC3.Meta.Claim
import ABC3.Interface.Falt1.Differentials
import ABC3.Interface.Falt1.GaloisCohomologyCompute

/-!
# [Falt1] Chapter I §4 —— Differentials and cohomology(5 項目)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)。物理 p.13-19(印字 p.266-272)。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。★§4 は 5/5。

Theorem 4.4 は `Skeleton/Falt1/ChapterI.lean`(既存、[LocProP] Lemma 2.1
が直接引用する箇所)にある——本ファイルには 4.1・4.2・4.3・4.5 を置く。
-/

namespace ABC3.Skeleton.Falt1

open ABC3.Interface.Falt1

/-- **`[Falt1] Theorem 4.1`**。

内容 (Falt1 p.13-14 = 印字 p.266-267): (i) The map `Ω_{R/V} → Ω_{R̄/V̄}`
induces almost isomorphisms `Ω_{R/V} ⊗_R R̄[1/p] ≅ Ω_{R̄/V̄}` and
`Ω_{R/V} ⊗_R (R̄[1/p]/R̄) ≅ Ω_{R̄/(R⊗_V V̄)}`. (ii) The sequence
`0 → Ω_{V̄/V} ⊗_V̄ R̄ → Ω_{R̄/R} → Ω_{R̄/(R⊗_V V̄)} → 0` is exact.

★**2026-09-05 に 260dpi で描画して逐語再確認し、転記を訂正した**。
以前の記述は `Ω̄_{R/V}`(`Ω` の上に bar)と書いていたが、原典の bar は
**`Ω` ではなく添字の `R`・`V` に付いている**——すなわち `Ω_{R̄/V̄}` であり、
`p` 進**完備化ではない**。この違いは形式化の難易度を大きく変える
(完備化の枠組みが不要で、通常の Kähler 微分と局所化 `R̄[1/p]` だけで
書ける)。(ii) の第 2 項も `Ω̄_{R̄/R}` ではなく `Ω_{R̄/R}` である。 -/
def theorem_4_1 (E : DifferentialsSetup) :=
  PProd.mk E.thm41i1_ker <| PProd.mk E.thm41i1_coker <| PProd.mk E.thm41i2_ker <|
  PProd.mk E.thm41i2_coker <| PProd.mk E.thm41ii_injective <|
  PProd.mk E.thm41ii_exact E.thm41ii_surjective

def theorem_4_1.nonvacuous := theorem_4_1 DifferentialsSetup.example

def theorem_4_1.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 13, item := "Theorem 4.1", sectionId := "falt1-thm-4-1" }

/-- **`[Falt1] Theorem 4.2`**。

内容 (Falt1 p.15、260dpi 目視、★簡略化——Interface/Falt1/
GaloisCohomologyCompute.lean 冒頭の注記参照): (i) H^i(Δ,R̂) ≅
Λ^i((R⊗_V V̄)^(-1)^d) ⊕ (rest), rest annihilated by p^{1/(p-1)},
Gal(K̄/K)-equivariant. (ii) étale base change induces almost
isomorphisms. (iii) Künneth morphisms are almost isomorphisms.
(iv) H^1(Δ,R̂)/(p-torsion) ≅ Ω_{R/V}⊗(R⊗V̄)^(-1) via a functorial
isomorphism. (v) cup products induce H^i(Δ,R̂)/(p-torsion) ≅
Ω^i_{R/V}⊗(R⊗V̄)^(-i). (vi) cup product duality induces almost
isomorphisms. -/
def theorem_4_2 (E : GaloisCohCompute42) :=
  (E.thm42i, E.thm42i_rest, E.thm42ii_ker_ann, E.thm42ii_coker_ann,
   E.thm42iii_ker_ann, E.thm42iii_coker_ann, E.thm42iv, E.thm42v,
   E.thm42vi_ker_ann, E.thm42vi_coker_ann)

def theorem_4_2.nonvacuous := theorem_4_2 GaloisCohCompute42.example

def theorem_4_2.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 15, item := "Theorem 4.2", sectionId := "falt1-thm-4-2" }

/-- **`[Falt1] Theorem 4.3`**。

内容 (Falt1 p.17、260dpi 目視): There exists a Γ-equivariant functorial
extension 0→ρ⁻¹R̂→E_ρ→Ω_{R/V}(dlog∞)⊗R̂(-1)→0. This extension is
obtained via pushout by R̂→ρ⁻¹R̂ from an extension
0→R̂→E→Ω_{R/V}(dlog∞)⊗R̂(-1)→0. -/
def theorem_4_3 (E : GaloisCohCompute43) :=
  And.intro E.thm43a_inj <| And.intro E.thm43a_exact <| And.intro E.thm43a_surj <|
  And.intro E.thm43b_inj <| And.intro E.thm43b_exact E.thm43b_surj

def theorem_4_3.nonvacuous := theorem_4_3 GaloisCohCompute43.example

def theorem_4_3.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 17, item := "Theorem 4.3", sectionId := "falt1-thm-4-3" }

/-- **`[Falt1] Theorem 4.5`**。

内容 (Falt1 p.19-20、260dpi 目視、★簡略化——Interface/Falt1/
GaloisCohomologyCompute.lean 冒頭の注記参照): (i) The spectral sequence
E_2^{a,b}=H^a(Δ,H^b(I,R̂)) ⇒ H^{a+b}(Δ,R̂) degenerates up to m-torsion.
(ii) The induced mappings are almost isomorphisms. (iii) Under the
isomorphism of 4.2(v), the filtration on H^i(Δ,R̂) corresponds to the
canonical filtration on Ω^i_{R/V}(dlog∞). -/
def theorem_4_5 (E : GaloisCohCompute45) :=
  PProd.mk E.thm45i_ker_ann <| PProd.mk E.thm45i_coker_ann <|
  PProd.mk E.thm45ii_ker_ann <| PProd.mk E.thm45ii_coker_ann E.thm45iii_holds

def theorem_4_5.nonvacuous := theorem_4_5 GaloisCohCompute45.example

def theorem_4_5.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 19, item := "Theorem 4.5", sectionId := "falt1-thm-4-5" }

end ABC3.Skeleton.Falt1
