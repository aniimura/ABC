/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DedekindExtend
import Mathlib.AlgebraicGeometry.SpreadingOut

/-!
# [GenEll] Remark 1.5.1 —— **局所延長を開近傍へ広げる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

## ★★★★★★★★段 2c の前半

`DedekindExtend.lean` で各素点 `v` の `(𝓞_F)_v`-点を得た（段 2a）。
★`Spec (𝓞_F)_v` は **Zariski 開ではない**ので、そのままでは貼り合わせられない。

★★そこで mathlib の `spread_out_of_isGermInjective′`
（`SpreadingOut.lean`、Stacks 0BX6）で**開近傍へ広げる**:

    Spec 𝒪_{X,v} ⟶ X′   ⟹   ∃ U ∋ v 開,  U ⟶ X′

★★★`Spec 𝓞_F` は**整**なので `IsGermInjective` は自動で付く
——mathlib の `instance (priority := 100) [IsIntegral X] : X.IsGermInjective`。
★底が `Spec ℤ`（終対象）なので、`S`-射であるという条件も自動。

## ★★★★茎と局所化の橋

`spread_out_of_isGermInjective′` は**茎**（`X.presheaf.stalk x`）で述べられている。
★`(𝓞_F)_v` と茎を繋ぐのは mathlib の
`Spec.stalkIso R x : (Spec R).presheaf.stalk x ≅ .of (Localization.AtPrime x.asIdeal)`。

## ★残っている段 2c の後半

★開近傍 `{U_v}` は `Spec 𝓞_F` を覆う——**空でない開は生成点を含む**ので、
閉点を全部覆えば自動的に全体を覆う。
★★あとは `Scheme.OpenCover.glueMorphisms` で貼る。整合は
`extend_unique` と同じ機構（生成点の稠密性）で出る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory IsDedekindDomain

variable (F : Type) [Field F] [NumberField F]

/-- ★素点 `v` を `Spec 𝓞_F` の点として見る。 -/
def placePoint (v : FinitePlace F) : specRingOfIntegers F := ⟨v.asIdeal, v.isPrime⟩

/-- ★★★★★★★★**局所延長は開近傍へ広がる**（Stacks 0BX6）。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

★★`Spec (𝓞_F)_v` は Zariski 開ではないので貼り合わせに使えない。
**開近傍 `U ∋ v` へ広げる**のが本補題である。

★★★機構は mathlib の `spread_out_of_isGermInjective′` 1 本——
`Spec 𝓞_F` は整なので `IsGermInjective` は自動、
底が `Spec ℤ`（終対象）なので `S`-射条件も自動。 -/
theorem exists_open_extend {X' : Scheme.{0}}
    (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    (v : FinitePlace F)
    (xv : Spec (CommRingCat.of (Localization.AtPrime v.asIdeal)) ⟶ X') :
    ∃ (U : (specRingOfIntegers F).Opens) (hxU : placePoint F v ∈ U)
      (g : U.toScheme ⟶ X'),
      Spec.map (Spec.stalkIso (CommRingCat.of (NumberField.RingOfIntegers F))
          (placePoint F v)).inv ≫ xv
        = U.fromSpecStalkOfMem _ hxU ≫ g := by
  obtain ⟨U, hxU, g, h1, -⟩ := spread_out_of_isGermInjective'
    (sX := specZIsTerminal.from (specRingOfIntegers F)) (sY := f')
    (x := placePoint F v)
    (Spec.map (Spec.stalkIso (CommRingCat.of (NumberField.RingOfIntegers F))
      (placePoint F v)).inv ≫ xv)
    (specZIsTerminal.hom_ext _ _)
  exact ⟨U, hxU, g, h1⟩

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** 残っているのは
開近傍 `{U_v}` から `Scheme.OpenCover` を作って `glueMorphisms` で貼る段である。 -/

def exists_open_extend.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(固有性から点を延ばす——局所延長を開近傍へ広げる)",
    sectionId := "genell-def-1-5" }

def exists_open_extend.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "spread_out_of_isGermInjective′(Stacks 0BX6)"
      (.inMathlib "AlgebraicGeometry.spread_out_of_isGermInjective'") 8,
    .citation "[mathlib]" "Spec.stalkIso(茎と AtPrime 局所化の同一視)"
      (.inMathlib "AlgebraicGeometry.Spec.stalkIso") 8,
    .citation "[ABC3]" "exists_unique_extend_atPrime(各素点での延長)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_unique_extend_atPrime") 8,
    .implicitStep
      ("★Spec 𝓞_F は整なので IsGermInjective は自動で付く" ++
       "(mathlib の instance (priority := 100) [IsIntegral X] : X.IsGermInjective)。" ++
       "★★底が Spec ℤ(終対象)なので S-射条件も自動") 8,
    .implicitStep
      ("★★★残る段: 開近傍 {U_v} から Scheme.OpenCover を作って glueMorphisms で貼る。" ++
       "★空でない開は生成点を含むので、" ++
       "閉点を全部覆えば自動的に全体を覆う。" ++
       "★★整合は extend_unique と同じ機構(生成点の稠密性)で出る") 8 ]

end ABC3.Found.GenEll
