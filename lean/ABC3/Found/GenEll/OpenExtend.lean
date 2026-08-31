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

/-! ## ★★★★★★空でない開は生成点を含む -/

/-- ★★★★**整域の `Spec` では、点を含む開集合は生成点（`⊥`）を含む**。

★機構: 任意の点 `x` は `⊥` の閉包（＝`zeroLocus ⊥` ＝全体）に入るので、
`x` の任意の開近傍は `{⊥}` と交わる。 -/
theorem generic_mem_of_mem {R : Type} [CommRing R] [IsDomain R]
    (U : Set (PrimeSpectrum R)) (hU : IsOpen U) {x : PrimeSpectrum R} (hx : x ∈ U) :
    (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum R) ∈ U := by
  have hcl : x ∈ closure ({(⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum R)} :
      Set (PrimeSpectrum R)) := by
    rw [PrimeSpectrum.closure_singleton, PrimeSpectrum.zeroLocus_bot]
    exact Set.mem_univ x
  obtain ⟨y, hyU, hy⟩ := mem_closure_iff.1 hcl U hU hx
  rw [Set.mem_singleton_iff] at hy
  exact hy ▸ hyU

omit [NumberField F] in
/-- ★★★★★★**`Spec 𝒪_F` では、点を含む開集合は生成点を含む**。

★★★これが「閉点を全部覆えば自動的に全体を覆う」ことの根拠である
——`Spec 𝒪_F` の点は生成点と閉点しかない。 -/
theorem generic_mem_opens (U : (specRingOfIntegers F).Opens)
    {x : specRingOfIntegers F} (hx : x ∈ U) :
    (⟨⊥, Ideal.isPrime_bot⟩ : specRingOfIntegers F) ∈ U := by
  have hx' : x ∈ (U : Set (specRingOfIntegers F)) := hx
  show (⟨⊥, Ideal.isPrime_bot⟩ : specRingOfIntegers F) ∈ (U : Set (specRingOfIntegers F))
  exact generic_mem_of_mem (R := NumberField.RingOfIntegers F) _ U.2 hx'

/-! ## ★★★★★★★★生成点を茎経由で見る -/

/-- ★★**茎 `𝒪_{Spec 𝒪_F, v}` から `F` への環準同型**（`Spec.stalkIso` 経由）。 -/
noncomputable def stalkToF (v : FinitePlace F) :
    (specRingOfIntegers F).presheaf.stalk (placePoint F v) ⟶ CommRingCat.of F :=
  (Spec.stalkIso (CommRingCat.of (NumberField.RingOfIntegers F)) (placePoint F v)).hom ≫
    CommRingCat.ofHom (algebraMap (Localization.AtPrime v.asIdeal) F)

/-- ★★★★★★**茎を経由した生成点は、素直な生成点と同じ**。

★機構は 3 本の合成:
`Spec.fromSpecStalk_eq`（アフィンでの `fromSpecStalk` の記述）、
`Spec.germ_stalkMapIso_hom`（芽と `stalkIso` の両立）、
`IsScalarTower.algebraMap_eq`（𝒪_F → (𝒪_F)_v → F）。

★★`rw` が芽（`germ`）の型で落ちるので、
`congrArg` と `calc` で項で繋いでいる（`tools/lean-idioms.md`）。 -/
theorem specMap_stalkToF_fromSpecStalk (v : FinitePlace F) :
    Spec.map (stalkToF F v) ≫ (specRingOfIntegers F).fromSpecStalk (placePoint F v)
      = Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) F)) := by
  set R := CommRingCat.of (NumberField.RingOfIntegers F) with hR
  set x := placePoint F v with hx
  set alg : CommRingCat.of (Localization.AtPrime v.asIdeal) ⟶ CommRingCat.of F :=
    CommRingCat.ofHom (algebraMap (Localization.AtPrime v.asIdeal) F) with halg
  rw [Spec.fromSpecStalk_eq, ← Spec.map_comp]
  congr 1
  have key := Spec.germ_stalkMapIso_hom (R := R) x
  calc ((Scheme.ΓSpecIso R).inv ≫ (Spec R).presheaf.germ ⊤ x trivial) ≫
        ((Spec.stalkIso R x).hom ≫ alg)
      = (Scheme.ΓSpecIso R).inv ≫
          ((Spec R).presheaf.germ ⊤ x trivial ≫ (Spec.stalkIso R x).hom) ≫ alg := by
        simp only [Category.assoc]
    _ = (Scheme.ΓSpecIso R).inv ≫
          ((Scheme.ΓSpecIso R).hom ≫ CommRingCat.ofHom (algebraMap _ _)) ≫ alg :=
        congrArg (fun z => (Scheme.ΓSpecIso R).inv ≫ z ≫ alg) key
    _ = CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F)
          (Localization.AtPrime v.asIdeal)) ≫ alg := by
        simp only [← Category.assoc, Iso.inv_hom_id, Category.id_comp]
    _ = CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) F) := by
        rw [halg, ← CommRingCat.ofHom_comp]
        congr 1
        exact (IsScalarTower.algebraMap_eq (NumberField.RingOfIntegers F)
          (Localization.AtPrime v.asIdeal) F).symm

/-- ★★★★**開近傍の中で見た生成点を `X` へ戻すと、素直な生成点**。 -/
theorem generic_comp_ι (v : FinitePlace F)
    (U : (specRingOfIntegers F).Opens) (hxU : placePoint F v ∈ U) :
    (Spec.map (stalkToF F v) ≫ U.fromSpecStalkOfMem _ hxU) ≫ U.ι
      = Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) F)) := by
  rw [Category.assoc, AlgebraicGeometry.Scheme.Opens.fromSpecStalkOfMem_ι,
    specMap_stalkToF_fromSpecStalk]

/-- ★★★★★★★★**開近傍へ広げた射を生成点へ制限すると、元の `F`-点に戻る**。

★★★これが貼り合わせの整合条件の根拠である
——どの `U_v` で見ても生成点では同じ `xK` になる。 -/
theorem generic_comp_open_extend {X' : Scheme.{0}}
    (v : FinitePlace F)
    (xv : Spec (CommRingCat.of (Localization.AtPrime v.asIdeal)) ⟶ X')
    (U : (specRingOfIntegers F).Opens) (hxU : placePoint F v ∈ U)
    (g : U.toScheme ⟶ X')
    (h : Spec.map (Spec.stalkIso (CommRingCat.of (NumberField.RingOfIntegers F))
          (placePoint F v)).inv ≫ xv = U.fromSpecStalkOfMem _ hxU ≫ g)
    (xK : Spec (CommRingCat.of F) ⟶ X')
    (hxv : Spec.map (CommRingCat.ofHom
        (algebraMap (Localization.AtPrime v.asIdeal) F)) ≫ xv = xK) :
    (Spec.map (stalkToF F v) ≫ U.fromSpecStalkOfMem _ hxU) ≫ g = xK := by
  rw [Category.assoc, ← h, ← Category.assoc, ← Spec.map_comp, stalkToF,
    ← Category.assoc, Iso.inv_hom_id, Category.id_comp]
  exact hxv

/-! ## ★★★★★★★★開近傍は全体を覆う -/

/-- ★**数体の整数環には素点が存在する**。

★体ではないので単元でない非零元があり、それを含む極大イデアルは `≠ ⊥`。 -/
theorem nonempty_finitePlace : Nonempty (FinitePlace F) := by
  obtain ⟨p, hp0, hp⟩ := Ring.exists_not_isUnit_of_not_isField
    (NumberField.RingOfIntegers.not_isField F)
  obtain ⟨m, hm, hpm⟩ := Ideal.exists_le_maximal (Ideal.span {p}) (Ideal.span_singleton_ne_top hp)
  refine ⟨⟨m, hm.isPrime, ?_⟩⟩
  intro h
  rw [h] at hpm
  have hmem : p ∈ (⊥ : Ideal (NumberField.RingOfIntegers F)) :=
    hpm (Ideal.mem_span_singleton_self p)
  exact hp0 (by simpa using hmem)

/-- ★★★★★★★★**各素点の開近傍は `Spec 𝒪_F` 全体を覆う**。

★`Spec 𝒪_F` の点は**生成点と閉点しかない**:
閉点 `x` は `x` 自身の近傍 `U_x` に入り、
生成点は空でない任意の開に入る（`generic_mem_opens`）。

★★★これで `Scheme.Opens.iSupOpenCover` が `Spec 𝒪_F` 自身の開被覆になる。 -/
theorem iSup_opens_eq_top (U : FinitePlace F → (specRingOfIntegers F).Opens)
    (hxU : ∀ v, placePoint F v ∈ U v) : ⨆ v, U v = ⊤ := by
  refine le_antisymm le_top (fun x _ => ?_)
  show x ∈ (⨆ v, U v : (specRingOfIntegers F).Opens)
  rw [TopologicalSpace.Opens.mem_iSup]
  by_cases h : x.asIdeal = ⊥
  · obtain ⟨v₀⟩ := nonempty_finitePlace F
    refine ⟨v₀, ?_⟩
    have hb := generic_mem_opens F (U v₀) (hxU v₀)
    have hxeq : x = (⟨⊥, Ideal.isPrime_bot⟩ : specRingOfIntegers F) := PrimeSpectrum.ext h
    rw [hxeq]
    exact hb
  · exact ⟨⟨x.asIdeal, x.isPrime, h⟩, hxU _⟩

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
