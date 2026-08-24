/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop53UntrPf
import ABC3.Found.FrdI.Prop32Frob
import ABC3.Found.FrdI.Prop44Core
import ABC3.Found.FrdI.Thm51Sec

/-!
# [FrdI] Proposition 5.5, (i) の「任意の対象へ」の段 —— pre-step の対で降ろす

原文 (FrdI p.104):
> tion (i) follows immediately for Frobenius-trivial A by considering base-identity

原文 (FrdI p.104):
> the case of arbitrary A then follows by considering "pairs of pre-steps" as in Theorem

★★**原文が畳んだ 1 文がこれ**である。`Proposition 5.5, (i)` は
**Frobenius-trivial な `A`** については「即座に従う」が、**任意の `A`** については
`Theorem 5.1, (i)` の「**pre-step の対**」を経由する、と書いてある。

## ★我々の木で何が起きていたか(記録)

`Prop55PfUnit.lean` の `otimes_pfRoot_eq_bot` は
仮定 `hA : IsFrobeniusTrivial P A` を**そのまま持っている**。
そのため、そこから作られた

* `otimes_pfRoot_eq_bot_of_root`(`Prop53UntrPf.lean`)
* `pfRoot_isOfUnitTrivialType` / `pfRoot_isOfModelType` / `pfRoot_ratFnData`
  (`Prop55Std.lean`)

はすべて **`hftr : ∀ X : C, IsFrobeniusTrivial P X`**(＝ `𝒞` が
Frobenius-trivial 型)を引きずっていた。★これは原文 `Proposition 5.5` の仮定
(Frobenius-isotropic ＋ Frobenius-normalized 型)より**強い**——
`𝔽_Φ` で見れば `n·α = α` すなわち `α = 0` を全対象に要求するからである。

## ★本ファイルが埋める段

**`𝒪^▷` は co-angular で linear な射に沿って単射**(`Proposition 1.11, (iv)`、
在庫 `otriPullHom_injective`)であり、**co-angular pre-step に沿っては全単射**
(`FrobenioidCore` の (iii)(c) `otriFwd` / `otriBwd`)である。

★★したがって `𝒞` が isotropic 型なら(そこでは pre-step は自動的に co-angular ——
在庫 `isCoAngular_of_isotropic`)、**pre-step の対**

  `A₀ ⟵ X ⟶ A`  (`A₀` は Frobenius-trivial、`F.baseSurj` ＋ `F.preStepSpan`)

を `𝒞^pf` へ送るだけで unit-trivial 性が `⟨A₀,1⟩` から `⟨A,1⟩` へ渡る。
★これで **`hftr` は落ちる**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★1. `𝒪^▷` の引き戻しは単射・全単射 -/

section OTriBij

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★★**`𝒪^▷` は co-angular pre-step に沿って全単射**。

★単射性は `Proposition 1.11, (iv)`(`otriPullHom_injective`)、
全射性は `FrobenioidCore` の **(iii)(c)** `otriFwd` / `otriBwd` である。 -/
theorem otriPullHom_bijective (F : FrobenioidCore P) {X A : C} (φ : X ⟶ A)
    (hc : IsCoAngular P φ) (hs : IsPreStep P φ) :
    Function.Bijective (otriPullHom P F φ hc hs.1) := by
  refine ⟨otriPullHom_injective P F φ hc hs.1, fun α => ?_⟩
  obtain ⟨β, hβ, -⟩ := F.otriFwd φ hc hs (α : End X) α.2
  refine ⟨⟨β, hβ.1⟩, ?_⟩
  obtain ⟨_a, -, huniq⟩ := F.otriBwd φ hc hs β hβ.1
  refine Subtype.ext ?_
  refine (huniq _ ⟨(otriPull P F φ hc hs.1 ⟨β, hβ.1⟩).2,
    otriPull_spec P F φ hc hs.1 ⟨β, hβ.1⟩⟩).trans
    (huniq (α : End X) ⟨α.2, hβ.2⟩).symm

/-- ★`𝒪^▷` の中の単元で unit-trivial 性を言い換える。 -/
theorem isUnitTrivial_iff_otri (A : C) :
    IsUnitTrivial P A ↔ ∀ x : OTri P A, IsUnit x → x = 1 := by
  constructor
  · intro h x hx
    have hmem : ((x : OTri P A) : End A) ∈ OTimes P A :=
      ⟨x.2, IsUnit.map (OTri P A).subtype hx⟩
    exact Subtype.ext (Submonoid.mem_bot.mp (h ▸ hmem))
  · intro h
    refine le_antisymm (fun u hu => ?_) bot_le
    rw [Submonoid.mem_bot]
    exact congrArg Subtype.val (h ⟨u, hu.1⟩ (isUnit_otri_of_otimes _ hu))

/-- ★★★**unit-trivial 性は `𝒪^▷` の単射に沿って「上」へ移る** ——
`φ : X ⟶ A` が co-angular かつ linear で `X` が unit-trivial なら `A` も unit-trivial。 -/
theorem isUnitTrivial_of_otriPull_inj (F : FrobenioidCore P) {X A : C} (φ : X ⟶ A)
    (hc : IsCoAngular P φ) (hl : IsLinear P φ) (hX : IsUnitTrivial P X) :
    IsUnitTrivial P A := by
  refine (isUnitTrivial_iff_otri P A).mpr (fun x hx => ?_)
  have h1 : otriPullHom P F φ hc hl x = 1 :=
    (isUnitTrivial_iff_otri P X).mp hX _ (hx.map (otriPullHom P F φ hc hl))
  have h2 : otriPullHom P F φ hc hl x = otriPullHom P F φ hc hl 1 := by
    rw [h1, map_one]
  exact otriPullHom_injective P F φ hc hl h2

/-- ★★★**unit-trivial 性は co-angular pre-step に沿って「下」へ移る** ——
`φ : X ⟶ A` が co-angular pre-step で `A` が unit-trivial なら `X` も unit-trivial。 -/
theorem isUnitTrivial_of_otriPull_surj (F : FrobenioidCore P) {X A : C} (φ : X ⟶ A)
    (hc : IsCoAngular P φ) (hs : IsPreStep P φ) (hA : IsUnitTrivial P A) :
    IsUnitTrivial P X := by
  refine (isUnitTrivial_iff_otri P X).mpr (fun x hx => ?_)
  set E := MulEquiv.ofBijective (otriPullHom P F φ hc hs.1)
    (otriPullHom_bijective P F φ hc hs) with hEdef
  have hyu : IsUnit (E.symm x) := hx.map (E.symm : ↥(OTri P X) →* ↥(OTri P A))
  have hy1 : E.symm x = 1 := (isUnitTrivial_iff_otri P A).mp hA _ hyu
  have hx' : x = E (E.symm x) := (E.apply_symm_apply x).symm
  rw [hx', hy1, map_one]

end OTriBij

/-! ## ★2. `𝒞 ⟶ 𝒞^pf` は pre-step を pre-step に送る -/

section ToPfRoot

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★`𝒞 → 𝒞^pf` は pre-step を pre-step に送る(次数と底の同型は保たれる)。 -/
theorem toPfRoot_isPreStep {A B : C} (φ : A ⟶ B) (hφ : IsPreStep P φ) :
    IsPreStep (pfRootPre P F) ((toPfRoot P F).map φ) := by
  refine ⟨?_, ?_⟩
  · show (pfRootPre P F).degFr ((toPfRoot P F).map φ) = 1
    rw [toPfRoot_degFr]
    exact hφ.1
  · show IsIso (rootBase (toRootHom (F := F) φ))
    rw [rootBase_toRootHom]
    exact hφ.2

end ToPfRoot

/-! ## ★3. 本題 —— `hftr` を落とす -/

section Arb

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★★★★★★★**[FrdI] Proposition 5.5, (i) の「任意の `A`」の段** ——
`𝒞` が isotropic 型・Frobenius-normalized 型・unit-trivial 型なら、
**`𝒞^pf` はすべての `⟨A,1⟩` で unit-trivial** である。

原文 (FrdI p.104):
> the case of arbitrary A then follows by considering "pairs of pre-steps" as in Theorem

★`F.baseSurj` で Frobenius-trivial な `A₀`(底は `A` と同型)を取り、
`F.preStepSpan` で **pre-step の対** `A₀ ⟵ X ⟶ A` を得る。
★`𝒞^pf` は isotropic なので、そこでは pre-step は自動的に co-angular であり、
`𝒪^▷` は `X ⟶ A₀` に沿って全単射、`X ⟶ A` に沿って単射である。
★★これで unit-trivial 性が `⟨A₀,1⟩ → ⟨X,1⟩ → ⟨A,1⟩` と渡る。 -/
theorem otimes_pfRoot_eq_bot_of_arb (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfn : ∀ X : C, IsFrobeniusNormalized P X)
    (hut : ∀ X : C, IsUnitTrivial P X)
    (Gpf : Frobenioid (pfRootPre P F)) (A : C) :
    IsUnitTrivial (pfRootPre P F) ((toPfRoot P F).obj A) := by
  obtain ⟨A₀, hA₀, ⟨e⟩⟩ := F.baseSurj ((P.toElem.obj A).base)
  obtain ⟨ζ, hdeg, hprop⟩ := hA₀
  obtain ⟨X, φ, ψ, hφ, hψ, -⟩ := F.preStepSpan A₀ A e.hom (by infer_instance)
  -- `𝒞^pf` へ送る
  have hφ' : IsPreStep (pfRootPre P F) ((toPfRoot P F).map φ) := toPfRoot_isPreStep φ hφ
  have hψ' : IsPreStep (pfRootPre P F) ((toPfRoot P F).map ψ) := toPfRoot_isPreStep ψ hψ
  have hcφ : IsCoAngular (pfRootPre P F) ((toPfRoot P F).map φ) :=
    isCoAngular_of_isotropic (pfRootPre P F) Gpf
      (pfRoot_isOfIsotropicType (F := F) hfi _) _ hφ'
  have hcψ : IsCoAngular (pfRootPre P F) ((toPfRoot P F).map ψ) :=
    isCoAngular_of_isotropic (pfRootPre P F) Gpf
      (pfRoot_isOfIsotropicType (F := F) hfi _) _ hψ'
  -- `⟨A₀,1⟩` は Frobenius-trivial な `A₀` から来るので unit-trivial
  have hA₀pf : IsUnitTrivial (pfRootPre P F) ((toPfRoot P F).obj A₀) :=
    otimes_pfRoot_eq_bot (F := F) hiso ⟨ζ, hdeg, hprop⟩ (hfn A₀) (hfn (rtObj P F A₀ 1))
      ζ hdeg hprop (hut A₀)
  -- `⟨X,1⟩` へ降ろし、`⟨A,1⟩` へ上げる
  have hXpf : IsUnitTrivial (pfRootPre P F) ((toPfRoot P F).obj X) :=
    isUnitTrivial_of_otriPull_surj (pfRootPre P F) Gpf.core _ hcφ hφ' hA₀pf
  exact isUnitTrivial_of_otriPull_inj (pfRootPre P F) Gpf.core _ hcψ hψ'.1 hXpf

/-- ★★★★★★**[FrdI] Proposition 5.3 (a) / 5.5, (iii)** ——
**完全化は unit-trivial 性を保つ**(★`hftr` なしの版)。 -/
theorem pfRoot_isOfUnitTrivialType_of_arb (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfn : ∀ X : C, IsFrobeniusNormalized P X)
    (hut : ∀ X : C, IsUnitTrivial P X)
    (Gpf : Frobenioid (pfRootPre P F)) :
    IsOfUnitTrivialType (pfRootPre P F) := by
  intro Z
  obtain ⟨A, n⟩ := Z
  refine (isUnitTrivial_pfRoot_iff A n).mpr ?_
  refine isUnitTrivial_pfCat_of_rtObj_one (F := F) (rtObj P F A (n * n)) ?_
  exact (isUnitTrivial_pfRoot_iff (rtObj P F A (n * n)) 1).mp
    (otimes_pfRoot_eq_bot_of_arb hfi hiso hfn hut Gpf (rtObj P F A (n * n)))

variable [IsConnected D]

/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** —— **`𝒞^pf` は model 型**
(★`hftr` なしの版)。 -/
theorem pfRoot_isOfModelType_of_arb (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfn : ∀ X : C, IsFrobeniusNormalized P X)
    (hut : ∀ X : C, IsUnitTrivial P X)
    (Gpf : Frobenioid (pfRootPre P F)) :
    IsOfModelType (PfRootObj P F) (pfRootPre P F) Gpf :=
  thm_5_1_iv Gpf (pfRoot_isOfIsotropicType (F := F) hfi)
    (pfRoot_isOfUnitTrivialType_of_arb hfi hiso hfn hut Gpf)

end Arb

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Proposition 5.5, (i)` の
「任意の `A` は pre-step の対で Frobenius-trivial な場合に帰着する」の段。 -/
def otimes_pfRoot_eq_bot_of_arb.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — 任意の A は pre-step の対で Frobenius-trivial に帰着",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★locator —— `hftr` なしの「完全化は unit-trivial 性を保つ」。 -/
def pfRoot_isOfUnitTrivialType_of_arb.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (iii) — 𝒞^pf は unit-trivial 型(Frobenius-trivial 型の仮定なし)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
