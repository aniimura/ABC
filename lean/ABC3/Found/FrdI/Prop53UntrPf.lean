/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfUnit

/-!
# `⟨A,n⟩` の End は `𝒞^pf` の End の共役(鎖 `prop55` の `p53-untrpf` の (a))

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.103。

原文 (FrdI p.103):
> the Frobenioid C the model Frobenioid [cf. Theorem 5.2, (ii)] associated to

## ★★`compRoot` は `compPf` の共役である

`Prop55PfUnit.lean` の `otimes_pfRoot_eq_bot` は `⟨A,1⟩` についてしか言っていない。
★すべての `⟨A,n⟩` へ広げるには、**`compRoot` の定義を開く**のが要点である。

`compRoot_eq_lift`(`Def31Pf.lean`)を `c = 1`、`PA = PB = PE = n·n`、
`ef = eg = er = n` で使うと、`X = Y = Z = ⟨A,n⟩` のとき

  `compRoot f g = Θ.hom (compPf (Θ.inv f) (Θ.inv g))`,
  `Θ := rtRootIso A A (n·n = n·n) (n·n = n·n)`

★★つまり **`Θ.inv` は単系の同型 `End ⟨A,n⟩ ≅ Hom^pf(A^{(n·n)}, A^{(n·n)})`** である。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `rootSelfIso` | ★`Θ`(根 `n·n` から根 `n` へ) |
| `rootSelfIso_inv_compRoot` | ★★**`Θ.inv` は合成を `compPf` へ移す** |
| `rootSelfIso_inv_id` | `Θ.inv` は単位を単位へ |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u v2 u2 w

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★**`Θ`** —— 根 `n·n` の `Hom^pf` から根 `n` の `Hom^pf` への同型。 -/
noncomputable def rootSelfIso (A : C) (n : ℕ+) :
    HomPf P F (rtObj P F A (n * n)) (rtObj P F A (n * n))
      ≅ HomPf P F (rtObj P F A n) (rtObj P F A n) :=
  rtRootIso P F A A (show n * n = n * n from rfl) (show n * n = n * n from rfl)

/-- ★★★**`Θ.inv` は `compRoot` を `compPf` へ移す**。 -/
theorem rootSelfIso_inv_compRoot (A : C) (n : ℕ+)
    (f g : HomRoot P F (⟨A, n⟩ : PfRootObj P F) ⟨A, n⟩) :
    (rootSelfIso (F := F) A n).inv (compRoot P F f g)
      = compPf P F ((rootSelfIso (F := F) A n).inv f) ((rootSelfIso (F := F) A n).inv g) := by
  have hlift := compRoot_eq_lift (P := P) (F := F) (X := (⟨A, n⟩ : PfRootObj P F))
    (Y := ⟨A, n⟩) (Z := ⟨A, n⟩) f g
    (c := 1) (PA := n * n) (PB := n * n) (PE := n * n)
    (hcA := (one_mul _).symm) (hcB := (one_mul _).symm) (hcE := (one_mul _).symm)
    (ef := n) (eg := n) (er := n)
    (hfA := rfl) (hfB := rfl) (hgA := rfl) (hgE := rfl) (hrA := rfl) (hrE := rfl)
  rw [hlift]
  exact Iso.hom_inv_id_apply _ _

/-- ★**`Θ.inv` は単位を単位へ**。 -/
theorem rootSelfIso_inv_id (A : C) (n : ℕ+) :
    (rootSelfIso (F := F) A n).inv (idRoot P F (⟨A, n⟩ : PfRootObj P F))
      = toHomPf (F := F) (𝟙 (rtObj P F A (n * n))) :=
  rtRootIso_inv_id P F A (show n * n = n * n from rfl)

/-! ## ★2. `OTri` の対応 -/

/-- ★**`rootBase` が恒等 ⟺ `pfBase` が恒等**(共役だから)。 -/
theorem rootBase_eq_id_iff (A : C) (n : ℕ+)
    (f : HomRoot P F (⟨A, n⟩ : PfRootObj P F) ⟨A, n⟩) :
    rootBase f = 𝟙 _ ↔ pfBase f = 𝟙 _ := by
  haveI : IsIso (P.Base (rtExt P F A n)) := (rtExt_frobType P F A n).2
  have hspec := rootBase_spec (P := P) (F := F) f
  constructor
  · intro h
    rw [h, Category.id_comp] at hspec
    exact (cancel_epi (P.Base (rtExt P F A n))).mp (by rw [Category.comp_id]; exact hspec.symm)
  · intro h
    rw [h, Category.comp_id] at hspec
    exact (cancel_mono (P.Base (rtExt P F A n))).mp (by rw [Category.id_comp]; exact hspec)

/-- ★**`Θ := (rootSelfIso).inv` は底の恒等性を保つ**。 -/
theorem pfBase_rootSelfIso_inv_eq_id_iff (A : C) (n : ℕ+)
    (f : HomRoot P F (⟨A, n⟩ : PfRootObj P F) ⟨A, n⟩) :
    pfBase ((rootSelfIso (F := F) A n).inv f) = 𝟙 _ ↔ pfBase f = 𝟙 _ := by
  haveI : IsIso (P.Base (rtLift P F A (show n * n = n * n from rfl))) :=
    (rtLift_frobType P F A _).2
  have h : pfBase f ≫ P.Base (rtLift P F A (show n * n = n * n from rfl))
      = P.Base (rtLift P F A (show n * n = n * n from rfl))
        ≫ pfBase ((rootSelfIso (F := F) A n).inv f) :=
    pfBase_rtRootIso_inv A A _ _ f
  constructor
  · intro hz
    rw [hz, Category.comp_id] at h
    exact (cancel_mono (P.Base (rtLift P F A (show n * n = n * n from rfl)))).mp
      (by rw [Category.id_comp]; exact h)
  · intro hz
    rw [hz, Category.id_comp] at h
    exact (cancel_epi (P.Base (rtLift P F A (show n * n = n * n from rfl)))).mp
      (by rw [Category.comp_id]; exact h.symm)

/-- ★`Θ` の行き先の対象 —— **`𝒞^pf` の対象として見た `A^{(n·n)}`**。

★`PfCat P F` は定義上 `C` そのものだが、`𝟙` や `End` の解決先を分けるために
**型を `PfCat P F` と書いた別名**が要る(`pfDown` と同じ事情)。 -/
noncomputable def rtObjPf (P : PreFrobenioid C Φ) (F : FrobenioidCore P) (A : C) (n : ℕ+) : PfCat P F :=
  (rtObj P F A (n * n) : C)

/-- ★`𝒞^pf` 側の `IsBaseIdentity` を `pfBase = 𝟙` に開く。 -/
theorem isBaseIdentity_pf_iff {X : PfCat P F} (g : End X) :
    IsBaseIdentity (pfPre P F) g ↔ pfBase g = 𝟙 _ :=
  ⟨fun h => h.trans ((pfPre P F).Base_id X), fun h => h.trans ((pfPre P F).Base_id X).symm⟩

/-- ★根つき側の `IsBaseIdentity` を `rootBase = 𝟙` に開く。 -/
theorem isBaseIdentity_pfRoot_iff {X : PfRootObj P F} (g : End X) :
    IsBaseIdentity (pfRootPre P F) g ↔ rootBase g = 𝟙 _ :=
  ⟨fun h => h.trans ((pfRootPre P F).Base_id X),
    fun h => h.trans ((pfRootPre P F).Base_id X).symm⟩

/-- ★★★**`OTri` は `Θ` で対応する**。 -/
theorem mem_otri_rootSelfIso_inv (A : C) (n : ℕ+)
    (f : HomRoot P F (⟨A, n⟩ : PfRootObj P F) ⟨A, n⟩) :
    f ∈ OTri (pfRootPre P F) (⟨A, n⟩ : PfRootObj P F)
      ↔ (rootSelfIso (F := F) A n).inv f ∈ OTri (pfPre P F) (rtObjPf P F A n) := by
  have hdeg : pfDeg ((rootSelfIso (F := F) A n).inv f) = pfDeg f :=
    pfDeg_rtRootIso_inv A A _ _ f
  constructor
  · rintro ⟨hbi, hl⟩
    refine ⟨(isBaseIdentity_pf_iff (X := rtObjPf P F A n) _).mpr ?_, ?_⟩
    · exact (pfBase_rootSelfIso_inv_eq_id_iff A n f).mpr
        ((rootBase_eq_id_iff A n f).mp ((isBaseIdentity_pfRoot_iff f).mp hbi))
    · show pfDeg ((rootSelfIso (F := F) A n).inv f) = 1
      rw [hdeg]
      exact hl
  · rintro ⟨hbi, hl⟩
    refine ⟨(isBaseIdentity_pfRoot_iff f).mpr ?_, ?_⟩
    · exact (rootBase_eq_id_iff A n f).mpr
        ((pfBase_rootSelfIso_inv_eq_id_iff A n f).mp ((isBaseIdentity_pf_iff (X := rtObjPf P F A n) _).mp hbi))
    · show pfDeg f = 1
      rw [← hdeg]
      exact hl

/-! ## ★3. `OTimes` の対応と unit-trivial 性の移送 -/

/-- ★★**`Θ` は単系の同型** —— `End ⟨A,n⟩ ≃* End_{𝒞^pf} A^{(n·n)}`。

★`End` の積は `x * y = y ≫ x` なので、`rootSelfIso_inv_compRoot` がそのまま `map_mul'`。 -/
noncomputable def endRootMulEquiv (A : C) (n : ℕ+) :
    End (⟨A, n⟩ : PfRootObj P F) ≃* End (rtObjPf P F A n) where
  toFun f := (rootSelfIso (F := F) A n).inv f
  invFun g := (rootSelfIso (F := F) A n).hom g
  left_inv f := (rootSelfIso (F := F) A n).inv_hom_id_apply f
  right_inv g := (rootSelfIso (F := F) A n).hom_inv_id_apply g
  map_mul' x y := rootSelfIso_inv_compRoot A n y x

@[simp] theorem endRootMulEquiv_apply (A : C) (n : ℕ+) (f : End (⟨A, n⟩ : PfRootObj P F)) :
    endRootMulEquiv (F := F) A n f = (rootSelfIso (F := F) A n).inv f := rfl

/-- ★★★**`OTimes` は `Θ` で対応する**。 -/
theorem mem_otimes_rootSelfIso_inv (A : C) (n : ℕ+) (f : End (⟨A, n⟩ : PfRootObj P F)) :
    f ∈ OTimes (pfRootPre P F) (⟨A, n⟩ : PfRootObj P F)
      ↔ endRootMulEquiv (F := F) A n f ∈ OTimes (pfPre P F) (rtObjPf P F A n) := by
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨(mem_otri_rootSelfIso_inv A n f).mp h1,
      h2.map (endRootMulEquiv (F := F) A n).toMonoidHom⟩
  · rintro ⟨h1, h2⟩
    refine ⟨(mem_otri_rootSelfIso_inv A n f).mpr h1, ?_⟩
    have hu := h2.map (endRootMulEquiv (F := F) A n).symm.toMonoidHom
    rwa [MulEquiv.coe_toMonoidHom, MulEquiv.symm_apply_apply] at hu

/-- ★★★★**根 `n` の unit-trivial 性は `𝒞^pf` の `A^{(n·n)}` の unit-trivial 性と同値**。 -/
theorem isUnitTrivial_pfRoot_iff (A : C) (n : ℕ+) :
    IsUnitTrivial (pfRootPre P F) (⟨A, n⟩ : PfRootObj P F)
      ↔ IsUnitTrivial (pfPre P F) (rtObjPf P F A n) := by
  constructor
  · intro h
    refine le_antisymm (fun g hg => ?_) bot_le
    rw [Submonoid.mem_bot]
    set f : End (⟨A, n⟩ : PfRootObj P F) := (endRootMulEquiv (F := F) A n).symm g with hf
    have hgf : endRootMulEquiv (F := F) A n f = g := (endRootMulEquiv (F := F) A n).apply_symm_apply g
    have hmem : f ∈ OTimes (pfRootPre P F) (⟨A, n⟩ : PfRootObj P F) :=
      (mem_otimes_rootSelfIso_inv A n f).mpr (hgf ▸ hg)
    have hf1 : f = 1 := Submonoid.mem_bot.mp (h ▸ hmem)
    rw [← hgf, hf1, map_one]
  · intro h
    refine le_antisymm (fun f hf => ?_) bot_le
    rw [Submonoid.mem_bot]
    have hΘ : endRootMulEquiv (F := F) A n f ∈ OTimes (pfPre P F) (rtObjPf P F A n) :=
      (mem_otimes_rootSelfIso_inv A n f).mp hf
    have h1 : endRootMulEquiv (F := F) A n f = 1 := Submonoid.mem_bot.mp (h ▸ hΘ)
    have h2 := congrArg (endRootMulEquiv (F := F) A n).symm h1
    rwa [MulEquiv.symm_apply_apply, map_one] at h2

/-! ## ★4. unit-trivial 性は次数 1 の同型で移る -/

universe v3 u3 w3

section Transport

variable {Φ' : MonoidOn.{v, u, w3} D} {C' : Type u3} [Category.{v3} C']
  {Q : PreFrobenioid C' Φ'}

/-- ★★★**unit-trivial 性は次数 `1` の同型に沿って移る**。

★共役 `φ ↦ u ≫ φ ≫ u⁻¹`(`Iso.conj`)は `End` の単系同型であり、
`Base`・`degFr` の 3 本の公理で `𝒪^▷`・`𝒪^×` を保つ。 -/
theorem isUnitTrivial_of_iso {A B : C'} (u : A ⟶ B) (hiso : IsIso u) (hd : Q.degFr u = 1)
    (h : IsUnitTrivial Q A) : IsUnitTrivial Q B := by
  haveI := hiso
  have hdinv : Q.degFr (inv u) = 1 := degFr_inv_eq_one (P := Q) u hd
  refine le_antisymm (fun φ hφ => ?_) bot_le
  rw [Submonoid.mem_bot]
  obtain ⟨⟨hbi, hl⟩, hu⟩ := hφ
  have hcv : (asIso u).symm.conj φ = u ≫ (show B ⟶ B from φ) ≫ inv u := Iso.conj_apply _ _
  have hmem : (asIso u).symm.conj φ ∈ OTimes Q A := by
    refine ⟨⟨?_, ?_⟩, hu.map (asIso u).symm.conj.toMonoidHom⟩
    · show Q.Base ((asIso u).symm.conj φ) = Q.Base (𝟙 A)
      rw [hcv, Q.Base_comp, Q.Base_comp,
        show Q.Base (show B ⟶ B from φ) = 𝟙 _ from hbi.trans (Q.Base_id B),
        Category.id_comp, ← Q.Base_comp, IsIso.hom_inv_id, Q.Base_id]
    · show Q.degFr ((asIso u).symm.conj φ) = 1
      rw [hcv, Q.degFr_comp, Q.degFr_comp, hd, hdinv,
        show Q.degFr (show B ⟶ B from φ) = 1 from hl, one_mul, one_mul]
  have h1 : (asIso u).symm.conj φ = 1 := Submonoid.mem_bot.mp (h ▸ hmem)
  have h2 := congrArg (asIso u).symm.conj.symm h1
  rwa [MulEquiv.symm_apply_apply, map_one] at h2

end Transport

/-! ## ★5. `Proposition 5.3` (a) —— 根 `n` の unit-trivial 性 -/

/-- ★★`𝒞^pf` の中で `A^{(1)} ≅ A`(次数 `1` の同型)。 -/
theorem isUnitTrivial_pfCat_of_rtObj_one (B : C)
    (h : IsUnitTrivial (pfPre P F) ((rtObj P F B 1 : C) : PfCat P F)) :
    IsUnitTrivial (pfPre P F) ((B : C) : PfCat P F) := by
  haveI := isIso_rtExt_one P F B
  haveI hw : IsIso ((toPfCat P F).map (inv (rtExt P F B 1))) :=
    ⟨(toPfCat P F).map (rtExt P F B 1),
      by rw [← Functor.map_comp, IsIso.inv_hom_id]; exact (toPfCat P F).map_id _,
      by rw [← Functor.map_comp, IsIso.hom_inv_id]; exact (toPfCat P F).map_id _⟩
  have hdeg : (pfPre P F).degFr ((toPfCat P F).map (inv (rtExt P F B 1))) = 1 := by
    show pfDeg (toHomPf (F := F) (inv (rtExt P F B 1))) = 1
    rw [pfDeg_toHomPf]
    exact degFr_inv_eq_one (rtExt P F B 1) (rtExt_degFr P F B 1)
  exact isUnitTrivial_of_iso (Q := pfPre P F) ((toPfCat P F).map (inv (rtExt P F B 1))) hw hdeg h

/-- ★★★★★**[FrdI] Proposition 5.3** (a) —— **完全化は unit-trivial 性を保つ**
(★原文が `⟨A,1⟩` でなく**すべての根 `⟨A,n⟩`** について要求する形)。

★★`otimes_pfRoot_eq_bot` は `⟨A,1⟩` しか言っていなかった。
`Θ : End ⟨A,n⟩ ≃* End_{𝒞^pf} A^{(n·n)}`(`endRootMulEquiv`)で `n` を `1` に落とす。

原文 (FrdI p.103):
> tion monoid Φbirat (respectively, Q · Φbirat = Φbirat ⊗Z Q = (Φbirat)pf). In particular, -/
theorem otimes_pfRoot_eq_bot_of_root (hiso : ∀ X : C, IsIsotropic P X)
    (hftr : ∀ X : C, IsFrobeniusTrivial P X)
    (hfn : ∀ X : C, IsFrobeniusNormalized P X)
    (ζ : ∀ X : C, ℕ+ →* End X)
    (hdeg : ∀ (X : C) (m : ℕ+), P.degFr ((ζ X m : End X) : X ⟶ X) = m)
    (hprop : ∀ (X : C) (m : ℕ+),
      IsBaseIdentity P (ζ X m) ∧ IsFrobeniusType P ((ζ X m : End X) : X ⟶ X))
    (hut : ∀ X : C, IsUnitTrivial P X) (A : C) (n : ℕ+) :
    IsUnitTrivial (pfRootPre P F) (⟨A, n⟩ : PfRootObj P F) := by
  refine (isUnitTrivial_pfRoot_iff A n).mpr ?_
  refine isUnitTrivial_pfCat_of_rtObj_one (F := F) (rtObj P F A (n * n)) ?_
  exact (isUnitTrivial_pfRoot_iff (rtObj P F A (n * n)) 1).mp
    (otimes_pfRoot_eq_bot hiso (hftr _) (hfn _) (hfn _) (ζ _) (hdeg _) (hprop _) (hut _))

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Proposition 5.3` の「`(𝒞^un-tr)^pf` は unit-trivial」をすべての
`⟨A,n⟩` へ広げる段(第 1 段)。 -/
def rootSelfIso_inv_compRoot.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — ⟨A,n⟩ の End は 𝒞^pf の End の共役",
    sectionId := "frdi-prop-5-3" }

/-- ★★★★★locator —— `Proposition 5.3` の (a):
**すべての根 `⟨A,n⟩` で `(𝒞^un-tr)^pf` は unit-trivial**。 -/
def otimes_pfRoot_eq_bot_of_root.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — (𝒞^un-tr)^pf も unit-trivial(すべての根 ⟨A,n⟩)",
    sectionId := "frdi-prop-5-3" }

/-- ★★locator —— unit-trivial 性は次数 1 の同型で移る(`Definition 1.2, (iv)` の一般論)。 -/
def isUnitTrivial_of_iso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 23,
    item := "Definition 1.2, (iv) — unit-trivial 性は次数 1 の同型で移る",
    sectionId := "frdi-def-1-2-iv" }

end ABC3.Found.FrdI
