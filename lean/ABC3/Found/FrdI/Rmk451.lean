/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def45
import ABC3.Found.FrdI.Def31
import ABC3.Found.FrdI.Remark311

/-!
# [FrdI] Remark 4.5.1 —— `𝒞^istr` は (rationally) standard 型を受け継ぐ

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.86。

原文 (FrdI p.86):
> We observe in passing that it is immediate from the definitions

原文 (FrdI p.86):
> that if C is of rationally standard type (respectively, of standard type), then so is

## ★原文は「immediate from the definitions」と書く

★★**実際には 5 条(standard 型)＋ 4 条(rationally standard 型)の移送**である。
本ファイルは `Definition 3.1, (i)` の 5 条を 1 つずつ移す。

| 条 | 内容 | 本ファイル |
|---|---|---|
| (a) | quasi-isotropic 型 | `istr_quasiIsotropicType` |
| (a) | Frobenius-isotropic 型 | `istr_frobIsotropicType` |
| (b) | group-like 型なら `𝒞^istr` に Frobenius-compact 対象 | ★残り |
| (c) | Frobenius-normalized 型 | `istr_frobNormalizedType` |
| (d) | `𝒟` が FSMFF-type | ★**そのまま**(`𝒟` は変わらない) |
| (e) | `Φ` が non-dilating | ★**そのまま**(`Φ` は変わらない) |

★(a) の 2 条は **`𝒞^istr` の全対象が isotropic** であること(`istr_isotropic`)から出る。
★(c) は `Istr P` が `𝒞` の充満部分圏で `istrPre` が `P` を制限したものだから、
`.hom` に落として `𝒞` の主張を当てるだけ —— ただし
**冪と `.hom` の交換**(`Functor.mapEnd` の `map_pow`)が要る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-! ## ★1. (a) —— isotropic 由来の 2 条 -/

/-- ★`𝒞^istr` は isotropic 型。 -/
theorem istr_isOfIsotropicType : IsOfIsotropicType (istrPre P F) :=
  fun A => istr_isotropic P F A

/-- ★★**(a) その 1** —— `𝒞^istr` は quasi-isotropic 型。 -/
theorem istr_quasiIsotropicType : IsOfQuasiIsotropicType (Istr P) (istrPre P F) :=
  isOfQuasiIsotropicType_of_isOfIsotropicType (istrPre P F) (istr_frobenioidCore P F)
    istr_isOfIsotropicType

/-- ★★**(a) その 2** —— `𝒞^istr` は Frobenius-isotropic 型。

★`𝟙` が Frobenius 型で終域が isotropic なので、条件は自明に満たされる。 -/
theorem istr_frobIsotropicType : IsOfFrobeniusIsotropicType (istrPre P F) := by
  intro A
  exact ⟨A, 𝟙 A, ⟨⟨isCoAngular_id (istrPre P F) A, isIsometric_of_isIso (istrPre P F) (𝟙 A)⟩,
    isBaseIsomorphism_of_isIso (istrPre P F) (𝟙 A)⟩, istr_isotropic P F A⟩

/-! ## ★2. (c) —— Frobenius-normalized 型の移送 -/

/-- ★`𝒞^istr` の `𝒪^▷` は `𝒞` の `𝒪^▷` に落ちる。 -/
theorem istr_mem_otri {A : Istr P} {α : End A} (h : α ∈ OTri (istrPre P F) A) :
    ((α : A ⟶ A).hom : A.obj ⟶ A.obj) ∈ OTri P A.obj := by
  refine ⟨?_, ?_⟩
  · show P.Base ((α : A ⟶ A).hom) = P.Base (𝟙 A.obj)
    have h1 : (istrPre P F).Base (α : A ⟶ A) = (istrPre P F).Base (𝟙 A) := h.1
    rw [istrPre_Base, istrPre_Base] at h1
    exact h1
  · show P.degFr ((α : A ⟶ A).hom) = 1
    have h2 : (istrPre P F).degFr (α : A ⟶ A) = 1 := h.2
    rw [istrPre_degFr] at h2
    exact h2

/-- ★★**(c)** —— Frobenius-normalized 型は `𝒞^istr` に移る。

★★**要点**: `Istr P` は充満部分圏なので合成は `𝒞` のものと同じだが、
`α ^ n` を `.hom` に落とすには**冪と `.hom` の交換**が要る。
`(isotropicProp P).ι.mapEnd A : End A →* End A.obj` の `map_pow` がそれである。 -/
theorem istr_frobNormalizedType (h : IsOfFrobeniusNormalizedType P) :
    IsOfFrobeniusNormalizedType (istrPre P F) := by
  intro A φ hφ α hα
  refine ObjectProperty.hom_ext (isotropicProp P) ?_
  have hb : IsBaseIdentity P ((φ : A ⟶ A).hom) := by
    show P.Base ((φ : A ⟶ A).hom) = P.Base (𝟙 A.obj)
    have h1 : (istrPre P F).Base (φ : A ⟶ A) = (istrPre P F).Base (𝟙 A) := hφ
    rw [istrPre_Base, istrPre_Base] at h1
    exact h1
  have hd : (istrPre P F).degFr φ = P.degFr ((φ : A ⟶ A).hom) := istrPre_degFr F φ
  have hcore := h A.obj ((φ : A ⟶ A).hom) hb _ (istr_mem_otri (F := F) hα)
  have hpow := map_pow ((isotropicProp P).ι.mapEnd A) α
    ((P.degFr ((φ : A ⟶ A).hom) : ℕ+) : ℕ)
  rw [hd]
  exact Eq.trans (congrArg (fun t : A.obj ⟶ A.obj => (φ : A ⟶ A).hom ≫ t) hpow) hcore

/-! ## ★3. (b) の第 1 段 —— group-like 型の降下 -/

/-- ★`IsGroupLike` は全射準同型に沿って移る。 -/
theorem isGroupLike_of_surjective {M N : Type w} [AddCommMonoid M] [AddCommMonoid N]
    (f : M →+ N) (hf : Function.Surjective f) (h : IsGroupLike M) : IsGroupLike N := by
  refine (isGroupLike_iff N).mpr ?_
  intro b
  obtain ⟨a, rfl⟩ := hf b
  exact ((isGroupLike_iff M).mp h a).map f

/-- ★★**`𝒞^istr` が group-like 型なら `𝒞` も group-like 型**。

★isotropic hull `A ⟶ A^istr` は pre-step、したがって**底同型**なので
`Φ(base A^istr) → Φ(base A)` は全単射(`MonoidOn.map_bijective_of_iso`)。
group-like 性はその全射に沿って移る。 -/
theorem istr_groupLikeType_down (h : IsOfGroupLikeType (istrPre P F)) :
    IsOfGroupLikeType P := by
  intro A
  obtain ⟨B, φ, _, hstep, hisoB, _⟩ := F.isotropicHullExists A
  haveI hb : IsIso (P.Base φ) := hstep.2
  exact isGroupLike_of_surjective (Φ.map (P.Base φ))
    (MonoidOn.map_bijective_of_iso Φ (asIso (P.Base φ))).2 (h ⟨B, hisoB⟩)

/-! ## ★4. (b) の第 2 段 —— Frobenius-compact 対象の移送

★`Istr (istrPre P F)` は `Istr P` の**全対象**からなる充満部分圏
(`istr_isotropic` により条件が空虚に真)なので、`End` も `OTimes` も一致する。
★`IsFrobeniusCompact` は `End` と `OTimes` と自己同型だけで書かれているので、
その 3 つの対応を作れば移送できる。 -/

/-- ★二重 `Istr` の対象(全対象が isotropic なので誰でも入る)。 -/
abbrev istr2 (A : Istr P) : Istr (istrPre P F) := ⟨A, istr_isotropic P F A⟩

/-- ★二重 `Istr` の自己射モノイドは `Istr P` のそれと同型。 -/
def istr2EndEquiv (A : Istr P) : End (istr2 (F := F) A) ≃* End A where
  toFun f := f.hom
  invFun g := ⟨g⟩
  left_inv _ := rfl
  right_inv _ := rfl
  map_mul' _ _ := rfl

/-- ★`OTimes` は二重 `Istr` と行き来する。 -/
theorem istr2_mem_otimes_iff (A : Istr P) (f : End (istr2 (F := F) A)) :
    f ∈ OTimes (istrPre (istrPre P F) (istr_frobenioidCore P F)) (istr2 (F := F) A)
      ↔ (istr2EndEquiv (F := F) A) f ∈ OTimes (istrPre P F) A := by
  constructor
  · rintro ⟨⟨hb, hl⟩, hu⟩
    refine ⟨⟨?_, ?_⟩, hu.map (istr2EndEquiv (F := F) A).toMonoidHom⟩
    · show (istrPre P F).Base (f.hom) = (istrPre P F).Base (𝟙 A)
      have h1 : (istrPre (istrPre P F) (istr_frobenioidCore P F)).Base f
          = (istrPre (istrPre P F) (istr_frobenioidCore P F)).Base (𝟙 (istr2 (F := F) A)) := hb
      rw [istrPre_Base _ f, istrPre_Base _ (𝟙 (istr2 (F := F) A))] at h1
      exact h1
    · show (istrPre P F).degFr (f.hom) = 1
      have h2 : (istrPre (istrPre P F) (istr_frobenioidCore P F)).degFr f = 1 := hl
      rw [istrPre_degFr _ f] at h2
      exact h2
  · rintro ⟨⟨hb, hl⟩, hu⟩
    refine ⟨⟨?_, ?_⟩, hu.map (istr2EndEquiv (F := F) A).symm.toMonoidHom⟩
    · show (istrPre (istrPre P F) (istr_frobenioidCore P F)).Base f
        = (istrPre (istrPre P F) (istr_frobenioidCore P F)).Base (𝟙 (istr2 (F := F) A))
      rw [istrPre_Base _ f, istrPre_Base _ (𝟙 (istr2 (F := F) A))]
      exact hb
    · show (istrPre (istrPre P F) (istr_frobenioidCore P F)).degFr f = 1
      rw [istrPre_degFr _ f]
      exact hl

/-- ★二重 `Istr` の自己同型は `Istr P` の自己同型と 1 対 1。 -/
def istr2IsoEquiv (A : Istr P) : ((istr2 (F := F) A) ≅ (istr2 (F := F) A)) ≃ (A ≅ A) where
  toFun θ := Iso.mk θ.hom.hom θ.inv.hom
    (congrArg InducedCategory.Hom.hom θ.hom_inv_id)
    (congrArg InducedCategory.Hom.hom θ.inv_hom_id)
  invFun θ := Iso.mk ⟨θ.hom⟩ ⟨θ.inv⟩
    (InducedCategory.Hom.ext θ.hom_inv_id)
    (InducedCategory.Hom.ext θ.inv_hom_id)
  left_inv _ := rfl
  right_inv _ := rfl

/-- ★共役は対応と交換する。 -/
theorem istr2_endConj (A : Istr P) (θ : (istr2 (F := F) A) ≅ (istr2 (F := F) A))
    (u : End (istr2 (F := F) A)) :
    (istr2EndEquiv (F := F) A) (endConj θ u)
      = endConj (istr2IsoEquiv (F := F) A θ) ((istr2EndEquiv (F := F) A) u) := rfl

/-- ★★**(b) 第 2 段** —— Frobenius-compact 対象は二重 `Istr` へ移る。 -/
theorem istr2_frobeniusCompact (A : Istr P) (h : IsFrobeniusCompact (istrPre P F) A) :
    IsFrobeniusCompact (istrPre (istrPre P F) (istr_frobenioidCore P F)) (istr2 (F := F) A) := by
  obtain ⟨hcomm, ⟨u₀, hu₀mem, hu₀pow⟩, hconj⟩ := h
  set e := istr2EndEquiv (F := F) A with he
  refine ⟨?_, ⟨e.symm u₀, ?_, ?_⟩, ?_⟩
  · intro x y hx hy
    refine e.injective ?_
    rw [map_mul, map_mul]
    exact hcomm _ _ ((istr2_mem_otimes_iff A x).mp hx) ((istr2_mem_otimes_iff A y).mp hy)
  · exact (istr2_mem_otimes_iff A (e.symm u₀)).mpr
      (by rw [MulEquiv.apply_symm_apply]; exact hu₀mem)
  · intro k hk
    refine hu₀pow k ?_
    have h1 := congrArg e hk
    rw [map_pow, MulEquiv.apply_symm_apply, map_one] at h1
    exact h1
  · intro θ c d H u hu
    have H' : ∀ v : End A, v ∈ OTimes (istrPre P F) A → ∃ k : ℕ+,
        ((endConj (istr2IsoEquiv (F := F) A θ) v) ^ (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)
          = (v ^ (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A) := by
      intro v hv
      obtain ⟨k, hk⟩ := H (e.symm v)
        ((istr2_mem_otimes_iff A (e.symm v)).mpr (by rw [MulEquiv.apply_symm_apply]; exact hv))
      refine ⟨k, ?_⟩
      have h2 := congrArg e hk
      rw [map_pow, map_pow, istr2_endConj, MulEquiv.apply_symm_apply] at h2
      exact h2
    obtain ⟨k, hk⟩ := hconj (istr2IsoEquiv (F := F) A θ) c d H' (e u)
      ((istr2_mem_otimes_iff A u).mp hu)
    refine ⟨k, e.injective ?_⟩
    rw [map_pow, map_pow, istr2_endConj]
    exact hk

/-! ## ★5. `Remark 4.5.1` の standard 型の側 -/

/-- ★★★★**[FrdI] Remark 4.5.1(standard 型の側)** ——
`𝒞` が standard 型なら `𝒞^istr` も standard 型。

★`Definition 3.1, (i)` の 5 条がすべて移った。 -/
theorem istr_standardType (h : IsOfStandardType D C P F) :
    IsOfStandardType D (Istr P) (istrPre P F) (istr_frobenioidCore P F) where
  quasiIsotropic := istr_quasiIsotropicType
  frobIsotropic := istr_frobIsotropicType
  groupLikeCompact := fun hgl => by
    obtain ⟨A, hA⟩ := h.groupLikeCompact (istr_groupLikeType_down hgl)
    exact ⟨istr2 (F := F) A, istr2_frobeniusCompact A hA⟩
  frobNormalized := istr_frobNormalizedType h.frobNormalized
  baseFSMFF := h.baseFSMFF
  phiNonDilating := h.phiNonDilating

/-! ## ★6. 残り —— rationally standard 型の側

★★**未実装**。`Definition 4.5, (iii)` の 4 条:

| 条 | 内容 | 状態 |
|---|---|---|
| (a) | birationally Frobenius-normalized 型 | 未 |
| (a) | rational 型 | 未 |
| (a) | standard 型 | ★**済**(`istr_standardType`) |
| (b) | `(𝒞^un-tr)^birat` に Frobenius-compact 対象 | 未 |

★いずれも `𝒞^istr` の birationalization / unit-trivialization を経由するので、
`Prop33UnTr.lean` / `Prop44Pre.lean` の語彙を `𝒞^istr` へ持ち上げる作業になる。
-/

/-! ## ★出典の紐付け(条つき) -/

def istr_frobNormalizedType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — Definition 3.1, (i) の (a)(c) の移送",
    sectionId := "frdi-remark-4-5-1" }

end ABC3.Found.FrdI
