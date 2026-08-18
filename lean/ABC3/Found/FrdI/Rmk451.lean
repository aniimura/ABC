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

/-! ## ★3. 残り —— (b) の Frobenius-compact 対象

★★**未実装**。(b) は

  `𝒞^istr` が group-like 型 ⟹ `(𝒞^istr)^istr` に Frobenius-compact 対象がある

であり、2 段の移送が要る:

1. `IsOfGroupLikeType (istrPre P F) → IsOfGroupLikeType P`
   —— 各 `A : C` の isotropic hull `A ⟶ A^istr` は**底同型**なので
   `Φ(base A) ≅ Φ(base A^istr)`、group-like 性はその同型に沿って移る。
2. `IsFrobeniusCompact (istrPre P F) A → IsFrobeniusCompact (istrPre (istrPre P F) _) A'`
   —— `Istr (istrPre P F)` は `Istr P` の**全対象**からなる充満部分圏
   (`istr_isotropic`)なので、`End` も `OTimes` も一致する。

★この 2 段が済めば `Definition 3.1, (i)` の 5 条が揃い、
`Remark 4.5.1` の **standard 型の半分**が閉じる。
★rationally standard 型の側((a) の 3 条 ＋ (b))はさらに
`(𝒞^istr)^un-tr,birat` の Frobenius-compact 対象を要する。
-/

/-! ## ★出典の紐付け(条つき) -/

def istr_frobNormalizedType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — Definition 3.1, (i) の (a)(c) の移送",
    sectionId := "frdi-remark-4-5-1" }

end ABC3.Found.FrdI
