import ABC3.Found.Arakelov.PicRestrictLocal
import ABC3.Found.Arakelov.PicAssoc

/-!
# Arakelov (B1) 第 16 ブロック —— **層化は局所自明性を保つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★これが B1 の最後の壁である

`Pic X` の乗法は `tensorModules L M = 層化 (L.val ⊗ M.val)` である。
★★これが再び可逆層になるには「**層化が局所自明性を保つ**」が要る。

## ★★4 部品(すべて mathlib に在った、2026-08-17 実測)

| 部品 | 在庫 |
|---|---|
| `Over.forget V` が cocontinuous | `Sites/Over.lean` |
| 制限が局所全射/単射を保つ | ★`Sites/PreservesLocallyBijective.lean` |
| 制限が層を保つ | ★`Functor.op_comp_isSheaf`(`Over.forget` は continuous) |
| 層どうしの局所全単射は同型 | ★`Sheaf.isLocallyBijective_iff_isIso` |

## ★機構

    η_P|_V : P|_V ⟶ (層化 P).val|_V
      ・局所全単射(制限が保つ)
      ・P|_V は層(≅ 𝟙_|_V、そして 𝟙_ の制限は層)
      ・(層化 P).val|_V も層(制限は層を保つ)
    ⟹ 同型
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable (X : Scheme.{u}) (V : X.Opens)

/-- ★**層の制限は層である**。 -/
theorem isSheaf_restrict (F : (X.Opens)ᵒᵖ ⥤ AddCommGrpCat.{u})
    (h : Presheaf.IsSheaf (Opens.grothendieckTopology X) F) :
    Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((Over.forget V).op ⋙ F) :=
  (Over.forget V).op_comp_isSheaf _ _ ⟨F, h⟩

/-- ★★★★★**制限した層化の単位は、局所自明な開集合の上で同型である**。

★これが壁の核心。★★機構は「層どうしの局所全単射は同型」。 -/
theorem isIso_restrict_sheafifyUnit (P : X.PresheafOfModules)
    (hP : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((Over.forget V).op ⋙ P.presheaf)) :
    IsIso (Functor.whiskerLeft (Over.forget V).op
      ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (sheafifyUnit X P))) := by
  haveI := isLocallySurjective_restrict X V
    ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (sheafifyUnit X P))
    (isLocallySurjective_sheafifyUnit X P)
  haveI := isLocallyInjective_restrict X V
    ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (sheafifyUnit X P))
    (isLocallyInjective_sheafifyUnit X P)
  have htgt : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((Over.forget V).op ⋙
        ((PresheafOfModules.sheafification (R := X.ringCatSheaf)
          (𝟙 X.ringCatSheaf.obj)).obj P).val.presheaf) :=
    isSheaf_restrict X V _
      (((PresheafOfModules.sheafification (R := X.ringCatSheaf)
        (𝟙 X.ringCatSheaf.obj)).obj P).isSheaf)
  let φ : (⟨_, hP⟩ : Sheaf ((Opens.grothendieckTopology X).over V) AddCommGrpCat.{u})
      ⟶ ⟨_, htgt⟩ :=
    ⟨Functor.whiskerLeft (Over.forget V).op
      ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (sheafifyUnit X P))⟩
  haveI hinj : Sheaf.IsLocallyInjective φ :=
    isLocallyInjective_restrict X V _ (isLocallyInjective_sheafifyUnit X P)
  haveI hsur : Sheaf.IsLocallySurjective φ :=
    isLocallySurjective_restrict X V _ (isLocallySurjective_sheafifyUnit X P)
  haveI : IsIso φ := (Sheaf.isLocallyBijective_iff_isIso φ).1 ⟨hinj, hsur⟩
  exact inferInstanceAs (IsIso ((sheafToPresheaf
    ((Opens.grothendieckTopology X).over V) AddCommGrpCat.{u}).map φ))

/-- ★★**制限した層化の単位は前層加群の射としても同型である**。

★機構は `toPresheaf` が同型を反映すること + 制限関手の台が `whiskerLeft` そのものであること
(どちらも 2026-08-17 に探りで確認)。 -/
theorem isIso_restrictMap_sheafifyUnit (P : X.PresheafOfModules)
    (hP : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((Over.forget V).op ⋙ P.presheaf)) :
    IsIso ((restrictPresheafFunctor X V).map (sheafifyUnit X P)) := by
  haveI := isIso_restrict_sheafifyUnit X V P hP
  refine (isIso_iff_of_reflects_iso ((restrictPresheafFunctor X V).map (sheafifyUnit X P))
    (PresheafOfModules.toPresheaf
      (((Over.forget V).op ⋙ X.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u}))).1 ?_
  exact ‹IsIso (Functor.whiskerLeft (Over.forget V).op
    ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (sheafifyUnit X P)))›

/-- ★★★★★★**層化は局所自明性を保つ**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★**これが B1 の最後の壁だった。**
これで `tensorModules L M = 層化 (L.val ⊗ M.val)` が再び可逆層になる。 -/
theorem isLocallyTrivial_sheafify (P : X.PresheafOfModules)
    (h : IsLocallyTrivial X P) :
    IsLocallyTrivial X ((PresheafOfModules.sheafification (R := X.ringCatSheaf)
      (𝟙 X.ringCatSheaf.obj)).obj P).val := by
  intro U
  obtain ⟨S, hS, hV⟩ := h U
  refine ⟨S, hS, ?_⟩
  intro V i hi
  obtain ⟨e⟩ := hV i hi
  -- ★`P|_V` は層である(`𝟙_|_V` と同型で、`𝟙_|_V` は層)
  have hunit : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((Over.forget V).op ⋙
        (PresheafOfModules.unit X.ringCatSheaf.obj).presheaf) :=
    isSheaf_restrict X V _
      ((sheafCompose (Opens.grothendieckTopology X)
        (forget₂ RingCat.{u} AddCommGrpCat.{u})).obj X.ringCatSheaf).property
  have hP : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((Over.forget V).op ⋙ P.presheaf) :=
    (Presheaf.isSheaf_of_iso_iff
      ((PresheafOfModules.toPresheaf _).mapIso e)).2 hunit
  haveI := isIso_restrictMap_sheafifyUnit X V P hP
  exact ⟨(asIso ((restrictPresheafFunctor X V).map (sheafifyUnit X P))).symm ≪≫ e⟩

/-! ## ★出典の紐付け(`.src`) -/

def isSheaf_restrict.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層の制限は層であること)",
    sectionId := "genell-def-1-1-i" }

def isLocallyTrivial_sheafify.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層化が局所自明性を保つこと)",
    sectionId := "genell-def-1-1-i" }

def isIso_restrict_sheafifyUnit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限した層化の単位が同型であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
