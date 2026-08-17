import ABC3.Found.Arakelov.PicSheafifyTrivial

/-!
# Arakelov (B1) 第 17 ブロック —— **層の側での並べ替え**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★これが `Pic X` の乗法の逆元を与える

    (A ⊗ B) ⊗ (A' ⊗ B') ≅ (A ⊗ A') ⊗ (B ⊗ B')

★★ここで `A ⊗ A' ≅ 𝟙` かつ `B ⊗ B' ≅ 𝟙` なら右辺は `𝟙` になる。
★したがって「可逆層のテンソル積の逆は、逆のテンソル積」である。

## ★機構

結合律(第 13 ブロック)と可換律(第 2 ブロック)を 5 回使う。
★★結合律には**外側 2 因子の局所階数 1** が要るので、
第 15・16 ブロック(局所自明性のテンソル閉性・層化での保存)がここで効く。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

variable {X : Scheme.{u}}

/-- ★局所自明なら局所階数 1 である(第 15 ブロックの言い換え、使いやすい形)。 -/
theorem rankOne_of_trivial {M : X.Modules} (h : IsLocallyTrivial X M.val) :
    IsLocallyRankOne X M.val := h.isLocallyRankOne

/-- ★★★★**層の側での並べ替え** `(A ⊗ B) ⊗ (A' ⊗ B') ≅ (A ⊗ A') ⊗ (B ⊗ B')`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `Pic X` の**逆元**を与える。 -/
noncomputable def tensorRearrange (A B A' B' : X.Modules)
    (hA : IsLocallyTrivial X A.val) (hB : IsLocallyTrivial X B.val)
    (hA' : IsLocallyTrivial X A'.val) (hB' : IsLocallyTrivial X B'.val) :
    tensorModules (tensorModules A B) (tensorModules A' B')
      ≅ tensorModules (tensorModules A A') (tensorModules B B') :=
  -- ★(A ⊗ B) ⊗ (A' ⊗ B') ≅ A ⊗ (B ⊗ (A' ⊗ B'))
  tensorModulesAssoc X A B (tensorModules A' B') (rankOne_of_trivial hA)
      (rankOne_of_trivial (isLocallyTrivial_tensorModules X A' B' hA' hB'))
  -- ★A ⊗ (B ⊗ (A' ⊗ B')) ≅ A ⊗ ((B ⊗ A') ⊗ B')
  ≪≫ tensorModulesIso (Iso.refl A)
      (tensorModulesAssoc X B A' B' (rankOne_of_trivial hB) (rankOne_of_trivial hB')).symm
  -- ★A ⊗ ((B ⊗ A') ⊗ B') ≅ A ⊗ ((A' ⊗ B) ⊗ B')
  ≪≫ tensorModulesIso (Iso.refl A)
      (tensorModulesIso (tensorModulesComm B A') (Iso.refl B'))
  -- ★A ⊗ ((A' ⊗ B) ⊗ B') ≅ A ⊗ (A' ⊗ (B ⊗ B'))
  ≪≫ tensorModulesIso (Iso.refl A)
      (tensorModulesAssoc X A' B B' (rankOne_of_trivial hA') (rankOne_of_trivial hB'))
  -- ★A ⊗ (A' ⊗ (B ⊗ B')) ≅ (A ⊗ A') ⊗ (B ⊗ B')
  ≪≫ (tensorModulesAssoc X A A' (tensorModules B B') (rankOne_of_trivial hA)
      (rankOne_of_trivial (isLocallyTrivial_tensorModules X B B' hB hB'))).symm

/-! ## ★★★★★可逆層(層の側) -/

/-- ★★★**可逆層(直線束)** —— 逆を持ち、局所自明な層加群。 -/
structure InvSheaf (X : Scheme.{u}) where
  /-- 台。 -/
  carrier : X.Modules
  /-- テンソル積についての逆。 -/
  inv : X.Modules
  /-- ★逆であること。 -/
  isInv : Nonempty (tensorModules carrier inv ≅ unitModules X)
  /-- ★局所自明。 -/
  trivial : IsLocallyTrivial X carrier.val
  /-- ★逆の側も局所自明。 -/
  invTrivial : IsLocallyTrivial X inv.val

namespace InvSheaf

/-- ★**構造層は可逆層である**(単位元)。 -/
noncomputable def one (X : Scheme.{u}) : InvSheaf X where
  carrier := unitModules X
  inv := unitModules X
  isInv := ⟨tensorUnitLeft (unitModules X)⟩
  trivial := isLocallyTrivial_unitModules X
  invTrivial := isLocallyTrivial_unitModules X

/-- ★★**逆を取る操作**(逆元)。 -/
noncomputable def symm (L : InvSheaf X) : InvSheaf X where
  carrier := L.inv
  inv := L.carrier
  isInv := L.isInv.map fun e => tensorModulesComm L.inv L.carrier ≪≫ e
  trivial := L.invTrivial
  invTrivial := L.trivial

/-- ★★★★**テンソル積**(乗法)。

★★★逆の存在は `tensorRearrange` で出る:

    (A ⊗ B) ⊗ (A' ⊗ B') ≅ (A ⊗ A') ⊗ (B ⊗ B') ≅ 𝟙 ⊗ 𝟙 ≅ 𝟙

★局所自明性は第 15・16 ブロック(前層テンソル・層化の両方が保つ)。 -/
noncomputable def mul (L M : InvSheaf X) : InvSheaf X where
  carrier := tensorModules L.carrier M.carrier
  inv := tensorModules L.inv M.inv
  isInv := by
    obtain ⟨eL⟩ := L.isInv
    obtain ⟨eM⟩ := M.isInv
    exact ⟨tensorRearrange L.carrier M.carrier L.inv M.inv
        L.trivial M.trivial L.invTrivial M.invTrivial
      ≪≫ tensorModulesIso eL eM ≪≫ tensorUnitLeft (unitModules X)⟩
  trivial := isLocallyTrivial_tensorModules X L.carrier M.carrier L.trivial M.trivial
  invTrivial := isLocallyTrivial_tensorModules X L.inv M.inv L.invTrivial M.invTrivial

@[simp] theorem mul_carrier (L M : InvSheaf X) :
    (L.mul M).carrier = tensorModules L.carrier M.carrier := rfl

@[simp] theorem one_carrier (X : Scheme.{u}) :
    (one X).carrier = unitModules X := rfl

/-- ★★**逆の一意性**——台が同型なら逆も同型である。

★モノイダル圏の標準の議論。★★結合律に局所階数 1 が要るので仮定を運ぶ。 -/
theorem invIso_of_carrierIso {L M : InvSheaf X} (e : L.carrier ≅ M.carrier) :
    Nonempty (L.inv ≅ M.inv) := by
  obtain ⟨eL⟩ := L.isInv
  obtain ⟨eM⟩ := M.isInv
  refine ⟨(tensorUnitLeft L.inv).symm
    ≪≫ tensorModulesIso eM.symm (Iso.refl L.inv)
    ≪≫ tensorModulesAssoc X M.carrier M.inv L.inv
        (rankOne_of_trivial M.trivial) (rankOne_of_trivial L.invTrivial)
    ≪≫ tensorModulesIso (Iso.refl M.carrier) (tensorModulesComm M.inv L.inv)
    ≪≫ (tensorModulesAssoc X M.carrier L.inv M.inv
        (rankOne_of_trivial M.trivial) (rankOne_of_trivial M.invTrivial)).symm
    ≪≫ tensorModulesIso (tensorModulesIso e.symm (Iso.refl L.inv)) (Iso.refl M.inv)
    ≪≫ tensorModulesIso eL (Iso.refl M.inv)
    ≪≫ tensorUnitLeft M.inv⟩

/-- ★可逆層の同型による同値関係(台が同型であること)。 -/
def setoid (X : Scheme.{u}) : Setoid (InvSheaf X) where
  r L M := Nonempty (L.carrier ≅ M.carrier)
  iseqv :=
    { refl := fun _ => ⟨Iso.refl _⟩
      symm := fun ⟨e⟩ => ⟨e.symm⟩
      trans := fun ⟨e⟩ ⟨f⟩ => ⟨e ≪≫ f⟩ }

end InvSheaf

/-- ★★★★★★**`Pic X`** —— 可逆層の同型類。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any -/
def PicSheaf (X : Scheme.{u}) : Type _ := Quotient (InvSheaf.setoid X)

namespace PicSheaf

/-- ★同型類を取る。 -/
def mk (L : InvSheaf X) : PicSheaf X := Quotient.mk _ L

theorem mk_eq_mk {L M : InvSheaf X} (e : L.carrier ≅ M.carrier) : mk L = mk M :=
  Quotient.sound ⟨e⟩

/-- ★★★★★★**`Pic X` は可換群である**。

★★★これが B1 の中核である。 -/
noncomputable instance : CommGroup (PicSheaf X) where
  mul := Quotient.map₂ InvSheaf.mul
    (by
      rintro L L' ⟨e⟩ M M' ⟨f⟩
      exact ⟨tensorModulesIso e f⟩)
  one := mk (InvSheaf.one X)
  inv := Quotient.map InvSheaf.symm
    (by
      rintro L L' ⟨e⟩
      exact InvSheaf.invIso_of_carrierIso e)
  mul_assoc := by
    rintro ⟨L⟩ ⟨M⟩ ⟨N⟩
    exact mk_eq_mk (tensorModulesAssoc X L.carrier M.carrier N.carrier
      (rankOne_of_trivial L.trivial) (rankOne_of_trivial N.trivial))
  one_mul := by
    rintro ⟨L⟩
    exact mk_eq_mk (tensorUnitLeft L.carrier)
  mul_one := by
    rintro ⟨L⟩
    exact mk_eq_mk (tensorUnitRight L.carrier)
  inv_mul_cancel := by
    rintro ⟨L⟩
    obtain ⟨e⟩ := L.symm.isInv
    exact mk_eq_mk e
  mul_comm := by
    rintro ⟨L⟩ ⟨M⟩
    exact mk_eq_mk (tensorModulesComm L.carrier M.carrier)

end PicSheaf

/-! ## ★出典の紐付け(`.src`) -/

def PicSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Pic X = 可逆層の同型類とその可換群構造)",
    sectionId := "genell-def-1-1-i" }

def InvSheaf.mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆層のテンソル積と逆元)",
    sectionId := "genell-def-1-1-i" }

def tensorRearrange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——テンソル積の並べ替え)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
