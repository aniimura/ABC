/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.AMetricMonoid
import ABC3.Meta.Claim

/-!
# 算術直線束の**等長同型**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★原文の `APic(X)` は**同型類**の群である

`AMetric X`（`AMetricMonoid.lean`）はテンソル積を持つが、群になるのは
**等長同型で割った後**である。★本ファイルはその同値関係を作る。

## ★★★★★★等長性は「基準ノルム `h` の等式」として書く

    `IsIsometry L M φ ≝ ∀ V e p, p ∈ V → L.h V (pullTriv φ V e) p = M.h V e p`

★★**切断のノルムで書くと穴が開く**——`Γ(L ⊗ M, ⊤) = Γ(L,⊤) ⊗ Γ(M,⊤)` の一般の元は
純テンソルの**和**であり、ノルムは加法的でないので純テンソルでの等式から一般へ行けない。
★★★`h` の側なら切断を経由しないので、テンソル積へ降ろす段で困らない。

★それでも「ノルムを保つ」は定理として出る（`normOf_of_isIsometry`）——
`trivValue` は線型なので `trivValue L V (pullTriv φ V e) s = trivValue M V e (φ s)` が
`φ.hom` の**自然性**から出るからである。

## ★残っている段（明示）

1. ★積が同値関係で降りること（`IsIsometry (L*M) (L'*M') (φ ⊗ᵢ ψ)`）
2. ★★結合律・単位律
3. ★★★逆元（双対計量）
4. ★★★★商そのものと群のインスタンス
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-- ★同型に沿って自明化を引き戻す。 -/
noncomputable def pullTriv {L M : X.PresheafOfModules} (φ : L ≅ M) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) :
    (restrictPresheafFunctor X V).obj L ≅ 𝟙_ (PresheafModulesOn X V) :=
  (restrictPresheafFunctor X V).mapIso φ ≪≫ e

set_option backward.isDefEq.respectTransparency false in
theorem trivValue_pullTriv {L M : X.PresheafOfModules} (φ : L ≅ M) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : L.obj (op ⊤)) :
    trivValue L V (pullTriv φ V e) s = trivValue M V e (φ.hom.app (op ⊤) s) := by
  have hnat := PresheafOfModules.naturality_apply φ.hom
    (homOfLE (le_top : V ≤ (⊤ : X.Opens))).op s
  have h1 : trivValue L V (pullTriv φ V e) s
      = e.hom.app (op (Over.mk (𝟙 V))) (φ.hom.app (op V) (secOn L V s)) := rfl
  have h2 : trivValue M V e (φ.hom.app (op ⊤) s)
      = e.hom.app (op (Over.mk (𝟙 V))) (secOn M V (φ.hom.app (op ⊤) s)) := rfl
  rw [h1, h2]
  exact congrArg _ hnat

theorem pullTriv_refl {L : X.PresheafOfModules} (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj L ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullTriv (Iso.refl L) V e = e := by
  ext1
  show ((restrictPresheafFunctor X V).mapIso (Iso.refl L)).hom ≫ e.hom = e.hom
  simp

theorem pullTriv_trans {L M N : X.PresheafOfModules} (φ : L ≅ M) (ψ : M ≅ N) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj N ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullTriv (φ ≪≫ ψ) V e = pullTriv φ V (pullTriv ψ V e) := by
  ext1
  show ((restrictPresheafFunctor X V).mapIso (φ ≪≫ ψ)).hom ≫ e.hom
    = ((restrictPresheafFunctor X V).mapIso φ).hom
      ≫ (((restrictPresheafFunctor X V).mapIso ψ).hom ≫ e.hom)
  simp

theorem pullTriv_symm_pullTriv {L M : X.PresheafOfModules} (φ : L ≅ M) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj L ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullTriv φ V (pullTriv φ.symm V e) = e := by
  rw [← pullTriv_trans, φ.self_symm_id, pullTriv_refl]

/-! ## ★★★★★★等長同型 -/

/-- ★★★★★★**等長同型**——基準ノルムが引き戻しで一致すること。 -/
def IsIsometry (L M : AMetric X) (φ : L.sheaf ≅ M.sheaf) : Prop :=
  ∀ (V : X.Opens) (e : (restrictPresheafFunctor X V).obj M.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X), p ⁻¹ᵁ V = ⊤ →
    L.metric.h V (pullTriv φ V e) p = M.metric.h V e p

/-- ★★★★★★★**等長同型は切断のノルムを保つ**。 -/
theorem normOf_of_isIsometry {L M : AMetric X} {φ : L.sheaf ≅ M.sheaf}
    (hφ : IsIsometry L M φ) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : L.sheaf.obj (op ⊤)) :
    M.normOf V e p hp (φ.hom.app (op ⊤) s) = L.normOf V (pullTriv φ V e) p hp s := by
  show trivSecNorm M.sheaf V e p hp _ * M.metric.h V e p
    = trivSecNorm L.sheaf V (pullTriv φ V e) p hp s * L.metric.h V (pullTriv φ V e) p
  rw [hφ V e p hp]
  congr 1
  show ‖evalOn p V hp (trivValue M.sheaf V e (φ.hom.app (op ⊤) s))‖ = _
  rw [← trivValue_pullTriv φ V e s]
  rfl

/-- ★★★★★★★★**チャートを選ばない形**——ノルムは同型で保たれる。 -/
theorem normOf_of_isIsometry' {L M : AMetric X} {φ : L.sheaf ≅ M.sheaf}
    (hφ : IsIsometry L M φ) (V V' : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (e' : (restrictPresheafFunctor X V').obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V'))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (hp' : p ⁻¹ᵁ V' = ⊤)
    (s : L.sheaf.obj (op ⊤)) :
    M.normOf V e p hp (φ.hom.app (op ⊤) s) = L.normOf V' e' p hp' s :=
  (normOf_of_isIsometry hφ V e p hp s).trans
    (L.normOf_chart_indep (pullTriv φ V e) e' p hp hp' s)

/-! ## ★★★★★同値関係 -/

theorem isIsometry_refl (L : AMetric X) : IsIsometry L L (Iso.refl L.sheaf) := by
  intro V e p _
  rw [pullTriv_refl]

theorem isIsometry_symm {L M : AMetric X} {φ : L.sheaf ≅ M.sheaf} (hφ : IsIsometry L M φ) :
    IsIsometry M L φ.symm := by
  intro V e p hp
  have h := hφ V (pullTriv φ.symm V e) p hp
  rw [pullTriv_symm_pullTriv] at h
  exact h.symm

theorem isIsometry_trans {L M N : AMetric X} {φ : L.sheaf ≅ M.sheaf} {ψ : M.sheaf ≅ N.sheaf}
    (hφ : IsIsometry L M φ) (hψ : IsIsometry M N ψ) : IsIsometry L N (φ ≪≫ ψ) := by
  intro V e p hp
  rw [pullTriv_trans, hφ V (pullTriv ψ V e) p hp, hψ V e p hp]

/-- ★★**等長同型でつながること**。 -/
def Isometric (L M : AMetric X) : Prop := ∃ φ : L.sheaf ≅ M.sheaf, IsIsometry L M φ

theorem isometric_refl (L : AMetric X) : Isometric L L := ⟨Iso.refl _, isIsometry_refl L⟩

theorem isometric_symm {L M : AMetric X} (h : Isometric L M) : Isometric M L :=
  h.elim fun φ hφ => ⟨φ.symm, isIsometry_symm hφ⟩

theorem isometric_trans {L M N : AMetric X} (h : Isometric L M) (h' : Isometric M N) :
    Isometric L N :=
  h.elim fun φ hφ => h'.elim fun ψ hψ => ⟨φ ≪≫ ψ, isIsometry_trans hφ hψ⟩

/-- ★★★**等長同型は同値関係である**。 -/
def isometricSetoid (X : Scheme.{0}) : Setoid (AMetric X) where
  r := Isometric
  iseqv := ⟨isometric_refl, isometric_symm, isometric_trans⟩

/-! ### ★出典の紐付け(`.src`) -/

def IsIsometry.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(算術直線束の等長同型)",
    sectionId := "genell-def-1-1-i" }

def normOf_of_isIsometry.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(等長同型は切断のノルムを保つこと)",
    sectionId := "genell-def-1-1-i" }

def isometricSetoid.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(等長同型は同値関係であること)",
    sectionId := "genell-def-1-1-i" }

def isometricSetoid.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "trivValue_pullTriv(引き戻した自明化の値)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivValue_pullTriv") 3,
    .implicitStep
      ("★等長性は**基準ノルム `h` の等式**として書いた——切断のノルムで書くと" ++
       "「純テンソルでしか確かめられないがノルムは加法的でない」という穴が開く。" ++
       "★★`h` の側なら切断を経由しないので、テンソル積へ降ろす段で困らない") 3,
    .implicitStep
      ("★★★残っているのは (1) 積が同値関係で降りること、" ++
       "(2) 結合律・単位律、(3) 逆元(双対計量)、そして商そのものである") 3 ]

end ABC3.Found.Arakelov
