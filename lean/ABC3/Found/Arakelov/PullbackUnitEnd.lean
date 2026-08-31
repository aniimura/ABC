/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PullbackTriv
import ABC3.Meta.Claim

/-!
# **遷移単元は引き戻しと可換である**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★★★★到達点

    `transUnit (f^*L) V' (pullTrivOfBase eW) (pullTrivOfBase eW')`
      `= ρ (transUnit L W eW eW')`      （`transUnit_pullTrivOfBase`）

★これが `Definition 1.1, (ii)` の**計量を運ぶ段の核**である
——`LocalMetric.tensor` における `transUnit_tensorTriv` に当たる。

## ★★★★★★★★なぜ `ρ` が現れるか —— **随伴の転置**

`PullbackTriv.lean` で 3 段まで削った:

| 削った段 | 内容 |
|---|---|
| `pullTrivOfBase_comp` | `bcIso` が両側で打ち消し、比は輸送の共役になる |
| `transUnit_eq_unitEnd` | 遷移単元は `unitEnd`（環準同型 `End(𝟙_) →+* Γ`） |
| `unitEnd_pullTransport_reduce` | 制限の層が剥がれる |

★残った 1 本 `pullUnitEnd f W Ψ = ρ (unitEnd W Ψ)` を、本ファイルが
**随伴の転置**で証明する:

1. `pullUnitEnd f W Ψ = freeYonedaEquiv (ftI'.hom ≫ 共役)`（`freeYonedaEquiv_termIso_comp`）
2. `ftI'.hom ≫ pullbackOnUnitIso.inv = pOF.inv ≫ (f|)^*(ftI.hom)`（`termIso_comp_pullUnitInv`）
3. `freeYonedaEquiv (pOF.inv ≫ κ) = κ.app (freeYonedaEquiv pOF.inv)` かつ
   `freeYonedaEquiv pOF.inv = 随伴の unit`（在庫の `freeYonedaEquivUnitGen`）
   ——これで `κ.app (unit) = (adj.homEquiv κ).app (生成元)` になる（`Adjunction.homEquiv_unit`）
4. `adj.homEquiv` の**左自然性**で `ftI.hom` と `Ψ` を外へ出す
   ——`pullUnitEnd f W Ψ = η.app (Ψ.app (unitOne W))`、
   ここで `η ≝ adj.homEquiv (pullbackOnUnitIso.hom)`
5. ★`Ψ = 𝟙` に取ると **`η.app (unitOne W) = 1`**（`etaOn_unitOne`）
6. ★★`η.app` は `Γ(Y,W)`-線型で、行き先は**係数制限**なので

       `η.app (u • 1) = u • η.app 1 = ρ(u) • 1 = ρ u`

★★★**`ρ` は `η` の線型性から出る**——これが答えである。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory TopologicalSpace MonoidalCategory Opposite

/-! ## ★★`freeYonedaEquiv` の道具 -/

/-- ★一般形の「生成元での値」。 -/
theorem freeYonedaEquiv_apply_gen {C : Type} [SmallCategory C] {R : Cᵒᵖ ⥤ RingCat}
    {M : PresheafOfModules R} {Z : C}
    (m : (PresheafOfModules.free R).obj (yoneda.obj Z) ⟶ M) :
    PresheafOfModules.freeYonedaEquiv m = (m.app (op Z)) (ModuleCat.freeMk (𝟙 Z)) := by
  have h := PresheafOfModules.freeYonedaEquiv_symm_app M Z (PresheafOfModules.freeYonedaEquiv m)
  rw [Equiv.symm_apply_apply] at h
  exact h.symm

/-- ★★**`unitEnd` は `freeYonedaEquiv` で書ける**（`unitMul_unitEnd` の証明の中身を外に出した）。 -/
theorem freeYonedaEquiv_termIso_comp {X : Scheme.{0}} (U : X.Opens)
    (φ : 𝟙_ (PresheafModulesOn X U) ⟶ 𝟙_ (PresheafModulesOn X U)) :
    PresheafOfModules.freeYonedaEquiv
        ((freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom ≫ φ)
      = unitEnd U φ := by
  rw [PresheafOfModules.freeYonedaEquiv_comp, freeYonedaEquiv_apply' U, termIso_hom_gen]
  rfl

/-- ★同上（`unitMul` の側）。 -/
theorem termIso_comp_unitMul {X : Scheme.{0}} (U : X.Opens) (c : (Γ(X, U) : Type)) :
    (freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom ≫ unitMul U c
      = PresheafOfModules.freeYonedaEquiv.symm c := by
  show _ ≫ ((freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).inv ≫ _) = _
  rw [← Category.assoc,
    (freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom_inv_id, Category.id_comp]

/-! ## ★★★★自由対象へ降ろす -/

/-- ★★`ftI'.hom ≫ pullbackOnUnitIso.inv` は自由対象の側で書ける。 -/
theorem termIso_comp_pullUnitInv {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens) :
    (freeYonedaTermIso (Over.mk (𝟙 ((Opens.map f.base).obj W)))
        (overTerminalUnique ((Opens.map f.base).obj W))).hom ≫ (pullbackOnUnitIso f W).inv
      = (pullbackOnFreeIsoGen f W (Over.mk (𝟙 W))).inv
        ≫ (pullbackPreOn f W).map
            (freeYonedaTermIso (Over.mk (𝟙 W)) (overTerminalUnique W)).hom := by
  simp only [pullbackOnUnitIso, Iso.trans_inv, Category.assoc,
    Functor.mapIso_inv, Iso.symm_inv]
  exact Iso.hom_inv_id_assoc _ _

/-- ★★★共役を自由対象の側へ移した形。 -/
theorem pullUnitEnd_keyComp {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    (Ψ : End (𝟙_ (PresheafModulesOn Y W))) :
    (freeYonedaTermIso (Over.mk (𝟙 ((Opens.map f.base).obj W)))
        (overTerminalUnique ((Opens.map f.base).obj W))).hom
      ≫ ((pullbackOnUnitIso f W).inv
        ≫ ((pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom))
    = (pullbackOnFreeIsoGen f W (Over.mk (𝟙 W))).inv
      ≫ ((pullbackPreOn f W).map (freeYonedaTermIso (Over.mk (𝟙 W)) (overTerminalUnique W)).hom
        ≫ ((pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom)) := by
  rw [← Category.assoc (freeYonedaTermIso (Over.mk (𝟙 ((Opens.map f.base).obj W)))
      (overTerminalUnique ((Opens.map f.base).obj W))).hom,
    ← Category.assoc (pullbackOnFreeIsoGen f W (Over.mk (𝟙 W))).inv]
  exact congrArg (fun m => m ≫ ((pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom))
    (termIso_comp_pullUnitInv f W)

/-- ★★★★`pullUnitEnd` を生成元での値として書く。 -/
theorem pullUnitEnd_eq_appGen {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    (Ψ : End (𝟙_ (PresheafModulesOn Y W))) :
    pullUnitEnd f W Ψ
      = (((pullbackPreOn f W).map (freeYonedaTermIso (Over.mk (𝟙 W)) (overTerminalUnique W)).hom
          ≫ ((pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom)).app
          (op (Over.mk (𝟙 ((Opens.map f.base).obj W)))))
        (PresheafOfModules.freeYonedaEquiv
          (pullbackOnFreeIsoGen f W (Over.mk (𝟙 W))).inv) := by
  refine (freeYonedaEquiv_termIso_comp _ _).symm.trans ?_
  refine (congrArg PresheafOfModules.freeYonedaEquiv (pullUnitEnd_keyComp f W Ψ)).trans ?_
  exact PresheafOfModules.freeYonedaEquiv_comp _ _

/-- ★★★★★★**生成元での値は随伴の転置で書ける**。

★機構は在庫の `freeYonedaEquivUnitGen`（`freeYonedaEquiv` の `pOF.inv` は随伴の `unit`）と
mathlib の `Adjunction.homEquiv_unit`。 -/
theorem appGen_eq_homEquiv {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    (κ : (pullbackPreOn f W).obj
        ((PresheafOfModules.free (((Over.forget W).op ⋙ Y.presheaf)
          ⋙ forget₂ CommRingCat RingCat)).obj (yoneda.obj (Over.mk (𝟙 W))))
      ⟶ 𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj W))) :
    (κ.app (op (Over.mk (𝟙 ((Opens.map f.base).obj W)))))
        (PresheafOfModules.freeYonedaEquiv
          (pullbackOnFreeIsoGen f W (Over.mk (𝟙 W))).inv)
      = (((PresheafOfModules.pullbackPushforwardAdjunction (phiOn f W)).homEquiv _ _ κ).app
          (op (Over.mk (𝟙 W)))) (ModuleCat.freeMk (𝟙 (Over.mk (𝟙 W)))) := by
  have hu : PresheafOfModules.freeYonedaEquiv
      (pullbackOnFreeIsoGen f W (Over.mk (𝟙 W))).inv
      = PresheafOfModules.freeYonedaEquiv
        ((PresheafOfModules.pullbackPushforwardAdjunction (phiOn f W)).unit.app
          ((PresheafOfModules.free (((Over.forget W).op ⋙ Y.presheaf)
            ⋙ forget₂ CommRingCat RingCat)).obj (yoneda.obj (Over.mk (𝟙 W))))) :=
    (freeYonedaEquivUnitGen (phiOn f W) (Over.mk (𝟙 W))).symm
  rw [Adjunction.homEquiv_unit]
  refine (congrArg (fun y => (κ.app (op (Over.mk (𝟙 ((Opens.map f.base).obj W))))) y) hu).trans ?_
  rw [freeYonedaEquiv_apply_gen]
  rfl

/-! ## ★★★★★★★★随伴の転置 `η` -/

/-- ★★★★**`pullbackOnUnitIso.hom` の随伴転置** `η : 𝟙_|_W ⟶ f_*(𝟙_|_{f⁻¹W})`。 -/
noncomputable def etaOn {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens) :=
  ((PresheafOfModules.pullbackPushforwardAdjunction (phiOn f W)).homEquiv
    (𝟙_ (PresheafModulesOn Y W))
    (𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj W)))) (pullbackOnUnitIso f W).hom

theorem pullUnitEnd_eq_homEquiv {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    (Ψ : End (𝟙_ (PresheafModulesOn Y W))) :
    pullUnitEnd f W Ψ
      = (((PresheafOfModules.pullbackPushforwardAdjunction (phiOn f W)).homEquiv _ _
          ((pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom)).app
          (op (Over.mk (𝟙 W)))) (unitOne W) := by
  refine (pullUnitEnd_eq_appGen f W Ψ).trans ?_
  refine (appGen_eq_homEquiv f W _).trans ?_
  have hnat : ((PresheafOfModules.pullbackPushforwardAdjunction (phiOn f W)).homEquiv
        ((PresheafOfModules.free (((Over.forget W).op ⋙ Y.presheaf)
          ⋙ forget₂ CommRingCat RingCat)).obj (yoneda.obj (Over.mk (𝟙 W))))
        (𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj W))))
      ((pullbackPreOn f W).map (freeYonedaTermIso (Over.mk (𝟙 W)) (overTerminalUnique W)).hom
        ≫ ((pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom))
      = (freeYonedaTermIso (Over.mk (𝟙 W)) (overTerminalUnique W)).hom
        ≫ ((PresheafOfModules.pullbackPushforwardAdjunction (phiOn f W)).homEquiv
            (𝟙_ (PresheafModulesOn Y W))
            (𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj W))))
          ((pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom) :=
    Adjunction.homEquiv_naturality_left _ _ _
  refine (congrArg (fun m => (m.app (op (Over.mk (𝟙 W))))
    (ModuleCat.freeMk (𝟙 (Over.mk (𝟙 W))))) hnat).trans ?_
  show ((((PresheafOfModules.pullbackPushforwardAdjunction (phiOn f W)).homEquiv _ _)
        ((pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom)).app (op (Over.mk (𝟙 W))))
      (((freeYonedaTermIso (Over.mk (𝟙 W)) (overTerminalUnique W)).hom.app (op (Over.mk (𝟙 W))))
        (ModuleCat.freeMk (𝟙 (Over.mk (𝟙 W))))) = _
  rw [termIso_hom_gen]

/-- ★★★★★**`pullUnitEnd` は `η` を `Ψ` の生成元での値に当てたものである**。 -/
theorem pullUnitEnd_eq_eta {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    (Ψ : End (𝟙_ (PresheafModulesOn Y W))) :
    pullUnitEnd f W Ψ
      = ((etaOn f W).app (op (Over.mk (𝟙 W))))
        ((Ψ.app (op (Over.mk (𝟙 W)))) (unitOne W)) := by
  refine (pullUnitEnd_eq_homEquiv f W Ψ).trans ?_
  have hnat : ((PresheafOfModules.pullbackPushforwardAdjunction (phiOn f W)).homEquiv
        (𝟙_ (PresheafModulesOn Y W))
        (𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj W))))
      ((pullbackPreOn f W).map Ψ ≫ (pullbackOnUnitIso f W).hom)
      = Ψ ≫ etaOn f W :=
    Adjunction.homEquiv_naturality_left _ _ _
  refine (congrArg (fun m => (m.app (op (Over.mk (𝟙 W)))) (unitOne W)) hnat).trans ?_
  rfl

/-- ★★★★★★**`η` は生成元を生成元へ送る**。

★機構は「`Ψ = 𝟙` に取る」だけである——`pullUnitEnd f W 𝟙 = unitEnd 𝟙 = 1`。 -/
theorem etaOn_unitOne {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens) :
    ((etaOn f W).app (op (Over.mk (𝟙 W)))) (unitOne W)
      = unitOne ((Opens.map f.base).obj W) := by
  have h := pullUnitEnd_eq_eta f W (𝟙 (𝟙_ (PresheafModulesOn Y W)))
  have h0 : pullUnitEnd f W (𝟙 (𝟙_ (PresheafModulesOn Y W))) = 1 := by
    show unitEnd _ ((pullbackOnUnitIso f W).inv
      ≫ (pullbackPreOn f W).map (𝟙 (𝟙_ (PresheafModulesOn Y W)))
      ≫ (pullbackOnUnitIso f W).hom) = 1
    rw [CategoryTheory.Functor.map_id, Category.id_comp, Iso.inv_hom_id]
    exact map_one (unitEnd _)
  exact h.symm.trans h0

/-! ## ★★★★★★★★★★到達点 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★**引き戻しが `End(𝟙_)` に誘導する写像は構造射である**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

    `pullUnitEnd f W Ψ = f.c.app (op W) (unitEnd W Ψ)`

★機構は 2 つだけ:
`η.app (unitOne W) = 1`（`etaOn_unitOne`）と、
`η.app` が `Γ(Y,W)`-線型でその行き先が**係数制限**であること。
★★後者から `η.app (u • 1) = u • 1 = ρ u` が出る——**これが `ρ` の出どころ**である。 -/
theorem pullUnitEnd_eq {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    (Ψ : End (𝟙_ (PresheafModulesOn Y W))) :
    pullUnitEnd f W Ψ = f.c.app (op W) (unitEnd W Ψ) := by
  refine (pullUnitEnd_eq_eta f W Ψ).trans ?_
  refine (congrArg (fun y => ((etaOn f W).app (op (Over.mk (𝟙 W)))) y)
    (smul_unitOne W ((Ψ.app (op (Over.mk (𝟙 W)))) (unitOne W))).symm).trans ?_
  rw [map_smul, etaOn_unitOne]
  have h2 := smul_unitOne ((Opens.map f.base).obj W)
    ((phiOn f W).app (op (Over.mk (𝟙 W))) ((Ψ.app (op (Over.mk (𝟙 W)))) (unitOne W)))
  exact (ModuleCat.restrictScalars.smul_def
    (M := (𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj W))).obj
      (op (Over.mk (𝟙 ((Opens.map f.base).obj W)))))
    ((phiOn f W).app (op (Over.mk (𝟙 W)))).hom
    ((Ψ.app (op (Over.mk (𝟙 W)))) (unitOne W))
    (unitOne ((Opens.map f.base).obj W))).trans h2

/-- ★★★★★★★★★★**遷移単元は引き戻しと可換である**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

    `transUnit (f^*L) V' (pullTrivOfBase eW) (pullTrivOfBase eW')`
      `= ρ (transUnit L W eW eW')`

★★これが `Definition 1.1, (ii)` の**計量を運ぶ段の核**である
——`LocalMetric.tensor` における `transUnit_tensorTriv` に当たる。 -/
theorem transUnit_pullTrivOfBase {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (eW eW' : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W) :
    transUnit ((pullbackPre f).obj L) V'
        (pullTrivOfBase f L W eW hV'W) (pullTrivOfBase f L W eW' hV'W)
      = X.presheaf.map (homOfLE hV'W).op (f.c.app (op W) (transUnit L W eW eW')) := by
  have e1 : transUnit ((pullbackPre f).obj L) V'
        (pullTrivOfBase f L W eW hV'W) (pullTrivOfBase f L W eW' hV'W)
      = unitEnd V' ((pullTrivOfBase f L W eW hV'W).symm
          ≪≫ pullTrivOfBase f L W eW' hV'W).hom := rfl
  have e2 : ((pullTrivOfBase f L W eW hV'W).symm
        ≪≫ pullTrivOfBase f L W eW' hV'W).hom
      = (pullUnitIsoOn f W hV'W).inv
        ≫ (pullTransportFun f W hV'W).map (eW.symm ≪≫ eW').hom
        ≫ (pullUnitIsoOn f W hV'W).hom :=
    congrArg Iso.hom (pullTrivOfBase_comp f L W eW eW' hV'W)
  rw [e1, e2, unitEnd_pullTransport_reduce, pullUnitEnd_eq]
  rfl

/-! ### ★出典の紐付け(`.src`) -/

def freeYonedaEquiv_termIso_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(unitEnd は freeYonedaEquiv で書けること)",
    sectionId := "genell-def-1-1-i" }

def etaOn.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(pullbackOnUnitIso の随伴転置 η)",
    sectionId := "genell-def-1-1-ii" }

def etaOn_unitOne.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(η は生成元を生成元へ送ること)",
    sectionId := "genell-def-1-1-ii" }

def pullUnitEnd_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻しが End(𝟙_) に誘導する写像は構造射であること)",
    sectionId := "genell-def-1-1-ii" }

def transUnit_pullTrivOfBase.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(遷移単元は引き戻しと可換であること)",
    sectionId := "genell-def-1-1-ii" }

def transUnit_pullTrivOfBase.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "pullTrivOfBase_comp(輸送した自明化の比は元の比の輸送の共役)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullTrivOfBase_comp") 3,
    .citation "[ABC3]" "unitEnd_pullTransport_reduce(制限の層が剥がれること)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.unitEnd_pullTransport_reduce") 3,
    .citation "[ABC3]" "freeYonedaEquivUnitGen(freeYonedaEquiv の pOF.inv は随伴の unit)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.freeYonedaEquivUnitGen") 3,
    .citation "[mathlib]" "Adjunction.homEquiv_unit / homEquiv_naturality_left"
      (.inMathlib "CategoryTheory.Adjunction.homEquiv_unit") 3,
    .implicitStep
      ("★★★★★**なぜ ρ が現れるか**: 随伴の転置 η ≝ adj.homEquiv (pullbackOnUnitIso.hom) は " ++
       "Γ(Y,W)-線型であり、その行き先 (f_* N).obj ⟨W,𝟙⟩ は**係数制限**なので、" ++
       "η.app (u • 1) = u • η.app 1 = ρ(u) • 1 = ρ u となる。" ++
       "★η.app 1 = 1 は「Ψ = 𝟙 に取る」だけで出る(etaOn_unitOne)") 3,
    .implicitStep
      ("★★残っている段: これで Definition 1.1 (ii) の計量を運ぶ段の**核**が揃った。" ++
       "あとは LocalMetric.tensor と同じ chartH_triv_indep / chartH_shrink / chartH_indep の " ++
       "3 段を書けば AMetric の引き戻しが出る") 3 ]

end ABC3.Found.Arakelov
