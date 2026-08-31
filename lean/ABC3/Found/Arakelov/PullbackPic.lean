/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PullbackTensor
import ABC3.Found.Arakelov.PicFreeTop
import ABC3.Found.Arakelov.AMetricPic
import ABC3.Meta.Claim

/-!
# **引き戻しは `APic` の群準同型である**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★★★★到達点

    `φ^* : APicM Y →* APicM X`

★原文が `evident` で畳んだ「引き戻し」は、`APic` が**群**である以上、
**群準同型**でなければ意味をなさない。★★本ファイルがそれを与える。

## ★★★★★★★★要る 3 つ

| 段 | 内容 | 場所 |
|---|---|---|
| 積 | `φ^*(L̄ ⊗ M̄) ≅ φ^*L̄ ⊗ φ^*M̄`（等長） | `PullbackTensor.lean`（§9-743） |
| 単位 | `φ^*Ō_Y ≅ Ō_X`（等長） | ★★**本ファイル** |
| 同値 | 等長同型は引き戻しで保たれる | ★★**本ファイル** |

★この 3 つが揃うと `AInv.pullback`（逆元つきの対象の引き戻し）が定まり、
`Quotient.map` で `APicM` の上に降りる。

## ★★★★★★★★単位の側 —— `f^*𝒪_Y ≅ 𝒪_X` が `1` を `1` に送ること

層の同型 `pullbackUnitPreIso` は在庫（`PicFreeTop.lean`）。★足すのは

    `(pullbackUnitPreIso f).hom (pullSec f 𝒪_Y ⊤ 1) = 1`

である。★★機構は `PullbackUnitEnd.lean` と同じ **`freeYonedaEquiv` の道**:

| 段 | 内容 | 道具 |
|---|---|---|
| 1 | `𝟙_ ≅ free (yoneda ⊤)` の下で `1` は生成元 | `freeYonedaTopIso_inv_one` |
| 2 | unit の自然性で自由対象の側へ移す | `unit.naturality` |
| 3 | **unit は生成元を生成元へ送る** | 在庫の `isoHomUnitGenGen` |
| 4 | `f⁻¹⊤ = ⊤` を通して戻す | `opensMap_top` ＋ `subst` |

★★★段 3 が在庫にあったのが効いた——`δ` と Beck–Chevalley の mate の**両方**に
使われている器具が、3 度目の出番を得た。

## ★★★★★★同値の側 —— またしても生成切断

`IsIsometry A B φ` から `IsIsometry (f^*A) (f^*B) (f^*φ)` を出すのに要るのは

    `trivGen (f^*A) V' (pullTriv (f^*φ) (pullTrivOfBase f B W e_W))`
      `= trivGen (f^*A) V' (pullTrivOfBase f A W (pullTriv φ W e_W))`

であり、★これは `pullSec` の**加群についての**自然性（`pullSec_naturality`）だけで出る。
★★`§9-741`・`§9-743` と同じ手口である——同型の等式ではなく生成切断の一致。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll

/-! ## ★★★★`free (yoneda ⊤)` の生成元と `1` -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★`free (yoneda ⊤) ≅ 𝟙_` は生成元を `1` に送る。 -/
theorem freeYonedaTopIso_hom_freeMk {Y : Scheme.{0}} (U : (Y.Opens)ᵒᵖ)
    (x : (yoneda.obj (⊤ : Y.Opens)).obj U) :
    (freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).hom.app U (ModuleCat.freeMk x)
      = (1 : ((Y.presheaf.obj U) : Type)) := by
  show (freeYonedaTopHom (R := Y.presheaf) (α := Y.Opens)).app U (ModuleCat.freeMk x) = _
  erw [ModuleCat.freeDesc_apply]
  rfl

set_option backward.isDefEq.respectTransparency false in
/-- ★★★逆向き——`1` は生成元へ戻る。 -/
theorem freeYonedaTopIso_inv_one {Y : Scheme.{0}} :
    (freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).inv.app (op (⊤ : Y.Opens))
        (1 : (Γ(Y, (⊤ : Y.Opens)) : Type))
      = ModuleCat.freeMk (𝟙 (⊤ : Y.Opens)) := by
  have hinj : Function.Injective
      ((freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).hom.app (op (⊤ : Y.Opens))).hom := by
    intro a b hab
    have h := congrArg
      ((freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).inv.app (op (⊤ : Y.Opens))).hom hab
    have ha := congrArg (fun (ψ : PresheafOfModules.freeObj
        (R := Y.presheaf ⋙ forget₂ CommRingCat RingCat) (yoneda.obj (⊤ : Y.Opens)) ⟶ _) =>
      ψ.app (op (⊤ : Y.Opens)) a) (freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).hom_inv_id
    have hb := congrArg (fun (ψ : PresheafOfModules.freeObj
        (R := Y.presheaf ⋙ forget₂ CommRingCat RingCat) (yoneda.obj (⊤ : Y.Opens)) ⟶ _) =>
      ψ.app (op (⊤ : Y.Opens)) b) (freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).hom_inv_id
    exact ha.symm.trans (h.trans hb)
  apply hinj
  rw [freeYonedaTopIso_hom_freeMk]
  exact congrArg (fun (ψ : 𝟙_ Y.PresheafOfModules ⟶ 𝟙_ Y.PresheafOfModules) =>
    ψ.app (op (⊤ : Y.Opens)) (1 : (Γ(Y, (⊤ : Y.Opens)) : Type)))
    (freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).inv_hom_id

/-! ## ★★★★★★★★`f^*𝒪_Y ≅ 𝒪_X` は `1` を `1` に送る -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★**引き戻した `1` は自由対象の側では unit の生成元での値である**。 -/
theorem pullSec_one_free {X Y : Scheme.{0}} (f : X ⟶ Y) :
    ((pullbackPre f).map (freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).inv).app
        (op ((Opens.map f.base).obj (⊤ : Y.Opens)))
        (pullSec f (𝟙_ Y.PresheafOfModules) ⊤ (1 : (Γ(Y, (⊤ : Y.Opens)) : Type)))
      = PresheafOfModules.freeYonedaEquiv
          (M := (PresheafOfModules.pushforward (pullbackPhi f)).obj
            ((pullbackPre f).obj (PresheafOfModules.freeObj
              (R := Y.presheaf ⋙ forget₂ CommRingCat RingCat) (yoneda.obj (⊤ : Y.Opens)))))
          (X := (⊤ : Y.Opens))
          ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app
            (PresheafOfModules.freeObj
              (R := Y.presheaf ⋙ forget₂ CommRingCat RingCat) (yoneda.obj (⊤ : Y.Opens)))) := by
  set Fr := PresheafOfModules.freeObj
    (R := Y.presheaf ⋙ forget₂ CommRingCat RingCat) (yoneda.obj (⊤ : Y.Opens)) with hFr
  set uF := (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app Fr with huF
  have hnat := congrArg
    (fun (ψ : 𝟙_ Y.PresheafOfModules ⟶ (PresheafOfModules.pushforward (pullbackPhi f)).obj
      ((pullbackPre f).obj Fr)) =>
        ψ.app (op (⊤ : Y.Opens)) (1 : (Γ(Y, (⊤ : Y.Opens)) : Type)))
    ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.naturality
      (freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).inv)
  have hsplit : ((((𝟭 (PresheafOfModules (Y.presheaf ⋙ forget₂ CommRingCat RingCat))).map
        (freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).inv ≫ uF)).app
        (op (⊤ : Y.Opens))) (1 : (Γ(Y, (⊤ : Y.Opens)) : Type))
      = uF.app (op (⊤ : Y.Opens))
        ((freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).inv.app (op (⊤ : Y.Opens))
          (1 : (Γ(Y, (⊤ : Y.Opens)) : Type))) := rfl
  refine hnat.symm.trans (hsplit.trans ?_)
  rw [freeYonedaTopIso_inv_one]
  exact (freeYonedaEquiv_apply_gen _).symm

set_option backward.isDefEq.respectTransparency false in
/-- ★前層加群の 3 つの合成を切断の水準で開く。★**`rfl`** である。 -/
theorem comp_app_apply₂ {C : Type} [SmallCategory C] {R : Cᵒᵖ ⥤ RingCat}
    {M N P Q : PresheafOfModules R} (a : M ⟶ N) (b : N ⟶ P) (c : P ⟶ Q)
    (Z : Cᵒᵖ) (x : (M.obj Z : Type)) :
    ((a ≫ b ≫ c).app Z) x = c.app Z (b.app Z (a.app Z x)) := rfl

set_option backward.isDefEq.respectTransparency false in
/-- ★★`f⁻¹⊤ = ⊤` を通して生成元を `1` へ送る最後の 2 段。 -/
theorem unitPre_last {X : Scheme.{0}} {U : X.Opens} (hU : U = ⊤)
    (hobj : PresheafOfModules.freeObj
        (R := X.presheaf ⋙ forget₂ CommRingCat RingCat) (yoneda.obj U)
      = PresheafOfModules.freeObj
        (R := X.presheaf ⋙ forget₂ CommRingCat RingCat) (yoneda.obj (⊤ : X.Opens)))
    (V : (X.Opens)ᵒᵖ) (y : (yoneda.obj U).obj V) :
    (((eqToIso hobj).hom
        ≫ (freeYonedaTopIso (R := X.presheaf) (α := X.Opens)).hom).app V)
        (ModuleCat.freeMk y) = (1 : ((X.presheaf.obj V) : Type)) := by
  subst hU
  exact freeYonedaTopIso_hom_freeMk V y

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★**`f^*𝒪_Y ≅ 𝒪_X` は引き戻した `1` を `1` に送る**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★機構は在庫の `isoHomUnitGenGen`（unit は生成元を生成元へ送る）である。 -/
theorem pullbackUnitPreIso_pullSec_one {X Y : Scheme.{0}} (f : X ⟶ Y) :
    (pullbackUnitPreIso f).hom.app (op ((Opens.map f.base).obj (⊤ : Y.Opens)))
        (pullSec f (𝟙_ Y.PresheafOfModules) ⊤ (1 : (Γ(Y, (⊤ : Y.Opens)) : Type)))
      = (1 : (Γ(X, (Opens.map f.base).obj (⊤ : Y.Opens)) : Type)) := by
  have hstep : (pullbackFreeYonedaIso f (⊤ : Y.Opens)).hom.app
        (op ((Opens.map f.base).obj (⊤ : Y.Opens)))
        (((pullbackPre f).map (freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).inv).app
          (op ((Opens.map f.base).obj (⊤ : Y.Opens)))
          (pullSec f (𝟙_ Y.PresheafOfModules) ⊤ (1 : (Γ(Y, (⊤ : Y.Opens)) : Type))))
      = ModuleCat.freeMk (𝟙 ((Opens.map f.base).obj (⊤ : Y.Opens))) := by
    rw [pullSec_one_free f]
    exact isoHomUnitGenGen (pullbackPhi f) (⊤ : Y.Opens)
  have hobj : PresheafOfModules.freeObj
      (R := X.presheaf ⋙ forget₂ CommRingCat RingCat)
        (yoneda.obj ((Opens.map f.base).obj (⊤ : Y.Opens)))
    = PresheafOfModules.freeObj
      (R := X.presheaf ⋙ forget₂ CommRingCat RingCat) (yoneda.obj (⊤ : X.Opens)) := by
    rw [opensMap_top]
  have hfun := congrArg
    (fun z => (((eqToIso hobj).hom
        ≫ (freeYonedaTopIso (R := X.presheaf) (α := X.Opens)).hom).app
      (op ((Opens.map f.base).obj (⊤ : Y.Opens)))) z) hstep
  simp only [pullbackUnitPreIso, Iso.trans_hom, Functor.mapIso_hom, Iso.symm_hom]
  refine (comp_app_apply₂ _ _ _ _ _).trans (hfun.trans ?_)
  exact unitPre_last (opensMap_top f) hobj _ _

/-! ## ★★★★★★★★単位元の等長性 -/

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★**輸送した基準の自明化と、引き戻した基準の自明化は同じ生成切断をもつ**。 -/
theorem trivGen_pullTriv_unit {X Y : Scheme.{0}} (f : X ⟶ Y) (V : X.Opens)
    (hV : V ≤ (Opens.map f.base).obj (⊤ : Y.Opens)) :
    trivGen ((pullbackPre f).obj (𝟙_ Y.PresheafOfModules)) V
        (pullTriv (pullbackUnitPreIso f) V (baseTriv X V))
      = trivGen ((pullbackPre f).obj (𝟙_ Y.PresheafOfModules)) V
          (pullTrivOfBase f (𝟙_ Y.PresheafOfModules) ⊤ (baseTriv Y ⊤) hV) := by
  have hL : trivGen ((pullbackPre f).obj (𝟙_ Y.PresheafOfModules)) V
      (pullTriv (pullbackUnitPreIso f) V (baseTriv X V))
      = (pullbackUnitPreIso f).inv.app (op V) (1 : (Γ(X, V) : Type)) :=
    trivGen_pullTriv (pullbackUnitPreIso f) V (baseTriv X V)
  have hR : trivGen ((pullbackPre f).obj (𝟙_ Y.PresheafOfModules)) V
      (pullTrivOfBase f (𝟙_ Y.PresheafOfModules) ⊤ (baseTriv Y ⊤) hV)
      = ((pullbackPre f).obj (𝟙_ Y.PresheafOfModules)).map (homOfLE hV).op
          (pullSec f (𝟙_ Y.PresheafOfModules) ⊤ (1 : (Γ(Y, (⊤ : Y.Opens)) : Type))) :=
    trivGen_pullTrivOfBase f (𝟙_ Y.PresheafOfModules) ⊤ (baseTriv Y ⊤) hV
  have hnat := PresheafOfModules.naturality_apply (pullbackUnitPreIso f).hom
    (homOfLE hV).op (pullSec f (𝟙_ Y.PresheafOfModules) ⊤ (1 : (Γ(Y, (⊤ : Y.Opens)) : Type)))
  have hone : (pullbackUnitPreIso f).hom.app (op V)
      (((pullbackPre f).obj (𝟙_ Y.PresheafOfModules)).map (homOfLE hV).op
        (pullSec f (𝟙_ Y.PresheafOfModules) ⊤ (1 : (Γ(Y, (⊤ : Y.Opens)) : Type))))
      = (1 : (Γ(X, V) : Type)) := by
    refine hnat.trans ?_
    show (𝟙_ X.PresheafOfModules).map (homOfLE hV).op
      ((pullbackUnitPreIso f).hom.app (op ((Opens.map f.base).obj (⊤ : Y.Opens)))
        (pullSec f (𝟙_ Y.PresheafOfModules) ⊤ (1 : (Γ(Y, (⊤ : Y.Opens)) : Type)))) = _
    rw [pullbackUnitPreIso_pullSec_one f]
    exact map_one (X.presheaf.map (homOfLE hV).op).hom
  rw [hL, hR, ← hone]
  exact congrArg (fun (ψ : (pullbackPre f).obj (𝟙_ Y.PresheafOfModules)
      ⟶ (pullbackPre f).obj (𝟙_ Y.PresheafOfModules)) => ψ.app (op V)
    (((pullbackPre f).obj (𝟙_ Y.PresheafOfModules)).map (homOfLE hV).op
      (pullSec f (𝟙_ Y.PresheafOfModules) ⊤ (1 : (Γ(Y, (⊤ : Y.Opens)) : Type)))))
    (pullbackUnitPreIso f).hom_inv_id

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★**`φ^*Ō_Y ≅ Ō_X`——等長**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★`V'' = V` でよい（`f⁻¹⊤ = ⊤` なので縮小が要らない）のが単位元の側の楽なところ。 -/
theorem isIsometry_pullback_one {X Y : Scheme.{0}} (f : X ⟶ Y) :
    IsIsometry (AMetricPullback f (AMetric.one Y)) (AMetric.one X) (pullbackUnitPreIso f) := by
  intro V e p hp
  have hV : V ≤ (Opens.map f.base).obj (⊤ : Y.Opens) := by
    rw [opensMap_top]; exact le_top
  have hchart : (AMetricPullback f (AMetric.one Y)).metric.h V
        (pullTriv (pullbackUnitPreIso f) V (baseTriv X V)) p
      = (AMetric.one X).metric.h V (baseTriv X V) p := by
    have h1 := LocalMetric.h_congr_trivGen (AMetricPullback f (AMetric.one Y)).metric V
      (pullTriv (pullbackUnitPreIso f) V (baseTriv X V))
      (pullTrivOfBase f (𝟙_ Y.PresheafOfModules) ⊤ (baseTriv Y ⊤) hV) p hp
      (trivGen_pullTriv_unit f V hV)
    have h2 : (AMetricPullback f (AMetric.one Y)).metric.h V
        (pullTrivOfBase f (𝟙_ Y.PresheafOfModules) ⊤ (baseTriv Y ⊤) hV) p
        = (structLocalMetric Y).h ⊤ (baseTriv Y ⊤) (p ≫ f) :=
      pullback_h_pullTrivOfBase f isLocallyTrivial_unit (structLocalMetric Y) ⊤
        (baseTriv Y ⊤) hV p hp
    have h3 : (structLocalMetric Y).h ⊤ (baseTriv Y ⊤) (p ≫ f) = 1 :=
      structLocalMetric_h_baseTriv Y ⊤ (p ≫ f) rfl
    have h4 : (AMetric.one X).metric.h V (baseTriv X V) p = 1 :=
      structLocalMetric_h_baseTriv X V p hp
    rw [h1, h2, h3, h4]
  have hL := (AMetricPullback f (AMetric.one Y)).metric.compat V
    (pullTriv (pullbackUnitPreIso f) V (baseTriv X V))
    (pullTriv (pullbackUnitPreIso f) V e) p hp
  have hR := (AMetric.one X).metric.compat V (baseTriv X V) e p hp
  have htu : transUnit (AMetricPullback f (AMetric.one Y)).sheaf V
        (pullTriv (pullbackUnitPreIso f) V (baseTriv X V))
        (pullTriv (pullbackUnitPreIso f) V e)
      = transUnit (AMetric.one X).sheaf V (baseTriv X V) e :=
    transUnit_pullTriv (pullbackUnitPreIso f) V (baseTriv X V) e
  rw [htu] at hL
  have hne : ‖evalOn p V hp (transUnit (AMetric.one X).sheaf V (baseTriv X V) e)‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p V hp (isUnit_transUnit _ V _ e))
  exact mul_right_cancel₀ hne (hL.trans (hchart.trans hR.symm))

/-! ## ★★★★★★★★等長同型は引き戻しで保たれる -/

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★**切断の引き戻しは加群の射について自然である**。

★これも `η_L` が**自然変換**であることそのものであり、ただで手に入る。 -/
theorem pullSec_naturality {X Y : Scheme.{0}} (f : X ⟶ Y) {A B : Y.PresheafOfModules}
    (ψ : A ⟶ B) (W : Y.Opens) (x : (A.obj (op W) : Type)) :
    pullSec f B W (ψ.app (op W) x)
      = ((pullbackPre f).map ψ).app (op ((Opens.map f.base).obj W)) (pullSec f A W x) :=
  congrArg
    (fun (χ : A ⟶ (PresheafOfModules.pushforward (pullbackPhi f)).obj
      ((pullbackPre f).obj B)) => χ.app (op W) x)
    ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.naturality ψ)

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★**同型で引き戻した自明化と、引き戻しの自明化は同じ生成切断をもつ**。 -/
theorem trivGen_pullTriv_pullTrivOfBase {X Y : Scheme.{0}} (f : X ⟶ Y)
    {A B : Y.PresheafOfModules} (φ : A ≅ B) (W : Y.Opens)
    (eW : (restrictPresheafFunctor Y W).obj B ≅ 𝟙_ (PresheafModulesOn Y W))
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W) :
    trivGen ((pullbackPre f).obj A) V'
        (pullTriv ((pullbackPre f).mapIso φ) V' (pullTrivOfBase f B W eW hV'W))
      = trivGen ((pullbackPre f).obj A) V' (pullTrivOfBase f A W (pullTriv φ W eW) hV'W) := by
  have hL : trivGen ((pullbackPre f).obj A) V'
      (pullTriv ((pullbackPre f).mapIso φ) V' (pullTrivOfBase f B W eW hV'W))
      = ((pullbackPre f).map φ.inv).app (op V')
        (((pullbackPre f).obj B).map (homOfLE hV'W).op (pullGen f B W eW)) := by
    rw [trivGen_pullTriv, trivGen_pullTrivOfBase]
    rfl
  have hR : trivGen ((pullbackPre f).obj A) V' (pullTrivOfBase f A W (pullTriv φ W eW) hV'W)
      = ((pullbackPre f).obj A).map (homOfLE hV'W).op
          (pullSec f A W (φ.inv.app (op W) (trivGen B W eW))) := by
    rw [trivGen_pullTrivOfBase]
    congr 1
  rw [hL, hR, pullSec_naturality f φ.inv W (trivGen B W eW)]
  exact PresheafOfModules.naturality_apply ((pullbackPre f).map φ.inv)
    (homOfLE hV'W).op (pullGen f B W eW)

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★**等長同型は引き戻しで保たれる**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★★これがないと引き戻しは**同値類の上に降りない**。 -/
theorem isIsometry_pullback {X Y : Scheme.{0}} (f : X ⟶ Y) {A B : AMetric Y}
    {φ : A.sheaf ≅ B.sheaf} (hφ : IsIsometry A B φ) :
    IsIsometry (AMetricPullback f A) (AMetricPullback f B) ((pullbackPre f).mapIso φ) := by
  intro V e p hp
  obtain ⟨c⟩ := nonempty_pullChart f B.triv V p hp
  have hpW : (p ≫ f) ⁻¹ᵁ c.W = ⊤ := comp_preimage_eq_top_of_le f c.hV'W c.hpV'
  set e₀ := pullTrivOfBase f B.sheaf c.W c.eW c.hV'W with he₀
  set e' := trivialOfLe ((pullbackPre f).obj B.sheaf) c.hV'V e with he'
  have hchart : (AMetricPullback f A).metric.h c.V' (pullTriv ((pullbackPre f).mapIso φ) c.V' e₀) p
      = (AMetricPullback f B).metric.h c.V' e₀ p := by
    have h1 := LocalMetric.h_congr_trivGen (AMetricPullback f A).metric c.V'
      (pullTriv ((pullbackPre f).mapIso φ) c.V' e₀)
      (pullTrivOfBase f A.sheaf c.W (pullTriv φ c.W c.eW) c.hV'W) p c.hpV'
      (trivGen_pullTriv_pullTrivOfBase f φ c.W c.eW c.hV'W)
    have h2 : (AMetricPullback f A).metric.h c.V'
        (pullTrivOfBase f A.sheaf c.W (pullTriv φ c.W c.eW) c.hV'W) p
        = A.metric.h c.W (pullTriv φ c.W c.eW) (p ≫ f) :=
      pullback_h_pullTrivOfBase f A.triv A.metric c.W (pullTriv φ c.W c.eW) c.hV'W p c.hpV'
    have h3 : (AMetricPullback f B).metric.h c.V' e₀ p = B.metric.h c.W c.eW (p ≫ f) :=
      pullback_h_pullTrivOfBase f B.triv B.metric c.W c.eW c.hV'W p c.hpV'
    rw [h1, h2, h3]
    exact hφ c.W c.eW (p ≫ f) hpW
  have hL := (AMetricPullback f A).metric.compat c.V'
    (pullTriv ((pullbackPre f).mapIso φ) c.V' e₀)
    (pullTriv ((pullbackPre f).mapIso φ) c.V' e') p c.hpV'
  have hR := (AMetricPullback f B).metric.compat c.V' e₀ e' p c.hpV'
  have htu : transUnit (AMetricPullback f A).sheaf c.V'
        (pullTriv ((pullbackPre f).mapIso φ) c.V' e₀)
        (pullTriv ((pullbackPre f).mapIso φ) c.V' e')
      = transUnit (AMetricPullback f B).sheaf c.V' e₀ e' :=
    transUnit_pullTriv ((pullbackPre f).mapIso φ) c.V' e₀ e'
  rw [htu] at hL
  have hne : ‖evalOn p c.V' c.hpV'
      (transUnit (AMetricPullback f B).sheaf c.V' e₀ e')‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p c.V' c.hpV' (isUnit_transUnit _ c.V' e₀ e'))
  have hkey : (AMetricPullback f A).metric.h c.V'
        (pullTriv ((pullbackPre f).mapIso φ) c.V' e') p
      = (AMetricPullback f B).metric.h c.V' e' p :=
    mul_right_cancel₀ hne (hL.trans (hchart.trans hR.symm))
  have hLres := (AMetricPullback f A).metric.restrict c.hV'V
    (pullTriv ((pullbackPre f).mapIso φ) V e) p c.hpV'
  have hRres := (AMetricPullback f B).metric.restrict c.hV'V e p c.hpV'
  exact hLres.symm.trans (hkey.trans hRres)

/-! ## ★★★★★★★★★★`APic` の群準同型 -/

theorem isometric_pullback {X Y : Scheme.{0}} (f : X ⟶ Y) {A B : AMetric Y}
    (h : Isometric A B) : Isometric (AMetricPullback f A) (AMetricPullback f B) :=
  h.elim fun φ hφ => ⟨(pullbackPre f).mapIso φ, isIsometry_pullback f hφ⟩

theorem isometric_pullback_one {X Y : Scheme.{0}} (f : X ⟶ Y) :
    Isometric (AMetricPullback f (1 : AMetric Y)) (1 : AMetric X) :=
  ⟨pullbackUnitPreIso f, isIsometry_pullback_one f⟩

/-- ★★★★★★★★**逆元つきの算術直線束の引き戻し**。

★逆であることは 3 段で出る: 積の等長性（§9-743）・等長性の保存・単位元の等長性。 -/
noncomputable def AInv.pullback {X Y : Scheme.{0}} (f : X ⟶ Y) (L : AInv Y) : AInv X where
  carrier := AMetricPullback f L.carrier
  inv := AMetricPullback f L.inv
  isInv :=
    isometric_trans (isometric_symm (isometric_pullback_mul f L.carrier L.inv))
      (isometric_trans (isometric_pullback f L.isInv) (isometric_pullback_one f))

/-- ★★★★★★★★★★**引き戻しは `APic` の群準同型である** `φ^* : APicM Y →* APicM X`。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★★原文が `evident` で畳んだ引き戻しは、`APic` が群である以上、
**群準同型**でなければ意味をなさない。★★★これがその中身である。 -/
noncomputable def APicMPullback {X Y : Scheme.{0}} (f : X ⟶ Y) : APicM Y →* APicM X where
  toFun := Quotient.map (AInv.pullback f) (fun _ _ h => isometric_pullback f h)
  map_one' := APicM.mk_eq_mk (isometric_pullback_one f)
  map_mul' := by
    refine Quotient.ind fun L => Quotient.ind fun M => ?_
    exact APicM.mk_eq_mk (isometric_pullback_mul f L.carrier M.carrier)

@[simp] theorem APicMPullback_mk {X Y : Scheme.{0}} (f : X ⟶ Y) (L : AInv Y) :
    APicMPullback f (APicM.mk L) = APicM.mk (AInv.pullback f L) := rfl

/-! ### ★出典の紐付け(`.src`) -/

def pullbackUnitPreIso_pullSec_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(f^*𝒪_Y ≅ 𝒪_X が引き戻した 1 を 1 に送ること)",
    sectionId := "genell-def-1-1-ii" }

def isIsometry_pullback_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(単位元の引き戻しが等長であること)",
    sectionId := "genell-def-1-1-ii" }

def isIsometry_pullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(等長同型が引き戻しで保たれること)",
    sectionId := "genell-def-1-1-ii" }

def AInv.pullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(逆元つきの算術直線束の引き戻し)",
    sectionId := "genell-def-1-1-ii" }

def APicMPullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻しが APic の群準同型であること)",
    sectionId := "genell-def-1-1-ii" }

def APicMPullback.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isometric_pullback_mul(引き戻しはテンソル積と可換、§9-743)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isometric_pullback_mul") 3,
    .citation "[ABC3]" "isoHomUnitGenGen(unit は生成元を生成元へ送る)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isoHomUnitGenGen") 3,
    .citation "[ABC3]" "pullbackUnitPreIso(前層の段の f^*𝒪_Y ≅ 𝒪_X)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullbackUnitPreIso") 3,
    .implicitStep
      ("★★手口は §9-741 / §9-743 と同じ——**同型の等式ではなく生成切断の一致**。" ++
       "単位元の側も等長性の保存も、要るのは pullSec の自然性(前層の射・加群の射の" ++
       "両方について)だけである") 3,
    .implicitStep
      ("★★★残っている段: 計量表示の ht_M̄(x) = deg_F(x_F^* M̄)。" ++
       "本ブロックで APicM Y →* APicM X が出たので x_F^* M̄ の類は取れるが、" ++
       "deg_F を APicM (Spec 𝓞_F) の上で読むには " ++
       "ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)([Szp] Prop 1.1、原文自身の引用)が要る") 3 ]

end ABC3.Found.Arakelov
