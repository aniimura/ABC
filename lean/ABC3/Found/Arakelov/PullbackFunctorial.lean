/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.MetricHeight
import ABC3.Meta.Claim

/-!
# **引き戻しは関手的である** `(φ ∘ ψ)^* ≅ ψ^* ∘ φ^*`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★★★★到達点

    `f^*(g^* M̄) ≅ (f ≫ g)^* M̄`   —— **等長**
    `APicMPullback f ∘ APicMPullback g = APicMPullback (f ≫ g)`

★原文の `ht_M̄(x)` が「`x` を与える**任意の**射」で定まるのは、
`deg_K(L|_{Spec 𝒪_K}) = deg_F(L)` があるからだが、
その式を**引き戻しの言葉で書く**には引き戻しの関手性が要る。
★★本ファイルがそれを無条件で与える。

## ★★★★★★★★層の段は mathlib にあった（実測 2026-08-28）

`PresheafOfModules.pullbackComp` が

    `pullback φ ⋙ pullback ψ ≅ pullback (φ ≫ whiskerLeft F.op ψ)`

を与える。★スキームの側で要るのは**添字が合うこと**だけであり、実測すると

| 主張 | 実測 |
|---|---|
| `Opens.map (f ≫ g).base = Opens.map g.base ⋙ Opens.map f.base` | ★**`rfl`** |
| `pullbackPhi (f ≫ g) = pullbackPhi g ≫ whiskerLeft _ (pullbackPhi f)` | ★**`rfl`** |

★★どちらも `rfl` だったので `pullbackPreComp` は 1 行である。

## ★★★★★★★★★★切断の段 —— `unit_conjugateEquiv_symm`

`pullSec` が合成と両立すること

    `(pullbackPreComp f g).hom (pullSec f (g^*L) (g⁻¹W) (pullSec g L W s))`
      `= pullSec (f ≫ g) L W s`

は、mathlib の **`unit_conjugateEquiv_symm`**（共役同型と unit の両立）そのものである。
★★★決め手は `PresheafOfModules.pushforwardComp` が **`Iso.refl`** であること
——押し出しは合成の上で**恒等**なので、共役の式の片側が消える。

## ★★★★★★あとはいつもの型 —— 生成切断で合わせる

`§9-741`・`§9-743`・`§9-745` と同じく、同型の等式ではなく生成切断の一致
（`trivGen_pullTriv_comp`）で `h` の一致（`LocalMetric.h_congr_trivGen`）を出す。

## ★★★★★★★★残っている仮定が 1 つに絞れた

    `htMetric_baseChange` は、[Szp] の同型の**族**が底変換と両立すること

      `∀ L, deg_K(ι^* L) = deg_F(L)`

だけを仮定して、原文の `deg_K(L|_{Spec 𝒪_K}) = deg_F(L)` を与える。
★因子の側ではこれは `degNormalized_baseChange` として**証明済み**である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace NumberField
open ABC3.Found.GenEll

/-! ## ★★★★★★層の段の関手性 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★**前層の引き戻しは関手的である** `g^* ⋙ f^* ≅ (f ≫ g)^*`。

★mathlib の `PresheafOfModules.pullbackComp` に添字を渡すだけである
——`Opens.map (f ≫ g).base = Opens.map g.base ⋙ Opens.map f.base` も
`pullbackPhi (f ≫ g) = pullbackPhi g ≫ whiskerLeft _ (pullbackPhi f)` も **`rfl`** だから。 -/
noncomputable def pullbackPreComp {X Y Z : Scheme.{0}} (f : X ⟶ Y) (g : Y ⟶ Z) :
    pullbackPre g ⋙ pullbackPre f ≅ pullbackPre (f ≫ g) :=
  PresheafOfModules.pullbackComp (pullbackPhi g) (pullbackPhi f)

/-! ## ★★★★★★★★★★切断の段 -/

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★**切断の引き戻しは合成と両立する**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★機構は mathlib の `unit_conjugateEquiv_symm`（共役同型と unit の両立）である。
★★決め手は `pushforwardComp` が **`Iso.refl`** であること——押し出しは合成の上で
恒等なので、共役の式の片側が消える。 -/
theorem pullSec_comp {X Y Z : Scheme.{0}} (f : X ⟶ Y) (g : Y ⟶ Z) (L : Z.PresheafOfModules)
    (W : Z.Opens) (s : (L.obj (op W) : Type)) :
    ((pullbackPreComp f g).hom.app L).app (op ((Opens.map (f ≫ g).base).obj W))
        (pullSec f ((pullbackPre g).obj L) ((Opens.map g.base).obj W) (pullSec g L W s))
      = pullSec (f ≫ g) L W s := by
  have h := unit_conjugateEquiv_symm
    (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi (f ≫ g)))
    ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi g)).comp
      (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)))
    (PresheafOfModules.pushforwardComp (pullbackPhi g) (pullbackPhi f)).inv L
  have h2 := congrArg (fun (χ : L ⟶ (PresheafOfModules.pushforward (pullbackPhi f)
      ⋙ PresheafOfModules.pushforward (pullbackPhi g)).obj
      ((pullbackPre (f ≫ g)).obj L)) => χ.app (op W) s) h
  exact h2.symm

/-! ## ★★★★★★生成切断で合わせる -/

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★**合成で輸送した自明化と、2 段で輸送した自明化は同じ生成切断をもつ**。 -/
theorem trivGen_pullTriv_comp {X Y Z : Scheme.{0}} (f : X ⟶ Y) (g : Y ⟶ Z)
    (L : Z.PresheafOfModules) (W : Z.Opens)
    (eW : (restrictPresheafFunctor Z W).obj L ≅ 𝟙_ (PresheafModulesOn Z W))
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map (f ≫ g).base).obj W) :
    trivGen ((pullbackPre f).obj ((pullbackPre g).obj L)) V'
        (pullTriv ((pullbackPreComp f g).app L) V'
          (pullTrivOfBase (f ≫ g) L W eW hV'W))
      = trivGen ((pullbackPre f).obj ((pullbackPre g).obj L)) V'
          (pullTrivOfBase f ((pullbackPre g).obj L) ((Opens.map g.base).obj W)
            (pullTrivOfBase g L W eW (le_refl _)) hV'W) := by
  set x := pullSec f ((pullbackPre g).obj L) ((Opens.map g.base).obj W)
    (pullSec g L W (trivGen L W eW)) with hxdef
  have hL : trivGen ((pullbackPre f).obj ((pullbackPre g).obj L)) V'
      (pullTriv ((pullbackPreComp f g).app L) V' (pullTrivOfBase (f ≫ g) L W eW hV'W))
      = ((pullbackPreComp f g).app L).inv.app (op V')
        (((pullbackPre (f ≫ g)).obj L).map (homOfLE hV'W).op
          (pullSec (f ≫ g) L W (trivGen L W eW))) := by
    rw [trivGen_pullTriv, trivGen_pullTrivOfBase]
    rfl
  have hR : trivGen ((pullbackPre f).obj ((pullbackPre g).obj L)) V'
      (pullTrivOfBase f ((pullbackPre g).obj L) ((Opens.map g.base).obj W)
        (pullTrivOfBase g L W eW (le_refl _)) hV'W)
      = ((pullbackPre f).obj ((pullbackPre g).obj L)).map (homOfLE hV'W).op x := by
    rw [trivGen_pullTrivOfBase]
    congr 1
    show pullSec f ((pullbackPre g).obj L) ((Opens.map g.base).obj W)
      (trivGen ((pullbackPre g).obj L) ((Opens.map g.base).obj W)
        (pullTrivOfBase g L W eW (le_refl _))) = _
    rw [pullGen_eq_trivGen]
    rfl
  have hx := PresheafOfModules.naturality_apply ((pullbackPreComp f g).app L).hom
    (homOfLE hV'W).op x
  have hcomp : ((pullbackPreComp f g).app L).hom.app
      (op ((Opens.map (f ≫ g).base).obj W)) x
      = pullSec (f ≫ g) L W (trivGen L W eW) := pullSec_comp f g L W (trivGen L W eW)
  rw [hcomp] at hx
  rw [hL, hR, ← hx]
  exact congrArg (fun (ψ : (pullbackPre g ⋙ pullbackPre f).obj L
      ⟶ (pullbackPre g ⋙ pullbackPre f).obj L) =>
    ψ.app (op V') (((pullbackPre f).obj ((pullbackPre g).obj L)).map (homOfLE hV'W).op x))
    ((pullbackPreComp f g).app L).hom_inv_id

/-! ## ★★★★★★★★★★計量の関手性 -/

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★**引き戻しは関手的である（等長）** `f^*(g^* M̄) ≅ (f ≫ g)^* M̄`。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★機構は 2 つの特徴づけ（`§9-742` の `pullback_h_pullTrivOfBase`）を継ぐだけ:

    `f^*(g^* m).h (pullTrivOfBase f (pullTrivOfBase g e_W))  = m.h W e_W ((p ≫ f) ≫ g)`
    `(f ≫ g)^* m.h (pullTrivOfBase (f ≫ g) e_W)             = m.h W e_W (p ≫ (f ≫ g))`

★★あとは `Category.assoc` である。 -/
theorem isIsometry_pullback_comp {X Y Z : Scheme.{0}} (f : X ⟶ Y) (g : Y ⟶ Z) (M : AMetric Z) :
    IsIsometry (AMetricPullback f (AMetricPullback g M)) (AMetricPullback (f ≫ g) M)
      ((pullbackPreComp f g).app M.sheaf) := by
  intro V e p hp
  obtain ⟨c⟩ := nonempty_pullChart (f ≫ g) M.triv V p hp
  have hpTop : p ⁻¹ᵁ ((Opens.map (f ≫ g).base).obj c.W) = ⊤ :=
    preimage_eq_top_of_le c.hV'W c.hpV'
  have hpfW : (p ≫ f) ⁻¹ᵁ ((Opens.map g.base).obj c.W) = ⊤ := hpTop
  set e₀ := pullTrivOfBase (f ≫ g) M.sheaf c.W c.eW c.hV'W with he₀
  set e' := trivialOfLe ((pullbackPre (f ≫ g)).obj M.sheaf) c.hV'V e with he'
  have hchart : (AMetricPullback f (AMetricPullback g M)).metric.h c.V'
        (pullTriv ((pullbackPreComp f g).app M.sheaf) c.V' e₀) p
      = (AMetricPullback (f ≫ g) M).metric.h c.V' e₀ p := by
    have h1 := LocalMetric.h_congr_trivGen (AMetricPullback f (AMetricPullback g M)).metric c.V'
      (pullTriv ((pullbackPreComp f g).app M.sheaf) c.V' e₀)
      (pullTrivOfBase f ((pullbackPre g).obj M.sheaf) ((Opens.map g.base).obj c.W)
        (pullTrivOfBase g M.sheaf c.W c.eW (le_refl _)) c.hV'W) p c.hpV'
      (trivGen_pullTriv_comp f g M.sheaf c.W c.eW c.hV'W)
    have h2 : (AMetricPullback f (AMetricPullback g M)).metric.h c.V'
        (pullTrivOfBase f ((pullbackPre g).obj M.sheaf) ((Opens.map g.base).obj c.W)
          (pullTrivOfBase g M.sheaf c.W c.eW (le_refl _)) c.hV'W) p
        = (AMetricPullback g M).metric.h ((Opens.map g.base).obj c.W)
            (pullTrivOfBase g M.sheaf c.W c.eW (le_refl _)) (p ≫ f) :=
      pullback_h_pullTrivOfBase f (AMetricPullback g M).triv (AMetricPullback g M).metric
        ((Opens.map g.base).obj c.W) (pullTrivOfBase g M.sheaf c.W c.eW (le_refl _))
        c.hV'W p c.hpV'
    have h3 : (AMetricPullback g M).metric.h ((Opens.map g.base).obj c.W)
        (pullTrivOfBase g M.sheaf c.W c.eW (le_refl _)) (p ≫ f)
        = M.metric.h c.W c.eW ((p ≫ f) ≫ g) :=
      pullback_h_pullTrivOfBase g M.triv M.metric c.W c.eW (le_refl _) (p ≫ f) hpfW
    have h4 : (AMetricPullback (f ≫ g) M).metric.h c.V' e₀ p
        = M.metric.h c.W c.eW (p ≫ (f ≫ g)) :=
      pullback_h_pullTrivOfBase (f ≫ g) M.triv M.metric c.W c.eW c.hV'W p c.hpV'
    rw [h1, h2, h3, h4, Category.assoc]
  have hL := (AMetricPullback f (AMetricPullback g M)).metric.compat c.V'
    (pullTriv ((pullbackPreComp f g).app M.sheaf) c.V' e₀)
    (pullTriv ((pullbackPreComp f g).app M.sheaf) c.V' e') p c.hpV'
  have hR := (AMetricPullback (f ≫ g) M).metric.compat c.V' e₀ e' p c.hpV'
  have htu : transUnit (AMetricPullback f (AMetricPullback g M)).sheaf c.V'
        (pullTriv ((pullbackPreComp f g).app M.sheaf) c.V' e₀)
        (pullTriv ((pullbackPreComp f g).app M.sheaf) c.V' e')
      = transUnit (AMetricPullback (f ≫ g) M).sheaf c.V' e₀ e' :=
    transUnit_pullTriv ((pullbackPreComp f g).app M.sheaf) c.V' e₀ e'
  rw [htu] at hL
  have hne : ‖evalOn p c.V' c.hpV'
      (transUnit (AMetricPullback (f ≫ g) M).sheaf c.V' e₀ e')‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p c.V' c.hpV' (isUnit_transUnit _ c.V' e₀ e'))
  have hkey : (AMetricPullback f (AMetricPullback g M)).metric.h c.V'
        (pullTriv ((pullbackPreComp f g).app M.sheaf) c.V' e') p
      = (AMetricPullback (f ≫ g) M).metric.h c.V' e' p :=
    mul_right_cancel₀ hne (hL.trans (hchart.trans hR.symm))
  have hLres := (AMetricPullback f (AMetricPullback g M)).metric.restrict c.hV'V
    (pullTriv ((pullbackPreComp f g).app M.sheaf) V e) p c.hpV'
  have hRres := (AMetricPullback (f ≫ g) M).metric.restrict c.hV'V e p c.hpV'
  exact hLres.symm.trans (hkey.trans hRres)

theorem isometric_pullback_comp {X Y Z : Scheme.{0}} (f : X ⟶ Y) (g : Y ⟶ Z) (M : AMetric Z) :
    Isometric (AMetricPullback f (AMetricPullback g M)) (AMetricPullback (f ≫ g) M) :=
  ⟨(pullbackPreComp f g).app M.sheaf, isIsometry_pullback_comp f g M⟩

theorem APicMPullback_mk_comp {X Y Z : Scheme.{0}} (f : X ⟶ Y) (g : Y ⟶ Z) (M : AInv Z) :
    APicMPullback f (APicMPullback g (APicM.mk M)) = APicMPullback (f ≫ g) (APicM.mk M) :=
  APicM.mk_eq_mk (isometric_pullback_comp f g M.carrier)

/-- ★★★★★★★★★★**`APic` の水準での関手性**。 -/
theorem APicMPullback_comp {X Y Z : Scheme.{0}} (f : X ⟶ Y) (g : Y ⟶ Z) (L : APicM Z) :
    APicMPullback f (APicMPullback g L) = APicMPullback (f ≫ g) L := by
  induction L using Quotient.ind with
  | _ M => exact APicMPullback_mk_comp f g M

/-! ## ★★★★★★★★高さの底変換不変性 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**高さは定義体の取り方に依らない**（計量表示）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★仮定は **[Szp] の同型の族が底変換と両立すること 1 つだけ**に絞られた
——引き戻しの関手性（本ファイル）が残りを吸収する。
★因子の側ではその両立は `degNormalized_baseChange` として**証明済み**である。 -/
theorem htMetric_baseChange {X : Scheme.{0}} (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (szpF : SzpData F) (szpK : SzpData K)
    (compat : ∀ L : APicM (Spec (CommRingCat.of (𝓞 F))),
      degMetric K szpK (APicMPullback
        (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L) = degMetric F szpF L)
    (M : AInv X) (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetric K szpK M (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF)
      = htMetric F szpF M xF := by
  show degMetric K szpK (APicMPullback
    (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) (APicM.mk M)) = _
  rw [← APicMPullback_comp (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) xF
    (APicM.mk M), compat]
  rfl

/-! ### ★出典の紐付け(`.src`) -/

def pullbackPreComp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(前層の引き戻しの関手性)",
    sectionId := "genell-def-1-1-ii" }

def pullSec_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(切断の引き戻しが合成と両立すること)",
    sectionId := "genell-def-1-1-ii" }

def isIsometry_pullback_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻しが関手的であること——等長)",
    sectionId := "genell-def-1-1-ii" }

def APicMPullback_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(APic の水準での引き戻しの関手性)",
    sectionId := "genell-def-1-1-ii" }

def htMetric_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(高さの底変換不変性——[Szp] の族の両立を仮定)",
    sectionId := "genell-def-1-1-ii" }

def isIsometry_pullback_comp.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "PresheafOfModules.pullbackComp(前層の引き戻しの合成)"
      (.inMathlib "PresheafOfModules.pullbackComp") 3,
    .citation "[mathlib]" "unit_conjugateEquiv_symm(共役同型と unit の両立)"
      (.inMathlib "CategoryTheory.unit_conjugateEquiv_symm") 3,
    .citation "[ABC3]" "pullback_h_pullTrivOfBase(引き戻した計量の特徴づけ、§9-742)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullback_h_pullTrivOfBase") 3,
    .implicitStep
      ("★実測(2026-08-28): Opens.map (f ≫ g).base = Opens.map g.base ⋙ Opens.map f.base も " ++
       "pullbackPhi (f ≫ g) = pullbackPhi g ≫ whiskerLeft _ (pullbackPhi f) も **rfl**。" ++
       "★★決め手は pushforwardComp が **Iso.refl** であること——" ++
       "押し出しは合成の上で恒等なので共役の式の片側が消える") 3,
    .implicitStep
      ("★★★残っている仮定が 1 つに絞れた: [Szp] の同型の**族**が底変換と両立すること。" ++
       "因子の側では degNormalized_baseChange として証明済みである") 3 ]

end ABC3.Found.Arakelov
