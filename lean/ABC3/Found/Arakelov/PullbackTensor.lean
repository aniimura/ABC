/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PullbackMetric
import ABC3.Found.Arakelov.AMetricIso
import ABC3.Found.Arakelov.PicDeltaIso
import ABC3.Found.Arakelov.PicPushMu
import ABC3.Meta.Claim

/-!
# **引き戻しはテンソル積と可換である**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★★★★到達点

    `φ^*(L̄ ⊗ M̄) ≅ φ^*L̄ ⊗ φ^*M̄`   —— **等長**（`isIsometry_pullback_mul`）

★層の水準の同型は在庫（`pullbackTensorIso`）。★★本ファイルが足すのは
**計量が一致すること**であり、これが `ht_M̄(x) = deg_F(x_F^* M̄)` の
**加法性**に要る段である。

## ★★★★★★★★機構 —— 特徴づけが全部やる

`§9-742` の `pullback_h_pullTrivOfBase`

    `(f^* m).h V' (pullTrivOfBase f L W e_W) p = m.h W e_W (p ≫ f)`

と `isTensorOf_tensor`

    `(m_A ⊗ m_B).h V (tensorTriv e_A e_B) p = m_A.h V e_A p · m_B.h V e_B p`

を突き合わせると、`p` の周りのチャートの上で**両辺とも**

    `m_A.h W e_A (p ≫ f) · m_B.h W e_B (p ≫ f)`

になる。★★★あとはチャートから任意の `(V, e)` へ `compat`・`restrict` で運ぶだけ。

## ★★★★★★★★★★再び —— **同型ではなく生成切断で合わせる**

チャートの上で 2 つの自明化

    `pullTriv δ (tensorTriv (pullTrivOfBase e_A) (pullTrivOfBase e_B))`  と
    `pullTrivOfBase (tensorTriv e_A e_B)`

を突き合わせる必要がある。★**同型として等しいことは示さない**——
`§9-741` と同じく**生成切断が一致すること**を示せば足りる
（`LocalMetric.h_congr_trivGen`：生成切断が同じなら `h` も同じ）。

★★その生成切断の一致は、結局

    `δ (pullSec f (A ⊗ B) W (a ⊗ b)) = pullSec f A W a ⊗ pullSec f B W b`

に落ちる（`delta_pullSec`）。★★★これは **`δ` の定義そのもの**から出る:
mathlib の `leftAdjointOplaxMonoidal_δ` は

    `δ = (adj.homEquiv).symm ((η ⊗ₘ η) ≫ μ)`

なので、`Adjunction.homEquiv_unit` で転置を戻せば
`η_{A⊗B} ≫ push δ = (η_A ⊗ₘ η_B) ≫ μ` であり、
切断に当てて在庫の `pushforward_μ_app_tmul`（`μ` は純テンソルの上で恒等）を使えばよい。

★★★★**ここでも「同型の水準 → 切断の水準」が効いている**（`PullbackGen.lean` と同じ手口）。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll

/-! ## ★★★★★★★★`δ` は切断の上で純テンソルを純テンソルに送る -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★**比較射 `δ` と切断の引き戻しは可換である**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

    `δ (pullSec f (A ⊗ B) W (a ⊗ b)) = pullSec f A W a ⊗ pullSec f B W b`

★機構は `δ` の**定義**（`leftAdjointOplaxMonoidal_δ`）を随伴で転置し、
在庫の `pushforward_μ_app_tmul`（`μ` は純テンソルの上で恒等）を当てるだけ。 -/
theorem delta_pullSec {X Y : Scheme.{0}} (f : X ⟶ Y) (A B : Y.PresheafOfModules)
    (W : Y.Opens) (a : (A.obj (op W) : Type)) (b : (B.obj (op W) : Type)) :
    (pullbackDelta f A B).app (op ((Opens.map f.base).obj W))
        (pullSec f (A ⊗ B) W (a ⊗ₜ[(Y.presheaf.obj (op W) : Type)] b))
      = (pullSec f A W a) ⊗ₜ[(X.presheaf.obj (op ((Opens.map f.base).obj W)) : Type)]
          (pullSec f B W b) := by
  have hu : (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).homEquiv _ _
        (pullbackDelta f A B)
      = (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app (A ⊗ B)
        ≫ (PresheafOfModules.pushforward (pullbackPhi f)).map (pullbackDelta f A B) :=
    Adjunction.homEquiv_unit _ _ _ (pullbackDelta f A B)
  have hd : (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).homEquiv _ _
        (pullbackDelta f A B)
      = (((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app A
          ⊗ₘ (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app B)
        ≫ Functor.LaxMonoidal.μ (PresheafOfModules.pushforward (pullbackPhi f)) _ _) := by
    rw [pullbackDelta, Adjunction.leftAdjointOplaxMonoidal_δ, Equiv.apply_symm_apply]
  have h := congrArg
    (fun (φ : (A ⊗ B) ⟶ (PresheafOfModules.pushforward (pullbackPhi f)).obj
      ((pullbackPre f).obj A ⊗ (pullbackPre f).obj B)) =>
        φ.app (op W) (a ⊗ₜ[(Y.presheaf.obj (op W) : Type)] b)) (hu.symm.trans hd)
  refine h.trans ?_
  exact pushforward_μ_app_tmul f ((pullbackPre f).obj A) ((pullbackPre f).obj B) W
    (pullSec f A W a) (pullSec f B W b)

/-! ## ★★★★★生成切断の計算 -/

/-- ★★★**テンソル自明化の生成切断は生成切断のテンソルである**。 -/
theorem trivGen_tensorTriv {X : Scheme.{0}} {A B : X.PresheafOfModules} {V : X.Opens}
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V)) :
    trivGen (A ⊗ B) V (tensorTriv eA eB)
      = (trivGen A V eA) ⊗ₜ[(Γ(X, V) : Type)] (trivGen B V eB) := by
  apply (trivEquiv (A ⊗ B) V (tensorTriv eA eB)).injective
  show (trivEquiv (A ⊗ B) V (tensorTriv eA eB)) ((trivEquiv (A ⊗ B) V (tensorTriv eA eB)).symm 1)
    = _
  rw [LinearEquiv.apply_symm_apply, trivEquiv_tensorTriv eA eB (trivGen A V eA) (trivGen B V eB)]
  show (1 : (Γ(X, V) : Type))
    = (trivEquiv A V eA) ((trivEquiv A V eA).symm 1)
      * (trivEquiv B V eB) ((trivEquiv B V eB).symm 1)
  rw [LinearEquiv.apply_symm_apply, LinearEquiv.apply_symm_apply, mul_one]

/-- ★★★**同型で引き戻した自明化の生成切断は、生成切断を逆に送ったものである**。 -/
theorem trivGen_pullTriv {X : Scheme.{0}} {L M : X.PresheafOfModules} (φ : L ≅ M) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) :
    trivGen L V (pullTriv φ V e) = φ.inv.app (op V) (trivGen M V e) := by
  apply (trivEquiv L V (pullTriv φ V e)).injective
  show (trivEquiv L V (pullTriv φ V e)) ((trivEquiv L V (pullTriv φ V e)).symm 1) = _
  rw [LinearEquiv.apply_symm_apply, trivEquiv_pullTriv φ V e (φ.inv.app (op V) (trivGen M V e))]
  have hiv : φ.hom.app (op V) (φ.inv.app (op V) (trivGen M V e)) = trivGen M V e :=
    congrArg (fun (ψ : M ⟶ M) => ψ.app (op V) (trivGen M V e)) φ.inv_hom_id
  rw [hiv]
  show (1 : (Γ(X, V) : Type)) = (trivEquiv M V e) ((trivEquiv M V e).symm 1)
  rw [LinearEquiv.apply_symm_apply]

set_option backward.isDefEq.respectTransparency false in
/-- ★前層のテンソル積の制限は純テンソルの上で成分ごとである。★**`rfl`** である。 -/
theorem tensor_map_tmul {X : Scheme.{0}} (A B : X.PresheafOfModules) {V W : X.Opens}
    (hWV : W ≤ V) (x : (A.obj (op V) : Type)) (y : (B.obj (op V) : Type)) :
    (A ⊗ B).map (homOfLE hWV).op (x ⊗ₜ[(Γ(X, V) : Type)] y)
      = (A.map (homOfLE hWV).op x) ⊗ₜ[(Γ(X, W) : Type)] (B.map (homOfLE hWV).op y) := rfl

/-! ## ★★★★★★★★2 つの自明化の生成切断が一致すること -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★**輸送したテンソル自明化と、テンソル自明化の輸送は、同じ生成切断をもつ**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★★**同型として等しいことは主張していない**——`§9-741` と同じく、
`transUnit`（したがって `h`）が使うのは生成切断だけなので、これで足りる。 -/
theorem trivGen_pullTriv_tensorTriv {X Y : Scheme.{0}} (f : X ⟶ Y) (A B : Y.PresheafOfModules)
    (W : Y.Opens)
    (eA : (restrictPresheafFunctor Y W).obj A ≅ 𝟙_ (PresheafModulesOn Y W))
    (eB : (restrictPresheafFunctor Y W).obj B ≅ 𝟙_ (PresheafModulesOn Y W))
    {V'' : X.Opens} (hV''W : V'' ≤ (Opens.map f.base).obj W) :
    trivGen ((pullbackPre f).obj (A ⊗ B)) V''
        (pullTriv (pullbackTensorIso f A B) V''
          (tensorTriv (pullTrivOfBase f A W eA hV''W) (pullTrivOfBase f B W eB hV''W)))
      = trivGen ((pullbackPre f).obj (A ⊗ B)) V''
          (pullTrivOfBase f (A ⊗ B) W (tensorTriv eA eB) hV''W) := by
  have hg : pullGen f (A ⊗ B) W (tensorTriv eA eB)
      = pullSec f (A ⊗ B) W ((trivGen A W eA) ⊗ₜ[(Γ(Y, W) : Type)] (trivGen B W eB)) := by
    show pullSec f (A ⊗ B) W (trivGen (A ⊗ B) W (tensorTriv eA eB)) = _
    rw [trivGen_tensorTriv]
  have hmid : (pullbackTensorIso f A B).hom.app (op V'')
        (((pullbackPre f).obj (A ⊗ B)).map (homOfLE hV''W).op
          (pullGen f (A ⊗ B) W (tensorTriv eA eB)))
      = (((pullbackPre f).obj A).map (homOfLE hV''W).op (pullGen f A W eA))
        ⊗ₜ[(Γ(X, V'') : Type)]
        (((pullbackPre f).obj B).map (homOfLE hV''W).op (pullGen f B W eB)) := by
    rw [PresheafOfModules.naturality_apply (pullbackTensorIso f A B).hom
      (homOfLE hV''W).op (pullGen f (A ⊗ B) W (tensorTriv eA eB))]
    show ((pullbackPre f).obj A ⊗ (pullbackPre f).obj B).map (homOfLE hV''W).op
      ((pullbackDelta f A B).app (op ((Opens.map f.base).obj W))
        (pullGen f (A ⊗ B) W (tensorTriv eA eB))) = _
    rw [hg, delta_pullSec f A B W (trivGen A W eA) (trivGen B W eB),
      tensor_map_tmul ((pullbackPre f).obj A) ((pullbackPre f).obj B) hV''W]
    rfl
  have hcancel : (pullbackTensorIso f A B).inv.app (op V'')
      ((pullbackTensorIso f A B).hom.app (op V'')
        (((pullbackPre f).obj (A ⊗ B)).map (homOfLE hV''W).op
          (pullGen f (A ⊗ B) W (tensorTriv eA eB))))
      = ((pullbackPre f).obj (A ⊗ B)).map (homOfLE hV''W).op
          (pullGen f (A ⊗ B) W (tensorTriv eA eB)) :=
    congrArg (fun (ψ : (pullbackPre f).obj (A ⊗ B) ⟶ (pullbackPre f).obj (A ⊗ B)) =>
      ψ.app (op V'') (((pullbackPre f).obj (A ⊗ B)).map (homOfLE hV''W).op
        (pullGen f (A ⊗ B) W (tensorTriv eA eB)))) (pullbackTensorIso f A B).hom_inv_id
  rw [trivGen_pullTriv, trivGen_tensorTriv, trivGen_pullTrivOfBase, trivGen_pullTrivOfBase,
    trivGen_pullTrivOfBase]
  exact (congrArg (fun z => (pullbackTensorIso f A B).inv.app (op V'') z) hmid.symm).trans hcancel

/-! ## ★★★★★★生成切断が同じなら計量の値も同じ -/

/-- ★★★★★★**生成切断が一致する 2 つの自明化は、同じ `h` の値をもつ**。

★★これが「同型の等式を避けて生成切断で済ませる」ことを可能にしている補題である
——`transUnit F V e e' = trivEquiv e' (trivGen e)` なので、
生成切断が同じなら遷移単元は `1` であり、`compat` がそのまま結論を与える。 -/
theorem LocalMetric.h_congr_trivGen {X : Scheme.{0}} {F : X.PresheafOfModules}
    (m : LocalMetric X F) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (hgen : trivGen F V e = trivGen F V e') : m.h V e p = m.h V e' p := by
  have h1 : transUnit F V e e' = 1 := by
    rw [transUnit_eq_trivEquiv_trivGen, hgen, ← transUnit_eq_trivEquiv_trivGen, transUnit_self]
  have h2 := m.compat V e e' p hp
  rw [h1, evalOn_one, norm_one, mul_one] at h2
  exact h2.symm

/-! ## ★★★★★★★★チャートの上での一致 -/

/-- ★★★★★★★★**チャートの上では両辺が一致する**。

★★機構は 2 つの特徴づけを突き合わせるだけ——両辺とも

    `m_A.h W e_A (p ≫ f) · m_B.h W e_B (p ≫ f)`

になる。 -/
theorem pullback_tensor_h_chart {X Y : Scheme.{0}} (f : X ⟶ Y) {A B : Y.PresheafOfModules}
    (hA : IsLocallyTrivial Y A) (hB : IsLocallyTrivial Y B)
    (mA : LocalMetric Y A) (mB : LocalMetric Y B) (W : Y.Opens)
    (eA : (restrictPresheafFunctor Y W).obj A ≅ 𝟙_ (PresheafModulesOn Y W))
    (eB : (restrictPresheafFunctor Y W).obj B ≅ 𝟙_ (PresheafModulesOn Y W))
    {V'' : X.Opens} (hV''W : V'' ≤ (Opens.map f.base).obj W)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V'' = ⊤) :
    (LocalMetric.pullback f (hA.tensor hB) (LocalMetric.tensor hA hB mA mB)).h V''
        (pullTriv (pullbackTensorIso f A B) V''
          (tensorTriv (pullTrivOfBase f A W eA hV''W) (pullTrivOfBase f B W eB hV''W))) p
      = (LocalMetric.tensor (isLocallyTrivial_pullbackPre f A hA)
          (isLocallyTrivial_pullbackPre f B hB)
          (LocalMetric.pullback f hA mA) (LocalMetric.pullback f hB mB)).h V''
          (tensorTriv (pullTrivOfBase f A W eA hV''W) (pullTrivOfBase f B W eB hV''W)) p := by
  have hpW : (p ≫ f) ⁻¹ᵁ W = ⊤ := comp_preimage_eq_top_of_le f hV''W hp
  rw [LocalMetric.h_congr_trivGen _ V'' _
      (pullTrivOfBase f (A ⊗ B) W (tensorTriv eA eB) hV''W) p hp
      (trivGen_pullTriv_tensorTriv f A B W eA eB hV''W),
    pullback_h_pullTrivOfBase f (hA.tensor hB) (LocalMetric.tensor hA hB mA mB) W
      (tensorTriv eA eB) hV''W p hp,
    isTensorOf_tensor hA hB mA mB W eA eB (p ≫ f) hpW,
    isTensorOf_tensor (isLocallyTrivial_pullbackPre f A hA)
      (isLocallyTrivial_pullbackPre f B hB)
      (LocalMetric.pullback f hA mA) (LocalMetric.pullback f hB mB) V''
      (pullTrivOfBase f A W eA hV''W) (pullTrivOfBase f B W eB hV''W) p hp,
    pullback_h_pullTrivOfBase f hA mA W eA hV''W p hp,
    pullback_h_pullTrivOfBase f hB mB W eB hV''W p hp]

/-! ## ★★★★★★★★★★到達点 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★**引き戻しはテンソル積と可換である（等長）**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

    `φ^*(L̄ ⊗ M̄) ≅ φ^*L̄ ⊗ φ^*M̄`   —— 等長

★★これが `ht_M̄(x) = deg_F(x_F^* M̄)` の**加法性**に要る段である。

★機構: `p` の周りに `Y` 側のチャート `(W, e_A, e_B)` を取り、
`V'' := V ⊓ f⁻¹W` へ降りて `pullback_tensor_h_chart` を使い、
`compat`（遷移単元は `pullTriv` で変わらない＝`transUnit_pullTriv`）と
`restrict` で任意の `(V, e)` へ戻す。 -/
theorem isIsometry_pullback_mul {X Y : Scheme.{0}} (f : X ⟶ Y) (L M : AMetric Y) :
    IsIsometry (AMetricPullback f (L * M)) (AMetricPullback f L * AMetricPullback f M)
      (pullbackTensorIso f L.sheaf M.sheaf) := by
  intro V e p hp
  obtain ⟨c⟩ := nonempty_tensorChart L.triv M.triv ⊤ (p ≫ f) rfl
  set V'' : X.Opens := V ⊓ (Opens.map f.base).obj c.W with hV''def
  have hV''V : V'' ≤ V := inf_le_left
  have hV''W : V'' ≤ (Opens.map f.base).obj c.W := inf_le_right
  have hpfW : p ⁻¹ᵁ ((Opens.map f.base).obj c.W) = ⊤ := c.hpW
  have hp'' : p ⁻¹ᵁ V'' = ⊤ := by
    show p ⁻¹ᵁ V ⊓ p ⁻¹ᵁ ((Opens.map f.base).obj c.W) = ⊤
    rw [hp, hpfW, inf_idem]
  set e₀ := tensorTriv (pullTrivOfBase f L.sheaf c.W c.eA hV''W)
    (pullTrivOfBase f M.sheaf c.W c.eB hV''W) with he₀
  set e' := trivialOfLe ((pullbackPre f).obj L.sheaf ⊗ (pullbackPre f).obj M.sheaf) hV''V e
    with he'
  have hchart := pullback_tensor_h_chart f L.triv M.triv L.metric M.metric c.W c.eA c.eB
    hV''W p hp''
  have hL := (AMetricPullback f (L * M)).metric.compat V''
    (pullTriv (pullbackTensorIso f L.sheaf M.sheaf) V'' e₀)
    (pullTriv (pullbackTensorIso f L.sheaf M.sheaf) V'' e') p hp''
  have hR := (AMetricPullback f L * AMetricPullback f M).metric.compat V'' e₀ e' p hp''
  have htu : transUnit (AMetricPullback f (L * M)).sheaf V''
        (pullTriv (pullbackTensorIso f L.sheaf M.sheaf) V'' e₀)
        (pullTriv (pullbackTensorIso f L.sheaf M.sheaf) V'' e')
      = transUnit (AMetricPullback f L * AMetricPullback f M).sheaf V'' e₀ e' :=
    transUnit_pullTriv (pullbackTensorIso f L.sheaf M.sheaf) V'' e₀ e'
  rw [htu] at hL
  have hne : ‖evalOn p V'' hp''
      (transUnit (AMetricPullback f L * AMetricPullback f M).sheaf V'' e₀ e')‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p V'' hp'' (isUnit_transUnit _ V'' e₀ e'))
  have hkey : (AMetricPullback f (L * M)).metric.h V''
        (pullTriv (pullbackTensorIso f L.sheaf M.sheaf) V'' e') p
      = (AMetricPullback f L * AMetricPullback f M).metric.h V'' e' p :=
    mul_right_cancel₀ hne (hL.trans (hchart.trans hR.symm))
  have hLres := (AMetricPullback f (L * M)).metric.restrict hV''V
    (pullTriv (pullbackTensorIso f L.sheaf M.sheaf) V e) p hp''
  have hRres := (AMetricPullback f L * AMetricPullback f M).metric.restrict hV''V e p hp''
  exact hLres.symm.trans (hkey.trans hRres)

/-- ★★★★★★★★★★**同値類の水準での形**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★これが `APicM Y → APicM X` が**群準同型**であることの中身である。 -/
theorem isometric_pullback_mul {X Y : Scheme.{0}} (f : X ⟶ Y) (L M : AMetric Y) :
    Isometric (AMetricPullback f (L * M)) (AMetricPullback f L * AMetricPullback f M) :=
  ⟨pullbackTensorIso f L.sheaf M.sheaf, isIsometry_pullback_mul f L M⟩

/-! ### ★出典の紐付け(`.src`) -/

def delta_pullSec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(比較射 δ と切断の引き戻しが可換であること)",
    sectionId := "genell-def-1-1-ii" }

def LocalMetric.h_congr_trivGen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(生成切断が一致すれば計量の値も一致すること)",
    sectionId := "genell-def-1-1-ii" }

def isIsometry_pullback_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻しがテンソル積と可換であること——等長)",
    sectionId := "genell-def-1-1-ii" }

def isometric_pullback_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻しの乗法性——同値類の水準)",
    sectionId := "genell-def-1-1-ii" }

def isIsometry_pullback_mul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "pullback_h_pullTrivOfBase(引き戻した計量の特徴づけ、§9-742)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullback_h_pullTrivOfBase") 3,
    .citation "[ABC3]" "isTensorOf_tensor(テンソル計量の特徴づけ)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isTensorOf_tensor") 3,
    .citation "[ABC3]" "pullbackTensorIso(層の水準で引き戻しがテンソル積を保つこと)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullbackTensorIso") 3,
    .citation "[ABC3]" "pushforward_μ_app_tmul(押し出しの μ が純テンソルの上で恒等)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pushforward_μ_app_tmul") 3,
    .citation "[mathlib]" "Adjunction.leftAdjointOplaxMonoidal_δ(δ の定義)"
      (.inMathlib "CategoryTheory.Adjunction.leftAdjointOplaxMonoidal_δ") 3,
    .implicitStep
      ("★★手口は §9-741 と同じ——**同型の等式ではなく生成切断の一致**で済ませる。" ++
       "h_congr_trivGen が「生成切断が同じなら h も同じ」を与えるので、" ++
       "pullTriv δ (tensorTriv (pullTrivOfBase …)) と pullTrivOfBase (tensorTriv …) が" ++
       "同型として等しいことは示さなくてよい") 3,
    .implicitStep
      ("★★★残っている段: 単位元の側 f^*Ō_Y ≅ Ō_X が等長であること。" ++
       "これと本定理を合わせると APicM Y → APicM X が群準同型になる") 3 ]

end ABC3.Found.Arakelov
