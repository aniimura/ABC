/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PullbackUnitEnd
import ABC3.Found.Arakelov.LocalMetric
import ABC3.Meta.Claim

/-!
# **生成切断の引き戻し**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★★★★本ファイルが外す障害

`Definition 1.1 (ii)` の計量を運ぶ段でチャート

    `(W, e_W)`（`W` は `Y` の開集合、`e_W : L|_W ≅ 𝟙_`）

を使うと、**`W` の取り方に依らないこと**を示す必要がある。
`§9-740` の `transUnit_pullTrivOfBase` は**同じ `W`** の上の話なので、
`W` を小さくする段（`W'' ≤ W`）が別に要る。

★素朴には `pullTrivOfBase` そのものが `W` の縮小と可換であること、すなわち

    `pullTrivOfBase f L W'' (e_W|_{W''}) = pullTrivOfBase f L W e_W`（`V''` 上）

を示すことになる。★★しかしこれは **`bcIso` の `W` についての自然性**を要し、
`bcMate → bcMate が同型 → bcIso` の三部作をもう一度やることになる
（`restrictOnFunctor` と `pullbackPreOn` が可換であることの Beck–Chevalley）。
★★★**実測: `rfl` ではない**（2026-08-28）。

## ★★★★★★★★★★代わりに通した道 —— **切断の水準へ落とす**

同型の水準の自然性ではなく、**随伴の unit を切断に適用したもの**を使う:

    `pullSec f L W : L(W) → (f^*L)(f⁻¹W)`,  `pullSec := η_L.app (op W)`

★これは**前層加群の射の成分**なので

| 要る性質 | 出どころ | 値段 |
|---|---|---|
| `pullSec (u • g) = ρ(u) • pullSec g` | `η_L` が加群の射（係数制限） | ★ただ |
| `pullSec (g|_{W''}) = (pullSec g)|_{f⁻¹W''}` | `η_L` が**前層の射**（自然性） | ★ただ |

★★同型の水準の自然性（高い）が、切断の水準の自然性（ただ）に化ける。
**これが本ファイルの手口である。**

## ★★★★★★★★★★橋 —— `trivEquiv_pullSec`

    `trivEquiv (f^*L) (f⁻¹W) (pullTrivOfBase f L W e_W) (pullSec f L W s)`
      `= ρ (trivEquiv L W e_W s)`

★「輸送した自明化は、引き戻した切断を、座標を引き戻して読む」。

★★機構は 4 段:

| 段 | 内容 | 道具 |
|---|---|---|
| 1 | `res_W η_L = η^{on}_{res_W L} ≫ push(bcMate.app L)` | `bcMate` の**定義**（転置）＋ `Adjunction.homEquiv_unit` |
| 2 | `bcIso.inv ∘ bcMate.app L = 𝟙` | `Iso.hom_inv_id` |
| 3 | `pOn(e_W.hom) ∘ η^{on}_{res L} = η^{on}_{𝟙_} ∘ e_W.hom` | `η^{on}` の**自然性** |
| 4 | `pullbackOnUnitIso.hom ∘ η^{on}_{𝟙_} = ρ` | `etaOn_app`（`§9-740` の `etaOn_unitOne`） |

★★★段 2 が効くのは `bcMate` が「`res_W η_L` の転置」として**定義されている**から。
つまり `bcIso` の中身を開かずに、`bcIso` の**作り方**だけで消える。

## ★★★★★★★★★★到達点

    `transUnit (f^*L) V'' (pullTrivOfBase f L W'' (e_W|_{W''})) e`
      `= transUnit (f^*L) V'' (pullTrivOfBase f L W e_W) e`

★★これで `W` の縮小が通る。★★★注意: **同型そのものの等式は示していない**
（`pullTrivOfBase_shrinkW` は未証明のまま）。示したのは
**生成切断が一致すること**であり、`transUnit` はそれだけで決まる。
★計量が要るのは `transUnit` だけなので、これで足りる。

## ★残っている段（明示）

★★`PullChart` と `pullChartH` の 3 段（`triv_indep` / `shrinkV` / `shrinkW`）、
そして `LocalMetric.pullback` と `AMetric` の引き戻し。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll

/-! ## ★★★★★★随伴の unit を切断に当てる -/

/-- ★★★★★★**切断の引き戻し** `L(W) → (f^*L)(f⁻¹W)`。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★引き戻し–押し出し随伴の **unit の `op W` 成分**そのものである。 -/
noncomputable def pullSec {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (g : (L.obj (op W) : Type)) :
    (((pullbackPre f).obj L).obj (op ((Opens.map f.base).obj W)) : Type) :=
  ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app L).app (op W) g

set_option backward.isDefEq.respectTransparency false in
/-- ★★★**切断の引き戻しは `ρ`-半線型である**。

★`η_L` は加群の射で、行き先は**係数制限**なので `u • x` が `ρ(u) • x` になる。 -/
theorem pullSec_smul {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (u : (Γ(Y, W) : Type)) (g : (L.obj (op W) : Type)) :
    pullSec f L W (u • g) = (f.c.app (op W)) u • pullSec f L W g := by
  have h := map_smul (((PresheafOfModules.pullbackPushforwardAdjunction
    (pullbackPhi f)).unit.app L).app (op W)).hom u g
  refine h.trans ?_
  exact ModuleCat.restrictScalars.smul_def
    (M := ((pullbackPre f).obj L).obj (op ((Opens.map f.base).obj W)))
    (f.c.app (op W)).hom u (pullSec f L W g)

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★**切断の引き戻しは制限と可換である**。

★★これが本ファイルの要である——`η_L` は**前層加群の射**なので、
これはその**自然性の四角形そのもの**であり、ただで手に入る。
★★★同型の水準で同じことを言うと `bcIso` の `W` についての自然性になり、
Beck–Chevalley の三部作をもう一度書く羽目になる。 -/
theorem pullSec_restrict {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    {W W'' : Y.Opens} (hW : W'' ≤ W) (g : (L.obj (op W) : Type))
    (hU : (Opens.map f.base).obj W'' ≤ (Opens.map f.base).obj W) :
    pullSec f L W'' (L.map (homOfLE hW).op g)
      = ((pullbackPre f).obj L).map (homOfLE hU).op (pullSec f L W g) :=
  PresheafOfModules.naturality_apply
    ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app L)
    (homOfLE hW).op g

/-! ## ★★★★★★★★段 1 —— `bcMate` の定義を開く -/

/-- ★★★★★★★★**`res_W η_L` は `η^{on}` と `bcMate` の合成である**。

★これは `bcMate` の**定義**（`res_W η_L` の随伴転置）と
`Adjunction.homEquiv_unit` から直ちに出る。 -/
theorem restrict_unit_eq {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) :
    (restrictPresheafFunctor Y W).map
        ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app L)
      = (adjOn f W).unit.app ((restrictPresheafFunctor Y W).obj L)
        ≫ (PresheafOfModules.pushforward (phiOn f W)).map ((bcMate f W).app L) := by
  have h : (adjOn f W).homEquiv _ _ ((bcMate f W).app L)
      = (restrictPresheafFunctor Y W).map
        ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app L) :=
    Equiv.apply_symm_apply _ _
  rw [← h, Adjunction.homEquiv_unit]
  rfl

/-- ★★★★同上を切断の水準で書いたもの。 -/
theorem pullSec_eq_bcMate {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (g : (L.obj (op W) : Type)) :
    pullSec f L W g
      = ((bcMate f W).app L).app (op (Over.mk (𝟙 ((Opens.map f.base).obj W))))
          (((adjOn f W).unit.app ((restrictPresheafFunctor Y W).obj L)).app
            (op (Over.mk (𝟙 W))) g) :=
  congrArg (fun φ => (PresheafOfModules.Hom.app φ (op (Over.mk (𝟙 W)))) g)
    (restrict_unit_eq f L W)

/-! ## ★★★★★★★段 2 —— `bcIso` を消す -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★**輸送した自明化を切断に当てた形**（`V' = f⁻¹W` の場合）。★**`rfl`** である。 -/
theorem trivEquiv_pullTrivOfBase_apply {X Y : Scheme.{0}} (f : X ⟶ Y)
    (L : Y.PresheafOfModules) (W : Y.Opens)
    (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    (x : (((pullbackPre f).obj L).obj (op ((Opens.map f.base).obj W)) : Type)) :
    trivEquiv ((pullbackPre f).obj L) ((Opens.map f.base).obj W)
        (pullTrivOfBase f L W eW (le_refl ((Opens.map f.base).obj W))) x
      = ((pullbackOnUnitIso f W).hom.app (op (Over.mk (𝟙 ((Opens.map f.base).obj W))))
          (((pullbackPreOn f W).map eW.hom).app
              (op (Over.mk (𝟙 ((Opens.map f.base).obj W))))
            ((bcIso f W L).inv.app (op (Over.mk (𝟙 ((Opens.map f.base).obj W)))) x))) := rfl

set_option backward.isDefEq.respectTransparency false in
/-- ★★★**`bcIso.inv` は `bcMate` を消す**——`bcIso` は `bcMate` の `asIso` だから。 -/
theorem bcIso_inv_bcMate {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (Z : (Over ((Opens.map f.base).obj W))ᵒᵖ)
    (y : ((((pullbackPreOn f W).obj ((restrictPresheafFunctor Y W).obj L))).obj Z : Type)) :
    (bcIso f W L).inv.app Z (((bcMate f W).app L).app Z y) = y :=
  congrArg (fun (φ : (pullbackPreOn f W).obj ((restrictPresheafFunctor Y W).obj L)
    ⟶ (pullbackPreOn f W).obj ((restrictPresheafFunctor Y W).obj L)) => φ.app Z y)
    (bcIso f W L).hom_inv_id

/-! ## ★★★★★★★段 3 —— `η^{on}` の自然性で `e_W` を外へ出す -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★**`η^{on}` の自然性**を `e_W.hom` に当てたもの。 -/
theorem unit_naturality_eW {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    (g : (L.obj (op W) : Type)) :
    ((pullbackPreOn f W).map eW.hom).app (op (Over.mk (𝟙 ((Opens.map f.base).obj W))))
        (((adjOn f W).unit.app ((restrictPresheafFunctor Y W).obj L)).app
          (op (Over.mk (𝟙 W))) g)
      = ((adjOn f W).unit.app (𝟙_ (PresheafModulesOn Y W))).app (op (Over.mk (𝟙 W)))
          (eW.hom.app (op (Over.mk (𝟙 W))) g) :=
  (congrArg (fun φ => (PresheafOfModules.Hom.app φ (op (Over.mk (𝟙 W)))) g)
    ((adjOn f W).unit.naturality eW.hom)).symm

/-! ## ★★★★★★★段 4 —— 単位対象の上では `ρ` そのもの -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★**`η` は単位対象の上で `ρ` である**。

★機構は `§9-740` の `etaOn_unitOne`（`η(1) = 1`）と半線型性:
`η(c • 1) = ρ(c) • 1 = ρ c`。 -/
theorem etaOn_app {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    (c : ((𝟙_ (PresheafModulesOn Y W)).obj (op (Over.mk (𝟙 W))) : Type)) :
    (etaOn f W).app (op (Over.mk (𝟙 W))) c = f.c.app (op W) c := by
  refine (congrArg (fun y => ((etaOn f W).app (op (Over.mk (𝟙 W)))) y)
    (smul_unitOne W c).symm).trans ?_
  rw [map_smul, etaOn_unitOne]
  exact (ModuleCat.restrictScalars.smul_def
    (M := (𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj W))).obj
      (op (Over.mk (𝟙 ((Opens.map f.base).obj W)))))
    ((phiOn f W).app (op (Over.mk (𝟙 W)))).hom c
    (unitOne ((Opens.map f.base).obj W))).trans
    (smul_unitOne ((Opens.map f.base).obj W) _)

set_option backward.isDefEq.respectTransparency false in
/-- ★★`etaOn` を unit と `pullbackOnUnitIso` の合成として切断に当てた形。 -/
theorem etaOn_app_eq_unit {X Y : Scheme.{0}} (f : X ⟶ Y) (W : Y.Opens)
    (c : ((𝟙_ (PresheafModulesOn Y W)).obj (op (Over.mk (𝟙 W))) : Type)) :
    (pullbackOnUnitIso f W).hom.app (op (Over.mk (𝟙 ((Opens.map f.base).obj W))))
        (((adjOn f W).unit.app (𝟙_ (PresheafModulesOn Y W))).app (op (Over.mk (𝟙 W))) c)
      = (etaOn f W).app (op (Over.mk (𝟙 W))) c := by
  have h : etaOn f W = (adjOn f W).unit.app (𝟙_ (PresheafModulesOn Y W))
      ≫ (PresheafOfModules.pushforward (phiOn f W)).map (pullbackOnUnitIso f W).hom :=
    Adjunction.homEquiv_unit _ _ _ (pullbackOnUnitIso f W).hom
  rw [h]
  rfl

/-! ## ★★★★★★★★★★橋 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★**輸送した自明化は、引き戻した切断を、座標を引き戻して読む**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

    `trivEquiv (f^*L) (f⁻¹W) (pullTrivOfBase f L W e_W) (pullSec f L W s)`
      `= ρ (trivEquiv L W e_W s)`

★★これが `Definition 1.1 (ii)` の計量を運ぶ段の**チャート独立性**を支える。 -/
theorem trivEquiv_pullSec {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    (s : (L.obj (op W) : Type)) :
    trivEquiv ((pullbackPre f).obj L) ((Opens.map f.base).obj W)
        (pullTrivOfBase f L W eW (le_refl ((Opens.map f.base).obj W)))
        (pullSec f L W s)
      = (f.c.app (op W)) (trivEquiv L W eW s) := by
  refine (trivEquiv_pullTrivOfBase_apply f L W eW (pullSec f L W s)).trans ?_
  have h1 : (bcIso f W L).inv.app (op (Over.mk (𝟙 ((Opens.map f.base).obj W))))
        (pullSec f L W s)
      = ((adjOn f W).unit.app ((restrictPresheafFunctor Y W).obj L)).app
          (op (Over.mk (𝟙 W))) s := by
    rw [pullSec_eq_bcMate]
    exact bcIso_inv_bcMate f L W _ _
  rw [h1, unit_naturality_eW, etaOn_app_eq_unit, etaOn_app]
  rfl

/-! ## ★★★★★自明化の生成切断 -/

/-- ★★**自明化 `e` の生成切断** `e⁻¹(1) ∈ F(V)`。 -/
noncomputable def trivGen {X : Scheme.{0}} (F : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    (F.obj (op V) : Type) := (trivEquiv F V e).symm 1

/-- ★★★**遷移単元は生成切断の座標である**。★**`rfl`** である。

★★これが効く理由: 遷移単元は `e` については**生成切断だけ**を通して見る。
したがって「同型そのものの一致」ではなく「生成切断の一致」で十分になる。 -/
theorem transUnit_eq_trivEquiv_trivGen {X : Scheme.{0}} (F : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    transUnit F V e e' = trivEquiv F V e' (trivGen F V e) := rfl

/-- ★★★**制限した自明化の生成切断は生成切断の制限である**。 -/
theorem trivGen_restrict {X : Scheme.{0}} (F : X.PresheafOfModules) {V W : X.Opens}
    (hWV : W ≤ V) (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    trivGen F W (trivialOfLe F hWV e) = F.map (homOfLE hWV).op (trivGen F V e) := by
  show (trivEquiv F W (trivialOfLe F hWV e)).symm 1
    = F.map (homOfLE hWV).op ((trivEquiv F V e).symm 1)
  apply (trivEquiv F W (trivialOfLe F hWV e)).injective
  rw [LinearEquiv.apply_symm_apply, trivEquiv_restrict F hWV e ((trivEquiv F V e).symm 1),
    LinearEquiv.apply_symm_apply]
  exact (map_one (X.presheaf.map (homOfLE hWV).op).hom).symm

/-! ## ★★★★★★★★★引き戻した生成切断 -/

/-- ★★★★**引き戻した生成切断**。 -/
noncomputable def pullGen {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W)) :
    (((pullbackPre f).obj L).obj (op ((Opens.map f.base).obj W)) : Type) :=
  pullSec f L W (trivGen L W eW)

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★**生成切断の引き戻しは、輸送した自明化の生成切断である**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★機構は `trivEquiv_pullSec` に `s = trivGen` を入れて `ρ(1) = 1` を使うだけ。 -/
theorem pullGen_eq_trivGen {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W)) :
    trivGen ((pullbackPre f).obj L) ((Opens.map f.base).obj W)
        (pullTrivOfBase f L W eW (le_refl ((Opens.map f.base).obj W)))
      = pullGen f L W eW := by
  apply (trivEquiv ((pullbackPre f).obj L) ((Opens.map f.base).obj W)
    (pullTrivOfBase f L W eW (le_refl ((Opens.map f.base).obj W)))).injective
  show (trivEquiv _ _ _) ((trivEquiv _ _ _).symm 1) = _
  rw [LinearEquiv.apply_symm_apply, pullGen, trivEquiv_pullSec f L W eW (trivGen L W eW)]
  show (1 : (Γ(X, (Opens.map f.base).obj W) : Type))
    = (f.c.app (op W)) ((trivEquiv L W eW) ((trivEquiv L W eW).symm 1))
  rw [LinearEquiv.apply_symm_apply]
  exact (map_one (f.c.app (op W)).hom).symm

set_option backward.isDefEq.respectTransparency false in
/-- ★★**輸送した自明化は `X` 側の縮小と可換である**。★**`rfl`** である（2026-08-28 実測）。

★`W` 側の縮小と違い、こちらは `restrictOnFunctor` の推移律だけなので `rfl` で通る。 -/
theorem pullTrivOfBase_shrinkV {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    {V' V'' : X.Opens} (hV''V' : V'' ≤ V')
    (hV'W : V' ≤ (Opens.map f.base).obj W) (hV''W : V'' ≤ (Opens.map f.base).obj W) :
    trivialOfLe ((pullbackPre f).obj L) hV''V' (pullTrivOfBase f L W eW hV'W)
      = pullTrivOfBase f L W eW hV''W := rfl

/-- ★★★★★**輸送した自明化の生成切断は、引き戻した生成切断の制限である**。 -/
theorem trivGen_pullTrivOfBase {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W) :
    trivGen ((pullbackPre f).obj L) V' (pullTrivOfBase f L W eW hV'W)
      = ((pullbackPre f).obj L).map (homOfLE hV'W).op (pullGen f L W eW) := by
  rw [← pullTrivOfBase_shrinkV f L W eW hV'W (le_refl ((Opens.map f.base).obj W)) hV'W,
    trivGen_restrict, pullGen_eq_trivGen]

/-- ★★★★**引き戻した生成切断は `W` の縮小と可換である**。

★機構は `trivGen_restrict`（`Y` 側）と `pullSec_restrict`（自然性）を継ぐだけ。 -/
theorem pullGen_restrict {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    {W W'' : Y.Opens} (hW : W'' ≤ W)
    (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    (hU : (Opens.map f.base).obj W'' ≤ (Opens.map f.base).obj W) :
    pullGen f L W'' (trivialOfLe L hW eW)
      = ((pullbackPre f).obj L).map (homOfLE hU).op (pullGen f L W eW) := by
  show pullSec f L W'' (trivGen L W'' (trivialOfLe L hW eW)) = _
  rw [trivGen_restrict]
  exact pullSec_restrict f L hW (trivGen L W eW) hU

/-- ★前層加群の制限は合成する（元の水準）。

★★`PresheafOfModules.map_comp` は `restrictScalarsComp'` を経由するので、
`map_id_apply'` と同じく**元の水準の形**を別に用意する。 -/
theorem PresheafOfModules.map_comp_apply' {X : Scheme.{0}} (M : X.PresheafOfModules)
    {A B C : X.Opens} (h1 : B ≤ A) (h2 : C ≤ B) (h3 : C ≤ A)
    (x : (M.obj (op A) : Type)) :
    M.map (homOfLE h2).op (M.map (homOfLE h1).op x) = M.map (homOfLE h3).op x :=
  (congrArg (fun (φ : M.obj (op A) ⟶ (ModuleCat.restrictScalars
      (RingCat.Hom.hom (X.ringCatSheaf.obj.map ((homOfLE h1).op ≫ (homOfLE h2).op)))).obj
        (M.obj (op C))) => φ.hom x)
    (M.map_comp (homOfLE h1).op (homOfLE h2).op)).symm

/-! ## ★★★★★★★★★★到達点 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★**遷移単元は `W` の縮小に依らない**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

    `transUnit (f^*L) V'' (pullTrivOfBase f L W'' (e_W|_{W''})) e`
      `= transUnit (f^*L) V'' (pullTrivOfBase f L W e_W) e`

★★★これが `Definition 1.1 (ii)` のチャート独立性の**残っていた半分**である
（もう半分＝同じ `W` の上での自明化独立性は `§9-740` の
`transUnit_pullTrivOfBase`）。

★★★★逸脱の記録: **同型そのものの等式（`pullTrivOfBase_shrinkW`）は示していない**。
それは `bcIso` の `W` についての自然性を要し、Beck–Chevalley の三部作の再演になる。
本ファイルが示したのは**生成切断の一致**であり、`transUnit` はそれだけで決まる
（`transUnit_eq_trivEquiv_trivGen`）。計量が使うのは `transUnit` だけなので足りる。 -/
theorem transUnit_pullTrivOfBase_shrinkW {X Y : Scheme.{0}} (f : X ⟶ Y)
    (L : Y.PresheafOfModules) {W W'' : Y.Opens} (hW : W'' ≤ W)
    (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    {V'' : X.Opens} (hV''W'' : V'' ≤ (Opens.map f.base).obj W'')
    (hV''W : V'' ≤ (Opens.map f.base).obj W)
    (e : (restrictPresheafFunctor X V'').obj ((pullbackPre f).obj L)
      ≅ 𝟙_ (PresheafModulesOn X V'')) :
    transUnit ((pullbackPre f).obj L) V''
        (pullTrivOfBase f L W'' (trivialOfLe L hW eW) hV''W'') e
      = transUnit ((pullbackPre f).obj L) V'' (pullTrivOfBase f L W eW hV''W) e := by
  have hU : (Opens.map f.base).obj W'' ≤ (Opens.map f.base).obj W := fun z hz => hW hz
  rw [transUnit_eq_trivEquiv_trivGen, transUnit_eq_trivEquiv_trivGen]
  congr 1
  rw [trivGen_pullTrivOfBase, trivGen_pullTrivOfBase, pullGen_restrict f L hW eW hU]
  exact PresheafOfModules.map_comp_apply' _ hU hV''W'' hV''W (pullGen f L W eW)

/-! ### ★出典の紐付け(`.src`) -/

def pullSec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(切断の引き戻し＝随伴の unit の成分)",
    sectionId := "genell-def-1-1-ii" }

def pullSec_restrict.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(切断の引き戻しが制限と可換であること)",
    sectionId := "genell-def-1-1-ii" }

def trivEquiv_pullSec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(輸送した自明化は引き戻した切断を座標の引き戻しで読むこと)",
    sectionId := "genell-def-1-1-ii" }

def pullGen_eq_trivGen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(生成切断の引き戻しは輸送した自明化の生成切断であること)",
    sectionId := "genell-def-1-1-ii" }

def transUnit_pullTrivOfBase_shrinkW.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(遷移単元が W の縮小に依らないこと)",
    sectionId := "genell-def-1-1-ii" }

def trivEquiv_pullSec.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "bcMate(Beck–Chevalley の mate は res_W (unit) の随伴転置)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.bcMate") 3,
    .citation "[ABC3]" "etaOn_unitOne(η は生成元を生成元へ送る、§9-740)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.etaOn_unitOne") 3,
    .citation "[mathlib]" "Adjunction.homEquiv_unit"
      (.inMathlib "CategoryTheory.Adjunction.homEquiv_unit") 3,
    .implicitStep
      ("★★★★手口: **同型の水準の自然性を切断の水準へ落とす**。" ++
       "bcIso の W についての自然性(Beck–Chevalley の三部作の再演)ではなく、" ++
       "随伴の unit を切断に当てた pullSec の自然性(前層加群の射なのでただ)を使う") 3 ]

def transUnit_pullTrivOfBase_shrinkW.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "trivEquiv_pullSec(座標の引き戻し)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivEquiv_pullSec") 3,
    .citation "[ABC3]" "transUnit_pullTrivOfBase(同じ W の上でのチャート独立性、§9-740)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.transUnit_pullTrivOfBase") 3,
    .implicitStep
      ("★★★★逸脱: 同型そのものの等式(pullTrivOfBase_shrinkW)は**示していない**——" ++
       "それは bcIso の W についての自然性を要し、Beck–Chevalley の三部作の再演になる。" ++
       "示したのは**生成切断の一致**であり、transUnit はそれだけで決まる" ++
       "(transUnit_eq_trivEquiv_trivGen)。計量が使うのは transUnit だけなので足りる") 3,
    .implicitStep
      ("★★残っている段: PullChart と pullChartH の 3 段" ++
       "(triv_indep / shrinkV / shrinkW)、そして LocalMetric.pullback と AMetric の引き戻し") 3 ]

end ABC3.Found.Arakelov
