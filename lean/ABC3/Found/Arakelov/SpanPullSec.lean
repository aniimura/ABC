/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PullbackPic
import ABC3.Meta.Claim

/-!
# `Γ(f^*L)` は**大域切断の引き戻しで生成される**（単位束と同型移送、`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★★★★これは何か——**span で書くと `T` と `Γ(X,⊤)` の摩擦が消える**

底変換の有限側は結局

    `Γ(f^*L) = Γ(X,⊤)-span { pullSec f L ⊤ s : s ∈ Γ(L) }`

に帰着する。★★これを**係数拡大の写像 `μ` の全射性**ではなく
**span の等式**として書くと、`T` と `Γ(X,⊤)` の同一視（`§9-788` の摩擦）が要らない。

## ★★★★★機構（3 段）

    span_pullSecT_unit   : ★単位束では成り立つ（`pullbackUnitPreIso_pullSec_one`、`§9-745`）
    pullSecT_naturality  : ★★`pullSec` は前層の射に沿って自然（随伴の unit の自然性）
    span_pullSecT_congr  : ★★★同型で移せる

★★★★残る段は**テンソル**である:
`A ⊗ A⁻¹ ≅ Ō` と在庫の `delta_pullSec`（`§9-743`）で
`span { pullSec a ⊗ pullSec b } = ⊤` を出し、
`§9-785` の `surjective_of_map_surjective` で `span { pullSec a } = ⊤` に落とす。

## ★測定の記録

単位束の場合は **`f` の全射性を要さない**。
★`f^*Ō ≅ Ō`（`pullbackUnitPreIso`、任意の `f`）で `Γ(f^*Ō) ≅ Γ(X,⊤)` となり、
引き戻した `1` が `1` に行くので span は全体になる。
★★したがってテンソルで回す道でも `f` の全射性は要らない。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll

/-- ★**大域切断の引き戻し**（戻り値の型を `op ⊤` で宣言し直したもの）。 -/
noncomputable def pullSecT {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (s : (L.obj (op ⊤) : Type)) :
    (((pullbackPre f).obj L).obj (op (⊤ : X.Opens)) : Type) :=
  pullSec f L ⊤ s

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**単位束では `Γ(f^*Ō)` は引き戻した大域切断で生成される**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★`f^*Ō ≅ Ō`（`pullbackUnitPreIso`）で `Γ(f^*Ō) ≅ Γ(X,⊤)` となり、
在庫の `pullbackUnitPreIso_pullSec_one`（`§9-745`）が
引き戻した `1` を `1` に送る。★★`1` を含む部分加群は全体である。 -/
theorem span_pullSecT_unit {X Y : Scheme.{0}} (f : X ⟶ Y) :
    Submodule.span (Γ(X, (⊤ : X.Opens)) : Type)
      (Set.range (pullSecT f (𝟙_ Y.PresheafOfModules))) = ⊤ := by
  set e : (((pullbackPre f).obj (𝟙_ Y.PresheafOfModules)).obj (op (⊤ : X.Opens)) : Type)
      ≃ₗ[(Γ(X, (⊤ : X.Opens)) : Type)] (Γ(X, (⊤ : X.Opens)) : Type) :=
    ((PresheafOfModules.evaluation X.ringCatSheaf.obj (op (⊤ : X.Opens))).mapIso
      (pullbackUnitPreIso f)).toLinearEquiv with he
  have hone : e (pullSecT f (𝟙_ Y.PresheafOfModules)
      (1 : (Γ(Y, (⊤ : Y.Opens)) : Type))) = (1 : (Γ(X, (⊤ : X.Opens)) : Type)) :=
    pullbackUnitPreIso_pullSec_one f
  rw [eq_top_iff]
  intro x _
  have h1 : x = (e x) • (pullSecT f (𝟙_ Y.PresheafOfModules)
      (1 : (Γ(Y, (⊤ : Y.Opens)) : Type))) := by
    apply e.injective
    rw [map_smul, hone, smul_eq_mul, mul_one]
  rw [h1]
  exact Submodule.smul_mem _ _ (Submodule.subset_span ⟨_, rfl⟩)

set_option maxHeartbeats 1000000 in
/-- ★★★**切断の引き戻しは前層の射に沿って自然である**（随伴の unit の自然性）。 -/
theorem pullSecT_naturality {X Y : Scheme.{0}} (f : X ⟶ Y) {A B : Y.PresheafOfModules}
    (φ : A ⟶ B) (s : (A.obj (op ⊤) : Type)) :
    pullSecT f B ((φ.app (op ⊤)) s)
      = (((pullbackPre f).map φ).app (op (⊤ : X.Opens))) (pullSecT f A s) := by
  have hnat := (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.naturality φ
  have h := congrArg (fun (m : _ ⟶ _) => (m.app (op (⊤ : Y.Opens))).hom s) hnat
  show (((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app B).app
      (op (⊤ : Y.Opens))).hom ((φ.app (op ⊤)) s) = _
  exact h

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**生成されるという性質は前層の同型で移る**。 -/
theorem span_pullSecT_congr {X Y : Scheme.{0}} (f : X ⟶ Y) {A B : Y.PresheafOfModules}
    (φ : A ≅ B)
    (h : Submodule.span (Γ(X, (⊤ : X.Opens)) : Type) (Set.range (pullSecT f A)) = ⊤) :
    Submodule.span (Γ(X, (⊤ : X.Opens)) : Type) (Set.range (pullSecT f B)) = ⊤ := by
  set ψ : (((pullbackPre f).obj A).obj (op (⊤ : X.Opens)) : Type)
      ≃ₗ[(Γ(X, (⊤ : X.Opens)) : Type)] (((pullbackPre f).obj B).obj (op (⊤ : X.Opens)) : Type) :=
    ((PresheafOfModules.evaluation X.ringCatSheaf.obj (op (⊤ : X.Opens))).mapIso
      ((pullbackPre f).mapIso φ)).toLinearEquiv with hψ
  set χ : (A.obj (op (⊤ : Y.Opens)) : Type) ≃ₗ[(Γ(Y, (⊤ : Y.Opens)) : Type)]
      (B.obj (op (⊤ : Y.Opens)) : Type) :=
    ((PresheafOfModules.evaluation Y.ringCatSheaf.obj (op (⊤ : Y.Opens))).mapIso φ).toLinearEquiv
    with hχ
  have hcomm : ∀ s, pullSecT f B (χ s) = ψ (pullSecT f A s) := by
    intro s
    exact pullSecT_naturality f φ.hom s
  have hrange : Set.range (pullSecT f B)
      = ((ψ : _ →ₗ[(Γ(X, (⊤ : X.Opens)) : Type)] _) : _ → _) '' (Set.range (pullSecT f A)) := by
    ext y
    constructor
    · rintro ⟨t, rfl⟩
      refine ⟨pullSecT f A (χ.symm t), ⟨_, rfl⟩, ?_⟩
      have hh := hcomm (χ.symm t)
      rw [χ.apply_symm_apply] at hh
      exact hh.symm
    · rintro ⟨z, ⟨s, rfl⟩, rfl⟩
      exact ⟨χ s, hcomm s⟩
  rw [hrange, ← Submodule.map_span, h, Submodule.map_top, LinearMap.range_eq_top]
  exact ψ.surjective

/-! ### ★出典の紐付け(`.src`) -/

def span_pullSecT_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(単位束では Γ(f^*Ō) は引き戻した大域切断で生成される)",
    sectionId := "genell-def-1-1-ii" }

def span_pullSecT_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(生成されるという性質は前層の同型で移る)",
    sectionId := "genell-def-1-1-ii" }

def span_pullSecT_congr.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "pullbackUnitPreIso_pullSec_one(f^*Ō ≅ Ō は 1 を 1 に送る、§9-745)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullbackUnitPreIso_pullSec_one") 3,
    .citation "[mathlib]" "Adjunction.unit の自然性"
      (.inMathlib "CategoryTheory.NatTrans.naturality") 3,
    .implicitStep
      ("★★残る段はテンソルである: A ⊗ A⁻¹ ≅ Ō と在庫の delta_pullSec(§9-743)で " ++
       "span { pullSec a ⊗ pullSec b } = ⊤ を出し、" ++
       "§9-785 の surjective_of_map_surjective で span { pullSec a } = ⊤ に落とす") 4,
    .implicitStep
      ("★測定: 単位束の場合は f の全射性を要さない。" ++
       "f^*Ō ≅ Ō(任意の f)で Γ(f^*Ō) ≅ Γ(X,⊤) となり、引き戻した 1 が 1 に行くからである") 3 ]

end ABC3.Found.Arakelov
