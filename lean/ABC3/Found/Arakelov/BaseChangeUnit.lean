/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.BaseChangeSetup
import ABC3.Found.Arakelov.PullbackPic
import ABC3.Found.Arakelov.DegArithPre
import ABC3.Meta.Claim

/-!
# 底変換 —— **単位束の場合** `𝓞_K ⊗_{𝓞_F} Γ(Ō) → Γ(f^*Ō)` は全射（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★これは何か

`§9-786` で係数拡大の写像

    `μ_L : T ⊗_S Γ(L) → Γ(f^*L)`   （`muBC`）

を作った。本ファイルは**単位束 `L = Ō` の場合に `μ` が全射である**ことを示す。

★★これが `§9-785` の道筋の**出発点**である:

    `μ_A ⊗ μ_{A⁻¹} ≅ μ_{Ō}`（全射）  ⟹  `μ_A` 全射  ⟹（可逆加群なので）`μ_A` 同型

## ★★★★★機構

    `ψ : Γ(f^*Ō) → T`,  `y ↦ Γ-Spec( (f^*Ō ≅ Ō)(y) )`

は**単射**（同型 2 つの合成）であり、在庫の `pullbackUnitPreIso_pullSec_one`（`§9-745`）が

    `ψ (pullSec f Ō ⊤ 1) = 1`

を与える。★したがって `ψ (μ (t ⊗ 1)) = t` であり、任意の `y` について
`t := ψ y` と取れば `μ (t ⊗ 1) = y`（`ψ` の単射性）。

## ★配管の記録

`pullSec f L ⊤` の行き先は `op ((Opens.map f.base).obj ⊤)` の成分で、
`op ⊤` とは `rfl` だが**インスタンス探索が通らない**。
★`pullSecTop`（`§9-786`）と `psiU` の両方で**戻り値・引数の型を `op ⊤` で宣言し直す**
ことでこの摩擦を避けている。★★`rw` が通らない箇所は
「ゴールの綴りで `have` を作る」（`lean-idioms.md`）で回避した。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll
open scoped TensorProduct

/-! ## ★★★★★★係数拡大の写像 `μ_L` -/

/-- ★★★★★★**`μ_L : T ⊗_S Γ(L) → Γ(f^*L)`**（切断の引き戻しの係数拡大）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★これが同型であることが底変換の**有限側**の内容である。 -/
noncomputable def muBC (S T : Type) [CommRing S] [CommRing T] [Algebra S T]
    (L : (Spec (CommRingCat.of S)).PresheafOfModules) :
    T ⊗[S] (gammaModPre (CommRingCat.of S) L : Type)
      →ₗ[T] (gammaModPre (CommRingCat.of T)
        ((pullbackPre (Spec.map (CommRingCat.ofHom (algebraMap S T)))).obj L) : Type) :=
  bcLift S T _ _ (pullSecTop (CommRingCat.ofHom (algebraMap S T)) L)
    (pullSecTop_add (CommRingCat.ofHom (algebraMap S T)) L)
    (fun a x => pullSecTop_smul (CommRingCat.ofHom (algebraMap S T)) L a x)

/-! ## ★★★★★★★単位束の場合 -/

/-- ★**`Γ(f^*Ō) → T`**——`f^*Ō ≅ Ō` と `Γ-Spec` 同型の合成。 -/
noncomputable def psiU (S T : Type) [CommRing S] [CommRing T] [Algebra S T]
    (y : (((pullbackPre (Spec.map (CommRingCat.ofHom (algebraMap S T)))).obj
        (𝟙_ (Spec (CommRingCat.of S)).PresheafOfModules)).obj
        (op (⊤ : (Spec (CommRingCat.of T)).Opens)) : Type)) : T :=
  (Scheme.ΓSpecIso (CommRingCat.of T)).hom.hom
    (((pullbackUnitPreIso (Spec.map (CommRingCat.ofHom (algebraMap S T)))).hom.app (op ⊤)) y)

/-- ★`ψ` は係数と両立する。 -/
theorem psiU_smul (S T : Type) [CommRing S] [CommRing T] [Algebra S T]
    (c : (Γ(Spec (CommRingCat.of T), (⊤ : (Spec (CommRingCat.of T)).Opens)) : Type))
    (y : (((pullbackPre (Spec.map (CommRingCat.ofHom (algebraMap S T)))).obj
        (𝟙_ (Spec (CommRingCat.of S)).PresheafOfModules)).obj
        (op (⊤ : (Spec (CommRingCat.of T)).Opens)) : Type)) :
    psiU S T (c • y) = (Scheme.ΓSpecIso (CommRingCat.of T)).hom.hom c * psiU S T y := by
  show (Scheme.ΓSpecIso (CommRingCat.of T)).hom.hom
      (((pullbackUnitPreIso (Spec.map (CommRingCat.ofHom (algebraMap S T)))).hom.app (op ⊤))
        (c • y)) = _
  rw [map_smul]
  exact map_mul (Scheme.ΓSpecIso (CommRingCat.of T)).hom.hom _ _

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★**`ψ` は引き戻した `1` を `1` に送る**（在庫の `pullbackUnitPreIso_pullSec_one`）。 -/
theorem psiU_pullSec_one (S T : Type) [CommRing S] [CommRing T] [Algebra S T] :
    psiU S T (pullSec (Spec.map (CommRingCat.ofHom (algebraMap S T)))
      (𝟙_ (Spec (CommRingCat.of S)).PresheafOfModules) ⊤
      (1 : (Γ(Spec (CommRingCat.of S), (⊤ : (Spec (CommRingCat.of S)).Opens)) : Type))) = 1 := by
  have h : ((pullbackUnitPreIso (Spec.map (CommRingCat.ofHom (algebraMap S T)))).hom.app
      (op (⊤ : (Spec (CommRingCat.of T)).Opens)))
      (pullSec (Spec.map (CommRingCat.ofHom (algebraMap S T)))
        (𝟙_ (Spec (CommRingCat.of S)).PresheafOfModules) ⊤
        (1 : (Γ(Spec (CommRingCat.of S), (⊤ : (Spec (CommRingCat.of S)).Opens)) : Type)))
      = (1 : (Γ(Spec (CommRingCat.of T), (⊤ : (Spec (CommRingCat.of T)).Opens)) : Type)) :=
    pullbackUnitPreIso_pullSec_one (Spec.map (CommRingCat.ofHom (algebraMap S T)))
  unfold psiU
  rw [h]
  exact map_one _

set_option maxHeartbeats 1000000 in
theorem psiU_pullSecTop_one (S T : Type) [CommRing S] [CommRing T] [Algebra S T] :
    psiU S T (pullSecTop (CommRingCat.ofHom (algebraMap S T))
      (𝟙_ (Spec (CommRingCat.of S)).PresheafOfModules)
      (1 : (Γ(Spec (CommRingCat.of S), (⊤ : (Spec (CommRingCat.of S)).Opens)) : Type))) = 1 :=
  psiU_pullSec_one S T

set_option maxHeartbeats 1000000 in
/-- ★★★`ψ` は単射である（同型 2 つの合成）。 -/
theorem psiU_injective (S T : Type) [CommRing S] [CommRing T] [Algebra S T] :
    Function.Injective (psiU S T) := by
  have h1 : Function.Injective (fun y =>
      ((pullbackUnitPreIso (Spec.map (CommRingCat.ofHom (algebraMap S T)))).hom.app
        (op (⊤ : (Spec (CommRingCat.of T)).Opens))) y) :=
    ((PresheafOfModules.evaluation (Spec (CommRingCat.of T)).ringCatSheaf.obj
      (op (⊤ : (Spec (CommRingCat.of T)).Opens))).mapIso
      (pullbackUnitPreIso (Spec.map (CommRingCat.ofHom (algebraMap S T))))).toLinearEquiv.injective
  have h2 : Function.Injective ((Scheme.ΓSpecIso (CommRingCat.of T)).hom.hom) :=
    (gammaSpecRingEquiv (CommRingCat.of T)).injective
  exact h2.comp h1

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★**単位束では `μ` は全射である**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★これが `§9-785` の道筋の出発点である——
`μ_A ⊗ μ_{A⁻¹} ≅ μ_{Ō}` が全射なので `μ_A` が全射、ゆえに同型。 -/
theorem muBC_unit_surjective (S T : Type) [CommRing S] [CommRing T] [Algebra S T] :
    Function.Surjective (muBC S T (𝟙_ (Spec (CommRingCat.of S)).PresheafOfModules)) := by
  intro y
  refine ⟨(psiU S T y) ⊗ₜ[S]
    (1 : (Γ(Spec (CommRingCat.of S), (⊤ : (Spec (CommRingCat.of S)).Opens)) : Type)), ?_⟩
  apply psiU_injective S T
  show psiU S T (((Scheme.ΓSpecIso (CommRingCat.of T)).inv.hom (psiU S T y)) •
      (pullSecTop (CommRingCat.ofHom (algebraMap S T))
        (𝟙_ (Spec (CommRingCat.of S)).PresheafOfModules)
        (1 : (Γ(Spec (CommRingCat.of S), (⊤ : (Spec (CommRingCat.of S)).Opens)) : Type)))) = _
  rw [psiU_smul, psiU_pullSecTop_one, mul_one]
  exact congrArg (fun (m : _ ⟶ _) => CommRingCat.Hom.hom m (psiU S T y))
    (Scheme.ΓSpecIso (CommRingCat.of T)).inv_hom_id

/-! ### ★出典の紐付け(`.src`) -/

def muBC.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(係数拡大の写像 μ_L : T ⊗_S Γ(L) → Γ(f^*L))",
    sectionId := "genell-def-1-1-ii" }

def muBC_unit_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(単位束では μ は全射である)",
    sectionId := "genell-def-1-1-ii" }

def muBC_unit_surjective.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "pullbackUnitPreIso_pullSec_one(f^*Ō ≅ Ō は引き戻した 1 を 1 に送る、§9-745)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullbackUnitPreIso_pullSec_one") 3,
    .citation "[ABC3]" "pullSecTop_smul(切断の引き戻しの半線形性、§9-786)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullSecTop_smul") 4,
    .implicitStep
      ("★★次の段: μ_{A ⊗ B} ≅ μ_A ⊗ μ_B(在庫の delta_pullSec、§9-743)を使い、" ++
       "A ⊗ A⁻¹ ≅ Ō で本定理に帰着させて μ_A の全射性を出す。" ++
       "★★★そのあと Module.Invertible.bijective_of_surjective で同型に上がる") 4 ]

end ABC3.Found.Arakelov
